! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copyright (c) 2015, Regents of the University of Colorado
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification, are 
! permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this list of 
!    conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice, this list
!    of conditions and the following disclaimer in the documentation and/or other 
!    materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be 
!    used to endorse or promote products derived from this software without specific prior
!    written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
! EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
! MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
! THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT 
! OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! History:
! Nov 2018- T. Michibata - Added CloudSat+MODIS Warmrain Diagnostics
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE MOD_COSP_DIAGNOSTICS
  USE COSP_KINDS, ONLY: wp
  USE COSP_PHYS_CONSTANTS,  ONLY: tmelt
  USE MOD_COSP_STATS,  ONLY: hist2D
  USE MOD_COSP_CONFIG, &
      ONLY: R_UNDEF,          R_GROUND,         &
            SGCLD_CLR,        SGCLD_ST,         &
            SGCLD_CUM,                          &
            CWP_THRESHOLD,    COT_THRESHOLD,    &
            CFODD_NDBZE,      CFODD_NICOD,      &
            CFODD_BNDRE,      CFODD_BNDZE,      &
            CFODD_NCLASS,                       &
            CFODD_DBZE_MIN,   CFODD_DBZE_MAX,   &
            CFODD_ICOD_MIN,   CFODD_ICOD_MAX,   &
            CFODD_DBZE_WIDTH, CFODD_ICOD_WIDTH, &
            CFODD_HISTDBZE,   CFODD_HISTICOD,   &
            WR_NREGIME
  IMPLICIT NONE
CONTAINS

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !-----------------------  SUBROUTINE COSP_DIAG_WARMRAIN  -----------------------
  ! (c) 2018, Research Institute for Applied Mechanics (RIAM), Kyushu Univ.
  ! All rights reserved.
  ! * Purpose:    1) Diagnose Contoured Frequency by Optical Depth Diagram (CFODD)
  !                  from CloudSat Radar and MODIS retrievals.
  !               2) Diagnose Warm-Rain Occurrence Frequency (nonprecip/drizzle/rain)
  !                  from CloudSat Radar.
  ! * History:    2018.01.25 (T.Michibata): Test run
  !               2018.05.17 (T.Michibata): update for COSP2
  !               2018.09.19 (T.Michibata): modified I/O
  !               2018.11.22 (T.Michibata): minor revisions
  ! * References: Suzuki et al. (JAS'10, JAS'11, JGR'13, GRL'13, JAS'15)
  !               Jing et al. (JGR'17); Jing and Suzuki (GRL'18)
  !               Kay et al. (JGR'18)
  !               Michibata et al. (ACP'14)
  ! * Contact:    Takuro Michibata (RIAM, Kyushu University, Japan).
  !               E-mail: michibata@riam.kyushu-u.ac.jp
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_DIAG_WARMRAIN( Npoints, Ncolumns, Nlevels,           & !! in
                                 temp,    zlev,     delz,              & !! in
                                 lwp,     liqcot,   liqreff, liqcfrc,  & !! in
                                 iwp,     icecot,   icereff, icecfrc,  & !! in
                                 fracout, dbze,                        & !! in
                                 cfodd_ntotal,                         & !! inout
                                 wr_occfreq_ntotal                     ) !! inout

    ! Inputs
    integer, intent(in) :: &
         Npoints,          & ! Number of horizontal gridpoints
         Ncolumns,         & ! Number of subcolumns
         Nlevels             ! Number of vertical layers
    real(wp), dimension(Npoints), intent(in) :: &
         lwp,              & ! MODIS LWP [kg m-2]
         liqcot,           & ! MODIS liq. COT
         liqreff,          & ! MODIS liq. Reff [m]
         liqcfrc             ! MODIS liq. cloud fraction
    real(wp), dimension(Npoints), intent(in) :: &
         iwp,              & ! MODIS IWP [kg m-2]
         icecot,           & ! MODIS ice COT
         icereff,          & ! MODIS ice Reff [m]
         icecfrc             ! MODIS ice cloud fraction
    real(wp), dimension(Npoints,Nlevels), intent(in) :: &
         zlev,             & ! altitude [m] for model level
         delz                ! delta Z [m] (= zlevm(k+1)-zlemv(k))
    real(wp), dimension(Npoints,1,Nlevels),intent(in) :: &
         temp                ! temperature [K]
    real(wp), dimension(Npoints,Ncolumns,Nlevels),intent(in) :: &
         fracout,          & ! SCOPS subcolumn retrieval
         dbze                ! Radar reflectivity [dBZe]

    ! Outputs
    real(wp),dimension(Npoints,CFODD_NDBZE,CFODD_NICOD,CFODD_NCLASS),intent(inout) :: &
         cfodd_ntotal        ! # of CFODD samples for each ICOD/dBZe bin
    real(wp),dimension(Npoints,WR_NREGIME),intent(inout) :: &
         wr_occfreq_ntotal   ! # of SLWC samples
    
    ! Local variables
    integer  :: i, j, k
    integer  :: ix, iy
    integer  :: kctop, kcbtm
    integer  :: icls
    integer  :: iregime
    real     :: cmxdbz
    real(wp) :: diagcgt   !! diagnosed cloud geometric thickness [m]
    real(wp) :: diagdbze  !! diagnosed dBZe
    real(wp) :: diagicod  !! diagnosed in-cloud optical depth
    real(wp) :: cbtmh     !! diagnosed in-cloud optical depth
    real(wp), dimension(Npoints,Ncolumns,Nlevels) :: icod  !! in-cloud optical depth (ICOD)
    logical  :: octop, ocbtm, oslwc


    !! initialize
    cfodd_ntotal(:,:,:,:)  = 0._wp
    wr_occfreq_ntotal(:,:) = 0._wp
    icod(:,:,:) = 0._wp

    do i = 1, Npoints
       !! check by MODIS retrieval
       if ( lwp(i)     .le.  CWP_THRESHOLD  .or.  &
          & liqcot(i)  .le.  COT_THRESHOLD  .or.  &
          & liqreff(i) .lt.  CFODD_BNDRE(1) .or.  &
          & liqreff(i) .gt.  CFODD_BNDRE(4) .or.  &
          & iwp(i)     .gt.  CWP_THRESHOLD  .or.  &       !! exclude ice cloud
          & icecot(i)  .gt.  COT_THRESHOLD  .or.  &       !! exclude ice cloud
          & icereff(i) .gt.  CFODD_BNDRE(1)       ) then  !! exclude ice cloud
          cycle
       else
          !! retrieve the CFODD array based on Reff
          icls = 0
          if (    liqreff(i) .ge. CFODD_BNDRE(1) .and. liqreff(i) .lt. CFODD_BNDRE(2) ) then
             icls = 1
          elseif( liqreff(i) .ge. CFODD_BNDRE(2) .and. liqreff(i) .lt. CFODD_BNDRE(3) ) then
             icls = 2
          elseif( liqreff(i) .ge. CFODD_BNDRE(3) .and. liqreff(i) .le. CFODD_BNDRE(4) ) then
             icls = 3
          endif
       endif

       !CDIR NOLOOPCHG
       do j = 1, Ncolumns, 1
          octop = .true.  !! initialize
          ocbtm = .true.  !! initialize
          kcbtm =     0   !! initialize
          kctop =     0   !! initialize
          !CDIR NOLOOPCHG
          do k = Nlevels, 1, -1  !! scan from cloud-bottom to cloud-top
             if ( dbze(i,j,k) .eq. R_GROUND .or. &
                  dbze(i,j,k) .eq. R_UNDEF       ) cycle
             if ( ocbtm                          .and. &
                & fracout(i,j,k) .ne. SGCLD_CLR  .and. &
                & dbze(i,j,k)    .ge. CFODD_DBZE_MIN   ) then
                ocbtm = .false.  !! cloud bottom detected
                kcbtm = k
                kctop = k
             endif
             if (       octop                    .and. &  !! scan cloud-top
                & .not. ocbtm                    .and. &  !! cloud-bottom already detected
                & fracout(i,j,k) .ne. SGCLD_CLR  .and. &  !! exclude clear sky
                & dbze(i,j,k)    .ge. CFODD_DBZE_MIN   ) then
                kctop = k  !! update
             endif
          enddo  !! k loop
          if( ocbtm ) cycle  !! cloud wasn't detected in this subcolumn

          !! check SLWC?
          if( temp(i,1,kctop) .lt. tmelt ) cycle  !! return to the j (subcolumn) loop
          oslwc = .true.
          cmxdbz = CFODD_DBZE_MIN  !! initialized

          !CDIR NOLOOPCHG
          do k = kcbtm, kctop, -1
             cmxdbz = max( cmxdbz, dbze(i,j,k) )  !! column maximum dBZe update
             if ( fracout(i,j,k) .eq. SGCLD_CLR  .or.  &
                & fracout(i,j,k) .eq. SGCLD_CUM  .or.  &
                & dbze   (i,j,k) .lt. CFODD_DBZE_MIN   ) then
                oslwc = .false.
             endif
          enddo
          if ( .not. oslwc ) cycle  !! return to the j (subcolumn) loop

          !! warm-rain occurrence frequency
          iregime = 0
          if( cmxdbz .lt. CFODD_BNDZE(1) ) then
             iregime = 1  !! non-precipitating
          elseif( cmxdbz .ge. CFODD_BNDZE(1) .and. cmxdbz .lt. CFODD_BNDZE(2) ) then
             iregime = 2  !! drizzling
          elseif( cmxdbz .ge. CFODD_BNDZE(2) ) then
             iregime = 3  !! raining
          endif
          wr_occfreq_ntotal(i,iregime) = wr_occfreq_ntotal(i,iregime) + 1._wp

          !! retrievals of ICOD and dBZe bin for CFODD plane
          diagcgt = zlev(i,kctop) - zlev(i,kcbtm)
          cbtmh   = zlev(i,kcbtm)
          !CDIR NOLOOPCHG
          do k = kcbtm, kctop, -1
             if( k .eq. kcbtm ) then
                diagicod = liqcot(i)
             else
                diagicod = liqcot(i) * ( 1._wp - ( (zlev(i,k)-cbtmh)/diagcgt)**(5._wp/3._wp) )
             endif
             icod(i,j,k) = min( diagicod, CFODD_ICOD_MAX )
          enddo

       enddo  ! j (Ncolumns)

       !! # of samples for CFODD (joint 2d-histogram dBZe vs ICOD)
       call hist2d( dbze(i,1:Ncolumns,1:Nlevels), icod(i,1:Ncolumns,1:Nlevels),  &
                  & Ncolumns*Nlevels,                                            &
                  & CFODD_HISTDBZE, CFODD_NDBZE, CFODD_HISTICOD, CFODD_NICOD,    &
                  & cfodd_ntotal( i, 1:CFODD_NDBZE, 1:CFODD_NICOD, icls )        )

    enddo     ! i (Npoints)

  RETURN
  END SUBROUTINE COSP_DIAG_WARMRAIN


END MODULE MOD_COSP_DIAGNOSTICS
