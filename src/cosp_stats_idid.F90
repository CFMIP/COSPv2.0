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
! May 2018 - T. Michibata - IDiD Initial version
! Sep 2018 - T. Michibata - Modified IDiD output to save memory
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE MOD_COSP_STATS_IDID
  USE COSP_KINDS, ONLY: wp
  USE MOD_COSP_CONFIG, &
      ONLY: R_UNDEF,     R_GROUND,                            &
            CFODD_NBINX, CFODD_NBINY, Nclass,                 &
            CFODD_BNDRE, CFODD_NREFF,                         &
            CFODD_XMIN,  CFODD_XMAX, CFODD_YMIN, CFODD_YMAX,  &
            CFODD_DELX,  CFODD_DELY,                          &
            PDFMAP_NPHASE
  IMPLICIT NONE
CONTAINS

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !--------------------------  SUBROUTINE COSP_IDID_WR  --------------------------
  ! (c) 2018, Research Institute for Applied Mechanics (RIAM), Kyushu Univ.
  ! All rights reserved.
  ! * Purpose:    1) Diagnose Contoured Frequency by Optical Depth Diagram (CFODD)
  !                  from CloudSat Radar and MODIS retrievals.
  !               2) Diagnose PDFs of Warm-Rain Phase (nonprecip/drizzle/rain)
  !                  from CloudSat Radar.
  ! * History:    2018.01.25 (T.Michibata): Test run
  !               2018.05.17 (T.Michibata): update for COSP2
  !               2018.09.19 (T.Michibata): modified I/O
  ! * References: Suzuki et al. (JAS'10, JAS'11, JGR'13, GRL'13, JAS'15)
  !               Jing et al. (JGR'17, GRL'18)
  !               Kay et al. (JGR'18)
  !               Michibata et al. (ACP'14)
  ! * Contact:    Takuro Michibata (RIAM, Kyushu University, Japan).
  !               E-mail: michibata@riam.kyushu-u.ac.jp
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_IDID_WR( Npoints, Ncolumns, Nlevels,           & !! in
                           temp,    zlev,     delz,              & !! in
                           lwp,     liqcot,   liqreff, liqcfrc,  & !! in
                           iwp,     icecot,   icereff, icecfrc,  & !! in
                           fracout, dbze,                        & !! in
                           ncfodd,                               & !! inout
                           cfodd,                                & !! inout
                           nslwc,                                & !! inout
                           pdf                                   ) !! inout

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
    real(wp),dimension(Npoints,1,Nlevels),intent(in) :: &
         temp                ! temperature [K]
    real(wp),dimension(Npoints,Ncolumns,Nlevels),intent(in) :: &
         fracout,          & ! SCOPS subcolumn retrieval
         dbze                ! Radar reflectivity [dBZe]

    ! Outputs
    real(wp),dimension(Npoints,CFODD_NBINX,CFODD_NBINY,Nclass),intent(inout) :: &
         ncfodd,           & ! # of CFODD samples for each ICOD/dBZe bin
         cfodd               ! CFODD histgram
    real(wp),dimension(Npoints,PDFMAP_NPHASE),intent(inout) :: &
         nslwc,            & ! # of SLWC samples
         pdf                 ! pdf of nonprecip/drizzle/precip occurrences
    
    ! Local variables
    integer  :: i, j, k
    integer  :: ix, iy
    integer  :: kctop, kcbtm
    integer  :: icls, isum
    integer  :: iphase
    real     :: cmxdbz
    real(wp) :: diagcgt   !! diagnosed cloud geometric thickness [m]
    real(wp) :: diagdbze  !! diagnosed dBZe
    real(wp) :: diagicod  !! diagnosed in-cloud optical depth
    real(wp) :: cbtmh     !! diagnosed in-cloud optical depth
    logical  :: octop, ocbtm, oslwc


    isum = CFODD_NREFF  !! for gathering all Reff samples
    !! initialize
    cfodd (:,:,:,:) = 0._wp
    ncfodd(:,:,:,:) = 0._wp
    pdf   (:,:)     = 0._wp
    nslwc (:,:)     = 0._wp

    do i = 1, Npoints
       !! check by MODIS retrieval
       if ( lwp(i)     .le.   0.0           .or.  &
          & liqcot(i)  .le.  0.30           .or.  &
          & liqreff(i) .lt.  CFODD_BNDRE(1) .or.  &
          & liqreff(i) .gt.  CFODD_BNDRE(4) .or.  &
          & iwp(i)     .gt.   1.D-3         .or.  &       !! exclude ice cloud
          & icecot(i)  .gt.  0.30           .or.  &       !! exclude ice cloud
          & icereff(i) .gt.   5.D-6               ) then  !! exclude ice cloud
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
             if ( ocbtm                  .and. &
                & fracout(i,j,k) .ne. 0  .and. &
                & dbze(i,j,k)    .ge. -30.0    ) then
                ocbtm = .false.  !! cloud bottom detected
                kcbtm = k
                kctop = k
             endif
             if (       octop            .and. &  !! scan cloud-top
                & .not. ocbtm            .and. &  !! cloud-bottom already detected
                & fracout(i,j,k) .ne. 0  .and. &  !! exclude clear sky
                & dbze(i,j,k)    .ge. -30.0    ) then
                kctop = k  !! update
             endif
          enddo  !! k loop
          if( ocbtm ) cycle  !! cloud wasn't detected in this subcolumn

          !! check SLWC?
          if( temp(i,1,kctop) .lt. 273.15 ) cycle  !! return to the j (subcolumn) loop
          oslwc = .true.
          cmxdbz = -100._wp  !! initialized

          !CDIR NOLOOPCHG
          do k = kcbtm, kctop, -1
             cmxdbz = max( cmxdbz, dbze(i,j,k) )  !! column maximum dBZe update
             if ( fracout(i,j,k) .eq. 0 .or.  &
                & fracout(i,j,k) .eq. 2 .or.  &
                & dbze   (i,j,k) .lt. -30.0   ) then
                oslwc = .false.
             endif
          enddo
          if ( .not. oslwc ) cycle  !! return to the j (subcolumn) loop

          !! warm rain fraction (PDFMAP)
          iphase = 0
          if( cmxdbz .lt. -15.0 ) then
             iphase = 1  !! non-precipitating
          elseif( cmxdbz .ge. -15.0 .and. cmxdbz .lt. 0.0 ) then
             iphase = 2  !! drizzling
          elseif( cmxdbz .ge. 0.0 ) then
             iphase = 3  !! raining
          endif
          nslwc(i,iphase) = nslwc(i,iphase) + 1._wp

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
             !! retrieve the CFODD array based on ICOD and dBZe
             diagdbze = min( dbze(i,j,k), CFODD_XMAX )
             diagicod = min( diagicod,    CFODD_YMAX )
             ix = int( ( diagdbze - CFODD_XMIN ) / CFODD_DELX ) + 1
             iy = int( ( diagicod - CFODD_YMIN ) / CFODD_DELY ) + 1
             ix = max( 1, min( ix, CFODD_NBINX ) )  !! for safety
             iy = max( 1, min( iy, CFODD_NBINY ) )  !! for safety
             ncfodd(i,ix,iy,icls) = ncfodd(i,ix,iy,icls) + 1._wp
             ncfodd(i,ix,iy,isum) = ncfodd(i,ix,iy,isum) + 1._wp
          enddo

       enddo  ! j (Ncolumns)
    enddo     ! i (Npoints)

    !! summary statistics of warm-rain PDFs (nonprecip/drizzling/raining)
    do iphase = 1, PDFMAP_NPHASE
       do i = 1, Npoints, 1
          if ( sum( nslwc(i,:) ) .ge. 1.0 ) then
             pdf(i,iphase) = nslwc(i,iphase) / sum( nslwc(i,:) ) * 100.0  !! PDF % in each phase
          else  !! in case nslwc==0
             pdf(i,iphase) = R_UNDEF
          endif
       enddo
    enddo
    !! summary statistics of CFODDs
    do icls = 1, Nclass
       do iy = 1, CFODD_NBINY
          do ix = 1, CFODD_NBINX
             do i = 1, Npoints
                if ( sum( ncfodd(i,:,iy,icls) ) .ge. 1.0 ) then
                   cfodd(i,ix,iy,icls) = real( ncfodd(i,ix,iy,icls) ) / real( sum( ncfodd(i,:,iy,icls) ) )
                   cfodd(i,ix,iy,icls) = cfodd(i,ix,iy,icls) * 100._wp  !! %
                else
                   cfodd(i,ix,iy,icls) = R_UNDEF
                endif
             enddo
          enddo
       enddo
    enddo

  RETURN
  END SUBROUTINE COSP_IDID_WR


END MODULE MOD_COSP_STATS_IDID
