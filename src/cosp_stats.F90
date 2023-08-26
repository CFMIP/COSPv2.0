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
! Jul 2007 - A. Bodas-Salcedo - Initial version
! Jul 2008 - A. Bodas-Salcedo - Added capability of producing outputs in standard grid
! Oct 2008 - J.-L. Dufresne   - Bug fixed. Assignment of Npoints,Nlevels,Nhydro,Ncolumns 
!                               in COSP_STATS
! Oct 2008 - H. Chepfer       - Added PARASOL reflectance arguments
! Jun 2010 - T. Yokohata, T. Nishimura and K. Ogochi - Added NEC SXs optimisations
! Jan 2013 - G. Cesana        - Added betaperp and temperature arguments 
!                             - Added phase 3D/3Dtemperature/Map output variables in diag_lidar 
! May 2015 - D. Swales        - Modified for cosp2.0 
! Nov 2018 - T. Michibata     - Added CloudSat+MODIS Warmrain Diagnostics
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE MOD_COSP_STATS
  USE COSP_KINDS, ONLY: wp
  USE MOD_COSP_CONFIG, &
      ONLY: R_UNDEF,          R_GROUND,         &
            SGCLD_CLR,        SGCLD_ST,         &
            SGCLD_CUM,                          &
            CWP_THRESHOLD,    COT_THRESHOLD,    &
            CFODD_NDBZE,      CFODD_NICOD,      &
            CFODD_BNDRE,      CFODD_BNDZE,      &
            CFODD_NCLASS,     COT_NCLASS,       &
            CFODD_DBZE_MIN,   CFODD_DBZE_MAX,   &
            CFODD_ICOD_MIN,   CFODD_ICOD_MAX,   &
            CFODD_DBZE_WIDTH, CFODD_ICOD_WIDTH, &
            CFODD_HISTDBZE,   CFODD_HISTICOD,   &
            WR_NREGIME, NOBSTYPE,               &
            SLWC_NCOT,                          &
            SLWC_COT_MAX, SLWC_HISTCOT          
            
  USE COSP_PHYS_CONSTANTS,  ONLY: tmelt

  IMPLICIT NONE
CONTAINS

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !---------- SUBROUTINE COSP_CHANGE_VERTICAL_GRID ----------------
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_CHANGE_VERTICAL_GRID(Npoints,Ncolumns,Nlevels,zfull,zhalf,y,Nglevels,newgrid_bot,newgrid_top,r,log_units)
   implicit none
   ! Input arguments
   integer,intent(in) :: Npoints  !# of grid points
   integer,intent(in) :: Nlevels  !# of levels
   integer,intent(in) :: Ncolumns !# of columns
   real(wp),dimension(Npoints,Nlevels),intent(in) :: zfull ! Height at model levels [m] (Bottom of model layer)
   real(wp),dimension(Npoints,Nlevels),intent(in) :: zhalf ! Height at half model levels [m] (Bottom of model layer)
   real(wp),dimension(Npoints,Ncolumns,Nlevels),intent(in) :: y     ! Variable to be changed to a different grid
   integer,intent(in) :: Nglevels  !# levels in the new grid
   real(wp),dimension(Nglevels),intent(in) :: newgrid_bot ! Lower boundary of new levels  [m]
   real(wp),dimension(Nglevels),intent(in) :: newgrid_top ! Upper boundary of new levels  [m]
   logical,optional,intent(in) :: log_units ! log units, need to convert to linear units
   ! Output
   real(wp),dimension(Npoints,Ncolumns,Nglevels),intent(out) :: r ! Variable on new grid

   ! Local variables
   integer :: i,j,k
   logical :: lunits
   integer :: l
   real(wp) :: w ! Weight
   real(wp) :: dbb, dtb, dbt, dtt ! Distances between edges of both grids
   integer :: Nw  ! Number of weights
   real(wp) :: wt  ! Sum of weights
   real(wp),dimension(Nlevels) :: oldgrid_bot,oldgrid_top ! Lower and upper boundaries of model grid
   real(wp) :: yp ! Local copy of y at a particular point.
              ! This allows for change of units.

   lunits=.false.
   if (present(log_units)) lunits=log_units

   r = 0._wp

   do i=1,Npoints
     ! Calculate tops and bottoms of new and old grids
     oldgrid_bot = zhalf(i,:)
     oldgrid_top(1:Nlevels-1) = oldgrid_bot(2:Nlevels)
     oldgrid_top(Nlevels) = zfull(i,Nlevels) +  zfull(i,Nlevels) - zhalf(i,Nlevels) ! Top level symmetric
     l = 0 ! Index of level in the old grid
     ! Loop over levels in the new grid
     do k = 1,Nglevels
       Nw = 0 ! Number of weigths
       wt = 0._wp ! Sum of weights
       ! Loop over levels in the old grid and accumulate total for weighted average
       do
         l = l + 1
         w = 0.0 ! Initialise weight to 0
         ! Distances between edges of both grids
         dbb = oldgrid_bot(l) - newgrid_bot(k)
         dtb = oldgrid_top(l) - newgrid_bot(k)
         dbt = oldgrid_bot(l) - newgrid_top(k)
         dtt = oldgrid_top(l) - newgrid_top(k)
         if (dbt >= 0.0) exit ! Do next level in the new grid
         if (dtb > 0.0) then
           if (dbb <= 0.0) then
             if (dtt <= 0) then
               w = dtb
             else
               w = newgrid_top(k) - newgrid_bot(k)
             endif
           else
             if (dtt <= 0) then
               w = oldgrid_top(l) - oldgrid_bot(l)
             else
               w = -dbt
             endif
           endif
           ! If layers overlap (w/=0), then accumulate
           if (w /= 0.0) then
             Nw = Nw + 1
             wt = wt + w
             do j=1,Ncolumns
               if (lunits) then
                 if (y(i,j,l) /= R_UNDEF) then
                   yp = 10._wp**(y(i,j,l)/10._wp)
                 else
                   yp = 0._wp
                 endif
               else
                 yp = y(i,j,l)
               endif
               r(i,j,k) = r(i,j,k) + w*yp
             enddo
           endif
         endif
       enddo
       l = l - 2
       if (l < 1) l = 0
       ! Calculate average in new grid
       if (Nw > 0) then
         do j=1,Ncolumns
           r(i,j,k) = r(i,j,k)/wt
         enddo
       endif
     enddo
   enddo

   ! Set points under surface to R_UNDEF, and change to dBZ if necessary
   do k=1,Nglevels
     do j=1,Ncolumns
       do i=1,Npoints
         if (newgrid_top(k) > zhalf(i,1)) then ! Level above model bottom level
           if (lunits) then
             if (r(i,j,k) <= 0.0) then
               r(i,j,k) = R_UNDEF
             else
               r(i,j,k) = 10._wp*log10(r(i,j,k))
             endif
           endif
         else ! Level below surface
           r(i,j,k) = R_GROUND
         endif
       enddo
     enddo
   enddo

END SUBROUTINE COSP_CHANGE_VERTICAL_GRID

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !------------- SUBROUTINE COSP_LIDAR_ONLY_CLOUD -----------------
  ! (c) 2008, Lawrence Livermore National Security Limited Liability Corporation.
  ! All rights reserved.
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_LIDAR_ONLY_CLOUD(Npoints, Ncolumns, Nlevels, beta_tot, beta_mol,       &
     Ze_tot, lidar_only_freq_cloud, tcc, radar_tcc, radar_tcc2)
    ! Inputs
    integer,intent(in) :: &
         Npoints,       & ! Number of horizontal gridpoints
         Ncolumns,      & ! Number of subcolumns
         Nlevels          ! Number of vertical layers
    real(wp),dimension(Npoints,Nlevels),intent(in) :: &
         beta_mol         ! Molecular backscatter
    real(wp),dimension(Npoints,Ncolumns,Nlevels),intent(in) :: &
         beta_tot,      & ! Total backscattered signal
         Ze_tot           ! Radar reflectivity
    ! Outputs
    real(wp),dimension(Npoints,Nlevels),intent(out) :: &
         lidar_only_freq_cloud
    real(wp),dimension(Npoints),intent(out) ::&
         tcc,       & !
         radar_tcc, & !
         radar_tcc2   !
    
    ! local variables
    real(wp) :: sc_ratio
    real(wp),parameter :: &
         s_cld=5.0, &
         s_att=0.01
    integer :: flag_sat,flag_cld,pr,i,j,flag_radarcld,flag_radarcld_no1km,j_1km
    
    lidar_only_freq_cloud = 0._wp
    tcc = 0._wp
    radar_tcc = 0._wp
    radar_tcc2 = 0._wp
    do pr=1,Npoints
       do i=1,Ncolumns
          flag_sat = 0
          flag_cld = 0
          flag_radarcld = 0 !+JEK
          flag_radarcld_no1km=0 !+JEK
          ! look for j_1km from bottom to top
          j = 1
          do while (Ze_tot(pr,i,j) .eq. R_GROUND)
             j = j+1
          enddo
          j_1km = j+1  !this is the vertical index of 1km above surface  
          
          do j=1,Nlevels
             sc_ratio = beta_tot(pr,i,j)/beta_mol(pr,j)
             if ((sc_ratio .le. s_att) .and. (flag_sat .eq. 0)) flag_sat = j
             if (Ze_tot(pr,i,j) .lt. -30.) then  !radar can't detect cloud
                if ( (sc_ratio .gt. s_cld) .or. (flag_sat .eq. j) ) then  !lidar sense cloud
                   lidar_only_freq_cloud(pr,j)=lidar_only_freq_cloud(pr,j)+1. !top->surf
                   flag_cld=1
                endif
             else  !radar sense cloud (z%Ze_tot(pr,i,j) .ge. -30.)
                flag_cld=1
                flag_radarcld=1
                if (j .gt. j_1km) flag_radarcld_no1km=1              
             endif
          enddo !levels
          if (flag_cld .eq. 1) tcc(pr)=tcc(pr)+1._wp
          if (flag_radarcld .eq. 1) radar_tcc(pr)=radar_tcc(pr)+1.
          if (flag_radarcld_no1km .eq. 1) radar_tcc2(pr)=radar_tcc2(pr)+1.        
       enddo !columns
    enddo !points
    lidar_only_freq_cloud=lidar_only_freq_cloud/Ncolumns
    tcc=tcc/Ncolumns
    radar_tcc=radar_tcc/Ncolumns
    radar_tcc2=radar_tcc2/Ncolumns
    
    ! Unit conversion
    where(lidar_only_freq_cloud /= R_UNDEF) &
            lidar_only_freq_cloud = lidar_only_freq_cloud*100._wp
    where(tcc /= R_UNDEF) tcc = tcc*100._wp
    where(radar_tcc /= R_UNDEF) radar_tcc = radar_tcc*100._wp
    where(radar_tcc2 /= R_UNDEF) radar_tcc2 = radar_tcc2*100._wp
    
  END SUBROUTINE COSP_LIDAR_ONLY_CLOUD
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !-----------------------  SUBROUTINE COSP_DIAG_WARMRAIN ------------------------
  ! (c) 2018, Research Institute for Applied Mechanics (RIAM), Kyushu Univ.
  ! All rights reserved.
  ! * Purpose:    1) Diagnose Contoured Frequency by Optical Depth Diagram (CFODD)
  !                  from CloudSat Radar and MODIS retrievals.
  !               2) Diagnose Warm-Rain Occurrence Frequency (nonprecip/drizzle/rain)
  !                  from CloudSat Radar.
  ! * History:    Jan 2018 (T.Michibata): Test run
  !               May 2018 (T.Michibata): update for COSP2
  !               Sep 2018 (T.Michibata): modified I/O
  !               Nov 2018 (T.Michibata): minor revisions
  !               May 2020 (T.Michibata and X.Jing): bug-fix for frac_out dimsize
  !               Apr 2022 (T.Michibata): bug-fix for non-sunlit columns
  ! * References: Michibata et al. (GMD'19, doi:10.5194/gmd-12-4297-2019)
  !               Michibata et al. (GRL'20, doi:10.1029/2020GL088340)
  !               Suzuki et al. (JAS'10, doi:10.1175/2010JAS3463.1)
  !               Suzuki et al. (JAS'15, doi:10.1175/JAS-D-14-0265.1)
  !               Jing et al. (JCLIM'19, doi:10.1175/jcli-d-18-0789.1)
  ! * Contact:    Takuro Michibata (RIAM, Kyushu University, Japan).
  !               E-mail: michibata@riam.kyushu-u.ac.jp
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_DIAG_WARMRAIN( Npoints, Ncolumns, Nlevels,           & !! in
                                 temp,    zlev,                        & !! in
                                 lchnk, sunlit,                        & !! in
                                 lwp,     liqcot,   liqreff, liqcfrc,  & !! in
                                 iwp,     icecot,   icereff, icecfrc,  & !! in
                                 fracout, dbze,                        & !! in
                                 tautot_liq, tautot_ice,               & !! in
                                 lidarcldflag,                         & !! in
                                 cfodd_ntotal,                         & !! inout
                                 wr_occfreq_ntotal,                    & !! inout
                                 lsmallcot, mice, lsmallreff,          & !! inout
                                 lbigreff, nmultilcld, nfracmulti,     & !! inout
                                 nhetcld, coldct,                      & !! inout
                                 coldct_cal,calice, obs_ntotal,        & !! inout
                                 slwccot_ntotal )    !! inout
    integer,parameter :: &
         Nphase = 6 ! Number of CALIPSO cloud layer phase types
                    ! [ice,liquid,undefined,false ice,false liquid,Percent of ice]

    ! Inputs
    integer, intent(in) :: &
         Npoints,          & ! Number of horizontal gridpoints
         Ncolumns,         & ! Number of subcolumns
         Nlevels,          & ! Number of vertical layers
         lchnk               ! Local chunk id for debugging by location
    integer, dimension(Npoints), intent(in) ::  &
         sunlit              ! cospgridIN%sunlit
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
         zlev                ! altitude [m] for model level, above ground
    real(wp), dimension(Npoints,1,Nlevels),intent(in) :: &
         temp                ! temperature [K]
    real(wp), dimension(Npoints,Ncolumns,Nlevels),intent(in) :: &
         fracout,          & ! SCOPS subcolumn retrieval
         dbze,             & ! Radar reflectivity [dBZe]
         tautot_liq,       & ! CALIPSO optical thickness liquid
         tautot_ice,       & ! CALIPSO optical thickness ice
         lidarcldflag        ! CALIPSO cloud flag 
    !real(wp), dimension(Npoints,Nlevels,Nphase),intent(in) :: &
    !     lidarcldphase       ! CALIPSO cloud fraction by phase [2=liquid]

    ! Outputs
    real(wp),dimension(Npoints,CFODD_NDBZE,CFODD_NICOD,CFODD_NCLASS),intent(inout) :: &
         cfodd_ntotal        ! # of CFODD samples for each ICOD/dBZe bin
    real(wp),dimension(Npoints,WR_NREGIME),intent(inout) :: &
         wr_occfreq_ntotal   ! # of SLWC samples
    real(wp),dimension(Npoints),intent(inout) :: &
         lsmallcot,         & ! # of liquid clouds that don't meet COT condition
         mice,              & ! # of ice clouds
         lsmallreff,        & ! # of liquid clouds that have too small reff to meet SLWC conditions
         lbigreff,          & ! # of liquid clouds that have too big reff to meet SLWC conditions
         coldct,            & ! # of subcolumns with cloud top temp < 273 K
         coldct_cal,        & ! # of subcolumns with cloud top temp < 273 K detected by CALIPSO (and not MODIS)
         calice,            &  ! # of columns where Calipso detected ice that was not detected by MODIS
         nfracmulti            ! # of subcolumns where fracout indicates multilayer cloud
    real(wp),dimension(Npoints,2),intent(inout) :: &
         nmultilcld,        & ! # of multilayer cloud subcolumns, excluded from SLWC counts, 1 = MODIS/CloudSat detected, 2 = CALIPSO/CloudSat detected
         nhetcld              ! # of heterogenous clouds (stratocumulus above/below cumulus) in continuous layer
     real(wp),dimension(Npoints,NOBSTYPE),intent(inout) :: obs_ntotal     ! # of Observations
     real(wp),dimension(Npoints,SLWC_NCOT,COT_NCLASS),intent(inout) :: slwccot_ntotal ! # of MODIS liquid COT samples for SLWCs only @ each ICOD bin
     !real(wp),dimension(Npoints,SLWC_NCOT,3),intent(inout) :: slwccot_ntotal     ! # of MODIS liquid COT samples for SLWCs only @ each ICOD bin, MODIS/CALIPSO


    ! Local variables
    integer  :: i, j, k
    integer  :: ix, iy, cs
    integer  :: kctop, kcbtm
    integer  :: icls
    integer  :: iregime
    real     :: cmxdbz
    real(wp) :: diagcgt   !! diagnosed cloud geometric thickness [m]
    real(wp) :: diagdbze  !! diagnosed dBZe
    real(wp) :: diagicod  !! diagnosed in-cloud optical depth
    real(wp) :: cbtmh     !! diagnosed in-cloud optical depth
    real(wp), dimension(Ncolumns) :: slwccot_cls !! masking COT by class
    real(wp), dimension(Ncolumns,Nlevels) :: dbze_cls, icod_cls
    integer,  dimension(Ncolumns,Nlevels) :: scolcls, scolcls2,scolcls3 !! masking by class 
    real(wp), dimension(Npoints,Ncolumns,Nlevels) :: icod, icod_cal  !! in-cloud optical depth (ICOD)
    logical  :: octop, ocbtm, oslwc, multilcld, hetcld, modis_ice, fracmulti, icoldct, ulmodis
    integer, dimension(Npoints,Ncolumns,Nlevels) :: fracout_int  !! fracout (decimal to integer)
    integer  :: obstype   !! 1 = all-sky; 2 = clear-sky; 3 = cloudy-sky
    real(wp),dimension(Npoints,Ncolumns) :: slwccot     ! MODIS liquid COT for SLWCs only
    real(wp),dimension(Npoints,Ncolumns) :: slwccot_cal ! CALIPSO liquid COT for SLWCs only
    logical, dimension(Npoints) :: modis_cond           ! MODIS column-level conditions for detecting SLWCs
    logical,dimension(Npoints,Ncolumns) :: modiscs_multi,modiscs_coldct  ! flag for multilayer and supercooled liq clouds detected by MODIS/CloudSat

    fracout_int(:,:,:) = NINT( fracout(:,:,:) )  !! assign an integer subpixel ID (0=clear-sky; 1=St; 2=Cu)
    modis_cond(:) = .false.
    
    !print*,"zlev(1,Nlevels) = ", zlev(1,Nlevels)
    !print*,"zlev(1,1) = ", zlev(1,1)

    !! initialize
    slwccot(:,:) = R_UNDEF
    slwccot_cal(:,:) = R_UNDEF
    icod(:,:,:) = R_UNDEF
    icod_cal(:,:,:) = R_UNDEF
    modiscs_multi(:,:) = .false.
    modiscs_coldct(:,:) = .false.
    do i = 1, Npoints
!       if ( lwp(i) .eq. R_UNDEF ) then  ! for non-sunlit columns
       if ( sunlit(i) .le. 0 ) then  ! for non-sunlit columns
          cfodd_ntotal(i,:,:,:) = R_UNDEF
          wr_occfreq_ntotal(i,:) = R_UNDEF
!          icod(i,:,:) = R_UNDEF
!         icod_cal(i,:,:) = R_UNDEF
          lsmallcot(i) = R_UNDEF
          mice(i) = R_UNDEF
          lsmallreff(i) = R_UNDEF
          lbigreff(i) = R_UNDEF
          obs_ntotal(i,:) = R_UNDEF
          nhetcld(i,:) = R_UNDEF
          nmultilcld(i,:) = R_UNDEF
          nfracmulti(i) = R_UNDEF
          coldct(i) = R_UNDEF
          coldct_cal(i) = R_UNDEF
!          slwccot(i,:) = R_UNDEF
!          slwccot_cal(i,:) = R_UNDEF
          slwccot_ntotal(i,:,:) = R_UNDEF
          !slwccot_ntotal_cal(i,:,:) = R_UNDEF
          calice(i) = R_UNDEF
       else
          cfodd_ntotal(i,:,:,:)  = 0._wp
          wr_occfreq_ntotal(i,:) = 0._wp
 !         icod(i,:,:) = 0._wp
 !         icod_cal(i,:,:) = 0._wp
          lsmallcot(i) = 0._wp
          mice(i) = 0._wp
          lsmallreff(i) = 0._wp
          lbigreff(i) = 0._wp
          obs_ntotal(i,:) = 0._wp
          nhetcld(i,:) = 0._wp
          nmultilcld(i,:) = 0._wp
          nfracmulti(i) = 0._wp
          coldct(i) = 0._wp
          coldct_cal(i) = 0._wp
!          slwccot(i,:) = 0._wp
!          slwccot_cal(i,:) = 0._wp
          slwccot_ntotal(i,:,:) = 0._wp
          !slwccot_ntotal_cal(i,:,:) = 
          calice(i) = 0._wp
       endif
    enddo

    do i = 1, Npoints
        !! Total Sampling Frequency
       do j = 1, Ncolumns, 1
!          if( lwp(i).eq.R_UNDEF ) cycle ! remove non-sunlit columns
          if ( sunlit(i) .le. 0 ) cycle ! remove non-sunlit columns
          obs_ntotal(i,1) = obs_ntotal(i,1) + 1._wp  ! all-sky (# of all samples)
          obstype = 2                      ! initial flag (2 = clear sky)
          !CDIR NOLOOPCHG
          do k = 1, Nlevels - 2  !masking out through 960m above ground, groundclutter
             if ( dbze   (i,j,k) .ge. CFODD_DBZE_MIN .and. &
                & dbze   (i,j,k) .le. CFODD_DBZE_MAX ) then
                !.and. &
                !& fracout(i,j,k) .ne. SGCLD_CLR            ) then
                obstype = 3  ! cloudy sky
             endif
          enddo !Nlevels
          
          obs_ntotal(i,obstype) = obs_ntotal(i,obstype) + 1._wp
       enddo !Ncolumns

       !! check by MODIS retrieval
!       if ( lwp(i)     .eq.  R_UNDEF        .or.  &       !! exclude non-sunlit
       if ( sunlit(i)  .le.  0              .or.  &       !! exclude non-sunlit
          & lwp(i)     .le.  CWP_THRESHOLD  .or.  &
          & liqcot(i)  .le.  COT_THRESHOLD  .or.  &
          & liqreff(i) .lt.  CFODD_BNDRE(1) .or.  &
          & liqreff(i) .gt.  CFODD_BNDRE(4) .or.  &
          & iwp(i)     .gt.  CWP_THRESHOLD  .or.  &       !! exclude ice cloud
          & icecot(i)  .gt.  COT_THRESHOLD  .or.  &       !! exclude ice cloud
          & icereff(i) .gt.  CFODD_BNDRE(1)       ) then  !! exclude ice cloud
          modis_cond(i) = .true.
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
          icoldct = .false. !! CMB
          !CDIR NOLOOPCHG
          do k = Nlevels-2, 1, -1  !! scan from cloud-bottom to cloud-top, masking out through 960m above ground
             if ( dbze(i,j,k) .eq. R_GROUND .or. &
                  dbze(i,j,k) .eq. R_UNDEF       ) cycle
             if ( ocbtm                              .and. &
                !& fracout_int(i,j,k) .ne. SGCLD_CLR  .and. &
                & dbze(i,j,k)        .ge. CFODD_DBZE_MIN   ) then
                ocbtm = .false.  !! cloud bottom detected
                kcbtm = k
                kctop = k
             endif
             if (       octop                        .and. &  !! scan cloud-top
                & .not. ocbtm                        .and. &  !! cloud-bottom already detected
                !& fracout_int(i,j,k) .ne. SGCLD_CLR  .and. &  !! exclude clear sky
                & dbze(i,j,k)        .ge. CFODD_DBZE_MIN   ) then
                kctop = k  !! update
             endif
          enddo  !! k loop

          if( ocbtm )  cycle  !! cloud wasn't detected in this subcolumn
          !! check SLWC?
          if( temp(i,1,kctop) .lt. tmelt ) then
              coldct(i) = coldct(i) + 1._wp 
              modiscs_coldct(i,j) = .true.
              icoldct = .true.
          endif
          !if( temp(i,1,kctop) .lt. tmelt ) cycle  !! return to the j (subcolumn) loop
          oslwc = .true.
          hetcld = .false.
          multilcld = .false.
          fracmulti = .false.
          cmxdbz = CFODD_DBZE_MIN  !! initialized

          !CDIR NOLOOPCHG
          do k = kcbtm, kctop, -1
             cmxdbz = max( cmxdbz, dbze(i,j,k) )  !! column maximum dBZe update
             !if ( fracout_int(i,j,k) .eq. SGCLD_CLR  .or.  &
             !   & fracout_int(i,j,k) .eq. SGCLD_CUM  .or.  &
             !   & dbze       (i,j,k) .lt. CFODD_DBZE_MIN   ) then
             !   oslwc = .false.
             !endif

             if ( fracout_int(i,j,k) .eq. SGCLD_CUM .and.  &
                & .not. multilcld ) then
                hetcld = .true.
                !oslwc = .false.
            endif
            if ( fracout_int(i,j,k) .eq. SGCLD_CLR) then
                fracmulti = .true.
            endif
            if ( dbze       (i,j,k) .lt. CFODD_DBZE_MIN ) then
                multilcld = .true.
                modiscs_multi(i,j) = .true.
                oslwc = .false.
            endif
          enddo
          
          if ( multilcld ) then
             nmultilcld(i,1) = nmultilcld(i,1) + 1._wp
          endif
          
          if ( hetcld ) then
             nhetcld(i,1) = nhetcld(i,1) + 1._wp
          endif
          
          if (fracmulti) then
             nfracmulti(i) = nfracmulti(i) + 1._wp
          endif
          
          if ( .not. oslwc ) cycle  !! return to the j (subcolumn) loop

          !! warm-rain occurrence frequency
          iregime = 0
          if( cmxdbz .lt. CFODD_BNDZE(1) .and. .not. icoldct ) then !.and. &
              !& .not. fracmulti ) then
             iregime = 1  !! non-precipitating
          elseif( (cmxdbz .ge. CFODD_BNDZE(1)) .and. (cmxdbz .lt. CFODD_BNDZE(2)) &
                  .and. .not. icoldct ) then !.and. .not. fracmulti ) then
             iregime = 2  !! drizzling
          elseif( cmxdbz .ge. CFODD_BNDZE(2) .and. .not. icoldct ) then!.and. &
                 !& .not. fracmulti ) then
             iregime = 3  !! raining
          elseif ( cmxdbz .lt. CFODD_BNDZE(1) .and. icoldct ) then
             iregime = 4  !! cold cloud top, non-precip
          elseif ( (cmxdbz .ge. CFODD_BNDZE(1)) .and. (cmxdbz .lt. CFODD_BNDZE(2)) &
                  & .and. icoldct ) then
             iregime = 5  !! cold cloud top, drizzling
          elseif ( cmxdbz .ge. CFODD_BNDZE(2) .and. icoldct ) then
             iregime = 6 !! cold cloud top, raining
          elseif ( cmxdbz .lt. CFODD_BNDZE(1) .and. fracmulti ) then
              iregime = 7
          elseif ( (cmxdbz .ge. CFODD_BNDZE(1)) .and. (cmxdbz .lt. CFODD_BNDZE(2)) &
                   .and. fracmulti ) then
              iregime = 8
          elseif ( cmxdbz .ge. CFODD_BNDZE(2) .and. fracmulti ) then
              iregime = 9
          endif
          wr_occfreq_ntotal(i,iregime) = wr_occfreq_ntotal(i,iregime) + 1._wp

          !! retrievals of ICOD and dBZe bin for CFODD plane
          diagcgt = zlev(i,kctop) - zlev(i,kcbtm)
          cbtmh   = zlev(i,kcbtm)
          !CDIR NOLOOPCHG
          do k = kcbtm, kctop, -1
             if( k .eq. kcbtm ) then
                diagicod = liqcot(i)
                slwccot(i,j) = min( liqcot(i), SLWC_COT_MAX )
             else
                diagicod = liqcot(i) * ( 1._wp - ( (zlev(i,k)-cbtmh)/diagcgt)**(5._wp/3._wp) )
             endif
             icod(i,j,k) = min( diagicod, CFODD_ICOD_MAX )
          enddo
          
           !! retrieve the CFODD array based on Reff
          icls = 0
          scolcls(:,:) = 0
          scolcls2(:,:) = 0 
          scolcls3(:,:) = 0 
          
          if (    liqreff(i) .ge. CFODD_BNDRE(1) .and. liqreff(i) .lt. CFODD_BNDRE(2) .and. &
               .not. icoldct ) then !.and. .not. fracmulti ) then
             icls = 1
          elseif( liqreff(i) .ge. CFODD_BNDRE(2) .and. liqreff(i) .lt. CFODD_BNDRE(3) .and. &
               .not. icoldct ) then !.and. .not. fracmulti ) then
             icls = 2
          elseif( liqreff(i) .ge. CFODD_BNDRE(3) .and. liqreff(i) .le. CFODD_BNDRE(4) .and. &
               .not. icoldct ) then ! .and. .not. fracmulti ) then
             icls = 3
          endif
          
          scolcls(j,1:Nlevels) = icls   ! save class assignment for each subcolumn for histogram
          
          icls = 0
          ! Generating CFODD for only SLWCs with max reflectivity < 20 dBZ, cot < 20 for linear
          ! regression, Suzuki et al. (2010)
          if ( liqreff(i) .ge. CFODD_BNDRE(1) .and. liqreff(i) .lt. CFODD_BNDRE(2) .and. &
               .not. icoldct .and. cmxdbz .lt. CFODD_DBZE_MAX .and. liqcot(i) < 20._wp) then
             icls = 4  ! small Reff size bin, only SLWCs with max reflectivity < 20 dBZ, cot < 20
          elseif( liqreff(i) .ge. CFODD_BNDRE(2) .and. liqreff(i) .lt. CFODD_BNDRE(3) .and. &
                 .not. icoldct .and. cmxdbz .lt. CFODD_DBZE_MAX .and. liqcot(i) < 20._wp) then
             icls = 5  ! medium Reff size bin, only SLWCs with max reflectivity < 20 dBZ, cot < 20
          elseif( liqreff(i) .ge. CFODD_BNDRE(3) .and. liqreff(i) .le. CFODD_BNDRE(4) .and. &
                  .not. icoldct .and. cmxdbz .lt. CFODD_DBZE_MAX .and. liqcot(i) < 20._wp) then
             icls = 6  ! large Reff size bin, only SLWCs with max reflectivity < 20 dBZ, cot < 20
          endif
          
          scolcls2(j,1:Nlevels) = icls   ! save class assignment for each subcolumn for histogram
          
          icls = 0
          ! Generating CFODD for only SLWCs with max reflectivity < 20 dBZ, 4 <= cot < 20 for linear
          ! regression, Suzuki et al. (2010)
          if ( liqreff(i) .ge. CFODD_BNDRE(1) .and. liqreff(i) .lt. CFODD_BNDRE(2) .and. &
               .not. icoldct .and. cmxdbz .lt. CFODD_DBZE_MAX .and. liqcot(i) .ge. 4._wp &
               .and. liqcot(i) < 20._wp ) then
             icls = 7  ! small Reff size bin, only SLWCs with max reflectivity < 20 dBZ, cot < 20
          elseif( liqreff(i) .ge. CFODD_BNDRE(2) .and. liqreff(i) .lt. CFODD_BNDRE(3) .and. &
                 .not. icoldct .and. cmxdbz .lt. CFODD_DBZE_MAX .and. liqcot(i) .ge. 4._wp  &
                 .and. liqcot(i) < 20._wp) then
             icls = 8  ! medium Reff size bin, only SLWCs with max reflectivity < 20 dBZ, cot < 20
          elseif( liqreff(i) .ge. CFODD_BNDRE(3) .and. liqreff(i) .le. CFODD_BNDRE(4) .and. &
                  .not. icoldct .and. cmxdbz .lt. CFODD_DBZE_MAX .and. liqcot(i) .ge. 4._wp &
                  .and. liqcot(i) < 20._wp) then
             icls = 9  ! large Reff size bin, only SLWCs with max reflectivity < 20 dBZ, cot < 20
          endif
          
          scolcls3(j,1:Nlevels) = icls   ! save class assignment for each subcolumn for histogram

       enddo  ! j (Ncolumns)
       
       ! Generate original CFODD
       do cs = 1,3,1
          !! initialize
          dbze_cls = dbze(i,1:Ncolumns,1:Nlevels)
          icod_cls = icod(i,1:Ncolumns,1:Nlevels)
          slwccot_cls = slwccot(i,1:Ncolumns)
          
          !! mask out subcolumns not in class i
          where (scolcls(1:Ncolumns,1:Nlevels) .ne. cs)
              dbze_cls(1:Ncolumns,1:Nlevels) = R_UNDEF
              icod_cls(1:Ncolumns,1:Nlevels) = R_UNDEF
          end where
          
          where (scolcls(1:Ncolumns,1) .ne. cs)
              slwccot_cls(1:Ncolumns) = R_UNDEF
          end where
          
          call hist2d( dbze_cls(1:Ncolumns,1:Nlevels), icod_cls(1:Ncolumns,1:Nlevels), &
                      & Ncolumns*Nlevels,                                          &
                      & CFODD_HISTDBZE, CFODD_NDBZE, CFODD_HISTICOD, CFODD_NICOD,  &
                      & cfodd_ntotal( i, 1:CFODD_NDBZE, 1:CFODD_NICOD, cs )      )
          
          slwccot_ntotal(i, 1:SLWC_NCOT, cs) = hist1d( Ncolumns,                 &
                     slwccot_cls(1:Ncolumns), SLWC_NCOT, SLWC_HISTCOT         )       
                      
       enddo ! cs (classes)
       
       ! Generate additional CFODD for just SLWCs with cmxdbz < 20
       do cs = 4,6,1
          !! initialize
          dbze_cls = dbze(i,1:Ncolumns,1:Nlevels)
          icod_cls = icod(i,1:Ncolumns,1:Nlevels)
          slwccot_cls = slwccot(i,1:Ncolumns)
          
          !! mask out subcolumns not in class i
          where (scolcls2(1:Ncolumns,1:Nlevels) .ne. cs)
              dbze_cls(1:Ncolumns,1:Nlevels) = R_UNDEF
              icod_cls(1:Ncolumns,1:Nlevels) = R_UNDEF
          end where
          
          where (scolcls2(1:Ncolumns,1) .ne. cs)
              slwccot_cls(1:Ncolumns) = R_UNDEF
          end where
          
          call hist2d( dbze_cls(1:Ncolumns,1:Nlevels), icod_cls(1:Ncolumns,1:Nlevels), &
                      & Ncolumns*Nlevels,                                          &
                      & CFODD_HISTDBZE, CFODD_NDBZE, CFODD_HISTICOD, CFODD_NICOD,  &
                      & cfodd_ntotal( i, 1:CFODD_NDBZE, 1:CFODD_NICOD, cs )      )
          
          slwccot_ntotal(i, 1:SLWC_NCOT, cs) = hist1d( Ncolumns,                 &
                     slwccot_cls(1:Ncolumns), SLWC_NCOT, SLWC_HISTCOT         )       
                      
       enddo ! cs (classes)
       
       ! Generate additional CFODD for just SLWCs with cmxdbz < 20 and 4 <= COT < 20
       do cs = 7,9,1
          !! initialize
          dbze_cls = dbze(i,1:Ncolumns,1:Nlevels)
          icod_cls = icod(i,1:Ncolumns,1:Nlevels)
          slwccot_cls = slwccot(i,1:Ncolumns)
          
          !! mask out subcolumns not in class i
          where (scolcls3(1:Ncolumns,1:Nlevels) .ne. cs)
              dbze_cls(1:Ncolumns,1:Nlevels) = R_UNDEF
              icod_cls(1:Ncolumns,1:Nlevels) = R_UNDEF
          end where
          
          where (scolcls3(1:Ncolumns,1) .ne. cs)
              slwccot_cls(1:Ncolumns) = R_UNDEF
          end where
          
          call hist2d( dbze_cls(1:Ncolumns,1:Nlevels), icod_cls(1:Ncolumns,1:Nlevels), &
                      & Ncolumns*Nlevels,                                          &
                      & CFODD_HISTDBZE, CFODD_NDBZE, CFODD_HISTICOD, CFODD_NICOD,  &
                      & cfodd_ntotal( i, 1:CFODD_NDBZE, 1:CFODD_NICOD, cs )      )
          
          slwccot_ntotal(i, 1:SLWC_NCOT, cs) = hist1d( Ncolumns,                 &
                     slwccot_cls(1:Ncolumns), SLWC_NCOT, SLWC_HISTCOT         )       
                      
       enddo ! cs (classes)
          
       !! # of samples for CFODD (joint 2d-histogram dBZe vs ICOD)
!       call hist2d( dbze(i,1:Ncolumns,1:Nlevels), icod(i,1:Ncolumns,1:Nlevels), &
!                  & Ncolumns*Nlevels,                                           &
!                  & CFODD_HISTDBZE, CFODD_NDBZE, CFODD_HISTICOD, CFODD_NICOD,   &
!                  & cfodd_ntotal( i, 1:CFODD_NDBZE, 1:CFODD_NICOD, icls )       )
    
!       slwccot_ntotal(i, 1:SLWC_NCOT, icls) = hist1d( Ncolumns,                     &
!                     slwccot(i,1:Ncolumns), SLWC_NCOT, SLWC_HISTCOT         )

    enddo     ! i (Npoints)

    ! Detection of SLWCs using CALIPSO/MODIS, only those that were not detectected by CloudSat/MODIS.
    ! In other words, if CloudSat/MODIS detected an SLWC anywhere in the column, this column will be excluded from 
    ! the following CALIPSO/MODIS analysis of SLWCs to avoid double-counting. The  subcolumn logic is
    ! based on a subcolumn CALIPSO cloud flag derived from the cospOUT%cal_lidarcldphase variable
    ! NOTE: CFODD is not possible from this section because it is CALIPSO/MODIS only
  
    do i  = 1, Npoints
        modis_ice = .false.
        ulmodis = .false.
        if( sunlit(i) .le. 0 .or.  & ! remove non-sunlit columns for MODIS detection of SLWCs 
            sum(wr_occfreq_ntotal(i,1:WR_NREGIME) ) .ge. 1.0_wp ) cycle !exclude columns where MODIS/CloudSat detected SLWCs already
           ! iwp(i)     .gt.  CWP_THRESHOLD  .or.  &       !! exclude columns with ice clouds detected by MODIS
           ! icecot(i)  .gt.  COT_THRESHOLD  .or.  &       !! exclude columns ice clouds detected by MODIS
           ! icereff(i) .gt.  CFODD_BNDRE(1)       ) cycle !! exclude coulmns with ice clouds detected by MODIS

        if( iwp(i)     .gt.  CWP_THRESHOLD  .or.  &       !! exclude columns with ice clouds detected by MODIS
            icecot(i)  .gt.  COT_THRESHOLD  .or.  &       !! exclude columns ice clouds detected by MODIS
            icereff(i) .gt.  CFODD_BNDRE(1)       ) then  !! exclude coulmns with ice clouds detected by MODIS
            modis_ice = .true.                           
        endif
        
        if ( modis_ice ) cycle 
        
        if( liqcot(i)  .le.  COT_THRESHOLD  .or.  &
            lwp(i) .le. CWP_THRESHOLD       .or. &
            liqreff(i) .lt.  CFODD_BNDRE(1) .or.  &
            liqreff(i) .gt.  CFODD_BNDRE(4)      ) then
            ulmodis = .true. !! if MODIS does not detect a cloud, check if CALIPSO detects
        endif

       !CDIR NOLOOPCHG
       do j = 1, Ncolumns, 1
          octop = .true.  !! initialize
          ocbtm = .true.  !! initialize
          kcbtm =     0   !! initialize
          kctop =     0   !! initialize
          if ( MAXVAL(tautot_ice(i,j,1:Nlevels)) .gt. COT_THRESHOLD & 
               .and. .not. modis_ice ) then !Calipso detected ice, but not excluding for consistency with obs
             calice(i) = calice(i) + 1._wp
          endif

        !  if ( MAXVAL(tautot_ice(i,j,1:Nlevels)) .gt. COT_THRESHOLD ) cycle !exclude if Calipso detected ice

        !  if ( MAXVAL(tautot_liq(i,j,1:Nlevels)) .le. COT_THRESHOLD ) cycle !exclude subcolumn if clearsky
      
                
          !CDIR NOLOOPCHG
          do k = Nlevels, 1, -1  !! scan from cloud-bottom to cloud-top
             if ( lidarcldflag(i,j,k) .eq. R_UNDEF ) cycle
             if ( ocbtm                              .and. &
                & lidarcldflag(i,j,k) .eq. 1._wp   ) then
                ocbtm = .false.  !! cloud bottom detected
                kcbtm = k
                kctop = k
             endif
             if (       octop    .and. &  !! scan cloud-top
                 & .not. ocbtm   .and. &  !! cloud-bottom already detected
                 & ( lidarcldflag(i,j,k) .eq. 1._wp ) ) then
                kctop = k  !! update
             endif
          enddo  !! k loop         

          if( ocbtm )  cycle  !! cloud wasn't detected in this subcolumn
          !! check SLWC?
          !! If cold cloud top and not detected already by MODIS/CloudSat, add to coldct_cal
          !! Supercooled liquid cloud fraction requires phase detection, which is only available through
          !! MODIS, so MODIS detection of cloud is required for the supercooled liquid frequency here
          if( temp(i,1,kctop) .lt. tmelt .and. .not. modiscs_coldct(i,j) .and. .not. ulmodis) then 
              coldct_cal(i) = coldct_cal(i) + 1._wp 
          endif

          if( temp(i,1,kctop) .lt. tmelt ) cycle  !! return to the j (subcolumn) loop
          oslwc = .true.
          multilcld = .false.

          !CDIR NOLOOPCHG
          do k = kcbtm, kctop, -1
             if ( lidarcldflag(i,j,k) .ne. 1._wp ) then
                multilcld = .true.
                oslwc = .false.
             endif
          enddo

          !If CALIPSO detected a multilayer cloud that MODIS/CloudSat did not, add to multilcld_cal
          if (multilcld .and. (.not. modiscs_multi(i,j) ) ) then
             nmultilcld(i,2) = nmultilcld(i,2) + 1._wp
          endif
          
          if ( .not. oslwc ) cycle  !! return to the j (subcolumn) loop

          !! warm-rain occurrence frequency
          !! If MODIS/CALIPSO detected the SLWC, add to npdfslwc_mcal
          !! and retrieve COT binned by cloud top Reff
          if (.not. ulmodis) then 
              iregime = 10
              wr_occfreq_ntotal(i,iregime) = wr_occfreq_ntotal(i,iregime) + 1._wp

              !! retrievals of ICOD and dBZe bin for CFODD plane
              diagcgt = zlev(i,kctop) - zlev(i,kcbtm)
              cbtmh   = zlev(i,kcbtm)
              !CDIR NOLOOPCHG
              do k = kcbtm, kctop, -1
                 if( k .eq. kcbtm ) then
                    diagicod = liqcot(i)
                    slwccot_cal(i,j) = min (liqcot(i), SLWC_COT_MAX)
                 else
                    diagicod = liqcot(i) * ( 1._wp - ( (zlev(i,k)-cbtmh)/diagcgt)**(5._wp/3._wp) )
                 endif
                 icod_cal(i,j,k) = min( diagicod, CFODD_ICOD_MAX )
              enddo ! k (Nlevels)
          
          !If only CALIPSO detected the SLWC, add to npdfslwc_calonly
          elseif (ulmodis) then
              iregime = 11
              wr_occfreq_ntotal(i,iregime) = wr_occfreq_ntotal(i,iregime) + 1._wp
          endif

       enddo  ! j (Ncolumns)
       
       !! retrieve the COT array based on Reff, if SLWCs detected by MODIS/CALIPSO
       icls = 0
       if (.not. ulmodis) then
           if ( liqreff(i) .ge. CFODD_BNDRE(1) .and. liqreff(i) .lt. CFODD_BNDRE(2) ) then
               icls = 10
           elseif ( liqreff(i) .ge. CFODD_BNDRE(2)  .and.  liqreff(i) .lt. CFODD_BNDRE(3)  ) then
               icls = 11
           elseif ( liqreff(i) .ge. CFODD_BNDRE(3) ) then
               icls = 12
           endif
       
       slwccot_ntotal(i, 1:SLWC_NCOT, icls) = hist1d( Ncolumns,                  &
                      slwccot_cal(i,1:Ncolumns), SLWC_NCOT, SLWC_HISTCOT         )
      endif   
    
    enddo     ! i (Npoints)

        
       
    !! CMB adding counts of other cloud types to assess frequency of single-layer warm phase clouds
    do i = 1, Npoints
!       if (lwp(i) .eq. R_UNDEF) cycle !remove non-sunlit columns
       if (sunlit(i) .le. 0 ) cycle !remove non-sunlit columns
             if ( lwp(i)     .gt.  CWP_THRESHOLD .and. &
                & liqcot(i)  .le.  COT_THRESHOLD  .and.  &
                & iwp(i)     .le.  CWP_THRESHOLD  .and.  &
                & icecot(i)  .le.  COT_THRESHOLD  .and.  &
                & icereff(i) .le.  CFODD_BNDRE(1)        ) then !liquid small COT
                lsmallcot(i) = lsmallcot(i) + 1._wp
             endif
             
             if ( iwp(i)     .gt. CWP_THRESHOLD   .or.  &
                & icecot(i)  .gt. COT_THRESHOLD   .or.  &
                & icereff(i) .gt. CFODD_BNDRE(1)        ) then !meets ice condition
                mice(i) = mice(i) + 1._wp
             endif

             if ( lwp(i)     .gt. CWP_THRESHOLD .and. &
                & liqcot(i)  .gt. COT_THRESHOLD   .and. &
                & liqreff(i) .lt. CFODD_BNDRE(1)  .and. &
                & iwp(i)     .le.  CWP_THRESHOLD  .and. &
                & icecot(i)  .le.  COT_THRESHOLD  .and. &
                & icereff(i) .le.  CFODD_BNDRE(1)        ) then !liquid small reff
                lsmallreff(i) = lsmallreff(i) + 1._wp
             endif
         
             if ( lwp(i)     .gt. CWP_THRESHOLD   .and. &
                & liqcot(i)  .gt. COT_THRESHOLD   .and. &
                & liqreff(i) .gt. CFODD_BNDRE(4)  .and. & 
                & iwp(i)     .le.  CWP_THRESHOLD  .and. &
                & icecot(i)  .le.  COT_THRESHOLD  .and. &
                & icereff(i) .le.  CFODD_BNDRE(1)        ) then !liquid big reff
                lbigreff(i) = lbigreff(i) + 1._wp
             endif
              
    enddo !i (Npoints)

  RETURN
  END SUBROUTINE COSP_DIAG_WARMRAIN

  ! ######################################################################################
  ! FUNCTION hist1D
  ! ######################################################################################
  function hist1d(Npoints,var,nbins,bins)
    ! Inputs
    integer,intent(in) :: &
         Npoints, & ! Number of points in input array
         Nbins      ! Number of bins for sorting
    real(wp),intent(in),dimension(Npoints) :: &
         var        ! Input variable to be sorted
    real(wp),intent(in),dimension(Nbins+1) :: &
         bins       ! Histogram bins [lowest,binTops]  
    ! Outputs
    real(wp),dimension(Nbins) :: &
         hist1d     ! Output histogram      
    ! Local variables
    integer :: ij
    
    do ij=2,Nbins+1  
       hist1D(ij-1) = count(var .ge. bins(ij-1) .and. var .lt. bins(ij))
       if (count(var .eq. R_GROUND) .ge. 1) hist1D(ij-1)=R_UNDEF
    enddo
    
  end function hist1D
  
  ! ######################################################################################
  ! SUBROUTINE hist2D
  ! ######################################################################################
  subroutine hist2D(var1,var2,npts,bin1,nbin1,bin2,nbin2,jointHist)
    implicit none
    
    ! INPUTS
    integer, intent(in) :: &
         npts,  & ! Number of data points to be sorted
         nbin1, & ! Number of bins in histogram direction 1 
         nbin2    ! Number of bins in histogram direction 2
    real(wp),intent(in),dimension(npts) :: &
         var1,  & ! Variable 1 to be sorted into bins
         var2     ! variable 2 to be sorted into bins
    real(wp),intent(in),dimension(nbin1+1) :: &
         bin1     ! Histogram bin 1 boundaries
    real(wp),intent(in),dimension(nbin2+1) :: &
         bin2     ! Histogram bin 2 boundaries
    ! OUTPUTS
    real(wp),intent(out),dimension(nbin1,nbin2) :: &
         jointHist
    
    ! LOCAL VARIABLES
    integer :: ij,ik
    
    do ij=2,nbin1+1
       do ik=2,nbin2+1
          jointHist(ij-1,ik-1)=count(var1 .ge. bin1(ij-1) .and. var1 .lt. bin1(ij) .and. &
               var2 .ge. bin2(ik-1) .and. var2 .lt. bin2(ik))        
       enddo
    enddo
  end subroutine hist2D
END MODULE MOD_COSP_STATS
