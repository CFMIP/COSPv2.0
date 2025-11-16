! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copyright (c) 2016, Regents of the University of Colorado
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
! History
! March 2016 - D. Swales - Original version
!   May 2018 - T. Michibata - Inline Diagnostic Driver (IDiD)
!   Sep 2018 - T. Michibata - modified IDiD output
!   Oct 2018 - T. Michibata - for MIROC6 interface
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 SUBROUTINE COSPEXE (            &
                      Npoints,   &
                      Nlevels,   &
                      lon,       &
                      lat,       &
                      p,         &
                      ph,        &
                      zlev,      &
                      zlev_half, &
                      T,         &
                      rh,        &
                      sh,        &
                      cca,       &
                      tca,       &
                      skt,       &
                      surfelev,  &
                      landmask,  &
                      mr_ozone,  &
                      sunlit,    &
                      mr_lsliq,  &
                      mr_lsice,  &
                      mr_ccliq,  &
                      mr_ccice,  &
#ifdef OPT_CHIMERRA
                      mr_lsrain, &
                      mr_lssnow, &
#ifdef OPT_GRPL
                      mr_lsgrpl, &
#endif
#endif
                      fl_lsrain, &
                      fl_lssnow, &
                      fl_lsgrpl, &
                      fl_ccrain, &
                      fl_ccsnow, &
                      Reff,      &
                      dtau_s,    &
                      dtau_c,    &
                      dem_s,     &
                      dem_c,     &
                      emsfc_lw,  &
#ifdef OPT_DPLRW
                      gwvel,     &
                      gcumf,     &
#endif
                      ocospchk,  &
                      ocospudf,  &
                      cosp_input_nl, &
                      cosp_output_nl, &
                      time )

  use cosp_kinds,          only: wp                         
  USE MOD_COSP_CONFIG,     ONLY: R_UNDEF,PARASOL_NREFL,LIDAR_NCAT,SR_BINS,                &
                                 N_HYDRO,RTTOV_MAX_CHANNELS,numMISRHgtBins,               &
                                 DBZE_BINS,LIDAR_NTEMP,calipso_histBsct,                  &
                                 CFODD_NDBZE,      CFODD_NICOD,                           &
                                 CFODD_BNDRE,      CFODD_NCLASS,                          &
                                 CFODD_DBZE_MIN,   CFODD_DBZE_MAX,                        &
                                 CFODD_ICOD_MIN,   CFODD_ICOD_MAX,                        &
                                 CFODD_DBZE_WIDTH, CFODD_ICOD_WIDTH,                      &
                                 WR_NREGIME,                                              &
                                 numMODISTauBins,numMODISPresBins,                        &
                                 numMODISReffIceBins,numMODISReffLiqBins,                 &
                                 numISCCPTauBins,numISCCPPresBins,numMISRTauBins,         &
                                 ntau,npres,modis_histTau,tau_binBounds,                  & !<- modified by Doppler
                                 modis_histTauEdges,tau_binEdges,                         &
                                 modis_histTauCenters,tau_binCenters,ntauV1p4,            &
                                 tau_binBoundsV1p4,tau_binEdgesV1p4, tau_binCentersV1p4,  &
                                 vgrid_zl,vgrid_zu,vgrid_z,   &  !Jing X.
                                 Nlvgrid,   &
                                 cloudsat_preclvl !<-added
  use cosp_phys_constants, only: amw,amd,amO3,amCO2,amCH4,amN2O,amCO
!  use mod_cosp_io,         only: nc_read_input_file !,nc_cmor_init   !do NOT use cmor. Jing X.
  USE mod_quickbeam_optics,only: size_distribution,hydro_class_init,quickbeam_optics
  use quickbeam,           only: radar_cfg
  use mod_cosp,            only: cosp_init,cosp_optical_inputs,cosp_column_inputs,        &
                                 cosp_outputs,cosp_cleanUp,cosp_simulator,linitialization
  USE mod_rng,             ONLY: rng_state, init_rng
  USE mod_scops,           ONLY: scops
  USE mod_prec_scops,      ONLY: prec_scops
  USE MOD_COSP_UTILS,      ONLY: cosp_precip_mxratio
  use cosp_optics,         ONLY: cosp_simulator_optics,lidar_optics,modis_optics,         &
                                 modis_optics_partition
#ifdef OPT_DPLRW
  use mod_cosp_config,     only: Nlvdplr,lvdplr_MIN,lvdplr_MAX,lvdplr_WID, &
                                 Nlvspwd,lvspwd_MIN,lvspwd_MAX,lvspwd_WID, &
                                 NlvdBZe,lvdBZe_MIN,lvdBZe_MAX,lvdBZe_WID, &
                                 Nlvtemp,lvtemp_MIN,lvtemp_MAX,lvtemp_WID, &
                                 N_ISCCP
#endif

  ! OUTPUT ONLY
!  use mod_cosp_io,         only: var1d,var2d,var3d,construct_cospOutList !,nc_cmor_init,    &
                                 !nc_cmor_associate_1d,nc_cmor_write_1d,nc_cmor_close,     &
                                 !nc_cmor_associate_2d,nc_cmor_write_2d  !do Not use cmor. Jing.
  
  implicit none

  ! Input/Output driver file control
  !character(len=64),parameter :: &
  !     cosp_input_namelist  = 'cosp_input_nl.v2.0.txt', &
  !     cosp_output_namelist = 'cosp_output_nl_v2.0.txt'

  ! Sample input data variables
  integer,intent(in) ::&
       Npoints, &    ! Number of gridpoints (lat*lon)
       Nlevels       ! Number of levels
  real(wp),dimension(Npoints),intent(in) :: &
       lon, &        ! longitude [degrees East]
       lat, &        ! latitude  [degrees North]
       skt, &        ! skin temperature [K]
       surfelev, &   ! Surface Elevation (m) !<-added
       landmask, &   ! landmask [0 - ocean, 1 - land]
       sunlit        ! 1 for day points, 0 for night points
 ! real(wp),dimension(Npoints) :: & !<-added
 !      surfelev      ! Surface Elevation (m) !<-added
  real(wp),dimension(Npoints,Nlevels),intent(inout) :: &
       p, &          ! Pressure at full model level [Pa]
!       ph, &         ! Pressure at half model levels [Pa]
       zlev, &       ! Height at model levels [m]
!       zlev_half, &  ! Height at half model levels [m] (Bottom of model layer)
       T, &          ! Temperature at model levels [K]
       rh, &         ! relative humidity
       sh, &         ! specific humidity of water vapor [kg/kg]
       cca, &        ! convective cloud fraction
       tca, &        ! total cloud fraction
       mr_ozone, &   ! ozone mass mixing ratio [kg/kg]
       mr_lsliq, &   ! mixing_ratio_large_scale_cloud_liquid [kg/kg]
       mr_lsice, &   ! mixing_ratio_large_scale_cloud_ice [kg/kg]
       mr_ccliq, &   ! mixing_ratio_convective_cloud_liquid [kg/kg]
       mr_ccice, &   ! mixing_ratio_convective_cloud_ice [kg/kg]
#ifdef OPT_CHIMERRA
       mr_lsrain, &  ! mixing_ratio_large_scale_rainwater [kg/kg]
       mr_lssnow, &  ! mixing_ratio_large_scale_snowwater [kg/kg]
#ifdef OPT_GRPL
       mr_lsgrpl, &  ! mixing_ratio_large_scale_grplwater [kg/kg]
#endif
#endif
       fl_lsrain, &  ! flux_large_scale_cloud_rain [kg m^-2 s^-1]
       fl_lssnow, &  ! flux_large_scale_cloud_snow [kg m^-2 s^-1]
       fl_lsgrpl, &  ! flux_large_scale_cloud_graupel [kg m^-2 s^-1]
       fl_ccrain, &  ! flux_convective_cloud_rain [kg m^-2 s^-1]
       fl_ccsnow, &  ! flux_convective_cloud_snow [kg m^-2 s^-1]
       dtau_s, &     ! mean 0.67 micron optical depth of stratiform clouds
       dtau_c, &     ! mean 0.67 micron optical depth of convective clouds
       dem_s, &      ! 10.5 micron longwave emissivity of stratiform clouds
       dem_c         ! 10.5 micron longwave emissivity of convective clouds
  real(wp),dimension(Npoints,Nlevels+1),intent(inout) :: &
       ph, &         ! Pressure at half model levels [Pa]
       zlev_half     ! Height at half model levels [m]
  real(wp),dimension(Npoints,Nlevels,N_HYDRO),intent(in) :: &
       Reff          ! effective radius of each hydrometeor [m]
#ifdef OPT_DPLRW
  real(wp),dimension(Npoints,Nlevels),intent(inout) :: &
       gwvel           ! w velocity m/s
  real(wp),dimension(Npoints,Nlevels+1),intent(inout) :: &
       gcumf           ! convective vertical velocity from cumulus scheme
#endif
  logical,intent(in) :: &
       ocospchk, &   ! Check COSP input variables?
       ocospudf      ! Change undefined variables?
  character(len=128),intent(in)  :: cosp_input_nl  ! fullpath for input namelist
  character(len=128),intent(in)  :: cosp_output_nl ! fullpath for output namelist
  !character(len=256),intent(in)  :: cosp_input_nl  ! fullpath for input namelist
  !character(len=256),intent(in)  :: cosp_output_nl ! fullpath for output namelist

  ! Test data
  integer :: &
       Nlon,Nlat,geomode
  real(wp) :: &
       emsfc_lw
  real(wp),dimension(Npoints)       :: &
       u_wind,v_wind
!!!  real(wp),dimension(:),allocatable,target:: &
!!!       lon,       & ! Longitude (deg)
!!!       lat,       & ! Latitude (deg)
!!!       skt,       & ! Skin temperature (K)
!!!       landmask,  & ! Land/sea mask (0/1)
!!!       u_wind,    & ! U-component of wind (m/s)
!!!       v_wind,    & ! V-component of wind (m/s)
!!!       sunlit       ! Sunlit flag
!!!  real(wp),dimension(:,:),allocatable,target :: &
!!!       p,         & ! Model pressure levels (pa)
!!!       ph,        & ! Moddel pressure @ half levels (pa)
!!!       zlev,      & ! Model level height (m)
!!!       zlev_half, & ! Model level height @ half-levels (m)
!!!       T,         & ! Temperature (K)
!!!       sh,        & ! Specific humidity (kg/kg)
!!!       rh,        & ! Relative humidity (1)
!!!       tca,       & ! Total cloud fraction (1)
!!!       cca,       & ! Convective cloud fraction (1) 
!!!       mr_lsliq,  & ! Mass mixing ratio for stratiform cloud liquid (kg/kg)
!!!       mr_lsice,  & ! Mass mixing ratio for stratiform cloud ice (kg/kg)
!!!       mr_ccliq,  & ! Mass mixing ratio for convective cloud liquid (kg/kg)
!!!       mr_ccice,  & ! Mass mixing ratio for convective cloud ice (kg/kg)
!!!       mr_ozone,  & ! Mass mixing ratio for ozone (kg/kg)
!!!       fl_lsrain, & ! Precipitation flux (rain) for stratiform cloud (kg/m^2/s)
!!!       fl_lssnow, & ! Precipitation flux (snow) for stratiform cloud (kg/m^2/s)
!!!       fl_lsgrpl, & ! Precipitation flux (groupel) for stratiform cloud (kg/m^2/s)
!!!       fl_ccrain, & ! Precipitation flux (rain) for convective cloud (kg/m^2/s)
!!!       fl_ccsnow, & ! Precipitation flux (snow) for convective cloud (kg/m^2/s)
!!!       dtau_s,    & ! 0.67micron optical depth (stratiform cloud) (1)
!!!       dtau_c,    & ! 0.67micron optical depth (convective cloud) (1)
!!!       dem_s,     & ! 11micron emissivity (stratiform cloud) 
!!!       dem_c        ! 11microm emissivity (convective cloud)
  real(wp),dimension(:,:,:),allocatable,target :: &
          frac_out   !Subcolumn cloud cover (0/1)
!!!       frac_out ,  & ! Subcolumn cloud cover (0/1)
!!!       Reff         ! Subcolumn effective radius

#ifndef OPT_GRPL
  ! temporary array for lsgrpl added by YN
  real(wp),dimension(Npoints,Nlevels) :: mr_lsgrpl
#endif

  ! Input namelist fields
  integer ::                      & !
!!!       Npoints,                   & ! Number of gridpoints
       Ncolumns,                  & ! Number of subcolumns
!!!       Nlevels,                   & ! Number of model vertical levels
       Npoints_it,                & ! Number of gridpoints to be processed in one 
                                    ! iteration
       Nvgrid,                    & ! Number of vertical levels for statistical outputs 
                                    ! (USE_VGRID=.true.)
       surface_radar,             & ! surface=1/spaceborne=0
       cloudsat_use_gas_abs,      & ! Include gaseous absorption (1=yes/0=no)
       cloudsat_do_ray,           & ! Calculate output Rayleigh (1=yes/0=no)
       lidar_ice_type,            & ! Ice particle shape in lidar calculations 
                                    ! (0=ice-spheres/1=ice-non-spherical)
       overlap,                   & ! Overlap type: 1=max, 2=rand, 3=max/rand
       isccp_topheight,           & ! ISCCP cloud top height
       isccp_topheight_direction, & ! ISCCP cloud top height direction
       rttov_platform,            & ! RTTOV: Satellite platform
       rttov_satellite,           & ! RTTOV: Satellite
       rttov_instrument,          & ! RTTOV: Instrument
       rttov_Nchannels              ! RTTOV: Number of channels to be computed
  real(wp) ::                     & !
       cloudsat_radar_freq,       & ! CloudSat radar frequency (GHz)
       cloudsat_k2,               & ! |K|^2, -1=use frequency dependent default
       rttov_ZenAng,              & ! RTTOV: Satellite Zenith Angle
       co2,                       & ! CO2 mixing ratio
       ch4,                       & ! CH4 mixing ratio
       n2o,                       & ! n2o mixing ratio
       co                           ! co mixing ratio
  logical ::                      & !
       use_vgrid,                 & ! Use fixed vertical grid for outputs?
       csat_vgrid,                & ! CloudSat vertical grid? 
       use_precipitation_fluxes     ! True if precipitation fluxes are input to the 
                                    ! algorithm 
  integer,dimension(RTTOV_MAX_CHANNELS) :: &
       rttov_Channels               ! RTTOV: Channel numbers
  real(wp),dimension(RTTOV_MAX_CHANNELS) :: &
       rttov_Surfem                 ! RTTOV: Surface emissivity
  character(len=64) :: &
       cloudsat_micro_scheme        ! Microphysical scheme used in cloudsat radar simulator
  logical ofirst
  save    ofirst
  data    ofirst / .true. /
! for dbze output. Jing X.
  real(wp),dimension(Npoints,Nlevels):: &
       dbze94_tmp,&           ! CloudSat radar reflectivity
       dbze94grdmax2d,&       ! CloudSat grid-maximum radar refrectivity
       fracout_tmp            ! CloudSat radar reflectivity
  real(wp),dimension(Npoints):: &
       dbze94grdmax1d         ! CloudSat grid-maximum radar refrectivity
  character(3) :: idcol       ! index # of subcolumns

  integer :: jfpar, ifile, ieof 
  save ifile, jfpar
  character :: hclas*3
  data         hclas / 'SFC' /
  save         hclas  !! T.Michibata
  integer  :: ncls
  integer :: i, j, k, n, tp
  real(wp) :: flux_limit, mxratio_limit
  data       flux_limit    / 1.D-20 /
  data       mxratio_limit / 1.D-20 / 


!!!  character(len=64) :: &
!!!       finput                       ! Input NetCDF file
!!!  character(len=512) :: &
!!!       dinput                       ! Directory where the input files are located
!!!  character(len=600) :: &
!!!       fileIN                       ! dinput+finput
!!!  namelist/COSP_INPUT/overlap,isccp_topheight,isccp_topheight_direction,npoints,         &
!!!                      npoints_it,ncolumns,nlevels,use_vgrid,Nlvgrid,csat_vgrid,dinput,   &
!!!                      finput,cloudsat_radar_freq,surface_radar,cloudsat_use_gas_abs,     &
!!!                      cloudsat_do_ray,cloudsat_k2,cloudsat_micro_scheme,lidar_ice_type,  &
!!!                      use_precipitation_fluxes,rttov_platform,rttov_satellite,           &
!!!                      rttov_Instrument,rttov_Nchannels,rttov_Channels,rttov_Surfem,      &
!!!                      rttov_ZenAng,co2,ch4,n2o,co

  namelist/COSP_INPUT/overlap,isccp_topheight,isccp_topheight_direction,         &
                      npoints_it,ncolumns,use_vgrid,Nvgrid,csat_vgrid,   & ! rename Nlvgrid -> Nvgrid by YN
                      cloudsat_radar_freq,surface_radar,cloudsat_use_gas_abs,     &
                      cloudsat_do_ray,cloudsat_k2,cloudsat_micro_scheme,lidar_ice_type,  &
                      use_precipitation_fluxes,rttov_platform,rttov_satellite,           &
                      rttov_Instrument,rttov_Nchannels,rttov_Channels,rttov_Surfem,      &
                      rttov_ZenAng,co2,ch4,n2o,co
  !add save. Jing X.
  save  overlap,isccp_topheight,isccp_topheight_direction,         &
                      npoints_it,ncolumns,use_vgrid,Nvgrid,csat_vgrid,   & ! rename Nlvgrid -> Nvgrid by YN
                      cloudsat_radar_freq,surface_radar,cloudsat_use_gas_abs,     &
                      cloudsat_do_ray,cloudsat_k2,cloudsat_micro_scheme,lidar_ice_type,  &
                      use_precipitation_fluxes,rttov_platform,rttov_satellite,           &
                      rttov_Instrument,rttov_Nchannels,rttov_Channels,rttov_Surfem,      &
                      rttov_ZenAng,co2,ch4,n2o,co
   

  ! Output namelist
  logical :: Lcfaddbze94,Ldbze94,Latb532,LcfadLidarsr532,Lclcalipso,Lclhcalipso,         &
             Lcllcalipso,Lclmcalipso,Lcltcalipso,LparasolRefl,Lclcalipsoliq,             &
             Lclcalipsoice,Lclcalipsoun,Lclcalipsotmp,Lclcalipsotmpliq,Lclcalipsotmpice, &
             Lclcalipsotmpun,Lclhcalipsoliq,Lcllcalipsoliq,Lclmcalipsoliq,               &
             Lcltcalipsoliq,Lclhcalipsoice,Lcllcalipsoice,Lclmcalipsoice,Lcltcalipsoice, &
             Lclhcalipsoun,Lcllcalipsoun,Lclmcalipsoun,Lcltcalipsoun,Lalbisccp,          &
             Lboxptopisccp,Lboxtauisccp,Lpctisccp,Lclisccp,Ltauisccp,Lcltisccp,          &
             Lmeantbisccp,Lmeantbclrisccp,LclMISR,Lclcalipso2,Lcltlidarradar,Lfracout,   &
             LlidarBetaMol532,Lcltmodis,Lclwmodis,Lclimodis,Lclhmodis,Lclmmodis,         &
             Lcllmodis,Ltautmodis,Ltauwmodis,Ltauimodis,Ltautlogmodis,Ltauwlogmodis,     &
             Ltauilogmodis,Lreffclwmodis,Lreffclimodis,Lpctmodis,Llwpmodis,Liwpmodis,    &
             Lclmodis,Ltbrttov

  namelist/COSP_OUTPUT/Lcfaddbze94,Ldbze94,Latb532,LcfadLidarsr532,Lclcalipso,           &
                       Lclhcalipso,Lcllcalipso,Lclmcalipso,Lcltcalipso,LparasolRefl,     &
                       Lclcalipsoliq,Lclcalipsoice,Lclcalipsoun,Lclcalipsotmp,           &
                       Lclcalipsotmpliq,Lclcalipsotmpice,Lclcalipsotmpun,Lclhcalipsoliq, &
                       Lcllcalipsoliq,Lclmcalipsoliq,Lcltcalipsoliq,Lclhcalipsoice,      &
                       Lcllcalipsoice,Lclmcalipsoice,Lcltcalipsoice,Lclhcalipsoun,       &
                       Lcllcalipsoun,Lclmcalipsoun,Lcltcalipsoun,Lalbisccp,Lboxptopisccp,&
                       Lboxtauisccp,Lpctisccp,Lclisccp,Ltauisccp,Lcltisccp,Lmeantbisccp, &
                       Lmeantbclrisccp,LclMISR,Lclcalipso2,Lcltlidarradar,Lfracout,      &
                       LlidarBetaMol532,Lcltmodis,Lclwmodis,Lclimodis,Lclhmodis,         &
                       Lclmmodis,Lcllmodis,Ltautmodis,Ltauwmodis,Ltauimodis,             &
                       Ltautlogmodis,Ltauwlogmodis,Ltauilogmodis,Lreffclwmodis,          &
                       Lreffclimodis,Lpctmodis,Llwpmodis,Liwpmodis,Lclmodis,Ltbrttov

   save                Lcfaddbze94,Ldbze94,Latb532,LcfadLidarsr532,Lclcalipso,           &
                       Lclhcalipso,Lcllcalipso,Lclmcalipso,Lcltcalipso,LparasolRefl,     &
                       Lclcalipsoliq,Lclcalipsoice,Lclcalipsoun,Lclcalipsotmp,           &
                       Lclcalipsotmpliq,Lclcalipsotmpice,Lclcalipsotmpun,Lclhcalipsoliq, &
                       Lcllcalipsoliq,Lclmcalipsoliq,Lcltcalipsoliq,Lclhcalipsoice,      &
                       Lcllcalipsoice,Lclmcalipsoice,Lcltcalipsoice,Lclhcalipsoun,       &
                       Lcllcalipsoun,Lclmcalipsoun,Lcltcalipsoun,Lalbisccp,Lboxptopisccp,&
                       Lboxtauisccp,Lpctisccp,Lclisccp,Ltauisccp,Lcltisccp,Lmeantbisccp, &
                       Lmeantbclrisccp,LclMISR,Lclcalipso2,Lcltlidarradar,Lfracout,      &
                       !Lcloudsat_tcc, Lcloudsat_tcc2,                                    &
                       LlidarBetaMol532,Lcltmodis,Lclwmodis,Lclimodis,Lclhmodis,         &
                       Lclmmodis,Lcllmodis,Ltautmodis,Ltauwmodis,Ltauimodis,             &
                       Ltautlogmodis,Ltauwlogmodis,Ltauilogmodis,Lreffclwmodis,          &
                       Lreffclimodis,Lpctmodis,Llwpmodis,Liwpmodis,Lclmodis,Ltbrttov     

  logical ::                Lwr_occfreq, Lcfodd
  namelist / COSP_OUTPUT /  Lwr_occfreq, Lcfodd
  save                      Lwr_occfreq, Lcfodd

  logical ::                Lptradarflag0,Lptradarflag1,Lptradarflag2,Lptradarflag3,Lptradarflag4,           & !<-added
                            Lptradarflag5,Lptradarflag6,Lptradarflag7,Lptradarflag8,Lptradarflag9,Lradarpia    !<-added
  namelist / COSP_OUTPUT /  Lptradarflag0,Lptradarflag1,Lptradarflag2,Lptradarflag3,Lptradarflag4,           & !<-added
                            Lptradarflag5,Lptradarflag6,Lptradarflag7,Lptradarflag8,Lptradarflag9,Lradarpia    !<-added
  save                      Lptradarflag0,Lptradarflag1,Lptradarflag2,Lptradarflag3,Lptradarflag4,           & !<-added
                            Lptradarflag5,Lptradarflag6,Lptradarflag7,Lptradarflag8,Lptradarflag9,Lradarpia    !<-added

  ! scops check: added by YN
  logical ::                Lscops
  namelist / COSP_OUTPUT /  Lscops
  save                      Lscops

#ifdef OPT_DPLRW
  logical                :: Ldplrw!, Ldplrw_LS, Ldplrw_CU
  namelist /COSP_OUTPUT/    Ldplrw!, Ldplrw_LS, Ldplrw_CU
  save                      Ldplrw!, Ldplrw_LS, Ldplrw_CU
#endif

! Working variable
  real(wp),dimension(Npoints,Nlevels) :: &
       zlev2, &             ! Height at model levels [m]
!       zlev_half2, &        ! Height at half model levels [m] (Bottom of model layer)
       mr_lsliq2, &         ! mixing_ratio_large_scale_cloud_liquid [kg/kg]
       mr_lsice2, &         ! mixing_ratio_large_scale_cloud_ice [kg/kg]
       mr_ccliq2, &         ! mixing_ratio_convective_cloud_liquid [kg/kg]
       mr_ccice2, &         ! mixing_ratio_convective_cloud_ice [kg/kg]
       fl_lsrain2, &        ! flux_large_scale_cloud_rain [kg m^-2 s^-1]
       fl_lssnow2, &        ! flux_large_scale_cloud_snow [kg m^-2 s^-1]
       fl_lsgrpl2, &        ! flux_large_scale_cloud_graupel [kg m^-2 s^-1]
       fl_ccrain2, &        ! flux_convective_cloud_rain [kg m^-2 s^-1]
       fl_ccsnow2           ! flux_convective_cloud_snow [kg m^-2 s^-1]  
  real(wp),dimension(Npoints,Nlevels+1) :: &
       zlev_half2           ! Height at half model levels [m]

  ! Local variables
  logical :: &
       lsingle   = .true.,  & ! True if using MMF_v3_single_moment CLOUDSAT microphysical scheme (default)
       ldouble   = .false., & ! True if using MMF_v3.5_two_moment CLOUDSAT microphysical scheme
       lisccp    = .false. ,& ! Local on/off switch for simulators (used by initialization)
       lmodis    = .false., & !
       lmisr     = .false., & !
       lcalipso  = .false., & !
       lcloudsat = .false., & !
       lrttov    = .false., & !
       lparasol  = .false.    !
  type(size_distribution) :: &
       sd                ! Hydrometeor description
  type(radar_cfg) :: &
       rcfg_cloudsat     ! Radar configuration
  type(cosp_outputs) :: &
       cospOUT           ! COSP simulator outputs
  type(cosp_optical_inputs) :: &
       cospIN            ! COSP optical (or derived?) fields needed by simulators
  type(cosp_column_inputs) :: &
       cospstateIN       ! COSP model fields needed by simulators
  integer :: iChunk,nChunks,start_idx,end_idx,nPtsPerIt
  real(wp),dimension(10) :: driver_time

  character(len=256),dimension(100) :: cosp_status
#ifdef OPT_DPLRW
  integer,parameter :: &
       I_LSCLIQ = 1, & ! Large-scale (stratiform) liquid
       I_LSCICE = 2, & ! Large-scale (stratiform) ice
       I_LSRAIN = 3, & ! Large-scale (stratiform) rain
       I_LSSNOW = 4, & ! Large-scale (stratiform) snow
       I_CVCLIQ = 5, & ! Convective liquid
       I_CVCICE = 6, & ! Convective ice
       I_CVRAIN = 7, & ! Convective rain
       I_CVSNOW = 8, & ! Convective snow
       I_LSGRPL = 9    ! Large-scale (stratiform) groupel
  
  ! ! Stratiform and convective clouds in frac_out (scops output).
  ! integer, parameter :: &
  !      I_LSC = 1, & ! Large-scale clouds
  !      I_CVC = 2    ! Convective clouds   
  
#endif
  ! Microphysical settings for the precipitation flux to mixing ratio conversion
  ! real(wp),parameter,dimension(N_HYDRO) :: &
  !                ! LSL   LSI      LSR       LSS   CVL  CVI      CVR       CVS       LSG
  !      N_ax    = (/-1., -1.,     8.e6,     3.e6, -1., -1.,     8.e6,     3.e6,     4.e6/),&
  !      N_bx    = (/-1., -1.,      0.0,      0.0, -1., -1.,      0.0,      0.0,      0.0/),&
  !      alpha_x = (/-1., -1.,      0.0,      0.0, -1., -1.,      0.0,      0.0,      0.0/),&
  !      c_x     = (/-1., -1.,    842.0,     4.84, -1., -1.,    842.0,     4.84,     94.5/),&
  !      d_x     = (/-1., -1.,      0.8,     0.25, -1., -1.,      0.8,     0.25,      0.5/),&
  !      g_x     = (/-1., -1.,      0.5,      0.5, -1., -1.,      0.5,      0.5,      0.5/),&
  !      a_x     = (/-1., -1.,    524.0,    52.36, -1., -1.,    524.0,    52.36,   209.44/),&
  !      b_x     = (/-1., -1.,      3.0,      3.0, -1., -1.,      3.0,      3.0,      3.0/),&
  !      gamma_1 = (/-1., -1., 17.83725, 8.284701, -1., -1., 17.83725, 8.284701, 11.63230/),&
  !      gamma_2 = (/-1., -1.,      6.0,      6.0, -1., -1.,      6.0,      6.0,      6.0/),&
  !      gamma_3 = (/-1., -1.,      2.0,      2.0, -1., -1.,      2.0,      2.0,      2.0/),&
  !      gamma_4 = (/-1., -1.,      6.0,      6.0, -1., -1.,      6.0,      6.0,      6.0/)       

  ! Fields used solely for output
  integer,parameter :: &
       n_out_list = 100,           & ! Number of possible output variables
       N3D        = 8,            & ! Number of 3D output variables
       N2D        = 14,           & ! Number of 2D output variables
       N1D        = 41              ! Number of 1D output variables  
  character(len=32),dimension(n_out_list) :: out_list  ! List of output variable names
  integer :: lon_axid,time_axid,height_axid,height_mlev_axid,grid_id,lonvar_id,       &
             latvar_id,column_axid,sza_axid,temp_axid,channel_axid,dbze_axid,sratio_axid,&
             MISR_CTH_axid,lat_axid,tau_axid,pressure2_axid 
!  type(var1d) :: v1d(N1D+1) ! Structures needed by output routines for 1D variables
!  type(var2d) :: v2d(N2D)   ! Structures needed by output routines for 2D variables
!  type(var3d) :: v3d(N3D)   ! Structures needed by output routines for 3D variables
  real(wp)    :: time     
  !double precision :: time,time_bnds(2),time_step,half_time_step

!!!  real(wp),dimension(Npoints)  :: &
!!!       clwmodis, &    ! MODIS liquid cloud fraction
!!!       climodis, &    ! MODIS ice cloud fraction
!!!       reffclwmodis, & ! MODIS liquid cloud effective radius
!!!       reffclimodis, & ! MODIS ice cloud effective radius
!!!       clwxreffclw, & ! clwmodis * reffclwmodis, for estimating liq effective radius
!!!       clixreffcli  ! climodis * reffclimodis, for estimating ice effective radius
!!!       cltcalipso, & ! lidar total cloud fraction
!!!       cltisccp             ! total cloud fraction by the ISCCP simulator
  real(wp),dimension(Npoints) :: &
       cllcalipso, & ! lidar low-level cloud fraction
       clmcalipso, & ! lidar mid-level cloud fraction
       clhcalipso, & ! lidar high-level cloud fraction
       cltcalipso, & ! lidar total cloud fraction
       cltlidarradar, & ! lidar and radar total cloud fraction
       cllcalipsoice, & ! lidar low-level ice cloud fraction
       clmcalipsoice, & ! lidar mid-level ice cloud fraction
       clhcalipsoice, & ! lidar high-level ice cloud fraction
       cltcalipsoice, & ! lidar total ice cloud fraction
       cllcalipsoliq, & ! lidar low-level liq cloud fraction
       clmcalipsoliq, & ! lidar mid-level liq cloud fraction
       clhcalipsoliq, & ! lidar high-level liq cloud fraction
       cltcalipsoliq, & ! lidar total liq cloud fraction
       cllcalipsoun, &  ! lidar low-level undefined-phase cloud fraction
       clmcalipsoun, &  ! lidar mid-level undefined-phase cloud fraction
       clhcalipsoun, &  ! lidar high-level undefined-phase cloud fraction
       cltcalipsoun, &  ! lidar total undefined-phase cloud fraction
       cltisccp, &   ! total cloud fraction by the ISCCP simulator
       pctisccp, &   ! mean cloud top pressure by the ISCCP simulator
       tauisccp, &   ! mean optical depth by the ISCCP simulator
       albisccp, &   ! mean cloud albedo by the ISCCP simulator
       meantbisccp, & ! mean all-sky 10.5 micron brightness temperature by the ISCCP simulator
       meantbclrisccp, & ! mean clear-sky 10.5 micron brightness temperature by the ISCCP simulator
       cltmodis, &    ! MODIS total cloud fraction
       clwmodis, &    ! MODIS liquid cloud fraction
       climodis, &    ! MODIS ice cloud fraction
       clhmodis, &    ! MODIS high-level cloud fraction
       clmmodis, &    ! MODIS mid-level cloud fraction
       cllmodis, &    ! MODIS low-level cloud fraction
       tautmodis, &   ! MODIS total cloud optical thickness
       tauwmodis, &   ! MODIS liquid cloud optical thickness
       tauimodis, &   ! MODIS ice cloud optical thickness
       tautlogmodis, & ! MODIS total cloud optical thickness (Log10 Mean)
       tauwlogmodis, & ! MODIS liquid cloud optical thickness (Log10 Mean)
       tauilogmodis, & ! MODIS ice cloud optical thickness (Log10 Mean)
       reffclwmodis, & ! MODIS liquid cloud effective radius
       reffclimodis, & ! MODIS ice cloud effective radius
       clwxreffclw, & ! clwmodis * reffclwmodis, for estimating liq effective radius
       clixreffcli, &  ! climodis * reffclimodis, for estimating ice effective radius
       pctmodis, &    ! MODIS cloud top pressure
       lwpmodis, &    ! MODIS cloud liquid water path
       iwpmodis       ! MODIS cloud ice water path

  real(wp),dimension(Npoints,7,7) :: &
       clisccp, &     ! cloud fraction as calculated by the ISCCP simulator
       clmodis        ! cloud fraction as calculated by the MODIS simulator
  real(8),dimension(Npoints,7,numMISRHgtBins) :: &
       clmisr         ! cloud fraction as calculated by the MISR simulator
   !
   ! Allocatable output variables
   real(wp),dimension(:,:),allocatable :: &
        clcalipso, &    ! CALIPSO cloud area fraction
        clcalipso2, &   ! CALIPSO cloud area fraction undetected by CloudSat
        clcalipsoice, & ! CALIPSO ice cloud area fraction 
        clcalipsoliq, & ! CALIPSO liq cloud area fraction 
        clcalipsoun,  & ! CALIPSO undefined-phase cloud area fraction 
        clcalipsotmp, & ! CALIPSO cloud area fraction wrt temperature
        clcalipsotmpice, & ! CALIPSO ice cloud area fraction wrt temperature
        clcalipsotmpliq, & ! CALIPSO liq cloud area fraction wrt temperature
        clcalipsotmpun,  & ! CALIPSO undefined-phase cloud area fraction wrt temperature
        lidarbetamol532, & ! Lidar Molecular Backscatter (532 nm)
        boxtauisccp, &  ! Optical Depth in Each Column calculated by ISCCP Simulator
        boxptopisccp, & ! Cloud Top Pressure in Each Column calculated by ISCCP Simulator
        parasolGrid_Refl     ! PARASOL-like mono-directional reflectance (column)
   real(wp),dimension(:,:,:),allocatable :: &
        cfadlidarsr532, &  ! Histogram of CALIPSO scattering ratio (CFAD)
        cfaddbze94, &      ! Histogram of CloudSat radar reflectivity (CFAD)
        dbze94, &          ! CloudSat radar reflectivity
        atb532, &          ! Lidar attenuated total backscatter (532 nm)
        fracout, &           ! Subcolumn output from SCOPS
        parasolPix_Refl    !PARASOL-like mono-directional reflectance (subcolumn)
!!!  real(wp),dimension(:),allocatable :: mgrid_z,mgrid_zu,mgrid_zl

  !! Inline Warmrain Diagnostics:
  real(wp),dimension(:),allocatable :: &
       idid_stats_1d  ! IDiD summary output 1D
  real(wp),dimension(:,:),allocatable :: &
       idid_stats_2d  ! IDiD summary output 2D
  real(wp),dimension(:,:,:),allocatable :: &
       idid_stats_3d  ! IDiD summary output 3D
  logical,dimension(:),allocatable :: &
       mm1
  logical,dimension(:,:),allocatable :: &
       mm2

  !! precip flag !<-added
  real(wp),dimension(:),allocatable :: &
       ptcloudsatflag0,  &
       ptcloudsatflag1,  &
       ptcloudsatflag2,  &
       ptcloudsatflag3,  &
       ptcloudsatflag4,  &
       ptcloudsatflag5,  &
       ptcloudsatflag6,  &
       ptcloudsatflag7,  &
       ptcloudsatflag8,  &
       ptcloudsatflag9,  &
       cloudsatpia 

  !! scops check: added by YN
  real(wp),dimension(:,:,:),allocatable :: &
       out_array
  integer :: id
  character(1) :: c1
 
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!! debug start
  !call errcnt(0)
!  call getjfp(jfpar)
!  call flush(jfpar)
!!! debug end

  call cpu_time(driver_time(1))
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Very small values of COSP input variables reset to zero
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  do k=1,Nlevels
    do i=1,Npoints
       if ( abs(fl_lsrain( i,k ) ) .lt. flux_limit ) then
          fl_lsrain2( i,k ) = 0.D0
       else
          fl_lsrain2( i,k ) = fl_lsrain( i,k )
       endif
       if ( abs(fl_lssnow( i,k ) ) .lt. flux_limit ) then
          fl_lssnow2( i,k ) = 0.D0
       else
          fl_lssnow2( i,k ) = fl_lssnow( i,k )
       endif
       if ( abs(fl_lsgrpl( i,k ) ) .lt. flux_limit ) then
          fl_lsgrpl2( i,k ) = 0.D0
       else
          fl_lsgrpl2( i,k ) = fl_lsgrpl( i,k )
       endif
       if ( abs(fl_ccrain( i,k ) ) .lt. flux_limit ) then
          fl_ccrain2( i,k ) = 0.D0
       else
          fl_ccrain2( i,k ) = fl_ccrain( i,k )
       endif
       if ( abs(fl_ccsnow( i,k ) ) .lt. flux_limit ) then
          fl_ccsnow2( i,k ) = 0.D0
       else
          fl_ccsnow2( i,k ) = fl_ccsnow( i,k )
       endif
       if ( abs(mr_lsliq( i,k ) ) .lt. mxratio_limit ) then
          mr_lsliq2( i,k ) = 0.D0
       else
          mr_lsliq2( i,k ) = mr_lsliq( i,k )
       endif
       if ( abs(mr_lsice( i,k ) ) .lt. mxratio_limit ) then
          mr_lsice2( i,k ) = 0.D0
       else
          mr_lsice2( i,k ) = mr_lsice( i,k )
       endif
       if ( abs(mr_ccliq( i,k ) ) .lt. mxratio_limit ) then
          mr_ccliq2( i,k ) = 0.D0
       else
          mr_ccliq2( i,k ) = mr_ccliq( i,k )
       endif
       if ( abs(mr_ccice( i,k ) ) .lt. mxratio_limit ) then
          mr_ccice2( i,k ) = 0.D0
       else
          mr_ccice2( i,k ) = mr_ccice( i,k )
       endif
    enddo
  enddo
#ifndef OPT_GRPL
  ! temporary array for lsgrpl added by YN
  mr_lsgrpl = 0._wp
#endif

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Check COSP input variables
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  call histin( lon, 'ALONCHK', 'Longitude', 'deg east', 'ASFC', hclas)
  call histin( lat, 'ALATCHK', 'Latitude',  'deg north', 'ASFC',hclas )
  call histin( p, 'GDPCHK', 'pressure at full level', 'Pa', 'ALEV', &
               hclas )
  call histin( ph, 'GDPMCHK', 'pressure at half level', 'Pa', 'AMLV', &
               hclas )
  call histin( zlev, 'GDZCHK', 'height at full level', 'm', 'ALEV', &
               hclas )
  call histin( zlev_half, 'GDZMCHK', 'height at half level',  'm', 'AMLV', &
               hclas )
  call histin( T, 'GDTCHK', 'temp at full level', 'K',   'ALEV', &
               hclas )
  call histin( rh, 'RHCHK', 'relative humidity', ' ', 'ALEV', &
               hclas )
  call histin( sh, 'SHCHK', 'specific humidity', ' ', 'ALEV', &
               hclas )
  call histin( cca, 'CCNVCHK', 'convective cloud fraction', ' ', 'ALEV', &
               hclas )
  call histin( tca, 'CTTLCHK', 'LSC+CNV cloud fraction', ' ', 'ALEV', &
               hclas )
  call histin( skt, 'GRTSCHK', 'skin temperature', 'K', 'ASFC', hclas)
  call histin( surfelev, 'SURFELEVCHK', 'surface temperature', 'K', 'ASFC', &
               hclas )
  call histin( landmask, 'MASKCHK', 'land mask (0-ocen)', '', 'ASFC', &
               hclas )
  call histin( mr_ozone, 'GDO3CHK', 'ozone mixing ratio', 'kg/kg', &
               'ALEV', hclas )
  call histin( sunlit, 'SUNLITCHK', 'solar insl. index', ' ', &
               'ASFC', hclas )
  call histin( mr_lsliq2, 'QCLIQLCHK', 'LSC liq mixing ratio', 'kg/kg', &
               'ALEV', hclas )
  call histin( mr_lsice2, 'QCICELCHK', 'LSC ice mixing ratio', 'kg/kg', &
               'ALEV', hclas )
#ifdef OPT_CHIMERRA
  call histin( mr_lsrain, 'QRAINLCHK', 'LSC rain mixing ratio', 'kg/kg', &
               'ALEV', hclas )
  call histin( mr_lssnow, 'QSNOWLCHK', 'LSC snow mixing ratio', 'kg/kg', &
               'ALEV', hclas )
  call histin( mr_lsgrpl, 'QGRPLLCHK', 'LSC snow mixing ratio', 'kg/kg', &
               'ALEV', hclas )
#endif
  call histin( mr_ccliq2, 'QCLIQCCHK', 'CNV liq mixing ratio', 'kg/kg', &
               'ALEV', hclas )
  call histin( mr_ccice2, 'QCICECCHK', 'CNV ice mixing ratio', 'kg/kg', &
               'ALEV', hclas )
  call histin( fl_lsrain2, 'FRANLFCHK', 'LSC rain flux', 'kg/m**2/s',  &
               'ALEV', hclas )
  call histin( fl_lssnow2, 'FSNWLFCHK', 'LSC snow flux', 'kg/m**2/s', &
               'ALEV', hclas )
  call histin( fl_ccrain2, 'FRANCFCHK', 'CNV rain flux', 'kg/m**2/s', &
               'ALEV', hclas )
  call histin( fl_ccsnow2, 'FSNWCFCHK', 'CNV snow flux', 'kg/m**2/s', &
               'ALEV', hclas )
  CALL HISTIN( Reff(:,:,1), 'REFF1CHK',  'eff radius of LS liq', &
        'm', 'ALEV', HCLAS )
  CALL HISTIN( Reff(:,:,2), 'REFF2CHK',  'eff radius of LS ice', &
        'm', 'ALEV', HCLAS )
  CALL HISTIN( Reff(:,:,3), 'REFF3CHK',  'eff radius of LS rain', &
        'm', 'ALEV', HCLAS )
  CALL HISTIN( Reff(:,:,4), 'REFF4CHK',  'eff radius of LS snow', &
        'm', 'ALEV', HCLAS )
  CALL HISTIN( Reff(:,:,5), 'REFF5CHK',  'eff radius of CV liq', &
        'm', 'ALEV', HCLAS )
  CALL HISTIN( Reff(:,:,6), 'REFF6CHK',  'eff radius of CV ice', &
        'm', 'ALEV', HCLAS )
  CALL HISTIN( Reff(:,:,7), 'REFF7CHK',  'eff radius of CV rain', &
        'm', 'ALEV', HCLAS )
  CALL HISTIN( Reff(:,:,8), 'REFF8CHK',  'eff radius of CV snow', &
        'm', 'ALEV', HCLAS )
  CALL HISTIN( Reff(:,:,9), 'REFF9CHK',  'eff radius of LS grpl', &
        'm', 'ALEV', HCLAS )
  CALL HISTIN( dtau_s, 'DTAUSCHK','0.67 micron OD strati', &
        '',     'ALEV', HCLAS )
  CALL HISTIN( dtau_c, 'DTAUCCHK','0.67 micron OD convec', &
        '',     'ALEV', HCLAS )
  CALL HISTIN( dem_s, 'DEMSCHK','10.5 micron emis strati', &
        '',     'ALEV', HCLAS )
  CALL HISTIN( dem_c, 'DEMCCHK','10.5 micron emis convec', &
        '',     'ALEV', HCLAS )
#ifdef OPT_DPLRW
  CALL HISTIN( gwvel, 'GDWVELCHK','grid vertical air morion', 'm/s', &
       & 'ALEV', HCLAS )
  CALL HISTIN( gcumf, 'GDCUMFCHK','grid cumulus mass flux', 'kg/m^2/s', &
       & 'AMLV', HCLAS )
#endif


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Read in namelist
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Input namelist (cosp setup)
!!!  open(10,file=cosp_input_namelist,status='unknown')
!!!  read(10,nml=cosp_input)
!!!  close(10)
!!!
!!!  ! Output namelist (logical flags to turn on/off outputs)
!!!  open(10,file=cosp_output_namelist,status='unknown')
!!!  read(10,nml=cosp_output)
!!!  close(10)

  call getjfp( jfpar )
  if ( ofirst ) then
     ofirst = .false.
    ! ifile = 38 !<-added
    !!call getjfp( jfpar )
     write( jfpar,* ) ' @@@ COSP v2.0.0: 2017/11/22 @@@'
     write( jfpar,* ) ' @@@ COSP IDiD (warm rain): 2018/10/29 @@@'
#ifdef OPT_DPLRW
     write( jfpar,* ) ' @@@ Doppler vertical velocity Developer 2022 AUG @@@'
#endif
    !input namelist
     write( jfpar,* ) ' @@@ #1 IFLOPN for COSP_INPUT_NL @@@'
     call iflopn( ifile, ieof, &
                  cosp_input_nl, 'Read', '(*)' )
     write( jfpar,* ) ifile, ieof, cosp_input_nl
     write( jfpar,* ) ' @@@ #1 COMPLETED! @@@'
     if ( ieof .ne. 0 ) then
        write( jfpar,* ) ' ### FILE NOT FOUND : ', 'cosp_input_nl.txt'
        call xabort( 1 )
     endif
     write( jfpar,* ) ' @@@ #2 REWIND COSP_INPUT_NL @@@'
     rewind( ifile )
     write( jfpar,* ) ' @@@ #2 COMPLETED ! @@@'
     write( jfpar,* ) ' @@@ #3 READ COSP_INPUT_NL @@@'
     read( ifile, nml=cosp_input )
     write( jfpar,* ) ' @@@ #3 COMPLETED ! @@@'
     write( jfpar,* ) ' @@@ #4 WRITE COSP_INPUT_NL @@@'
     write( jfpar, cosp_input)
     write( jfpar,* ) ' @@@ #4 COMPLETED ! @@@'
     close( ifile )
    !output namelist 
     write( jfpar,* ) ' @@@ #5 IFLOPN for COSP_OUTPUT_NL @@@'
     call iflopn( ifile, ieof, &
                  cosp_output_nl, 'read', '(*)' )
     write( jfpar,* ) ' @@@ #5 COMPLETED ! @@@'
     if ( ieof .ne. 0 ) then
        write( jfpar,* ) ' ### FILE NOT FOUND : ', 'cosp_input_nl.txt'
        call xabort( 1 )
     endif
     rewind( ifile )
     read( ifile, nml=cosp_output )
     write( jfpar, cosp_output)
     close( ifile )
     !call read_cosp_output_nl(cosp_output_nl,N_OUT_LIST,cfg)
     call flush( jfpar )
  endif

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Read in sample input data.
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !!!allocate(lon(Npoints),lat(Npoints),p(Npoints,Nlevels),ph(Npoints,Nlevels),             &
  !!!         zlev(Npoints,Nlevels),zlev_half(Npoints,Nlevels),T(Npoints,Nlevels),          &
  !!!         sh(Npoints,Nlevels),rh(Npoints,Nlevels),tca(Npoints,Nlevels),                 &
  !!!         cca(Npoints,Nlevels),mr_lsliq(Npoints,Nlevels),mr_lsice(Npoints,Nlevels),     &
  !!!         mr_ccliq(Npoints,Nlevels),mr_ccice(Npoints,Nlevels),                          &
  !!!         fl_lsrain(Npoints,Nlevels),fl_lssnow(Npoints,Nlevels),                        &
  !!!         fl_lsgrpl(Npoints,Nlevels),fl_ccrain(Npoints,Nlevels),                        &
  !!!         fl_ccsnow(Npoints,Nlevels),Reff(Npoints,Nlevels,N_HYDRO),                     &
  !!!         dtau_s(Npoints,Nlevels),dtau_c(Npoints,Nlevels),dem_s(Npoints,Nlevels),       &
  !!!         dem_c(Npoints,Nlevels),skt(Npoints),landmask(Npoints),                        &
  !!!         mr_ozone(Npoints,Nlevels),u_wind(Npoints),v_wind(Npoints),sunlit(Npoints),    &
  !!!         frac_out(Npoints,Ncolumns,Nlevels))

!  allocate(frac_out(Npoints,Ncolumns,Nlevels))  ! T.Michibata removed (2021.06.08)
  do i=1,Npoints
     u_wind(i)=0._wp
     v_wind(i)=0._wp
  enddo

!  Time information
!  time           = 8*1._wp/8._wp  ! First time step
!  is this useful?
!###  cosptime       = time/86400.D0
!###  time_step      = 3._wp/24._wp
!###  half_time_step = 0.5_wp*time_step
!###  time_bnds      = (/cosptime-half_time_step,cosptime+half_time_step/)

 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 ! Modification of undef variables: -999.D0 --> -1.0E30
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  call undef_m2n(p, Npoints, Nlevels, 1)
  call undef_m2n(ph, Npoints, Nlevels+1, 1)
  call undef_m2n(zlev, Npoints, Nlevels, 1)
  do j=1,Npoints
     if ( zlev(j,1) .lt. 0.0 .and. zlev(j,1) .gt. -1.0E30+1.0 ) then
        do k=2,Nlevels
           zlev2(j,k) = zlev(j,k) - zlev(j,1) + 15.0
        enddo
        zlev2(j,1) = 15.0
     else
        do k=1,Nlevels
           zlev2(j,k) = zlev(j,k)
        enddo
     endif
  enddo
  call undef_m2n(zlev_half, Npoints, Nlevels+1, 1)
  do j=1,Npoints
     if ( zlev_half(j,1) .lt. 0.0 .and. zlev_half(j,1) .gt. -1.0E30+1.0 ) then
        do k=2,Nlevels+1
           zlev_half2(j,k) = zlev_half(j,k) - zlev_half(j,1)
        enddo
        zlev_half2(j,1) = 0.0
     else
        do k=1,Nlevels+1
           zlev_half2(j,k) = zlev_half(j,k)
        enddo
     endif
  enddo

#ifdef OPT_DPLRW
  where (cca > 0._wp .and. cca <= 1._wp)
     gcumf(:,1:Nlevels) = gcumf(:,1:Nlevels)/cca(:,:)
  elsewhere
     gcumf = -999._wp  ! MIROC_UNDEF
  end where
#endif

  call undef_m2n(T, Npoints, Nlevels, 1)
  call undef_m2n(rh, Npoints, Nlevels, 1)
  call undef_m2n(sh, Npoints, Nlevels, 1)
  call undef_m2n(cca, Npoints, Nlevels, 1)
  call undef_m2n(tca, Npoints, Nlevels, 1)
  call undef_m2n(skt, Npoints, 1, 1)
  call undef_m2n(surfelev, Npoints, 1, 1)
  call undef_m2n(landmask, Npoints, 1, 1)
  call undef_m2n(mr_ozone, Npoints, Nlevels, 1)
  call undef_m2n(sunlit, Npoints, 1, 1)
  call undef_m2n(mr_lsliq2, Npoints, Nlevels, 1)
  call undef_m2n(mr_lsice2, Npoints, Nlevels, 1)
  call undef_m2n(mr_ccliq2, Npoints, Nlevels, 1)
  call undef_m2n(mr_ccice2, Npoints, Nlevels, 1)
  call undef_m2n(fl_lsrain2, Npoints, Nlevels, 1)
  call undef_m2n(fl_lssnow2, Npoints, Nlevels, 1)
  call undef_m2n(fl_ccrain2, Npoints, Nlevels, 1)
  call undef_m2n(fl_ccsnow2, Npoints, Nlevels, 1)
  call undef_m2n(Reff, Npoints, Nlevels, 9)
  call undef_m2n(dtau_s, Npoints, Nlevels, 1)
  call undef_m2n(dtau_c, Npoints, Nlevels, 1)
  call undef_m2n(dem_s, Npoints, Nlevels, 1)
  call undef_m2n(dem_c, Npoints, Nlevels, 1)
#ifdef OPT_DPLRW
  call undef_m2n(gwvel, Npoints, Nlevels,   1)
  call undef_m2n(gcumf, Npoints, Nlevels+1, 1)
#endif

  !!!fileIN = trim(dinput)//trim(finput)
  !!!call nc_read_input_file(fileIN,Npoints,Nlevels,N_HYDRO,lon,lat,p,ph,zlev,zlev_half,    &
  !!!                        T,sh,rh,tca,cca,mr_lsliq,mr_lsice,mr_ccliq,mr_ccice,fl_lsrain, &
  !!!                        fl_lssnow,fl_lsgrpl,fl_ccrain,fl_ccsnow,Reff,dtau_s,dtau_c,    &
  !!!                        dem_s,dem_c,skt,landmask,mr_ozone,u_wind,v_wind,sunlit,        &
  !!!                        emsfc_lw,geomode,Nlon,Nlat)
  call cpu_time(driver_time(2))
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Which simulators need to be run? Look at which outputs are requested.
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (Lpctisccp .or. Lclisccp .or. Lboxptopisccp .or.  Lboxtauisccp .or. Ltauisccp .or. &
       Lcltisccp .or. Lmeantbisccp .or. Lmeantbclrisccp .or. Lalbisccp) Lisccp = .true.
  if (LclMISR) Lmisr = .true.
  if (Lcltmodis .or. Lclwmodis .or. Lclimodis .or. Lclhmodis .or. Lclmmodis .or.         &
       Lcllmodis .or. Ltautmodis .or. Ltauwmodis .or. Ltauimodis .or. Ltautlogmodis .or. &
       Ltauwlogmodis .or. Ltauilogmodis .or. Lreffclwmodis .or. Lreffclimodis .or.       &
       Lpctmodis .or. Llwpmodis .or. Liwpmodis .or. Lclmodis) Lmodis = .true.
  if (Lclcalipso2 .or. Lclcalipso .or.  Lclhcalipso .or. Lcllcalipso .or. Lclmcalipso    &
       .or. Lcltcalipso .or. Lcltlidarradar .or. Lclcalipsoliq .or. Lclcalipsoice .or.   &
       Lclcalipsoun .or. Lclcalipsotmp .or. Lclcalipsotmpliq .or. Lclcalipsotmpice .or.  &
       Lclcalipsotmpun .or. Lcltcalipsoliq .or. Lcltcalipsoice .or. Lcltcalipsoun .or.   &
       Lclhcalipsoliq .or. Lclhcalipsoice .or. Lclhcalipsoun .or. Lclmcalipsoliq .or.    &
       Lclmcalipsoice .or. Lclmcalipsoun .or. Lcllcalipsoliq .or. Lcllcalipsoice .or.    &
       Lcllcalipsoun .or. LlidarBetaMol532 .or. LcfadLidarsr532 .or. Lcltlidarradar .or. &
       Lcltlidarradar) lcalipso = .true.
  if (LcfadDbze94 .or. Ldbze94 .or. Lcltlidarradar .or. Lptradarflag0 .or. Lptradarflag1 & !<-added
       .or. Lptradarflag2 .or. Lptradarflag3 .or. Lptradarflag4 .or. Lptradarflag5 .or.  & !<-added
       Lptradarflag6 .or. Lptradarflag7 .or. Lptradarflag8 .or. Lptradarflag9 .or.       & !<-added
       Lradarpia) Lcloudsat = .true. !<-added
  if (Lparasolrefl) Lparasol = .true.
  if (Ltbrttov) Lrttov = .true.
#ifdef OPT_DPLRW
  !if (Ldplrw_LS .or. Ldplrw_CU) Ldplrw = .true.
  if (Ldplrw) Lcloudsat = .true.
#endif
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !                 IF IMPLEMTING COSP IN GCM, HERE IS WHERE TO START!!!
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Initialize COSP
  !*This only needs to be done the first time that COSP is called.*
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 if (.not. allocated(vgrid_zl) .or. .not. allocated(vgrid_zu) .or. .not. allocated(vgrid_z)) then
 !if (linitialization) then

  ! Initialize quickbeam_optics, also if two-moment radar microphysics scheme is wanted...
  if (cloudsat_micro_scheme == 'MMF_v3.5_two_moment')  then
     ldouble = .true. 
     lsingle = .false.
  endif  
  
  ! Initialize the distributional parameters for hydrometeors in radar simulator
  call hydro_class_init(lsingle,ldouble,sd)

  ! Initialize COSP simulator
  !call COSP_INIT(Lisccp, Lmodis, Lmisr, Lcloudsat, Lcalipso, Lparasol, Lrttov,           &
  !     Npoints, Nlevels, cloudsat_radar_freq, cloudsat_k2, cloudsat_use_gas_abs,         &
  !     cloudsat_do_ray, isccp_topheight, isccp_topheight_direction, surface_radar,       &
  !     rcfg_cloudsat, rttov_Nchannels, rttov_Channels, rttov_platform,                   &
  !     rttov_satellite, rttov_instrument, use_vgrid, csat_vgrid, Nlvgrid,                &
  !     cloudsat_micro_scheme, cospOUT)
  call COSP_INIT(Lisccp, Lmodis, Lmisr, Lcloudsat, Lcalipso, Lparasol, Lrttov,           &
       !Lwr_occfreq, Lcfodd,                                                              &
       Npoints, Nlevels, cloudsat_radar_freq, cloudsat_k2, cloudsat_use_gas_abs,         &
       cloudsat_do_ray, isccp_topheight, isccp_topheight_direction, surface_radar,       &
       rcfg_cloudsat, rttov_Nchannels, rttov_Channels, rttov_platform,                   &
       rttov_satellite, rttov_instrument, use_vgrid, csat_vgrid, Nvgrid,      &
!#ifdef OPT_DPLRW
!       Ldplrw, &
!#endif
       cloudsat_micro_scheme, cospOUT)

!#ifdef OPT_DPLRW
!  write(jfpar,*) '@@@@@ CFAD dims CHECK in COSP2 @@@@@'
!  write(jfpar,'(a,3(x,ES11.3e2),x,I4)') 'dplrLS: MIN,MAX,WID,Num = ',dplrLS_MIN,dplrLS_MAX,dplrLS_WID,NdplrLS
!  write(jfpar,'(a,3(x,ES11.3e2),x,I4)') 'dplrCV: MIN,MAX,WID,Num = ',dplrCU_MIN,dplrCU_MAX,dplrCU_WID,NdplrCU
!  write(jfpar,'(a,3(x,ES11.3e2),x,I4)') 'lvtemp: MIN,MAX,WID,Num = ',lvtemp_MIN,lvtemp_MAX,lvtemp_WID,Nlvtemp
!  write(jfpar,'(a,3(x,ES11.3e2),x,I4)') 'lvdBZe: MIN,MAX,WID,Num = ',lvdBZe_MIN,lvdBZe_MAX,lvdBZe_WID,NlvdBZe
!#endif
  call flush(jfpar)
  call cpu_time(driver_time(3))
  end if 
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Construct output derived type.
  ! *NOTE* The "construct/destroy" subroutines are local to this module and should be
  !        modified for your configuration. E.g. it may be overkill to query each field.
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  call construct_cosp_outputs(Lpctisccp, Lclisccp, Lboxptopisccp, Lboxtauisccp,          &
       Ltauisccp, Lcltisccp, Lmeantbisccp, Lmeantbclrisccp, Lalbisccp, LclMISR,          &
       Lcltmodis, Lclwmodis, Lclimodis, Lclhmodis, Lclmmodis, Lcllmodis, Ltautmodis,     &
       Ltauwmodis, Ltauimodis, Ltautlogmodis, Ltauwlogmodis, Ltauilogmodis,              &
       Lreffclwmodis, Lreffclimodis, Lpctmodis, Llwpmodis, Liwpmodis, Lclmodis, Latb532, &
       LlidarBetaMol532, LcfadLidarsr532, Lclcalipso2, Lclcalipso, Lclhcalipso,          &
       Lcllcalipso, Lclmcalipso, Lcltcalipso, Lcltlidarradar, Lclcalipsoliq,             &
       Lclcalipsoice, Lclcalipsoun, Lclcalipsotmp, Lclcalipsotmpliq, Lclcalipsotmpice,   &
       Lclcalipsotmpun, Lcltcalipsoliq, Lcltcalipsoice, Lcltcalipsoun, Lclhcalipsoliq,   &
       Lclhcalipsoice, Lclhcalipsoun, Lclmcalipsoliq, Lclmcalipsoice, Lclmcalipsoun,     &
       Lcllcalipsoliq, Lcllcalipsoice, Lcllcalipsoun, LcfadDbze94, Ldbze94,              &
       Lparasolrefl, Ltbrttov,                                                           &
       Lptradarflag0,Lptradarflag1,Lptradarflag2,Lptradarflag3,Lptradarflag4,            &
       Lptradarflag5,Lptradarflag6,Lptradarflag7,Lptradarflag8,Lptradarflag9,Lradarpia,  &
       Lwr_occfreq, Lcfodd,                                                              &
       Lscops,  & ! added by YN
#ifdef OPT_DPLRW
       Ldplrw, &!Ldplrw_LS, Ldplrw_CU, &
#endif
       Npoints, Ncolumns, Nlevels, rttov_Nchannels, cospOUT)

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Break COSP up into pieces and loop over each COSP 'chunk'.
  ! nChunks = # Points to Process (nPoints) / # Points per COSP iteration (nPoints_it)
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  nChunks = nPoints/nPoints_it+1

!  if (nPoints .eq. nPoints_it) nChunks = 1
  if (mod(nPoints,nPoints_it)==0) nChunks = nChunks - 1  ! modified by YN
  do iChunk=1,nChunks
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! Determine indices for "chunking" (again, if necessary)
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     if (nChunks .eq. 1) then
        start_idx = 1
        end_idx   = nPoints
        nPtsPerIt = nPoints
     else
        start_idx = (iChunk-1)*nPoints_it+1
        end_idx   = iChunk*nPoints_it
        if (end_idx .gt. nPoints) end_idx=nPoints
        nPtsPerIt = end_idx-start_idx+1
     endif

     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! Construct COSP input types
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     if (iChunk .eq. 1) then
        call construct_cospIN(Nptsperit,nColumns,nLevels, &
                              Lmodis, Lmisr, Lisccp, Lcalipso, Lcloudsat,cospIN)
        call construct_cospstateIN(Nptsperit,nLevels,rttov_nChannels,cospstateIN)
     endif
     if (iChunk .eq. nChunks) then
        call destroy_cospIN(cospIN)
        call destroy_cospstateIN(cospstateIN)
        call construct_cospIN(Nptsperit,nColumns,nLevels, &
                              Lmodis, Lmisr, Lisccp, Lcalipso, Lcloudsat,cospIN)
        call construct_cospstateIN(Nptsperit,nLevels,rttov_nChannels,cospstateIN)    
     endif
     call cpu_time(driver_time(4))
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! Populate input types with model fields.
     ! Here the 3D sample model fields (temperature,pressure,etc...) are ordered from the
     ! surface-2-TOA, whereas COSP expects all fields to be ordered from TOA-2-SFC. So the
     ! vertical fields are flipped prior to storing to COSP input type.
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     cospIN%emsfc_lw         = emsfc_lw
     cospIN%rcfg_cloudsat    = rcfg_cloudsat
     cospstateIN%hgt_matrix  = zlev2(start_idx:end_idx,Nlevels:1:-1) ! km
     cospstateIN%sunlit      = sunlit(start_idx:end_idx)            ! 0-1
     cospstateIN%skt         = skt(start_idx:end_idx)               ! K
     cospstateIN%surfelev    = surfelev(start_idx:end_idx)          ! m !<-added
     cospstateIN%land        = landmask(start_idx:end_idx)          ! 0-1 (*note* model specific)
     cospstateIN%qv          = sh(start_idx:end_idx,Nlevels:1:-1)   ! kg/kg
     cospstateIN%at          = T(start_idx:end_idx,Nlevels:1:-1)    ! K
     cospstateIN%pfull       = p(start_idx:end_idx,Nlevels:1:-1)    ! Pa 
     ! Pressure at interface (nlevels+1). Set uppermost interface to 0.
!     cospstateIN%phalf(:,2:Nlevels+1) = ph(start_idx:end_idx,Nlevels:1:-1)   ! Pa  
!     cospstateIN%phalf(:,1)           = 0._wp
     cospstateIN%phalf(:,1:Nlevels+1)  = ph(start_idx:end_idx,Nlevels+1:1:-1)
     ! Height at interface (nlevels+1). Set lowermost interface to 0.
     cospstateIN%hgt_matrix_half(:,1:Nlevels) = zlev_half2(start_idx:end_idx,Nlevels+1:2:-1) ! km
     cospstateIN%hgt_matrix_half(:,Nlevels+1) = 0._wp
     
     ! RTTOV inputs (by default, COSP is distributed not using RTTOV, but the infrastructure is there.)
     cospstateIN%u_sfc                = u_wind(start_idx:end_idx)            ! m/s
     cospstateIN%v_sfc                = v_wind(start_idx:end_idx)            ! m/s
     cospstateIN%emis_sfc             = rttov_surfem
     cospstateIN%zenang               = rttov_zenang
     cospstateIN%lat                  = lat(start_idx:end_idx)
     cospstateIN%lon                  = lon(start_idx:end_idx)
     cospstateIN%month                = 2 ! This is needed by RTTOV only for the surface emissivity calculation.
     cospstateIN%co2                  = co2*(amd/amCO2)*1e6
     cospstateIN%ch4                  = ch4*(amd/amCH4)*1e6  
     cospstateIN%n2o                  = n2o*(amd/amN2O)*1e6
     cospstateIN%co                   = co*(amd/amCO)*1e6
     cospstateIN%o3                   = mr_ozone(start_idx:end_idx,Nlevels:1:-1)*(amd/amO3)*1e6 ! microns

#ifdef OPT_DPLRW
     cospstateIN%gwvel  = gwvel(start_idx:end_idx,Nlevels:1:-1)
     cospstateIN%gcumf  = gcumf(start_idx:end_idx,Nlevels+1:1:-1)
#endif
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! Generate subcolumns and compute optical inputs.
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     call subsample_and_optics(nPtsPerIt,nLevels,nColumns,N_HYDRO,overlap,                     &
          use_precipitation_fluxes, Lcalipso, Lcloudsat, Lmodis, lidar_ice_type, sd,           &
          tca(start_idx:end_idx,Nlevels:1:-1),cca(start_idx:end_idx,Nlevels:1:-1),             &
          fl_lsrain2(start_idx:end_idx,Nlevels:1:-1),fl_lssnow2(start_idx:end_idx,Nlevels:1:-1), &
          fl_lsgrpl2(start_idx:end_idx,Nlevels:1:-1),fl_ccrain2(start_idx:end_idx,Nlevels:1:-1), &
          fl_ccsnow2(start_idx:end_idx,Nlevels:1:-1),mr_lsliq2(start_idx:end_idx,Nlevels:1:-1),  &
          mr_lsice2(start_idx:end_idx,Nlevels:1:-1),mr_ccliq2(start_idx:end_idx,Nlevels:1:-1),   &
          mr_ccice2(start_idx:end_idx,Nlevels:1:-1),Reff(start_idx:end_idx,Nlevels:1:-1,:),     &
          dtau_c(start_idx:end_idx,nLevels:1:-1),dtau_s(start_idx:end_idx,nLevels:1:-1),       &
          dem_c(start_idx:end_idx,nLevels:1:-1),dem_s(start_idx:end_idx,nLevels:1:-1),         &
#ifdef OPT_CHIMERRA
          mr_lsrain(start_idx:end_idx,Nlevels:1:-1),mr_lssnow(start_idx:end_idx,Nlevels:1:-1),  &
          mr_lsgrpl(start_idx:end_idx,Nlevels:1:-1), &
#endif
          cospstateIN,cospIN)
     call cpu_time(driver_time(6))

     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! Call COSP
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     cosp_status = COSP_SIMULATOR(cospIN, cospstateIN, cospOUT,start_idx,end_idx,.false.)
     call cpu_time(driver_time(7))
  enddo
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Output
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Create list of output varibles.
!  call construct_cospOutList(Lpctisccp,Lclisccp,Lboxptopisccp, Lboxtauisccp,Ltauisccp,   &
!                             Lcltisccp,Lmeantbisccp,Lmeantbclrisccp,Lalbisccp,LclMISR,   &
!                             Lcltmodis,Lclwmodis,Lclimodis,Lclhmodis,Lclmmodis,          &
!                             Lcllmodis,Ltautmodis,Ltauwmodis,Ltauimodis,Ltautlogmodis,   &
!                             Ltauwlogmodis,Ltauilogmodis,Lreffclwmodis,Lreffclimodis,    &
!                             Lpctmodis,Llwpmodis,Liwpmodis,Lclmodis,Latb532,             &
!                             LlidarBetaMol532,LcfadLidarsr532,Lclcalipso2,Lclcalipso,    &
!                             Lclhcalipso,Lcllcalipso,Lclmcalipso,Lcltcalipso,            &
!                             Lcltlidarradar,Lclcalipsoliq,Lclcalipsoice,Lclcalipsoun,    &
!                             Lclcalipsotmp,Lclcalipsotmpliq,Lclcalipsotmpice,            &
!                             Lclcalipsotmpun,Lcltcalipsoliq,Lcltcalipsoice,              &
!                             Lcltcalipsoun,Lclhcalipsoliq,Lclhcalipsoice,Lclhcalipsoun,  &
!                             Lclmcalipsoliq,Lclmcalipsoice,Lclmcalipsoun,Lcllcalipsoliq, &
!                             Lcllcalipsoice,Lcllcalipsoun,LcfadDbze94,Ldbze94,           &
!                             Lparasolrefl,Ltbrttov,N_OUT_LIST,out_list)  

  ! Time information needed by cmor output.
  !time           = 8*1._wp/8._wp
  !time_step      = 3._wp/24._wp
  !half_time_step = 0.5_wp*time_step
  !time_bnds      = (/time-half_time_step,time+half_time_step/)

  ! Model grid info for cmor output
!!!  allocate(mgrid_z(Nlevels),mgrid_zl(Nlevels),mgrid_zu(Nlevels))
!!!  mgrid_z  = zlev(1,:)
!!!  mgrid_zl = zlev_half(1,:)
!!!  mgrid_zu(1:Nlevels-1) = zlev_half(1,2:Nlevels)
!!!  mgrid_zu(Nlevels)     = zlev(1,Nlevels) + (zlev(1,Nlevels) - mgrid_zl(Nlevels))
  
!!!  if (geomode .eq. 1) then
!!!     call nc_cmor_init('../cmor/cosp_cmor_nl_1D.txt','replace',nPoints,nColumns,nLevels, &
!!!                       rttov_nChannels,nLvgrid_local,lon,lat,mgrid_zl,mgrid_zu,mgrid_z,cospOUT,&
!!!                       geomode,Nlon,Nlat,N1D+1,N2D,N3D,N_OUT_LIST,out_list,lon_axid,     &
!!!                       lat_axid,time_axid,height_axid,height_mlev_axid,grid_id,lonvar_id,&
!!!                       latvar_id,column_axid,sza_axid,temp_axid,channel_axid,dbze_axid,  &
!!!                       sratio_axid,MISR_CTH_axid,tau_axid,pressure2_axid,v1d(1:N1D+1),   &
!!!                       v2d,v3d)
!!!     call nc_cmor_associate_1d(grid_id,height_axid,height_mlev_axid,                     &
!!!                               column_axid,sza_axid,temp_axid,channel_axid,dbze_axid,    &
!!!                               sratio_axid,MISR_CTH_axid,tau_axid,pressure2_axid,        &
!!!                               nPoints,nColumns,nLevels,rttov_nChannels,nLvgrid_local,   &
!!!                               cospOUT,N1D+1,N2D,N3D,v1d,v2d,v3d)
!!!     call nc_cmor_write_1d(nPoints,lon,lat,time_bnds,lonvar_id,latvar_id,N1D+1,N2D,N3D,  &
!!!                           v1d(1:N1D+1),v2d,v3d)
!!!  endif
!!!  if (geomode .gt. 1) then
!!!     call nc_cmor_init('../cmor/cosp_cmor_nl_2D.txt','replace',nPoints,nColumns,nLevels, &
!!!                       rttov_nChannels,nLvgrid_local,lon,lat,mgrid_zl,mgrid_zu,mgrid_z,cospOUT,&
!!!                       geomode,Nlon,Nlat,N1D,N2D,N3D,N_OUT_LIST,out_list,lon_axid,       &
!!!                       lat_axid,time_axid,height_axid,height_mlev_axid,grid_id,lonvar_id,&
!!!                       latvar_id,column_axid,sza_axid,temp_axid,channel_axid,dbze_axid,  &
!!!                       sratio_axid,MISR_CTH_axid,tau_axid,pressure2_axid,v1d(1:N1D),v2d, &
!!!                       v3d)
!!!     call nc_cmor_associate_2d(lon_axid,lat_axid,time_axid,height_axid,height_mlev_axid, &
!!!                               column_axid,sza_axid,temp_axid,channel_axid,dbze_axid,    &
!!!                               sratio_axid,MISR_CTH_axid,tau_axid,pressure2_axid,Nlon,   &
!!!                               Nlat,nColumns,nLevels,rttov_nChannels,nLvgrid_local,      &
!!!                               cospOUT,N1D,N2D,N3D,v1d(1:N1D),v2d,v3d)
!!!     call nc_cmor_write_2d(time_bnds,geomode,Nlon,Nlat,N1D,N2D,N3D,v1d(1:N1D),v2d,v3d)
!!!  endif
!!!  deallocate(mgrid_z,mgrid_zu,mgrid_zl)
!!!  call nc_cmor_close()
!!!  print*,'Time to write to output:  ',driver_time(8)-driver_time(7)

     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
     ! Allocate output variables
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
!!!     allocate( cfadlidarsr532(Npoints,SR_BINS,Nlvgrid), cfaddbze94(Npoints,DBZE_BINS,Nlvgrid) )
     allocate( clcalipso(Npoints,Nlvgrid), clcalipso2(Npoints,Nlvgrid),              &
           clcalipsoice(Npoints,Nlvgrid), parasolGrid_Refl(Npoints,PARASOL_NREFL),   &
           parasolPix_Refl(Npoints,Ncolumns,PARASOL_NREFL),                          &
           clcalipsoliq(Npoints,Nlvgrid), clcalipsoun(Npoints,Nlvgrid),              &
           clcalipsotmp(Npoints,LIDAR_NTEMP), clcalipsotmpice(Npoints,LIDAR_NTEMP),  &
           clcalipsotmpliq(Npoints,LIDAR_NTEMP), clcalipsotmpun(Npoints,LIDAR_NTEMP),&
           lidarbetamol532(Npoints,Nlevels), boxtauisccp(Npoints,Ncolumns),          &
           boxptopisccp(Npoints,Ncolumns), cfadlidarsr532(Npoints,SR_BINS,Nlvgrid),  &
           cfaddbze94(Npoints,DBZE_BINS,Nlvgrid), dbze94(Npoints,Ncolumns,Nlevels),  &
           atb532(Npoints,Ncolumns,Nlevels), fracout(Npoints,Ncolumns,Nlevels) )   

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  ! HISTIN COSP variables
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  !
  ! 1D variables
  if (Lcllcalipso) then
     cllcalipso(:) = cospOUT%calipso_cldlayer(:,1)
     call undef_n2m( cllcalipso, Npoints, 1, 1 )
     call histin( cllcalipso, 'cllcalipso', &
          & 'Lidar Low-level Cloud Fraction', '%', 'ASFC', hclas )
  endif
  if (Lclmcalipso) then
     clmcalipso(:) = cospOUT%calipso_cldlayer(:,2)
     call undef_n2m( clmcalipso, Npoints, 1, 1 )
     call histin( clmcalipso, 'clmcalipso', &
          & 'Lidar Mid-level Cloud Fraction', '%', 'ASFC', hclas )
  endif
  if (Lclhcalipso) then
     clhcalipso(:) = cospOUT%calipso_cldlayer(:,3)
     call undef_n2m( clhcalipso, Npoints, 1, 1 )
     call histin( clhcalipso, 'clhcalipso', &
          & 'Lidar High-level Cloud Fraction', '%', 'ASFC', hclas )
  endif
  if (Lcltcalipso) then
     cltcalipso(:) = cospOUT%calipso_cldlayer(:,4)
     call undef_n2m( cltcalipso, Npoints, 1, 1 )
     call histin( cltcalipso, 'cltcalipso', &
          & 'Lidar Total Cloud Fraction', '%', 'ASFC', hclas )
  endif
  if (Lcltlidarradar) then
     cltlidarradar = cospOUT%radar_lidar_tcc
     call undef_n2m( cltlidarradar, Npoints, 1, 1 )
     call histin( cltlidarradar, 'cltlidarradar', &
          & 'Lidar and Radar Total Cloud Fraction', '%', 'ASFC', hclas )
  endif
  if (Lcllcalipsoice) then
     cllcalipsoice(:) = cospOUT%calipso_cldlayerphase(:,1,1)
     call undef_n2m( cllcalipsoice, Npoints, 1, 1 )
     call histin( cllcalipsoice, 'cllcalipsoice', &
          & 'Lidar Low-level Ice Cloud Fraction', '%', 'ASFC', hclas )
  endif
  if (Lclmcalipsoice) then
     clmcalipsoice(:) = cospOUT%calipso_cldlayerphase(:,2,1)
     call undef_n2m( clmcalipsoice, Npoints, 1, 1 )
     call histin( clmcalipsoice, 'clmcalipsoice', &
          & 'Lidar Mid-level Ice Cloud Fraction', '%', 'ASFC', hclas )
  endif
  if (Lclhcalipsoice) then
     clhcalipsoice(:) = cospOUT%calipso_cldlayerphase(:,3,1)
     call undef_n2m( clhcalipsoice, Npoints, 1, 1 )
     call histin( clhcalipsoice, 'clhcalipsoice', &
          & 'Lidar High-level Ice Cloud Fraction', '%', 'ASFC', hclas )
  endif
  if (Lcltcalipsoice) then
     cltcalipsoice(:) = cospOUT%calipso_cldlayerphase(:,4,1)
     call undef_n2m( cltcalipsoice, Npoints, 1, 1 )
     call histin( cltcalipsoice, 'cltcalipsoice', &
          & 'Lidar Total Ice Cloud Fraction', '%', 'ASFC', hclas )
  endif
  if (Lcllcalipsoliq) then
     cllcalipsoliq(:) = cospOUT%calipso_cldlayerphase(:,1,2)
     call undef_n2m( cllcalipsoliq, Npoints, 1, 1 )
     call histin( cllcalipsoliq, 'cllcalipsoliq', &
          & 'Lidar Low-level Liq Cloud Fraction', '%', 'ASFC', hclas )
  endif
  if (Lclmcalipsoliq) then
     clmcalipsoliq(:) = cospOUT%calipso_cldlayerphase(:,2,2)
     call undef_n2m( clmcalipsoliq, Npoints, 1, 1 )
     call histin( clmcalipsoliq, 'clmcalipsoliq', &
          & 'Lidar Mid-level Liq Cloud Fraction', '%', 'ASFC', hclas )
  endif
  if (Lclhcalipsoliq) then
     clhcalipsoliq(:) = cospOUT%calipso_cldlayerphase(:,3,2)
     call undef_n2m( clhcalipsoliq, Npoints, 1, 1 )
     call histin( clhcalipsoliq, 'clhcalipsoliq', &
          & 'Lidar High-level Liq Cloud Fraction', '%', 'ASFC', hclas )
  endif
  if (Lcltcalipsoliq) then
     cltcalipsoliq(:) = cospOUT%calipso_cldlayerphase(:,4,2)
     call undef_n2m( cltcalipsoliq, Npoints, 1, 1 )
     call histin( cltcalipsoliq, 'cltcalipsoliq', &
          & 'Lidar Total Liq Cloud Fraction', '%', 'ASFC', hclas )
  endif
  if (Lcllcalipsoun) then
     cllcalipsoun(:) = cospOUT%calipso_cldlayerphase(:,1,3)
     call undef_n2m( cllcalipsoun, Npoints, 1, 1 )
     call histin( cllcalipsoun, 'cllcalipsoun', &
          & 'Lidar Low-level Undefined-phase Cloud Fraction', '%', 'ASFC', hclas )
  endif
  if (Lclmcalipsoun) then
     clmcalipsoun(:) = cospOUT%calipso_cldlayerphase(:,2,3)
     call undef_n2m( clmcalipsoun, Npoints, 1, 1 )
     call histin( clmcalipsoun, 'clmcalipsoun', &
          & 'Lidar Mid-level Undefined-phase Cloud Fraction', '%', 'ASFC', hclas )
  endif
  if (Lclhcalipsoun) then
     clhcalipsoun(:) = cospOUT%calipso_cldlayerphase(:,3,3)
     call undef_n2m( clhcalipsoun, Npoints, 1, 1 )
     call histin( clhcalipsoun, 'clhcalipsoun', &
          & 'Lidar High-level Undefined-phase Cloud Fraction', '%', 'ASFC', hclas )
  endif
  if (Lcltcalipsoun) then
     cltcalipsoun(:) = cospOUT%calipso_cldlayerphase(:,4,3)
     call undef_n2m( cltcalipsoun, Npoints, 1, 1 )
     call histin( cltcalipsoun, 'cltcalipsoun', &
          & 'Lidar Total Undefined-phase Cloud Fraction', '%', 'ASFC', hclas )
  endif
  if (Lcltisccp) then
     cltisccp = cospOUT%isccp_totalcldarea
     call undef_n2m( cltisccp, Npoints, 1, 1 )
     call histin( cltisccp, 'cltisccp', &
          & 'Total Cloud Fraction by ISCCP Simulator', '%', 'ASFC', hclas )
  endif
  if (Lpctisccp) then
     pctisccp = cospOUT%isccp_meanptop
     call undef_n2m( pctisccp, Npoints, 1, 1 )
     call histin( pctisccp, 'pctisccp', &
          & 'Mean Cloud Top Pressure by ISCCP Simulator', 'Pa', 'ASFC', hclas )
  endif
  if (Ltauisccp) then
     tauisccp = cospOUT%isccp_meantaucld
     call undef_n2m( tauisccp, Npoints, 1, 1 )
     call histin( tauisccp, 'tauisccp', &
          & 'Mean Optical Depth by ISCCP Simulator', '1', 'ASFC', hclas )
  endif
  if (Lalbisccp) then
     albisccp = cospOUT%isccp_meanalbedocld
     call undef_n2m( albisccp, Npoints, 1, 1 )
     call histin( albisccp, 'albisccp', &
          & 'Mean Cloud Albedo by ISCCP Simulator', '1', 'ASFC', hclas )
  endif
  if (Lmeantbisccp) then
     meantbisccp = cospOUT%isccp_meantb
     call undef_n2m( meantbisccp, Npoints, 1, 1 )
     call histin( meantbisccp, 'meantbisccp', &
          & 'Mean all-sky 10.5 micron BT by ISCCP Simulator', 'K', 'ASFC', hclas )
  endif
  if (Lmeantbclrisccp) then
     meantbclrisccp = cospOUT%isccp_meantbclr
     call undef_n2m( meantbclrisccp, Npoints, 1, 1 )
     call histin( meantbclrisccp, 'meantbclrisccp', &
          & 'Mean clear-sky 10.5 micron BT by ISCCP Simulator', 'K', 'ASFC', hclas )
  endif
  if (Lcltmodis) then
     cltmodis = cospOUT%modis_Cloud_Fraction_Total_Mean
     call undef_n2m( cltmodis, Npoints, 1, 1 )
     call histin( cltmodis, 'cltmodis', &
          & 'MODIS total cloud fraction', '%', 'ASFC', hclas )
  endif
  if (Lclwmodis) then
     clwmodis = cospOUT%modis_Cloud_Fraction_Water_Mean
     call undef_n2m( clwmodis, Npoints, 1, 1 )
     call histin( clwmodis, 'clwmodis', &
          & 'MODIS liquid cloud fraction', '%', 'ASFC', hclas )
  endif
  if (Lclimodis) then
     climodis = cospOUT%modis_Cloud_Fraction_Ice_Mean
     call undef_n2m( climodis, Npoints, 1, 1 )
     call histin( climodis, 'climodis', &
          & 'MODIS ice cloud fraction', '%', 'ASFC', hclas )
  endif
  if (Lclhmodis) then
     clhmodis = cospOUT%modis_Cloud_Fraction_High_Mean
     call undef_n2m( clhmodis, Npoints, 1, 1 )
     call histin( clhmodis, 'clhmodis', &
          & 'MODIS high-level cloud fraction', '%', 'ASFC', hclas )
  endif
  if (Lclmmodis) then
     clmmodis = cospOUT%modis_Cloud_Fraction_Mid_Mean
     call undef_n2m( clmmodis, Npoints, 1, 1 )
     call histin( clmmodis, 'clmmodis', &
          & 'MODIS mid-level cloud fraction', '%', 'ASFC', hclas )
  endif
  if (Lcllmodis) then
     cllmodis = cospOUT%modis_Cloud_Fraction_Low_Mean
     call undef_n2m( cllmodis, Npoints, 1, 1 )
     call histin( cllmodis, 'cllmodis', &
          & 'MODIS low-level cloud fraction', '%', 'ASFC', hclas )
  endif
  if (Ltautmodis) then
     tautmodis = cospOUT%modis_Optical_Thickness_Total_Mean
     call undef_n2m( tautmodis, Npoints, 1, 1 )
     call histin( tautmodis, 'tautmodis', &
          & 'MODIS Total Cloud Optical Thickness', '1', 'ASFC', hclas )
  endif
  if (Ltauwmodis) then
     tauwmodis = cospOUT%modis_Optical_Thickness_Water_Mean
     call undef_n2m( tauwmodis, Npoints, 1, 1 )
     call histin( tauwmodis, 'tauwmodis', &
          & 'MODIS Liquid Cloud Optical Thickness', '1', 'ASFC', hclas )
  endif
  if (Ltauimodis) then
     tauimodis = cospOUT%modis_Optical_Thickness_Ice_Mean
     call undef_n2m( tauimodis, Npoints, 1, 1 )
     call histin( tauimodis, 'tauimodis', &
          & 'MODIS Ice Cloud Optical Thickness', '1', 'ASFC', hclas )
  endif
  if (Ltautlogmodis) then
     tautlogmodis = cospOUT%modis_Optical_Thickness_Total_LogMean
     call undef_n2m( tautlogmodis, Npoints, 1, 1 )
     call histin( tautlogmodis, 'tautlogmodis', &
          & 'MODIS Total Cld Opt Thcknss (Log10 Mean)', '1', 'ASFC', hclas )
  endif
  if (Ltauwlogmodis) then
     tauwlogmodis = cospOUT%modis_Optical_Thickness_Water_LogMean
     call undef_n2m( tauwlogmodis, Npoints, 1, 1 )
     call histin( tauwlogmodis, 'tauwlogmodis', &
          & 'MODIS Liq Cld Opt Thcknss (Log10 Mean)', '1', 'ASFC', hclas )
  endif
  if (Ltauilogmodis) then
     tauilogmodis = cospOUT%modis_Optical_Thickness_Ice_LogMean
     call undef_n2m( tauilogmodis, Npoints, 1, 1 )
     call histin( tauilogmodis, 'tauilogmodis', &
          & 'MODIS Ice Cld Opt Thcknss (Log10 Mean)', '1', 'ASFC', hclas )
  endif
  if (Lreffclwmodis) then
     reffclwmodis = cospOUT%modis_Cloud_Particle_Size_Water_Mean
     call undef_n2m( reffclwmodis, Npoints, 1, 1 )
     call histin( reffclwmodis, 'reffclwmodis', &
          & 'MODIS Liq Cld effective radius', 'm', 'ASFC', hclas )
  endif
  if (Lreffclimodis) then
     reffclimodis = cospOUT%modis_Cloud_Particle_Size_Ice_Mean
     call undef_n2m( reffclimodis, Npoints, 1, 1 )
     call histin( reffclimodis, 'reffclimodis', &
          & 'MODIS Ice Cld effective radius', 'm', 'ASFC', hclas )
  endif
  if (Lpctmodis) then
     pctmodis = cospOUT%modis_Cloud_Top_Pressure_Total_Mean
     call undef_n2m( pctmodis, Npoints, 1, 1 )
     call histin( pctmodis, 'pctmodis', &
          & 'MODIS cloud top pressure', 'Pa', 'ASFC', hclas )
  endif
  if (Llwpmodis) then
     lwpmodis = cospOUT%modis_Liquid_Water_Path_Mean
     call undef_n2m( lwpmodis, Npoints, 1, 1 )
     call histin( lwpmodis, 'lwpmodis', &
          & 'MODIS cloud liq water path', 'kg m-2', 'ASFC', hclas )
  endif
  if (Liwpmodis) then
     iwpmodis = cospOUT%modis_Ice_Water_Path_Mean
     call undef_n2m( iwpmodis, Npoints, 1, 1 )
     call histin( iwpmodis, 'iwpmodis', &
          & 'MODIS cloud ice water path', 'kg m-2', 'ASFC', hclas )
  endif
  if (Lclwmodis .and. Lreffclwmodis) then
     do j=1,Npoints
        if ( clwmodis(j) .gt. -998.D0 .and. &
           & reffclwmodis(j) .gt. -998.D0     ) then
           clwxreffclw(j) = clwmodis(j) * reffclwmodis(j)
        else
           clwxreffclw(j) = -999.D0
        endif
     enddo
     call undef_n2m( clwxreffclw, Npoints, 1, 1 )
     call histin( clwxreffclw, 'clwxreffclw', &
          & 'clwmodis * reffclwmodis', 'm', 'ASFC', hclas )
  endif
  if (Lclimodis .and. Lreffclimodis) then
     do j=1,Npoints
        if ( climodis(j) .gt. -998.D0 .and. &
           & reffclimodis(j) .gt. -998.D0     ) then
           clixreffcli(j) = climodis(j) * reffclimodis(j)
        else
           clixreffcli(j) = -999.D0
        endif
     enddo
     call undef_n2m( clixreffcli, Npoints, 1, 1 )
     call histin( clixreffcli, 'clixreffcli', &
          & 'climodis * reffclimodis', 'm', 'ASFC', hclas )
  endif
!!!  if (Lcltcalipso) then
!!!     cltcalipso(:) = cospOUT%cldlayer(:,4)
!!!     call undef_n2m( cltcalipso, Npoints, 1, 1 )
!!!     call histin( cltcalipso, 'cltcalipso', &
!!!          & 'Lidar Total Cloud Fraction', '%', 'ASFC', hclas )
!!!  endif
!!!  if (Lcltisccp) then
!!!     cltisccp = cospOUT%totalcldarea
!!!     call undef_n2m( cltisccp, Npoints, 1, 1 )
!!!     call histin( cltisccp, 'cltisccp', &
!!!          & 'Total Cloud Fraction by ISCCP Simulator', '%', 'ASFC', hclas )
!!!  endif
  !
  ! 2D variables
  if (Lclcalipso) then
     clcalipso = cospOUT%calipso_lidarcld(:,Nlvgrid:1:-1)
     ncls = Nlvgrid
     call undef_n2m( clcalipso, Npoints, Nlvgrid, 1 )
     call histnn( clcalipso, 'clcalipso', &
          & 'CALIPSO cloud area fraction', '%', 'A', hclas, ncls )
  endif
  if (Lclcalipso2) then
     clcalipso2 = cospOUT%lidar_only_freq_cloud(:,Nlvgrid:1:-1)
     ncls = Nlvgrid
     call undef_n2m( clcalipso2, Npoints, Nlvgrid, 1 )
     call histnn( clcalipso2, 'clcalipso2', &
          & 'CALIPSO cld frac undetected by CloudSat', '%', 'A', hclas, ncls )
  endif
  if (Lclcalipsoice) then
     clcalipsoice(:,:) = cospOUT%calipso_lidarcldphase(:,Nlvgrid:1:-1,1)
     ncls = Nlvgrid
     call undef_n2m( clcalipsoice, Npoints, Nlvgrid, 1 )
     call histnn( clcalipsoice, 'clcalipsoice', &
          & 'CALIPSO ice cld frac', '%', 'A', hclas, ncls )
  endif
  if (Lclcalipsoliq) then
     clcalipsoliq(:,:) = cospOUT%calipso_lidarcldphase(:,Nlvgrid:1:-1,2)
     ncls = Nlvgrid
     call undef_n2m( clcalipsoliq, Npoints, Nlvgrid, 1 )
     call histnn( clcalipsoliq, 'clcalipsoliq', &
          & 'CALIPSO liq cld frac', '%', 'A', hclas, ncls )
  endif
  if (Lclcalipsoun) then
     clcalipsoun(:,:) = cospOUT%calipso_lidarcldphase(:,Nlvgrid:1:-1,3)
     ncls = Nlvgrid
     call undef_n2m( clcalipsoun, Npoints, Nlvgrid, 1 )
     call histnn( clcalipsoun, 'clcalipsoun', &
          & 'CALIPSO undef-phase cld frac', '%', 'A', hclas, ncls )
  endif
  if (Lclcalipsotmp) then
     clcalipsotmp(:,:) = cospOUT%calipso_lidarcldtmp(:,:,1)
     ncls = LIDAR_NTEMP
     call undef_n2m( clcalipsotmp, Npoints, LIDAR_NTEMP, 1 )
     call histnn( clcalipsotmp, 'clcalipsotmp', &
          & 'CALIPSO cld frac wrt temperature', '%', 'A', hclas, ncls )
  endif
  if (Lclcalipsotmpice) then
     clcalipsotmpice(:,:) = cospOUT%calipso_lidarcldtmp(:,:,2)
     ncls = LIDAR_NTEMP
     call undef_n2m( clcalipsotmpice, Npoints, LIDAR_NTEMP, 1 )
     call histnn( clcalipsotmpice, 'clcalipsotmpice', &
          & 'CALIPSO ice cld frac wrt temperature', '%', 'A', hclas, ncls )
  endif
  if (Lclcalipsotmpliq) then
     clcalipsotmpliq(:,:) = cospOUT%calipso_lidarcldtmp(:,:,3)
     ncls = LIDAR_NTEMP
     call undef_n2m( clcalipsotmpliq, Npoints, LIDAR_NTEMP, 1 )
     call histnn( clcalipsotmpliq, 'clcalipsotmpliq', &
          & 'CALIPSO liq cld frac wrt temperature', '%', 'A', hclas, ncls )
  endif
  if (Lclcalipsotmpun) then
     clcalipsotmpun(:,:) = cospOUT%calipso_lidarcldtmp(:,:,4)
     ncls = LIDAR_NTEMP
     call undef_n2m( clcalipsotmpun, Npoints, LIDAR_NTEMP, 1 )
     call histnn( clcalipsotmpun, 'clcalipsotmpun', &
          & 'CALIPSO undef-phase cld frac wrt temperature', '%', 'A', hclas, ncls )
  endif
  if (LlidarBetaMol532) then
     lidarbetamol532 = cospOUT%calipso_beta_mol
     ncls = Nlevels
     call undef_n2m( lidarbetamol532, Npoints, Nlevels, 1 )
     call histnn( lidarbetamol532, 'lidarbetamol532', &
          & 'Lidar Molecular Backscatter (532 nm)', 'm-1 sr-1', 'A', hclas, ncls )
  endif
  if (Lboxtauisccp) then
     boxtauisccp = cospOUT%isccp_boxtau
     ncls = Ncolumns
     call undef_n2m( boxtauisccp, Npoints, Ncolumns, 1 )
     call histnn( boxtauisccp, 'boxtauisccp', &
          & 'Optical Depth in Each Column', '1', 'A', hclas, ncls )
  endif
  if (Lboxptopisccp) then
     boxptopisccp = cospOUT%isccp_boxptop
     ncls = Ncolumns
     call undef_n2m( boxptopisccp, Npoints, Ncolumns, 1 )
     call histnn( boxptopisccp, 'boxptopisccp', &
          & 'Cloud Top Pressure in Each Column', 'Pa', 'A', hclas, ncls )
  endif
  if (LparasolRefl) then
     parasolGrid_Refl = cospOUT%parasolGrid_refl
     ncls = PARASOL_NREFL
     call undef_n2m( parasolGrid_Refl, Npoints, PARASOL_NREFL, 1 )
     call histnn( parasolGrid_Refl, 'parasolRefl', &
          & 'PARASOL-like mono-directional reflectance', '1', 'A', hclas, ncls )
  endif
  ! <--- added by YN
  if (Lscops) then
     allocate(out_array(npoints,nlevels,4))

     out_array(:,:,1) = cospOUT%frac_ls
     call undef_n2m(out_array,Npoints,Nlevels,1)
     call histin(out_array(:,:,1),'frac_ls','frac_ls scops',' ','ALEV',hclas)
     
     out_array(:,:,2) = cospOUT%frac_cv
     call undef_n2m(out_array(:,:,2),Npoints,Nlevels,1)
     call histin(out_array(:,:,2),'frac_cv','frac_cv scops',' ','ALEV',hclas)

     out_array(:,:,3) = cospOUT%prec_ls
     call undef_n2m(out_array(:,:,3),Npoints,Nlevels,1)
     call histin(out_array(:,:,3),'prec_ls','prec_ls prec scops',' ','ALEV',hclas)

     out_array(:,:,4) = cospOUT%prec_cv
     call undef_n2m(out_array(:,:,4),Npoints,Nlevels,1)
     call histin(out_array(:,:,4),'prec_cv','prec_cv prec scops',' ','ALEV',hclas)

     deallocate(out_array)
  end if

  !
  ! 3D variables
  if (Lclisccp) then
     clisccp = cospOUT%isccp_fq
     ncls = 7 * 7
     call undef_n2m( clisccp, Npoints, 7, 7 )
     call histnn( clisccp, 'clisccp', &
          & 'Cloud Fraction by ISCCP Simulator', '%', 'A', hclas, ncls )
  endif
  if (Lclmodis) then
     clmodis = cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure
     ncls = 7 * 7
     call undef_n2m( clmodis, Npoints, 7, 7 )
     call histnn( clmodis, 'clmodis', &
          & 'Cloud Fraction by MODIS Simulator', '%', 'A', hclas, ncls )
  endif
  if (LclMISR) then
     clmisr = cospOUT%MISR_fq
     ncls = 7 * numMISRHgtBins
     call undef_n2m( clmisr, Npoints, 7, numMISRHgtBins )
     call histnn( clmisr, 'clmisr', &
          & 'Cloud Fraction by MISR Simulator', '%', 'A', hclas, ncls )
  endif
  if (LcfadLidarsr532) then
     cfadlidarsr532 = cospOUT%calipso_cfad_sr(:,:,Nlvgrid:1:-1)
     ncls = SR_BINS * Nlvgrid
     call undef_n2m( cfadlidarsr532, Npoints, SR_BINS, Nlvgrid )
     call histnn( cfadlidarsr532, 'cfadlidarsr532', &
          & 'Histogram of CALIPSO scattering ratio', '1', 'A', hclas, ncls )
  endif
  if (LcfadDbze94) then
     cfaddbze94 = cospOUT%cloudsat_cfad_ze(:,:,Nlvgrid:1:-1)
     ncls = DBZE_BINS * Nlvgrid
     call undef_n2m( cfaddbze94, Npoints, DBZE_BINS, Nlvgrid )
     call histnn( cfaddbze94, 'cfaddbze94', &
          & 'Histogram of CloudSat reflectivity', '1', 'A', hclas, ncls )
  endif
  !!if (Ldbze94) then
  !!   dbze94 = cospOUT%cloudsat_Ze_tot
  !!   ncls = Ncolumns * Nlevels
  !!   call undef_n2m( dbze94, Npoints, Ncolumns, Nlevels )
  !!   call histnn( dbze94, 'dbze94', &
  !!        & 'CloudSat reflectivity', '1', 'A', hclas, ncls )
  !!endif

  !-- write out each subcolumn as a single data file. Jing X.
  if (Ldbze94) then
    do i=1,Ncolumns
     write(idcol,'(i3)') i
     dbze94_tmp = cospOUT%cloudsat_Ze_tot(1:Npoints,i,1:Nlevels)
     ncls = Nlevels
     call undef_n2m( dbze94_tmp, Npoints, Nlevels, 1)
     call histin( dbze94_tmp, 'dbze'//trim(adjustl(idcol)), &
          & 'CloudSat reflectivity', '1', 'ALEV', hclas)
    end do
    !! Grid maximum dBZe (T.Michibata)
!    dbze94grdmax2d(:,:) = -999.d0
!    dbze94grdmax1d(:)   = -999.d0
!    do j = 1, Npoints
!       do k = 1, Nlevels
!          do i = 1, Ncolumns
!             if ( cospOUT%cloudsat_Ze_tot(j,i,k) .ge. -30 .and. &
!                & cospOUT%cloudsat_Ze_tot(j,i,k) .le.  20       ) then
!!                dbze94grdmax2d(j,k) = &
!!                & maxval( cospOUT%cloudsat_Ze_tot(j,:,k) )
!!                dbze94grdmax1d(j) = max( dbze94grdmax1d(j), &
!!                & maxval( cospOUT%cloudsat_Ze_tot(j,:,k) ) )
!                dbze94grdmax2d(j,k) = max( dbze94grdmax2d(j,k), &
!                                    & cospOUT%cloudsat_Ze_tot(j,i,k) )
!                dbze94grdmax1d(j) = max( dbze94grdmax1d(j), &
!                                  &   cospOUT%cloudsat_Ze_tot(j,i,k) )
!             endif
!          enddo
!       enddo
!    enddo
!
!<--- X.Jing
!    allocate(mm1(Npoints))
!    allocate(mm2(Npoints,Nlevels))
!    do j = 1, Npoints
!       do k = 1, Nlevels
!          mm1=(cospOUT%cloudsat_Ze_tot(j,:,k) .ge. -30 .and. cospOUT%cloudsat_Ze_tot(j,:,k) .le.  20)
!          dbze94grdmax2d(j,k) = maxval(cospOUT%cloudsat_Ze_tot(j,:,k),mask=mm1)
!          if(dbze94grdmax2d(j,k) .gt. 20 .or. dbze94grdmax2d(j,k) .lt. -30) dbze94grdmax2d(j,k) = -999.d0
!       enddo
!       mm2=(cospOUT%cloudsat_Ze_tot(j,:,:) .ge. -30 .and. cospOUT%cloudsat_Ze_tot(j,:,:) .le.  20)
!       dbze94grdmax1d(j) = maxval(cospOUT%cloudsat_Ze_tot(j,:,:),mask=mm2)
!       if(dbze94grdmax1d(j) .gt. 20 .or. dbze94grdmax1d(j) .lt. -30) dbze94grdmax1d(j) = -999.d0
!    enddo
!    deallocate(mm1)
!    deallocate(mm2)
!    call histin( dbze94grdmax1d, 'dbze94grdmax1d', &
!         & 'CloudSat maximum reflectivity (IJ)', '1', 'ASFC', hclas)
!    call histin( dbze94grdmax2d, 'dbze94grdmax2d', &
!         & 'CloudSat maximum reflectivity (IJ,K)', '1', 'ALEV', hclas)
!---> X.Jing
  endif

  if (Latb532) then
     atb532 = cospOUT%calipso_beta_tot
     ncls = Ncolumns * Nlevels
     call undef_n2m( atb532, Npoints, Ncolumns, Nlevels )
     call histnn( atb532, 'atb532', &
          & 'Lidar attenuated backscatter', 'm-1 sr-1', 'A', hclas, ncls )
  endif
  !!if (Lfracout) then
  !!   fracout = cospIN%frac_out
  !!   ncls = Ncolumns * Nlevels
  !!   call undef_n2m( fracout, Npoints, Ncolumns, Nlevels )
  !!   call histnn( fracout, 'fracout', &
  !!        & 'Subcolumn output from SCOPS', '1', 'A', hclas, ncls )
  !!endif

#ifdef OPT_DPLRW
  if (Ldplrw) then
     do id=1,9
        write(c1,'(I1)') id

        ! allocate(out_array(npoints, Nlvgrid, 3))
        ! out_array = reshape(cospOUT%dplrw_sumZ(:,:,id,:),(/ Npoints, Nlvgrid, 3 /) )
        ! call undef_n2m(out_array,Npoints,Nlvgrid,3)
        ! ncls = Nlvgrid*3
        ! call histnn(out_array,'dplrw_sumZ_ISC'//c1,'dplrw sumZ','#','A',hclas,ncls)
        ! deallocate(out_array)

        ! allocate(out_array(npoints, Nlvtemp, 3))
        ! out_array = reshape(cospOUT%dplrw_sumT(:,:,id,:),(/ Npoints, Nlvtemp, 3 /) )
        ! call undef_n2m(out_array,Npoints,Nlvtemp,3)
        ! ncls = Nlvtemp*3
        ! call histnn(out_array,'dplrw_sumT_ISC'//c1,'dplrw sumT','#','A',hclas,ncls)
        ! deallocate(out_array)

        !---
        allocate(out_array(npoints, Nlvdplr*Nlvgrid, 3))
        out_array = reshape(cospOUT%dplrw_Z(:,:,:,id,:),(/ Npoints, Nlvdplr*Nlvgrid, 3 /) )
        call undef_n2m(out_array,npoints,Nlvdplr*Nlvgrid,3)
        ncls = Nlvgrid * Nlvdplr * 3
        call histnn(out_array,'dplrw_Z_ISC'//c1,'dplrw Z count','#','A',hclas,ncls)
        deallocate(out_array)

        allocate(out_array(npoints, Nlvspwd*Nlvgrid, 3))
        out_array = reshape(cospOUT%spwid_Z(:,:,:,id,:),(/ Npoints, Nlvspwd*Nlvgrid, 3 /) )
        call undef_n2m(out_array,npoints,Nlvspwd*Nlvgrid,3)
        ncls = Nlvgrid * Nlvspwd * 3
        call histnn(out_array,'spwid_Z_ISC'//c1,'spwid Z count','#','A',hclas,ncls)
        deallocate(out_array)

        allocate(out_array(npoints, NlvdBZe*Nlvgrid, 3))
        out_array = reshape(cospOUT%Zef94_Z(:,:,:,id,:),(/ Npoints, NlvdBZe*Nlvgrid, 3 /) )
        call undef_n2m(out_array,npoints,NlvdBZe*Nlvgrid,3)
        ncls = Nlvgrid * NlvdBZe * 3
        call histnn(out_array,'Zef94_Z_ISC'//c1,'Zef94 Z count','#','A',hclas,ncls)
        deallocate(out_array)

        !---
        allocate(out_array(npoints, Nlvdplr*Nlvtemp, 3))
        out_array = reshape(cospOUT%dplrw_T(:,:,:,id,:),(/ Npoints, Nlvdplr*Nlvtemp, 3 /) )
        call undef_n2m(out_array,npoints,Nlvdplr*Nlvtemp,3)
        ncls = Nlvtemp * Nlvdplr * 3
        call histnn(out_array,'dplrw_T_ISC'//c1,'dplrw T count','#','A',hclas,ncls)
        deallocate(out_array)

        allocate(out_array(npoints, Nlvspwd*Nlvtemp, 3))
        out_array = reshape(cospOUT%spwid_T(:,:,:,id,:),(/ Npoints, Nlvspwd*Nlvtemp, 3 /) )
        call undef_n2m(out_array,npoints,Nlvspwd*Nlvtemp,3)
        ncls = Nlvtemp * Nlvspwd * 3
        call histnn(out_array,'spwid_T_ISC'//c1,'spwid T count','#','A',hclas,ncls)
        deallocate(out_array)

        allocate(out_array(npoints, NlvdBZe*Nlvtemp, 3))
        out_array = reshape(cospOUT%Zef94_T(:,:,:,id,:),(/ Npoints, NlvdBZe*Nlvtemp, 3 /) )
        call undef_n2m(out_array,npoints,NlvdBZe*Nlvtemp,3)
        ncls = Nlvtemp * NlvdBZe * 3
        call histnn(out_array,'Zef94_T_ISC'//c1,'Zef94 T count','#','A',hclas,ncls)
        deallocate(out_array)

        !---
        allocate(out_array(npoints, Nlvdplr*NlvdBZe, 3))
        out_array = reshape(cospOUT%ZefVd_2(:,:,:,id,:),(/ Npoints, Nlvdplr*NlvdBZe, 3 /) )
        call undef_n2m(out_array,npoints,Nlvdplr*NlvdBZe,3)
        ncls = Nlvdplr * NlvdBZe * 3
        call histnn(out_array,'ZefVd_2_ISC'//c1,'ZefVd 2 count','#','A',hclas,ncls)
        deallocate(out_array)
     end do
  end if
#endif

  !-- write out each subcolumn as a single data file. Jing X.
  if (Lfracout) then
    do i=1,Ncolumns
     write(idcol,'(i3)') i
     fracout_tmp = cospIN%frac_out(1:Npoints,i,Nlevels:1:-1)
     ncls = Nlevels
     call undef_n2m(fracout_tmp, Npoints, Nlevels, 1)
     call histin( fracout_tmp, 'cfrac'//trim(adjustl(idcol)), &
          & 'subcolumn output from SCOPS', '1', 'ALEV', hclas)
    end do
  endif

  !--- Inline Diagnostics OUTPUTS (T.Michibata, 2018.11.22):
  if (Lwr_occfreq) then
     !! Occurrence Frequency Map GENERATOR
     allocate( idid_stats_1d(Npoints) )
        !! # of non-precipitating clouds (dBZmax < -15)
     idid_stats_1d(:) = cospOUT%wr_occfreq_ntotal(:,1)
     call undef_n2m( idid_stats_1d, Npoints, 1, 1 )
     call histin( idid_stats_1d, 'npdfcld', &
         & '# of non-precip clouds', '#', 'ASFC', hclas )
        !! # of drizzling clouds (-15 < dBZmax < 0)
     idid_stats_1d(:) = cospOUT%wr_occfreq_ntotal(:,2)
     call undef_n2m( idid_stats_1d, Npoints, 1, 1 )
     call histin( idid_stats_1d, 'npdfdrz', &
         & '# of drizzling clouds', '#', 'ASFC', hclas )
        !! # of precipitating clouds (0 < dBZmax)
     idid_stats_1d(:) = cospOUT%wr_occfreq_ntotal(:,3)
     call undef_n2m( idid_stats_1d, Npoints, 1, 1 )
     call histin( idid_stats_1d, 'npdfrain', &
         & '# of precipitating clouds', '#', 'ASFC', hclas )
     deallocate( idid_stats_1d )
  endif

  if (Lcfodd) then
     !! CFODD GENERATOR
     allocate( idid_stats_3d(Npoints,CFODD_NDBZE,CFODD_NICOD) )
     ncls = CFODD_NDBZE * CFODD_NICOD
        !! # of CFODD for Re < 12 um
     idid_stats_3d(:,:,:) = cospOUT%cfodd_ntotal(:,:,:,1)
     call undef_n2m( idid_stats_3d, Npoints, CFODD_NDBZE, CFODD_NICOD )
     call histnn( idid_stats_3d, 'ncfodd1', &
          & '# of CFODD (Re < 12 um)', '#', 'A', hclas, ncls )
        !! # of CFODD for 12 um < Re < 18 um
     idid_stats_3d(:,:,:) = cospOUT%cfodd_ntotal(:,:,:,2)
     call undef_n2m( idid_stats_3d, Npoints, CFODD_NDBZE, CFODD_NICOD )
     call histnn( idid_stats_3d, 'ncfodd2', &
          & '# of CFODD (12 um < Re < 18 um)', '#', 'A', hclas, ncls )
        !! # of CFODD for 18 um < Re
     idid_stats_3d(:,:,:) = cospOUT%cfodd_ntotal(:,:,:,3)
     call undef_n2m( idid_stats_3d, Npoints, CFODD_NDBZE, CFODD_NICOD )
     call histnn( idid_stats_3d, 'ncfodd3', &
          & '# of CFODD (18 um < Re)', '#', 'A', hclas, ncls )
!--> removed (2019.04.14)
        !! # of CFODD for all Re range
!     idid_stats_3d(:,:,:) = cospOUT%cfodd_ntotal(:,:,:,4)
!     call undef_n2m( idid_stats_3d, Npoints, CFODD_NDBZE, CFODD_NICOD )
!     call histnn( idid_stats_3d, 'ncfodd4', &
!          & '# of CFODD (all Re range)', '#', 'A', hclas, ncls )
!<-- removed (2019.04.14)
     deallocate( idid_stats_3d )
  endif

!<-added
  if (Lptradarflag0) then
     allocate( ptcloudsatflag0(Npoints) )
     ptcloudsatflag0(:) = cospOUT%cloudsat_precip_cover(:,1)
     call undef_n2m( ptcloudsatflag0, Npoints, 1, 1 )
     call histin( ptcloudsatflag0, 'ptcloudsatflag0', &
          & 'flag0: non-precip', '%', 'ASFC', hclas )
     deallocate( ptcloudsatflag0 )
  endif

  if (Lptradarflag1) then
     allocate( ptcloudsatflag1(Npoints) )
     ptcloudsatflag1(:) = cospOUT%cloudsat_precip_cover(:,2)
     call undef_n2m( ptcloudsatflag1, Npoints, 1, 1 )
     call histin( ptcloudsatflag1, 'ptcloudsatflag1', &
          & 'flag1: rain possible', '%', 'ASFC', hclas )
     deallocate( ptcloudsatflag1 )
  endif

  if (Lptradarflag2) then
     allocate( ptcloudsatflag2(Npoints) )
     ptcloudsatflag2(:) = cospOUT%cloudsat_precip_cover(:,3)
     call undef_n2m( ptcloudsatflag2, Npoints, 1, 1 )
     call histin( ptcloudsatflag2, 'ptcloudsatflag2', &
          & 'flag2: rain probable', '%', 'ASFC', hclas )
     deallocate( ptcloudsatflag2 )
  endif

  if (Lptradarflag3) then
     allocate( ptcloudsatflag3(Npoints) )
     ptcloudsatflag3(:) = cospOUT%cloudsat_precip_cover(:,4)
     call undef_n2m( ptcloudsatflag3, Npoints, 1, 1 )
     call histin( ptcloudsatflag3, 'ptcloudsatflag3', &
          & 'flag3: rain certain', '%', 'ASFC', hclas )
     deallocate( ptcloudsatflag3 )
  endif

  if (Lptradarflag4) then
     allocate( ptcloudsatflag4(Npoints) )
     ptcloudsatflag4(:) = cospOUT%cloudsat_precip_cover(:,5)
     call undef_n2m( ptcloudsatflag4, Npoints, 1, 1 )
     call histin( ptcloudsatflag4, 'ptcloudsatflag4', &
          & 'flag4: snow possible', '%', 'ASFC', hclas )
     deallocate( ptcloudsatflag4 )
  endif

  if (Lptradarflag5) then
     allocate( ptcloudsatflag5(Npoints) )
     ptcloudsatflag5(:) = cospOUT%cloudsat_precip_cover(:,6)
     call undef_n2m( ptcloudsatflag5, Npoints, 1, 1 )
     call histin( ptcloudsatflag5, 'ptcloudsatflag5', &
          & 'flag5: snow certain', '%', 'ASFC', hclas )
     deallocate( ptcloudsatflag5 )
  endif

  if (Lptradarflag6) then
     allocate( ptcloudsatflag6(Npoints) )
     ptcloudsatflag6(:) = cospOUT%cloudsat_precip_cover(:,7)
     call undef_n2m( ptcloudsatflag6, Npoints, 1, 1 )
     call histin( ptcloudsatflag6, 'ptcloudsatflag6', &
          & 'flag6: mixed possible', '%', 'ASFC', hclas )
     deallocate( ptcloudsatflag6 )
  endif

  if (Lptradarflag7) then
     allocate( ptcloudsatflag7(Npoints) )
     ptcloudsatflag7(:) = cospOUT%cloudsat_precip_cover(:,8)
     call undef_n2m( ptcloudsatflag7, Npoints, 1, 1 )
     call histin( ptcloudsatflag7, 'ptcloudsatflag7', &
          & 'flag7: mixed certain', '%', 'ASFC', hclas )
     deallocate( ptcloudsatflag7 )
  endif

  if (Lptradarflag8) then
     allocate( ptcloudsatflag8(Npoints) )
     ptcloudsatflag8(:) = cospOUT%cloudsat_precip_cover(:,9)
     call undef_n2m( ptcloudsatflag8, Npoints, 1, 1 )
     call histin( ptcloudsatflag8, 'ptcloudsatflag8', &
          & 'flag8: heavy rain', '%', 'ASFC', hclas )
     deallocate( ptcloudsatflag8 )
  endif

  if (Lptradarflag9) then
     allocate( ptcloudsatflag9(Npoints) )
     ptcloudsatflag9(:) = cospOUT%cloudsat_precip_cover(:,10)
     call undef_n2m( ptcloudsatflag9, Npoints, 1, 1 )
     call histin( ptcloudsatflag9, 'ptcloudsatflag9', &
          & 'flag9: uncertain', '%', 'ASFC', hclas )
     deallocate( ptcloudsatflag9 )
  endif

  if (Lradarpia) then
     allocate( cloudsatpia(Npoints) )
     cloudsatpia(:) = cospOUT%cloudsat_pia(:)
     call undef_n2m( cloudsatpia, Npoints, 1, 1 )
     call histin( cloudsatpia, 'cloudsatpia', &
          & 'Cloudsat path integrated attenuation', '1', 'ASFC', hclas )
     deallocate( cloudsatpia )
  endif
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Modification of undef variables: -1.0E30 --> -999.D0 
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  call undef_n2m(p, Npoints, Nlevels, 1)
  call undef_n2m(ph, Npoints, Nlevels+1, 1)
  call undef_n2m(zlev, Npoints, Nlevels, 1)
  call undef_n2m(zlev_half, Npoints, Nlevels+1, 1)
  call undef_n2m(T, Npoints, Nlevels, 1)
  call undef_n2m(rh, Npoints, Nlevels, 1)
  call undef_n2m(sh, Npoints, Nlevels, 1)
  call undef_n2m(cca, Npoints, Nlevels, 1)
  call undef_n2m(tca, Npoints, Nlevels, 1)
  call undef_n2m(skt, Npoints, 1, 1) 
  call undef_n2m(surfelev, Npoints, 1, 1)
  call undef_n2m(landmask, Npoints, 1, 1)
  call undef_n2m(mr_ozone, Npoints, Nlevels, 1)
  call undef_n2m(sunlit, Npoints, 1, 1)
  call undef_n2m(Reff, Npoints, Nlevels, 9)
  call undef_n2m(dtau_s, Npoints, Nlevels, 1)
  call undef_n2m(dtau_c, Npoints, Nlevels, 1)
  call undef_n2m(dem_s, Npoints, Nlevels, 1)
  call undef_n2m(dem_c, Npoints, Nlevels, 1)
#ifdef OPT_DPLRW
  call undef_n2m(gwvel, Npoints, Nlevels, 1)
  call undef_n2m(gcumf, Npoints, Nlevels, 1)
#endif
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Free up memory
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  call destroy_cosp_outputs(cospOUT)
  call destroy_cospIN(cospIN)
  call destroy_cospstateIN(cospstateIN)
  call cosp_cleanUp()
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  ! Free up local memory
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  !deallocate(u_wind,v_wind,cfadlidarsr532,cfaddbze94)
  deallocate( clcalipso, clcalipso2,              &
              clcalipsoice, parasolGrid_Refl,   &
              parasolPix_Refl,                          &
              clcalipsoliq, clcalipsoun,              &
              clcalipsotmp, clcalipsotmpice,  &
              clcalipsotmpliq, clcalipsotmpun,&
              lidarbetamol532, boxtauisccp,          &
              boxptopisccp, cfadlidarsr532,  &
              cfaddbze94, dbze94,  &
              atb532, fracout )
  call cpu_time(driver_time(8))

  !print*,'Time to read in data:     ',driver_time(2)-driver_time(1)
  !print*,'Time to initialize:       ',driver_time(3)-driver_time(2)
  !print*,'Time to construct types:  ',driver_time(4)-driver_time(3)
  !print*,'Time to compute optics:   ',driver_time(6)-driver_time(4)
  !print*,'Time to run COSP:         ',driver_time(7)-driver_time(6)
  !print*,'Time to write to output:  ',driver_time(8)-driver_time(7)
  !print*,'Total time:               ',driver_time(8)-driver_time(1)

return
END SUBROUTINE COSPEXE

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  ! SUBROUTINE subsample_and_optics
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine subsample_and_optics(nPoints, nLevels, nColumns, nHydro, overlap,              &
       use_precipitation_fluxes, Lcalipso, Lcloudsat, Lmodis, &
       lidar_ice_type, sd, tca, cca, fl_lsrainIN, fl_lssnowIN,    &
       fl_lsgrplIN, fl_ccrainIN, fl_ccsnowIN, mr_lsliq, mr_lsice, mr_ccliq, mr_ccice,       &
       reffIN, dtau_c, dtau_s, dem_c, dem_s, &
#ifdef OPT_CHIMERRA
       mr_lsrain, mr_lssnow, mr_lsgrpl, &
#endif
       cospstateIN, cospIN)
  USE COSP_KINDS,   ONLY: wp,dp
  USE mod_quickbeam_optics,only: size_distribution,quickbeam_optics
  use quickbeam,           only: radar_cfg
  USE MOD_COSP_CONFIG,     ONLY: R_UNDEF,N_HYDRO, Nlvgrid,              & !<-added
                                 vgrid_zl, vgrid_zu, vgrid_z, cloudsat_preclvl, use_vgrid
  USE mod_quickbeam_optics,only: size_distribution,quickbeam_optics,gases !<-added
  use mod_cosp,            only: cosp_init,cosp_optical_inputs,cosp_column_inputs,        &
                                 cosp_outputs,cosp_cleanUp,cosp_simulator
  USE mod_rng,             ONLY: rng_state, init_rng
  USE mod_scops,           ONLY: scops
  USE mod_prec_scops,      ONLY: prec_scops
  USE MOD_COSP_UTILS,      ONLY: cosp_precip_mxratio
  use cosp_optics,         ONLY: cosp_simulator_optics,lidar_optics,modis_optics,         &
                                 modis_optics_partition
  use quickbeam,           only: radar_cfg
  use mod_cosp_stats,      ONLY: COSP_CHANGE_VERTICAL_GRID !<-added

  ! Indices to address arrays of LS and CONV hydrometeors
  integer,parameter :: &
       I_LSCLIQ = 1, & ! Large-scale (stratiform) liquid
       I_LSCICE = 2, & ! Large-scale (stratiform) ice
       I_LSRAIN = 3, & ! Large-scale (stratiform) rain
       I_LSSNOW = 4, & ! Large-scale (stratiform) snow
       I_CVCLIQ = 5, & ! Convective liquid
       I_CVCICE = 6, & ! Convective ice
       I_CVRAIN = 7, & ! Convective rain
       I_CVSNOW = 8, & ! Convective snow
       I_LSGRPL = 9    ! Large-scale (stratiform) groupel
  
  ! Stratiform and convective clouds in frac_out (scops output).
  integer, parameter :: &
       I_LSC = 1, & ! Large-scale clouds
       I_CVC = 2    ! Convective clouds    
  
 ! real(wp),parameter,dimension(N_HYDRO) :: &
 !                ! LSL   LSI      LSR       LSS   CVL  CVI      CVR       CVS       LSG
 !      N_ax    = (/-1., -1.,     8.e6,     3.e6, -1., -1.,     8.e6,     3.e6,     4.e6/),&
 !      N_bx    = (/-1., -1.,      0.0,      0.0, -1., -1.,      0.0,      0.0,      0.0/),&
 !      alpha_x = (/-1., -1.,      0.0,      0.0, -1., -1.,      0.0,      0.0,      0.0/),&
 !      c_x     = (/-1., -1.,    842.0,     4.84, -1., -1.,    842.0,     4.84,     94.5/),&
 !      d_x     = (/-1., -1.,      0.8,     0.25, -1., -1.,      0.8,     0.25,      0.5/),&
 !      g_x     = (/-1., -1.,      0.5,      0.5, -1., -1.,      0.5,      0.5,      0.5/),&
 !      a_x     = (/-1., -1.,    524.0,    52.36, -1., -1.,    524.0,    52.36,   209.44/),&
 !      b_x     = (/-1., -1.,      3.0,      3.0, -1., -1.,      3.0,      3.0,      3.0/),&
 !      gamma_1 = (/-1., -1., 17.83725, 8.284701, -1., -1., 17.83725, 8.284701, 11.63230/),&
 !      gamma_2 = (/-1., -1.,      6.0,      6.0, -1., -1.,      6.0,      6.0,      6.0/),&
 !      gamma_3 = (/-1., -1.,      2.0,      2.0, -1., -1.,      2.0,      2.0,      2.0/),&
 !      gamma_4 = (/-1., -1.,      6.0,      6.0, -1., -1.,      6.0,      6.0,      6.0/)       

    ! Inputs
    integer,intent(in) :: nPoints, nLevels, nColumns, nHydro, overlap, lidar_ice_type
    real(wp),intent(in),dimension(nPoints,nLevels) :: tca,cca,mr_lsliq,mr_lsice,mr_ccliq,   &
         mr_ccice,dtau_c,dtau_s,dem_c,dem_s,fl_lsrainIN,fl_lssnowIN,fl_lsgrplIN,fl_ccrainIN,&
         fl_ccsnowIN
#ifdef OPT_CHIMERRA
    real(wp),intent(in),dimension(nPoints,nLevels) :: mr_lsrain,mr_lssnow,mr_lsgrpl
#endif
    real(wp),intent(in),dimension(nPoints,nLevels,nHydro) :: reffIN
    type(size_distribution),intent(inout) :: sd
    logical,intent(in) :: use_precipitation_fluxes
    logical,intent(in) :: Lcalipso, Lcloudsat, Lmodis
!<-added
  type(radar_cfg) :: &
       rcfg_cloudsat     ! Radar configuration

    ! Outputs
    type(cosp_optical_inputs),intent(inout) :: cospIN
    type(cosp_column_inputs),intent(inout)  :: cospstateIN

    ! Local variables
    type(rng_state),allocatable,dimension(:) :: rngs  ! Seeds for random number generator
    integer,dimension(:),allocatable :: seed
    integer,dimension(:),allocatable :: cloudsat_preclvl_index !<-added
    integer :: i,j,k,tp
    real(wp) :: zstep !<-added
    real(wp),dimension(:,:), allocatable :: &
         ls_p_rate, cv_p_rate, frac_ls, frac_cv, prec_ls, prec_cv,g_vol
    real(wp),dimension(:,:,:),  allocatable :: &
         frac_prec, MODIS_cloudWater, MODIS_cloudIce, fracPrecipIce, fracPrecipIce_statGrid, & !<-added
         MODIS_watersize,MODIS_iceSize, MODIS_opticalThicknessLiq,MODIS_opticalThicknessIce
    real(wp),dimension(:,:,:,:),allocatable :: &
         mr_hydro, Reff, Np
    real(wp),dimension(nPoints,nLevels) :: &
         column_frac_out, column_prec_out, fl_lsrain, fl_lssnow, fl_lsgrpl, fl_ccrain, fl_ccsnow
    logical :: cmpGases=.true.
    real(wp),dimension(:),allocatable :: mmar
    logical,dimension(:),allocatable :: mm

    if (Ncolumns .gt. 1) then
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Generate subcolumns for clouds (SCOPS) and precipitation type (PREC_SCOPS)
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! RNG used for subcolumng eneration
       allocate(rngs(nPoints),seed(nPoints))
       seed(:)=0
       seed = int(cospstateIN%phalf(:,Nlevels+1))  ! In case of NPoints=1
       ! *NOTE* Chunking will change the seed
       if (NPoints .gt. 1) seed=int((cospstateIN%phalf(:,Nlevels+1)-minval(cospstateIN%phalf(:,Nlevels+1)))/      &
            (maxval(cospstateIN%phalf(:,Nlevels+1))-minval(cospstateIN%phalf(:,Nlevels+1)))*100000) + 1
       call init_rng(rngs, seed)
      
       ! Call scops
       call scops(NPoints,Nlevels,Ncolumns,rngs,tca,cca,overlap,cospIN%frac_out,0)
       deallocate(seed,rngs)
       
       ! Sum up precipitation rates
       allocate(ls_p_rate(nPoints,nLevels),cv_p_rate(nPoints,Nlevels))
       if(use_precipitation_fluxes) then
          ls_p_rate(:,1:nLevels) = fl_lsrainIN + fl_lssnowIN + fl_lsgrplIN
          cv_p_rate(:,1:nLevels) = fl_ccrainIN + fl_ccsnowIN
       else
          ls_p_rate(:,1:nLevels) = 0 ! mixing_ratio(rain) + mixing_ratio(snow) + mixing_ratio (groupel)
          cv_p_rate(:,1:nLevels) = 0 ! mixing_ratio(rain) + mixing_ratio(snow)
       endif
       
       ! Call PREC_SCOPS
       allocate(frac_prec(nPoints,nColumns,nLevels))
       call prec_scops(nPoints,nLevels,nColumns,ls_p_rate,cv_p_rate,cospIN%frac_out,frac_prec)
       cospIN%frac_prec = frac_prec
       deallocate(ls_p_rate,cv_p_rate)

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Compute fraction in each gridbox for precipitation  and cloud type.
       !%%%%%%%%%%%%%%%%0%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Allocate
       allocate(frac_ls(nPoints,nLevels),prec_ls(nPoints,nLevels),                       &
                frac_cv(nPoints,nLevels),prec_cv(nPoints,nLevels))
       
       ! Initialize
       frac_ls(1:nPoints,1:nLevels) = 0._wp
       prec_ls(1:nPoints,1:nLevels) = 0._wp
       frac_cv(1:nPoints,1:nLevels) = 0._wp
       prec_cv(1:nPoints,1:nLevels) = 0._wp
!       do j=1,nPoints
!          do k=1,nLevels
!             do i=1,nColumns
!                if (cospIN%frac_out(j,i,k)  .eq. 1)  frac_ls(j,k) = frac_ls(j,k)+1._wp
!                if (cospIN%frac_out(j,i,k)  .eq. 2)  frac_cv(j,k) = frac_cv(j,k)+1._wp
!                if (frac_prec(j,i,k) .eq. 1)  prec_ls(j,k) = prec_ls(j,k)+1._wp
!                if (frac_prec(j,i,k) .eq. 2)  prec_cv(j,k) = prec_cv(j,k)+1._wp
!                if (frac_prec(j,i,k) .eq. 3)  prec_cv(j,k) = prec_cv(j,k)+1._wp
!                if (frac_prec(j,i,k) .eq. 3)  prec_ls(j,k) = prec_ls(j,k)+1._wp
!             enddo
!             frac_ls(j,k)=frac_ls(j,k)/nColumns
!             frac_cv(j,k)=frac_cv(j,k)/nColumns
!             prec_ls(j,k)=prec_ls(j,k)/nColumns
!             prec_cv(j,k)=prec_cv(j,k)/nColumns
!          enddo
!       enddo       
       allocate(mm(nColumns))
       allocate(mmar(nColumns))
       mmar(:) = 1._wp
       do j=1,nPoints
          do k=1,nLevels
             mm=(cospIN%frac_out(j,:,k) .eq. 1)
             frac_ls(j,k) = sum(mmar,mask=mm)/nColumns
             mm=(cospIN%frac_out(j,:,k) .eq. 2)
             frac_cv(j,k) = sum(mmar,mask=mm)/nColumns
             mm=(frac_prec(j,:,k) .eq. 1 .or. frac_prec(j,:,k) .eq. 3)
             prec_ls(j,k) = sum(mmar,mask=mm)/nColumns
             mm=(frac_prec(j,:,k) .eq. 2 .or. frac_prec(j,:,k) .eq. 3)
             prec_cv(j,k) = sum(mmar,mask=mm)/nColumns
          enddo
       enddo
       deallocate(mm)
       deallocate(mmar)

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Assign gridmean mixing-ratios (mr_XXXXX), effective radius (ReffIN) and number
       ! concentration (not defined) to appropriate sub-column. Here we are using scops. 
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       allocate(mr_hydro(nPoints,nColumns,nLevels,nHydro),                               &
                Reff(nPoints,nColumns,nLevels,nHydro),                                   &
                Np(nPoints,nColumns,nLevels,nHydro))

       ! Initialize
       mr_hydro(:,:,:,:) = 0._wp
       Reff(:,:,:,:)     = 0._wp
       Np(:,:,:,:)       = 0._wp
       do k=1,nColumns
          !@@@@@ Subcolumn cloud fraction
          column_frac_out = cospIN%frac_out(:,k,:)

          ! LS clouds
          where (column_frac_out == I_LSC)
             mr_hydro(:,k,:,I_LSCLIQ) = mr_lsliq
             mr_hydro(:,k,:,I_LSCICE) = mr_lsice
             Reff(:,k,:,I_LSCLIQ)     = ReffIN(:,:,I_LSCLIQ)
             Reff(:,k,:,I_LSCICE)     = ReffIN(:,:,I_LSCICE)
          ! CONV clouds   
          elsewhere (column_frac_out == I_CVC)
             mr_hydro(:,k,:,I_CVCLIQ) = mr_ccliq
             mr_hydro(:,k,:,I_CVCICE) = mr_ccice
             Reff(:,k,:,I_CVCLIQ)     = ReffIN(:,:,I_CVCLIQ)
             Reff(:,k,:,I_CVCICE)     = ReffIN(:,:,I_CVCICE)
          end where
          
          !@@@@@ Subcolumn precipitation
          column_prec_out = frac_prec(:,k,:)

          ! LS Precipitation
          where ((column_prec_out == 1) .or. (column_prec_out == 3))
             Reff(:,k,:,I_LSRAIN) = ReffIN(:,:,I_LSRAIN)
             Reff(:,k,:,I_LSSNOW) = ReffIN(:,:,I_LSSNOW)
             Reff(:,k,:,I_LSGRPL) = ReffIN(:,:,I_LSGRPL)
             ! CONV precipitation   
          elsewhere ((column_prec_out == 2) .or. (column_prec_out == 3))
             Reff(:,k,:,I_CVRAIN) = ReffIN(:,:,I_CVRAIN)
             Reff(:,k,:,I_CVSNOW) = ReffIN(:,:,I_CVSNOW)
          end where
       enddo
       
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Convert the subcolumn mixing ratio and precipitation fluxes from gridbox mean
       ! values to fraction-based values. 
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Initialize
       fl_lsrain(:,:) = 0._wp
       fl_lssnow(:,:) = 0._wp
       fl_lsgrpl(:,:) = 0._wp
       fl_ccrain(:,:) = 0._wp
       fl_ccsnow(:,:) = 0._wp
       do k=1,nLevels
          do j=1,nPoints
             ! In-cloud mixing ratios.
             if (frac_ls(j,k) .ne. 0.) then
                mr_hydro(j,:,k,I_LSCLIQ) = mr_hydro(j,:,k,I_LSCLIQ)/frac_ls(j,k)
                mr_hydro(j,:,k,I_LSCICE) = mr_hydro(j,:,k,I_LSCICE)/frac_ls(j,k)
             endif
             if (frac_cv(j,k) .ne. 0.) then
                mr_hydro(j,:,k,I_CVCLIQ) = mr_hydro(j,:,k,I_CVCLIQ)/frac_cv(j,k)
                mr_hydro(j,:,k,I_CVCICE) = mr_hydro(j,:,k,I_CVCICE)/frac_cv(j,k)
             endif
             ! Precipitation
             if (use_precipitation_fluxes) then
                if (prec_ls(j,k) .ne. 0.) then
                   fl_lsrain(j,k) = fl_lsrainIN(j,k)/prec_ls(j,k)
                   fl_lssnow(j,k) = fl_lssnowIN(j,k)/prec_ls(j,k)
                   fl_lsgrpl(j,k) = fl_lsgrplIN(j,k)/prec_ls(j,k)
                endif
                if (prec_cv(j,k) .ne. 0.) then
                   fl_ccrain(j,k) = fl_ccrainIN(j,k)/prec_cv(j,k)
                   fl_ccsnow(j,k) = fl_ccsnowIN(j,k)/prec_cv(j,k)
                endif
             else
                if (prec_ls(j,k) .ne. 0.) then
                   mr_hydro(j,:,k,I_LSRAIN) = mr_hydro(j,:,k,I_LSRAIN)/prec_ls(j,k)
                   mr_hydro(j,:,k,I_LSSNOW) = mr_hydro(j,:,k,I_LSSNOW)/prec_ls(j,k)
                   mr_hydro(j,:,k,I_LSGRPL) = mr_hydro(j,:,k,I_LSGRPL)/prec_ls(j,k)
                endif
                if (prec_cv(j,k) .ne. 0.) then
                   mr_hydro(j,:,k,I_CVRAIN) = mr_hydro(j,:,k,I_CVRAIN)/prec_cv(j,k)
                   mr_hydro(j,:,k,I_CVSNOW) = mr_hydro(j,:,k,I_CVSNOW)/prec_cv(j,k)
                endif
             endif
          enddo
       enddo

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Convert precipitation fluxes to mixing ratios
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       if (use_precipitation_fluxes) then
          ! LS rain
          ! call cosp_precip_mxratio(nPoints, nLevels, nColumns, cospstateIN%pfull,        &
          !      cospstateIN%at, frac_prec, 1._wp, n_ax(I_LSRAIN), n_bx(I_LSRAIN),         &
          !      alpha_x(I_LSRAIN), c_x(I_LSRAIN),   d_x(I_LSRAIN),   g_x(I_LSRAIN),       &
          !      a_x(I_LSRAIN),   b_x(I_LSRAIN),   gamma_1(I_LSRAIN), gamma_2(I_LSRAIN),   &
          !      gamma_3(I_LSRAIN), gamma_4(I_LSRAIN), fl_lsrain,                          &
          !      mr_hydro(:,:,:,I_LSRAIN), Reff(:,:,:,I_LSRAIN))
          ! ! LS snow
          ! call cosp_precip_mxratio(nPoints, nLevels, nColumns, cospstateIN%pfull,        &
          !      cospstateIN%at, frac_prec, 1._wp,  n_ax(I_LSSNOW),  n_bx(I_LSSNOW),       &
          !      alpha_x(I_LSSNOW), c_x(I_LSSNOW),  d_x(I_LSSNOW),  g_x(I_LSSNOW),         &
          !      a_x(I_LSSNOW),   b_x(I_LSSNOW),   gamma_1(I_LSSNOW),  gamma_2(I_LSSNOW),  &
          !      gamma_3(I_LSSNOW), gamma_4(I_LSSNOW), fl_lssnow,                          &
          !      mr_hydro(:,:,:,I_LSSNOW), Reff(:,:,:,I_LSSNOW))
          ! ! CV rain
          ! call cosp_precip_mxratio(nPoints, nLevels, nColumns, cospstateIN%pfull,        &
          !      cospstateIN%at, frac_prec, 2._wp, n_ax(I_CVRAIN),  n_bx(I_CVRAIN),        &
          !      alpha_x(I_CVRAIN), c_x(I_CVRAIN),   d_x(I_CVRAIN),   g_x(I_CVRAIN),       &
          !      a_x(I_CVRAIN),   b_x(I_CVRAIN),   gamma_1(I_CVRAIN), gamma_2(I_CVRAIN),   &
          !      gamma_3(I_CVRAIN), gamma_4(I_CVRAIN), fl_ccrain,                          &
          !      mr_hydro(:,:,:,I_CVRAIN), Reff(:,:,:,I_CVRAIN))
          ! ! CV snow
          ! call cosp_precip_mxratio(nPoints, nLevels, nColumns, cospstateIN%pfull,        &
          !      cospstateIN%at, frac_prec, 2._wp, n_ax(I_CVSNOW),  n_bx(I_CVSNOW),        &
          !      alpha_x(I_CVSNOW),  c_x(I_CVSNOW),   d_x(I_CVSNOW),   g_x(I_CVSNOW),      &
          !      a_x(I_CVSNOW),   b_x(I_CVSNOW),   gamma_1(I_CVSNOW), gamma_2(I_CVSNOW),   &
          !      gamma_3(I_CVSNOW), gamma_4(I_CVSNOW), fl_ccsnow,                          &
          !      mr_hydro(:,:,:,I_CVSNOW), Reff(:,:,:,I_CVSNOW))
          ! ! LS groupel.
          ! call cosp_precip_mxratio(nPoints, nLevels, nColumns, cospstateIN%pfull,        &
          !      cospstateIN%at, frac_prec, 1._wp, n_ax(I_LSGRPL),  n_bx(I_LSGRPL),        &
          !      alpha_x(I_LSGRPL), c_x(I_LSGRPL),   d_x(I_LSGRPL),   g_x(I_LSGRPL),       &
          !      a_x(I_LSGRPL),   b_x(I_LSGRPL),   gamma_1(I_LSGRPL),  gamma_2(I_LSGRPL),  &
          !      gamma_3(I_LSGRPL), gamma_4(I_LSGRPL), fl_lsgrpl,                          &
          !      mr_hydro(:,:,:,I_LSGRPL), Reff(:,:,:,I_LSGRPL))

          ! modified by YN: cosp_precip_mxratio
          ! CV rain
          call cosp_precip_mxratio(Npoints, Nlevels, Ncolumns, &
               cospstateIN%pfull, cospstateIN%at, frac_prec, 2._wp, I_CVRAIN, sd, &
               fl_ccrain, mr_hydro(:,:,:,I_CVRAIN), Reff(:,:,:,I_CVRAIN) )

          ! CV snow
          call cosp_precip_mxratio(Npoints, Nlevels, Ncolumns, &
               cospstateIN%pfull, cospstateIN%at, frac_prec, 2._wp, I_CVSNOW, sd, &
               fl_ccsnow, mr_hydro(:,:,:,I_CVSNOW), Reff(:,:,:,I_CVSNOW) )

          deallocate(frac_prec)

!!!  T.Michibata (2019.03.18) added --->
!!!  if you would like to use model-native Qr and Qs prognosed in the CHIMERRA for inputs into CloudSat simulator, 
!!!  validate the following description (strongly recommended):
#ifdef OPT_CHIMERRA
          do k=1,nColumns
             column_prec_out = cospIN%frac_prec(:,k,:)

             where ((column_prec_out == 1) .or. (column_prec_out == 3))
                ! LS rain
                mr_hydro(:,k,:,I_LSRAIN) = mr_lsrain(:,:)/prec_ls(:,:)
                ! LS snow
                mr_hydro(:,k,:,I_LSSNOW) = mr_lssnow(:,:)/prec_ls(:,:)
                ! LS grpl
                mr_hydro(:,k,:,I_LSGRPL) = mr_lsgrpl(:,:)/prec_ls(:,:)
             end where
          enddo
#endif
!!!  <--
       endif
       deallocate(frac_ls,prec_ls,frac_cv,prec_cv)

    else
       cospIN%frac_out(:,:,:) = 1  
       allocate(mr_hydro(nPoints,1,nLevels,nHydro),Reff(nPoints,1,nLevels,nHydro),       &
                Np(nPoints,1,nLevels,nHydro))
       mr_hydro(:,1,:,I_LSCLIQ) = mr_lsliq
       mr_hydro(:,1,:,I_LSCICE) = mr_lsice
       mr_hydro(:,1,:,I_CVCLIQ) = mr_ccliq
       mr_hydro(:,1,:,I_CVCICE) = mr_ccice
       Reff(:,1,:,:)            = ReffIN
    endif
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 11 micron emissivity
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,dem_c,dem_s,       &
                               cospIN%emiss_11)
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 0.67 micron optical depth
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,dtau_c,dtau_s,     &
                               cospIN%tau_067)
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! LIDAR Polarized optics
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (Lcalipso) then
!       call lidar_optics(nPoints, nColumns, nLevels, 4, lidar_ice_type,                       &
!            mr_hydro(:,:,:,I_LSCLIQ),  mr_hydro(:,:,:,I_LSCICE), mr_hydro(:,:,:,I_CVCLIQ),    &
!            mr_hydro(:,:,:,I_CVCICE), ReffIN(:,:,I_LSCLIQ), ReffIN(:,:,I_LSCICE),             &
!            ReffIN(:,:,I_CVCLIQ), ReffIN(:,:,I_CVCICE), cospstateIN%pfull, cospstateIN%phalf, &
!            cospstateIN%at, cospIN%beta_mol, cospIN%betatot, cospIN%tau_mol, cospIN%tautot,   &
!            cospIN%tautot_S_liq, cospIN%tautot_S_ice, cospIN%betatot_ice, cospIN%betatot_liq, &
!            cospIN%tautot_ice, cospIN%tautot_liq)
!--> for snow effects (T.Michibata, 2021/01/03 added)
       call lidar_optics(nPoints, nColumns, nLevels, 5, lidar_ice_type,                       &
            mr_hydro(:,:,:,I_LSCLIQ),  mr_hydro(:,:,:,I_LSCICE), mr_hydro(:,:,:,I_CVCLIQ),    &
            mr_hydro(:,:,:,I_CVCICE), ReffIN(:,:,I_LSCLIQ), ReffIN(:,:,I_LSCICE),             &
            ReffIN(:,:,I_CVCLIQ), ReffIN(:,:,I_CVCICE),                                       &
            mr_hydro(:,:,:,I_LSSNOW), ReffIN(:,:,I_LSSNOW),                                   & ! added
            cospstateIN%pfull, cospstateIN%phalf,                                             &
            cospstateIN%at, cospIN%beta_mol, cospIN%betatot, cospIN%tau_mol, cospIN%tautot,   &
            cospIN%tautot_S_liq, cospIN%tautot_S_ice, cospIN%betatot_ice, cospIN%betatot_liq, &
            cospIN%tautot_ice, cospIN%tautot_liq)
!<-- for snow effects (T.Michibata, 2021/01/03 added)
    endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! CLOUDSAT RADAR OPTICS
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (Lcloudsat) then

       ! Compute gaseous absorption (assume identical for each subcolun) !<-added(a lot)
       !allocate(g_vol(nPoints,nLevels))
       !g_vol(:,:)=0._wp
       if (allocated(g_vol)) then
          continue
       else
         allocate(g_vol(nPoints,nLevels))
       endif
       g_vol(:,:)=0._wp
       do i=1,nPoints
          do j=1,nLevels
             if (rcfg_cloudsat%use_gas_abs == 1 .or. (rcfg_cloudsat%use_gas_abs == 2 .and. j .eq. 1)) then
                g_vol(i,j) = gases(cospstateIN%pfull(i,j), cospstateIN%at(i,j),cospstateIN%qv(i,j),rcfg_cloudsat%freq)
             endif
             cospIN%g_vol_cloudsat(i,:,j)=g_vol(i,j)
          end do
       end do

       ! Loop over all subcolumns
       !allocate(fracPrecipIce(nPoints,nColumns,nLevels))
       !fracPrecipIce(:,:,:) = 0._wp
       if (allocated(fracPrecipIce)) then
          continue
       else
         allocate(fracPrecipIce(nPoints,nColumns,nLevels))
       endif
       fracPrecipIce(:,:,:) = 0._wp
       do k=1,nColumns
          ! Try to avoid memory copy... Pass mr_hydro, Reff and Np directly into quickbeam_optics
          if (k .eq. 1) then
             cmpGases = .true.
             !call quickbeam_optics(sd, rcfg_cloudsat, nPoints, nLevels, R_UNDEF,  &
             call quickbeam_optics(sd, cospIN%rcfg_cloudsat, nPoints, nLevels, R_UNDEF,  &
                  mr_hydro(:,k,:,1:nHydro)*1000._wp, Reff(:,k,:,1:nHydro)*1.e6_wp,&
                  Np(:,k,:,1:nHydro), cospstateIN%pfull, cospstateIN%at,          &
                  cospstateIN%qv, cmpGases, cospIN%z_vol_cloudsat(1:nPoints,k,:), &
                  cospIN%kr_vol_cloudsat(1:nPoints,k,:),                          &
                  cospIN%g_vol_cloudsat(1:nPoints,k,:),                           &
#ifdef OPT_DPLRW
                  cospIN%vfall(1:nPoints,k,:,1:nHydro),                      &
                  cospIN%vfsqu(1:nPoints,k,:,1:nHydro),                      &
                  cospIN%zehyd(1:nPoints,k,:,1:nHydro),                      &
#endif
                  g_vol_out=g_vol)
          else
             cmpGases = .false.
             !call quickbeam_optics(sd, rcfg_cloudsat, nPoints, nLevels, R_UNDEF,  &
             call quickbeam_optics(sd, cospIN%rcfg_cloudsat, nPoints, nLevels, R_UNDEF,  &
                  mr_hydro(:,k,:,1:nHydro)*1000._wp, Reff(:,k,:,1:nHydro)*1.e6_wp,&
                  Np(:,k,:,1:nHydro), cospstateIN%pfull, cospstateIN%at,          &
                  cospstateIN%qv, cmpGases, cospIN%z_vol_cloudsat(1:nPoints,k,:), &
                  cospIN%kr_vol_cloudsat(1:nPoints,k,:),                          &
                  cospIN%g_vol_cloudsat(1:nPoints,k,:),                           &
#ifdef OPT_DPLRW
                  cospIN%vfall(1:nPoints,k,:,1:nHydro),                      &
                  cospIN%vfsqu(1:nPoints,k,:,1:nHydro),                      &
                  cospIN%zehyd(1:nPoints,k,:,1:nHydro),                      &
#endif
                  g_vol_in=g_vol)
          endif

          ! At each model level, what fraction of the precipitation is frozen?
          where(mr_hydro(:,k,:,I_LSRAIN) .gt. 0 .or. mr_hydro(:,k,:,I_LSSNOW) .gt. 0 .or. &
                mr_hydro(:,k,:,I_CVRAIN) .gt. 0 .or. mr_hydro(:,k,:,I_CVSNOW) .gt. 0 .or. &
                mr_hydro(:,k,:,I_LSGRPL) .gt. 0)
             fracPrecipIce(:,k,:) = (mr_hydro(:,k,:,I_LSSNOW) + mr_hydro(:,k,:,I_CVSNOW) + &
                  mr_hydro(:,k,:,I_LSGRPL)) / &
                  (mr_hydro(:,k,:,I_LSSNOW) + mr_hydro(:,k,:,I_CVSNOW) + mr_hydro(:,k,:,I_LSGRPL) + &
                  mr_hydro(:,k,:,I_LSRAIN)  + mr_hydro(:,k,:,I_CVRAIN))
          elsewhere
             fracPrecipIce(:,k,:) = 0._wp
          endwhere
       enddo  ! do k=1,nColumn

       ! Regrid frozen fraction to Cloudsat/Calipso statistical grid
       if (allocated(fracPrecipIce_statGrid)) then
          continue
       else
         allocate(fracPrecipIce_statGrid(nPoints,nColumns,Nlvgrid))
       endif
       fracPrecipIce_statGrid(:,:,:) = 0._wp

       ! Find proper layer above de surface elevation to compute precip flags in Cloudsat/Calipso statistical grid
       if (allocated(cloudsat_preclvl_index)) then
          continue
       else
         allocate(cloudsat_preclvl_index(nPoints))
       endif

       if (use_vgrid) then  ! DEBUG OUT by YN 221213
          call cosp_change_vertical_grid(Npoints, Ncolumns, Nlevels, cospstateIN%hgt_matrix(:,Nlevels:1:-1), &
               cospstateIN%hgt_matrix_half(:,Nlevels:1:-1), fracPrecipIce(:,:,Nlevels:1:-1), Nlvgrid,  &
               vgrid_zl(Nlvgrid:1:-1), vgrid_zu(Nlvgrid:1:-1), fracPrecipIce_statGrid(:,:,Nlvgrid:1:-1))

          cloudsat_preclvl_index(:) = 0._wp
          ! Compute the zstep distance between two atmopsheric layers
          zstep = vgrid_zl(1)-vgrid_zl(2)
          ! Computing altitude index for precip flags calculation (one layer above surfelev layer)
          cloudsat_preclvl_index(:) = cloudsat_preclvl - floor( cospstateIN%surfelev(:)/zstep )
       else
          cloudsat_preclvl_index(:) = nLevels - 1
       endif

       ! For near-surface diagnostics, we only need the frozen fraction at one layer.
        do i=1,nPoints
           cospIN%fracPrecipIce(i,:) = fracPrecipIce_statGrid(i,:,cloudsat_preclvl_index(i))
        enddo

        deallocate(fracPrecipIce,fracPrecipIce_statGrid,cloudsat_preclvl_index)  ! added by YN
     endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! MODIS optics
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (Lmodis) then
       allocate(MODIS_cloudWater(nPoints,nColumns,nLevels),                                &
                MODIS_cloudIce(nPoints,nColumns,nLevels),                                  &
                MODIS_waterSize(nPoints,nColumns,nLevels),                                 &
                MODIS_iceSize(nPoints,nColumns,nLevels),                                   &
                MODIS_opticalThicknessLiq(nPoints,nColumns,nLevels),                       &
                MODIS_opticalThicknessIce(nPoints,nColumns,nLevels))
       ! Cloud water
       call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,                &
            mr_hydro(:,:,:,I_CVCLIQ),mr_hydro(:,:,:,I_LSCLIQ),MODIS_cloudWater)
       ! Cloud ice
       call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,                &
            mr_hydro(:,:,:,I_CVCICE),mr_hydro(:,:,:,I_LSCICE),MODIS_cloudIce)  
       ! Water droplet size
       call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,                &
            Reff(:,:,:,I_CVCLIQ),Reff(:,:,:,I_LSCLIQ),MODIS_waterSize)
       ! Ice crystal size
       call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,                &
            Reff(:,:,:,I_CVCICE),Reff(:,:,:,I_LSCICE),MODIS_iceSize)
       
       ! Partition optical thickness into liquid and ice parts
       call modis_optics_partition(nPoints, nLevels, nColumns, MODIS_cloudWater,           &
            MODIS_cloudIce, MODIS_waterSize, MODIS_iceSize, cospIN%tau_067,                &
            MODIS_opticalThicknessLiq, MODIS_opticalThicknessIce)
       
       ! Compute assymetry parameter and single scattering albedo 
       call modis_optics(nPoints, nLevels, nColumns, MODIS_opticalThicknessLiq,            &
            MODIS_waterSize*1.0e6_wp, MODIS_opticalThicknessIce,                           &
            MODIS_iceSize*1.0e6_wp, cospIN%fracLiq, cospIN%asym, cospIN%ss_alb)

       ! Deallocate memory     ! debug modified 2022AUG by YN
       deallocate(MODIS_cloudWater,MODIS_cloudIce,MODIS_WaterSize,MODIS_iceSize,           &
            MODIS_opticalThicknessLiq,MODIS_opticalThicknessIce)
    end if
    ! debug modified 2022AUG by YN
    deallocate(mr_hydro,Np,Reff)

  end subroutine subsample_and_optics
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE construct_cospIN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine construct_cospIN(npoints,ncolumns,nlevels, &
             Lmodis, Lmisr, Lisccp, Lcalipso, Lcloudsat, y)
  USE mod_cosp, ONLY: cosp_optical_inputs
  USE MOD_COSP_CONFIG,     ONLY: PARASOL_NREFL
#ifdef OPT_DPLRW
  USE MOD_COSP_CONFIG,     ONLY: N_HYDRO
#endif

    ! Inputs
    integer,intent(in) :: &
         npoints,  & ! Number of horizontal gridpoints
         ncolumns, & ! Number of subcolumns
         nlevels     ! Number of vertical levels
    logical :: &
         Lmodis    , & 
         Lmisr     , &
         Lisccp    , & 
         Lcalipso  , &
         Lcloudsat 
    ! Outputs 
    type(cosp_optical_inputs),intent(out) :: y
    
    ! Dimensions
    y%Npoints  = Npoints
    y%Ncolumns = Ncolumns
    y%Nlevels  = Nlevels
!    y%Npart    = 4
    y%Npart    = 5  ! added for snow effects (T.Michibata)
    y%Nrefl    = PARASOL_NREFL
    allocate(y%frac_out (npoints,       ncolumns,nlevels))
    allocate(y%frac_prec(npoints,       ncolumns,nlevels)) ! added by YN

    if (Lmodis .or. Lmisr .or. Lisccp) then
       allocate(y%tau_067(npoints,        ncolumns,nlevels),&
                y%emiss_11(npoints,       ncolumns,nlevels))
    endif
    if (Lcalipso) then
       allocate(y%betatot(npoints,        ncolumns,nlevels),&
                y%betatot_ice(npoints,    ncolumns,nlevels),&
                y%betatot_liq(npoints,    ncolumns,nlevels),&
                y%tautot(npoints,         ncolumns,nlevels),&
                y%tautot_ice(npoints,     ncolumns,nlevels),&
                y%tautot_liq(npoints,     ncolumns,nlevels),&
                y%beta_mol(npoints,                nlevels),&
                y%tau_mol(npoints,                 nlevels),&
                y%tautot_S_ice(npoints,   ncolumns        ),&
                y%tautot_S_liq(npoints,   ncolumns))

    endif
    if (Lcloudsat) then
       allocate(y%z_vol_cloudsat(npoints, ncolumns,nlevels),&
                y%kr_vol_cloudsat(npoints,ncolumns,nlevels),&
                y%g_vol_cloudsat(npoints, ncolumns,nlevels),&
                y%fracPrecipIce(npoints,  ncolumns)         )  !<-added
#ifdef OPT_DPLRW
       allocate(y%vfall(npoints, ncolumns, nlevels, N_HYDRO),  &
                y%vfsqu(npoints, ncolumns, nlevels, N_HYDRO),  &
                y%zehyd(npoints, ncolumns, nlevels, N_HYDRO) )
#endif
    endif
    if (Lmodis) then
       allocate(y%fracLiq(npoints,        ncolumns,nlevels),&
                y%asym(npoints,           ncolumns,nlevels),&
                y%ss_alb(npoints,         ncolumns,nlevels))
    endif
    
  end subroutine construct_cospIN
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE construct_cospstateIN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
  subroutine construct_cospstateIN(npoints,nlevels,nchan,y)
  USE mod_cosp, ONLY: cosp_column_inputs
#ifdef OPT_DPLRW
  USE MOD_COSP_CONFIG,     ONLY: N_HYDRO
#endif
    ! Inputs
    integer,intent(in) :: &
         npoints, & ! Number of horizontal gridpoints
         nlevels, & ! Number of vertical levels
         nchan      ! Number of channels
    ! Outputs
    type(cosp_column_inputs),intent(out) :: y         
    
    allocate(y%sunlit(npoints),y%skt(npoints),y%land(npoints),y%at(npoints,nlevels),     &
             y%pfull(npoints,nlevels),y%phalf(npoints,nlevels+1),y%qv(npoints,nlevels),  &
             y%o3(npoints,nlevels),y%hgt_matrix(npoints,nlevels),y%u_sfc(npoints),       &
             y%v_sfc(npoints),y%lat(npoints),y%lon(nPoints),y%emis_sfc(nchan),           &
             y%cloudIce(nPoints,nLevels),y%cloudLiq(nPoints,nLevels),y%surfelev(npoints),&
             y%fl_snow(nPoints,nLevels),y%fl_rain(nPoints,nLevels),y%seaice(npoints),    &
             y%tca(nPoints,nLevels),y%hgt_matrix_half(npoints,nlevels+1)                 &
#ifdef OPT_DPLRW
             , y%gwvel(npoints,nlevels), y%gcumf(npoints,nlevels+1)                      &
#endif
             )

  end subroutine construct_cospstateIN

  ! ######################################################################################
  ! SUBROUTINE construct_cosp_outputs
  !
  ! This subroutine allocates output fields based on input logical flag switches.
  ! ######################################################################################  
  subroutine construct_cosp_outputs(Lpctisccp,Lclisccp,&
                                    Lboxptopisccp,Lboxtauisccp,Ltauisccp,Lcltisccp,      &
                                    Lmeantbisccp,Lmeantbclrisccp,Lalbisccp,LclMISR,      &
                                    Lcltmodis,Lclwmodis,Lclimodis,Lclhmodis,Lclmmodis,   &
                                    Lcllmodis,Ltautmodis,Ltauwmodis,Ltauimodis,          &
                                    Ltautlogmodis,Ltauwlogmodis,Ltauilogmodis,           &
                                    Lreffclwmodis,Lreffclimodis,Lpctmodis,Llwpmodis,     &
                                    Liwpmodis,Lclmodis,Latb532,LlidarBetaMol532,         &
                                    LcfadLidarsr532,Lclcalipso2,                         &
                                    Lclcalipso,Lclhcalipso,Lcllcalipso,Lclmcalipso,      &
                                    Lcltcalipso,Lcltlidarradar,Lclcalipsoliq,            &
                                    Lclcalipsoice,Lclcalipsoun,Lclcalipsotmp,            &
                                    Lclcalipsotmpliq,Lclcalipsotmpice,Lclcalipsotmpun,   &
                                    Lcltcalipsoliq,Lcltcalipsoice,Lcltcalipsoun,         &
                                    Lclhcalipsoliq,Lclhcalipsoice,Lclhcalipsoun,         &
                                    Lclmcalipsoliq,Lclmcalipsoice,Lclmcalipsoun,         &
                                    Lcllcalipsoliq,Lcllcalipsoice,Lcllcalipsoun,         & 
                                    LcfadDbze94,Ldbze94,Lparasolrefl,Ltbrttov,           &
                                    Lptradarflag0,Lptradarflag1,Lptradarflag2,           & !<-added
                                    Lptradarflag3,Lptradarflag4,Lptradarflag5,           & !<-added
                                    Lptradarflag6,Lptradarflag7,Lptradarflag8,           & !<-added
                                    Lptradarflag9,Lradarpia,Lwr_occfreq, Lcfodd,         & !<-added
                                    Lscops,                                              & !<-added by YN
#ifdef OPT_DPLRW
                                    Ldplrw, &
#endif
                                    Npoints,Ncolumns,Nlevels,Nchan,x)

  USE mod_cosp, ONLY: cosp_outputs
  USE MOD_COSP_CONFIG,     ONLY: R_UNDEF,PARASOL_NREFL,LIDAR_NCAT,SR_BINS,                &
                                 N_HYDRO,RTTOV_MAX_CHANNELS,numMISRHgtBins,               &
                                 DBZE_BINS,LIDAR_NTEMP,calipso_histBsct,                  &
                                 WR_NREGIME,    CFODD_NCLASS,                             &
                                 CFODD_NDBZE,   CFODD_NICOD,                              &
                                 numMODISTauBins,numMODISPresBins,                        &
                                 numMODISReffIceBins,numMODISReffLiqBins,                 &
                                 numISCCPTauBins,numISCCPPresBins,numMISRTauBins,         &
                                 ntau,npres,modis_histTau,tau_binBounds,                  & ! <- modified by Doppler
                                 modis_histTauEdges,tau_binEdges,                         &
                                 modis_histTauCenters,tau_binCenters,ntauV1p4,            &
                                 tau_binBoundsV1p4,tau_binEdgesV1p4, tau_binCentersV1p4,  &
                                 Nlvgrid
#ifdef OPT_DPLRW
  use mod_cosp_config,     only: Nlvdplr, Nlvspwd, NlvdBZe, Nlvtemp, N_ISCCP
#endif

     ! Inputs
     logical,intent(in) :: &
         Lpctisccp,        & ! ISCCP mean cloud top pressure
         Lclisccp,         & ! ISCCP cloud area fraction
         Lboxptopisccp,    & ! ISCCP CTP in each column
         Lboxtauisccp,     & ! ISCCP optical epth in each column
         Ltauisccp,        & ! ISCCP mean optical depth
         Lcltisccp,        & ! ISCCP total cloud fraction
         Lmeantbisccp,     & ! ISCCP mean all-sky 10.5micron brightness temperature
         Lmeantbclrisccp,  & ! ISCCP mean clear-sky 10.5micron brightness temperature
         Lalbisccp,        & ! ISCCP mean cloud albedo         
         LclMISR,          & ! MISR cloud fraction
         Lcltmodis,        & ! MODIS total cloud fraction
         Lclwmodis,        & ! MODIS liquid cloud fraction
         Lclimodis,        & ! MODIS ice cloud fraction
         Lclhmodis,        & ! MODIS high-level cloud fraction
         Lclmmodis,        & ! MODIS mid-level cloud fraction
         Lcllmodis,        & ! MODIS low-level cloud fraction
         Ltautmodis,       & ! MODIS total cloud optical thicknes
         Ltauwmodis,       & ! MODIS liquid optical thickness
         Ltauimodis,       & ! MODIS ice optical thickness
         Ltautlogmodis,    & ! MODIS total cloud optical thickness (log10 mean)
         Ltauwlogmodis,    & ! MODIS liquid optical thickness (log10 mean)
         Ltauilogmodis,    & ! MODIS ice optical thickness (log10 mean)
         Lreffclwmodis,    & ! MODIS liquid cloud particle size
         Lreffclimodis,    & ! MODIS ice particle size
         Lpctmodis,        & ! MODIS cloud top pressure
         Llwpmodis,        & ! MODIS cloud liquid water path
         Liwpmodis,        & ! MODIS cloud ice water path
         Lclmodis,         & ! MODIS cloud area fraction
         Latb532,          & ! CALIPSO attenuated total backscatter (532nm)
         LlidarBetaMol532, & ! CALIPSO molecular backscatter (532nm)         
         LcfadLidarsr532,  & ! CALIPSO scattering ratio CFAD
         Lclcalipso2,      & ! CALIPSO cloud fraction undetected by cloudsat
         Lclcalipso,       & ! CALIPSO cloud area fraction
         Lclhcalipso,      & ! CALIPSO high-level cloud fraction
         Lcllcalipso,      & ! CALIPSO low-level cloud fraction
         Lclmcalipso,      & ! CALIPSO mid-level cloud fraction
         Lcltcalipso,      & ! CALIPSO total cloud fraction
         Lcltlidarradar,   & ! CALIPSO-CLOUDSAT total cloud fraction
         Lclcalipsoliq,    & ! CALIPSO liquid cloud area fraction
         Lclcalipsoice,    & ! CALIPSO ice cloud area fraction 
         Lclcalipsoun,     & ! CALIPSO undetected cloud area fraction
         Lclcalipsotmp,    & ! CALIPSO undetected cloud area fraction
         Lclcalipsotmpliq, & ! CALIPSO liquid cloud area fraction
         Lclcalipsotmpice, & ! CALIPSO ice cloud area fraction
         Lclcalipsotmpun,  & ! CALIPSO undetected cloud area fraction
         Lcltcalipsoliq,   & ! CALIPSO liquid total cloud fraction
         Lcltcalipsoice,   & ! CALIPSO ice total cloud fraction
         Lcltcalipsoun,    & ! CALIPSO undetected total cloud fraction
         Lclhcalipsoliq,   & ! CALIPSO high-level liquid cloud fraction
         Lclhcalipsoice,   & ! CALIPSO high-level ice cloud fraction
         Lclhcalipsoun,    & ! CALIPSO high-level undetected cloud fraction
         Lclmcalipsoliq,   & ! CALIPSO mid-level liquid cloud fraction
         Lclmcalipsoice,   & ! CALIPSO mid-level ice cloud fraction
         Lclmcalipsoun,    & ! CALIPSO mid-level undetected cloud fraction
         Lcllcalipsoliq,   & ! CALIPSO low-level liquid cloud fraction
         Lcllcalipsoice,   & ! CALIPSO low-level ice cloud fraction
         Lcllcalipsoun,    & ! CALIPSO low-level undetected cloud fraction
         LcfadDbze94,      & ! CLOUDSAT radar reflectivity CFAD
         Ldbze94,          & ! CLOUDSAT radar reflectivity
         LparasolRefl,     & ! PARASOL reflectance
         Ltbrttov            ! RTTOV mean clear-sky brightness temperature
     logical,intent(in) :: & !<-added
         Lptradarflag0,    & ! CLOUDSAT 
         Lptradarflag1,    & ! CLOUDSAT 
         Lptradarflag2,    & ! CLOUDSAT 
         Lptradarflag3,    & ! CLOUDSAT 
         Lptradarflag4,    & ! CLOUDSAT 
         Lptradarflag5,    & ! CLOUDSAT 
         Lptradarflag6,    & ! CLOUDSAT 
         Lptradarflag7,    & ! CLOUDSAT 
         Lptradarflag8,    & ! CLOUDSAT 
         Lptradarflag9,    & ! CLOUDSAT
         Lradarpia,        & ! CLOUDSAT
         Lwr_occfreq,      & ! Warmrain Diagnostics from CloudSat+MODIS sims.
         Lcfodd              ! Warmrain Diagnostics from CloudSat+MODIS sims.
     ! added by YN
     logical,intent(in) :: Lscops
#ifdef OPT_DPLRW
     logical,intent(in)    :: Ldplrw              ! doppler velocity
     !logical,intent(inout) :: Ldplrw_LS, Ldplrw_CU   ! output control
     !integer :: tp
#endif

     integer,intent(in) :: &
          Npoints,         & ! Number of sampled points
          Ncolumns,        & ! Number of subgrid columns
          Nlevels,         & ! Number of model levels
!          Nlvgrid,         & ! Number of levels in L3 stats computation
          Nchan              ! Number of RTTOV channels  
          
     ! Outputs
     type(cosp_outputs),intent(out) :: &
          x           ! COSP output structure  
   
     ! ISCCP simulator outputs
    if (Lboxtauisccp)    allocate(x%isccp_boxtau(Npoints,Ncolumns)) 
    if (Lboxptopisccp)   allocate(x%isccp_boxptop(Npoints,Ncolumns))
    if (Lclisccp)        allocate(x%isccp_fq(Npoints,numISCCPTauBins,numISCCPPresBins))
    if (Lcltisccp)       allocate(x%isccp_totalcldarea(Npoints))
    if (Lpctisccp)       allocate(x%isccp_meanptop(Npoints))
    if (Ltauisccp)       allocate(x%isccp_meantaucld(Npoints))
    if (Lmeantbisccp)    allocate(x%isccp_meantb(Npoints))
    if (Lmeantbclrisccp) allocate(x%isccp_meantbclr(Npoints))
    if (Lalbisccp)       allocate(x%isccp_meanalbedocld(Npoints))
    
    ! MISR simulator
    if (LclMISR) then 
       allocate(x%misr_fq(Npoints,numMISRTauBins,numMISRHgtBins))
       ! *NOTE* These 3 fields are not output, but were part of the v1.4.0 cosp_misr, so
       !        they are still computed. Should probably have a logical to control these
       !        outputs.
       allocate(x%misr_dist_model_layertops(Npoints,numMISRHgtBins))
       allocate(x%misr_meanztop(Npoints))
       allocate(x%misr_cldarea(Npoints))    
    endif
    
    ! MODIS simulator
    if (Lcltmodis)     allocate(x%modis_Cloud_Fraction_Total_Mean(Npoints))
    if (Lclwmodis)     allocate(x%modis_Cloud_Fraction_Water_Mean(Npoints))
    if (Lclimodis)     allocate(x%modis_Cloud_Fraction_Ice_Mean(Npoints))
    if (Lclhmodis)     allocate(x%modis_Cloud_Fraction_High_Mean(Npoints))
    if (Lclmmodis)     allocate(x%modis_Cloud_Fraction_Mid_Mean(Npoints))
    if (Lcllmodis)     allocate(x%modis_Cloud_Fraction_Low_Mean(Npoints))
    if (Ltautmodis)    allocate(x%modis_Optical_Thickness_Total_Mean(Npoints))
    if (Ltauwmodis)    allocate(x%modis_Optical_Thickness_Water_Mean(Npoints))
    if (Ltauimodis)    allocate(x%modis_Optical_Thickness_Ice_Mean(Npoints))
    if (Ltautlogmodis) allocate(x%modis_Optical_Thickness_Total_LogMean(Npoints))
    if (Ltauwlogmodis) allocate(x%modis_Optical_Thickness_Water_LogMean(Npoints))
    if (Ltauilogmodis) allocate(x%modis_Optical_Thickness_Ice_LogMean(Npoints))
    if (Lreffclwmodis) allocate(x%modis_Cloud_Particle_Size_Water_Mean(Npoints))
    if (Lreffclimodis) allocate(x%modis_Cloud_Particle_Size_Ice_Mean(Npoints))
    if (Lpctmodis)     allocate(x%modis_Cloud_Top_Pressure_Total_Mean(Npoints))
    if (Llwpmodis)     allocate(x%modis_Liquid_Water_Path_Mean(Npoints))
    if (Liwpmodis)     allocate(x%modis_Ice_Water_Path_Mean(Npoints))
    if (Lclmodis) then
        allocate(x%modis_Optical_Thickness_vs_Cloud_Top_Pressure(nPoints,numModisTauBins,numMODISPresBins))
        allocate(x%modis_Optical_thickness_vs_ReffLIQ(nPoints,numMODISTauBins,numMODISReffLiqBins))   
        allocate(x%modis_Optical_Thickness_vs_ReffICE(nPoints,numMODISTauBins,numMODISReffIceBins))
    endif
    
    ! LIDAR simulator
    if (LlidarBetaMol532) allocate(x%calipso_beta_mol(Npoints,Nlevels))
    if (Latb532)          allocate(x%calipso_beta_tot(Npoints,Ncolumns,Nlevels))
    if (LcfadLidarsr532)  then
        allocate(x%calipso_srbval(SR_BINS+1))
        allocate(x%calipso_cfad_sr(Npoints,SR_BINS,Nlvgrid))
        allocate(x%calipso_betaperp_tot(Npoints,Ncolumns,Nlevels))  
    endif
    if (Lclcalipso)       allocate(x%calipso_lidarcld(Npoints,Nlvgrid))
    if (Lclhcalipso .or. Lclmcalipso .or. Lcllcalipso .or. Lcltcalipso) then
        allocate(x%calipso_cldlayer(Npoints,LIDAR_NCAT))
    endif   
    if (Lclcalipsoice .or. Lclcalipsoliq .or. Lclcalipsoun) then
        allocate(x%calipso_lidarcldphase(Npoints,Nlvgrid,6))
    endif
    if (Lclcalipsotmp .or. Lclcalipsotmpliq .or. Lclcalipsoice .or. Lclcalipsotmpun .or. Lclcalipsotmpice) then
        allocate(x%calipso_lidarcldtmp(Npoints,LIDAR_NTEMP,5))
    endif
    if (Lcllcalipsoice .or. Lclmcalipsoice .or. Lclhcalipsoice .or.                   &
        Lcltcalipsoice .or. Lcllcalipsoliq .or. Lclmcalipsoliq .or.                   &
        Lclhcalipsoliq .or. Lcltcalipsoliq .or. Lcllcalipsoun  .or.                   &
        Lclmcalipsoun  .or. Lclhcalipsoun  .or. Lcltcalipsoun) then
        allocate(x%calipso_cldlayerphase(Npoints,LIDAR_NCAT,6))     
    endif
    
  !  if (Lclopaquecalipso .or. Lclthincalipso .or. Lclzopaquecalipso) then
  !      allocate(x%calipso_cldtype(Npoints,LIDAR_NTYPE))
  !  endif
  !  if (Lclopaquetemp .or. Lclthintemp .or. Lclzopaquetemp) then
  !      allocate(x%calipso_cldtypetemp(Npoints,LIDAR_NTYPE))
  !  endif
  !  if (Lclopaquemeanz .or. Lclthinmeanz) then
  !      allocate(x%calipso_cldtypemeanz(Npoints,2))
  !  endif
  !  if (Lclopaquemeanzse .or. Lclthinmeanzse .or. Lclzopaquecalipsose) then
  !      allocate(x%calipso_cldtypemeanzse(Npoints,3))
  !  endif
  !  if (Lclthinemis) then
  !      allocate(x%calipso_cldthinemis(Npoints))
  !  endif
  !  if (Lclcalipsoopaque .or. Lclcalipsothin .or. Lclcalipsozopaque .or. Lclcalipsoopacity) then
  !      allocate(x%calipso_lidarcldtype(Npoints,Nlvgrid,LIDAR_NTYPE+1))
  !  endif 

    ! These 2 outputs are part of the calipso output type, but are not controlled by an 
    ! logical switch in the output namelist, so if all other fields are on, then allocate
    if (LlidarBetaMol532 .or. Latb532        .or. LcfadLidarsr532 .or. Lclcalipso  .or.  &
        Lclcalipsoice    .or. Lclcalipsoliq  .or. Lclcalipsoun    .or. Lclcalipso2 .or.  &
        Lclhcalipso      .or. Lclmcalipso    .or. Lcllcalipso     .or. Lcltcalipso .or.  &
        Lclcalipsotmp    .or. Lclcalipsoice  .or. Lclcalipsotmpun .or.                   &
        Lclcalipsotmpliq .or. Lcllcalipsoice .or. Lclmcalipsoice  .or.                   &
        Lclhcalipsoice   .or. Lcltcalipsoice .or. Lcllcalipsoliq  .or.                   &
        Lclmcalipsoliq   .or. Lclhcalipsoliq .or. Lcltcalipsoliq  .or.                   &
        Lcllcalipsoun    .or. Lclmcalipsoun  .or. Lclhcalipsoun   .or. Lcltcalipsoun) then
       allocate(x%calipso_tau_tot(Npoints,Ncolumns,Nlevels))       
       allocate(x%calipso_temp_tot(Npoints,Nlevels))               
    endif 
      
    ! PARASOL
    if (Lparasolrefl) then
        allocate(x%parasolPix_refl(Npoints,Ncolumns,PARASOL_NREFL))
        allocate(x%parasolGrid_refl(Npoints,PARASOL_NREFL))
    endif 

    ! Cloudsat simulator(added flags:Y.Imura 2021.06.16)
    if (Ldbze94)        allocate(x%cloudsat_Ze_tot(Npoints,Ncolumns,Nlevels))
    if (LcfadDbze94)    allocate(x%cloudsat_cfad_ze(Npoints,DBZE_BINS,Nlvgrid))
    if (Lptradarflag0 .or. Lptradarflag1 .or. Lptradarflag2 .or. Lptradarflag3 .or. & !<-added
        Lptradarflag4 .or. Lptradarflag5 .or. Lptradarflag6 .or. Lptradarflag7 .or. & !<-added
        Lptradarflag8 .or. Lptradarflag9) then !<-added
       allocate(x%cloudsat_precip_cover(Npoints,DBZE_BINS)) !<-added
    endif 
    if (Lradarpia) allocate(x%cloudsat_pia(Npoints)) !<-added


    ! Combined CALIPSO/CLOUDSAT fields
    if (Lclcalipso2)    allocate(x%lidar_only_freq_cloud(Npoints,Nlvgrid))
    if (Lcltlidarradar) allocate(x%radar_lidar_tcc(Npoints))
        
    ! RTTOV
    if (Ltbrttov) allocate(x%rttov_tbs(Npoints,Nchan))

    ! Joint MODIS/CloudSat Statistics (T.Michibata, 2018.10.29)
    if (Lwr_occfreq) allocate(x%wr_occfreq_ntotal(Npoints,WR_NREGIME))
    if (Lcfodd)      allocate(x%cfodd_ntotal(Npoints,CFODD_NDBZE,CFODD_NICOD,CFODD_NCLASS))

    ! added by YN
    if (Lscops) then
       allocate(x%frac_ls(Npoints,Nlevels))
       allocate(x%frac_cv(Npoints,Nlevels))
       allocate(x%prec_ls(Npoints,Nlevels))
       allocate(x%prec_cv(Npoints,Nlevels))
    end if

#ifdef OPT_DPLRW
    if (Ldplrw) then
       !allocate(x%dplrw_sumZ(Npoints, Nlvgrid, N_ISCCP, 0:2))
       !allocate(x%dplrw_sumT(Npoints, Nlvtemp, N_ISCCP, 0:2))
       
       allocate(x%dplrw_Z(Npoints, Nlvdplr, Nlvgrid, N_ISCCP, 0:2))
       allocate(x%spwid_Z(Npoints, Nlvspwd, Nlvgrid, N_ISCCP, 0:2))
       allocate(x%Zef94_Z(Npoints, NlvdBZe, Nlvgrid, N_ISCCP, 0:2))

       allocate(x%dplrw_T(Npoints, Nlvdplr, Nlvtemp, N_ISCCP, 0:2))
       allocate(x%spwid_T(Npoints, Nlvspwd, Nlvtemp, N_ISCCP, 0:2))
       allocate(x%Zef94_T(Npoints, NlvdBZe, Nlvtemp, N_ISCCP, 0:2))

       allocate(x%ZefVd_2(Npoints, Nlvdplr, NlvdBZe, N_ISCCP, 0:2))

       allocate(x%gcumw(npoints,nlevels))
    end if
#endif

  end subroutine construct_cosp_outputs
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE destroy_cospIN     
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine destroy_cospIN(y)
  USE mod_cosp, ONLY: cosp_optical_inputs
    type(cosp_optical_inputs),intent(inout) :: y

    if (allocated(y%tau_067))         deallocate(y%tau_067)
    if (allocated(y%emiss_11))        deallocate(y%emiss_11)
    if (allocated(y%frac_out))        deallocate(y%frac_out)
    if (allocated(y%frac_prec))       deallocate(y%frac_prec) ! added by YN
    if (allocated(y%beta_mol))        deallocate(y%beta_mol)
    if (allocated(y%tau_mol))         deallocate(y%tau_mol)
    if (allocated(y%betatot))         deallocate(y%betatot)
    if (allocated(y%betatot_ice))     deallocate(y%betatot_ice)
    if (allocated(y%betatot_liq))     deallocate(y%betatot_liq)
    if (allocated(y%tautot))          deallocate(y%tautot)
    if (allocated(y%tautot_ice))      deallocate(y%tautot_ice)
    if (allocated(y%tautot_liq))      deallocate(y%tautot_liq)
    if (allocated(y%tautot_S_liq))    deallocate(y%tautot_S_liq)
    if (allocated(y%tautot_S_ice))    deallocate(y%tautot_S_ice)
    if (allocated(y%z_vol_cloudsat))  deallocate(y%z_vol_cloudsat)
    if (allocated(y%kr_vol_cloudsat)) deallocate(y%kr_vol_cloudsat)
    if (allocated(y%g_vol_cloudsat))  deallocate(y%g_vol_cloudsat)
    if (allocated(y%asym))            deallocate(y%asym)
    if (allocated(y%ss_alb))          deallocate(y%ss_alb)
    if (allocated(y%fracLiq))         deallocate(y%fracLiq)
    if (allocated(y%fracPrecipIce))   deallocate(y%fracPrecipIce) !<-added
#ifdef OPT_DPLRW
    if (allocated(y%vfall))      deallocate(y%vfall)
    if (allocated(y%vfsqu))      deallocate(y%vfsqu)
    if (allocated(y%zehyd))      deallocate(y%zehyd)
#endif

  end subroutine destroy_cospIN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE destroy_cospstateIN     
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  subroutine destroy_cospstateIN(y)
  USE mod_cosp, ONLY: cosp_column_inputs
    type(cosp_column_inputs),intent(inout) :: y

    if (allocated(y%sunlit))          deallocate(y%sunlit)
    if (allocated(y%skt))             deallocate(y%skt)
    if (allocated(y%land))            deallocate(y%land)
    if (allocated(y%at))              deallocate(y%at)
    if (allocated(y%pfull))           deallocate(y%pfull)
    if (allocated(y%phalf))           deallocate(y%phalf)
    if (allocated(y%qv))              deallocate(y%qv)
    if (allocated(y%o3))              deallocate(y%o3)
    if (allocated(y%hgt_matrix))      deallocate(y%hgt_matrix)
    if (allocated(y%u_sfc))           deallocate(y%u_sfc)
    if (allocated(y%v_sfc))           deallocate(y%v_sfc)
    if (allocated(y%lat))             deallocate(y%lat)
    if (allocated(y%lon))             deallocate(y%lon)
    if (allocated(y%emis_sfc))        deallocate(y%emis_sfc)
    if (allocated(y%cloudIce))        deallocate(y%cloudIce)
    if (allocated(y%cloudLiq))        deallocate(y%cloudLiq)
    if (allocated(y%seaice))          deallocate(y%seaice)
    if (allocated(y%fl_rain))         deallocate(y%fl_rain)
    if (allocated(y%fl_snow))         deallocate(y%fl_snow)
    if (allocated(y%tca))             deallocate(y%tca)
    if (allocated(y%hgt_matrix_half)) deallocate(y%hgt_matrix_half)    
    if (allocated(y%surfelev))        deallocate(y%surfelev) !<-added
#ifdef OPT_DPLRW
    if (allocated(y%gwvel))           deallocate(y%gwvel)
    if (allocated(y%gcumf))           deallocate(y%gcumf)
#endif

  end subroutine destroy_cospstateIN
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE destroy_cosp_outputs
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  subroutine destroy_cosp_outputs(y)
  USE mod_cosp, ONLY: cosp_outputs
     type(cosp_outputs),intent(inout) :: y
#ifdef OPT_DPLRW
     integer :: iD, nD
#endif

     ! Deallocate and nullify
     if (associated(y%calipso_beta_mol))          then
        deallocate(y%calipso_beta_mol)
        nullify(y%calipso_beta_mol)
     endif
     if (associated(y%calipso_temp_tot))          then
        deallocate(y%calipso_temp_tot)
        nullify(y%calipso_temp_tot)     
     endif
     if (associated(y%calipso_betaperp_tot))      then
        deallocate(y%calipso_betaperp_tot)
        nullify(y%calipso_betaperp_tot)     
     endif
     if (associated(y%calipso_beta_tot))          then
        deallocate(y%calipso_beta_tot)    
        nullify(y%calipso_beta_tot)     
     endif
     if (associated(y%calipso_tau_tot))           then
        deallocate(y%calipso_tau_tot) 
        nullify(y%calipso_tau_tot)     
     endif
     if (associated(y%calipso_lidarcldphase))     then
        deallocate(y%calipso_lidarcldphase)
        nullify(y%calipso_lidarcldphase)     
     endif
     if (associated(y%calipso_cldlayerphase))     then
        deallocate(y%calipso_cldlayerphase)
        nullify(y%calipso_cldlayerphase)     
     endif
     if (associated(y%calipso_lidarcldtmp))       then
        deallocate(y%calipso_lidarcldtmp)
        nullify(y%calipso_lidarcldtmp)     
     endif
     if (associated(y%calipso_cldlayer))          then
        deallocate(y%calipso_cldlayer)
        nullify(y%calipso_cldlayer)     
     endif
   !   if (associated(y%calipso_cldtype))          then
   !     deallocate(y%calipso_cldtype)
   !     nullify(y%calipso_cldtype)
   !  endif
   !  if (associated(y%calipso_cldtypetemp))      then
   !     deallocate(y%calipso_cldtypetemp)
   !     nullify(y%calipso_cldtypetemp)
   !  endif
   !  if (associated(y%calipso_cldtypemeanz))     then
   !     deallocate(y%calipso_cldtypemeanz)
   !     nullify(y%calipso_cldtypemeanz)
   !  endif
   !  if (associated(y%calipso_cldtypemeanzse))   then
   !     deallocate(y%calipso_cldtypemeanzse)
   !     nullify(y%calipso_cldtypemeanzse)
   !  endif
   !  if (associated(y%calipso_cldthinemis))      then
   !     deallocate(y%calipso_cldthinemis)
   !     nullify(y%calipso_cldthinemis)
   !  endif
     if (associated(y%calipso_lidarcld))         then
        deallocate(y%calipso_lidarcld)
        nullify(y%calipso_lidarcld)     
     endif
     if (associated(y%calipso_srbval))            then
        deallocate(y%calipso_srbval)
        nullify(y%calipso_srbval)     
     endif
     if (associated(y%calipso_cfad_sr))          then
        deallocate(y%calipso_cfad_sr)
        nullify(y%calipso_cfad_sr)     
     endif
     if (associated(y%parasolPix_refl))           then
        deallocate(y%parasolPix_refl)
        nullify(y%parasolPix_refl)     
     endif
     if (associated(y%parasolGrid_refl))          then
        deallocate(y%parasolGrid_refl) 
        nullify(y%parasolGrid_refl)     
     endif
     if (associated(y%cloudsat_Ze_tot))           then
        deallocate(y%cloudsat_Ze_tot) 
        nullify(y%cloudsat_Ze_tot)  
     endif
     if (associated(y%cloudsat_cfad_ze))          then
        deallocate(y%cloudsat_cfad_ze)
        nullify(y%cloudsat_cfad_ze)     
     endif
     if (associated(y%radar_lidar_tcc))           then
        deallocate(y%radar_lidar_tcc) 
        nullify(y%radar_lidar_tcc)  
     endif
     if (associated(y%lidar_only_freq_cloud))     then
        deallocate(y%lidar_only_freq_cloud)
        nullify(y%lidar_only_freq_cloud)     
     endif
     if (associated(y%isccp_totalcldarea))        then
        deallocate(y%isccp_totalcldarea) 
        nullify(y%isccp_totalcldarea)  
     endif
     if (associated(y%isccp_meantb))              then
        deallocate(y%isccp_meantb) 
        nullify(y%isccp_meantb)     
     endif
     if (associated(y%isccp_meantbclr))           then
        deallocate(y%isccp_meantbclr)
        nullify(y%isccp_meantbclr)  
     endif
     if (associated(y%isccp_meanptop))            then
        deallocate(y%isccp_meanptop)
        nullify(y%isccp_meanptop)     
     endif
     if (associated(y%isccp_meantaucld))          then
        deallocate(y%isccp_meantaucld) 
        nullify(y%isccp_meantaucld)       
     endif
     if (associated(y%isccp_meanalbedocld))       then
        deallocate(y%isccp_meanalbedocld)
        nullify(y%isccp_meanalbedocld)     
     endif
     if (associated(y%isccp_boxtau))              then
        deallocate(y%isccp_boxtau)
        nullify(y%isccp_boxtau)       
     endif
     if (associated(y%isccp_boxptop))             then
        deallocate(y%isccp_boxptop)
        nullify(y%isccp_boxptop)     
     endif
     if (associated(y%isccp_fq))                  then
        deallocate(y%isccp_fq)
        nullify(y%isccp_fq)       
     endif
     if (associated(y%misr_fq))                   then
        deallocate(y%misr_fq) 
        nullify(y%misr_fq)     
     endif
     if (associated(y%misr_dist_model_layertops)) then
        deallocate(y%misr_dist_model_layertops)
        nullify(y%misr_dist_model_layertops)       
     endif
     if (associated(y%misr_meanztop))             then
        deallocate(y%misr_meanztop)
        nullify(y%misr_meanztop)     
     endif
     if (associated(y%misr_cldarea))              then
        deallocate(y%misr_cldarea)
        nullify(y%misr_cldarea)      
     endif
     if (associated(y%rttov_tbs))                 then
        deallocate(y%rttov_tbs)
        nullify(y%rttov_tbs)     
     endif
     if (associated(y%modis_Cloud_Fraction_Total_Mean))                      then
        deallocate(y%modis_Cloud_Fraction_Total_Mean)       
        nullify(y%modis_Cloud_Fraction_Total_Mean)       
     endif
     if (associated(y%modis_Cloud_Fraction_Ice_Mean))                        then
        deallocate(y%modis_Cloud_Fraction_Ice_Mean)     
        nullify(y%modis_Cloud_Fraction_Ice_Mean)     
     endif
     if (associated(y%modis_Cloud_Fraction_Water_Mean))                      then
        deallocate(y%modis_Cloud_Fraction_Water_Mean)           
        nullify(y%modis_Cloud_Fraction_Water_Mean)           
     endif
     if (associated(y%modis_Cloud_Fraction_High_Mean))                       then
        deallocate(y%modis_Cloud_Fraction_High_Mean)     
        nullify(y%modis_Cloud_Fraction_High_Mean)     
     endif
     if (associated(y%modis_Cloud_Fraction_Mid_Mean))                        then
        deallocate(y%modis_Cloud_Fraction_Mid_Mean)       
        nullify(y%modis_Cloud_Fraction_Mid_Mean)       
     endif
     if (associated(y%modis_Cloud_Fraction_Low_Mean))                        then
        deallocate(y%modis_Cloud_Fraction_Low_Mean)     
        nullify(y%modis_Cloud_Fraction_Low_Mean)     
     endif
     if (associated(y%modis_Optical_Thickness_Total_Mean))                   then
        deallocate(y%modis_Optical_Thickness_Total_Mean)  
        nullify(y%modis_Optical_Thickness_Total_Mean)  
     endif
     if (associated(y%modis_Optical_Thickness_Water_Mean))                   then
        deallocate(y%modis_Optical_Thickness_Water_Mean)     
        nullify(y%modis_Optical_Thickness_Water_Mean)     
     endif
     if (associated(y%modis_Optical_Thickness_Ice_Mean))                     then
        deallocate(y%modis_Optical_Thickness_Ice_Mean)       
        nullify(y%modis_Optical_Thickness_Ice_Mean)       
     endif
     if (associated(y%modis_Optical_Thickness_Total_LogMean))                then
        deallocate(y%modis_Optical_Thickness_Total_LogMean)    
        nullify(y%modis_Optical_Thickness_Total_LogMean)    
     endif
     if (associated(y%modis_Optical_Thickness_Water_LogMean))                then
        deallocate(y%modis_Optical_Thickness_Water_LogMean)     
        nullify(y%modis_Optical_Thickness_Water_LogMean)     
     endif
     if (associated(y%modis_Optical_Thickness_Ice_LogMean))                  then
        deallocate(y%modis_Optical_Thickness_Ice_LogMean)     
        nullify(y%modis_Optical_Thickness_Ice_LogMean)     
     endif
     if (associated(y%modis_Cloud_Particle_Size_Water_Mean))                 then
        deallocate(y%modis_Cloud_Particle_Size_Water_Mean)       
        nullify(y%modis_Cloud_Particle_Size_Water_Mean)       
     endif
     if (associated(y%modis_Cloud_Particle_Size_Ice_Mean))                   then
        deallocate(y%modis_Cloud_Particle_Size_Ice_Mean)     
        nullify(y%modis_Cloud_Particle_Size_Ice_Mean)     
     endif
     if (associated(y%modis_Cloud_Top_Pressure_Total_Mean))                  then
        deallocate(y%modis_Cloud_Top_Pressure_Total_Mean)           
        nullify(y%modis_Cloud_Top_Pressure_Total_Mean)           
     endif
     if (associated(y%modis_Liquid_Water_Path_Mean))                         then
        deallocate(y%modis_Liquid_Water_Path_Mean)     
        nullify(y%modis_Liquid_Water_Path_Mean)     
     endif
     if (associated(y%modis_Ice_Water_Path_Mean))                            then
        deallocate(y%modis_Ice_Water_Path_Mean)       
        nullify(y%modis_Ice_Water_Path_Mean)       
     endif
     if (associated(y%modis_Optical_Thickness_vs_Cloud_Top_Pressure))        then
        deallocate(y%modis_Optical_Thickness_vs_Cloud_Top_Pressure)     
        nullify(y%modis_Optical_Thickness_vs_Cloud_Top_Pressure)     
     endif
     if (associated(y%modis_Optical_thickness_vs_ReffLIQ))                   then
        deallocate(y%modis_Optical_thickness_vs_ReffLIQ)
        nullify(y%modis_Optical_thickness_vs_ReffLIQ)
     endif
     if (associated(y%modis_Optical_thickness_vs_ReffICE))                   then
        deallocate(y%modis_Optical_thickness_vs_ReffICE)
        nullify(y%modis_Optical_thickness_vs_ReffICE)
     endif
     !! IDiD WR (T.Michibata, 2018.10.29)
     if (associated(y%cfodd_ntotal)) then
        deallocate(y%cfodd_ntotal)
        nullify(y%cfodd_ntotal)
     endif
     if (associated(y%wr_occfreq_ntotal)) then
        deallocate(y%wr_occfreq_ntotal)
        nullify(y%wr_occfreq_ntotal)
     endif
     !! (Y.Imura, 2021.06.16) !<-added
     if (associated(y%cloudsat_precip_cover)) then
        deallocate(y%cloudsat_precip_cover)
        nullify(y%cloudsat_precip_cover)
     endif
     if (associated(y%cloudsat_pia)) then
        deallocate(y%cloudsat_pia)
        nullify(y%cloudsat_pia)
     endif
     ! added by YN
     if (associated(y%frac_ls)) then
        deallocate(y%frac_ls) ; nullify(y%frac_ls)
     end if
     if (associated(y%frac_cv)) then
        deallocate(y%frac_cv) ; nullify(y%frac_cv)
     end if
     if (associated(y%prec_ls)) then
        deallocate(y%prec_ls) ; nullify(y%prec_ls)
     end if
     if (associated(y%prec_cv)) then
        deallocate(y%prec_cv) ; nullify(y%prec_cv)
     end if
#ifdef OPT_DPLRW
     ! if (associated(y%dplrw_sumZ)) then
     !    deallocate(y%dplrw_sumZ) ; nullify(y%dplrw_sumZ)
     ! end if
     ! if (associated(y%dplrw_sumT)) then
     !    deallocate(y%dplrw_sumT) ; nullify(y%dplrw_sumT)
     ! end if

     if (associated(y%dplrw_Z)) then
        deallocate(y%dplrw_Z) ; nullify(y%dplrw_Z)
     end if
     if (associated(y%spwid_Z)) then
        deallocate(y%spwid_Z) ; nullify(y%spwid_Z)
     end if
     if (associated(y%Zef94_Z)) then
        deallocate(y%Zef94_Z) ; nullify(y%Zef94_Z)
     end if

     if (associated(y%dplrw_T)) then
        deallocate(y%dplrw_T) ; nullify(y%dplrw_T)
     end if
     if (associated(y%spwid_T)) then
        deallocate(y%spwid_T) ; nullify(y%spwid_T)
     end if
     if (associated(y%Zef94_T)) then
        deallocate(y%Zef94_T) ; nullify(y%Zef94_T)
     end if

     if (associated(y%ZefVd_2)) then
        deallocate(y%ZefVd_2) ; nullify(y%ZefVd_2)
     end if

     if (associated(y%gcumw)) then
        deallocate(y%gcumw) ; nullify(y%gcumw)
     end if
#endif

   end subroutine destroy_cosp_outputs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------- SUBROUTINE UNDEF_M2N ----------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE UNDEF_M2N(x, nx, ny, nz)
!  USE MOD_COSP_CONSTANTS
  IMPLICIT NONE
  !
  integer,intent(in) :: nx, ny, nz
  real(8),dimension(nx,ny,nz),intent(inout) :: x
  real(8), parameter :: MIROC_UNDEF = -999.D0
  real(8), parameter :: R_UNDEF = -1.0E30
  !
  where( x /= MIROC_UNDEF )
     x = x
  elsewhere
     x = R_UNDEF
  end where
  !
END SUBROUTINE UNDEF_M2N

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------- SUBROUTINE UNDEF_N2M ----------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE UNDEF_N2M(x, nx, ny, nz)
!  USE MOD_COSP_CONSTANTS
  IMPLICIT NONE
  !
  integer,intent(in) :: nx, ny, nz
  real(8),dimension(nx,ny,nz),intent(inout) :: x
  real(8), parameter :: MIROC_UNDEF = -999.D0
  real(8), parameter :: R_UNDEF = -1.0E30
  !
!  where( x /= R_UNDEF )
  where( x > R_UNDEF * 0.1D0 )
     x = x
  elsewhere
     x = MIROC_UNDEF
  end where
  !
END SUBROUTINE UNDEF_N2M

