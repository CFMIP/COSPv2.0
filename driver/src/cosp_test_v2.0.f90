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
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
program cosp_test_v2
  use cosp_kinds,          only: wp
  use mod_cosp_config,     only: RTTOV_MAX_CHANNELS,N_HYDRO,numMODISTauBins,modis_histTau,&
                                 modis_histTauEdges,modis_histTauCenters,ntau,ntauV1p4,   &
                                 tau_binBounds,tau_binEdges,tau_binCenters,R_UNDEF,       &
                                 tau_binBoundsV1p4,tau_binEdgesV1p4,tau_binCentersV1p4
  use cosp_phys_constants, only: amw,amd,amO3,amCO2,amCH4,amN2O,amCO
  use mod_cosp_io,         only: nc_read_input_file,nc_cmor_init
  USE mod_quickbeam_optics,only: size_distribution,hydro_class_init,quickbeam_optics
  use quickbeam,           only: radar_cfg
  use mod_cosp,            only: cosp_init,cosp_optical_inputs,cosp_column_inputs,        &
                                 construct_cosp_outputs,cosp_outputs,construct_cospIN,    &
                                 construct_cospstateIN,destroy_cosp_outputs,              &
                                 destroy_cospIN,destroy_cospstateIN,cosp_simulator
  USE mod_rng,             ONLY: rng_state, init_rng
  USE mod_scops,           ONLY: scops
  USE mod_prec_scops,      ONLY: prec_scops
  USE MOD_COSP_UTILS,      ONLY: cosp_precip_mxratio
  use cosp_optics,         ONLY: cosp_simulator_optics,lidar_optics,modis_optics,         &
                                 modis_optics_partition,num_trial_res,clara_optics
  ! OUTPUT ONLY
  use mod_cosp_io,         only: var1d,var2d,var3d,construct_cospOutList,nc_cmor_init,    &
                                 nc_cmor_associate_1d,nc_cmor_write_1d,nc_cmor_close,     &
                                 nc_cmor_associate_2d,nc_cmor_write_2d
  
  implicit none

  ! Input/Output driver file control
  character(len=64),parameter :: &
       cosp_input_namelist  = 'cosp_input_nl.v2.0.txt', &
       cosp_output_namelist = 'cosp_output_nl_v2.0.txt'

  ! Test data
  integer :: &
       Nlon,Nlat,geomode
  real(wp) :: &
       emsfc_lw
  real(wp),dimension(:),allocatable,target:: &
       lon,       & ! Longitude (deg)
       lat,       & ! Latitude (deg)
       skt,       & ! Skin temperature (K)
       landmask,  & ! Land/sea mask (0/1)
       u_wind,    & ! U-component of wind (m/s)
       v_wind,    & ! V-component of wind (m/s)
       sunlit       ! Sunlit flag
  real(wp),dimension(:,:),allocatable,target :: &
       p,         & ! Model pressure levels (pa)
       ph,        & ! Moddel pressure @ half levels (pa)
       zlev,      & ! Model level height (m)
       zlev_half, & ! Model level height @ half-levels (m)
       T,         & ! Temperature (K)
       sh,        & ! Specific humidity (kg/kg)
       rh,        & ! Relative humidity (1)
       tca,       & ! Total cloud fraction (1)
       cca,       & ! Convective cloud fraction (1) 
       mr_lsliq,  & ! Mass mixing ratio for stratiform cloud liquid (kg/kg)
       mr_lsice,  & ! Mass mixing ratio for stratiform cloud ice (kg/kg)
       mr_ccliq,  & ! Mass mixing ratio for convective cloud liquid (kg/kg)
       mr_ccice,  & ! Mass mixing ratio for convective cloud ice (kg/kg)
       mr_ozone,  & ! Mass mixing ratio for ozone (kg/kg)
       fl_lsrain, & ! Precipitation flux (rain) for stratiform cloud (kg/m^2/s)
       fl_lssnow, & ! Precipitation flux (snow) for stratiform cloud (kg/m^2/s)
       fl_lsgrpl, & ! Precipitation flux (groupel) for stratiform cloud (kg/m^2/s)
       fl_ccrain, & ! Precipitation flux (rain) for convective cloud (kg/m^2/s)
       fl_ccsnow, & ! Precipitation flux (snow) for convective cloud (kg/m^2/s)
       dtau_s,    & ! 0.67micron optical depth (stratiform cloud) (1)
       dtau_c,    & ! 0.67micron optical depth (convective cloud) (1)
       dem_s,     & ! 11micron emissivity (stratiform cloud) 
       dem_c        ! 11microm emissivity (convective cloud)
  real(wp),dimension(:,:,:),allocatable,target :: &
       frac_out,  & ! Subcolumn cloud cover (0/1)
       frac_prec, & ! Subcolumn precipitation fraction
       Reff         ! Subcolumn effective radius

  ! Input namelist fields
  integer ::                      & !
       Npoints,                   & ! Number of gridpoints
       Ncolumns,                  & ! Number of subcolumns
       Nlevels,                   & ! Number of model vertical levels
       Npoints_it,                & ! Number of gridpoints to be processed in one 
                                    ! iteration
       Nlvgrid,                   & ! Number of vertical levels for statistical outputs 
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
       rttov_Nchannels,           & ! RTTOV: Number of channels to be computed
       CLARA_Tb_subvis,           & ! Method for Tb (sub-visible clouds)
       CLARA_Tb_semitrans,        & ! Method for Tb (semi-transparent clouds)
       CLARA_Tb_opaque,           & ! Method for Tb (opaque clouds)
       claraRTTOV_platform,       & ! Satellite platform
       claraRTTOV_satellite,      & ! Satellite ID
       claraRTTOV_sensor,         & ! Sensor ID
       claraRTTOV_nchan             ! Number of channels
  integer,dimension(3) :: &
       claraRTTOV_channels          ! Channel numbers
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
       use_precipitation_fluxes,  & ! True if precipitation fluxes are input to the 
                                    ! algorithm 
       CLARA_RTTOVclr,            & ! True => Use RTTOV for cloudy free scenes 
       CLARA_retSize,             & ! True => use TOA reflectance minimization to determine cloud particle size.
       claraRTTOV_addrefrac,          & !
       claraRTTOV_use_q2m,            & !
       claraRTTOV_clw_data,           & !
       claraRTTOV_addsolar,           & !
       claraRTTOV_addclouds,          & !
       claraRTTOV_addaerosol,         & !
       claraRTTOV_use_cld_opts_param, & !
       claraRTTOV_ozone_data,         & !
       claraRTTOV_co2,                & !
       claraRTTOV_n2o,                & !
       claraRTTOV_ch4,                & !
       claraRTTOV_co,                 & !
       claraRTTOV_addinterp,          & !
       claraRTTOV_calcemis,           & ! 
       claraRTTOV_calcrefl       
  integer,dimension(RTTOV_MAX_CHANNELS) :: &
       rttov_Channels                 ! RTTOV: Channel numbers
  real(wp),dimension(RTTOV_MAX_CHANNELS) :: &
       rttov_Surfem                   ! RTTOV: Surface emissivity
  character(len=64) :: &
       cloudsat_micro_scheme,       & ! Microphysical scheme used in cloudsat radar simulator
       claraRTTOV_coefdir             ! Location of RTTOV coefficient files
  character(len=64) :: &
       finput                       ! Input NetCDF file
  character(len=512) :: &
       dinput                       ! Directory where the input files are located
  character(len=600) :: &
       fileIN                       ! dinput+finput
  namelist/COSP_INPUT/overlap,isccp_topheight,isccp_topheight_direction,npoints,         &
                      npoints_it,ncolumns,nlevels,use_vgrid,Nlvgrid,csat_vgrid,dinput,   &
                      finput,cloudsat_radar_freq,surface_radar,cloudsat_use_gas_abs,     &
                      cloudsat_do_ray,cloudsat_k2,cloudsat_micro_scheme,lidar_ice_type,  &
                      use_precipitation_fluxes,rttov_platform,rttov_satellite,           &
                      rttov_Instrument,rttov_Nchannels,rttov_Channels,rttov_Surfem,      &
                      rttov_ZenAng,co2,ch4,n2o,co,                                       &
                      CLARA_Tb_subvis,CLARA_Tb_semitrans,CLARA_Tb_opaque,CLARA_RTTOVclr, &
                      claraRTTOV_coefdir,claraRTTOV_platform,claraRTTOV_satellite,       &
                      claraRTTOV_sensor,claraRTTOV_nchan,claraRTTOV_channels,            &
                      claraRTTOV_addrefrac,claraRTTOV_use_q2m,claraRTTOV_clw_data,       &
                      claraRTTOV_addsolar,claraRTTOV_addclouds,claraRTTOV_addaerosol,    &
                      claraRTTOV_use_cld_opts_param,claraRTTOV_ozone_data,claraRTTOV_co2,&
                      claraRTTOV_n2o,claraRTTOV_ch4,claraRTTOV_co,claraRTTOV_addinterp,  &
                      claraRTTOV_calcemis,claraRTTOV_calcrefl,CLARA_retSize
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
             Lclmodis,Ltbrttov,Ltauclara,Lpctclara,Ltctclara,Lhctclara,Lreffclwclara,    &
             Lreffcliclara,Lcltclara,Llwpclara,Liwpclara,Lclclara
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
                       Lreffclimodis,Lpctmodis,Llwpmodis,Liwpmodis,Lclmodis,Ltbrttov,    &
                       Ltauclara,Lpctclara,Ltctclara,Lhctclara,Lreffclwclara,            &
                       Lreffcliclara,Lcltclara,Llwpclara,Liwpclara,Lclclara

  ! Local variables
  character(len=32),parameter :: &
       cospvID = 'COSP v1.5'        ! COSP version I
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

   ! Microphysical settings for the precipitation flux to mixing ratio conversion
  real(wp),parameter,dimension(N_HYDRO) :: &
                 ! LSL   LSI      LSR       LSS   CVL  CVI      CVR       CVS       LSG
       N_ax    = (/-1., -1.,     8.e6,     3.e6, -1., -1.,     8.e6,     3.e6,     4.e6/),&
       N_bx    = (/-1., -1.,      0.0,      0.0, -1., -1.,      0.0,      0.0,      0.0/),&
       alpha_x = (/-1., -1.,      0.0,      0.0, -1., -1.,      0.0,      0.0,      0.0/),&
       c_x     = (/-1., -1.,    842.0,     4.84, -1., -1.,    842.0,     4.84,     94.5/),&
       d_x     = (/-1., -1.,      0.8,     0.25, -1., -1.,      0.8,     0.25,      0.5/),&
       g_x     = (/-1., -1.,      0.5,      0.5, -1., -1.,      0.5,      0.5,      0.5/),&
       a_x     = (/-1., -1.,    524.0,    52.36, -1., -1.,    524.0,    52.36,   209.44/),&
       b_x     = (/-1., -1.,      3.0,      3.0, -1., -1.,      3.0,      3.0,      3.0/),&
       gamma_1 = (/-1., -1., 17.83725, 8.284701, -1., -1., 17.83725, 8.284701, 11.63230/),&
       gamma_2 = (/-1., -1.,      6.0,      6.0, -1., -1.,      6.0,      6.0,      6.0/),&
       gamma_3 = (/-1., -1.,      2.0,      2.0, -1., -1.,      2.0,      2.0,      2.0/),&
       gamma_4 = (/-1., -1.,      6.0,      6.0, -1., -1.,      6.0,      6.0,      6.0/)       

  ! Fields used solely for output
  integer,parameter :: &
       n_out_list = 73,           & ! Number of possible output variables
       N3D        = 9,            & ! Number of 3D output variables
       N2D        = 14,           & ! Number of 2D output variables
       N1D        = 49              ! Number of 1D output variables  
  character(len=32),dimension(n_out_list) :: out_list  ! List of output variable names
  integer :: lon_axid,time_axid,height_axid,height_mlev_axid,grid_id,lonvar_id,       &
             latvar_id,column_axid,sza_axid,temp_axid,channel_axid,dbze_axid,sratio_axid,&
             MISR_CTH_axid,lat_axid,tau_axid,pressure2_axid 
  type(var1d) :: v1d(N1D+1) ! Structures needed by output routines for 1D variables
  type(var2d) :: v2d(N2D)   ! Structures needed by output routines for 2D variables
  type(var3d) :: v3d(N3D)   ! Structures needed by output routines for 3D variables
  double precision :: time,time_bnds(2),time_step,half_time_step
  real(wp),dimension(:),allocatable :: mgrid_z,mgrid_zu,mgrid_zl
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  call cpu_time(driver_time(1))
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Read in namelists
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Input namelist (cosp setup)
  open(10,file=cosp_input_namelist,status='unknown')
  read(10,nml=cosp_input)
  close(10)

  ! Output namelist (logical flags to turn on/off outputs)
  open(10,file=cosp_output_namelist,status='unknown')
  read(10,nml=cosp_output)
  close(10)

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Read in sample input data.
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  allocate(lon(Npoints),lat(Npoints),p(Npoints,Nlevels),ph(Npoints,Nlevels),             &
           zlev(Npoints,Nlevels),zlev_half(Npoints,Nlevels),T(Npoints,Nlevels),          &
           sh(Npoints,Nlevels),rh(Npoints,Nlevels),tca(Npoints,Nlevels),                 &
           cca(Npoints,Nlevels),mr_lsliq(Npoints,Nlevels),mr_lsice(Npoints,Nlevels),     &
           mr_ccliq(Npoints,Nlevels),mr_ccice(Npoints,Nlevels),                          &
           fl_lsrain(Npoints,Nlevels),fl_lssnow(Npoints,Nlevels),                        &
           fl_lsgrpl(Npoints,Nlevels),fl_ccrain(Npoints,Nlevels),                        &
           fl_ccsnow(Npoints,Nlevels),Reff(Npoints,Nlevels,N_HYDRO),                     &
           dtau_s(Npoints,Nlevels),dtau_c(Npoints,Nlevels),dem_s(Npoints,Nlevels),       &
           dem_c(Npoints,Nlevels),skt(Npoints),landmask(Npoints),                        &
           mr_ozone(Npoints,Nlevels),u_wind(Npoints),v_wind(Npoints),sunlit(Npoints),    &
           frac_out(Npoints,Ncolumns,Nlevels))

  fileIN = trim(dinput)//trim(finput)
  call nc_read_input_file(fileIN,Npoints,Nlevels,N_HYDRO,lon,lat,p,ph,zlev,zlev_half,    &
                          T,sh,rh,tca,cca,mr_lsliq,mr_lsice,mr_ccliq,mr_ccice,fl_lsrain, &
                          fl_lssnow,fl_lsgrpl,fl_ccrain,fl_ccsnow,Reff,dtau_s,dtau_c,    &
                          dem_s,dem_c,skt,landmask,mr_ozone,u_wind,v_wind,sunlit,        &
                          emsfc_lw,geomode,Nlon,Nlat)
  call cpu_time(driver_time(2))

  ! Which simulators need to be run? Look at which outputs are requested.
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
  if (LcfadDbze94 .or. Ldbze94 .or. Lcltlidarradar) Lcloudsat = .true.
  if (Lparasolrefl) Lparasol = .true.
  if (Ltbrttov) Lrttov = .true.
  
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

  ! Initialize quickbeam_optics, also if two-moment radar microphysics scheme is wanted...
  if (cloudsat_micro_scheme == 'MMF_v3.5_two_moment')  then
     ldouble = .true. 
     lsingle = .false.
  endif  

  ! Initialize the distributional parameters for hydrometeors in radar simulator
  call hydro_class_init(R_UNDEF,lsingle,ldouble,sd)
  
  ! Initialize COSP simulator
  call COSP_INIT(Lisccp,Lmodis,Lmisr,Lcloudsat,Lcalipso,Lparasol,Lrttov,                 &
       Npoints,Nlevels,cloudsat_radar_freq,cloudsat_k2,cloudsat_use_gas_abs,  &
                 cloudsat_do_ray,isccp_topheight,isccp_topheight_direction,surface_radar,&
                 rcfg_cloudsat,rttov_Nchannels,rttov_Channels,rttov_platform,           &
                 rttov_satellite,rttov_instrument,use_vgrid,csat_vgrid,Nlvgrid,         &
                 cloudsat_micro_scheme,                                                 &
                 CLARA_Tb_subvis,CLARA_Tb_semitrans,CLARA_Tb_opaque,CLARA_RTTOVclr,     &
                 claraRTTOV_coefdir,claraRTTOV_platform,claraRTTOV_satellite,           &
                 claraRTTOV_sensor,claraRTTOV_nchan,claraRTTOV_channels,                &
                 claraRTTOV_addrefrac,claraRTTOV_use_q2m,claraRTTOV_clw_data,           &
                 claraRTTOV_addsolar,claraRTTOV_addclouds,claraRTTOV_addaerosol,        &
                 claraRTTOV_use_cld_opts_param,claraRTTOV_ozone_data,claraRTTOV_co2,    &
                 claraRTTOV_n2o,claraRTTOV_ch4,claraRTTOV_co,claraRTTOV_addinterp,      &
                 claraRTTOV_calcemis,claraRTTOV_calcrefl,CLARA_retSize)
  call cpu_time(driver_time(3))

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Construct output derived type.
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  call construct_cosp_outputs(Lpctisccp,Lclisccp,Lboxptopisccp, Lboxtauisccp,Ltauisccp,  &
                              Lcltisccp,Lmeantbisccp,Lmeantbclrisccp,Lalbisccp,LclMISR,  &
                              Lcltmodis,Lclwmodis,Lclimodis,Lclhmodis,Lclmmodis,         &
                              Lcllmodis,Ltautmodis,Ltauwmodis,Ltauimodis,Ltautlogmodis,  &
                              Ltauwlogmodis,Ltauilogmodis,Lreffclwmodis,Lreffclimodis,   &
                              Lpctmodis,Llwpmodis,Liwpmodis,Lclmodis,Latb532,            &
                              LlidarBetaMol532,LcfadLidarsr532,Lclcalipso2,Lclcalipso,   &
                              Lclhcalipso,Lcllcalipso,Lclmcalipso,Lcltcalipso,           &
                              Lcltlidarradar,Lclcalipsoliq,Lclcalipsoice,Lclcalipsoun,   &
                              Lclcalipsotmp,Lclcalipsotmpliq,Lclcalipsotmpice,           &
                              Lclcalipsotmpun,Lcltcalipsoliq,Lcltcalipsoice,             &
                              Lcltcalipsoun,Lclhcalipsoliq,Lclhcalipsoice,Lclhcalipsoun, &
                              Lclmcalipsoliq,Lclmcalipsoice,Lclmcalipsoun,Lcllcalipsoliq,&
                              Lcllcalipsoice,Lcllcalipsoun,LcfadDbze94,Ldbze94,          &
                              Lparasolrefl,Ltbrttov,Ltauclara,Lpctclara,Ltctclara,       &
                              Lhctclara,Lreffclwclara,Lreffcliclara,Lcltclara,Llwpclara, &
                              Liwpclara,Lclclara,Npoints,Ncolumns,Nlevels,Nlvgrid,       &
                              rttov_Nchannels,cospOUT)

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Break COSP up into pieces and loop over each COSP 'chunk'.
  ! nChunks = # Points to Process (nPoints) / # Points per COSP iteration (nPoints_it)
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  nChunks = nPoints/nPoints_it+1
  if (nPoints .eq. nPoints_it) nChunks = 1
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
        call construct_cospIN(Nptsperit,nColumns,nLevels,cospIN)
        call construct_cospstateIN(Nptsperit,nLevels,rttov_nChannels,cospstateIN)
     endif
     if (iChunk .eq. nChunks) then
        call destroy_cospIN(cospIN)
        call destroy_cospstateIN(cospstateIN)
        call construct_cospIN(Nptsperit,nColumns,nLevels,cospIN)
        call construct_cospstateIN(Nptsperit,nLevels,rttov_nChannels,cospstateIN)    
     endif
     call cpu_time(driver_time(4))

     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! Populate input types with model fields.
     ! Here the 3D sample model fields (temperature,pressure,etc...) are ordered from the
     ! surface-2-TOA, whereas COSP expects all fields to be ordered from TOA-2-SFC. So the
     ! vertical fields are flipped prior to storing to COSP input type.
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     cospstateIN%hgt_matrix           = zlev(start_idx:end_idx,Nlevels:1:-1)
     cospstateIN%sunlit               = sunlit(start_idx:end_idx)
     cospstateIN%skt                  = skt(start_idx:end_idx)
     cospstateIN%land                 = landmask(start_idx:end_idx)
     cospstateIN%qv                   = sh(start_idx:end_idx,Nlevels:1:-1) 
     cospstateIN%at                   = T(start_idx:end_idx,Nlevels:1:-1) 
     cospstateIN%pfull                = p(start_idx:end_idx,Nlevels:1:-1) 
     cospstateIN%o3                   = mr_ozone(start_idx:end_idx,Nlevels:1:-1)*(amd/amO3)*1e6
     cospstateIN%u_sfc                = u_wind(start_idx:end_idx)
     cospstateIN%v_sfc                = v_wind(start_idx:end_idx)
     cospstateIN%emis_sfc             = rttov_surfem
     cospstateIN%zenang               = rttov_zenang
     cospstateIN%lat                  = lat(start_idx:end_idx)
     cospstateIN%lon                  = lon(start_idx:end_idx)
     cospstateIN%month                = 2 ! This is needed by RTTOV only for the surface emissivity calculation.
     cospstateIN%co2                  = co2*(amd/amCO2)*1e6
     cospstateIN%ch4                  = ch4*(amd/amCH4)*1e6  
     cospstateIN%n2o                  = n2o*(amd/amN2O)*1e6
     cospstateIN%co                   = co*(amd/amCO)*1e6
     cospstateIN%phalf(:,1)           = 0._wp
     cospstateIN%phalf(:,2:Nlevels+1) = ph(start_idx:end_idx,Nlevels:1:-1)      
     cospstateIN%hgt_matrix_half(:,1:Nlevels) = zlev_half(start_idx:end_idx,Nlevels:1:-1)
     cospstateIN%hgt_matrix_half(:,Nlevels+1) = 0._wp  
     cospIN%tautot_S_liq              = 0._wp
     cospIN%tautot_S_ice              = 0._wp
     cospIN%emsfc_lw                  = emsfc_lw
     cospIN%rcfg_cloudsat             = rcfg_cloudsat
     
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! Compute subcolumns and optical inputs to COSP in the same manner as in v1.4, call
     ! subsample_and_optics.
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     call subsample_and_optics(nPtsPerIt,nLevels,nColumns,N_HYDRO,overlap,                     &
          use_precipitation_fluxes,lidar_ice_type,sd,                                          &
          tca(start_idx:end_idx,Nlevels:1:-1),cca(start_idx:end_idx,Nlevels:1:-1),             &
          fl_lsrain(start_idx:end_idx,Nlevels:1:-1),fl_lssnow(start_idx:end_idx,Nlevels:1:-1), &
          fl_lsgrpl(start_idx:end_idx,Nlevels:1:-1),fl_ccrain(start_idx:end_idx,Nlevels:1:-1), &
          fl_ccsnow(start_idx:end_idx,Nlevels:1:-1),mr_lsliq(start_idx:end_idx,Nlevels:1:-1),  &
          mr_lsice(start_idx:end_idx,Nlevels:1:-1),mr_ccliq(start_idx:end_idx,Nlevels:1:-1),   &
          mr_ccice(start_idx:end_idx,Nlevels:1:-1),Reff(start_idx:end_idx,Nlevels:1:-1,:),     &
          dtau_c(start_idx:end_idx,nLevels:1:-1),dtau_s(start_idx:end_idx,nLevels:1:-1),       &
          dem_c(start_idx:end_idx,nLevels:1:-1),dem_s(start_idx:end_idx,nLevels:1:-1),         &
          cospstateIN,cospIN)

     call cpu_time(driver_time(6))
     
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! Call COSP
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     cosp_status = COSP_SIMULATOR(cospIN, cospstateIN, cospOUT,start_idx,end_idx,.false.)
     
     call cpu_time(driver_time(7))
  enddo
  print*,'Time to read in data:     ',driver_time(2)-driver_time(1)
  print*,'Time to initialize:       ',driver_time(3)-driver_time(2)
  print*,'Time to construct types:  ',driver_time(4)-driver_time(3)
  print*,'Time to compute optics:   ',driver_time(6)-driver_time(4)
  print*,'Time to run COSP:         ',driver_time(7)-driver_time(6)
  print*,'Total time:               ',driver_time(7)-driver_time(1)

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Output
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Create list of output varibles.
  call construct_cospOutList(Lpctisccp,Lclisccp,Lboxptopisccp, Lboxtauisccp,Ltauisccp,   &
                             Lcltisccp,Lmeantbisccp,Lmeantbclrisccp,Lalbisccp,LclMISR,   &
                             Lcltmodis,Lclwmodis,Lclimodis,Lclhmodis,Lclmmodis,          &
                             Lcllmodis,Ltautmodis,Ltauwmodis,Ltauimodis,Ltautlogmodis,   &
                             Ltauwlogmodis,Ltauilogmodis,Lreffclwmodis,Lreffclimodis,    &
                             Lpctmodis,Llwpmodis,Liwpmodis,Lclmodis,Latb532,             &
                             LlidarBetaMol532,LcfadLidarsr532,Lclcalipso2,Lclcalipso,    &
                             Lclhcalipso,Lcllcalipso,Lclmcalipso,Lcltcalipso,            &
                             Lcltlidarradar,Lclcalipsoliq,Lclcalipsoice,Lclcalipsoun,    &
                             Lclcalipsotmp,Lclcalipsotmpliq,Lclcalipsotmpice,            &
                             Lclcalipsotmpun,Lcltcalipsoliq,Lcltcalipsoice,              &
                             Lcltcalipsoun,Lclhcalipsoliq,Lclhcalipsoice,Lclhcalipsoun,  &
                             Lclmcalipsoliq,Lclmcalipsoice,Lclmcalipsoun,Lcllcalipsoliq, &
                             Lcllcalipsoice,Lcllcalipsoun,LcfadDbze94,Ldbze94,           &
                             Lparasolrefl,Ltbrttov,Ltauclara,Lpctclara,Ltctclara,        &
                             Lhctclara,Lreffclwclara,Lreffcliclara,Lcltclara,Llwpclara,  &
                             Liwpclara,Lclclara,N_OUT_LIST,out_list)  

  ! Time information for cmor output.
  time           = 8*1._wp/8._wp
  time_step      = 3._wp/24._wp
  half_time_step = 0.5_wp*time_step
  time_bnds      = (/time-half_time_step,time+half_time_step/)

  ! Model grid info for cmor output
  allocate(mgrid_z(Nlevels),mgrid_zl(Nlevels),mgrid_zu(Nlevels))
  mgrid_z             = cospstateIN%hgt_matrix(1,:)
  mgrid_zl            = cospstateIN%hgt_matrix_half(1,:)
  mgrid_zu(2:Nlevels) = cospstateIN%hgt_matrix_half(1,1:Nlevels-1)
  mgrid_zu(1)         = cospstateIN%hgt_matrix(1,1)+(cospstateIN%hgt_matrix(1,1)-mgrid_zl(1))
  
  if (geomode .eq. 1) then
     call nc_cmor_init('../cmor/cosp_cmor_nl_1D.txt','replace',nPoints,nColumns,nLevels, &
                       rttov_nChannels,nLvgrid,lon,lat,mgrid_zl,mgrid_zu,mgrid_z,cospOUT,&
                       geomode,Nlon,Nlat,N1D+1,N2D,N3D,N_OUT_LIST,out_list,lon_axid,     &
                       lat_axid,time_axid,height_axid,height_mlev_axid,grid_id,lonvar_id,&
                       latvar_id,column_axid,sza_axid,temp_axid,channel_axid,dbze_axid,  &
                       sratio_axid,MISR_CTH_axid,tau_axid,pressure2_axid,v1d(1:N1D+1),   &
                       v2d,v3d)
     call nc_cmor_associate_1d(grid_id,time_axid,height_axid,height_mlev_axid,           &
                               column_axid,sza_axid,temp_axid,channel_axid,dbze_axid,    &
                               sratio_axid,MISR_CTH_axid,tau_axid,pressure2_axid,Nlon,   &
                               Nlat,nPoints,nColumns,nLevels,rttov_nChannels,nLvgrid,    &
                               cospOUT,N1D+1,N2D,N3D,v1d,v2d,v3d)
     call nc_cmor_write_1d(nPoints,lon,lat,time_bnds,lonvar_id,latvar_id,N1D+1,N2D,N3D,  &
                           v1d(1:N1D+1),v2d,v3d)
  endif
  if (geomode .gt. 1) then
     call nc_cmor_init('../cmor/cosp_cmor_nl_2D.txt','replace',nPoints,nColumns,nLevels, &
                       rttov_nChannels,nLvgrid,lon,lat,mgrid_zl,mgrid_zu,mgrid_z,cospOUT,&
                       geomode,Nlon,Nlat,N1D,N2D,N3D,N_OUT_LIST,out_list,lon_axid,       &
                       lat_axid,time_axid,height_axid,height_mlev_axid,grid_id,lonvar_id,&
                       latvar_id,column_axid,sza_axid,temp_axid,channel_axid,dbze_axid,  &
                       sratio_axid,MISR_CTH_axid,tau_axid,pressure2_axid,v1d(1:N1D),v2d, &
                       v3d)
     call nc_cmor_associate_2d(lon_axid,lat_axid,time_axid,height_axid,height_mlev_axid, &
                               column_axid,sza_axid,temp_axid,channel_axid,dbze_axid,    &
                               sratio_axid,MISR_CTH_axid,tau_axid,pressure2_axid,Nlon,   &
                               Nlat,nPoints,nColumns,nLevels,rttov_nChannels,nLvgrid,    &
                               cospOUT,N1D,N2D,N3D,v1d(1:N1D),v2d,v3d)
     call nc_cmor_write_2d(time_bnds,geomode,Nlon,Nlat,N1D,N2D,N3D,v1d(1:N1D),v2d,v3d)
  endif
  deallocate(mgrid_z,mgrid_zu,mgrid_zl)
  call nc_cmor_close()
  call cpu_time(driver_time(8))
  print*,'Time to write to output:  ',driver_time(8)-driver_time(7)

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Free up memory
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  call destroy_cosp_outputs(cospOUT)
  call destroy_cospIN(cospIN)
  call destroy_cospstateIN(cospstateIN)
contains
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  ! SUBROUTINE subsample_and_optics
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  subroutine subsample_and_optics(nPoints,nLevels,nColumns,nHydro,overlap,                  &
       use_precipitation_fluxes,lidar_ice_type,sd,tca,cca,fl_lsrainIN,fl_lssnowIN,          &
       fl_lsgrplIN,fl_ccrainIN,fl_ccsnowIN,mr_lsliq,mr_lsice,mr_ccliq,mr_ccice,reffIN,dtau_c,&
       dtau_s,dem_c,dem_s,cospstateIN,cospIN)
    ! Inputs
    integer,intent(in) :: nPoints,nLevels,nColumns,nHydro,overlap,lidar_ice_type
    real(wp),intent(in),dimension(nPoints,nLevels) :: tca,cca,mr_lsliq,mr_lsice,mr_ccliq,   &
         mr_ccice,dtau_c,dtau_s,dem_c,dem_s,fl_lsrainIN,fl_lssnowIN,fl_lsgrplIN,     &
         fl_ccrainIN,fl_ccsnowIN
    real(wp),intent(in),dimension(nPoints,nLevels,nHydro) :: reffIN
    logical,intent(in) :: use_precipitation_fluxes
    type(size_distribution),intent(inout) :: sd
    
    ! Outputs
    type(cosp_optical_inputs),intent(inout) :: cospIN
    type(cosp_column_inputs),intent(inout)  :: cospstateIN

    ! Local variables
    integer :: i,j,k
    type(rng_state),allocatable,dimension(:) :: rngs  ! Seeds for random number generator
    integer,dimension(:),allocatable :: seed
    real(wp),dimension(:,:),    allocatable :: ls_p_rate,cv_p_rate,frac_ls,frac_cv,prec_ls,prec_cv
    real(wp),dimension(nPoints,nLevels)     :: column_frac_out,column_prec_out,fl_lsrain, &
         fl_lssnow,fl_lsgrpl,fl_ccrain,fl_ccsnow
    real(wp),dimension(:,:,:),  allocatable :: frac_prec,hm_matrix,re_matrix,     &
                                             Np_matrix,MODIS_cloudWater,MODIS_cloudIce,  &
                                             MODIS_watersize,MODIS_iceSize,              &
                                             MODIS_opticalThicknessLiq,                  &
                                             MODIS_opticalThicknessIce
    real(wp),dimension(:,:,:,:),allocatable :: mr_hydro,Reff,Np
    logical :: cmpGases=.true.

    if (Ncolumns .gt. 1) then
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Generate subcolumns for clouds (SCOPS) and precipitation type (PREC_SCOPS)
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! RNG used for subcolumn generation
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
       deallocate(ls_p_rate,cv_p_rate)
       
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Compute precipitation fraction in each gridbox
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Allocate
       allocate(frac_ls(nPoints,nLevels),prec_ls(nPoints,nLevels),                       &
                frac_cv(nPoints,nLevels),prec_cv(nPoints,nLevels))
       
       ! Initialize
       frac_ls(1:nPoints,1:nLevels) = 0._wp
       prec_ls(1:nPoints,1:nLevels) = 0._wp
       frac_cv(1:nPoints,1:nLevels) = 0._wp
       prec_cv(1:nPoints,1:nLevels) = 0._wp
       do j=1,nPoints
          do k=1,nLevels
             do i=1,nColumns
                if (cospIN%frac_out(j,i,nLevels+1-k)  .eq. 1)  frac_ls(j,k) = frac_ls(j,k)+1._wp
                if (cospIN%frac_out(j,i,nLevels+1-k)  .eq. 2)  frac_cv(j,k) = frac_cv(j,k)+1._wp
                if (frac_prec(j,i,nLevels+1-k) .eq. 1)  prec_ls(j,k) = prec_ls(j,k)+1._wp
                if (frac_prec(j,i,nLevels+1-k) .eq. 2)  prec_cv(j,k) = prec_cv(j,k)+1._wp
                if (frac_prec(j,i,nLevels+1-k) .eq. 3)  prec_cv(j,k) = prec_cv(j,k)+1._wp
                if (frac_prec(j,i,nLevels+1-k) .eq. 3)  prec_ls(j,k) = prec_ls(j,k)+1._wp
             enddo
             frac_ls(j,k)=frac_ls(j,k)/nColumns
             frac_cv(j,k)=frac_cv(j,k)/nColumns
             prec_ls(j,k)=prec_ls(j,k)/nColumns
             prec_cv(j,k)=prec_cv(j,k)/nColumns
          enddo
       enddo

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Compute mixing ratios, effective radii and precipitation fluxes for clouds
       ! and precipitation
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       allocate(mr_hydro(nPoints,nColumns,nLevels,nHydro),                               &
                Reff(nPoints,nColumns,nLevels,nHydro),                                   &
                Np(nPoints,nColumns,nLevels,nHydro))

       ! Initialize
       mr_hydro(:,:,:,:) = 0._wp
       Reff(:,:,:,:)    = 0._wp
       Np(:,:,:,:)       = 0._wp
       do k=1,nColumns
          ! Subcolumn cloud fraction
          column_frac_out = cospIN%frac_out(:,k,nLevels:1:-1)
               
          ! LS clouds
          where (column_frac_out == I_LSC)
             mr_hydro(:,k,:,I_LSCLIQ) = mr_lsliq(1:nPoints,Nlevels:1:-1)
             mr_hydro(:,k,:,I_LSCICE) = mr_lsice(1:nPoints,Nlevels:1:-1)
             Reff(:,k,:,I_LSCLIQ)     = ReffIN(:,Nlevels:1:-1,I_LSCLIQ)
             Reff(:,k,:,I_LSCICE)     = ReffIN(:,Nlevels:1:-1,I_LSCICE)
             Np(:,k,:,I_LSCLIQ)       = 0._wp ! Should be inputs
             Np(:,k,:,I_LSCICE)       = 0._wp ! Should be inputs
             ! CONV clouds   
          elsewhere (column_frac_out == I_CVC)
             mr_hydro(:,k,:,I_CVCLIQ) = mr_ccliq(1:nPoints,Nlevels:1:-1)
             mr_hydro(:,k,:,I_CVCICE) = mr_ccice(1:nPoints,Nlevels:1:-1)
             Reff(:,k,:,I_CVCLIQ)     = ReffIN(:,Nlevels:1:-1,I_CVCLIQ)
             Reff(:,k,:,I_CVCICE)     = ReffIN(:,Nlevels:1:-1,I_CVCICE)
             Np(:,k,:,I_CVCLIQ)       = 0._wp ! Should be inputs
             Np(:,k,:,I_CVCICE)       = 0._wp ! Should be inputs
          end where
          
          ! Subcolumn precipitation
          column_prec_out = frac_prec(:,k,nLevels:1:-1)

          ! LS Precipitation
          where ((column_prec_out == 1) .or. (column_prec_out == 3) )
             Reff(:,k,:,I_LSRAIN) = ReffIN(:,Nlevels:1:-1,I_LSRAIN)
             Reff(:,k,:,I_LSSNOW) = ReffIN(:,Nlevels:1:-1,I_LSSNOW)
             Reff(:,k,:,I_LSGRPL) = ReffIN(:,Nlevels:1:-1,I_LSGRPL)
             Np(:,k,:,I_LSRAIN)   = 0._wp! Should be inputs
             Np(:,k,:,I_LSSNOW)   = 0._wp! Should be inputs
             Np(:,k,:,I_LSGRPL)   = 0._wp! Should be inputs
             ! CONV precipitation   
          elsewhere ((column_prec_out == 2) .or. (column_prec_out == 3))
             Reff(:,k,:,I_CVRAIN) = ReffIN(:,Nlevels:1:-1,I_CVRAIN)
             Reff(:,k,:,I_CVSNOW) = ReffIN(:,Nlevels:1:-1,I_CVSNOW)
             Np(:,k,:,I_CVRAIN)   = 0._wp! Should be inputs
             Np(:,k,:,I_CVSNOW)   = 0._wp! Should be inputs
          end where
       enddo
       
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Convert the mixing ratio and precipitation fluxes from gridbox mean to
       ! the fraction-based values
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       do k=1,nLevels
          do j=1,nPoints
             ! Clouds
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
                   fl_lsrain(j,k) = fl_lsrainIN(j,Nlevels-k+1)/prec_ls(j,k)
                   fl_lssnow(j,k) = fl_lssnowIN(j,Nlevels-k+1)/prec_ls(j,k)
                   fl_lsgrpl(j,k) = fl_lsgrplIN(j,Nlevels-k+1)/prec_ls(j,k)
                endif
                if (prec_cv(j,k) .ne. 0.) then
                   fl_ccrain(j,k) = fl_ccrainIN(j,Nlevels-k+1)/prec_cv(j,k)
                   fl_ccsnow(j,k) = fl_ccsnowIN(j,Nlevels-k+1)/prec_cv(j,k)
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
       deallocate(frac_ls,prec_ls,frac_cv,prec_cv)

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Convert precipitation fluxes to mixing ratios
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       if (use_precipitation_fluxes) then

       call cosp_precip_mxratio(nPoints,nLevels, nColumns,                            &
               cospstateIN%pfull(:,Nlevels:1:-1),                    &
               cospstateIN%at(:,Nlevels:1:-1),                       &
               frac_prec(:,:,Nlevels:1:-1), 1._wp,                   &
               n_ax(I_LSRAIN), n_bx(I_LSRAIN),   alpha_x(I_LSRAIN),  &
               c_x(I_LSRAIN),   d_x(I_LSRAIN),   g_x(I_LSRAIN),      &
               a_x(I_LSRAIN),   b_x(I_LSRAIN),   gamma_1(I_LSRAIN),  &
               gamma_2(I_LSRAIN),gamma_3(I_LSRAIN),gamma_4(I_LSRAIN),&
               fl_lsrain,mr_hydro(:,:,:,I_LSRAIN),Reff(:,:,:,I_LSRAIN))
          call cosp_precip_mxratio(nPoints,nLevels,nColumns,                             &
               cospstateIN%pfull(:,Nlevels:1:-1),                    &
               cospstateIN%at(:,Nlevels:1:-1),                       &
               frac_prec(:,:,Nlevels:1:-1), 1._wp,                   &
               n_ax(I_LSSNOW),  n_bx(I_LSSNOW),  alpha_x(I_LSSNOW),  &
               c_x(I_LSSNOW),   d_x(I_LSSNOW),   g_x(I_LSSNOW),      &
               a_x(I_LSSNOW),   b_x(I_LSSNOW),   gamma_1(I_LSSNOW),  &
               gamma_2(I_LSSNOW),gamma_3(I_LSSNOW),gamma_4(I_LSSNOW),&
               fl_lssnow,mr_hydro(:,:,:,I_LSSNOW),Reff(:,:,:,I_LSSNOW))
          call cosp_precip_mxratio(nPoints,nLevels,nColumns,                             &
               cospstateIN%pfull(:,Nlevels:1:-1),                    &
               cospstateIN%at(:,Nlevels:1:-1),                       &
               frac_prec(:,:,Nlevels:1:-1), 2._wp,                   &
               n_ax(I_CVRAIN),  n_bx(I_CVRAIN),  alpha_x(I_CVRAIN),  &
               c_x(I_CVRAIN),   d_x(I_CVRAIN),   g_x(I_CVRAIN),      &
               a_x(I_CVRAIN),   b_x(I_CVRAIN),   gamma_1(I_CVRAIN),  &
               gamma_2(I_CVRAIN),gamma_3(I_CVRAIN),gamma_4(I_CVRAIN),&
               fl_ccrain,mr_hydro(:,:,:,I_CVRAIN),Reff(:,:,:,I_CVRAIN))
          call cosp_precip_mxratio(nPoints,nLevels,nColumns,                             &
               cospstateIN%pfull(:,Nlevels:1:-1),                    &
               cospstateIN%at(:,Nlevels:1:-1),                       &
               frac_prec(:,:,Nlevels:1:-1), 2._wp,                   &      
               n_ax(I_CVSNOW),  n_bx(I_CVSNOW),  alpha_x(I_CVSNOW),  &
               c_x(I_CVSNOW),   d_x(I_CVSNOW),   g_x(I_CVSNOW),      &
               a_x(I_CVSNOW),   b_x(I_CVSNOW),   gamma_1(I_CVSNOW),  &
               gamma_2(I_CVSNOW),gamma_3(I_CVSNOW),gamma_4(I_CVSNOW),&
               fl_ccsnow,mr_hydro(:,:,:,I_CVSNOW),Reff(:,:,:,I_CVSNOW))
          call cosp_precip_mxratio(nPoints,nLevels,nColumns,                             &
               cospstateIN%pfull(:,Nlevels:1:-1),                    &
               cospstateIN%at(:,Nlevels:1:-1),                       &
               frac_prec(:,:,Nlevels:1:-1), 1._wp,                   &       
               n_ax(I_LSGRPL),  n_bx(I_LSGRPL),  alpha_x(I_LSGRPL),  &
               c_x(I_LSGRPL),   d_x(I_LSGRPL),   g_x(I_LSGRPL),      &
               a_x(I_LSGRPL),   b_x(I_LSGRPL),   gamma_1(I_LSGRPL),  &
               gamma_2(I_LSGRPL),gamma_3(I_LSGRPL),gamma_4(I_LSGRPL),&
               fl_lsgrpl,mr_hydro(:,:,:,I_LSGRPL),Reff(:,:,:,I_LSGRPL))
          deallocate(frac_prec)
       endif

    else
       cospIN%frac_out(:,:,:) = 1  
       allocate(mr_hydro(nPoints, 1,nLevels,nHydro),Reff(nPoints,1,nLevels,nHydro),  &
            Np(nPoints,1,nLevels,nHydro))
       mr_hydro(:,1,:,I_LSCLIQ) = mr_lsliq(1:nPoints,1:nLevels:-1)
       mr_hydro(:,1,:,I_LSCICE) = mr_lsice(1:nPoints,1:nLevels:-1)
       mr_hydro(:,1,:,I_CVCLIQ) = mr_ccliq(1:nPoints,1:nLevels:-1)
       mr_hydro(:,1,:,I_CVCICE) = mr_ccice(1:nPoints,1:nLevels:-1)
       Reff(:,1,:,:)            = ReffIN(:,1:nLevels:-1,:)
       Np(:,1,:,:)              = 0._wp! Should be inputs
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
    call lidar_optics(nPoints,nColumns,nLevels,4,lidar_ice_type,                           &
                      mr_hydro(:,:,nLevels:1:-1,I_LSCLIQ),                                 &
                      mr_hydro(:,:,nLevels:1:-1,I_LSCICE),                                 &
                      mr_hydro(:,:,nLevels:1:-1,I_CVCLIQ),                                 &
                      mr_hydro(:,:,nLevels:1:-1,I_CVCICE),                                 &
                      ReffIN(:,:,I_LSCLIQ),ReffIN(:,:,I_LSCICE),     &
                      ReffIN(:,:,I_CVCLIQ),ReffIN(:,:,I_CVCICE),     & 
                      cospstateIN%pfull,cospstateIN%phalf,cospstateIN%at,                  &
                      cospIN%beta_mol,cospIN%betatot,cospIN%taupart,                       &
                      cospIN%tau_mol,cospIN%tautot,cospIN%tautot_S_liq,                    &
                      cospIN%tautot_S_ice, betatot_ice = cospIN%betatot_ice,               &
                      betatot_liq=cospIN%betatot_liq,tautot_ice=cospIN%tautot_ice,         &
                      tautot_liq = cospIN%tautot_liq)
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! CLOUDSAT RADAR OPTICS
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    allocate(hm_matrix(nHydro,nPoints,nLevels),re_matrix(nHydro,nPoints,nLevels),          &
             Np_matrix(nHydro,nPoints,nLevels))           
    
    ! Loop over all subcolumns
    do k=1,nColumns
       do i=1,nHydro
          hm_matrix(i,1:nPoints,nLevels:1:-1) = mr_hydro(:,k,:,i)*1000._wp
          re_matrix(i,1:nPoints,nLevels:1:-1) = Reff(:,k,:,i)*1.e6_wp  
          Np_matrix(i,1:nPoints,nLevels:1:-1) = Np(:,k,:,i)       
       enddo

       call quickbeam_optics(sd,cospIN%rcfg_cloudsat,nPoints,nLevels,R_UNDEF,hm_matrix,    &
            re_matrix,Np_matrix,cospstateIN%pfull,cospstateIN%at,cospstateIN%qv,cmpGases,  &
            cospIN%z_vol_cloudsat(1:nPoints,k,:),cospIN%kr_vol_cloudsat(1:nPoints,k,:),    &
            cospIN%g_vol_cloudsat(1:nPoints,k,:))
    enddo
    deallocate(hm_matrix,re_matrix,Np_matrix)
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! MODIS optics
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    allocate(MODIS_cloudWater(nPoints,nColumns,nLevels),                                   &
             MODIS_cloudIce(nPoints,nColumns,nLevels),                                     &
             MODIS_waterSize(nPoints,nColumns,nLevels),                                    &
             MODIS_iceSize(nPoints,nColumns,nLevels),                                      &
             MODIS_opticalThicknessLiq(nPoints,nColumns,nLevels),                          &
             MODIS_opticalThicknessIce(nPoints,nColumns,nLevels))
    ! Cloud water
    call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,                   &
                               mr_hydro(:,:,nLevels:1:-1,I_CVCLIQ),                        &
                               mr_hydro(:,:,nLevels:1:-1,I_LSCLIQ),MODIS_cloudWater)
    ! Cloud ice
    call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,                   &
                               mr_hydro(:,:,nLevels:1:-1,I_CVCICE),                        &
                               mr_hydro(:,:,nLevels:1:-1,I_LSCICE),MODIS_cloudIce)  
    ! Water droplet size
    call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,                   &
                               Reff(:,:,nLevels:1:-1,I_CVCLIQ),                            &
                               Reff(:,:,nLevels:1:-1,I_LSCLIQ),MODIS_waterSize)
    ! Ice crystal size
    call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,                   &
                                Reff(:,:,nLevels:1:-1,I_CVCICE),                           &
                                Reff(:,:,nLevels:1:-1,I_LSCICE),MODIS_iceSize)
    
    ! Partition optical thickness into liquid and ice parts
    call modis_optics_partition(nPoints,nLevels,nColumns,MODIS_cloudWater,MODIS_cloudIce,  &
                                MODIS_waterSize,MODIS_iceSize,cospIN%tau_067,              &
                                MODIS_opticalThicknessLiq,MODIS_opticalThicknessIce)
    
    ! Compute assymetry parameter and single scattering albedo 
    call modis_optics(nPoints,nLevels,nColumns,num_trial_res,MODIS_opticalThicknessLiq,    &
                      MODIS_waterSize*1.0e6_wp,MODIS_opticalThicknessIce,                  &
                      MODIS_iceSize*1.0e6_wp,cospIN%fracLiq, cospIN%asym, cospIN%ss_alb)
    
    ! Deallocate memory
    deallocate(MODIS_cloudWater,MODIS_cloudIce,MODIS_WaterSize,MODIS_iceSize,              &
               MODIS_opticalThicknessLiq,MODIS_opticalThicknessIce,mr_hydro,               &
               Np,Reff)

  end subroutine subsample_and_optics
end program cosp_test_v2

