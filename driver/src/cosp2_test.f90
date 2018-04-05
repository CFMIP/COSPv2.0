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
! April 2018 - R. Guzman - Added OPAQ and Ground LIDar (GLID) diagnostics
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
program cosp2_test
  use cosp_kinds,          only: wp                         
  USE MOD_COSP_CONFIG,     ONLY: R_UNDEF,PARASOL_NREFL,LIDAR_NCAT,LIDAR_NTYPE,SR_BINS,    & !OPAQ
                                 N_HYDRO,RTTOV_MAX_CHANNELS,numMISRHgtBins,               &
                                 DBZE_BINS,LIDAR_NTEMP,calipso_histBsct,                  &
                                 numMODISTauBins,numMODISPresBins,                        &
                                 numMODISReffIceBins,numMODISReffLiqBins,                 &
                                 numISCCPTauBins,numISCCPPresBins,numMISRTauBins,         &
                                 ntau,modis_histTau,tau_binBounds,                        &
                                 modis_histTauEdges,tau_binEdges,                         &
                                 modis_histTauCenters,tau_binCenters,ntauV1p4,            &
                                 tau_binBoundsV1p4,tau_binEdgesV1p4, tau_binCentersV1p4,  &
                                 groundlidar_histBsct,                                    & !GLID
                                 Nlvgrid_local  => Nlvgrid,                               & !OPAQ
                                 vgrid_z_local  => vgrid_z                                  !OPAQ

  use cosp_phys_constants, only: amw,amd,amO3,amCO2,amCH4,amN2O,amCO
  use mod_cosp_io,         only: nc_read_input_file,write_cosp2_output
  USE mod_quickbeam_optics,only: size_distribution,hydro_class_init,quickbeam_optics,     &
                                 quickbeam_optics_init,gases
  use quickbeam,           only: radar_cfg
  use mod_cosp,            only: cosp_init,cosp_optical_inputs,cosp_column_inputs,        &
                                 cosp_outputs,cosp_cleanUp,cosp_simulator
  USE mod_rng,             ONLY: rng_state, init_rng
  USE mod_scops,           ONLY: scops
  USE mod_prec_scops,      ONLY: prec_scops
  USE MOD_COSP_UTILS,      ONLY: cosp_precip_mxratio
  use cosp_optics,         ONLY: cosp_simulator_optics,lidar_optics,modis_optics,         &
                                 modis_optics_partition,groundlidar_optics                  !GLID

  implicit none

  ! Input/Output driver file control
  character(len=64),parameter :: &
       cosp_input_namelist  = 'cosp2_input_nl.txt', &
       cosp_output_namelist = 'cosp2_output_nl.txt'

  ! Test data
  integer :: &
       Nlon,Nlat,geomode
  real(wp) :: &
       emsfc_lw
  real(wp),dimension(:),allocatable,target:: &
       lon,       & ! Longitude (deg)
       lat,       & ! Latitude (deg)
       skt,       & ! Skin temperature (K)
       surfelev,  & ! Surface Elevation (m) !TIBO2
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
       rttov_Nchannels              ! RTTOV: Number of channels to be computed
  real(wp),dimension(:),allocatable :: &                                      !OPAQ
       vgrid_z                      ! mid-level altitude of the vertical grid !OPAQ
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
  character(len=64) :: &
       finput                       ! Input NetCDF file
  character(len=256) :: &
       foutput
  character(len=512) :: &
       dinput                       ! Directory where the input files are located
  character(len=600) :: &
       fileIN                       ! dinput+finput
  namelist/COSP_INPUT/overlap, isccp_topheight, isccp_topheight_direction, npoints,      &
       npoints_it, ncolumns, nlevels, use_vgrid, Nlvgrid, csat_vgrid, dinput, finput,    &
       foutput, cloudsat_radar_freq, surface_radar, cloudsat_use_gas_abs,cloudsat_do_ray,&
       cloudsat_k2, cloudsat_micro_scheme, lidar_ice_type, use_precipitation_fluxes,     &
       rttov_platform, rttov_satellite, rttov_Instrument, rttov_Nchannels,               &
       rttov_Channels, rttov_Surfem, rttov_ZenAng, co2, ch4, n2o, co

  ! Output namelist
  logical :: Lcfaddbze94,Ldbze94,Latb532,LcfadLidarsr532,Lclcalipso,Lclhcalipso,         &
             Lcllcalipso,Lclmcalipso,Lcltcalipso,LparasolRefl,Lclcalipsoliq,             &
             Lclcalipsoice,Lclcalipsoun,Lclcalipsotmp,Lclcalipsotmpliq,Lclcalipsotmpice, &
             Lclcalipsotmpun,Lclhcalipsoliq,Lcllcalipsoliq,Lclmcalipsoliq,Lcltcalipsoliq,& !OPAQ
             Lclhcalipsoice,Lcllcalipsoice,Lclmcalipsoice,Lcltcalipsoice,Lclhcalipsoun,  & !OPAQ
             Lcllcalipsoun,Lclmcalipsoun,Lcltcalipsoun,Lclopaquecalipso,Lclthincalipso,  & !OPAQ
             Lclzopaquecalipso,Lclcalipsoopaque,Lclcalipsothin,Lclcalipsozopaque,        & !OPAQ
             Lclcalipsoopacity,Lclopaquetemp,Lclthintemp,Lclzopaquetemp,Lclopaquemeanz,  & !TIBO
             Lclthinmeanz,Lclthinemis,Lclopaquemeanzse,Lclthinmeanzse,Lclzopaquecalipsose,& !TIBO2
             LlidarBetaMol532gr,LcfadLidarsr532gr,Latb532gr,Lclgroundlidar,              & !GLID
             Lclhgroundlidar,Lcllgroundlidar,Lclmgroundlidar,Lcltgroundlidar,            & !GLID
             Lalbisccp,Lboxptopisccp,Lboxtauisccp,Lpctisccp,Lclisccp,Ltauisccp,Lcltisccp,&
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
                       Lcllcalipsoun,Lclmcalipsoun,Lcltcalipsoun,Lclopaquecalipso,       & !OPAQ
                       Lclthincalipso,Lclzopaquecalipso,Lclcalipsoopaque,Lclcalipsothin, & !OPAQ
                       Lclcalipsozopaque,Lclcalipsoopacity,Lclopaquetemp,Lclthintemp,    & !TIBO
                       Lclzopaquetemp,Lclopaquemeanz,Lclthinmeanz,Lclthinemis,           & !TIBO
                       Lclopaquemeanzse,Lclthinmeanzse,Lclzopaquecalipsose,              & !TIBO2
                       LlidarBetaMol532gr,LcfadLidarsr532gr,Latb532gr,Lclgroundlidar,    & !GLID
                       Lclhgroundlidar,Lcllgroundlidar,Lclmgroundlidar,Lcltgroundlidar,  & !GLID
                       Lalbisccp,Lboxptopisccp,                                          &
                       Lboxtauisccp,Lpctisccp,Lclisccp,Ltauisccp,Lcltisccp,Lmeantbisccp, &
                       Lmeantbclrisccp,LclMISR,Lclcalipso2,Lcltlidarradar,Lfracout,      &
                       LlidarBetaMol532,Lcltmodis,Lclwmodis,Lclimodis,Lclhmodis,         &
                       Lclmmodis,Lcllmodis,Ltautmodis,Ltauwmodis,Ltauimodis,             &
                       Ltautlogmodis,Ltauwlogmodis,Ltauilogmodis,Lreffclwmodis,          &
                       Lreffclimodis,Lpctmodis,Llwpmodis,Liwpmodis,Lclmodis,Ltbrttov

  ! Local variables
  logical :: &
       lsingle   = .true.,  & ! True if using MMF_v3_single_moment CLOUDSAT microphysical scheme (default)
       ldouble   = .false., & ! True if using MMF_v3.5_two_moment CLOUDSAT microphysical scheme
       lisccp    = .false. ,& ! Local on/off switch for simulators (used by initialization)
       lmodis    = .false., & !
       lmisr     = .false., & !
       lcalipso  = .false., & !
       lgroundlidar = .false., & !GLID
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
       n_out_list = 86,           & ! Number of possible output variables !OPAQ !N1D+1 !GLID !TIBO !TIBO2
       N3D        = 9,            & ! Number of 3D output variables       !OPAQ !GLID
       N2D        = 20,           & ! Number of 2D output variables       !OPAQ !GLID
       N1D        = 56              ! Number of 1D output variables       !OPAQ !GLID !TIBO !TIBO2
  character(len=32),dimension(n_out_list) :: out_list  ! List of output variable names
  integer :: lon_axid,time_axid,height_axid,height_mlev_axid,grid_id,lonvar_id,       &
             latvar_id,column_axid,sza_axid,temp_axid,channel_axid,dbze_axid,sratio_axid,&
             MISR_CTH_axid,lat_axid,tau_axid,pressure2_axid 
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
           frac_out(Npoints,Ncolumns,Nlevels),surfelev(Npoints)) !TIBO2

  fileIN = trim(dinput)//trim(finput)
  call nc_read_input_file(fileIN,Npoints,Nlevels,N_HYDRO,lon,lat,p,ph,zlev,zlev_half,    &
                          T,sh,rh,tca,cca,mr_lsliq,mr_lsice,mr_ccliq,mr_ccice,fl_lsrain, &
                          fl_lssnow,fl_lsgrpl,fl_ccrain,fl_ccsnow,Reff,dtau_s,dtau_c,    &
                          dem_s,dem_c,skt,landmask,mr_ozone,u_wind,v_wind,sunlit,        &
                          emsfc_lw,geomode,Nlon,Nlat,surfelev) !TIBO2
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
       Lcltlidarradar .or. Lclopaquecalipso .or. Lclthincalipso .or. Lclzopaquecalipso   & !OPAQ
       .or. Lclcalipsoopaque .or. Lclcalipsothin .or. Lclcalipsozopaque .or.             & !OPAQ
       Lclcalipsoopacity .or. Lclopaquetemp .or. Lclthintemp .or. Lclzopaquetemp .or.    & !OPAQ
       Lclopaquemeanz .or. Lclthinmeanz .or. Lclthinemis .or. Lclopaquemeanzse .or.      & !TIBO
       Lclthinmeanzse .or. Lclzopaquecalipsose) Lcalipso = .true.              !OPAQ !TIBO !TIBO2

  if (LlidarBetaMol532gr .or. LcfadLidarsr532gr .or. Latb532gr .or. Lclgroundlidar .or.  & !GLID
       Lclhgroundlidar .or. Lcllgroundlidar .or. Lclmgroundlidar .or. Lcltgroundlidar)   & !GLID
       Lgroundlidar = .true.                                                               !GLID

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
  call quickbeam_optics_init()

  ! Initialize the distributional parameters for hydrometeors in radar simulator
  call hydro_class_init(lsingle,ldouble,sd)

  ! Initialize COSP simulator
  call COSP_INIT(Lisccp, Lmodis, Lmisr, Lcloudsat, Lcalipso, Lgroundlidar, Lparasol, Lrttov, &
       Npoints, Nlevels, cloudsat_radar_freq, cloudsat_k2, cloudsat_use_gas_abs,         &
       cloudsat_do_ray, isccp_topheight, isccp_topheight_direction, surface_radar,       &
       rcfg_cloudsat, rttov_Nchannels, rttov_Channels, rttov_platform,                   &
       rttov_satellite, rttov_instrument, use_vgrid, csat_vgrid, Nlvgrid,                &
       vgrid_z_local, cloudsat_micro_scheme, cospOUT)                                      !OPAQ
  call cpu_time(driver_time(3))
  
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
       Latb532gr, LlidarBetaMol532, LlidarBetaMol532gr, LcfadLidarsr532, LcfadLidarsr532gr, & !GLID
       Lclcalipso2, Lclcalipso, Lclgroundlidar, Lclhcalipso, Lcllcalipso, Lclmcalipso,   & !GLID
       Lcltcalipso, Lclhgroundlidar, Lcllgroundlidar, Lclmgroundlidar, Lcltgroundlidar,  & !GLID
       Lcltlidarradar, Lclcalipsoliq,                                                    &
       Lclcalipsoice, Lclcalipsoun, Lclcalipsotmp, Lclcalipsotmpliq, Lclcalipsotmpice,   &
       Lclcalipsotmpun, Lcltcalipsoliq, Lcltcalipsoice, Lcltcalipsoun, Lclhcalipsoliq,   &
       Lclhcalipsoice, Lclhcalipsoun, Lclmcalipsoliq, Lclmcalipsoice, Lclmcalipsoun,     &
       Lcllcalipsoliq, Lcllcalipsoice, Lcllcalipsoun, Lclopaquecalipso, Lclthincalipso,  & !OPAQ
       Lclzopaquecalipso, Lclcalipsoopaque, Lclcalipsothin, Lclcalipsozopaque,           & !OPAQ
       Lclcalipsoopacity, Lclopaquetemp, Lclthintemp, Lclzopaquetemp, Lclopaquemeanz,    & !TIBO
       Lclthinmeanz, Lclthinemis, Lclopaquemeanzse, Lclthinmeanzse, Lclzopaquecalipsose, & !TIBO !TIBO2
       LcfadDbze94, Ldbze94, Lparasolrefl,                                               &
       Ltbrttov, Npoints, Ncolumns, Nlevels, Nlvgrid_local, rttov_Nchannels, cospOUT)

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
     cospIN%emsfc_lw         = emsfc_lw
     cospIN%rcfg_cloudsat    = rcfg_cloudsat
     cospstateIN%hgt_matrix  = zlev(start_idx:end_idx,Nlevels:1:-1) ! km
     cospstateIN%sunlit      = sunlit(start_idx:end_idx)            ! 0-1
     cospstateIN%skt         = skt(start_idx:end_idx)               ! K
     cospstateIN%surfelev    = surfelev(start_idx:end_idx)          ! m !TIBO2
     cospstateIN%land        = landmask(start_idx:end_idx)          ! 0-1 (*note* model specific)
     cospstateIN%qv          = sh(start_idx:end_idx,Nlevels:1:-1)   ! kg/kg
     cospstateIN%at          = T(start_idx:end_idx,Nlevels:1:-1)    ! K
     cospstateIN%pfull       = p(start_idx:end_idx,Nlevels:1:-1)    ! Pa 
     ! Pressure at interface (nlevels+1). Set uppermost interface to 0.
     cospstateIN%phalf(:,2:Nlevels+1) = ph(start_idx:end_idx,Nlevels:1:-1)   ! Pa  
     cospstateIN%phalf(:,1)           = 0._wp
     ! Height at interface (nlevels+1). Set lowermost interface to 0.
     cospstateIN%hgt_matrix_half(:,1:Nlevels) = zlev_half(start_idx:end_idx,Nlevels:1:-1) ! km
     cospstateIN%hgt_matrix_half(:,Nlevels+1) = 0._wp
     
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! Generate subcolumns and compute optical inputs.
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
  call write_cosp2_output(Npoints, Ncolumns, Nlevels, zlev(1,Nlevels:1:-1), lon, lat, cospOUT, foutput)

  call cpu_time(driver_time(8))
  print*,'Time to write to output:  ',driver_time(8)-driver_time(7)

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Free up memory
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  call destroy_cosp_outputs(cospOUT)
  call destroy_cospIN(cospIN)
  call destroy_cospstateIN(cospstateIN)
  call cosp_cleanUp()
contains
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  ! SUBROUTINE subsample_and_optics
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  subroutine subsample_and_optics(nPoints, nLevels, nColumns, nHydro, overlap,              &
       use_precipitation_fluxes, lidar_ice_type, sd, tca, cca, fl_lsrainIN, fl_lssnowIN,    &
       fl_lsgrplIN, fl_ccrainIN, fl_ccsnowIN, mr_lsliq, mr_lsice, mr_ccliq, mr_ccice,       &
       reffIN, dtau_c, dtau_s, dem_c, dem_s, cospstateIN, cospIN)
    ! Inputs
    integer,intent(in) :: nPoints, nLevels, nColumns, nHydro, overlap, lidar_ice_type
    real(wp),intent(in),dimension(nPoints,nLevels) :: tca,cca,mr_lsliq,mr_lsice,mr_ccliq,   &
         mr_ccice,dtau_c,dtau_s,dem_c,dem_s,fl_lsrainIN,fl_lssnowIN,fl_lsgrplIN,fl_ccrainIN,&
         fl_ccsnowIN
    real(wp),intent(in),dimension(nPoints,nLevels,nHydro) :: reffIN
    logical,intent(in) :: use_precipitation_fluxes
    type(size_distribution),intent(inout) :: sd
    
    ! Outputs
    type(cosp_optical_inputs),intent(inout) :: cospIN
    type(cosp_column_inputs),intent(inout)  :: cospstateIN

    ! Local variables
    type(rng_state),allocatable,dimension(:) :: rngs  ! Seeds for random number generator
    integer,dimension(:),allocatable :: seed
    integer :: i,j,k
    real(wp),dimension(:,:), allocatable :: &
         ls_p_rate, cv_p_rate, frac_ls, frac_cv, prec_ls, prec_cv,g_vol
    real(wp),dimension(:,:,:),  allocatable :: &
         frac_prec, MODIS_cloudWater, MODIS_cloudIce,       &
         MODIS_watersize,MODIS_iceSize, MODIS_opticalThicknessLiq,MODIS_opticalThicknessIce
    real(wp),dimension(:,:,:,:),allocatable :: &
         mr_hydro, Reff, Np
    real(wp),dimension(nPoints,nLevels) :: &
         column_frac_out, column_prec_out, fl_lsrain, fl_lssnow, fl_lsgrpl, fl_ccrain, fl_ccsnow
!    logical :: cmpGases=.true.

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
       ! Compute fraction in each gridbox for precipitation  and cloud type.
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
                if (cospIN%frac_out(j,i,k)  .eq. 1)  frac_ls(j,k) = frac_ls(j,k)+1._wp
                if (cospIN%frac_out(j,i,k)  .eq. 2)  frac_cv(j,k) = frac_cv(j,k)+1._wp
                if (frac_prec(j,i,k) .eq. 1)  prec_ls(j,k) = prec_ls(j,k)+1._wp
                if (frac_prec(j,i,k) .eq. 2)  prec_cv(j,k) = prec_cv(j,k)+1._wp
                if (frac_prec(j,i,k) .eq. 3)  prec_cv(j,k) = prec_cv(j,k)+1._wp
                if (frac_prec(j,i,k) .eq. 3)  prec_ls(j,k) = prec_ls(j,k)+1._wp
             enddo
             frac_ls(j,k)=frac_ls(j,k)/nColumns
             frac_cv(j,k)=frac_cv(j,k)/nColumns
             prec_ls(j,k)=prec_ls(j,k)/nColumns
             prec_cv(j,k)=prec_cv(j,k)/nColumns
          enddo
       enddo       

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
          ! Subcolumn cloud fraction
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
          
          ! Subcolumn precipitation
          column_prec_out = frac_prec(:,k,:)

          ! LS Precipitation
          where ((column_prec_out == 1) .or. (column_prec_out == 3) )
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
       deallocate(frac_ls,prec_ls,frac_cv,prec_cv)

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Convert precipitation fluxes to mixing ratios
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       if (use_precipitation_fluxes) then
          ! LS rain
          call cosp_precip_mxratio(nPoints, nLevels, nColumns, cospstateIN%pfull,        &
               cospstateIN%at, frac_prec, 1._wp, n_ax(I_LSRAIN), n_bx(I_LSRAIN),         &
               alpha_x(I_LSRAIN), c_x(I_LSRAIN),   d_x(I_LSRAIN),   g_x(I_LSRAIN),       &
               a_x(I_LSRAIN),   b_x(I_LSRAIN),   gamma_1(I_LSRAIN), gamma_2(I_LSRAIN),   &
               gamma_3(I_LSRAIN), gamma_4(I_LSRAIN), fl_lsrain,                          &
               mr_hydro(:,:,:,I_LSRAIN), Reff(:,:,:,I_LSRAIN))
          ! LS snow
          call cosp_precip_mxratio(nPoints, nLevels, nColumns, cospstateIN%pfull,        &
               cospstateIN%at, frac_prec, 1._wp,  n_ax(I_LSSNOW),  n_bx(I_LSSNOW),       &
               alpha_x(I_LSSNOW), c_x(I_LSSNOW),  d_x(I_LSSNOW),  g_x(I_LSSNOW),         &
               a_x(I_LSSNOW),   b_x(I_LSSNOW),   gamma_1(I_LSSNOW),  gamma_2(I_LSSNOW),  &
               gamma_3(I_LSSNOW), gamma_4(I_LSSNOW), fl_lssnow,                          &
               mr_hydro(:,:,:,I_LSSNOW), Reff(:,:,:,I_LSSNOW))
          ! CV rain
          call cosp_precip_mxratio(nPoints, nLevels, nColumns, cospstateIN%pfull,        &
               cospstateIN%at, frac_prec, 2._wp, n_ax(I_CVRAIN),  n_bx(I_CVRAIN),        &
               alpha_x(I_CVRAIN), c_x(I_CVRAIN),   d_x(I_CVRAIN),   g_x(I_CVRAIN),       &
               a_x(I_CVRAIN),   b_x(I_CVRAIN),   gamma_1(I_CVRAIN), gamma_2(I_CVRAIN),   &
               gamma_3(I_CVRAIN), gamma_4(I_CVRAIN), fl_ccrain,                          &
               mr_hydro(:,:,:,I_CVRAIN), Reff(:,:,:,I_CVRAIN))
          ! CV snow
          call cosp_precip_mxratio(nPoints, nLevels, nColumns, cospstateIN%pfull,        &
               cospstateIN%at, frac_prec, 2._wp, n_ax(I_CVSNOW),  n_bx(I_CVSNOW),        &
               alpha_x(I_CVSNOW),  c_x(I_CVSNOW),   d_x(I_CVSNOW),   g_x(I_CVSNOW),      &
               a_x(I_CVSNOW),   b_x(I_CVSNOW),   gamma_1(I_CVSNOW), gamma_2(I_CVSNOW),   &
               gamma_3(I_CVSNOW), gamma_4(I_CVSNOW), fl_ccsnow,                          &
               mr_hydro(:,:,:,I_CVSNOW), Reff(:,:,:,I_CVSNOW))
          ! LS groupel.
          call cosp_precip_mxratio(nPoints, nLevels, nColumns, cospstateIN%pfull,        &
               cospstateIN%at, frac_prec, 1._wp, n_ax(I_LSGRPL),  n_bx(I_LSGRPL),        &
               alpha_x(I_LSGRPL), c_x(I_LSGRPL),   d_x(I_LSGRPL),   g_x(I_LSGRPL),       &
               a_x(I_LSGRPL),   b_x(I_LSGRPL),   gamma_1(I_LSGRPL),  gamma_2(I_LSGRPL),  &
               gamma_3(I_LSGRPL), gamma_4(I_LSGRPL), fl_lsgrpl,                          &
               mr_hydro(:,:,:,I_LSGRPL), Reff(:,:,:,I_LSGRPL))
          deallocate(frac_prec)
       endif

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
       call lidar_optics(nPoints, nColumns, nLevels, 4, lidar_ice_type,                       &
            mr_hydro(:,:,:,I_LSCLIQ),  mr_hydro(:,:,:,I_LSCICE), mr_hydro(:,:,:,I_CVCLIQ),    &
            mr_hydro(:,:,:,I_CVCICE), ReffIN(:,:,I_LSCLIQ), ReffIN(:,:,I_LSCICE),             &
            ReffIN(:,:,I_CVCLIQ), ReffIN(:,:,I_CVCICE), cospstateIN%pfull, cospstateIN%phalf, &
            cospstateIN%at, cospIN%beta_mol, cospIN%betatot, cospIN%tau_mol, cospIN%tautot,   &
            cospIN%tautot_S_liq, cospIN%tautot_S_ice, cospIN%betatot_ice, cospIN%betatot_liq, &
            cospIN%tautot_ice, cospIN%tautot_liq)
    endif

    if (Lgroundlidar) then
       call groundlidar_optics(nPoints,nColumns,nLevels,4,lidar_ice_type,            & !GLID
                      mr_hydro(:,:,:,I_LSCLIQ),                           & !GLID
                      mr_hydro(:,:,:,I_LSCICE),                           & !GLID
                      mr_hydro(:,:,:,I_CVCLIQ),                           & !GLID
                      mr_hydro(:,:,:,I_CVCICE),                           & !GLID
                      ReffIN(:,:,I_LSCLIQ),ReffIN(:,:,I_LSCICE),                     & !GLID
                      ReffIN(:,:,I_CVCLIQ),ReffIN(:,:,I_CVCICE),                     & !GLID
                      cospstateIN%pfull,cospstateIN%phalf,cospstateIN%at,            & !GLID
                      cospIN%beta_mol_gr,cospIN%betatot_gr,                          & !GLID
                      cospIN%tau_mol_gr,cospIN%tautot_gr)                              !GLID
    endif
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! CLOUDSAT RADAR OPTICS
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (lcloudsat) then

       ! Compute gaseous absorption (assume identical for each subcolun)
       allocate(g_vol(nPoints,nLevels))
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
       do k=1,nColumns
          call quickbeam_optics(sd, rcfg_cloudsat, nPoints, nLevels, R_UNDEF,  &
               mr_hydro(:,k,:,1:nHydro)*1000._wp, Reff(:,k,:,1:nHydro)*1.e6_wp,&
               Np(:,k,:,1:nHydro), cospstateIN%pfull, cospstateIN%at,          &
               cospstateIN%qv, cospIN%z_vol_cloudsat(1:nPoints,k,:),           &
               cospIN%kr_vol_cloudsat(1:nPoints,k,:))
       enddo
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
       
       ! Deallocate memory
       deallocate(MODIS_cloudWater,MODIS_cloudIce,MODIS_WaterSize,MODIS_iceSize,           &
            MODIS_opticalThicknessLiq,MODIS_opticalThicknessIce,mr_hydro,                  &
            Np,Reff)
    endif
  end subroutine subsample_and_optics
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE construct_cospIN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine construct_cospIN(npoints,ncolumns,nlevels,y)
    ! Inputs
    integer,intent(in) :: &
         npoints,  & ! Number of horizontal gridpoints
         ncolumns, & ! Number of subcolumns
         nlevels     ! Number of vertical levels
    ! Outputs 
    type(cosp_optical_inputs),intent(out) :: y
    
    ! Dimensions
    y%Npoints  = Npoints
    y%Ncolumns = Ncolumns
    y%Nlevels  = Nlevels
    y%Npart    = 4
    y%Nrefl    = PARASOL_NREFL
    allocate(y%frac_out(npoints,       ncolumns,nlevels))

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
                y%tautot_S_liq(npoints,   ncolumns        ))
    endif

    if (Lgroundlidar) then                                    !GLID
       allocate(y%beta_mol_gr(npoints,             nlevels),& !GLID
                y%betatot_gr(npoints,     ncolumns,nlevels),& !GLID
                y%tau_mol_gr(npoints,              nlevels),& !GLID
                y%tautot_gr(npoints,      ncolumns,nlevels))  !GLID
    endif                                                     !GLID
    if (Lcloudsat) then
       allocate(y%z_vol_cloudsat(npoints, ncolumns,nlevels),&
                y%kr_vol_cloudsat(npoints,ncolumns,nlevels),&
                y%g_vol_cloudsat(npoints, ncolumns,nlevels))
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
             y%cloudIce(nPoints,nLevels),y%cloudLiq(nPoints,nLevels),y%surfelev(npoints),& !TIBO2
             y%fl_snow(nPoints,nLevels),y%fl_rain(nPoints,nLevels),y%seaice(npoints),    &
             y%tca(nPoints,nLevels),y%hgt_matrix_half(npoints,nlevels+1))

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
                                    Liwpmodis,Lclmodis,Latb532,Latb532gr,                & !GLID
                                    LlidarBetaMol532,LlidarBetaMol532gr,LcfadLidarsr532, & !GLID
                                    LcfadLidarsr532gr,Lclcalipso2,                       & !GLID
                                    Lclcalipso,Lclgroundlidar,Lclhcalipso,Lcllcalipso,   & !GLID
                                    Lclmcalipso,Lcltcalipso,Lclhgroundlidar,             & !GLID
                                    Lcllgroundlidar,Lclmgroundlidar,Lcltgroundlidar,     & !GLID
                                    Lcltlidarradar,Lclcalipsoliq,                        & !GLID
                                    Lclcalipsoice,Lclcalipsoun,Lclcalipsotmp,            &
                                    Lclcalipsotmpliq,Lclcalipsotmpice,Lclcalipsotmpun,   &
                                    Lcltcalipsoliq,Lcltcalipsoice,Lcltcalipsoun,         &
                                    Lclhcalipsoliq,Lclhcalipsoice,Lclhcalipsoun,         &
                                    Lclmcalipsoliq,Lclmcalipsoice,Lclmcalipsoun,         &
                                    Lcllcalipsoliq,Lcllcalipsoice,Lcllcalipsoun,         & 
                                    Lclopaquecalipso,Lclthincalipso,Lclzopaquecalipso,   & !OPAQ
                                    Lclcalipsoopaque,Lclcalipsothin,Lclcalipsozopaque,   & !OPAQ
                                    Lclcalipsoopacity,Lclopaquetemp,Lclthintemp,         & !TIBO
                                    Lclzopaquetemp,Lclopaquemeanz,Lclthinmeanz,          & !TIBO
                                    Lclthinemis,Lclopaquemeanzse,Lclthinmeanzse,         & !TIBO2
                                    Lclzopaquecalipsose,LcfadDbze94,Ldbze94,Lparasolrefl,& !TIBO2
                                    Ltbrttov,Npoints,Ncolumns,Nlevels,Nlvgrid,Nchan,x)
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
         Latb532gr,        & ! GROUND LIDAR attenuated total backscatter (532nm) !GLID 
         LlidarBetaMol532, & ! CALIPSO molecular backscatter (532nm)         
         LlidarBetaMol532gr,&! GROUND LIDAR molecular backscatter (532nm)      !GLID    
         LcfadLidarsr532,  & ! CALIPSO scattering ratio CFAD
         LcfadLidarsr532gr,& ! GROUND LIDAR scattering ratio CFAD              !GLID
         Lclcalipso2,      & ! CALIPSO cloud fraction undetected by cloudsat
         Lclcalipso,       & ! CALIPSO cloud area fraction
         Lclgroundlidar,   & ! GROUND LIDAR cloud area fraction                !GLID
         Lclhcalipso,      & ! CALIPSO high-level cloud fraction
         Lcllcalipso,      & ! CALIPSO low-level cloud fraction
         Lclmcalipso,      & ! CALIPSO mid-level cloud fraction
         Lcltcalipso,      & ! CALIPSO total cloud fraction
         Lclhgroundlidar,  & ! GROUND LIDAR high-level cloud fraction          !GLID
         Lcllgroundlidar,  & ! GROUND LIDAR low-level cloud fraction           !GLID
         Lclmgroundlidar,  & ! GROUND LIDAR mid-level cloud fraction           !GLID
         Lcltgroundlidar,  & ! GROUND LIDAR total cloud fraction               !GLID
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
         Lclopaquecalipso, & ! CALIPSO opaque cloud cover (2D Map)                   !OPAQ
         Lclthincalipso,   & ! CALIPSO thin cloud cover (2D Map)                     !OPAQ
         Lclzopaquecalipso,& ! CALIPSO z_opaque altitude (opaque clouds only, 2D Map)!OPAQ
         Lclcalipsoopaque, & ! CALIPSO opaque cloud profiles 3D fraction             !OPAQ
         Lclcalipsothin,   & ! CALIPSO thin cloud profiles 3D fraction               !OPAQ
         Lclcalipsozopaque,& ! CALIPSO z_opaque 3D fraction                          !OPAQ
         Lclcalipsoopacity,& ! CALIPSO opacity 3D fraction                           !OPAQ
         Lclopaquetemp,    & ! CALIPSO opaque cloud temperature                   !TIBO
         Lclthintemp,      & ! CALIPSO thin cloud temperature                     !TIBO
         Lclzopaquetemp,   & ! CALIPSO z_opaque temperature                       !TIBO
         Lclopaquemeanz,   & ! CALIPSO opaque cloud altitude                      !TIBO
         Lclthinmeanz,     & ! CALIPSO thin cloud altitude                        !TIBO
         Lclthinemis,      & ! CALIPSO thin cloud emissivity                      !TIBO
         Lclopaquemeanzse,   & ! CALIPSO opaque cloud altitude with respect to SE !TIBO2
         Lclthinmeanzse,     & ! CALIPSO thin cloud altitude with respect to SE   !TIBO2
         Lclzopaquecalipsose,& ! CALIPSO z_opaque altitude with respect to SE     !TIBO2
         LcfadDbze94,      & ! CLOUDSAT radar reflectivity CFAD
         Ldbze94,          & ! CLOUDSAT radar reflectivity
         LparasolRefl,     & ! PARASOL reflectance
         Ltbrttov            ! RTTOV mean clear-sky brightness temperature
     
     integer,intent(in) :: &
          Npoints,         & ! Number of sampled points
          Ncolumns,        & ! Number of subgrid columns
          Nlevels,         & ! Number of model levels
          Nlvgrid,         & ! Number of levels in L3 stats computation
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
    if (Lclopaquecalipso .or. Lclthincalipso .or. Lclzopaquecalipso) then !OPAQ
        allocate(x%calipso_cldtype(Npoints,LIDAR_NTYPE))                  !OPAQ
    endif                                                                 !OPAQ
    if (Lclopaquetemp .or. Lclthintemp .or. Lclzopaquetemp) then          !TIBO
        allocate(x%calipso_cldtypetemp(Npoints,LIDAR_NTYPE))              !TIBO
    endif                                                                 !TIBO
    if (Lclopaquemeanz .or. Lclthinmeanz) then                            !TIBO
        allocate(x%calipso_cldtypemeanz(Npoints,2))                       !TIBO
    endif                                                                 !TIBO
    if (Lclopaquemeanzse .or. Lclthinmeanzse .or. Lclzopaquecalipsose) then !TIBO2
        allocate(x%calipso_cldtypemeanzse(Npoints,3))                       !TIBO2
    endif                                                                   !TIBO2
    if (Lclthinemis) then                                                 !TIBO
        allocate(x%calipso_cldthinemis(Npoints))                          !TIBO
    endif                                                                 !TIBO
    if (Lclcalipsoopaque .or. Lclcalipsothin .or. Lclcalipsozopaque .or. Lclcalipsoopacity) then !OPAQ
        allocate(x%calipso_lidarcldtype(Npoints,Nlvgrid,LIDAR_NTYPE+1))   !OPAQ
    endif                                                                 !OPAQ
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

    ! GROUND LIDAR simulator
    if (LlidarBetaMol532gr) allocate(x%groundlidar_beta_mol(Npoints,Nlevels))          !GLID
    if (Latb532gr)          allocate(x%groundlidar_beta_tot(Npoints,Ncolumns,Nlevels)) !GLID
    if (LcfadLidarsr532gr) then                                                        !GLID
        allocate(x%groundlidar_srbval(SR_BINS+1))                                      !GLID
        allocate(x%groundlidar_cfad_sr(Npoints,SR_BINS,Nlvgrid))                       !GLID
    endif                                                                              !GLID
    if (Lclgroundlidar)     allocate(x%groundlidar_lidarcld(Npoints,Nlvgrid))          !GLID
    if (Lclhgroundlidar .or. Lclmgroundlidar .or. Lcllgroundlidar .or. Lcltgroundlidar) then  !GLID
        allocate(x%groundlidar_cldlayer(Npoints,LIDAR_NCAT))                           !GLID
    endif                                                                              !GLID
      
    ! PARASOL
    if (Lparasolrefl) then
        allocate(x%parasolPix_refl(Npoints,Ncolumns,PARASOL_NREFL))
        allocate(x%parasolGrid_refl(Npoints,PARASOL_NREFL))
    endif 

    ! Cloudsat simulator
    if (Ldbze94)        allocate(x%cloudsat_Ze_tot(Npoints,Ncolumns,Nlevels))
    if (LcfadDbze94)    allocate(x%cloudsat_cfad_ze(Npoints,DBZE_BINS,Nlvgrid))

    ! Combined CALIPSO/CLOUDSAT fields
    if (Lclcalipso2)    allocate(x%lidar_only_freq_cloud(Npoints,Nlvgrid))
    if (Lcltlidarradar) allocate(x%radar_lidar_tcc(Npoints))
        
    ! RTTOV
    if (Ltbrttov) allocate(x%rttov_tbs(Npoints,Nchan))
 
  end subroutine construct_cosp_outputs
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE destroy_cospIN     
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine destroy_cospIN(y)
    type(cosp_optical_inputs),intent(inout) :: y

    if (allocated(y%tau_067))         deallocate(y%tau_067)
    if (allocated(y%emiss_11))        deallocate(y%emiss_11)
    if (allocated(y%frac_out))        deallocate(y%frac_out)
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
    if (allocated(y%beta_mol_gr))     deallocate(y%beta_mol_gr) !GLID
    if (allocated(y%betatot_gr))      deallocate(y%betatot_gr)  !GLID
    if (allocated(y%tau_mol_gr))      deallocate(y%tau_mol_gr)  !GLID
    if (allocated(y%tautot_gr))       deallocate(y%tautot_gr)   !GLID

  end subroutine destroy_cospIN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE destroy_cospstateIN     
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  subroutine destroy_cospstateIN(y)
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
    if (allocated(y%surfelev))        deallocate(y%surfelev)        !TIBO2
    
  end subroutine destroy_cospstateIN
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE destroy_cosp_outputs
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  subroutine destroy_cosp_outputs(y)
     type(cosp_outputs),intent(inout) :: y

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
     if (associated(y%calipso_lidarcldtype))     then !OPAQ
        deallocate(y%calipso_lidarcldtype)            !OPAQ
        nullify(y%calipso_lidarcldtype)               !OPAQ
     endif                                            !OPAQ
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
     if (associated(y%calipso_cldtype))          then !OPAQ
        deallocate(y%calipso_cldtype)                 !OPAQ
        nullify(y%calipso_cldtype)                    !OPAQ
     endif                                            !OPAQ
     if (associated(y%calipso_cldtypetemp))      then !TIBO
        deallocate(y%calipso_cldtypetemp)             !TIBO
        nullify(y%calipso_cldtypetemp)                !TIBO
     endif                                            !TIBO
     if (associated(y%calipso_cldtypemeanz))     then !TIBO
        deallocate(y%calipso_cldtypemeanz)            !TIBO
        nullify(y%calipso_cldtypemeanz)               !TIBO
     endif                                            !TIBO
     if (associated(y%calipso_cldtypemeanzse))   then !TIBO2
        deallocate(y%calipso_cldtypemeanzse)          !TIBO2
        nullify(y%calipso_cldtypemeanzse)             !TIBO2
     endif                                            !TIBO2
     if (associated(y%calipso_cldthinemis))      then !TIBO
        deallocate(y%calipso_cldthinemis)             !TIBO
        nullify(y%calipso_cldthinemis)                !TIBO
     endif                                            !TIBO
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
     if (associated(y%groundlidar_beta_mol))     then !GLID
        deallocate(y%groundlidar_beta_mol)            !GLID
        nullify(y%groundlidar_beta_mol)               !GLID
     endif                                            !GLID
     if (associated(y%groundlidar_beta_tot))     then !GLID
        deallocate(y%groundlidar_beta_tot)            !GLID
        nullify(y%groundlidar_beta_tot)               !GLID
     endif                                            !GLID
     if (associated(y%groundlidar_cldlayer))     then !GLID
        deallocate(y%groundlidar_cldlayer)            !GLID
        nullify(y%groundlidar_cldlayer)               !GLID
     endif                                            !GLID
     if (associated(y%groundlidar_lidarcld))     then !GLID
        deallocate(y%groundlidar_lidarcld)            !GLID
        nullify(y%groundlidar_lidarcld)               !GLID
     endif                                            !GLID
     if (associated(y%groundlidar_cfad_sr))      then !GLID
        deallocate(y%groundlidar_cfad_sr)             !GLID
        nullify(y%groundlidar_cfad_sr)                !GLID
     endif                                            !GLID
     if (associated(y%groundlidar_srbval))       then !GLID
        deallocate(y%groundlidar_srbval)              !GLID
        nullify(y%groundlidar_srbval)                 !GLID
     endif                                            !GLID
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
        
   end subroutine destroy_cosp_outputs
  
 end program cosp2_test

