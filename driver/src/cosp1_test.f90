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
program cosp1_test
  use cosp_kinds,          only: wp                         
  USE MOD_COSP_CONFIG,     ONLY: R_UNDEF,PARASOL_NREFL,LIDAR_NCAT,SR_BINS,                &
                                 N_HYDRO,RTTOV_MAX_CHANNELS,numMISRHgtBins,               &
                                 cloudsat_DBZE_BINS,LIDAR_NTEMP,calipso_histBsct,                  &
                                 numMODISTauBins,numMODISPresBins,                        &
                                 numMODISReffIceBins,numMODISReffLiqBins,                 &
                                 numISCCPTauBins,numISCCPPresBins,numMISRTauBins,         &
                                 ntau,modis_histTau,tau_binBounds,                        &
                                 modis_histTauEdges,tau_binEdges,                         &
                                 modis_histTauCenters,tau_binCenters,ntauV1p4,            &
                                 tau_binBoundsV1p4,tau_binEdgesV1p4, tau_binCentersV1p4,  &
                                 Nlvgrid_local  => Nlvgrid
  use mod_cosp1_io,         only: nc_read_input_file, write_cosp1_output, read_cosp_output_nl

  USE MOD_COSP_INTERFACE_v1p4, ONLY: cosp                   => cosp_interface_v1p4,      &
                                     cosp_gridbox,construct_cosp_vgrid,                  &
                                     construct_cosp_gridbox,                             &
                                     free_cosp_gridbox      => destroy_cosp_gridbox,     &
                                     free_cosp_sgradar      => destroy_cosp_sgradar,     &
                                     free_cosp_radarstats   => destroy_cosp_radarstats,  &
                                     free_cosp_sglidar      => destroy_cosp_sglidar,     &
                                     free_cosp_lidarstats   => destroy_cosp_lidarstats,  &
                                     free_cosp_isccp        => destroy_cosp_isccp,       &
                                     free_cosp_misr         => destroy_cosp_misr,        &
                                     free_cosp_rttov        => destroy_cosp_rttov,       &
                                     free_cosp_modis        => destroy_cosp_modis,       &
                                     free_cosp_vgrid        => destroy_cosp_vgrid,       &
                                     free_cosp_subgrid      => destroy_cosp_subgrid,     &
                                     construct_cosp_subgrid,cosp_config,cosp_subgrid,    &
                                     cosp_sglidar,cosp_lidarstats,                       &
                                     construct_cosp_lidarstats,construct_cosp_sglidar,   &
                                     cosp_isccp,construct_cosp_isccp,cosp_misr,          &
                                     construct_cosp_misr,cosp_rttov,construct_cosp_rttov,&
                                     cosp_sgradar,cosp_radarstats,                       &
                                     construct_cosp_radarstats,construct_cosp_sgradar,   &
                                     cosp_modis,construct_cosp_modis,                    &
                                     cosp_vgrid,I_CVCLIQ,I_LSCLIQ,I_CVCICE,I_LSCICE,     &
                                     I_LSRAIN,I_LSSNOW,I_LSGRPL,I_CVRAIN,I_CVSNOW
  implicit none

  ! Input/Output driver file control
  character(len=64),parameter :: &
       cosp_input_namelist  = 'cosp1_input_nl.txt', &
       cosp_output_namelist = 'cosp1_output_nl.txt'

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
  namelist/COSP_INPUT/overlap,isccp_topheight,isccp_topheight_direction,npoints,         &
                      npoints_it,ncolumns,nlevels,use_vgrid,Nlvgrid,csat_vgrid,dinput,   &
                      finput,foutput,cloudsat_radar_freq,surface_radar,cloudsat_use_gas_abs,     &
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

  integer :: iChunk,nChunks,start_idx,end_idx,nPtsPerIt
  real(wp),dimension(10) :: driver_time

  character(len=256),dimension(100) :: cosp_status

  ! Fields used solely for output
  integer,parameter :: &
       n_out_list = 63,           & ! Number of possible output variables
       N3D        = 8,            & ! Number of 3D output variables
       N2D        = 14,           & ! Number of 2D output variables
       N1D        = 40              ! Number of 1D output variables  
  character(len=32),dimension(n_out_list) :: out_list  ! List of output variable names
  integer :: lon_axid,time_axid,height_axid,height_mlev_axid,grid_id,lonvar_id,       &
             latvar_id,column_axid,sza_axid,temp_axid,channel_axid,dbze_axid,sratio_axid,&
             MISR_CTH_axid,lat_axid,tau_axid,pressure2_axid ,k
  double precision :: time,time_bnds(2),time_step,half_time_step
  real(wp),dimension(:),allocatable :: mgrid_z,mgrid_zu,mgrid_zl

  ! COSP1 types
  type(cosp_gridbox)    :: gbx     ! Gridbox information. Input for COSP
  type(cosp_subgrid)    :: sgx     ! Subgrid outputs
  type(cosp_config)     :: cfg     ! Configuration options
  type(cosp_vgrid)      :: vgrid   ! Information on vertical grid of stats
  type(cosp_sgradar)    :: sgradar ! Output from radar simulator
  type(cosp_sglidar)    :: sglidar ! Output from lidar simulator
  type(cosp_isccp)      :: isccp   ! Output from ISCCP simulator
  type(cosp_modis)      :: modis   ! Output from MODIS simulator
  type(cosp_misr)       :: misr    ! Output from MISR simulator
  type(cosp_rttov)      :: rttov   ! Output from RTTOV 
  type(cosp_radarstats) :: stradar ! Summary statistics from radar simulator
  type(cosp_lidarstats) :: stlidar ! Summary statistics from lidar simulator
  integer,parameter ::      &
       use_mie_tables=0,    &
       melt_lay=0,          & ! melting layer model off=0, on=1
       Nprmts_max_hydro=12, & ! Max number of parameters for hydrometeor size distributions
       Naero=1,             & ! Number of aerosol species (Not used)
       Nprmts_max_aero=1      ! Max number of parameters for aerosol size distributions (Not used)
  logical,parameter ::      &
       use_reff=.true.        ! True if you want effective radius to be used by radar simulator (always used by lidar)
  real(wp) :: toffset_step

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
  call read_cosp_output_nl(cosp_output_namelist,63,cfg)

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

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Allocate memory for gridbox type
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  call construct_cosp_gridbox(time, time_bnds, cloudsat_radar_freq,                      &
       surface_radar, use_mie_tables, cloudsat_use_gas_abs, cloudsat_do_ray, melt_lay,   &
       cloudsat_k2, Npoints, Nlevels, Ncolumns, N_HYDRO,              &
       Nprmts_max_hydro, Naero, Nprmts_max_aero, Npoints_it,lidar_ice_type,              &
       isccp_topheight, isccp_topheight_direction, overlap, emsfc_lw,                    &
       use_precipitation_fluxes, use_reff, rttov_platform, rttov_satellite,              &
       rttov_instrument, rttov_Nchannels, rttov_ZenAng,rttov_channels(1:rttov_Nchannels),&
       rttov_surfem(1:rttov_Nchannels), co2, ch4, n2o, co,gbx)

  gbx%longitude = lon
  gbx%latitude  = lat
  ! Toffset. This assumes that time is the mid-point of the interval.
  do k=1,Npoints
     gbx%toffset(k) = -half_time_step + toffset_step*(k-0.5)
  enddo
  gbx%p         = p
  gbx%ph        = ph
  gbx%zlev      = zlev
  gbx%zlev_half = zlev_half
  gbx%T         = T
  gbx%q         = rh
  gbx%sh        = sh
  gbx%cca       = cca
  gbx%tca       = tca
  gbx%psfc      = ph(:,1)
  gbx%skt       = skt
  gbx%land      = landmask
  gbx%mr_ozone  = mr_ozone
  gbx%u_wind    = u_wind
  gbx%v_wind    = v_wind
  gbx%sunlit    = sunlit
  gbx%rain_ls   = fl_lsrain
  gbx%snow_ls   = fl_lssnow
  gbx%grpl_ls   = fl_lsgrpl
  gbx%rain_cv   = fl_ccrain
  gbx%snow_cv   = fl_ccsnow
  gbx%dtau_s    = dtau_s
  gbx%dtau_c    = dtau_c
  gbx%dem_s     = dem_s
  gbx%dem_c     = dem_c
  gbx%Reff      = Reff
  gbx%Reff(:,:,I_LSRAIN) = 0._wp
  gbx%mr_hydro(:,:,I_LSCLIQ) = mr_lsliq
  gbx%mr_hydro(:,:,I_LSCICE) = mr_lsice
  gbx%mr_hydro(:,:,I_CVCLIQ) = mr_ccliq
  gbx%mr_hydro(:,:,I_CVCICE) = mr_ccice
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  ! Define new vertical grid
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  call construct_cosp_vgrid(gbx, Nlvgrid, use_vgrid, csat_vgrid, vgrid)
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  ! Subgrid information
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  call construct_cosp_subgrid(Npoints, Ncolumns, Nlevels, sgx)
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  ! Allocate memory for other types.
  ! *NOTE* These construct subroutines are different than the original construct 
  !        subroutines provided with cospv1.4.0. The old subroutines required the 
  !        derived type cosp_config to be provided, which was used to determine which
  !        simulators were to be run. For simulators that were not run, a minimum
  !        amount of space was allocated for that simulator.
  !        In COSPv2.0, which simulators are run is determined by looking at which
  !        output fields are allocated (i.e. if the output field for the modis tau vs.
  !        cloud-top height joint histogram is allocated, we know that the ISCCP and
  !        MODIS simulators need to be run). This change in v2.0 makes the way that
  !        the simulators outputs were allocated in compatable, so these subroutines
  !        needed to be modified, albeit only slightly.
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  if (cfg%Lradar_sim) call construct_cosp_sgradar(Npoints,Ncolumns,Nlevels,N_HYDRO,sgradar)
  if (cfg%Lradar_sim) call construct_cosp_radarstats(Npoints,Ncolumns,vgrid%Nlvgrid,N_HYDRO,stradar)  
  if (cfg%Llidar_sim) call construct_cosp_sglidar(Npoints,Ncolumns,Nlevels,N_HYDRO,PARASOL_NREFL,sglidar)
  if (cfg%Llidar_sim) call construct_cosp_lidarstats(Npoints,Ncolumns,vgrid%Nlvgrid,N_HYDRO,PARASOL_NREFL,stlidar)
  if (cfg%Lisccp_sim) call construct_cosp_isccp(Npoints,Ncolumns,Nlevels,isccp)
  if (cfg%Lmodis_sim) call construct_cosp_modis(Npoints,modis)
  if (cfg%Lmisr_sim)  call construct_cosp_misr(Npoints,misr)
  if (cfg%Lrttov_sim) call construct_cosp_rttov(Npoints,rttov_Nchannels,rttov)
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  ! Call simulator
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  call cpu_time(driver_time(6))
  call cosp(overlap,Ncolumns,cfg,vgrid,gbx,sgx,sgradar,sglidar,isccp,misr,modis,rttov,&
       stradar,stlidar)
  call cpu_time(driver_time(7))

  print*,'Time to read in data:     ',driver_time(2)-driver_time(1)
  print*,'Time to run COSP:         ',driver_time(7)-driver_time(6)
  print*,'Total time:               ',driver_time(7)-driver_time(1)

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Output
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  call write_cosp1_output(Npoints, Ncolumns, Nlevels, gbx%zlev(1,Nlevels:1:-1),lon, lat,  cfg, vgrid, gbx, sgx,    &
       sgradar, sglidar, isccp, misr, modis, rttov, stradar, stlidar, foutput)
  call cpu_time(driver_time(8))
  print*,'Time to write to output:  ',driver_time(8)-driver_time(7)

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Free up memory
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  call free_cosp_gridbox(gbx)
  call free_cosp_subgrid(sgx)
  call free_cosp_vgrid(vgrid)
  if (cfg%Lradar_sim) call free_cosp_sgradar(sgradar)
  if (cfg%Lradar_sim) call free_cosp_radarstats(stradar)
  if (cfg%Llidar_sim) call free_cosp_sglidar(sglidar)
  if (cfg%Llidar_sim) call free_cosp_lidarstats(stlidar)
  if (cfg%Lisccp_sim) call free_cosp_isccp(isccp)
  if (cfg%Lmisr_sim)  call free_cosp_misr(misr)
  if (cfg%Lmodis_sim) call free_cosp_modis(modis)
  if (cfg%Lrttov_sim) call free_cosp_rttov(rttov)
  
 end program cosp1_test

