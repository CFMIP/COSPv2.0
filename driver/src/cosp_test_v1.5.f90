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
! History
! May 2015 - D. Swales - Original version
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PROGRAM COSPTEST_v1p5
  USE MOD_COSP_INTERFACE_v1p5,    ONLY: cosp_interface => cosp_interface_v1p5,           &
                                        cosp_interface_init,                             &
                                        cosp_config,cosp_gridbox,cosp_subgrid,& 
                                        construct_cosp_gridbox,  destroy_cosp_gridbox,   &
                                        construct_cosp_subgrid,  destroy_cosp_subgrid,   &
                                        construct_cosp_config,I_CVCLIQ,I_LSCLIQ,         &
                                        I_CVCICE,I_LSCICE,I_LSRAIN,I_LSSNOW,I_LSGRPL,    &
                                        I_CVRAIN,I_CVSNOW
  USE MOD_COSP_IO,                ONLY: nc_read_input_file,nc_cmor_init,                 &
                                        nc_cmor_associate_1d,nc_cmor_write_1d,           &
                                        nc_cmor_associate_2d,nc_cmor_write_2d,           &
                                        nc_cmor_close,nc_cmor_write_1d_v1p5,var1d,var2d,var3d  
  USE COSP_KINDS,                 ONLY: wp,dp
  USE MOD_COSP,                   ONLY: linitialization,cosp_column_inputs,cosp_outputs, &
                                        construct_cosp_outputs,destroy_cosp_outputs
  USE MOD_COSP_CONFIG,            ONLY: RTTOV_MAX_CHANNELS,N_HYDRO,PARASOL_NREFL
  USE COSP_PHYS_CONSTANTS,        ONLY: amw,amd,amO3,amCO2,amCH4,amN2O,amCO
  
  IMPLICIT NONE
  
  ! Parameters
  character(len=64),parameter :: &
       cosp_input_namelist  = 'cosp_input_nl_1D.v1p5.txt', &
       cosp_output_namelist = 'cosp_output_nl.txt'
  character(len=32),parameter :: &
       cospvID = 'COSP v1.5'        ! COSP version ID
  integer,parameter :: &
       N_MAX_INPUT_FILES = 10000, &
       N_OUT_LIST = 63,           & ! Number of possible output variables
       N3D        = 8,            & ! Number of 3D output variables
       N2D        = 14,           & ! Number of 2D output variables
       N1D        = 40              ! Number of 1D output variables  

  ! Local variables
  type(cosp_gridbox)       :: gbx        ! Gridbox information. Input for COSP
  type(cosp_subgrid)       :: sgx        ! Subgrid outputs
  type(cosp_config)        :: cfg        ! Configuration options
  type(cosp_column_inputs) :: cospgridIN ! Model state needed by cospv1.5
  type(cosp_outputs)       :: cospOUT   ! COSP simulator outputs (flat)
  
  ! Sample input data variables
  integer :: &
       Nlon,Nlat,geomode,k
  real(wp) :: &
       emsfc_lw
  real(wp),dimension(:),    allocatable        :: &
       lon,lat,skt,landmask,u_wind,v_wind,sunlit
  real(wp),dimension(:,:),  allocatable,target :: &
       p,ph,zlev,zlev_half,T,sh,rh,tca,cca,mr_lsliq,mr_lsice,mr_ccliq,mr_ccice,fl_lsrain,&
       fl_lssnow,fl_lsgrpl,fl_ccrain,fl_ccsnow,dtau_s,dtau_c,dem_s,dem_c,mr_ozone
  real(wp),dimension(:,:,:),allocatable :: &
       frac_out,frac_prec,Reff
  character(len=600) :: &
       dfinput ! Input file
  
  ! Variables for hydrometeor description
  double precision :: time,time_bnds(2),time_step,half_time_step
  
  ! Output stuff
  integer :: N1,lon_axid,time_axid,height_axid,height_mlev_axid,grid_id,lonvar_id,       &
             latvar_id,column_axid,sza_axid,temp_axid,channel_axid,dbze_axid,sratio_axid,&
             MISR_CTH_axid,lat_axid,tau_axid,pressure2_axid 
  type(var1d) :: v1d(N1D+1) ! Structures needed by output routines for 1D variables
  type(var2d) :: v2d(N2D)   ! Structures needed by output routines for 2D variables
  type(var3d) :: v3d(N3D)   ! Structures needed by output routines for 3D variables
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Input Namelist fields
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  integer ::                      & !
       Npoints,                   & ! Number of gridpoints
       Ncolumns,                  & ! Number of subcolumns
       Nlevels,                   & ! Number of model vertical levels
       Npoints_it,                & ! Number of gridpoints to be processed in one 
                                    ! iteration
       Nlvgrid,                   & ! Number of vertical levels for statistical outputs 
                                    ! (USE_VGRID=.true.)
       surface_radar,             & ! surface=1/spaceborne=0
       use_mie_tables,            & ! Use a precomputed lookup-table (1=yes/0=no)
       use_gas_abs,               & ! Include gaseous absorption (1=yes/0=no)
       do_ray,                    & ! Calculate output Rayleigh (1=yes/0=no)
       melt_lay,                  & ! Melting layer model (1=on/0=off)
       Nprmts_max_hydro,          & ! Max number of parameters for hydrometeor size 
                                    ! distributions
       Naero,                     & ! Number of aerosol species (Not used)
       Nprmts_max_aero,           & ! Max number of parameters for aerosol size 
                                    ! distributions (Not used)
       lidar_ice_type,            & ! Ice particle shape in lidar calculations 
                                    ! (0=ice-spheres/1=ice-non-spherical)
       overlap,                   & ! Overlap type: 1=max, 2=rand, 3=max/rand
       isccp_topheight,           & ! ISCCP cloud top height
       isccp_topheight_direction, & ! ISCCP cloud top height direction
       platform,                  & ! RTTOV: Satellite platform
       satellite,                 & ! RTTOV: Satellite
       instrument,                & ! RTTOV: Instrument
       Nchannels ,                & ! RTTOV: Number of channels to be computed
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
       radar_freq,                & ! CloudSat radar frequency (GHz)
       k2,                        & ! |K|^2, -1=use frequency dependent default
       ZenAng,                    & ! RTTOV: Satellite Zenith Angle
       co2,                       & ! CO2 mixing ratio
       ch4,                       & ! CH4 mixing ratio
       n2o,                       & ! n2o mixing ratio
       co                           ! co mixing ratio
  logical ::                      & !
       use_vgrid,                 & ! Use fixed vertical grid for outputs?
       csat_vgrid,                & ! CloudSat vertical grid? 
       use_precipitation_fluxes,  & ! True if precipitation fluxes are input to the 
                                    ! algorithm 
       use_reff,                  & ! True if you want effective radius to be used by 
                                    ! radar simulator (always used by lidar)
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
       claraRTTOV_calcrefl              !
  integer,dimension(RTTOV_MAX_CHANNELS) :: &
       Channels                     ! RTTOV: Channel numbers
  real(wp),dimension(RTTOV_MAX_CHANNELS) :: &
       Surfem                       ! RTTOV: Surface emissivity
  character(len=64) :: &
       radar_micro_scheme,        & ! Microphysical scheme used in radar simulator
       claraRTTOV_coefdir           ! Location of RTTOV coefficient files
  character(len=64),dimension(N_MAX_INPUT_FILES) :: &
       finput                       ! List input NetCDF files
  character(len=512) :: &
       dinput,                    & ! Directory where the input files are located
       cmor_nl                      ! CMOR namelist
  namelist/COSP_INPUT/cmor_nl,overlap,isccp_topheight,isccp_topheight_direction,         &
                      npoints,npoints_it,ncolumns,nlevels,use_vgrid,Nlvgrid,csat_vgrid,  &
                      dinput,finput,radar_freq,surface_radar,use_mie_tables,             &
                      use_gas_abs,do_ray,melt_lay,k2,radar_micro_scheme,                 &
                      Nprmts_max_hydro,Naero,Nprmts_max_aero,lidar_ice_type,             &
                      use_precipitation_fluxes,use_reff,platform,satellite,              &
                      Instrument,Nchannels,Channels,Surfem,ZenAng,co2,ch4,n2o,co,        &
                      CLARA_Tb_subvis,CLARA_Tb_semitrans,CLARA_Tb_opaque,CLARA_RTTOVclr, &
                      claraRTTOV_coefdir,claraRTTOV_platform,claraRTTOV_satellite,       &
                      claraRTTOV_sensor,claraRTTOV_nchan,claraRTTOV_channels,            &
                      claraRTTOV_addrefrac,claraRTTOV_use_q2m,claraRTTOV_clw_data,       &
                      claraRTTOV_addsolar,claraRTTOV_addclouds,claraRTTOV_addaerosol,    &
                      claraRTTOV_use_cld_opts_param,claraRTTOV_ozone_data,claraRTTOV_co2,&
                      claraRTTOV_n2o,claraRTTOV_ch4,claraRTTOV_co,claraRTTOV_addinterp,  &
                      claraRTTOV_calcemis,claraRTTOV_calcrefl,CLARA_retSize
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Read in namelist
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Input namelist
  finput(:) = ''
  open(10,file=trim(cosp_input_namelist),status='old')
  read(10,nml=COSP_INPUT)
  close(10)

  ! Output namelist contains only logical switches that go directly into the derived type
  ! for the configuration settings, cosp_config 
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Read in sample input data
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Allocate space
  allocate(lon(Npoints),lat(Npoints),p(Npoints,Nlevels),ph(Npoints,Nlevels),             &
           zlev(Npoints,Nlevels),zlev_half(Npoints,Nlevels+1),T(Npoints,Nlevels),        &
           sh(Npoints,Nlevels),rh(Npoints,Nlevels),tca(Npoints,Nlevels),                 &
           cca(Npoints,Nlevels),mr_lsliq(Npoints,Nlevels),mr_lsice(Npoints,Nlevels),     &
           mr_ccliq(Npoints,Nlevels),mr_ccice(Npoints,Nlevels),                          &
           fl_lsrain(Npoints,Nlevels),fl_lssnow(Npoints,Nlevels),                        &
           fl_lsgrpl(Npoints,Nlevels),fl_ccrain(Npoints,Nlevels),                        &
           fl_ccsnow(Npoints,Nlevels),Reff(Npoints,Nlevels,N_HYDRO),                     &
           dtau_s(Npoints,Nlevels),dtau_c(Npoints,Nlevels),dem_s(Npoints,Nlevels),       &
           dem_c(Npoints,Nlevels),skt(Npoints),landmask(Npoints),                        &
           mr_ozone(Npoints,Nlevels),u_wind(Npoints),v_wind(Npoints),sunlit(Npoints),    &
           frac_out(Npoints,Ncolumns,Nlevels),                                           &
           frac_prec(Npoints,Ncolumns,Nlevels))
  
  ! Read in data
  dfinput = trim(dinput)//trim(finput(1))
  call nc_read_input_file(dfinput,Npoints,Nlevels,N_HYDRO,lon,lat,p,ph,zlev,zlev_half,   &
                          T,sh,rh,tca,cca,mr_lsliq,mr_lsice,mr_ccliq,mr_ccice,fl_lsrain, &
                          fl_lssnow,fl_lsgrpl,fl_ccrain,fl_ccsnow,Reff,dtau_s,dtau_c,    &
                          dem_s,dem_c,skt,landmask,mr_ozone,u_wind,v_wind,sunlit,        &
                          emsfc_lw,geomode,Nlon,Nlat)                      
  
  ! Time information
  time           = 8*1._wp/8._wp  ! First time step
  time_step      = 3._wp/24._wp!3.D0/24.D0
  half_time_step = 0.5_wp*time_step
  time_bnds      = (/time-half_time_step,time+half_time_step/) 
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Construct input derived types to use an input into COSP.
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  ! Obtain configuration options from namelist and populate configuration derived type
  call construct_cosp_config(cosp_output_namelist,N_OUT_LIST,cfg)

  ! Gridbox
  call construct_cosp_gridbox(time,time_bnds,Npoints,Nlevels,Ncolumns,                   &
                              N_HYDRO,Nprmts_max_hydro,Naero,Nprmts_max_aero,Npoints_it, &
                              emsfc_lw,Nchannels,ZenAng,surfem(1:Nchannels),co2,ch4,n2o, &
                              co,lon,lat,p,ph,zlev,zlev_half,T,sh,cca,tca,skt,landmask,  &
                              mr_ozone,u_wind,v_wind,sunlit,fl_lsrain,fl_lssnow,         &
                              fl_lsgrpl,fl_ccrain,fl_ccsnow,dtau_s,dtau_c,dem_s,dem_c,   &
                              Reff,mr_lsliq,mr_lsice,mr_ccliq,mr_ccice,gbx)

  ! Subgrid
  call construct_cosp_subgrid(Npoints, Ncolumns, Nlevels, sgx)

  ! Call initialization
  linitialization=.TRUE.
  if (linitialization) then
    call cosp_interface_init(Npoints,Nlevels,Npoints_it,overlap,use_precipitation_fluxes,&
                             radar_freq,radar_micro_scheme,k2,use_gas_abs,do_ray,        &
                             isccp_topheight,isccp_topheight_direction,                  &
                             zlev(:,Nlevels:1:-1),zlev_half(:,Nlevels:1:-1),             &
                             surface_radar,Nchannels,Channels,platform,satellite,        &
                             instrument,lidar_ice_type,use_vgrid,Nlvgrid,csat_vgrid,     &
                             cospvID,CLARA_Tb_subvis,CLARA_Tb_semitrans,CLARA_Tb_opaque, &
                             CLARA_RTTOVclr,claraRTTOV_coefdir,claraRTTOV_platform,      &
                             claraRTTOV_satellite,claraRTTOV_sensor,claraRTTOV_nchan,    &
                             claraRTTOV_channels,claraRTTOV_addrefrac,claraRTTOV_use_q2m,&
                             claraRTTOV_clw_data,claraRTTOV_addsolar,                    &
                             claraRTTOV_addclouds,claraRTTOV_addaerosol,                 &
                             claraRTTOV_use_cld_opts_param,claraRTTOV_ozone_data,        &
                             claraRTTOV_co2,claraRTTOV_n2o,claraRTTOV_ch4,claraRTTOV_co, &
                             claraRTTOV_addinterp,claraRTTOV_calcemis,                   &
                             claraRTTOV_calcrefl,CLARA_retSize)
  endif                        
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Construct output derived types.
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  call construct_cosp_outputs(cfg%Lpctisccp,cfg%Lclisccp,cfg%Lboxptopisccp,              &
                              cfg%Lboxtauisccp,cfg%Ltauisccp,cfg%Lcltisccp,              &
                              cfg%Lmeantbisccp,cfg%Lmeantbclrisccp,cfg%Lalbisccp,        &
                              cfg%LclMISR,cfg%Lcltmodis,cfg%Lclwmodis,cfg%Lclimodis,     &
                              cfg%Lclhmodis,cfg%Lclmmodis,cfg%Lcllmodis,cfg%Ltautmodis,  &
                              cfg%Ltauwmodis,cfg%Ltauimodis,cfg%Ltautlogmodis,           &
                              cfg%Ltauwlogmodis,cfg%Ltauilogmodis,cfg%Lreffclwmodis,     &
                              cfg%Lreffclimodis,cfg%Lpctmodis,cfg%Llwpmodis,             &
                              cfg%Liwpmodis,cfg%Lclmodis,cfg%Latb532,                    &
                              cfg%LlidarBetaMol532,cfg%LcfadLidarsr532,cfg%Lclcalipso2,  &
                              cfg%Lclcalipso,cfg%Lclhcalipso,cfg%Lcllcalipso,            &
                              cfg%Lclmcalipso,cfg%Lcltcalipso,cfg%Lcltlidarradar,        &
                              cfg%Lclcalipsoliq,cfg%Lclcalipsoice,cfg%Lclcalipsoun,      &
                              cfg%Lclcalipsotmp,cfg%Lclcalipsotmpliq,                    &
                              cfg%Lclcalipsotmpice,cfg%Lclcalipsotmpun,                  &
                              cfg%Lcltcalipsoliq,cfg%Lcltcalipsoice,cfg%Lcltcalipsoun,   &
                              cfg%Lclhcalipsoliq,cfg%Lclhcalipsoice,cfg%Lclhcalipsoun,   &
                              cfg%Lclmcalipsoliq,cfg%Lclmcalipsoice,cfg%Lclmcalipsoun,   &
                              cfg%Lcllcalipsoliq,cfg%Lcllcalipsoice,cfg%Lcllcalipsoun,   &
                              cfg%LcfadDbze94,cfg%Ldbze94,cfg%Lparasolrefl,cfg%Ltbrttov, &
                              Npoints,Ncolumns,Nlevels,Nlvgrid,Nchannels,cospOUT)
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Call COSP
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! Call main cosp engine
  call cosp_interface(Npoints,gbx,sgx,cospOUT)
 
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Write outputs to CMOR-compliant netCDF format.
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (cfg%Lwrite_output) then
     N1 = N1D
     if (geomode == 1) N1 = N1D+1
     call nc_cmor_init(cmor_nl,'replace',cfg,gbx,sgx,cospOUT,                     &
                       geomode,Nlon,Nlat,N1,N2D,N3D,N_OUT_LIST,lon_axid,lat_axid,        &
                       time_axid,height_axid,                                            &
                       height_mlev_axid,grid_id,lonvar_id,latvar_id, column_axid,        &
                       sza_axid,temp_axid,channel_axid,dbze_axid,sratio_axid,            &
                       MISR_CTH_axid,tau_axid,pressure2_axid,v1d(1:N1),v2d,v3d)
     if (geomode == 1) then
        call nc_cmor_associate_1d(grid_id,time_axid,height_axid,height_mlev_axid,        &
                                  column_axid,sza_axid,temp_axid,channel_axid,dbze_axid, &
                                  sratio_axid,MISR_CTH_axid,tau_axid,pressure2_axid,Nlon,&
                                  Nlat,gbx,sgx,                                    &
                                  cospOUT,N1D,N2D,N3D,v1d(1:N1),v2d,v3d)
        call nc_cmor_write_1d_v1p5(gbx,time_bnds,lonvar_id,latvar_id,N1,N2D,N3D,         &
                                   v1d(1:N1),v2d,v3d)
     elseif (geomode >  1) then
        call nc_cmor_associate_2d(lon_axid,lat_axid,time_axid,height_axid,               &
                                  height_mlev_axid,column_axid,sza_axid,temp_axid,       &
                                  channel_axid,dbze_axid,sratio_axid,MISR_CTH_axid,      &
                                  tau_axid,pressure2_axid,Nlon,Nlat,gbx,sgx,       &
                                  cospOUT,N1D,N2D,N3D,v1d(1:N1),v2d,v3d)
        call nc_cmor_write_2d(time,time_bnds,geomode,Nlon,Nlat,N1,N2D,N3D,v1d(1:N1),v2d,&
                              v3d)
     endif
     call nc_cmor_close()
  endif
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Deallocate local arrays from memory
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  deallocate(lon,lat,p,ph,zlev,zlev_half,T,sh,rh,tca,cca, mr_lsliq,mr_lsice,mr_ccliq,    &
             mr_ccice,fl_lsrain,fl_lssnow,fl_lsgrpl,fl_ccrain,fl_ccsnow,Reff,dtau_s,     &
             dtau_c,dem_s,dem_c,skt,landmask,mr_ozone,u_wind,v_wind,sunlit,frac_out)
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Free up memory
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			 
  call destroy_cosp_gridbox(gbx)
  call destroy_cosp_subgrid(sgx)
  call destroy_cosp_outputs(cospOUT)
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! END PROGRAM
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END PROGRAM COSPTEST_v1p5
