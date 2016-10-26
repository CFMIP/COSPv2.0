! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copyright (c) 2008, British Crown Copyright, The Met Office
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
! Feb 2008 - A. Bodas-Salcedo - Initial version
! Dec 2010 - A. Bodas-Salcedo - Added capability for processing multiple files
! May 2015 - D. Swales        - Modified for COSPv2.0
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PROGRAM COSPTEST_trunk
  USE COSP_KINDS,              ONLY: wp,dp
  USE MOD_COSP_CONFIG,         ONLY: RTTOV_MAX_CHANNELS,N_HYDRO,PARASOL_NREFL
  USE MOD_COSP_IO,             ONLY: nc_read_input_file,nc_cmor_init,                    &
                                     nc_cmor_associate_1d,nc_cmor_write_1d,              &
                                     nc_cmor_associate_2d,nc_cmor_write_2d,              &
                                     nc_cmor_close,var1d,var2d,var3d,read_cosp_output_nl,&
                                     cosp_error
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
  IMPLICIT NONE
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Parameters
  character(len=64),parameter :: &
       cosp_input_namelist  = 'cosp_input_nl.v1.4.txt', &
       cosp_output_namelist = 'cosp_output_nl_v1.4.txt'
  integer,parameter :: &
       N_MAX_INPUT_FILES = 10000, &
       N_OUT_LIST = 65,           & ! Number of possible output variables
       N3D        = 10,           & ! Number of 3D output variables
       N2D        = 14,           & ! Number of 2D output variables
       N1D        = 40              ! Number of 1D output variables

  ! Local variables
  type(cosp_gridbox)    :: gbx     ! Gridbox information. Input for COSP
  type(cosp_subgrid)    :: sgx     ! Subgrid outputs
  type(cosp_config)     :: cfg     ! Configuration options
  type(cosp_vgrid)      :: vgrid   ! Information on vertical grid of stats
  type(cosp_sgradar)    :: sgradar ! Output from radar simulator
  type(cosp_sgradar)    :: armsgradar ! Output from radar simulator
  type(cosp_sglidar)    :: sglidar ! Output from lidar simulator
  type(cosp_isccp)      :: isccp   ! Output from ISCCP simulator
  type(cosp_modis)      :: modis   ! Output from MODIS simulator
  type(cosp_misr)       :: misr    ! Output from MISR simulator
  type(cosp_rttov)      :: rttov   ! Output from RTTOV 
  type(cosp_radarstats) :: stradar ! Summary statistics from radar simulator
  type(cosp_radarstats) :: armstradar ! Summary statistics from radar simulator
  type(cosp_lidarstats) :: stlidar ! Summary statistics from lidar simulator
  
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
       Reff
  character(len=600) :: &
       dfinput ! Input file
  real(wp),dimension(10) :: driver_time
  real(wp),dimension(:),allocatable :: mgrid_zl,mgrid_zu,mgrid_z

  ! Variables for hydrometeor description
  double precision :: time,time_bnds(2),time_step,half_time_step
  
  ! Output stuff
  integer :: N1,lon_axid,time_axid,height_axid,height_mlev_axid,grid_id,lonvar_id,       &
             latvar_id,column_axid,sza_axid,temp_axid,channel_axid,dbze_axid,sratio_axid,&
             MISR_CTH_axid,lat_axid,tau_axid,pressure2_axid,i,count_rate,count_max,t0,t1,&
             t2,t3,Nfiles
  real(wp) :: toffset_step
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
       arm_surface_radar,         & ! surface=1/spaceborne=0
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
       Nchannels                    ! RTTOV: Number of channels to be computed
  real(wp) ::                     & !
       radar_freq,                & ! CloudSat radar frequency (GHz)
       arm_radar_freq,            & ! arm radar frequency (GHz)
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
       use_reff                     ! True if you want effective radius to be used by 
                                    ! radar simulator (always used by lidar)
  integer,dimension(RTTOV_MAX_CHANNELS) :: &
       Channels                     ! RTTOV: Channel numbers
  real(wp),dimension(RTTOV_MAX_CHANNELS) :: &
       Surfem                       ! RTTOV: Surface emissivity
  character(len=64),dimension(N_MAX_INPUT_FILES) :: &
       finput                       ! List input NetCDF files
  character(len=512) :: &
       dinput,                    & ! Directory where the input files are located
       cmor_nl                      ! CMOR namelist
  
  namelist/COSP_INPUT/cmor_nl,overlap,isccp_topheight,isccp_topheight_direction, &
                      npoints,npoints_it,ncolumns,nlevels,use_vgrid,nlvgrid,     &
                      csat_vgrid,dinput,finput,radar_freq,arm_radar_freq,        &
                      surface_radar,arm_surface_radar,use_mie_tables,            &
                      use_gas_abs,do_ray,melt_lay,k2,                            &
                      Nprmts_max_hydro,Naero,Nprmts_max_aero,lidar_ice_type,     &
                      use_precipitation_fluxes,use_reff,platform,satellite,      &
                      Instrument,Nchannels,Channels,Surfem,ZenAng,co2,ch4,n2o,co
  ! Start timer
  call cpu_time(driver_time(1))
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Read COSP namelists
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  finput(:) = ''
  open(10,file=cosp_input_namelist,status='old')
  read(10,nml=cosp_input)
  close(10)
  call read_cosp_output_nl(cosp_output_namelist,N_OUT_LIST,cfg)
    
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Find number of input files
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  i=1
  do while (i <= N_MAX_INPUT_FILES)
     if (len_trim(finput(i)) < 1) exit
     i=i+1
  enddo
  Nfiles = i-1
  if (Nfiles < 1) call cosp_error('cosp_test','Number of files < 1')
  if (Nfiles > N_MAX_INPUT_FILES) then
     call cosp_error('cosp_test','Number of files > N_MAX_INPUT_FILES')
  endif
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Allocate local arrays
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
           mr_ozone(Npoints,Nlevels),u_wind(Npoints),v_wind(Npoints),sunlit(Npoints))

  ! Time information
  time           = 8*1._wp/8._wp  ! First time step
  time_step      = 3._wp/24._wp
  half_time_step = 0.5_wp*time_step
  time_bnds      = (/time-half_time_step,time+half_time_step/) 
  call cpu_time(driver_time(2))

  do i=1,Nfiles
     dfinput=trim(dinput)//trim(finput(1))
     time_bnds = (/time-half_time_step,time+half_time_step/) ! This may need to be adjusted, 
                                                             ! depending on the approx_interval in the MIP table
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! Read input geophysical variables from NetCDF file
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     call nc_read_input_file(dfinput,Npoints,Nlevels,N_HYDRO,lon,lat,p,ph,zlev,zlev_half,&
                             T,sh,rh,tca,cca,mr_lsliq,mr_lsice,mr_ccliq,mr_ccice,        &
                             fl_lsrain,fl_lssnow,fl_lsgrpl,fl_ccrain,fl_ccsnow,Reff,     &
                             dtau_s,dtau_c,dem_s,dem_c,skt,landmask,mr_ozone,u_wind,     &
                             v_wind,sunlit,emsfc_lw,geomode,Nlon,Nlat)
     ! geomode = 2 for (lon,lat) mode.
     ! geomode = 3 for (lat,lon) mode.
     ! In those modes it returns Nlon and Nlat with the correct values

     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! Allocate memory for gridbox type
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
     call construct_cosp_gridbox(time,time_bnds,radar_freq,arm_radar_freq,surface_radar, &
                                 arm_surface_radar,use_mie_tables,                       &
                                 use_gas_abs,do_ray,melt_lay,k2,Npoints,Nlevels,Ncolumns,&
                                 N_HYDRO,Nprmts_max_hydro,Naero,Nprmts_max_aero,         &
                                 Npoints_it,lidar_ice_type,isccp_topheight,              &
                                 isccp_topheight_direction,overlap,emsfc_lw,             &
                                 use_precipitation_fluxes,use_reff,Platform,Satellite,   &
                                 Instrument,Nchannels,ZenAng, channels(1:Nchannels),     &
                                 surfem(1:Nchannels),co2,ch4,n2o,co,gbx)
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
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
     ! Define new vertical grid
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
     call construct_cosp_vgrid(gbx,Nlvgrid,use_vgrid,csat_vgrid,vgrid)
  
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
     ! Subgrid information
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
     call construct_cosp_subgrid(Npoints, Ncolumns, Nlevels, sgx)

     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
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
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
     if (cfg%Lradar_sim) call construct_cosp_sgradar(Npoints,Ncolumns,Nlevels,N_HYDRO,sgradar)
     if (cfg%Lradar_sim) call construct_cosp_radarstats(Npoints,Ncolumns,vgrid%Nlvgrid,N_HYDRO,stradar)  
     if (cfg%Larmradar_sim) call construct_cosp_sgradar(Npoints,Ncolumns,Nlevels,N_HYDRO,armsgradar)
     if (cfg%Larmradar_sim) call construct_cosp_radarstats(Npoints,Ncolumns,vgrid%Nlvgrid,N_HYDRO,armstradar)  
     if (cfg%Llidar_sim) call construct_cosp_sglidar(Npoints,Ncolumns,Nlevels,N_HYDRO,PARASOL_NREFL,sglidar)
     if (cfg%Llidar_sim) call construct_cosp_lidarstats(Npoints,Ncolumns,vgrid%Nlvgrid,N_HYDRO,PARASOL_NREFL,stlidar)
     if (cfg%Lisccp_sim) call construct_cosp_isccp(Npoints,Ncolumns,Nlevels,isccp)
     if (cfg%Lmodis_sim) call construct_cosp_modis(Npoints,modis)
     if (cfg%Lmisr_sim)  call construct_cosp_misr(Npoints,misr)
     if (cfg%Lrttov_sim) call construct_cosp_rttov(Npoints,Nchannels,rttov)

     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
     ! Call simulator
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
     call cosp(overlap,Ncolumns,cfg,vgrid,gbx,sgx,sgradar,armsgradar,sglidar,isccp,misr,&
               modis,rttov,stradar,armstradar,stlidar)
     write(6,*) 'after cosp'
     call cpu_time(driver_time(3))
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! Write outputs to CMOR-compliant netCDF format.
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     if (cfg%Lwrite_output) then

       ! Model grid info for cmor output
        allocate(mgrid_z(Nlevels),mgrid_zl(Nlevels),mgrid_zu(Nlevels))
        mgrid_z             = zlev(1,Nlevels:1:-1)
        mgrid_zl            = zlev_half(1,Nlevels:1:-1)
        mgrid_zu(2:Nlevels) = zlev_half(1,Nlevels:2:-1)
        mgrid_zu(1)         = zlev(1,Nlevels)+(zlev(1,Nlevels)-mgrid_zl(Nlevels))
        
        N1 = N1D
        if (geomode == 1) N1 = N1D+1
        if (i .eq. 1) then
           call nc_cmor_init(cmor_nl,'replace',cfg,vgrid,gbx,sgx,sglidar,isccp,misr,     &
                             modis,rttov,sgradar,stradar,armsgradar,armstradar,stlidar,  &
                             geomode,Nlon,Nlat,N1, &
                             N2D,N3D,N_OUT_LIST,mgrid_zl,mgrid_zu,mgrid_z,lon_axid,      &
                             lat_axid,time_axid,height_axid, &
                             height_mlev_axid,grid_id,lonvar_id,latvar_id,column_axid,   &
                             sza_axid,temp_axid,channel_axid,dbze_axid,sratio_axid,      &
                             MISR_CTH_axid,tau_axid,pressure2_axid,v1d(1:N1),v2d,v3d)
        endif
        if (geomode == 1) then
           call nc_cmor_associate_1d(grid_id,time_axid,height_axid,height_mlev_axid,     &
                                     column_axid,sza_axid,temp_axid,channel_axid,        &
                                     dbze_axid,sratio_axid,MISR_CTH_axid,tau_axid,       &
                                     pressure2_axid,Nlon,Nlat,vgrid,gbx,sgx,sglidar,     &
                                     isccp,misr,modis,rttov,sgradar,armsgradar,          &
                                     stradar,armstradar,stlidar,N1D,N2D,N3D,v1d(1:N1),   &
                                     v2d,v3d)
           call nc_cmor_write_1d(gbx,time_bnds,lonvar_id,latvar_id,N1,N2D,N3D,v1d(1:N1),v2d, &
                                 v3d)
        elseif (geomode >  1) then
           call nc_cmor_associate_2d(lon_axid,lat_axid,time_axid,height_axid,            &
                                     height_mlev_axid,column_axid,sza_axid,temp_axid,    &
                                     channel_axid,dbze_axid,sratio_axid,MISR_CTH_axid,   &
                                     tau_axid,pressure2_axid,Nlon,Nlat,vgrid,gbx,sgx,    &
                                     sglidar,isccp,misr,modis,rttov,sgradar,armsgradar,  &
                                     stradar,armstradar,stlidar,N1D,N2D,N3D,v1d(1:N1),   &
                                     v2d,v3d)
                                     
           call nc_cmor_write_2d(time,time_bnds,geomode,Nlon,Nlat,N1,N2D,N3D,v1d(1:N1),  &
                                 v2d,v3d)
        endif
        if (i .eq. Nfiles) then
           call nc_cmor_close()
        endif
        deallocate(mgrid_zl,mgrid_zu,mgrid_z)
     endif
     call cpu_time(driver_time(4))

     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
     ! Deallocate memory in derived types
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
     call free_cosp_gridbox(gbx)
     call free_cosp_subgrid(sgx)
     call free_cosp_vgrid(vgrid)
     if (cfg%Lradar_sim) call free_cosp_sgradar(sgradar)
     if (cfg%Lradar_sim) call free_cosp_radarstats(stradar)
     if (cfg%Larmradar_sim) call free_cosp_sgradar(armsgradar)
     if (cfg%Larmradar_sim) call free_cosp_radarstats(armstradar)
     if (cfg%Llidar_sim) call free_cosp_sglidar(sglidar)
     if (cfg%Llidar_sim) call free_cosp_lidarstats(stlidar)
     if (cfg%Lisccp_sim) call free_cosp_isccp(isccp)
     if (cfg%Lmisr_sim)  call free_cosp_misr(misr)
     if (cfg%Lmodis_sim) call free_cosp_modis(modis)
     if (cfg%Lrttov_sim) call free_cosp_rttov(rttov)

  enddo

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  ! Free up local memory
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  deallocate(lon,lat,p,ph,zlev,zlev_half,T,sh,rh,tca,cca, mr_lsliq,mr_lsice,mr_ccliq,    &
             mr_ccice,fl_lsrain,fl_lssnow,fl_lsgrpl,fl_ccrain,fl_ccsnow,Reff,dtau_s,     &
             dtau_c,dem_s,dem_c,skt,landmask,mr_ozone,u_wind,v_wind,sunlit)

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  ! Driver timing
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  call cpu_time(driver_time(5))
  print*,'Time to read in data:     ',driver_time(2)-driver_time(1)
  print*,'Time to run COSP:         ',driver_time(3)-driver_time(2)
  print*,'Time to write output:     ',driver_time(4)-driver_time(3)
  print*,'Total time to run driver: ',driver_time(5)-driver_time(1)
END PROGRAM COSPTEST_trunk
