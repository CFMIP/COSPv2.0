! (c) British Crown Copyright 2008, the Met Office.
! All rights reserved.
! $Revision$, $Date$
! $URL$
! 
! Redistribution and use in source and binary forms, with or without modification, are permitted 
! provided that the following conditions are met:
! 
!     * Redistributions of source code must retain the above copyright notice, this list 
!       of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright notice, this list
!       of conditions and the following disclaimer in the documentation and/or other materials 
!       provided with the distribution.
!     * Neither the name of the Met Office nor the names of its contributors may be used 
!       to endorse or promote products derived from this software without specific prior written 
!       permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
! IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
! FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR 
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
! IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
! OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!
! History:
! Feb 2008 - A. Bodas-Salcedo - Initial version
! Dec 2010 - A. Bodas-Salcedo - Added capability for processing multiple files
!

#include "cosp_defs.h"
PROGRAM COSPTEST
!  USE MOD_COSP_CONSTANTS
!  USE MOD_COSP_TYPES
!  USE MOD_COSP
!  USE MOD_COSP_IO

  USE COSP_KINDS,              ONLY: wp,dp
  USE MOD_COSP_CONFIG,         ONLY: RTTOV_MAX_CHANNELS,N_HYDRO,PARASOL_NREFL
  USE MOD_COSP_IO, 	       ONLY: nc_read_input_file,nc_cmor_init,                    &
                                     nc_cmor_associate_1d,nc_cmor_write_1d_v1p4,         &
                                     nc_cmor_associate_2d,nc_cmor_write_2d,              &
                                     nc_cmor_close,var1d,var2d,var3d,read_cosp_output_nl,&
                                     cosp_error
  USE MOD_COSP_INTERFACE_v1p4, ONLY: cosp                   => cosp_interface_v1p4,      &
                                     cosp_gridbox           => cosp_gridbox_v1p4,        &
                                     construct_cosp_vgrid   => construct_cosp_vgrid_v1p4,&
                                     construct_cosp_gridbox => construct_cosp_gridbox_v1p4,&
                                     free_cosp_gridbox      => destroy_cosp_gridbox_v1p4,&
                                     free_cosp_sgradar      => destroy_cosp_sgradar,     &
                                     free_cosp_radarstats   => destroy_cosp_radarstats,  &
                                     free_cosp_sglidar      => destroy_cosp_sglidar,     &
                                     free_cosp_lidarstats   => destroy_cosp_lidarstats,  &
                                     free_cosp_isccp        => destroy_cosp_isccp,       &
                                     free_cosp_misr         => destroy_cosp_misr,        &
                                     free_cosp_rttov        => destroy_cosp_rttov,       &
                                     free_cosp_modis        => destroy_cosp_modis,       &
                                     free_cosp_vgrid        => destroy_cosp_vgrid,       &
                                     cosp_sglidar,cosp_lidarstats,                       &
                                     construct_cosp_lidarstats,construct_cosp_sglidar,   &
                                     cosp_isccp,construct_cosp_isccp,cosp_misr,          &
                                     construct_cosp_misr,cosp_rttov,construct_cosp_rttov,&
                                     cosp_sgradar,cosp_radarstats,                       & 
                                     construct_cosp_radarstats,construct_cosp_sgradar,   &    
                                     cosp_modis,construct_cosp_modis,                    &
                                     cosp_vgrid
                                        
  USE MOD_COSP_INTERFACE_v1p5, ONLY: free_cosp_subgrid => destroy_cosp_subgrid,          &
                                     construct_cosp_subgrid,cosp_config,cosp_subgrid,    &
                                     I_CVCLIQ,I_LSCLIQ,I_CVCICE,I_LSCICE,I_LSRAIN,       &
                                     I_LSSNOW,I_LSGRPL,I_CVRAIN,I_CVSNOW


  IMPLICIT NONE
  ! Parameters
  integer,parameter ::   &
       N_SIMULATORS = 7,&
       N_OUT_LIST = 63,  & ! Number of possible output variables
       N3D        = 8,   & ! Number of 3D output variables
       N2D        = 14,  & ! Number of 2D output variables
       N1D        = 40     ! Number of 1D output variables
  character*32, dimension(N_SIMULATORS) :: SIM_NAME = (/'Radar','Lidar','ISCCP','MISR ','MODIS','RTTOV','Stats'/)
  integer,dimension(N_SIMULATORS) :: tsim
  data tsim/N_SIMULATORS*0.0/

  ! Local variables
  character(len=64)  :: cosp_input_nl='cosp_input_nl_1D.v1p4.txt'
!   character(len=64)  :: cosp_output_nl='cfmip2/cosp_output_cfmip2_short_offline.txt'
  character(len=64)  :: cosp_output_nl='cosp_output_nl.txt'
  character(len=512) :: dinput ! Directory with input files
  integer,parameter :: N_MAX_INPUT_FILES = 10000 ! Maximum number of input files
  character(len=64),dimension(N_MAX_INPUT_FILES) :: finput ! File names
  character(len=600) :: dfinput ! Input file
  character(len=512) :: cmor_nl
  character(len=8)  :: wmode ! Writing mode 'replace' or 'append'
  integer :: overlap   !  overlap type: 1=max, 2=rand, 3=max/rand
  integer :: isccp_topheight,isccp_topheight_direction
  integer :: Ncolumns     ! Number of subcolumns in SCOPS
  integer :: Npoints      ! Number of gridpoints
  integer :: Nlevels      ! Number of levels
  integer :: Nlr          ! Number of levels in statistical outputs
  integer :: Npoints_it   ! Max number of gridpoints to be processed in one iteration
  integer :: Nfiles       ! Number of files to be processed
  integer,parameter :: ntsteps=5 
  integer :: i,k
  type(cosp_config) :: cfg   ! Configuration options
  type(cosp_gridbox) :: gbx ! Gridbox information. Input for COSP
  type(cosp_subgrid) :: sgx     ! Subgrid outputs
  type(cosp_sgradar) :: sgradar ! Output from radar simulator
  type(cosp_sglidar) :: sglidar ! Output from lidar simulator
  type(cosp_isccp)   :: isccp   ! Output from ISCCP simulator
  type(cosp_modis)   :: modis   ! Output from MODIS simulator
  type(cosp_misr)    :: misr    ! Output from MISR simulator
  type(cosp_rttov)   :: rttov   ! Output from RTTOV 
  type(cosp_vgrid)   :: vgrid   ! Information on vertical grid of stats
  type(cosp_radarstats) :: stradar ! Summary statistics from radar simulator
  type(cosp_lidarstats) :: stlidar ! Summary statistics from lidar simulator
  type(var1d) :: v1d(N1D+1) ! Structures needed by output routines for 1D variables
  type(var2d) :: v2d(N2D) ! Structures needed by output routines for 2D variables
  type(var3d) :: v3d(N3D) ! Structures needed by output routines for 3D variables
  integer ::  grid_id,latvar_id,lonvar_id,lon_axid,lat_axid,time_axid,height_axid,height_mlev_axid,column_axid,sza_axid, &
              temp_axid,channel_axid,dbze_axid,sratio_axid,MISR_CTH_axid,tau_axid,pressure2_axid
  real,dimension(:),allocatable :: lon,lat
  real,dimension(:,:),allocatable,target :: p,ph,zlev,zlev_half,T,sh,rh,tca,cca, &
                    mr_lsliq,mr_lsice,mr_ccliq,mr_ccice,fl_lsrain,fl_lssnow,fl_lsgrpl, &
                    fl_ccrain,fl_ccsnow,dtau_s,dtau_c,dem_s,dem_c,mr_ozone
  real,dimension(:,:,:),allocatable :: Reff
  real,dimension(:),allocatable :: skt,landmask,u_wind,v_wind,sunlit
  integer :: t0,t1,t2,t3,count_rate,count_max
  integer :: Nlon,Nlat,geomode
  real :: radar_freq,k2,ZenAng,co2,ch4,n2o,co,emsfc_lw
  integer,dimension(RTTOV_MAX_CHANNELS) :: Channels
  real,dimension(RTTOV_MAX_CHANNELS) :: Surfem
  integer :: surface_radar,use_mie_tables,use_gas_abs,do_ray,melt_lay
  integer :: Nprmts_max_hydro,Naero,Nprmts_max_aero,lidar_ice_type
  integer :: platform,satellite,Instrument,Nchannels,N1
  logical :: use_vgrid,csat_vgrid,use_precipitation_fluxes,use_reff
  double precision :: time,time_bnds(2),time_step
  real :: toffset_step,half_time_step
  namelist/COSP_INPUT/cmor_nl,overlap,isccp_topheight,isccp_topheight_direction, &
              npoints,npoints_it,ncolumns,nlevels,use_vgrid,nlr,csat_vgrid,dinput,finput, &
              radar_freq,surface_radar,use_mie_tables, &
              use_gas_abs,do_ray,melt_lay,k2,Nprmts_max_hydro,Naero,Nprmts_max_aero, &
              lidar_ice_type,use_precipitation_fluxes,use_reff, &
              platform,satellite,Instrument,Nchannels, &
              Channels,Surfem,ZenAng,co2,ch4,n2o,co

  !---------------- End of declaration of variables --------------
   

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Read COSP namelists
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  finput(:) = ''
  open(10,file=cosp_input_nl,status='old')
  read(10,nml=cosp_input)
  close(10)
  call read_cosp_output_nl(cosp_output_nl,N_OUT_LIST,cfg)
  
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Find number of input files
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  i=1
  do while (i <= N_MAX_INPUT_FILES)
     if (len_trim(finput(i)) < 1) exit
     i=i+1
  enddo
  Nfiles = i-1
  if (Nfiles < 1) call cosp_error('cosp_test','Number of files < 1')
  if (Nfiles > N_MAX_INPUT_FILES) call cosp_error('cosp_test','Number of files > N_MAX_INPUT_FILES')

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Allocate local arrays
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  allocate(lon(Npoints),lat(Npoints), &
           p(Npoints,Nlevels),ph(Npoints,Nlevels), &
           zlev(Npoints,Nlevels),zlev_half(Npoints,Nlevels),T(Npoints,Nlevels), &
           sh(Npoints,Nlevels),rh(Npoints,Nlevels), &
           tca(Npoints,Nlevels),cca(Npoints,Nlevels),mr_lsliq(Npoints,Nlevels), &
           mr_lsice(Npoints,Nlevels),mr_ccliq(Npoints,Nlevels),mr_ccice(Npoints,Nlevels), &
           fl_lsrain(Npoints,Nlevels),fl_lssnow(Npoints,Nlevels),fl_lsgrpl(Npoints,Nlevels),fl_ccrain(Npoints,Nlevels), &
           fl_ccsnow(Npoints,Nlevels),Reff(Npoints,Nlevels,N_HYDRO),dtau_s(Npoints,Nlevels),dtau_c(Npoints,Nlevels), &
           dem_s(Npoints,Nlevels),dem_c(Npoints,Nlevels),skt(Npoints),landmask(Npoints), &
           mr_ozone(Npoints,Nlevels),u_wind(Npoints),v_wind(Npoints),sunlit(Npoints))

  call system_clock(t0,count_rate,count_max) !!! Only for testing purposes

  ! Example that processes ntsteps. It always uses the same input data
  time_step      = 3.D0/24.D0
  time           = 8*1.D0/8.D0 ! First time step
  toffset_step   = time_step/Npoints
  half_time_step = 0.5*time_step
  do i=1,Nfiles
        dfinput=trim(dinput)//trim(finput(i))
        time_bnds = (/time-half_time_step,time+half_time_step/) ! This may need to be adjusted, 
                                                                ! depending on the approx_interval in the MIP table
        print *, 'Processing file: ', trim(dfinput)
        print *, 'Time: ',time
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Read input geophysical variables from NetCDF file
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! input : surface to top
        call nc_read_input_file(dfinput,Npoints,Nlevels,N_HYDRO,lon,lat,p,ph,zlev,zlev_half,T,sh,rh,tca,cca, &
                mr_lsliq,mr_lsice,mr_ccliq,mr_ccice,fl_lsrain,fl_lssnow,fl_lsgrpl,fl_ccrain,fl_ccsnow,Reff, &
                dtau_s,dtau_c,dem_s,dem_c,skt,landmask,mr_ozone,u_wind,v_wind,sunlit, &
                emsfc_lw,geomode,Nlon,Nlat)
                ! geomode = 2 for (lon,lat) mode.
                ! geomode = 3 for (lat,lon) mode.
                ! In those modes it returns Nlon and Nlat with the correct values

        call system_clock(t1,count_rate,count_max)
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Allocate memory for gridbox type
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        print *, 'Allocating memory for gridbox type...'
!ds        call construct_cosp_gridbox(time,time_bnds,radar_freq,surface_radar,use_mie_tables,use_gas_abs, &
!ds                                    do_ray,melt_lay,k2, &
!ds                                    Npoints,Nlevels,Ncolumns,N_HYDRO,Nprmts_max_hydro,Naero,Nprmts_max_aero,Npoints_it, &
!ds                                    lidar_ice_type,isccp_topheight,isccp_topheight_direction,overlap,emsfc_lw, &
!ds                                    use_precipitation_fluxes,use_reff, &
!ds                                    Platform,Satellite,Instrument,Nchannels,ZenAng, &
!ds                                    channels(1:Nchannels),surfem(1:Nchannels),co2,ch4,n2o,co,gbx)
!ds
!ds        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!ds        ! Here code to populate input structure
!ds        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!ds        print *, 'Populating input structure...'
!ds        gbx%longitude = lon
!ds        gbx%latitude = lat
!ds        ! Toffset. This assumes that time is the mid-point of the interval.
!ds        do k=1,Npoints
!ds          gbx%toffset(k) = -half_time_step + toffset_step*(k-0.5)
!ds        enddo
!ds        gbx%p = p
!ds        gbx%ph = ph
!ds        gbx%zlev = zlev
!ds        gbx%zlev_half = zlev_half
!ds        gbx%T = T
!ds        gbx%q = rh
!ds        gbx%sh = sh
!ds        gbx%cca = cca
!ds        gbx%tca = tca
!ds        gbx%psfc = ph(:,1)
!ds        gbx%skt  = skt
!ds        gbx%land = landmask
!ds        gbx%mr_ozone  = mr_ozone
!ds        gbx%u_wind  = u_wind
!ds        gbx%v_wind  = v_wind
!ds        gbx%sunlit  = sunlit
!ds
!ds        gbx%mr_hydro(:,:,I_LSCLIQ) = mr_lsliq
!ds        gbx%mr_hydro(:,:,I_LSCICE) = mr_lsice
!ds        gbx%mr_hydro(:,:,I_CVCLIQ) = mr_ccliq
!ds        gbx%mr_hydro(:,:,I_CVCICE) = mr_ccice
!ds        gbx%rain_ls = fl_lsrain
!ds        gbx%snow_ls = fl_lssnow
!ds        gbx%grpl_ls = fl_lsgrpl
!ds        gbx%rain_cv = fl_ccrain
!ds        gbx%snow_cv = fl_ccsnow
!ds
!ds        gbx%Reff = Reff
!ds        gbx%Reff(:,:,I_LSRAIN) = 0.0
!ds
!ds        ! ISCCP simulator
!ds        gbx%dtau_s   = dtau_s
!ds        gbx%dtau_c   = dtau_c
!ds        gbx%dem_s    = dem_s
!ds        gbx%dem_c    = dem_c

        call construct_cosp_gridbox(time,time_bnds,radar_freq,surface_radar, use_mie_tables,&
                                    use_gas_abs,do_ray,melt_lay,k2,Npoints,Nlevels,Ncolumns,&
                                    N_HYDRO,Nprmts_max_hydro,Naero,Nprmts_max_aero,         &
                                    Npoints_it,lidar_ice_type,isccp_topheight,              &
                                    isccp_topheight_direction,overlap,emsfc_lw,             &
                                    use_precipitation_fluxes,use_reff,Platform,Satellite,   &
                                    Instrument,Nchannels,ZenAng, channels(1:Nchannels),     &
                                    surfem(1:Nchannels),co2,ch4,n2o,co,lon,lat,p,ph,zlev,   &
                                    zlev_half,T,rh,sh,cca,tca,skt,landmask,mr_ozone,u_wind, &
                                    v_wind,sunlit,fl_lsrain,fl_lssnow,fl_lsgrpl,fl_ccrain,  &
                                    fl_ccsnow,dtau_s,dtau_c,dem_s,dem_c,Reff,mr_lsliq,      &
                                    mr_lsice,mr_ccliq,mr_ccice,gbx)
        
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Define new vertical grid
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        print *, 'Defining new vertical grid...'
        call construct_cosp_vgrid(gbx,Nlr,use_vgrid,csat_vgrid,vgrid)

        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Allocate memory for other types
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        print *, 'Allocating memory for other types...'
        call construct_cosp_subgrid(Npoints, Ncolumns, Nlevels, sgx)
!ds        call construct_cosp_sgradar(cfg,Npoints,Ncolumns,Nlevels,N_HYDRO,sgradar)
!ds        call construct_cosp_radarstats(cfg,Npoints,Ncolumns,vgrid%Nlvgrid,N_HYDRO,stradar)
!ds        call construct_cosp_sglidar(cfg,Npoints,Ncolumns,Nlevels,N_HYDRO,PARASOL_NREFL,sglidar)
!ds        call construct_cosp_lidarstats(cfg,Npoints,Ncolumns,vgrid%Nlvgrid,N_HYDRO,PARASOL_NREFL,stlidar)
!ds        call construct_cosp_isccp(cfg,Npoints,Ncolumns,Nlevels,isccp)
!ds        call construct_cosp_modis(cfg,Npoints,modis)
!ds        call construct_cosp_misr(cfg,Npoints,misr)
!ds        call construct_cosp_rttov(cfg,Npoints,Nchannels,rttov)
        call construct_cosp_sgradar(Npoints,Ncolumns,Nlevels,N_HYDRO,sgradar)
        call construct_cosp_radarstats(Npoints,Ncolumns,vgrid%Nlvgrid,N_HYDRO,stradar)  
        call construct_cosp_sglidar(Npoints,Ncolumns,Nlevels,N_HYDRO,PARASOL_NREFL,sglidar)
        call construct_cosp_lidarstats(Npoints,Ncolumns,vgrid%Nlvgrid,N_HYDRO,PARASOL_NREFL,stlidar) 
        call construct_cosp_isccp(Npoints,Ncolumns,Nlevels,isccp)
        call construct_cosp_modis(Npoints,modis)
        call construct_cosp_misr(Npoints,misr)
        call construct_cosp_rttov(Npoints,Nchannels,rttov)
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Call simulator
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        print *, 'Calling simulator...'
        call cosp(overlap,Ncolumns,cfg,vgrid,gbx,sgx,sgradar,sglidar,isccp,misr,modis,rttov,stradar,stlidar)

!ds#ifdef RTTOV
!ds        call cosp(overlap,Ncolumns,cfg,vgrid,gbx,sgx,sgradar,sglidar,isccp,misr,modis,rttov,stradar,stlidar)
!ds#else
!ds        call cosp(overlap,Ncolumns,cfg,vgrid,gbx,sgx,sgradar,sglidar,isccp,misr,modis,stradar,stlidar)
!ds#endif
        call system_clock(t2,count_rate,count_max)

        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Write outputs to CMOR-compliant NetCDF
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Initialise CMOR interface
        gbx%time = time
        ! Write one time step to file
        if (cfg%Lwrite_output) then
            N1 = N1D
            if (geomode == 1) N1 = N1D+1
            print *, 'Writing outputs...'
            if (i .eq. 1) then
               call nc_cmor_init(cmor_nl,'replace',cfg,vgrid,gbx,sgx,sglidar,isccp,misr,     &
                                 modis,rttov,sgradar,stradar,stlidar,geomode,Nlon,Nlat,N1,   &
                                 N2D,N3D,N_OUT_LIST,lon_axid,lat_axid,time_axid,height_axid, &
                                 height_mlev_axid,grid_id,lonvar_id,latvar_id,column_axid,   &
                                 sza_axid,temp_axid,channel_axid,dbze_axid,sratio_axid,      &
                                 MISR_CTH_axid,tau_axid,pressure2_axid,v1d(1:N1),v2d,v3d)
            endif
            if (geomode == 1) then
               call nc_cmor_associate_1d(grid_id,time_axid,height_axid,height_mlev_axid,     &
                                         column_axid,sza_axid,temp_axid,channel_axid,        &
                                         dbze_axid,sratio_axid,MISR_CTH_axid,tau_axid,       &
                                         pressure2_axid,Nlon,Nlat,vgrid,gbx,sgx,sglidar,     &
                                         isccp,misr,modis,rttov,sgradar,stradar,stlidar,     &
                                         N1D,N2D,N3D,v1d(1:N1),v2d,v3d)
               call nc_cmor_write_1d_v1p4(gbx,time_bnds,lonvar_id,latvar_id,N1,N2D,N3D,      &
                                          v1d(1:N1),v2d,v3d)
            elseif (geomode >  1) then
               call nc_cmor_associate_2d(lon_axid,lat_axid,time_axid,height_axid,            &
                                         height_mlev_axid,column_axid,sza_axid,temp_axid,    &
                                         channel_axid,dbze_axid,sratio_axid,MISR_CTH_axid,   &
                                         tau_axid,pressure2_axid,Nlon,Nlat,vgrid,gbx,sgx,    &
                                         sglidar,isccp,misr,modis,rttov,sgradar,stradar,     &
                                         stlidar,N1D,N2D,N3D,v1d(1:N1),v2d,v3d)
               call nc_cmor_write_2d(time,time_bnds,geomode,Nlon,Nlat,N1,N2D,N3D,v1d(1:N1),  &
                                     v2d,v3d)
            endif
!ds            if (i == 1) call nc_cmor_init(cmor_nl,'replace',cfg,vgrid,gbx,sgx,sgradar,sglidar, &
!ds                                  isccp,misr,modis,rttov,stradar,stlidar,geomode,Nlon,Nlat,N1,N2D,N3D, &
!ds                                  lon_axid,lat_axid,time_axid,height_axid,height_mlev_axid,grid_id,lonvar_id,latvar_id, &
!ds                                  column_axid,sza_axid,temp_axid,channel_axid,dbze_axid, &
!ds                                  sratio_axid,MISR_CTH_axid,tau_axid,pressure2_axid, &
!ds                                  v1d(1:N1),v2d,v3d)
!ds            if (geomode == 1) then
!ds               print *, 'Associate'
!ds               call nc_cmor_associate_1d(grid_id,time_axid,height_axid,height_mlev_axid,column_axid,sza_axid, &
!ds                         temp_axid,channel_axid,dbze_axid,sratio_axid,MISR_CTH_axid,tau_axid,pressure2_axid, &
!ds                         Nlon,Nlat,vgrid,gbx,sgx,sgradar,sglidar,isccp,misr,modis,rttov,stradar,stlidar, &
!ds                         v1d(1:N1),v2d,v3d)
!ds               print *, 'Write'
!ds               call nc_cmor_write_1d(gbx,time_bnds,lonvar_id,latvar_id,N1,N2D,N3D,v1d(1:N1),v2d,v3d)
!ds            elseif (geomode >  1) then
!ds               call nc_cmor_associate_2d(lon_axid,lat_axid,time_axid,height_axid,height_mlev_axid,column_axid,sza_axid, &
!ds                         temp_axid,channel_axid,dbze_axid,sratio_axid,MISR_CTH_axid,tau_axid,pressure2_axid, &
!ds                         Nlon,Nlat,vgrid,gbx,sgx,sgradar,sglidar,isccp,misr,modis,rttov,stradar,stlidar, &
!ds                         v1d(1:N1),v2d,v3d)
!ds               call nc_cmor_write_2d(time,time_bnds,geomode,Nlon,Nlat,N1,N2D,N3D,v1d(1:N1),v2d,v3d)
!ds
!ds            endif
            if (i == Nfiles) then
              call nc_cmor_close()
            endif
        endif

        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Deallocate memory in derived types
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        print *, 'Deallocating memory...'
        call free_cosp_gridbox(gbx)
        call free_cosp_subgrid(sgx)
        call free_cosp_sgradar(sgradar)
        call free_cosp_radarstats(stradar)
        call free_cosp_sglidar(sglidar)
        call free_cosp_lidarstats(stlidar)
        call free_cosp_isccp(isccp)
        call free_cosp_misr(misr)
        call free_cosp_modis(modis)
        call free_cosp_rttov(rttov)
        call free_cosp_vgrid(vgrid)
        ! Update time
        time = time + time_step
  enddo
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Deallocate memory in local arrays
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  deallocate(lon,lat,p,ph,zlev,zlev_half,T,sh,rh,tca,cca, mr_lsliq,mr_lsice,mr_ccliq,mr_ccice, &
           fl_lsrain,fl_lssnow,fl_lsgrpl,fl_ccrain,fl_ccsnow,Reff,dtau_s,dtau_c,dem_s,dem_c,skt, &
           landmask,mr_ozone,u_wind,v_wind,sunlit)

  ! Time in s. Only for testing purposes
  call system_clock(t3,count_rate,count_max)
  do i=1,N_SIMULATORS
      print *,'=== '//trim(SIM_NAME(i))//': ', float(tsim(i))/count_rate
  enddo
  print *,'=== COSP: ', (t2-t1)*1.0/count_rate
  print *,'=== TOTAL: ', (t3-t0)*1.0/count_rate
    
END PROGRAM COSPTEST
