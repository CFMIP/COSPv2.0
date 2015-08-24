! (c) British Crown Copyright 2008, the Met Office.
! All rights reserved.
! $Revision: 80 $, $Date: 2013-09-12 10:14:54 -0600 (Thu, 12 Sep 2013) $
! $URL: https://cfmip-obs-sim.googlecode.com/svn/devel/branches/dustinswales/cosp_test.F90 $
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
  USE COSP_KINDS,      ONLY: wp
  USE MOD_COSP_CONFIG
  USE MOD_COSP_TYPES,  ONLY:cosp_config,cosp_vgrid,cosp_subgrid,cosp_gridbox
  USE MOD_COSP_INTERFACE
  USE MOD_COSP_IO
  IMPLICIT NONE

  ! Parameters
  integer,     parameter :: ntsteps=5 
  integer,     parameter :: N_MAX_INPUT_FILES = 10000 ! Maximum number of input files
  character*64,parameter :: cosp_input_nl='cosp_input_nl.txt'
  character*64,parameter :: cosp_output_nl='cosp_output_nl.txt'

  ! Local variables
  character*64,dimension(N_MAX_INPUT_FILES) :: finput ! File names
  character*64 :: radar_micro_scheme
  character*512 :: dinput,             & ! Directory with input files
                   cmor_nl               !
  character*600 :: dfinput               ! Input file
  character*8   :: wmode                 ! Writing mode 'replace' or 'append'
  integer       :: overlap, isccp_topheight,isccp_topheight_direction,Ncolumns, Npoints,            & 
                   Nlevels,Nlr, Npoints_it,Nfiles,i,k,grid_id,latvar_id,lonvar_id,lon_axid,         &
                   lat_axid,time_axid, height_axid,height_mlev_axid,column_axid,sza_axid,           &
                   temp_axid,channel_axid,dbze_axid,sratio_axid,MISR_CTH_axid,tau_axid,             &
                   pressure2_axid,t0,t1,t2,t3,count_rate,count_max,Nlon,Nlat,geomode,               &
                   surface_radar,use_mie_tables,use_gas_abs,do_ray,melt_lay,Nprmts_max_hydro,       &
                   Naero,Nprmts_max_aero,lidar_ice_type,platform,satellite,Instrument,Nchannels,N1, &
                   modis_tau_inputs
  integer,dimension(RTTOV_MAX_CHANNELS)  :: Channels

  type(cosp_config)     :: cfg        ! Configuration options
  type(cosp_gridbox)    :: gbx        ! Gridbox information. Input for COSP
  type(cosp_subgrid)    :: sgx        ! Subgrid outputs
  type(cosp_vgrid)      :: vgrid      ! Information on vertical grid of stats
  type(cosp_lidarstats) :: stlidar    ! Summary statistics from lidar simulator
  type(cosp_outputs)    :: cospOUT
  type(var1d)           :: v1d(N1D+1) ! Structures needed by output routines for 1D variables
  type(var2d)           :: v2d(N2D)   ! Structures needed by output routines for 2D variables
  type(var3d)           :: v3d(N3D)   ! Structures needed by output routines for 3D variables
  real(wp)                                     :: radar_freq,k2,ZenAng,co2,ch4,n2o,co,      &
                                                  emsfc_lw,toffset_step,half_time_step
  real(wp),dimension(:),    allocatable        :: lon,lat,skt,landmask,u_wind,v_wind,sunlit
  real(wp),dimension(:,:),  allocatable,target :: p,ph,zlev,zlev_half,T,sh,rh,tca,cca,     &
                                                  mr_lsliq,mr_lsice,mr_ccliq,mr_ccice,     &
                                                  fl_lsrain,fl_lssnow,fl_lsgrpl,fl_ccrain, &
                                                  fl_ccsnow,dtau_s,dtau_c,dem_s,dem_c,mr_ozone
  real(wp),dimension(:,:,:),allocatable        :: Reff
  real(wp),dimension(RTTOV_MAX_CHANNELS)       :: Surfem
  logical                                      :: use_vgrid,csat_vgrid,use_precipitation_fluxes,use_reff
  double precision                             :: time,time_bnds(2),time_step

  namelist/COSP_INPUT/cmor_nl,overlap,isccp_topheight,isccp_topheight_direction, &
              npoints,npoints_it,ncolumns,nlevels,use_vgrid,nlr,csat_vgrid,dinput,finput, &
              radar_freq,surface_radar,use_mie_tables, &
              use_gas_abs,do_ray,melt_lay,k2,radar_micro_scheme,Nprmts_max_hydro,Naero,Nprmts_max_aero, &
              lidar_ice_type,use_precipitation_fluxes,use_reff, &
              platform,satellite,Instrument,Nchannels, &
              Channels,Surfem,ZenAng,co2,ch4,n2o,co,MODIS_TAU_INPUTS

  !---------------- End of declaration of variables --------------


  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Read COSP namelists
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  finput(:) = ''
  open(10,file=cosp_input_nl,status='old')
  read(10,nml=cosp_input)
  close(10)
  call read_cosp_output_nl(cosp_output_nl,cfg)
  
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
  time_step      = 3._wp/24._wp!3.D0/24.D0
  time           = 8*1._wp/8._wp!8*1.D0/8.D0 ! First time step
  toffset_step   = time_step/Npoints
  half_time_step = 0.5_wp*time_step
  print*,'pre-loop'
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
        call construct_cosp_gridbox(time,time_bnds,radar_freq,surface_radar,use_mie_tables,use_gas_abs, &
                                    do_ray,melt_lay,k2, &
                                    Npoints,Nlevels,Ncolumns,N_HYDRO,Nprmts_max_hydro,Naero,Nprmts_max_aero,Npoints_it, &
                                    lidar_ice_type,isccp_topheight,isccp_topheight_direction,overlap,emsfc_lw, &
                                    use_precipitation_fluxes,use_reff, &
                                    Platform,Satellite,Instrument,Nchannels,ZenAng, &
                                    channels(1:Nchannels),surfem(1:Nchannels),co2,ch4,n2o,co,modis_tau_inputs,gbx)

        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Here code to populate input structure
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        print *, 'Populating input structure...'
        gbx%longitude = lon
        gbx%latitude = lat
        ! Toffset. This assumes that time is the mid-point of the interval.
        do k=1,Npoints
          gbx%toffset(k) = -half_time_step + toffset_step*(k-0.5)
        enddo
        gbx%p = p
        gbx%ph = ph
        gbx%zlev = zlev
        gbx%zlev_half = zlev_half
        gbx%T = T
        gbx%q = rh
        gbx%sh = sh
        gbx%cca = cca
        gbx%tca = tca
        gbx%psfc = ph(:,1)
        gbx%skt  = skt
        gbx%land = landmask
        gbx%mr_ozone  = mr_ozone
        gbx%u_wind  = u_wind
        gbx%v_wind  = v_wind
        gbx%sunlit  = sunlit

        gbx%mr_hydro(:,:,I_LSCLIQ) = mr_lsliq
        gbx%mr_hydro(:,:,I_LSCICE) = mr_lsice
        gbx%mr_hydro(:,:,I_CVCLIQ) = mr_ccliq
        gbx%mr_hydro(:,:,I_CVCICE) = mr_ccice
        gbx%rain_ls = fl_lsrain
        gbx%snow_ls = fl_lssnow
        gbx%grpl_ls = fl_lsgrpl
        gbx%rain_cv = fl_ccrain
        gbx%snow_cv = fl_ccsnow

        gbx%Reff = Reff
        gbx%Reff(:,:,I_LSRAIN) = 0.0

        ! ISCCP simulator
        gbx%dtau_s   = dtau_s
        gbx%dtau_c   = dtau_c
        gbx%dem_s    = dem_s
        gbx%dem_c    = dem_c

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
        call construct_cosp_sgcloudsat(cfg,Npoints,Ncolumns,Nlevels,N_HYDRO,cospOUT%cloudsat%cloudsat_sgx)
        call construct_cosp_cloudsatstats(cfg,Npoints,Ncolumns,vgrid%Nlvgrid,N_HYDRO,cospOUT%cloudsat%cloudsat_gbx)
        call construct_cosp_sglidar(cfg,Npoints,Ncolumns,Nlevels,N_HYDRO,cospOUT%sglidar)
        call construct_cosp_lidarstats(cfg,Npoints,Ncolumns,vgrid%Nlvgrid,N_HYDRO,PARASOL_NREFL,cospOUT%stlidar)
        call construct_cosp_isccp(cfg,Npoints,Ncolumns,Nlevels,cospOUT%isccp)
        call construct_cosp_modis(cfg,Npoints,cospOUT%modis)
        call construct_cosp_misr(cfg,Npoints,Ncolumns,Nlevels,cospOUT%misr)
        call construct_cosp_rttov(cfg,Npoints,Nchannels,cospOUT%rttov)
        call construct_cosp_parasol(cfg,Npoints,Ncolumns,PARASOL_NREFL,cospOUT%parasol)

        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Call simulator
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        print *, 'Calling simulator...'
        
        call cosp(overlap,Ncolumns,cfg,vgrid,gbx,sgx,cospOUT)

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
            if (i == 1) call nc_cmor_init(cmor_nl,'replace',cfg,vgrid,gbx,sgx,cospOUT%sglidar, &
                                  cospOUT%isccp,cospOUT%misr,cospOUT%modis,cospOUT%rttov,cospOUT%cloudsat,cospOUT%stlidar,cospOUT%parasol,geomode,Nlon,Nlat,N1,N2D,N3D, &
                                  lon_axid,lat_axid,time_axid,height_axid,height_mlev_axid,grid_id,lonvar_id,latvar_id, &
                                  column_axid,sza_axid,temp_axid,channel_axid,dbze_axid, &
                                  sratio_axid,MISR_CTH_axid,tau_axid,pressure2_axid, &
                                  v1d(1:N1),v2d,v3d)
            if (geomode == 1) then
               print *, 'Associate'
               call nc_cmor_associate_1d(grid_id,time_axid,height_axid,height_mlev_axid,column_axid,sza_axid, &
                         temp_axid,channel_axid,dbze_axid,sratio_axid,MISR_CTH_axid,tau_axid,pressure2_axid, &
                         Nlon,Nlat,vgrid,gbx,sgx,cospOUT%sglidar,cospOUT%isccp,cospOUT%misr,cospOUT%modis,cospOUT%rttov,cospOUT%cloudsat,cospOUT%stlidar,cospOUT%parasol, &
                         v1d(1:N1),v2d,v3d)
               print *, 'Write'
               call nc_cmor_write_1d_new(gbx,time_bnds,lonvar_id,latvar_id,N1,N2D,N3D,v1d(1:N1),v2d,v3d)
            elseif (geomode >  1) then
               call nc_cmor_associate_2d(lon_axid,lat_axid,time_axid,height_axid,height_mlev_axid,column_axid,sza_axid, &
                         temp_axid,channel_axid,dbze_axid,sratio_axid,MISR_CTH_axid,tau_axid,pressure2_axid, &
                         Nlon,Nlat,vgrid,gbx,sgx,cospOUT%sglidar,cospOUT%isccp,cospOUT%misr,cospOUT%modis,cospOUT%rttov,cospOUT%cloudsat,cospOUT%stlidar,cospOUT%parasol, &
                         v1d(1:N1),v2d,v3d)
               call nc_cmor_write_2d(time,time_bnds,geomode,Nlon,Nlat,N1,N2D,N3D,v1d(1:N1),v2d,v3d)

            endif
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
        call free_cosp_vgrid(vgrid)
        call free_cosp_sgcloudsat(cospOUT%cloudsat%cloudsat_sgx)
        call free_cosp_cloudsatstats(cospOUT%cloudsat%cloudsat_gbx)
        call free_cosp_sglidar(cospOUT%sglidar)
        call free_cosp_lidarstats(cospOUT%stlidar)
        call free_cosp_isccp(cospOUT%isccp)
        call free_cosp_misr(cospOUT%misr)
        call free_cosp_modis(cospOUT%modis)
        call free_cosp_rttov(cospOUT%rttov)
        call free_cosp_parasol(cospOUT%parasol)
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
