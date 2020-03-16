! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copyright (c) 2020, Regents of the University of Colorado
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
! April 2018 - R. Guzman - Added OPAQ diagnostics and Ground LIDar (GLID) simulator
! April 2018 - R. Guzman - Added ATLID simulator
!   Nov 2018 - T. Michibata - Added CloudSat+MODIS Warmrain Diagnostics
! March 2020 - D. Swales - Modified to use CESM2 subcolumn inputs.
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
program cosp2_cesm2_test
  use cosp_kinds,      only: wp                         
  USE MOD_COSP_CONFIG, only: R_UNDEF,PARASOL_NREFL,LIDAR_NCAT,LIDAR_NTYPE,SR_BINS,    &
       numMISRHgtBins, numMISRTauBins, cloudsat_DBZE_BINS,LIDAR_NTEMP, CFODD_NDBZE,   &
       CFODD_NICOD, CFODD_BNDRE, CFODD_NCLASS, CFODD_DBZE_MIN, CFODD_DBZE_MAX,        &
       CFODD_ICOD_MIN, CFODD_ICOD_MAX, CFODD_DBZE_WIDTH, CFODD_ICOD_WIDTH, WR_NREGIME,&
       numMODISTauBins, numMODISPresBins, numMODISReffIceBins, numMODISReffLiqBins,   &
       numISCCPTauBins, numISCCPPresBins
  use mod_cosp_io,     only: nc_read_input_file,write_cosp2_output
  use quickbeam,       only: radar_cfg
  use mod_cosp,        only: cosp_init,cosp_optical_inputs,cosp_column_inputs,        &
       cosp_outputs,cosp_cleanUp,cosp_simulator
  use netcdf
  
  implicit none

  ! Input/Output driver file control
  character(len=64),parameter :: &
       cosp_input_namelist  = 'cosp2_cesm2_input_nl.txt', &
       cosp_output_namelist = 'cosp2_cesm2_output_nl.txt'

  ! Test data
  ! 1D [nPoints]
  real(wp),dimension(:),allocatable,target:: &
       lon,                 & ! Longitude (deg)
       lat,                 & ! Latitude (deg)
       skt,                 & ! Skin temperature (K)
       surfelev,            & ! Surface Elevation (m)
       landmask,            & ! Land/sea mask (0/1) LANDFRAC
       sunlit                 ! Sunlit flag
  ! 2D [nPoints, nLevels]
  real(wp),dimension(:,:),allocatable,target :: &
       p,                   & ! Model pressure levels (pa)
       ph,                  & ! Moddel pressure @ half levels (pa)
       zlev,                & ! Model level height (m)
       zlev_half,           & ! Model level height @ half-levels (m)
       T,                   & ! Temperature (K)
       q                      ! Specific humidity (kg/kg)
  real(wp),dimension(:,:,:),allocatable,target :: &
       cld_frac,            & ! Subcolumn cloud-fraction
       tau_067,             & ! ISCCP/MISR/MODIS: 0.67micron optical depth
       emiss_11,            & ! MISR: 11micron emissivity
       MODIS_fracliq,       & ! MODIS: Fraction of tau due to liquid
       MODIS_asym,          & ! MODIS: Asymmetry parater
       MODIS_ssa,           & ! MODIS: Single-scattering albedo
       calipso_betatot,     & ! CALIPSO(LIDAR): Lidar backscatter coefficient
       calipso_betatot_ice, & ! CALIPSO(LIDAR): Lidar backscatter coefficient (ice)
       calipso_betatot_liq, & ! CALIPSO(LIDAR): Lidar backscatter coefficient (liquid)
       calipso_tautot,      & ! CALIPSO(LIDAR): Vertically integreated od 
       calipso_tautot_ice,  & ! CALIPSO(LIDAR): Vertically integreated od (ice)
       calipso_tautot_liq,  & ! CALIPSO(LIDAR): Vertically integreated od (liquid)
       cloudsat_z_vol,      & ! CLOUDSAT(RADAR): Effective reflectivity factor
       cloudsat_kr_vol,     & ! CLOUDSAT(RADAR): Attenuation coefficient (hydrometeors)
       cloudsat_g_vol         ! CLOUDSAT(RADAR): Attenuation coefficient (gases only)

  ! Input namelist fields
  integer ::                      & !
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
       isccp_topheight_direction    ! ISCCP cloud top height direction

  real(wp),dimension(:),allocatable :: & 
       vgrid_z                      ! mid-level altitude of the vertical grid
  real(wp) ::                     & !
       cloudsat_radar_freq,       & ! CloudSat radar frequency (GHz)
       cloudsat_k2                  ! |K|^2, -1=use frequency dependent default

  logical ::                      & !
       use_vgrid,                 & ! Use fixed vertical grid for outputs?
       csat_vgrid,                & ! CloudSat vertical grid? 
       use_precipitation_fluxes     ! True if precipitation fluxes are input to the 
                                    ! algorithm 
  character(len=64) :: &
       cloudsat_micro_scheme        ! Microphysical scheme used in cloudsat radar simulator
  character(len=256) :: &
       finput                       ! Input NetCDF file
  character(len=256) :: &
       foutput
  character(len=512) :: &
       dinput                       ! Directory where the input files are located
  character(len=600) :: &
       fileIN                       ! dinput+finput
  namelist/COSP_INPUT/overlap, isccp_topheight, isccp_topheight_direction,      &
       npoints_it, use_vgrid, Nlvgrid, csat_vgrid, dinput, finput,    &
       foutput, cloudsat_radar_freq, surface_radar, cloudsat_use_gas_abs,cloudsat_do_ray,&
       cloudsat_k2, cloudsat_micro_scheme, lidar_ice_type, use_precipitation_fluxes

  ! Output namelist
  logical :: Lcfaddbze94,Ldbze94,Latb532,LcfadLidarsr532,Lclcalipso,Lclhcalipso,         &
             Lcllcalipso,Lclmcalipso,Lcltcalipso,LparasolRefl,Lclcalipsoliq,             &
             Lclcalipsoice,Lclcalipsoun,Lclcalipsotmp,Lclcalipsotmpliq,Lclcalipsotmpice, &
             Lclcalipsotmpun,Lclhcalipsoliq,Lcllcalipsoliq,Lclmcalipsoliq,Lcltcalipsoliq,& 
             Lclhcalipsoice,Lcllcalipsoice,Lclmcalipsoice,Lcltcalipsoice,Lclhcalipsoun,  & 
             Lcllcalipsoun,Lclmcalipsoun,Lcltcalipsoun,Lclopaquecalipso,Lclthincalipso,  & 
             Lclzopaquecalipso,Lclcalipsoopaque,Lclcalipsothin,Lclcalipsozopaque,        & 
             Lclcalipsoopacity,Lclopaquetemp,Lclthintemp,Lclzopaquetemp,Lclopaquemeanz,  &
             Lclthinmeanz,Lclthinemis,Lclopaquemeanzse,Lclthinmeanzse,Lclzopaquecalipsose,&
             LlidarBetaMol532gr,LcfadLidarsr532gr,Latb532gr,LclgrLidar532,LclhgrLidar532,&
             LcllgrLidar532,LclmgrLidar532,LcltgrLidar532,LlidarBetaMol355,              &
             LcfadLidarsr355,Latb355,Lclatlid,Lclhatlid,Lcllatlid,Lclmatlid,Lcltatlid,   &
             Lalbisccp,Lboxptopisccp,Lboxtauisccp,Lpctisccp,Lclisccp,Ltauisccp,Lcltisccp,&
             Lmeantbisccp,Lmeantbclrisccp,LclMISR,Lclcalipso2,Lcltlidarradar,            &
             Lcloudsat_tcc, Lcloudsat_tcc2,Lfracout,   &
             LlidarBetaMol532,Lcltmodis,Lclwmodis,Lclimodis,Lclhmodis,Lclmmodis,         &
             Lcllmodis,Ltautmodis,Ltauwmodis,Ltauimodis,Ltautlogmodis,Ltauwlogmodis,     &
             Ltauilogmodis,Lreffclwmodis,Lreffclimodis,Lpctmodis,Llwpmodis,Liwpmodis,    &
             Lclmodis,Lptradarflag0,Lptradarflag1,Lptradarflag2,Lptradarflag3,  &
             Lptradarflag4,Lptradarflag5,Lptradarflag6,Lptradarflag7,Lptradarflag8,      &
             Lptradarflag9,Lradarpia,                                                    &
             Lwr_occfreq, Lcfodd
  namelist/COSP_OUTPUT/Lcfaddbze94,Ldbze94,Latb532,LcfadLidarsr532,Lclcalipso,           &
                       Lclhcalipso,Lcllcalipso,Lclmcalipso,Lcltcalipso,LparasolRefl,     &
                       Lclcalipsoliq,Lclcalipsoice,Lclcalipsoun,Lclcalipsotmp,           &
                       Lclcalipsotmpliq,Lclcalipsotmpice,Lclcalipsotmpun,Lclhcalipsoliq, &
                       Lcllcalipsoliq,Lclmcalipsoliq,Lcltcalipsoliq,Lclhcalipsoice,      &
                       Lcllcalipsoice,Lclmcalipsoice,Lcltcalipsoice,Lclhcalipsoun,       &
                       Lcllcalipsoun,Lclmcalipsoun,Lcltcalipsoun,Lclopaquecalipso,       &
                       Lclthincalipso,Lclzopaquecalipso,Lclcalipsoopaque,Lclcalipsothin, &
                       Lclcalipsozopaque,Lclcalipsoopacity,Lclopaquetemp,Lclthintemp,    &
                       Lclzopaquetemp,Lclopaquemeanz,Lclthinmeanz,Lclthinemis,           &
                       Lclopaquemeanzse,Lclthinmeanzse,Lclzopaquecalipsose,              &
                       LlidarBetaMol532gr,LcfadLidarsr532gr,Latb532gr,LclgrLidar532,     &
                       LclhgrLidar532,LcllgrLidar532,LclmgrLidar532,LcltgrLidar532,      &
                       LlidarBetaMol355,LcfadLidarsr355,Latb355,Lclatlid,                &
                       Lclhatlid,Lcllatlid,Lclmatlid,Lcltatlid,Lalbisccp,Lboxptopisccp,  &
                       Lboxtauisccp,Lpctisccp,Lclisccp,Ltauisccp,Lcltisccp,Lmeantbisccp, &
                       Lmeantbclrisccp,LclMISR,Lclcalipso2,Lcltlidarradar,               &
                       Lcloudsat_tcc, Lcloudsat_tcc2, Lfracout,                          &
                       LlidarBetaMol532,Lcltmodis,Lclwmodis,Lclimodis,Lclhmodis,         &
                       Lclmmodis,Lcllmodis,Ltautmodis,Ltauwmodis,Ltauimodis,             &
                       Ltautlogmodis,Ltauwlogmodis,Ltauilogmodis,Lreffclwmodis,          &
                       Lreffclimodis,Lpctmodis,Llwpmodis,Liwpmodis,Lclmodis,             &
                       Lptradarflag0,Lptradarflag1,Lptradarflag2,Lptradarflag3,          &
                       Lptradarflag4,Lptradarflag5,Lptradarflag6,Lptradarflag7,          &
                       Lptradarflag8,Lptradarflag9,Lradarpia,                            &
                       Lwr_occfreq, Lcfodd

  ! Local variables
  logical :: &
       lisccp      = .false. ,& ! Local on/off switch for simulators (used by initialization)
       lmodis      = .false., & !
       lmisr       = .false., & !
       lcalipso    = .false., & !
       lgrLidar532 = .false., & !
       latlid      = .false., & !
       lcloudsat   = .false., & !
       lparasol    = .false.    !
  type(radar_cfg) :: &
       rcfg_cloudsat     ! Radar configuration
  type(cosp_outputs) :: &
       cospOUT           ! COSP simulator outputs
  type(cosp_optical_inputs) :: &
       cospIN            ! COSP optical (or derived?) fields needed by simulators
  type(cosp_column_inputs) :: &
       cospstateIN       ! COSP model fields needed by simulators
  integer :: iChunk,nChunks,start_idx,end_idx,nPtsPerIt,ij, status,ncid,  ndims, nvars, ngatts, recdim, &
       dimID(4), varID(40),nLon, nLat, nColumns, nLevels, nPoints
  real(wp),dimension(10) :: driver_time
  character(len=256),dimension(100) :: cosp_status
  character(len=256) :: varName
  real(wp),dimension(:,:),    allocatable :: temp2D
  real(wp),dimension(:,:,:),  allocatable :: temp3D
  real(wp),dimension(:,:,:,:),allocatable :: temp4D
  ! Debug fields
  integer, parameter :: nLon_d=2
  
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
       Lcltlidarradar .or. Lclopaquecalipso .or. Lclthincalipso .or. Lclzopaquecalipso   &
       .or. Lclcalipsoopaque .or. Lclcalipsothin .or. Lclcalipsozopaque .or.             &
       Lclcalipsoopacity .or. Lclopaquetemp .or. Lclthintemp .or. Lclzopaquetemp .or.    &
       Lclopaquemeanz .or. Lclthinmeanz .or. Lclthinemis .or. Lclopaquemeanzse .or.      & 
       Lclthinmeanzse .or. Lclzopaquecalipsose) Lcalipso = .true. 

  if (LlidarBetaMol532gr .or. LcfadLidarsr532gr .or. Latb532gr .or. LclgrLidar532 .or.   & 
       LclhgrLidar532 .or. LcllgrLidar532 .or. LclmgrLidar532 .or. LcltgrLidar532)       & 
       LgrLidar532 = .true.

  if (LlidarBetaMol355 .or. LcfadLidarsr355 .or. Latb355 .or. Lclatlid .or.              & 
       Lclhatlid .or. Lcllatlid .or. Lclmatlid .or. Lcltatlid)                           & 
       Latlid = .true. 

  if (LcfadDbze94 .or. Ldbze94 .or. Lcltlidarradar) Lcloudsat = .true.

  if (LcfadDbze94 .or. Ldbze94 .or. Lcltlidarradar .or. Lptradarflag0 .or. Lptradarflag1 &
       .or. Lptradarflag2 .or. Lptradarflag3 .or. Lptradarflag4 .or. Lptradarflag5 .or.  &
       Lptradarflag6 .or. Lptradarflag7 .or. Lptradarflag8 .or. Lptradarflag9 .or.       &
       Lradarpia) Lcloudsat = .true.
  if (Lparasolrefl) Lparasol = .true.

  if (lisccp)    print*,'Running ISCCP simulator'
  if (lmisr)     print*,'Running MISR simulator'
  if (lmodis)    print*,'Running MODIS simulator'
  if (lcloudsat) print*,'Running Cloudsat RADAR simulator'
  if (lcalipso)  print*,'Running CALIPSO LIDAR simulator'
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Read in sample input data.
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fileIN = trim(dinput)//trim(finput)
  print*,fileIN
  status = nf90_open(fileIN, nf90_nowrite, ncid)
  status = nf90_inq_dimid(ncid,'lev',dimID(1))
  status = nf90_inquire_dimension(ncid, dimID(1), varName, len = nLevels) 
  status = nf90_inq_dimid(ncid,'cosp_scol',dimID(2))
  status = nf90_inquire_dimension(ncid, dimID(2), varName, len = nColumns) 
  status = nf90_inq_dimid(ncid,'lon',dimID(3))
  status = nf90_inquire_dimension(ncid, dimID(3), varName, len = nLon) 
  status = nf90_inq_dimid(ncid,'lat',dimID(4))
  status = nf90_inquire_dimension(ncid, dimID(4), varName, len = nLat)

  nLon = nLon_d
  
  ! Total number of points is nLon x nLat
  nPoints = nLon*nLat
  print*,nLon,nLat,nPoints,nLevels

  !
  ! Read in fields, squash along lon/lat dimension
  !
  
  ! 1D
  allocate(temp2D(nLon,nLat),lon(nPoints), lat(nPoints), skt(nPoints), &
       surfelev(nPoints), landmask(nPoints), sunlit(nPoints))
  status = nf90_inq_varid(ncid,"lon",varID(1))
  status = nf90_get_var(ncid,varID(1),temp2D,count=(/nLon,nLat/))
  lon    = reshape(temp2D,(/nPoints/))
  status = nf90_inq_varid(ncid,"lat",varID(2))
  status = nf90_get_var(ncid,varID(2),temp2D,count=(/nLon,nLat/))
  lat    = reshape(temp2D,(/nPoints/))
  status = nf90_inq_varid(ncid,"TS_COSP",varID(3))
  status = nf90_get_var(ncid,varID(3),temp2D,count=(/nLon,nLat/))
  skt    = reshape(temp2D,(/nPoints/))
  ! Surface elevation
  ! Land/sea mask
  ! Sunlit
  ! 2D
  allocate(temp3D(nLon,nLat,nLevels), p(nPoints,nLevels), &
       ph(nPoints,nLevels), zlev(nPoints,nLevels),      &
       zlev_half(nPoints,nLevels+1), T(nPoints,nLevels),  &
       q(nPoints,nLevels))
  status    = nf90_inq_varid(ncid,"P_COSP",varID(7))
  status    = nf90_get_var(ncid,varID(7),temp3D,count=(/nLon,nLat,nLevels/))
  p         = reshape(temp3D,(/nPoints, nLevels/))
  status    = nf90_inq_varid(ncid,"PH_COSP",varID(8))  
  status    = nf90_get_var(ncid,varID(8),temp3D,count=(/nLon,nLat,nLevels/))
  ph        = reshape(temp3D,(/nPoints, nLevels/))
  status    = nf90_inq_varid(ncid,"ZLEV_COSP",varID(9))  
  status    = nf90_get_var(ncid,varID(9),temp3D,count=(/nLon,nLat,nLevels/))
  zlev      = reshape(temp3D,(/nPoints, nLevels/))
  status    = nf90_inq_varid(ncid,"ZLEV_HALF_COSP",varID(10))  
  status    = nf90_get_var(ncid,varID(10),temp3D,count=(/nLon,nLat,nLevels/))
  zlev_half = reshape(temp3D,(/nPoints, nLevels/))
  status    = nf90_inq_varid(ncid,"T_COSP",varID(11))  
  status    = nf90_get_var(ncid,varID(11),temp3D,count=(/nLon,nLat,nLevels/))
  T         = reshape(temp3D,(/nPoints, nLevels/))
  status    = nf90_inq_varid(ncid,"RH_COSP",varID(12))  
  status    = nf90_get_var(ncid,varID(12),temp3D,count=(/nLon,nLat,nLevels/))
  q         = reshape(temp3D,(/nPoints, nLevels/))

  ! 3D, only when simulator is requested
  print*,'Reading in subcolumn inputs...'     
  allocate(temp4D(nLon, nLat, nColumns, nLevels))
  status  = nf90_inq_varid(ncid,"SCOPS_OUT",varID(14))  
  status  = nf90_get_var(ncid,varID(14),temp4D,count=(/nLon,nLat,nColumns,nLevels/))
  cld_frac = reshape(temp4D,(/nPoints, nColumns, nLevels/))
  
  if (lmisr .or. lmodis .or. lisccp) then
     print*,'   Reading in MISR/MODIS/ISCCP subcolumn inputs...'     
     allocate(tau_067(nPoints,nColumns,nLevels))
     status  = nf90_inq_varid(ncid,"TAU_067",varID(15))  
     status  = nf90_get_var(ncid,varID(15),temp4D,count=(/nLon,nLat,nColumns,nLevels/))
     tau_067 = reshape(temp4D,(/nPoints, nColumns, nLevels/))
  endif

  if (lmisr) then
     print*,'   Reading in MISR subcolumn inputs...'
     allocate(emiss_11(nPoints,nColumns,nLevels))
     status   = nf90_inq_varid(ncid,"EMISS_11",varID(16))  
     status   = nf90_get_var(ncid,varID(16),temp4D,count=(/nLon,nLat,nColumns,nLevels/))
     emiss_11 = reshape(temp4D,(/nPoints, nColumns, nLevels/))
  endif
  
  if (lmodis) then
     print*,'   Reading in MODIS subcolumn inputs...'
     allocate(MODIS_fracliq(nPoints,nColumns,nLevels),       &
              MODIS_asym(nPoints,nColumns,nLevels),          &
              MODIS_ssa(nPoints,nColumns,nLevels))
     status        = nf90_inq_varid(ncid,"MODIS_fracliq",varID(17))  
     status        = nf90_get_var(ncid,varID(17),temp4D,count=(/nLon,nLat,nColumns,nLevels/))
     MODIS_fracliq = reshape(temp4D,(/nPoints, nColumns, nLevels/))
     status        = nf90_inq_varid(ncid,"MODIS_asym",varID(18))  
     status        = nf90_get_var(ncid,varID(18),temp4D,count=(/nLon,nLat,nColumns,nLevels/))
     MODIS_asym    = reshape(temp4D,(/nPoints, nColumns, nLevels/))
     status        = nf90_inq_varid(ncid,"MODIS_ssa",varID(19))  
     status        = nf90_get_var(ncid,varID(19),temp4D,count=(/nLon,nLat,nColumns,nLevels/))
     MODIS_ssa     = reshape(temp4D,(/nPoints, nColumns, nLevels/))
  endif

  if (lcalipso) then
     print*,'   Reading in Calipso subcolumn inputs...'
     allocate(calipso_betatot(nPoints,nColumns,nLevels),     &
              calipso_betatot_ice(nPoints,nColumns,nLevels), &
              calipso_betatot_liq(nPoints,nColumns,nLevels), &
              calipso_tautot(nPoints,nColumns,nLevels),      &
              calipso_tautot_ice(nPoints,nColumns,nLevels),  &
              calipso_tautot_liq(nPoints,nColumns,nLevels))
     status              = nf90_inq_varid(ncid,"CAL_betatot",varID(20))  
     status              = nf90_get_var(ncid,varID(20),temp4D,count=(/nLon,nLat,nColumns,nLevels/))
     calipso_betatot     = reshape(temp4D,(/nPoints, nColumns, nLevels/))
     status              = nf90_inq_varid(ncid,"CAL_betatot_ice",varID(21))  
     status              = nf90_get_var(ncid,varID(21),temp4D,count=(/nLon,nLat,nColumns,nLevels/))
     calipso_betatot_ice = reshape(temp4D,(/nPoints, nColumns, nLevels/))
     status              = nf90_inq_varid(ncid,"CAL_betatot_liq",varID(22))  
     status              = nf90_get_var(ncid,varID(22),temp4D,count=(/nLon,nLat,nColumns,nLevels/))
     calipso_betatot_liq = reshape(temp4D,(/nPoints, nColumns, nLevels/))
     status              = nf90_inq_varid(ncid,"CAL_tautot",varID(23))  
     status              = nf90_get_var(ncid,varID(23),temp4D,count=(/nLon,nLat,nColumns,nLevels/))
     calipso_tautot      = reshape(temp4D,(/nPoints, nColumns, nLevels/))
     status              = nf90_inq_varid(ncid,"CAL_tautot_ice",varID(24))  
     status              = nf90_get_var(ncid,varID(24),temp4D,count=(/nLon,nLat,nColumns,nLevels/))
     calipso_tautot_ice  = reshape(temp4D,(/nPoints, nColumns, nLevels/))
     status              = nf90_inq_varid(ncid,"CAL_tautot_liq",varID(25))  
     status              = nf90_get_var(ncid,varID(25),temp4D,count=(/nLon,nLat,nColumns,nLevels/))
     calipso_tautot_liq  = reshape(temp4D,(/nPoints, nColumns, nLevels/))
  endif
  
  if (lcloudsat) then
     print*,'   Reading in Cloudsat subcolumn inputs...'
     allocate(cloudsat_z_vol(nPoints,nColumns,nLevels),      &
              cloudsat_kr_vol(nPoints,nColumns,nLevels),     &
              cloudsat_g_vol(nPoints,nColumns,nLevels))
     status              = nf90_inq_varid(ncid,"CS_z_vol",varID(26))  
     status              = nf90_get_var(ncid,varID(26),temp4D,count=(/nLon,nLat,nColumns,nLevels/))
     cloudsat_z_vol      = reshape(temp4D,(/nPoints, nColumns, nLevels/))
     status              = nf90_inq_varid(ncid,"CS_kr_vol",varID(27))  
     status              = nf90_get_var(ncid,varID(27),temp4D,count=(/nLon,nLat,nColumns,nLevels/))
     cloudsat_kr_vol     = reshape(temp4D,(/nPoints, nColumns, nLevels/))
     status              = nf90_inq_varid(ncid,"CS_g_vol",varID(28))  
     status              = nf90_get_var(ncid,varID(28),temp4D,count=(/nLon,nLat,nColumns,nLevels/))
     cloudsat_g_vol      = reshape(temp4D,(/nPoints, nColumns, nLevels/))
  endif
  call cpu_time(driver_time(2))
       
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
       Latb532gr, Latb355, LlidarBetaMol532, LlidarBetaMol532gr, LlidarBetaMol355,       & 
       LcfadLidarsr532, LcfadLidarsr532gr, LcfadLidarsr355, Lclcalipso2,                 & 
       Lclcalipso, LclgrLidar532, Lclatlid, Lclhcalipso, Lcllcalipso, Lclmcalipso,       & 
       Lcltcalipso, LclhgrLidar532, LcllgrLidar532, LclmgrLidar532, LcltgrLidar532,      & 
       Lclhatlid, Lcllatlid, Lclmatlid, Lcltatlid, Lcltlidarradar,  Lcloudsat_tcc,       &
       Lcloudsat_tcc2, Lclcalipsoliq,                                                    & 
       Lclcalipsoice, Lclcalipsoun, Lclcalipsotmp, Lclcalipsotmpliq, Lclcalipsotmpice,   &
       Lclcalipsotmpun, Lcltcalipsoliq, Lcltcalipsoice, Lcltcalipsoun, Lclhcalipsoliq,   &
       Lclhcalipsoice, Lclhcalipsoun, Lclmcalipsoliq, Lclmcalipsoice, Lclmcalipsoun,     &
       Lcllcalipsoliq, Lcllcalipsoice, Lcllcalipsoun, Lclopaquecalipso, Lclthincalipso,  & 
       Lclzopaquecalipso, Lclcalipsoopaque, Lclcalipsothin, Lclcalipsozopaque,           & 
       Lclcalipsoopacity, Lclopaquetemp, Lclthintemp, Lclzopaquetemp, Lclopaquemeanz,    & 
       Lclthinmeanz, Lclthinemis, Lclopaquemeanzse, Lclthinmeanzse, Lclzopaquecalipsose, &
       LcfadDbze94, Ldbze94, Lparasolrefl,                                               &
       Lptradarflag0,Lptradarflag1,Lptradarflag2,Lptradarflag3,Lptradarflag4,            &
       Lptradarflag5,Lptradarflag6,Lptradarflag7,Lptradarflag8,Lptradarflag9,Lradarpia,  &
       Lwr_occfreq, Lcfodd,                                                              &
       Npoints, Ncolumns, Nlevels, Nlvgrid, cospOUT)

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Initialize COSP
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  call COSP_INIT(Lisccp, Lmodis, Lmisr, Lcloudsat, Lcalipso, LgrLidar532, Latlid,        &
       Lparasol, .false., cloudsat_radar_freq, cloudsat_k2, cloudsat_use_gas_abs,        &
       cloudsat_do_ray, isccp_topheight, isccp_topheight_direction, surface_radar,       &
       rcfg_cloudsat, use_vgrid, csat_vgrid, Nlvgrid, Nlevels, cloudsat_micro_scheme)
  call cpu_time(driver_time(3))

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
     print*,'Chunk ',iChunk,' of ',nChunks
     print*,'      ',start_idx,end_idx
     
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! Construct COSP input types
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     if (iChunk .eq. 1) then
        call construct_cospIN(Nptsperit,nColumns,nLevels,cospIN)
        call construct_cospstateIN(Nptsperit,nLevels,cospstateIN)
     endif
     if (iChunk .eq. nChunks) then
        call destroy_cospIN(cospIN)
        call destroy_cospstateIN(cospstateIN)
        call construct_cospIN(Nptsperit,nColumns,nLevels,cospIN)
        call construct_cospstateIN(Nptsperit,nLevels,cospstateIN)    
     endif
     call cpu_time(driver_time(4))

     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! Populate COSP2 input types with CESM2 model fields.
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
     cospIN%emsfc_lw                = 1
     cospIN%frac_out                = cld_frac
     cospstateIN%hgt_matrix         = zlev(start_idx:end_idx,:)
     cospstateIN%hgt_matrix_half    = zlev_half(start_idx:end_idx,:)
     cospstateIN%sunlit             = sunlit(start_idx:end_idx)
     cospstateIN%skt                = skt(start_idx:end_idx)
     cospstateIN%land               = landmask(start_idx:end_idx)
     cospstateIN%qv                 = q(start_idx:end_idx,:)
     cospstateIN%at                 = T(start_idx:end_idx,:)
     cospstateIN%pfull              = p(start_idx:end_idx,:)
     cospstateIN%phalf(:,1:nLevels) = ph(start_idx:end_idx,:)
     cospstateIN%phalf(:,nLevels)   = 0._wp
     if (lmisr .or. lisccp .or. lmodis) then
        cospIN%tau_067              = tau_067(start_idx:end_idx,:,:)
     endif
     if (lmisr) then
        cospIN%emiss_11             = emiss_11(start_idx:end_idx,:,:)
     endif
     if (lmodis) then
        cospIN%fracLiq              = MODIS_fracliq(start_idx:end_idx,:,:)
        cospIN%asym                 = MODIS_asym(start_idx:end_idx,:,:)
        cospIN%ss_alb               = MODIS_ssa(start_idx:end_idx,:,:)
     endif
     if (lcalipso) then
        cospIN%betatot_calipso      = calipso_betatot(start_idx:end_idx,:,:)
        cospIN%betatot_ice_calipso  = calipso_betatot_ice(start_idx:end_idx,:,:)
        cospIN%betatot_liq_calipso  = calipso_betatot_liq(start_idx:end_idx,:,:)
        cospIN%tautot_calipso       = calipso_tautot(start_idx:end_idx,:,:)
        cospIN%tautot_ice_calipso   = calipso_tautot(start_idx:end_idx,:,:)
        cospIN%tautot_liq_calipso   = calipso_tautot(start_idx:end_idx,:,:)
     endif
     if (lcloudsat) then
        cospIN%rcfg_cloudsat        = rcfg_cloudsat
        cospstateIN%surfelev        = surfelev(start_idx:end_idx)
        cospIN%z_vol_cloudsat       = cloudsat_z_vol(start_idx:end_idx,:,:)
        cospIN%kr_vol_cloudsat      = cloudsat_kr_vol(start_idx:end_idx,:,:)
        cospIN%g_vol_cloudsat       = cloudsat_g_vol(start_idx:end_idx,:,:)
     endif
     call cpu_time(driver_time(5))
 
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! Call COSP
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     cosp_status = COSP_SIMULATOR(cospIN, cospstateIN, cospOUT, start_idx, end_idx,.false.)
     do ij=1,size(cosp_status,1)
        if (cosp_status(ij) .ne. '') print*,trim(cosp_status(ij))
     end do
     
     call cpu_time(driver_time(7))
  enddo
  print*,'Time to read in data:     ',driver_time(2)-driver_time(1)
  print*,'Time to initialize:       ',driver_time(3)-driver_time(2)
  print*,'Time to construct types:  ',driver_time(4)-driver_time(3)
  print*,'Time to copy inputs:      ',driver_time(5)-driver_time(4)
  print*,'Time to run COSP:         ',driver_time(7)-driver_time(6)
  print*,'Total time:               ',driver_time(7)-driver_time(1)

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Output
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  call write_cosp2_output(Npoints, Ncolumns, Nlevels, zlev(1,:), Nlvgrid, lon, lat, cospOUT, foutput)

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
       allocate(y%betatot_calipso(npoints,        ncolumns,nlevels),&
                y%betatot_ice_calipso(npoints,    ncolumns,nlevels),&
                y%betatot_liq_calipso(npoints,    ncolumns,nlevels),&
                y%tautot_calipso(npoints,         ncolumns,nlevels),&
                y%tautot_ice_calipso(npoints,     ncolumns,nlevels),&
                y%tautot_liq_calipso(npoints,     ncolumns,nlevels),&
                y%beta_mol_calipso(npoints,                nlevels),&
                y%tau_mol_calipso(npoints,                 nlevels),&
                y%tautot_S_ice(npoints,   ncolumns        ),&
                y%tautot_S_liq(npoints,   ncolumns        ))
    endif

    if (LgrLidar532) then
       allocate(y%beta_mol_grLidar532(npoints,          nlevels),& 
                y%betatot_grLidar532(npoints,  ncolumns,nlevels),& 
                y%tau_mol_grLidar532(npoints,           nlevels),& 
                y%tautot_grLidar532(npoints,   ncolumns,nlevels)) 
    endif

    if (Latlid) then
       allocate(y%beta_mol_atlid(npoints,             nlevels),& 
                y%betatot_atlid(npoints,     ncolumns,nlevels),& 
                y%tau_mol_atlid(npoints,              nlevels),& 
                y%tautot_atlid(npoints,      ncolumns,nlevels))
    endif 

    if (Lcloudsat) then
       allocate(y%z_vol_cloudsat(npoints,  ncolumns,nlevels),&
                y%kr_vol_cloudsat(npoints, ncolumns,nlevels),&
                y%g_vol_cloudsat(npoints,  ncolumns,nlevels),&
                y%fracPrecipIce(npoints,   ncolumns))
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
  subroutine construct_cospstateIN(npoints,nlevels,y)
    ! Inputs
    integer,intent(in) :: &
         npoints, & ! Number of horizontal gridpoints
         nlevels    ! Number of vertical levels
    ! Outputs
    type(cosp_column_inputs),intent(out) :: y         
    
    allocate(y%sunlit(npoints),y%skt(npoints),y%land(npoints),y%at(npoints,nlevels),     &
             y%pfull(npoints,nlevels),y%phalf(npoints,nlevels+1),y%qv(npoints,nlevels),  &
             y%o3(npoints,nlevels),y%hgt_matrix(npoints,nlevels),y%u_sfc(npoints),       &
             y%v_sfc(npoints),y%lat(npoints),y%lon(nPoints),                             &
             y%cloudIce(nPoints,nLevels),y%cloudLiq(nPoints,nLevels),y%surfelev(npoints),&
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
                                    Liwpmodis,Lclmodis,Latb532,Latb532gr,Latb355,        &
                                    LlidarBetaMol532,LlidarBetaMol532gr,LlidarBetaMol355,&
                                    LcfadLidarsr532,LcfadLidarsr532gr,LcfadLidarsr355,   & 
                                    Lclcalipso2,Lclcalipso,LclgrLidar532,Lclatlid,      & 
                                    Lclhcalipso,Lcllcalipso,Lclmcalipso,Lcltcalipso,     &
                                    LclhgrLidar532,LcllgrLidar532,LclmgrLidar532,     & 
                                    LcltgrLidar532,Lclhatlid,Lcllatlid,Lclmatlid,       & 
                                    Lcltatlid,Lcltlidarradar,Lcloudsat_tcc,            &
                                    Lcloudsat_tcc2,Lclcalipsoliq,              &
                                    Lclcalipsoice,Lclcalipsoun,Lclcalipsotmp,            &
                                    Lclcalipsotmpliq,Lclcalipsotmpice,Lclcalipsotmpun,   &
                                    Lcltcalipsoliq,Lcltcalipsoice,Lcltcalipsoun,         &
                                    Lclhcalipsoliq,Lclhcalipsoice,Lclhcalipsoun,         &
                                    Lclmcalipsoliq,Lclmcalipsoice,Lclmcalipsoun,         &
                                    Lcllcalipsoliq,Lcllcalipsoice,Lcllcalipsoun,         & 
                                    Lclopaquecalipso,Lclthincalipso,Lclzopaquecalipso,   & 
                                    Lclcalipsoopaque,Lclcalipsothin,Lclcalipsozopaque,   & 
                                    Lclcalipsoopacity,Lclopaquetemp,Lclthintemp,         & 
                                    Lclzopaquetemp,Lclopaquemeanz,Lclthinmeanz,          & 
                                    Lclthinemis,Lclopaquemeanzse,Lclthinmeanzse,         &
                                    Lclzopaquecalipsose,LcfadDbze94,Ldbze94,Lparasolrefl,&
                                    Lptradarflag0,Lptradarflag1,Lptradarflag2,           &
                                    Lptradarflag3,Lptradarflag4,Lptradarflag5,           &
                                    Lptradarflag6,Lptradarflag7,Lptradarflag8,           &
                                    Lptradarflag9,Lradarpia,Lwr_occfreq,Lcfodd,          &
                                    Npoints,Ncolumns,Nlevels,Nlvgrid,x)
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
         Latb532gr,        & ! GROUND LIDAR @ 532NM attenuated total backscatter (532nm)
         Latb355,          & ! ATLID attenuated total backscatter (355nm)
         LlidarBetaMol532, & ! CALIPSO molecular backscatter (532nm)         
         LlidarBetaMol532gr,&! GROUND LIDAR @ 532NM molecular backscatter (532nm)
         LlidarBetaMol355, & ! ATLID molecular backscatter (355nm) 
         LcfadLidarsr532,  & ! CALIPSO scattering ratio CFAD
         LcfadLidarsr532gr,& ! GROUND LIDAR @ 532NM scattering ratio CFAD  
         LcfadLidarsr355,  & ! ATLID scattering ratio CFAD 
         Lclcalipso2,      & ! CALIPSO cloud fraction undetected by cloudsat
         Lclcalipso,       & ! CALIPSO cloud area fraction
         LclgrLidar532,   & ! GROUND LIDAR @ 532NM cloud area fraction 
         Lclatlid,         & ! ATLID cloud area fraction 
         Lclhcalipso,      & ! CALIPSO high-level cloud fraction
         Lcllcalipso,      & ! CALIPSO low-level cloud fraction
         Lclmcalipso,      & ! CALIPSO mid-level cloud fraction
         Lcltcalipso,      & ! CALIPSO total cloud fraction
         LclhgrLidar532,  & ! GROUND LIDAR @ 532NM high-level cloud fraction 
         LcllgrLidar532,  & ! GROUND LIDAR @ 532NM low-level cloud fraction 
         LclmgrLidar532,  & ! GROUND LIDAR @ 532NM mid-level cloud fraction
         LcltgrLidar532,  & ! GROUND LIDAR @ 532NM total cloud fraction
         Lclhatlid,        & ! ATLID high-level cloud fraction
         Lcllatlid,        & ! ATLID low-level cloud fraction  
         Lclmatlid,        & ! ATLID mid-level cloud fraction 
         Lcltatlid,        & ! ATLID total cloud fraction
         Lcltlidarradar,   & ! CALIPSO-CLOUDSAT total cloud fraction
         Lcloudsat_tcc,    & !
         Lcloudsat_tcc2,   & !
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
         Lclopaquecalipso, & ! CALIPSO opaque cloud cover (2D Map)
         Lclthincalipso,   & ! CALIPSO thin cloud cover (2D Map)
         Lclzopaquecalipso,& ! CALIPSO z_opaque altitude (opaque clouds only, 2D Map)
         Lclcalipsoopaque, & ! CALIPSO opaque cloud profiles 3D fraction 
         Lclcalipsothin,   & ! CALIPSO thin cloud profiles 3D fraction 
         Lclcalipsozopaque,& ! CALIPSO z_opaque 3D fraction 
         Lclcalipsoopacity,& ! CALIPSO opacity 3D fraction 
         Lclopaquetemp,    & ! CALIPSO opaque cloud temperature 
         Lclthintemp,      & ! CALIPSO thin cloud temperature
         Lclzopaquetemp,   & ! CALIPSO z_opaque temperature  
         Lclopaquemeanz,   & ! CALIPSO opaque cloud altitude  
         Lclthinmeanz,     & ! CALIPSO thin cloud altitude 
         Lclthinemis,      & ! CALIPSO thin cloud emissivity
         Lclopaquemeanzse,   & ! CALIPSO opaque cloud altitude with respect to SE 
         Lclthinmeanzse,     & ! CALIPSO thin cloud altitude with respect to SE
         Lclzopaquecalipsose,& ! CALIPSO z_opaque altitude with respect to SE
         LcfadDbze94,      & ! CLOUDSAT radar reflectivity CFAD
         Ldbze94,          & ! CLOUDSAT radar reflectivity
         LparasolRefl,     & ! PARASOL reflectance
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
         Lwr_occfreq,      & ! CloudSat+MODIS joint diagnostics
         Lcfodd              ! CloudSat+MODIS joint diagnostics
         
     integer,intent(in) :: &
          Npoints,         & ! Number of sampled points
          Ncolumns,        & ! Number of subgrid columns
          Nlevels,         & ! Number of model levels
          Nlvgrid            ! Number of levels in L3 stats computation
          
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
    if (Lclopaquecalipso .or. Lclthincalipso .or. Lclzopaquecalipso) then
        allocate(x%calipso_cldtype(Npoints,LIDAR_NTYPE))
    endif 
    if (Lclopaquetemp .or. Lclthintemp .or. Lclzopaquetemp) then 
        allocate(x%calipso_cldtypetemp(Npoints,LIDAR_NTYPE))  
    endif
    if (Lclopaquemeanz .or. Lclthinmeanz) then 
        allocate(x%calipso_cldtypemeanz(Npoints,2))
    endif 
    if (Lclopaquemeanzse .or. Lclthinmeanzse .or. Lclzopaquecalipsose) then 
        allocate(x%calipso_cldtypemeanzse(Npoints,3)) 
    endif 
    if (Lclthinemis) then 
        allocate(x%calipso_cldthinemis(Npoints))
    endif
    if (Lclcalipsoopaque .or. Lclcalipsothin .or. Lclcalipsozopaque .or. Lclcalipsoopacity) then 
        allocate(x%calipso_lidarcldtype(Npoints,Nlvgrid,LIDAR_NTYPE+1))
    endif
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

    ! GROUND LIDAR @ 532NM simulator
    if (LlidarBetaMol532gr) allocate(x%grLidar532_beta_mol(Npoints,Nlevels))
    if (Latb532gr)          allocate(x%grLidar532_beta_tot(Npoints,Ncolumns,Nlevels))
    if (LcfadLidarsr532gr) then 
        allocate(x%grLidar532_srbval(SR_BINS+1)) 
        allocate(x%grLidar532_cfad_sr(Npoints,SR_BINS,Nlvgrid))
    endif
    if (LclgrLidar532)     allocate(x%grLidar532_lidarcld(Npoints,Nlvgrid)) 
    if (LclhgrLidar532 .or. LclmgrLidar532 .or. LcllgrLidar532 .or. LcltgrLidar532) then
        allocate(x%grLidar532_cldlayer(Npoints,LIDAR_NCAT))  
    endif
      
    ! ATLID simulator
    if (LlidarBetaMol355) allocate(x%atlid_beta_mol(Npoints,Nlevels))
    if (Latb355)          allocate(x%atlid_beta_tot(Npoints,Ncolumns,Nlevels))
    if (LcfadLidarsr355) then
        allocate(x%atlid_srbval(SR_BINS+1)) 
        allocate(x%atlid_cfad_sr(Npoints,SR_BINS,Nlvgrid))
    endif 
    if (Lclatlid)     allocate(x%atlid_lidarcld(Npoints,Nlvgrid)) 
    if (Lclhatlid .or. Lclmatlid .or. Lcllatlid .or. Lcltatlid) then
        allocate(x%atlid_cldlayer(Npoints,LIDAR_NCAT)) 
    endif
      
    ! PARASOL
    if (Lparasolrefl) then
        allocate(x%parasolPix_refl(Npoints,Ncolumns,PARASOL_NREFL))
        allocate(x%parasolGrid_refl(Npoints,PARASOL_NREFL))
    endif 

    ! Cloudsat simulator
    if (Ldbze94)        allocate(x%cloudsat_Ze_tot(Npoints,Ncolumns,Nlevels))
    if (LcfadDbze94)    allocate(x%cloudsat_cfad_ze(Npoints,cloudsat_DBZE_BINS,Nlvgrid))
    if (Lptradarflag0 .or. Lptradarflag1 .or. Lptradarflag2 .or. Lptradarflag3 .or. &
        Lptradarflag4 .or. Lptradarflag5 .or. Lptradarflag6 .or. Lptradarflag7 .or. &
        Lptradarflag8 .or. Lptradarflag9) then
       allocate(x%cloudsat_precip_cover(Npoints,cloudsat_DBZE_BINS))
    endif
    if (Lradarpia) allocate(x%cloudsat_pia(Npoints))

    ! Combined CALIPSO/CLOUDSAT fields
    if (Lclcalipso2)    allocate(x%lidar_only_freq_cloud(Npoints,Nlvgrid))
    if (Lcltlidarradar) allocate(x%cloudsat_calipso_tcc(Npoints))
    if (Lcloudsat_tcc) allocate(x%cloudsat_tcc(Npoints))
    if (Lcloudsat_tcc2) allocate(x%cloudsat_tcc2(Npoints))
            
    ! Joint MODIS/CloudSat Statistics
    if (Lwr_occfreq)  allocate(x%wr_occfreq_ntotal(Npoints,WR_NREGIME))
    if (Lcfodd)       allocate(x%cfodd_ntotal(Npoints,CFODD_NDBZE,CFODD_NICOD,CFODD_NCLASS))

  end subroutine construct_cosp_outputs
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE destroy_cospIN     
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine destroy_cospIN(y)
    type(cosp_optical_inputs),intent(inout) :: y
    if (allocated(y%tau_067))             deallocate(y%tau_067)
    if (allocated(y%emiss_11))            deallocate(y%emiss_11)
    if (allocated(y%frac_out))            deallocate(y%frac_out)
    if (allocated(y%beta_mol_calipso))    deallocate(y%beta_mol_calipso)
    if (allocated(y%tau_mol_calipso))     deallocate(y%tau_mol_calipso)
    if (allocated(y%betatot_calipso))     deallocate(y%betatot_calipso)
    if (allocated(y%betatot_ice_calipso)) deallocate(y%betatot_ice_calipso)
    if (allocated(y%betatot_liq_calipso)) deallocate(y%betatot_liq_calipso)
    if (allocated(y%tautot_calipso))      deallocate(y%tautot_calipso)
    if (allocated(y%tautot_ice_calipso))  deallocate(y%tautot_ice_calipso)
    if (allocated(y%tautot_liq_calipso))  deallocate(y%tautot_liq_calipso)
    if (allocated(y%tautot_S_liq))        deallocate(y%tautot_S_liq)
    if (allocated(y%tautot_S_ice))        deallocate(y%tautot_S_ice)
    if (allocated(y%z_vol_cloudsat))      deallocate(y%z_vol_cloudsat)
    if (allocated(y%kr_vol_cloudsat))     deallocate(y%kr_vol_cloudsat)
    if (allocated(y%g_vol_cloudsat))      deallocate(y%g_vol_cloudsat)
    if (allocated(y%asym))                deallocate(y%asym)
    if (allocated(y%ss_alb))              deallocate(y%ss_alb)
    if (allocated(y%fracLiq))             deallocate(y%fracLiq)
    if (allocated(y%beta_mol_grLidar532)) deallocate(y%beta_mol_grLidar532)
    if (allocated(y%betatot_grLidar532))  deallocate(y%betatot_grLidar532)
    if (allocated(y%tau_mol_grLidar532))  deallocate(y%tau_mol_grLidar532)
    if (allocated(y%tautot_grLidar532))   deallocate(y%tautot_grLidar532)
    if (allocated(y%beta_mol_atlid))      deallocate(y%beta_mol_atlid) 
    if (allocated(y%betatot_atlid))       deallocate(y%betatot_atlid) 
    if (allocated(y%tau_mol_atlid))       deallocate(y%tau_mol_atlid) 
    if (allocated(y%tautot_atlid))        deallocate(y%tautot_atlid)
    if (allocated(y%fracPrecipIce))      deallocate(y%fracPrecipIce)
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
    if (allocated(y%surfelev))        deallocate(y%surfelev)
    
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
     if (associated(y%calipso_lidarcldtype))     then
        deallocate(y%calipso_lidarcldtype) 
        nullify(y%calipso_lidarcldtype) 
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
     if (associated(y%calipso_cldtype))          then
        deallocate(y%calipso_cldtype)
        nullify(y%calipso_cldtype)  
     endif  
     if (associated(y%calipso_cldtypetemp))      then
        deallocate(y%calipso_cldtypetemp) 
        nullify(y%calipso_cldtypetemp) 
     endif 
     if (associated(y%calipso_cldtypemeanz))     then
        deallocate(y%calipso_cldtypemeanz)
        nullify(y%calipso_cldtypemeanz)  
     endif  
     if (associated(y%calipso_cldtypemeanzse))   then 
        deallocate(y%calipso_cldtypemeanzse) 
        nullify(y%calipso_cldtypemeanzse)  
     endif 
     if (associated(y%calipso_cldthinemis))      then
        deallocate(y%calipso_cldthinemis)
        nullify(y%calipso_cldthinemis)
     endif 
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
     if (associated(y%grLidar532_beta_mol))     then
        deallocate(y%grLidar532_beta_mol)
        nullify(y%grLidar532_beta_mol)  
     endif
     if (associated(y%grLidar532_beta_tot))     then 
        deallocate(y%grLidar532_beta_tot) 
        nullify(y%grLidar532_beta_tot)
     endif 
     if (associated(y%grLidar532_cldlayer))     then 
        deallocate(y%grLidar532_cldlayer) 
        nullify(y%grLidar532_cldlayer) 
     endif 
     if (associated(y%grLidar532_lidarcld))     then 
        deallocate(y%grLidar532_lidarcld) 
        nullify(y%grLidar532_lidarcld) 
     endif 
     if (associated(y%grLidar532_cfad_sr))      then 
        deallocate(y%grLidar532_cfad_sr)  
        nullify(y%grLidar532_cfad_sr) 
     endif 
     if (associated(y%grLidar532_srbval))       then
        deallocate(y%grLidar532_srbval) 
        nullify(y%grLidar532_srbval) 
     endif  
     if (associated(y%atlid_beta_mol))           then
        deallocate(y%atlid_beta_mol) 
        nullify(y%atlid_beta_mol) 
     endif 
     if (associated(y%atlid_beta_tot))           then
        deallocate(y%atlid_beta_tot) 
        nullify(y%atlid_beta_tot) 
     endif
     if (associated(y%atlid_cldlayer))           then 
        deallocate(y%atlid_cldlayer)  
        nullify(y%atlid_cldlayer) 
     endif 
     if (associated(y%atlid_lidarcld))           then 
        deallocate(y%atlid_lidarcld)  
        nullify(y%atlid_lidarcld) 
     endif  
     if (associated(y%atlid_cfad_sr))            then
        deallocate(y%atlid_cfad_sr) 
        nullify(y%atlid_cfad_sr) 
     endif
     if (associated(y%atlid_srbval))             then 
        deallocate(y%atlid_srbval) 
        nullify(y%atlid_srbval) 
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
     if (associated(y%cloudsat_precip_cover))     then
        deallocate(y%cloudsat_precip_cover)
        nullify(y%cloudsat_precip_cover)
     endif
     if (associated(y%cloudsat_pia))              then
        deallocate(y%cloudsat_pia)
        nullify(y%cloudsat_pia)
     endif
     if (associated(y%cloudsat_tcc))           then
        deallocate(y%cloudsat_tcc) 
        nullify(y%cloudsat_tcc)  
     endif
     if (associated(y%cloudsat_tcc2))           then
        deallocate(y%cloudsat_tcc2) 
        nullify(y%cloudsat_tcc2)  
     endif
     if (associated(y%cloudsat_calipso_tcc))           then
        deallocate(y%cloudsat_calipso_tcc) 
        nullify(y%cloudsat_calipso_tcc)  
     endif
     if (associated(y%cloudsat_tcc))           then
        deallocate(y%cloudsat_tcc) 
        nullify(y%cloudsat_tcc)  
     endif
     if (associated(y%cloudsat_tcc2))           then
        deallocate(y%cloudsat_tcc2) 
        nullify(y%cloudsat_tcc2)  
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
     if (associated(y%cfodd_ntotal)) then
        deallocate(y%cfodd_ntotal)
        nullify(y%cfodd_ntotal)
     endif
     if (associated(y%wr_occfreq_ntotal)) then
        deallocate(y%wr_occfreq_ntotal)
        nullify(y%wr_occfreq_ntotal)
     endif

   end subroutine destroy_cosp_outputs
  
 end program cosp2_cesm2_test

