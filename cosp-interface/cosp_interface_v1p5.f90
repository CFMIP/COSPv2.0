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
MODULE MOD_COSP_INTERFACE_v1p5
  USE COSP_KINDS,           ONLY: wp
  USE MOD_COSP,             ONLY: cosp_simulator,cosp_optical_inputs,cosp_column_inputs, &
                                  cosp_init,cosp_outputs,                                &
                                  construct_cospstateIN, destroy_cospstateIN,            &
                                  construct_cospIN,      destroy_cospIN
  USE COSP_PHYS_CONSTANTS,  ONLY: amw,amd,amO3,amCO2,amCH4,amN2O,amCO
  USE MOD_COSP_CONFIG,      ONLY: PARASOL_NREFL,RTTOV_MAX_CHANNELS,N_HYDRO,R_UNDEF,      &
                                  COSP_VERSION,numMODISTauBins,modis_histTau,            &
                                  modis_histTauEdges,modis_histTauCenters,ntau,          &
                                  tau_binBounds,tau_binEdges,tau_binCenters
  USE cosp_optics,          ONLY: cosp_simulator_optics,lidar_optics,                &
                                  modis_optics_partition, num_trial_res,modis_optics
  USE MOD_COSP_UTILS,       ONLY: cosp_precip_mxratio
  USE quickbeam,            ONLY: radar_cfg,save_scale_LUTs
  USE mod_quickbeam_optics, ONLY: quickbeam_optics,size_distribution,hydro_class_init
  USE mod_rng,              ONLY: rng_state, init_rng
  USE mod_scops,            ONLY: scops
  USE mod_prec_scops,       ONLY: prec_scops
  
  implicit none
  
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
  
  ! Stratiform and convective clouds in frac_out.
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
       
       
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Initialization variables
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  type(radar_cfg)                          :: rcfg_cloudsat  ! Radar configuration
  type(size_distribution)                  :: sd             ! Hydrometeor description
  integer ::         & 
       overlap,      & !
       Npoints_it,   & ! Number of points per iteration
       lidar_ice_type  ! LIDAR ice type
  logical :: use_precipitation_fluxes
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !TYPE COSP_CONFIG
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  TYPE COSP_CONFIG
     logical :: &
          Lstats,           & ! Control for L3 stats output
          Lwrite_output,    & ! Control for output
          Ltoffset,         & ! Time difference between each profile and the value 
                              ! recorded in varaible time.
          Lradar_sim,       & ! Radar simulator on/off switch 
          Llidar_sim,       & ! LIDAR simulator on/off switch 
          Lisccp_sim,       & ! ISCCP simulator on/off switch
          Lmodis_sim,       & ! MODIS simulatoe on/off switch
          Lmisr_sim,        & ! MISR simulator on/off switch 
          Lrttov_sim,       & ! RTTOV simulator on/off switch 
          Lparasol_sim,     & ! PARASOL simulator on/off switch 
          Lpctisccp,        & ! ISCCP mean cloud top pressure
          Lclisccp,         & ! ISCCP cloud area fraction
          Lboxptopisccp,    & ! ISCCP CTP in each column
          Lboxtauisccp,     & ! ISCCP optical epth in each column
          Ltauisccp,        & ! ISCCP mean optical depth
          Lcltisccp,        & ! ISCCP total cloud fraction
          Lmeantbisccp,     & ! ISCCP mean all-sky 10.5micron brightness temperature
          Lmeantbclrisccp,  & ! ISCCP mean clear-sky 10.5micron brightness temperature
          Lalbisccp,        & ! ISCCP mean cloud albedo
          LcfadDbze94,      & ! CLOUDSAT radar reflectivity CFAD
          Ldbze94,          & ! CLOUDSAT radar reflectivity
          LparasolRefl,     & ! PARASOL reflectance
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
          Llwpmodis,        & ! MODIS cloud ice water path
          Liwpmodis,        & ! MODIS cloud liquid water path
          Lclmodis,         & ! MODIS cloud area fraction
          LclMISR,          & ! MISR cloud fraction
          Lfracout,         & ! SCOPS Subcolumn output
          Ltbrttov            ! RTTOV mean clear-sky brightness temperature
     character(len=32),dimension(:),allocatable :: out_list
  END TYPE COSP_CONFIG
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !                              TYPE COSP_SUBGRID
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  TYPE COSP_SUBGRID
     integer ::      &
          Npoints,   & ! Number of gridpoints
          Ncolumns,  & ! Number of columns
          Nlevels,   & ! Number of levels
          Nhydro       ! Number of hydrometeors
     real(wp),dimension(:,:,:),pointer :: &
          prec_frac, & ! Subgrid precip array (Npoints,Ncolumns,Nlevels)
          frac_out     ! Subgrid cloud array  (Npoints,Ncolumns,Nlevels)
  END TYPE COSP_SUBGRID
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !                              TYPE COSP_GRIDBOX
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  TYPE COSP_GRIDBOX
     integer,pointer :: &
          Npoints,          & ! Number of gridpoints
          Nlevels,          & ! Number of levels
          Ncolumns,         & ! Number of columns
          Nhydro,           & ! Number of hydrometeors
          Nprmts_max_hydro, & ! Max number of parameters for hydrometeor size distribution
          Naero,            & ! Number of aerosol species
          Nprmts_max_aero,  & ! Max number of parameters for aerosol size distributions
          Npoints_it          ! Max number of gridpoints to be processed in one iteration
     
     ! Time [days]
     double precision,pointer :: time
     double precision,dimension(:),pointer :: time_bnds!(2)

     integer,pointer :: &
          nsizes ! Number of discrete drop sizes (um) used to represent the distribution
     
     ! Geolocation and point information (Npoints)
     real(wp),dimension(:),pointer :: &
          toffset,   & ! Time offset of esch point from the value in time 
          longitude, & ! Longitude [degrees East]                       
          latitude,  & ! Latitude [deg North]                          
          land,      & ! Landmask [0 - Ocean, 1 - Land]              
          psfc,      & ! Surface pressure [Pa]                      
          sunlit,    & ! 1 for day points, 0 for nightime            
          skt,       & ! Skin temperature (K)                      
          u_wind,    & ! Eastward wind [m s-1]                   
          v_wind       ! Northward wind [m s-1]                    
     
     ! Gridbox information (Npoints,Nlevels)
     real(wp),dimension(:,:),pointer :: &
          zlev,      & ! Height of model levels [m]                           
          zlev_half, & ! Height at half model levels [m] (Bottom of layer)   
          dlev,      & ! Depth of model levels  [m]                         
          p,         & ! Pressure at full model levels [Pa]      
          ph,        & ! Pressure at half model levels [Pa]             
          T,         & ! Temperature at model levels [K]                 
          sh,        & ! Specific humidity to water [kg/kg]             
          dtau_s,    & ! mean 0.67 micron optical depth of stratiform clouds  
          dtau_c,    & ! mean 0.67 micron optical depth of convective clouds 
          dem_s,     & ! 10.5 micron longwave emissivity of stratiform clouds 
          dem_c,     & ! 10.5 micron longwave emissivity of convective clouds 
          mr_ozone     ! Ozone mass mixing ratio [kg/kg]                    
     
     ! TOTAL and CONV cloud fraction for SCOPS
     real(wp),dimension(:,:),pointer :: &
          tca,       & ! Total cloud fraction
          cca          ! Convective cloud fraction
     
     ! Precipitation fluxes on model levels
     real(wp),dimension(:,:),pointer :: &
          rain_ls,   & ! Large-scale precipitation flux of rain [kg/m2.s]
          rain_cv,   & ! Convective precipitation flux of rain [kg/m2.s]
          snow_ls,   & ! Large-scale precipitation flux of snow [kg/m2.s]
          snow_cv,   & ! Convective precipitation flux of snow [kg/m2.s]
          grpl_ls      ! large-scale precipitation flux of graupel [kg/m2.s]
     
     ! Hydrometeors concentration and distribution parameters
     real(wp),dimension(:,:,:),allocatable :: &
          mr_hydro         ! Mixing ratio of each hydrometeor 
                           ! (Npoints,Nlevels,Nhydro) [kg/kg]
     real(wp),dimension(:,:),pointer :: &
          dist_prmts_hydro ! Distributional parameters for hydrometeors 
                           ! (Nprmts_max_hydro,Nhydro)
     real(wp),dimension(:,:,:),allocatable :: &
          Reff             ! Effective radius [m]. 
                           ! (Npoints,Nlevels,Nhydro)
     real(wp),dimension(:,:,:),allocatable :: &
          Np               ! Total Number Concentration [#/kg]. 
                           ! (Npoints,Nlevels,Nhydro)
 
     ! Aerosols concentration and distribution parameters
     real(wp),dimension(:,:,:),pointer :: &
          conc_aero       ! Aerosol concentration for each species 
                          ! (Npoints,Nlevels,Naero)
     integer,dimension(:),pointer :: &
          dist_type_aero  ! Particle size distribution type for each aerosol species 
                          ! (Naero)
     real(wp),dimension(:,:,:,:),pointer :: &
          dist_prmts_aero ! Distributional parameters for aerosols 
                          ! (Npoints,Nlevels,Nprmts_max_aero,Naero)

     ! ISCCP simulator inputs
     real(wp),pointer :: &
          isccp_emsfc_lw      ! 10.5 micron emissivity of surface (fraction)
     
     ! RTTOV inputs/options
     integer,pointer :: &
          Nchan     ! Number of channels to be computed
     real(wp),dimension(:), pointer :: &
          Surfem    ! Surface emissivity
     real(wp),pointer :: &
          ZenAng, & ! Satellite Zenith Angles
          co2,    & ! CO2 mixing ratio
          ch4,    & ! CH4 mixing ratio
          n2o,    & ! N2O mixing ratio
          co        ! CO mixing ratio
  END TYPE COSP_GRIDBOX
  
contains
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !                      SUBROUTINE cosp_interface_init
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine cosp_interface_init(Npoints,Nlevels,Npointsit,cloud_overlap,lusePrecip,     &
                                 cloudsat_radar_freq,cloudsat_micro_scheme,              &
                                 cloudsat_k2,cloudsat_use_gas_abs,                       &
                                 cloudsat_do_ray,isccp_top_height,                       &
                                 isccp_top_height_direction,hgt_matrix,hgt_matrix_half,  &
                                 surface_radar,rttov_Nchannels,rttov_Channels,           &
                                 rttov_platform,rttov_satellite,rttov_instrument,        &
                                 lidarIceType,lusevgrid,Nvgrid,luseCSATvgrid,cospvID)
    
    ! Inputs
    integer,intent(in) ::            & !
         Npoints,                    & ! Number of horizontal points
         Nlevels,                    & ! Number of vertical points
         Npointsit,                  & ! Number of points per iteration
         Nvgrid,                     & ! Number of vertical layers grid (used by CALIPSO/CLOUDSAT)
         cloud_overlap,              & ! Cloud overlap type
         cloudsat_use_gas_abs,       & !
         cloudsat_do_ray,            & !
         isccp_top_height,           & ! ISCCP cloud top height convention
         isccp_top_height_direction, & ! ISCCP cloud top height direction
         lidarIceType,               & ! LIDAR ice type
         surface_radar,              & !
         rttov_Nchannels,            & ! RTTOV: Number of channels
         rttov_platform,             & ! RTTOV: Satellite platform
         rttov_satellite,            & ! RTTOV: Satellite 
         rttov_instrument              ! RTTOV: Instrument
    integer,intent(in),dimension(RTTOV_MAX_CHANNELS) :: &
         rttov_Channels                ! RTTOV: Channel numbers         
    real(wp),intent(in) ::           & !
         cloudsat_radar_freq,        & ! Cloudsat radar frequency
         cloudsat_k2                   !
    real(wp),intent(in),dimension(Npoints,Nlevels) :: & 
         hgt_matrix     
    real(wp),intent(in),dimension(Npoints,Nlevels+1) :: &
         hgt_matrix_half     
    logical,intent(in) :: &
         lusePrecip,      & ! True if precipitation fluxes are input to the algorithm
         lusevgrid,       & ! True if using new grid for L3 CALIPSO and CLOUDSAT 
         luseCSATvgrid      ! True to use CLOUDSAT vertical grid spacing of 480m
    character(len=64),intent(in) :: &
         cloudsat_micro_scheme ! Microphysical scheme used in radar simulator
    character(len=32),intent(in) :: &
         cospvID

    ! Local variables
    logical :: &
       lsingle=.true., & ! True if using MMF_v3_single_moment CLOUDSAT microphysical scheme (default)
       ldouble=.false.    ! True if using MMF_v3.5_two_moment CLOUDSAT microphysical scheme    
       
    ! Optical depth axis used in cospV1.4.0 joint-histogram
    integer :: l,k
    integer,parameter :: &
       ntauV1p4 = 6
    real(wp),parameter,dimension(ntauV1p4+1) :: &
       tau_binBoundsV1p4 = (/0.3, 1.3, 3.6, 9.4, 23., 60., 10000./)
    real(wp),parameter,dimension(2,ntauV1p4) :: &
       tau_binEdgesV1p4 = reshape(source =(/tau_binBoundsV1p4(1),((tau_binBoundsV1p4(k),l=1,2),   &
                                         k=2,ntauV1p4),100000._wp/),shape = (/2,ntauV1p4/)) 
    real(wp),parameter,dimension(ntauV1p4) :: &
       tau_binCentersV1p4 = (tau_binEdgesV1p4(1,:)+tau_binEdgesV1p4(2,:))/2._wp  
     
    ! Initialization  
    overlap                  = cloud_overlap
    use_precipitation_fluxes = lusePrecip
    lidar_ice_type           = lidarIceType
    Npoints_it               = Npointsit
    COSP_VERSION             = cospvID
    
    ! If COSPv1.4.0 is being run we use a different optical depth axis for the joint-histograms
    if (COSP_VERSION == 'COSP v1.4') then
       print*,'Using cospV1.4.0 optical depth axis for MODIS tau/ctp joint-histogram'
       allocate(modis_histTau(ntauV1p4+1),modis_histTauEdges(2,ntauV1p4),                &
                modis_histTauCenters(ntauV1p4))
       numMODISTauBins      = ntauV1p4+1            ! Number of tau bins for joint-histogram
!ds       numMODISTauBins      = ntauV1p4            ! Number of tau bins for joint-histogram
       modis_histTau        = tau_binBoundsV1p4   ! Joint-histogram boundaries (optical depth)
       modis_histTauEdges   = tau_binEdgesV1p4    ! Joint-histogram bin edges (optical depth)
       modis_histTauCenters = tau_binCentersV1p4  ! Joint-histogram bin centers (optical depth)       
    else
       allocate(modis_histTau(ntau),modis_histTauEdges(2,ntau),                        &
!ds       allocate(modis_histTau(ntau+1),modis_histTauEdges(2,ntau),                        &
                modis_histTauCenters(ntau))
       numMODISTauBins      = ntau
       modis_histTau        = tau_binBounds
       modis_histTauEdges   = tau_binEdges
       modis_histTauCenters = tau_binCenters
    endif   
    
    ! If two-moment radar microphysics scheme is wanted
    if (cloudsat_micro_scheme == 'MMF_v3.5_two_moment')  then
       ldouble = .true. 
       lsingle = .false.
    endif

    ! Initialize the distributional parameters for hydrometeors in radar simulator
    call hydro_class_init(R_UNDEF,lsingle,ldouble,sd)

    ! Initialize COSP simulators
    call COSP_INIT(Npoints,Nlevels,cloudsat_radar_freq,cloudsat_k2,cloudsat_use_gas_abs,  &
                   cloudsat_do_ray,isccp_top_height,isccp_top_height_direction,           &
                   hgt_matrix,hgt_matrix_half,surface_radar,rcfg_cloudsat,rttov_Nchannels,&
                   rttov_Channels,rttov_platform,rttov_satellite,rttov_instrument,        &
                   lusevgrid,luseCSATvgrid,Nvgrid,cloudsat_micro_scheme)
  end subroutine cosp_interface_init
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !                            SUBROUTINE COSP_INTERFACE (v1.5)
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine cosp_interface_v1p5(Npoints,gbx,sgx,cospOUT)
    
    ! Inputs
    type(cosp_gridbox), intent(in)    :: gbx        ! Model grid description
    type(cosp_subgrid), intent(inout) :: sgx        ! Subgrid description    
    integer,            intent(in)    :: Npoints    ! Number of points to process
    
    ! Outputs
    type(cosp_outputs), intent(OUT) :: cospOUT   ! Nested COSP simulator output
    
    ! Local variables
    type(cosp_optical_inputs) :: &
         cospIN        ! Optical (or derived?) fields needed by simulators
    type(cosp_column_inputs) :: &
         cospstateIN   ! Model fields needed by simulators
    integer :: &
         num_chunks, & ! Number of iterations to make
         start_idx,  & ! Starting index when looping over points
         end_idx,    & ! Ending index when looping over points
         Nptsperit     ! Number of points for current iteration
    integer :: i

    ! Number of calls to COSP to make
    num_chunks = Npoints/Npoints_it+1
    print*,'Number of COSP iterations:',num_chunks
    
    ! Figure out which simulators to run given what output fields are allocated
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Loop over each chunk
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    DO i = 1, num_chunks    
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
       ! Compute indices to process
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       if (num_chunks .eq. 1) then
          start_idx = 1
          end_idx   = Npoints
          Nptsperit = Npoints
       else
          start_idx = (i-1)*Npoints_it+1
          end_idx   = i*Npoints_it
          if (end_idx .gt. Npoints) end_idx=Npoints
          Nptsperit = end_idx-start_idx+1
       endif
       print*,'COSP iteration number ',i,'    Processing ',start_idx,' to ',end_idx
       
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Allocate space
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       if (i .eq. 1) then
          call construct_cospIN(Nptsperit,gbx%ncolumns,gbx%nlevels,cospIN)
          call construct_cospstateIN(Nptsperit,gbx%nlevels,gbx%nchan,cospstateIN)
       endif
       if (i .eq. num_chunks) then
          call destroy_cospIN(cospIN)
          call destroy_cospstateIN(cospstateIN)
          call construct_cospIN(Nptsperit,gbx%ncolumns,gbx%nlevels,cospIN)
          call construct_cospstateIN(Nptsperit,gbx%nlevels,gbx%nchan,cospstateIN)    
       endif
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
       ! Generate subcolumns and compute optical properties needed by simulators  
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       call subsample_and_optics_v1p4(gbx,sgx,Nptsperit,start_idx,end_idx,         &
                                      cospIN,cospstateIN)
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! COSP engine 
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       CALL COSP_SIMULATOR(cospIN, cospstateIN, cospOUT, start_idx,end_idx)   
    END DO
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Deallocate space
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    call destroy_cospIN(cospIN)
    call destroy_cospstateIN(cospstateIN)

  end subroutine cosp_interface_v1p5
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !                      SUBROUTINE SUBSAMPLE_AND_OPTICS (cosp1.4)
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine subsample_and_optics_v1p4(gbx,sgx,npoints,start_idx,end_idx,cospIN,cospgridIN)

    ! Inputs
    type(cosp_gridbox),intent(in)    :: gbx   ! Grid box description
    type(cosp_subgrid),intent(inout) :: sgx   ! Sub-grid scale description
    integer,intent(in) :: &
         npoints,     & ! Number of points
         start_idx,   & ! Starting index for subsetting input data.
         end_idx        ! Ending index for subsetting input data.
    ! Outputs
    type(cosp_optical_inputs),intent(inout) :: &
         cospIN         ! Optical (or derived) fields needed by simulators
    type(cosp_column_inputs),intent(inout) :: &
         cospgridIN     ! Model fields needed by simulators
    
    ! Local variables
    integer :: i,j,k,ij
    real(wp),dimension(npoints,gbx%Nlevels) :: column_frac_out,column_prec_out
    real(wp),dimension(:,:),    allocatable :: frac_ls,frac_cv,prec_ls,prec_cv,ls_p_rate,&
                                               cv_p_rate
    real(wp),dimension(:,:,:),allocatable :: frac_out,frac_prec,hm_matrix,re_matrix,     &
                                             Np_matrix,MODIS_cloudWater,MODIS_cloudIce,  &
                                             MODIS_watersize,MODIS_iceSize,              &
                                             MODIS_opticalThicknessLiq,                  &
                                             MODIS_opticalThicknessIce
    real(wp),dimension(:,:,:,:),allocatable :: mr_hydro,Reff,Np
    type(rng_state),allocatable,dimension(:) :: rngs  ! Seeds for random number generator
    integer,dimension(:),allocatable :: seed
    logical :: cmpGases=.true.
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Initialize COSP inputs
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cospIN%tautot_S_liq                 = 0._wp
    cospIN%tautot_S_ice                 = 0._wp
    cospIN%emsfc_lw                     = gbx%isccp_emsfc_lw
    cospIN%frac_out                     = sgx%frac_out(start_idx:end_idx,:,gbx%Nlevels:1:-1)
    cospIN%rcfg_cloudsat                = rcfg_cloudsat
    cospgridIN%hgt_matrix               = gbx%zlev(start_idx:end_idx,gbx%Nlevels:1:-1)
    cospgridIN%hgt_matrix_half          = gbx%zlev_half(start_idx:end_idx,gbx%Nlevels:1:-1)
    cospgridIN%sunlit                   = gbx%sunlit(start_idx:end_idx)
    cospgridIN%skt                      = gbx%skt(start_idx:end_idx)
    cospgridIN%land                     = gbx%land(start_idx:end_idx)
    cospgridIN%qv                       = gbx%sh(start_idx:end_idx,gbx%Nlevels:1:-1) 
    cospgridIN%at                       = gbx%T(start_idx:end_idx,gbx%Nlevels:1:-1) 
    cospgridIN%pfull                    = gbx%p(start_idx:end_idx,gbx%Nlevels:1:-1) 
    cospgridIN%o3                       = gbx%mr_ozone(start_idx:end_idx,gbx%Nlevels:1:-1)*(amd/amO3)*1e6
    cospgridIN%u_sfc                    = gbx%u_wind(start_idx:end_idx)
    cospgridIN%v_sfc                    = gbx%v_wind(start_idx:end_idx)
    cospgridIN%t_sfc                    = gbx%t(start_idx:end_idx,1)
    cospgridIN%emis_sfc                 = gbx%surfem
    cospgridIN%lat                      = gbx%latitude(start_idx:end_idx)
    cospgridIN%co2                      = gbx%co2*(amd/amCO2)*1e6
    cospgridIN%ch4                      = gbx%ch4*(amd/amCH4)*1e6  
    cospgridIN%n2o                      = gbx%n2o*(amd/amN2O)*1e6
    cospgridIN%co                       = gbx%co*(amd/amCO)*1e6
    cospgridIN%zenang                   = gbx%zenang
    cospgridIN%phalf(:,1)               = 0._wp
    cospgridIN%phalf(:,2:gbx%Nlevels+1) = gbx%ph(start_idx:end_idx,gbx%Nlevels:1:-1)    
    
    if (gbx%Ncolumns .gt. 1) then
       
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Random number generator
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       allocate(rngs(Npoints),seed(Npoints))
       seed(:)=0
       seed = int(gbx%psfc)  ! In case of Npoints=1
       if (Npoints .gt. 1) seed=int((gbx%psfc(start_idx:end_idx)-minval(gbx%psfc))/      &
            (maxval(gbx%psfc)-minval(gbx%psfc))*100000) + 1
       call init_rng(rngs, seed)  

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Generate subcolumns for clouds (SCOPS) and precipitation type (PREC_SCOPS)
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Call SCOPS
       if (gbx%Ncolumns .gt. 1) then
          call scops(npoints,gbx%Nlevels,gbx%Ncolumns,seed,rngs,                         &
                     gbx%tca(start_idx:end_idx,gbx%Nlevels:1:-1),                        &
                     gbx%cca(start_idx:end_idx,gbx%Nlevels:1:-1),overlap,                &
                     sgx%frac_out(start_idx:end_idx,:,:),0)
          deallocate(seed,rngs)
       else
          sgx%frac_out(start_idx:end_idx,:,:) = 1  
       endif
       
       ! Sum up precipitation rates
       allocate(ls_p_rate(npoints,gbx%Nlevels),cv_p_rate(npoints,gbx%Nlevels))
       if(use_precipitation_fluxes) then
          ls_p_rate(:,gbx%Nlevels:1:-1) = gbx%rain_ls(start_idx:end_idx,1:gbx%Nlevels) + &
               gbx%snow_ls(start_idx:end_idx,1:gbx%Nlevels) + &
               gbx%grpl_ls(start_idx:end_idx,1:gbx%Nlevels)
          cv_p_rate(:,gbx%Nlevels:1:-1) = gbx%rain_cv(start_idx:end_idx,1:gbx%Nlevels) + &
               gbx%snow_cv(start_idx:end_idx,1:gbx%Nlevels)
       else
          ls_p_rate(:,gbx%Nlevels:1:-1) = &
               gbx%mr_hydro(start_idx:end_idx,1:gbx%Nlevels,I_LSRAIN) +                  &
               gbx%mr_hydro(start_idx:end_idx,1:gbx%Nlevels,I_LSSNOW) +                  &
               gbx%mr_hydro(start_idx:end_idx,1:gbx%Nlevels,I_LSGRPL)
          cv_p_rate(:,gbx%Nlevels:1:-1) =                                                &
               gbx%mr_hydro(start_idx:end_idx,1:gbx%Nlevels,I_CVRAIN) +                  &
               gbx%mr_hydro(start_idx:end_idx,1:gbx%Nlevels,I_CVSNOW)
       endif
       
       ! Call PREC_SCOPS
       call prec_scops(npoints,gbx%Nlevels,gbx%Ncolumns,ls_p_rate,cv_p_rate,             &
                       sgx%frac_out(start_idx:end_idx,:,:),                              &
                       sgx%prec_frac(start_idx:end_idx,:,:))
       deallocate(ls_p_rate,cv_p_rate)

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Compute precipitation fraction in each gridbox
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Allocate
       allocate(frac_ls(npoints,gbx%Nlevels),prec_ls(npoints,gbx%Nlevels),               &
                frac_cv(npoints,gbx%Nlevels),prec_cv(npoints,gbx%Nlevels))

       ! Initialize
       frac_ls(1:npoints,1:gbx%Nlevels) = 0._wp
       prec_ls(1:npoints,1:gbx%Nlevels) = 0._wp
       frac_cv(1:npoints,1:gbx%Nlevels) = 0._wp
       prec_cv(1:npoints,1:gbx%Nlevels) = 0._wp
       do j=1,npoints,1
          do k=1,gbx%Nlevels,1
             do i=1,gbx%Ncolumns,1
                if (sgx%frac_out(start_idx+j-1,i,gbx%Nlevels+1-k) == I_LSC)              &
                     frac_ls(j,k) = frac_ls(j,k)+1._wp
                if (sgx%frac_out(start_idx+j-1,i,gbx%Nlevels+1-k) == I_CVC)              &
                     frac_cv(j,k) = frac_cv(j,k)+1._wp
                if (sgx%prec_frac(start_idx+j-1,i,gbx%Nlevels+1-k) .eq. 1)               &
                     prec_ls(j,k) = prec_ls(j,k)+1._wp
                if (sgx%prec_frac(start_idx+j-1,i,gbx%Nlevels+1-k) .eq. 2)               &
                     prec_cv(j,k) = prec_cv(j,k)+1._wp
                if (sgx%prec_frac(start_idx+j-1,i,gbx%Nlevels+1-k) .eq. 3)               &
                     prec_cv(j,k) = prec_cv(j,k)+1._wp
                if (sgx%prec_frac(start_idx+j-1,i,gbx%Nlevels+1-k) .eq. 3)               &
                     prec_ls(j,k) = prec_ls(j,k)+1._wp
             enddo
             frac_ls(j,k)=frac_ls(j,k)/gbx%Ncolumns
             frac_cv(j,k)=frac_cv(j,k)/gbx%Ncolumns
             prec_ls(j,k)=prec_ls(j,k)/gbx%Ncolumns
             prec_cv(j,k)=prec_cv(j,k)/gbx%Ncolumns
          enddo
       enddo

       ! Flip SCOPS output from TOA-to-SFC to SFC-to-TOA
       sgx%frac_out(start_idx:end_idx,:,1:gbx%Nlevels)  =                                &
            sgx%frac_out(start_idx:end_idx,:,gbx%Nlevels:1:-1)
       sgx%prec_frac(start_idx:end_idx,:,1:gbx%Nlevels) =                                &
            sgx%prec_frac(start_idx:end_idx,:,gbx%Nlevels:1:-1)
       
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Compute mixing ratios, effective radii and precipitation fluxes for clouds
       ! and precipitation
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       allocate(mr_hydro(npoints, gbx%Ncolumns, gbx%Nlevels, gbx%Nhydro),                &
                Reff(    npoints, gbx%Ncolumns, gbx%Nlevels, gbx%Nhydro),                &
                Np(      npoints, gbx%Ncolumns, gbx%Nlevels, gbx%Nhydro))
       mr_hydro(:,:,:,:) = 0._wp
       Reff(:,:,:,:)     = 0._wp
       Np(:,:,:,:)       = 0._wp
       do k=1,gbx%Ncolumns
          ! Subcolumn cloud fraction
          column_frac_out = sgx%frac_out(start_idx:end_idx,k,:)

          ! LS clouds
          where (column_frac_out == I_LSC)
             mr_hydro(:,k,:,I_LSCLIQ) = gbx%mr_hydro(start_idx:end_idx,:,I_LSCLIQ)
             mr_hydro(:,k,:,I_LSCICE) = gbx%mr_hydro(start_idx:end_idx,:,I_LSCICE)
             Reff(:,k,:,I_LSCLIQ)     = gbx%Reff(start_idx:end_idx,:,I_LSCLIQ)
             Reff(:,k,:,I_LSCICE)     = gbx%Reff(start_idx:end_idx,:,I_LSCICE)
             Np(:,k,:,I_LSCLIQ)       = gbx%Np(start_idx:end_idx,:,I_LSCLIQ)
             Np(:,k,:,I_LSCICE)       = gbx%Np(start_idx:end_idx,:,I_LSCICE)
             ! CONV clouds   
          elsewhere (column_frac_out == I_CVC)
             mr_hydro(:,k,:,I_CVCLIQ) = gbx%mr_hydro(start_idx:end_idx,:,I_CVCLIQ)
             mr_hydro(:,k,:,I_CVCICE) = gbx%mr_hydro(start_idx:end_idx,:,I_CVCICE)
             Reff(:,k,:,I_CVCLIQ)     = gbx%Reff(start_idx:end_idx,:,I_CVCLIQ)
             Reff(:,k,:,I_CVCICE)     = gbx%Reff(start_idx:end_idx,:,I_CVCICE)
             Np(:,k,:,I_CVCLIQ)       = gbx%Np(start_idx:end_idx,:,I_CVCLIQ)
             Np(:,k,:,I_CVCICE)       = gbx%Np(start_idx:end_idx,:,I_CVCICE)
          end where
          
          ! Subcolumn precipitation
          column_prec_out = sgx%prec_frac(start_idx:end_idx,k,:)
          
          ! LS Precipitation
          where ((column_prec_out == 1) .or. (column_prec_out == 3) )
             Reff(:,k,:,I_LSRAIN) = gbx%Reff(start_idx:end_idx,:,I_LSRAIN)
             Reff(:,k,:,I_LSSNOW) = gbx%Reff(start_idx:end_idx,:,I_LSSNOW)
             Reff(:,k,:,I_LSGRPL) = gbx%Reff(start_idx:end_idx,:,I_LSGRPL)
             Np(:,k,:,I_LSRAIN)   = gbx%Np(start_idx:end_idx,:,I_LSRAIN)
             Np(:,k,:,I_LSSNOW)   = gbx%Np(start_idx:end_idx,:,I_LSSNOW)
             Np(:,k,:,I_LSGRPL)   = gbx%Np(start_idx:end_idx,:,I_LSGRPL)
          ! CONV precipitation   
          elsewhere ((column_prec_out == 2) .or. (column_prec_out == 3))
             Reff(:,k,:,I_CVRAIN) = gbx%Reff(start_idx:end_idx,:,I_CVRAIN)
             Reff(:,k,:,I_CVSNOW) = gbx%Reff(start_idx:end_idx,:,I_CVSNOW)
             Np(:,k,:,I_CVRAIN)   = gbx%Np(start_idx:end_idx,:,I_CVRAIN)
             Np(:,k,:,I_CVSNOW)   = gbx%Np(start_idx:end_idx,:,I_CVSNOW)
          end where
       enddo
       
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Convert the mixing ratio and precipitation fluxes from gridbox mean to
       ! the fraction-based values
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       do k=1,gbx%Nlevels
          do j=1,npoints
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
                   gbx%rain_ls(start_idx+j-1,k) = gbx%rain_ls(start_idx+j-1,k)/prec_ls(j,k)
                   gbx%snow_ls(start_idx+j-1,k) = gbx%snow_ls(start_idx+j-1,k)/prec_ls(j,k)
                   gbx%grpl_ls(start_idx+j-1,k) = gbx%grpl_ls(start_idx+j-1,k)/prec_ls(j,k)
                endif
                if (prec_cv(j,k) .ne. 0.) then
                   gbx%rain_cv(start_idx+j-1,k) = gbx%rain_cv(start_idx+j-1,k)/prec_cv(j,k)
                   gbx%snow_cv(start_idx+j-1,k) = gbx%snow_cv(start_idx+j-1,k)/prec_cv(j,k)
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
          call cosp_precip_mxratio(npoints, gbx%Nlevels, gbx%Ncolumns,                   &
                                   gbx%p(start_idx:end_idx,:),gbx%T(start_idx:end_idx,:),&
                                   sgx%prec_frac(start_idx:end_idx,:,:), 1._wp,          &
                                   n_ax(I_LSRAIN), n_bx(I_LSRAIN),   alpha_x(I_LSRAIN),  &
                                   c_x(I_LSRAIN),   d_x(I_LSRAIN),   g_x(I_LSRAIN),      &
                                   a_x(I_LSRAIN),   b_x(I_LSRAIN),   gamma_1(I_LSRAIN),  &
                                   gamma_2(I_LSRAIN),gamma_3(I_LSRAIN),gamma_4(I_LSRAIN),&
                                   gbx%rain_ls(start_idx:end_idx,:),                     &
                                   mr_hydro(:,:,:,I_LSRAIN),Reff(:,:,:,I_LSRAIN))
          call cosp_precip_mxratio(npoints, gbx%Nlevels, gbx%Ncolumns,                   &
                                   gbx%p(start_idx:end_idx,:),gbx%T(start_idx:end_idx,:),&
                                   sgx%prec_frac(start_idx:end_idx,:,:), 1._wp,          &          
                                   n_ax(I_LSSNOW),  n_bx(I_LSSNOW),  alpha_x(I_LSSNOW),  &
                                   c_x(I_LSSNOW),   d_x(I_LSSNOW),   g_x(I_LSSNOW),      &
                                   a_x(I_LSSNOW),   b_x(I_LSSNOW),   gamma_1(I_LSSNOW),  &
                                   gamma_2(I_LSSNOW),gamma_3(I_LSSNOW),gamma_4(I_LSSNOW),&
                                   gbx%snow_ls(start_idx:end_idx,:),                     &
                                   mr_hydro(:,:,:,I_LSSNOW),Reff(:,:,:,I_LSSNOW))
          call cosp_precip_mxratio(npoints, gbx%Nlevels, gbx%Ncolumns,                   &
                                   gbx%p(start_idx:end_idx,:),gbx%T(start_idx:end_idx,:),&
                                   sgx%prec_frac(start_idx:end_idx,:,:), 2._wp,          &
                                   n_ax(I_CVRAIN),  n_bx(I_CVRAIN),  alpha_x(I_CVRAIN),  &
                                   c_x(I_CVRAIN),   d_x(I_CVRAIN),   g_x(I_CVRAIN),      &
                                   a_x(I_CVRAIN),   b_x(I_CVRAIN),   gamma_1(I_CVRAIN),  &
                                   gamma_2(I_CVRAIN),gamma_3(I_CVRAIN),gamma_4(I_CVRAIN),&
                                   gbx%rain_cv(start_idx:end_idx,:),                     &
                                   mr_hydro(:,:,:,I_CVRAIN),Reff(:,:,:,I_CVRAIN))
          call cosp_precip_mxratio(npoints, gbx%Nlevels, gbx%Ncolumns,                   &
                                   gbx%p(start_idx:end_idx,:),gbx%T(start_idx:end_idx,:),&
                                   sgx%prec_frac(start_idx:end_idx,:,:), 2._wp,          &          
                                   n_ax(I_CVSNOW),  n_bx(I_CVSNOW),  alpha_x(I_CVSNOW),  &
                                   c_x(I_CVSNOW),   d_x(I_CVSNOW),   g_x(I_CVSNOW),      &
                                   a_x(I_CVSNOW),   b_x(I_CVSNOW),   gamma_1(I_CVSNOW),  &
                                   gamma_2(I_CVSNOW),gamma_3(I_CVSNOW),gamma_4(I_CVSNOW),&
                                   gbx%snow_cv(start_idx:end_idx,:),                     &
                                   mr_hydro(:,:,:,I_CVSNOW),Reff(:,:,:,I_CVSNOW))
          call cosp_precip_mxratio(npoints, gbx%Nlevels, gbx%Ncolumns,                   &
                                   gbx%p(start_idx:end_idx,:),gbx%T(start_idx:end_idx,:),&
                                   sgx%prec_frac(start_idx:end_idx,:,:), 1._wp,          &         
                                   n_ax(I_LSGRPL),  n_bx(I_LSGRPL),  alpha_x(I_LSGRPL),  &
                                   c_x(I_LSGRPL),   d_x(I_LSGRPL),   g_x(I_LSGRPL),      &
                                   a_x(I_LSGRPL),   b_x(I_LSGRPL),   gamma_1(I_LSGRPL),  &
                                   gamma_2(I_LSGRPL),gamma_3(I_LSGRPL),gamma_4(I_LSGRPL),&
                                   gbx%grpl_ls(start_idx:end_idx,:),                     &
                                   mr_hydro(:,:,:,I_LSGRPL),Reff(:,:,:,I_LSGRPL))
       endif
    else
       allocate(mr_hydro(npoints, 1, gbx%Nlevels, gbx%Nhydro),                           &
                Reff(npoints,     1, gbx%Nlevels, gbx%Nhydro),                           &
                Np(npoints,       1, gbx%Nlevels, gbx%Nhydro))
       mr_hydro(:,1,:,:) = gbx%mr_hydro(start_idx:end_idx,:,:)
       Reff(:,1,:,:)     = gbx%Reff(start_idx:end_idx,:,:)
       Np(:,1,:,:)       = gbx%Np(start_idx:end_idx,:,:)
       where(gbx%dtau_s(start_idx:end_idx,:) .gt. 0)
          sgx%frac_out(start_idx:end_idx,1,:) = 1
       endwhere
    endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 11 micron emissivity
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    call cosp_simulator_optics(npoints,gbx%Ncolumns,gbx%Nlevels,                         &
                               sgx%frac_out(start_idx:end_idx,:,gbx%Nlevels:1:-1),       &
                               gbx%dem_c(start_idx:end_idx,gbx%Nlevels:1:-1),            &
                               gbx%dem_s(start_idx:end_idx,gbx%Nlevels:1:-1),            &
                               cospIN%emiss_11)
 
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 0.67 micron optical depth
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    call cosp_simulator_optics(npoints,gbx%Ncolumns,gbx%Nlevels,                         &
                               sgx%frac_out(start_idx:end_idx,:,gbx%Nlevels:1:-1),       &
                               gbx%dtau_c(start_idx:end_idx,gbx%Nlevels:1:-1),           &
                               gbx%dtau_s(start_idx:end_idx,gbx%Nlevels:1:-1),           &
                               cospIN%tau_067)
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! LIDAR Polarized optics
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    call lidar_optics(npoints,gbx%Ncolumns,gbx%Nlevels,4,lidar_ice_type,                 &
                          mr_hydro(:,:,cospIN%Nlevels:1:-1,I_LSCLIQ),                    &
                          mr_hydro(:,:,cospIN%Nlevels:1:-1,I_LSCICE),                    &
                          mr_hydro(:,:,cospIN%Nlevels:1:-1,I_CVCLIQ),                    &
                          mr_hydro(:,:,cospIN%Nlevels:1:-1,I_CVCICE),                    &
                          gbx%Reff(start_idx:end_idx,cospIN%Nlevels:1:-1,I_LSCLIQ),      &
                          gbx%Reff(start_idx:end_idx,cospIN%Nlevels:1:-1,I_LSCICE),      &
                          gbx%Reff(start_idx:end_idx,cospIN%Nlevels:1:-1,I_CVCLIQ),      &
                          gbx%Reff(start_idx:end_idx,cospIN%Nlevels:1:-1,I_CVCICE),      & 
                          cospgridIN%pfull,cospgridIN%phalf,cospgridIN%at,               &
                          cospIN%beta_mol,cospIN%betatot,cospIN%taupart,                 &
                          cospIN%tau_mol,cospIN%tautot,cospIN%tautot_S_liq,              &
                          cospIN%tautot_S_ice, betatot_ice = cospIN%betatot_ice,         &
                          betatot_liq=cospIN%betatot_liq,tautot_ice=cospIN%tautot_ice,   &
                          tautot_liq = cospIN%tautot_liq)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! CLOUDSAT RADAR OPTICS
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    ! Allocate memory
    allocate(hm_matrix(N_HYDRO,npoints,gbx%Nlevels),                                     &
             re_matrix(N_HYDRO,npoints,gbx%Nlevels),                                     &
             Np_matrix(N_HYDRO,npoints,gbx%Nlevels))           

    do ij=1,gbx%Ncolumns
       do i=1,N_HYDRO
          hm_matrix(i,1:npoints,gbx%Nlevels:1:-1) = mr_hydro(:,ij,:,i)*1000._wp 
          re_matrix(i,1:npoints,gbx%Nlevels:1:-1) = Reff(:,ij,:,i)*1.e6_wp  
          Np_matrix(i,1:npoints,gbx%Nlevels:1:-1) = Np(:,ij,:,i)       
       enddo
       call quickbeam_optics(sd, rcfg_cloudsat,npoints,gbx%Nlevels, R_UNDEF, hm_matrix,  &
                             re_matrix, Np_matrix,                                       &
                             gbx%p(start_idx:end_idx,gbx%Nlevels:1:-1),                  & 
                             gbx%T(start_idx:end_idx,gbx%Nlevels:1:-1),                  &
                             gbx%sh(start_idx:end_idx,gbx%Nlevels:1:-1),cmpGases,        &
                             cospIN%z_vol_cloudsat(1:npoints,ij,:),                      &
                             cospIN%kr_vol_cloudsat(1:npoints,ij,:),                     &
                             cospIN%g_vol_cloudsat(1:npoints,ij,:))
    enddo
    
    ! Deallocate memory
    deallocate(hm_matrix,re_matrix,Np_matrix)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! MODIS optics
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Allocate memory
    allocate(MODIS_cloudWater(npoints,gbx%Ncolumns,gbx%Nlevels),                         &
             MODIS_cloudIce(npoints,gbx%Ncolumns,gbx%Nlevels),                           &
             MODIS_waterSize(npoints,gbx%Ncolumns,gbx%Nlevels),                          &
             MODIS_iceSize(npoints,gbx%Ncolumns,gbx%Nlevels),                            &
             MODIS_opticalThicknessLiq(npoints,gbx%Ncolumns,gbx%Nlevels),                &
             MODIS_opticalThicknessIce(npoints,gbx%Ncolumns,gbx%Nlevels))
    ! Cloud water
    call cosp_simulator_optics(npoints,gbx%Ncolumns,gbx%Nlevels,                         &
                               sgx%frac_out(start_idx:end_idx,:,:),                      &
                               mr_hydro(:,:,:,I_CVCLIQ),mr_hydro(:,:,:,I_LSCLIQ),        &
                               MODIS_cloudWater(:, :, gbx%Nlevels:1:-1))   
    ! Cloud ice
    call cosp_simulator_optics(npoints,gbx%Ncolumns,gbx%Nlevels,                         &
                               sgx%frac_out(start_idx:end_idx,:,:),                      &
                               mr_hydro(:,:,:,I_CVCICE), mr_hydro(:,:,:,I_LSCICE),       &
                               MODIS_cloudIce(:, :, gbx%Nlevels:1:-1))  
    ! Water droplet size
    call cosp_simulator_optics(npoints,gbx%Ncolumns,gbx%Nlevels,                         &
                               sgx%frac_out(start_idx:end_idx,:,:),reff(:,:,:,I_CVCLIQ), &
                               reff(:,:,:,I_LSCLIQ),                                     &
                               MODIS_waterSize(:, :, gbx%Nlevels:1:-1))
    ! Ice crystal size
    call cosp_simulator_optics(npoints,gbx%Ncolumns,gbx%Nlevels,                         &
                               sgx%frac_out(start_idx:end_idx,:,:),reff(:,:,:,I_CVCICE), &
                               reff(:,:,:,I_LSCICE),                                     &
                               MODIS_iceSize(:, :, gbx%Nlevels:1:-1))
    ! Partition optical thickness into liquid and ice parts
    call modis_optics_partition(npoints,gbx%Nlevels,gbx%Ncolumns,                        &
                                MODIS_cloudWater,MODIS_cloudIce,MODIS_waterSize,         &
                                MODIS_iceSize,cospIN%tau_067,MODIS_opticalThicknessLiq,  &
                                MODIS_opticalThicknessIce)
    ! Compute assymetry parameter and single scattering albedo 
    call modis_optics(npoints,gbx%Nlevels,gbx%Ncolumns,num_trial_res,                    &
                      MODIS_opticalThicknessLiq, MODIS_waterSize*1.0e6_wp,               &
                      MODIS_opticalThicknessIce, MODIS_iceSize*1.0e6_wp,                 &
                      cospIN%fracLiq, cospIN%asym, cospIN%ss_alb)
    
    ! Deallocate memory
    deallocate(MODIS_cloudWater,MODIS_cloudIce,MODIS_WaterSize,MODIS_iceSize,            &
               MODIS_opticalThicknessLiq,MODIS_opticalThicknessIce,mr_hydro,             &
               Reff,Np)
    
  end subroutine subsample_and_optics_v1p4
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !                          SUBROUTINE construct_cosp_config	
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE construct_cosp_config(cosp_nl,N_OUT_LIST,cfg)
    character(len=*),intent(in) :: cosp_nl
    integer,intent(in) :: N_OUT_LIST
    type(cosp_config),intent(out) :: cfg
    
    ! Local variables
    integer :: i
    logical :: Lradar_sim,Llidar_sim,Lisccp_sim,Lmodis_sim,Lmisr_sim,Lrttov_sim,        &
               Lparasol_sim,Lalbisccp,Latb532,Lboxptopisccp,Lboxtauisccp,LcfadDbze94,   &
               LcfadLidarsr532,Lclcalipso2,Lclcalipso,Lclhcalipso,Lclisccp,Lcllcalipso, &
               Lclmcalipso,Lcltcalipso,Lcltlidarradar,Lpctisccp,Ldbze94,Ltauisccp,      &
               Lcltisccp,Lclcalipsoliq,Lclcalipsoice,Lclcalipsoun,clcalipsotmp,         &
               Lclcalipsotmpliq,Lclcalipsotmpice,Lclcalipsotmpun,Lcltcalipsoliq,        &
               Lcltcalipsoice,Lcltcalipsoun,Lclhcalipsoliq,Lclhcalipsoice,Lclhcalipsoun,&
               Lclmcalipsoliq,Lclmcalipsoice,Lclmcalipsoun,Lcllcalipsoliq,              &
               Lcllcalipsoice,Lcllcalipsoun,Ltoffset,LparasolRefl,LclMISR,Lmeantbisccp, &
               Lmeantbclrisccp,Lfracout,LlidarBetaMol532,Ltbrttov,Lcltmodis,Lclwmodis,  &
               Lclimodis,Lclhmodis,Lclmmodis,Lcllmodis,Ltautmodis,Ltauwmodis,Ltauimodis,&
               Ltautlogmodis,Ltauwlogmodis,Ltauilogmodis,Lreffclwmodis,Lreffclimodis,   &
               Lpctmodis,Llwpmodis,Liwpmodis,Lclmodis,Lclcalipsotmp
    
    namelist/COSP_OUTPUT/Lradar_sim,Llidar_sim,Lisccp_sim,Lmodis_sim,Lmisr_sim,         &
                         Lrttov_sim,Lparasol_sim,Lalbisccp,Latb532,Lboxptopisccp,       &
                         Lboxtauisccp,LcfadDbze94,LcfadLidarsr532,Lclcalipso2,          &
                         Lclcalipso,Lclhcalipso,Lclisccp,Lcllcalipso,Lclmcalipso,       &
                         Lcltcalipso,Lcltlidarradar,Lpctisccp,Ldbze94,Ltauisccp,        &
                         Lclcalipsoliq,Lclcalipsoice,Lclcalipsoun,Lclcalipsotmp,        &
                         Lclcalipsotmpliq,Lclcalipsotmpice,Lclcalipsotmpun,             &
                         Lcltcalipsoliq,Lcltcalipsoice,Lcltcalipsoun,Lclhcalipsoliq,    &
                         Lclhcalipsoice,Lclhcalipsoun,Lclmcalipsoliq,Lclmcalipsoice,    &
                         Lclmcalipsoun,Lcllcalipsoliq,Lcllcalipsoice,Lcllcalipsoun,     &
                         Lcltisccp,Ltoffset,LparasolRefl,LclMISR,Lmeantbisccp,          &
                         Lmeantbclrisccp,Lfracout,LlidarBetaMol532,Ltbrttov,Lcltmodis,  &
                         Lclwmodis,Lclimodis,Lclhmodis,Lclmmodis,Lcllmodis,Ltautmodis,  &
                         Ltauwmodis,Ltauimodis,Ltautlogmodis,Ltauwlogmodis,             &
                         Ltauilogmodis,Lreffclwmodis,Lreffclimodis,Lpctmodis,Llwpmodis, &
                         Liwpmodis,Lclmodis
    
    allocate(cfg%out_list(N_OUT_LIST))
    do i=1,N_OUT_LIST
       cfg%out_list(i)=''
    enddo
    open(10,file=cosp_nl,status='old')
    read(10,nml=cosp_output)
    close(10)
    
    ! Deal with dependencies
    if (.not.Lradar_sim) then
       LcfadDbze94    = .false.
       Lclcalipso2    = .false.
       Lcltlidarradar = .false. ! Needs radar & lidar
       Ldbze94        = .false.
       Lclcalipso2    = .false. ! Needs radar & lidar
    endif
    if (.not.Llidar_sim) then
       Latb532 = .false.
       LcfadLidarsr532  = .false.
       Lclcalipso2      = .false. ! Needs radar & lidar
       Lclcalipso       = .false.
       Lclhcalipso      = .false.
       Lcllcalipso      = .false.
       Lclmcalipso      = .false.
       Lcltcalipso      = .false.
       Lcltlidarradar   = .false.
       LparasolRefl     = .false.
       LlidarBetaMol532 = .false.
       Lcltlidarradar   = .false. ! Needs radar & lidar
       Lclcalipsoliq    = .false.
       Lclcalipsoice    = .false.
       Lclcalipsoun     = .false.
       Lclcalipsotmp    = .false.
       Lclcalipsotmpun  = .false.
       Lclcalipsotmpliq = .false.
       Lclcalipsotmpice = .false.
       Lclhcalipsoliq   = .false.
       Lcllcalipsoliq   = .false.
       Lclmcalipsoliq   = .false.
       Lcltcalipsoliq   = .false.
       Lclhcalipsoice   = .false.
       Lcllcalipsoice   = .false.
       Lclmcalipsoice   = .false.
       Lcltcalipsoice   = .false.
       Lclhcalipsoun    = .false.
       Lcllcalipsoun    = .false.
       Lclmcalipsoun    = .false.
       Lcltcalipsoun    = .false.
    endif
    if (.not.Lisccp_sim) then
       Lalbisccp        = .false.
       Lboxptopisccp    = .false.
       Lboxtauisccp     = .false.
       Lclisccp         = .false.
       Lpctisccp        = .false.
       Ltauisccp        = .false.
       Lcltisccp        = .false.
       Lmeantbisccp     = .false.
       Lmeantbclrisccp  = .false.
    endif
    if (.not.Lmisr_sim) then
       LclMISR          = .false.
    endif
    if (.not.Lrttov_sim) then
       Ltbrttov         = .false.
    endif
    if ((.not.Lradar_sim) .and. (.not.Llidar_sim) .and. (.not.Lisccp_sim) .and.         &
         (.not.Lmisr_sim)) then
       Lfracout         = .false.
    endif
    if (.not.Lmodis_sim) then
       Lcltmodis        = .false.
       Lclwmodis        = .false.
       Lclimodis        = .false.
       Lclhmodis        = .false.
       Lclmmodis        = .false.
       Lcllmodis        = .false.
       Ltautmodis       = .false.
       Ltauwmodis       = .false.
       Ltauimodis       = .false.
       Ltautlogmodis    = .false.
       Ltauwlogmodis    = .false.
       Ltauilogmodis    = .false.
       Lreffclwmodis    = .false.
       Lreffclimodis    = .false.
       Lpctmodis        = .false.
       Llwpmodis        = .false.
       Liwpmodis        = .false.
       Lclmodis         = .false.
    endif
    if (Lmodis_sim) then
       Lisccp_sim       = .true.
    endif
    
    cfg%Lstats = .false.
    if ((Lradar_sim).or.(Llidar_sim).or.(Lisccp_sim)) then
       cfg%Lstats = .true.
    endif
    
    ! Copy instrument flags to cfg structure
    cfg%Lradar_sim   = Lradar_sim
    cfg%Llidar_sim   = Llidar_sim
    cfg%Lisccp_sim   = Lisccp_sim
    cfg%Lmodis_sim   = Lmodis_sim
    cfg%Lmisr_sim    = Lmisr_sim
    cfg%Lrttov_sim   = Lrttov_sim
    cfg%Lparasol_sim = Lparasol_sim
    
    ! Flag to control output to file
    cfg%Lwrite_output = .false.
    if (cfg%Lstats.or.cfg%Lmisr_sim.or.cfg%Lrttov_sim) then
       cfg%Lwrite_output = .true.
    endif
    
    ! Output diagnostics
    i = 1
    if (Lalbisccp)         cfg%out_list(1)  = 'albisccp'
    if (Latb532)           cfg%out_list(2)  = 'atb532'
    if (Lboxptopisccp)     cfg%out_list(3)  = 'boxptopisccp'
    if (Lboxtauisccp)      cfg%out_list(4)  = 'boxtauisccp'
    if (LcfadDbze94)       cfg%out_list(5)  = 'cfadDbze94'
    if (LcfadLidarsr532)   cfg%out_list(6)  = 'cfadLidarsr532'
    if (Lclcalipso2)       cfg%out_list(7)  = 'clcalipso2'
    if (Lclcalipso)        cfg%out_list(8)  = 'clcalipso'
    if (Lclhcalipso)       cfg%out_list(9)  = 'clhcalipso'
    if (Lclisccp)          cfg%out_list(10) = 'clisccp'
    if (Lcllcalipso)       cfg%out_list(11) = 'cllcalipso'
    if (Lclmcalipso)       cfg%out_list(12) = 'clmcalipso'
    if (Lcltcalipso)       cfg%out_list(13) = 'cltcalipso'
    if (Lcllcalipsoice)    cfg%out_list(14) = 'cllcalipsoice'
    if (Lclmcalipsoice)    cfg%out_list(15) = 'clmcalipsoice'
    if (Lclhcalipsoice)    cfg%out_list(16) = 'clhcalipsoice'
    if (Lcltcalipsoice)    cfg%out_list(17) = 'cltcalipsoice'
    if (Lcllcalipsoliq)    cfg%out_list(18) = 'cllcalipsoliq'
    if (Lclmcalipsoliq)    cfg%out_list(19) = 'clmcalipsoliq'
    if (Lclhcalipsoliq)    cfg%out_list(20) = 'clhcalipsoliq'
    if (Lcltcalipsoliq)    cfg%out_list(21) = 'cltcalipsoliq'
    if (Lcllcalipsoun)     cfg%out_list(22) = 'cllcalipsoun'
    if (Lclmcalipsoun)     cfg%out_list(23) = 'clmcalipsoun'
    if (Lclhcalipsoun)     cfg%out_list(24) = 'clhcalipsoun'
    if (Lcltcalipsoun)     cfg%out_list(25) = 'cltcalipsoun'
    if (Lclcalipsoice)     cfg%out_list(26) = 'clcalipsoice'
    if (Lclcalipsoliq)     cfg%out_list(27) = 'clcalipsoliq'
    if (Lclcalipsoun)      cfg%out_list(28) = 'clcalipsoun'
    if (Lclcalipsotmp)     cfg%out_list(29) = 'clcalipsotmp'
    if (Lclcalipsotmpice)  cfg%out_list(30) = 'clcalipsotmpice'
    if (Lclcalipsotmpliq)  cfg%out_list(31) = 'clcalipsotmpliq'
    if (Lclcalipsotmpun)   cfg%out_list(32) = 'clcalipsotmpun'
    if (Lcltlidarradar)    cfg%out_list(33) = 'cltlidarradar'
    if (Lpctisccp)         cfg%out_list(34) = 'pctisccp'
    if (Ldbze94)           cfg%out_list(35) = 'dbze94'
    if (Ltauisccp)         cfg%out_list(36) = 'tauisccp'
    if (Lcltisccp)         cfg%out_list(37) = 'cltisccp'
    if (Ltoffset)          cfg%out_list(38) = 'toffset'
    if (LparasolRefl)      cfg%out_list(39) = 'parasolRefl'
    if (LclMISR)           cfg%out_list(40) = 'clMISR'
    if (Lmeantbisccp)      cfg%out_list(41) = 'meantbisccp'
    if (Lmeantbclrisccp)   cfg%out_list(42) = 'meantbclrisccp'
    if (Lfracout)          cfg%out_list(43) = 'fracout'
    if (LlidarBetaMol532)  cfg%out_list(44) = 'lidarBetaMol532'
    if (Ltbrttov)          cfg%out_list(45) = 'tbrttov'
    if (Lcltmodis)         cfg%out_list(46) = 'cltmodis'
    if (Lclwmodis)         cfg%out_list(47) = 'clwmodis'
    if (Lclimodis)         cfg%out_list(48) = 'climodis'
    if (Lclhmodis)         cfg%out_list(49) = 'clhmodis'
    if (Lclmmodis)         cfg%out_list(50) = 'clmmodis'
    if (Lcllmodis)         cfg%out_list(51) = 'cllmodis'
    if (Ltautmodis)        cfg%out_list(52) = 'tautmodis'
    if (Ltauwmodis)        cfg%out_list(53) = 'tauwmodis'
    if (Ltauimodis)        cfg%out_list(54) = 'tauimodis'
    if (Ltautlogmodis)     cfg%out_list(55) = 'tautlogmodis'
    if (Ltauwlogmodis)     cfg%out_list(56) = 'tauwlogmodis'
    if (Ltauilogmodis)     cfg%out_list(57) = 'tauilogmodis'
    if (Lreffclwmodis)     cfg%out_list(58) = 'reffclwmodis'
    if (Lreffclimodis)     cfg%out_list(59) = 'reffclimodis'
    if (Lpctmodis)         cfg%out_list(60) = 'pctmodis'
    if (Llwpmodis)         cfg%out_list(61) = 'lwpmodis'
    if (Liwpmodis)         cfg%out_list(62) = 'iwpmodis'
    if (Lclmodis)          cfg%out_list(63) = 'clmodis'
    
    ! Copy diagnostic flags to cfg structure
    ! ISCCP simulator  
    cfg%Lalbisccp        = Lalbisccp
    cfg%Latb532          = Latb532
    cfg%Lboxptopisccp    = Lboxptopisccp
    cfg%Lboxtauisccp     = Lboxtauisccp
    cfg%Lmeantbisccp     = Lmeantbisccp
    cfg%Lmeantbclrisccp  = Lmeantbclrisccp
    cfg%Lclisccp         = Lclisccp
    cfg%Lpctisccp        = Lpctisccp
    cfg%Ltauisccp        = Ltauisccp
    cfg%Lcltisccp        = Lcltisccp
    
    ! CloudSat simulator  
    cfg%Ldbze94          = Ldbze94
    cfg%LcfadDbze94      = LcfadDbze94
    
    ! CALIPSO/PARASOL simulator  
    cfg%LcfadLidarsr532  = LcfadLidarsr532
    cfg%Lclcalipso2      = Lclcalipso2
    cfg%Lclcalipso       = Lclcalipso
    cfg%Lclhcalipso      = Lclhcalipso
    cfg%Lcllcalipso      = Lcllcalipso
    cfg%Lclmcalipso      = Lclmcalipso
    cfg%Lcltcalipso      = Lcltcalipso
    cfg%Lclhcalipsoice   = Lclhcalipsoice
    cfg%Lcllcalipsoice   = Lcllcalipsoice
    cfg%Lclmcalipsoice   = Lclmcalipsoice
    cfg%Lcltcalipsoice   = Lcltcalipsoice
    cfg%Lclhcalipsoliq   = Lclhcalipsoliq
    cfg%Lcllcalipsoliq   = Lcllcalipsoliq
    cfg%Lclmcalipsoliq   = Lclmcalipsoliq
    cfg%Lcltcalipsoliq   = Lcltcalipsoliq
    cfg%Lclhcalipsoun    = Lclhcalipsoun
    cfg%Lcllcalipsoun    = Lcllcalipsoun
    cfg%Lclmcalipsoun    = Lclmcalipsoun
    cfg%Lcltcalipsoun    = Lcltcalipsoun
    cfg%Lclcalipsoice    = Lclcalipsoice
    cfg%Lclcalipsoliq    = Lclcalipsoliq
    cfg%Lclcalipsoun     = Lclcalipsoun
    cfg%Lclcalipsotmp    = Lclcalipsotmp
    cfg%Lclcalipsotmpice = Lclcalipsotmpice
    cfg%Lclcalipsotmpliq = Lclcalipsotmpliq
    cfg%Lclcalipsotmpun  = Lclcalipsotmpun
    cfg%Lcltlidarradar   = Lcltlidarradar
    cfg%LparasolRefl     = LparasolRefl
    
    ! MISR simulator  
    cfg%LclMISR          = LclMISR
    
    ! Other
    cfg%Ltoffset         = Ltoffset
    cfg%Lfracout         = Lfracout
    cfg%LlidarBetaMol532 = LlidarBetaMol532
    
    ! RTTOV
    cfg%Ltbrttov         = Ltbrttov
    
    ! MODIS simulator  
    cfg%Lcltmodis        = Lcltmodis
    cfg%Lclwmodis        = Lclwmodis
    cfg%Lclimodis        = Lclimodis
    cfg%Lclhmodis        = Lclhmodis
    cfg%Lclmmodis        = Lclmmodis
    cfg%Lcllmodis        = Lcllmodis
    cfg%Ltautmodis       = Ltautmodis
    cfg%Ltauwmodis       = Ltauwmodis
    cfg%Ltauimodis       = Ltauimodis
    cfg%Ltautlogmodis    = Ltautlogmodis
    cfg%Ltauwlogmodis    = Ltauwlogmodis
    cfg%Ltauilogmodis    = Ltauilogmodis
    cfg%Lreffclwmodis    = Lreffclwmodis
    cfg%Lreffclimodis    = Lreffclimodis
    cfg%Lpctmodis        = Lpctmodis
    cfg%Llwpmodis        = Llwpmodis
    cfg%Liwpmodis        = Liwpmodis
    cfg%Lclmodis         = Lclmodis
    
  END SUBROUTINE construct_cosp_config

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !                        SUBROUTINE construct_cosp_subgrid
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_COSP_SUBGRID(Npoints,Ncolumns,Nlevels,y)
    ! Inputs
    integer,intent(in) :: &
         Npoints,  & ! Number of gridpoints
         Ncolumns, & ! Number of columns
         Nlevels     ! Number of levels
    ! Outputs
    type(cosp_subgrid),intent(out) :: y
    
    ! Dimensions
    y%Npoints  = Npoints
    y%Ncolumns = Ncolumns
    y%Nlevels  = Nlevels
    
    ! Allocate
    allocate(y%frac_out(Npoints,Ncolumns,Nlevels))
    if (Ncolumns > 1) then
       allocate(y%prec_frac(Npoints,Ncolumns,Nlevels))
    else ! CRM mode, not needed
       allocate(y%prec_frac(1,1,1))
    endif
    
    ! Initialize
    y%prec_frac = 0._wp
    y%frac_out  = 0._wp
  END SUBROUTINE CONSTRUCT_COSP_SUBGRID
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !                         SUBROUTINE construct_cosp_gridbox
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_COSP_GRIDBOX(time,time_bnds,                                      &
                                    Npoints,Nlevels,Ncolumns,Nhydro,Nprmts_max_hydro,    &
                                    Naero,Nprmts_max_aero,Npoints_it,isccp_emsfc_lw,     &
                                    Nchan,ZenAng,SurfEm,co2,ch4,n2o,co,lon,lat,p,        &
                                    ph,zlev,zlev_half,T,sh,cca,tca,skt,landmask,         &
                                    mr_ozone,u_wind,v_wind,sunlit,fl_lsrain,fl_lssnow,   &
                                    fl_lsgrpl,fl_ccrain,fl_ccsnow,dtau_s,dtau_c,dem_s,   &
                                    dem_c,Reff,mr_lsliq,mr_lsice,mr_ccliq,mr_ccice,y)
    ! Inputs
    double precision,intent(in),target :: &
         time,          & ! Time since start of run [days] 
         time_bnds(2)     ! Time boundaries
    integer,intent(in),target :: &
         Npoints,           & ! Number of gridpoints
         Nlevels,           & ! Number of levels
         Ncolumns,          & ! Number of columns
         Nhydro,            & ! Number of hydrometeors
         Nprmts_max_hydro,  & ! Max number of parameters for hydrometeor size 
                              ! distributions
         Naero,             & ! Number of aerosol species
         Nprmts_max_aero,   & ! Max number of parameters for aerosol size distributions
         Npoints_it,        & ! Number of gridpoints processed in one iteration
         Nchan                ! RTTOV number of channels
    
    real(wp),intent(in),target :: &
         isccp_emsfc_lw,   & ! 11microm surface emissivity
         co2,              & ! CO2 
         ch4,              & ! CH4
         n2o,              & ! N2O
         co,               & ! CO
         ZenAng              ! RTTOV zenith abgle
    real(wp),intent(in),dimension(Nchan),target :: &
         SurfEm
    real(wp),intent(in),dimension(Npoints),target :: &
         lon,        &
         lat,        &
         u_wind,     &
         v_wind,     &
         skt,        &
         landmask,   &
         sunlit
    
    real(wp),intent(in),dimension(Npoints,Nlevels),target :: &
         p,          &
         ph,         &
         zlev,       &
         T,          & 
         sh,         &
         cca,        &
         tca,        &
         mr_ozone,   &
         fl_lsrain,  &
         fl_lssnow,  &
         fl_lsgrpl,  &
         fl_ccrain,  &
         fl_ccsnow,  &
         dtau_c,     &
         dtau_s,     &
         dem_c,      &
         dem_s!,      &
 !        mr_lsice,   &
!         mr_lsliq,   &
!         mr_ccliq,   &
!         mr_ccice
    real(wp),intent(in),dimension(Npoints,Nlevels+1),target ::&
         zlev_half
    real(wp),intent(in),dimension(Npoints,Nlevels,N_HYDRO) :: &
         Reff
    real(wp),intent(in),dimension(Npoints,Nlevels) :: &
         mr_lsice,   &
         mr_lsliq,   &
         mr_ccliq,   &
         mr_ccice    
    
    ! Outputs
    type(cosp_gridbox),intent(out) :: y
    
    ! Local variables
    integer :: k

    ! Dimensions and scalars
    y%Npoints           => Npoints
    y%Nlevels           => Nlevels
    y%Ncolumns          => Ncolumns
    y%Nhydro            => Nhydro
    y%Nprmts_max_hydro  => Nprmts_max_hydro
    y%Naero             => Naero
    y%Nprmts_max_aero   => Nprmts_max_aero
    y%Npoints_it        => Npoints_it
    y%isccp_emsfc_lw    => isccp_emsfc_lw
    y%time              => time
    y%time_bnds         => time_bnds
   
    ! Allocate surface information and geolocation
     allocate(y%toffset(Npoints))

    ! Hydrometeors concentration and distribution parameters
    allocate(y%mr_hydro(Npoints,Nlevels,Nhydro),y%Reff(Npoints,Nlevels,Nhydro),          &
             y%Np(Npoints,Nlevels,Nhydro),y%dist_prmts_hydro(Nprmts_max_hydro,Nhydro))

    ! Aerosols concentration and distribution parameters
    allocate(y%conc_aero(Npoints,Nlevels,Naero), y%dist_type_aero(Naero),                &
             y%dist_prmts_aero(Npoints,Nlevels,Nprmts_max_aero,Naero))
    
    ! Point towards inputs where applicable and define where necessary
    y%longitude              => lon
    y%latitude               => lat
    y%p                      => p
    y%ph                     => ph
    y%zlev                   => zlev
    y%zlev_half              => zlev_half
    y%T                      => T
    y%sh                     => sh
    y%cca                    => cca
    y%tca                    => tca
    y%psfc                   => ph(:,1)
    y%skt                    => skt
    y%land                   => landmask
    y%mr_ozone               => mr_ozone
    y%u_wind                 => u_wind
    y%v_wind                 => v_wind
    y%sunlit                 => sunlit
    y%rain_ls                => fl_lsrain
    y%snow_ls                => fl_lssnow
    y%grpl_ls                => fl_lsgrpl
    y%rain_cv                => fl_ccrain
    y%snow_cv                => fl_ccsnow
    y%dtau_s                 => dtau_s
    y%dtau_c                 => dtau_c
    y%dem_s                  => dem_s
    y%dem_c                  => dem_c   
    y%dlev                   = 0._wp               
    y%dist_prmts_hydro       = 0._wp
    y%conc_aero              = 0._wp
    y%dist_type_aero         = 0   
    y%dist_prmts_aero        = 0._wp 
    y%Np                     = 0._wp     
    y%Reff                   = Reff
    y%Reff(:,:,I_LSRAIN)     = 0.0   
    y%mr_hydro               = 0._wp     
    y%mr_hydro(:,:,I_LSCLIQ) = mr_lsliq
    y%mr_hydro(:,:,I_LSCICE) = mr_lsice
    y%mr_hydro(:,:,I_CVCLIQ) = mr_ccliq
    y%mr_hydro(:,:,I_CVCICE) = mr_ccice
    
    ! RTTOV parameters
    y%surfem => surfem
    y%Nchan  => Nchan
    y%ZenAng => ZenAng
    y%co2    => co2
    y%ch4    => ch4
    y%n2o    => n2o
    y%co     => co  
      
    ! Toffset. This assumes that time is the mid-point of the interval.
    y%toffset = -0.5_wp*3._wp/24._wp + 3._wp/24._wp*([1:Npoints]-0.5)/Npoints
    
  END SUBROUTINE CONSTRUCT_COSP_GRIDBOX
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !                            SUBROUTINE destroy_cosp_subgrid
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine destroy_cosp_subgrid(y)
    type(cosp_subgrid),intent(inout) :: y   
    deallocate(y%prec_frac, y%frac_out)
  end subroutine destroy_cosp_subgrid
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !                            SUBROUTINE destroy_cosp_gridbox
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine destroy_cosp_gridbox(y,dglobal,save_LUT)
    type(cosp_gridbox),intent(inout)       :: y
    logical,           intent(in),optional :: dglobal
    logical,           intent(in),optional :: save_LUT
    logical                                :: local_save_LUT
    
    if (present(save_LUT)) then
       local_save_LUT = save_LUT
    else
       local_save_LUT = .false.!RADAR_SIM_UPDATE_scale_LUTs_flag
    endif
    
    ! Save any updates to radar simulator LUT
    if (local_save_LUT) call save_scale_LUTs(rcfg_cloudsat)
    
    deallocate(y%mr_hydro, y%dist_prmts_hydro, y%conc_aero, y%dist_type_aero,            &
               y%dist_prmts_aero, y%Reff,y%Np,y%toffset)              
  end subroutine destroy_cosp_gridbox
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !                                    END MODULE
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module MOD_COSP_INTERFACE_v1p5
