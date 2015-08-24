#include "cosp_defs.h"
MODULE MOD_COSP_INTERFACE_v1p4
  use COSP_KINDS,                 ONLY: wp,dp
  use MOD_COSP,                   ONLY: cosp_outputs,construct_cosp_outputs,             &
                                        destroy_cosp_outputs
  USE MOD_COSP_CONFIG,            ONLY: PARASOL_NREFL,R_UNDEF,                           &
                                        SR_BINS,LIDAR_NCAT,LIDAR_NTEMP,DBZE_BINS,        &
                                        numMISRHgtBins,N_HYDRO,numMODISReffLiqBins,      &
                                        numMODISReffIceBins,numMODISTauBins,             &
                                        numMODISPresBins
  USE MOD_COSP_INTERFACE_v1p5,    ONLY: cosp_config,cosp_subgrid,cosp_gridbox,&
                                        cosp_interface_init,cosp_interface_v1p5,         &
                                        construct_cosp_gridbox,destroy_cosp_gridbox,     &
                                        I_LSCLIQ,I_LSCICE,I_CVCICE,I_CVCLIQ,I_LSRAIN
  USE quickbeam,                  ONLY: maxhclass,nRe_types,nd,mt_ntt,Re_BIN_LENGTH,     &
                                        Re_MAX_BIN

  implicit none
  
  character(len=120),parameter :: &
       RADAR_SIM_LUT_DIRECTORY = './'
  logical,parameter :: &
       RADAR_SIM_LOAD_scale_LUTs_flag   = .false., &
       RADAR_SIM_UPDATE_scale_LUTs_flag = .false.     
       
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !                                  TYPE COSP_VGRID
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  TYPE COSP_VGRID
     logical ::  &
          use_vgrid,  & ! Logical flag that indicates change of grid
          csat_vgrid    ! Flag for Cloudsat grid
     integer :: &
          Npoints,    & ! Number of sampled points
          Ncolumns,   & ! Number of subgrid columns
          Nlevels,    & ! Number of model levels
          Nlvgrid       ! Number of levels of new grid
     real(wp), dimension(:), pointer :: &
          z,          & ! Height of new level              (Nlvgrid)
          zl,         & ! Lower boundaries of new levels   (Nlvgrid)
          zu,         & ! Upper boundaries of new levels   (Nlvgrid)
          mz,         & ! Height of model levels           (Nlevels)
          mzl,        & ! Lower boundaries of model levels (Nlevels)
          mzu           ! Upper boundaries of model levels (Nlevels)
  END TYPE COSP_VGRID
  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !  TYPE class_param
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  type class_param
     ! Variables used to store hydrometeor "default" properties
     real(dp),dimension(maxhclass) :: p1,p2,p3,dmin,dmax,apm,bpm,rho
     integer, dimension(maxhclass) :: dtype,col,cp,phase
     
     ! Radar properties
     real(dp) :: freq,k2
     integer  :: nhclass           ! number of hydrometeor classes in use
     integer  :: use_gas_abs, do_ray
     
     ! Defines location of radar relative to hgt_matrix.   
     logical :: radar_at_layer_one ! If true radar is assume to be at the edge 
                                   ! of the first layer, if the first layer is the
                                   ! surface than a ground-based radar.   If the
                                   ! first layer is the top-of-atmosphere, then
                                   ! a space borne radar. 
     
     ! Variables used to store Z scale factors
     character(len=240)                             :: scale_LUT_file_name
     logical                                        :: load_scale_LUTs, update_scale_LUTs
     logical, dimension(maxhclass,nRe_types)        :: N_scale_flag
     logical, dimension(maxhclass,mt_ntt,nRe_types) :: Z_scale_flag,Z_scale_added_flag
     real(dp),dimension(maxhclass,mt_ntt,nRe_types) :: Ze_scaled,Zr_scaled,kr_scaled
     real(dp),dimension(maxhclass,nd,nRe_types)     :: fc, rho_eff     
  end type class_param
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! TYPE cosp_gridbox_v1p4
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  TYPE COSP_GRIDBOX_v1p4
     integer :: &
          Npoints,          & ! Number of gridpoints
          Nlevels,          & ! Number of levels
          Ncolumns,         & ! Number of columns
          Nhydro,           & ! Number of hydrometeors
          Nprmts_max_hydro, & ! Max number of parameters for hydrometeor size distribution
          Naero,            & ! Number of aerosol species
          Nprmts_max_aero,  & ! Max number of parameters for aerosol size distributions
          Npoints_it          ! Max number of gridpoints to be processed in one iteration
     
     ! Time [days]
     double precision :: time
     double precision :: time_bnds(2)
     
     ! Radar ancillary info
     real(wp) :: &
          radar_freq,    & ! Radar frequency [GHz]
          k2               ! |K|^2, -1=use frequency dependent default
     integer :: surface_radar,  & ! surface=1, spaceborne=0
          use_mie_tables, & ! use a precomputed loopup table? yes=1,no=0
          use_gas_abs,    & ! include gaseous absorption? yes=1,no=0
          do_ray,         & ! calculate/output Rayleigh refl=1, not=0
          melt_lay          ! melting layer model off=0, on=1
     
     
     ! Structures used by radar simulator that need to be set only ONCE per 
     ! radar configuration (e.g. freq, pointing direction) ... added by roj Feb 2008
     type(class_param) :: &
          hp     ! Structure used by radar simulator to store Ze and N scaling constants 
                 ! and other information
     integer :: &
          nsizes ! Number of discrete drop sizes (um) used to represent the distribution
     
     ! Lidar
     integer :: &
          lidar_ice_type ! Ice particle shape hypothesis in lidar calculations
                         ! (ice_type=0 for spheres, ice_type=1 for non spherical particles)
    
     ! Radar
     logical :: &
          use_precipitation_fluxes, & ! True if precipitation fluxes are input to the 
                                      ! algorithm 
          use_reff                    ! True if Reff is to be used by radar (memory not 
                                      ! allocated)       
     
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
          q,         & ! Relative humidity to water (%)                       
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
     real(wp),dimension(:,:,:),pointer :: &
          mr_hydro         ! Mixing ratio of each hydrometeor 
                           ! (Npoints,Nlevels,Nhydro) [kg/kg]
     real(wp),dimension(:,:),pointer :: &
          dist_prmts_hydro ! Distributional parameters for hydrometeors 
                           ! (Nprmts_max_hydro,Nhydro)
     real(wp),dimension(:,:,:),pointer :: &
          Reff             ! Effective radius [m]. 
                           ! (Npoints,Nlevels,Nhydro)
     real(wp),dimension(:,:,:),pointer :: &
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
     integer :: &
          ! ISCCP_TOP_HEIGHT
          ! 1 = adjust top height using both a computed infrared brightness temperature and
          !     the visible optical depth to adjust cloud top pressure. Note that this 
          !     calculation is most appropriate to compare to ISCCP data during sunlit 
          !     hours.
          ! 2 = do not adjust top height, that is cloud top pressure is the actual cloud 
          !     top pressure in the model.
          ! 3 = adjust top height using only the computed infrared brightness temperature. 
          !     Note that this calculation is most appropriate to compare to ISCCP IR only 
          !     algortihm (i.e. you can compare to nighttime ISCCP data with this option)
          isccp_top_height, &
          ! ISCCP_TOP_HEIGHT_DIRECTION
          ! Direction for finding atmosphere pressure level with interpolated temperature 
          ! equal to the radiance determined cloud-top temperature
          ! 1 = find the *lowest* altitude (highest pressure) level with interpolated 
          !     temperature equal to the radiance determined cloud-top temperature
          ! 2 = find the *highest* altitude (lowest pressure) level with interpolated 
          !     temperature equal to the radiance determined cloud-top temperature
          !     ONLY APPLICABLE IF top_height EQUALS 1 or 3
          ! 1 = default setting, and matches all versions of ISCCP simulator with versions 
          !     numbers 3.5.1 and lower; 2 = experimental setting  
          isccp_top_height_direction, &
          ! Overlap type (1=max, 2=rand, 3=max/rand)
          isccp_overlap 
     real(wp) :: &
          isccp_emsfc_lw      ! 10.5 micron emissivity of surface (fraction)
     
     ! RTTOV inputs/options
     integer :: &
          plat,   & ! Satellite platform
          sat,    & ! Satellite
          inst,   & ! Instrument
          Nchan     ! Number of channels to be computed
     integer, dimension(:), pointer :: &
          Ichan     ! Channel numbers
     real(wp),dimension(:), pointer :: &
          Surfem    ! Surface emissivity
     real(wp) :: &
          ZenAng, & ! Satellite Zenith Angles
          co2,    & ! CO2 mixing ratio
          ch4,    & ! CH4 mixing ratio
          n2o,    & ! N2O mixing ratio
          co        ! CO mixing ratio
  END TYPE COSP_GRIDBOX_v1p4
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! TYPE cosp_modis
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  type cosp_modis
     integer,pointer ::                    & !
          Npoints                            ! Number of gridpoints
     real(wp),pointer,dimension(:) ::      & !  
          Cloud_Fraction_Total_Mean,       & ! L3 MODIS retrieved cloud fraction (total) 
          Cloud_Fraction_Water_Mean,       & ! L3 MODIS retrieved cloud fraction (liq) 
          Cloud_Fraction_Ice_Mean,         & ! L3 MODIS retrieved cloud fraction (ice) 
          Cloud_Fraction_High_Mean,        & ! L3 MODIS retrieved cloud fraction (high) 
          Cloud_Fraction_Mid_Mean,         & ! L3 MODIS retrieved cloud fraction (middle) 
          Cloud_Fraction_Low_Mean,         & ! L3 MODIS retrieved cloud fraction (low ) 
          Optical_Thickness_Total_Mean,    & ! L3 MODIS retrieved optical thickness (tot)
          Optical_Thickness_Water_Mean,    & ! L3 MODIS retrieved optical thickness (liq)
          Optical_Thickness_Ice_Mean,      & ! L3 MODIS retrieved optical thickness (ice)
          Optical_Thickness_Total_LogMean, & ! L3 MODIS retrieved log10 optical thickness 
          Optical_Thickness_Water_LogMean, & ! L3 MODIS retrieved log10 optical thickness 
          Optical_Thickness_Ice_LogMean,   & ! L3 MODIS retrieved log10 optical thickness
          Cloud_Particle_Size_Water_Mean,  & ! L3 MODIS retrieved particle size (liquid)
          Cloud_Particle_Size_Ice_Mean,    & ! L3 MODIS retrieved particle size (ice)
          Cloud_Top_Pressure_Total_Mean,   & ! L3 MODIS retrieved cloud top pressure
          Liquid_Water_Path_Mean,          & ! L3 MODIS retrieved liquid water path
          Ice_Water_Path_Mean                ! L3 MODIS retrieved ice water path
     real(wp),pointer,dimension(:,:,:) ::  &
          Optical_Thickness_vs_Cloud_Top_Pressure,  & ! Tau/Pressure joint histogram
          Optical_Thickness_vs_ReffICE,             & ! Tau/ReffICE joint histogram
          Optical_Thickness_vs_ReffLIQ                ! Tau/ReffLIQ joint histogram

  end type cosp_modis  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! TYPE cosp_misr	
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  TYPE COSP_MISR
     integer,pointer :: &
        Npoints,       & ! Number of gridpoints
        Ntau,          & ! Number of tau intervals
        Nlevels          ! Number of cth levels  
     real(wp),dimension(:,:,:),pointer ::   & !
        fq_MISR          ! Fraction of the model grid box covered by each of the MISR 
          				 ! cloud types
     real(wp),dimension(:,:),pointer ::   & !
        MISR_dist_model_layertops !  
     real(wp),dimension(:),pointer ::   & !
        MISR_meanztop, & ! Mean MISR cloud top height
        MISR_cldarea     ! Mean MISR cloud cover area
  END TYPE COSP_MISR  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !									Type cosp_rttov
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  TYPE COSP_RTTOV
     ! Dimensions
     integer,pointer :: &
        Npoints,  & ! Number of gridpoints
        Nchan       ! Number of channels
     
     ! Brightness temperatures (Npoints,Nchan)
     real(wp),pointer :: tbs(:,:)
  END TYPE COSP_RTTOV
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 !									TYPE cosp_isccp
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  TYPE COSP_ISCCP
     integer,pointer  ::&
        Npoints,      & ! Number of gridpoints.
        Ncolumns,     & ! Number of columns.
        Nlevels         ! Number of levels.
     real(wp),dimension(:),pointer :: &
        totalcldarea, & ! The fraction of model grid box columns with cloud somewhere in 
          				  ! them.
        meantb,       & ! Mean all-sky 10.5 micron brightness temperature.
        meantbclr,    & ! Mean clear-sky 10.5 micron brightness temperature.
        meanptop,     & ! Mean cloud top pressure (mb).
        meantaucld,   & ! Mean optical thickness.
        meanalbedocld   ! Mean cloud albedo.
     real(wp),dimension(:,:),pointer ::&
        boxtau,       & ! Optical thickness in each column   .
        boxptop         ! Cloud top pressure (mb) in each column.
     real(wp),dimension(:,:,:),pointer :: &
        fq_isccp        ! The fraction of the model grid box covered by each of the 49 
          			    ! ISCCP D level cloud types.
  END TYPE COSP_ISCCP
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! TYPE cosp_sglidar
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  type cosp_sglidar
     integer,pointer :: &
          Npoints,         & ! Number of sampled points
          Ncolumns,        & ! Number of subgrid columns
          Nlevels,         & ! Number of model levels
          Nhydro,          & ! Number of hydrometeors
          Nrefl              ! Number of parasol reflectances
     real(wp),dimension(:,:),pointer :: &
          beta_mol,      & ! Molecular backscatter
          temp_tot
     real(wp),dimension(:,:,:),pointer :: &
          betaperp_tot,  & ! Total backscattered signal
          beta_tot,      & ! Total backscattered signal
          tau_tot,       & ! Optical thickness integrated from top to level z
          refl             ! PARASOL reflectances 
  end type cosp_sglidar
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! TYPE cosp_lidarstats
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  type cosp_lidarstats
     integer,pointer :: &
          Npoints,         & ! Number of sampled points
          Ncolumns,        & ! Number of subgrid columns
          Nlevels,         & ! Number of model levels
          Nhydro,          & ! Number of hydrometeors
          Nrefl              ! Number of parasol reflectances
     real(wp), dimension(:,:,:),pointer :: &
          lidarcldphase,   & ! 3D "lidar" phase cloud fraction 
          cldlayerphase,   & ! low, mid, high-level lidar phase cloud cover
          lidarcldtmp,     & ! 3D "lidar" phase cloud temperature
          cfad_sr            ! CFAD of scattering ratio
     real(wp), dimension(:,:),pointer :: &
          lidarcld,        & ! 3D "lidar" cloud fraction 
          cldlayer,        & ! low, mid, high-level, total lidar cloud cover
          parasolrefl
     real(wp), dimension(:),pointer :: &
          srbval             ! SR bins in cfad_sr
  end type cosp_lidarstats  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! TYPE cosp_sgradar
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  type cosp_sgradar
     ! Dimensions
     integer,pointer :: &
          Npoints,            & ! Number of gridpoints
          Ncolumns,           & ! Number of columns
          Nlevels,            & ! Number of levels
          Nhydro                ! Number of hydrometeors
     real(wp),dimension(:,:),pointer :: &
          att_gas               ! 2-way attenuation by gases [dBZ] (Npoints,Nlevels)
     real(wp),dimension(:,:,:),pointer :: &
          Ze_tot                ! Effective reflectivity factor (Npoints,Ncolumns,Nlevels)
  end type cosp_sgradar
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! TYPE cosp_radarstats
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  type cosp_radarstats
     integer,pointer  :: &
          Npoints,            & ! Number of sampled points
          Ncolumns,           & ! Number of subgrid columns
          Nlevels,            & ! Number of model levels
          Nhydro                ! Number of hydrometeors
     real(wp), dimension(:,:,:), pointer :: &
          cfad_ze               ! Ze CFAD(Npoints,dBZe_bins,Nlevels)
     real(wp),dimension(:),pointer :: &
          radar_lidar_tcc       ! Radar&lidar total cloud amount, grid-box scale (Npoints)
     real(wp), dimension(:,:),pointer :: &
          lidar_only_freq_cloud !(Npoints,Nlevels)
  end type cosp_radarstats
      
contains
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !                            SUBROUTINE COSP_INTERFACE (v1.4)
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine cosp_interface_v1p4(overlap,Ncolumns,cfg,vgrid,gbxIN,sgx,sgradar,sglidar,   &
                                 isccp,misr,modis,rttov,stradar,stlidar)
    ! Inputs 
    integer,                intent(in)    :: overlap  ! Overlap type in SCOPS: 1=max, 
                                                      ! 2=rand, 3=max/rand
    integer,                intent(in)    :: Ncolumns ! Number of columns
    type(cosp_config),      intent(in)    :: cfg      ! Configuration options
    type(cosp_vgrid),target,intent(in)    :: vgrid    ! Information on vertical grid of 
                                                      ! stats
    type(cosp_subgrid),     intent(inout) :: sgx      ! Subgrid info
    type(cosp_sgradar),     intent(inout) :: sgradar  ! Output from radar simulator (pixel)
    type(cosp_sglidar),     intent(inout) :: sglidar  ! Output from lidar simulator (pixel)
    type(cosp_isccp),       intent(inout) :: isccp    ! Output from ISCCP simulator
    type(cosp_misr),        intent(inout) :: misr     ! Output from MISR simulator
    type(cosp_modis),       intent(inout) :: modis    ! Output from MODIS simulator
    type(cosp_rttov),       intent(inout) :: rttov    ! Output from RTTOV
    type(cosp_radarstats),  intent(inout) :: stradar  ! Summary statistics from cloudsat
                                                      ! simulator (gridbox)
    type(cosp_lidarstats),  intent(inout) :: stlidar  ! Output from LIDAR simulator (gridbox)
    type(cosp_gridbox_v1p4),intent(inout),target :: gbxIN ! COSP gridbox type from v1.4
                                                          ! Shares memory with new type
                                                          ! for v1.5.
    
    ! Inputs to cosp_interface_v1p5
    type(cosp_gridbox) :: gbxOUT   ! COSP gridbox type for v1.5
    
    ! Outputs from cosp_interface_v1p5
    type(cosp_outputs),target :: cospOUT  ! NEW derived type output that contains all 
    					                  ! simulator information
     ! Local variables
     character(len=32) :: &
        cospvID = 'COSP v1.4' ! COSP version ID				                  
#ifdef MMF_V3_SINGLE_MOMENT    					  
     character(len=64) :: &
        cloudsat_micro_scheme = 'MMF_v3_single_moment'
#endif
#ifdef MMF_V3p5_TWO_MOMENT
     character(len=64) :: &
        cloudsat_micro_scheme = 'MMF_v3.5_two_moment'
#endif 
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Call cosp_init
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    call cosp_interface_init(gbxIN%Npoints,gbxIN%Nlevels,gbxIN%Npoints_it,overlap,       &
                             gbxIN%use_precipitation_fluxes,gbxIN%radar_freq,            &
                             cloudsat_micro_scheme,gbxIN%k2,gbxIN%use_gas_abs,           &
                             gbxIN%do_ray,gbxIN%isccp_top_height,                        &
                             gbxIN%isccp_top_height_direction,                           &
                             gbxIN%zlev(:,gbxIN%Nlevels:1:-1),                           &
                             gbxIN%zlev_half(:,gbxIN%Nlevels:1:-1),gbxIN%surface_radar,  &
                             gbxIN%Nchan,gbxIN%Ichan,gbxIN%Plat,gbxIN%Sat,gbxIN%Inst,    &
                             gbxIN%lidar_ice_type,vgrid%use_vgrid,vgrid%Nlvgrid,         &
                             vgrid%csat_vgrid,cospvID)
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Call construct_cosp_outputs
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    call construct_cosp_outputs(cfg%Lpctisccp,cfg%Lclisccp,cfg%Lboxptopisccp,            &
                                cfg%Lboxtauisccp,cfg%Ltauisccp,cfg%Lcltisccp,            &
                                cfg%Lmeantbisccp,cfg%Lmeantbclrisccp,cfg%Lalbisccp,      &
                                cfg%LclMISR,cfg%Lcltmodis,cfg%Lclwmodis,cfg%Lclimodis,   &
                                cfg%Lclhmodis,cfg%Lclmmodis,cfg%Lcllmodis,cfg%Ltautmodis,&
                                cfg%Ltauwmodis,cfg%Ltauimodis,cfg%Ltautlogmodis,         &
                                cfg%Ltauwlogmodis,cfg%Ltauilogmodis,cfg%Lreffclwmodis,   &
                                cfg%Lreffclimodis,cfg%Lpctmodis,cfg%Llwpmodis,           &
                                cfg%Liwpmodis,cfg%Lclmodis,cfg%Latb532,                  &
                                cfg%LlidarBetaMol532,cfg%LcfadLidarsr532,cfg%Lclcalipso2,&
                                cfg%Lclcalipso,cfg%Lclhcalipso,cfg%Lcllcalipso,          &
                                cfg%Lclmcalipso,cfg%Lcltcalipso,cfg%Lcltlidarradar,      &
                                cfg%Lclcalipsoliq,cfg%Lclcalipsoice,cfg%Lclcalipsoun,    &
                                cfg%Lclcalipsotmp,cfg%Lclcalipsotmpliq,                  &
                                cfg%Lclcalipsotmpice,cfg%Lclcalipsotmpun,                &
                                cfg%Lcltcalipsoliq,cfg%Lcltcalipsoice,cfg%Lcltcalipsoun, &
                                cfg%Lclhcalipsoliq,cfg%Lclhcalipsoice,cfg%Lclhcalipsoun, &
                                cfg%Lclmcalipsoliq,cfg%Lclmcalipsoice,cfg%Lclmcalipsoun, &
                                cfg%Lcllcalipsoliq,cfg%Lcllcalipsoice,cfg%Lcllcalipsoun, &
                                cfg%LcfadDbze94,cfg%Ldbze94,cfg%Lparasolrefl,            &
                                cfg%Ltbrttov,gbxIN%Npoints,gbxIN%Ncolumns,gbxIN%Nlevels, &
                                vgrid%Nlvgrid,gbxIN%Nchan,cospOUT)    

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Copy input data from old derived type to new derived type
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gbxOUT%Npoints          => gbxIN%Npoints
    gbxOUT%Nlevels          => gbxIN%Nlevels
    gbxOUT%Ncolumns         => gbxIN%Ncolumns
    gbxOUT%Nhydro           => gbxIN%Nhydro
    gbxOUT%Nprmts_max_hydro => gbxIN%Nprmts_max_hydro
    gbxOUT%Naero            => gbxIN%Naero
    gbxOUT%Nprmts_max_aero  => gbxIN%Nprmts_max_aero
    gbxOUT%Npoints_it       => gbxIN%Npoints_it
    gbxOUT%time             => gbxIN%time
    gbxOUT%time_bnds        => gbxIN%time_bnds
    gbxOUT%nsizes           => gbxIN%nsizes
    gbxOUT%toffset          => gbxIN%toffset
    gbxOUT%longitude        => gbxIN%longitude
    gbxOUT%latitude         => gbxIN%latitude
    gbxOUT%zlev             => gbxIN%zlev
    gbxOUT%zlev_half        => gbxIN%zlev_half
    gbxOUT%dlev             => gbxIN%dlev
    gbxOUT%p                => gbxIN%p
    gbxOUT%ph               => gbxIN%ph
    gbxOUT%T                => gbxIN%T
    gbxOUT%sh               => gbxIN%sh
    gbxOUT%dtau_s           => gbxIN%dtau_s
    gbxOUT%dtau_c           => gbxIN%dtau_c
    gbxOUT%dem_s            => gbxIN%dem_s
    gbxOUT%dem_c            => gbxIN%dem_c
    gbxOUT%mr_ozone         => gbxIN%mr_ozone
    gbxOUT%land             => gbxIN%land
    gbxOUT%psfc             => gbxIN%psfc
    gbxOUT%sunlit           => gbxIN%sunlit
    gbxOUT%skt              => gbxIN%skt
    gbxOUT%u_wind           => gbxIN%u_wind
    gbxOUT%v_wind           => gbxIN%v_wind
    gbxOUT%tca              => gbxIN%tca
    gbxOUT%cca              => gbxIN%cca
    gbxOUT%rain_ls          => gbxIN%rain_ls
    gbxOUT%rain_cv          => gbxIN%rain_cv
    gbxOUT%snow_ls          => gbxIN%snow_ls
    gbxOUT%snow_cv          => gbxIN%snow_cv
    gbxOUT%grpl_ls          => gbxIN%grpl_ls
    gbxOUT%dist_prmts_hydro => gbxIN%dist_prmts_hydro
    gbxOUT%conc_aero        => gbxIN%conc_aero
    gbxOUT%dist_type_aero   => gbxIN%dist_type_aero
    gbxOUT%dist_prmts_aero  => gbxIN%dist_prmts_aero
    gbxOUT%isccp_emsfc_lw   => gbxIN%isccp_emsfc_lw
    gbxOUT%Nchan            => gbxIN%Nchan
    gbxOUT%Surfem           => gbxIN%Surfem
    gbxOUT%ZenAng           => gbxIN%ZenAng
    gbxOUT%co2              => gbxIN%co2
    gbxOUT%ch4              => gbxIN%ch4
    gbxOUT%n2o              => gbxIN%n2o
    gbxOUT%co               => gbxIN%co 
    allocate(gbxOUT%mr_hydro(gbxIN%Npoints,gbxIN%Nlevels,gbxIN%Nhydro),                  &
             gbxOUT%Reff(gbxIN%Npoints,gbxIN%Nlevels,gbxIN%Nhydro),                      &
             gbxOUT%Np(gbxIN%Npoints,gbxIN%Nlevels,gbxIN%Nhydro))
    gbxOUT%mr_hydro         =  gbxIN%mr_hydro
    gbxOUT%Reff             =  gbxIN%Reff
    gbxOUT%Np               =  gbxIN%Np

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Call cosp_interface_v1p5
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    call cosp_interface_V1p5(gbxIN%Npoints,gbxOUT,sgx,cospOUT)
 
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Copy output from v1.5 to v1.4
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    ! MISR
    if (cfg%Lmisr_sim) then
       if (cfg%LclMISR) misr%fq_MISR  => cospOUT%misr_fq
       ! *NOTE* These 3 fields are not output, but were part of the v1.4.0 cosp_misr, so
       !        they are still computed. Should probably have a logical to control these
       !        outputs.       
       misr%MISR_meanztop             => cospOUT%misr_meanztop
       misr%MISR_cldarea              => cospOUT%misr_cldarea
       misr%MISR_dist_model_layertops => cospOUT%misr_dist_model_layertops
    endif
    
    ! ISCCP
    if (cfg%Lisccp_sim) then
       if (cfg%Lboxtauisccp)    isccp%boxtau        => cospOUT%isccp_boxtau
       if (cfg%Lboxptopisccp)   isccp%boxptop       => cospOUT%isccp_boxptop
       if (cfg%Lclisccp)        isccp%fq_isccp      => cospOUT%isccp_fq
       if (cfg%Lcltisccp)       isccp%totalcldarea  => cospOUT%isccp_totalcldarea
       if (cfg%Lmeantbisccp)    isccp%meantb        => cospOUT%isccp_meantb
       if (cfg%Lmeantbclrisccp) isccp%meantbclr     => cospOUT%isccp_meantbclr
       if (cfg%Lpctisccp)       isccp%meanptop      => cospOUT%isccp_meanptop
       if (cfg%Ltauisccp)       isccp%meantaucld    => cospOUT%isccp_meantaucld
       if (cfg%Lalbisccp)       isccp%meanalbedocld => cospOUT%isccp_meanalbedocld
   endif

    ! MODIS
    if (cfg%Lmodis_sim) then
       if (cfg%Lcltmodis)     modis%Cloud_Fraction_Total_Mean =>                         &
                          cospOUT%modis_Cloud_Fraction_Total_Mean
       if (cfg%Lclwmodis)     modis%Cloud_Fraction_Water_Mean =>                         &
                          cospOUT%modis_Cloud_Fraction_Water_Mean
       if (cfg%Lclimodis)     modis%Cloud_Fraction_Ice_Mean =>                           &
                          cospOUT%modis_Cloud_Fraction_Ice_Mean
       if (cfg%Lclhmodis)     modis%Cloud_Fraction_High_Mean =>                          &
                          cospOUT%modis_Cloud_Fraction_High_Mean
       if (cfg%Lclmmodis)     modis%Cloud_Fraction_Mid_Mean =>                           &
                          cospOUT%modis_Cloud_Fraction_Mid_Mean
       if (cfg%Lcllmodis)     modis%Cloud_Fraction_Low_Mean =>                           &
                          cospOUT%modis_Cloud_Fraction_Low_Mean
       if (cfg%Ltautmodis)    modis%Optical_Thickness_Total_Mean =>                      &
                          cospOUT%modis_Optical_Thickness_Total_Mean
       if (cfg%Ltauwmodis)    modis%Optical_Thickness_Water_Mean =>                      &
                          cospOUT%modis_Optical_Thickness_Water_Mean
       if (cfg%Ltauimodis)    modis%Optical_Thickness_Ice_Mean =>                        &
                          cospOUT%modis_Optical_Thickness_Ice_Mean
       if (cfg%Ltautlogmodis) modis%Optical_Thickness_Total_LogMean =>                   &
                          cospOUT%modis_Optical_Thickness_Total_LogMean
       if (cfg%Ltauwlogmodis) modis%Optical_Thickness_Water_LogMean =>                   &
                          cospOUT%modis_Optical_Thickness_Water_LogMean
       if (cfg%Ltauilogmodis) modis%Optical_Thickness_Ice_LogMean =>                     &
                          cospOUT%modis_Optical_Thickness_Ice_LogMean
       if (cfg%Lreffclwmodis) modis%Cloud_Particle_Size_Water_Mean =>                    &
                          cospOUT%modis_Cloud_Particle_Size_Water_Mean
       if (cfg%Lreffclimodis) modis%Cloud_Particle_Size_Ice_Mean =>                      &
                          cospOUT%modis_Cloud_Particle_Size_Ice_Mean
       if (cfg%Lpctmodis)     modis%Cloud_Top_Pressure_Total_Mean =>                     &
                          cospOUT%modis_Cloud_Top_Pressure_Total_Mean
       if (cfg%Llwpmodis)     modis%Liquid_Water_Path_Mean =>                            &
                          cospOUT%modis_Liquid_Water_Path_Mean
       if (cfg%Liwpmodis)     modis%Ice_Water_Path_Mean =>                               &
                          cospOUT%modis_Ice_Water_Path_Mean
       if (cfg%Lclmodis) then
          modis%Optical_Thickness_vs_Cloud_Top_Pressure =>                               &
             cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure
          modis%Optical_Thickness_vs_ReffICE => cospOUT%modis_Optical_Thickness_vs_ReffICE
          modis%Optical_Thickness_vs_ReffLIQ => cospOUT%modis_Optical_Thickness_vs_ReffLIQ
       endif
    endif

    ! PARASOL
    if (cfg%Lparasol_sim) then
       if (cfg%Lparasolrefl) sglidar%refl        => cospOUT%parasolPix_refl
       if (cfg%Lparasolrefl) stlidar%parasolrefl => cospOUT%parasolGrid_refl
    endif

    ! RTTOV
    if (cfg%Lrttov_sim) rttov%tbs => cospOUT%rttov_tbs  

    ! CALIPSO
    if (cfg%Llidar_sim) then
       ! *NOTE* In COSPv1.5 all outputs are ordered from TOA-2-SFC, but in COSPv1.4 this is
       !        not true. To maintain the outputs of v1.4, the affected fields are flipped.

       if (cfg%LlidarBetaMol532) then
          cospOUT%calipso_beta_mol = cospOUT%calipso_beta_mol(:,sglidar%Nlevels:1:-1)
          sglidar%beta_mol         => cospOUT%calipso_beta_mol
       endif
       if (cfg%Latb532) then
          cospOUT%calipso_beta_tot = cospOUT%calipso_beta_tot(:,:,sglidar%Nlevels:1:-1)
          sglidar%beta_tot         => cospOUT%calipso_beta_tot
       endif
       if (cfg%LcfadLidarsr532)  then
          cospOUT%calipso_cfad_sr       = cospOUT%calipso_cfad_sr(:,:,stlidar%Nlevels:1:-1)
          cospOUT%calipso_betaperp_tot  = cospOUT%calipso_betaperp_tot(:,:,sglidar%Nlevels:1:-1)
          stlidar%srbval                => cospOUT%calipso_srbval
          stlidar%cfad_sr               => cospOUT%calipso_cfad_sr
          sglidar%betaperp_tot          => cospOUT%calipso_betaperp_tot
       endif   
       if (cfg%Lclcalipso) then
          cospOUT%calipso_lidarcld = cospOUT%calipso_lidarcld(:,stlidar%Nlevels:1:-1)
          stlidar%lidarcld         => cospOUT%calipso_lidarcld
       endif       
       if (cfg%Lclhcalipso .or. cfg%Lclmcalipso .or. cfg%Lcllcalipso .or. cfg%Lcltcalipso) then
          stlidar%cldlayer => cospOUT%calipso_cldlayer
       endif
       if (cfg%Lclcalipsoice .or. cfg%Lclcalipsoliq .or. cfg%Lclcalipsoun) then
          stlidar%lidarcldphase => cospOUT%calipso_lidarcldphase
       endif
       if (cfg%Lcllcalipsoice .or. cfg%Lclmcalipsoice .or. cfg%Lclhcalipsoice .or.                   &
           cfg%Lcltcalipsoice .or. cfg%Lcllcalipsoliq .or. cfg%Lclmcalipsoliq .or.                   &
           cfg%Lclhcalipsoliq .or. cfg%Lcltcalipsoliq .or. cfg%Lcllcalipsoun  .or.                   &
           cfg%Lclmcalipsoun  .or. cfg%Lclhcalipsoun  .or. cfg%Lcltcalipsoun) then       
           cospOUT%calipso_lidarcldphase = cospOUT%calipso_lidarcldphase(:,stlidar%Nlevels:1:-1,:) 
           stlidar%cldlayerphase         => cospOUT%calipso_cldlayerphase
       endif
       if (cfg%Lclcalipsotmp .or. cfg%Lclcalipsotmpliq .or. cfg%Lclcalipsoice .or. cfg%Lclcalipsotmpun) then
          stlidar%lidarcldtmp => cospOUT%calipso_lidarcldtmp
       endif
       ! Fields present, but not controlled by logical switch
       cospOUT%calipso_temp_tot = cospOUT%calipso_temp_tot(:,sglidar%Nlevels:1:-1)
       cospOUT%calipso_tau_tot  = cospOUT%calipso_tau_tot(:,:,sglidar%Nlevels:1:-1)
       sglidar%temp_tot => cospOUT%calipso_temp_tot
       sglidar%tau_tot  => cospOUT%calipso_tau_tot
    endif

    ! Cloudsat             
    if (cfg%Lradar_sim) then
       ! *NOTE* In COSPv1.5 all outputs are ordered from TOA-2-SFC, but in COSPv1.4 this is
       !        not true. To maintain the outputs of v1.4, the affected fields are flipped.    
       if (cfg%Ldbze94) then
          cospOUT%cloudsat_Ze_tot = cospOUT%cloudsat_Ze_tot(:,:,sgradar%Nlevels:1:-1) 
          sgradar%Ze_tot                => cospOUT%cloudsat_Ze_tot  
       endif
       if (cfg%LcfadDbze94) then 
          cospOUT%cloudsat_cfad_ze      = cospOUT%cloudsat_cfad_ze(:,:,stradar%Nlevels:1:-1)
          stradar%cfad_ze               => cospOUT%cloudsat_cfad_ze              
       endif
 
    endif

    ! Combined instrument products
    if (cfg%Lclcalipso2) then
       cospOUT%lidar_only_freq_cloud = cospOUT%lidar_only_freq_cloud(:,stradar%Nlevels:1:-1)
       stradar%lidar_only_freq_cloud => cospOUT%lidar_only_freq_cloud    
    endif
    if (cfg%Lcltlidarradar) stradar%radar_lidar_tcc => cospOUT%radar_lidar_tcc      

    
    ! *NOTE* In COSPv1.5 all outputs are ordered from TOA-2-SFC, but in COSPv1.4 this is
    !        not true. To maintain the outputs of v1.4, the affected fields are flipped.
    sgx%frac_out                  = sgx%frac_out(:,:,sgx%Nlevels:1:-1)

  end subroutine cosp_interface_v1p4
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE construct_cosp_gridbox_v1p4
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_COSP_GRIDBOX_v1p4(time,time_bnds,radar_freq,surface_radar,         &
                                         use_mie_tables,use_gas_abs,do_ray,melt_lay,k2,   &
                                         Npoints,Nlevels,Ncolumns,Nhydro,Nprmts_max_hydro,&
                                         Naero,Nprmts_max_aero,Npoints_it,lidar_ice_type, &
                                         isccp_top_height,isccp_top_height_direction,     &
                                         isccp_overlap,isccp_emsfc_lw,                    &
                                         use_precipitation_fluxes,use_reff,Plat,Sat,Inst, &
                                         Nchan,ZenAng,Ichan,SurfEm,co2,ch4,n2o,co,        &
                                         y,                                               &
                                         lon,lat,p, ph,zlev,zlev_half,T,rh,sh,cca,tca,skt,&
                                         landmask,mr_ozone,u_wind,v_wind,sunlit,fl_lsrain,&
                                         fl_lssnow,fl_lsgrpl,fl_ccrain,fl_ccsnow,dtau_s,  &
                                         dtau_c,dem_s,dem_c,Reff,mr_lsliq,mr_lsice,       &
                                         mr_ccliq,mr_ccice,load_LUT)
    
    ! Inputs
    double precision,intent(in) :: &
         time,          & ! Time since start of run [days] 
         time_bnds(2)     ! Time boundaries
    integer,intent(in) :: &
         surface_radar,     & ! surface=1,spaceborne=0
         use_mie_tables,    & ! use a precomputed lookup table? yes=1,no=0,2=use first
                              ! column everywhere
         use_gas_abs,       & ! include gaseous absorption? yes=1,no=0
         do_ray,            & ! calculate/output Rayleigh refl=1, not=0
         melt_lay,          & ! melting layer model off=0, on=1
         Npoints,           & ! Number of gridpoints
         Nlevels,           & ! Number of levels
         Ncolumns,          & ! Number of columns
         Nhydro,            & ! Number of hydrometeors
         Nprmts_max_hydro,  & ! Max number of parameters for hydrometeor size 
                              ! distributions
         Naero,             & ! Number of aerosol species
         Nprmts_max_aero,   & ! Max number of parameters for aerosol size distributions
         Npoints_it,        & ! Number of gridpoints processed in one iteration
         lidar_ice_type,    & ! Ice particle shape in lidar calculations (0=ice-spheres ;
                              ! 1=ice-non-spherical)
         isccp_top_height , & !
         isccp_top_height_direction, & !
         isccp_overlap,     & !
         Plat,              & ! RTTOV satellite platform
         Sat,               & ! RTTOV satellite
         Inst,              & ! RTTOV instrument
         Nchan                ! RTTOV number of channels
    integer,intent(in),dimension(Nchan) :: &
         Ichan
    real(wp),intent(in) :: &
         radar_freq,       & ! Radar frequency [GHz]
         k2,               & ! |K|^2, -1=use frequency dependent default
         isccp_emsfc_lw,   & ! 11microm surface emissivity
         co2,              & ! CO2 
         ch4,              & ! CH4
         n2o,              & ! N2O
         co,               & ! CO
         ZenAng              ! RTTOV zenith abgle
    real(wp),intent(in),dimension(Nchan) :: &
         SurfEm
    logical,intent(in) :: &
         use_precipitation_fluxes,&
         use_reff
    logical,intent(in),optional :: load_LUT
    
    ! RTTOV inputs
    real(wp),intent(in),dimension(Npoints),optional :: &
         lon,        &
         lat,        &
         u_wind,     &
         v_wind,     &
         skt,        &
         landmask,   &
         sunlit   
    real(wp),intent(in),dimension(Npoints,Nlevels),optional :: &
         p,          &
         ph,         &
         zlev,       &
         T,          & 
         rh,         &
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
         dem_s,      &
         mr_lsice,   &
         mr_lsliq,   &
         mr_ccliq,   &
         mr_ccice
    real(wp),intent(in),dimension(Npoints,Nlevels+1),optional ::&
         zlev_half
    real(wp),intent(in),dimension(Npoints,Nlevels,N_HYDRO),optional :: &
         Reff    
    
    ! Outputs
    type(cosp_gridbox_v1p4),intent(out) :: y
    
    ! local variables
    integer :: k
    character(len=240) :: LUT_file_name
    logical :: local_load_LUT,rttovInputs
    
    if (present(load_LUT)) then
       local_load_LUT = load_LUT
    else
       local_load_LUT = RADAR_SIM_LOAD_scale_LUTs_flag
    endif
    
    ! Check to see if RTTOV inputs are provided (default is no RTTOV)
    rttovInputs = .false.
    if (present(lon)       .and. present(lat)       .and.                                &
        present(ph)        .and. present(zlev)      .and. present(zlev_half) .and.       &
        present(rh)        .and. present(sh)        .and. present(cca)       .and.       &
        present(p)         .and. present(T)         .and. present(tca)       .and.       &
        present(skt)       .and. present(landmask)  .and. present(mr_ozone)  .and.       &
        present(u_wind)    .and. present(v_wind)    .and. present(sunlit)    .and.       &
        present(fl_lsrain) .and. present(fl_lssnow) .and. present(fl_lsgrpl) .and.       &
        present(fl_ccrain) .and. present(fl_ccsnow) .and. present(dtau_s)    .and.       &
        present(dtau_c)    .and. present(dem_s)     .and. present(dem_c)     .and.       &
        present(Reff)      .and. present(mr_lsliq)  .and. present(mr_lsice)  .and.       &
        present(mr_ccliq)  .and. present(mr_ccice)) then 
        rttovInputs = .true.
    endif
    
    ! Dimensions and scalars
    y%radar_freq       = radar_freq
    y%surface_radar    = surface_radar
    y%use_mie_tables   = use_mie_tables
    y%use_gas_abs      = use_gas_abs
    y%do_ray           = do_ray
    y%melt_lay         = melt_lay
    y%k2               = k2
    y%Npoints          = Npoints
    y%Nlevels          = Nlevels
    y%Ncolumns         = Ncolumns
    y%Nhydro           = Nhydro
    y%Nprmts_max_hydro = Nprmts_max_hydro
    y%Naero            = Naero
    y%Nprmts_max_aero  = Nprmts_max_aero
    y%Npoints_it       = Npoints_it
    y%lidar_ice_type   = lidar_ice_type
    y%isccp_top_height = isccp_top_height
    y%isccp_top_height_direction = isccp_top_height_direction
    y%isccp_overlap    = isccp_overlap
    y%isccp_emsfc_lw   = isccp_emsfc_lw
    y%use_precipitation_fluxes = use_precipitation_fluxes
    y%use_reff = use_reff
    y%time      = time
    y%time_bnds = time_bnds
    
    ! RTTOV parameters
    y%Plat   = Plat
    y%Sat    = Sat
    y%Inst   = Inst
    y%Nchan  = Nchan
    y%ZenAng = ZenAng
    y%co2    = co2
    y%ch4    = ch4
    y%n2o    = n2o
    y%co     = co
    
    ! Gridbox information (Npoints,Nlevels)
    allocate(y%zlev(Npoints,Nlevels),y%zlev_half(Npoints,Nlevels),                       &
             y%dlev(Npoints,Nlevels),y%p(Npoints,Nlevels),y%ph(Npoints,Nlevels),         &
             y%T(Npoints,Nlevels),y%q(Npoints,Nlevels), y%sh(Npoints,Nlevels),           &
             y%dtau_s(Npoints,Nlevels),y%dtau_c(Npoints,Nlevels),                        &
             y%dem_s(Npoints,Nlevels),y%dem_c(Npoints,Nlevels),y%tca(Npoints,Nlevels),   &
             y%cca(Npoints,Nlevels),y%rain_ls(Npoints,Nlevels),                          &
             y%rain_cv(Npoints,Nlevels),y%grpl_ls(Npoints,Nlevels),                      &
             y%snow_ls(Npoints,Nlevels),y%snow_cv(Npoints,Nlevels),                      &
             y%mr_ozone(Npoints,Nlevels))
    
    ! Surface information and geolocation (Npoints)
    allocate(y%toffset(Npoints),y%longitude(Npoints),y%latitude(Npoints),y%psfc(Npoints),&
             y%land(Npoints),y%sunlit(Npoints),y%skt(Npoints),y%u_wind(Npoints),         &
             y%v_wind(Npoints))
    
    ! Hydrometeors concentration and distribution parameters
    allocate(y%mr_hydro(Npoints,Nlevels,Nhydro),y%Reff(Npoints,Nlevels,Nhydro),          &
             y%dist_prmts_hydro(Nprmts_max_hydro,Nhydro),y%Np(Npoints,Nlevels,Nhydro)) 

    ! Aerosols concentration and distribution parameters
    allocate(y%conc_aero(Npoints,Nlevels,Naero), y%dist_type_aero(Naero), &
             y%dist_prmts_aero(Npoints,Nlevels,Nprmts_max_aero,Naero))
    
    ! RTTOV channels and sfc. emissivity
    allocate(y%ichan(Nchan),y%surfem(Nchan))
    
    ! Initialize    
    y%zlev      = 0.0
    y%zlev_half = 0.0
    y%dlev      = 0.0
    y%p         = 0.0
    y%ph        = 0.0
    y%T         = 0.0
    y%q         = 0.0
    y%sh        = 0.0
    y%dtau_s    = 0.0
    y%dtau_c    = 0.0
    y%dem_s     = 0.0
    y%dem_c     = 0.0
    y%tca       = 0.0
    y%cca       = 0.0
    y%rain_ls   = 0.0
    y%rain_cv   = 0.0
    y%grpl_ls   = 0.0
    y%snow_ls   = 0.0
    y%snow_cv   = 0.0
    y%Reff      = 0.0
    y%Np        = 0.0 
    y%mr_ozone  = 0.0
    y%u_wind    = 0.0
    y%v_wind    = 0.0
    y%toffset   = 0.0
    y%longitude = 0.0
    y%latitude  = 0.0
    y%psfc      = 0.0
    y%land      = 0.0
    y%sunlit    = 0.0
    y%skt       = 0.0
    y%mr_hydro  = 0.0
    y%dist_prmts_hydro = 0.0 
    y%conc_aero        = 0.0 
    y%dist_type_aero   = 0   
    y%dist_prmts_aero  = 0.0 

    ! Toffset. This assumes that time is the mid-point of the interval.
    do k=1,Npoints
       y%toffset(k) = -0.5_wp*3._wp/24._wp + 3._wp/24._wp*(k-0.5)/Npoints
    enddo
    
    ! Setup RTTOV inputs (if provided)
    if (rttovInputs) then
       y%longitude = lon
       y%latitude  = lat
       y%ichan     = ichan
       y%surfem    = surfem
       y%p         = p
       y%ph        = ph
       y%zlev      = zlev
       y%zlev_half = zlev_half
       y%T         = T 
       y%q         = rh
       y%sh        = sh
       y%cca       = cca
       y%tca       = tca
       y%psfc      = ph(:,1)
       y%skt       = skt
       y%land      = landmask
       y%mr_ozone  = mr_ozone
       y%u_wind    = u_wind
       y%v_wind    = v_wind
       y%sunlit    = sunlit
       y%rain_ls   = fl_lsrain
       y%snow_ls   = fl_lssnow
       y%grpl_ls   = fl_lsgrpl
       y%rain_cv   = fl_ccrain
       y%snow_cv   = fl_ccsnow
       y%dtau_s    = dtau_s
       y%dtau_c    = dtau_c
       y%dem_s     = dem_s
       y%dem_c     = dem_c
       y%mr_hydro(:,:,I_LSCLIQ) = mr_lsliq
       y%mr_hydro(:,:,I_LSCICE) = mr_lsice
       y%mr_hydro(:,:,I_CVCLIQ) = mr_ccliq
       y%mr_hydro(:,:,I_CVCICE) = mr_ccice
       y%Reff = Reff
       y%Reff(:,:,I_LSRAIN) = 0.0
    endif      
    
  END SUBROUTINE CONSTRUCT_COSP_GRIDBOX_v1p4
    
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE destroy_cosp_gridbox_v1p4
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE destroy_cosp_gridbox_v1p4(y,save_LUT)
    
    type(cosp_gridbox_v1p4),intent(inout) :: y
    logical,intent(in),optional :: save_LUT
    
    logical :: local_save_LUT
    if (present(save_LUT)) then
       local_save_LUT = save_LUT
    else
       local_save_LUT = RADAR_SIM_UPDATE_scale_LUTs_flag
    endif
    
    ! save any updates to radar simulator LUT
    if (local_save_LUT) call save_scale_LUTs_v1p4(y%hp)
    
    deallocate(y%zlev,y%zlev_half,y%dlev,y%p,y%ph,y%T,y%q,y%sh,y%dtau_s,y%dtau_c,y%dem_s,&
               y%dem_c,y%toffset,y%longitude,y%latitude,y%psfc,y%land,y%tca,y%cca,       &
               y%mr_hydro,y%dist_prmts_hydro,y%conc_aero,y%dist_type_aero,               &
               y%dist_prmts_aero,y%rain_ls,y%rain_cv,y%snow_ls,y%snow_cv,y%grpl_ls,      &
               y%sunlit,y%skt,y%Reff,y%Np,y%ichan,y%surfem,y%mr_ozone,y%u_wind,y%v_wind)
    
  END SUBROUTINE destroy_cosp_gridbox_v1p4
  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE save_scale_LUTs_v1p4
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine save_scale_LUTs_v1p4(hp)
    type(class_param), intent(inout) :: hp
    logical                          :: LUT_file_exists
    integer                          :: i,j,k,ind
    
    inquire(file=trim(hp%scale_LUT_file_name) // '_radar_Z_scale_LUT.dat', &
         exist=LUT_file_exists)
    
    OPEN(unit=12,file=trim(hp%scale_LUT_file_name) // '_radar_Z_scale_LUT.dat',&
         form='unformatted',err= 99,access='DIRECT',recl=28)
    
    write(*,*) 'Creating or Updating radar LUT file: ', &
         trim(hp%scale_LUT_file_name) // '_radar_Z_scale_LUT.dat'
    
    do i=1,maxhclass
       do j=1,mt_ntt
          do k=1,nRe_types
             ind = i+(j-1)*maxhclass+(k-1)*(nRe_types*mt_ntt)
             if(.not.LUT_file_exists .or. hp%Z_scale_added_flag(i,j,k)) then
                hp%Z_scale_added_flag(i,j,k)=.false.
                write(12,rec=ind) hp%Z_scale_flag(i,j,k), &
                     hp%Ze_scaled(i,j,k), &
                     hp%Zr_scaled(i,j,k), &
                     hp%kr_scaled(i,j,k)
             endif
          enddo
       enddo
    enddo
    close(unit=12)
    return 
    
99  write(*,*) 'Error: Unable to create/update radar LUT file: ', &
         trim(hp%scale_LUT_file_name) // '_radar_Z_scale_LUT.dat'
    return  
  end subroutine save_scale_LUTs_v1p4

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !SUBROUTINE construct_cosp_vgrid_v1p4
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_COSP_VGRID_v1p4(gbx,Nlvgrid,use_vgrid,cloudsat,x)
    type(cosp_gridbox_v1p4),intent(in) :: gbx ! Gridbox information
    integer,intent(in) :: Nlvgrid  ! Number of new levels    
    logical,intent(in) :: use_vgrid! Logical flag that controls the output on a different grid
    logical,intent(in) :: cloudsat ! TRUE if a CloudSat like grid (480m) is requested
    type(cosp_vgrid),intent(out) :: x
    
    ! Local variables
    integer :: i
    real :: zstep
    
    x%use_vgrid  = use_vgrid
    x%csat_vgrid = cloudsat
    
    ! Dimensions
    x%Npoints  = gbx%Npoints
    x%Ncolumns = gbx%Ncolumns
    x%Nlevels  = gbx%Nlevels
    
    ! --- Allocate arrays ---
    if (use_vgrid) then
       x%Nlvgrid = Nlvgrid
    else 
       x%Nlvgrid = gbx%Nlevels
    endif
    allocate(x%z(x%Nlvgrid),x%zl(x%Nlvgrid),x%zu(x%Nlvgrid))
    allocate(x%mz(x%Nlevels),x%mzl(x%Nlevels),x%mzu(x%Nlevels))
    
    ! --- Model vertical levels ---
    ! Use height levels of first model gridbox
    x%mz  = gbx%zlev(1,:)
    x%mzl = gbx%zlev_half(1,:)
    x%mzu(1:x%Nlevels-1) = gbx%zlev_half(1,2:x%Nlevels)
    x%mzu(x%Nlevels) = gbx%zlev(1,x%Nlevels) + (gbx%zlev(1,x%Nlevels) - x%mzl(x%Nlevels))
    
    if (use_vgrid) then
       ! --- Initialise to zero ---
       x%z  = 0.0
       x%zl = 0.0
       x%zu = 0.0
       if (cloudsat) then ! --- CloudSat grid requested ---
          zstep = 480.0
       else
          ! Other grid requested. Constant vertical spacing with top at 20 km
          zstep = 20000.0/x%Nlvgrid
       endif
       do i=1,x%Nlvgrid
          x%zl(i) = (i-1)*zstep
          x%zu(i) = i*zstep
       enddo
       x%z = (x%zl + x%zu)/2.0
    else
       x%z  = x%mz
       x%zl = x%mzl
       x%zu = x%mzu
    endif
    
  END SUBROUTINE CONSTRUCT_COSP_VGRID_v1p4
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE construct_cosp_sgradar
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine construct_cosp_sgradar(Npoints,Ncolumns,Nlevels,Nhydro,x)
    integer,target,     intent(in)  :: Npoints  ! Number of sampled points
    integer,target,     intent(in)  :: Ncolumns ! Number of subgrid columns
    integer,target,     intent(in)  :: Nlevels  ! Number of model levels
    integer,target,     intent(in)  :: Nhydro   ! Number of hydrometeors
    type(cosp_sgradar), intent(out) :: x

    ! Dimensions
    x%Npoints  => Npoints
    x%Ncolumns => Ncolumns
    x%Nlevels  => Nlevels
    x%Nhydro   => Nhydro

    ! Allocate
    allocate(x%att_gas(Npoints,Nlevels),x%Ze_tot(Npoints,Ncolumns,Nlevels))

    ! Initialize
    x%att_gas = 0._wp
    x%Ze_tot  = 0._wp

  end subroutine construct_cosp_sgradar
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE construct_cosp_radarstats
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine construct_cosp_radarstats(Npoints,Ncolumns,Nlevels,Nhydro,x)
    integer,target,       intent(in)  :: Npoints  ! Number of sampled points
    integer,target,       intent(in)  :: Ncolumns ! Number of subgrid columns
    integer,target,       intent(in)  :: Nlevels  ! Number of model levels
    integer,target,       intent(in)  :: Nhydro   ! Number of hydrometeors
    type(cosp_radarstats),intent(out) :: x

    ! Dimensions
    x%Npoints  => Npoints
    x%Ncolumns => Ncolumns
    x%Nlevels  => Nlevels
    x%Nhydro   => Nhydro

    ! Allocate
    allocate(x%cfad_ze(Npoints,DBZE_BINS,Nlevels),x%lidar_only_freq_cloud(Npoints,Nlevels), &
             x%radar_lidar_tcc(Npoints))
    
    ! Initialize
    x%cfad_ze               = 0._wp
    x%lidar_only_freq_cloud = 0._wp
    x%radar_lidar_tcc       = 0._wp    
    
  end subroutine construct_cosp_radarstats
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE destroy_cosp_sgradar
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine destroy_cosp_sgradar(x)
    type(cosp_sgradar),intent(inout) :: x

    deallocate(x%att_gas,x%Ze_tot)

  end subroutine destroy_cosp_sgradar
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE destroy_cosp_radarstats
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine destroy_cosp_radarstats(x)
    type(cosp_radarstats),intent(inout) :: x

    deallocate(x%cfad_ze,x%lidar_only_freq_cloud,x%radar_lidar_tcc)

  end subroutine destroy_cosp_radarstats
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE construct_cosp_sglidar
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine construct_cosp_sglidar(Npoints,Ncolumns,Nlevels,Nhydro,Nrefl,x)
    ! Inputs
    integer,intent(in),target :: &
         Npoints,  & ! Number of sampled points
         Ncolumns, & ! Number of subgrid columns
         Nlevels,  & ! Number of model levels
         Nhydro,   & ! Number of hydrometeors
         Nrefl       ! Number of parasol reflectances ! parasol
    ! Outputs
    type(cosp_sglidar),intent(out) :: x

    ! Dimensions
    x%Npoints  => Npoints
    x%Ncolumns => Ncolumns
    x%Nlevels  => Nlevels
    x%Nhydro   => Nhydro
    x%Nrefl    => Nrefl

    ! Allocate
    allocate(x%beta_mol(x%Npoints,x%Nlevels), x%beta_tot(x%Npoints,x%Ncolumns,x%Nlevels), &
             x%tau_tot(x%Npoints,x%Ncolumns,x%Nlevels),x%refl(x%Npoints,x%Ncolumns,x%Nrefl), &
             x%temp_tot(x%Npoints,x%Nlevels),x%betaperp_tot(x%Npoints,x%Ncolumns,x%Nlevels))

    ! Initialize
    x%beta_mol     = 0._wp
    x%beta_tot     = 0._wp
    x%tau_tot      = 0._wp
    x%refl         = 0._wp
    x%temp_tot     = 0._wp
    x%betaperp_tot = 0._wp
  end subroutine construct_cosp_sglidar
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE construct_cosp_lidarstats
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine construct_cosp_lidarstats(Npoints,Ncolumns,Nlevels,Nhydro,Nrefl,x)
    ! Inputs
    integer,intent(in),target :: &
         Npoints,  & ! Number of sampled points
         Ncolumns, & ! Number of subgrid columns
         Nlevels,  & ! Number of model levels
         Nhydro,   & ! Number of hydrometeors
         Nrefl       ! Number of parasol reflectances
    ! Outputs
    type(cosp_lidarstats),intent(out) :: x
    ! Local variables
    integer :: i,j,k,l,m

    ! Dimensions
    x%Npoints  => Npoints
    x%Ncolumns => Ncolumns
    x%Nlevels  => Nlevels
    x%Nhydro   => Nhydro
    x%Nrefl    => Nrefl

    ! Allocate
    allocate(x%srbval(SR_BINS),x%cfad_sr(x%Npoints,SR_BINS,x%Nlevels), &
         x%lidarcld(x%Npoints,x%Nlevels), x%cldlayer(x%Npoints,LIDAR_NCAT),&
         x%parasolrefl(x%Npoints,x%Nrefl),x%lidarcldphase(x%Npoints,x%Nlevels,6),&
         x%lidarcldtmp(x%Npoints,LIDAR_NTEMP,5),x%cldlayerphase(x%Npoints,LIDAR_NCAT,6))

    ! Initialize
    x%srbval        = 0._wp
    x%cfad_sr       = 0._wp
    x%lidarcld      = 0._wp
    x%cldlayer      = 0._wp
    x%parasolrefl   = 0._wp
    x%lidarcldphase = 0._wp
    x%cldlayerphase = 0._wp
    x%lidarcldtmp   = 0._wp

  end subroutine construct_cosp_lidarstats

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE destroy_cosp_lidarstats
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine destroy_cosp_lidarstats(x)
    type(cosp_lidarstats),intent(inout) :: x

    deallocate(x%srbval,x%cfad_sr,x%lidarcld,x%cldlayer,x%parasolrefl,x%cldlayerphase,   &
               x%lidarcldtmp,x%lidarcldphase)

  end subroutine destroy_cosp_lidarstats
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE destroy_cosp_sglidar
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine destroy_cosp_sglidar(x)
    type(cosp_sglidar),intent(inout) :: x

    deallocate(x%beta_mol,x%beta_tot,x%tau_tot,x%refl,x%temp_tot,x%betaperp_tot)
  end subroutine destroy_cosp_sglidar
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !                           SUBROUTINE construct_cosp_isccp
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_COSP_ISCCP(Npoints,Ncolumns,Nlevels,x)
    integer,target,   intent(in)  :: Npoints  ! Number of sampled points
    integer,target,   intent(in)  :: Ncolumns ! Number of subgrid columns
    integer,target,   intent(in)  :: Nlevels  ! Number of model levels
    type(cosp_isccp), intent(out) :: x        ! Output

    x%Npoints  => Npoints
    x%Ncolumns => Ncolumns
    x%Nlevels  => Nlevels
    x%Npoints  => Npoints
    x%Ncolumns => Ncolumns
    x%Nlevels  => Nlevels

    ! Allocate 
    allocate(x%fq_isccp(Npoints,7,7),x%totalcldarea(Npoints),x%meanptop(Npoints),        &
             x%meantaucld(Npoints),x%meantb(Npoints),x%meantbclr(Npoints),               &
             x%meanalbedocld(Npoints),x%boxtau(Npoints,Ncolumns),                        &
             x%boxptop(Npoints,Ncolumns))

    ! Initialize
    x%fq_isccp     = 0._wp
    x%totalcldarea = 0._wp
    x%meanptop     = 0._wp
    x%meantaucld   = 0._wp
    x%meantb       = 0._wp
    x%meantbclr    = 0._wp
    x%meanalbedocld= 0._wp
    x%boxtau       = 0._wp
    x%boxptop      = 0._wp

  END SUBROUTINE CONSTRUCT_COSP_ISCCP

 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 !                          SUBROUTINE destroy_cosp_isccp
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE destroy_cosp_isccp(x)
    type(cosp_isccp),intent(inout) :: x
    
    deallocate(x%fq_isccp,x%totalcldarea,x%meanptop,x%meantaucld,x%meantb,x%meantbclr,   &
               x%meanalbedocld,x%boxtau,x%boxptop)
  END SUBROUTINE destroy_cosp_isccp

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !  					SUBROUTINE construct_cosp_misr
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_COSP_MISR(Npoints,x)
    integer,          intent(in),target   :: Npoints  ! Number of gridpoints
    type(cosp_misr),  intent(out)         :: x

    ! Local variables
    integer,target :: &
         Ntau=7,Ncth=numMISRHgtBins
    
    x%Npoints => Npoints
    x%Ntau    => Ntau
    x%Nlevels => Ncth

    ! Allocate
    allocate(x%fq_MISR(x%Npoints,x%Ntau,x%Nlevels),x%MISR_meanztop(x%Npoints),           &
             x%MISR_cldarea(x%Npoints),x%MISR_dist_model_layertops(x%Npoints,x%Nlevels))

    ! Initialize
    x%fq_MISR                   = 0._wp
    x%MISR_meanztop             = 0._wp
    x%MISR_cldarea              = 0._wp
    x%MISR_dist_model_layertops = 0._wp
   
  END SUBROUTINE CONSTRUCT_COSP_MISR
 
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 !                           SUBROUTINE destroy_cosp_misr
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE destroy_cosp_misr(x)
    type(cosp_misr),intent(inout) :: x

    deallocate(x%fq_MISR,x%MISR_meanztop,x%MISR_cldarea,x%MISR_dist_model_layertops)

  END SUBROUTINE destroy_cosp_misr
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE construct_cosp_modis
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_COSP_MODIS(nPoints, x)
    integer,target,   intent(in)  :: Npoints  ! Number of sampled points
    type(cosp_MODIS), intent(out) :: x
    
    x%nPoints  => nPoints
    
    ! Allocate gridmean variables
    allocate(x%Cloud_Fraction_Total_Mean(Npoints),x%Cloud_Fraction_Water_Mean(Npoints),  &
             x%Cloud_Fraction_Ice_Mean(Npoints),x%Cloud_Fraction_High_Mean(Npoints),     &
             x%Cloud_Fraction_Mid_Mean(Npoints),x%Cloud_Fraction_Low_Mean(Npoints),      &
             x%Optical_Thickness_Total_Mean(Npoints),                                    &
             x%Optical_Thickness_Water_Mean(Npoints),                                    &
             x%Optical_Thickness_Ice_Mean(Npoints),                                      &
             x%Optical_Thickness_Total_LogMean(Npoints),                                 &
             x%Optical_Thickness_Water_LogMean(Npoints),                                 &
             x%Optical_Thickness_Ice_LogMean(Npoints),                                   &
             x%Cloud_Particle_Size_Water_Mean(Npoints),                                  &
             x%Cloud_Particle_Size_Ice_Mean(Npoints),                                    &
             x%Cloud_Top_Pressure_Total_Mean(Npoints),x%Liquid_Water_Path_Mean(Npoints), &
             x%Ice_Water_Path_Mean(Npoints),                                             &
             x%Optical_Thickness_vs_Cloud_Top_Pressure(nPoints,numMODISTauBins+1,numMODISPresBins),&
             x%Optical_Thickness_vs_ReffICE(nPoints,numModisTauBins+1,numMODISReffIceBins),&
             x%Optical_Thickness_vs_ReffLIQ(nPoints,numModisTauBins+1,numMODISReffLiqBins))
    x%Optical_Thickness_vs_Cloud_Top_Pressure(:, :, :) = R_UNDEF
    x%Optical_Thickness_vs_ReffICE(:,:,:)              = R_UNDEF
    x%Optical_Thickness_vs_ReffLIQ(:,:,:)              = R_UNDEF

  END SUBROUTINE CONSTRUCT_COSP_MODIS
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE destroy_cosp_modis
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE destroy_cosp_modis(x)
    type(cosp_MODIS),intent(inout) :: x
    
    ! Free space used by cosp_modis variable.     
    if(associated(x%Cloud_Fraction_Total_Mean))  deallocate(x%Cloud_Fraction_Total_Mean) 
    if(associated(x%Cloud_Fraction_Water_Mean))  deallocate(x%Cloud_Fraction_Water_Mean) 
    if(associated(x%Cloud_Fraction_Ice_Mean))    deallocate(x%Cloud_Fraction_Ice_Mean) 
    if(associated(x%Cloud_Fraction_High_Mean))   deallocate(x%Cloud_Fraction_High_Mean) 
    if(associated(x%Cloud_Fraction_Mid_Mean))    deallocate(x%Cloud_Fraction_Mid_Mean) 
    if(associated(x%Cloud_Fraction_Low_Mean))    deallocate(x%Cloud_Fraction_Low_Mean) 
    if(associated(x%Liquid_Water_Path_Mean))     deallocate(x%Liquid_Water_Path_Mean) 
    if(associated(x%Ice_Water_Path_Mean))        deallocate(x%Ice_Water_Path_Mean)
    if(associated(x%Optical_Thickness_Total_Mean))                                       &
         deallocate(x%Optical_Thickness_Total_Mean) 
    if(associated(x%Optical_Thickness_Water_Mean))                                       &
         deallocate(x%Optical_Thickness_Water_Mean) 
    if(associated(x%Optical_Thickness_Ice_Mean))                                         &
         deallocate(x%Optical_Thickness_Ice_Mean) 
    if(associated(x%Optical_Thickness_Total_LogMean))                                    &
         deallocate(x%Optical_Thickness_Total_LogMean) 
    if(associated(x%Optical_Thickness_Water_LogMean))                                    &
         deallocate(x%Optical_Thickness_Water_LogMean) 
    if(associated(x%Optical_Thickness_Ice_LogMean))                                      &
         deallocate(x%Optical_Thickness_Ice_LogMean) 
    if(associated(x%Cloud_Particle_Size_Water_Mean))                                     &
         deallocate(x%Cloud_Particle_Size_Water_Mean) 
    if(associated(x%Cloud_Particle_Size_Ice_Mean))                                       &
         deallocate(x%Cloud_Particle_Size_Ice_Mean) 
    if(associated(x%Cloud_Top_Pressure_Total_Mean))                                      &
         deallocate(x%Cloud_Top_Pressure_Total_Mean) 
    if(associated(x%Optical_Thickness_vs_Cloud_Top_Pressure))                            &
         deallocate(x%Optical_Thickness_vs_Cloud_Top_Pressure) 
    if(associated(x%Optical_Thickness_vs_ReffICE))                                       &
         deallocate(x%Optical_Thickness_vs_ReffICE) 
    if(associated(x%Optical_Thickness_vs_ReffLIQ))                                       &
         deallocate(x%Optical_Thickness_vs_ReffLIQ) 
  END SUBROUTINE destroy_cosp_modis  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !           					 SUBROUTINE construct_cosp_rttov
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_COSP_RTTOV(Npoints,Nchan,x)
    integer,          intent(in)  :: Npoints  ! Number of sampled points
    integer,          intent(in)  :: Nchan    ! Number of channels
    type(cosp_rttov), intent(out) :: x
    
    ! Local variables
    integer :: i,j
   
    ! Allocate
    allocate(x%tbs(Npoints,Nchan))
    
    ! Initialize
    x%tbs     = 0.0
  END SUBROUTINE CONSTRUCT_COSP_RTTOV
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !                             SUBROUTINE destroy_cosp_rttov
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE destroy_cosp_rttov(x)
    type(cosp_rttov),intent(inout) :: x
    
    ! Deallocate
    deallocate(x%tbs)
  END SUBROUTINE destroy_cosp_rttov
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !                            SUBROUTINE destroy_cosp_
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine destroy_cosp_vgrid(x)
    type(cosp_vgrid),intent(inout) :: x
    deallocate(x%z, x%zl, x%zu, x%mz, x%mzl, x%mzu)
  end subroutine destroy_cosp_vgrid


    
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !                                    END MODULE
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module MOD_COSP_INTERFACE_v1p4
