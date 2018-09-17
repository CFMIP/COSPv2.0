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
#include "cosp_defs.h"
MODULE MOD_COSP_INTERFACE_v1p4
  use COSP_KINDS,          only: wp,dp
  use cosp_phys_constants, only: amw,amd,amO3,amCO2,amCH4,amN2O,amCO
  use MOD_COSP,            only: cosp_init,cosp_outputs,cosp_optical_inputs,              &
                                 cosp_column_inputs,cosp_simulator,linitialization
  use mod_cosp_config,     only: RTTOV_MAX_CHANNELS,N_HYDRO,numMODISTauBins,modis_histTau,&
                                 modis_histTauEdges,modis_histTauCenters,ntau,ntauV1p4,   &
                                 tau_binBounds,tau_binEdges,tau_binCenters,R_UNDEF,       &
                                 tau_binBoundsV1p4,tau_binEdgesV1p4,tau_binCentersV1p4,   &
                                 numMISRHgtBins,SR_BINS,LIDAR_NCAT,LIDAR_NTEMP,cloudsat_DBZE_BINS, &
                                 numMODISReffIceBins, numMODISPresBins,PARASOL_NREFL,     &
                                 numMODISReffLiqBins,vgrid_zl,vgrid_zu,vgrid_z,           &
                                 numISCCPTauBins,numISCCPPresBins,numMISRTauBins
  use mod_quickbeam_optics,only: size_distribution,hydro_class_init,quickbeam_optics_init,&
                                 quickbeam_optics,gases
  use cosp_optics,         only: cosp_simulator_optics,lidar_optics,modis_optics,         &
                                 modis_optics_partition
  use quickbeam,           only: maxhclass,nRe_types,nd,mt_ntt,radar_cfg
  use mod_rng,             only: rng_state, init_rng
  use mod_scops,           only: scops
  use mod_prec_scops,      only: prec_scops
  use mod_cosp_utils,      only: cosp_precip_mxratio

  implicit none
  
  ! Outputs from COSP2
  type(cosp_outputs) :: cospOUT  ! NEW derived type output that contains all 
                                        ! simulator information
  
  character(len=120),parameter :: &
       RADAR_SIM_LUT_DIRECTORY = './'
  logical,parameter :: &
       RADAR_SIM_LOAD_scale_LUTs_flag   = .false., &
       RADAR_SIM_UPDATE_scale_LUTs_flag = .false.
  
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
  
  ! Initialization fields
  type(size_distribution) :: &
       sd                ! Hydrometeor description
  type(radar_cfg) :: &
       rcfg_cloudsat     ! Radar configuration
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! TYPE COSP_CONFIG
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
  ! TYPE cosp_vgrid
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
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! TYPE COSP_SUBGRID
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  TYPE COSP_SUBGRID
     integer ::      &
          Npoints,   & ! Number of gridpoints
          Ncolumns,  & ! Number of columns
          Nlevels,   & ! Number of levels
          Nhydro       ! Number of hydrometeor types
     real(wp),dimension(:,:,:),pointer :: &
          prec_frac, & ! Subgrid precip array (Npoints,Ncolumns,Nlevels)
          frac_out     ! Subgrid cloud array  (Npoints,Ncolumns,Nlevels)
  END TYPE COSP_SUBGRID  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! TYPE class_param
  ! With the reorganizing of COSPv2.0, this derived type
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
  ! TYPE cosp_gridbox
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  TYPE cosp_gridbox
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
  END TYPE cosp_gridbox
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
  ! TYPE cosp_rttov
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
 ! TYPE cosp_isccp
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
          cfad_ze               ! Ze CFAD(Npoints,cloudsat_dBZe_bins,Nlevels)
     real(wp),dimension(:),pointer :: &
          radar_lidar_tcc       ! Radar&lidar total cloud amount, grid-box scale (Npoints)
     real(wp), dimension(:,:),pointer :: &
          lidar_only_freq_cloud !(Npoints,Nlevels)
  end type cosp_radarstats
      
contains
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !                            SUBROUTINE COSP_INTERFACE (v1.4)
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine cosp_interface_v1p4(overlap,Ncolumns,cfg,vgrid,gbx,sgx,sgradar,sglidar,   &
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
    type(cosp_gridbox),intent(inout),target :: gbx ! COSP gridbox type from v1.4
                                                          ! Shares memory with new type
 
    ! Outputs from COSP2
    !type(cosp_outputs),target :: cospOUT  ! NEW derived type output that contains all 
    					                  ! simulator information
    ! Local variables
    integer :: i
    integer :: &
         num_chunks, & ! Number of iterations to make
         start_idx,  & ! Starting index when looping over points
         end_idx,    & ! Ending index when looping over points
         Nptsperit     ! Number of points for current iteration				                  
    logical :: &
         lsingle=.true., & ! True if using MMF_v3_single_moment CLOUDSAT microphysical scheme (default)
         ldouble=.false.   ! True if using MMF_v3.5_two_moment CLOUDSAT microphysical scheme  
    type(cosp_optical_inputs) :: &
         cospIN            ! COSP optical (or derived?) fields needed by simulators
    type(cosp_column_inputs) :: &
         cospstateIN       ! COSP model fields needed by simulators
    character(len=256),dimension(100) :: cosp_status

#ifdef MMF_V3_SINGLE_MOMENT    					  
    character(len=64) :: &
         cloudsat_micro_scheme = 'MMF_v3_single_moment'
#endif
#ifdef MMF_V3p5_TWO_MOMENT
    character(len=64) :: &
         cloudsat_micro_scheme = 'MMF_v3.5_two_moment'
#endif 
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Initialize COSP
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Initialize MODIS optical-depth bin boundaries for joint-histogram. (defined in cosp_config.F90)
    if (.not. allocated(modis_histTau)) then
       allocate(modis_histTau(ntauV1p4+1),modis_histTauEdges(2,ntauV1p4),modis_histTauCenters(ntauV1p4))
       numMODIStauBins      = ntauV1p4+1
       modis_histTau        = tau_binBoundsV1p4
       modis_histTauEdges   = tau_binEdgesV1p4
       modis_histTauCenters = tau_binCentersV1p4
    endif

    if (.not. allocated(vgrid_zl) .or. .not. allocated(vgrid_zu) .or. .not. allocated(vgrid_z)) then
       
       ! Initialize quickbeam_optics, also if two-moment radar microphysics scheme is wanted...
       if (cloudsat_micro_scheme == 'MMF_v3.5_two_moment')  then
          ldouble = .true. 
          lsingle = .false.
       endif
       call quickbeam_optics_init()

       ! Initialize the distributional parameters for hydrometeors in radar simulator
       call hydro_class_init(lsingle,ldouble,sd)
       
       ! Initialize COSP simulator
       call COSP_INIT(cfg%Lisccp_sim,cfg%Lmodis_sim,cfg%Lmisr_sim,cfg%Lradar_sim,        &
            cfg%Llidar_sim,cfg%Lparasol_sim,cfg%Lrttov_sim,gbx%Npoints,gbx%Nlevels,      &
            gbx%radar_freq,gbx%k2,gbx%use_gas_abs,gbx%do_ray,gbx%isccp_top_height,       &
            gbx%isccp_top_height_direction,gbx%surface_radar,rcfg_cloudsat,gbx%Nchan,    &
            gbx%Ichan,gbx%plat,gbx%sat,gbx%inst,vgrid%use_vgrid,vgrid%csat_vgrid,        &
            vgrid%Nlvgrid,cloudsat_micro_scheme)
    endif
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Construct output type for cosp
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
                                cfg%Ltbrttov,gbx%Npoints,gbx%Ncolumns,gbx%Nlevels,       &
                                vgrid%Nlvgrid,gbx%Nchan,cospOUT)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Break COSP into chunks, only applicable when gbx%Npoints_it > gbx%Npoints
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    num_chunks = gbx%Npoints/gbx%Npoints_it+1
    do i=1,num_chunks
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Determine indices for "chunking" (again, if necessary)
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       if (num_chunks .eq. 1) then
          start_idx = 1
          end_idx   = gbx%Npoints
          Nptsperit = gbx%Npoints
       else
          start_idx = (i-1)*gbx%Npoints_it+1
          end_idx   = i*gbx%Npoints_it
          if (end_idx .gt. gbx%Npoints) end_idx=gbx%Npoints
          Nptsperit = end_idx-start_idx+1
       endif

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
       ! Generate subcolumns and compute optical inputs to COSP.
       ! This subroutine essentially contains all of the pieces of code that were removed
       ! from the simulators during the v2.0 reconstruction.
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       call subsample_and_optics(overlap,gbx,sgx,cfg,Nptsperit,start_idx,end_idx,cospIN,     &
                                 cospstateIN)
       
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Call COSPv2.0
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       cosp_status = COSP_SIMULATOR(cospIN, cospstateIN, cospOUT, start_idx,end_idx,.false.) 
    enddo
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Free up memory
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    call destroy_cospIN(cospIN)
    call destroy_cospstateIN(cospstateIN)
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Copy new output to old output types.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! MISR
    if (cfg%Lmisr_sim) then
       if (cfg%LclMISR) misr%fq_MISR  = cospOUT%misr_fq
       ! *NOTE* These 3 fields are not output, but were part of the v1.4.0 cosp_misr, so
       !        they are still computed. Should probably have a logical to control these
       !        outputs in cosp_config. In the meantime, only assign v1.4.0 outputs to
       !        v2.0 outputs IF a MISR diagnostic was requested.
       if (cfg%LclMISR) misr%MISR_meanztop             = cospOUT%misr_meanztop
       if (cfg%LclMISR) misr%MISR_cldarea              = cospOUT%misr_cldarea
       if (cfg%LclMISR) misr%MISR_dist_model_layertops = cospOUT%misr_dist_model_layertops
    endif
    
    ! ISCCP
    if (cfg%Lisccp_sim) then
       if (cfg%Lboxtauisccp)    isccp%boxtau        = cospOUT%isccp_boxtau
       if (cfg%Lboxptopisccp)   isccp%boxptop       = cospOUT%isccp_boxptop
       if (cfg%Lclisccp)        isccp%fq_isccp      = cospOUT%isccp_fq
       if (cfg%Lcltisccp)       isccp%totalcldarea  = cospOUT%isccp_totalcldarea
       if (cfg%Lmeantbisccp)    isccp%meantb        = cospOUT%isccp_meantb
       if (cfg%Lmeantbclrisccp) isccp%meantbclr     = cospOUT%isccp_meantbclr
       if (cfg%Lpctisccp)       isccp%meanptop      = cospOUT%isccp_meanptop
       if (cfg%Ltauisccp)       isccp%meantaucld    = cospOUT%isccp_meantaucld
       if (cfg%Lalbisccp)       isccp%meanalbedocld = cospOUT%isccp_meanalbedocld
   endif

    ! MODIS
    if (cfg%Lmodis_sim) then
       if (cfg%Lcltmodis)     modis%Cloud_Fraction_Total_Mean =                         &
                          cospOUT%modis_Cloud_Fraction_Total_Mean
       if (cfg%Lclwmodis)     modis%Cloud_Fraction_Water_Mean =                         &
                          cospOUT%modis_Cloud_Fraction_Water_Mean
       if (cfg%Lclimodis)     modis%Cloud_Fraction_Ice_Mean =                           &
                          cospOUT%modis_Cloud_Fraction_Ice_Mean
       if (cfg%Lclhmodis)     modis%Cloud_Fraction_High_Mean =                          &
                          cospOUT%modis_Cloud_Fraction_High_Mean
       if (cfg%Lclmmodis)     modis%Cloud_Fraction_Mid_Mean =                           &
                          cospOUT%modis_Cloud_Fraction_Mid_Mean
       if (cfg%Lcllmodis)     modis%Cloud_Fraction_Low_Mean =                           &
                          cospOUT%modis_Cloud_Fraction_Low_Mean
       if (cfg%Ltautmodis)    modis%Optical_Thickness_Total_Mean =                      &
                          cospOUT%modis_Optical_Thickness_Total_Mean
       if (cfg%Ltauwmodis)    modis%Optical_Thickness_Water_Mean =                      &
                          cospOUT%modis_Optical_Thickness_Water_Mean
       if (cfg%Ltauimodis)    modis%Optical_Thickness_Ice_Mean =                        &
                          cospOUT%modis_Optical_Thickness_Ice_Mean
       if (cfg%Ltautlogmodis) modis%Optical_Thickness_Total_LogMean =                   &
                          cospOUT%modis_Optical_Thickness_Total_LogMean
       if (cfg%Ltauwlogmodis) modis%Optical_Thickness_Water_LogMean =                   &
                          cospOUT%modis_Optical_Thickness_Water_LogMean
       if (cfg%Ltauilogmodis) modis%Optical_Thickness_Ice_LogMean =                     &
                          cospOUT%modis_Optical_Thickness_Ice_LogMean
       if (cfg%Lreffclwmodis) modis%Cloud_Particle_Size_Water_Mean =                    &
                          cospOUT%modis_Cloud_Particle_Size_Water_Mean
       if (cfg%Lreffclimodis) modis%Cloud_Particle_Size_Ice_Mean =                      &
                          cospOUT%modis_Cloud_Particle_Size_Ice_Mean
       if (cfg%Lpctmodis)     modis%Cloud_Top_Pressure_Total_Mean =                     &
                          cospOUT%modis_Cloud_Top_Pressure_Total_Mean
       if (cfg%Llwpmodis)     modis%Liquid_Water_Path_Mean =                            &
                          cospOUT%modis_Liquid_Water_Path_Mean
       if (cfg%Liwpmodis)     modis%Ice_Water_Path_Mean =                               &
                          cospOUT%modis_Ice_Water_Path_Mean
       if (cfg%Lclmodis) then
          modis%Optical_Thickness_vs_Cloud_Top_Pressure =                               &
             cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure
          modis%Optical_Thickness_vs_ReffICE = cospOUT%modis_Optical_Thickness_vs_ReffICE
          modis%Optical_Thickness_vs_ReffLIQ = cospOUT%modis_Optical_Thickness_vs_ReffLIQ
       endif
    endif

    ! PARASOL
    if (cfg%Lparasol_sim) then
       if (cfg%Lparasolrefl) sglidar%refl        = cospOUT%parasolPix_refl
       if (cfg%Lparasolrefl) stlidar%parasolrefl = cospOUT%parasolGrid_refl
    endif

    ! RTTOV
    if (cfg%Lrttov_sim) then
       if (cfg%Ltbrttov) rttov%tbs = cospOUT%rttov_tbs  
    endif
    
    ! CALIPSO
    if (cfg%Llidar_sim) then
       ! *NOTE* In COSPv2.0 all outputs are ordered from TOA-2-SFC, but in COSPv1.4 this is
       !        not true. To maintain the outputs of v1.4, the affected fields are flipped.

       if (cfg%LlidarBetaMol532) then
          sglidar%beta_mol         = cospOUT%calipso_beta_mol
       endif
       if (cfg%Latb532) then
          sglidar%beta_tot         = cospOUT%calipso_beta_tot
       endif
       if (cfg%LcfadLidarsr532)  then
          stlidar%srbval       = cospOUT%calipso_srbval
          stlidar%cfad_sr      = cospOUT%calipso_cfad_sr(:,:,vgrid%Nlvgrid:1:-1)
          sglidar%betaperp_tot = cospOUT%calipso_betaperp_tot
       endif

       if (cfg%Lclcalipso) then
          stlidar%lidarcld = cospOUT%calipso_lidarcld(:,stlidar%Nlevels:1:-1)
       endif       
       if (cfg%Lclhcalipso .or. cfg%Lclmcalipso .or. cfg%Lcllcalipso .or. cfg%Lcltcalipso) then
          stlidar%cldlayer = cospOUT%calipso_cldlayer
       endif
       if (cfg%Lclcalipsoice .or. cfg%Lclcalipsoliq .or. cfg%Lclcalipsoun) then
          stlidar%lidarcldphase = cospOUT%calipso_lidarcldphase(:,vgrid%Nlvgrid:1:-1,:)
       endif
       if (cfg%Lcllcalipsoice .or. cfg%Lclmcalipsoice .or. cfg%Lclhcalipsoice .or.                   &
           cfg%Lcltcalipsoice .or. cfg%Lcllcalipsoliq .or. cfg%Lclmcalipsoliq .or.                   &
           cfg%Lclhcalipsoliq .or. cfg%Lcltcalipsoliq .or. cfg%Lcllcalipsoun  .or.                   &
           cfg%Lclmcalipsoun  .or. cfg%Lclhcalipsoun  .or. cfg%Lcltcalipsoun) then       
          stlidar%cldlayerphase         = cospOUT%calipso_cldlayerphase
       endif
       if (cfg%Lclcalipsotmp .or. cfg%Lclcalipsotmpliq .or. cfg%Lclcalipsoice .or. cfg%Lclcalipsotmpun) then
          stlidar%lidarcldtmp = cospOUT%calipso_lidarcldtmp
       endif
       ! Fields present, but not controlled by logical switch
       if (any([cfg%Lclcalipsoliq,cfg%Lclcalipsoice,cfg%Lclcalipsoun,cfg%Lclcalipsotmp,          &
            cfg%Lclcalipsotmpliq,cfg%Lclcalipsotmpice,cfg%Lclcalipsotmpun,cfg%Lclhcalipsoliq,&
            cfg%Lcllcalipsoliq,cfg%Lclmcalipsoliq,cfg%Lcltcalipsoliq,cfg%Lclhcalipsoice,&
            cfg%Lcllcalipsoice,cfg%Lclmcalipsoice,cfg%Lcltcalipsoice,cfg%Lclhcalipsoun,&
            cfg%Lcllcalipsoun,cfg%Lclmcalipsoun,cfg%Lcltcalipsoun])) then
          sglidar%temp_tot = cospOUT%calipso_temp_tot(:,sglidar%Nlevels:1:-1)
          sglidar%tau_tot  = cospOUT%calipso_tau_tot(:,:,sglidar%Nlevels:1:-1)
       endif       
    endif
    
    ! Cloudsat             
    if (cfg%Lradar_sim) then
       ! *NOTE* In COSP2 all outputs are ordered from TOA-2-SFC, but in COSPv1.4 this is
       !        not true. To maintain the outputs of v1.4, the affected fields are flipped.    
       if (cfg%Ldbze94) then
          sgradar%Ze_tot = cospOUT%cloudsat_Ze_tot
       endif
       if (cfg%LcfadDbze94) then 
          stradar%cfad_ze = cospOUT%cloudsat_cfad_ze(:,:,stradar%Nlevels:1:-1)              
       endif
    endif

    ! Combined instrument products
    if (cfg%Lclcalipso2) then
       stradar%lidar_only_freq_cloud = cospOUT%lidar_only_freq_cloud(:,stradar%Nlevels:1:-1)    
    endif
    if (cfg%Lcltlidarradar) stradar%radar_lidar_tcc = cospOUT%radar_lidar_tcc      
    
    ! Subcolumns
    sgx%frac_out = sgx%frac_out(:,:,sgx%Nlevels:1:-1)
    
  end subroutine cosp_interface_v1p4
   
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! SUBROUTINE subsample_and_optics
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine subsample_and_optics(overlap,gbx,sgx,cfg,npoints,start_idx,end_idx,cospIN,cospgridIN)

    ! Inputs
    integer, intent(in) :: overlap  ! Overlap type in SCOPS: 1=max, 2=rand, 3=max/rand
    type(cosp_gridbox),intent(in)    :: gbx   ! Grid box description
    type(cosp_config),intent(in)     :: cfg   ! Configuration information
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
                                               cv_p_rate,g_vol
    real(wp),dimension(:,:,:),allocatable :: hm_matrix,re_matrix,                        &
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
    cospIN%tautot_S_liq                                 = 0._wp
    cospIN%tautot_S_ice                                 = 0._wp
    cospIN%emsfc_lw                                     = gbx%isccp_emsfc_lw
    cospIN%rcfg_cloudsat                                = rcfg_cloudsat
    cospgridIN%hgt_matrix(1:nPoints,1:gbx%Nlevels)      = gbx%zlev(start_idx:end_idx,gbx%Nlevels:1:-1)
    cospgridIN%hgt_matrix_half(1:nPoints,1:gbx%Nlevels) = gbx%zlev_half(start_idx:end_idx,gbx%Nlevels:1:-1)
    cospgridIN%sunlit(1:nPoints)                        = gbx%sunlit(start_idx:end_idx)
    cospgridIN%skt(1:nPoints)                           = gbx%skt(start_idx:end_idx)
    cospgridIN%land(1:nPoints)                          = gbx%land(start_idx:end_idx)
    cospgridIN%qv(1:nPoints,1:gbx%Nlevels)              = gbx%sh(start_idx:end_idx,gbx%Nlevels:1:-1) 
    cospgridIN%at(1:nPoints,1:gbx%Nlevels)              = gbx%T(start_idx:end_idx,gbx%Nlevels:1:-1) 
    cospgridIN%pfull(1:nPoints,1:gbx%Nlevels)           = gbx%p(start_idx:end_idx,gbx%Nlevels:1:-1) 
    cospgridIN%o3(1:nPoints,1:gbx%Nlevels)              = gbx%mr_ozone(start_idx:end_idx,gbx%Nlevels:1:-1)*(amd/amO3)*1e6
    cospgridIN%u_sfc(1:nPoints)                         = gbx%u_wind(start_idx:end_idx)
    cospgridIN%v_sfc(1:nPoints)                         = gbx%v_wind(start_idx:end_idx)
    cospgridIN%emis_sfc                                 = gbx%surfem
    cospgridIN%lat(1:nPoints)                           = gbx%latitude(start_idx:end_idx)
    cospgridIN%lon(1:nPoints)                           = gbx%longitude(start_idx:end_idx)
    cospgridIN%month                                    = 2 ! This is needed by RTTOV only for the surface emissivity calculation.
    cospgridIN%co2                                      = gbx%co2*(amd/amCO2)*1e6
    cospgridIN%ch4                                      = gbx%ch4*(amd/amCH4)*1e6  
    cospgridIN%n2o                                      = gbx%n2o*(amd/amN2O)*1e6
    cospgridIN%co                                       = gbx%co*(amd/amCO)*1e6
    cospgridIN%zenang                                   = gbx%zenang
    cospgridIN%phalf(:,1)                               = 0._wp
    cospgridIN%phalf(:,2:gbx%Nlevels+1)                 = gbx%ph(start_idx:end_idx,gbx%Nlevels:1:-1)    
    if (gbx%Ncolumns .gt. 1) then
       
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Random number generator
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       allocate(rngs(Npoints),seed(Npoints))
       seed(:)=0
       seed = int(gbx%psfc)  ! In case of Npoints=1
       if (Npoints .gt. 1) seed=int((gbx%psfc(start_idx:end_idx)-minval(gbx%psfc(start_idx:end_idx)))/      &
            (maxval(gbx%psfc(start_idx:end_idx))-minval(gbx%psfc(start_idx:end_idx)))*100000) + 1
       call init_rng(rngs, seed)  

       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Generate subcolumns for clouds (SCOPS) and precipitation type (PREC_SCOPS)
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Call SCOPS
       if (gbx%Ncolumns .gt. 1) then
          call scops(npoints,gbx%Nlevels,gbx%Ncolumns,rngs,                              &
                     gbx%tca(start_idx:end_idx,gbx%Nlevels:1:-1),                        &
                     gbx%cca(start_idx:end_idx,gbx%Nlevels:1:-1),overlap,                &
                     sgx%frac_out(start_idx:end_idx,:,:),0)
          deallocate(seed,rngs)
       else
          sgx%frac_out(start_idx:end_idx,:,:) = 1  
       endif
       cospIN%frac_out=sgx%frac_out(start_idx:end_idx,:,:)
       
       ! Sum up precipitation rates
       allocate(ls_p_rate(npoints,gbx%Nlevels),cv_p_rate(npoints,gbx%Nlevels))
       if(gbx%use_precipitation_fluxes) then
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
             if (gbx%use_precipitation_fluxes) then
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
       if (gbx%use_precipitation_fluxes) then
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
    if (cfg%Llidar_sim) then
       call lidar_optics(npoints,gbx%Ncolumns,gbx%Nlevels,4,gbx%lidar_ice_type,          &
                         mr_hydro(:,:,cospIN%Nlevels:1:-1,I_LSCLIQ),                     &
                         mr_hydro(:,:,cospIN%Nlevels:1:-1,I_LSCICE),                     &
                         mr_hydro(:,:,cospIN%Nlevels:1:-1,I_CVCLIQ),                     &
                         mr_hydro(:,:,cospIN%Nlevels:1:-1,I_CVCICE),                     &
                         gbx%Reff(start_idx:end_idx,cospIN%Nlevels:1:-1,I_LSCLIQ),       &
                         gbx%Reff(start_idx:end_idx,cospIN%Nlevels:1:-1,I_LSCICE),       &
                         gbx%Reff(start_idx:end_idx,cospIN%Nlevels:1:-1,I_CVCLIQ),       &
                         gbx%Reff(start_idx:end_idx,cospIN%Nlevels:1:-1,I_CVCICE),       & 
                         cospgridIN%pfull,cospgridIN%phalf,cospgridIN%at,                &
                         cospIN%beta_mol,cospIN%betatot,cospIN%tau_mol,cospIN%tautot,    &
                         cospIN%tautot_S_liq,cospIN%tautot_S_ice,                        &
                         betatot_ice = cospIN%betatot_ice,                               &
                         betatot_liq=cospIN%betatot_liq,tautot_ice=cospIN%tautot_ice,    &
                         tautot_liq = cospIN%tautot_liq)
    endif
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! CLOUDSAT RADAR OPTICS
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (cfg%Lradar_sim) then
       allocate(g_vol(nPoints,gbx%Nlevels))
       ! Compute gaseous absorption (assume identical for each subcolun)
       g_vol(:,:)=0._wp
       do i=1,nPoints
          do j=1,cospIN%nLevels
             if (rcfg_cloudsat%use_gas_abs == 1 .or. (rcfg_cloudsat%use_gas_abs == 2 .and. j .eq. cospIN%Nlevels)) then
                g_vol(i,j) = gases(cospgridIN%pfull(i,j), cospgridIN%at(i,j),cospgridIN%qv(i,j),rcfg_cloudsat%freq)
             endif
             cospIN%g_vol_cloudsat(i,:,j)=g_vol(i,j)
          end do
       end do
       
       ! Loop over all subcolumns
       do ij=1,gbx%nColumns
          call quickbeam_optics(sd, rcfg_cloudsat, nPoints, cospIN%nLevels, R_UNDEF,      &
               mr_hydro(:,ij,cospIN%Nlevels:1:-1,1:N_HYDRO)*1000._wp,                     &
               Reff(:,ij,cospIN%Nlevels:1:-1,1:N_HYDRO)*1.e6_wp,                          &
               Np(:,ij,cospIN%Nlevels:1:-1,1:N_HYDRO),cospgridIN%pfull, cospgridIN%at,    &
               cospgridIN%qv, cospIN%z_vol_cloudsat(1:nPoints,ij,:),                      &
               cospIN%kr_vol_cloudsat(1:nPoints,ij,:))
       enddo
    end if

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! MODIS optics
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (cfg%Lmodis_sim) then
       ! Allocate memory
       allocate(MODIS_cloudWater(npoints,gbx%Ncolumns,gbx%Nlevels),                      &
                MODIS_cloudIce(npoints,gbx%Ncolumns,gbx%Nlevels),                        &
                MODIS_waterSize(npoints,gbx%Ncolumns,gbx%Nlevels),                       &
                MODIS_iceSize(npoints,gbx%Ncolumns,gbx%Nlevels),                         &
                MODIS_opticalThicknessLiq(npoints,gbx%Ncolumns,gbx%Nlevels),             &
                MODIS_opticalThicknessIce(npoints,gbx%Ncolumns,gbx%Nlevels))
       ! Cloud water
       call cosp_simulator_optics(npoints,gbx%Ncolumns,gbx%Nlevels,                      &
            sgx%frac_out(start_idx:end_idx,:,:),mr_hydro(:,:,:,I_CVCLIQ),                &
            mr_hydro(:,:,:,I_LSCLIQ),MODIS_cloudWater(:, :, gbx%Nlevels:1:-1))   
       ! Cloud ice
       call cosp_simulator_optics(npoints,gbx%Ncolumns,gbx%Nlevels,                      &
            sgx%frac_out(start_idx:end_idx,:,:),mr_hydro(:,:,:,I_CVCICE),                &
            mr_hydro(:,:,:,I_LSCICE),MODIS_cloudIce(:, :, gbx%Nlevels:1:-1))  
       ! Water droplet size
       call cosp_simulator_optics(npoints,gbx%Ncolumns,gbx%Nlevels,                      &
            sgx%frac_out(start_idx:end_idx,:,:),reff(:,:,:,I_CVCLIQ),                    &
            reff(:,:,:,I_LSCLIQ),MODIS_waterSize(:, :, gbx%Nlevels:1:-1))
       ! Ice crystal size
       call cosp_simulator_optics(npoints,gbx%Ncolumns,gbx%Nlevels,                      &
            sgx%frac_out(start_idx:end_idx,:,:),reff(:,:,:,I_CVCICE),                    &
            reff(:,:,:,I_LSCICE),MODIS_iceSize(:, :, gbx%Nlevels:1:-1))
       ! Partition optical thickness into liquid and ice parts
       call modis_optics_partition(npoints,gbx%Nlevels,gbx%Ncolumns,MODIS_cloudWater,    &
            MODIS_cloudIce,MODIS_waterSize,MODIS_iceSize,cospIN%tau_067,                 &
            MODIS_opticalThicknessLiq, MODIS_opticalThicknessIce)
       ! Compute assymetry parameter and single scattering albedo 
       call modis_optics(npoints,gbx%Nlevels,gbx%Ncolumns,MODIS_opticalThicknessLiq,     &
            MODIS_waterSize*1.0e6_wp,MODIS_opticalThicknessIce,MODIS_iceSize*1.0e6_wp,   &
            cospIN%fracLiq, cospIN%asym, cospIN%ss_alb)
       
       ! Deallocate memory
       deallocate(MODIS_cloudWater,MODIS_cloudIce,MODIS_WaterSize,MODIS_iceSize,         &
                  MODIS_opticalThicknessLiq,MODIS_opticalThicknessIce,mr_hydro,          &
                  Reff,Np)
    end if
    
  end subroutine subsample_and_optics
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE construct_cosp_gridbox
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_cosp_gridbox(time,time_bnds,radar_freq,surface_radar,         &
                                         use_mie_tables,use_gas_abs,do_ray,melt_lay,k2,   &
                                         Npoints,Nlevels,Ncolumns,Nhydro,Nprmts_max_hydro,&
                                         Naero,Nprmts_max_aero,Npoints_it,lidar_ice_type, &
                                         isccp_top_height,isccp_top_height_direction,     &
                                         isccp_overlap,isccp_emsfc_lw,                    &
                                         use_precipitation_fluxes,use_reff,Plat,Sat,Inst, &
                                         Nchan,ZenAng,Ichan,SurfEm,co2,ch4,n2o,co,        &
                                         y,load_LUT)
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

    ! Outputs
    type(cosp_gridbox),intent(out) :: y
    
    ! local variables
    logical :: local_load_LUT
    
    if (present(load_LUT)) then
       local_load_LUT = load_LUT
    else
       local_load_LUT = RADAR_SIM_LOAD_scale_LUTs_flag
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
    y%ichan  = ichan
    y%surfem = surfem
    
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
    
  END SUBROUTINE CONSTRUCT_cosp_gridbox
    
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE destroy_cosp_gridbox
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE destroy_cosp_gridbox(y,save_LUT)
    
    type(cosp_gridbox),intent(inout) :: y
    logical,intent(in),optional :: save_LUT
    
    logical :: local_save_LUT
    if (present(save_LUT)) then
       local_save_LUT = save_LUT
    else
       local_save_LUT = RADAR_SIM_UPDATE_scale_LUTs_flag
    endif
    
    ! save any updates to radar simulator LUT
    if (local_save_LUT) call save_scale_LUTs(y%hp)
    
    deallocate(y%zlev,y%zlev_half,y%dlev,y%p,y%ph,y%T,y%q,y%sh,y%dtau_s,y%dtau_c,y%dem_s,&
               y%dem_c,y%toffset,y%longitude,y%latitude,y%psfc,y%land,y%tca,y%cca,       &
               y%mr_hydro,y%dist_prmts_hydro,y%conc_aero,y%dist_type_aero,               &
               y%dist_prmts_aero,y%rain_ls,y%rain_cv,y%snow_ls,y%snow_cv,y%grpl_ls,      &
               y%sunlit,y%skt,y%Reff,y%Np,y%ichan,y%surfem,y%mr_ozone,y%u_wind,y%v_wind)
    
  END SUBROUTINE destroy_cosp_gridbox
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE construct_cosp_subgrid
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
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE save_scale_LUTs
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine save_scale_LUTs(hp)
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
  end subroutine save_scale_LUTs

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !SUBROUTINE construct_cosp_vgrid
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_COSP_VGRID(gbx,Nlvgrid,use_vgrid,cloudsat,x)
    type(cosp_gridbox),intent(in) :: gbx ! Gridbox information
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
    
  END SUBROUTINE CONSTRUCT_COSP_VGRID
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
    allocate(x%cfad_ze(Npoints,cloudsat_DBZE_BINS,Nlevels),x%lidar_only_freq_cloud(Npoints,Nlevels), &
             x%radar_lidar_tcc(Npoints))
    
    ! Initialize
    x%cfad_ze               = 0._wp
    x%lidar_only_freq_cloud = 0._wp
    x%radar_lidar_tcc       = 0._wp    
    
  end subroutine construct_cosp_radarstats
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE destroy_cosp_subgrid
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine destroy_cosp_subgrid(y)
    type(cosp_subgrid),intent(inout) :: y   
    deallocate(y%prec_frac, y%frac_out)
  end subroutine destroy_cosp_subgrid
  
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

    if (associated(x%fq_MISR))                   deallocate(x%fq_MISR)
    if (associated(x%MISR_meanztop))             deallocate(x%MISR_meanztop)
    if (associated(x%MISR_cldarea))              deallocate(x%MISR_cldarea)
    if (associated(x%MISR_dist_model_layertops)) deallocate(x%MISR_dist_model_layertops)

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
             x%Optical_Thickness_vs_Cloud_Top_Pressure(nPoints,ntauV1p4+1,numMODISPresBins),&
             x%Optical_Thickness_vs_ReffICE(nPoints,ntauV1p4+1,numMODISReffIceBins),&
             x%Optical_Thickness_vs_ReffLIQ(nPoints,ntauV1p4+1,numMODISReffLiqBins))
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
    
    allocate(y%tau_067(npoints,        ncolumns,nlevels),&
             y%emiss_11(npoints,       ncolumns,nlevels),&
             y%frac_out(npoints,       ncolumns,nlevels),&
             y%betatot(npoints,        ncolumns,nlevels),&
             y%betatot_ice(npoints,    ncolumns,nlevels),&
             y%fracLiq(npoints,        ncolumns,nlevels),&
             y%betatot_liq(npoints,    ncolumns,nlevels),&
             y%tautot(npoints,         ncolumns,nlevels),&
             y%tautot_ice(npoints,     ncolumns,nlevels),&
             y%tautot_liq(npoints,     ncolumns,nlevels),&
             y%z_vol_cloudsat(npoints, ncolumns,nlevels),&
             y%kr_vol_cloudsat(npoints,ncolumns,nlevels),&
             y%g_vol_cloudsat(npoints, ncolumns,nlevels),&
             y%asym(npoints,           ncolumns,nlevels),&
             y%ss_alb(npoints,         ncolumns,nlevels),&
             y%beta_mol(npoints,                nlevels),&
             y%tau_mol(npoints,                 nlevels),&
             y%tautot_S_ice(npoints,   ncolumns        ),&
             y%tautot_S_liq(npoints,   ncolumns))
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
             y%cloudIce(nPoints,nLevels),y%cloudLiq(nPoints,nLevels),                    &
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
                                    Liwpmodis,Lclmodis,Latb532,LlidarBetaMol532,         &
                                    LcfadLidarsr532,Lclcalipso2,                         &
                                    Lclcalipso,Lclhcalipso,Lcllcalipso,Lclmcalipso,      &
                                    Lcltcalipso,Lcltlidarradar,Lclcalipsoliq,            &
                                    Lclcalipsoice,Lclcalipsoun,Lclcalipsotmp,            &
                                    Lclcalipsotmpliq,Lclcalipsotmpice,Lclcalipsotmpun,   &
                                    Lcltcalipsoliq,Lcltcalipsoice,Lcltcalipsoun,         &
                                    Lclhcalipsoliq,Lclhcalipsoice,Lclhcalipsoun,         &
                                    Lclmcalipsoliq,Lclmcalipsoice,Lclmcalipsoun,         &
                                    Lcllcalipsoliq,Lcllcalipsoice,Lcllcalipsoun,         & 
                                    LcfadDbze94,Ldbze94,Lparasolrefl,Ltbrttov, &
                                    Npoints,Ncolumns,Nlevels,Nlvgrid,Nchan,x)
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
      
    ! PARASOL
    if (Lparasolrefl) then
        allocate(x%parasolPix_refl(Npoints,Ncolumns,PARASOL_NREFL))
        allocate(x%parasolGrid_refl(Npoints,PARASOL_NREFL))
    endif 

    ! Cloudsat simulator
    if (Ldbze94)        allocate(x%cloudsat_Ze_tot(Npoints,Ncolumns,Nlevels))
    if (LcfadDbze94)    allocate(x%cloudsat_cfad_ze(Npoints,cloudsat_DBZE_BINS,Nlvgrid))

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

    
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !                                    END MODULE
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module MOD_COSP_INTERFACE_v1p4
