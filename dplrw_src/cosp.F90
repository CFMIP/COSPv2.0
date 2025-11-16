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
! History:
! May 2015- D. Swales    - Original version
! Oct 2018- T. Michibata - Inline Diagnostic Driver (IDiD) for MIROC6 interface
! 
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE MOD_COSP
  USE COSP_KINDS,                  ONLY: wp
  USE MOD_COSP_CONFIG,             ONLY: R_UNDEF,PARASOL_NREFL,LIDAR_NCAT,SR_BINS,       &
                                         N_HYDRO,RTTOV_MAX_CHANNELS,numMISRHgtBins,      &
                                         DBZE_BINS,LIDAR_NTEMP,calipso_histBsct,         &
                                         use_vgrid,Nlvgrid,vgrid_zu,vgrid_zl,vgrid_z,dz, &
                                         WR_NREGIME,    CFODD_NCLASS,                    &
                                         CFODD_NDBZE,   CFODD_NICOD,                     &
                                         numMODISTauBins,numMODISPresBins,               &
                                         numMODISReffIceBins,numMODISReffLiqBins,        &
                                         numISCCPTauBins,numISCCPPresBins,numMISRTauBins,&
                                         ntau,npres,modis_histTau,tau_binBounds,         & !<- modified by Doppler
                                         modis_histTauEdges,tau_binEdges,nCloudsatPrecipClass,& !<-added
                                         modis_histTauCenters,tau_binCenters,            &
                                         cloudsat_preclvl!,grLidar532_histBsct,atlid_histBsct !<-added 
  USE MOD_COSP_MODIS_INTERFACE,    ONLY: cosp_modis_init,     modis_IN
  USE MOD_COSP_RTTOV_INTERFACE,    ONLY: cosp_rttov_init,     rttov_IN
  USE MOD_COSP_MISR_INTERFACE,     ONLY: cosp_misr_init,      misr_IN
  USE MOD_COSP_ISCCP_INTERFACE,    ONLY: cosp_isccp_init,     isccp_IN
  USE MOD_COSP_CALIPSO_INTERFACE,  ONLY: cosp_calipso_init,   calipso_IN
  USE MOD_COSP_PARASOL_INTERFACE,  ONLY: cosp_parasol_init,   parasol_in
  USE MOD_COSP_CLOUDSAT_INTERFACE, ONLY: cosp_cloudsat_init,  cloudsat_IN
  USE quickbeam,                   ONLY: quickbeam_subcolumn, quickbeam_column, radar_cfg
  USE MOD_ICARUS,                  ONLY: icarus_subcolumn,    icarus_column
  USE MOD_MISR_SIMULATOR,          ONLY: misr_subcolumn,      misr_column
  USE MOD_LIDAR_SIMULATOR,         ONLY: lidar_subcolumn,     lidar_column
  USE MOD_MODIS_SIM,               ONLY: modis_subcolumn,     modis_column
  USE MOD_PARASOL,                 ONLY: parasol_subcolumn,   parasol_column
  use mod_cosp_rttov,              ONLY: rttov_column
  USE MOD_COSP_STATS,              ONLY: COSP_LIDAR_ONLY_CLOUD,COSP_CHANGE_VERTICAL_GRID, &
                                         COSP_DIAG_WARMRAIN
  !USE MOD_COSP_DIAGNOSTICS,         ONLY: COSP_DIAG_WARMRAIN, COSP_DIAG_PHASE, &
  !                                        COSP_DIAG_PARTICLE_TYPES

#ifdef OPT_DPLRW
  use quickbeam,                   ONLY: quickbeam_dplrw
  use MOD_COSP_CONFIG,             ONLY: Nlvtemp,NlvdBZe,Nlvdplr,Nlvspwd
#endif

  IMPLICIT NONE
  
  logical :: linitialization ! Initialization flag
  save linitialization           !Jing X.
  data linitialization /.true./  !Jing X.
  
  ! ######################################################################################
  ! TYPE cosp_column_inputs
  ! ######################################################################################
  type cosp_column_inputs
     integer :: &
          Npoints,             & ! Number of gridpoints.
          Ncolumns,            & ! Number of columns.
          Nlevels                ! Number of levels.
         
     integer,allocatable,dimension(:) :: &
          sunlit                 ! Sunlit flag                            (0-1)

     real(wp),allocatable,dimension(:,:) :: &
          at,                  & ! Temperature                            (K)
          pfull,               & ! Pressure                               (Pa)
          phalf,               & ! Pressure at half-levels                (Pa)
          qv,                  & ! Specific humidity                      (kg/kg)
          hgt_matrix,          & ! Height of hydrometeors                 (km)
          hgt_matrix_half        ! Height of hydrometeors at half levels  (km)
#ifdef OPT_DPLRW
     real(wp),allocatable,dimension(:,:) :: &
          gwvel,               & ! w velocity dynamical diag (full level) [m/s]
          gcumf                  ! convective mass flux from cumulus scheme (half level) [kg m^-2 s^-1]
!     logical :: Ldplrw_Taxis
#endif

     real(wp),allocatable,dimension(:) :: &
          land,                & ! Land/Sea mask                          (0-1)
          skt,                 & ! Surface temperature                    (K)
          surfelev               ! Surface Elevation                      (m)  <-added
     ! Fields used ONLY by RTTOV
     integer :: &
          month                  ! Month for surface emissivty atlas      (1-12)
     real(wp) :: &
          zenang,              & ! Satellite zenith angle for RTTOV       (deg)
          co2,                 & ! CO2                                    (kg/kg)
          ch4,                 & ! Methane                                (kg/kg)
          n2o,                 & ! N2O                                    (kg/kg)
          co                     ! CO                                     (kg/kg)
     real(wp),allocatable,dimension(:) :: &
          emis_sfc,            & ! Surface emissivity                     (1)
          u_sfc,               & ! Surface u-wind                         (m/s)
          v_sfc,               & ! Surface v-wind                         (m/s)
          seaice,              & ! Sea-ice fraction                       (0-1)
          lat,                 & ! Latitude                              (deg)
          lon                    ! Longitude                              (deg)
     real(wp),allocatable,dimension(:,:) :: &
          o3,                  & ! Ozone                                  (kg/kg)
          tca,                 & ! Total column cloud fraction            (0-1)
          cloudIce,            & ! Cloud ice water mixing ratio           (kg/kg)
          cloudLiq,            & ! Cloud liquid water mixing ratio        (kg/kg)
          fl_rain,             & ! Precipitation (rain) flux              (kg/m2/s)
          fl_snow                ! Precipitation (snow) flux              (kg/m2/s)
  end type cosp_column_inputs
  
  ! ######################################################################################
  ! TYPE cosp_optical_inputs
  ! ######################################################################################  
  type cosp_optical_inputs
     integer :: &
          Npoints,             & ! Number of gridpoints.
          Ncolumns,            & ! Number of columns.
          Nlevels,             & ! Number of levels.
          Npart,               & ! Number of cloud meteors for LIDAR simulator.
          Nrefl                  ! Number of reflectances for PARASOL simulator
     real(wp) :: &
          emsfc_lw               ! 11 micron surface emissivity
     real(wp),allocatable,dimension(:,:,:) :: &
          frac_out,            & ! Cloud fraction
          frac_prec,           & ! prec SCOPS  ! added by YN
          tau_067,             & ! Optical depth
          fracLiq,             & ! Cloud fraction
          emiss_11,            & ! Emissivity
          asym,                & ! Assymetry parameter
          ss_alb,              & ! Single-scattering albedo
          betatot,             & ! Backscatter coefficient for polarized optics (total)
          betatot_ice,         & ! Backscatter coefficient for polarized optics (ice)
          betatot_liq,         & ! Backscatter coefficient for polarized optics (liquid)
          tautot,              & ! Optical thickess integrated from top (total)
          tautot_ice,          & ! Optical thickess integrated from top (ice)
          tautot_liq,          & ! Optical thickess integrated from top (liquid)
          z_vol_cloudsat,      & ! Effective reflectivity factor (mm^6/m^3)
          kr_vol_cloudsat,     & ! Attenuation coefficient hydro (dB/km) 
          g_vol_cloudsat         ! Attenuation coefficient gases (dB/km)
#ifdef OPT_DPLRW
     real(wp),allocatable,dimension(:,:,:,:) :: &
          vfall,           & ! subcolumn fall velocity
          vfsqu,           & ! subcolumn spectrum width (square mean of fall velocity)
          zehyd              ! radar reflectivity for each hydrometeor
#endif

     real(wp),allocatable,dimension(:,:) :: &
          beta_mol,            & ! Molecular backscatter coefficient
          tau_mol,             & ! Molecular optical depth
          tautot_S_liq,        & ! Liquid water optical thickness, from TOA to SFC
          tautot_S_ice,        & ! Ice water optical thickness, from TOA to SFC
          fracPrecipIce          ! Fraction of precipitation which is frozen (1).  <-added
     type(radar_cfg) :: &
          rcfg_cloudsat         ! Radar comfiguration information (CLOUDSAT)
  end type cosp_optical_inputs
  
  ! ######################################################################################
  ! TYPE cosp_outputs
  ! ######################################################################################
  type cosp_outputs

     ! CALIPSO outputs
     real(wp),dimension(:,:,:),pointer :: &
          calipso_betaperp_tot => null(),  & ! Total backscattered signal
          calipso_beta_tot => null(),      & ! Total backscattered signal
          calipso_tau_tot => null(),       & ! Optical thickness integrated from top to level z
          calipso_lidarcldphase => null(), & ! 3D "lidar" phase cloud fraction 
         ! calipso_lidarcldtype => null(),  & ! 3D "lidar" OPAQ type fraction 
          calipso_cldlayerphase => null(), & ! low, mid, high-level lidar phase cloud cover
          calipso_lidarcldtmp => null(),   & ! 3D "lidar" phase cloud temperature
          calipso_cfad_sr => null()          ! CFAD of scattering ratio
     real(wp), dimension(:,:),pointer :: &
          calipso_lidarcld => null(),      & ! 3D "lidar" cloud fraction 
          calipso_cldlayer => null(),      & ! low, mid, high-level, total lidar cloud cover
          calipso_beta_mol => null(),      & ! Molecular backscatter
          calipso_temp_tot => null()
     real(wp), dimension(:),pointer :: &
          calipso_srbval => null()           ! SR bins in cfad_sr
     
     ! PARASOL outputs
     real(wp),dimension(:,:,:),pointer :: &
          parasolPix_refl => null()            ! PARASOL reflectances (subcolumn)    
     real(wp),dimension(:,:),pointer :: &
          parasolGrid_refl => null()           ! PARASOOL reflectances (column)

     ! CLOUDSAT outputs
     real(wp),dimension(:,:,:),pointer :: &
          cloudsat_Ze_tot => null(),         & ! Effective reflectivity factor (Npoints,Ncolumns,Nlevels)     
          cloudsat_cfad_ze => null()           ! Ze CFAD(Npoints,dBZe_bins,Nlevels)
     real(wp), dimension(:,:),pointer :: &
          lidar_only_freq_cloud => null(),   & ! (Npoints,Nlevels)
          cloudsat_precip_cover => null()      ! Radar total cloud amount by CloudSat precip flag (Npoints,dBZe_bins) <-added
     real(wp),dimension(:),pointer :: &
          radar_lidar_tcc => null(), &         ! Radar&lidar totalcloud amount, grid-box scale (Npoints)
         ! cloudsat_tcc2 => null(),            &
         ! radar_lidar_tcc => null(),          & ! Radar&lidar total cloud amount, grid-box scale (Npoints)
          cloudsat_pia => null()               ! Radar path integrated attenuation (Npoints) <-added 
#ifdef OPT_DPLRW
     real(wp),dimension(:,:),pointer :: gcumw => null()
!     real(wp),dimension(:,:),pointer :: &
!          cldfracT_lid => null(), cldfracT_rad => null(), cldfracT_tot => null()
     
     !real(wp),dimension(:,:,:,:),pointer :: &
     !     dplrw_sumZ => null(), dplrw_sumT => null()
     
     real(wp),dimension(:,:,:,:,:),pointer :: &
          dplrw_Z => null(), spwid_Z => null(), Zef94_Z => null(), &
          dplrw_T => null(), spwid_T => null(), Zef94_T => null(), &
          ZefVd_2 => null()
#endif

     ! ISCCP outputs       
     real(wp),dimension(:),pointer :: &
          isccp_totalcldarea => null(), & ! The fraction of model grid box columns with cloud 
           		                  ! somewhere in them. (%)
          isccp_meantb => null(),       & ! Mean all-sky 10.5 micron brightness temperature. (K)
          isccp_meantbclr => null(),    & ! Mean clear-sky 10.5 micron brightness temperature. (K)
          isccp_meanptop => null(),     & ! Mean cloud top pressure (mb).
          isccp_meantaucld => null(),   & ! Mean optical thickness. (1)
          isccp_meanalbedocld => null()   ! Mean cloud albedo. (1)
     real(wp),dimension(:,:),pointer ::&
          isccp_boxtau => null(),       & ! Optical thickness in each column. (1)
          isccp_boxptop => null()         ! Cloud top pressure in each column. (mb)
     real(wp),dimension(:,:,:),pointer :: &
          isccp_fq  => null()             ! The fraction of the model grid box covered by each of
                                          ! the 49 ISCCP D level cloud types. (%)
     
     ! MISR outputs  			    
     real(wp),dimension(:,:,:),pointer ::   & !
          misr_fq => null()          ! Fraction of the model grid box covered by each of the MISR 
                           ! cloud types
     real(wp),dimension(:,:),pointer ::   & !
          misr_dist_model_layertops => null() !  
     real(wp),dimension(:),pointer ::   & !
          misr_meanztop => null(), & ! Mean MISR cloud top height
          misr_cldarea => null()     ! Mean MISR cloud cover area         			    

     ! MODIS outptus		    
     real(wp),pointer,dimension(:) ::      & !  
          modis_Cloud_Fraction_Total_Mean => null(),       & ! L3 MODIS retrieved cloud fraction (total) 
          modis_Cloud_Fraction_Water_Mean => null(),       & ! L3 MODIS retrieved cloud fraction (liq) 
          modis_Cloud_Fraction_Ice_Mean => null(),         & ! L3 MODIS retrieved cloud fraction (ice) 
          modis_Cloud_Fraction_High_Mean => null(),        & ! L3 MODIS retrieved cloud fraction (high) 
          modis_Cloud_Fraction_Mid_Mean => null(),         & ! L3 MODIS retrieved cloud fraction (middle) 
          modis_Cloud_Fraction_Low_Mean => null(),         & ! L3 MODIS retrieved cloud fraction (low ) 
          modis_Optical_Thickness_Total_Mean => null(),    & ! L3 MODIS retrieved optical thickness (tot)
          modis_Optical_Thickness_Water_Mean => null(),    & ! L3 MODIS retrieved optical thickness (liq)
          modis_Optical_Thickness_Ice_Mean => null(),      & ! L3 MODIS retrieved optical thickness (ice)
          modis_Optical_Thickness_Total_LogMean => null(), & ! L3 MODIS retrieved log10 optical thickness 
          modis_Optical_Thickness_Water_LogMean => null(), & ! L3 MODIS retrieved log10 optical thickness 
          modis_Optical_Thickness_Ice_LogMean => null(),   & ! L3 MODIS retrieved log10 optical thickness
          modis_Cloud_Particle_Size_Water_Mean => null(),  & ! L3 MODIS retrieved particle size (liquid)
          modis_Cloud_Particle_Size_Ice_Mean => null(),    & ! L3 MODIS retrieved particle size (ice)
          modis_Cloud_Top_Pressure_Total_Mean => null(),   & ! L3 MODIS retrieved cloud top pressure
          modis_Liquid_Water_Path_Mean => null(),          & ! L3 MODIS retrieved liquid water path
          modis_Ice_Water_Path_Mean => null()                ! L3 MODIS retrieved ice water path
     real(wp),pointer,dimension(:,:,:) ::  &
          modis_Optical_Thickness_vs_Cloud_Top_Pressure => null(), & ! Tau/Pressure joint histogram          			    
          modis_Optical_Thickness_vs_ReffICE => null(),            & ! Tau/ReffICE joint histogram
          modis_Optical_Thickness_vs_ReffLIQ => null()               ! Tau/ReffLIQ joint histogram

     ! RTTOV outputs
     real(wp),pointer :: &
          rttov_tbs(:,:) => null() ! Brightness Temperature	    

     ! IDiD output
     real(wp),dimension(:,:,:,:),pointer :: &
          cfodd_ntotal => null()       ! # of total samples for CFODD (Npoints,CFODD_NDBZE,CFODD_NICOD,CFODD_NCLASS)
     real(wp),dimension(:,:),    pointer :: &
          wr_occfreq_ntotal => null()  ! # of total samples for warm-rain occurrence frequency of
                                       ! nonprecip/drizzle/precip (Npoints,WR_NREGIME)
     
     ! SCOPS check ! added by YN
     real(wp),pointer :: &
          frac_ls(:,:) => null(), frac_cv(:,:) => null(), &
          prec_ls(:,:) => null(), prec_cv(:,:) => null()

  end type cosp_outputs

CONTAINS
  ! ######################################################################################
  ! FUNCTION cosp_simulator
  ! ######################################################################################
  function COSP_SIMULATOR(cospIN,cospgridIN,cospOUT,start_idx,stop_idx,debug)
    type(cosp_optical_inputs),intent(in),target :: cospIN     ! Optical inputs to COSP simulator
    type(cosp_column_inputs), intent(in),target :: cospgridIN ! Host model inputs to COSP
    
    ! Inputs into the simulators
    type(isccp_IN)    :: isccpIN    ! Input to the ISCCP simulator
    type(misr_IN)     :: misrIN     ! Input to the LIDAR simulator
    type(calipso_IN)  :: calipsoIN  ! Input to the LIDAR simulator
    type(parasol_IN)  :: parasolIN  ! Input to the PARASOL simulator
    type(cloudsat_IN) :: cloudsatIN ! Input to the CLOUDSAT radar simulator
    type(modis_IN)    :: modisIN    ! Input to the MODIS simulator
    type(rttov_IN)    :: rttovIN    ! Input to the RTTOV simulator
    integer,optional  :: start_idx,stop_idx
    logical,optional  :: debug
    
    ! Outputs from the simulators (nested simulator output structure)
    type(cosp_outputs), intent(inout) :: cospOUT
    character(len=256),dimension(100) :: cosp_simulator
    
    ! Local variables
    integer :: &
         i,icol,ij,ik,nError
    integer :: k
    integer,target :: &
         Npoints
    logical :: &
         Lisccp_subcolumn,    & ! On/Off switch for subcolumn ISCCP simulator
         Lmisr_subcolumn,     & ! On/Off switch for subcolumn MISR simulator
         Lcalipso_subcolumn,  & ! On/Off switch for subcolumn CALIPSO simulator
         Lparasol_subcolumn,  & ! On/Off switch for subcolumn PARASOL simulator
         Lcloudsat_subcolumn, & ! On/Off switch for subcolumn CLOUDSAT simulator
         Lmodis_subcolumn,    & ! On/Off switch for subcolumn MODIS simulator
         Lrttov_subcolumn,    & ! On/Off switch for subcolumn RTTOV simulator
         Lisccp_column,       & ! On/Off switch for column ISCCP simulator
         Lmisr_column,        & ! On/Off switch for column MISR simulator
         Lcalipso_column,     & ! On/Off switch for column CALIPSO simulator
         Lparasol_column,     & ! On/Off switch for column PARASOL simulator
         Lcloudsat_column,    & ! On/Off switch for column CLOUDSAT simulator
         Lmodis_column,       & ! On/Off switch for column MODIS simulator
         Lrttov_column,       & ! On/Off switch for column RTTOV simulator (not used)      
         Lradar_lidar_tcc,    & ! On/Off switch from joint Calipso/Cloudsat product
         !Lcloudsat_tcc,       & !
         !Lcloudsat_tcc2,      & ! 
         Llidar_only_freq_cloud, & ! On/Off switch from joint Calipso/Cloudsat product
         Lcloudsat_modis_wr     ! On/Off switch from joint CloudSat/MODIS warm rain product

    logical :: Lscops ! added by YN

#ifdef OPT_DPLRW
    logical :: &
         Ldplrw!            &   ! On/Off switch for column dplrw product
!         Ldplrw_LS,         &
!         Ldplrw_CU
#endif
    
    logical :: &
         ok_lidar_cfad    = .false., &
         !ok_lidar_cfad_grLidar532 = .false., &
         !ok_lidar_cfad_atlid = .false., &
         lrttov_cleanUp   = .false.
    !! CloudSat+MODIS Joint Product I/O control
    logical, save :: ofirst  = .true.
    logical :: &
         Lwarmrain = .false.

    integer, dimension(:,:),allocatable  :: &
         modisRetrievedPhase,isccpLEVMATCH
    real(wp), dimension(:),  allocatable  :: &
         modisCfTotal,modisCfLiquid,modisMeanIceWaterPath, isccp_meantbclr,     &                         
         modisCfIce, modisCfHigh, modisCfMid, modisCfLow,modisMeanTauTotal,     &       
         modisMeanTauLiquid, modisMeanTauIce, modisMeanLogTauTotal,             &       
         modisMeanLogTauLiquid, modisMeanLogTauIce, modisMeanSizeLiquid,        &        
         modisMeanSizeIce, modisMeanCloudTopPressure, modisMeanLiquidWaterPath, &
         radar_lidar_tcc!,cloudsat_tcc, cloudsat_tcc2
    REAL(WP), dimension(:,:),allocatable  :: &
         modisRetrievedCloudTopPressure,modisRetrievedTau,modisRetrievedSize,   &
         misr_boxtau,misr_boxztop,misr_dist_model_layertops,isccp_boxtau,       &
         isccp_boxttop,isccp_boxptop,calipso_beta_mol,lidar_only_freq_cloud
    REAL(WP), dimension(:,:,:),allocatable :: &
         modisJointHistogram,modisJointHistogramIce,modisJointHistogramLiq,     &
         calipso_beta_tot,calipso_betaperp_tot,cloudsatDBZe,parasolPix_refl,cloudsatZe_non !<-added

    real(wp),dimension(:),allocatable,target :: &
         out1D_1,out1D_2,out1D_3,out1D_4,out1D_5,out1D_6
   ! Jing X.-->
    real(wp),dimension(:,:),allocatable,target :: &
         out2D_1,out2D_2
    real(wp),dimension(:,:,:),allocatable,target :: &
         out3D_1,out3D_2,out3D_3,out3D_4,out3D_5,out3D_6,out3D_7,out3D_8
   !   <--
    !! IDiD (T.Michibata)
    real(wp),dimension(:,:,:),allocatable :: &
       t_in,betamol_in,tmpFlip,betamolFlip,pnormFlip,pnorm_perpFlip,ze_totFlip
    real(wp),dimension(:,:,:),allocatable :: &
       tempI, frac_outI, ze_totI, ze_nonI
    real(wp), allocatable :: &
          zlev   (:,:),           & ! altitude (used only when use_vgrid=.true.)
          delz   (:,:),           & ! delta Z
          cfodd_ntotal (:,:,:,:), & ! # of total samples for CFODD (Npoints,CFODD_NDBZE,CFODD_NICOD,CFODD_NCLASS)
          wr_occfreq_ntotal(:,:)    ! # of warm-rain (nonprecip/drizzle/precip) (Npoints,WR_NREGIME)

    ! for SCOPS check ! added by YN
    real(wp),dimension(:),allocatable :: mmar
    logical, dimension(:),allocatable :: mmlg

#ifdef OPT_DPLRW
    ! Indices to address arrays of LS and CONV hydrometeors
!    integer,parameter :: &
!         I_LSCLIQ = 1, & ! Large-scale (stratiform) liquid
!         I_LSCICE = 2, & ! Large-scale (stratiform) ice
!         I_LSRAIN = 3, & ! Large-scale (stratiform) rain
!         I_LSSNOW = 4, & ! Large-scale (stratiform) snow
!         I_CVCLIQ = 5, & ! Convective liquid
!         I_CVCICE = 6, & ! Convective ice
!         I_CVRAIN = 7, & ! Convective rain
!         I_CVSNOW = 8, & ! Convective snow
!         I_LSGRPL = 9    ! Large-scale (stratiform) groupel
    ! Stratiform and convective clouds in frac_out (scops output).
!    integer, parameter :: &
!         I_LSC = 1, & ! Large-scale clouds
!         I_CVC = 2    ! Convective clouds
    ! output control
!    real(wp),dimension(:,:),allocatable,target :: &
!         out2D_3, out2D_4
!    real(wp),dimension(:,:,:),allocatable,target :: &
!         out3D_9, out3D_0, out3D_a, out3D_b, out3D_c, out3D_d, out3D_e, out3D_f
#endif
! add by NEC
    integer ncount

    integer jfpar
    call getjfp(jfpar)

    ! Initialize error reporting for output
    cosp_simulator(:)=''

!    write( jfpar,* ) ' ### test init by T.Michibata. (test No.1)'
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 1) Determine if using full inputs or subset
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (present(start_idx) .and. present(stop_idx)) then
       ij=start_idx
       ik=stop_idx
    else
       ij=1
       ik=cospIN%Npoints
    endif
    Npoints = ik-ij+1
    
!    write( jfpar,* ) ' ### test init by T.Michibata. (test No.2)'
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 2a) Determine which simulators to run and which statistics to compute
    !    - If any of the subcolumn fields are allocated, then run the subcolumn simulators. 
    !    - If any of the column fields are allocated, then compute the statistics for that
    !      simulator, but only save the requested fields.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Start with all simulators and joint-diagnostics off
    Lisccp_subcolumn    = .false.
    Lmisr_subcolumn     = .false.
    Lcalipso_subcolumn  = .false.
    Lparasol_subcolumn  = .false.
    Lcloudsat_subcolumn = .false.
    Lmodis_subcolumn    = .false.
    Lrttov_subcolumn    = .false.
    Lisccp_column       = .false.
    Lmisr_column        = .false.
    Lcalipso_column     = .false.
    Lparasol_column     = .false.
    Lcloudsat_column    = .false.
    Lmodis_column       = .false.
    Lrttov_column       = .false.
    Lradar_lidar_tcc    = .false.
    Llidar_only_freq_cloud = .false.
    !Lcloudsat_tcc       = .false.
    !Lcloudsat_tcc2      = .false.
    Lcloudsat_modis_wr  = .false.
    Lscops              = .false. ! added by YN
#ifdef OPT_DPLRW
    Ldplrw              = .false.
!    Ldplrw_LS           = .false.
!    Ldplrw_CU           = .false.
#endif

    ! CLOUDSAT subcolumn
    if (associated(cospOUT%cloudsat_Ze_tot)) Lcloudsat_subcolumn = .true.

    ! MODIS subcolumn
    if (associated(cospOUT%modis_Cloud_Fraction_Water_Mean)                .or.          &
        associated(cospOUT%modis_Cloud_Fraction_Total_Mean)                .or.          &
        associated(cospOUT%modis_Cloud_Fraction_Ice_Mean)                  .or.          &
        associated(cospOUT%modis_Cloud_Fraction_High_Mean)                 .or.          &
        associated(cospOUT%modis_Cloud_Fraction_Mid_Mean)                  .or.          &
        associated(cospOUT%modis_Cloud_Fraction_Low_Mean)                  .or.          &
        associated(cospOUT%modis_Optical_Thickness_Total_Mean)             .or.          &
        associated(cospOUT%modis_Optical_Thickness_Water_Mean)             .or.          &
        associated(cospOUT%modis_Optical_Thickness_Ice_Mean)               .or.          &
        associated(cospOUT%modis_Optical_Thickness_Total_LogMean)          .or.          &
        associated(cospOUT%modis_Optical_Thickness_Water_LogMean)          .or.          &
        associated(cospOUT%modis_Optical_Thickness_Ice_LogMean)            .or.          &
        associated(cospOUT%modis_Cloud_Particle_Size_Water_Mean)           .or.          &
        associated(cospOUT%modis_Cloud_Particle_Size_Ice_Mean)             .or.          &
        associated(cospOUT%modis_Cloud_Top_Pressure_Total_Mean)            .or.          &
        associated(cospOUT%modis_Liquid_Water_Path_Mean)                   .or.          &
        associated(cospOUT%modis_Ice_Water_Path_Mean)                      .or.          &
        associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure))               &
       Lmodis_subcolumn    = .true.

    ! ISCCP subcolumn
    if (associated(cospOUT%isccp_boxtau)                                   .or.          &
        associated(cospOUT%isccp_boxptop))                                               &                 
       Lisccp_subcolumn    = .true.

    ! MISR subcolumn
    if (associated(cospOUT%misr_dist_model_layertops))                                   &
       Lmisr_subcolumn     = .true.

    ! CALIPOSO subcolumn
    if (associated(cospOUT%calipso_tau_tot)                                .or.          &
        associated(cospOUT%calipso_beta_mol)                               .or.          &
        associated(cospOUT%calipso_temp_tot)                               .or.          &
        associated(cospOUT%calipso_betaperp_tot)                           .or.          &
        associated(cospOUT%calipso_beta_tot))                                            &
       Lcalipso_subcolumn  = .true.

    ! PARASOL subcolumn
    if (associated(cospOUT%parasolPix_refl))                                             &
       Lparasol_subcolumn  = .true.

    ! RTTOV column
    if (associated(cospOUT%rttov_tbs))                                                   &
       Lrttov_column    = .true.

    ! Set flag to deallocate rttov types (only done on final call to simulator)
    if (size(cospOUT%isccp_meantb) .eq. stop_idx) lrttov_cleanUp = .true.    
    
    ! ISCCP column
    if (associated(cospOUT%isccp_fq)                                       .or.          &
        associated(cospOUT%isccp_meanalbedocld)                            .or.          &
        associated(cospOUT%isccp_meanptop)                                 .or.          &
        associated(cospOUT%isccp_meantaucld)                               .or.          &
        associated(cospOUT%isccp_totalcldarea)                             .or.          &
        associated(cospOUT%isccp_meantb)) then
       Lisccp_column    = .true.             
       Lisccp_subcolumn = .true.
    endif

    ! MISR column
    if (associated(cospOUT%misr_cldarea)                                   .or.          &
        associated(cospOUT%misr_meanztop)                                  .or.          &
        associated(cospOUT%misr_fq)) then
       Lmisr_column    = .true.
       Lmisr_subcolumn = .true.
    endif

    ! CALIPSO column
    if (associated(cospOUT%calipso_cfad_sr)                                .or.          &
        associated(cospOUT%calipso_lidarcld)                               .or.          &
        associated(cospOUT%calipso_lidarcldphase)                          .or.          &
       ! associated(cospOUT%calipso_lidarcldtype)                           .or.          &
        associated(cospOUT%calipso_cldlayer)                               .or.          &
        associated(cospOUT%calipso_cldlayerphase)                          .or.          &
        associated(cospOUT%calipso_lidarcldtmp)) then
       Lcalipso_column    = .true.
       Lcalipso_subcolumn = .true.
    endif

    ! PARASOL column
    if (associated(cospOUT%parasolGrid_refl)) then
       Lparasol_column    = .true.
       Lparasol_subcolumn = .true.
    endif

    ! CLOUDSAT column
    if (associated(cospOUT%cloudsat_cfad_ze) .or. associated(cospOUT%cloudsat_precip_cover) .or. &
        associated(cospOUT%cloudsat_pia)) then
       Lcloudsat_column    = .true.
       Lcloudsat_subcolumn = .true.
    endif

#ifdef OPT_DPLRW
!    if (associated(cospOUT%dplrw_LSZ)) Ldplrw_LS = .true.
!    if (associated(cospOUT%dplrw_CUZ)) Ldplrw_CU = .true.

    if (associated(cospOUT%dplrw_Z)) Ldplrw = .true.
    if (Ldplrw) Lisccp_subcolumn = .true.
    if (Ldplrw) Lcloudsat_subcolumn = .true.
#endif

    ! MODIS column
    if (associated(cospOUT%modis_Cloud_Fraction_Total_Mean)                .or.          &
        associated(cospOUT%modis_Cloud_Fraction_Water_Mean)                .or.          &
        associated(cospOUT%modis_Cloud_Fraction_Ice_Mean)                  .or.          &
        associated(cospOUT%modis_Cloud_Fraction_High_Mean)                 .or.          &
        associated(cospOUT%modis_Cloud_Fraction_Mid_Mean)                  .or.          &
        associated(cospOUT%modis_Cloud_Fraction_Low_Mean)                  .or.          &
        associated(cospOUT%modis_Optical_Thickness_Total_Mean)             .or.          &
        associated(cospOUT%modis_Optical_Thickness_Water_Mean)             .or.          &
        associated(cospOUT%modis_Optical_Thickness_Ice_Mean)               .or.          &
        associated(cospOUT%modis_Optical_Thickness_Total_LogMean)          .or.          &
        associated(cospOUT%modis_Optical_Thickness_Water_LogMean)          .or.          &
        associated(cospOUT%modis_Optical_Thickness_Ice_LogMean)            .or.          &
        associated(cospOUT%modis_Cloud_Particle_Size_Water_Mean)           .or.          & 
        associated(cospOUT%modis_Cloud_Particle_Size_Ice_Mean)             .or.          &
        associated(cospOUT%modis_Cloud_Top_Pressure_Total_Mean)            .or.          &
        associated(cospOUT%modis_Liquid_Water_Path_Mean)                   .or.          &
        associated(cospOUT%modis_Ice_Water_Path_Mean)                      .or.          &
        associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure)) then
       Lmodis_column    = .true.
       Lmodis_subcolumn = .true.
    endif

    ! Joint simulator products
    if (associated(cospOUT%lidar_only_freq_cloud) .or. associated(cospOUT%radar_lidar_tcc)) then  !.or. &
       ! associated(cospOUT%cloudsat_tcc) .or. associated(cospOUT%cloudsat_tcc2)) then 
       Lcalipso_column     = .true.
       Lcalipso_subcolumn  = .true.
       Lcloudsat_column    = .true.
       Lcloudsat_subcolumn = .true.
       Lradar_lidar_tcc    = .true.
       Llidar_only_freq_cloud = .true.
       !Lcloudsat_tcc       = .true.
       !Lcloudsat_tcc2      = .true.
    endif

    ! CloudSat+MODIS joint simulator product
    if ( associated(cospOUT%cfodd_ntotal) .or. associated(cospOUT%wr_occfreq_ntotal) ) then
       Lmodis_column       = .true.
       Lmodis_subcolumn    = .true.
       Lcloudsat_column    = .true.
       Lcloudsat_subcolumn = .true.
       Lcloudsat_modis_wr  = .true. ! WR: warm rain product
    endif

    ! SCOPS check ! added by YN
    if (associated(cospOUT%frac_ls)) Lscops = .true.

!    write( jfpar,* ) ' ### test init by T.Michibata. (test No.2b)'
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 2b) Error Checking
    !     Enforce bounds on input fields. If input field is out-of-bounds, report error 
    !     and turn off simulator
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    call cosp_errorCheck(cospgridIN,cospIN,Lisccp_subcolumn,Lisccp_column,               &
                         Lmisr_subcolumn,Lmisr_column,Lmodis_subcolumn,Lmodis_column,    &
                         Lcloudsat_subcolumn,Lcloudsat_column,Lcalipso_subcolumn,        &
                         Lcalipso_column,Lrttov_subcolumn,Lrttov_column,                 &
                         Lparasol_subcolumn,Lparasol_column,Lradar_lidar_tcc,            &
                         Llidar_only_freq_cloud,Lcloudsat_modis_wr,                       &
                         !,Lcloudsat_tcc,Lcloudsat_tcc2,            &
                         cospOUT,cosp_simulator,nError)

!    write( jfpar,* ) ' ### test init by T.Michibata. (test No.3)'
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 3) Populate instrument simulator inputs
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (Lisccp_subcolumn .or. Lmodis_subcolumn) then
       isccpIN%Npoints  => Npoints
       isccpIN%Ncolumns => cospIN%Ncolumns
       isccpIN%Nlevels  => cospIN%Nlevels
       isccpIN%emsfc_lw => cospIN%emsfc_lw
       isccpIN%skt      => cospgridIN%skt
       isccpIN%qv       => cospgridIN%qv
       isccpIN%at       => cospgridIN%at
       isccpIN%frac_out => cospIN%frac_out
       isccpIN%dtau     => cospIN%tau_067
       isccpIN%dem      => cospIN%emiss_11
       isccpIN%phalf    => cospgridIN%phalf
       isccpIN%sunlit   => cospgridIN%sunlit
       isccpIN%pfull    => cospgridIN%pfull
    endif
    
    if (Lmisr_subcolumn) then
       misrIN%Npoints  => Npoints
       misrIN%Ncolumns => cospIN%Ncolumns
       misrIN%Nlevels  => cospIN%Nlevels
       misrIN%dtau     => cospIN%tau_067
       misrIN%sunlit   => cospgridIN%sunlit
       misrIN%zfull    => cospgridIN%hgt_matrix
       misrIN%at       => cospgridIN%at
    endif
    
    if (Lcalipso_subcolumn) then
       calipsoIN%Npoints     => Npoints
       calipsoIN%Ncolumns    => cospIN%Ncolumns
       calipsoIN%Nlevels     => cospIN%Nlevels
       calipsoIN%beta_mol    => cospIN%beta_mol
       calipsoIN%betatot     => cospIN%betatot
       calipsoIN%betatot_liq => cospIN%betatot_liq
       calipsoIN%betatot_ice => cospIN%betatot_ice
       calipsoIN%tau_mol     => cospIN%tau_mol
       calipsoIN%tautot      => cospIN%tautot
       calipsoIN%tautot_liq  => cospIN%tautot_liq
       calipsoIN%tautot_ice  => cospIN%tautot_ice
    endif
    
    if (Lparasol_subcolumn) then
       parasolIN%Npoints      => Npoints
       parasolIN%Nlevels      => cospIN%Nlevels
       parasolIN%Ncolumns     => cospIN%Ncolumns
       parasolIN%Nrefl        => cospIN%Nrefl
       parasolIN%tautot_S_liq => cospIN%tautot_S_liq
       parasolIN%tautot_S_ice => cospIN%tautot_S_ice
    endif

    if (Lcloudsat_subcolumn) then
       cloudsatIN%Npoints    => Npoints
       cloudsatIN%Nlevels    => cospIN%Nlevels
       cloudsatIN%Ncolumns   => cospIN%Ncolumns
       cloudsatIN%z_vol      => cospIN%z_vol_cloudsat
       cloudsatIN%kr_vol     => cospIN%kr_vol_cloudsat
       cloudsatIN%g_vol      => cospIN%g_vol_cloudsat
       cloudsatIN%rcfg       => cospIN%rcfg_cloudsat
       cloudsatIN%hgt_matrix => cospgridIN%hgt_matrix
    endif
    
    if (Lmodis_subcolumn) then
       modisIN%Ncolumns  => cospIN%Ncolumns
       modisIN%Nlevels   => cospIN%Nlevels
       modisIN%Npoints   => Npoints
       modisIN%liqFrac   => cospIN%fracLiq
       modisIN%tau       => cospIN%tau_067
       modisIN%g         => cospIN%asym
       modisIN%w0        => cospIN%ss_alb
       modisIN%Nsunlit   = count(cospgridIN%sunlit > 0)
       if (modisIN%Nsunlit .gt. 0) then
          allocate(modisIN%sunlit(modisIN%Nsunlit),modisIN%pres(modisIN%Nsunlit,cospIN%Nlevels+1))
!          modisIN%sunlit    = pack((/ (i, i = 2, Npoints ) /),mask = cospgridIN%sunlit > 0)
! mod by NEC
          ncount = 1
          do i = 1, Npoints
             if (cospgridIN%sunlit(i) > 0) then
               if (ncount <= size(modisIN%sunlit)) then
                 modisIN%sunlit(ncount) = i
                 ncount = ncount + 1
               end if
             end if
          enddo
          modisIN%pres      = cospgridIN%phalf(int(modisIN%sunlit(:)),:)
       endif
       if (count(cospgridIN%sunlit <= 0) .gt. 0) then
          allocate(modisIN%notSunlit(count(cospgridIN%sunlit <= 0)))
!          modisIN%notSunlit = pack((/ (i, i = 1, Npoints ) /),mask = .not. cospgridIN%sunlit > 0)
! mod by NEC
          ncount = 1
          do i = 1, Npoints
            if ( .not. (cospgridIN%sunlit(i) > 0) ) then
              if (ncount <= size(modisIN%notSunlit)) then
                modisIN%notSunlit(ncount) = i
                ncount = ncount + 1
              endif
            endif
          enddo
       endif
    endif
    
    if (Lrttov_column) then
       rttovIN%nPoints    => Npoints
       rttovIN%nLevels    => cospIN%nLevels
       rttovIN%nSubCols   => cospIN%nColumns
       rttovIN%zenang     => cospgridIN%zenang
       rttovIN%co2        => cospgridIN%co2
       rttovIN%ch4        => cospgridIN%ch4
       rttovIN%n2o        => cospgridIN%n2o
       rttovIN%co         => cospgridIN%co
       rttovIN%surfem     => cospgridIN%emis_sfc
       rttovIN%h_surf     => cospgridIN%hgt_matrix_half(:,cospIN%Nlevels+1)
       rttovIN%u_surf     => cospgridIN%u_sfc
       rttovIN%v_surf     => cospgridIN%v_sfc
       rttovIN%t_skin     => cospgridIN%skt
       rttovIN%p_surf     => cospgridIN%phalf(:,cospIN%Nlevels+1)
       rttovIN%q2m        => cospgridIN%qv(:,cospIN%Nlevels)
       rttovIN%t2m        => cospgridIN%at(:,cospIN%Nlevels)
       rttovIN%lsmask     => cospgridIN%land
       rttovIN%latitude   => cospgridIN%lat
       rttovIN%longitude  => cospgridIN%lon
       rttovIN%seaice     => cospgridIN%seaice
       rttovIN%p          => cospgridIN%pfull
       rttovIN%ph         => cospgridIN%phalf
       rttovIN%t          => cospgridIN%at
       rttovIN%q          => cospgridIN%qv
       rttovIN%o3         => cospgridIN%o3
       ! Below only needed for all-sky RTTOV calculation
       rttovIN%month      => cospgridIN%month
       rttovIN%tca        => cospgridIN%tca
       rttovIN%cldIce     => cospgridIN%cloudIce
       rttovIN%cldLiq     => cospgridIN%cloudLiq
       rttovIN%fl_rain    => cospgridIN%fl_rain
       rttovIN%fl_snow    => cospgridIN%fl_snow
    endif

!    write( jfpar,* ) ' ### test init by T.Michibata. (test No.4)'
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 4) Call subcolumn simulators
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! ISCCP (icarus) subcolumn simulator
    if (Lisccp_subcolumn .or. Lmodis_subcolumn) then
       ! Allocate space for local variables
       allocate(isccpLEVMATCH(Npoints,isccpIN%Ncolumns),                                 &
                isccp_boxttop(Npoints,isccpIN%Ncolumns),                                 &
                isccp_boxptop(Npoints,isccpIN%Ncolumns),                                 &
                isccp_boxtau(Npoints,isccpIN%Ncolumns), isccp_meantbclr(Npoints))
       ! Call simulator
       call icarus_subcolumn(isccpIN%npoints,isccpIN%ncolumns,isccpIN%nlevels,           &
                             isccpIN%sunlit,isccpIN%dtau,isccpIN%dem,isccpIN%skt,        &
                             isccpIN%emsfc_lw,isccpIN%qv,isccpIN%at,isccpIN%pfull,       &
                             isccpIN%phalf,isccpIN%frac_out,isccpLEVMATCH,               &
                             isccp_boxtau(:,:),isccp_boxptop(:,:),                       &
                             isccp_boxttop(:,:),isccp_meantbclr(:))
       ! Store output (if requested)
       if (associated(cospOUT%isccp_boxtau)) then
          cospOUT%isccp_boxtau(ij:ik,:)  = isccp_boxtau
       endif
       if (associated(cospOUT%isccp_boxptop)) then
          cospOUT%isccp_boxptop(ij:ik,:) = isccp_boxptop
       endif
       if (associated(cospOUT%isccp_meantbclr)) then
          cospOUT%isccp_meantbclr(ij:ik) = isccp_meantbclr
       endif
   endif

   ! MISR subcolumn simulator
    if (Lmisr_subcolumn) then
       ! Allocate space for local variables
       allocate(misr_boxztop(Npoints,misrIN%Ncolumns),                                   &
                misr_boxtau(Npoints,misrIN%Ncolumns),                                    &
                misr_dist_model_layertops(Npoints,numMISRHgtBins))
       ! Call simulator
       call misr_subcolumn(misrIN%Npoints,misrIN%Ncolumns,misrIN%Nlevels,misrIN%dtau,    &
                           misrIN%zfull,misrIN%at,misrIN%sunlit,misr_boxtau,             &
                           misr_dist_model_layertops,misr_boxztop)
       ! Store output (if requested)
       if (associated(cospOUT%misr_dist_model_layertops)) then
          cospOUT%misr_dist_model_layertops(ij:ik,:) = misr_dist_model_layertops
       endif
    endif

    ! Calipso subcolumn simulator
    if (Lcalipso_subcolumn) then
       ! Allocate space for local variables
       allocate(calipso_beta_mol(calipsoIN%Npoints,calipsoIN%Nlevels),                   &
                calipso_beta_tot(calipsoIN%Npoints,calipsoIN%Ncolumns,calipsoIN%Nlevels),&
                calipso_betaperp_tot(calipsoIN%Npoints,calipsoIN%Ncolumns,calipsoIN%Nlevels))
       ! Call simulator
      ! call lidar_subcolumn(calipsoIN%npoints, calipsoIN%ncolumns, calipsoIN%nlevels, .false., &
      !                      calipsoIN%beta_mol, calipsoIN%tau_mol, calipsoIN%betatot, calipsoIN%tautot,  &
      !                      calipso_beta_mol(:,:), calipso_beta_tot(:,:,:), calipsoIN%betatot_ice,       &
      !                      calipsoIN%tautot_ice, calipsoIN%betatot_liq, calipsoIN%tautot_liq,           &
      !                      calipso_betaperp_tot(:,:,:)) 
      call lidar_subcolumn(calipsoIN%npoints,calipsoIN%ncolumns,calipsoIN%nlevels,      &
                            calipsoIN%beta_mol,calipsoIN%tau_mol,                       &
                            calipsoIN%betatot,calipsoIN%tautot,calipsoIN%betatot_ice,   &
                            calipsoIN%tautot_ice,calipsoIN%betatot_liq,                 &
                            calipsoIN%tautot_liq,calipso_beta_mol(:,:),                 &
                            calipso_beta_tot(:,:,:),calipso_betaperp_tot(:,:,:))
       ! Store output (if requested)
       if (associated(cospOUT%calipso_beta_mol))                                         &
            cospOUT%calipso_beta_mol(ij:ik,calipsoIN%Nlevels:1:-1) = calipso_beta_mol
       if (associated(cospOUT%calipso_beta_tot))                                         &
            cospOUT%calipso_beta_tot(ij:ik,:,calipsoIN%Nlevels:1:-1) = calipso_beta_tot
       if (associated(cospOUT%calipso_betaperp_tot))                                     &
            cospOUT%calipso_betaperp_tot(ij:ik,:,:) = calipso_betaperp_tot

    endif

    ! PARASOL subcolumn simulator
    if (Lparasol_subcolumn) then
       ! Allocate space for local variables
       allocate(parasolPix_refl(parasolIN%Npoints,parasolIN%Ncolumns,PARASOL_NREFL))
       ! Call simulator
       do icol=1,parasolIN%Ncolumns
          call parasol_subcolumn(parasolIN%npoints, PARASOL_NREFL,                       &
                                 parasolIN%tautot_S_liq(1:parasolIN%Npoints,icol),       &
                                 parasolIN%tautot_S_ice(1:parasolIN%Npoints,icol),       &
                                 parasolPix_refl(:,icol,1:PARASOL_NREFL))
          ! Store output (if requested)
          if (associated(cospOUT%parasolPix_refl)) then
             cospOUT%parasolPix_refl(ij:ik,icol,1:PARASOL_NREFL) =                          &
                  parasolPix_refl(:,icol,1:PARASOL_NREFL)
          endif
       enddo
    endif    

    ! Cloudsat (quickbeam) subcolumn simulator
    if (Lcloudsat_subcolumn) then
       ! Allocate space for local variables
        allocate(cloudsatDBZe(cloudsatIN%Npoints,cloudsatIN%Ncolumns,cloudsatIN%Nlevels))
       if (allocated(cloudsatZe_non)) then
         continue
       else
         allocate(cloudsatZe_non(cloudsatIN%Npoints,cloudsatIN%Ncolumns,cloudsatIN%Nlevels))  !<-added
       endif 
       do icol=1,cloudsatIN%ncolumns
          call quickbeam_subcolumn(cloudsatIN%rcfg,cloudsatIN%Npoints,cloudsatIN%Nlevels,&
                                   cloudsatIN%hgt_matrix/1000._wp,                       &
                                   cloudsatIN%z_vol(:,icol,:),                           &
                                   cloudsatIN%kr_vol(:,icol,:),                          &
                                   cloudsatIN%g_vol(:,1,:),cloudsatDBze(:,icol,:),cloudsatZe_non(:,icol,:)) !<-added
       enddo
       ! Store output (if requested)
       if (associated(cospOUT%cloudsat_Ze_tot)) then
          cospOUT%cloudsat_Ze_tot(ij:ik,:,:) = cloudsatDBZe(:,:,1:cloudsatIN%Nlevels) !<-added(change)
       endif
    endif

    if (Lmodis_subcolumn) then
       if(modisiN%nSunlit > 0) then 
          ! Allocate space for local variables
          allocate(modisRetrievedTau(modisIN%nSunlit,modisIN%nColumns),                  &
                   modisRetrievedSize(modisIN%nSunlit,modisIN%nColumns),                 &
                   modisRetrievedPhase(modisIN%nSunlit,modisIN%nColumns),                &
                   modisRetrievedCloudTopPressure(modisIN%nSunlit,modisIN%nColumns))
          ! Call simulator
          do i = 1, modisIN%nSunlit
             call modis_subcolumn(modisIN%Ncolumns,modisIN%Nlevels,modisIN%pres(i,:),    &
                                  modisIN%tau(int(modisIN%sunlit(i)),:,:),               &
                                  modisIN%liqFrac(int(modisIN%sunlit(i)),:,:),           &
                                  modisIN%g(int(modisIN%sunlit(i)),:,:),                 &
                                  modisIN%w0(int(modisIN%sunlit(i)),:,:),                &
                                  isccp_boxptop(int(modisIN%sunlit(i)),:),               &
                                  modisRetrievedPhase(i,:),                              &
                                  modisRetrievedCloudTopPressure(i,:),                   &
                                  modisRetrievedTau(i,:),modisRetrievedSize(i,:))
          end do
       endif
    endif

!    write( jfpar,* ) ' ### test init by T.Michibata. (test No.5)'
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 5) Call column simulators
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    ! ISCCP
    if (Lisccp_column) then
       ! Check to see which outputs are requested. If not requested, use a local dummy array
       if(.not. associated(cospOUT%isccp_meanalbedocld)) then
          allocate(out1D_1(Npoints))
          !cospOUT%isccp_meanalbedocld(ij:ik) => out1D_1
          cospOUT%isccp_meanalbedocld => out1D_1 !Jing X.
       endif
       if(.not. associated(cospOUT%isccp_meanptop)) then
          allocate(out1D_2(Npoints))
          !cospOUT%isccp_meanptop(ij:ik) => out1D_2
          cospOUT%isccp_meanptop => out1D_2  !Jing X.
       endif
       if(.not. associated(cospOUT%isccp_meantaucld)) then
          allocate(out1D_3(Npoints))
          !cospOUT%isccp_meantaucld(ij:ik) => out1D_3
          cospOUT%isccp_meantaucld => out1D_3  !Jing X.
       endif   
       if(.not. associated(cospOUT%isccp_totalcldarea)) then
          allocate(out1D_4(Npoints))
          !cospOUT%isccp_totalcldarea(ij:ik) => out1D_4
          cospOUT%isccp_totalcldarea => out1D_4  !Jing X.
       endif
       if(.not. associated(cospOUT%isccp_meantb)) then
          allocate(out1D_5(Npoints))
          !cospOUT%isccp_meantb(ij:ik) => out1D_5    
          cospOUT%isccp_meantb => out1D_5  !Jing X.    
       endif
       if(.not. associated(cospOUT%isccp_fq)) then
          !allocate(out1D_6(Npoints*numISCCPTauBins*numISCCPPresBins))       
          !cospOUT%isccp_fq(ij:ik,1:numISCCPTauBins,1:numISCCPPresBins) => out1D_6
          allocate(out3D_1(Npoints,numISCCPTauBins,numISCCPPresBins))
          cospOUT%isccp_fq => out3D_1  !Jing X.
       endif   
                                
       ! Call simulator
       call icarus_column(isccpIN%npoints, isccpIN%ncolumns,isccp_boxtau(:,:),           &
                          isccp_boxptop(:,:)/100._wp, isccpIN%sunlit,isccp_boxttop,      &
                          cospOUT%isccp_fq(ij:ik,:,:),                                   &
                          cospOUT%isccp_meanalbedocld(ij:ik),                            &
                          cospOUT%isccp_meanptop(ij:ik),cospOUT%isccp_meantaucld(ij:ik), &
                          cospOUT%isccp_totalcldarea(ij:ik),cospOUT%isccp_meantb(ij:ik))
       cospOUT%isccp_fq(ij:ik,:,:) = cospOUT%isccp_fq(ij:ik,:,7:1:-1)
       
       ! Check if there is any value slightly greater than 1
       where ((cospOUT%isccp_totalcldarea > 1.0-1.e-5) .and.                             &
              (cospOUT%isccp_totalcldarea < 1.0+1.e-5))
              cospOUT%isccp_totalcldarea = 1.0
       endwhere
       
       ! Clear up memory (if necessary)
       if (allocated(isccp_boxttop))   deallocate(isccp_boxttop)
       ! modified by Doppler simulator
       !if (allocated(isccp_boxptop))   deallocate(isccp_boxptop)
       !if (allocated(isccp_boxtau))    deallocate(isccp_boxtau)
       if (allocated(isccp_meantbclr)) deallocate(isccp_meantbclr)
       if (allocated(isccpLEVMATCH))   deallocate(isccpLEVMATCH)
       if (allocated(out1D_1)) then
          deallocate(out1D_1)
          nullify(cospOUT%isccp_meanalbedocld)
       endif
       if (allocated(out1D_2)) then
          deallocate(out1D_2)
          nullify(cospOUT%isccp_meanptop)
       endif
       if (allocated(out1D_3)) then
          deallocate(out1D_3)
          nullify(cospOUT%isccp_meantaucld)
       endif
       if (allocated(out1D_4)) then
          deallocate(out1D_4)
          nullify(cospOUT%isccp_totalcldarea)
       endif
       if (allocated(out1D_5)) then
          deallocate(out1D_5)
          nullify(cospOUT%isccp_meantb)
       endif
       !if (allocated(out1D_6)) then
       !   deallocate(out1D_6)
       !   nullify(cospOUT%isccp_fq)
       !endif
       if (allocated(out3D_1)) then !Jing X.
          deallocate(out3D_1)
          nullify(cospOUT%isccp_fq)
       endif
    endif
    
    ! MISR
    if (Lmisr_column) then
       ! Check to see which outputs are requested. If not requested, use a local dummy array
       if (.not. associated(cospOUT%misr_cldarea)) then
          allocate(out1D_1(Npoints))
          !cospOUT%misr_cldarea(ij:ik) => out1D_1              
          cospOUT%misr_cldarea => out1D_1    !Jing X.             
       endif
       if (.not. associated(cospOUT%misr_meanztop)) then 
          allocate(out1D_2(Npoints))
          !cospOUT%misr_meanztop(ij:ik) => out1D_2
          cospOUT%misr_meanztop => out1D_2  !Jing X.
       endif
       if (.not. associated(cospOUT%misr_fq)) then
          !allocate(out1D_3(Npoints*numMISRTauBins*numMISRHgtBins))
          !cospOUT%misr_fq(ij:ik,1:numMISRTauBins,1:numMISRHgtBins) => out1D_3     
          allocate(out3D_1(Npoints,numMISRTauBins,numMISRHgtBins))  !Jing X.
          cospOUT%misr_fq => out3D_1     
        endif   
    
       ! Call simulator
        call misr_column(misrIN%Npoints,misrIN%Ncolumns,misr_boxztop,misrIN%sunlit,&
                         misr_boxtau,cospOUT%misr_cldarea(ij:ik),                  &
                         cospOUT%misr_meanztop(ij:ik),cospOUT%misr_fq(ij:ik,:,:))              

       ! Clear up memory
       if (allocated(misr_boxtau))               deallocate(misr_boxtau)
       if (allocated(misr_boxztop))              deallocate(misr_boxztop)
       if (allocated(misr_dist_model_layertops)) deallocate(misr_dist_model_layertops)   
       if (allocated(out1D_1)) then
          deallocate(out1D_1)
          nullify(cospOUT%misr_cldarea)
       endif
       if (allocated(out1D_2)) then
          deallocate(out1D_2)
          nullify(cospOUT%misr_meanztop)
       endif
       !if (allocated(out1D_3)) then
       !   deallocate(out1D_3)
       !   nullify(cospOUT%misr_fq)
       !endif
     ! Jing X. -->
       if (allocated(out3D_1)) then
          deallocate(out3D_1)
          nullify(cospOUT%misr_fq)
       endif
     ! <--
    endif
    
    ! CALIPSO LIDAR Simulator
    if (Lcalipso_column) then
       ! Check to see which outputs are requested. If not requested, use a local dummy array
       if (.not. associated(cospOUT%calipso_cfad_sr)) then
          !allocate(out1D_1(Npoints*SR_BINS*Nlvgrid))
          !cospOUT%calipso_cfad_sr(ij:ik,1:SR_BINS,1:Nlvgrid) => out1D_1
          allocate(out3D_1(Npoints,SR_BINS,Nlvgrid))  !Jing X.
          cospOUT%calipso_cfad_sr => out3D_1
       endif
       if (.not. associated(cospOUT%calipso_lidarcld)) then
          !allocate(out1D_2(Npoints*Nlvgrid))
          !cospOUT%calipso_lidarcld(ij:ik,1:Nlvgrid) => out1D_2
          allocate(out2D_1(Npoints,Nlvgrid)) !Jing X.
          cospOUT%calipso_lidarcld => out2D_1
       endif   
       if (.not. associated(cospOUT%calipso_lidarcldphase)) then
          !allocate(out1D_3(Npoints*Nlvgrid*6))
          !cospOUT%calipso_lidarcldphase(ij:ik,1:Nlvgrid,1:6) => out1D_3
          allocate(out3D_2(Npoints,Nlvgrid,6)) !Jing X.
          cospOUT%calipso_lidarcldphase => out3D_2
       endif
       if (.not. associated(cospOUT%calipso_cldlayer)) then
          !allocate(out1D_4(Npoints*LIDAR_NCAT))
          !cospOUT%calipso_cldlayer(ij:ik,1:LIDAR_NCAT) => out1D_4
          allocate(out2D_2(Npoints,LIDAR_NCAT)) !Jing X.
          cospOUT%calipso_cldlayer => out2D_2
       endif
       if (.not. associated(cospOUT%calipso_cldlayerphase)) then
          !allocate(out1D_5(Npoints*LIDAR_NCAT*6))
          !cospOUT%calipso_cldlayerphase(ij:ik,1:LIDAR_NCAT,1:6) => out1D_5
          allocate(out3D_3(Npoints,LIDAR_NCAT,6)) !Jing X.
          cospOUT%calipso_cldlayerphase => out3D_3
       endif   
       if (.not. associated(cospOUT%calipso_lidarcldtmp)) then
          !allocate(out1D_6(Npoints*40*5))
          !cospOUT%calipso_lidarcldtmp(ij:ik,1:40,1:5) => out1D_6
          allocate(out3D_4(Npoints,40,5)) !Jing X.
          cospOUT%calipso_lidarcldtmp => out3D_4
       endif  
      ! if (.not. associated(cospOUT%calipso_lidarcldtype)) then
      !    allocate(out1D_7(Npoints*Nlvgrid*4))
      !    cospOUT%calipso_lidarcldtype(ij:ik,1:Nlvgrid,1:4) => out1D_7
      ! endif
      ! if (.not. associated(cospOUT%calipso_cldtype)) then
      !    allocate(out1D_8(Npoints*LIDAR_NTYPE))
      !    cospOUT%calipso_cldtype(ij:ik,1:LIDAR_NTYPE) => out1D_8
      ! endif
      ! if (.not. associated(cospOUT%calipso_cldtypetemp)) then
      !    allocate(out1D_9(Npoints*LIDAR_NTYPE))
      !    cospOUT%calipso_cldtypetemp(ij:ik,1:LIDAR_NTYPE) => out1D_9
      ! endif
      ! if (.not. associated(cospOUT%calipso_cldtypemeanz)) then
      !    allocate(out1D_10(Npoints*2))
      !    cospOUT%calipso_cldtypemeanz(ij:ik,1:2) => out1D_10
      ! endif
      ! if (.not. associated(cospOUT%calipso_cldtypemeanzse)) then
      !    allocate(out1D_12(Npoints*3))
      !    cospOUT%calipso_cldtypemeanzse(ij:ik,1:3) => out1D_12
      ! endif
      ! if (.not. associated(cospOUT%calipso_cldthinemis)) then
      !    allocate(out1D_11(Npoints))
      !    cospOUT%calipso_cldthinemis(ij:ik) => out1D_11
      ! endif 
       
       ! Call simulator
       ok_lidar_cfad=.true.
       call lidar_column(calipsoIN%Npoints,calipsoIN%Ncolumns,calipsoIN%Nlevels,        &
                         Nlvgrid,SR_BINS,cospgridIN%at(:,:),                            &
                         calipso_beta_tot(:,:,:),calipso_betaperp_tot(:,:,:),           &
                         calipso_beta_mol(:,:),                                         &
                         cospgridIN%phalf(:,2:calipsoIN%Nlevels),ok_lidar_cfad,         &
                         LIDAR_NCAT,cospOUT%calipso_cfad_sr(ij:ik,:,:),                 &
                         cospOUT%calipso_lidarcld(ij:ik,:),                             &
                         cospOUT%calipso_lidarcldphase(ij:ik,:,:),                      &
                         cospOUT%calipso_cldlayer(ij:ik,:),                             &
                         cospgridIN%hgt_matrix,cospgridIN%hgt_matrix_half,              &
                         cospOUT%calipso_cldlayerphase(ij:ik,:,:),                      &
                         cospOUT%calipso_lidarcldtmp(ij:ik,:,:))
          !  cospOUT%calipso_lidarcldtype(ij:ik,:,:), cospOUT%calipso_cldtype(ij:ik,:),  &
          !  cospOUT%calipso_cldtypetemp(ij:ik,:), cospOUT%calipso_cldtypemeanz(ij:ik,:), &
          !  cospOUT%calipso_cldtypemeanzse(ij:ik,:), cospOUT%calipso_cldthinemis(ij:ik), &
          !  cospOUT%calipso_cldlayerphase(ij:ik,:,:), cospOUT%calipso_lidarcldtmp(ij:ik,:,:))

       if (associated(cospOUT%calipso_srbval)) cospOUT%calipso_srbval = calipso_histBsct
 
       ! Free up memory (if necessary)
       !if (allocated(out1D_1)) then
       !   deallocate(out1D_1)
       if (allocated(out3D_1)) then  !Jing X.
          deallocate(out3D_1)
          nullify(cospOUT%calipso_cfad_sr)
       endif
       !if (allocated(out1D_2)) then
       !   deallocate(out1D_2)
       if (allocated(out2D_1)) then  !Jing X.
          deallocate(out2D_1)
          nullify(cospOUT%calipso_lidarcld)
       endif
       !if (allocated(out1D_3)) then
       !   deallocate(out1D_3)
       if (allocated(out3D_2)) then  !Jing X.
          deallocate(out3D_2)
          nullify(cospOUT%calipso_lidarcldphase)
       endif
       !if (allocated(out1D_4)) then
       !   deallocate(out1D_4)
       if (allocated(out2D_2)) then  !Jing X.
          deallocate(out2D_2)
          nullify(cospOUT%calipso_cldlayer)
       endif
       !if (allocated(out1D_5)) then
       !   deallocate(out1D_5)
       if (allocated(out3D_3)) then  !Jing X.
          deallocate(out3D_3)
          nullify(cospOUT%calipso_cldlayerphase)
       endif
       !if (allocated(out1D_6)) then
       !   deallocate(out1D_6)
       if (allocated(out3D_4)) then   !Jing X.
          deallocate(out3D_4)
          nullify(cospOUT%calipso_lidarcldtmp)
       endif       
    !   if (allocated(out1D_7)) then
    !      deallocate(out1D_7)
    !      nullify(cospOUT%calipso_lidarcldtype)
    !   endif
    !   if (allocated(out1D_8)) then
    !      deallocate(out1D_8)
    !      nullify(cospOUT%calipso_cldtype)
    !   endif
    !   if (allocated(out1D_9)) then
    !      deallocate(out1D_9)
    !      nullify(cospOUT%calipso_cldtypetemp)
    !   endif
    !   if (allocated(out1D_10)) then
    !      deallocate(out1D_10)
    !      nullify(cospOUT%calipso_cldtypemeanz)
    !   endif
    !   if (allocated(out1D_12)) then
    !      deallocate(out1D_12)
    !      nullify(cospOUT%calipso_cldtypemeanzse)
    !   endif
    !   if (allocated(out1D_11)) then
    !      deallocate(out1D_11)
    !      nullify(cospOUT%calipso_cldthinemis)
    !   endif 
    endif

    ! PARASOL
    if (Lparasol_column) then
       call parasol_column(parasolIN%Npoints,PARASOL_NREFL,parasolIN%Ncolumns,           &
                            cospgridIN%land(:),parasolPix_refl(:,:,:),                   &
                            cospOUT%parasolGrid_refl(ij:ik,:))
       if (allocated(parasolPix_refl)) deallocate(parasolPix_refl)
    endif

    ! CLOUDSAT
    if (Lcloudsat_column) then
       ! Check to see which outputs are requested. If not requested, use a local dummy array
       if (.not. associated(cospOUT%cloudsat_cfad_ze)) then
          !allocate(out1D_1(Npoints*DBZE_BINS*Nlvgrid))
          !cospOUT%cloudsat_cfad_ze(ij:ik,1:DBZE_BINS,1:Nlvgrid) => out1D_1
          allocate(out3D_1(Npoints,DBZE_BINS,Nlvgrid))  !Jing X.
          cospOUT%cloudsat_cfad_ze => out3D_1
       endif

       if (.not. associated(cospOUT%cloudsat_pia)) then !<-added
          allocate(out1D_2(Npoints)) !<-added
          cospOUT%cloudsat_pia(ij:ik) => out1D_2 !<-added
       endif !<-added 
       if (.not. associated(cospOUT%cloudsat_precip_cover)) then !<-added
          allocate(out1D_3(Npoints*nCloudsatPrecipClass)) !<-added
          cospOUT%cloudsat_precip_cover(ij:ik,1:nCloudsatPrecipClass) => out1D_3 !<-added
       endif !<-added

       ! Call simulator
       call quickbeam_column(cloudsatIN%Npoints,cloudsatIN%Ncolumns,cloudsatIN%Nlevels,  &
            Nlvgrid, DBZE_BINS, 'cloudsat', cloudsatDBZe, cloudsatZe_non,                & !<-added
            cospgridIN%land(:), cospgridIN%surfelev(:), cospgridIN%at(:,cospIN%Nlevels), & !<-added
            cospIN%fracPrecipIce,cospgridIN%hgt_matrix,cospgridIN%hgt_matrix_half,       & !<-added
            cospOUT%cloudsat_cfad_ze(ij:ik,:,:),cospOUT%cloudsat_precip_cover,  & !<-added
            cospOUT%cloudsat_pia ) !<-added
       ! Free up memory  (if necessary)
       !if (allocated(out1D_1)) then
       !   deallocate(out1D_1)
       if (allocated(out3D_1)) then  !Jing X.
          deallocate(out3D_1)
          nullify(cospOUT%cloudsat_cfad_ze)
       endif
       if (allocated(out1D_2)) then !<-added
          deallocate(out1D_2) !<-added
          nullify(cospOUT%cloudsat_pia) !<-added
       endif !<-added
       if (allocated(out1D_3)) then !<-added
          deallocate(out1D_3) !<-added
          nullify(cospOUT%cloudsat_precip_cover) !<-added
       endif !<-added

    endif !<-added

#ifdef OPT_DPLRW
    if (Ldplrw) then
       call quickbeam_dplrw(Npoints, ij, ik, cospIN%Ncolumns, cospIN%Nlevels, cloudsatIN%rcfg, &
                            cospgridIN%hgt_matrix, cospgridIN%hgt_matrix_half, cospgridIN%at, cospgridIN%pfull, &
                            cospgridIN%gwvel, cospgridIN%gcumf, &
                            cospIN%vfall, cospIN%vfsqu, cospIN%zehyd, &
                            cloudsatIN%g_vol, cloudsatIN%kr_vol, &
                            cospOUT%gcumw(ij:ik,:), &
                            !cospOUT%dplrw_sumZ(ij:ik,:,:,:), cospOUT%dplrw_sumT(ij:ik,:,:,:), &
                            cospOUT%dplrw_Z(ij:ik,:,:,:,:), cospOUT%spwid_Z(ij:ik,:,:,:,:), cospOUT%Zef94_Z(ij:ik,:,:,:,:), &
                            cospOUT%dplrw_T(ij:ik,:,:,:,:), cospOUT%spwid_T(ij:ik,:,:,:,:), cospOUT%Zef94_T(ij:ik,:,:,:,:), &
                            cospOUT%ZefVd_2(ij:ik,:,:,:,:), &
                            isccp_boxptop(:,:), isccp_boxtau(:,:) )
    end if
#endif
    ! modified by Doppler simulator
    if (allocated(isccp_boxptop))   deallocate(isccp_boxptop)
    if (allocated(isccp_boxtau))    deallocate(isccp_boxtau)

    ! MODIS
    if (Lmodis_column) then
       if(modisiN%nSunlit > 0) then 
          ! Allocate space for local variables
          allocate(modisCftotal(modisIN%nSunlit), modisCfLiquid(modisIN%nSunlit),        &
                   modisCfIce(modisIN%nSunlit),modisCfHigh(modisIN%nSunlit),             &
                   modisCfMid(modisIN%nSunlit),modisCfLow(modisIN%nSunlit),              &
                   modisMeanTauTotal(modisIN%nSunlit),                                   &
                   modisMeanTauLiquid(modisIN%nSunlit),modisMeanTauIce(modisIN%nSunlit), &
                   modisMeanLogTauTotal(modisIN%nSunlit),                                &       
                   modisMeanLogTauLiquid(modisIN%nSunlit),                               &
                   modisMeanLogTauIce(modisIN%nSunlit),                                  &
                   modisMeanSizeLiquid(modisIN%nSunlit),                                 &
                   modisMeanSizeIce(modisIN%nSunlit),                                    &
                   modisMeanCloudTopPressure(modisIN%nSunlit),                           &
                   modisMeanLiquidWaterPath(modisIN%nSunlit),                            &
                   modisMeanIceWaterPath(modisIN%nSunlit),                               &
                   modisJointHistogram(modisIN%nSunlit,numMODISTauBins,numMODISPresBins),&
                   modisJointHistogramIce(modisIN%nSunlit,numModisTauBins,numMODISReffIceBins),&
                   modisJointHistogramLiq(modisIN%nSunlit,numModisTauBins,numMODISReffLiqBins))
          ! Call simulator
          call modis_column(modisIN%nSunlit, modisIN%Ncolumns,modisRetrievedPhase,       &
                             modisRetrievedCloudTopPressure,modisRetrievedTau,           &
                             modisRetrievedSize, modisCfTotal, modisCfLiquid, modisCfIce,&
                             modisCfHigh, modisCfMid, modisCfLow, modisMeanTauTotal,     &
                             modisMeanTauLiquid, modisMeanTauIce, modisMeanLogTauTotal,  &
                             modisMeanLogTauLiquid, modisMeanLogTauIce,                  &
                             modisMeanSizeLiquid, modisMeanSizeIce,                      &
                             modisMeanCloudTopPressure, modisMeanLiquidWaterPath,        &
                             modisMeanIceWaterPath, modisJointHistogram,                 &
                             modisJointHistogramIce,modisJointHistogramLiq)
          ! Store data (if requested)
          if (associated(cospOUT%modis_Cloud_Fraction_Total_Mean)) then
             cospOUT%modis_Cloud_Fraction_Total_Mean(ij+int(modisIN%sunlit(:))-1)   =    &
                  modisCfTotal
          endif
          if (associated(cospOUT%modis_Cloud_Fraction_Water_Mean)) then
             cospOUT%modis_Cloud_Fraction_Water_Mean(ij+int(modisIN%sunlit(:))-1)   =    &
                  modisCfLiquid
          endif
          if (associated(cospOUT%modis_Cloud_Fraction_Ice_Mean)) then
             cospOUT%modis_Cloud_Fraction_Ice_Mean(ij+int(modisIN%sunlit(:))-1)     =    &
                  modisCfIce
          endif
          if (associated(cospOUT%modis_Cloud_Fraction_High_Mean)) then
             cospOUT%modis_Cloud_Fraction_High_Mean(ij+int(modisIN%sunlit(:))-1)    =    &
                  modisCfHigh
          endif
          if (associated(cospOUT%modis_Cloud_Fraction_Mid_Mean)) then
             cospOUT%modis_Cloud_Fraction_Mid_Mean(ij+int(modisIN%sunlit(:))-1)     =    &
                  modisCfMid
          endif
          if (associated(cospOUT%modis_Cloud_Fraction_Low_Mean)) then
             cospOUT%modis_Cloud_Fraction_Low_Mean(ij+int(modisIN%sunlit(:))-1)     =    &
                  modisCfLow
          endif
          if (associated(cospOUT%modis_Optical_Thickness_Total_Mean)) then
             cospOUT%modis_Optical_Thickness_Total_Mean(ij+int(modisIN%sunlit(:))-1) =   &
                  modisMeanTauTotal
          endif
          if (associated(cospOUT%modis_Optical_Thickness_Water_Mean)) then
             cospOUT%modis_Optical_Thickness_Water_Mean(ij+int(modisIN%sunlit(:))-1) =   &
                  modisMeanTauLiquid
          endif
          if (associated(cospOUT%modis_Optical_Thickness_Ice_Mean)) then
             cospOUT%modis_Optical_Thickness_Ice_Mean(ij+int(modisIN%sunlit(:))-1)  =    &
                  modisMeanTauIce
          endif
          if (associated(cospOUT%modis_Optical_Thickness_Total_LogMean)) then
             cospOUT%modis_Optical_Thickness_Total_LogMean(ij+int(modisIN%sunlit(:))-1)= &
                  modisMeanLogTauTotal
          endif
          if (associated(cospOUT%modis_Optical_Thickness_Water_LogMean)) then
             cospOUT%modis_Optical_Thickness_Water_LogMean(ij+int(modisIN%sunlit(:))-1) = &
                  modisMeanLogTauLiquid
          endif
          if (associated(cospOUT%modis_Optical_Thickness_Ice_LogMean)) then
             cospOUT%modis_Optical_Thickness_Ice_LogMean(ij+int(modisIN%sunlit(:))-1) =  &
                  modisMeanLogTauIce
          endif        
          if (associated(cospOUT%modis_Cloud_Particle_Size_Water_Mean)) then 
             cospOUT%modis_Cloud_Particle_Size_Water_Mean(ij+int(modisIN%sunlit(:))-1) = &
                  modisMeanSizeLiquid
          endif        
          if (associated(cospOUT%modis_Cloud_Particle_Size_Ice_Mean)) then 
             cospOUT%modis_Cloud_Particle_Size_Ice_Mean(ij+int(modisIN%sunlit(:))-1) =   &
                  modisMeanSizeIce
          endif        
          if (associated(cospOUT%modis_Cloud_Top_Pressure_Total_Mean)) then
             cospOUT%modis_Cloud_Top_Pressure_Total_Mean(ij+int(modisIN%sunlit(:))-1) =  &
                  modisMeanCloudTopPressure
          endif        
          if (associated(cospOUT%modis_Liquid_Water_Path_Mean)) then 
             cospOUT%modis_Liquid_Water_Path_Mean(ij+int(modisIN%sunlit(:))-1)      =    &
                  modisMeanLiquidWaterPath
          endif        
          if (associated(cospOUT%modis_Ice_Water_Path_Mean)) then
              cospOUT%modis_Ice_Water_Path_Mean(ij+int(modisIN%sunlit(:))-1)         =   &
                  modisMeanIceWaterPath
          endif        
          if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure)) then
             cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure(ij+            &
                  int(modisIN%sunlit(:))-1, 1:numModisTauBins, :) = modisJointHistogram(:, :, :)           
             ! Reorder pressure bins in joint histogram to go from surface to TOA 
             cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure(ij:ik,:,:) = &
                  cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure(ij:ik,:,numMODISPresBins:1:-1)
          endif
          if (associated(cospOUT%modis_Optical_Thickness_vs_ReffIce)) then
             cospOUT%modis_Optical_Thickness_vs_ReffIce(ij+int(modisIN%sunlit(:))-1, 1:numMODISTauBins,:) = &
                modisJointHistogramIce(:,:,:)
          endif
          if (associated(cospOUT%modis_Optical_Thickness_vs_ReffLiq)) then
             cospOUT%modis_Optical_Thickness_vs_ReffLiq(ij+int(modisIN%sunlit(:))-1, 1:numMODISTauBins,:) = &
                modisJointHistogramLiq(:,:,:)
          endif
                    
          if(modisIN%nSunlit < modisIN%Npoints) then 
             ! Where it's night and we haven't done the retrievals the values are undefined
             if (associated(cospOUT%modis_Cloud_Fraction_Total_Mean))                    &
                cospOUT%modis_Cloud_Fraction_Total_Mean(ij+int(modisIN%notSunlit(:))-1) = R_UNDEF
             if (associated(cospOUT%modis_Cloud_Fraction_Water_Mean))                    &
                cospOUT%modis_Cloud_Fraction_Water_Mean(ij+int(modisIN%notSunlit(:))-1) = R_UNDEF
             if (associated(cospOUT%modis_Cloud_Fraction_Ice_Mean))                      &
                cospOUT%modis_Cloud_Fraction_Ice_Mean(ij+int(modisIN%notSunlit(:))-1) = R_UNDEF
             if (associated(cospOUT%modis_Cloud_Fraction_High_Mean))                     &
                cospOUT%modis_Cloud_Fraction_High_Mean(ij+int(modisIN%notSunlit(:))-1) = R_UNDEF
             if (associated(cospOUT%modis_Cloud_Fraction_Mid_Mean))                      &
                cospOUT%modis_Cloud_Fraction_Mid_Mean(ij+int(modisIN%notSunlit(:))-1) = R_UNDEF
             if (associated(cospOUT%modis_Cloud_Fraction_Low_Mean))                      &
                cospOUT%modis_Cloud_Fraction_Low_Mean(ij+int(modisIN%notSunlit(:))-1) = R_UNDEF
             if (associated(cospOUT%modis_Optical_Thickness_Total_Mean))                 &
                cospOUT%modis_Optical_Thickness_Total_Mean(ij+int(modisIN%notSunlit(:))-1) = R_UNDEF
             if (associated(cospOUT%modis_Optical_Thickness_Water_Mean))                 &
                cospOUT%modis_Optical_Thickness_Water_Mean(ij+int(modisIN%notSunlit(:))-1) = R_UNDEF
             if (associated(cospOUT%modis_Optical_Thickness_Ice_Mean))                   &
                cospOUT%modis_Optical_Thickness_Ice_Mean(ij+int(modisIN%notSunlit(:))-1) = R_UNDEF
             if (associated(cospOUT%modis_Optical_Thickness_Total_LogMean))              &
                cospOUT%modis_Optical_Thickness_Total_LogMean(ij+int(modisIN%notSunlit(:))-1) = R_UNDEF
             if (associated(cospOUT%modis_Optical_Thickness_Water_LogMean))              &
                cospOUT%modis_Optical_Thickness_Water_LogMean(ij+int(modisIN%notSunlit(:))-1) = R_UNDEF
             if (associated(cospOUT%modis_Optical_Thickness_Ice_LogMean))                &
                cospOUT%modis_Optical_Thickness_Ice_LogMean(ij+int(modisIN%notSunlit(:))-1) = R_UNDEF
             if (associated(cospOUT%modis_Cloud_Particle_Size_Water_Mean))               &
                cospOUT%modis_Cloud_Particle_Size_Water_Mean(ij+int(modisIN%notSunlit(:))-1) = R_UNDEF
             if (associated(cospOUT%modis_Cloud_Particle_Size_Ice_Mean))                 &
                cospOUT%modis_Cloud_Particle_Size_Ice_Mean(ij+int(modisIN%notSunlit(:))-1) = R_UNDEF
             if (associated(cospOUT%modis_Cloud_Top_Pressure_Total_Mean))                &
                cospOUT%modis_Cloud_Top_Pressure_Total_Mean(ij+int(modisIN%notSunlit(:))-1) = R_UNDEF
             if (associated(cospOUT%modis_Liquid_Water_Path_Mean))                       &
                cospOUT%modis_Liquid_Water_Path_Mean(ij+int(modisIN%notSunlit(:))-1) = R_UNDEF
             if (associated(cospOUT%modis_Ice_Water_Path_Mean))                          &
                cospOUT%modis_Ice_Water_Path_Mean(ij+int(modisIN%notSunlit(:))-1) = R_UNDEF
             if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure))      &
                cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure(ij+int(modisIN%notSunlit(:))-1, :, :) = R_UNDEF
          end if
       else
          ! It's nightime everywhere - everything is undefined
          if (associated(cospOUT%modis_Cloud_Fraction_Total_Mean))                       &
             cospOUT%modis_Cloud_Fraction_Total_Mean(ij:ik) = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_Water_Mean))                       &
             cospOUT%modis_Cloud_Fraction_Water_Mean(ij:ik) = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_Ice_Mean))                         &
             cospOUT%modis_Cloud_Fraction_Ice_Mean(ij:ik) = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_High_Mean))                        &
             cospOUT%modis_Cloud_Fraction_High_Mean(ij:ik) = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_Mid_Mean))                         &
             cospOUT%modis_Cloud_Fraction_Mid_Mean(ij:ik) = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_Low_Mean))                         &
             cospOUT%modis_Cloud_Fraction_Low_Mean(ij:ik) = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Total_Mean))                    &
             cospOUT%modis_Optical_Thickness_Total_Mean(ij:ik) = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Water_Mean))                    &
             cospOUT%modis_Optical_Thickness_Water_Mean(ij:ik) = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Ice_Mean))                      &
             cospOUT%modis_Optical_Thickness_Ice_Mean(ij:ik) = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Total_LogMean))                 &          
             cospOUT%modis_Optical_Thickness_Total_LogMean(ij:ik) = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Water_LogMean))                 &          
             cospOUT%modis_Optical_Thickness_Water_LogMean(ij:ik) = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Ice_LogMean))                   &          
             cospOUT%modis_Optical_Thickness_Ice_LogMean(ij:ik) = R_UNDEF  
          if (associated(cospOUT%modis_Cloud_Particle_Size_Water_Mean))                  &
             cospOUT%modis_Cloud_Particle_Size_Water_Mean(ij:ik) = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Particle_Size_Ice_Mean))                    &
              cospOUT%modis_Cloud_Particle_Size_Ice_Mean(ij:ik) = R_UNDEF 
          if (associated(cospOUT%modis_Cloud_Top_Pressure_Total_Mean))                   &
             cospOUT%modis_Cloud_Top_Pressure_Total_Mean(ij:ik) = R_UNDEF  
          if (associated(cospOUT%modis_Liquid_Water_Path_Mean))                          &
             cospOUT%modis_Liquid_Water_Path_Mean(ij:ik) = R_UNDEF
          if (associated(cospOUT%modis_Ice_Water_Path_Mean))                             &
             cospOUT%modis_Ice_Water_Path_Mean(ij:ik) = R_UNDEF 
          if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure))         & 
             cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure(ij:ik, :, :) = R_UNDEF
       endif
       ! Free up memory (if necessary)
       if (allocated(modisRetrievedTau))               deallocate(modisRetrievedTau)
       if (allocated(modisRetrievedSize))              deallocate(modisRetrievedSize)
       if (allocated(modisRetrievedPhase))             deallocate(modisRetrievedPhase)
       if (allocated(modisRetrievedCloudTopPressure))  deallocate(modisRetrievedCloudTopPressure)
       if (allocated(modisCftotal))                    deallocate(modisCftotal)
       if (allocated(modisCfLiquid))                   deallocate(modisCfLiquid)
       if (allocated(modisCfIce))                      deallocate(modisCfIce)
       if (allocated(modisCfHigh))                     deallocate(modisCfHigh)
       if (allocated(modisCfMid))                      deallocate(modisCfMid)
       if (allocated(modisCfLow))                      deallocate(modisCfLow)
       if (allocated(modisMeanTauTotal))               deallocate(modisMeanTauTotal)
       if (allocated(modisMeanTauLiquid))              deallocate(modisMeanTauLiquid)
       if (allocated(modisMeanTauIce))                 deallocate(modisMeanTauIce)
       if (allocated(modisMeanLogTauTotal))            deallocate(modisMeanLogTauTotal)
       if (allocated(modisMeanLogTauLiquid))           deallocate(modisMeanLogTauLiquid)
       if (allocated(modisMeanLogTauIce))              deallocate(modisMeanLogTauIce)
       if (allocated(modisMeanSizeLiquid))             deallocate(modisMeanSizeLiquid)
       if (allocated(modisMeanSizeIce))                deallocate(modisMeanSizeIce)
       if (allocated(modisMeanCloudTopPressure))       deallocate(modisMeanCloudTopPressure)
       if (allocated(modisMeanLiquidWaterPath))        deallocate(modisMeanLiquidWaterPath)
       if (allocated(modisMeanIceWaterPath))           deallocate(modisMeanIceWaterPath)
       if (allocated(modisJointHistogram))             deallocate(modisJointHistogram)
       if (allocated(modisJointHistogramIce))          deallocate(modisJointHistogramIce)
       if (allocated(modisJointHistogramLiq))          deallocate(modisJointHistogramLiq)
       if (allocated(isccp_boxttop))                   deallocate(isccp_boxttop)
       if (allocated(isccp_boxptop))                   deallocate(isccp_boxptop)
       if (allocated(isccp_boxtau))                    deallocate(isccp_boxtau)
       if (allocated(isccp_meantbclr))                 deallocate(isccp_meantbclr)
       if (allocated(isccpLEVMATCH))                   deallocate(isccpLEVMATCH)
    endif

    ! RTTOV
    if (lrttov_column) then
       call rttov_column(rttovIN%nPoints,rttovIN%nLevels,rttovIN%nSubCols,rttovIN%q,    &
                         rttovIN%p,rttovIN%t,rttovIN%o3,rttovIN%ph,rttovIN%h_surf,      &
                         rttovIN%u_surf,rttovIN%v_surf,rttovIN%p_surf,rttovIN%t_skin,   &
                         rttovIN%t2m,rttovIN%q2m,rttovIN%lsmask,rttovIN%longitude,      &
                         rttovIN%latitude,rttovIN%seaice,rttovIN%co2,rttovIN%ch4,       &
                         rttovIN%n2o,rttovIN%co,rttovIN%zenang,lrttov_cleanUp,          &
                         cospOUT%rttov_tbs(ij:ik,:),cosp_simulator(nError+1),           &
                         ! Optional arguments for surface emissivity calculation
                         month=rttovIN%month)
                         ! Optional arguments to rttov for all-sky calculation
                         ! rttovIN%month, rttovIN%tca,rttovIN%cldIce,rttovIN%cldLiq,     &
                         ! rttovIN%fl_rain,rttovIN%fl_snow)
    endif

    ! SCOPS fractions ! added by YN
    if (Lscops) then
       allocate(mmar(cospIN%ncolumns)) ; allocate(mmlg(cospIN%ncolumns))
       mmar(:) = 1._wp
       do k=1,cospIN%nlevels
          do i=1,npoints
             mmlg = (cospIN%frac_out(i,:,k) .eq. 1)
             cospOUT%frac_ls(ij+i-1,k) = sum(mmar,mask=mmlg)/cospIN%ncolumns

             mmlg = (cospIN%frac_out(i,:,k) .eq. 2)
             cospOUT%frac_cv(ij+i-1,k) = sum(mmar,mask=mmlg)/cospIN%ncolumns

             mmlg = ( (cospIN%frac_prec(i,:,k) .eq. 1) .or. (cospIN%frac_prec(i,:,k) .eq. 3) )
             cospOUT%prec_ls(ij+i-1,k) = sum(mmar,mask=mmlg)/cospIN%ncolumns

             mmlg = ( (cospIN%frac_prec(i,:,k) .eq. 2) .or. (cospIN%frac_prec(i,:,k) .eq. 3) )
             cospOUT%prec_cv(ij+i-1,k) = sum(mmar,mask=mmlg)/cospIN%ncolumns
          end do
       end do
       deallocate(mmar,mmlg)
       
    end if

!    write( jfpar,* ) ' ### test init by T.Michibata. (test No.6)'
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 6) Compute multi-instrument products
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! CLOUDSAT/CALIPSO products
    if (Lradar_lidar_tcc .or. Llidar_only_freq_cloud) then ! .or. Lcloudsat_tcc .or. Lcloudsat_tcc2) then
      
       if (use_vgrid) then
          allocate(lidar_only_freq_cloud(cloudsatIN%Npoints,Nlvgrid),                    &
               radar_lidar_tcc(cloudsatIN%Npoints)) !, cloudsat_tcc(cloudsatIN%Npoints),    &
               !cloudsat_tcc2(cloudsatIN%Npoints))
          allocate(betamol_in(cloudsatIN%Npoints,1,cloudsatIN%Nlevels),                  &
                   betamolFlip(cloudsatIN%Npoints,1,Nlvgrid),                            &
                   pnormFlip(cloudsatIN%Npoints,cloudsatIN%Ncolumns,Nlvgrid),            &
                   Ze_totFlip(cloudsatIN%Npoints,cloudsatIN%Ncolumns,Nlvgrid))

          betamol_in(:,1,:) = calipso_beta_mol(:,cloudsatIN%Nlevels:1:-1)
          call cosp_change_vertical_grid(cloudsatIN%Npoints,1,cloudsatIN%Nlevels,        &
               cospgridIN%hgt_matrix(:,cloudsatIN%Nlevels:1:-1),                         &
               cospgridIN%hgt_matrix_half(:,cloudsatIN%Nlevels:1:-1),betamol_in,         &
               Nlvgrid,vgrid_zl(Nlvgrid:1:-1),vgrid_zu(Nlvgrid:1:-1),                    &
               betamolFlip(:,1,Nlvgrid:1:-1))
    
          call cosp_change_vertical_grid(cloudsatIN%Npoints,cloudsatIN%Ncolumns,         &
               cloudsatIN%Nlevels,cospgridIN%hgt_matrix(:,cloudsatIN%Nlevels:1:-1),      &
               cospgridIN%hgt_matrix_half(:,cloudsatIN%Nlevels:1:-1),                    &
               calipso_beta_tot(:,:,cloudsatIN%Nlevels:1:-1),Nlvgrid,                    &
               vgrid_zl(Nlvgrid:1:-1),vgrid_zu(Nlvgrid:1:-1),pnormFlip(:,:,Nlvgrid:1:-1))
          
          call cosp_change_vertical_grid(cloudsatIN%Npoints,cloudsatIN%Ncolumns,         &
               cloudsatIN%Nlevels,cospgridIN%hgt_matrix(:,cloudsatIN%Nlevels:1:-1),      &
               cospgridIN%hgt_matrix_half(:,cloudsatIN%Nlevels:1:-1),                    &
               cloudsatDBZe(:,:,cloudsatIN%Nlevels:1:-1),Nlvgrid,vgrid_zl(Nlvgrid:1:-1), &
               vgrid_zu(Nlvgrid:1:-1),Ze_totFlip(:,:,Nlvgrid:1:-1),log_units=.true.)    

          call cosp_lidar_only_cloud(cloudsatIN%Npoints,cloudsatIN%Ncolumns,             &
                                     Nlvgrid,pnormFlip,betamolFlip,Ze_totFlip,           &
                                     lidar_only_freq_cloud,radar_lidar_tcc)  !,              &
                                     !cloudsat_tcc, cloudsat_tcc2)
          
          deallocate(betamol_in,betamolFlip,pnormFlip,ze_totFlip)
       else
          allocate(lidar_only_freq_cloud(cloudsatIN%Npoints,cloudsatIN%Nlevels),         &
               radar_lidar_tcc(cloudsatIN%Npoints)) !, cloudsat_tcc(cloudsatIN%Npoints),    &
               !cloudsat_tcc2(cloudsatIN%Npoints))
          call cosp_lidar_only_cloud(cloudsatIN%Npoints,cloudsatIN%Ncolumns,             &
               cospIN%Nlevels,calipso_beta_tot(:,:,cloudsatIN%Nlevels:1:-1),             &
               calipso_beta_mol(:,cloudsatIN%Nlevels:1:-1),                              &
               cloudsatDBZe(:,:,cloudsatIN%Nlevels:1:-1),lidar_only_freq_cloud,          &
               radar_lidar_tcc) !,cloudsat_tcc, cloudsat_tcc2)
       endif
       
       ! Store, when necessary
       if (associated(cospOUT%lidar_only_freq_cloud)) then
          cospOUT%lidar_only_freq_cloud(ij:ik,:) = lidar_only_freq_cloud
       endif
       if (associated(cospOUT%radar_lidar_tcc)) then
          cospOUT%radar_lidar_tcc(ij:ik) = radar_lidar_tcc
       endif
       !if (associated(cospOUT%cloudsat_tcc)) then
       !   cospOUT%cloudsat_tcc(ij:ik) = cloudsat_tcc
       !endif
       !if (associated(cospOUT%cloudsat_tcc2)) then
       !   cospOUT%cloudsat_tcc2(ij:ik) = cloudsat_tcc2
       !endif
    endif

    ! CLOUDSAT/MODIS products (CFODDs and Occurrence Frequency of Warm Clouds)
    ! CloudSat/MODIS joint products (CFODDs and Occurrence Frequency of Warm Clouds)
    if (Lcloudsat_modis_wr) then
       allocate( cfodd_ntotal(cloudsatIN%Npoints, CFODD_NDBZE, CFODD_NICOD, CFODD_NCLASS) )
       allocate( wr_occfreq_ntotal(cloudsatIN%Npoints, WR_NREGIME) )

       if ( use_vgrid ) then
          !! interporation for fixed vertical grid:
          allocate( zlev(cloudsatIN%Npoints,Nlvgrid),                         &
                    delz(cloudsatIN%Npoints,Nlvgrid),                         &
                    t_in(cloudsatIN%Npoints,1,cloudsatIN%Nlevels),            &
                    tempI(cloudsatIN%Npoints,1,Nlvgrid),                      &
                    Ze_totI(cloudsatIN%Npoints,cloudsatIN%Ncolumns,Nlvgrid),  &
                    frac_outI(cloudsatIN%Npoints,cloudsatIN%Ncolumns,Nlvgrid) )
          do k = 1, Nlvgrid
             zlev(:,k) = vgrid_zu(k)
             delz(:,k) = dz(k)
          enddo
          t_in(:,1,:) = cospgridIN%at(:,:)
          call cosp_change_vertical_grid (                                    &
               cloudsatIN%Npoints, 1, cloudsatIN%Nlevels,                     &
               cospgridIN%hgt_matrix(:,cloudsatIN%Nlevels:1:-1),              &
               cospgridIN%hgt_matrix_half(:,cloudsatIN%Nlevels:1:-1),         &
               t_in(:,:,cloudsatIN%Nlevels:1:-1), Nlvgrid,                    &
               vgrid_zl(Nlvgrid:1:-1), vgrid_zu(Nlvgrid:1:-1),                &
               tempI(:,:,Nlvgrid:1:-1)                                        )
          call cosp_change_vertical_grid (                                    &
               cloudsatIN%Npoints, cloudsatIN%Ncolumns, cloudsatIN%Nlevels,   &
               cospgridIN%hgt_matrix(:,cloudsatIN%Nlevels:1:-1),              &
               cospgridIN%hgt_matrix_half(:,cloudsatIN%Nlevels:1:-1),         &
               cloudsatDBZe(:,:,cloudsatIN%Nlevels:1:-1), Nlvgrid,            &
               vgrid_zl(Nlvgrid:1:-1), vgrid_zu(Nlvgrid:1:-1),                &
               Ze_totI(:,:,Nlvgrid:1:-1), log_units=.true.                    )
          call cosp_change_vertical_grid (                                    &
               cloudsatIN%Npoints, cloudsatIN%Ncolumns, cloudsatIN%Nlevels,   &
               cospgridIN%hgt_matrix(:,cloudsatIN%Nlevels:1:-1),              &
               cospgridIN%hgt_matrix_half(:,cloudsatIN%Nlevels:1:-1),         &
               cospIN%frac_out(:,:,cloudsatIN%Nlevels:1:-1), Nlvgrid,         &
               vgrid_zl(Nlvgrid:1:-1), vgrid_zu(Nlvgrid:1:-1),                &
               frac_outI(:,:,Nlvgrid:1:-1)                                    )
          call cosp_diag_warmrain(                                            &
               cloudsatIN%Npoints, cloudsatIN%Ncolumns, Nlvgrid,              & !! in
               tempI, zlev, delz,                                             & !! in
               cospOUT%modis_Liquid_Water_Path_Mean,                          & !! in
               cospOUT%modis_Optical_Thickness_Water_Mean,                    & !! in
               cospOUT%modis_Cloud_Particle_Size_Water_Mean,                  & !! in
               cospOUT%modis_Cloud_Fraction_Water_Mean,                       & !! in
               cospOUT%modis_Ice_Water_Path_Mean,                             & !! in
               cospOUT%modis_Optical_Thickness_Ice_Mean,                      & !! in
               cospOUT%modis_Cloud_Particle_Size_Ice_Mean,                    & !! in
               cospOUT%modis_Cloud_Fraction_Ice_Mean,                         & !! in
               frac_outI,                                                     & !! in
               Ze_totI,                                                       & !! in
               cfodd_ntotal, wr_occfreq_ntotal                                ) !! inout
          deallocate( zlev, delz, t_in, tempI, frac_outI, Ze_totI )
       else  ! do not use vgrid interporation ---------------------------------------!
          !! original model grid
          allocate( delz(cloudsatIN%Npoints,cospIN%Nlevels) )
          do k = 1, cospIN%Nlevels-1
             delz(:,k) = cospgridIN%hgt_matrix_half(:,k+1) &
                         - cospgridIN%hgt_matrix_half(:,k)
          enddo
          delz(:,cospIN%Nlevels) = 2.0*( cospgridIN%hgt_matrix(:,cospIN%Nlevels) &
                                  - cospgridIN%hgt_matrix_half(:,cospIN%Nlevels) )
          call cosp_diag_warmrain(                                            &
               cloudsatIN%Npoints, cloudsatIN%Ncolumns, cospIN%Nlevels,       & !! in
               cospgridIN%at, cospgridIN%hgt_matrix, delz,                    & !! in
               cospOUT%modis_Liquid_Water_Path_Mean,                          & !! in
               cospOUT%modis_Optical_Thickness_Water_Mean,                    & !! in
               cospOUT%modis_Cloud_Particle_Size_Water_Mean,                  & !! in
               cospOUT%modis_Cloud_Fraction_Water_Mean,                       & !! in
               cospOUT%modis_Ice_Water_Path_Mean,                             & !! in
               cospOUT%modis_Optical_Thickness_Ice_Mean,                      & !! in
               cospOUT%modis_Cloud_Particle_Size_Ice_Mean,                    & !! in
               cospOUT%modis_Cloud_Fraction_Ice_Mean,                         & !! in
               cospIN%frac_out,                                               & !! in
               cloudsatDBZe,                                                  & !! in
               cfodd_ntotal, wr_occfreq_ntotal                                ) !! inout
          deallocate( delz )
       endif  !! use_vgrid or not

       ! Store, when necessary
       if ( associated(cospOUT%cfodd_ntotal) ) then
          cospOUT%cfodd_ntotal(ij:ik,:,:,:) = cfodd_ntotal
       endif
       if ( associated(cospOUT%wr_occfreq_ntotal) ) then
          cospOUT%wr_occfreq_ntotal(ij:ik,:) = wr_occfreq_ntotal
       endif
    endif

!    write( jfpar,* ) ' ### test init by T.Michibata. (test No.7)'
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 7) Cleanup
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (Lisccp_subcolumn .or. Lmodis_subcolumn) then
       nullify(isccpIN%Npoints,isccpIN%Ncolumns,isccpIN%Nlevels,isccpIN%emsfc_lw,        &
               isccpIN%skt,isccpIN%qv,isccpIN%at,isccpIN%frac_out,isccpIN%dtau,          &
               isccpIN%dem,isccpIN%phalf,isccpIN%sunlit,isccpIN%pfull)
    endif
    
    if (Lmisr_subcolumn) then
       nullify(misrIN%Npoints,misrIN%Ncolumns,misrIN%Nlevels,misrIN%dtau,misrIN%sunlit,  &
               misrIN%zfull,misrIN%at)
    endif

    if (Lcalipso_subcolumn) then
       nullify(calipsoIN%Npoints,calipsoIN%Ncolumns,calipsoIN%Nlevels,calipsoIN%beta_mol,&
               calipsoIN%betatot,calipsoIN%betatot_liq,calipsoIN%betatot_ice,            &
               calipsoIN%tau_mol,calipsoIN%tautot,calipsoIN%tautot_liq,calipsoIN%tautot_ice)
    endif
    
    if (Lparasol_subcolumn) then
       nullify(parasolIN%Npoints,parasolIN%Nlevels,parasolIN%Ncolumns,parasolIN%Nrefl,   &
            parasolIN%tautot_S_liq,parasolIN%tautot_S_ice)
    endif
    
    if (Lcloudsat_subcolumn) then
       nullify(cloudsatIN%Npoints,cloudsatIN%Nlevels,cloudsatIN%Ncolumns,cloudsatIN%rcfg,&
               cloudsatIN%kr_vol,cloudsatIN%g_vol,cloudsatIN%z_vol,cloudsatIN%hgt_matrix)
    endif

    if (Lmodis_subcolumn) then
       nullify(modisIN%Npoints,modisIN%Ncolumns,modisIN%Nlevels,modisIN%tau,modisIN%g,   &
               modisIN%liqFrac,modisIN%w0)
       if (allocated(modisIN%sunlit))    deallocate(modisIN%sunlit)
       if (allocated(modisIN%notSunlit)) deallocate(modisIN%notSunlit)
       if (allocated(modisIN%pres))      deallocate(modisIN%pres)
    endif

    if (allocated(calipso_beta_tot))      deallocate(calipso_beta_tot)
    if (allocated(calipso_beta_mol))      deallocate(calipso_beta_mol)
    if (allocated(calipso_betaperp_tot))  deallocate(calipso_betaperp_tot)
    if (allocated(cloudsatDBZe))          deallocate(cloudsatDBZe)
    if (allocated(lidar_only_freq_cloud)) deallocate(lidar_only_freq_cloud)
    if (allocated(radar_lidar_tcc))       deallocate(radar_lidar_tcc)
    !if (allocated(cloudsat_tcc))          deallocate(cloudsat_tcc)
    !if (allocated(cloudsat_tcc2))         deallocate(cloudsat_tcc2)
    !! IDiD (T.Michibata, 2018.10.29)
    if( allocated(cfodd_ntotal)  )        deallocate(cfodd_ntotal)
    if( allocated(wr_occfreq_ntotal) )    deallocate(wr_occfreq_ntotal)

    call flush(jfpar)
  end function COSP_SIMULATOR

  ! ######################################################################################
  ! SUBROUTINE cosp_init
  ! ######################################################################################
  !SUBROUTINE COSP_INIT(Lisccp,Lmodis,Lmisr,Lcloudsat,Lcalipso,Lparasol,Lrttov,           &
  !                     Npoints,Nlevels,cloudsat_radar_freq,cloudsat_k2,                  &
  !                     cloudsat_use_gas_abs,cloudsat_do_ray,isccp_top_height,            &
  !                     isccp_top_height_direction,surface_radar,rcfg,rttov_Nchannels,    &
  !                     rttov_Channels,rttov_platform,rttov_satellite,rttov_instrument,   &
  !                     lusevgrid,luseCSATvgrid,Nvgrid,cloudsat_micro_scheme,cospOUT)
  SUBROUTINE COSP_INIT(Lisccp,Lmodis,Lmisr,Lcloudsat,Lcalipso,Lparasol,Lrttov,           &
                       Npoints,Nlevels,cloudsat_radar_freq,cloudsat_k2,                  &
                       cloudsat_use_gas_abs,cloudsat_do_ray,isccp_top_height,            &
                       isccp_top_height_direction,surface_radar,rcfg,rttov_Nchannels,    &
                       rttov_Channels,rttov_platform,rttov_satellite,rttov_instrument,   &
                       lusevgrid,luseCSATvgrid,Nvgrid,                                   &
!#ifdef OPT_DPLRW
!                       Ldplrw, &
!#endif
                       cloudsat_micro_scheme,cospOUT)
    
    ! INPUTS
    logical,intent(in) :: Lisccp,Lmodis,Lmisr,Lcloudsat,Lcalipso,Lparasol,Lrttov
    integer,intent(in)  :: &
         cloudsat_use_gas_abs,       & ! 
         cloudsat_do_ray,            & !
         isccp_top_height,           & !
         isccp_top_height_direction, & !
         Npoints,                    & !
         Nlevels,                    & !
         Nvgrid,                     & ! Number of levels for new L3 grid
         surface_radar,              & ! 
         rttov_Nchannels,            & ! Number of RTTOV channels
         rttov_platform,             & ! RTTOV platform
         rttov_satellite,            & ! RTTOV satellite
         rttov_instrument              ! RTTOV instrument
    integer,intent(in),dimension(RTTOV_MAX_CHANNELS) :: &
         rttov_channels                ! RTTOV channels    
    real(wp),intent(in) :: &
         cloudsat_radar_freq,        & !
         cloudsat_k2                   !   
    logical,intent(in) :: &
         lusevgrid,                  & ! Switch to use different vertical grid
         luseCSATvgrid                 ! Switch to use CLOUDSAT grid spacing for new  
                                       ! vertical grid
!#ifdef OPT_DPLRW
!    logical,intent(in) :: Ldplrw
!q#endif
    character(len=64),intent(in) :: &
       cloudsat_micro_scheme           ! Microphysical scheme used by CLOUDSAT
    type(cosp_outputs),intent(inout) :: cospOUT
    
    ! OUTPUTS
    type(radar_cfg) :: rcfg
 
    ! Local variables
    integer  :: i
    real(wp) :: zstep
    
    integer :: jfpar 

    call getjfp(jfpar)
    !write( jfpar,* ) ' ### test init by Jing X. (test No.1)'
!    write( jfpar,* ) ' ### test init by T.Michibata. (test INIT)'

    ! Initialize MODIS optical-depth bin boundaries for joint-histogram. (defined in cosp_config.F90)
    if (.not. allocated(modis_histTau)) then
       allocate(modis_histTau(ntau+1),modis_histTauEdges(2,ntau),modis_histTauCenters(ntau))
       numMODISTauBins      = ntau
       modis_histTau        = tau_binBounds
       modis_histTauEdges   = tau_binEdges
       modis_histTauCenters = tau_binCenters
    endif

    !write( jfpar,* ) ' ### test init by Jing X. (test No.1.2)'
    
    ! Set up vertical grid used by CALIPSO and CLOUDSAT L3
    use_vgrid = lusevgrid
    
    if (use_vgrid) then
      Nlvgrid  = Nvgrid
       allocate(vgrid_zl(Nlvgrid),vgrid_zu(Nlvgrid),vgrid_z(Nlvgrid),dz(Nlvgrid))
       ! CloudSat grid requested
       if (luseCSATvgrid)       zstep = 480._wp
       ! Other grid requested. Constant vertical spacing with top at 20 km
       if (.not. luseCSATvgrid) zstep = 20000._wp/Nvgrid
       do i=1,Nvgrid
          vgrid_zl(Nlvgrid-i+1) = (i-1)*zstep
          vgrid_zu(Nlvgrid-i+1) = i*zstep
       enddo
       vgrid_z = (vgrid_zl+vgrid_zu)/2._wp
       dz = zstep
    else
      ! cancelled. Jing X.
       Nlvgrid = Nlevels
       allocate(vgrid_zl(Nlvgrid),vgrid_zu(Nlvgrid),vgrid_z(Nlvgrid),dz(Nlvgrid))
    endif

    !write( jfpar,* ) ' ### test init by Jing X. (test No.2)'
    ! Initialize simulators
    if (Lisccp) call cosp_isccp_init(isccp_top_height,isccp_top_height_direction)
    !write( jfpar,* ) ' ### test init by Jing X. (test No.3)'
    if (Lmodis) call cosp_modis_init()
    !write( jfpar,* ) ' ### test init by Jing X. (test No.4)'
    if (Lmisr)  call cosp_misr_init()
    !write( jfpar,* ) ' ### test init by Jing X. (test No.5)'
    !if (Lrttov) call cosp_rttov_init(rttov_Nchannels,rttov_platform,rttov_satellite,     &
    !     rttov_instrument,rttov_channels)
    if (Lrttov) call cosp_rttov_init()
    !write( jfpar,* ) ' ### test init by Jing X. (test No.6)'
    if (Lcloudsat) call cosp_cloudsat_init(cloudsat_radar_freq,cloudsat_k2,              &
         cloudsat_use_gas_abs,cloudsat_do_ray,R_UNDEF,N_HYDRO, surface_radar,            &
         rcfg,cloudsat_micro_scheme)
    !write( jfpar,* ) ' ### test init by Jing X. (test No.7)'
    if (Lcalipso) call cosp_calipso_init()
    !write( jfpar,* ) ' ### test init by Jing X. (test No.8)'
    if (Lparasol) call cosp_parasol_init()
    !write( jfpar,* ) ' ### test init by Jing X. (test No.9)'

!    if ( Lwarmrain ) then
!       nullify(cospOUT%cfodd_ntotal, cospOUT%wr_occfreq_ntotal)  !! IDiD Warm-Rain
!       write( jfpar,* ) ' ### test init IDiD by T. Michibata (test No.10): true'
!    else
!       write( jfpar,* ) ' ### test init IDiD by T. Michibata (test No.10): false'
!    endif

    linitialization = .FALSE.
  END SUBROUTINE COSP_INIT

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_cleanUp
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine cosp_cleanUp()
    deallocate(vgrid_zl,vgrid_zu,vgrid_z,dz)
  end subroutine cosp_cleanUp
   
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_errorCheck
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine cosp_errorCheck(cospgridIN,cospIN,Lisccp_subcolumn,Lisccp_column,Lmisr_subcolumn,Lmisr_column,    &
                             Lmodis_subcolumn,Lmodis_column,Lcloudsat_subcolumn,Lcloudsat_column,Lcalipso_subcolumn,  &
                             Lcalipso_column,Lrttov_subcolumn,Lrttov_column,Lparasol_subcolumn,Lparasol_column,    &
                             !,Lcloudsat_tcc, Lcloudsat_tcc2,
                             Lradar_lidar_tcc,Llidar_only_freq_cloud,Lcloudsat_modis_wr, &
                             cospOUT, errorMessage, nError)
  ! Inputs
  type(cosp_column_inputs),intent(in) :: &
     cospgridIN       ! Model grid inputs to COSP
  type(cosp_optical_inputs),intent(in) :: &
     cospIN           ! Derived (optical) input to COSP

  ! Outputs   
  logical,intent(inout) :: &
      Lisccp_subcolumn,    & ! ISCCP subcolumn simulator on/off switch
      Lisccp_column,       & ! ISCCP column simulator on/off switch
      Lmisr_subcolumn,     & ! MISR subcolumn simulator on/off switch
      Lmisr_column,        & ! MISR column simulator on/off switch
      Lmodis_subcolumn,    & ! MODIS subcolumn simulator on/off switch
      Lmodis_column,       & ! MODIS column simulator on/off switch
      Lcloudsat_subcolumn, & ! CLOUDSAT subcolumn simulator on/off switch
      Lcloudsat_column,    & ! CLOUDSAT column simulator on/off switch
      Lcalipso_subcolumn,  & ! CALIPSO subcolumn simulator on/off switch
      Lcalipso_column,     & ! CALIPSO column simulator on/off switch
      Lparasol_subcolumn,  & ! PARASOL subcolumn simulator on/off switch
      Lparasol_column,     & ! PARASOL column simulator on/off switch
      Lrttov_subcolumn,    & ! RTTOV subcolumn simulator on/off switch
      Lrttov_column,       & ! RTTOV column simulator on/off switch      
      !Lcloudsat_tcc,       & !
      !Lcloudsat_tcc2,      & ! 
      Lradar_lidar_tcc,    & ! On/Off switch for joint Calipso/Cloudsat product
      Llidar_only_freq_cloud, & ! On/Off switch for joint Calipso/Cloudsat product
      Lcloudsat_modis_wr     ! On/Off switch for joint CloudSat/MODIS warm rain product

  type(cosp_outputs),intent(inout) :: &
       cospOUT                ! COSP Outputs
  character(len=256),dimension(100) :: errorMessage
  integer,intent(out) :: nError
  
  ! Local variables
  character(len=100) :: parasolErrorMessage

  nError = 0
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! PART 1: Check input array values for out-of-bounds values. When an out-of-bound value
  !         is encountered, COSP outputs that are dependent on that input are filled with
  !         an undefined value (set in cosp_config.f90) and if necessary, that simulator 
  !         is turned off.
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (any(cospgridIN%sunlit .lt. 0)) then
     nError=nError+1
     errorMessage(nError) = 'ERROR: COSP input variable: cospgridIN%sunlit contains values out of range (0 or 1)'
     Lisccp_subcolumn = .false.
     Lisccp_column    = .false.
     Lmisr_subcolumn  = .false.
     Lmisr_column     = .false.
     Lmodis_subcolumn = .false.
     Lmodis_column    = .false.
     Lcloudsat_modis_wr = .false.
     if (associated(cospOUT%isccp_totalcldarea))  cospOUT%isccp_totalcldarea(:)  = R_UNDEF
     if (associated(cospOUT%isccp_meantb))        cospOUT%isccp_meantb(:)        = R_UNDEF
     if (associated(cospOUT%isccp_meantbclr))     cospOUT%isccp_meantbclr(:)     = R_UNDEF
     if (associated(cospOUT%isccp_meanptop))      cospOUT%isccp_meanptop(:)      = R_UNDEF
     if (associated(cospOUT%isccp_meantaucld))    cospOUT%isccp_meantaucld(:)    = R_UNDEF
     if (associated(cospOUT%isccp_meanalbedocld)) cospOUT%isccp_meanalbedocld(:) = R_UNDEF
     if (associated(cospOUT%isccp_boxtau))        cospOUT%isccp_boxtau(:,:)      = R_UNDEF
     if (associated(cospOUT%isccp_boxptop))       cospOUT%isccp_boxptop(:,:)     = R_UNDEF
     if (associated(cospOUT%isccp_fq))            cospOUT%isccp_fq(:,:,:)        = R_UNDEF 
     if (associated(cospOUT%misr_fq))                   cospOUT%misr_fq(:,:,:)                 = R_UNDEF
     if (associated(cospOUT%misr_dist_model_layertops)) cospOUT%misr_dist_model_layertops(:,:) = R_UNDEF
     if (associated(cospOUT%misr_meanztop))             cospOUT%misr_meanztop(:)               = R_UNDEF
     if (associated(cospOUT%misr_cldarea))              cospOUT%misr_cldarea(:)                = R_UNDEF
     if (associated(cospOUT%modis_Cloud_Fraction_Total_Mean))                          &
          cospOUT%modis_Cloud_Fraction_Total_Mean(:)                   = R_UNDEF
     if (associated(cospOUT%modis_Cloud_Fraction_Water_Mean))                          &
          cospOUT%modis_Cloud_Fraction_Water_Mean(:)                   = R_UNDEF
     if (associated(cospOUT%modis_Cloud_Fraction_Ice_Mean))                            &
          cospOUT%modis_Cloud_Fraction_Ice_Mean(:)                     = R_UNDEF
     if (associated(cospOUT%modis_Cloud_Fraction_High_Mean))                           &
          cospOUT%modis_Cloud_Fraction_High_Mean(:)                    = R_UNDEF
     if (associated(cospOUT%modis_Cloud_Fraction_Mid_Mean))                            &
          cospOUT%modis_Cloud_Fraction_Mid_Mean(:)                     = R_UNDEF
     if (associated(cospOUT%modis_Cloud_Fraction_Low_Mean))                            &
          cospOUT%modis_Cloud_Fraction_Low_Mean(:)                     = R_UNDEF
     if (associated(cospOUT%modis_Optical_Thickness_Total_Mean))                       &
          cospOUT%modis_Optical_Thickness_Total_Mean(:)                = R_UNDEF
     if (associated(cospOUT%modis_Optical_Thickness_Water_Mean))                       &
          cospOUT%modis_Optical_Thickness_Water_Mean(:)                = R_UNDEF
     if (associated(cospOUT%modis_Optical_Thickness_Ice_Mean))                         &
          cospOUT%modis_Optical_Thickness_Ice_Mean(:)                  = R_UNDEF
     if (associated(cospOUT%modis_Optical_Thickness_Total_LogMean))                    &
          cospOUT%modis_Optical_Thickness_Total_LogMean(:)             = R_UNDEF
     if (associated(cospOUT%modis_Optical_Thickness_Water_LogMean))                    &
          cospOUT%modis_Optical_Thickness_Water_LogMean(:)             = R_UNDEF
     if (associated(cospOUT%modis_Optical_Thickness_Ice_LogMean))                      &
          cospOUT%modis_Optical_Thickness_Ice_LogMean(:)               = R_UNDEF
     if (associated(cospOUT%modis_Cloud_Particle_Size_Water_Mean))                     &
          cospOUT%modis_Cloud_Particle_Size_Water_Mean(:)              = R_UNDEF
     if (associated(cospOUT%modis_Cloud_Particle_Size_Ice_Mean))                       &
          cospOUT%modis_Cloud_Particle_Size_Ice_Mean(:)                = R_UNDEF
     if (associated(cospOUT%modis_Cloud_Top_Pressure_Total_Mean))                      &
          cospOUT%modis_Cloud_Top_Pressure_Total_Mean(:)               = R_UNDEF
     if (associated(cospOUT%modis_Liquid_Water_Path_Mean))                             &
          cospOUT%modis_Liquid_Water_Path_Mean(:)                      = R_UNDEF
     if (associated(cospOUT%modis_Ice_Water_Path_Mean))                                &
          cospOUT%modis_Ice_Water_Path_Mean(:)                         = R_UNDEF
     if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure))            &
          cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure(:,:,:) = R_UNDEF
     if (associated(cospOUT%modis_Optical_Thickness_vs_ReffICE))                       &
          cospOUT%modis_Optical_Thickness_vs_ReffICE(:,:,:)            = R_UNDEF
     if (associated(cospOUT%modis_Optical_Thickness_vs_ReffLIQ))                       &
          cospOUT%modis_Optical_Thickness_vs_ReffLIQ(:,:,:)            = R_UNDEF
  endif
  if (any(cospgridIN%at .lt. 0)) then   
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospgridIN%at contains values out of range (at<0), expected units (K)'
       Lisccp_subcolumn = .false.
       Lisccp_column    = .false.
       Lmisr_subcolumn  = .false.
       Lmisr_column     = .false.
       Lrttov_subcolumn = .false.
       Lcalipso_column  = .false.
       Lcloudsat_column = .false.
       Lradar_lidar_tcc = .false.
       Llidar_only_freq_cloud = .false.
      ! Lcloudsat_tcc    = .false.
      ! Lcloudsat_tcc2   = .false.
       Lcloudsat_modis_wr = .false.
       if (associated(cospOUT%rttov_tbs)) cospOUT%rttov_tbs(:,:)         = R_UNDEF       
       if (associated(cospOUT%isccp_totalcldarea))  cospOUT%isccp_totalcldarea(:)  = R_UNDEF
       if (associated(cospOUT%isccp_meantb))        cospOUT%isccp_meantb(:)        = R_UNDEF
       if (associated(cospOUT%isccp_meantbclr))     cospOUT%isccp_meantbclr(:)     = R_UNDEF
       if (associated(cospOUT%isccp_meanptop))      cospOUT%isccp_meanptop(:)      = R_UNDEF
       if (associated(cospOUT%isccp_meantaucld))    cospOUT%isccp_meantaucld(:)    = R_UNDEF
       if (associated(cospOUT%isccp_meanalbedocld)) cospOUT%isccp_meanalbedocld(:) = R_UNDEF
       if (associated(cospOUT%isccp_boxtau))        cospOUT%isccp_boxtau(:,:)      = R_UNDEF
       if (associated(cospOUT%isccp_boxptop))       cospOUT%isccp_boxptop(:,:)     = R_UNDEF
       if (associated(cospOUT%isccp_fq))            cospOUT%isccp_fq(:,:,:)        = R_UNDEF
       if (associated(cospOUT%misr_fq))                   cospOUT%misr_fq(:,:,:)                 = R_UNDEF
       if (associated(cospOUT%misr_dist_model_layertops)) cospOUT%misr_dist_model_layertops(:,:) = R_UNDEF
       if (associated(cospOUT%misr_meanztop))             cospOUT%misr_meanztop(:)               = R_UNDEF
       if (associated(cospOUT%misr_cldarea))              cospOUT%misr_cldarea(:)                = R_UNDEF
       if (associated(cospOUT%calipso_cfad_sr))       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF
       if (associated(cospOUT%calipso_lidarcld))      cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldphase)) cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_cldlayerphase)) cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldtmp))   cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF           
     !  if (associated(cospOUT%calipso_lidarcldtype)) cospOUT%calipso_lidarcldtype(:,:,:)   = R_UNDEF
     !  if (associated(cospOUT%calipso_cldtype)) cospOUT%calipso_cldtype(:,:)               = R_UNDEF
     !  if (associated(cospOUT%calipso_cldtypetemp)) cospOUT%calipso_cldtypetemp(:,:)       = R_UNDEF
     !  if (associated(cospOUT%calipso_cldtypemeanz)) cospOUT%calipso_cldtypemeanz(:,:)     = R_UNDEF
     !  if (associated(cospOUT%calipso_cldtypemeanzse)) cospOUT%calipso_cldtypemeanzse(:,:) = R_UNDEF
     !  if (associated(cospOUT%calipso_cldthinemis)) cospOUT%calipso_cldthinemis(:)         = R_UNDEF 
       if (associated(cospOUT%cloudsat_cfad_ze))      cospOUT%cloudsat_cfad_ze(:,:,:)      = R_UNDEF
       if (associated(cospOUT%lidar_only_freq_cloud)) cospOUT%lidar_only_freq_cloud(:,:)   = R_UNDEF
       if (associated(cospOUT%radar_lidar_tcc))       cospOUT%radar_lidar_tcc(:)           = R_UNDEF 
     !  if (associated(cospOUT%cloudsat_tcc)) cospOUT%cloudsat_tcc(:)   = R_UNDEF
     !  if (associated(cospOUT%cloudsat_tcc2)) cospOUT%cloudsat_tcc2(:) = R_UNDEF
       if (associated(cospOUT%cfodd_ntotal)) cospOUT%cfodd_ntotal(:,:,:,:) = R_UNDEF !<-added
       if (associated(cospOUT%wr_occfreq_ntotal)) cospOUT%wr_occfreq_ntotal(:,:) = R_UNDEF !<-added
    endif
    if (any(cospgridIN%pfull .lt. 0)) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospgridIN%pfull contains values out of range'
       Lisccp_subcolumn = .false.
       Lisccp_column    = .false.     
       Lrttov_subcolumn = .false.
       if (associated(cospOUT%rttov_tbs))           cospOUT%rttov_tbs(:,:)         = R_UNDEF       
       if (associated(cospOUT%isccp_totalcldarea))  cospOUT%isccp_totalcldarea(:)  = R_UNDEF
       if (associated(cospOUT%isccp_meantb))        cospOUT%isccp_meantb(:)        = R_UNDEF
       if (associated(cospOUT%isccp_meantbclr))     cospOUT%isccp_meantbclr(:)     = R_UNDEF
       if (associated(cospOUT%isccp_meanptop))      cospOUT%isccp_meanptop(:)      = R_UNDEF
       if (associated(cospOUT%isccp_meantaucld))    cospOUT%isccp_meantaucld(:)    = R_UNDEF
       if (associated(cospOUT%isccp_meanalbedocld)) cospOUT%isccp_meanalbedocld(:) = R_UNDEF
       if (associated(cospOUT%isccp_boxtau))        cospOUT%isccp_boxtau(:,:)      = R_UNDEF
       if (associated(cospOUT%isccp_boxptop))       cospOUT%isccp_boxptop(:,:)     = R_UNDEF
       if (associated(cospOUT%isccp_fq))            cospOUT%isccp_fq(:,:,:)        = R_UNDEF      
    endif
    if (any(cospgridIN%phalf .lt. 0)) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospgridIN%phalf contains values out of range'
       Lisccp_subcolumn = .false.
       Lisccp_column    = .false.     
       Lrttov_subcolumn = .false.
       Lmodis_subcolumn = .false.
       Lmodis_column    = .false.
       Lcalipso_column  = .false.
       if (associated(cospOUT%rttov_tbs))           cospOUT%rttov_tbs(:,:)         = R_UNDEF       
       if (associated(cospOUT%isccp_totalcldarea))  cospOUT%isccp_totalcldarea(:)  = R_UNDEF
       if (associated(cospOUT%isccp_meantb))        cospOUT%isccp_meantb(:)        = R_UNDEF
       if (associated(cospOUT%isccp_meantbclr))     cospOUT%isccp_meantbclr(:)     = R_UNDEF
       if (associated(cospOUT%isccp_meanptop))      cospOUT%isccp_meanptop(:)      = R_UNDEF
       if (associated(cospOUT%isccp_meantaucld))    cospOUT%isccp_meantaucld(:)    = R_UNDEF
       if (associated(cospOUT%isccp_meanalbedocld)) cospOUT%isccp_meanalbedocld(:) = R_UNDEF
       if (associated(cospOUT%isccp_boxtau))        cospOUT%isccp_boxtau(:,:)      = R_UNDEF
       if (associated(cospOUT%isccp_boxptop))       cospOUT%isccp_boxptop(:,:)     = R_UNDEF
       if (associated(cospOUT%isccp_fq))            cospOUT%isccp_fq(:,:,:)        = R_UNDEF 
       if (associated(cospOUT%modis_Cloud_Fraction_Total_Mean))                          &
            cospOUT%modis_Cloud_Fraction_Total_Mean(:)                   = R_UNDEF
       if (associated(cospOUT%modis_Cloud_Fraction_Water_Mean))                          &
            cospOUT%modis_Cloud_Fraction_Water_Mean(:)                   = R_UNDEF
       if (associated(cospOUT%modis_Cloud_Fraction_Ice_Mean))                            &
            cospOUT%modis_Cloud_Fraction_Ice_Mean(:)                     = R_UNDEF
       if (associated(cospOUT%modis_Cloud_Fraction_High_Mean))                           &
            cospOUT%modis_Cloud_Fraction_High_Mean(:)                    = R_UNDEF
       if (associated(cospOUT%modis_Cloud_Fraction_Mid_Mean))                            &
            cospOUT%modis_Cloud_Fraction_Mid_Mean(:)                     = R_UNDEF
       if (associated(cospOUT%modis_Cloud_Fraction_Low_Mean))                            &
            cospOUT%modis_Cloud_Fraction_Low_Mean(:)                     = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_Total_Mean))                       &
            cospOUT%modis_Optical_Thickness_Total_Mean(:)                = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_Water_Mean))                       &
            cospOUT%modis_Optical_Thickness_Water_Mean(:)                = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_Ice_Mean))                         &
            cospOUT%modis_Optical_Thickness_Ice_Mean(:)                  = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_Total_LogMean))                    &
            cospOUT%modis_Optical_Thickness_Total_LogMean(:)             = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_Water_LogMean))                    &
            cospOUT%modis_Optical_Thickness_Water_LogMean(:)             = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_Ice_LogMean))                      &
            cospOUT%modis_Optical_Thickness_Ice_LogMean(:)               = R_UNDEF
       if (associated(cospOUT%modis_Cloud_Particle_Size_Water_Mean))                     &
            cospOUT%modis_Cloud_Particle_Size_Water_Mean(:)              = R_UNDEF
       if (associated(cospOUT%modis_Cloud_Particle_Size_Ice_Mean))                       &
            cospOUT%modis_Cloud_Particle_Size_Ice_Mean(:)                = R_UNDEF
       if (associated(cospOUT%modis_Cloud_Top_Pressure_Total_Mean))                      &
            cospOUT%modis_Cloud_Top_Pressure_Total_Mean(:)               = R_UNDEF
       if (associated(cospOUT%modis_Liquid_Water_Path_Mean))                             &
            cospOUT%modis_Liquid_Water_Path_Mean(:)                      = R_UNDEF
       if (associated(cospOUT%modis_Ice_Water_Path_Mean))                                &
            cospOUT%modis_Ice_Water_Path_Mean(:)                         = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure))            &
            cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure(:,:,:) = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_vs_ReffICE))                       &
            cospOUT%modis_Optical_Thickness_vs_ReffICE(:,:,:)            = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_vs_ReffLIQ))                       &
            cospOUT%modis_Optical_Thickness_vs_ReffLIQ(:,:,:)            = R_UNDEF        
       if (associated(cospOUT%calipso_cfad_sr))       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF
       if (associated(cospOUT%calipso_lidarcld))      cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldphase)) cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_cldlayerphase)) cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldtmp))   cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF      
    endif
    if (any(cospgridIN%qv .lt. 0)) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospgridIN%qv contains values out of range'
       Lisccp_subcolumn = .false.
       Lisccp_column    = .false.     
       Lrttov_subcolumn = .false.
       if (associated(cospOUT%rttov_tbs))           cospOUT%rttov_tbs(:,:)         = R_UNDEF       
       if (associated(cospOUT%isccp_totalcldarea))  cospOUT%isccp_totalcldarea(:)  = R_UNDEF
       if (associated(cospOUT%isccp_meantb))        cospOUT%isccp_meantb(:)        = R_UNDEF
       if (associated(cospOUT%isccp_meantbclr))     cospOUT%isccp_meantbclr(:)     = R_UNDEF
       if (associated(cospOUT%isccp_meanptop))      cospOUT%isccp_meanptop(:)      = R_UNDEF
       if (associated(cospOUT%isccp_meantaucld))    cospOUT%isccp_meantaucld(:)    = R_UNDEF
       if (associated(cospOUT%isccp_meanalbedocld)) cospOUT%isccp_meanalbedocld(:) = R_UNDEF
       if (associated(cospOUT%isccp_boxtau))        cospOUT%isccp_boxtau(:,:)      = R_UNDEF
       if (associated(cospOUT%isccp_boxptop))       cospOUT%isccp_boxptop(:,:)     = R_UNDEF
       if (associated(cospOUT%isccp_fq))            cospOUT%isccp_fq(:,:,:)        = R_UNDEF                
    endif
    if (any(cospgridIN%hgt_matrix .lt. -300)) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospgridIN%hgt_matrix contains values out of range'
       Lmisr_subcolumn     = .false.
       Lmisr_column        = .false.
       Lcloudsat_subcolumn = .false.
       Lcloudsat_column    = .false.
       Lcalipso_column     = .false.
       Lradar_lidar_tcc = .false.
       !Lcloudsat_tcc    = .false.
       !Lcloudsat_tcc2   = .false.
       Llidar_only_freq_cloud = .false.
       Lcloudsat_modis_wr = .false.
       if (associated(cospOUT%misr_fq))                   cospOUT%misr_fq(:,:,:)                 = R_UNDEF
       if (associated(cospOUT%misr_dist_model_layertops)) cospOUT%misr_dist_model_layertops(:,:) = R_UNDEF
       if (associated(cospOUT%misr_meanztop))             cospOUT%misr_meanztop(:)               = R_UNDEF
       if (associated(cospOUT%misr_cldarea))              cospOUT%misr_cldarea(:)                = R_UNDEF
       if (associated(cospOUT%calipso_cfad_sr))           cospOUT%calipso_cfad_sr(:,:,:)         = R_UNDEF
       if (associated(cospOUT%calipso_lidarcld))          cospOUT%calipso_lidarcld(:,:)          = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldphase))     cospOUT%calipso_lidarcldphase(:,:,:)   = R_UNDEF
       if (associated(cospOUT%calipso_cldlayer))          cospOUT%calipso_cldlayer(:,:)          = R_UNDEF
       if (associated(cospOUT%calipso_cldlayerphase))     cospOUT%calipso_cldlayerphase(:,:,:)   = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldtmp))       cospOUT%calipso_lidarcldtmp(:,:,:)     = R_UNDEF
   !    if (associated(cospOUT%calipso_lidarcldtype))      cospOUT%calipso_lidarcldtype(:,:,:)    = R_UNDEF
   !    if (associated(cospOUT%calipso_cldtype))           cospOUT%calipso_cldtype(:,:)           = R_UNDEF
   !    if (associated(cospOUT%calipso_cldtypetemp))       cospOUT%calipso_cldtypetemp(:,:)       = R_UNDEF
   !    if (associated(cospOUT%calipso_cldtypemeanz))      cospOUT%calipso_cldtypemeanz(:,:)      = R_UNDEF
   !    if (associated(cospOUT%calipso_cldtypemeanzse))    cospOUT%calipso_cldtypemeanzse(:,:)    = R_UNDEF
   !    if (associated(cospOUT%calipso_cldthinemis))       cospOUT%calipso_cldthinemis(:)         = R_UNDEF            
       if (associated(cospOUT%cloudsat_cfad_ze))          cospOUT%cloudsat_cfad_ze(:,:,:)        = R_UNDEF
       if (associated(cospOUT%cloudsat_Ze_tot))           cospOUT%cloudsat_Ze_tot(:,:,:)         = R_UNDEF
       if (associated(cospOUT%lidar_only_freq_cloud))     cospOUT%lidar_only_freq_cloud(:,:)     = R_UNDEF
       if (associated(cospOUT%radar_lidar_tcc))           cospOUT%radar_lidar_tcc(:)             = R_UNDEF
   !    if (associated(cospOUT%cloudsat_tcc))              cospOUT%cloudsat_tcc(:)                = R_UNDEF
   !    if (associated(cospOUT%cloudsat_tcc2))             cospOUT%cloudsat_tcc2(:)               = R_UNDEF
       if (associated(cospOUT%cfodd_ntotal)) cospOUT%cfodd_ntotal(:,:,:,:) = R_UNDEF !<-added
       if (associated(cospOUT%wr_occfreq_ntotal)) cospOUT%wr_occfreq_ntotal(:,:) = R_UNDEF !<-added
    endif
    if (any(cospgridIN%hgt_matrix_half .lt. -300)) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospgridIN%hgt_matrix_half contains values out of range'
       Lrttov_subcolumn = .false.
       Lcloudsat_column = .false.
       Lcalipso_column  = .false.
       Lradar_lidar_tcc = .false.
       !Lcloudsat_tcc    = .false.
       !Lcloudsat_tcc2   = .false.
       Llidar_only_freq_cloud = .false.
       Lcloudsat_modis_wr = .false.
       if (associated(cospOUT%rttov_tbs))             cospOUT%rttov_tbs(:,:)               = R_UNDEF       
       if (associated(cospOUT%calipso_cfad_sr))       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF
       if (associated(cospOUT%calipso_lidarcld))      cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldphase)) cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_cldlayerphase)) cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldtmp))   cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF            
       if (associated(cospOUT%cloudsat_cfad_ze))      cospOUT%cloudsat_cfad_ze(:,:,:)      = R_UNDEF
       !if (associated(cospOUT%cloudsat_tcc))          cospOUT%cloudsat_tcc(:)              = R_UNDEF
       !if (associated(cospOUT%cloudsat_tcc2))         cospOUT%cloudsat_tcc2(:)             = R_UNDEF
       if (associated(cospOUT%lidar_only_freq_cloud)) cospOUT%lidar_only_freq_cloud(:,:)   = R_UNDEF
       if (associated(cospOUT%radar_lidar_tcc))       cospOUT%radar_lidar_tcc(:)           = R_UNDEF                 
    endif
    if (any(cospgridIN%land .lt. 0)) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospgridIN%land contains values out of range'
       Lrttov_subcolumn = .false.
       Lcalipso_column  = .false.       
       Lparasol_column  = .false.
       if (associated(cospOUT%rttov_tbs))             cospOUT%rttov_tbs(:,:)               = R_UNDEF       
       if (associated(cospOUT%calipso_cfad_sr))       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF
       if (associated(cospOUT%calipso_lidarcld))      cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldphase)) cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_cldlayerphase)) cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldtmp))   cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF
       if (associated(cospOUT%parasolGrid_refl))      cospOUT%parasolGrid_refl(:,:)        = R_UNDEF
    endif
    if (any(cospgridIN%skt .lt. 0)) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospgridIN%skt contains values out of range'
       Lisccp_subcolumn = .false.
       Lisccp_column    = .false.     
       Lrttov_subcolumn = .false.
       if (associated(cospOUT%rttov_tbs))           cospOUT%rttov_tbs(:,:)         = R_UNDEF       
       if (associated(cospOUT%isccp_totalcldarea))  cospOUT%isccp_totalcldarea(:)  = R_UNDEF
       if (associated(cospOUT%isccp_meantb))        cospOUT%isccp_meantb(:)        = R_UNDEF
       if (associated(cospOUT%isccp_meantbclr))     cospOUT%isccp_meantbclr(:)     = R_UNDEF
       if (associated(cospOUT%isccp_meanptop))      cospOUT%isccp_meanptop(:)      = R_UNDEF
       if (associated(cospOUT%isccp_meantaucld))    cospOUT%isccp_meantaucld(:)    = R_UNDEF
       if (associated(cospOUT%isccp_meanalbedocld)) cospOUT%isccp_meanalbedocld(:) = R_UNDEF
       if (associated(cospOUT%isccp_boxtau))        cospOUT%isccp_boxtau(:,:)      = R_UNDEF
       if (associated(cospOUT%isccp_boxptop))       cospOUT%isccp_boxptop(:,:)     = R_UNDEF
       if (associated(cospOUT%isccp_fq))            cospOUT%isccp_fq(:,:,:)        = R_UNDEF     
    endif

	! RTTOV Inputs
    if (cospgridIN%zenang .lt. -90. .OR. cospgridIN%zenang .gt. 90) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospgridIN%zenang contains values out of range'
       Lrttov_subcolumn = .false.
       if (associated(cospOUT%rttov_tbs)) cospOUT%rttov_tbs(:,:) = R_UNDEF       
    endif
    if (cospgridIN%co2 .lt. 0) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospgridIN%co2 contains values out of range'
       Lrttov_subcolumn = .false.
       if (associated(cospOUT%rttov_tbs)) cospOUT%rttov_tbs(:,:) = R_UNDEF       
    endif
    if (cospgridIN%ch4 .lt. 0) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospgridIN%ch4 contains values out of range'
       Lrttov_subcolumn = .false.
       if (associated(cospOUT%rttov_tbs)) cospOUT%rttov_tbs(:,:) = R_UNDEF       
    endif
    if (cospgridIN%n2o .lt. 0) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospgridIN%n2o contains values out of range'
       Lrttov_subcolumn = .false.
       if (associated(cospOUT%rttov_tbs)) cospOUT%rttov_tbs(:,:) = R_UNDEF       
    endif
    if (cospgridIN%co.lt. 0) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospgridIN%co contains values out of range'
       Lrttov_subcolumn = .false.
       if (associated(cospOUT%rttov_tbs)) cospOUT%rttov_tbs(:,:) = R_UNDEF       
    endif
    if (any(cospgridIN%o3 .lt. 0)) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospgridIN%o3 contains values out of range'
       Lrttov_subcolumn = .false.
       if (associated(cospOUT%rttov_tbs)) cospOUT%rttov_tbs(:,:) = R_UNDEF       
    endif
    if (any(cospgridIN%emis_sfc .lt. 0. .OR. cospgridIN%emis_sfc .gt. 1)) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospgridIN%emis_sfc contains values out of range'
       Lrttov_subcolumn = .false.
       if (associated(cospOUT%rttov_tbs)) cospOUT%rttov_tbs(:,:) = R_UNDEF       
    endif
    if (any(cospgridIN%u_sfc .lt. -100. .OR. cospgridIN%u_sfc .gt. 100.)) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospIN%u_sfc contains values out of range'
       if (associated(cospOUT%rttov_tbs)) cospOUT%rttov_tbs(:,:) = R_UNDEF       
       Lrttov_subcolumn = .false.
    endif
    if (any(cospgridIN%v_sfc .lt. -100. .OR. cospgridIN%v_sfc .gt. 100.)) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospIN%v_sfc contains values out of range'
       Lrttov_subcolumn = .false.
       if (associated(cospOUT%rttov_tbs)) cospOUT%rttov_tbs(:,:) = R_UNDEF       
    endif
    if (any(cospgridIN%lat .lt. -90 .OR. cospgridIN%lat .gt. 90)) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospIN%lat contains values out of range'
       Lrttov_subcolumn = .false.
       if (associated(cospOUT%rttov_tbs)) cospOUT%rttov_tbs(:,:) = R_UNDEF       
    endif

    ! COSP_INPUTS
    if (cospIN%emsfc_lw .lt. 0. .OR. cospIN%emsfc_lw .gt. 1.) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospIN%emsfc_lw contains values out of range'
       Lisccp_subcolumn = .false.
       Lisccp_column    = .false.
       if (associated(cospOUT%isccp_totalcldarea))  cospOUT%isccp_totalcldarea(:)  = R_UNDEF
       if (associated(cospOUT%isccp_meantb))        cospOUT%isccp_meantb(:)        = R_UNDEF
       if (associated(cospOUT%isccp_meantbclr))     cospOUT%isccp_meantbclr(:)     = R_UNDEF
       if (associated(cospOUT%isccp_meanptop))      cospOUT%isccp_meanptop(:)      = R_UNDEF
       if (associated(cospOUT%isccp_meantaucld))    cospOUT%isccp_meantaucld(:)    = R_UNDEF
       if (associated(cospOUT%isccp_meanalbedocld)) cospOUT%isccp_meanalbedocld(:) = R_UNDEF
       if (associated(cospOUT%isccp_boxtau))        cospOUT%isccp_boxtau(:,:)      = R_UNDEF
       if (associated(cospOUT%isccp_boxptop))       cospOUT%isccp_boxptop(:,:)     = R_UNDEF
       if (associated(cospOUT%isccp_fq))            cospOUT%isccp_fq(:,:,:)        = R_UNDEF  
       
    endif
    if (any(cospIN%tau_067 .lt. 0)) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospIN%tau_067 contains values out of range'
       Lisccp_subcolumn = .false.
       Lisccp_column    = .false.
       Lmisr_subcolumn  = .false.
       Lmisr_column     = .false.
       Lmodis_subcolumn = .false.
       Lmodis_column    = .false.
       
       if (associated(cospOUT%isccp_totalcldarea))  cospOUT%isccp_totalcldarea(:)  = R_UNDEF
       if (associated(cospOUT%isccp_meantb))        cospOUT%isccp_meantb(:)        = R_UNDEF
       if (associated(cospOUT%isccp_meantbclr))     cospOUT%isccp_meantbclr(:)     = R_UNDEF
       if (associated(cospOUT%isccp_meanptop))      cospOUT%isccp_meanptop(:)      = R_UNDEF
       if (associated(cospOUT%isccp_meantaucld))    cospOUT%isccp_meantaucld(:)    = R_UNDEF
       if (associated(cospOUT%isccp_meanalbedocld)) cospOUT%isccp_meanalbedocld(:) = R_UNDEF
       if (associated(cospOUT%isccp_boxtau))        cospOUT%isccp_boxtau(:,:)      = R_UNDEF
       if (associated(cospOUT%isccp_boxptop))       cospOUT%isccp_boxptop(:,:)     = R_UNDEF
       if (associated(cospOUT%isccp_fq))            cospOUT%isccp_fq(:,:,:)        = R_UNDEF
       if (associated(cospOUT%misr_fq))                   cospOUT%misr_fq(:,:,:)                 = R_UNDEF
       if (associated(cospOUT%misr_dist_model_layertops)) cospOUT%misr_dist_model_layertops(:,:) = R_UNDEF
       if (associated(cospOUT%misr_meanztop))             cospOUT%misr_meanztop(:)               = R_UNDEF
       if (associated(cospOUT%misr_cldarea))              cospOUT%misr_cldarea(:)                = R_UNDEF
       if (associated(cospOUT%modis_Cloud_Fraction_Total_Mean))                          &
            cospOUT%modis_Cloud_Fraction_Total_Mean(:)                   = R_UNDEF
       if (associated(cospOUT%modis_Cloud_Fraction_Water_Mean))                          &
            cospOUT%modis_Cloud_Fraction_Water_Mean(:)                   = R_UNDEF
       if (associated(cospOUT%modis_Cloud_Fraction_Ice_Mean))                            &
            cospOUT%modis_Cloud_Fraction_Ice_Mean(:)                     = R_UNDEF
       if (associated(cospOUT%modis_Cloud_Fraction_High_Mean))                           &
            cospOUT%modis_Cloud_Fraction_High_Mean(:)                    = R_UNDEF
       if (associated(cospOUT%modis_Cloud_Fraction_Mid_Mean))                            &
            cospOUT%modis_Cloud_Fraction_Mid_Mean(:)                     = R_UNDEF
       if (associated(cospOUT%modis_Cloud_Fraction_Low_Mean))                            &
            cospOUT%modis_Cloud_Fraction_Low_Mean(:)                     = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_Total_Mean))                       &
            cospOUT%modis_Optical_Thickness_Total_Mean(:)                = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_Water_Mean))                       &
            cospOUT%modis_Optical_Thickness_Water_Mean(:)                = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_Ice_Mean))                         &
            cospOUT%modis_Optical_Thickness_Ice_Mean(:)                  = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_Total_LogMean))                    &
            cospOUT%modis_Optical_Thickness_Total_LogMean(:)             = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_Water_LogMean))                    &
            cospOUT%modis_Optical_Thickness_Water_LogMean(:)             = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_Ice_LogMean))                      &
            cospOUT%modis_Optical_Thickness_Ice_LogMean(:)               = R_UNDEF
       if (associated(cospOUT%modis_Cloud_Particle_Size_Water_Mean))                     &
            cospOUT%modis_Cloud_Particle_Size_Water_Mean(:)              = R_UNDEF
       if (associated(cospOUT%modis_Cloud_Particle_Size_Ice_Mean))                       &
            cospOUT%modis_Cloud_Particle_Size_Ice_Mean(:)                = R_UNDEF
       if (associated(cospOUT%modis_Cloud_Top_Pressure_Total_Mean))                      &
            cospOUT%modis_Cloud_Top_Pressure_Total_Mean(:)               = R_UNDEF
       if (associated(cospOUT%modis_Liquid_Water_Path_Mean))                             &
            cospOUT%modis_Liquid_Water_Path_Mean(:)                      = R_UNDEF
       if (associated(cospOUT%modis_Ice_Water_Path_Mean))                                &
            cospOUT%modis_Ice_Water_Path_Mean(:)                         = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure))            &
            cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure(:,:,:) = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_vs_ReffICE))                       &
            cospOUT%modis_Optical_Thickness_vs_ReffICE(:,:,:)            = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_vs_ReffLIQ))                       &
            cospOUT%modis_Optical_Thickness_vs_ReffLIQ(:,:,:)            = R_UNDEF        
       
    endif
    if (any(cospIN%emiss_11 .lt. 0. .OR. cospIN%emiss_11 .gt. 1)) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospIN%emiss_11 contains values out of range'
       Lisccp_subcolumn = .false.
       Lisccp_column    = .false.
       if (associated(cospOUT%isccp_totalcldarea))  cospOUT%isccp_totalcldarea(:)  = R_UNDEF
       if (associated(cospOUT%isccp_meantb))        cospOUT%isccp_meantb(:)        = R_UNDEF
       if (associated(cospOUT%isccp_meantbclr))     cospOUT%isccp_meantbclr(:)     = R_UNDEF
       if (associated(cospOUT%isccp_meanptop))      cospOUT%isccp_meanptop(:)      = R_UNDEF
       if (associated(cospOUT%isccp_meantaucld))    cospOUT%isccp_meantaucld(:)    = R_UNDEF
       if (associated(cospOUT%isccp_meanalbedocld)) cospOUT%isccp_meanalbedocld(:) = R_UNDEF
       if (associated(cospOUT%isccp_boxtau))        cospOUT%isccp_boxtau(:,:)      = R_UNDEF
       if (associated(cospOUT%isccp_boxptop))       cospOUT%isccp_boxptop(:,:)     = R_UNDEF
       if (associated(cospOUT%isccp_fq))            cospOUT%isccp_fq(:,:,:)        = R_UNDEF
         
    endif
    if (any(cospIN%asym .lt. -1. .OR. cospIN%asym .gt. 1)) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospIN%asym contains values out of range'
       Lmodis_subcolumn = .false.
       Lmodis_column    = .false.
       if (associated(cospOUT%modis_Cloud_Fraction_Total_Mean))                          &
            cospOUT%modis_Cloud_Fraction_Total_Mean(:)                   = R_UNDEF
       if (associated(cospOUT%modis_Cloud_Fraction_Water_Mean))                          &
            cospOUT%modis_Cloud_Fraction_Water_Mean(:)                   = R_UNDEF
       if (associated(cospOUT%modis_Cloud_Fraction_Ice_Mean))                            &
            cospOUT%modis_Cloud_Fraction_Ice_Mean(:)                     = R_UNDEF
       if (associated(cospOUT%modis_Cloud_Fraction_High_Mean))                           &
            cospOUT%modis_Cloud_Fraction_High_Mean(:)                    = R_UNDEF
       if (associated(cospOUT%modis_Cloud_Fraction_Mid_Mean))                            &
            cospOUT%modis_Cloud_Fraction_Mid_Mean(:)                     = R_UNDEF
       if (associated(cospOUT%modis_Cloud_Fraction_Low_Mean))                            &
            cospOUT%modis_Cloud_Fraction_Low_Mean(:)                     = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_Total_Mean))                       &
            cospOUT%modis_Optical_Thickness_Total_Mean(:)                = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_Water_Mean))                       &
            cospOUT%modis_Optical_Thickness_Water_Mean(:)                = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_Ice_Mean))                         &
            cospOUT%modis_Optical_Thickness_Ice_Mean(:)                  = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_Total_LogMean))                    &
            cospOUT%modis_Optical_Thickness_Total_LogMean(:)             = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_Water_LogMean))                    &
            cospOUT%modis_Optical_Thickness_Water_LogMean(:)             = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_Ice_LogMean))                      &
            cospOUT%modis_Optical_Thickness_Ice_LogMean(:)               = R_UNDEF
       if (associated(cospOUT%modis_Cloud_Particle_Size_Water_Mean))                     &
            cospOUT%modis_Cloud_Particle_Size_Water_Mean(:)              = R_UNDEF
       if (associated(cospOUT%modis_Cloud_Particle_Size_Ice_Mean))                       &
            cospOUT%modis_Cloud_Particle_Size_Ice_Mean(:)                = R_UNDEF
       if (associated(cospOUT%modis_Cloud_Top_Pressure_Total_Mean))                      &
            cospOUT%modis_Cloud_Top_Pressure_Total_Mean(:)               = R_UNDEF
       if (associated(cospOUT%modis_Liquid_Water_Path_Mean))                             &
            cospOUT%modis_Liquid_Water_Path_Mean(:)                      = R_UNDEF
       if (associated(cospOUT%modis_Ice_Water_Path_Mean))                                &
            cospOUT%modis_Ice_Water_Path_Mean(:)                         = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure))            &
            cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure(:,:,:) = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_vs_ReffICE))                       &
            cospOUT%modis_Optical_Thickness_vs_ReffICE(:,:,:)            = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_vs_ReffLIQ))                       &
            cospOUT%modis_Optical_Thickness_vs_ReffLIQ(:,:,:)            = R_UNDEF             
    endif
    if (any(cospIN%ss_alb .lt. 0 .OR. cospIN%ss_alb .gt. 1)) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospIN%ss_alb contains values out of range'
       Lmodis_subcolumn = .false.
       Lmodis_column    = .false.
       if (associated(cospOUT%modis_Cloud_Fraction_Total_Mean))                          &
            cospOUT%modis_Cloud_Fraction_Total_Mean(:)                   = R_UNDEF
       if (associated(cospOUT%modis_Cloud_Fraction_Water_Mean))                          &
            cospOUT%modis_Cloud_Fraction_Water_Mean(:)                   = R_UNDEF
       if (associated(cospOUT%modis_Cloud_Fraction_Ice_Mean))                            &
            cospOUT%modis_Cloud_Fraction_Ice_Mean(:)                     = R_UNDEF
       if (associated(cospOUT%modis_Cloud_Fraction_High_Mean))                           &
            cospOUT%modis_Cloud_Fraction_High_Mean(:)                    = R_UNDEF
       if (associated(cospOUT%modis_Cloud_Fraction_Mid_Mean))                            &
            cospOUT%modis_Cloud_Fraction_Mid_Mean(:)                     = R_UNDEF
       if (associated(cospOUT%modis_Cloud_Fraction_Low_Mean))                            &
            cospOUT%modis_Cloud_Fraction_Low_Mean(:)                     = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_Total_Mean))                       &
            cospOUT%modis_Optical_Thickness_Total_Mean(:)                = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_Water_Mean))                       &
            cospOUT%modis_Optical_Thickness_Water_Mean(:)                = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_Ice_Mean))                         &
            cospOUT%modis_Optical_Thickness_Ice_Mean(:)                  = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_Total_LogMean))                    &
            cospOUT%modis_Optical_Thickness_Total_LogMean(:)             = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_Water_LogMean))                    &
            cospOUT%modis_Optical_Thickness_Water_LogMean(:)             = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_Ice_LogMean))                      &
            cospOUT%modis_Optical_Thickness_Ice_LogMean(:)               = R_UNDEF
       if (associated(cospOUT%modis_Cloud_Particle_Size_Water_Mean))                     &
            cospOUT%modis_Cloud_Particle_Size_Water_Mean(:)              = R_UNDEF
       if (associated(cospOUT%modis_Cloud_Particle_Size_Ice_Mean))                       &
            cospOUT%modis_Cloud_Particle_Size_Ice_Mean(:)                = R_UNDEF
       if (associated(cospOUT%modis_Cloud_Top_Pressure_Total_Mean))                      &
            cospOUT%modis_Cloud_Top_Pressure_Total_Mean(:)               = R_UNDEF
       if (associated(cospOUT%modis_Liquid_Water_Path_Mean))                             &
            cospOUT%modis_Liquid_Water_Path_Mean(:)                      = R_UNDEF
       if (associated(cospOUT%modis_Ice_Water_Path_Mean))                                &
            cospOUT%modis_Ice_Water_Path_Mean(:)                         = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure))            &
            cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure(:,:,:) = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_vs_ReffICE))                       &
            cospOUT%modis_Optical_Thickness_vs_ReffICE(:,:,:)            = R_UNDEF
       if (associated(cospOUT%modis_Optical_Thickness_vs_ReffLIQ))                       &
            cospOUT%modis_Optical_Thickness_vs_ReffLIQ(:,:,:)            = R_UNDEF                 
    endif
    if (any(cospIN%betatot .lt. 0)) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospIN%betatot contains values out of range'
       Lcalipso_subcolumn = .false.
       Lcalipso_column    = .false.
       if (associated(cospOUT%calipso_cfad_sr))       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF
       if (associated(cospOUT%calipso_lidarcld))      cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldphase)) cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_cldlayerphase)) cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldtmp))   cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF
       if (associated(cospOUT%calipso_srbval))        cospOUT%calipso_srbval(:)            = R_UNDEF
    endif
    if (any(cospIN%betatot_liq .lt. 0)) then
       nError=nError+1
       errorMessage(nError) = ('ERROR: COSP input variable: cospIN%betatot_liq contains values out of range')
       Lcalipso_subcolumn = .false.
       Lcalipso_column    = .false.
       if (associated(cospOUT%calipso_cfad_sr))       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF
       if (associated(cospOUT%calipso_lidarcld))      cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldphase)) cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_cldlayerphase)) cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldtmp))   cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF
       if (associated(cospOUT%calipso_srbval))        cospOUT%calipso_srbval(:)            = R_UNDEF       
    endif
    if (any(cospIN%betatot_ice .lt. 0)) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospIN%betatot_ice contains values out of range'
       Lcalipso_subcolumn = .false.
       Lcalipso_column    = .false.
       if (associated(cospOUT%calipso_cfad_sr))       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF
       if (associated(cospOUT%calipso_lidarcld))      cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldphase)) cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_cldlayerphase)) cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldtmp))   cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF
       if (associated(cospOUT%calipso_srbval))        cospOUT%calipso_srbval(:)            = R_UNDEF
    endif 
    if (any(cospIN%beta_mol .lt. 0)) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospIN%beta_mol contains values out of range'
       Lcalipso_subcolumn = .false.
       Lcalipso_column    = .false.
       Lcloudsat_column   = .false.
       Lradar_lidar_tcc = .false.
       Llidar_only_freq_cloud = .false.
       if (associated(cospOUT%calipso_cfad_sr))       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF
       if (associated(cospOUT%calipso_lidarcld))      cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldphase)) cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_cldlayerphase)) cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldtmp))   cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF
       if (associated(cospOUT%calipso_srbval))        cospOUT%calipso_srbval(:)            = R_UNDEF
       if (associated(cospOUT%cloudsat_cfad_ze))      cospOUT%cloudsat_cfad_ze(:,:,:)      = R_UNDEF
       if (associated(cospOUT%lidar_only_freq_cloud)) cospOUT%lidar_only_freq_cloud(:,:)   = R_UNDEF
       if (associated(cospOUT%radar_lidar_tcc))       cospOUT%radar_lidar_tcc(:)           = R_UNDEF          
    endif    
    if (any(cospIN%tautot .lt. 0)) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospIN%tautot contains values out of range'
       Lcalipso_subcolumn = .false.
       Lcalipso_column    = .false.
       if (associated(cospOUT%calipso_cfad_sr))       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF
       if (associated(cospOUT%calipso_lidarcld))      cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldphase)) cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_cldlayerphase)) cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldtmp))   cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF
       if (associated(cospOUT%calipso_srbval))        cospOUT%calipso_srbval(:)            = R_UNDEF      
    endif
    if (any(cospIN%tautot_liq .lt. 0)) then
       nError=nError+1
       errorMessage(nError) = ('ERROR: COSP input variable: cospIN%tautot_liq contains values out of range')
       Lcalipso_subcolumn = .false.
       Lcalipso_column    = .false.
       if (associated(cospOUT%calipso_cfad_sr))       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF
       if (associated(cospOUT%calipso_lidarcld))      cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldphase)) cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_cldlayerphase)) cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldtmp))   cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF
       if (associated(cospOUT%calipso_srbval))        cospOUT%calipso_srbval(:)            = R_UNDEF       
    endif
    if (any(cospIN%tautot_ice .lt. 0)) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospIN%tautot_ice contains values out of range'
       Lcalipso_subcolumn = .false.
       Lcalipso_column    = .false.
       if (associated(cospOUT%calipso_cfad_sr))       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF
       if (associated(cospOUT%calipso_lidarcld))      cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldphase)) cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_cldlayerphase)) cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldtmp))   cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF
       if (associated(cospOUT%calipso_srbval))        cospOUT%calipso_srbval(:)            = R_UNDEF        
    endif
    if (any(cospIN%tau_mol .lt. 0)) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospIN%tau_mol contains values out of range'
       Lcalipso_subcolumn = .false.
       Lcalipso_column    = .false.
       if (associated(cospOUT%calipso_cfad_sr))       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF
       if (associated(cospOUT%calipso_lidarcld))      cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldphase)) cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_cldlayerphase)) cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldtmp))   cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF
       if (associated(cospOUT%calipso_srbval))        cospOUT%calipso_srbval(:)            = R_UNDEF        
   !    if (associated(cospOUT%calipso_lidarcldtype))  cospOUT%calipso_lidarcldtype(:,:,:)  = R_UNDEF
   !    if (associated(cospOUT%calipso_cldtype))       cospOUT%calipso_cldtype(:,:)         = R_UNDEF 
   !    if (associated(cospOUT%calipso_cldtypetemp))   cospOUT%calipso_cldtypetemp(:,:)     = R_UNDEF
   !    if (associated(cospOUT%calipso_cldtypemeanz))  cospOUT%calipso_cldtypemeanz(:,:)    = R_UNDEF
   !    if (associated(cospOUT%calipso_cldtypemeanzse)) cospOUT%calipso_cldtypemeanzse(:,:) = R_UNDEF
   !    if (associated(cospOUT%calipso_cldthinemis))   cospOUT%calipso_cldthinemis(:)       = R_UNDEF    
    endif    
    if (any(cospIN%tautot_S_liq .lt. 0)) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospIN%tautot_S_liq contains values out of range'
       Lparasol_subcolumn = .false.
       Lparasol_column    = .false.
       if (associated(cospOUT%parasolPix_refl))  cospOUT%parasolPix_refl(:,:,:) = R_UNDEF      
       if (associated(cospOUT%parasolGrid_refl)) cospOUT%parasolGrid_refl(:,:)  = R_UNDEF
    endif
    if (any(cospIN%tautot_S_ice .lt. 0)) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospIN%tautot_S_ice contains values out of range'
       Lparasol_subcolumn = .false.
       Lparasol_column    = .false.
       if (associated(cospOUT%parasolPix_refl))  cospOUT%parasolPix_refl(:,:,:) = R_UNDEF      
       if (associated(cospOUT%parasolGrid_refl)) cospOUT%parasolGrid_refl(:,:)  = R_UNDEF       
    endif    
    if (any(cospIN%z_vol_cloudsat .lt. 0)) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospIN%z_vol_cloudsat contains values out of range'
       Lcloudsat_subcolumn = .false.
       Lcloudsat_column    = .false.
       Lradar_lidar_tcc    = .false.
       !Lcloudsat_tcc = .false.
       !Lcloudsat_tcc2 = .false.
       Llidar_only_freq_cloud = .false.
       Lcloudsat_modis_wr = .false.
       if (associated(cospOUT%cloudsat_cfad_ze))          cospOUT%cloudsat_cfad_ze(:,:,:)        = R_UNDEF
       if (associated(cospOUT%cloudsat_Ze_tot))           cospOUT%cloudsat_Ze_tot(:,:,:)         = R_UNDEF
       !if (associated(cospOUT%cloudsat_tcc)) cospOUT%cloudsat_tcc(:) = R_UNDEF
       !if (associated(cospOUT%cloudsat_tcc2)) cospOUT%cloudsat_tcc2(:) = R_UNDEF
       if (associated(cospOUT%lidar_only_freq_cloud))     cospOUT%lidar_only_freq_cloud(:,:)     = R_UNDEF
       if (associated(cospOUT%radar_lidar_tcc))           cospOUT%radar_lidar_tcc(:)             = R_UNDEF     
       if (associated(cospOUT%cfodd_ntotal))              cospOUT%cfodd_ntotal(:,:,:,:)          = R_UNDEF
       if (associated(cospOUT%wr_occfreq_ntotal))         cospOUT%wr_occfreq_ntotal(:,:)         = R_UNDEF
    endif
    if (any(cospIN%kr_vol_cloudsat .lt. 0)) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospIN%kr_vol_cloudsat contains values out of range'
       Lcloudsat_subcolumn = .false.
       Lcloudsat_column    = .false.
       Lradar_lidar_tcc    = .false.
       !Lcloudsat_tcc = .false.
       !Lcloudsat_tcc2 = .false.
       Llidar_only_freq_cloud = .false.
       Lcloudsat_modis_wr = .false.
       if (associated(cospOUT%cloudsat_cfad_ze))          cospOUT%cloudsat_cfad_ze(:,:,:)        = R_UNDEF
       if (associated(cospOUT%cloudsat_Ze_tot))           cospOUT%cloudsat_Ze_tot(:,:,:)         = R_UNDEF
       !if (associated(cospOUT%cloudsat_tcc)) cospOUT%cloudsat_tcc(:) = R_UNDEF
       !if (associated(cospOUT%cloudsat_tcc2)) cospOUT%cloudsat_tcc2(:) = R_UNDEF
       if (associated(cospOUT%lidar_only_freq_cloud))     cospOUT%lidar_only_freq_cloud(:,:)     = R_UNDEF
       if (associated(cospOUT%radar_lidar_tcc))           cospOUT%radar_lidar_tcc(:)             = R_UNDEF      
       if (associated(cospOUT%cfodd_ntotal))              cospOUT%cfodd_ntotal(:,:,:,:)          = R_UNDEF
       if (associated(cospOUT%wr_occfreq_ntotal))         cospOUT%wr_occfreq_ntotal(:,:)         = R_UNDEF
    endif    
    if (any(cospIN%g_vol_cloudsat .lt. 0)) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospIN%g_vol_cloudsat contains values out of range'
       Lcloudsat_subcolumn = .false.
       Lcloudsat_column    = .false.
       Lradar_lidar_tcc    = .false.
       !Lcloudsat_tcc = .false.
       !Lcloudsat_tcc2 = .false.
       Llidar_only_freq_cloud = .false.
       Lcloudsat_modis_wr = .false.
       if (associated(cospOUT%cloudsat_cfad_ze))          cospOUT%cloudsat_cfad_ze(:,:,:)        = R_UNDEF
       if (associated(cospOUT%cloudsat_Ze_tot))           cospOUT%cloudsat_Ze_tot(:,:,:)         = R_UNDEF
       !if (associated(cospOUT%cloudsat_tcc)) cospOUT%cloudsat_tcc(:) = R_UNDEF
       !if (associated(cospOUT%cloudsat_tcc2)) cospOUT%cloudsat_tcc2(:) = R_UNDEF 
       if (associated(cospOUT%lidar_only_freq_cloud))     cospOUT%lidar_only_freq_cloud(:,:)     = R_UNDEF
       if (associated(cospOUT%radar_lidar_tcc))           cospOUT%radar_lidar_tcc(:)             = R_UNDEF
       if (associated(cospOUT%cfodd_ntotal))              cospOUT%cfodd_ntotal(:,:,:,:)          = R_UNDEF
       if (associated(cospOUT%wr_occfreq_ntotal))         cospOUT%wr_occfreq_ntotal(:,:)         = R_UNDEF
    endif   
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Part 2: Check input fields array size for consistency. This needs to be done for each
  !         simulator
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! ISCCP
  if (size(cospIN%frac_out,1)  .ne. cospIN%Npoints .OR. &
      size(cospIN%tau_067,1)   .ne. cospIN%Npoints .OR. &
      size(cospIN%emiss_11,1)  .ne. cospIN%Npoints .OR. &
      size(cospgridIN%skt)     .ne. cospIN%Npoints .OR. &
      size(cospgridIN%qv,1)    .ne. cospIN%Npoints .OR. &
      size(cospgridIN%at,1)    .ne. cospIN%Npoints .OR. &
      size(cospgridIN%phalf,1) .ne. cospIN%Npoints .OR. &
      size(cospgridIN%sunlit)  .ne. cospIN%Npoints .OR. &
      size(cospgridIN%pfull,1) .ne. cospIN%Npoints) then
      Lisccp_subcolumn = .false.
      Lisccp_column    = .false.
      nError=nError+1
      errorMessage(nError) = 'ERROR(isccp_simulator): The number of points in the input fields are inconsistent'
  endif
  if (size(cospIN%frac_out,2) .ne. cospIN%Ncolumns .OR. &
      size(cospIN%tau_067,2)  .ne. cospIN%Ncolumns .OR. &
      size(cospIN%emiss_11,2) .ne. cospIN%Ncolumns) then
      Lisccp_subcolumn = .false.
      Lisccp_column    = .false.
      nError=nError+1
      errorMessage(nError) = 'ERROR(isccp_simulator): The number of sub-columns in the input fields are inconsistent'
  endif
  if (size(cospIN%frac_out,3)  .ne. cospIN%Nlevels .OR. &
      size(cospIN%tau_067,3)   .ne. cospIN%Nlevels .OR. &
      size(cospIN%emiss_11,3)  .ne. cospIN%Nlevels .OR. &
      size(cospgridIN%qv,2)    .ne. cospIN%Nlevels .OR. &
      size(cospgridIN%at,2)    .ne. cospIN%Nlevels .OR. &
      size(cospgridIN%pfull,2) .ne. cospIN%Nlevels .OR. &    
      size(cospgridIN%phalf,2) .ne. cospIN%Nlevels+1) then
      Lisccp_subcolumn = .false.
      Lisccp_column    = .false.
      nError=nError+1
      errorMessage(nError) = 'ERROR(isccp_simulator): The number of levels in the input fields are inconsistent'
  endif
      
  ! MISR
  if (size(cospIN%tau_067,1)        .ne. cospIN%Npoints .OR. &
      size(cospgridIN%sunlit)       .ne. cospIN%Npoints .OR. & 
      size(cospgridIN%hgt_matrix,1) .ne. cospIN%Npoints .OR. &
      size(cospgridIN%at,1)         .ne. cospIN%Npoints) then
      Lmisr_subcolumn = .false.
      Lmisr_column    = .false.
      nError=nError+1
      errorMessage(nError) = 'ERROR(misr_simulator): The number of points in the input fields are inconsistent'
  endif
  if (size(cospIN%tau_067,2) .ne. cospIN%Ncolumns) then
      Lmisr_subcolumn = .false.
      Lmisr_column    = .false.
      nError=nError+1
      errorMessage(nError) = 'ERROR(misr_simulator): The number of sub-columns in the input fields are inconsistent'
  endif
  if (size(cospIN%tau_067,3)        .ne. cospIN%Nlevels .OR. &
      size(cospgridIN%hgt_matrix,2) .ne. cospIN%Nlevels .OR. &
      size(cospgridIN%at,2)         .ne. cospIN%Nlevels) then
      Lmisr_subcolumn = .false.
      Lmisr_column    = .false.
      nError=nError+1
      errorMessage(nError) = 'ERROR(misr_simulator): The number of levels in the input fields are inconsistent'
  endif    

  ! MODIS
  if (size(cospIN%fracLiq,1) .ne. cospIN%Npoints .OR. &
      size(cospIN%tau_067,1) .ne. cospIN%Npoints .OR. &
      size(cospIN%asym,1)    .ne. cospIN%Npoints .OR. &
      size(cospIN%ss_alb,1)  .ne. cospIN%Npoints) then
      Lmodis_subcolumn = .false.
      Lmodis_column    = .false.
      nError=nError+1
      errorMessage(nError) = 'ERROR(modis_simulator): The number of points in the input fields are inconsistent'
  endif
  if (size(cospIN%fracLiq,2) .ne. cospIN%Ncolumns .OR. &
      size(cospIN%tau_067,2) .ne. cospIN%Ncolumns .OR. &
      size(cospIN%asym,2)    .ne. cospIN%Ncolumns .OR. &
      size(cospIN%ss_alb,2)  .ne. cospIN%Ncolumns) then
      Lmodis_subcolumn = .false.
      Lmodis_column    = .false.
      nError=nError+1
      errorMessage(nError) = 'ERROR(modis_simulator): The number of sub-columns in the input fields are inconsistent'
  endif        
  if (size(cospIN%fracLiq,3) .ne. cospIN%Nlevels .OR. &
      size(cospIN%tau_067,3) .ne. cospIN%Nlevels .OR. &
      size(cospIN%asym,3)    .ne. cospIN%Nlevels .OR. &
      size(cospIN%ss_alb,3)  .ne. cospIN%Nlevels) then
      Lmodis_subcolumn = .false.
      Lmodis_column    = .false.
      nError=nError+1
      errorMessage(nError) = 'ERROR(modis_simulator): The number of levels in the input fields are inconsistent'
  endif  
  
  ! CLOUDSAT    
  if (size(cospIN%z_vol_cloudsat,1)   .ne. cospIN%Npoints .OR. &
      size(cospIN%kr_vol_cloudsat,1)  .ne. cospIN%Npoints .OR. &
      size(cospIN%g_vol_cloudsat,1)   .ne. cospIN%Npoints .OR. &
      size(cospgridIN%hgt_matrix,1)   .ne. cospIN%Npoints) then
      Lcloudsat_subcolumn = .false.
      Lcloudsat_column    = .false.
      nError=nError+1
      errorMessage(nError) = 'ERROR(cloudsat_simulator): The number of points in the input fields are inconsistent'
  endif
  if (size(cospIN%z_vol_cloudsat,2)  .ne. cospIN%Ncolumns .OR. &
      size(cospIN%kr_vol_cloudsat,2) .ne. cospIN%Ncolumns .OR. &
      size(cospIN%g_vol_cloudsat,2)  .ne. cospIN%Ncolumns) then
      Lcloudsat_subcolumn = .false.
      Lcloudsat_column    = .false.
      nError=nError+1
      errorMessage(nError) = 'ERROR(cloudsat_simulator): The number of sub-columns in the input fields are inconsistent'
  endif       
  if (size(cospIN%z_vol_cloudsat,3)  .ne. cospIN%Nlevels .OR. &
      size(cospIN%kr_vol_cloudsat,3) .ne. cospIN%Nlevels .OR. &
      size(cospIN%g_vol_cloudsat,3)  .ne. cospIN%Nlevels .OR. &
      size(cospgridIN%hgt_matrix,2)  .ne. cospIN%Nlevels) then
      Lcloudsat_subcolumn = .false.
      Lcloudsat_column    = .false.
      nError=nError+1
      errorMessage(nError) = 'ERROR(cloudsat_simulator): The number of levels in the input fields are inconsistent'
  endif

  ! CALIPSO
  if (size(cospIN%beta_mol,1)    .ne. cospIN%Npoints .OR. &
      size(cospIN%betatot,1)     .ne. cospIN%Npoints .OR. &
      size(cospIN%betatot_liq,1) .ne. cospIN%Npoints .OR. &
      size(cospIN%betatot_ice,1) .ne. cospIN%Npoints .OR. &
      size(cospIN%tau_mol,1)     .ne. cospIN%Npoints .OR. &
      size(cospIN%tautot,1)      .ne. cospIN%Npoints .OR. &
      size(cospIN%tautot_liq,1)  .ne. cospIN%Npoints .OR. &
      size(cospIN%tautot_ice,1)  .ne. cospIN%Npoints) then
      Lcalipso_subcolumn = .false.
      Lcalipso_column    = .false.
      nError=nError+1
      errorMessage(nError) = 'ERROR(calipso_simulator): The number of points in the input fields are inconsistent'
  endif          
   if (size(cospIN%betatot,2)     .ne. cospIN%Ncolumns .OR. &
       size(cospIN%betatot_liq,2) .ne. cospIN%Ncolumns .OR. &
       size(cospIN%betatot_ice,2) .ne. cospIN%Ncolumns .OR. &
       size(cospIN%tautot,2)      .ne. cospIN%Ncolumns .OR. &
       size(cospIN%tautot_liq,2)  .ne. cospIN%Ncolumns .OR. &
       size(cospIN%tautot_ice,2)  .ne. cospIN%Ncolumns) then
       Lcalipso_subcolumn = .false.
       Lcalipso_column    = .false.
      nError=nError+1
      errorMessage(nError) = 'ERROR(calipso_simulator): The number of sub-columns in the input fields are inconsistent'
  endif       
  if (size(cospIN%beta_mol,2)    .ne. cospIN%Nlevels .OR. &
      size(cospIN%betatot,3)     .ne. cospIN%Nlevels .OR. &
      size(cospIN%betatot_liq,3) .ne. cospIN%Nlevels .OR. &
      size(cospIN%betatot_ice,3) .ne. cospIN%Nlevels .OR. &
      size(cospIN%tau_mol,2)     .ne. cospIN%Nlevels .OR. &
      size(cospIN%tautot,3)      .ne. cospIN%Nlevels .OR. &
      size(cospIN%tautot_liq,3)  .ne. cospIN%Nlevels .OR. &
      size(cospIN%tautot_ice,3)  .ne. cospIN%Nlevels) then
      Lcalipso_subcolumn = .false.
      Lcalipso_column    = .false.
      nError=nError+1
      errorMessage(nError) = 'ERROR(calipso_simulator): The number of levels in the input fields are inconsistent'
  endif 
  
  ! PARASOL
  if (size(cospIN%tautot_S_liq,1) .ne. cospIN%Npoints .OR. &
      size(cospIN%tautot_S_ice,1) .ne. cospIN%Npoints) then
      Lparasol_subcolumn = .false.
      Lparasol_column    = .false.
      nError=nError+1
      errorMessage(nError) = 'ERROR(parasol_simulator): The number of points in the input fields are inconsistent'
  endif
  if (size(cospIN%tautot_S_liq,2) .ne. cospIN%Ncolumns .OR. &
      size(cospIN%tautot_S_ice,2) .ne. cospIN%Ncolumns) then
      Lparasol_subcolumn = .false.
      Lparasol_column    = .false.
      nError=nError+1
      errorMessage(nError) = 'ERROR(parasol_simulator): The number of levels in the input fields are inconsistent'
  endif  
  
  ! RTTOV
  if (size(cospgridIN%pfull,1)           .ne. cospIN%Npoints .OR. &
      size(cospgridIN%at,1)              .ne. cospIN%Npoints .OR. &
      size(cospgridIN%qv,1)              .ne. cospIN%Npoints .OR. &
      size(cospgridIN%hgt_matrix_half,1) .ne. cospIN%Npoints .OR. &
      size(cospgridIN%u_sfc)             .ne. cospIN%Npoints .OR. &
      size(cospgridIN%v_sfc)             .ne. cospIN%Npoints .OR. &
      size(cospgridIN%skt)               .ne. cospIN%Npoints .OR. &
      size(cospgridIN%phalf,1)           .ne. cospIN%Npoints .OR. &
      size(cospgridIN%qv,1)              .ne. cospIN%Npoints .OR. &
      size(cospgridIN%land)              .ne. cospIN%Npoints .OR. &
      size(cospgridIN%lat)               .ne. cospIN%Npoints) then
      Lrttov_subcolumn = .false.
      Lrttov_column    = .false.
      nError=nError+1
      errorMessage(nError) = 'ERROR(rttov_simulator): The number of points in the input fields are inconsistent'
  endif      
  if (size(cospgridIN%pfull,2)           .ne. cospIN%Nlevels   .OR. &
      size(cospgridIN%at,2)              .ne. cospIN%Nlevels   .OR. &
      size(cospgridIN%qv,2)              .ne. cospIN%Nlevels   .OR. &
      size(cospgridIN%hgt_matrix_half,2) .ne. cospIN%Nlevels+1 .OR. &
      size(cospgridIN%phalf,2)           .ne. cospIN%Nlevels+1 .OR. &
      size(cospgridIN%qv,2)              .ne. cospIN%Nlevels) then
      Lrttov_subcolumn = .false.
      Lrttov_column    = .false.
      nError=nError+1
      errorMessage(nError) = 'ERROR(rttov_simulator): The number of levels in the input fields are inconsistent'
  endif
  
  end subroutine cosp_errorCheck
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! END MODULE
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
END MODULE MOD_COSP
