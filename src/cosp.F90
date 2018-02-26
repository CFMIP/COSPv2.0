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
! May 2015- D. Swales - Original version
! Jul 2017- R. Guzman - Added OPAQ and GLID diagnostics
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


MODULE MOD_COSP
  USE COSP_KINDS,                  ONLY: wp
  USE MOD_COSP_CONFIG,             ONLY: R_UNDEF,PARASOL_NREFL,LIDAR_NCAT,LIDAR_NTYPE, SR_BINS,& !OPAQ
                                         N_HYDRO,RTTOV_MAX_CHANNELS,numMISRHgtBins,      &
                                         DBZE_BINS,LIDAR_NTEMP,calipso_histBsct,         &
                                         use_vgrid,Nlvgrid,vgrid_zu,vgrid_zl,vgrid_z,    &
                                         numMODISTauBins,numMODISPresBins,               &
                                         numMODISReffIceBins,numMODISReffLiqBins,        &
                                         numISCCPTauBins,numISCCPPresBins,numMISRTauBins,&
                                         ntau,modis_histTau,tau_binBounds,               &
                                         modis_histTauEdges,tau_binEdges,                &
                                         modis_histTauCenters,tau_binCenters
  
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
  USE MOD_LIDAR_SIMULATOR,         ONLY: lidar_subcolumn,     lidar_column,   & !GLID
                                         lidar_subcolumn_gr,  lidar_column_gr   !GLID
  USE MOD_MODIS_SIM,               ONLY: modis_subcolumn,     modis_column
  USE MOD_PARASOL,                 ONLY: parasol_subcolumn,   parasol_column, ntetas
  use mod_cosp_rttov,              ONLY: rttov_column
  USE MOD_COSP_STATS,              ONLY: COSP_LIDAR_ONLY_CLOUD,COSP_CHANGE_VERTICAL_GRID
  
  IMPLICIT NONE
  
  logical :: linitialization ! Initialization flag
  
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

     real(wp),allocatable,dimension(:) :: &
          land,                & ! Land/Sea mask                          (0-1)
          skt,                 & ! Surface temperature                    (K)    !TIBO2
          surfelev               ! Surface Elevation                      (m)    !TIBO2
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
     real(wp),allocatable,dimension(:,:,:,:) :: &
          taupart,             & !GLID
          taupart_gr             !GLID
     real(wp),allocatable,dimension(:,:,:) :: &
          frac_out,            & ! Cloud fraction
          tau_067,             & ! Optical depth
          fracLiq,             & ! Cloud fraction
          emiss_11,            & ! Emissivity
          asym,                & ! Assymetry parameter
          ss_alb,              & ! Single-scattering albedo
          betatot,             & ! Backscatter coefficient for polarized optics (total)
          betatot_gr,          & !GLID
          betatot_ice,         & ! Backscatter coefficient for polarized optics (ice)
          betatot_liq,         & ! Backscatter coefficient for polarized optics (liquid)
          tautot,              & ! Optical thickess integrated from top (total)
          tautot_gr,           & !GLID
          tautot_ice,          & ! Optical thickess integrated from top (ice)
          tautot_liq,          & ! Optical thickess integrated from top (liquid)
          z_vol_cloudsat,      & ! Effective reflectivity factor (mm^6/m^3)
          kr_vol_cloudsat,     & ! Attenuation coefficient hydro (dB/km) 
          g_vol_cloudsat         ! Attenuation coefficient gases (dB/km)
     real(wp),allocatable,dimension(:,:) :: &
          beta_mol,            & ! Molecular backscatter coefficient
          beta_mol_gr,         & !GLID
          tau_mol,             & ! Molecular optical depth
          tau_mol_gr,          & !GLID
          tautot_S_liq,        & ! Liquid water optical thickness, from TOA to SFC
          tautot_S_ice           ! Ice water optical thickness, from TOA to SFC 
     type(radar_cfg) :: &
          rcfg_cloudsat         ! Radar comfiguration information (CLOUDSAT)
  end type cosp_optical_inputs
  
  ! ######################################################################################
  ! TYPE cosp_outputs
  ! ######################################################################################
  type cosp_outputs

     ! CALIPSO outputs
     real(wp),dimension(:,:,:),pointer :: &
          calipso_betaperp_tot,  & ! Total backscattered signal
          calipso_beta_tot,      & ! Total backscattered signal
          calipso_beta_tot_gr,   & ! Total backscattered signal                                 !GLID
          calipso_tau_tot,       & ! Optical thickness integrated from top to level z
          calipso_lidarcldphase, & ! 3D "lidar" phase cloud fraction 
          calipso_lidarcldtype,  & ! 3D "lidar" OPAQ type fraction, opaque and thin clouds + z_opaque and opacity fraction (Ntype+1=4) !OPAQ
          calipso_cldlayerphase, & ! low, mid, high-level lidar phase cloud cover
          calipso_lidarcldtmp,   & ! 3D "lidar" phase cloud temperature
          calipso_cfad_sr,       & ! CFAD of scattering ratio
          calipso_cfad_sr_gr       ! CFAD of GROUND LIDAR scattering ratio                      !GLID
     real(wp), dimension(:,:),pointer :: &
          calipso_lidarcld,      & ! 3D "lidar" cloud fraction 
          calipso_lidarcld_gr,   & ! 3D GROUND "lidar" cloud fraction                           !GLID
          calipso_cldlayer,      & ! low, mid, high-level, total lidar cloud cover
          calipso_cldlayer_gr,   & ! low, mid, high-level, total GROUND lidar cloud cover       !GLID
          calipso_cldtype,       & ! opaque and thin lidar cloud cover + z_opaque altitude      !OPAQ
          calipso_cldtypetemp,   & ! opaque and thin cloud temperature                          !TIBO
          calipso_cldtypemeanz,  & ! opaque and thin cloud altitude                             !TIBO
          calipso_cldtypemeanzse,& ! opaque, thin cloud and z_opaque altitude with respect to SE !TIBO2
          calipso_beta_mol,      & ! Molecular backscatter
          calipso_beta_mol_gr,   & ! Molecular backscatter                                      !GLID
          calipso_temp_tot
     real(wp), dimension(:),pointer :: &
          calipso_cldthinemis,   & ! thin cloud emissivity                                      !TIBO
          calipso_srbval           ! SR bins in cfad_sr
     
     ! PARASOL outputs
     real(wp),dimension(:,:,:),pointer :: &
          parasolPix_refl            ! PARASOL reflectances (subcolumn)    
     real(wp),dimension(:,:),pointer :: &
          parasolGrid_refl           ! PARASOOL reflectances (column)

     ! CLOUDSAT outputs
     real(wp),dimension(:,:,:),pointer :: &
          cloudsat_Ze_tot,         & ! Effective reflectivity factor (Npoints,Ncolumns,Nlevels)     
          cloudsat_cfad_ze           ! Ze CFAD(Npoints,dBZe_bins,Nlevels)
     real(wp), dimension(:,:),pointer :: &
          lidar_only_freq_cloud      ! (Npoints,Nlevels)
     real(wp),dimension(:),pointer :: &
          radar_lidar_tcc            ! Radar&lidar total cloud amount, grid-box scale (Npoints)
          
     ! ISCCP outputs       
     real(wp),dimension(:),pointer :: &
          isccp_totalcldarea, & ! The fraction of model grid box columns with cloud 
           				        ! somewhere in them. (%)
          isccp_meantb,       & ! Mean all-sky 10.5 micron brightness temperature. (K)
          isccp_meantbclr,    & ! Mean clear-sky 10.5 micron brightness temperature. (K)
          isccp_meanptop,     & ! Mean cloud top pressure (mb).
          isccp_meantaucld,   & ! Mean optical thickness. (1)
          isccp_meanalbedocld   ! Mean cloud albedo. (1)
     real(wp),dimension(:,:),pointer ::&
          isccp_boxtau,       & ! Optical thickness in each column. (1)
          isccp_boxptop         ! Cloud top pressure in each column. (mb)
     real(wp),dimension(:,:,:),pointer :: &
          isccp_fq              ! The fraction of the model grid box covered by each of
                                ! the 49 ISCCP D level cloud types. (%)
     
     ! MISR outptus  			    
     real(wp),dimension(:,:,:),pointer ::   & !
          misr_fq          ! Fraction of the model grid box covered by each of the MISR 
                           ! cloud types
     real(wp),dimension(:,:),pointer ::   & !
          misr_dist_model_layertops !  
     real(wp),dimension(:),pointer ::   & !
          misr_meanztop, & ! Mean MISR cloud top height
          misr_cldarea     ! Mean MISR cloud cover area         			    

     ! MODIS outptus		    
     real(wp),pointer,dimension(:) ::      & !  
          modis_Cloud_Fraction_Total_Mean,       & ! L3 MODIS retrieved cloud fraction (total) 
          modis_Cloud_Fraction_Water_Mean,       & ! L3 MODIS retrieved cloud fraction (liq) 
          modis_Cloud_Fraction_Ice_Mean,         & ! L3 MODIS retrieved cloud fraction (ice) 
          modis_Cloud_Fraction_High_Mean,        & ! L3 MODIS retrieved cloud fraction (high) 
          modis_Cloud_Fraction_Mid_Mean,         & ! L3 MODIS retrieved cloud fraction (middle) 
          modis_Cloud_Fraction_Low_Mean,         & ! L3 MODIS retrieved cloud fraction (low ) 
          modis_Optical_Thickness_Total_Mean,    & ! L3 MODIS retrieved optical thickness (tot)
          modis_Optical_Thickness_Water_Mean,    & ! L3 MODIS retrieved optical thickness (liq)
          modis_Optical_Thickness_Ice_Mean,      & ! L3 MODIS retrieved optical thickness (ice)
          modis_Optical_Thickness_Total_LogMean, & ! L3 MODIS retrieved log10 optical thickness 
          modis_Optical_Thickness_Water_LogMean, & ! L3 MODIS retrieved log10 optical thickness 
          modis_Optical_Thickness_Ice_LogMean,   & ! L3 MODIS retrieved log10 optical thickness
          modis_Cloud_Particle_Size_Water_Mean,  & ! L3 MODIS retrieved particle size (liquid)
          modis_Cloud_Particle_Size_Ice_Mean,    & ! L3 MODIS retrieved particle size (ice)
          modis_Cloud_Top_Pressure_Total_Mean,   & ! L3 MODIS retrieved cloud top pressure
          modis_Liquid_Water_Path_Mean,          & ! L3 MODIS retrieved liquid water path
          modis_Ice_Water_Path_Mean                ! L3 MODIS retrieved ice water path
     real(wp),pointer,dimension(:,:,:) ::  &
          modis_Optical_Thickness_vs_Cloud_Top_Pressure, & ! Tau/Pressure joint histogram          			    
          modis_Optical_Thickness_vs_ReffICE,            & ! Tau/ReffICE joint histogram
          modis_Optical_Thickness_vs_ReffLIQ               ! Tau/ReffLIQ joint histogram

     ! RTTOV outputs
     real(wp),pointer :: &
          rttov_tbs(:,:) ! Brightness Temperature	    
     
  end type cosp_outputs

CONTAINS
  ! ######################################################################################
  ! FUNCTION cosp_simulator
  ! ######################################################################################
!  SUBROUTINE COSP_SIMULATOR(cospIN,cospgridIN,cospOUT,start_idx,stop_idx)
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
         isim,t0,t1,i,icol,nloop,rmod,nprof,nlevels,istart,istop,maxlim,ij,ik,i1,i2,nError
    integer,target :: &
         Npoints
    logical :: &
         Lisccp_subcolumn,    & ! On/Off switch for subcolumn ISCCP simulator
         Lmisr_subcolumn,     & ! On/Off switch for subcolumn MISR simulator
         Lcalipso_subcolumn,  & ! On/Off switch for subcolumn CALIPSO simulator
         Lcalipso_subcolumn_gr,&! On/Off switch for subcolumn GROUND LIDAR simulator !GLID
         Lparasol_subcolumn,  & ! On/Off switch for subcolumn PARASOL simulator
         Lcloudsat_subcolumn, & ! On/Off switch for subcolumn CLOUDSAT simulator
         Lmodis_subcolumn,    & ! On/Off switch for subcolumn MODIS simulator
         Lrttov_subcolumn,    & ! On/Off switch for subcolumn RTTOV simulator
         Lisccp_column,       & ! On/Off switch for column ISCCP simulator
         Lmisr_column,        & ! On/Off switch for column MISR simulator
         Lcalipso_column,     & ! On/Off switch for column CALIPSO simulator
         Lcalipso_column_gr,  & ! On/Off switch for column GROUND LIDAR simulator !GLID
         Lparasol_column,     & ! On/Off switch for column PARASOL simulator
         Lcloudsat_column,    & ! On/Off switch for column CLOUDSAT simulator
         Lmodis_column,       & ! On/Off switch for column MODIS simulator
         Lrttov_column,       & ! On/Off switch for column RTTOV simulator (not used)      
         Lradar_lidar_tcc,    & ! On/Off switch from joint Calipso/Cloudsat product
         Llidar_only_freq_cloud  ! On/Off switch from joint Calipso/Cloudsat product
    logical :: &
         ok_lidar_cfad    = .false., &
         ok_lidar_cfad_gr = .false., & !GLID
         lrttov_cleanUp   = .false.
    
    integer, dimension(:,:),allocatable  :: &
         modisRetrievedPhase,isccpLEVMATCH
    real(wp), dimension(:),  allocatable  :: &
         modisCfTotal,modisCfLiquid,modisMeanIceWaterPath, isccp_meantbclr,     &                         
         modisCfIce, modisCfHigh, modisCfMid, modisCfLow,modisMeanTauTotal,     &       
         modisMeanTauLiquid, modisMeanTauIce, modisMeanLogTauTotal,             &       
         modisMeanLogTauLiquid, modisMeanLogTauIce, modisMeanSizeLiquid,        &        
         modisMeanSizeIce, modisMeanCloudTopPressure, modisMeanLiquidWaterPath, &
         radar_lidar_tcc
    REAL(WP), dimension(:,:),allocatable  :: &
         modisRetrievedCloudTopPressure,modisRetrievedTau,modisRetrievedSize,   &
         misr_boxtau,misr_boxztop,misr_dist_model_layertops,isccp_boxtau,       &
         isccp_boxttop,isccp_boxptop,calipso_beta_mol,lidar_only_freq_cloud,    &
         calipso_beta_mol_gr !GLID
    REAL(WP), dimension(:,:,:),allocatable :: &
         modisJointHistogram,modisJointHistogramIce,modisJointHistogramLiq,     &
         calipso_beta_tot,calipso_betaperp_tot, cloudsatDBZe,parasolPix_refl,   &
         calipso_beta_tot_gr !GLID
    real(wp),dimension(:),allocatable,target :: &
         out1D_1,out1D_2,out1D_3,out1D_4,out1D_5,out1D_6,out1D_7,out1D_8,       & !OPAQ
         out1D_9,out1D_10,out1D_11,out1D_12                                       !TIBO !TIBO2
    real(wp),dimension(:,:,:),allocatable :: &
       t_in,betamol_in,tmpFlip,betamolFlip,pnormFlip,pnorm_perpFlip,ze_totFlip
    real(wp),dimension(20) :: cosp_time

    call cpu_time(cosp_time(1))
    ! Initialize error reporting for output
    cosp_simulator(:)=''

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
    Lcalipso_subcolumn_gr = .false. !GLID
    Lparasol_subcolumn  = .false.
    Lcloudsat_subcolumn = .false.
    Lmodis_subcolumn    = .false.
    Lrttov_subcolumn    = .false.
    Lisccp_column       = .false.
    Lmisr_column        = .false.
    Lcalipso_column     = .false.
    Lcalipso_column_gr  = .false. !GLID
    Lparasol_column     = .false.
    Lcloudsat_column    = .false.
    Lmodis_column       = .false.
    Lrttov_column       = .false.
    Lradar_lidar_tcc    = .false.
    Llidar_only_freq_cloud = .false.

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

    ! GROUND LIDAR subcolumn !GLID
    if (associated(cospOUT%calipso_beta_mol_gr)                            .or.          & !GLID
        associated(cospOUT%calipso_beta_tot_gr))                                         & !GLID
       Lcalipso_subcolumn_gr  = .true. !GLID

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
        associated(cospOUT%calipso_lidarcldtype)                           .or.          & !OPAQ
        associated(cospOUT%calipso_cldlayer)                               .or.          &
        associated(cospOUT%calipso_cldtype)                                .or.          & !OPAQ
        associated(cospOUT%calipso_cldtypetemp)                            .or.          & !TIBO
        associated(cospOUT%calipso_cldtypemeanz)                           .or.          & !TIBO
        associated(cospOUT%calipso_cldtypemeanzse)                         .or.          & !TIBO2
        associated(cospOUT%calipso_cldthinemis)                            .or.          & !TIBO
        associated(cospOUT%calipso_cldlayerphase)                          .or.          &
        associated(cospOUT%calipso_lidarcldtmp)) then
       Lcalipso_column    = .true.
       Lcalipso_subcolumn = .true.
    endif

    ! GROUND LIDAR column !GLID
    if (associated(cospOUT%calipso_cfad_sr_gr)                             .or.          & !GLID
        associated(cospOUT%calipso_lidarcld_gr)                            .or.          & !GLID
        associated(cospOUT%calipso_cldlayer_gr)) then                                      !GLID
       Lcalipso_column_gr    = .true. !GLID
       Lcalipso_subcolumn_gr = .true. !GLID
    endif !GLID

    ! PARASOL column
    if (associated(cospOUT%parasolGrid_refl)) then
       Lparasol_column    = .true.
       Lparasol_subcolumn = .true.
    endif

    ! CLOUDSAT column
    if (associated(cospOUT%cloudsat_cfad_ze)) then
       Lcloudsat_column    = .true.
       Lcloudsat_subcolumn = .true.
    endif

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
    if (associated(cospOUT%lidar_only_freq_cloud) .or. associated(cospOUT%radar_lidar_tcc)) then
       Lcalipso_column     = .true.
       Lcalipso_subcolumn  = .true.
       Lcloudsat_column    = .true.
       Lcloudsat_subcolumn = .true.
       Lradar_lidar_tcc    = .true.
       Llidar_only_freq_cloud = .true.
    endif
    
    call cpu_time(cosp_time(2))
    if (debug) print*,'   Time to check outputs to see which simualtor to run:  ',cosp_time(2)-cosp_time(1)
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 2b) Error Checking
    !     Enforce bounds on input fields. If input field is out-of-bounds, report error 
    !     and turn off simulator
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    call cosp_errorCheck(cospgridIN,cospIN,Lisccp_subcolumn,Lisccp_column,               &
                         Lmisr_subcolumn,Lmisr_column,Lmodis_subcolumn,Lmodis_column,    &
                         Lcloudsat_subcolumn,Lcloudsat_column,Lcalipso_subcolumn,        & !GLID
                         Lcalipso_subcolumn_gr,Lcalipso_column,Lcalipso_column_gr,       & !GLID
                         Lrttov_subcolumn,Lrttov_column,                                 & !GLID
                         Lparasol_subcolumn,Lparasol_column,Lradar_lidar_tcc,            &
                         Llidar_only_freq_cloud,cospOUT,cosp_simulator,nError)
    call cpu_time(cosp_time(3))
    if (debug) print*,'   Time for cosp_errorCheck:                             ',cosp_time(3)-cosp_time(2)

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

    if (Lcalipso_subcolumn_gr) then !GLID
       calipsoIN%beta_mol_gr => cospIN%beta_mol_gr !GLID
       calipsoIN%betatot_gr  => cospIN%betatot_gr  !GLID
       calipsoIN%tau_mol_gr  => cospIN%tau_mol_gr  !GLID
       calipsoIN%tautot_gr   => cospIN%tautot_gr   !GLID
    endif !GLID
    
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
          modisIN%sunlit    = pack((/ (i, i = 1, Npoints ) /),mask = cospgridIN%sunlit > 0)
          modisIN%pres      = cospgridIN%phalf(int(modisIN%sunlit(:)),:)
       endif
       if (count(cospgridIN%sunlit <= 0) .gt. 0) then
          allocate(modisIN%notSunlit(count(cospgridIN%sunlit <= 0)))
          modisIN%notSunlit = pack((/ (i, i = 1, Npoints ) /),mask = .not. cospgridIN%sunlit > 0)
       endif
       !allocate(modisIN%sunlit(modisIN%Nsunlit),                                         &
       !         modisIN%notSunlit(count(cospgridIN%sunlit <= 0)),                        &       
       !         modisIN%pres(modisIN%Nsunlit,cospIN%Nlevels+1))             
       !modisIN%sunlit    = pack((/ (i, i = 1, Npoints ) /),                              &
       !     mask = cospgridIN%sunlit > 0)
       !modisIN%notSunlit = pack((/ (i, i = 1, Npoints ) /),                              &
       !     mask = .not. cospgridIN%sunlit > 0)
       !modisIN%pres      = cospgridIN%phalf(int(modisIN%sunlit(:)),:)
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
    call cpu_time(cosp_time(3))
    if (debug) print*,'   Time to populate simulator imputs:                    ',cosp_time(3)-cosp_time(2)

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
   call cpu_time(cosp_time(4))
   if (debug) print*,'   Time to run isccp_subcolumn:                          ',cosp_time(4)-cosp_time(3)

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
   call cpu_time(cosp_time(5))
   if (debug) print*,'   Time to run misr_subcolumn:                           ',cosp_time(5)-cosp_time(4)

    ! Calipso subcolumn simulator
    if (Lcalipso_subcolumn) then
       ! Allocate space for local variables
       allocate(calipso_beta_mol(calipsoIN%Npoints,calipsoIN%Nlevels),                   &
                calipso_beta_tot(calipsoIN%Npoints,calipsoIN%Ncolumns,calipsoIN%Nlevels),&
                calipso_betaperp_tot(calipsoIN%Npoints,calipsoIN%Ncolumns,calipsoIN%Nlevels))
       ! Call simulator
       call lidar_subcolumn(calipsoIN%npoints,calipsoIN%ncolumns,calipsoIN%nlevels,      &
                            calipsoIN%beta_mol,calipsoIN%tau_mol,                        &
                            calipsoIN%betatot,calipsoIN%tautot,calipsoIN%betatot_ice,    &
                            calipsoIN%tautot_ice,calipsoIN%betatot_liq,                  &
                            calipsoIN%tautot_liq,calipso_beta_mol(:,:),                  &
                            calipso_beta_tot(:,:,:),calipso_betaperp_tot(:,:,:))
       ! Store output (if requested)
       if (associated(cospOUT%calipso_beta_mol))                                         &
            cospOUT%calipso_beta_mol(ij:ik,calipsoIN%Nlevels:1:-1) = calipso_beta_mol
       if (associated(cospOUT%calipso_beta_tot))                                         &
            cospOUT%calipso_beta_tot(ij:ik,:,calipsoIN%Nlevels:1:-1) = calipso_beta_tot
       if (associated(cospOUT%calipso_betaperp_tot))                                     &
            cospOUT%calipso_betaperp_tot(ij:ik,:,:) = calipso_betaperp_tot

    endif

    ! GROUND LIDAR subcolumn simulator                                                       !GLID
    if (Lcalipso_subcolumn_gr) then                                                          !GLID
       ! Allocate space for local variables                                                  !GLID
       allocate(calipso_beta_mol_gr(calipsoIN%Npoints,calipsoIN%Nlevels),                &   !GLID
                calipso_beta_tot_gr(calipsoIN%Npoints,calipsoIN%Ncolumns,calipsoIN%Nlevels)) !GLID
       ! Call simulator                                                                      !GLID
       call lidar_subcolumn_gr(calipsoIN%npoints,calipsoIN%ncolumns,calipsoIN%nlevels,   &   !GLID
                               calipsoIN%beta_mol_gr,calipsoIN%tau_mol_gr,               &   !GLID
                               calipsoIN%betatot_gr,calipsoIN%tautot_gr,                 &   !GLID
                               calipso_beta_mol_gr(:,:),calipso_beta_tot_gr(:,:,:))          !GLID
       ! Store output (if requested)                                                         !GLID
       if (associated(cospOUT%calipso_beta_mol_gr))                                      &   !GLID
            cospOUT%calipso_beta_mol_gr(ij:ik,calipsoIN%Nlevels:1:-1) = calipso_beta_mol_gr  !GLID
    endif                                                                                    !GLID

   call cpu_time(cosp_time(6))
   if (debug) print*,'   Time to run lidar_subcolumn:                          ',cosp_time(6)-cosp_time(5)

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
    call cpu_time(cosp_time(7))
    if (debug) print*,'   Time to run parasol_subcolumn:                        ',cosp_time(7)-cosp_time(6)

    ! Cloudsat (quickbeam) subcolumn simulator
    if (Lcloudsat_subcolumn) then
       ! Allocate space for local variables
       allocate(cloudsatDBZe(cloudsatIN%Npoints,cloudsatIN%Ncolumns,cloudsatIN%Nlevels))
       do icol=1,cloudsatIN%ncolumns
          call quickbeam_subcolumn(cloudsatIN%rcfg,cloudsatIN%Npoints,cloudsatIN%Nlevels,&
                                   cloudsatIN%hgt_matrix/1000._wp,                       &
                                   cloudsatIN%z_vol(:,icol,:),                           &
                                   cloudsatIN%kr_vol(:,icol,:),                          &
                                   cloudsatIN%g_vol(:,1,:),cloudsatDBze(:,icol,:))
       enddo
       ! Store output (if requested)
       if (associated(cospOUT%cloudsat_Ze_tot)) then
          cospOUT%cloudsat_Ze_tot(ij:ik,:,:) = cloudsatDBZe(:,:,cloudsatIN%Nlevels:1:-1)
       endif
    endif
    call cpu_time(cosp_time(8))
    if (debug) print*,'   Time to run radar_subcolumn:                          ',cosp_time(8)-cosp_time(7)
   
    if (Lmodis_subcolumn) then
       if(modisiN%nSunlit > 0) then 
          ! Allocate space for local variables
          allocate(modisRetrievedTau(modisIN%nSunlit,modisIN%nColumns),                  &
                   modisRetrievedSize(modisIN%nSunlit,modisIN%nColumns),                 &
                   modisRetrievedPhase(modisIN%nSunlit,modisIN%nColumns),                &
                   modisRetrievedCloudTopPressure(modisIN%nSunlit,modisIN%nColumns))
          ! Call simulator
          do i = 1, modisIN%nSunlit
             call modis_subcolumn(modisIN%Ncolumns,modisIN%Nlevels,                      &
                                  modisIN%pres(int(modisIN%sunlit(i)),:),                &
                                  modisIN%tau(int(modisIN%sunlit(i)),:,:),               &
                                  modisIN%liqFrac(int(modisIN%sunlit(i)),:,:),           &
                                  modisIN%g(int(modisIN%sunlit(i)),:,:),                 &
                                  modisIN%w0(int(modisIN%sunlit(i)),:,:),                &
                                  isccp_boxtau(int(modisIN%sunlit(i)),:),                &
                                  isccp_boxptop(int(modisIN%sunlit(i)),:),               &
                                  modisRetrievedPhase(i,:),                              &
                                  modisRetrievedCloudTopPressure(i,:),                   &
                                  modisRetrievedTau(i,:),modisRetrievedSize(i,:))
          end do
       endif
    endif
    call cpu_time(cosp_time(9))
    if (debug) print*,'   Time to run modis_subcolum:                           ',cosp_time(9)-cosp_time(8)


    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 5) Call column simulators
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    ! ISCCP
    if (Lisccp_column) then
       ! Check to see which outputs are requested. If not requested, use a local dummy array
       if(.not. associated(cospOUT%isccp_meanalbedocld)) then
          allocate(out1D_1(Npoints))
          cospOUT%isccp_meanalbedocld(ij:ik) => out1D_1
       endif
       if(.not. associated(cospOUT%isccp_meanptop)) then
          allocate(out1D_2(Npoints))
          cospOUT%isccp_meanptop(ij:ik) => out1D_2
       endif
       if(.not. associated(cospOUT%isccp_meantaucld)) then
          allocate(out1D_3(Npoints))
          cospOUT%isccp_meantaucld(ij:ik) => out1D_3
       endif   
       if(.not. associated(cospOUT%isccp_totalcldarea)) then
          allocate(out1D_4(Npoints))
          cospOUT%isccp_totalcldarea(ij:ik) => out1D_4
       endif
       if(.not. associated(cospOUT%isccp_meantb)) then
          allocate(out1D_5(Npoints))
          cospOUT%isccp_meantb(ij:ik) => out1D_5    
       endif
       if(.not. associated(cospOUT%isccp_fq)) then
          allocate(out1D_6(Npoints*numISCCPTauBins*numISCCPPresBins))       
          cospOUT%isccp_fq(ij:ik,1:numISCCPTauBins,1:numISCCPPresBins) => out1D_6
       endif   
                                
       ! Call simulator
       call icarus_column(isccpIN%npoints, isccpIN%ncolumns, isccpIN%nlevels,            &
                          isccp_boxtau(:,:),isccp_boxptop(:,:)/100._wp,                  &
                          isccpIN%sunlit,isccpIN%pfull,isccpIN%phalf,isccpIN%qv,         &
                          isccpIN%at,isccpIN%skt,isccpIN%emsfc_lw,isccp_boxttop,         &
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
       if (allocated(isccp_boxptop))   deallocate(isccp_boxptop)
       if (allocated(isccp_boxtau))    deallocate(isccp_boxtau)
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
       if (allocated(out1D_6)) then
          deallocate(out1D_6)
          nullify(cospOUT%isccp_fq)
       endif
    endif
    call cpu_time(cosp_time(10))
    if (debug) print*,'   Time to run isccp_column:                             ',cosp_time(10)-cosp_time(9)

    ! MISR
    if (Lmisr_column) then
       ! Check to see which outputs are requested. If not requested, use a local dummy array
       if (.not. associated(cospOUT%misr_cldarea)) then
          allocate(out1D_1(Npoints))
          cospOUT%misr_cldarea(ij:ik) => out1D_1                 
       endif
       if (.not. associated(cospOUT%misr_meanztop)) then 
          allocate(out1D_2(Npoints))
          cospOUT%misr_meanztop(ij:ik) => out1D_2
       endif
       if (.not. associated(cospOUT%misr_fq)) then
          allocate(out1D_3(Npoints*numMISRTauBins*numMISRHgtBins))
          cospOUT%misr_fq(ij:ik,1:numMISRTauBins,1:numMISRHgtBins) => out1D_3     
        endif   
    
       ! Call simulator
       call misr_column(misrIN%Npoints,misrIN%Ncolumns,misrIN%Nlevels,misr_boxztop,     &
                        misrIN%sunlit,misr_boxtau,cospOUT%misr_cldarea(ij:ik),          &
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
       if (allocated(out1D_3)) then
          deallocate(out1D_3)
          nullify(cospOUT%misr_fq)
       endif
    endif
    call cpu_time(cosp_time(11))
    if (debug) print*,'   Time to run misr_column:                              ',cosp_time(11)-cosp_time(10)
    
    ! CALIPSO LIDAR Simulator
    if (Lcalipso_column) then
       ! Check to see which outputs are requested. If not requested, use a local dummy array
       if (.not. associated(cospOUT%calipso_cfad_sr)) then
          allocate(out1D_1(Npoints*SR_BINS*Nlvgrid))
          cospOUT%calipso_cfad_sr(ij:ik,1:SR_BINS,1:Nlvgrid) => out1D_1
       endif
       if (.not. associated(cospOUT%calipso_lidarcld)) then
          allocate(out1D_2(Npoints*Nlvgrid))
          cospOUT%calipso_lidarcld(ij:ik,1:Nlvgrid) => out1D_2
       endif   
       if (.not. associated(cospOUT%calipso_lidarcldphase)) then
          allocate(out1D_3(Npoints*Nlvgrid*6))
          cospOUT%calipso_lidarcldphase(ij:ik,1:Nlvgrid,1:6) => out1D_3
       endif
       if (.not. associated(cospOUT%calipso_cldlayer)) then
          allocate(out1D_4(Npoints*LIDAR_NCAT))
          cospOUT%calipso_cldlayer(ij:ik,1:LIDAR_NCAT) => out1D_4
       endif
       if (.not. associated(cospOUT%calipso_cldlayerphase)) then
          allocate(out1D_5(Npoints*LIDAR_NCAT*6))
          cospOUT%calipso_cldlayerphase(ij:ik,1:LIDAR_NCAT,1:6) => out1D_5
       endif   
       if (.not. associated(cospOUT%calipso_lidarcldtmp)) then
          allocate(out1D_6(Npoints*40*5))
          cospOUT%calipso_lidarcldtmp(ij:ik,1:40,1:5) => out1D_6
       endif   
       if (.not. associated(cospOUT%calipso_lidarcldtype)) then        !OPAQ
          allocate(out1D_7(Npoints*Nlvgrid*4))                         !OPAQ
          cospOUT%calipso_lidarcldtype(ij:ik,1:Nlvgrid,1:4) => out1D_7 !OPAQ
       endif                                                           !OPAQ
       if (.not. associated(cospOUT%calipso_cldtype)) then             !OPAQ
          allocate(out1D_8(Npoints*LIDAR_NTYPE))                       !OPAQ
          cospOUT%calipso_cldtype(ij:ik,1:LIDAR_NTYPE) => out1D_8      !OPAQ
       endif                                                           !OPAQ
       if (.not. associated(cospOUT%calipso_cldtypetemp)) then         !TIBO
          allocate(out1D_9(Npoints*LIDAR_NTYPE))                       !TIBO
          cospOUT%calipso_cldtypetemp(ij:ik,1:LIDAR_NTYPE) => out1D_9  !TIBO
       endif                                                           !TIBO
       if (.not. associated(cospOUT%calipso_cldtypemeanz)) then        !TIBO
          allocate(out1D_10(Npoints*2))                                !TIBO
          cospOUT%calipso_cldtypemeanz(ij:ik,1:2) => out1D_10          !TIBO
       endif                                                           !TIBO
       if (.not. associated(cospOUT%calipso_cldtypemeanzse)) then      !TIBO2
          allocate(out1D_12(Npoints*3))                                !TIBO2
          cospOUT%calipso_cldtypemeanzse(ij:ik,1:3) => out1D_12        !TIBO2
       endif                                                           !TIBO2
       if (.not. associated(cospOUT%calipso_cldthinemis)) then         !TIBO
          allocate(out1D_11(Npoints))                                  !TIBO
          cospOUT%calipso_cldthinemis(ij:ik) => out1D_11               !TIBO
       endif                                                           !TIBO
       
       ! Call simulator
       ok_lidar_cfad=.true.
       call lidar_column(calipsoIN%Npoints,calipsoIN%Ncolumns,calipsoIN%Nlevels,         &
                         Nlvgrid,SR_BINS,cospgridIN%at(:,:),                             &
                         calipso_beta_tot(:,:,:),calipso_betaperp_tot(:,:,:),            &
                         calipso_beta_mol(:,:),cospgridIN%land,cospgridIN%surfelev,      & !TIBO2
                         cospgridIN%phalf(:,2:calipsoIN%Nlevels),ok_lidar_cfad,          &
                         LIDAR_NCAT,LIDAR_NTYPE,cospOUT%calipso_cfad_sr(ij:ik,:,:),      & !OPAQ
                         cospOUT%calipso_lidarcld(ij:ik,:),                              &
                         cospOUT%calipso_lidarcldphase(ij:ik,:,:),                       &
                         cospOUT%calipso_lidarcldtype(ij:ik,:,:),                        & !OPAQ
                         cospOUT%calipso_cldlayer(ij:ik,:),                              &
                         cospOUT%calipso_cldtype(ij:ik,:),                               & !OPAQ
                         cospOUT%calipso_cldtypetemp(ij:ik,:),                           & !TIBO
                         cospOUT%calipso_cldtypemeanz(ij:ik,:),                          & !TIBO
                         cospOUT%calipso_cldtypemeanzse(ij:ik,:),                        & !TIBO2
                         cospOUT%calipso_cldthinemis(ij:ik),                             & !TIBO
                         cospgridIN%hgt_matrix,cospgridIN%hgt_matrix_half,               &
                         cospOUT%calipso_cldlayerphase(ij:ik,:,:),                       &
                         cospOUT%calipso_lidarcldtmp(ij:ik,:,:),                         & !OPAQ
                         vgrid_z(:))                                                       !OPAQ

       if (associated(cospOUT%calipso_srbval)) cospOUT%calipso_srbval = calipso_histBsct

       ! Free up memory (if necessary)
       if (allocated(out1D_1)) then
          deallocate(out1D_1)
          nullify(cospOUT%calipso_cfad_sr)
       endif
       if (allocated(out1D_2)) then
          deallocate(out1D_2)
          nullify(cospOUT%calipso_lidarcld)
       endif
       if (allocated(out1D_3)) then
          deallocate(out1D_3)
          nullify(cospOUT%calipso_lidarcldphase)
       endif
       if (allocated(out1D_4)) then
          deallocate(out1D_4)
          nullify(cospOUT%calipso_cldlayer)
       endif
       if (allocated(out1D_5)) then
          deallocate(out1D_5)
          nullify(cospOUT%calipso_cldlayerphase)
       endif
       if (allocated(out1D_6)) then
          deallocate(out1D_6)
          nullify(cospOUT%calipso_lidarcldtmp)
       endif
       if (allocated(out1D_7)) then             !OPAQ
          deallocate(out1D_7)                   !OPAQ
          nullify(cospOUT%calipso_lidarcldtype) !OPAQ
       endif                                    !OPAQ
       if (allocated(out1D_8)) then             !OPAQ
          deallocate(out1D_8)                   !OPAQ
          nullify(cospOUT%calipso_cldtype)      !OPAQ
       endif                                    !OPAQ
       if (allocated(out1D_9)) then             !TIBO
          deallocate(out1D_9)                   !TIBO
          nullify(cospOUT%calipso_cldtypetemp)  !TIBO
       endif                                    !TIBO
       if (allocated(out1D_10)) then            !TIBO
          deallocate(out1D_10)                  !TIBO
          nullify(cospOUT%calipso_cldtypemeanz) !TIBO
       endif                                    !TIBO
       if (allocated(out1D_12)) then              !TIBO2
          deallocate(out1D_12)                    !TIBO2
          nullify(cospOUT%calipso_cldtypemeanzse) !TIBO2
       endif                                      !TIBO2
       if (allocated(out1D_11)) then            !TIBO
          deallocate(out1D_11)                  !TIBO
          nullify(cospOUT%calipso_cldthinemis)  !TIBO
       endif                                    !TIBO

    endif

    ! GROUND LIDAR Simulator     !GLID
    if (Lcalipso_column_gr) then !GLID
       ! Check to see which outputs are requested. If not requested, use a local dummy array !GLID
       if (.not. associated(cospOUT%calipso_cfad_sr_gr)) then              !GLID  
          allocate(out1D_1(Npoints*SR_BINS*Nlvgrid))                       !GLID
          cospOUT%calipso_cfad_sr_gr(ij:ik,1:SR_BINS,1:Nlvgrid) => out1D_1 !GLID
       endif                                                               !GLID
       if (.not. associated(cospOUT%calipso_lidarcld_gr)) then    !GLID
          allocate(out1D_2(Npoints*Nlvgrid))                      !GLID
          cospOUT%calipso_lidarcld_gr(ij:ik,1:Nlvgrid) => out1D_2 !GLID
       endif                                                      !GLID
       if (.not. associated(cospOUT%calipso_cldlayer_gr)) then       !GLID
          allocate(out1D_3(Npoints*LIDAR_NCAT))                      !GLID
          cospOUT%calipso_cldlayer_gr(ij:ik,1:LIDAR_NCAT) => out1D_3 !GLID
       endif                                                         !GLID
       
       ! Call simulator        !GLID
       ok_lidar_cfad_gr=.true. !GLID
       call lidar_column_gr(calipsoIN%Npoints,calipsoIN%Ncolumns,calipsoIN%Nlevels,     & !GLID
                            Nlvgrid,SR_BINS,calipso_beta_tot_gr(:,:,:),                 & !GLID
                            calipso_beta_mol_gr(:,:),cospgridIN%land,                   & !GLID
                            cospgridIN%phalf(:,2:calipsoIN%Nlevels),                    & !GLID
                            ok_lidar_cfad_gr,LIDAR_NCAT,                                & !GLID
                            cospOUT%calipso_cfad_sr_gr(ij:ik,:,:),                      & !GLID
                            cospOUT%calipso_lidarcld_gr(ij:ik,:),                       & !GLID
                            cospOUT%calipso_cldlayer_gr(ij:ik,:),                       & !GLID
                            cospgridIN%hgt_matrix,cospgridIN%hgt_matrix_half)             !GLID

       ! Free up memory (if necessary)         !GLID
       if (allocated(out1D_1)) then            !GLID
          deallocate(out1D_1)                  !GLID
          nullify(cospOUT%calipso_cfad_sr_gr)  !GLID
       endif                                   !GLID
       if (allocated(out1D_2)) then            !GLID
          deallocate(out1D_2)                  !GLID
          nullify(cospOUT%calipso_lidarcld_gr) !GLID
       endif                                   !GLID
       if (allocated(out1D_3)) then            !GLID
          deallocate(out1D_3)                  !GLID
          nullify(cospOUT%calipso_cldlayer_gr) !GLID
       endif                                   !GLID

    endif !GLID

    call cpu_time(cosp_time(12))
    if (debug) print*,'   Time to run lidar_column:                             ',cosp_time(12)-cosp_time(11)

    ! PARASOL
    if (Lparasol_column) then
       call parasol_column(parasolIN%Npoints,PARASOL_NREFL,parasolIN%Ncolumns,           &
                            cospgridIN%land(:),parasolPix_refl(:,:,:),                   &
                            cospOUT%parasolGrid_refl(ij:ik,:))
       if (allocated(parasolPix_refl)) deallocate(parasolPix_refl)
    endif
    call cpu_time(cosp_time(13))
    if (debug) print*,'   Time to run parasol_column:                           ',cosp_time(13)-cosp_time(12)
    
    ! CLOUDSAT
    if (Lcloudsat_column) then
       ! Check to see which outputs are requested. If not requested, use a local dummy array
       if (.not. associated(cospOUT%cloudsat_cfad_ze)) then
          allocate(out1D_1(Npoints*DBZE_BINS*Nlvgrid))
          cospOUT%cloudsat_cfad_ze(ij:ik,1:DBZE_BINS,1:Nlvgrid) => out1D_1
       endif

       ! Call simulator
       call quickbeam_column(cloudsatIN%Npoints,cloudsatIN%Ncolumns,cloudsatIN%Nlevels,  &
                             Nlvgrid,cloudsatDBZe,cospgridIN%hgt_matrix,                 &
                             cospgridIN%hgt_matrix_half,cospOUT%cloudsat_cfad_ze(ij:ik,:,:))
       ! Free up memory  (if necessary)
       if (allocated(out1D_1)) then
          deallocate(out1D_1)
          nullify(cospOUT%cloudsat_cfad_ze)
       endif
    endif
    call cpu_time(cosp_time(14))
    if (debug) print*,'   Time to run radar_column:                             ',cosp_time(14)-cosp_time(13)

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
    call cpu_time(cosp_time(15))
    if (debug) print*,'   Time to run modis_column:                             ',cosp_time(15)-cosp_time(14)
    
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
    call cpu_time(cosp_time(16))
    if (debug) print*,'   Time to run rttov_column:                             ',cosp_time(16)-cosp_time(15)
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 6) Compute multi-instrument products
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! CLOUDSAT/CALIPSO products
    if (Lradar_lidar_tcc .or. Llidar_only_freq_cloud) then
      
       if (use_vgrid) then
          allocate(lidar_only_freq_cloud(cloudsatIN%Npoints,Nlvgrid),                    &
               radar_lidar_tcc(cloudsatIN%Npoints))
          allocate(t_in(cloudsatIN%Npoints,1,cloudsatIN%Nlevels),                        &
                   betamol_in(cloudsatIN%Npoints,1,cloudsatIN%Nlevels),                  &
                   tmpFlip(cloudsatIN%Npoints,1,Nlvgrid),                                &
                   betamolFlip(cloudsatIN%Npoints,1,Nlvgrid),                            &
                   pnormFlip(cloudsatIN%Npoints,cloudsatIN%Ncolumns,Nlvgrid),            &
                   pnorm_perpFlip(cloudsatIN%Npoints,cloudsatIN%Ncolumns,Nlvgrid),       &
                   Ze_totFlip(cloudsatIN%Npoints,cloudsatIN%Ncolumns,Nlvgrid))

          t_in(:,1,:)=cospgridIN%at(:,cloudsatIN%Nlevels:1:-1)
          call cosp_change_vertical_grid(cloudsatIN%Npoints,1,cloudsatIN%Nlevels,        &
               cospgridIN%hgt_matrix(:,cloudsatIN%Nlevels:1:-1),                         &
               cospgridIN%hgt_matrix_half(:,cloudsatIN%Nlevels:1:-1),t_in,Nlvgrid,       &
               vgrid_zl(Nlvgrid:1:-1),vgrid_zu(Nlvgrid:1:-1),tmpFlip(:,1,Nlvgrid:1:-1))
          
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
               calipso_betaperp_tot(:,:,cloudsatIN%Nlevels:1:-1),Nlvgrid,                &
               vgrid_zl(Nlvgrid:1:-1),vgrid_zu(Nlvgrid:1:-1),                            &
               pnorm_perpFlip(:,:,Nlvgrid:1:-1))
          
          call cosp_change_vertical_grid(cloudsatIN%Npoints,cloudsatIN%Ncolumns,         &
               cloudsatIN%Nlevels,cospgridIN%hgt_matrix(:,cloudsatIN%Nlevels:1:-1),      &
               cospgridIN%hgt_matrix_half(:,cloudsatIN%Nlevels:1:-1),                    &
               cloudsatDBZe(:,:,cloudsatIN%Nlevels:1:-1),Nlvgrid,vgrid_zl(Nlvgrid:1:-1), &
               vgrid_zu(Nlvgrid:1:-1),Ze_totFlip(:,:,Nlvgrid:1:-1),log_units=.true.)    

          call cosp_lidar_only_cloud(cloudsatIN%Npoints,cloudsatIN%Ncolumns,             &
                                     Nlvgrid,tmpFlip,pnormFlip,pnorm_perpFlip,           &
                                     betamolFlip,Ze_totFlip,                             &
                                     lidar_only_freq_cloud,radar_lidar_tcc)                                            
          deallocate(t_in,betamol_in,tmpFlip,betamolFlip,pnormFlip,pnorm_perpFlip,       &
                     ze_totFlip)
       else
          allocate(lidar_only_freq_cloud(cloudsatIN%Npoints,cloudsatIN%Nlevels),         &
               radar_lidar_tcc(cloudsatIN%Npoints))
          call cosp_lidar_only_cloud(cloudsatIN%Npoints,cloudsatIN%Ncolumns,             &
               cospIN%Nlevels,cospgridIN%at(:,cloudsatIN%Nlevels:1:-1),                  &
               calipso_beta_tot(:,:,cloudsatIN%Nlevels:1:-1),                            &
               calipso_betaperp_tot(:,:,cloudsatIN%Nlevels:1:-1),                        &
               calipso_beta_mol(:,cloudsatIN%Nlevels:1:-1),                              &
               cloudsatDBZe(:,:,cloudsatIN%Nlevels:1:-1),lidar_only_freq_cloud,          &
               radar_lidar_tcc)
       endif
       
       ! Store, when necessary
       if (associated(cospOUT%lidar_only_freq_cloud)) then
          cospOUT%lidar_only_freq_cloud(ij:ik,:) = lidar_only_freq_cloud
       endif
       if (associated(cospOUT%radar_lidar_tcc)) then
          cospOUT%radar_lidar_tcc(ij:ik) = radar_lidar_tcc
       endif

    endif
    call cpu_time(cosp_time(17))
    if (debug) print*,'   Time to create joint products:                        ',cosp_time(17)-cosp_time(16)
    
    ! Cleanup 
    if (allocated(calipso_beta_tot))      deallocate(calipso_beta_tot)
    if (allocated(calipso_beta_tot_gr))   deallocate(calipso_beta_tot_gr) !GLID
    if (allocated(calipso_beta_mol))      deallocate(calipso_beta_mol)
    if (allocated(calipso_beta_mol_gr))   deallocate(calipso_beta_mol_gr) !GLID
    if (allocated(calipso_betaperp_tot))  deallocate(calipso_betaperp_tot)
    if (allocated(cloudsatDBZe))          deallocate(cloudsatDBZe)
    if (allocated(lidar_only_freq_cloud)) deallocate(lidar_only_freq_cloud)
    if (allocated(radar_lidar_tcc))       deallocate(radar_lidar_tcc)

  end function COSP_SIMULATOR
  ! ######################################################################################
  ! SUBROUTINE cosp_init
  ! ######################################################################################
  SUBROUTINE COSP_INIT(Lisccp,Lmodis,Lmisr,Lcloudsat,Lcalipso,Lparasol,Lrttov,            &
                       Npoints,Nlevels,cloudsat_radar_freq,cloudsat_k2,                  &
                       cloudsat_use_gas_abs,cloudsat_do_ray,isccp_top_height,            &
                       isccp_top_height_direction,surface_radar,rcfg,rttov_Nchannels,    &
                       rttov_Channels,rttov_platform,rttov_satellite,rttov_instrument,   &
                       lusevgrid,luseCSATvgrid,Nvgrid,vgrid_z,cloudsat_micro_scheme,cospOUT) !OPAQ
    
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
    character(len=64),intent(in) :: &
       cloudsat_micro_scheme           ! Microphysical scheme used by CLOUDSAT
    type(cosp_outputs),intent(inout) :: cospOUT
    real(wp),intent(inout),dimension(:),allocatable :: &                         !OPAQ
         vgrid_z                       ! mid-level altitude of the vertical grid !OPAQ

    ! OUTPUTS
    type(radar_cfg) :: rcfg
  
    ! Local variables
    integer  :: i
    real(wp) :: zstep

    ! Initialize MODIS optical-depth bin boundaries for joint-histogram. (defined in cosp_config.F90)
    if (.not. allocated(modis_histTau)) then
       allocate(modis_histTau(ntau+1),modis_histTauEdges(2,ntau),modis_histTauCenters(ntau))
       numMODISTauBins      = ntau
       modis_histTau        = tau_binBounds
       modis_histTauEdges   = tau_binEdges
       modis_histTauCenters = tau_binCenters
    endif
    
    ! Set up vertical grid used by CALIPSO and CLOUDSAT L3
    use_vgrid = lusevgrid
    
    if (use_vgrid) then
      Nlvgrid  = Nvgrid
       allocate(vgrid_zl(Nlvgrid),vgrid_zu(Nlvgrid),vgrid_z(Nlvgrid))
       ! CloudSat grid requested
       if (luseCSATvgrid)       zstep = 480._wp
       ! Other grid requested. Constant vertical spacing with top at 20 km
       if (.not. luseCSATvgrid) zstep = 20000._wp/Nvgrid
       do i=1,Nvgrid
          vgrid_zl(Nlvgrid-i+1) = (i-1)*zstep
          vgrid_zu(Nlvgrid-i+1) = i*zstep
       enddo
       vgrid_z = (vgrid_zl+vgrid_zu)/2._wp
    else
       Nlvgrid = Nlevels
       allocate(vgrid_zl(Nlvgrid),vgrid_zu(Nlvgrid),vgrid_z(Nlvgrid))
    endif

    ! Initialize simulators
    if (Lisccp) call cosp_isccp_init(isccp_top_height,isccp_top_height_direction)
    if (Lmodis) call cosp_modis_init()
    if (Lmisr)  call cosp_misr_init()
    !if (Lrttov) call cosp_rttov_init(rttov_Nchannels,rttov_platform,rttov_satellite,     &
    !     rttov_instrument,rttov_channels)
    if (Lrttov) call cosp_rttov_init()
    if (Lcloudsat) call cosp_cloudsat_init(cloudsat_radar_freq,cloudsat_k2,              &
         cloudsat_use_gas_abs,cloudsat_do_ray,R_UNDEF,N_HYDRO,Npoints,Nlevels,           &
         surface_radar,rcfg,cloudsat_micro_scheme)
    if (Lcalipso) call cosp_calipso_init()
    if (Lparasol) call cosp_parasol_init()

    ! Set all output diagnostics as disassociated.
    nullify(cospOUT%calipso_betaperp_tot,cospOUT%calipso_beta_tot,                       &
            cospOUT%calipso_tau_tot,cospOUT%calipso_lidarcldphase,                       &
            cospOUT%calipso_cldlayerphase,cospOUT%calipso_lidarcldtmp,                   &
            cospOUT%calipso_cfad_sr,cospOUT%calipso_lidarcld,cospOUT%calipso_cldlayer,   &
            cospOUT%calipso_cfad_sr_gr,cospOUT%calipso_lidarcld_gr,                      & !GLID
            cospOUT%calipso_lidarcldtype,cospOUT%calipso_cldtype,                        & !OPAQ
            cospOUT%calipso_cldtypetemp,cospOUT%calipso_cldtypemeanz,                    & !TIBO
            cospOUT%calipso_cldtypemeanzse,cospOUT%calipso_cldthinemis,                  & !TIBO !TIBO2
            cospOUT%calipso_beta_mol,cospOUT%calipso_temp_tot,cospOUT%calipso_srbval,    &
            cospOUT%calipso_beta_mol_gr,cospOUT%calipso_cldlayer_gr,                     & !GLID
            cospOUT%parasolPix_refl,cospOUT%parasolGrid_refl,cospOUT%cloudsat_Ze_tot,    &   
            cospOUT%cloudsat_cfad_ze,cospOUT%lidar_only_freq_cloud,                      &
            cospOUT%radar_lidar_tcc,cospOUT%isccp_totalcldarea,cospOUT%isccp_meantb,     &
            cospOUT%isccp_meantbclr,cospOUT%isccp_meanptop,cospOUT%isccp_meantaucld,     &
            cospOUT%isccp_meanalbedocld,cospOUT%isccp_boxtau,cospOUT%isccp_boxptop,      &
            cospOUT%isccp_fq,cospOUT%misr_fq,cospOUT%misr_dist_model_layertops,          &
            cospOUT%misr_meanztop,cospOUT%misr_cldarea,                                  &
            cospOUT%modis_Cloud_Fraction_Total_Mean,                                     &
            cospOUT%modis_Cloud_Fraction_Water_Mean,                                     &
            cospOUT%modis_Cloud_Fraction_Ice_Mean,                                       &
            cospOUT%modis_Cloud_Fraction_High_Mean,                                      &
            cospOUT%modis_Cloud_Fraction_Mid_Mean,                                       &
            cospOUT%modis_Cloud_Fraction_Low_Mean,                                       &
            cospOUT%modis_Optical_Thickness_Total_Mean,                                  &
            cospOUT%modis_Optical_Thickness_Water_Mean,                                  &
            cospOUT%modis_Optical_Thickness_Ice_Mean,                                    &
            cospOUT%modis_Optical_Thickness_Total_LogMean,                               &
            cospOUT%modis_Optical_Thickness_Water_LogMean,                               &
            cospOUT%modis_Optical_Thickness_Ice_LogMean,                                 &
            cospOUT%modis_Cloud_Particle_Size_Water_Mean,                                &
            cospOUT%modis_Cloud_Particle_Size_Ice_Mean,                                  &
            cospOUT%modis_Cloud_Top_Pressure_Total_Mean,                                 &
            cospOUT%modis_Liquid_Water_Path_Mean,cospOUT%modis_Ice_Water_Path_Mean,      &
            cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure,                       &	    
            cospOUT%modis_Optical_Thickness_vs_ReffICE,                                  &
            cospOUT%modis_Optical_Thickness_vs_ReffLIQ,cospOUT%rttov_tbs)	    
     
    
    linitialization = .FALSE.
  END SUBROUTINE COSP_INIT
  
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
    
    allocate(y%tau_067(npoints,ncolumns,nlevels),    y%emiss_11(npoints,ncolumns,nlevels),&
             y%frac_out(npoints,ncolumns,nlevels),   y%beta_mol(npoints,nlevels),        &
             y%tau_mol(npoints,nlevels),             y%betatot(npoints,ncolumns,nlevels),&
             y%betatot_ice(npoints,ncolumns,nlevels),y%fracLiq(npoints,nColumns,nlevels),&
             y%betatot_liq(npoints,ncolumns,nlevels),y%tautot(npoints,ncolumns,nlevels), &
             y%tautot_ice(npoints,ncolumns,nlevels), y%tautot_S_ice(npoints,nlevels),    &
             y%tautot_liq(npoints,ncolumns,nlevels), y%tautot_S_liq(npoints,nlevels),    &
             y%taupart(npoints,ncolumns,nlevels,4),                                      &
             y%z_vol_cloudsat(npoints,Ncolumns,nlevels),                                 &
             y%kr_vol_cloudsat(npoints,Ncolumns,nlevels),                                &
             y%g_vol_cloudsat(npoints,Ncolumns,nlevels),                                 &
             y%asym(npoints,nColumns,nlevels),   y%ss_alb(npoints,nColumns,nlevels), & !GLID
             y%beta_mol_gr(npoints,nlevels), y%betatot_gr(npoints,ncolumns,nlevels), & !GLID
             y%taupart_gr(npoints,ncolumns,nlevels,4), y%tau_mol_gr(npoints,nlevels),& !GLID
             y%tautot_gr(npoints,ncolumns,nlevels)) !GLID
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
             y%cloudIce(nPoints,nLevels),y%cloudLiq(nPoints,nLevels),y%surfelev(npoints),& !TIBO2
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
                                    LlidarBetaMol532gr,LcfadLidarsr532,                  & !GLID
                                    LcfadLidarsr532gr,Lclcalipso2,                       & !GLID
                                    Lclcalipso,Lclcalipsogr,Lclhcalipso,Lcllcalipso,     & !GLID
                                    Lclmcalipso,Lcltcalipso,Lclhcalipsogr,               & !GLID
                                    Lcllcalipsogr,Lclmcalipsogr,Lcltcalipsogr,           & !GLID
                                    Lcltlidarradar,Lclcalipsoliq,                        & !GLID
                                    Lclcalipsoice,Lclcalipsoun,Lclcalipsotmp,            &
                                    Lclcalipsotmpliq,Lclcalipsotmpice,Lclcalipsotmpun,   &
                                    Lcltcalipsoliq,Lcltcalipsoice,Lcltcalipsoun,         &
                                    Lclhcalipsoliq,Lclhcalipsoice,Lclhcalipsoun,         &
                                    Lclmcalipsoliq,Lclmcalipsoice,Lclmcalipsoun,         &
                                    Lcllcalipsoliq,Lcllcalipsoice,Lcllcalipsoun,         & 
                                    Lclopaquecalipso,Lclthincalipso,Lclzopaquecalipso,   & !OPAQ
                                    Lclcalipsoopaque,Lclcalipsothin,Lclcalipsozopaque,   & !OPAQ
                                    Lclcalipsoopacity,Lclopaquetemp,Lclthintemp,         & !TIBO
                                    Lclzopaquetemp,Lclopaquemeanz,Lclthinmeanz,          & !TIBO
                                    Lclthinemis,Lclopaquemeanzse,Lclthinmeanzse,         & !TIBO2
                                    Lclzopaquecalipsose,LcfadDbze94,Ldbze94,Lparasolrefl,& !TIBO !TIBO2
                                    Ltbrttov,Npoints,Ncolumns,Nlevels,Nlvgrid,Nchan,x)
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
         LlidarBetaMol532gr,& ! GROUND LIDAR molecular backscatter (532nm)    !GLID    
         LcfadLidarsr532,  & ! CALIPSO scattering ratio CFAD
         LcfadLidarsr532gr, & ! GROUND LIDAR scattering ratio CFAD            !GLID
         Lclcalipso2,      & ! CALIPSO cloud fraction undetected by cloudsat
         Lclcalipso,       & ! CALIPSO cloud area fraction
         Lclcalipsogr,     & ! GROUND LIDAR cloud area fraction                !GLID
         Lclhcalipso,      & ! CALIPSO high-level cloud fraction
         Lcllcalipso,      & ! CALIPSO low-level cloud fraction
         Lclmcalipso,      & ! CALIPSO mid-level cloud fraction
         Lcltcalipso,      & ! CALIPSO total cloud fraction
         Lclhcalipsogr,    & ! GROUND LIDAR high-level cloud fraction          !GLID
         Lcllcalipsogr,    & ! GROUND LIDAR low-level cloud fraction           !GLID
         Lclmcalipsogr,    & ! GROUND LIDAR mid-level cloud fraction           !GLID
         Lcltcalipsogr,    & ! GROUND LIDAR total cloud fraction               !GLID
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
         Lclopaquecalipso, & ! CALIPSO opaque cloud cover (2D Map)                   !OPAQ
         Lclthincalipso,   & ! CALIPSO thin cloud cover (2D Map)                     !OPAQ
         Lclzopaquecalipso,& ! CALIPSO z_opaque altitude for opaque clouds only (2D Map)!OPAQ
         Lclcalipsoopaque, & ! CALIPSO opaque cloud profiles 3D fraction             !OPAQ
         Lclcalipsothin,   & ! CALIPSO thin cloud profiles 3D fraction               !OPAQ
         Lclcalipsozopaque,& ! CALIPSO z_opaque 3D fraction                          !OPAQ
         Lclcalipsoopacity,& ! CALIPSO opacity 3D fraction                           !OPAQ
         Lclopaquetemp,    & ! CALIPSO opaque cloud temperature                   !TIBO
         Lclthintemp,      & ! CALIPSO thin cloud temperature                     !TIBO
         Lclzopaquetemp,   & ! CALIPSO z_opaque temperature                       !TIBO
         Lclopaquemeanz,   & ! CALIPSO opaque cloud altitude                      !TIBO
         Lclthinmeanz,     & ! CALIPSO thin cloud altitude                        !TIBO
         Lclthinemis,      & ! CALIPSO thin cloud emissivity                      !TIBO
         Lclopaquemeanzse,   & ! CALIPSO opaque cloud altitude with respect to SE !TIBO2
         Lclthinmeanzse,     & ! CALIPSO thin cloud altitude with respect to SE   !TIBO2
         Lclzopaquecalipsose,& ! CALIPSO z_opaque altitude with respect to SE     !TIBO2
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
    if (LlidarBetaMol532gr) allocate(x%calipso_beta_mol_gr(Npoints,Nlevels))      !GLID
    if (Latb532)          allocate(x%calipso_beta_tot(Npoints,Ncolumns,Nlevels))
    if (LcfadLidarsr532)  then
        allocate(x%calipso_srbval(SR_BINS+1))
        allocate(x%calipso_cfad_sr(Npoints,SR_BINS,Nlvgrid))
        allocate(x%calipso_betaperp_tot(Npoints,Ncolumns,Nlevels))  
    endif
    if (LcfadLidarsr532gr) then                                                   !GLID
        allocate(x%calipso_srbval(SR_BINS+1))                                     !GLID
        allocate(x%calipso_cfad_sr_gr(Npoints,SR_BINS,Nlvgrid))                   !GLID
    endif                                                                         !GLID
    if (Lclcalipso)       allocate(x%calipso_lidarcld(Npoints,Nlvgrid))
    if (Lclcalipsogr)     allocate(x%calipso_lidarcld_gr(Npoints,Nlvgrid))        !GLID
    if (Lclhcalipso .or. Lclmcalipso .or. Lcllcalipso .or. Lcltcalipso) then
        allocate(x%calipso_cldlayer(Npoints,LIDAR_NCAT))
    endif   
    if (Lclhcalipsogr .or. Lclmcalipsogr .or. Lcllcalipsogr .or. Lcltcalipsogr) then     !GLID
        allocate(x%calipso_cldlayer_gr(Npoints,LIDAR_NCAT))                              !GLID
    endif                                                                                !GLID
    if (Lclcalipsoice .or. Lclcalipsoliq .or. Lclcalipsoun) then
        allocate(x%calipso_lidarcldphase(Npoints,Nlvgrid,6))
    endif
    if (Lclcalipsotmp .or. Lclcalipsotmpliq .or. Lclcalipsoice .or. Lclcalipsotmpun) then
        allocate(x%calipso_lidarcldtmp(Npoints,LIDAR_NTEMP,5))
    endif
    if (Lcllcalipsoice .or. Lclmcalipsoice .or. Lclhcalipsoice .or.                   &
        Lcltcalipsoice .or. Lcllcalipsoliq .or. Lclmcalipsoliq .or.                   &
        Lclhcalipsoliq .or. Lcltcalipsoliq .or. Lcllcalipsoun  .or.                   &
        Lclmcalipsoun  .or. Lclhcalipsoun  .or. Lcltcalipsoun) then
        allocate(x%calipso_cldlayerphase(Npoints,LIDAR_NCAT,6))     
    endif
    if (Lclopaquecalipso .or. Lclthincalipso .or. Lclzopaquecalipso) then !OPAQ
        allocate(x%calipso_cldtype(Npoints,LIDAR_NTYPE))                  !OPAQ
    endif                                                                 !OPAQ
    if (Lclopaquetemp .or. Lclthintemp .or. Lclzopaquetemp) then          !TIBO
        allocate(x%calipso_cldtypetemp(Npoints,LIDAR_NTYPE))              !TIBO
    endif                                                                 !TIBO
    if (Lclopaquemeanz .or. Lclthinmeanz) then                            !TIBO
        allocate(x%calipso_cldtypemeanz(Npoints,2))                       !TIBO
    endif                                                                 !TIBO
    if (Lclopaquemeanzse .or. Lclthinmeanzse .or. Lclzopaquecalipsose) then !TIBO2
        allocate(x%calipso_cldtypemeanzse(Npoints,3))                       !TIBO2
    endif                                                                   !TIBO2
    if (Lclthinemis) then                                                 !TIBO
        allocate(x%calipso_cldthinemis(Npoints))                          !TIBO
    endif                                                                 !TIBO
    if (Lclcalipsoopaque .or. Lclcalipsothin .or. Lclcalipsozopaque .or. Lclcalipsoopacity) then !OPAQ
        allocate(x%calipso_lidarcldtype(Npoints,Nlvgrid,LIDAR_NTYPE+1))   !OPAQ
    endif                                                                 !OPAQ

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
    if (LcfadDbze94)    allocate(x%cloudsat_cfad_ze(Npoints,DBZE_BINS,Nlvgrid))

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
    deallocate(y%tau_067,y%emiss_11,y%frac_out,y%beta_mol,y%tau_mol,y%betatot,           &
               y%betatot_ice,y%betatot_liq,y%tautot,y%tautot_ice,y%tautot_liq,           &
               y%tautot_S_liq,y%tautot_S_ice,y%z_vol_cloudsat,y%kr_vol_cloudsat,         &
               y%g_vol_cloudsat,y%asym,y%ss_alb,y%fracLiq,y%taupart,                     & !GLID
               y%beta_mol_gr,y%betatot_gr,y%taupart_gr,y%tau_mol_gr,y%tautot_gr)           !GLID

  end subroutine destroy_cospIN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE destroy_cospstateIN     
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  subroutine destroy_cospstateIN(y)
    type(cosp_column_inputs),intent(inout) :: y
    deallocate(y%sunlit,y%skt,y%land,y%at,y%pfull ,y%phalf,y%qv,y%o3,y%hgt_matrix,       &
               y%u_sfc,y%v_sfc,y%lat,y%lon,y%emis_sfc,y%cloudIce,y%cloudLiq,y%seaice,    &
               y%fl_rain,y%fl_snow,y%tca,y%hgt_matrix_half,y%surfelev)                     !TIBO2

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
     if (associated(y%calipso_beta_mol_gr))       then !GLID
        deallocate(y%calipso_beta_mol_gr)              !GLID
        nullify(y%calipso_beta_mol_gr)                 !GLID
     endif                                             !GLID
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
     if (associated(y%calipso_lidarcldtype))     then !OPAQ
        deallocate(y%calipso_lidarcldtype)            !OPAQ
        nullify(y%calipso_lidarcldtype)               !OPAQ
     endif                                            !OPAQ
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
     if (associated(y%calipso_cldlayer_gr))      then !GLID
        deallocate(y%calipso_cldlayer_gr)             !GLID
        nullify(y%calipso_cldlayer_gr)                !GLID
     endif                                            !GLID
     if (associated(y%calipso_cldtype))          then !OPAQ
        deallocate(y%calipso_cldtype)                 !OPAQ
        nullify(y%calipso_cldtype)                    !OPAQ
     endif                                            !OPAQ
     if (associated(y%calipso_cldtypetemp))      then !TIBO
        deallocate(y%calipso_cldtypetemp)             !TIBO
        nullify(y%calipso_cldtypetemp)                !TIBO
     endif                                            !TIBO
     if (associated(y%calipso_cldtypemeanz))     then !TIBO
        deallocate(y%calipso_cldtypemeanz)            !TIBO
        nullify(y%calipso_cldtypemeanz)               !TIBO
     endif                                            !TIBO
     if (associated(y%calipso_cldtypemeanzse))   then !TIBO2
        deallocate(y%calipso_cldtypemeanzse)          !TIBO2
        nullify(y%calipso_cldtypemeanzse)             !TIBO2
     endif                                            !TIBO2
     if (associated(y%calipso_cldthinemis))      then !TIBO
        deallocate(y%calipso_cldthinemis)             !TIBO
        nullify(y%calipso_cldthinemis)                !TIBO
     endif                                            !TIBO
     if (associated(y%calipso_lidarcld))         then
        deallocate(y%calipso_lidarcld)
        nullify(y%calipso_lidarcld)     
     endif
     if (associated(y%calipso_lidarcld_gr))      then !GLID
        deallocate(y%calipso_lidarcld_gr)             !GLID
        nullify(y%calipso_lidarcld_gr)                !GLID
     endif                                            !GLID
     if (associated(y%calipso_srbval))           then
        deallocate(y%calipso_srbval)
        nullify(y%calipso_srbval)     
     endif
     if (associated(y%calipso_cfad_sr))          then
        deallocate(y%calipso_cfad_sr)
        nullify(y%calipso_cfad_sr)     
     endif
     if (associated(y%calipso_cfad_sr_gr))       then !GLID
        deallocate(y%calipso_cfad_sr_gr)              !GLID
        nullify(y%calipso_cfad_sr_gr)                 !GLID
     endif                                            !GLID
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
   ! SUBROUTINE cosp_cleanUp
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   subroutine cosp_cleanUp()
     deallocate(vgrid_zl,vgrid_zu,vgrid_z)
   end subroutine cosp_cleanUp
   
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_errorCheck
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine cosp_errorCheck(cospgridIN,cospIN,Lisccp_subcolumn,Lisccp_column,     & !GLID
                             Lmisr_subcolumn,Lmisr_column,Lmodis_subcolumn,        & !GLID
                             Lmodis_column,Lcloudsat_subcolumn,Lcloudsat_column,   & !GLID
                             Lcalipso_subcolumn,Lcalipso_subcolumn_gr,             & !GLID
                             Lcalipso_column,Lcalipso_column_gr,Lrttov_subcolumn,  & !GLID
                             Lrttov_column,Lparasol_subcolumn,Lparasol_column,     & !GLID
                             Lradar_lidar_tcc,Llidar_only_freq_cloud,cospOUT,errorMessage,nError)
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
      Lcalipso_subcolumn_gr, & ! GROUND LIDAR subcolumn simulator on/off switch !GLID
      Lcalipso_column,     & ! CALIPSO column simulator on/off switch
      Lcalipso_column_gr,  & ! GROUND LIDAR column simulator on/off switch      !GLID
      Lparasol_subcolumn,  & ! PARASOL subcolumn simulator on/off switch
      Lparasol_column,     & ! PARASOL column simulator on/off switch
      Lrttov_subcolumn,    & ! RTTOV subcolumn simulator on/off switch
      Lrttov_column,       & ! RTTOV column simulator on/off switch       
      Lradar_lidar_tcc,    & ! On/Off switch for joint Calipso/Cloudsat product
      Llidar_only_freq_cloud ! On/Off switch for joint Calipso/Cloudsat product
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
       Lcalipso_column_gr  = .false. !GLID
       Lcloudsat_column = .false.
       Lradar_lidar_tcc = .false.
       Llidar_only_freq_cloud = .false.
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
       if (associated(cospOUT%calipso_cfad_sr_gr))    cospOUT%calipso_cfad_sr_gr(:,:,:)    = R_UNDEF !GLID
       if (associated(cospOUT%calipso_lidarcld))      cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_lidarcld_gr))   cospOUT%calipso_lidarcld_gr(:,:)     = R_UNDEF !GLID
       if (associated(cospOUT%calipso_lidarcldphase)) cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldtype))  cospOUT%calipso_lidarcldtype(:,:,:)  = R_UNDEF !OPAQ
       if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_cldlayer_gr))   cospOUT%calipso_cldlayer_gr(:,:)     = R_UNDEF !GLID
       if (associated(cospOUT%calipso_cldtype))       cospOUT%calipso_cldtype(:,:)         = R_UNDEF !OPAQ
       if (associated(cospOUT%calipso_cldtypetemp))   cospOUT%calipso_cldtypetemp(:,:)     = R_UNDEF !TIBO
       if (associated(cospOUT%calipso_cldtypemeanz))  cospOUT%calipso_cldtypemeanz(:,:)    = R_UNDEF !TIBO
       if (associated(cospOUT%calipso_cldtypemeanzse)) cospOUT%calipso_cldtypemeanzse(:,:) = R_UNDEF !TIBO2
       if (associated(cospOUT%calipso_cldthinemis))   cospOUT%calipso_cldthinemis(:)       = R_UNDEF !TIBO
       if (associated(cospOUT%calipso_cldlayerphase)) cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldtmp))   cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF            
       if (associated(cospOUT%cloudsat_cfad_ze))      cospOUT%cloudsat_cfad_ze(:,:,:)      = R_UNDEF
       if (associated(cospOUT%lidar_only_freq_cloud)) cospOUT%lidar_only_freq_cloud(:,:)   = R_UNDEF
       if (associated(cospOUT%radar_lidar_tcc))       cospOUT%radar_lidar_tcc(:)           = R_UNDEF       
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
       Lcalipso_column_gr  = .false. !GLID
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
       if (associated(cospOUT%calipso_cfad_sr_gr))    cospOUT%calipso_cfad_sr_gr(:,:,:)    = R_UNDEF !GLID
       if (associated(cospOUT%calipso_lidarcld))      cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_lidarcld_gr))   cospOUT%calipso_lidarcld_gr(:,:)     = R_UNDEF !GLID
       if (associated(cospOUT%calipso_lidarcldphase)) cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldtype))  cospOUT%calipso_lidarcldtype(:,:,:)  = R_UNDEF !OPAQ
       if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_cldlayer_gr))   cospOUT%calipso_cldlayer_gr(:,:)     = R_UNDEF !GLID
       if (associated(cospOUT%calipso_cldtype))       cospOUT%calipso_cldtype(:,:)         = R_UNDEF !OPAQ
       if (associated(cospOUT%calipso_cldtypetemp))   cospOUT%calipso_cldtypetemp(:,:)     = R_UNDEF !TIBO
       if (associated(cospOUT%calipso_cldtypemeanz))  cospOUT%calipso_cldtypemeanz(:,:)    = R_UNDEF !TIBO
       if (associated(cospOUT%calipso_cldtypemeanzse)) cospOUT%calipso_cldtypemeanzse(:,:) = R_UNDEF !TIBO2
       if (associated(cospOUT%calipso_cldthinemis))   cospOUT%calipso_cldthinemis(:)       = R_UNDEF !TIBO
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
       Lcalipso_column_gr  = .false. !GLID
       Lradar_lidar_tcc = .false.
       Llidar_only_freq_cloud = .false.
       if (associated(cospOUT%misr_fq))                   cospOUT%misr_fq(:,:,:)                 = R_UNDEF
       if (associated(cospOUT%misr_dist_model_layertops)) cospOUT%misr_dist_model_layertops(:,:) = R_UNDEF
       if (associated(cospOUT%misr_meanztop))             cospOUT%misr_meanztop(:)               = R_UNDEF
       if (associated(cospOUT%misr_cldarea))              cospOUT%misr_cldarea(:)                = R_UNDEF
       if (associated(cospOUT%calipso_cfad_sr))           cospOUT%calipso_cfad_sr(:,:,:)         = R_UNDEF
       if (associated(cospOUT%calipso_cfad_sr_gr))        cospOUT%calipso_cfad_sr_gr(:,:,:)      = R_UNDEF !GLID
       if (associated(cospOUT%calipso_lidarcld))          cospOUT%calipso_lidarcld(:,:)          = R_UNDEF
       if (associated(cospOUT%calipso_lidarcld_gr))       cospOUT%calipso_lidarcld_gr(:,:)       = R_UNDEF !GLID
       if (associated(cospOUT%calipso_lidarcldphase))     cospOUT%calipso_lidarcldphase(:,:,:)   = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldtype))      cospOUT%calipso_lidarcldtype(:,:,:)    = R_UNDEF !OPAQ
       if (associated(cospOUT%calipso_cldlayer))          cospOUT%calipso_cldlayer(:,:)          = R_UNDEF
       if (associated(cospOUT%calipso_cldlayer_gr))       cospOUT%calipso_cldlayer_gr(:,:)       = R_UNDEF !GLID
       if (associated(cospOUT%calipso_cldtype))           cospOUT%calipso_cldtype(:,:)           = R_UNDEF !OPAQ
       if (associated(cospOUT%calipso_cldtypetemp))       cospOUT%calipso_cldtypetemp(:,:)     = R_UNDEF !TIBO
       if (associated(cospOUT%calipso_cldtypemeanz))      cospOUT%calipso_cldtypemeanz(:,:)    = R_UNDEF !TIBO
       if (associated(cospOUT%calipso_cldtypemeanzse))    cospOUT%calipso_cldtypemeanzse(:,:)  = R_UNDEF !TIBO2
       if (associated(cospOUT%calipso_cldthinemis))       cospOUT%calipso_cldthinemis(:)       = R_UNDEF !TIBO
       if (associated(cospOUT%calipso_cldlayerphase))     cospOUT%calipso_cldlayerphase(:,:,:)   = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldtmp))       cospOUT%calipso_lidarcldtmp(:,:,:)     = R_UNDEF            
       if (associated(cospOUT%cloudsat_cfad_ze))          cospOUT%cloudsat_cfad_ze(:,:,:)        = R_UNDEF
       if (associated(cospOUT%cloudsat_Ze_tot))           cospOUT%cloudsat_Ze_tot(:,:,:)         = R_UNDEF
       if (associated(cospOUT%lidar_only_freq_cloud))     cospOUT%lidar_only_freq_cloud(:,:)     = R_UNDEF
       if (associated(cospOUT%radar_lidar_tcc))           cospOUT%radar_lidar_tcc(:)             = R_UNDEF       
    endif
    if (any(cospgridIN%hgt_matrix_half .lt. -300)) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospgridIN%hgt_matrix_half contains values out of range'
       Lrttov_subcolumn = .false.
       Lcloudsat_column = .false.
       Lcalipso_column  = .false.
       Lcalipso_column_gr  = .false. !GLID
       Lradar_lidar_tcc = .false.
       Llidar_only_freq_cloud = .false.
       if (associated(cospOUT%rttov_tbs))             cospOUT%rttov_tbs(:,:)               = R_UNDEF       
       if (associated(cospOUT%calipso_cfad_sr))       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF
       if (associated(cospOUT%calipso_cfad_sr_gr))    cospOUT%calipso_cfad_sr_gr(:,:,:)    = R_UNDEF !GLID
       if (associated(cospOUT%calipso_lidarcld))      cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_lidarcld_gr))   cospOUT%calipso_lidarcld_gr(:,:)     = R_UNDEF !GLID
       if (associated(cospOUT%calipso_lidarcldphase)) cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldtype))  cospOUT%calipso_lidarcldtype(:,:,:)  = R_UNDEF !OPAQ
       if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_cldlayer_gr))   cospOUT%calipso_cldlayer_gr(:,:)     = R_UNDEF !GLID
       if (associated(cospOUT%calipso_cldtype))       cospOUT%calipso_cldtype(:,:)         = R_UNDEF !OPAQ
       if (associated(cospOUT%calipso_cldtypetemp))   cospOUT%calipso_cldtypetemp(:,:)     = R_UNDEF !TIBO
       if (associated(cospOUT%calipso_cldtypemeanz))  cospOUT%calipso_cldtypemeanz(:,:)    = R_UNDEF !TIBO
       if (associated(cospOUT%calipso_cldtypemeanzse)) cospOUT%calipso_cldtypemeanzse(:,:) = R_UNDEF !TIBO2
       if (associated(cospOUT%calipso_cldthinemis))   cospOUT%calipso_cldthinemis(:)       = R_UNDEF !TIBO
       if (associated(cospOUT%calipso_cldlayerphase)) cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldtmp))   cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF            
       if (associated(cospOUT%cloudsat_cfad_ze))      cospOUT%cloudsat_cfad_ze(:,:,:)      = R_UNDEF
       if (associated(cospOUT%lidar_only_freq_cloud)) cospOUT%lidar_only_freq_cloud(:,:)   = R_UNDEF
       if (associated(cospOUT%radar_lidar_tcc))       cospOUT%radar_lidar_tcc(:)           = R_UNDEF                 
    endif
    if (any(cospgridIN%land .lt. 0)) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospgridIN%land contains values out of range'
       Lrttov_subcolumn = .false.
       Lcalipso_column  = .false.       
       Lcalipso_column_gr  = .false. !GLID       
       Lparasol_column  = .false.
       if (associated(cospOUT%rttov_tbs))             cospOUT%rttov_tbs(:,:)               = R_UNDEF       
       if (associated(cospOUT%calipso_cfad_sr))       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF
       if (associated(cospOUT%calipso_cfad_sr_gr))    cospOUT%calipso_cfad_sr_gr(:,:,:)    = R_UNDEF !GLID
       if (associated(cospOUT%calipso_lidarcld))      cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_lidarcld_gr))   cospOUT%calipso_lidarcld_gr(:,:)     = R_UNDEF !GLID
       if (associated(cospOUT%calipso_lidarcldphase)) cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldtype))  cospOUT%calipso_lidarcldtype(:,:,:)  = R_UNDEF !OPAQ
       if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_cldlayer_gr))   cospOUT%calipso_cldlayer_gr(:,:)     = R_UNDEF !GLID
       if (associated(cospOUT%calipso_cldtype))       cospOUT%calipso_cldtype(:,:)         = R_UNDEF !OPAQ
       if (associated(cospOUT%calipso_cldtypetemp))   cospOUT%calipso_cldtypetemp(:,:)     = R_UNDEF !TIBO
       if (associated(cospOUT%calipso_cldtypemeanz))  cospOUT%calipso_cldtypemeanz(:,:)    = R_UNDEF !TIBO
       if (associated(cospOUT%calipso_cldtypemeanzse)) cospOUT%calipso_cldtypemeanzse(:,:) = R_UNDEF !TIBO2
       if (associated(cospOUT%calipso_cldthinemis))   cospOUT%calipso_cldthinemis(:)       = R_UNDEF !TIBO
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
       if (associated(cospOUT%calipso_lidarcldtype))  cospOUT%calipso_lidarcldtype(:,:,:)  = R_UNDEF !OPAQ
       if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_cldtype))       cospOUT%calipso_cldtype(:,:)         = R_UNDEF !OPAQ
       if (associated(cospOUT%calipso_cldtypetemp))   cospOUT%calipso_cldtypetemp(:,:)     = R_UNDEF !TIBO
       if (associated(cospOUT%calipso_cldtypemeanz))  cospOUT%calipso_cldtypemeanz(:,:)    = R_UNDEF !TIBO
       if (associated(cospOUT%calipso_cldtypemeanzse)) cospOUT%calipso_cldtypemeanzse(:,:) = R_UNDEF !TIBO2
       if (associated(cospOUT%calipso_cldthinemis))   cospOUT%calipso_cldthinemis(:)       = R_UNDEF !TIBO
       if (associated(cospOUT%calipso_cldlayerphase)) cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldtmp))   cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF
       if (associated(cospOUT%calipso_srbval))        cospOUT%calipso_srbval(:)            = R_UNDEF
    endif

    if (any(cospIN%betatot_gr .lt. 0)) then !GLID
       nError=nError+1                      !GLID
       errorMessage(nError) = 'ERROR: COSP input variable: cospIN%betatot_gr contains values out of range' !GLID
       Lcalipso_subcolumn_gr = .false. !GLID
       Lcalipso_column_gr    = .false. !GLID
       if (associated(cospOUT%calipso_cfad_sr_gr))    cospOUT%calipso_cfad_sr_gr(:,:,:)    = R_UNDEF !GLID
       if (associated(cospOUT%calipso_lidarcld_gr))   cospOUT%calipso_lidarcld_gr(:,:)     = R_UNDEF !GLID
       if (associated(cospOUT%calipso_cldlayer_gr))   cospOUT%calipso_cldlayer_gr(:,:)     = R_UNDEF !GLID
    endif !GLID

    if (any(cospIN%betatot_liq .lt. 0)) then
       nError=nError+1
       errorMessage(nError) = ('ERROR: COSP input variable: cospIN%betatot_liq contains values out of range')
       Lcalipso_subcolumn = .false.
       Lcalipso_column    = .false.
       if (associated(cospOUT%calipso_cfad_sr))       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF
       if (associated(cospOUT%calipso_lidarcld))      cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldphase)) cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldtype))  cospOUT%calipso_lidarcldtype(:,:,:)  = R_UNDEF !OPAQ
       if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_cldtype))       cospOUT%calipso_cldtype(:,:)         = R_UNDEF !OPAQ
       if (associated(cospOUT%calipso_cldtypetemp))   cospOUT%calipso_cldtypetemp(:,:)     = R_UNDEF !TIBO
       if (associated(cospOUT%calipso_cldtypemeanz))  cospOUT%calipso_cldtypemeanz(:,:)    = R_UNDEF !TIBO
       if (associated(cospOUT%calipso_cldtypemeanzse)) cospOUT%calipso_cldtypemeanzse(:,:) = R_UNDEF !TIBO2
       if (associated(cospOUT%calipso_cldthinemis))   cospOUT%calipso_cldthinemis(:)       = R_UNDEF !TIBO
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
       if (associated(cospOUT%calipso_lidarcldtype))  cospOUT%calipso_lidarcldtype(:,:,:)  = R_UNDEF !OPAQ
       if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_cldtype))       cospOUT%calipso_cldtype(:,:)         = R_UNDEF !OPAQ
       if (associated(cospOUT%calipso_cldtypetemp))   cospOUT%calipso_cldtypetemp(:,:)     = R_UNDEF !TIBO
       if (associated(cospOUT%calipso_cldtypemeanz))  cospOUT%calipso_cldtypemeanz(:,:)    = R_UNDEF !TIBO
       if (associated(cospOUT%calipso_cldtypemeanzse)) cospOUT%calipso_cldtypemeanzse(:,:) = R_UNDEF !TIBO2
       if (associated(cospOUT%calipso_cldthinemis))   cospOUT%calipso_cldthinemis(:)       = R_UNDEF !TIBO
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
       if (associated(cospOUT%calipso_lidarcldtype))  cospOUT%calipso_lidarcldtype(:,:,:)  = R_UNDEF !OPAQ
       if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_cldtype))       cospOUT%calipso_cldtype(:,:)         = R_UNDEF !OPAQ
       if (associated(cospOUT%calipso_cldtypetemp))   cospOUT%calipso_cldtypetemp(:,:)     = R_UNDEF !TIBO
       if (associated(cospOUT%calipso_cldtypemeanz))  cospOUT%calipso_cldtypemeanz(:,:)    = R_UNDEF !TIBO
       if (associated(cospOUT%calipso_cldtypemeanzse)) cospOUT%calipso_cldtypemeanzse(:,:) = R_UNDEF !TIBO2
       if (associated(cospOUT%calipso_cldthinemis))   cospOUT%calipso_cldthinemis(:)       = R_UNDEF !TIBO
       if (associated(cospOUT%calipso_cldlayerphase)) cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldtmp))   cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF
       if (associated(cospOUT%calipso_srbval))        cospOUT%calipso_srbval(:)            = R_UNDEF
       if (associated(cospOUT%cloudsat_cfad_ze))      cospOUT%cloudsat_cfad_ze(:,:,:)      = R_UNDEF
       if (associated(cospOUT%lidar_only_freq_cloud)) cospOUT%lidar_only_freq_cloud(:,:)   = R_UNDEF
       if (associated(cospOUT%radar_lidar_tcc))       cospOUT%radar_lidar_tcc(:)           = R_UNDEF          
    endif    

    if (any(cospIN%beta_mol_gr .lt. 0)) then !GLID
       nError=nError+1 !GLID
       errorMessage(nError) = 'ERROR: COSP input variable: cospIN%beta_mol_gr contains values out of range' !GLID
       Lcalipso_subcolumn_gr = .false. !GLID
       Lcalipso_column_gr    = .false. !GLID
       if (associated(cospOUT%calipso_cfad_sr_gr))    cospOUT%calipso_cfad_sr_gr(:,:,:)    = R_UNDEF !GLID
       if (associated(cospOUT%calipso_lidarcld_gr))   cospOUT%calipso_lidarcld_gr(:,:)     = R_UNDEF !GLID
       if (associated(cospOUT%calipso_cldlayer_gr))   cospOUT%calipso_cldlayer_gr(:,:)     = R_UNDEF !GLID
    endif !GLID

    if (any(cospIN%tautot .lt. 0)) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospIN%tautot contains values out of range'
       Lcalipso_subcolumn = .false.
       Lcalipso_column    = .false.
       if (associated(cospOUT%calipso_cfad_sr))       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF
       if (associated(cospOUT%calipso_lidarcld))      cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldphase)) cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldtype))  cospOUT%calipso_lidarcldtype(:,:,:)  = R_UNDEF !OPAQ
       if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_cldtype))       cospOUT%calipso_cldtype(:,:)         = R_UNDEF !OPAQ
       if (associated(cospOUT%calipso_cldtypetemp))   cospOUT%calipso_cldtypetemp(:,:)     = R_UNDEF !TIBO
       if (associated(cospOUT%calipso_cldtypemeanz))  cospOUT%calipso_cldtypemeanz(:,:)    = R_UNDEF !TIBO
       if (associated(cospOUT%calipso_cldtypemeanzse)) cospOUT%calipso_cldtypemeanzse(:,:) = R_UNDEF !TIBO2
       if (associated(cospOUT%calipso_cldthinemis))   cospOUT%calipso_cldthinemis(:)       = R_UNDEF !TIBO
       if (associated(cospOUT%calipso_cldlayerphase)) cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldtmp))   cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF
       if (associated(cospOUT%calipso_srbval))        cospOUT%calipso_srbval(:)            = R_UNDEF      
    endif

    if (any(cospIN%tautot_gr .lt. 0)) then !GLID
       nError=nError+1 !GLID
       errorMessage(nError) = 'ERROR: COSP input variable: cospIN%tautot_gr contains values out of range' !GLID
       Lcalipso_subcolumn_gr = .false. !GLID
       Lcalipso_column_gr    = .false. !GLID
       if (associated(cospOUT%calipso_cfad_sr_gr))    cospOUT%calipso_cfad_sr_gr(:,:,:)    = R_UNDEF !GLID
       if (associated(cospOUT%calipso_lidarcld_gr))   cospOUT%calipso_lidarcld_gr(:,:)     = R_UNDEF !GLID
       if (associated(cospOUT%calipso_cldlayer_gr))   cospOUT%calipso_cldlayer_gr(:,:)     = R_UNDEF !GLID
    endif !GLID

    if (any(cospIN%tautot_liq .lt. 0)) then
       nError=nError+1
       errorMessage(nError) = ('ERROR: COSP input variable: cospIN%tautot_liq contains values out of range')
       Lcalipso_subcolumn = .false.
       Lcalipso_column    = .false.
       if (associated(cospOUT%calipso_cfad_sr))       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF
       if (associated(cospOUT%calipso_lidarcld))      cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldphase)) cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldtype))  cospOUT%calipso_lidarcldtype(:,:,:)  = R_UNDEF !OPAQ
       if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_cldtype))       cospOUT%calipso_cldtype(:,:)         = R_UNDEF !OPAQ
       if (associated(cospOUT%calipso_cldtypetemp))   cospOUT%calipso_cldtypetemp(:,:)     = R_UNDEF !TIBO
       if (associated(cospOUT%calipso_cldtypemeanz))  cospOUT%calipso_cldtypemeanz(:,:)    = R_UNDEF !TIBO
       if (associated(cospOUT%calipso_cldtypemeanzse)) cospOUT%calipso_cldtypemeanzse(:,:) = R_UNDEF !TIBO2
       if (associated(cospOUT%calipso_cldthinemis))   cospOUT%calipso_cldthinemis(:)       = R_UNDEF !TIBO
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
       if (associated(cospOUT%calipso_lidarcldtype))  cospOUT%calipso_lidarcldtype(:,:,:)  = R_UNDEF !OPAQ
       if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_cldtype))       cospOUT%calipso_cldtype(:,:)         = R_UNDEF !OPAQ
       if (associated(cospOUT%calipso_cldtypetemp))   cospOUT%calipso_cldtypetemp(:,:)     = R_UNDEF !TIBO
       if (associated(cospOUT%calipso_cldtypemeanz))  cospOUT%calipso_cldtypemeanz(:,:)    = R_UNDEF !TIBO
       if (associated(cospOUT%calipso_cldtypemeanzse)) cospOUT%calipso_cldtypemeanzse(:,:) = R_UNDEF !TIBO2
       if (associated(cospOUT%calipso_cldthinemis))   cospOUT%calipso_cldthinemis(:)       = R_UNDEF !TIBO
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
       if (associated(cospOUT%calipso_lidarcldtype))  cospOUT%calipso_lidarcldtype(:,:,:)  = R_UNDEF !OPAQ
       if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
       if (associated(cospOUT%calipso_cldtype))       cospOUT%calipso_cldtype(:,:)         = R_UNDEF !OPAQ
       if (associated(cospOUT%calipso_cldtypetemp))   cospOUT%calipso_cldtypetemp(:,:)     = R_UNDEF !TIBO
       if (associated(cospOUT%calipso_cldtypemeanz))  cospOUT%calipso_cldtypemeanz(:,:)    = R_UNDEF !TIBO
       if (associated(cospOUT%calipso_cldtypemeanzse)) cospOUT%calipso_cldtypemeanzse(:,:) = R_UNDEF !TIBO2
       if (associated(cospOUT%calipso_cldthinemis))   cospOUT%calipso_cldthinemis(:)       = R_UNDEF !TIBO
       if (associated(cospOUT%calipso_cldlayerphase)) cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
       if (associated(cospOUT%calipso_lidarcldtmp))   cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF
       if (associated(cospOUT%calipso_srbval))        cospOUT%calipso_srbval(:)            = R_UNDEF         
    endif   

    if (any(cospIN%tau_mol_gr .lt. 0)) then !GLID
       nError=nError+1 !GLID
       errorMessage(nError) = 'ERROR: COSP input variable: cospIN%tau_mol_gr contains values out of range' !GLID
       Lcalipso_subcolumn_gr = .false. !GLID
       Lcalipso_column_gr    = .false. !GLID
       if (associated(cospOUT%calipso_cfad_sr_gr))    cospOUT%calipso_cfad_sr_gr(:,:,:)    = R_UNDEF !GLID
       if (associated(cospOUT%calipso_lidarcld_gr))   cospOUT%calipso_lidarcld_gr(:,:)     = R_UNDEF !GLID
       if (associated(cospOUT%calipso_cldlayer_gr))   cospOUT%calipso_cldlayer_gr(:,:)     = R_UNDEF !GLID
    endif !GLID
 
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
       Llidar_only_freq_cloud = .false.
       if (associated(cospOUT%cloudsat_cfad_ze))          cospOUT%cloudsat_cfad_ze(:,:,:)        = R_UNDEF
       if (associated(cospOUT%cloudsat_Ze_tot))           cospOUT%cloudsat_Ze_tot(:,:,:)         = R_UNDEF
       if (associated(cospOUT%lidar_only_freq_cloud))     cospOUT%lidar_only_freq_cloud(:,:)     = R_UNDEF
       if (associated(cospOUT%radar_lidar_tcc))           cospOUT%radar_lidar_tcc(:)             = R_UNDEF     
    endif
    if (any(cospIN%kr_vol_cloudsat .lt. 0)) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospIN%kr_vol_cloudsat contains values out of range'
       Lcloudsat_subcolumn = .false.
       Lcloudsat_column    = .false.
       Lradar_lidar_tcc    = .false.
       Llidar_only_freq_cloud = .false.
       if (associated(cospOUT%cloudsat_cfad_ze))          cospOUT%cloudsat_cfad_ze(:,:,:)        = R_UNDEF
       if (associated(cospOUT%cloudsat_Ze_tot))           cospOUT%cloudsat_Ze_tot(:,:,:)         = R_UNDEF
       if (associated(cospOUT%lidar_only_freq_cloud))     cospOUT%lidar_only_freq_cloud(:,:)     = R_UNDEF
       if (associated(cospOUT%radar_lidar_tcc))           cospOUT%radar_lidar_tcc(:)             = R_UNDEF      
    endif    
    if (any(cospIN%g_vol_cloudsat .lt. 0)) then
       nError=nError+1
       errorMessage(nError) = 'ERROR: COSP input variable: cospIN%g_vol_cloudsat contains values out of range'
       Lcloudsat_subcolumn = .false.
       Lcloudsat_column    = .false.
       Lradar_lidar_tcc    = .false.
       Llidar_only_freq_cloud = .false.
       if (associated(cospOUT%cloudsat_cfad_ze))          cospOUT%cloudsat_cfad_ze(:,:,:)        = R_UNDEF
       if (associated(cospOUT%cloudsat_Ze_tot))           cospOUT%cloudsat_Ze_tot(:,:,:)         = R_UNDEF
       if (associated(cospOUT%lidar_only_freq_cloud))     cospOUT%lidar_only_freq_cloud(:,:)     = R_UNDEF
       if (associated(cospOUT%radar_lidar_tcc))           cospOUT%radar_lidar_tcc(:)             = R_UNDEF
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
  
  ! GROUND LIDAR !GLID
  if (size(cospIN%beta_mol_gr,1)    .ne. cospIN%Npoints .OR. & !GLID
      size(cospIN%betatot_gr,1)     .ne. cospIN%Npoints .OR. & !GLID
      size(cospIN%tau_mol_gr,1)     .ne. cospIN%Npoints .OR. & !GLID
      size(cospIN%tautot_gr,1)      .ne. cospIN%Npoints) then  !GLID
      Lcalipso_subcolumn_gr = .false. !GLID
      Lcalipso_column_gr    = .false. !GLID
      nError=nError+1                 !GLID
      errorMessage(nError) = 'ERROR(calipso_simulator): The number of points in the input fields are inconsistent' !GLID
  endif !GLID
  if (size(cospIN%betatot_gr,2)    .ne. cospIN%Ncolumns .OR. & !GLID
      size(cospIN%tautot_gr,2)     .ne. cospIN%Ncolumns) then  !GLID
      Lcalipso_subcolumn_gr = .false. !GLID
      Lcalipso_column_gr    = .false. !GLID
      nError=nError+1               !GLID
      errorMessage(nError) = 'ERROR(calipso_simulator): The number of sub-columns in the input fields are inconsistent' !GLID
  endif !GLID
  if (size(cospIN%beta_mol,2)    .ne. cospIN%Nlevels .OR. & !GLID
      size(cospIN%betatot,3)     .ne. cospIN%Nlevels .OR. & !GLID
      size(cospIN%tau_mol,2)     .ne. cospIN%Nlevels .OR. & !GLID
      size(cospIN%tautot,3)      .ne. cospIN%Nlevels) then  !GLID
      Lcalipso_subcolumn_gr = .false. !GLID
      Lcalipso_column_gr    = .false. !GLID
      nError=nError+1 !GLID
      errorMessage(nError) = 'ERROR(calipso_simulator): The number of levels in the input fields are inconsistent' !GLID
  endif !GLID

  ! PARASOL
  if (size(cospIN%tautot_S_liq,1) .ne. cospIN%Npoints .OR. &
      size(cospIN%tautot_S_ice,1) .ne. cospIN%Npoints) then
      Lparasol_subcolumn = .false.
      Lparasol_column    = .false.
      nError=nError+1
      errorMessage(nError) = 'ERROR(parasol_simulator): The number of points in the input fields are inconsistent'
  endif
  if (size(cospIN%tautot_S_liq,2) .ne. cospIN%Nlevels .OR. &
      size(cospIN%tautot_S_ice,2) .ne. cospIN%Nlevels) then
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
             
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Part 3: Instrument simulator specific error checking. This section contains error
  !         checking that was originally contained within the simulators.
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! PARASOL Simulator
  if (Lparasol_subcolumn) then
      IF (PARASOL_NREFL .gt. ntetas ) THEN
         write(parasolErrorMessage,"(a50,i2,a4,i2)") 'ERROR(lidar_simulator): nrefl should be less then ',ntetas,' not',PARASOL_NREFL
         nError=nError+1
         errorMessage(nError) = parasolErrorMessage
         Lparasol_subcolumn = .false.
         Lparasol_column    = .false.
      endif
  endif  
    
  end subroutine cosp_errorCheck
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! END MODULE
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
END MODULE MOD_COSP
