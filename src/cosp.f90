! (c) British Crown Copyright 2008, the Met Office.
! All rights reserved.
! 
! Redistribution and use iurce and binary forms, with or without modification, are permitted 
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

MODULE MOD_COSP
  USE COSP_KINDS,                  ONLY: wp
  USE MOD_COSP_CONFIG,             ONLY: R_UNDEF,PARASOL_NREFL,LIDAR_NCAT,SR_BINS,       &
                                         N_HYDRO,RTTOV_MAX_CHANNELS,numMISRHgtBins,      &
                                         DBZE_BINS,LIDAR_NTEMP,calipso_histBsct,         &
                                         use_vgrid,Nlvgrid,vgrid_zu,vgrid_zl,vgrid_z,    &
                                         mgrid_zl,mgrid_zu,mgrid_z,numMODISTauBins,      &
                                         numMODISPresBins,numMODISReffIceBins,           &
                                         numMODISReffLiqBins,numISCCPTauBins,            &
                                         numISCCPPresBins,numMISRTauBins
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
  USE MOD_PARASOL,                 ONLY: parasol_subcolumn,   parasol_column, ntetas
  use mod_cosp_rttov,              ONLY: rttov_subcolumn
  USE mod_quickbeam_optics,        ONLY: size_distribution,quickbeam_optics_init
  USE MOD_COSP_STATS,              ONLY: COSP_LIDAR_ONLY_CLOUD,COSP_CHANGE_VERTICAL_GRID
  use mod_cosp_error,              ONLY: errorMessage
  
  IMPLICIT NONE
  
  ! ######################################################################################
  ! TYPE cosp_column_inputs
  ! ######################################################################################
  type cosp_column_inputs
     integer :: &
          Npoints,             & ! Number of gridpoints.
          Ncolumns,            & ! Number of columns.
          Nlevels                ! Number of levels.
     integer,allocatable,dimension(:) :: &
          sunlit                 ! Sunlit flag

     real(wp),allocatable,dimension(:,:) :: &
          at,                  & ! Temperature
          pfull,               & ! Pressure
          phalf,               & ! Pressure at half-levels
          qv,                  & ! Specific humidity
          hgt_matrix,          & ! Height of hydrometeors (km)
          hgt_matrix_half        ! Height of hydrometeors at half levels (km)
     real(wp),allocatable,dimension(:) :: &
          land,                & ! Land/Sea mask
          skt                    ! Surface temperature
     ! RTTOV fields
     real(wp) :: &
          zenang,              & ! Satellite zenith angle for RTTOV       
          co2,                 & ! CO2
          ch4,                 & ! Methane
          n2o,                 & ! N2O
          co                     ! CO           
     real(wp),allocatable,dimension(:) :: &
          emis_sfc,            & ! Surface emissivity
          u_sfc,               & ! Surface u-wind
          v_sfc,               & ! Surface v-wind
          t_sfc,               & ! Skin temperature
          lat                    ! Latitude 
     real(wp),allocatable,dimension(:,:) :: &
          o3                     ! Ozone
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
          taupart
     real(wp),allocatable,dimension(:,:,:) :: &
          frac_out,            & ! Cloud fraction
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
          z_vol,               & ! Effective reflectivity factor (mm^6/m^3)
          kr_vol,              & ! Attenuation coefficient hydro (dB/km) 
          g_vol                  ! Attenuation coefficient gases (dB/km)
     real(wp),allocatable,dimension(:,:) :: &
          beta_mol,            & ! Molecular backscatter coefficient
          tau_mol,             & ! Molecular optical depth
          tautot_S_liq,        & ! Liquid water optical thickness, from TOA to SFC
          tautot_S_ice           ! Ice water optical thickness, from TOA to SFC 
     type(radar_cfg) :: rcfg
  end type cosp_optical_inputs
  
  ! ######################################################################################
  ! TYPE cosp_outputs
  ! ######################################################################################
  type cosp_outputs

     ! CALIPSO outputs
     real(wp),dimension(:,:,:),pointer :: &
          calipso_betaperp_tot,  & ! Total backscattered signal
          calipso_beta_tot,      & ! Total backscattered signal
          calipso_tau_tot,       & ! Optical thickness integrated from top to level z
          calipso_lidarcldphase, & ! 3D "lidar" phase cloud fraction 
          calipso_cldlayerphase, & ! low, mid, high-level lidar phase cloud cover
          calipso_lidarcldtmp,   & ! 3D "lidar" phase cloud temperature
          calipso_cfad_sr          ! CFAD of scattering ratio
     real(wp), dimension(:,:),pointer :: &
          calipso_lidarcld,      & ! 3D "lidar" cloud fraction 
          calipso_cldlayer,      & ! low, mid, high-level, total lidar cloud cover
          calipso_beta_mol,      & ! Molecular backscatter
          calipso_temp_tot
     real(wp), dimension(:),pointer :: &
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
  ! SUBROUTINE cosp_simulator
  ! ######################################################################################
  SUBROUTINE COSP_SIMULATOR(cospIN,cospgridIN,cospOUT,start_idx,stop_idx)
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
    
    ! Outputs from the simulators (nested simulator output structure)
    type(cosp_outputs), intent(inout) :: cospOUT
    
    ! Local variables
    integer :: &
         isim,t0,t1,i,icol,nloop,rmod,nprof,nlevels,istart,istop,maxlim,ij,ik,i1,i2
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
         Lrttov_column          ! On/Off switch for column RTTOV simulator (not used)      
    logical :: ok_lidar_cfad = .false.
    
    integer, dimension(:,:),allocatable  :: &
         modisRetrievedPhase,isccpLEVMATCH
    real(wp), dimension(:),  allocatable  :: &
         modisCfTotal,modisCfLiquid,                                            &                         
         modisCfIce, modisCfHigh, modisCfMid, modisCfLow,modisMeanTauTotal,     &       
         modisMeanTauLiquid, modisMeanTauIce, modisMeanLogTauTotal,             &       
         modisMeanLogTauLiquid, modisMeanLogTauIce, modisMeanSizeLiquid,        &        
         modisMeanSizeIce, modisMeanCloudTopPressure, modisMeanLiquidWaterPath, &
         modisMeanIceWaterPath      
    REAL(WP), dimension(:,:),allocatable  :: &
         cloudsatZe_non,cloudsatZe_ray,cloudsatH_atten_to_vol,cloudsatG_atten_to_vol,    &
         cloudsatDBZe, modisRetrievedCloudTopPressure, modisRetrievedTau,                &
         modisRetrievedSize,boxttop,boxtau,boxztop
    REAL(WP), dimension(:,:,:),allocatable :: &
         modisJointHistogram,modisJointHistogramIce,modisJointHistogramLiq
    real(wp),dimension(:),allocatable,target :: &
         out1D_1,out1D_2,out1D_3,out1D_4,out1D_5,out1D_6
    real(wp),dimension(:,:,:),allocatable :: &
       t_in,betamol_in,tmpFlip,betamolFlip,pnormFlip,pnorm_perpFlip,ze_totFlip
    
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
    ! RTTOV subcolumn
    if (associated(cospOUT%rttov_tbs))                                                   &
       Lrttov_subcolumn    = .true.
    ! ISCCP column
    if (associated(cospOUT%isccp_fq)                                       .or.          &
        associated(cospOUT%isccp_meanalbedocld)                            .or.          &
        associated(cospOUT%isccp_meanptop)                                 .or.          &
        associated(cospOUT%isccp_meantaucld)                               .or.          &
        associated(cospOUT%isccp_totalcldarea)                             .or.          &
        associated(cospOUT%isccp_meantb))                                                &
       Lisccp_column    = .true.
    ! MISR column
    if (associated(cospOUT%misr_cldarea)                                   .or.          &
        associated(cospOUT%misr_meanztop)                                  .or.          &
        associated(cospOUT%misr_fq))                                                     &
       Lmisr_column    = .true.
    ! CALIPSO column
    if (associated(cospOUT%calipso_cfad_sr)                                .or.          &
        associated(cospOUT%calipso_lidarcld)                               .or.          &
        associated(cospOUT%calipso_lidarcldphase)                          .or.          &
        associated(cospOUT%calipso_cldlayer)                               .or.          &
        associated(cospOUT%calipso_cldlayerphase)                          .or.          &
        associated(cospOUT%calipso_lidarcldtmp))                                         &
       Lcalipso_column    = .true.
    ! PARASOL column
    if (associated(cospOUT%parasolGrid_refl))                                            &
       Lparasol_column    = .true.
    ! CLOUDSAT column
    if (associated(cospOUT%cloudsat_cfad_ze)                               .or.          &
        associated(cospOUT%lidar_only_freq_cloud)                          .or.          &
        associated(cospOUT%radar_lidar_tcc))                                             &
       Lcloudsat_column    = .true.
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
        associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure))               &
       Lmodis_column    = .true.

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 2b) Error Checking
    !     Enforce bounds on input fields. If input field is out-of-bounds, report error 
    !     and do not run simulator
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    call cosp_errorCheck(cospgridIN,cospIN,Lisccp_subcolumn,Lisccp_column,               &
                         Lmisr_subcolumn,Lmisr_column,Lmodis_subcolumn,Lmodis_column,    &
                         Lcloudsat_subcolumn,Lcloudsat_column,Lcalipso_subcolumn,        &
                         Lcalipso_column,Lrttov_subcolumn,Lrttov_column,                 &
                         Lparasol_subcolumn,Lparasol_column,cospOUT)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 3) Populate instrument simulator inputs
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (Lisccp_subcolumn) then
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
       
       ! Local variables used by the ISCCP simulator
       allocate(isccpLEVMATCH(Npoints,isccpIN%Ncolumns),boxttop(Npoints,isccpIN%Ncolumns))
    endif
    
    if (Lmisr_subcolumn) then
       misrIN%Npoints  => Npoints
       misrIN%Ncolumns => cospIN%Ncolumns
       misrIN%Nlevels  => cospIN%Nlevels
       misrIN%dtau     => cospIN%tau_067
       misrIN%sunlit   => cospgridIN%sunlit
       misrIN%zfull    => cospgridIN%hgt_matrix
       misrIN%at       => cospgridIN%at
       allocate(boxztop(Npoints,misrIN%Ncolumns),boxtau(Npoints,misrIN%Ncolumns))
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
       cloudsatIN%z_vol      => cospIN%z_vol
       cloudsatIN%kr_vol     => cospIN%kr_vol
       cloudsatIN%g_vol      => cospIN%g_vol
       cloudsatIN%rcfg       => cospIN%rcfg
       cloudsatIN%hgt_matrix => cospgridIN%hgt_matrix
       
       ! Local variables used by the CLOUDSAT simulator
       allocate(cloudsatZe_non(cloudsatIN%Npoints,cloudsatIN%Nlevels),                   &
                cloudsatZe_ray(cloudsatIN%Npoints,cloudsatIN%Nlevels),                   &
                cloudsatH_atten_to_vol(cloudsatIN%Npoints,cloudsatIN%Nlevels),           &
                cloudsatG_atten_to_vol(cloudsatIN%Npoints,cloudsatIN%Nlevels),           &
                cloudsatDBZe(cloudsatIN%Npoints,cloudsatIN%Nlevels))
    endif
    
    if (Lrttov_subcolumn) then
       rttovIN%nlevels    => cospIN%Nlevels
       rttovIN%zenang     => cospgridIN%zenang
       rttovIN%co2        => cospgridIN%co2
       rttovIN%ch4        => cospgridIN%ch4
       rttovIN%n2o        => cospgridIN%n2o
       rttovIN%co         => cospgridIN%co
       rttovIN%o3         => cospgridIN%o3
       rttovIN%p          => cospgridIN%pfull
       rttovIN%t          => cospgridIN%at
       rttovIN%q          => cospgridIN%qv
       rttovIN%h_surf     => cospgridIN%hgt_matrix_half(:,1)
       rttovIN%u_surf     => cospgridIN%u_sfc
       rttovIN%v_surf     => cospgridIN%v_sfc
       rttovIN%t_skin     => cospgridIN%skt
       rttovIN%p_surf     => cospgridIN%phalf(:,cospIN%Nlevels+1)
       rttovIN%t_surf     => cospgridIN%t_sfc
       rttovIN%q_surf     => cospgridIN%qv(:,cospIN%Nlevels)
       rttovIN%lsmask     => cospgridIN%land
       rttovIN%latitude   => cospgridIN%lat
       rttovIN%surfem     => cospgridIN%emis_sfc
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
       allocate(modisIN%sunlit(modisIN%Nsunlit),                                         &
                modisIN%notSunlit(count(cospgridIN%sunlit <= 0)),                        &       
                modisIN%pres(modisIN%Nsunlit,cospIN%Nlevels+1))             
       modisIN%sunlit    = pack((/ (i, i = 1, Npoints ) /),                              &
            mask = cospgridIN%sunlit > 0)
       modisIN%notSunlit = pack((/ (i, i = 1, Npoints ) /),                              &
            mask = .not. cospgridIN%sunlit < 0)
       modisIN%pres      = cospgridIN%phalf(int(modisIN%sunlit(:)),:)
       ! Local variables used by the MODIS simulator
       allocate(modisRetrievedTau(modisIN%nSunlit,modisIN%nColumns),                     &
                modisRetrievedSize(modisIN%nSunlit,modisIN%nColumns),                    &
                modisRetrievedPhase(modisIN%nSunlit,modisIN%nColumns),                   &
                modisRetrievedCloudTopPressure(modisIN%nSunlit,modisIN%nColumns),        &
                modisCftotal(modisIN%nSunlit), modisCfLiquid(modisIN%nSunlit),           &
                modisCfIce(modisIN%nSunlit),modisCfHigh(modisIN%nSunlit),                &
                modisCfMid(modisIN%nSunlit),modisCfLow(modisIN%nSunlit),                 &
                modisMeanTauTotal(modisIN%nSunlit),modisMeanTauLiquid(modisIN%nSunlit),  &
                modisMeanTauIce(modisIN%nSunlit),modisMeanLogTauTotal(modisIN%nSunlit),  &       
                modisMeanLogTauLiquid(modisIN%nSunlit),                                  &
                modisMeanLogTauIce(modisIN%nSunlit),modisMeanSizeLiquid(modisIN%nSunlit),&
                modisMeanSizeIce(modisIN%nSunlit),                                       &
                modisMeanCloudTopPressure(modisIN%nSunlit),                              &
                modisMeanLiquidWaterPath(modisIN%nSunlit),                               &
                modisMeanIceWaterPath(modisIN%nSunlit),                                  &
                modisJointHistogram(modisIN%nSunlit,numMODISTauBins,numMODISPresBins),&
                modisJointHistogramIce(modisIN%nSunlit,numModisTauBins,numMODISReffIceBins),&
                modisJointHistogramLiq(modisIN%nSunlit,numModisTauBins,numMODISReffLiqBins))

    endif
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 4) Call subcolumn simulators
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (Lisccp_subcolumn) then
       call icarus_subcolumn(isccpIN%npoints,isccpIN%ncolumns,isccpIN%nlevels,           &
                             isccpIN%sunlit,isccpIN%dtau,isccpIN%dem,isccpIN%skt,        &
                             isccpIN%emsfc_lw,isccpIN%qv,isccpIN%at,isccpIN%pfull,       &
                             isccpIN%phalf,isccpIN%frac_out,isccpLEVMATCH,               &
                             cospOUT%isccp_boxtau(ij:ik,:),                              &
                             cospOUT%isccp_boxptop(ij:ik,:), boxttop(ij:ik,:),           &
                             cospOUT%isccp_meantbclr(ij:ik))
    endif
    
    if (Lmisr_subcolumn) then
       call misr_subcolumn(misrIN%Npoints,misrIN%Ncolumns,misrIN%Nlevels,misrIN%dtau,    &
                       misrIN%zfull,misrIN%at,misrIN%sunlit,boxtau,                      &
                       cospOUT%misr_dist_model_layertops(ij:ik,:),boxztop)
       
    endif
    if (Lcalipso_subcolumn) then
       call lidar_subcolumn(calipsoIN%npoints,calipsoIN%ncolumns,calipsoIN%nlevels,          &
                        calipsoIN%beta_mol,calipsoIN%tau_mol,                            &
                        calipsoIN%betatot,calipsoIN%tautot,calipsoIN%betatot_ice,        &
                        calipsoIN%tautot_ice,calipsoIN%betatot_liq,calipsoIN%tautot_liq, &
                        cospOUT%calipso_beta_mol(ij:ik,:),                               &
                        cospOUT%calipso_beta_tot(ij:ik,:,:),                             &
                        cospOUT%calipso_betaperp_tot(ij:ik,:,:))
    endif
  
    if (Lparasol_subcolumn) then
       do icol=1,parasolIN%Ncolumns
          call parasol_subcolumn(parasolIN%npoints, PARASOL_NREFL,                       &
                             parasolIN%tautot_S_liq(1:parasolIN%Npoints,icol),           &
                             parasolIN%tautot_S_ice(1:parasolIN%Npoints,icol),           &
                             cospOUT%parasolPix_refl(ij:ik,icol,1:PARASOL_NREFL))
       ENDDO
    endif
    
    if (Lcloudsat_subcolumn) then
       do icol=1,cloudsatIN%ncolumns
          call quickbeam_subcolumn(cloudsatIN%rcfg,cloudsatIN%Npoints,cloudsatIN%Nlevels,&
                               cloudsatIN%hgt_matrix/1000._wp,                           &
                               cloudsatIN%z_vol(:,icol,:),cloudsatIN%kr_vol(:,icol,:),   &
                               cloudsatIN%g_vol(:,icol,:),cloudsatH_atten_to_vol,        &
                               cloudsatG_atten_to_vol,cloudsatDBze,cloudsatZe_non,       &
                               cloudsatZe_ray)
          
          ! Store caluculated dBZe values for later output/processing
          cospOUT%cloudsat_Ze_tot(ij:ik,icol,:) = cloudsatDBZe
       enddo
    endif
    
    if (Lmodis_subcolumn) then
       if(modisiN%nSunlit > 0) then 
          do i = 1, modisIN%nSunlit
             call modis_subcolumn(modisIN%Ncolumns,modisIN%Nlevels, modisIN%pres(i,:),   &
                              modisIN%tau(i,:,:),modisIN%liqFrac(i,:,:),modisIN%g(i,:,:),&
                              modisIN%w0(i,:,:),                                         &
                              cospOUT%isccp_boxtau(ij+int(modisIN%sunlit(i))-1,:),       &
                              cospOUT%isccp_boxptop(ij+int(modisIN%sunlit(i))-1,:),      &
                              modisRetrievedPhase(i,:),                                  &
                              modisRetrievedCloudTopPressure(i,:),modisRetrievedTau(i,:),&
                              modisRetrievedSize(i,:))
          end do
       endif
    endif
    
    if (Lrttov_subcolumn) then
        call rttov_subcolumn(rttovIN%surfem,npoints,rttovIN%Nlevels,rttovIN%zenang,      &
                         rttovIN%p/100._wp,rttovIN%t,                                    &
                         (rttovIN%q/(rttovIN%q+0.622_wp*(1._wp - rttovIN%q)))*1e6,       &
                         rttovIN%o3,rttovIN%co2,rttovIN%ch4,rttovIN%n2o,rttovIN%co,      &
                         rttovIN%h_surf,rttovIN%u_surf,rttovIN%v_surf,rttovIN%t_skin,    &
                         rttovIN%p_surf/100._wp,rttovIN%t_surf,rttovIN%q_surf,           &
                         rttovIN%lsmask,rttovIN%latitude,cospOUT%rttov_tbs(ij:ik,:))     

    endif

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
                           cospOUT%isccp_boxtau(ij:ik,:) ,                               &
                           cospOUT%isccp_boxptop(ij:ik,:)/100._wp,                       &
                           isccpIN%sunlit,isccpIN%pfull,isccpIN%phalf,isccpIN%qv,        &
                           isccpIN%at,isccpIN%skt,isccpIN%emsfc_lw,boxttop,              &
                           cospOUT%isccp_fq(ij:ik,:,:),                                  &
                           cospOUT%isccp_meanalbedocld(ij:ik),                           &
                           cospOUT%isccp_meanptop(ij:ik),cospOUT%isccp_meantaucld(ij:ik),&
                           cospOUT%isccp_totalcldarea(ij:ik),cospOUT%isccp_meantb(ij:ik))
       cospOUT%isccp_fq(ij:ik,:,:) = cospOUT%isccp_fq(ij:ik,:,7:1:-1)
       
       ! Check if there is any value slightly greater than 1
       where ((cospOUT%isccp_totalcldarea > 1.0-1.e-5) .and.                             &
              (cospOUT%isccp_totalcldarea < 1.0+1.e-5))
              cospOUT%isccp_totalcldarea = 1.0
       endwhere
       
       ! Clear up memory (if necessary)
       deallocate(boxttop)
       if (allocated(out1D_1)) deallocate(out1D_1)
       if (allocated(out1D_2)) deallocate(out1D_2)
       if (allocated(out1D_3)) deallocate(out1D_3)
       if (allocated(out1D_4)) deallocate(out1D_4)
       if (allocated(out1D_5)) deallocate(out1D_5)
       if (allocated(out1D_6)) deallocate(out1D_6)
    endif
    
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
       call misr_column(misrIN%Npoints,misrIN%Ncolumns,misrIN%Nlevels,boxztop,           &
                         misrIN%sunlit,boxtau,cospOUT%misr_cldarea(ij:ik),               &
                         cospOUT%misr_meanztop(ij:ik),cospOUT%misr_fq(ij:ik,:,:))              

       ! Clear up memory
       deallocate(boxtau,boxztop)   
       if (allocated(out1D_1)) deallocate(out1D_1)
       if (allocated(out1D_2)) deallocate(out1D_2)
       if (allocated(out1D_3)) deallocate(out1D_3)
    endif
    
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

       ! Call simulator
       ok_lidar_cfad=.true. 
       call lidar_column(calipsoIN%Npoints,calipsoIN%Ncolumns,calipsoIN%Nlevels,         &
                         Nlvgrid,SR_BINS,cospgridIN%at(:,:),                             &
                         cospOUT%calipso_beta_tot(ij:ik,:,:),                            &
                         cospOUT%calipso_betaperp_tot(ij:ik,:,:),                        &
                         cospOUT%calipso_beta_mol(ij:ik,:),cospgridIN%land,              &
                         cospgridIN%phalf(:,2:calipsoIN%Nlevels),ok_lidar_cfad,          &
                         LIDAR_NCAT,cospOUT%calipso_cfad_sr(ij:ik,:,:),                  &
                         cospOUT%calipso_lidarcld(ij:ik,:),                              &
                         cospOUT%calipso_lidarcldphase(ij:ik,:,:),                       &
                         cospOUT%calipso_cldlayer(ij:ik,:),                              &
                         cospgridIN%hgt_matrix,cospgridIN%hgt_matrix_half,               &
                         cospOUT%calipso_cldlayerphase(ij:ik,:,:),                       &
                         cospOUT%calipso_lidarcldtmp(ij:ik,:,:))                                      
       cospOUT%calipso_srbval = calipso_histBsct     
      
      ! Free up memory (if necessary)
      if (allocated(out1D_1)) deallocate(out1D_1)
      if (allocated(out1D_2)) deallocate(out1D_2)
      if (allocated(out1D_3)) deallocate(out1D_3)
      if (allocated(out1D_4)) deallocate(out1D_4)
      if (allocated(out1D_5)) deallocate(out1D_5)
      if (allocated(out1D_6)) deallocate(out1D_6)
    endif

    ! PARASOL
    if (Lparasol_column) then
       call parasol_column(parasolIN%Npoints,PARASOL_NREFL,parasolIN%Ncolumns,           &
                            cospgridIN%land(:),cospOUT%parasolPix_refl(ij:ik,:,:),       &
                            cospOUT%parasolGrid_refl(ij:ik,:))
    endif
    
    ! CLOUDSAT
    if (Lcloudsat_column) then
       ! Check to see which outputs are requested. If not requested, use a local dummy array
       if (.not. associated(cospOUT%cloudsat_cfad_ze)) then
          allocate(out1D_1(Npoints*DBZE_BINS*Nlvgrid))
          cospOUT%cloudsat_cfad_ze(ij:ik,1:DBZE_BINS,1:Nlvgrid) => out1D_1
       endif

       ! Call simulator
       call quickbeam_column(cloudsatIN%Npoints,cloudsatIN%Ncolumns,cloudsatIN%Nlevels,  &
                              Nlvgrid,cospOUT%cloudsat_Ze_tot(ij:ik,:,:),                &
                              cospgridIN%hgt_matrix,cospgridIN%hgt_matrix_half,          &
                              cospOUT%cloudsat_cfad_ze(ij:ik,:,:))
       ! Free up memory  (if necessary)
       if (allocated(out1D_1)) deallocate(out1D_1)
    endif
    
    ! MODIS
    if (Lmodis_column) then
       if(modisiN%nSunlit > 0) then 
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
          ! Put results into COSP output structure
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
                  int(modisIN%sunlit(:))-1, 2:numModisTauBins+1, :) = modisJointHistogram(:, :, :)           
             ! Reorder pressure bins in joint histogram to go from surface to TOA 
             cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure(ij:ik,:,:) = &
                  cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure(ij:ik,:,numMODISPresBins:1:-1)
          endif
          if (associated(cospOUT%modis_Optical_Thickness_vs_ReffIce)) then
             cospOUT%modis_Optical_Thickness_vs_ReffIce(ij+int(modisIN%sunlit(:))-1, 2:numMODISTauBins+1,:) = &
                modisJointHistogramIce(:,:,:)
          endif
          if (associated(cospOUT%modis_Optical_Thickness_vs_ReffLiq)) then
             cospOUT%modis_Optical_Thickness_vs_ReffLiq(ij+int(modisIN%sunlit(:))-1, 2:numMODISTauBins+1,:) = &
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
       deallocate(modisRetrievedTau,modisRetrievedSize,modisRetrievedPhase,              &
                  modisRetrievedCloudTopPressure,modisCftotal,modisCfLiquid,             &
                  modisCfIce,modisCfHigh,modisCfMid,modisCfLow,modisMeanTauTotal,        &
                  modisMeanTauLiquid,modisMeanTauIce,modisMeanLogTauTotal,               &       
                  modisMeanLogTauLiquid,modisMeanLogTauIce,modisMeanSizeLiquid,          & 
                  modisMeanSizeIce,modisMeanCloudTopPressure,modisMeanLiquidWaterPath,   &
                  modisMeanIceWaterPath,modisJointHistogram,modisJointHistogramIce,      &
                  modisJointHistogramLiq)       
    endif
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 6) Compute multi-instrument products
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! CLOUDSAT/CALIPSO products
    if (associated(cospOUT%radar_lidar_tcc) .or. associated(cospOUT%lidar_only_freq_cloud)) then    
       if (use_vgrid) then
          allocate(t_in(cloudsatIN%Npoints,1,cloudsatIN%Nlevels),                        &
                   betamol_in(cloudsatIN%Npoints,1,cloudsatIN%Nlevels),                  &
                   tmpFlip(cloudsatIN%Npoints,1,Nlvgrid),                                &
                   betamolFlip(cloudsatIN%Npoints,1,Nlvgrid),                            &
                   pnormFlip(cloudsatIN%Npoints,cloudsatIN%Ncolumns,Nlvgrid),            &
                   pnorm_perpFlip(cloudsatIN%Npoints,cloudsatIN%Ncolumns,Nlvgrid),       &
                   Ze_totFlip(cloudsatIN%Npoints,cloudsatIN%Ncolumns,Nlvgrid))
          t_in(:,1,:)=cospgridIN%at(:,:)
          call cosp_change_vertical_grid(cloudsatIN%Npoints,1,cloudsatIN%Nlevels,        &
                                         cospgridIN%hgt_matrix,                          &
                                         cospgridIN%hgt_matrix_half,t_in,Nlvgrid,        &
                                         vgrid_zl,vgrid_zu,tmpFlip)
          betamol_in(:,1,:) = cospOUT%calipso_beta_mol(ij:ik,:)
          call cosp_change_vertical_grid(cloudsatIN%Npoints,1,cloudsatIN%Nlevels,        &
                                         cospgridIN%hgt_matrix,                          &
                                         cospgridIN%hgt_matrix_half,betamol_in,          &
                                         Nlvgrid,vgrid_zl,vgrid_zu,betamolFlip)
          call cosp_change_vertical_grid(cloudsatIN%Npoints,cloudsatIN%Ncolumns,         &
                                         cloudsatIN%Nlevels,cospgridIN%hgt_matrix,       &
                                         cospgridIN%hgt_matrix_half,                     &
                                         cospOUT%calipso_beta_tot(ij:ik,:,:),            &
                                         Nlvgrid,vgrid_zl,vgrid_zu,pnormFlip)
          call cosp_change_vertical_grid(cloudsatIN%Npoints,cloudsatIN%Ncolumns,         &
                                         cloudsatIN%Nlevels,cospgridIN%hgt_matrix,       &
                                         cospgridIN%hgt_matrix_half,                     &
                                         cospOUT%calipso_betaperp_tot(ij:ik,:,:),        &
                                         Nlvgrid,vgrid_zl,vgrid_zu,pnorm_perpFlip)      
          call cosp_change_vertical_grid(cloudsatIN%Npoints,cloudsatIN%Ncolumns,         &
                                         cloudsatIN%Nlevels,cospgridIN%hgt_matrix,       &
                                         cospgridIN%hgt_matrix_half,                     &
                                         cospOUT%cloudsat_Ze_tot(ij:ik,:,:),             &
                                         Nlvgrid,vgrid_zl,vgrid_zu,                      &
                                         Ze_totFlip,log_units=.true.)                               
                                                           
          call cosp_lidar_only_cloud(cloudsatIN%Npoints,cloudsatIN%Ncolumns,             &
                                     Nlvgrid,tmpFlip,pnormFlip,pnorm_perpFlip,           &
                                     betamolFlip,Ze_totFlip,                             &
                                     cospOUT%lidar_only_freq_cloud(ij:ik,:),             &
                                     cospOUT%radar_lidar_tcc(ij:ik))                                            
          deallocate(t_in,betamol_in,tmpFlip,betamolFlip,pnormFlip,pnorm_perpFlip,       &
                     ze_totFlip)
       else
          call cosp_lidar_only_cloud(cloudsatIN%Npoints,cloudsatIN%Ncolumns,             &
                                     cospIN%Nlevels,cospgridIN%at,                       &
                                     cospOUT%calipso_beta_tot(ij:ik,:,:),                &
                                     cospOUT%calipso_betaperp_tot(ij:ik,:,:),            &
                                     cospOUT%calipso_beta_mol(ij:ik,:),                  &
                                     cospOUT%cloudsat_Ze_tot(ij:ik,:,:),                 &
                                     cospOUT%lidar_only_freq_cloud(ij:ik,:),             &
                                     cospOUT%radar_lidar_tcc(ij:ik))
       endif
    endif

  end SUBROUTINE COSP_SIMULATOR
  ! ######################################################################################
  ! SUBROUTINE cosp_init
  ! ######################################################################################
  SUBROUTINE COSP_INIT(Npoints,Nlevels,cloudsat_radar_freq,cloudsat_k2,                  &
                       cloudsat_use_gas_abs,cloudsat_do_ray,isccp_top_height,            &
                       isccp_top_height_direction,hgt_matrix,hgt_matrix_half,            &
                       surface_radar,rcfg,rttov_Nchannels,rttov_Channels,rttov_platform, &
                       rttov_satellite,rttov_instrument,lusevgrid,luseCSATvgrid,Nvgrid,  &
                       cloudsat_micro_scheme)
    
    ! INPUTS
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
    real(wp),intent(in),dimension(Npoints,Nlevels) :: &
         hgt_matrix                    !
    real(wp),intent(in),dimension(Npoints,Nlevels+1) :: &
         hgt_matrix_half     
    logical,intent(in) :: &
         lusevgrid,                  & ! Switch to use different vertical grid
         luseCSATvgrid                 ! Switch to use CLOUDSAT grid spacing for new  
                                       ! vertical grid
    character(len=64),intent(in) :: &
       cloudsat_micro_scheme           ! Microphysical scheme used by CLOUDSAT
    
    ! OUTPUTS
    type(radar_cfg) :: rcfg
  
    ! Local variables
    integer  :: i
    real(wp) :: zstep

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
    endif
    
    ! Define model grid 
    allocate(mgrid_z(Nlevels),mgrid_zl(Nlevels),mgrid_zu(Nlevels))
    mgrid_z             = hgt_matrix(1,:)
    mgrid_zl            = hgt_matrix_half(1,:)
    mgrid_zu(2:Nlevels) = hgt_matrix_half(1,:)
    mgrid_zu(1)         = hgt_matrix(1,1)+(hgt_matrix(1,1)-mgrid_zl(1))    
    
    ! Initialize ISCCP
    call cosp_isccp_init(isccp_top_height,isccp_top_height_direction)
    ! Initialize MODIS
    call cosp_modis_init()
    ! Initialize MISR
    call cosp_misr_init()
    ! Initialize RTTOV
    call cosp_rttov_init(rttov_Nchannels,rttov_platform,rttov_satellite,rttov_instrument,&
         rttov_channels)
    ! Initialize quickbeam_optics
    call quickbeam_optics_init()
    ! Initialize radar
    call cosp_cloudsat_init(cloudsat_radar_freq,cloudsat_k2,cloudsat_use_gas_abs,        &
         cloudsat_do_ray,R_UNDEF,N_HYDRO,Npoints,Nlevels,hgt_matrix,surface_radar,rcfg,  &
         cloudsat_micro_scheme)
    ! Initialize lidar
    call cosp_calipso_init()
    ! Initialize PARASOL
    call cosp_parasol_init()
    
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
             y%taupart(npoints,ncolumns,nlevels,4),  y%z_vol(npoints,Ncolumns,nlevels),  &
             y%kr_vol(npoints,Ncolumns,nlevels),     y%g_vol(npoints,Ncolumns,nlevels),  &
             y%asym(npoints,nColumns,nlevels),       y%ss_alb(npoints,nColumns,nlevels))
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
             y%v_sfc(npoints),y%t_sfc(npoints),y%lat(npoints),y%emis_sfc(nchan),         &
             y%hgt_matrix_half(npoints,nlevels))
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
        allocate(x%modis_Optical_Thickness_vs_Cloud_Top_Pressure(nPoints,numModisTauBins+1,numMODISPresBins))
        allocate(x%modis_Optical_thickness_vs_ReffLIQ(nPoints,numMODISTauBins+1,numMODISReffLiqBins))   
        allocate(x%modis_Optical_Thickness_vs_ReffICE(nPoints,numMODISTauBins+1,numMODISReffIceBins))
    endif
    
    ! LIDAR simulator
    if (LlidarBetaMol532) allocate(x%calipso_beta_mol(Npoints,Nlevels))
    if (Latb532)          allocate(x%calipso_beta_tot(Npoints,Ncolumns,Nlevels))
    if (LcfadLidarsr532)  then
        allocate(x%calipso_srbval(SR_BINS))
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
    if (Lclcalipsotmp .or. Lclcalipsotmpliq .or. Lclcalipsoice .or. Lclcalipsotmpun) then
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
               y%tautot_S_liq,y%tautot_S_ice,y%z_vol,y%kr_vol,y%g_vol,y%asym,y%ss_alb,   &
               y%fracLiq,y%taupart)

  end subroutine destroy_cospIN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE destroy_cospstateIN     
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  subroutine destroy_cospstateIN(y)
    type(cosp_column_inputs),intent(inout) :: y
    deallocate(y%sunlit,y%skt,y%land,y%at,y%pfull,y%phalf,y%qv,y%o3,y%hgt_matrix,       &
         y%u_sfc,y%v_sfc,y%t_sfc,y%lat,y%emis_sfc,y%hgt_matrix_half)
  end subroutine destroy_cospstateIN
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE destroy_cosp_outputs
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  subroutine destroy_cosp_outputs(y)
     type(cosp_outputs),intent(inout) :: y

     ! Deallocate
     if (associated(y%calipso_beta_mol))          deallocate(y%calipso_beta_mol)     
     if (associated(y%calipso_temp_tot))          deallocate(y%calipso_temp_tot)     
     if (associated(y%calipso_betaperp_tot))      deallocate(y%calipso_betaperp_tot)     
     if (associated(y%calipso_beta_tot))          deallocate(y%calipso_beta_tot)     
     if (associated(y%calipso_tau_tot))           deallocate(y%calipso_tau_tot)     
     if (associated(y%calipso_lidarcldphase))     deallocate(y%calipso_lidarcldphase)     
     if (associated(y%calipso_cldlayerphase))     deallocate(y%calipso_cldlayerphase)     
     if (associated(y%calipso_lidarcldtmp))       deallocate(y%calipso_lidarcldtmp)     
     if (associated(y%calipso_cldlayer))          deallocate(y%calipso_cldlayer)     
     if (associated(y%calipso_lidarcld))          deallocate(y%calipso_lidarcld)     
     if (associated(y%calipso_srbval))            deallocate(y%calipso_srbval)     
     if (associated(y%calipso_cfad_sr))           deallocate(y%calipso_cfad_sr)     
     if (associated(y%parasolPix_refl))           deallocate(y%parasolPix_refl)     
     if (associated(y%parasolGrid_refl))          deallocate(y%parasolGrid_refl)     
     if (associated(y%cloudsat_Ze_tot))           deallocate(y%cloudsat_Ze_tot)  
     if (associated(y%cloudsat_cfad_ze))          deallocate(y%cloudsat_cfad_ze)     
     if (associated(y%radar_lidar_tcc))           deallocate(y%radar_lidar_tcc)  
     if (associated(y%lidar_only_freq_cloud))     deallocate(y%lidar_only_freq_cloud)     
     if (associated(y%isccp_totalcldarea))        deallocate(y%isccp_totalcldarea)  
     if (associated(y%isccp_meantb))              deallocate(y%isccp_meantb)     
     if (associated(y%isccp_meantbclr))           deallocate(y%isccp_meantbclr)  
     if (associated(y%isccp_meanptop))            deallocate(y%isccp_meanptop)     
     if (associated(y%isccp_meantaucld))          deallocate(y%isccp_meantaucld)       
     if (associated(y%isccp_meanalbedocld))       deallocate(y%isccp_meanalbedocld)     
     if (associated(y%isccp_boxtau))              deallocate(y%isccp_boxtau)       
     if (associated(y%isccp_boxptop))             deallocate(y%isccp_boxptop)     
     if (associated(y%isccp_fq))                  deallocate(y%isccp_fq)       
     if (associated(y%misr_fq))                   deallocate(y%misr_fq)     
     if (associated(y%misr_dist_model_layertops)) deallocate(y%misr_dist_model_layertops)       
     if (associated(y%misr_meanztop))             deallocate(y%misr_meanztop)     
     if (associated(y%misr_cldarea))              deallocate(y%misr_cldarea)      
     if (associated(y%rttov_tbs))                 deallocate(y%rttov_tbs)     
     if (associated(y%modis_Cloud_Fraction_Total_Mean))                                  &
        deallocate(y%modis_Cloud_Fraction_Total_Mean)       
     if (associated(y%modis_Cloud_Fraction_Ice_Mean))                                    &
        deallocate(y%modis_Cloud_Fraction_Ice_Mean)     
     if (associated(y%modis_Cloud_Fraction_Water_Mean))                                  &
        deallocate(y%modis_Cloud_Fraction_Water_Mean)           
     if (associated(y%modis_Cloud_Fraction_High_Mean))                                   &
        deallocate(y%modis_Cloud_Fraction_High_Mean)     
     if (associated(y%modis_Cloud_Fraction_Mid_Mean))                                    &
        deallocate(y%modis_Cloud_Fraction_Mid_Mean)       
     if (associated(y%modis_Cloud_Fraction_Low_Mean))                                    &
        deallocate(y%modis_Cloud_Fraction_Low_Mean)     
     if (associated(y%modis_Optical_Thickness_Total_Mean))                               &
        deallocate(y%modis_Optical_Thickness_Total_Mean)  
     if (associated(y%modis_Optical_Thickness_Water_Mean))                               &
        deallocate(y%modis_Optical_Thickness_Water_Mean)     
     if (associated(y%modis_Optical_Thickness_Ice_Mean))                                 &
        deallocate(y%modis_Optical_Thickness_Ice_Mean)       
     if (associated(y%modis_Optical_Thickness_Total_LogMean))                            &
        deallocate(y%modis_Optical_Thickness_Total_LogMean)    
     if (associated(y%modis_Optical_Thickness_Water_LogMean))                            &
        deallocate(y%modis_Optical_Thickness_Water_LogMean)     
     if (associated(y%modis_Optical_Thickness_Ice_LogMean))                              &
        deallocate(y%modis_Optical_Thickness_Ice_LogMean)     
     if (associated(y%modis_Cloud_Particle_Size_Water_Mean))                             &
        deallocate(y%modis_Cloud_Particle_Size_Water_Mean)       
     if (associated(y%modis_Cloud_Particle_Size_Ice_Mean))                               &
        deallocate(y%modis_Cloud_Particle_Size_Ice_Mean)     
     if (associated(y%modis_Cloud_Top_Pressure_Total_Mean))                              &
        deallocate(y%modis_Cloud_Top_Pressure_Total_Mean)           
     if (associated(y%modis_Liquid_Water_Path_Mean))                                     &
        deallocate(y%modis_Liquid_Water_Path_Mean)     
     if (associated(y%modis_Ice_Water_Path_Mean))                                        &
        deallocate(y%modis_Ice_Water_Path_Mean)       
     if (associated(y%modis_Optical_Thickness_vs_Cloud_Top_Pressure))                    &
        deallocate(y%modis_Optical_Thickness_vs_Cloud_Top_Pressure)     
     if (associated(y%modis_Optical_thickness_vs_ReffLIQ))                               &
        deallocate(y%modis_Optical_thickness_vs_ReffLIQ)
     if (associated(y%modis_Optical_thickness_vs_ReffICE))                               &
        deallocate(y%modis_Optical_thickness_vs_ReffICE)
        
  end subroutine destroy_cosp_outputs
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_errorCheck
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine cosp_errorCheck(cospgridIN,cospIN,Lisccp_subcolumn,Lisccp_column,Lmisr_subcolumn,Lmisr_column,    &
                             Lmodis_subcolumn,Lmodis_column,Lcloudsat_subcolumn,Lcloudsat_column,Lcalipso_subcolumn,  &
                             Lcalipso_column,Lrttov_subcolumn,Lrttov_column,Lparasol_subcolumn,Lparasol_column,    &
                             cospOUT)
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
      Lrttov_column          ! RTTOV column simulator on/off switch       
  type(cosp_outputs),intent(inout) :: &
      cospOUT                ! COSP Outputs   
  
  ! Local variables
  character(len=100) :: parasolErrorMessage    
                   
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! PART 1: Check input array values for out-of-bounds values. When an out-of-bound value
  !         is encountered, COSP outputs that are dependent on that input are filled with
  !         an undefined value (set in cosp_config.f90) and if necessary, that simulator 
  !         is turned off.
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (any(cospgridIN%sunlit .lt. 0)) then
       call errorMessage('ERROR: COSP input variable: cospgridIN%sunlit contains values out of range')
       Lisccp_subcolumn = .false.
       Lisccp_column    = .false.
       Lmisr_subcolumn  = .false.
       Lmisr_column     = .false.
       Lmodis_subcolumn = .false.
       Lmodis_column    = .false.
       cospOUT%isccp_totalcldarea(:)                                = R_UNDEF
       cospOUT%isccp_meantb(:)                                      = R_UNDEF
       cospOUT%isccp_meantbclr(:)                                   = R_UNDEF
       cospOUT%isccp_meanptop(:)                                    = R_UNDEF
       cospOUT%isccp_meantaucld(:)                                  = R_UNDEF
       cospOUT%isccp_meanalbedocld(:)                               = R_UNDEF
       cospOUT%isccp_boxtau(:,:)                                    = R_UNDEF
       cospOUT%isccp_boxptop(:,:)                                   = R_UNDEF
       cospOUT%isccp_fq(:,:,:)                                      = R_UNDEF
       cospOUT%misr_fq(:,:,:)                                       = R_UNDEF
       cospOUT%misr_dist_model_layertops(:,:)                       = R_UNDEF
       cospOUT%misr_meanztop(:)                                     = R_UNDEF
       cospOUT%misr_cldarea(:)                                      = R_UNDEF
       cospOUT%modis_Cloud_Fraction_Total_Mean(:)                   = R_UNDEF
       cospOUT%modis_Cloud_Fraction_Water_Mean(:)                   = R_UNDEF
       cospOUT%modis_Cloud_Fraction_Ice_Mean(:)                     = R_UNDEF
       cospOUT%modis_Cloud_Fraction_High_Mean(:)                    = R_UNDEF
       cospOUT%modis_Cloud_Fraction_Mid_Mean(:)                     = R_UNDEF
       cospOUT%modis_Cloud_Fraction_Low_Mean(:)                     = R_UNDEF
       cospOUT%modis_Optical_Thickness_Total_Mean(:)                = R_UNDEF
       cospOUT%modis_Optical_Thickness_Water_Mean(:)                = R_UNDEF
       cospOUT%modis_Optical_Thickness_Ice_Mean(:)                  = R_UNDEF
       cospOUT%modis_Optical_Thickness_Total_LogMean(:)             = R_UNDEF
       cospOUT%modis_Optical_Thickness_Water_LogMean(:)             = R_UNDEF
       cospOUT%modis_Optical_Thickness_Ice_LogMean(:)               = R_UNDEF
       cospOUT%modis_Cloud_Particle_Size_Water_Mean(:)              = R_UNDEF
       cospOUT%modis_Cloud_Particle_Size_Ice_Mean(:)                = R_UNDEF
       cospOUT%modis_Cloud_Top_Pressure_Total_Mean(:)               = R_UNDEF
       cospOUT%modis_Liquid_Water_Path_Mean(:)                      = R_UNDEF
       cospOUT%modis_Ice_Water_Path_Mean(:)                         = R_UNDEF
       cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure(:,:,:) = R_UNDEF
       cospOUT%modis_Optical_Thickness_vs_ReffICE(:,:,:)            = R_UNDEF
       cospOUT%modis_Optical_Thickness_vs_ReffLIQ(:,:,:)            = R_UNDEF
    endif
    if (any(cospgridIN%at .lt. 0)) then   
       call errorMessage('ERROR: COSP input variable: cospgridIN%at contains values out of range')
       Lisccp_subcolumn = .false.
       Lisccp_column    = .false.
       Lmisr_subcolumn  = .false.
       Lmisr_column     = .false.
       Lrttov_subcolumn = .false.
       Lcalipso_column  = .false.
       Lcloudsat_column = .false.
       cospOUT%rttov_tbs(:,:) = R_UNDEF
       cospOUT%isccp_totalcldarea(:)          = R_UNDEF
       cospOUT%isccp_meantb(:)                = R_UNDEF
       cospOUT%isccp_meantbclr(:)             = R_UNDEF
       cospOUT%isccp_meanptop(:)              = R_UNDEF
       cospOUT%isccp_meantaucld(:)            = R_UNDEF
       cospOUT%isccp_meanalbedocld(:)         = R_UNDEF
       cospOUT%isccp_boxtau(:,:)              = R_UNDEF    
       cospOUT%isccp_boxptop(:,:)             = R_UNDEF     
       cospOUT%isccp_fq(:,:,:)                = R_UNDEF       
       cospOUT%misr_fq(:,:,:)                 = R_UNDEF 
       cospOUT%misr_dist_model_layertops(:,:) = R_UNDEF 
       cospOUT%misr_meanztop(:)               = R_UNDEF 
       cospOUT%misr_cldarea(:)                = R_UNDEF
       cospOUT%calipso_cfad_sr(:,:,:)         = R_UNDEF
       cospOUT%calipso_lidarcld(:,:)          = R_UNDEF
       cospOUT%calipso_lidarcldphase(:,:,:)   = R_UNDEF
       cospOUT%calipso_cldlayer(:,:)          = R_UNDEF
       cospOUT%calipso_cldlayerphase(:,:,:)   = R_UNDEF
       cospOUT%calipso_lidarcldtmp(:,:,:)     = R_UNDEF            
       cospOUT%cloudsat_cfad_ze(:,:,:)        = R_UNDEF
       cospOUT%lidar_only_freq_cloud(:,:)     = R_UNDEF
       cospOUT%radar_lidar_tcc(:)             = R_UNDEF       
    endif
    if (any(cospgridIN%pfull .lt. 0)) then
       call errorMessage('ERROR: COSP input variable: cospgridIN%pfull contains values out of range')
       Lisccp_subcolumn = .false.
       Lisccp_column    = .false.     
       Lrttov_subcolumn = .false.
       cospOUT%rttov_tbs(:,:) = R_UNDEF
       cospOUT%isccp_totalcldarea(:)          = R_UNDEF
       cospOUT%isccp_meantb(:)                = R_UNDEF
       cospOUT%isccp_meantbclr(:)             = R_UNDEF
       cospOUT%isccp_meanptop(:)              = R_UNDEF
       cospOUT%isccp_meantaucld(:)            = R_UNDEF
       cospOUT%isccp_meanalbedocld(:)         = R_UNDEF
       cospOUT%isccp_boxtau(:,:)              = R_UNDEF    
       cospOUT%isccp_boxptop(:,:)             = R_UNDEF     
       cospOUT%isccp_fq(:,:,:)                = R_UNDEF         
    endif
    if (any(cospgridIN%phalf .lt. 0)) then
       call errorMessage('ERROR: COSP input variable: cospgridIN%phalf contains values out of range')
       Lisccp_subcolumn = .false.
       Lisccp_column    = .false.     
       Lrttov_subcolumn = .false.
       Lmodis_subcolumn = .false.
       Lmodis_column    = .false.
       Lcalipso_column  = .false.
       cospOUT%rttov_tbs(:,:) = R_UNDEF
       cospOUT%isccp_totalcldarea(:)                                = R_UNDEF
       cospOUT%isccp_meantb(:)                                      = R_UNDEF
       cospOUT%isccp_meantbclr(:)                                   = R_UNDEF
       cospOUT%isccp_meanptop(:)                                    = R_UNDEF
       cospOUT%isccp_meantaucld(:)                                  = R_UNDEF
       cospOUT%isccp_meanalbedocld(:)                               = R_UNDEF
       cospOUT%isccp_boxtau(:,:)                                    = R_UNDEF    
       cospOUT%isccp_boxptop(:,:)                                   = R_UNDEF     
       cospOUT%isccp_fq(:,:,:)                                      = R_UNDEF        
       cospOUT%modis_Cloud_Fraction_Total_Mean(:)                   = R_UNDEF
       cospOUT%modis_Cloud_Fraction_Water_Mean(:)                   = R_UNDEF
       cospOUT%modis_Cloud_Fraction_Ice_Mean(:)                     = R_UNDEF
       cospOUT%modis_Cloud_Fraction_High_Mean(:)                    = R_UNDEF
       cospOUT%modis_Cloud_Fraction_Mid_Mean(:)                     = R_UNDEF
       cospOUT%modis_Cloud_Fraction_Low_Mean(:)                     = R_UNDEF
       cospOUT%modis_Optical_Thickness_Total_Mean(:)                = R_UNDEF
       cospOUT%modis_Optical_Thickness_Water_Mean(:)                = R_UNDEF
       cospOUT%modis_Optical_Thickness_Ice_Mean(:)                  = R_UNDEF
       cospOUT%modis_Optical_Thickness_Total_LogMean(:)             = R_UNDEF 
       cospOUT%modis_Optical_Thickness_Water_LogMean(:)             = R_UNDEF
       cospOUT%modis_Optical_Thickness_Ice_LogMean(:)               = R_UNDEF
       cospOUT%modis_Cloud_Particle_Size_Water_Mean(:)              = R_UNDEF
       cospOUT%modis_Cloud_Particle_Size_Ice_Mean(:)                = R_UNDEF
       cospOUT%modis_Cloud_Top_Pressure_Total_Mean(:)               = R_UNDEF
       cospOUT%modis_Liquid_Water_Path_Mean(:)                      = R_UNDEF
       cospOUT%modis_Ice_Water_Path_Mean(:)                         = R_UNDEF 
       cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure(:,:,:) = R_UNDEF
       cospOUT%modis_Optical_Thickness_vs_ReffICE(:,:,:)            = R_UNDEF
       cospOUT%modis_Optical_Thickness_vs_ReffLIQ(:,:,:)            = R_UNDEF
       cospOUT%calipso_cfad_sr(:,:,:)                               = R_UNDEF
       cospOUT%calipso_lidarcld(:,:)                                = R_UNDEF
       cospOUT%calipso_lidarcldphase(:,:,:)                         = R_UNDEF
       cospOUT%calipso_cldlayer(:,:)                                = R_UNDEF
       cospOUT%calipso_cldlayerphase(:,:,:)                         = R_UNDEF
       cospOUT%calipso_lidarcldtmp(:,:,:)                           = R_UNDEF
    endif
    if (any(cospgridIN%qv .lt. 0)) then
       call errorMessage('ERROR: COSP input variable: cospgridIN%qv contains values out of range')
       Lisccp_subcolumn = .false.
       Lisccp_column    = .false.     
       Lrttov_subcolumn = .false.
       cospOUT%rttov_tbs(:,:)         = R_UNDEF
       cospOUT%isccp_totalcldarea(:)  = R_UNDEF
       cospOUT%isccp_meantb(:)        = R_UNDEF
       cospOUT%isccp_meantbclr(:)     = R_UNDEF
       cospOUT%isccp_meanptop(:)      = R_UNDEF
       cospOUT%isccp_meantaucld(:)    = R_UNDEF
       cospOUT%isccp_meanalbedocld(:) = R_UNDEF
       cospOUT%isccp_boxtau(:,:)      = R_UNDEF    
       cospOUT%isccp_boxptop(:,:)     = R_UNDEF     
       cospOUT%isccp_fq(:,:,:)        = R_UNDEF         
    endif
    if (any(cospgridIN%hgt_matrix .lt. 0)) then
       call errorMessage('ERROR: COSP input variable: cospgridIN%hgt_matrix contains values out of range')
       Lmisr_subcolumn     = .false.
       Lmisr_column        = .false.
       Lcloudsat_subcolumn = .false.
       Lcloudsat_column    = .false.
       Lcalipso_column     = .false.
       cospOUT%misr_fq(:,:,:)                 = R_UNDEF 
       cospOUT%misr_dist_model_layertops(:,:) = R_UNDEF 
       cospOUT%misr_meanztop(:)               = R_UNDEF 
       cospOUT%misr_cldarea(:)                = R_UNDEF       
       cospOUT%calipso_cfad_sr(:,:,:)         = R_UNDEF
       cospOUT%calipso_lidarcld(:,:)          = R_UNDEF
       cospOUT%calipso_lidarcldphase(:,:,:)   = R_UNDEF
       cospOUT%calipso_cldlayer(:,:)          = R_UNDEF
       cospOUT%calipso_cldlayerphase(:,:,:)   = R_UNDEF
       cospOUT%calipso_lidarcldtmp(:,:,:)     = R_UNDEF  
       cospOUT%cloudsat_Ze_tot(:,:,:)         = R_UNDEF
       cospOUT%cloudsat_cfad_ze(:,:,:)        = R_UNDEF
       cospOUT%lidar_only_freq_cloud(:,:)     = R_UNDEF
       cospOUT%radar_lidar_tcc(:)             = R_UNDEF                  
    endif
    if (any(cospgridIN%hgt_matrix_half .lt. 0)) then
       call errorMessage('ERROR: COSP input variable: cospgridIN%hgt_matrix_half contains values out of range')
       Lrttov_subcolumn = .false.
       Lcloudsat_column = .false.
       Lcalipso_column  = .false.       
       cospOUT%rttov_tbs(:,:)               = R_UNDEF
       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF
       cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
       cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
       cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
       cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
       cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF  
       cospOUT%cloudsat_cfad_ze(:,:,:)      = R_UNDEF
       cospOUT%lidar_only_freq_cloud(:,:)   = R_UNDEF
       cospOUT%radar_lidar_tcc(:)           = R_UNDEF                
    endif
    if (any(cospgridIN%land .lt. 0)) then
       call errorMessage('ERROR: COSP input variable: cospgridIN%land contains values out of range')
       Lrttov_subcolumn = .false.
       Lcalipso_column  = .false.       
       Lparasol_column  = .false.
       cospOUT%rttov_tbs(:,:)               = R_UNDEF
       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF
       cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
       cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
       cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
       cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
       cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF        
       cospOUT%parasolGrid_refl(:,:)        = R_UNDEF
    endif
    if (any(cospgridIN%skt .lt. 0)) then
       call errorMessage('ERROR: COSP input variable: cospgridIN%skt contains values out of range')
       Lisccp_subcolumn = .false.
       Lisccp_column    = .false.     
       Lrttov_subcolumn = .false.
       cospOUT%rttov_tbs(:,:)         = R_UNDEF
       cospOUT%isccp_totalcldarea(:)  = R_UNDEF
       cospOUT%isccp_meantb(:)        = R_UNDEF
       cospOUT%isccp_meantbclr(:)     = R_UNDEF
       cospOUT%isccp_meanptop(:)      = R_UNDEF
       cospOUT%isccp_meantaucld(:)    = R_UNDEF
       cospOUT%isccp_meanalbedocld(:) = R_UNDEF
       cospOUT%isccp_boxtau(:,:)      = R_UNDEF    
       cospOUT%isccp_boxptop(:,:)     = R_UNDEF     
       cospOUT%isccp_fq(:,:,:)        = R_UNDEF        
    endif

	! RTTOV Inputs
    if (cospgridIN%zenang .lt. -90. .OR. cospgridIN%zenang .gt. 90) then
       call errorMessage('ERROR: COSP input variable: cospgridIN%zenang contains values out of range')
       Lrttov_subcolumn = .false.
       cospOUT%rttov_tbs(:,:) = R_UNDEF
    endif
    if (cospgridIN%co2 .lt. 0) then
       call errorMessage('ERROR: COSP input variable: cospgridIN%co2 contains values out of range')
       Lrttov_subcolumn = .false.
       cospOUT%rttov_tbs(:,:) = R_UNDEF
    endif
    if (cospgridIN%ch4 .lt. 0) then
       call errorMessage('ERROR: COSP input variable: cospgridIN%ch4 contains values out of range')
       Lrttov_subcolumn = .false.
       cospOUT%rttov_tbs(:,:) = R_UNDEF
    endif
    if (cospgridIN%n2o .lt. 0) then
       call errorMessage('ERROR: COSP input variable: cospgridIN%n2o contains values out of range')
       Lrttov_subcolumn = .false.
       cospOUT%rttov_tbs(:,:) = R_UNDEF
    endif
    if (cospgridIN%co.lt. 0) then
       call errorMessage('ERROR: COSP input variable: cospgridIN%co contains values out of range')
       Lrttov_subcolumn = .false.
       cospOUT%rttov_tbs(:,:) = R_UNDEF
    endif
    if (any(cospgridIN%o3 .lt. 0)) then
       call errorMessage('ERROR: COSP input variable: cospgridIN%o3 contains values out of range')
       Lrttov_subcolumn = .false.
       cospOUT%rttov_tbs(:,:) = R_UNDEF
    endif
    if (any(cospgridIN%emis_sfc .lt. 0. .OR. cospgridIN%emis_sfc .gt. 1)) then
       call errorMessage('ERROR: COSP input variable: cospgridIN%emis_sfc contains values out of range')
       Lrttov_subcolumn = .false.
       cospOUT%rttov_tbs(:,:) = R_UNDEF
    endif
    if (any(cospgridIN%u_sfc .lt. -100. .OR. cospgridIN%u_sfc .gt. 100.)) then
       call errorMessage('ERROR: COSP input variable: cospIN%u_sfc contains values out of range')
       Lrttov_subcolumn = .false.
       cospOUT%rttov_tbs(:,:) = R_UNDEF
    endif
    if (any(cospgridIN%v_sfc .lt. -100. .OR. cospgridIN%v_sfc .gt. 100.)) then
       call errorMessage('ERROR: COSP input variable: cospIN%v_sfc contains values out of range')
       Lrttov_subcolumn = .false.
       cospOUT%rttov_tbs(:,:) = R_UNDEF
    endif
    if (any(cospgridIN%t_sfc .lt. 150 .OR. cospgridIN%t_sfc .gt. 350.)) then
       call errorMessage('ERROR: COSP input variable: cospIN%t_sfc contains values out of range')
       Lrttov_subcolumn = .false.
       cospOUT%rttov_tbs(:,:) = R_UNDEF
    endif
    if (any(cospgridIN%lat .lt. -90 .OR. cospgridIN%lat .gt. 90)) then
       call errorMessage('ERROR: COSP input variable: cospIN%lat contains values out of range')
       Lrttov_subcolumn = .false.
       cospOUT%rttov_tbs(:,:) = R_UNDEF
    endif

    ! COSP_INPUTS
    if (cospIN%emsfc_lw .lt. 0. .OR. cospIN%emsfc_lw .gt. 1.) then
       call errorMessage('ERROR: COSP input variable: cospIN%emsfc_lw contains values out of range')
       Lisccp_subcolumn = .false.
       Lisccp_column    = .false.
       cospOUT%isccp_totalcldarea(:)  = R_UNDEF
       cospOUT%isccp_meantb(:)        = R_UNDEF
       cospOUT%isccp_meantbclr(:)     = R_UNDEF
       cospOUT%isccp_meanptop(:)      = R_UNDEF
       cospOUT%isccp_meantaucld(:)    = R_UNDEF
       cospOUT%isccp_meanalbedocld(:) = R_UNDEF
       cospOUT%isccp_boxtau(:,:)      = R_UNDEF    
       cospOUT%isccp_boxptop(:,:)     = R_UNDEF     
       cospOUT%isccp_fq(:,:,:)        = R_UNDEF       
    endif
    if (any(cospIN%tau_067 .lt. 0)) then
       call errorMessage('ERROR: COSP input variable: cospIN%tau_067 contains values out of range')
       Lisccp_subcolumn = .false.
       Lisccp_column    = .false.
       Lmisr_subcolumn  = .false.
       Lmisr_column     = .false.
       Lmodis_subcolumn = .false.
       Lmodis_column    = .false.
       cospOUT%isccp_totalcldarea(:)                                = R_UNDEF
       cospOUT%isccp_meantb(:)                                      = R_UNDEF
       cospOUT%isccp_meantbclr(:)                                   = R_UNDEF
       cospOUT%isccp_meanptop(:)                                    = R_UNDEF
       cospOUT%isccp_meantaucld(:)                                  = R_UNDEF
       cospOUT%isccp_meanalbedocld(:)                               = R_UNDEF
       cospOUT%isccp_boxtau(:,:)                                    = R_UNDEF    
       cospOUT%isccp_boxptop(:,:)                                   = R_UNDEF     
       cospOUT%isccp_fq(:,:,:)                                      = R_UNDEF       
       cospOUT%misr_fq(:,:,:)                                       = R_UNDEF 
       cospOUT%misr_dist_model_layertops(:,:)                       = R_UNDEF
       cospOUT%misr_meanztop(:)                                     = R_UNDEF 
       cospOUT%misr_cldarea(:)                                      = R_UNDEF
       cospOUT%modis_Cloud_Fraction_Total_Mean(:)                   = R_UNDEF
       cospOUT%modis_Cloud_Fraction_Water_Mean(:)                   = R_UNDEF
       cospOUT%modis_Cloud_Fraction_Ice_Mean(:)                     = R_UNDEF
       cospOUT%modis_Cloud_Fraction_High_Mean(:)                    = R_UNDEF
       cospOUT%modis_Cloud_Fraction_Mid_Mean(:)                     = R_UNDEF
       cospOUT%modis_Cloud_Fraction_Low_Mean(:)                     = R_UNDEF
       cospOUT%modis_Optical_Thickness_Total_Mean(:)                = R_UNDEF
       cospOUT%modis_Optical_Thickness_Water_Mean(:)                = R_UNDEF
       cospOUT%modis_Optical_Thickness_Ice_Mean(:)                  = R_UNDEF
       cospOUT%modis_Optical_Thickness_Total_LogMean(:)             = R_UNDEF
       cospOUT%modis_Optical_Thickness_Water_LogMean(:)             = R_UNDEF
       cospOUT%modis_Optical_Thickness_Ice_LogMean(:)               = R_UNDEF
       cospOUT%modis_Cloud_Particle_Size_Water_Mean(:)              = R_UNDEF
       cospOUT%modis_Cloud_Particle_Size_Ice_Mean(:)                = R_UNDEF
       cospOUT%modis_Cloud_Top_Pressure_Total_Mean(:)               = R_UNDEF
       cospOUT%modis_Liquid_Water_Path_Mean(:)                      = R_UNDEF
       cospOUT%modis_Ice_Water_Path_Mean(:)                         = R_UNDEF
       cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure(:,:,:) = R_UNDEF
       cospOUT%modis_Optical_Thickness_vs_ReffICE(:,:,:)            = R_UNDEF
       cospOUT%modis_Optical_Thickness_vs_ReffLIQ(:,:,:)            = R_UNDEF
    endif
    if (any(cospIN%emiss_11 .lt. 0. .OR. cospIN%emiss_11 .gt. 1)) then
       call errorMessage('ERROR: COSP input variable: cospIN%emiss_11 contains values out of range')
       Lisccp_subcolumn = .false.
       Lisccp_column    = .false.
       cospOUT%isccp_totalcldarea(:)  = R_UNDEF
       cospOUT%isccp_meantb(:)        = R_UNDEF
       cospOUT%isccp_meantbclr(:)     = R_UNDEF
       cospOUT%isccp_meanptop(:)      = R_UNDEF
       cospOUT%isccp_meantaucld(:)    = R_UNDEF
       cospOUT%isccp_meanalbedocld(:) = R_UNDEF
       cospOUT%isccp_boxtau(:,:)      = R_UNDEF    
       cospOUT%isccp_boxptop(:,:)     = R_UNDEF     
       cospOUT%isccp_fq(:,:,:)        = R_UNDEF          
    endif
    if (any(cospIN%asym .lt. -1. .OR. cospIN%asym .gt. 1)) then
       call errorMessage('ERROR: COSP input variable: cospIN%asym contains values out of range')
       Lmodis_subcolumn = .false.
       Lmodis_column    = .false.
       cospOUT%modis_Cloud_Fraction_Total_Mean(:)                   = R_UNDEF
       cospOUT%modis_Cloud_Fraction_Water_Mean(:)                   = R_UNDEF
       cospOUT%modis_Cloud_Fraction_Ice_Mean(:)                     = R_UNDEF
       cospOUT%modis_Cloud_Fraction_High_Mean(:)                    = R_UNDEF
       cospOUT%modis_Cloud_Fraction_Mid_Mean(:)                     = R_UNDEF
       cospOUT%modis_Cloud_Fraction_Low_Mean(:)                     = R_UNDEF
       cospOUT%modis_Optical_Thickness_Total_Mean(:)                = R_UNDEF
       cospOUT%modis_Optical_Thickness_Water_Mean(:)                = R_UNDEF
       cospOUT%modis_Optical_Thickness_Ice_Mean(:)                  = R_UNDEF
       cospOUT%modis_Optical_Thickness_Total_LogMean(:)             = R_UNDEF 
       cospOUT%modis_Optical_Thickness_Water_LogMean(:)             = R_UNDEF
       cospOUT%modis_Optical_Thickness_Ice_LogMean(:)               = R_UNDEF
       cospOUT%modis_Cloud_Particle_Size_Water_Mean(:)              = R_UNDEF
       cospOUT%modis_Cloud_Particle_Size_Ice_Mean(:)                = R_UNDEF
       cospOUT%modis_Cloud_Top_Pressure_Total_Mean(:)               = R_UNDEF
       cospOUT%modis_Liquid_Water_Path_Mean(:)                      = R_UNDEF
       cospOUT%modis_Ice_Water_Path_Mean(:)                         = R_UNDEF 
       cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure(:,:,:) = R_UNDEF
       cospOUT%modis_Optical_Thickness_vs_ReffICE(:,:,:)            = R_UNDEF   
       cospOUT%modis_Optical_Thickness_vs_ReffLIQ(:,:,:)            = R_UNDEF        
    endif
    if (any(cospIN%ss_alb .lt. 0 .OR. cospIN%ss_alb .gt. 1)) then
       call errorMessage('ERROR: COSP input variable: cospIN%ss_alb contains values out of range')
       Lmodis_subcolumn = .false.
       Lmodis_column    = .false.
       cospOUT%modis_Cloud_Fraction_Total_Mean(:)                   = R_UNDEF
       cospOUT%modis_Cloud_Fraction_Water_Mean(:)                   = R_UNDEF
       cospOUT%modis_Cloud_Fraction_Ice_Mean(:)                     = R_UNDEF
       cospOUT%modis_Cloud_Fraction_High_Mean(:)                    = R_UNDEF
       cospOUT%modis_Cloud_Fraction_Mid_Mean(:)                     = R_UNDEF
       cospOUT%modis_Cloud_Fraction_Low_Mean(:)                     = R_UNDEF
       cospOUT%modis_Optical_Thickness_Total_Mean(:)                = R_UNDEF
       cospOUT%modis_Optical_Thickness_Water_Mean(:)                = R_UNDEF
       cospOUT%modis_Optical_Thickness_Ice_Mean(:)                  = R_UNDEF
       cospOUT%modis_Optical_Thickness_Total_LogMean(:)             = R_UNDEF 
       cospOUT%modis_Optical_Thickness_Water_LogMean(:)             = R_UNDEF
       cospOUT%modis_Optical_Thickness_Ice_LogMean(:)               = R_UNDEF
       cospOUT%modis_Cloud_Particle_Size_Water_Mean(:)              = R_UNDEF
       cospOUT%modis_Cloud_Particle_Size_Ice_Mean(:)                = R_UNDEF
       cospOUT%modis_Cloud_Top_Pressure_Total_Mean(:)               = R_UNDEF
       cospOUT%modis_Liquid_Water_Path_Mean(:)                      = R_UNDEF
       cospOUT%modis_Ice_Water_Path_Mean(:)                         = R_UNDEF 
       cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure(:,:,:) = R_UNDEF
       cospOUT%modis_Optical_Thickness_vs_ReffICE(:,:,:)            = R_UNDEF   
       cospOUT%modis_Optical_Thickness_vs_ReffLIQ(:,:,:)            = R_UNDEF          
    endif
    if (any(cospIN%betatot .lt. 0)) then
       call errorMessage('ERROR: COSP input variable: cospIN%betatot contains values out of range')
       Lcalipso_subcolumn = .false.
       Lcalipso_column    = .false.
       cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
       cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
       cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF
       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF    
       cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
       cospOUT%calipso_cldlayer(:,:)        = R_UNDEF   
       cospOUT%calipso_srbval(:)            = R_UNDEF 
    endif
    if (any(cospIN%betatot_liq .lt. 0)) then
       call errorMessage('ERROR: COSP input variable: cospIN%betatot_liq contains values out of range')
       Lcalipso_subcolumn = .false.
       Lcalipso_column    = .false.
       cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
       cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
       cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF
       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF    
       cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
       cospOUT%calipso_cldlayer(:,:)        = R_UNDEF   
       cospOUT%calipso_srbval(:)            = R_UNDEF
    endif
    if (any(cospIN%betatot_ice .lt. 0)) then
       call errorMessage('ERROR: COSP input variable: cospIN%betatot_ice contains values out of range')
       Lcalipso_subcolumn = .false.
       Lcalipso_column    = .false.
       cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
       cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
       cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF
       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF    
       cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
       cospOUT%calipso_cldlayer(:,:)        = R_UNDEF   
       cospOUT%calipso_srbval(:)            = R_UNDEF
    endif 
    if (any(cospIN%beta_mol .lt. 0)) then
       call errorMessage('ERROR: COSP input variable: cospIN%beta_mol contains values out of range')
       Lcalipso_subcolumn = .false.
       Lcalipso_column    = .false.
       Lcloudsat_column   = .false.
       cospOUT%calipso_lidarcldphase(:,:,:)       = R_UNDEF
       cospOUT%calipso_cldlayerphase(:,:,:)       = R_UNDEF
       cospOUT%calipso_lidarcldtmp(:,:,:)         = R_UNDEF
       cospOUT%calipso_cfad_sr(:,:,:)             = R_UNDEF    
       cospOUT%calipso_lidarcld(:,:)              = R_UNDEF
       cospOUT%calipso_cldlayer(:,:)              = R_UNDEF   
       cospOUT%calipso_srbval(:)                  = R_UNDEF   
       cospOUT%cloudsat_cfad_ze(:,:,:)            = R_UNDEF
       cospOUT%lidar_only_freq_cloud(:,:)         = R_UNDEF
       cospOUT%radar_lidar_tcc(:)                 = R_UNDEF            
    endif    
    if (any(cospIN%tautot .lt. 0)) then
       call errorMessage('ERROR: COSP input variable: cospIN%tautot contains values out of range')
       Lcalipso_subcolumn = .false.
       Lcalipso_column    = .false.
       cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
       cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
       cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF
       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF    
       cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
       cospOUT%calipso_cldlayer(:,:)        = R_UNDEF   
       cospOUT%calipso_srbval(:)            = R_UNDEF        
    endif
    if (any(cospIN%tautot_liq .lt. 0)) then
       call errorMessage('ERROR: COSP input variable: cospIN%tautot_liq contains values out of range')
       Lcalipso_subcolumn = .false.
       Lcalipso_column    = .false.
       cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
       cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
       cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF
       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF    
       cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
       cospOUT%calipso_cldlayer(:,:)        = R_UNDEF   
       cospOUT%calipso_srbval(:)            = R_UNDEF          
    endif
    if (any(cospIN%tautot_ice .lt. 0)) then
       call errorMessage('ERROR: COSP input variable: cospIN%tautot_ice contains values out of range')
       Lcalipso_subcolumn = .false.
       Lcalipso_column    = .false.
       cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
       cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
       cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF
       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF    
       cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
       cospOUT%calipso_cldlayer(:,:)        = R_UNDEF   
       cospOUT%calipso_srbval(:)            = R_UNDEF          
    endif
    if (any(cospIN%tau_mol .lt. 0)) then
       call errorMessage('ERROR: COSP input variable: cospIN%tau_mol contains values out of range')
       Lcalipso_subcolumn = .false.
       Lcalipso_column    = .false.
       cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
       cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
       cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF
       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF    
       cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
       cospOUT%calipso_cldlayer(:,:)        = R_UNDEF   
       cospOUT%calipso_srbval(:)            = R_UNDEF           
    endif    
    if (any(cospIN%tautot_S_liq .lt. 0)) then
       call errorMessage('ERROR: COSP input variable: cospIN%tautot_S_liq contains values out of range')
       Lparasol_subcolumn = .false.
       Lparasol_column    = .false.
       cospOUT%parasolPix_refl(:,:,:) = R_UNDEF      
       cospOUT%parasolGrid_refl(:,:)  = R_UNDEF
    endif
    if (any(cospIN%tautot_S_ice .lt. 0)) then
       call errorMessage('ERROR: COSP input variable: cospIN%tautot_S_ice contains values out of range')
       Lparasol_subcolumn = .false.
       Lparasol_column    = .false.
       cospOUT%parasolPix_refl(:,:,:) = R_UNDEF      
       cospOUT%parasolGrid_refl(:,:)  = R_UNDEF       
    endif    
    if (any(cospIN%z_vol .lt. 0)) then
       call errorMessage('ERROR: COSP input variable: cospIN%z_vol contains values out of range')
       Lcloudsat_subcolumn = .false.
       Lcloudsat_column    = .false.
       cospOUT%cloudsat_Ze_tot(:,:,:)     = R_UNDEF
       cospOUT%cloudsat_cfad_ze(:,:,:)    = R_UNDEF
       cospOUT%lidar_only_freq_cloud(:,:) = R_UNDEF
       cospOUT%radar_lidar_tcc(:)         = R_UNDEF         
    endif
    if (any(cospIN%kr_vol .lt. 0)) then
       call errorMessage('ERROR: COSP input variable: cospIN%kr_vol contains values out of range')
       Lcloudsat_subcolumn = .false.
       Lcloudsat_column    = .false.
       cospOUT%cloudsat_Ze_tot(:,:,:)     = R_UNDEF
       cospOUT%cloudsat_cfad_ze(:,:,:)    = R_UNDEF
       cospOUT%lidar_only_freq_cloud(:,:) = R_UNDEF
       cospOUT%radar_lidar_tcc(:)         = R_UNDEF       
    endif    
    if (any(cospIN%g_vol .lt. 0)) then
       call errorMessage('ERROR: COSP input variable: cospIN%g_vol contains values out of range')
       Lcloudsat_subcolumn = .false.
       Lcloudsat_column    = .false.
       cospOUT%cloudsat_Ze_tot(:,:,:)     = R_UNDEF
       cospOUT%cloudsat_cfad_ze(:,:,:)    = R_UNDEF
       cospOUT%lidar_only_freq_cloud(:,:) = R_UNDEF
       cospOUT%radar_lidar_tcc(:)         = R_UNDEF        
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
      call errorMessage('ERROR(isccp_simulator): The number of points in the input fields are inconsistent')
  endif
  if (size(cospIN%frac_out,2) .ne. cospIN%Ncolumns .OR. &
      size(cospIN%tau_067,2)  .ne. cospIN%Ncolumns .OR. &
      size(cospIN%emiss_11,2) .ne. cospIN%Ncolumns) then
      Lisccp_subcolumn = .false.
      Lisccp_column    = .false.
      call errorMessage('ERROR(isccp_simulator): The number of sub-columns in the input fields are inconsistent')
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
      call errorMessage('ERROR(isccp_simulator): The number of levels in the input fields are inconsistent')
  endif
      
  ! MISR
  if (size(cospIN%tau_067,1)        .ne. cospIN%Npoints .OR. &
      size(cospgridIN%sunlit)       .ne. cospIN%Npoints .OR. & 
      size(cospgridIN%hgt_matrix,1) .ne. cospIN%Npoints .OR. &
      size(cospgridIN%at,1)         .ne. cospIN%Npoints) then
      Lmisr_subcolumn = .false.
      Lmisr_column    = .false.
      call errorMessage('ERROR(misr_simulator): The number of points in the input fields are inconsistent')
  endif
  if (size(cospIN%tau_067,2) .ne. cospIN%Ncolumns) then
      Lmisr_subcolumn = .false.
      Lmisr_column    = .false.
      call errorMessage('ERROR(misr_simulator): The number of sub-columns in the input fields are inconsistent')
  endif
  if (size(cospIN%tau_067,3)        .ne. cospIN%Nlevels .OR. &
      size(cospgridIN%hgt_matrix,2) .ne. cospIN%Nlevels .OR. &
      size(cospgridIN%at,2)         .ne. cospIN%Nlevels) then
      Lmisr_subcolumn = .false.
      Lmisr_column    = .false.
      call errorMessage('ERROR(misr_simulator): The number of levels in the input fields are inconsistent')
  endif    

  ! MODIS
  if (size(cospIN%fracLiq,1) .ne. cospIN%Npoints .OR. &
      size(cospIN%tau_067,1) .ne. cospIN%Npoints .OR. &
      size(cospIN%asym,1)    .ne. cospIN%Npoints .OR. &
      size(cospIN%ss_alb,1)  .ne. cospIN%Npoints) then
      Lmodis_subcolumn = .false.
      Lmodis_column    = .false.
      call errorMessage('ERROR(modis_simulator): The number of points in the input fields are inconsistent')
  endif
  if (size(cospIN%fracLiq,2) .ne. cospIN%Ncolumns .OR. &
      size(cospIN%tau_067,2) .ne. cospIN%Ncolumns .OR. &
      size(cospIN%asym,2)    .ne. cospIN%Ncolumns .OR. &
      size(cospIN%ss_alb,2)  .ne. cospIN%Ncolumns) then
      Lmodis_subcolumn = .false.
      Lmodis_column    = .false.
      call errorMessage('ERROR(modis_simulator): The number of sub-columns in the input fields are inconsistent')
  endif        
  if (size(cospIN%fracLiq,3) .ne. cospIN%Nlevels .OR. &
      size(cospIN%tau_067,3) .ne. cospIN%Nlevels .OR. &
      size(cospIN%asym,3)    .ne. cospIN%Nlevels .OR. &
      size(cospIN%ss_alb,3)  .ne. cospIN%Nlevels) then
      Lmodis_subcolumn = .false.
      Lmodis_column    = .false.
      call errorMessage('ERROR(modis_simulator): The number of levels in the input fields are inconsistent')
  endif  
  
  ! CLOUDSAT    
  if (size(cospIN%z_vol,1)          .ne. cospIN%Npoints .OR. &
      size(cospIN%kr_vol,1)         .ne. cospIN%Npoints .OR. &
      size(cospIN%g_vol,1)          .ne. cospIN%Npoints .OR. &
      size(cospgridIN%hgt_matrix,1) .ne. cospIN%Npoints) then
      Lcloudsat_subcolumn = .false.
      Lcloudsat_column    = .false.
      call errorMessage('ERROR(cloudsat_simulator): The number of points in the input fields are inconsistent')
  endif
  if (size(cospIN%z_vol,2)  .ne. cospIN%Ncolumns .OR. &
      size(cospIN%kr_vol,2) .ne. cospIN%Ncolumns .OR. &
      size(cospIN%g_vol,2)  .ne. cospIN%Ncolumns) then
      Lcloudsat_subcolumn = .false.
      Lcloudsat_column    = .false.
      call errorMessage('ERROR(cloudsat_simulator): The number of sub-columns in the input fields are inconsistent')
  endif       
  if (size(cospIN%z_vol,3)          .ne. cospIN%Nlevels .OR. &
      size(cospIN%kr_vol,3)         .ne. cospIN%Nlevels .OR. &
      size(cospIN%g_vol,3)          .ne. cospIN%Nlevels .OR. &
      size(cospgridIN%hgt_matrix,2) .ne. cospIN%Nlevels) then
      Lcloudsat_subcolumn = .false.
      Lcloudsat_column    = .false.
      call errorMessage('ERROR(cloudsat_simulator): The number of levels in the input fields are inconsistent')
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
      call errorMessage('ERROR(calipso_simulator): The number of points in the input fields are inconsistent')
  endif          
   if (size(cospIN%betatot,2)     .ne. cospIN%Ncolumns .OR. &
       size(cospIN%betatot_liq,2) .ne. cospIN%Ncolumns .OR. &
       size(cospIN%betatot_ice,2) .ne. cospIN%Ncolumns .OR. &
       size(cospIN%tautot,2)      .ne. cospIN%Ncolumns .OR. &
       size(cospIN%tautot_liq,2)  .ne. cospIN%Ncolumns .OR. &
       size(cospIN%tautot_ice,2)  .ne. cospIN%Ncolumns) then
       Lcalipso_subcolumn = .false.
       Lcalipso_column    = .false.
      call errorMessage('ERROR(calipso_simulator): The number of sub-columns in the input fields are inconsistent')
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
      call errorMessage('ERROR(calipso_simulator): The number of levels in the input fields are inconsistent')
  endif 
  
  ! PARASOL
  if (size(cospIN%tautot_S_liq,1) .ne. cospIN%Npoints .OR. &
      size(cospIN%tautot_S_ice,1) .ne. cospIN%Npoints) then
      Lparasol_subcolumn = .false.
      Lparasol_column    = .false.
      call errorMessage('ERROR(parasol_simulator): The number of points in the input fields are inconsistent')
  endif
  if (size(cospIN%tautot_S_liq,2) .ne. cospIN%Nlevels .OR. &
      size(cospIN%tautot_S_ice,2) .ne. cospIN%Nlevels) then
      Lparasol_subcolumn = .false.
      Lparasol_column    = .false.
      call errorMessage('ERROR(parasol_simulator): The number of levels in the input fields are inconsistent')
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
      size(cospgridIN%t_sfc)             .ne. cospIN%Npoints .OR. &
      size(cospgridIN%qv,1)              .ne. cospIN%Npoints .OR. &
      size(cospgridIN%land)              .ne. cospIN%Npoints .OR. &
      size(cospgridIN%lat)               .ne. cospIN%Npoints) then
      Lrttov_subcolumn = .false.
      Lrttov_column    = .false.
      call errorMessage('ERROR(rttov_simulator): The number of points in the input fields are inconsistent')
  endif      
  if (size(cospgridIN%pfull,2)           .ne. cospIN%Nlevels   .OR. &
      size(cospgridIN%at,2)              .ne. cospIN%Nlevels   .OR. &
      size(cospgridIN%qv,2)              .ne. cospIN%Nlevels   .OR. &
      size(cospgridIN%hgt_matrix_half,2) .ne. cospIN%Nlevels   .OR. &
      size(cospgridIN%phalf,2)           .ne. cospIN%Nlevels+1 .OR. &
      size(cospgridIN%qv,2)              .ne. cospIN%Nlevels) then
      Lrttov_subcolumn = .false.
      Lrttov_column    = .false.
      call errorMessage('ERROR(rttov_simulator): The number of levels in the input fields are inconsistent')
  endif       
             
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Part 3: Instrument simulator specific error checking. This section contains error
  !         checking that was originally contained within the simulators.
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! PARASOL Simulator
  if (Lparasol_subcolumn) then
      IF (PARASOL_NREFL .gt. ntetas ) THEN
      write(parasolErrorMessage,"(a50,i2,a4,i2)") 'ERROR(lidar_simulator): nrefl should be less then ',ntetas,' not',PARASOL_NREFL
         call errorMessage(parasolErrorMessage)
      ENDIF
      Lparasol_subcolumn = .false.
      Lparasol_column    = .false.
  endif  
    
  end subroutine cosp_errorCheck
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! END MODULE
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
END MODULE MOD_COSP
