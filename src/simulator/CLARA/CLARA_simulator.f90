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
! November 2015- D. Swales - Original version
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module mod_clara_sim
  USE COSP_KINDS,            ONLY: wp
  USE COSP_PHYS_CONSTANTS,   ONLY: amw,amd,grav,avo,rd
  USE clara_rttov_interface, ONLY: get_Tb_rttov_clear,rttov_simple,get_Tb_rttov_clouds
  USE mod_modis_sim,         ONLY: compute_toa_reflectace,two_stream_reflectance,        &
                                   interpolate_to_min
  USE mod_cosp_config,       ONLY: numCLARAtauBins,numCLARApresBins,clara_histTau,       &
                                   clara_histPres
  use MOD_COSP_STATS,        ONLY: hist2D

  implicit none

  integer,parameter :: &
     CLARA_phaseIsNone         = 0, & ! No retrieved phase 
     CLARA_phaseIsLiquid       = 1, & ! Liquid phase
     CLARA_phaseIsIce          = 2, & ! Ice phase
     CLARA_phaseIsUndetermined = 3, & ! Undetermined
     num_trial_res             = 15   ! Number of reference radii to use for size retrieval.
  real(wp),parameter :: &
     clara_t0 =  296._wp,      & ! *Not the best description*  
     pstd     = 1013250._wp,   & ! Mean sea-level pressure (Pa)
     tb4_rad  = 1331.922_wp      ! Conversion factor for radiance to brightness temperature
            
  ! Retrieval parameters (set up during initialization)
  real(wp) :: &
     CLARA_upperCloudTauLim,   & ! Optical depth into cloud AVHRR can see for phase retrieval
     CLARA_minOpticalThickness,& ! Lower limit of optical sensitivity for AVHRR
     CLARA_phaseThresh,        & ! Fraction of total extincton needs to be in a single 
                                 ! category to make phase determination work     
     CLARA_STlimit,            & ! Optical depth limit for opaque clouds
     CLARA_re_water_min,       & ! Minimum effective radius (liquid)
     CLARA_re_water_max,       & ! Maximum effective radius (liquid)
     CLARA_re_ice_min,         & ! Minimum effective radius (ice)
     CLARA_re_ice_max            ! Minimum effective radius (ice)     
  logical :: & 
     CLARA_RTTOVclr,           & ! Flag to use RTTOV for clear-sky brightness temperature
     CLARA_retSize               ! If true, use TOA reflectance minimization (same as MODIS) 
                                 ! for particle size retrieval. If false, use cloud top 
                                 ! weighted-effective radius retrieval method.
  integer :: &
     nChannels,                & ! Number of channels used by RTTOV
     CLARA_Tb_subvis,          & ! Options for how to handle sub-visible clouds.
                                 ! 1) Use RTTOV with scattering
                                 ! 2) Use clear-sky RTTOV retrieval
     CLARA_Tb_semitrans,       & ! Options for how to handle semi-transparent clouds.
                                 ! 1) Use RTTOV with scattering (Only realistic option)
     CLARA_Tb_opaque             ! Options for how to handle opaque clouds.
                                 ! 1) Use RTTOV, with scattering.
                                 ! 2) Use RTTOV, assuming black-body w/o scattering
  real(wp),dimension(:),pointer :: &
     RTTOV_satLat,             & ! Latitude array for satellite zenith angles
     RTTOV_satzen                ! Satellite zenith angles and a function of latitude for
                                 ! satellite being used.
  real(wp),dimension(num_trial_res) :: &
     trial_re_w,               & ! Near-IR optical params vs size for retrieval scheme (liquid)
     trial_re_i                  ! Near-IR optical params vs size for retrieval scheme (ice)
  real(wp),dimension(num_trial_res) :: &
     g_w,                      & ! Assymettry parameter for size retrieval (liquid)
     g_i,                      & ! Assymettry parameter for size retrieval (ice)
     w0_w,                     & ! Single-scattering albedo for size retrieval (liquid)
     w0_i                        ! Single-scattering albedo for size retrieval (ice)

  contains
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE clara_subcolumn
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine clara_subcolumn(nPoints,nLevels,nSubCols,pfull,phalf,qv,at,skt,t2m,p2m,     &
                             emsfc_lw,lsmask,lat,cldLiqRTTOV,cldIceRTTOV,tau,            &
                             tautotliq,tautotice,reffLiq,reffIce,g,w0,tauLiqFrac,        &
                             clara_tau,clara_ctp,clara_ctt,clara_cth,clara_size,clara_phase)
    ! Inputs
    integer,intent(in) :: &
      nPoints,      & ! Number of gridpoints
      nLevels,      & ! Number of vertical levels
      nSubCols        ! Number of subcolumns
    real(wp),intent(in) :: &
      emsfc_lw        ! 10.5 micron emissivity of surface (fraction) 
    real(wp),intent(in), dimension(nPoints) :: &
      skt,          & ! Skin temperature
      t2m,          & ! Temperature at 2 meters.
      p2m,          & ! Pressure at 2 meters.
      lsmask,       & ! Land/Sea mask
      lat             ! Latitude
    real(wp),intent(in), dimension(nPoints,nLevels) ::  &
      at,           & ! Temperature 
      pfull,        & ! Pressure at model full levels
      qv              ! Specific humidity
    real(wp),intent(in), dimension(nPoints,nLevels+1) :: &
      phalf           ! Pressure at model half-levels
    real(wp),dimension(nPoints,nSubCols,nLevels),intent(in) :: &
      tau,          & ! TOA-2-SFC level subcolumn optical thickness @ 0.67 microns (total).
      tautotliq,    & ! TOA-2-SFC level subcolumn optical thickness @ 0.67 microns (liquid).
      tautotice,    & ! TOA-2-SFC level subcolumn optical thickness @ 0.67 microns (ice).
      reffLiq,      & ! TOA-2-SFC model cloud particle size profile (liquid)
      reffIce,      & ! TOA-2-SFC model cloud particle size profile (ice)
      g,            & ! Subcolumn assymetry parameter @ 3.7 microns.
      w0,           & ! Subcolumn single-scattering albedo @ 3.7 microns.
      tauLiqFrac      ! Fractional contribution of optical depth due to liquid water
    real,dimension(nPoints,nSubCols,6,nLevels),intent(in) :: &  
      cldLiqRTTOV,  & ! Subcolumn liquid cloud concentration (needed by RTTOV)
      cldIceRTTOV     ! Subcolumn ice cloud concentration (needed by RTTOV)
      
    ! Outputs  
    real(wp),dimension(nPoints,nSubcols),intent(out) :: &
       clara_tau,     & ! Retrieved CLARA optical-depth
       clara_ctp,     & ! Retrieved CLARA cloud-top pressure
       clara_ctt,     & ! Retrieved CLARA cloud-top temperature
       clara_cth,     & ! Retrieved CLARA cloud-top height
       clara_size       ! Retrieved CLARA particle size
    integer(wp),dimension(nPoints,nSubCols),intent(out) :: &
       clara_phase      ! Retrieved CLARA cloud phase
       
    ! Local variables
    integer                                          :: ij,ik,il,nchan
    integer,dimension(2)                             :: xi
    integer,dimension(nPoints)                       :: rt_sfc,rt_h2otype,ninv
    integer,dimension(nPoints,nLevels)               :: inversion 
    integer,dimension(nPoints,nSubCols)              :: clara_cldFlag
    real(wp)                                         :: intLiqFrac,obs_Refl_nir
    real(wp),dimension(nPoints)                      :: rt_satzen,irradSfc,clara_tbClr
    real(wp),dimension(nPoints,nSubCols,6,nLevels-1) :: rt_cloud,cldLiqRTTOVl,cldIceRTTOVl
    real(wp),dimension(nPoints,nSubCols,nLevels-1)   :: rt_cfrac
    real(wp),dimension(nPoints,nSubCols,nLevels  )   :: cfrac,tautot 
    real(wp),dimension(nPoints,nLevels+1)            :: Thalf 
    real(wp),dimension(nPoints,nChannels)            :: rt_emis, TbTempClr
    real(wp),dimension(nPoints,nSubCols)             :: clara_tbAll,test1
    real(wp),dimension(nChannels)                    :: TbTempALL    
    real(wp),dimension(nPoints,nLevels)              :: irrad,emis_wv,atC
    real(wp),dimension(num_trial_res)                :: predicted_Refl_nir
    logical,dimension(nPoints,nSubCols)              :: clara_cldMask
        
    ! ####################################################################################    
    ! Need to compute some quantities needed by RTTOV.
    ! ####################################################################################
    
    ! Cumulative optical depth profile
    do ij=1,nLevels
       tautot(1:nPoints,1:nSubCols,ij) = sum(tau(1:nPoints,1:nSubCols,1:ij),3)
    enddo   
       
	! Determine surface type (RTTOV has a different convention than most GCMs, so we 
	! need to flip it here. RTTOV requires to convention land=0,ocean=1)
    rt_sfc(1:nPoints) = 1
    where(lsmask(1:nPoints) .eq. 1) rt_sfc(1:nPoints) = 0

	! Set water type flag. Currently just set to 1.
    rt_h2otype(1:nPoints) = 1

    ! Compute surface emissivity for RTTOV channels. *NOT YET IMPLEMENTED, so use 10.5 
    ! micron emissivity*
    rt_emis(1:nPoints,1:nChannels) = emsfc_lw
       
	! Satellite viewing angle for RTTOV. Need to interpolate points to correct satellite
	! viewing angle. This is satellite dependent and set up during initialization.
    do il=1,nPoints
       ! Interpolation bounds
       xi(1) = minloc(abs((RTTOV_satLat+90) - (lat(il)+90)), 1)
       if (lat(il) .lt. RTTOV_satLat(xi(1))) xi(1)=xi(1)-1
       xi(2) = xi(1)+1
       ! Interpolation
       if (lat(il) .gt. minval(RTTOV_satLat) .and. lat(il) .lt. maxval(RTTOV_satLat)) then
          rt_satzen(il) = RTTOV_satZen(xi(1))+(lat(il)-RTTOV_satLat(xi(1)))*&
                          (RTTOV_satZen(xi(2))-RTTOV_satZen(xi(1)))/&
                          (RTTOV_satLat(xi(2))-RTTOV_satLat(xi(1)))
       endif
       ! Out of bounds, so use southernmost/northernmost satellite zenith angle.
       if (lat(il) .lt. minval(RTTOV_satLat)) rt_satzen(il) = RTTOV_satZen(1)
       if (lat(il) .gt. maxval(RTTOV_satLat)) rt_satzen(il) = RTTOV_satZen(size(RTTOV_satZen,1))
    enddo

    ! Compute layer cloud condensate amounts for both liquid and ice
    cldLiqRTTOVl = 0.5*(cldLiqRTTOV(1:nPoints,1:nSubCols,:,2:nLevels) + &
                        cldLiqRTTOV(1:nPoints,1:nSubCols,:,1:nLevels-1))
    cldIceRTTOVl = 0.5*(cldIceRTTOV(1:nPoints,1:nSubCols,:,2:nLevels) + &
                        cldIceRTTOV(1:nPoints,1:nSubCols,:,1:nLevels-1)) 

    ! Use liquid condensate for RTTOV cloud-types 1-5
    rt_cloud(1:nPoints,1:nSubCols,:,1:nLevels-1) = cldLiqRTTOVl
    ! Use ice condensate fot RTTOV cloud-type 6.
    where(sum(rt_cloud,3) .eq. 0) 
       rt_cloud(1:nPoints,1:nSubCols,6,1:nLevels-1) = sum(cldIceRTTOVl,3)
    end where
    
    ! Compute LAYER cloud fraction amounts. If any of the 6 RTTOV cloud-types contain
    ! condensate, then cloudy subcolumn.
    rt_cfrac(1:nPoints,1:nSubCols,1:nLevels-1) = 0._wp
    where(sum(cldIceRTTOVl+cldLiqRTTOVl,3) .ne. 0) rt_cfrac = 0.999_wp

    ! Compute LEVEL cloud fraction amounts. Anywhere (in the 6 cloud type for both liquid 
    ! and ice) where a cloud is identified, assign cloud-fraction to subcolumns
    cfrac(1:nPoints,1:nSubCols,1:nLevels) = 0._wp
    where(sum(cldIceRTTOV+cldLiqRTTOV,3) .ne. 0) cfrac = 1._wp
    
    ! RTTOV requires temperature at the interfaces, so here we compute it.
    thalf = expand_profile(nPoints,nLevels,at,pfull,phalf)    
    
    ! We need to identify temperature inversions within each profile.
    do ij=1,nPoints  
       inversion(ij,1:nLevels) = find_temperature_inversion(nLevels,thalf(ij,1:nLevels+1),&
                                                            phalf(ij,1:nLevels+1))
    enddo 
    ! ####################################################################################
    ! ####################################################################################
    !                                BEGIN RETRIEVAL
    ! ####################################################################################
    ! ####################################################################################
    
    ! ####################################################################################
    ! Compute total column optical depth. 
    ! ####################################################################################
    clara_tau(1:nPoints,1:nSubCols) = sum(tau(1:nPoints,1:nSubCols,1:nLevels),3)

    ! ####################################################################################
    ! Compute cloud-mask.
    ! Cloudy scenes are only detectable above a certain optical depth, so use as
    ! ####################################################################################
    clara_cldMask(1:nPoints,1:nSubCols) = .false.
    where(clara_tau(1:nPoints, 1:nSubCols) .gt. CLARA_minOpticalThickness)               &
       clara_cldMask(1:nPoints,1:nSubCols) = .true.
    
    ! ####################################################################################
    ! Compute cloud thickness flags (0=clear, 1=sub-visible, 2=semi-transparent,3=opaque)
    ! ####################################################################################          
    where(clara_tau(1:nPoints, 1:nSubCols) .lt. 0.001)                                   &
       clara_cldFlag(1:nPoints, 1:nSubCols) = 0
    where(clara_tau(1:nPoints, 1:nSubCols) .lt. CLARA_minOpticalThickness)               &
        clara_cldFlag(1:nPoints, 1:nSubCols) = 1
    where(clara_tau(1:nPoints, 1:nSubCols) .ge. CLARA_minOpticalThickness .and. &
          clara_tau(1:nPoints, 1:nSubCols) .le. CLARA_STlimit)                           &
       clara_cldFlag(1:nPoints, 1:nSubCols) = 2
    where(clara_tau(1:nPoints, 1:nSubCols) .gt. CLARA_STlimit)                           &
       clara_cldFlag(1:nPoints, 1:nSubCols) = 3
    
    ! ####################################################################################
    ! Compute clear-sky brightness temperature
    ! a) Compute water vapor continuum emissivity this treatment follows Schwarkzopf 
    !    and Ramasamy JGR 1999,vol 104, pages 9467-9499. The emissivity is 
    !    calculated at a wavenumber of 955 cm-1, or 10.47 microns.           
    ! b) Use RTTOV
    ! ####################################################################################

    call clearSkyTb_CDK(nPoints,nLevels,phalf,pfull,qv,at,skt,emsfc_lw,irrad,irradSfc,   &
                        emis_wv,atC,clara_tbClr)
    if (CLARA_RTTOVclr) then

      ! Cloud free radiances from RTTOV	   
      call get_Tb_rttov_clear(nPoints,nLevels,nChannels,                                 &
                              pfull(1:nPoints,1:nLevels)/100._wp,                        &
                              at(1:nPoints,1:nLevels),qv(1:nPoints,1:nLevels),           &
                              t2m(1:nPoints),p2m(1:nPoints)/100._wp,skt(1:nPoints),      &
                              rt_h2otype(1:nPoints),rt_satzen(1:nPoints),                &
                              reshape(rt_emis,[nPoints*nChannels]),rt_sfc(1:nPoints),    &
                              tbTempClr(1:nPoints,1:nChannels))
      ! Use channel 2 for IR (10.8 microns) brightness temperature. *FIXME* The IR channel
      ! index should not be hardcoded in here. What if the number of RTTOV channels 
      ! changes? 
      clara_tbClr(1:nPoints) = tbTempClr(1:nPoints,2)      
    endif    
    
    ! ####################################################################################
    ! Loop over all points and subcolumns.
    ! ####################################################################################
    clara_size(1:nPoints,1:nSubcols) = 0._wp
    do ij=1,nPoints
       do ik=1,nSubCols 
          if (CLARA_cldMask(ij,ik)) then 
          ! ##############################################################################
          ! Compute cloud phase
          ! ##############################################################################          
          if (clara_cldMask(ij,ik)) then
             intLiqFrac = weight_by_extinction(nLevels,tau(ij,ik,1:nLevels),          &
                                       tauLiqFrac(ij,ik,1:nLevels),CLARA_upperCloudTauLim)
             if(intLiqFrac .ge. CLARA_phaseThresh) then 
                clara_phase(ij,ik) = CLARA_phaseIsLiquid
             else if (intLiqFrac .le. 1._wp - CLARA_phaseThresh) then 
                clara_phase(ij,ik) = CLARA_phaseIsIce
             else 
                clara_phase(ij,ik) = CLARA_phaseIsUndetermined
             end if       
          else
             clara_phase(ij,ik) = CLARA_phaseIsNone   
          endif       
          
          ! ##############################################################################
          ! Compute cloud particle size
          ! At the moment there are two ways to retrieve phase and is set using the     
          ! namelist parameter, CLARA_retSize.
          ! ##############################################################################
          obs_Refl_nir = compute_toa_reflectace(nLevels,tau(ij,ik,1:nLevels),         &
                                                g(ij,ik,1:nLevels), w0(ij,ik,1:nLevels))
          ! Compute predicted reflectance
          if(any(clara_phase(ij,ik) == (/ CLARA_phaseIsLiquid, CLARA_phaseIsUndetermined,&
                                          CLARA_phaseIsIce /))) then 
             if (CLARA_retSize) then
                ! Retrieve phase using TOA reflectance matching
                if (clara_phase(ij,ik) == CLARA_phaseIsLiquid .OR. clara_phase(ij,ik) ==    &
                    CLARA_phaseIsUndetermined) then
                   predicted_Refl_nir(1:num_trial_res) = two_stream_reflectance(clara_tau(ij,ik), &
                        g_w(1:num_trial_res), w0_w(1:num_trial_res))
                   clara_size(ij,ik) = 1.0e-06_wp*interpolate_to_min(trial_re_w(1:num_trial_res), &
                        predicted_Refl_nir(1:num_trial_res), obs_Refl_nir)
                else
                   predicted_Refl_nir(1:num_trial_res) = two_stream_reflectance(clara_tau(ij,ik), &
                        g_i(1:num_trial_res), w0_i(1:num_trial_res))
                   clara_size(ij,ik) = 1.0e-06_wp*interpolate_to_min(trial_re_i(1:num_trial_res), &
                        predicted_Refl_nir(1:num_trial_res), obs_Refl_nir)
                endif
                if (clara_size(ij,ik) .lt. 0) clara_size(ij,ik)=-999._wp
             else        
                ! Retrieve particle size using cloud-top weighted effective radius method
                call ctw_reff(nLevels,tau(ij,ik,1:nLevels),tautotliq(ij,ik,1:nLevels),   &
                              tautotice(ij,ik,1:nLevels),cfrac(ij,ik,1:nLevels),         &
                              clara_phase(ij,ik),reffLiq(ij,ik,1:nLevels),               &
                              reffIce(ij,ik,1:nLevels),clara_size(ij,ik))              
             endif
             
          else 
             clara_size(ij,ik) = -999._wp
          endif        

          ! ##############################################################################
          ! Compute all-sky brightness temperature and cloud-top temperature
          ! There are 3 different types of clouds we need to account for and several 
          ! methods to compute the brightness temperature depending on this type:
          !   1) Sub-visible clouds
          !      a) Treat as black-cloud, using RTTOV_simple
          !      b) Full RTTOV calculation
          !      c) Treat scene as clear
          !   2) Semi-transparent clouds
          !      a) Treat as black-cloud, using RTTOV_simple
          !      b) Full RTTOV calculation 
          !      c) Treat as black-cloud, using simple-scheme (CDK)
          !   3) Opaque clouds
          !      a) Treat as black-cloud, using RTTOV_simple
          !      b) Full RTTOV calculation
          !      c) Treat as black-cloud, using simple-scheme (CDK)
          ! ##############################################################################
 
           ! #############################################################################
           ! 1) Sub-visible clouds
           ! #############################################################################
           if (clara_cldFlag(ij,ik) .eq. 1) then
              ! a) Brightness temperature (using RTTOV w/o scattering)
              if (CLARA_Tb_subvis .eq. 1) then
                 call rttov_simple(nLevels,nChannels,tautot(ij,ik,:),                    &
                                   phalf(ij,1:nLevels+1)/100._wp,thalf(ij,1:nLevels+1),  &
                                   pfull(ij,1:nLevels)/100._wp,at(ij,1:nLevels),         &
                                   qv(ij,1:nLevels),t2m(ij),p2m(ij)/100._wp,skt(ij),     &
                                   rt_emis(ij,1:nChannels),rt_satzen(ij),                &
                                   CLARA_minOpticalThickness,rt_h2otype(ij),rt_sfc(ij),  &
                                   CLARA_upperCloudTauLim,                               &
                                   tbTempALL(1:nChannels))                        
              endif
              ! b) Brightness temperature (using RTTOV w/ scattering)
              if (CLARA_Tb_subvis .eq. 2) then  
                 call get_Tb_rttov_clouds(1,nLevels,nChannels,pfull(ij,1:nLevels)/100.,  &
                                          at(ij,1:nLevels),qv(ij,1:nLevels),t2m(ij),     &
                                          p2m(ij)/100._wp,skt(ij),                       &
                                          rt_emis(ij,1:nChannels),rt_satzen(ij),         &
                                          rt_h2otype(ij),rt_sfc(ij),                     &
                                          rt_cloud(ij,ik,:,1:nLevels-1),                 &
                                          rt_cfrac(ij,ik,1:nLevels-1),                   &
                                          tbTempALL(1:nChannels))    
                 clara_TbALL(ij,ik) = tbTempALL(2)                            
              endif 
              ! c) Brightness temperature (use RTTOV clear-sky)
              if (CLARA_Tb_subvis .eq. 3) clara_TbALL(ij,ik) = clara_tbClr(ij)
              
              ! Cloud-top pressure
              clara_ctp(ij,ik) = -999._wp
              clara_ctt(ij,ik) = -999._wp  
           endif

           ! #############################################################################
           ! 2) Semi-transparent clouds
           ! #############################################################################
           if (clara_cldFlag(ij,ik) .eq. 2) then
              ! a) Brightness temperature (using RTTOV w/o scattering)
              if (CLARA_Tb_semitrans .eq. 1) then
                 call rttov_simple(nLevels,nChannels,tautot(ij,ik,:),                    &
                                   phalf(ij,1:nLevels+1)/100._wp,thalf(ij,1:nLevels+1),  &
                                   pfull(ij,1:nLevels)/100._wp,at(ij,1:nLevels),         &
                                   qv(ij,1:nLevels),t2m(ij),p2m(ij)/100._wp,skt(ij),     &
                                   rt_emis(ij,1:nChannels),rt_satzen(ij),                &
                                   CLARA_minOpticalThickness,rt_h2otype(ij),rt_sfc(ij),  &
                                   CLARA_upperCloudTauLim,                               &
                                   tbTempALL(1:nChannels))  
              endif
              ! b) Brightness temperature (using RTTOV w/ scattering)
              if (CLARA_Tb_semitrans .eq. 2) then              
                 call get_Tb_rttov_clouds(1,nLevels,nChannels,pfull(ij,1:nLevels)/100.,  &
                                          at(ij,1:nLevels),qv(ij,1:nLevels),t2m(ij),     &
                                          p2m(ij)/100._wp,skt(ij),                       &
                                          rt_emis(ij,1:nChannels),rt_satzen(ij),         &
                                          rt_h2otype(ij),rt_sfc(ij),                     &
                                          rt_cloud(ij,ik,:,1:nLevels-1),                 &
                                          rt_cfrac(ij,ik,1:nLevels-1),                   &
                                          tbTempALL(1:nChannels))    
              endif  
              ! c) Brightness temperature (using simple CDK method)
              if (CLARA_Tb_semitrans .eq. 3) then
                 call allskyTb_CDK(nLevels,irrad(ij,1:nLevels),irradSfc(ij),             &
                                   emis_wv(ij,1:nLevels),emsfc_lw,tautot(ij,ik,:),       &
                                   cfrac(ij,ik,1:nLevels),TbTempALL(2))              
              endif
              clara_TbALL(ij,ik) = tbTempALL(2)

              ! Cloud-top pressure
              call clara_CTTH_ST_simple(nLevels, tautot(ij,ik,1:nLevels),                &
                                        phalf(ij,1:nLevels+1),thalf(ij,1:nLevels+1),     &
                                        clara_ctp(ij,ik),clara_ctt(ij,ik))                                           
           endif

           ! #############################################################################
           ! 3) Opaque clouds
           ! #############################################################################
           if (clara_cldFlag(ij,ik) .eq. 3) then
              ! a) Brightness temperature (using RTTOV w/o scattering)
              if (CLARA_Tb_opaque .eq. 1) then
                 call rttov_simple(nLevels,nChannels,tautot(ij,ik,:),                    &
                                   phalf(ij,1:nLevels+1)/100._wp,thalf(ij,1:nLevels+1),  &
                                   pfull(ij,1:nLevels)/100._wp,at(ij,1:nLevels),         &
                                   qv(ij,1:nLevels),t2m(ij),p2m(ij)/100._wp,skt(ij),     &
                                   rt_emis(ij,1:nChannels),rt_satzen(ij),                &
                                   CLARA_minOpticalThickness,rt_h2otype(ij),rt_sfc(ij),  &
                                   CLARA_upperCloudTauLim,                               &
                                   tbTempALL(1:nChannels))              
              endif
              ! b) Brightness temperature (using RTTOV w/ scattering)
              if (CLARA_Tb_opaque .eq. 2) then
                 call get_Tb_rttov_clouds(1,nLevels,nChannels,pfull(ij,1:nLevels)/100.,  &
                                          at(ij,1:nLevels),qv(ij,1:nLevels),t2m(ij),     &
                                          p2m(ij)/100._wp,skt(ij),                       &
                                          rt_emis(ij,1:nChannels),rt_satzen(ij),         &
                                          rt_h2otype(ij),rt_sfc(ij),                     &
                                          rt_cloud(ij,ik,:,1:nLevels-1),                 &
                                          rt_cfrac(ij,ik,1:nLevels-1),                   &
                                          tbTempALL(1:nChannels))    
              endif
              ! c) Brightness temperature (using simple CDK method)
              if (CLARA_Tb_opaque .eq. 3) then
                 call allskyTb_CDK(nLevels,irrad(ij,1:nLevels),irradSfc(ij),             &
                                   emis_wv(ij,1:nLevels),emsfc_lw,tautot(ij,ik,:),       &
                                   cfrac(ij,ik,1:nLevels),TbTempALL(2))
              endif
              clara_TbALL(ij,ik) = tbTempALL(2)      
              
              ! Cloud-top pressure
              call clara_CTTH_opaque(nLevels,clara_TbAll(ij,ik),inversion(ij,1:nLevels), &
                                     thalf(ij,1:nLevels+1),phalf(ij,1:nLevels+1),        &
                                     clara_ctp(ij,ik),clara_ctt(ij,ik))              
           endif
           ! Compute cloud-top height
           clara_cth(ij,ik) = cloud_top_height(nLevels,clara_ctp(ij,ik),              &
                                               at(ij,1:nLevels),phalf(ij,1:nLevels+1))           
           !WRITE(*,"(3i3,f15.8)") ij,ik,clara_phase(ij,ik),clara_size(ij,ik)
           endif
       enddo
    enddo
  end subroutine clara_subcolumn

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE clara_column
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine clara_column(nPoints,nSubCols,claraSC_tau,claraSC_ctp,claraSC_ctt,          &
                          claraSC_cth,claraSC_size,claraSC_phase,clara_tau,clara_ctp,    &
                          clara_cth,clara_ctt,clara_sizeLIQ,clara_sizeICE,clara_cfrac,   &
                          clara_IWP,clara_LWP,clara_fq)
      ! Inputs
      integer,intent(in) :: &
         nPoints,       & ! Number of points
         nSubCols         ! Number of subcolumns
      real(wp),intent(in),dimension(nPoints,nSubCols) :: &
         claraSC_tau,   & ! Optical depth
         claraSC_ctp,   & ! Clara cloud-top pressure (Pa) 
         claraSC_ctt,   & ! Clara cloud-top temperature
         claraSC_cth,   & ! Clara cloud-top height (m)
         claraSC_size     ! Clara retrieved particle size                        
      integer,intent(in),dimension(nPoints,nSubCols) :: &
         claraSC_phase    ! Clara cloud-phase

      ! Outputs
      real(wp),dimension(nPoints),intent(out) :: &
         clara_tau,     & ! CLARA retrieved optical depth
         clara_ctp,     & ! CLARA retrieved cloud-top pressure
         clara_ctt,     & ! CLARA retrieved cloud-top temperature
         clara_cth,     & ! CLARA retrieved cloud-top height
         clara_sizeLIQ, & ! CLARA retrieved particle size (liquid)
         clara_sizeICE, & ! CLARA retrieved particle size (ice)
         clara_cfrac,   & ! CLARA retrieved cloud fraction 
         clara_IWP,     & ! CLARA retrieved ice-water path
         clara_LWP        ! CLARA retrieved liquid-water path
      real(wp),dimension(nPoints,numCLARAtauBins,numCLARApresBins),intent(out) :: &
         clara_fq         ! Joint-histogram of cloud-top pressure and optical-depth
       
      ! Local variables
      integer :: ij,ik
      logical,dimension(nPoints,nSubCols)  :: retMask,cloudMask,waterCloudMask,iceCloudMask
      real(wp),dimension(nPoints)          :: clara_cfracLIQ,clara_cfracICE
      real(wp),dimension(nPoints,nSubCols) :: tauWRK,ctpWRK
      
      ! Parameters
      real(wp), parameter :: &
         LWP_conversion = 2._wp/3._wp * 1000._wp, &     
         ice_density    = 0.93_wp ! Liquid density is 1.           

    ! Initialize
    clara_tau(1:nPoints)     = -999._wp
    clara_ctp(1:nPoints)     = -999._wp
    clara_ctt(1:nPoints)     = -999._wp
    clara_cth(1:nPoints)     = -999._wp
    clara_sizeLIQ(1:nPoints) = -999._wp
    clara_sizeICE(1:nPoints) = -999._wp
    clara_cfrac(1:nPoints)   = -999._wp
    clara_IWP(1:nPoints)     = -999._wp
    clara_LWP(1:nPoints)     = -999._wp

    ! ####################################################################################
    ! Include only those pixels with successful retrievals in the statistics 
    ! ####################################################################################
    retMask(1:nPoints,1:nSubCols)   = claraSC_size(1:nPoints,1:nSubCols) > 0.
    cloudMask(1:nPoints,1:nSubCols) = claraSC_phase(1:nPoints,1:nSubCols) .ne.           &
                                   CLARA_phaseIsNone .and. retMask(1:nPoints,1:nSubCols)      
    waterCloudMask(1:nPoints,1:nSubCols) = claraSC_phase(1:nPoints,1:nSubCols) .eq.      &
                                   CLARA_phaseIsLiquid .and. retMask(1:nPoints,1:nSubCols)
    iceCloudMask(1:nPoints,1:nSubCols) = claraSC_phase(1:nPoints,1:nSubCols) .eq.        &
                                   CLARA_phaseIsIce .and. retMask(1:nPoints,1:nSubCols)   
                                      
    ! ####################################################################################
    ! Compute column cloud fraction (actually pixel count in each column that is cloudy for now)
    ! ####################################################################################
    clara_cfrac(1:nPoints)    = real(count(cloudMask,      dim = 2))    
    clara_cfracLIQ(1:nPoints) = real(count(waterCloudMask, dim = 2))
    clara_cfracICE(1:nPoints) = real(count(iceCloudMask,   dim = 2)) 
    ! Don't want to divide by 0, even though the sums will be 0 where the pixel counts are 0.  
    where(clara_cfrac(1:nPoints)    .eq. 0) clara_cfrac(1:nPoints)    = -1._wp
    where(clara_cfracLIQ(1:nPoints) .eq. 0) clara_cfracLIQ(1:nPoints) = -1._wp
    where(clara_cfracICE(1:nPoints) .eq. 0) clara_cfracICE(1:nPoints) = -1._wp

    ! ####################################################################################
    ! Column optical-depth
    ! ####################################################################################
    clara_tau(1:nPoints) = sum(claraSC_tau,mask=cloudMask,dim=2)/clara_cfrac(1:nPoints) 

    ! ####################################################################################
    ! Column liquid and ice water paths
    ! ####################################################################################
    clara_LWP(1:nPoints) = LWP_conversion*sum(claraSC_size*claraSC_tau, &
         mask=waterCloudMask,dim=2)/clara_cfracLIQ(1:nPoints)
    clara_IWP(1:nPoints) = LWP_conversion * ice_density*sum(claraSC_size*claraSC_tau,&
         mask=iceCloudMask,dim = 2)/clara_cfracICE(1:nPoints)
    
    ! ####################################################################################
    ! Column liquid and ice particle sizes
    ! ####################################################################################
    clara_sizeLIQ(1:nPoints) = sum(claraSC_size, mask = waterCloudMask, dim = 2) /       &
                                   clara_cfracLIQ(1:nPoints)
    clara_sizeICE(1:nPoints) = sum(claraSC_size, mask = iceCloudMask, dim = 2) /         &
                                   clara_cfracICE(1:nPoints)                          

    ! ####################################################################################
    ! Column cloud-top pressure, height and temperature
    ! ####################################################################################
    clara_ctp(1:nPoints) = sum(claraSC_ctp, mask = cloudMask, dim = 2) / &
                               max(1, count(cloudMask, dim = 2))
    clara_ctt(1:nPoints) = sum(claraSC_ctt, mask = cloudMask, dim = 2) / &
                               max(1, count(cloudMask, dim = 2))                      
    clara_cth(1:nPoints) = sum(claraSC_cth, mask = cloudMask, dim = 2) / &
                               max(1, count(cloudMask, dim = 2))

    ! ####################################################################################
    ! Now compute cloud fraction from number of cloud subcolumns. 
    ! ####################################################################################                         
    clara_cfrac(1:nPoints)    = max(0._wp, clara_cfrac(1:nPoints)/nSubcols)
    clara_cfracLIQ(1:nPoints) = max(0._wp, clara_cfracLIQ(1:nPoints)/nSubcols)
    clara_cfracICE(1:nPoints) = max(0._wp, clara_cfracICE(1:nPoints)/nSubcols)
 
    ! ####################################################################################
    ! Compute optical-depth vs. cloud-top pressure joint histogram
    ! ####################################################################################
    tauWRK(1:nPoints,1:nSubCols) = claraSC_tau(1:nPoints,1:nSubCols)
    ctpWRK(1:nPoints,1:nSubCols) = claraSC_ctp(1:nPoints,1:nSubCols)/100._wp
    do ij=1,nPoints
       ! Fill clear and optically thin subcolumns with fill
       where(.not. cloudMask(ij,1:nSubCols)) 
          tauWRK(ij,1:nSubCols) = -999._wp
          ctpWRK(ij,1:nSubCols) = -999._wp
       endwhere
       ! Joint histogram of tau/CTP
       call hist2D(tauWRK(ij,1:nSubCols),ctpWRK(ij,1:nSubCols),nSubCols,&
                   clara_histTau,numCLARAtauBins,&
                   clara_histPres,numCLARApresBins,&
                   clara_fq(ij,1:numCLARAtauBins,1:numCLARApresBins))    
    enddo     
    clara_fq(1:nPoints,1:numCLARAtauBins,1:numCLARApresBins)=clara_fq/nSubCols   
         
  end subroutine clara_column
  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE clara_CTTH_opaque
  !
  ! Find the cloud top of opaque cloud, Find this layer by looking from the BOTTOM UP 
  ! (surface->troposphere) as is done in PPS. This simple approach does not work for 
  ! semi-transparent clouds due to the contamination from the lower layers. Clara
  ! semi-transparent clouds are retrieved using clara_CTTH_ST
  !
  ! Originally coded by Salomon Eliasson (Salomon.Eliasson@smhi.se)
  ! Modified for COSPv2.0 by Dustin Swales (dustin.swales@noaa.gov)
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine clara_CTTH_opaque(nLevels,Tb,inv,temp,pres,ctp,ctt)
     ! Inputs
     integer,intent(in) :: &
        nLevels ! Number of vertical levels
     integer,dimension(nLevels),intent(in) :: &
        inv     ! Indices for levels with inversions
     real(wp),intent(in) :: &
        Tb      ! Brightness temperature (all-sky)
     real(wp),dimension(nLevels+1),intent(in) :: &
        temp, & ! Temperature profile   
        pres    ! Pressure profile
        
     ! Outputs
     real(wp),intent(out) :: ctt,ctp
     
     ! Local variables
     integer :: trop_lev,ii,ij,flagged
     real(wp) :: T_inv,temp_T1,temp_T2,T1,T2,P1,P2,fractionalDeltaT
     integer,dimension(nLevels+1) :: multi
     logical,dimension(nLevels+1) :: adjustedToInversion
     
     ! Parameters
     real(wp), parameter :: &
        DistanceToInversion = 0.5   ! PPS 2014 places a cloud in an inversion if Tb is 
                                    ! within 0.5K of the inversion temperature.

     ! Initialize
     multi(1:nLevels+1)               = 0
     adjustedToInversion(1:nLevels+1) = .false.

     ! Loop over the levels to find the model level that CONTAINS the all-sky temperature.
     ! Loop from the model surface to the top of the atmosphere. Keep an eye out for 
     ! ambiguous clouds, i.e., a solution that could give several different cloud heights

     ! Find the tropopause level (highest inversion)
     trop_lev = MINVAL(inv,mask=inv .gt. 0)
     
     ! Get the temperature of the lowest inversion
     ii    = 1            
     T_inv = temp(inv(ii)) 

     do ij = nLevels,1,-1
       if (ij .ge. trop_lev) then
          ! Find between which model levels the all sky brightness temperature is found by
          ! looping from the Earth's surface up through the model levels within the 
          ! troposphere. Include the bottom and top of the troposphere. SE: Don't use 
          ! thalf to find temperature in the rest of the profile.
          temp_T2 = temp(ij+1)
          temp_T1 = temp(ij)

          ! Find the nearest inversion above the current layer
          if (ij .lt. inv(ii)) then
             ! There was no cloud close enough to the current inversion. Find the next 
             ! inversion higher up in the atmosphere 
             ii    = ii+1
             T_inv = temp(inv(ii))
          end if

          if ((ABS(Tb - temp(ij)) .lt. DistanceToInversion) .and.                        &
              (ij .eq. inv(ii)+1)) then ! and is actually an inversion
             ! The +1 to inversion gives the lower (colder) boundary of the inversion layer

             ! If Tb is within .5K of the temperature profile and that level is the 
             ! location of a temperature profile, allocate the cloud there.
             multi(ij)               = ij
             adjustedToInversion(ij) =.true.

          else
             ! Normal way of finding the cloud top height by comparing the all-sky  
             ! temperature to the model layers
             ! *) Non-inversion. (Tb is colder than the lower level but warmer than the 
             !    layer above)
             if (((Tb .le. temp_T1) .AND. (Tb .gt. temp_T2)) .or. &
             ! *) Inversion. (Tb is warmer than the lower level but colder than the layer
             !    above)
                 ((Tb .ge. temp_T1) .AND. (Tb .lt. temp_T2))) then
                multi(ij) = ij
             end if ! --- Tb_IR vs. model level temperatures
          end if ! --- If close to an inversion
       end if ! ---- if within troposphere
    enddo !--- end loop over levels     

    call get_boundaries_and_flags(nLevels,multi,adjustedToInversion,temp,pres,Tb,        &
                                  flagged,T1,T2,P1,P2)

    ! Cloud-top temperature
    ! Get the "exact" cloud location relative to T1 and T2
    if (flagged .eq. 2 )  then ! Too warm
       ! T2 corresponds to the is the warmest temperature in the profile in this case
       FractionalDeltaT = 0
       ctt              = T2
    else if  (flagged .EQ. 5 ) then     ! adjustedToInversion
       ! T1 is the local minimum of the inversion
       FractionalDeltaT = 1
       ctt              = T1
    else if (T1 .EQ. T2) then
       ! Unlikely event
       !stop "Does this ever happen?"
       FractionalDeltaT = 0
       ctt              = T1
    else
       ! The normal approach              
       FractionalDeltaT = ABS(T2 - Tb)/ABS(T2 - T1)
       ctt              = Tb
    end if
    
    ! Cloud-top pressure
    if (P2 .eq. 0) then
       ctp  = P1
    else
       ctp = EXP( LOG(P2) + FractionalDeltaT*( LOG( P1/P2 ) ) )    
    end if
                                                            
  end subroutine clara_CTTH_opaque
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! END SUBROUTINE clara_CTTH_opaque
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE clara_CTTH_ST_simple
  !
  ! Simple approach using only the documented biases for semi-transparent clouds. The 
  ! biases are from the "PPS algoritms scientific and validation report for thecloud 
  ! product processes of the NWC/PPS (27March 2015)"
  ! 
  ! Originally coded by Salomon Eliasson (Salomon.Eliasson@smhi.se)
  ! Modified for COSPv2.0 by Dustin Swales (dustin.swales@noaa.gov)
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  subroutine clara_CTTH_ST_simple(nlev, tau, P, T, ctp, ctt)
     ! Inputs 
     integer, intent(in) :: &
        nlev   ! Number of vertical levels
     real(wp),intent(in),dimension(nlev) :: &
        tau    ! Optical thickness
     real(wp),intent(in),dimension(nlev+1) :: &
        P,   & ! Pressure
        T      ! Temperature   
     ! Outputs    
     real(wp),intent(out) :: &
        ctp(1), & ! Cloud-top pressure
        ctt(1)    ! Cloud-top temperature

    ! Parameters low, middle, high
    real(wp),dimension(3),parameter :: &
       PPS_bias=(/154,947,-1405/)
    real(wp),parameter :: &
       lw_cl_lim = 3000._wp, & ! fixme: I'm just guessing these height. yet to find the real limits
       hi_cl_lim = 7000._wp    !  

    ! Local variables
    real(wp),dimension(nlev)   :: tau_int,rho
    real(wp),dimension(nlev+1) :: Z
    real(wp)                   :: dz, height,newctp(1),newctt(1),tau2,tau1,tau_ct_ratio,ctt_,ctp_
    integer                    :: L2, L1, inl, inl2

    ! Initialize
    newctp       = 0._wp
    newctt       = 0._wp
    Z(1:nlev+1)  = 0._wp
    ctp_         = -999._wp
    ctt_         = -999._wp
    ctp          = -999._wp
    ctt          = -999._wp
    tau2         = 0._wp
    tau1         = 0._wp
    tau_int      = 0._wp
    tau_CT_ratio = 0._wp
    L2 = 0
    L1 = 0    
    
    ! Compute density
    rho = P/(rd*T)

    ! Get the ctth using the modis approach
    ! Loop from the top of the atmosphere down.
    do inl = 1, nlev

       ! gather tau from the top down. copy it over
       if (tau(inl) .le. 0) cycle 
       tau_int(inl) = tau(inl)

       ! Check if the sum of the tau_int (so far) is higher than
       ! the boundary we want. If it is, find the model layers
       ! between which the integrated tau = tau_eqCT

       ! this level is the upper bounds of tau (tau(L2)<tau1<tau(L1))
       tau1 = sum( tau_int(1:inl) )

       if (tau1 .GE. CLARA_upperCloudTaulim) then

          ! L1 = the layer boundary below where tau_eqCT is exceeded
          L1 = inl +1 

          ! Look BACK UP the profile to find the closest level for
          !  which the integrated optical depth is less than
          !  tau_eqCT. Once this is found, add +1 to inl so that we
          !  pick <this> model layers' LOWER INTERFACE

          do inl2 = L1-2,1,-1
             if (tau_int(inl2)>0) then
                !doing this once
                tau2 = tau1 - tau_int(L1-1)
                L2 = inl2+1
                exit
             end if
          end do ! end looking back

          ! If the first level from the top down has tau >tau_eq
          ! then set L2 to L1 -1 (enough cloud is contained in
          ! this layer alone)
          if (L2 .eq. 0) then
             tau2 = 0.
             L2=L1-1
          end if

          ! Find the CTH between these layers using linear
          ! interpolation of log(pressure) according to 
          ! where tau_eqCT falls compared to the layer boundaries

          if ( (tau2 .gt. CLARA_upperCloudTauLim) .or. (tau2 .gt. tau1) ) then
             !stop "This has to be false: (tau2 .gt. tau_eq)&
             !     & .or. (tau2 .gt. tau1)"
          end if

          tau_CT_ratio = (CLARA_upperCloudTauLim-tau2)/(tau1-tau2)

          ! The pressures and Temperatures are at the layer interfaces
          ctp_ = EXP( LOG(P(L1)) + tau_CT_ratio*( LOG( P(L1)/P(L2) ) ) )
          ctt_ = T(L2) + tau_CT_ratio*( T(L1)-T(L2) )

          ctp = ctp_
          ctt = ctt_

       end if ! tau_int > tau_eqCT

       if (L2 .gt. 0) then
          exit
       end if
    end do ! loop levels
    if ( (sum(tau_int) .ge. CLARA_minOpticalThickness) .and. (sum(tau_int) .lt. CLARA_upperCloudTauLim) ) then
       ! Then there is a cloud that is semi-transparent and
       ! therefore we'll retrieve as far down as there are
       ! clouds and call that the equivalent cloud top.

       do inl = nlev,1,-1
          ! find the bottom of the cloud by looping up from
          !  the surface

          if (tau_int(inl)>0) then
             ctp = P(inl+1)
             ctt = T(inl+1)
             tau1=sum(tau_int)
             tau2=0.0
             tau_CT_ratio = 1
             exit
          end if
       end do

    end if ! if cloud is transparent
    ! END MODISish retrieval

    ! Make a height vector (P and T are at the intefaces, rho is in the middle)
    do inl = nlev,1,-1
       dz = ( P(inl+1)-P(inl) ) / ( rho(inl) * grav)
       Z(inl) = Z(inl+1)+dz
    end do

    ! Find the new ctp from the ground up
    do inl = nlev,1,-1
       if (Z(inl) < lw_cl_lim) then
          height = Z(inl) + PPS_bias(1)
       elseif (Z(inl) >= lw_cl_lim .and. Z(inl) <= hi_cl_lim) then
          height = Z(inl) + PPS_bias(2)
       else
          height = Z(inl) + PPS_bias(3)
       end if

       ! Kill the loop when I have the right height
       if (Z(inl) > height) then
          newctp = P(inl)
          newctt = T(inl)
          EXIT
       end if
    end do

  end subroutine clara_CTTH_ST_simple  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! END SUBROUTINE clara_CTTH_ST_simple
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE get_boundaries_and_flags
  !
  ! 1) Assign flags to the cloud top height allocations, to get some sort of idea of the 
  !    level of uncertainty
  ! 2) Get the temperature and pressure at the level interfaces surrounding the level that
  !    we think the cloud is located
  !
  ! Originally coded by Salomon Eliasson (Salomon.Eliasson@smhi.se)
  ! Modified for COSPv2.0 by Dustin Swales (dustin.swales@noaa.gov)  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  subroutine get_boundaries_and_flags(nLevels,multi,adjustedToInversion,thalf,phalf,Tb,  &
                                      flagged,T1,T2,P1,P2)

      ! Inputs
      integer,intent(in) :: &
          nLevels ! Number of vertical levels
      integer,dimension(nLevels+1) :: &
          multi
      logical,dimension(nLevels+1),intent(in) :: &
          adjustedToInversion
      real(wp),intent(in) :: &
          Tb       ! Brightness temperature
      real(wp),dimension(nLevels+1),intent(in) :: &
          phalf, & ! Pressure at model half levels
          thalf    ! Temperature at model half levels
   
      ! Outputs
      real(wp),intent(out) :: T1,T2,P1,P2
      integer,intent(out)  :: flagged
   
      ! Local variables
      integer,dimension(nLevels+1) :: cloud_layers
      real(wp),dimension(nLevels+1) :: mindiff
      integer :: ii,jj,ij,N_clouds,clevel
      logical :: adjacent_clouds
      real :: minmindiff,eps

      ! Count how many model layers where it would be valid to place the cloud by using Tb
      N_clouds = count(multi .gt. 0)

      ! Compact the vector of indices where clouds where found  
      jj=1
      do ij = nLevels+1,1,-1
         if (multi(ij) .gt. 0) then
            cloud_layers(jj) = multi(ij)
            jj = jj+1
         end if
      end do

      ! Assign flags
      if (N_clouds .gt. 0) then
         ! PICKING MODEL LEVEL
         ! Picking the solution at the lowest altitude in accordance with section 4.1.2.1 
         ! in CTTH-PGE03 v3.1
         clevel = cloud_layers(1) ! level middle (including half-levels at tropopause and surface)

         if (adjustedToInversion(clevel)) then
            ! Flag if cloud level is found by being close enough to the inversion
            flagged = 5

         elseif (N_clouds .EQ. 1) then
            ! The cloud can be safely allocated to just one model level
            flagged = 1

         ELSEIF (N_clouds .GT. 1) THEN
            ! The equivalent cloud can be allocated to several adjacent layers based on Tb and t11. 

            ! Test to see if the multiple cloud layers are both/all adjacent to each other (flagged=3) or not (flagged=4)
            adjacent_clouds = .true.
            do ii = 1,N_clouds-1
               adjacent_clouds = adjacent_clouds .and. & 
                    ((cloud_layers(ii)-cloud_layers(ii+1) .eq. 1) .or. &! levels adjacent
                    (cloud_layers(ii)-cloud_layers(ii+1) .eq. cloud_layers(ii+1))) ! no more clouds
            end do

            if (adjacent_clouds) then
               ! The cloud can be two or more adjacent layers.
               flagged = 3 ! somewhat ambiguous 

            else
               ! The cloud can be placed at quite different parts of the model column. As 
               ! long as this technique is used, I'll just pick the lowest layer and be 
               ! honest about the ambiguity. 
               flagged = 4 ! very ambiguous
            end if
         end if !-flags for Nclouds>0
      else
         ! No cloud was found, although one should be present. This is because Tb is 
         ! actually warmer than the warmest model layer within the tropopause (sometimes 
         ! the warmest temperature is in the model atmosphere). (FIXME: How can this be?)

         mindiff(:) = ABS(thalf-Tb)
         minmindiff = MINVAL(mindiff)
         eps = EPSILON(mindiff)
         do ij = nLevels+1,1,-1
            if (abs(mindiff(ij) - minmindiff) .LT. eps) then
               ! exit when you have found the layer
               exit
            end if
         end do

         clevel = ij
         flagged = 2 ! Tb is too warm
      endif

      ! Get the interfaces
      T2 = thalf(clevel)
      P2 = phalf(clevel)
      IF (clevel .EQ. nLevels+1) THEN
         ! Then it is the surface layer
         T1=T2
         P1=P2
      ELSE
         T1 = thalf(clevel+1)
         P1 = phalf(clevel+1)
      END IF

  end subroutine get_boundaries_and_flags
  ! ######################################################################################  
  ! END SUBROUTINE get_boundaries_and_flags
  ! ###################################################################################### 
   
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! FUNCTION find_temperature_inversion
  ! The purpose of this is to find all of the temperature inversions in each model grid. 
  ! These are used in some of the cloud tests and the tropopause is also detected using 
  ! this function. 
  !
  ! - The inversions are placed at the local minima in the profile. 
  ! - As a last step the highest inversion, assumed to be the tropopause is
  !   bumped up a level.
  ! - The tropopause must exist between 500-50hPa if none is found,
  !   the tropopause is place at the level closest to the 50hPa level
  !
  ! Originally coded by Salomon Eliasson (Salomon.Eliasson@smhi.se)
  ! Modified for COSPv2.0 by Dustin Swales (dustin.swales@noaa.gov)
  !  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function find_temperature_inversion(nLevels,temp,pres)
     ! Inputs
     integer,intent(in) :: &
        nLevels
     real(wp),dimension(nLevels+1),intent(in) :: &
        temp,   & ! Temperature at model full levels(K)
        pres      ! Pressure at model half levels (hPa)   
     
     ! Outputs
     integer,dimension(nLevels) :: &
        find_temperature_inversion     ! Layer indices for inversions   
     
     ! Local Variables
     logical  :: previousNotInversion,tropfail
     integer  :: ij,layer,ninv
     real(wp) :: T_comp,nearestToMin
     integer,dimension(nLevels) :: inv
     
     ! Parameters
     real(wp), parameter :: &
        stratosphereGradient = 2.5,   & ! [K] 
        MaxTropPres          = 50000, & ! [Pa]
        MinTropPres          = 5000     ! [Pa]  

    ! Loop from the surface up to the tropopause and look for dt/dz > 0
    inv(1:nLevels)        = 0
    previousNotInversion  = .true.
    ninv                  = 0
    layer                 = 0
    T_comp                = temp(nLevels+1)
    tropfail              = .false.
    do ij=nLevels,1,-1
       if (temp(ij) .gt. T_comp .and. previousNotInversion .and. .not. pres(ij) .lt. MinTropPres) then
          ninv                 = ninv+1
          layer                = ij
          inv(ninv)            = layer
          previousNotInversion = .FALSE.          
       elseif (temp(ij) .lt. T_comp) then
          previousNotInversion = .TRUE.
       end if       
       T_comp = temp(ij)
    enddo

    ! Tropopause was not found. Assign to the highest level allowed.
    if (ninv .eq. 0) then
       tropfail=.true.
    elseif (pres(inv(ninv)) .ge. MaxTropPres) then 
       tropfail=.true.
    endif
    if (tropfail) then
       ninv = ninv + 1
       ! Need to find the model level closest to the lowest pressure allowed and call that
       ! the tropopause (although it is not an inversion)
       nearestToMin = minval(abs(pres-minTropPres))

       do ij=1,nLevels
          if ((abs(pres(ij)-minTropPres) - nearestToMin) .lt. EPSILON(nearestToMin)) then
             exit
          end if
       enddo
       layer     = ij+1
       inv(ninv) = layer
    end if
    
    ! Place tropopause (last entry) outside the troposphere. Therefore, decrease model 
    ! layer by one
    layer     = layer-1
    inv(ninv) = layer
    if (temp(layer-1) - temp(layer) .lt. stratosphereGradient ) then
        ! Ad hoc fix for clouds just above the found tropopause (The last "layer"). Check 
        ! that the temperature gradient is high enough to really be in the stratosphere. 
        ! Sometimes clouds are in the tropopause region with slightly increasing 
        ! temperature. Bump up the tropopause if need be.
        layer     = layer-1
        inv(ninv) = layer
    end if
    find_temperature_inversion = inv
  end function find_temperature_inversion
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! END FUNCTION find_temperature_inversion
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! FUNCTION expand_profiles
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function expand_profile(nPoints,nLevels,temp,pfull,phalf)
     ! Inputs
     integer,intent(in) :: &
        nPoints, & ! Number of points
        nLevels    ! Number of vertical levels
     real(wp),dimension(nPoints,nLevels),intent(in) :: &
        temp,    & ! Corrected temperature profile   
        pfull      ! Pressure at model full levels
     real(wp),dimension(nPoints,nLevels+1),intent(in) :: &
        phalf      ! Pressure at model half levels   

     ! Outputs
     real(wp),dimension(nPoints,nLevels+1) :: &
        expand_profile  ! Corrected temperature profile at model half-levels
     
     ! Local variables
     integer :: lt,ilev
     real(wp) :: ratio

     ! Initialize
     ratio = 0._wp
     expand_profile(1:nPoints,1:nLevels+1) = 0._wp

    ! Loop to get the linearly interpolated temperature at the layer interfaces. The 
    ! interfaces might not be exactly in between two midpoints so linearly interpolate 
    ! according to how the pressure interfaces relate to the pressure midpoints.
    do lt = 1, nPoints
       do ilev = 2, nLevels
          ratio = log(phalf(lt,ilev)/pfull(lt,ilev-1))/&
                  log(pfull(lt,ilev)/pfull(lt,ilev-1))
          expand_profile(lt,ilev) = temp(lt,ilev-1) + ratio*(temp(lt,ilev)-temp(lt,ilev-1));
       end do

       ! Special for the highest model layer and lowest model layer interfaces
       expand_profile(lt,1) = temp(lt,1) - (expand_profile(lt,2)-temp(lt,1));
       expand_profile(lt,nLevels+1) = expand_profile(lt,nLevels);
    end do  
  
  end function expand_profile
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! END FUNCTION expand_profile
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE clearSkyTb_CDK
  ! Compute water vapor continuum emissivity this treatment follows Schwarkzopf 
  ! and Ramasamy JGR 1999,vol 104, pages 9467-9499. The emissivity is calculated at a 
  ! wavenumber of 955 cm-1, or 10.47 microns. 
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  subroutine clearSkyTb_CDK(nPoints,nLevels,phalf,pfull,qv,at,skt,emsfc_lw,irradiance,   &
                            irradianceSfc,dem_wv,atCor,Tb)
    ! Inputs
    integer,intent(in) :: &
       nPoints,nLevels
    real(wp),intent(in) :: &
       emsfc_lw
    real(wp),intent(in),dimension(nPoints) :: &
       skt   
    real(wp),intent(in),dimension(nPoints,nLevels) :: &
       pfull,qv,at
    real(wp),intent(in),dimension(nPoints,nLevels+1) :: &
       phalf   

    ! Outputs
    real(wp), intent(out),dimension(nPoints) :: &
       Tb,           & ! Clear-sky brightness temperature
       irradianceSfc   ! Surface irradiance
    real(wp),intent(out),dimension(nPoints,nLevels) :: &
       dem_wv,       & ! Water vapor emissivity
       atCor,        & ! Corrected temperature profile
       irradiance      ! Irradiance
       
    ! Need to add in more outputs. Anything in brightness_temperature.F90 in the "inter"
    ! derived type needs to be output!!!   

    ! Local variables
    integer :: il
    real(wp),dimension(nPoints) :: dp,atmrho,rvh2o,wk,rhoave,rh2os,rfrgn,tmpexp,tau_wv,  &
                                   fluxtoa_clr,trans_above_clear,corrected_irradiance
    
    ! Compute water vapor emissivity.
    do il=1,nLevels
        dp(1:nPoints)     = 10._wp*(phalf(1:nPoints,il+1)-phalf(1:nPoints,il))
        atmrho(1:nPoints) = dp(1:nPoints)/(grav*100._wp)        
        rvh2o(1:nPoints)  = qv(1:npoints,il)*amd/amw
        wk(1:nPoints)     = rvh2o(1:npoints)*avo*atmrho/amd
        rhoave(1:nPoints) = ((pfull(1:npoints,il)*10._wp)/pstd)*(clara_t0/at(1:nPoints,il))
        rh2os(1:nPoints)  = rvh2o(1:nPoints)*rhoave(1:nPoints)
        rfrgn(1:nPoints)  = rhoave(1:nPoints)-rh2os(1:nPoints)
        tmpexp(1:nPoints) = EXP(-0.02_wp*(at(1:nPoints,il)-clara_t0)) 
        tau_wv(1:nPoints) = wk(1:nPoints)*1.e-20*((0.0224697_wp*rh2os(1:nPoints)*     &
                              tmpexp(1:nPoints)) + (3.41817e-07*rfrgn(1:nPoints)))*0.98_wp
 
        ! Water vapor emissivity 
        dem_wv(1:nPoints,il) = 1._wp - exp( -1._wp * tau_wv(1:npoints))
    enddo
        
    fluxtoa_clr(1:nPoints)       = 0._wp
    trans_above_clear(1:nPoints) = 1._wp
    do il=1,nLevels
        ! Black body emission from model vertical layer temperature.
        irradiance(1:nPoints,il) = 1._wp/(EXP(tb4_rad/at(1:nPoints,il))-1.)    

        ! Increase TOA flux by flux emitted from layer times total transmittance in layers above          
        fluxtoa_clr(1:nPoints) = fluxtoa_clr(1:nPoints) +  &
            trans_above_clear(1:nPoints)*dem_wv(1:nPoints,il)*irradiance(1:nPoints,il)                

        ! Brightness temperature decreased by absorption by layer water vapor
        corrected_irradiance(1:nPoints) = irradiance(1:nPoints,il)*(1._wp-dem_wv(1:nPoints,il))

        ! Corrected temperature profile
        atCor(1:nPoints,il) = tb4_rad*(1./LOG(1.+1./corrected_irradiance(1:nPoints)))

        ! Update trans_layers_above with transmissivity from this layer for next time around loop
        trans_above_clear(1:nPoints) = trans_above_clear(1:nPoints)*(1._wp-dem_wv(1:nPoints,il))
    enddo

    ! Add in surface emission
    irradianceSfc(1:nPoints)  = 1._wp/( exp(tb4_rad/skt(1:nPoints)) - 1._wp )
    fluxtoa_clr(1:npoints) = fluxtoa_clr(1:nPoints)+emsfc_lw*irradianceSfc(1:nPoints)*   &
                             trans_above_clear(1:nPoints)

    ! Clear Sky brightness temperature
    Tb(1:npoints) = tb4_rad/(log(1._wp+(1._wp/fluxtoa_clr(1:nPoints))))  
    
  end subroutine clearSkyTb_CDK
  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE allSkyTb_CDK  
  !
  ! Function to calculate the allsky temperature using the same method as the ISCCP 
  ! simulator 
  !
  ! Originally coded by Salomon Eliasson (Salomon.Eliasson@smhi.se)
  ! Modified for COSPv2.0 by Dustin Swales (dustin.swales@noaa.gov)
  !  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  SUBROUTINE  allskyTb_CDK(nLevels,irrad,irradSfc,emis_wv,emsfc_lw,tau,cfrac,Tb)

    ! Inputs
    integer,intent(in) :: &
       nLevels     ! Number of vertical levels
    real(wp),intent(in),dimension(nLevels) :: &
       irrad,    & ! Irradiance
       emis_wv,  & ! Water vapor emissivity
       tau,      & ! Cloud optical depth (cumulative from TOA-2-SFC)   
       cfrac       ! Cloud fraction (subcolumn cloud mask, where 1=cloudy,0=clear)
    real(wp),intent(in) :: &
       irradSfc, & ! Surface irradiance    
       emsfc_lw    ! Surface emissivity   
       
    ! Outputs
    real(wp),intent(out) :: &
       Tb ! All-sky brightness temperature (K)
       
    ! Local variables
    real(wp) :: fluxtoa_all, trans_above_all,EMS
    real(wp),dimension(nLevels) :: cloud_emis
    integer  :: inl

    ! Initialize
    fluxtoa_all     = 0._wp
    trans_above_all = 1._wp
    EMS             = 0._wp
    Tb              = 0._wp
    cloud_emis      = 1._wp - exp(-0.5_wp*tau)
    DO inl = 1, nLevels

       ! EMS=emissivity in each sub-grid column in each grid influenced by water vapour 
       ! and clouds (emis originally calculated from tau_vis *0.5), and is looped over 
       ! each level and used to calculate the cumulative 10.5 micron TOA flux, hence 
       ! brightness temperature for each sub grid column.
       IF (cfrac(inl) .ne. 0) THEN      
          ! Cloudy
          EMS = 1._wp - (1._wp - emis_wv(inl))*(1._wp - cloud_emis(inl))
       ELSE  
          ! Cloud free
          EMS = 1._wp - (1._wp - emis_wv(inl))
       END IF
       
       trans_above_all = trans_above_all* (1._wp - EMS)
       fluxtoa_all     = fluxtoa_all + trans_above_all*EMS*irrad(inl)
    END DO

    ! Surface contribution to TOA all-sky brightness temperature
    fluxtoa_all = fluxtoa_all + trans_above_all*emsfc_lw*irradSfc

    ! and here's the all sky brightness temperature by this approach   
    Tb = tb4_rad*(1._wp/LOG(1._wp + 1._wp/fluxtoa_all))

  END SUBROUTINE allskyTb_CDK  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! END SUBROUTINE allskyTb_CDK
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! FUNCTION get_g_nir
  ! Compute asymmetry parameter at channel 383 using provided radius and phase.
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elemental function get_g_nir(re,phase)
    ! INPUTS
    integer, intent(in) :: phase
    real(wp),intent(in) :: re
    
    ! Outputs
    real(wp) :: get_g_nir
    
    ! Local variables
    integer :: xi(2)
    integer,parameter :: &
       nreIce = 9, & ! Number of ice radii in LUT
       nreLiq = 8    ! Number of liquid radii in LUT
    real(wp),dimension(nreIce),parameter :: &
       re_ice = (/5.00,7.07,10.00,14.15,20.00,28.28,40.00,56.58,80.00/),                 &
! djs2015: Use LUT without delta approximation (email from Jan Fokke on dec 11,2015)
!       g0_ice = (/0.750507,0.779607,0.810919,0.843880,0.875778,0.903978,0.922744,        &
!                  0.921151,0.904428/)
       g0_ice = (/0.756000, 0.784000, 0.814000, 0.846000, 0.878000, 0.907000, 0.932000,  &
                  0.948000, 0.956000/)
    real(wp),dimension(nreLiq),parameter :: &   
       re_liq = (/3.00,4.25, 6.00, 8.50,12.00,17.00,24.00,34.00/),                       &
       g0_liq = (/0.806269,0.789406,0.767975,0.779810,0.813295,0.843947,0.868019,        &
                  0.889109/)

    ! Liquid    
    if (phase .eq. 1) then
       ! Find interpolation bounds
       xi = [maxloc(re_liq-re,re_liq-re .le. 0),maxloc(re_liq-re,re_liq-re .le. 0)+1]
       if (minval(abs(re_liq-re)) .eq. 0) xi=[xi(1),xi(1)]
       ! Interpolate
       if (re .ge. CLARA_re_water_min .and. re .le. CLARA_re_water_max) then
          get_g_nir = g0_liq(xi(1))+(g0_liq(xi(2))-g0_liq(xi(1)))*(re-re_liq(xi(1)))/    &
                                    (re_liq(xi(2))-re_liq(xi(1))) 
       endif                             
       if (xi(1) .eq. xi(2)) get_g_nir = g0_liq(xi(1))             ! On table node                                           
       if (re .lt. CLARA_re_water_min)  get_g_nir = g0_liq(1)      ! Re too small
       if (re .gt. CLARA_re_water_max)  get_g_nir = g0_liq(nreLiq) ! Re too big
    endif
    ! Ice
    if (phase .eq. 2) then
       ! Find interpolation bounds
       xi = [maxloc(re_ice-re,re_ice-re .le. 0),maxloc(re_ice-re,re_ice-re .le. 0)+1]
       if (minval(abs(re_ice-re)) .eq. 0) xi=[xi(1),xi(1)]
       ! Interpolate
       if (re .ge. CLARA_re_ice_min .and. re .le. CLARA_re_ice_max) then
          get_g_nir = g0_ice(xi(1))+(g0_ice(xi(2))-g0_ice(xi(1)))*(re-re_ice(xi(1)))/    &
                                    (re_ice(xi(2))-re_ice(xi(1)))  
       endif                             
       if (xi(1) .eq. xi(2)) get_g_nir = g0_ice(xi(1))           ! On table node                                           
       if (re .lt. CLARA_re_ice_min)  get_g_nir = g0_ice(1)      ! Re too small
       if (re .gt. CLARA_re_ice_max)  get_g_nir = g0_ice(nreIce) ! Re too big
    endif    
    
  end function get_g_nir
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! END FUNCTION get_g_nir
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! FUNCTION get_ssa_nir
  ! Compute single-scattering albedo at channel 383 for provided radius and phase.
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elemental function get_ssa_nir(re,phase)
    ! INPUTS
    integer, intent(in) :: phase
    real(wp),intent(in) :: re
    
    ! Outputs
    real(wp) :: get_ssa_nir
    
    ! Local variables
    integer :: xi(2)
    integer,parameter :: &
       nreIce = 9, & ! Number of ice radii in LUT
       nreLiq = 8    ! Number of liquid radii in LUT
    real(wp),dimension(9),parameter :: &
       re_ice = (/5.00,7.07,10.00,14.15,20.00,28.28,40.00,56.58,80.00/),                 &
! djs2015: Use LUT without delta approximation (email from Jan Fokke on dec 11,2015)
!       w0_ice = (/0.882276,0.844964,0.800201,0.749888,0.696207,0.640230,0.575341,        &
!                  0.476380,0.367772/)
       w0_ice = (/0.882300, 0.845000, 0.800300, 0.750300, 0.698000, 0.647700, 0.606200,  &
                  0.577700, 0.558300/)                     
    real(wp),dimension(8),parameter :: &   
       re_liq = (/3.00,4.25, 6.00, 8.50,12.00,17.00,24.00,34.00/),                       &
       w0_liq = (/0.976040,0.965643,0.945708,0.919468,0.890786,0.857294,0.817537,        &
                  0.770678/)
    
    ! Liquid    
    if (phase .eq. 1) then
       ! Find interpolation bounds
       xi = [maxloc(re_liq-re,re_liq-re .le. 0),maxloc(re_liq-re,re_liq-re .le. 0)+1]
       if (minval(abs(re_liq-re)) .eq. 0) xi=[xi(1),xi(1)]
       ! Interpolate
       if (re .ge. CLARA_re_water_min .and. re .le. CLARA_re_water_max) then
          get_ssa_nir = w0_liq(xi(1))+(w0_liq(xi(2))-w0_liq(xi(1)))*(re-re_liq(xi(1)))/  &
                                      (re_liq(xi(2))-re_liq(xi(1))) 
       endif                             
       if (xi(1) .eq. xi(2)) get_ssa_nir = w0_liq(xi(1))             ! On table node                                           
       if (re .lt. CLARA_re_water_min)  get_ssa_nir = w0_liq(1)      ! Re too small
       if (re .gt. CLARA_re_water_max)  get_ssa_nir = w0_liq(nreLiq) ! Re too big
    endif
    ! Ice
    if (phase .eq. 2) then
       ! Find interpolation bounds
       xi = [maxloc(re_ice-re,re_ice-re .le. 0),maxloc(re_ice-re,re_ice-re .le. 0)+1]
       if (minval(abs(re_ice-re)) .eq. 0) xi=[xi(1),xi(1)]
       ! Interpolate
       if (re .ge. CLARA_re_ice_min .and. re .le. CLARA_re_ice_max) then
          get_ssa_nir = w0_ice(xi(1))+(w0_ice(xi(2))-w0_ice(xi(1)))*(re-re_ice(xi(1)))/  &
                                      (re_ice(xi(2))-re_ice(xi(1)))  
       endif                             
       if (xi(1) .eq. xi(2)) get_ssa_nir = w0_ice(xi(1))           ! On table node                                           
       if (re .lt. CLARA_re_ice_min)  get_ssa_nir = w0_ice(1)      ! Re too small
       if (re .gt. CLARA_re_ice_max)  get_ssa_nir = w0_ice(nreIce) ! Re too big
    endif        
    
  end function get_ssa_nir
  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! FUNCTION weight_by_extinction
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function weight_by_extinction(nLevels,tauIncrement, f, tauLimit) 
    ! INPUTS
    integer, intent(in)                    :: nLevels
    real(wp),intent(in),dimension(nLevels) :: tauIncrement, f
    real(wp),intent(in)                    :: tauLimit
    ! OUTPUTS
    real(wp)                               :: weight_by_extinction
    ! LOCAL VARIABLES
    real(wp)                               :: deltaX, totalTau, totalProduct
    integer                                :: i 
    
    ! Find the extinction-weighted value of f(tau), assuming constant f within each layer
    totalTau = 0._wp; totalProduct = 0._wp
    do i = 1, size(tauIncrement)
      if(totalTau + tauIncrement(i) > tauLimit) then 
        deltaX       = tauLimit - totalTau
        totalTau     = totalTau     + deltaX
        totalProduct = totalProduct + deltaX * f(i) 
      else
        totalTau     = totalTau     + tauIncrement(i) 
        totalProduct = totalProduct + tauIncrement(i) * f(i) 
      end if 
      if(totalTau >= tauLimit) exit
    end do 
    weight_by_extinction = totalProduct/totalTau
  end function weight_by_extinction

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! FUNCTION cloud_top height
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  pure function cloud_top_height(nLevels,ctp,temp,pres)
     ! Inputs
     integer,intent(in) :: &
        nLevels ! Number of vertical levels
     real(wp),intent(in) :: &
        ctp     ! Cloud-top pressure   
     real(wp),intent(in),dimension(nLevels) :: &
        temp    ! Temperature profile
     real(wp),intent(in),dimension(nLevels+1) :: &
        pres    ! Pressure at layer interfaces   
     ! Outputs
     real(wp) :: cloud_top_height
     
     ! Local variables
     integer :: ij
     real(wp) :: z0
     real(wp),dimension(nLevels) :: z

     z0=0._wp
     cloud_top_height=0._wp
     do ij=nLevels,1,-1
        z(ij) = z0-(temp(ij)*rd/grav)*alog10(pres(ij)/pres(ij+1))
        z0    = z(ij)
        if (ctp .gt. pres(ij)) exit
     enddo
     if (ij .ne. nLevels) cloud_top_height = z(ij+1)+(z(ij)-z(ij+1))*(pres(ij+1)-ctp)/(pres(ij+1)-pres(ij))
     if (ij .eq. nLevels) cloud_top_height = z(ij)*(pres(ij+1)-ctp)/(pres(ij+1)-pres(ij))

  end function cloud_top_height
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE ctw_reff
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine ctw_reff(nLevels,tau,tauL,tauI,cfrac,phase,reffL,reffI,retReff)
     ! Inputs
     integer, intent(in) :: &
        nLevels,    & ! Number of vertical levels
        phase
     real(wp),intent(in),dimension(nLevels) :: &
        tau,   & ! Optical depth (total)
        tauI,  & ! Optical depth (ice)
        tauL,  & ! Optical depth (liquid) 
        cfrac, & ! Subcolumn cloud fraction
        reffL, & ! Model effective radius (liquid)
        reffI    ! Model effective radius (ice)
     ! Outputs
     real(wp),intent(out) :: &
        retReff  ! Retrieved effective radius
       
     ! Local variables
     integer :: ij,ik,il
     real(wp),dimension(nLevels) :: tau_ac_prof,frac_in_upc
     real(wp) :: num,den
     
     ! Initialize
     tau_ac_prof(1:nLevels)  = 0._wp
     frac_in_upc(1:nLevels)  = 0._wp

     ! Find where in the cloud layer where tau > CLARR_upperCloudTauLim
     do il = 2, nLevels 
        ! First find the layers that belong to the upper layers of the cloud
        if (sum(cfrac) > 0 ) then
            tau_ac_prof(il)  = tau_ac_prof(il-1)  + tau (il)
             if ( tau_ac_prof(il) < CLARA_upperCloudTauLim ) then
                ! The whole layer is in the upper part of the cloud
                frac_in_upc(il) = 1.0
             else
                if ( tau_ac_prof(il-1) < CLARA_upperCloudTauLim ) then
                   ! The boundary is found in this layer
                   frac_in_upc(il) = (1.0 - tau_ac_prof(il-1)) / &
                        (tauL(il) + tauI(il))
                endif
             endif
        else
            ! Pass the value from the above layer
            tau_ac_prof(il)  = tau_ac_prof(il-1)
        endif
     enddo  
        
     ! Fill sub grid variables if there is a cloud
     if (tau_ac_prof(nLevels) > 0.0 ) THEN
        ! Calculate COT, REFF and CWP
        num = 0._wp
        den = 0._wp
        if ( phase .eq. 1) THEN
           DO il = 1, nLevels
              num = num + frac_in_upc(il) * reffL(il) * tauL(il)
              den = den + frac_in_upc(il) * tauL(il)
           END DO
           retReff = num/den
        else ! 2 == ice
            DO il = 1, nLevels
               num = num + frac_in_upc(il) * reffI(il) * tauI(il)
               den = den + frac_in_upc(il) * tauI(il)
            END DO
            retReff = num/den
        endif     
     endif
     
  end subroutine ctw_reff

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! END MODULE mod_clara_sim
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module mod_clara_sim