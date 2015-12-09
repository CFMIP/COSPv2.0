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
MODULE MOD_COSP_CLARA_INTERFACE
  USE COSP_KINDS,            ONLY: wp
  USE MOD_CLARA_SIM,         ONLY: CLARA_upperCloudTauLim,CLARA_minOpticalThickness,     &
                                   CLARA_STlimit,CLARA_RTTOVclr,CLARA_Tb_subvis,         &
                                   CLARA_Tb_semitrans,CLARA_Tb_opaque,RTTOV_satZen,      &
                                   RTTOV_satLat,nChannels,CLARA_re_water_min,            &
                                   CLARA_re_water_max,CLARA_re_ice_min,CLARA_re_ice_max, &
                                   trial_re_w,trial_re_i,g_w,g_i,w0_i,w0_w,num_trial_res,&
                                   get_g_nir,get_ssa_nir,CLARA_phaseThresh
  USE clara_rttov_interface, ONLY: jplm,jpim,get_rttov_coeffs,R
  
  implicit none
    
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Type clara_IN
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  type clara_IN
     integer,pointer ::      &
          Npoints,           & ! Number of horizontal gridpoints
          Ncolumns,          & ! Number of subcolumns
          Nlevels              ! Number of vertical levels  
     real(wp),pointer ::     &
          emsfc_lw             ! 10.5 micron surface emissivity.
     real(wp),pointer ::     &
          skt(:),            & ! Skin temperature (nPoints)
          at2m(:),           & ! Model temperature at 2meters.
          p2m(:),            & ! Model pressure at 2meters.
          lsmask(:),         & ! Model land/sea mask
          latitude(:)          ! Latitude
     real(wp),pointer ::     &
          pfull(:,:),        & ! Model pressure at full levels (Npoints,Nlevels)
          phalf(:,:),        & ! Model pressure at half levels (Npoints,Nlevels+1)
          at(:,:),           & ! Model temperature at full levels
          qv(:,:)              ! Model specific humidity at full levels
     real(wp),pointer ::     &
          tautot(:,:,:),     & ! TOA-2-SFC integrated subcolumn total optical thickness @ 0.67 microns.
          tautotliq(:,:,:),  & ! TOA-2-SFC integrated subcolumn liquid optical thickness @ 0.67 microns.
          tautotice(:,:,:),  & ! TOA-2-SFC integrated subcolumn ice optical thickness @ 0.67 microns.
          g(:,:,:),          & ! Subcolumn assymetry parameter  
          w0(:,:,:),         & ! Subcolumn single-scattering albedo
          liqFrac(:,:,:)       ! Fractional contribution to optical depth from liquid water
     real(wp),pointer ::      &    
          cldLiqRTTOV(:,:,:,:), & ! RTTOV cloud liquid concentration.
          cldIceRTTOV(:,:,:,:)    ! RTTOV cloud ice concentration.
  end type clara_IN
  
  contains
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_clara_init
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE cosp_clara_init(CLARA_RTTOVclrIN,CLARA_Tb_subvisIN,CLARA_Tb_semitransIN,    &
                             CLARA_Tb_opaqueIN)
     ! Inputs
     integer,intent(in) :: &
        CLARA_Tb_subvisIN,CLARA_Tb_semitransIN,CLARA_Tb_opaqueIN   
     logical,intent(in) :: &
        CLARA_RTTOVclrIN    
     ! Local variables
     integer :: i
     
     CLARA_upperCloudTauLim     = 1._wp  ! How many optical depths does AVHRR see into the cloud?
     CLARA_minOpticalThickness  = 0.2_wp ! Lower limit of optical sensitivity of AVHRR
     CLARA_STlimit              = 10._wp ! Optical depth limit for opaque cloud @ 3.7microns
     CLARA_RTTOVclr             = CLARA_RTTOVclrIN! Use RTTOV for clear-sky brightness temperature
     CLARA_Tb_subvis            = CLARA_Tb_subvisIN     ! RTTOV clear-sky
     CLARA_Tb_semitrans         = CLARA_Tb_semitransIN      ! RTTOV w/ scattering
     CLARA_Tb_opaque            = CLARA_Tb_opaqueIN      ! RTTOV w/o scattering (black-body)
     CLARA_re_water_min         = 3._wp  ! Minimum effective radius (liquid)
     CLARA_re_water_max         = 34._wp ! Maximum effective radius (liquid)
     CLARA_re_ice_min           = 5._wp  ! Minimum effective radius (ice)
     CLARA_re_ice_max           = 80._wp ! Minimum effective radius (ice)
     CLARA_phaseThresh          = 0.7_wp ! Fraction of phase needed for retrieval (i.e. 
                                         ! 70% of the cloud must be one phase or the other)
    ! Precompute near-IR optical params vs size for retrieval scheme    
    trial_re_w(1:num_trial_res) = CLARA_re_water_min + (CLARA_re_water_max -             &
                                  CLARA_re_water_min) /                                  &
                                  (num_trial_res-1) * (/(i, i=0, num_trial_res-1)/)
    trial_re_i(1:num_trial_res) = CLARA_re_ice_min   + (CLARA_re_ice_max -               &
                                  CLARA_re_ice_min) /                                    &
                                  (num_trial_res-1) * (/(i, i=0, num_trial_res-1)/)
    
    ! Initialize estimates for size retrieval
    g_w(1:num_trial_res)  = get_g_nir(trial_re_w(1:num_trial_res),1)
    w0_w(1:num_trial_res) = get_ssa_nir(trial_re_w(1:num_trial_res),1)
    g_i(1:num_trial_res)  = get_g_nir(trial_re_i(1:num_trial_res),2)
    w0_i(1:num_trial_res) = get_ssa_nir(trial_re_i(1:num_trial_res),2)        

  END SUBROUTINE cosp_clara_init
  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_clara_rttov_init
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine cosp_clara_rttov_init(rttovDir,addrefrac,use_q2m,clw_data,addsolar,addclouds,&
                                   addaerosl,user_cld_opt_param,ozone_data,co2_data,     &
                                   n2o_data,ch4_data,co_data,addinterp,calcemis,calcrefl,&
                                   id_platform,id_satellite,sensor,nchan,channels,dbg)
    ! Inputs
    character(len=*),intent(in) :: &
       rttovDir              ! RTTOV source code root directory
    logical(kind=jplm),intent(in) :: &
       addrefrac,use_q2m,clw_data,addsolar,addclouds,addaerosl,user_cld_opt_param,       & 
       ozone_data,co2_data,n2o_data,ch4_data,co_data,addinterp,calcemis,calcrefl
    integer(kind=jpim),intent(in) ::  &
       id_platform,id_satellite,sensor,nchan
    integer(kind=jpim),dimension(nchan),intent(in) :: &
       channels
    integer,intent(in) :: &
       dbg                   ! Debug flag, 0 = no debug; -1 = use earlier

    ! Local variables
    integer            :: ii
    integer(kind=jpim) :: error
    
    ! Parameters
    ! NOAA-14 satellite zenith angle
    integer,parameter :: nlatNOAA14 = 160
    real(wp),target,dimension(nlatNOAA14) :: & 
       satzenN14 = (/48.6953, 52.16074, 49.21415, 43.10466, 35.28501, 26.86696,          &
                    18.07936,  9.29647,  4.27182,  5.44160,  7.06823,  8.47386,  9.73470,& 
                    10.92017, 12.04863, 13.14443, 14.16183, 15.13755, 16.08984, 17.02158,& 
                    17.93092, 18.7867,  19.62353, 20.43138, 21.22349, 21.99617, 22.73001,& 
                    23.4451,  24.13814, 24.8405,  25.49737, 26.12563, 26.73702, 27.34391,& 
                    27.93138, 28.50525, 29.04329, 29.58218, 30.09974, 30.6055,  31.08346,& 
                    31.54523, 31.99565, 32.43642, 32.83337, 33.25727, 33.64569, 34.02335,& 
                    34.38442, 34.73474, 35.0783,  35.38177, 35.68592, 35.95063, 36.25259,& 
                    36.5273,  36.77239, 37.02624, 37.24519, 37.4687,  37.65719, 37.83723,& 
                    38.01817, 38.17473, 38.31023, 38.44192, 38.54707, 38.67994, 38.78153,& 
                    38.8672,  38.92857, 39.02708, 39.0919,  39.14176, 39.17746, 39.22013,& 
                    39.26681, 39.29843, 39.26291, 39.15234, 39.17629, 39.2086,  39.17621,& 
                    39.20073, 39.16491, 39.16323, 39.16189, 39.10093, 39.06634, 39.02916,& 
                    38.96496, 38.89199, 38.50542, 37.45097, 37.38972, 37.46801, 37.49936,& 
                    37.47264, 37.5317,  37.50639, 37.43009, 37.33318, 37.20377, 36.96389,& 
                    36.7553,  36.46271, 36.23306, 35.9711,  35.68075, 35.34879, 35.04118,& 
                    34.72102, 34.36323, 34.01497, 33.64992, 33.24688, 32.88052, 32.4819, &
                    32.06557, 31.61389, 31.13428, 30.64936, 30.13377, 29.61869, 29.08395,& 
                    28.51819, 27.93199, 27.33422, 26.71016, 26.12987, 25.53804, 24.98637,& 
                    24.39119, 23.83387, 23.26536, 22.6984,  22.04817, 21.39887, 20.65174,& 
                    19.79424, 18.83796, 17.77637, 16.7094,  15.60699, 14.53348, 13.47033,& 
                    12.30018, 11.09956,  9.81493, 8.27127,   6.12766,  3.5408,   6.17367,& 
                    14.8725,  23.13026, 30.81659, 37.50245, 43.33991, 48.55217, 53.36691/)
   real(wp),target,dimension(nlatNOAA14) :: &
      latN14 = (/-89.14152, -88.02943, -86.91077, -85.79063, -84.66992, -83.54895,       &
                 -82.42782, -81.30659, -80.18531, -79.06398, -77.94263, -76.82124,       &
                 -75.69984, -74.57843, -73.45701, -72.33558, -71.21413, -70.09269,       &
                 -68.97124, -67.84978, -66.72832, -65.60686, -64.4854,  -63.36393,       &
                 -62.24246, -61.12099, -59.99952, -58.87804, -57.75657, -56.63509,       &
                 -55.51361, -54.39214, -53.27066, -52.14917, -51.02769, -49.90621,       & 
                 -48.78473, -47.66325, -46.54176, -45.42028, -44.29879, -43.17731,       &
                 -42.05582, -40.93434, -39.81285, -38.69136, -37.56988, -36.44839,       &
                 -35.3269,  -34.20542, -33.08393, -31.96244, -30.84096, -29.71947,       &
                 -28.59798, -27.47649, -26.355,   -25.23351, -24.11202, -22.99054,       &
                 -21.86905, -20.74756, -19.62607, -18.50458, -17.38309, -16.2616,        &
                 -15.14011, -14.01862, -12.89713, -11.77564, -10.65415,  -9.532663,      & 
                  -8.411174, -7.289684, -6.168194, -5.046704, -3.925215, -2.803725,      &
                  -1.682235, -0.5607449, 0.5607449, 1.682235,  2.803725,  3.925215,      &
                   5.046704,  6.168194,  7.289684,  8.411174,  9.532663, 10.65415,       &
                  11.77564,  12.89713,  14.01862,  15.14011,  16.2616,   17.38309,       &
                  18.50458,  19.62607,  20.74756,  21.86905,  22.99054,  24.11202,       &
                  25.23351,  26.355,    27.47649,  28.59798,  29.71947,  30.84096,       &
                  31.96244,  33.08393,  34.20542,  35.3269,   36.44839,  37.56988,       &
                  38.69136,  39.81285,  40.93434,  42.05582,  43.17731,  44.29879,       &
                  45.42028,  46.54176,  47.66325,  48.78473,  49.90621,  51.02769,       &
                  52.14917,  53.27066,  54.39214,  55.51361,  56.63509,  57.75657,       &
                  58.87804,  59.99952,  61.12099,  62.24246,  63.36393,  64.4854,        &
                  65.60686,  66.72832,  67.84978,  68.97124,  70.09269,  71.21413,       &
                  72.33558,  73.45701,  74.57843,  75.69984,  76.82124,  77.94263,       &
                  79.06398,  80.18531,  81.30659,  82.42782,  83.54895,  84.66992,       &
                  85.79063,  86.91077,  88.02943,  89.14152/)
    
    ! Number of channels used by RTTOV in CLARA. This is used directly by the simulator.  
    Nchannels = nchan
      
    ! Allocations
    allocate(R%id%channels(nchan),R%other_opts%calcemis(nchan),&
             R%other_opts%calcrefl(nchan))    

    ! Populate RTTOV options structure
    R%opts%config%apply_reg_limits     = .true.
    R%opts%rt_all%addrefrac            = addrefrac  
    R%opts%rt_all%use_q2m              = use_q2m   
    R%opts%rt_mw%clw_data              = clw_data
    R%opts%rt_ir%addsolar              = addsolar
    ! Set this to true so that the coefficient tables are read in. 
    R%opts%rt_ir%addclouds             = .true.!addclouds  
    R%opts%rt_ir%addaerosl             = addaerosl 
    R%opts%rt_ir%user_cld_opt_param    = user_cld_opt_param
    R%opts%rt_ir%ozone_data            = ozone_data
    R%opts%rt_ir%co2_data              = co2_data  
    R%opts%rt_ir%n2o_data              = n2o_data  
    R%opts%rt_ir%ch4_data              = ch4_data  
    R%opts%rt_ir%co_data               = co_data   
    R%opts%interpolation%addinterp     = addinterp         
    R%id%rttov_dir                     = rttovDir
    R%id%platform                      = id_platform 
    R%id%satellite                     = id_satellite
    R%id%sensor                        = sensor     
    R%id%nchannels                     = nchannels
    R%id%channels(1:nchan)             = channels         
    R%other_opts%calcemis(1:nchan)     = calcemis
    R%other_opts%calcrefl(1:nchan)     = calcrefl
    R%other_opts%dbg                   = dbg

    ! Get RTTOV coefficients 
    call get_rttov_coeffs (error)
    R%opts%rt_ir%addclouds             = addclouds  

    ! Get satellite zenith angle for satellite beging used for CLARA.
    if (id_satellite .eq. 14) then 
       RTTOV_satLat => latN14
       RTTOV_satZen => satzenN14
    endif
    if (id_satellite .eq. 15) then 
       RTTOV_satLat => latN14
       RTTOV_satZen => satzenN14
    endif
    if (id_satellite .eq. 16) then 
       RTTOV_satLat => latN14
       RTTOV_satZen => satzenN14
    endif    
     if (id_satellite .eq. 17) then 
       RTTOV_satLat => latN14
       RTTOV_satZen => satzenN14
    endif           
    if (id_satellite .eq. 18) then 
       RTTOV_satLat => latN14
       RTTOV_satZen => satzenN14
    endif
        
  end subroutine cosp_clara_rttov_init

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! END MODULE cosp_clara_interface
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end module MOD_COSP_CLARA_INTERFACE