! ###########################################################################
!
! Copyright 2014, Regents of the University of Colorado. All right reserved.
! Use and duplication is permitted under the terms of the
! BSD 3-Clause License: https://opensource.org/licenses/BSD-3-Clause
!
! ###########################################################################

MODULE MOD_COSP_RTTOV
  USE COSP_KINDS,      ONLY: wp
  USE MOD_COSP_CONFIG, ONLY: RTTOV_MAX_CHANNELS
  IMPLICIT NONE
  
  ! Fields set during initialization
  integer :: &
     nch_in,  & ! Number of RTTOV channels
     plat_in, & ! RTTOV platform
     sat_in,  & ! RTTOV satellite
     sens_in    ! RTTOV instrument
  integer,dimension(RTTOV_MAX_CHANNELS) :: &
     ichan_in   ! RTTOV channel indices
   
CONTAINS
  SUBROUTINE RTTOV_subcolumn(surfem_in, prf_num_in, nlevels_in, zenang_in, p_in,t_in, q_in,  &
                         o3_in, co2_in, ch4_in, n2o_in, co_in, h_surf, u_surf, v_surf,   &
                         t_skin, p_surf, t_surf, q_surf, lsmask, latitude, tbs)   
    ! INPUTS
    integer,intent(in) :: &
         prf_num_in,      & ! Number of profiles to simulate
         nlevels_in         ! Number of pressure levels
    real(wp),intent(in) :: &
         zenang_in,       & ! Satellite zenith angle
         co2_in,          & ! Carbon dioxide 
         ch4_in,          & ! Methane 
         n2o_in,          & ! n2o 
         co_in              ! Carbon monoxide
    real(wp),intent(in),dimension(nch_in) :: &
         surfem_in          ! Surface emissivities for the channels
    real(wp),intent(in),dimension(prf_num_in) :: &
         h_surf,          & ! Surface height
         u_surf,          & ! U component of surface wind
         v_surf,          & ! V component of surface wind
         t_skin,          & ! Surface skin temperature
         p_surf,          & ! Surface pressure
         t_surf,          & ! 1.5 m Temperature
         q_surf,          & ! 1.5 m Specific humidity
         lsmask,          & ! land-sea mask
         latitude           ! Latitude
    real(wp),intent(in),dimension(prf_num_in,nlevels_in) :: &
         p_in,            & ! Pressure profiles  
         t_in,            & ! Temperature profiles
         q_in,            & ! Humidity profiles
         o3_in              ! Ozone profiles
    
    ! OUTPUTS
    real(wp),intent(inout),dimension(prf_num_in,nch_in) :: &
         tbs                ! Tbs (in the right format)
  end subroutine RTTOV_subcolumn

END MODULE MOD_COSP_RTTOV
