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
