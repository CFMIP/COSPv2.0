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
MODULE MOD_COSP_RTTOV_INTERFACE
  USE COSP_KINDS,               ONLY: wp
  USE MOD_COSP_CONFIG,          ONLY: RTTOV_MAX_CHANNELS
  USE MOD_COSP_RTTOV,           ONLY: nch_in,plat_in,sat_in,sens_in,ichan_in
   
  IMPLICIT NONE

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !									TYPE rttov_in
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  type rttov_in
     integer,pointer :: &
          nprof,        & ! Number of profiles to simulate
          nlevels         ! Number of pressure levels
     real(wp),pointer :: &
          zenang,       & ! Satellite zenith angle
          co2,          & ! Carbon dioxide 
          ch4,          & ! Methane 
          n2o,          & ! n2o 
          co              ! Carbon monoxide
     real(wp),dimension(:),pointer :: &
          surfem          ! Surface emissivities for the channels
     real(wp),dimension(:),pointer :: &
          h_surf,       & ! Surface height
          u_surf,       & ! U component of surface wind
          v_surf,       & ! V component of surface wind
          t_skin,       & ! Surface skin temperature
          p_surf,       & ! Surface pressure
          t_surf,       & ! 1.5 m Temperature
          q_surf,       & ! 1.5 m Specific humidity
          lsmask,       & ! land-sea mask
          latitude        ! Latitude
     real(wp),dimension(:,:),pointer :: &
          p,            & ! Pressure profiles  
          t,            & ! Temperature profiles
          q,            & ! Humidity profiles
          o3              ! Ozone profiles
  end type rttov_in
CONTAINS

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !								 SUBROUTINE cosp_rttov_init
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_RTTOV_INIT(Nchannels,platform,satellite,instrument,channels)
    integer,intent(in) :: & 
       Nchannels, & ! Number of channels
       platform,  & ! Satellite platform
       satellite, & ! Satellite
       instrument   ! Instrument
     integer,intent(in),dimension(RTTOV_MAX_CHANNELS) :: &
       channels     ! RTTOV channels  
    
    ! Initialize
    nch_in   = Nchannels
    plat_in  = platform
    sat_in   = satellite
    sens_in  = instrument
    ichan_in = channels
    
  END SUBROUTINE COSP_RTTOV_INIT
  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! 									END MODULE
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE MOD_COSP_RTTOV_INTERFACE
