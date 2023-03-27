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
! Apr 2015 - D. Swales - Modified for RTTOVv11.3
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE MOD_COSP_RTTOV_INTERFACE
  USE COSP_KINDS, ONLY: wp
!  USE MOD_COSP_CONFIG,  ONLY: RTTOV_MAX_CHANNELS,rttovDir
  USE MOD_COSP_RTTOV,   ONLY: rttov_IN
  IMPLICIT NONE
  
CONTAINS

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_rttov_init
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_RTTOV_INIT(NchanIN,nlevels,Lrttov_bt,Lrttov_rad,Lrttov_refl,                 &
                             Lrttov_cld,Lrttov_aer,Lrttov_cldparam,Lrttov_aerparam,    &
                             rttov_input_namelist)
    integer,intent(in) :: &
         NchanIN,    &
         nlevels
!    integer,intent(in),dimension(RTTOV_MAX_CHANNELS) :: &
!         channelsIN     ! RTTOV channels
    logical,intent(in)   :: &
         Lrttov_bt,        &
         Lrttov_rad,       &
         Lrttov_refl,      &
         Lrttov_cld,       &
         Lrttov_aer,       &
         Lrttov_cldparam,  &
         Lrttov_aerparam
    
    ! JKS testing using a RTTOV input namelist here 
    ! (default cosp_rttov namelist is set in cosp.F90)
    character(len=256),intent(in) :: rttov_input_namelist
     
     print*,'Running COSP_RTTOV_INIT from STUB files.', &
         'To run RTTOV, compile COSP after setting environmental variable "RTTOV"'
    
  END SUBROUTINE COSP_RTTOV_INIT
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_rttov_simulate
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_RTTOV_SIMULATE(rttovIN,lCleanup,                                 & ! Inputs
                                 bt_total,bt_clear,                                & ! Brightness Temp Outputs
                                 rad_total,rad_clear,rad_cloudy,                   & ! Radiance Outputs
                                 refl_total,refl_clear,                            & ! Reflectance Outputs
                                 error)        
  
      type(rttov_in),intent(in) :: &
          rttovIN
      logical,intent(in) :: &
          lCleanup   ! Flag to determine whether to deallocate RTTOV types          
      real(wp),dimension(rttovIN%nPoints,rttovIN%nChannels) :: & ! Can I do this? I guess so!
          bt_total,                          &        ! All-sky
          bt_clear,                          &        ! Clear-sky
          rad_total,                         &        ! All-sky
          rad_clear,                         &        ! Clear-sky
          rad_cloudy,                        &        ! Cloudy-sky
          refl_total,                        &        ! All-sky
          refl_clear                                  ! Clear-sky
      character(len=128) :: &
          error     ! Error messages (only populated if error encountered)  
  
  
      print*,'Running COSP_RTTOV_SIMULATE from STUB files.', &
             'To run RTTOV, compile COSP after setting environmental variable "RTTOV"'
  ! How do I want the interface to function? How should it to be consistent with the rest of COSP?
  
  END SUBROUTINE COSP_RTTOV_SIMULATE
  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! END MODULE
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE MOD_COSP_RTTOV_INTERFACE
