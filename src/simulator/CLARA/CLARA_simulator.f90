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
  USE COSP_KINDS,  ONLY: wp
  implicit none

  ! Retrieval parameters (set up during initialization)
  real(wp) :: &
     CLARA_upperCloudTauLim,   & ! Optical depth into cloud AVHRR can see for phase retrieval
     CLARA_minOpticalThickness,& ! Lower limit of optical sensitivity for AVHRR
     CLARA_STlimit,            & ! Optical depth limit for opaque clouds
     CLARA_
  integer,parameter :: &
     CLARA_phaseIsNone = 0       ! No retrieved phase

  contains
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE clara_subcolumn
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine clara_subcolumn(nPoints,nLevels,nSubCols,tautot,tautotliq,tautotice,g,w0)
    ! Inputs
    integer,intent(in) :: &
      nPoints,      & ! Number of gridpoints
      nLevels,      & ! Number of vertical levels
      nSubCols        ! Number of subcolumns
    real(wp),dimension(nPoints,nSubCols,nLevels),intent(in) :: &
      tautot,       & ! TOA-2-SFC integrated total subcolumn optical thickness @ 0.67 microns.
      tautotliq,    & ! TOA-2-SFC integrated liquid subcolumn optical thickness @ 0.67 microns.
      tautotice,    & ! TOA-2-SFC integrated ice subcolumn optical thickness @ 0.67 microns.
      g,            & ! Subcolumn assymetry parameter @ 3.7 microns.
      w0              ! Subcolumn single-scattering albedo @ 3.7 microns.
    ! Outputs
    real(wp),dimension(nPoints,nSubcols) :: &
       clara_tau      ! Retrieved CLARA optical-depth
    integer(wp),dimension(nPoints,nSubCols) :: &
       clara_phase,   & ! Retrieved CLARA cloud phase
       clara_cldType, & ! Retrieved CLARA cloud type. (1-Stratus Continental,2-Stratus Maritime,
                        ! 3-Cumulus Continental Clean, 4-Cumulus Continental Polluted,
                        ! 5-Cumulus Maritime, 6-Cirrus)
       clara_cldFlag    ! Retrieved CLARA cloud flag. (0-Clear, 1-Subvisible, 
                        ! 2-Semi-transparent, 3-Opaque)                  
    logical,dimension(nPoints,nSubCols) :: &   
       clara_cldMask  ! Retrieved CLARA cloud-mask
       
    ! Local variables
    integer                     :: ij,ik,il
    real(wp),dimension(nLevels) :: frac,tauI,tauL
    real(wp) :: num,den
    
    ! ####################################################################################
    ! Compute total column optical depth. 
    ! The input optical-depth, tautot, is the TOA-2-SFC integrated optical depth, so the 
    ! retrieved total-column optical depth is simply the optical depth in the lowest layer.
    ! ####################################################################################
    clara_tau(1:nPoints,1:nSubCols) = tautot(1:nPoints,1:nSubCols,nLevels)

    ! ####################################################################################
    ! Compute cloud-mask.
    ! Cloudy scenes are only detectable above a certain optical depth, so use as
    ! ####################################################################################
    clara_cldMask(1:nPoints,1:nSubCols) = .false.
    where(clara_tau(1:nPoints, 1:nSubCols) .gt. CLARA_minOpticalThickness) clara_cldMask(1:nPoints,1:nSubCols) = .true.
    
    ! ##############################################################################
    ! Compute cloud thickness flags
    ! ##############################################################################          
    where(clara_tau(1:nPoints, 1:nSubCols) .lt. 0.001)                     clara_cldFlag = 0
    where(clara_tau(1:nPoints, 1:nSubCols) .lt. CLARA_minOpticalThickness) clara_cldFlag = 1
    where(clara_tau(1:nPoints, 1:nSubCols) .ge. CLARA_minOpticalThickness .and. &
          clara_tau(1:nPoints, 1:nSubCols) .le. CLARA_STlimit)             clara_cldFlag = 2
    where(clara_tau(1:nPoints, 1:nSubCols) .gt. CLARA_STlimit)             clara_cldFlag = 3
 
    do ij=1,nPoints
       do ik=1,nSubCols
          ! ##############################################################################
	      ! Compute cloud phase
          ! ##############################################################################
          frac(1:nLevels) = 0._wp
          tauL(1:nLevels) = 0._wp
          tauI(1:nLevels) = 0._wp
          num             = 0._wp
          den             = 0._wp
          if (clara_cldMask(ij,ik)) then
             do il=2,nLevels
                ! Determine fraction of each layer that contributes to optical depth use for
                ! weighting.
                ! a) Above cloud layer => include this whole layer.
                if (tautot(ij,ik,il) .lt. CLARA_upperCloudTauLim) frac(il)=1._wp
                ! b) In cloud layer => determine fraction of layer to include.
                if (tautot(ij,ik,il) .gt. CLARA_upperCloudTauLim) then
                   frac(il)=(1._wp - tautot(ij,ik,il-1)) / tautot(ij,ik,il)
                endif
                ! c) Below cloud layer => omit this layer.
                if (frac(il) .lt. 0) frac(il)=0._wp
             
                ! Weight optical depth in each layer.
                num = num + frac(il) * (tautotliq(ij,ik,il)+2._wp*tautotice(ij,ik,il))
                den = den + frac(il) * (tautotliq(ij,ik,il)+tautotice(ij,ik,il))             
             enddo
             ! Determine cloud-top phase for this sub-grid column
             clara_phase(ij,ik) = nint(num/den)
          else
             clara_phase(ij,ik) = CLARA_phaseIsNone   
          endif
          ! ##############################################################################
          ! Compute cloud type
          ! I'm not sure on how to do this at the moment....
          ! ##############################################################################
          if (clara_cldMask(ij,ik)) then
          
          endif
          
          print*,ij,ik,clara_phase(ij,ik),clara_cldFlag(ij,ik)
       enddo
    enddo


     
  
  end subroutine clara_subcolumn

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE clara_column
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine clara_column()

  end subroutine clara_column
  
  
  
  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! END MODULE mod_clara_sim
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module mod_clara_sim