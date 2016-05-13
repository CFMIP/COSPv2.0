! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copyright (c) 2009, Centre National de la Recherche Scientifique
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
! December 2008, S. Bony,  H. Chepfer and J-L. Dufresne : 
! - optimization for vectorization
! Version 2.0 (October 2008)
! Version 2.1 (December 2008)
! May 2015 - D. Swales - Modified for COSPv2.0
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module mod_parasol
  USE COSP_KINDS,          ONLY: wp
  USE COSP_MATH_CONSTANTS, ONLY: pi
  use mod_cosp_config,     ONLY: R_UNDEF
  implicit none

  ! LUT Parameters
  INTEGER,PARAMETER :: &
       ntetas = 5, & ! Number of angles in LUT
       nbtau  = 7    ! Number of optical depths in LUT

  ! Optical depth
  REAL(WP),parameter,dimension(nbtau) :: &
       tau = (/0., 1., 5., 10., 20., 50., 100./)
  REAL(WP),parameter,dimension(ntetas) :: &
       tetas = (/0., 20., 40., 60., 80./)
  ! LUTs
  REAL(WP),parameter,dimension(ntetas,nbtau) :: &
       ! LUT for spherical liquid particles
       rlumA = reshape(source=(/ 0.03,     0.03,     0.03,     0.03,     0.03,           &
                                 0.090886, 0.072185, 0.058410, 0.052498, 0.034730,       &
                                 0.283965, 0.252596, 0.224707, 0.175844, 0.064488,       &
                                 0.480587, 0.436401, 0.367451, 0.252916, 0.081667,       &
                                 0.695235, 0.631352, 0.509180, 0.326551, 0.098215,       &
                                 0.908229, 0.823924, 0.648152, 0.398581, 0.114411,       &
                                 1.0,      0.909013, 0.709554, 0.430405, 0.121567/),     &
                                 shape=(/ntetas,nbtau/)), & 
	   ! LUT for ice particles         			     
       rlumB = reshape(source=(/ 0.03,     0.03,     0.03,     0.03,     0.03,           &
                                 0.092170, 0.087082, 0.083325, 0.084935, 0.054157,       &
                                 0.311941, 0.304293, 0.285193, 0.233450, 0.089911,       &
                                 0.511298, 0.490879, 0.430266, 0.312280, 0.107854,       &
                                 0.712079, 0.673565, 0.563747, 0.382376, 0.124127,       &
                                 0.898243, 0.842026, 0.685773, 0.446371, 0.139004,       &
                                 0.976646, 0.912966, 0.737154, 0.473317, 0.145269/),     &
                                 shape=(/ntetas,nbtau/))  
  
contains
  SUBROUTINE parasol_subcolumn(npoints,nrefl,tautot_S_liq,tautot_S_ice,refl)
    ! ##########################################################################
    ! Purpose: To compute Parasol reflectance signal from model-simulated profiles 
    !          of cloud water and cloud fraction in each sub-column of each model 
    !          gridbox.
    !
    !
    ! December 2008, S. Bony,  H. Chepfer and J-L. Dufresne : 
    ! - optimization for vectorization
    !
    ! Version 2.0 (October 2008)
    ! Version 2.1 (December 2008)
    ! ##########################################################################
    
    ! INPUTS
    INTEGER,intent(in) :: &
         npoints,              & ! Number of horizontal gridpoints
         nrefl                   ! Number of angles for which the reflectance is computed
    REAL(WP),intent(inout),dimension(npoints) :: &
         tautot_S_liq,         & ! Liquid water optical thickness, from TOA to SFC
         tautot_S_ice            ! Ice water optical thickness, from TOA to SFC
    ! OUTPUTS
    REAL(WP),intent(inout),dimension(npoints,nrefl) :: &
         refl                    ! Parasol reflectances
    
    ! LOCAL VARIABLES
    REAL(WP),dimension(npoints) :: &
         tautot_S,             & ! Cloud optical thickness, from TOA to surface
         frac_taucol_liq,      & !
         frac_taucol_ice         !
    
    ! Look up table variables:
    INTEGER                            :: ny,it 
    REAL(WP),dimension(ntetas)         :: r_norm
    REAL(WP),dimension(ntetas,nbtau-1) :: aa,ab,ba,bb
    REAL(WP),dimension(npoints,5)      :: rlumA_mod,rlumB_mod
    
    !--------------------------------------------------------------------------------
    ! Lum_norm=f(tetaS,tau_cloud) derived from adding-doubling calculations
    !        valid ONLY ABOVE OCEAN (albedo_sfce=5%)
    !        valid only in one viewing direction (theta_v=30�, phi_s-phi_v=320�)
    !        based on adding-doubling radiative transfer computation
    !        for tau values (0 to 100) and for tetas values (0 to 80)
    !        for 2 scattering phase functions: liquid spherical, ice non spherical
    
    ! Initialize
    rlumA_mod(1:npoints,1:5) = 0._wp
    rlumB_mod(1:npoints,1:5) = 0._wp

    r_norm(1:ntetas)=1._wp/ cos(pi/180._wp*tetas(1:ntetas))
    
    tautot_S_liq(1:npoints) = max(tautot_S_liq(1:npoints),tau(1))
    tautot_S_ice(1:npoints) = max(tautot_S_ice(1:npoints),tau(1))
    tautot_S(1:npoints)     = tautot_S_ice(1:npoints) + tautot_S_liq(1:npoints)

    ! Relative fraction of the opt. thick due to liquid or ice clouds
    WHERE (tautot_S(1:npoints) .gt. 0.)
       frac_taucol_liq(1:npoints) = tautot_S_liq(1:npoints) / tautot_S(1:npoints)
       frac_taucol_ice(1:npoints) = tautot_S_ice(1:npoints) / tautot_S(1:npoints)
    ELSEWHERE
       frac_taucol_liq(1:npoints) = 1._wp
       frac_taucol_ice(1:npoints) = 0._wp
    END WHERE
    tautot_S(1:npoints)=MIN(tautot_S(1:npoints),tau(nbtau))
    
    ! Linear interpolation    
    DO ny=1,nbtau-1
       ! Microphysics A (liquid clouds) 
       aA(1:ntetas,ny) = (rlumA(1:ntetas,ny+1)-rlumA(1:ntetas,ny))/(tau(ny+1)-tau(ny))
       bA(1:ntetas,ny) = rlumA(1:ntetas,ny) - aA(1:ntetas,ny)*tau(ny)
       ! Microphysics B (ice clouds)
       aB(1:ntetas,ny) = (rlumB(1:ntetas,ny+1)-rlumB(1:ntetas,ny))/(tau(ny+1)-tau(ny))
       bB(1:ntetas,ny) = rlumB(1:ntetas,ny) - aB(1:ntetas,ny)*tau(ny)
    ENDDO
    
    DO it=1,ntetas
       DO ny=1,nbtau-1
          WHERE (tautot_S(1:npoints) .ge. tau(ny).and. &
                 tautot_S(1:npoints) .le. tau(ny+1))
             rlumA_mod(1:npoints,it) = aA(it,ny)*tautot_S(1:npoints) + bA(it,ny)
             rlumB_mod(1:npoints,it) = aB(it,ny)*tautot_S(1:npoints) + bB(it,ny)
          END WHERE
       END DO
    END DO
    
    DO it=1,ntetas
       refl(1:npoints,it) = frac_taucol_liq(1:npoints) * rlumA_mod(1:npoints,it) &
            + frac_taucol_ice(1:npoints) * rlumB_mod(1:npoints,it)
       ! Normalized radiance -> reflectance: 
       refl(1:npoints,it) = refl(1:npoints,it) * r_norm(it)
    ENDDO
    
    RETURN
  END SUBROUTINE parasol_subcolumn
  ! ######################################################################################
  ! SUBROUTINE parasol_gridbox
  ! ######################################################################################
  subroutine parasol_column(npoints,nrefl,ncol,land,refl,parasolrefl)

    ! Inputs
    integer,intent(in) :: &
         npoints, & ! Number of horizontal grid points
         ncol,    & ! Number of subcolumns
         nrefl      ! Number of solar zenith angles for parasol reflectances
    real(wp),intent(in),dimension(npoints) :: &
         land       ! Landmask [0 - Ocean, 1 - Land]
    real(wp),intent(in),dimension(npoints,ncol,nrefl) :: &
         refl       ! Subgrid parasol reflectance ! parasol

    ! Outputs
    real(wp),intent(out),dimension(npoints,nrefl) :: &
         parasolrefl   ! Grid-averaged parasol reflectance

    ! Local variables
    integer :: k,ic

    ! Compute grid-box averaged Parasol reflectances
    parasolrefl(:,:) = 0._wp
    do k = 1, nrefl
       do ic = 1, ncol
          parasolrefl(:,k) = parasolrefl(:,k) + refl(:,ic,k)
       enddo
    enddo
    
    do k = 1, nrefl
       parasolrefl(:,k) = parasolrefl(:,k) / float(ncol)
       ! if land=1 -> parasolrefl=R_UNDEF
       ! if land=0 -> parasolrefl=parasolrefl
       parasolrefl(:,k) = parasolrefl(:,k) * MAX(1._wp-land(:),0.0) &
            + (1._wp - MAX(1._wp-land(:),0.0))*R_UNDEF
    enddo
  end subroutine parasol_column

end module mod_parasol
