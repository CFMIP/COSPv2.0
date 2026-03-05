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
! Jul 2007 - A. Bodas-Salcedo - Initial version
! May 2015 - Dustin Swales    - Modified for COSPv2.0
! 
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE MOD_COSP_UTILS
  USE COSP_KINDS, ONLY: wp
  USE MOD_COSP_CONFIG
  USE mod_quickbeam_optics, only: size_distribution
  use mod_cosp_error, only: errorMessage
  IMPLICIT NONE

CONTAINS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------- SUBROUTINE COSP_PRECIP_MXRATIO --------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! SUBROUTINE COSP_PRECIP_MXRATIO(Npoints,Nlevels,Ncolumns,p,T,prec_frac,prec_type, &
!                           n_ax,n_bx,alpha_x,c_x,d_x,g_x,a_x,b_x,gamma1,gamma2,gamma3,gamma4, &
!                           flux,mxratio,reff)

!     ! Input arguments, (IN)
!     integer,intent(in) :: Npoints,Nlevels,Ncolumns
!     real(wp),intent(in),dimension(Npoints,Nlevels) :: p,T,flux
!     real(wp),intent(in),dimension(Npoints,Ncolumns,Nlevels) :: prec_frac
!     real(wp),intent(in) :: n_ax,n_bx,alpha_x,c_x,d_x,g_x,a_x,b_x,gamma1,gamma2,gamma3,gamma4,prec_type
!     ! Input arguments, (OUT)
!     real(wp),intent(out),dimension(Npoints,Ncolumns,Nlevels) :: mxratio
!     real(wp),intent(inout),dimension(Npoints,Ncolumns,Nlevels) :: reff
!     ! Local variables
!     integer :: i,j,k
!     real(wp) :: sigma,one_over_xip1,xi,rho0,rho,lambda_x,gamma_4_3_2,delta
    
!     mxratio = 0.0

!     if (n_ax >= 0.0) then ! N_ax is used to control which hydrometeors need to be computed
!         xi      = d_x/(alpha_x + b_x - n_bx + 1._wp)
!         rho0    = 1.29_wp
!         sigma   = (gamma2/(gamma1*c_x))*(n_ax*a_x*gamma2)**xi
!         one_over_xip1 = 1._wp/(xi + 1._wp)
!         gamma_4_3_2 = 0.5_wp*gamma4/gamma3
!         delta = (alpha_x + b_x + d_x - n_bx + 1._wp)
!         do k=1,Nlevels
!             do j=1,Ncolumns
!                 do i=1,Npoints
!                     if ((prec_frac(i,j,k)==prec_type).or.(prec_frac(i,j,k)==3.)) then
!                         rho = p(i,k)/(287.05_wp*T(i,k))
!                         mxratio(i,j,k)=(flux(i,k)*((rho/rho0)**g_x)*sigma)**one_over_xip1
!                         mxratio(i,j,k)=mxratio(i,j,k)/rho
!                         ! Compute effective radius
!                         if ((reff(i,j,k) <= 0._wp).and.(flux(i,k) /= 0._wp)) then
!                            lambda_x = (a_x*c_x*((rho0/rho)**g_x)*n_ax*gamma1/flux(i,k))**(1._wp/delta)
!                            reff(i,j,k) = gamma_4_3_2/lambda_x
!                         endif
!                     endif
!                 enddo
!             enddo
!         enddo
!     endif
! END SUBROUTINE COSP_PRECIP_MXRATIO

  SUBROUTINE COSP_PRECIP_MXRATIO(Npoints,Nlevels,Ncolumns, &
                                 p,T,prec_frac,prec_type, index, sd, &
                                 flux, mxratio, reff)
    ! Input arguments, (IN)
    integer,intent(in) :: Npoints,Nlevels,Ncolumns,index
    real(wp),intent(in),dimension(Npoints,Nlevels) :: p,T,flux
    real(wp),intent(in),dimension(Npoints,Ncolumns,Nlevels) :: prec_frac
    real(wp),intent(in) :: prec_type
    type(size_distribution),intent(in) :: sd
    ! Input arguments, (OUT)
    real(wp),intent(out),dimension(Npoints,Ncolumns,Nlevels) :: mxratio
    real(wp),intent(inout),dimension(Npoints,Ncolumns,Nlevels) :: reff
    ! Local variables
    real(wp) :: pi = 3.14159265358979323846264338327950288419717_wp
    integer  :: i,j,k
    real(wp) :: sigma,one_over_xip1,xi,rho0,rho,lambda_x,gamma_4_3_2,delta
    real(wp) :: apm,bpm,mu,D0,Dm  ! modified gamma DSD
    real(wp) :: lambda,N0         ! exponential DSD
    real(wp) :: avf,bvf,vscs_fct  ! fall velocity params
    
    !!! single-moment bulk microphysics is assumed !!!

    if (sd%apm(index) > 0) then
       apm = sd%apm(index)
       bpm = sd%bpm(index)
    else if (sd%rho(index) > 0) then
       apm = sd%rho(index)*pi/6._wp
       bpm = 3._wp
    else
       call errorMessage('!!! unexpected sd, at COSP_PRECIP_MXRATIO !!!')
    end if

    ! fall velocity params
    select case(sd%ftype(index))
    case(1)
       ! vf = a * D**b
       avf = sd%f1(index)
       bvf = sd%f2(index)
    case(2)
       ! PL08 formulation
       call errorMessage('!!! not implemented yet, PL08 formulation, at COSP_PRECIP_MXRATIO')
    end select

    ! size distribution params
    select case(sd%dtype(index))
    case(1)
       Dm = sd%p2(index)
       mu = sd%p3(index)
       if (Dm < 0 .or. mu < 0) then
          call errorMessage('!!! mod.gamma DSD requires mean D and mu, at COSP_PRECIP_MXRATIO')
       end if
    case(2)
       N0 = sd%p1(index)
       if (N0 < 0) then
          call errorMessage('!!! exp. DSD requires N0, at COSP_PRECIP_MXRATIO')
       end if
    case default
       call errorMessage('!!! not implemented yet, sd%dtype=3,4,5, at COSP_PRECIP_MXRATIO')
    end select

    do k=1,Nlevels
       do j=1,Ncolumns
          do i=1,Npoints
             if ((prec_frac(i,j,k)/=prec_type).and.(prec_frac(i,j,k)/=3.)) cycle
             rho0 = 1.29_wp
             rho  = p(i,k)/(287.05_wp*T(i,k))
             vscs_fct = merge(sqrt(rho0/rho),1._wp,sd%fvscs(index)==1)

             select case(sd%dtype(index))
             case(1) ! modified gamma
                if (flux(i,k) > 1.D-20) then
                   D0 = Dm*gamma(mu)/gamma(mu+1._wp) * 1.D-6  ! um -> m
                   mxratio(i,j,k) = flux(i,k)/vscs_fct/avf/D0**bvf * gamma(mu+bpm)/gamma(mu+bvf+bpm) / rho
                   if (reff(i,j,k) < 1.D-20) then
                      reff(i,j,k) = D0/2._wp*gamma(mu+3._wp)/gamma(mu+2._wp)
                   end if
                else
                   mxratio(i,j,k) = 0._wp ; reff(i,j,k) = 0._wp
                end if

             case(2) ! exponential
                if (flux(i,k) > 1.D-20) then
                   lambda = (vscs_fct*apm*avf*N0*gamma(1._wp+bpm+bvf)/flux(i,k))**(1._wp/(1._wp+bpm+bvf))
                   mxratio(i,j,k) = flux(i,k)/vscs_fct/avf*lambda**bvf*gamma(1._wp+bpm)/gamma(1._wp+bpm+bvf)
                   if (reff(i,j,k) < 1.D-20) then
                      reff(i,j,k) = 0.5_wp/lambda*gamma(4._wp)/gamma(3._wp)
                   end if
                else
                   mxratio(i,j,k) = 0._wp ; reff(i,j,k) = 0._wp
                end if

             case default
                call errorMessage('!!! not implemented yet, sd%dtype=3,4,5, at COSP_PRECIP_MXRATIO')
             end select
             
          end do
       end do
    end do
    
  END SUBROUTINE COSP_PRECIP_MXRATIO

END MODULE MOD_COSP_UTILS
