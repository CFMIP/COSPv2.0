! (c) British Crown Copyright 2008, the Met Office.
! All rights reserved.
! $Revision: 23 $, $Date: 2011-03-31 07:41:37 -0600 (Thu, 31 Mar 2011) $
! $URL: https://cfmip-obs-sim.googlecode.com/svn/devel/branches/dustinswales/cosp_utils.F90 $
! 
! Redistribution and use in source and binary forms, with or without modification, are permitted 
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

!
! History:
! Jul 2007 - A. Bodas-Salcedo - Initial version
!

MODULE MOD_COSP_UTILS
  USE COSP_KINDS, ONLY: wp
  USE MOD_COSP_CONFIG
  IMPLICIT NONE

  INTERFACE COSP_CHECK_INPUT
    MODULE PROCEDURE COSP_CHECK_INPUT_1D,COSP_CHECK_INPUT_2D,COSP_CHECK_INPUT_3D
  END INTERFACE
CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------- SUBROUTINE COSP_PRECIP_MXRATIO --------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_PRECIP_MXRATIO(Npoints,Nlevels,Ncolumns,p,T,prec_frac,prec_type, &
                          n_ax,n_bx,alpha_x,c_x,d_x,g_x,a_x,b_x,gamma1,gamma2,gamma3,gamma4, &
                          flux,mxratio,reff)

    ! Input arguments, (IN)
    integer,intent(in) :: Npoints,Nlevels,Ncolumns
    real(wp),intent(in),dimension(Npoints,Nlevels) :: p,T,flux
    real(wp),intent(in),dimension(Npoints,Ncolumns,Nlevels) :: prec_frac
    real(wp),intent(in) :: n_ax,n_bx,alpha_x,c_x,d_x,g_x,a_x,b_x,gamma1,gamma2,gamma3,gamma4,prec_type
    ! Input arguments, (OUT)
    real(wp),intent(out),dimension(Npoints,Ncolumns,Nlevels) :: mxratio
    real(wp),intent(inout),dimension(Npoints,Ncolumns,Nlevels) :: reff
    ! Local variables
    integer :: i,j,k
    real(wp) :: sigma,one_over_xip1,xi,rho0,rho,lambda_x,gamma_4_3_2,delta
    
    mxratio = 0.0

    if (n_ax >= 0.0) then ! N_ax is used to control which hydrometeors need to be computed
        xi      = d_x/(alpha_x + b_x - n_bx + 1.0)
        rho0    = 1.29
        sigma   = (gamma2/(gamma1*c_x))*(n_ax*a_x*gamma2)**xi
        one_over_xip1 = 1.0/(xi + 1.0)
        gamma_4_3_2 = 0.5*gamma4/gamma3
        delta = (alpha_x + b_x + d_x - n_bx + 1.0)
        
        do k=1,Nlevels
            do j=1,Ncolumns
                do i=1,Npoints
                    if ((prec_frac(i,j,k)==prec_type).or.(prec_frac(i,j,k)==3.)) then
                        rho = p(i,k)/(287.05*T(i,k))
                        mxratio(i,j,k)=(flux(i,k)*((rho/rho0)**g_x)*sigma)**one_over_xip1
                        mxratio(i,j,k)=mxratio(i,j,k)/rho
                        ! Compute effective radius
                        if ((reff(i,j,k) <= 0.0).and.(flux(i,k) /= 0.0)) then
                           lambda_x = (a_x*c_x*((rho0/rho)**g_x)*n_ax*gamma1/flux(i,k))**(1./delta)
                           reff(i,j,k) = gamma_4_3_2/lambda_x
                        endif
                    endif
                enddo
            enddo
        enddo
    endif
END SUBROUTINE COSP_PRECIP_MXRATIO

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------- SUBROUTINES COSP_CHECK_INPUT_1D ---------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_CHECK_INPUT_1D(vname,x,min_val,max_val)
    character(len=*) :: vname
    real(wp),intent(inout) :: x(:)
    real(wp),intent(in),optional :: min_val,max_val
    logical :: l_min,l_max
    character(len=128) :: pro_name='COSP_CHECK_INPUT_1D'
    
    l_min=.false.
    l_max=.false.
    
    if (present(min_val)) then
!       if (x < min_val) x = min_val
      if (any(x < min_val)) then 
      l_min = .true.
        where (x < min_val)
          x = min_val
        end where
      endif
    endif    
    if (present(max_val)) then
!       if (x > max_val) x = max_val
      if (any(x > max_val)) then 
        l_max = .true.
        where (x > max_val)
          x = max_val
        end where  
      endif    
    endif    
    
    if (l_min) print *,'----- WARNING: '//trim(pro_name)//': minimum value of '//trim(vname)//' set to: ',min_val
    if (l_max) print *,'----- WARNING: '//trim(pro_name)//': maximum value of '//trim(vname)//' set to: ',max_val
  END SUBROUTINE COSP_CHECK_INPUT_1D
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------- SUBROUTINES COSP_CHECK_INPUT_2D ---------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_CHECK_INPUT_2D(vname,x,min_val,max_val)
    character(len=*) :: vname
    real(wp),intent(inout) :: x(:,:)
    real(wp),intent(in),optional :: min_val,max_val
    logical :: l_min,l_max
    character(len=128) :: pro_name='COSP_CHECK_INPUT_2D'
    
    l_min=.false.
    l_max=.false.
    
    if (present(min_val)) then
!       if (x < min_val) x = min_val
      if (any(x < min_val)) then 
      l_min = .true.
        where (x < min_val)
          x = min_val
        end where
      endif
    endif    
    if (present(max_val)) then
!       if (x > max_val) x = max_val
      if (any(x > max_val)) then 
        l_max = .true.
        where (x > max_val)
          x = max_val
        end where  
      endif    
    endif    
    
    if (l_min) print *,'----- WARNING: '//trim(pro_name)//': minimum value of '//trim(vname)//' set to: ',min_val
    if (l_max) print *,'----- WARNING: '//trim(pro_name)//': maximum value of '//trim(vname)//' set to: ',max_val
  END SUBROUTINE COSP_CHECK_INPUT_2D
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------- SUBROUTINES COSP_CHECK_INPUT_3D ---------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_CHECK_INPUT_3D(vname,x,min_val,max_val)
    character(len=*) :: vname
    real(wp),intent(inout) :: x(:,:,:)
    real(wp),intent(in),optional :: min_val,max_val
    logical :: l_min,l_max
    character(len=128) :: pro_name='COSP_CHECK_INPUT_3D'
    
    l_min=.false.
    l_max=.false.
    
    if (present(min_val)) then
!       if (x < min_val) x = min_val
      if (any(x < min_val)) then 
      l_min = .true.
        where (x < min_val)
          x = min_val
        end where
      endif
    endif    
    if (present(max_val)) then
!       if (x > max_val) x = max_val
      if (any(x > max_val)) then 
        l_max = .true.
        where (x > max_val)
          x = max_val
        end where  
      endif    
    endif    
    
    if (l_min) print *,'----- WARNING: '//trim(pro_name)//': minimum value of '//trim(vname)//' set to: ',min_val
    if (l_max) print *,'----- WARNING: '//trim(pro_name)//': maximum value of '//trim(vname)//' set to: ',max_val
  END SUBROUTINE COSP_CHECK_INPUT_3D


END MODULE MOD_COSP_UTILS
