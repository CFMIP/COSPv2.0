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
! May 2007: ActSim code of M. Chiriaco and H. Chepfer rewritten by S. Bony
!
! May 2008, H. Chepfer:
! - Units of pressure inputs: Pa 
! - Non Spherical particles : LS Ice NS coefficients, CONV Ice NS coefficients
! - New input: ice_type (0=ice-spheres ; 1=ice-non-spherical)
!
! June 2008, A. Bodas-Salcedo:
! - Ported to Fortran 90 and optimisation changes
!
! August 2008, J-L Dufresne:
! - Optimisation changes (sum instructions suppressed)
!
! October 2008, S. Bony,  H. Chepfer and J-L. Dufresne :  
! - Interface with COSP v2.0:
!      cloud fraction removed from inputs
!      in-cloud condensed water now in input (instead of grid-averaged value)
!      depolarisation diagnostic removed
!      parasol (polder) reflectances (for 5 different solar zenith angles) added
!
! December 2008, S. Bony,  H. Chepfer and J-L. Dufresne : 
! - Modification of the integration of the lidar equation.
! - change the cloud detection threshold
!
! April 2008, A. Bodas-Salcedo:
! - Bug fix in computation of pmol and pnorm of upper layer
!
! April 2008, J-L. Dufresne
! - Bug fix in computation of pmol and pnorm, thanks to Masaki Satoh: a factor 2 
! was missing. This affects the ATB values but not the cloud fraction. 
!
! January 2013, G. Cesana and H. Chepfer:
! - Add the perpendicular component of the backscattered signal (pnorm_perp_tot) in the arguments
! - Add the temperature for each levels (temp) in the arguments
! - Add the computation of the perpendicular component of the backscattered lidar signal 
! Reference: Cesana G. and H. Chepfer (2013): Evaluation of the cloud water phase
! in a climate model using CALIPSO-GOCCP, J. Geophys. Res., doi: 10.1002/jgrd.50376
!
! May 2015 - D. Swales - Modified for COSPv2.0
!
! Apr 2018 - R. Guzman - Added Ground LIDar (GLID) subroutines
! Reference
!
! GLID: Chiriaco et al. (2018): ReOBS: a new approach to synthetize long-term 
! multi-variable dataset and application to the SIRTA supersite. ESSD (in press)
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module mod_groundlidar_simulator
  USE COSP_KINDS,         ONLY: wp
  USE MOD_COSP_CONFIG,    ONLY: SR_BINS,S_CLD,S_ATT,S_CLD_ATT,R_UNDEF,groundlidar_histBsct,  &
                                use_vgrid,vgrid_zl,vgrid_zu,vgrid_z
  USE MOD_COSP_STATS,     ONLY: COSP_CHANGE_VERTICAL_GRID,hist1d
  implicit none
  
  ! Polynomial coefficients (Alpha, Beta, Gamma) which allow to compute the 
  ! ATBperpendicular as a function of the ATB for ice or liquid cloud particles 
  ! derived from CALIPSO-GOCCP observations at 120m vertical grid 
  ! (Cesana and Chepfer, JGR, 2013).
  !
  ! Relationship between ATBice and ATBperp,ice for ice particles:
  !                ATBperp,ice = Alpha*ATBice 
  ! Relationship between ATBice and ATBperp,ice for liquid particles:
  !          ATBperp,ice = Beta*ATBice^2 + Gamma*ATBice
!  real(wp) :: &
!       alpha,beta,gamma    

contains
  ! ######################################################################################
  ! The subroutines below compute the attenuated backscatter signal and the lidar 
  ! backscatter coefficients using eq (1) from doi:0094-8276/08/2008GL034207
  ! ######################################################################################
  subroutine cmp_backsignal(nlev,npoints,beta,tau,pnorm)
    ! INPUTS
    integer, intent(in) :: nlev,npoints
    real(wp),intent(in),dimension(npoints,nlev) :: beta,tau

    ! OUTPUTS
    real(wp),intent(out),dimension(npoints,nlev) :: pnorm

    ! Internal Variables
    real(wp), dimension(npoints) :: tautot_lay
    integer :: k

    ! Uppermost layer 
    pnorm(:,1) = beta(:,1) / (2._wp*tau(:,1)) * (1._wp-exp(-2._wp*tau(:,1)))

    ! Other layers
    do k=2,nlev
       tautot_lay(:) = tau(:,k)-tau(:,k-1) 
       WHERE ( EXP(-2._wp*tau(:,k-1)) .gt. 0. )
          WHERE (tautot_lay(:) .gt. 0.)
             pnorm(:,k) = beta(:,k)*EXP(-2._wp*tau(:,k-1)) /&
                  (2._wp*tautot_lay(:))*(1._wp-EXP(-2._wp*tautot_lay(:)))
          ELSEWHERE
             ! This must never happen, but just in case, to avoid div. by 0
             pnorm(:,k) = beta(:,k) * EXP(-2._wp*tau(:,k-1))
          END WHERE
       ELSEWHERE
          pnorm(:,k) = 0._wp!beta(:,k)
       END WHERE
    END DO
  end subroutine cmp_backsignal

! BEGINNING OF GLID CHANGES
  ! ######################################################################################
  ! SUBROUTINE groundlidar_subcolumn  FROM THE GROUND
  ! ######################################################################################
  subroutine groundlidar_subcolumn(npoints,ncolumns,nlev,beta_mol_gr,tau_mol_gr,betatot_gr, &
                                   tautot_gr, pmol_gr, pnorm_gr)

    ! INPUTS
    INTEGER,intent(in) :: & 
         npoints,      & ! Number of gridpoints
         ncolumns,     & ! Number of subcolumns
         nlev            ! Number of levels
    REAL(WP),intent(in),dimension(npoints,nlev) :: &
         beta_mol_gr,  & ! Molecular backscatter coefficient
         tau_mol_gr      ! Molecular optical depth

    REAL(WP),intent(in),dimension(npoints,ncolumns,nlev)       :: &
         betatot_gr,   & ! 
         tautot_gr       ! Optical thickess integrated from top

    ! OUTPUTS
    REAL(WP),intent(out),dimension(npoints,nlev) :: &
         pmol_gr         ! Molecular attenuated backscatter lidar signal power(m^-1.sr^-1)
    REAL(WP),intent(out),dimension(npoints,ncolumns,nlev) :: &
         pnorm_gr        ! Molecular backscatter signal power (m^-1.sr^-1)

    ! LOCAL VARIABLES
    INTEGER :: k,icol

! we flip the profiles in calling the subroutines so the computation usually made from
! TOA to SFC is done the other way around for the ground lidar, from SFC to TOA
    ! ####################################################################################
    ! *) Molecular signal
    ! ####################################################################################
    call cmp_backsignal(nlev,npoints,beta_mol_gr(1:npoints,nlev:1:-1),&
                        tau_mol_gr(1:npoints,nlev:1:-1),pmol_gr(1:npoints,nlev:1:-1))

    do icol=1,ncolumns
       ! #################################################################################
       ! *) Total Backscatter signal
       ! #################################################################################
       call cmp_backsignal(nlev,npoints,betatot_gr(1:npoints,icol,nlev:1:-1),&
            tautot_gr(1:npoints,icol,nlev:1:-1),pnorm_gr(1:npoints,icol,nlev:1:-1))

    enddo

  end subroutine groundlidar_subcolumn

  ! ######################################################################################
  ! SUBROUTINE groundlidar_column  FROM THE GROUND
  ! ######################################################################################
  subroutine groundlidar_column(npoints,ncol,nlevels,llm,max_bin, pnorm_gr,              &
                                pmol_gr, pplay, ok_lidar_cfad_gr, ncat, cfad2_gr,  &
                                lidarcld_gr, cldlayer_gr, zlev, zlev_half)

    ! Inputs
    integer,intent(in) :: &
         npoints, & ! Number of horizontal grid points
         ncol,    & ! Number of subcolumns
         nlevels, & ! Number of vertical layers (OLD grid)
         llm,     & ! Number of vertical layers (NEW grid)
         max_bin, & ! Number of bins for SR CFADs
         ncat       ! Number of cloud layer types (low,mid,high,total)
    real(wp),intent(in),dimension(npoints,ncol,Nlevels) :: &
         pnorm_gr   ! Lidar ATB
    real(wp),intent(in),dimension(npoints,Nlevels) :: &
         pmol_gr, & ! Molecular ATB
         pplay      ! Pressure on model levels (Pa)
    logical,intent(in) :: &
         ok_lidar_cfad_gr ! True if GROUND lidar CFAD diagnostics need to be computed
    real(wp),intent(in),dimension(npoints,nlevels) :: &
         zlev        ! Model full levels
    real(wp),intent(in),dimension(npoints,nlevels+1) :: &
         zlev_half   ! Model half levels

    ! Outputs
    real(wp),intent(inout),dimension(npoints,llm) :: &
         lidarcld_gr   ! 3D "lidar" cloud fraction
    real(wp),intent(inout),dimension(npoints,ncat) :: &
         cldlayer_gr   ! "lidar" cloud layer fraction (low, mid, high, total)
    real(wp),intent(inout),dimension(npoints,max_bin,llm) :: &
         cfad2_gr      ! CFADs of GROUND lidar SR

    ! Local Variables
    integer :: ic,i,j
    real(wp),dimension(npoints,ncol,llm) :: &
         x3d_gr
    real(wp),dimension(npoints,llm) :: &
         x3d_c_gr, pnorm_c_gr
    real(wp)  :: &
         xmax
    real(wp),dimension(npoints,1,Nlevels) :: ph_in,betamol_in
    real(wp),dimension(npoints,ncol,llm)  :: pnormFlip_gr
    real(wp),dimension(npoints,1,llm)     :: pplayFlip, betamolFlip

    ! Vertically regrid input data
    if (use_vgrid) then 
       ph_in(:,1,:) = pplay(:,nlevels:1:-1)
       call cosp_change_vertical_grid(Npoints,1,Nlevels,zlev(:,nlevels:1:-1),zlev_half(:,nlevels:1:-1),&
            ph_in,llm,vgrid_zl(llm:1:-1),vgrid_zu(llm:1:-1),pplayFlip(:,1,llm:1:-1))
       betamol_in(:,1,:) = pmol_gr(:,nlevels:1:-1)
       call cosp_change_vertical_grid(Npoints,1,Nlevels,zlev(:,nlevels:1:-1),zlev_half(:,nlevels:1:-1),&
            betamol_in,llm,vgrid_zl(llm:1:-1),vgrid_zu(llm:1:-1),betamolFlip(:,1,llm:1:-1))
       call cosp_change_vertical_grid(Npoints,Ncol,Nlevels,zlev(:,nlevels:1:-1),zlev_half(:,nlevels:1:-1),&
            pnorm_gr(:,:,nlevels:1:-1),llm,vgrid_zl(llm:1:-1),vgrid_zu(llm:1:-1),pnormFlip_gr(:,:,llm:1:-1))
    endif

    ! Initialization (The histogram bins, are set up during initialization and the
    ! maximum value is used as the upper bounds.)
    xmax = maxval(groundlidar_histBsct)

    ! Compute GROUND LIDAR scattering ratio
    if (use_vgrid) then
       do ic = 1, ncol
          pnorm_c_gr = pnormFlip_gr(:,ic,:)
          where ((pnorm_c_gr .lt. xmax) .and. (betamolFlip(:,1,:) .lt. xmax) .and.       &
                (betamolFlip(:,1,:) .gt. 0.0 ))
             x3d_c_gr = pnorm_c_gr/betamolFlip(:,1,:)
          elsewhere
             x3d_c_gr = R_UNDEF
          end where
          x3d_gr(:,ic,:) = x3d_c_gr
       enddo
       ! Diagnose cloud fractions for subcolumn GROUND lidar scattering ratios
       CALL COSP_CLDFRAC_GR(npoints,ncol,llm,ncat,x3d_gr,pnormFlip_gr,pplayFlip,       &
                            S_att,S_cld,S_cld_att,R_UNDEF,lidarcld_gr,cldlayer_gr)
    else
       do ic = 1, ncol
          pnorm_c_gr = pnorm_gr(:,ic,:)
          where ((pnorm_c_gr.lt.xmax) .and. (pmol_gr.lt.xmax) .and. (pmol_gr.gt. 0.0 ))
             x3d_c_gr = pnorm_c_gr/pmol_gr
          elsewhere
             x3d_c_gr = R_UNDEF
          end where
          x3d_gr(:,ic,:) = x3d_c_gr
       enddo
       ! Diagnose cloud fractions for subcolumn GROUND lidar scattering ratios
       CALL COSP_CLDFRAC_GR(npoints,ncol,nlevels,ncat,x3d_gr,pnorm_gr,pplay,      &
                            S_att,S_cld,S_cld_att,R_UNDEF,lidarcld_gr,cldlayer_gr)
    endif

    ! GROUND CFADs
    if (ok_lidar_cfad_gr) then
       ! CFADs of subgrid-scale GROUND lidar scattering ratios
       do i=1,Npoints
          do j=1,llm
             cfad2_gr(i,:,j) = hist1D(ncol,x3d_gr(i,:,j),SR_BINS,groundlidar_histBsct)
          enddo
       enddo
       where(cfad2_gr .ne. R_UNDEF) cfad2_gr=cfad2_gr/ncol
    endif 

    ! Unit conversions
    where(lidarcld_gr /= R_UNDEF)      lidarcld_gr      = lidarcld_gr*100._wp
    where(cldlayer_gr /= R_UNDEF)      cldlayer_gr      = cldlayer_gr*100._wp

  end subroutine groundlidar_column

    ! ####################################################################################
    ! SUBROUTINE cosp_cldfrac_gr
    ! Conventions: Ncat must be equal to 4
    ! ####################################################################################
    SUBROUTINE COSP_CLDFRAC_GR(Npoints,Ncolumns,Nlevels,Ncat,x_gr,ATB_gr,pplay,        &
                               S_att,S_cld,S_cld_att,undef,lidarcld_gr,cldlayer_gr)

	! Inputs
    integer,intent(in) :: &
       Npoints,  & ! Number of gridpoints
       Ncolumns, & ! Number of subcolumns
       Nlevels,  & ! Number of vertical levels
       Ncat        ! Number of cloud layer types
    real(wp),intent(in) :: &
       S_att,    & !
       S_cld,    & !
       S_cld_att,& ! New threshold for undefine cloud phase detection
       undef       ! Undefined value
    real(wp),intent(in),dimension(Npoints,Ncolumns,Nlevels) :: &
       x_gr,        & ! 
       ATB_gr         ! 3D attenuated backscatter
    real(wp),intent(in),dimension(Npoints,Nlevels) :: &
       pplay       ! Pressure

	! Outputs
    real(wp),intent(out),dimension(Npoints,Nlevels) :: &
       lidarcld_gr      ! 3D cloud fraction from GROUND
    real(wp),intent(out),dimension(Npoints,Ncat) :: &
       cldlayer_gr      ! Low, middle, high, total cloud fractions
    
    ! Local variables
    integer  :: &
       ip, k, iz, ic, ncol, nlev, i
    real(wp) :: &
       p1
    real(wp),dimension(Npoints,Nlevels) :: &
       nsub
    real(wp),dimension(Npoints,Ncolumns,Ncat) :: &
       cldlay,nsublay   
    real(wp),dimension(Npoints,Ncat) :: &
       nsublayer
    real(wp),dimension(Npoints,Ncolumns,Nlevels) :: &   
       cldy, & ! 
       srok    !

    ! ####################################################################################
	! 1) Initialize    
    ! ####################################################################################
    lidarcld_gr           = 0._wp
    nsub                  = 0._wp
    cldlay                = 0._wp
    nsublay               = 0._wp

    ! ####################################################################################
    ! 2) Cloud detection
    ! ####################################################################################
    do k=1,Nlevels
       ! Cloud detection at subgrid-scale:
       where ((x_gr(:,:,k) .gt. S_cld) .and. (x_gr(:,:,k) .ne. undef) )
          cldy(:,:,k)=1._wp
       elsewhere
          cldy(:,:,k)=0._wp
       endwhere
       
       ! Number of usefull sub-columns:
       where ((x_gr(:,:,k) .gt. S_att) .and. (x_gr(:,:,k) .ne. undef) )
          srok(:,:,k)=1._wp
       elsewhere
          srok(:,:,k)=0._wp
       endwhere
    enddo    

    ! ####################################################################################
    ! 3) Grid-box 3D cloud fraction and layered cloud fractions(ISCCP pressure categories)
    ! ####################################################################################
    do k=1,Nlevels
       do ic = 1, Ncolumns
          do ip = 1, Npoints

             iz=1
             p1 = pplay(ip,k)
             if ( p1.gt.0. .and. p1.lt.(440._wp*100._wp)) then ! high clouds
                iz=3
             else if(p1.ge.(440._wp*100._wp) .and. p1.lt.(680._wp*100._wp)) then ! mid clouds
                iz=2
             endif
             
             cldlay(ip,ic,iz) = MAX(cldlay(ip,ic,iz),cldy(ip,ic,k))
             cldlay(ip,ic,4)  = MAX(cldlay(ip,ic,4),cldy(ip,ic,k))
             lidarcld_gr(ip,k)   = lidarcld_gr(ip,k) + cldy(ip,ic,k)
             
             nsublay(ip,ic,iz) = MAX(nsublay(ip,ic,iz),srok(ip,ic,k))
             nsublay(ip,ic,4)  = MAX(nsublay(ip,ic,4),srok(ip,ic,k))
             nsub(ip,k)        = nsub(ip,k) + srok(ip,ic,k)
             
          enddo
       enddo
    enddo   
    
    ! Grid-box 3D cloud fraction
    where ( nsub(:,:).gt.0.0 )
       lidarcld_gr(:,:) = lidarcld_gr(:,:)/nsub(:,:)
    elsewhere
       lidarcld_gr(:,:) = undef
    endwhere
    
    ! Layered cloud fractions
    cldlayer_gr  = 0._wp
    nsublayer = 0._wp
    do iz = 1, Ncat
       do ic = 1, Ncolumns
          cldlayer_gr(:,iz)  = cldlayer_gr(:,iz)  + cldlay(:,ic,iz)
          nsublayer(:,iz) = nsublayer(:,iz) + nsublay(:,ic,iz)
       enddo
    enddo
    where (nsublayer(:,:) .gt. 0.0)
       cldlayer_gr(:,:) = cldlayer_gr(:,:)/nsublayer(:,:)
    elsewhere
       cldlayer_gr(:,:) = undef
    endwhere

    RETURN
  END SUBROUTINE COSP_CLDFRAC_GR
! END OF GLID CHANGES

end module mod_groundlidar_simulator
