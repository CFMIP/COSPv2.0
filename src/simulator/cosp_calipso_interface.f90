MODULE MOD_COSP_CALIPSO_INTERFACE
  USE COSP_KINDS,              ONLY: wp
  USE MOD_LIDAR_SIMULATOR,     ONLY: alpha,beta,gamma
  IMPLICIT NONE
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! TYPE calipso_in
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  type calipso_IN
     integer,pointer ::       &
          Npoints,      & ! Number of gridpoints.
          Ncolumns,     & ! Number of columns.
          Nlevels         ! Number of levels.

     real(wp),dimension(:,:),pointer :: &
          beta_mol,     & ! Molecular backscatter coefficient
          tau_mol         ! Molecular optical depth
     real(wp),dimension(:,:,:),pointer :: &
          betatot,      & ! 
          tautot,       & ! Optical thickess integrated from top
          betatot_ice,  & ! Backscatter coefficient for ice particles
          betatot_liq,  & ! Backscatter coefficient for liquid particles
          tautot_ice,   & ! Total optical thickness of ice
          tautot_liq      ! Total optical thickness of liq
     real(wp),dimension(:,:,:,:),pointer :: &
          taupart
  end type calipso_IN

CONTAINS
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_calipso_init
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine cosp_calipso_init() 
    
    ! Polynomial coefficients (Alpha, Beta, Gamma) which allow to compute the 
    ! ATBperpendicular as a function of the ATB for ice or liquid cloud particles 
    ! derived from CALIPSO-GOCCP observations at 120m vertical grid 
    ! (Cesana and Chepfer, JGR, 2013).
    !
    ! Relationship between ATBice and ATBperp,ice for ice particles:
    !                ATBperp,ice = Alpha*ATBice 
    ! Relationship between ATBice and ATBperp,ice for liquid particles:
    !          ATBperp,ice = Beta*ATBice^2 + Gamma*ATBice
    Alpha = 0.2904_wp
    Beta  = 0.4099_wp
    Gamma = 0.009_wp    
    
  end subroutine cosp_calipso_init

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !			 	    END MODULE
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE MOD_COSP_CALIPSO_INTERFACE
