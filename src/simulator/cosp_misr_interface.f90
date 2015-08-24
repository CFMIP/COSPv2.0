MODULE MOD_COSP_MISR_INTERFACE
  USE COSP_KINDS,              ONLY: wp
 
  IMPLICIT NONE

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! 								TYPE misr_in
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  type misr_IN
     integer,pointer  ::  &
          Npoints,        & ! Number of gridpoints.
          Ncolumns,       & ! Number of columns.
          Nlevels           ! Number of levels.
     integer,pointer ::   &
          sunlit(:)         ! Sunlit points (npoints).
     real(wp),pointer ::  &
          zfull(:,:),     & ! Height of full model levels (i.e. midpoints). (npoints,nlev)
          at(:,:)           ! Temperature. (npoints,nlev)
     real(wp),pointer ::  &                           
          dtau(:,:,:)       ! Optical depth. (npoints,ncolumns,nlev)
          
  end type misr_IN

CONTAINS

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !  							SUBROUTINE cosp_misr_init
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_MISR_INIT()

  END SUBROUTINE COSP_MISR_INIT

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! 							    	END MODULE
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE MOD_COSP_MISR_INTERFACE
