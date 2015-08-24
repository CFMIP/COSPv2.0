MODULE MOD_COSP_ISCCP_INTERFACE
  USE COSP_KINDS,      ONLY: wp
  USE mod_icarus,      ONLY: isccp_top_height,isccp_top_height_direction
  IMPLICIT NONE
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !                                  TYPE isccp_in
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Derived input type for ISCCP simulator
  type isccp_IN
     integer,pointer  ::       &
          Npoints,             & ! Number of gridpoints.
          Ncolumns,            & ! Number of columns.
          Nlevels,             & ! Number of levels.
          top_height,          & !
          top_height_direction   !
     integer,pointer ::        &
          sunlit(:)              ! Sunlit points (npoints)
     real(wp),pointer ::       &
          emsfc_lw
     real(wp),pointer ::       &
          skt(:)                 ! Surface temperature (npoints)
     real(wp),pointer ::       &
          at(:,:),             & ! Temperature (npoint,nlev)
          pfull(:,:),          & ! Pressure (npoints,nlev)
          qv(:,:)                ! Specific humidity (npoints,nlev)
     real(wp),pointer  ::       &          
          phalf(:,:)             ! Pressure at half levels (npoints,nlev+1)
     real(wp),pointer ::       &
          frac_out(:,:,:),     & ! Cloud fraction (npoints,ncolumns,nlevels)
          dtau(:,:,:),         & ! Optical depth (npoints,ncolumns,nlevels)
          dem(:,:,:)             ! Emissivity (npoints,ncolumns,nlevels)
  end type isccp_IN

CONTAINS

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !  							SUBROUTINE cosp_isccp_init
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_ISCCP_INIT(top_height,top_height_direction)
     integer,intent(in) :: &
         top_height, &
         top_height_direction

    ! Cloud-top height determination
    isccp_top_height           = top_height
    isccp_top_height_direction = top_height_direction

  END SUBROUTINE COSP_ISCCP_INIT

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !                                    END MODULE
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE MOD_COSP_ISCCP_INTERFACE
