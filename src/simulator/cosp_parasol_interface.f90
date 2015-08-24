MODULE MOD_COSP_PARASOL_INTERFACE
  USE COSP_KINDS,  ONLY: WP
  implicit none

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !									TYPE cosp_parasol
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  TYPE PARASOL_SGX
     ! Dimensions
     integer :: &
          Npoints,  & ! Number of gridpoints
          Ncolumns, & ! Number of columns
          Nrefl       ! Number of parasol reflectances
     
     ! Arrays with dimensions (Npoints,Ncolumns,Nrefl)
     real(wp),dimension(:,:,:),pointer :: &
          refl        ! parasol reflectances

  END TYPE PARASOL_SGX
  TYPE PARASOL_GBX
     integer :: &
          Npoints,  & ! Number of gridpoints
          Ncolumns, & ! Number of columns
          Nrefl       ! Number of parasol reflectances
     real(wp), dimension(:,:),pointer :: &
          parasolrefl ! Mean parasol reflectance

  END TYPE PARASOL_GBX
  TYPE COSP_PARASOL 
     type(PARASOL_SGX) :: PARASOL_SGX
     type(PARASOL_GBX) :: PARASOL_GBX
  END TYPE COSP_PARASOL
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !										TYPE parasol_in
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  TYPE parasol_IN
     integer,pointer :: &
        Npoints,       & ! Number of horizontal gridpoints
        Nlevels,       & ! Number of vertical levels
        Ncolumns,      & ! Number of columns
        Nrefl            ! Number of angles for which the reflectance is computed
     real(wp),dimension(:,:),pointer ::   &
        tautot_S_liq,  & ! Liquid water optical thickness, from TOA to SFC
        tautot_S_ice     ! Ice water optical thickness, from TOA to SFC
  END TYPE parasol_IN
  
contains

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !                           SUBROUTINE cosp_parasol_init
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_PARASOL_INIT()
    
  end subroutine COSP_PARASOL_INIT

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! 								    END MODULE
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module MOD_COSP_PARASOL_INTERFACE
