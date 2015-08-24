MODULE MOD_COSP_CLOUDSAT_INTERFACE
  USE MOD_COSP_CONFIG, ONLY: DBZE_BINS,CFAD_ZE_MIN,CFAD_ZE_WIDTH,SR_BINS,DBZE_MAX,       &
                             DBZE_MIN
  USE COSP_KINDS,      ONLY: wp
  USE quickbeam,       ONLY: quickbeam_init,radar_cfg,Re_MAX_BIN,Re_BIN_LENGTH,          &
                             radar_at_layer_one
  IMPLICIT NONE
         
  ! Directory where LUTs will be stored
  character(len=120) :: RADAR_SIM_LUT_DIRECTORY = './'
  logical :: RADAR_SIM_LOAD_scale_LUTs_flag   = .false.
  logical :: RADAR_SIM_UPDATE_scale_LUTs_flag = .false.

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! TYPE cloudsat_IN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  type cloudsat_IN
     integer,pointer ::            &
          Npoints,         & ! Number of horizontal grid-points
          Nlevels,         & ! Number of vertical levels
          Ncolumns           ! Number of subcolumns
     real(wp),pointer ::   &
          hgt_matrix(:,:),   & ! Height of hydrometeors (km)
          z_vol(:,:,:),      & ! Effective reflectivity factor (mm^6/m^3)
          kr_vol(:,:,:),     & ! Attenuation coefficient hydro (dB/km)
          g_vol(:,:,:),      & ! Attenuation coefficient gases (dB/km)
          g_to_vol_in(:,:)     ! Gaseous atteunation, radar to vol (dB)
     type(radar_cfg),pointer :: rcfg   ! Radar simulator configuration
  end type cloudsat_IN

CONTAINS

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !                              SUBROUTINE cosp_cloudsat_in
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_CLOUDSAT_INIT(radar_freq,k2,use_gas_abs,do_ray,undef,nhydro,Npoints,   &
                                Nlevels,hgt_matrix,surface_radar,rcfg,                   &
                                cloudsat_micro_scheme,load_LUT)
    ! INPUTS
    real(wp),intent(in) :: &
         radar_freq,  & ! Radar frequency (GHz)
         k2,          & ! |K|^2, the dielectric constant
         undef          ! Undefined
    integer,intent(in) :: &
         use_gas_abs, & ! 1 = do gaseous abs calcs, 0=no gasesous absorbtion calculated,
                        ! 2 = calculate absorption for first profile on all profiles
         do_ray,      & !
         nhydro,      & !
         Npoints,     & !
         Nlevels,     & !
         surface_radar
    real(wp),intent(in),dimension(Npoints,Nlevels) :: &
         hgt_matrix     !
    logical,intent(in),optional :: &
         load_LUT       !
    character(len=64),intent(in) :: &
       cloudsat_micro_scheme
    
    ! OUTPUTS
    type(radar_cfg) :: &
         rcfg           !
    
    ! LOCAL VARIABLES
    character(len=240) :: LUT_file_name
    logical       :: local_load_LUT,hgt_descending
    integer       :: i,j
    
    if (present(load_LUT)) then
       local_load_LUT = load_LUT
    else
       local_load_LUT = RADAR_SIM_LOAD_scale_LUTs_flag
    endif
    
    write(*,*) 'RADAR_SIM microphysics scheme is set to: ',&
                trim(cloudsat_micro_scheme)
    
    ! LUT file name
    LUT_file_name = trim(RADAR_SIM_LUT_DIRECTORY) // &
         trim(cloudsat_micro_scheme)
    
    ! Initialize for NEW radar-configurarion derived type (radar_cfg)
    rcfg%freq                = radar_freq
    rcfg%k2                  = k2
    rcfg%use_gas_abs         = use_gas_abs
    rcfg%do_ray              = do_ray
    rcfg%nhclass             = nhydro
    rcfg%load_scale_LUTs     = local_load_LUT
    rcfg%update_scale_LUTs   = .false.
    rcfg%scale_LUT_file_name = LUT_file_name
    rcfg%N_scale_flag        = .false.
    rcfg%fc                  = undef
    rcfg%rho_eff             = undef
    rcfg%Z_scale_flag        = .false.
    rcfg%Ze_scaled           = 0._wp
    rcfg%Zr_scaled           = 0._wp
    rcfg%kr_scaled           = 0._wp
    
    ! Set up Re bin "structure" for z_scaling
    rcfg%base_list(1)=0
    do j=1,Re_MAX_BIN
       rcfg%step_list(j)=0.1_wp+0.1_wp*((j-1)**1.5)
       if(rcfg%step_list(j)>Re_BIN_LENGTH) then
          rcfg%step_list(j)=Re_BIN_LENGTH
       endif
       if(j>1) then
          rcfg%base_list(j)=rcfg%base_list(j-1)+floor(Re_BIN_LENGTH/rcfg%step_list(j-1))
       endif
    enddo
    
    ! Set flag denoting position of radar relative to hgt_matrix orientation
    hgt_descending = hgt_matrix(1,1) > hgt_matrix(1,size(hgt_matrix,2))
    if ((surface_radar == 1 .and. hgt_descending) .or.     &
         (surface_radar == 0 .and. (.not. hgt_descending))) then
       radar_at_layer_one = .false.
    else
       radar_at_layer_one = .true.
    endif

  END SUBROUTINE COSP_CLOUDSAT_INIT

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !                              	  END MODULE
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE MOD_COSP_CLOUDSAT_INTERFACE
