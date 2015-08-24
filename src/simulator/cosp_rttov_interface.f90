MODULE MOD_COSP_RTTOV_INTERFACE
  USE COSP_KINDS,               ONLY: wp
  USE MOD_COSP_CONFIG,          ONLY: RTTOV_MAX_CHANNELS
  USE MOD_COSP_RTTOV,           ONLY: nch_in,plat_in,sat_in,sens_in,ichan_in
   
  IMPLICIT NONE

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !									TYPE rttov_in
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  type rttov_in
     integer,pointer :: &
          nprof,        & ! Number of profiles to simulate
          nlevels         ! Number of pressure levels
     real(wp),pointer :: &
          zenang,       & ! Satellite zenith angle
          co2,          & ! Carbon dioxide 
          ch4,          & ! Methane 
          n2o,          & ! n2o 
          co              ! Carbon monoxide
     real(wp),dimension(:),pointer :: &
          surfem          ! Surface emissivities for the channels
     real(wp),dimension(:),pointer :: &
          h_surf,       & ! Surface height
          u_surf,       & ! U component of surface wind
          v_surf,       & ! V component of surface wind
          t_skin,       & ! Surface skin temperature
          p_surf,       & ! Surface pressure
          t_surf,       & ! 1.5 m Temperature
          q_surf,       & ! 1.5 m Specific humidity
          lsmask,       & ! land-sea mask
          latitude        ! Latitude
     real(wp),dimension(:,:),pointer :: &
          p,            & ! Pressure profiles  
          t,            & ! Temperature profiles
          q,            & ! Humidity profiles
          o3              ! Ozone profiles
  end type rttov_in
CONTAINS

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !								 SUBROUTINE cosp_rttov_init
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_RTTOV_INIT(Nchannels,platform,satellite,instrument,channels)
    integer,intent(in) :: & 
       Nchannels, & ! Number of channels
       platform,  & ! Satellite platform
       satellite, & ! Satellite
       instrument   ! Instrument
     integer,intent(in),dimension(RTTOV_MAX_CHANNELS) :: &
       channels     ! RTTOV channels  
    
    ! Initialize
    nch_in   = Nchannels
    plat_in  = platform
    sat_in   = satellite
    sens_in  = instrument
    ichan_in = channels
    
  END SUBROUTINE COSP_RTTOV_INIT
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! 									END MODULE
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE MOD_COSP_RTTOV_INTERFACE
