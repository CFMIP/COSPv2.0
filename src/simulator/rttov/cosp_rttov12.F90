! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copyright (c) 2016, Regents of the University of Colorado
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
! March 2016 - M. Johnston - Original version
! April 2016 - D. Swales   - Modified for use in COSPv2.0
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module mod_cosp_rttov
  use rttov_const,         only : errorstatus_success, errorstatus_fatal
  use rttov_types,         only : rttov_options, rttov_coefs, rttov_profile,                &
                                  rttov_transmission, rttov_radiance,rttov_chanprof,        &
                                  rttov_emissivity,rttov_profile_cloud,rttov_scatt_coef,  &
                                  rttov_options_scatt
  USE rttov_const,         ONLY : surftype_sea, surftype_land, surftype_seaice, &
                                  sensor_id_mw, sensor_id_po
  use rttov_unix_env,      only : rttov_exit

  USE mod_rttov_emis_atlas, ONLY : rttov_emis_atlas_data, &
                                   atlas_type_ir, atlas_type_mw

  use cosp_kinds,          only : wp
  use mod_cosp_config,     only : RTTOV_MAX_CHANNELS,N_HYDRO,rttovDir
  use cosp_phys_constants, only : mdry=>amd,mO3=>amO3,mco2=>amCO2,mCH4=>amCH4,           &
                                  mn2o=>amN2O,mco=>amCO
  implicit none
#include "rttov_direct.interface"
!#include "rttov_alloc_prof.interface"
!#include "rttov_alloc_rad.interface"
!#include "rttov_alloc_transmission.interface"
#include "rttov_alloc_direct.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_user_options_checkinput.interface"
#include "rttov_read_coefs.interface"
#include "rttov_get_emis.interface"
#include "rttov_boundaryconditions.interface"

  ! Module parameters
  integer, parameter :: maxlim =  10000
  real(wp),parameter :: eps    =  0.622

  integer :: atlas_type

  ! Initialization parameters
  integer :: &
       platform,   & ! RTTOV platform
       sensor,     & ! RTTOV instrument
       satellite,  & ! RTTOV satellite
       nChannels     ! Number of channels
  integer,dimension(RTTOV_MAX_CHANNELS) :: &
       iChannel      ! RTTOV channel numbers

  ! Scattering coefficients (read in once during initialization)
  TYPE(rttov_coefs),SAVE :: &
       coef_rttov
  TYPE(rttov_scatt_coef),SAVE :: &
       coef_scatt
  ! RTTOV setup and options (set during initialization)
  TYPE(rttov_options),SAVE :: &
       opts     ! defaults to everything optional switched off
  TYPE(rttov_options_scatt),SAVE :: &
       opts_scatt

  TYPE(rttov_emis_atlas_data),SAVE :: emis_atlas(12)
contains

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE rttov_column
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine rttov_column(nPoints,nLevels,nSubCols,q,p,t,o3,ph,h_surf,u_surf,v_surf,     &
                          p_surf,t_skin,t2m,q2m,lsmask,lon,lat,seaice,co2,ch4,n2o,co,    &
                          zenang,lCleanup,                                               &
                          ! Outputs
                          Tb,error,                                                      &
                          ! Optional arguments for surface emissivity calculation.
                          surfem,month,                                                  &
                          ! Optional arguments for all-sky calculation.
                          tca,ciw,clw,rain,snow)
    ! Inputs
    integer,intent(in) :: &
         nPoints, & ! Number of gridpoints
         nLevels, & ! Number of vertical levels
         nSubCols   ! Number of subcolumns
    real(wp),intent(in) :: &
         co2,     & ! CO2 mixing ratio (kg/kg)
         ch4,     & ! CH4 mixing ratio (kg/kg)
         n2o,     & ! N2O mixing ratio (kg/kg)
         co,      & ! CO mixing ratio (kg/kg)
         zenang     ! Satellite zenith angle
    real(wp),dimension(nPoints),intent(in) :: &
         h_surf,  & ! Surface height (m)
         u_surf,  & ! Surface u-wind (m/s)
         v_surf,  & ! Surface v-wind (m/s)
         p_surf,  & ! Surface pressure (Pa)
         t_skin,  & ! Skin temperature (K)
         t2m,     & ! 2-meter temperature (K)
         q2m,     & ! 2-meter specific humidity (kg/kg)
         lsmask,  & ! Land/sea mask
         lon,     & ! Longitude (deg)
         lat,     & ! Latitude (deg)
         seaice     ! Seaiec fraction (0-1)
    real(wp),dimension(nPoints,nLevels),intent(in) :: &
         q,       & ! Specific humidity (kg/kg)
         p,       & ! Pressure(Pa)
         t,       & ! Temperature (K)
         o3         ! Ozone
    real(wp),dimension(nPoints,nLevels+1),intent(in) :: &
         ph         ! Pressure @ half-levels (Pa)
    logical,intent(in) :: &
         lCleanup   ! Flag to determine whether to deallocate RTTOV types

    ! Optional inputs (Needed for surface emissivity calculation)
    integer,optional :: &
         month      ! Month (needed to determine table to load)
    real(wp),dimension(nChannels),optional :: &
         surfem     ! Surface emissivity for each RTTOV channel

    ! Optional inputs (Needed for all-sky calculation)
    real(wp),dimension(nPoints,nLevels),optional :: &
         tca       ! Total column cloud amount (0-1)
    real(wp),dimension(nPoints,nSubCols,nLevels),optional :: &
         ciw,    & ! Cloud ice
         clw,    & ! Cloud liquid
         rain,   & ! Precipitation flux (kg/m2/s)
         snow      ! Precipitation flux (kg/m2/s)

    ! Outputs
    real(wp),dimension(nPoints,nChannels) :: &
         Tb        ! RTTOV brightness temperature.
    character(len=128) :: &
         error     ! Error messages (only populated if error encountered)
    
    ! Local variables
    integer :: &
         nloop,rmod,il,istart,istop,za,i,j,subcol,errorstatus,npts_it
    integer,dimension(60) :: &
         alloc_status
    real(wp),dimension(nPoints) :: &
         sh_surf
    real(wp),dimension(nPoints,nLevels) :: &
         sh,totalice
    real(wp),dimension(nSubCols,nPoints,nChannels) :: &
         Tbs ! Subcolumn brightness temperature
    logical :: &
         use_totalice, mmr_snowrain, cfrac
    logical :: &
         lallSky, & ! Control for type of brightness temperature calculation
                    ! (False(default) => clear-sky brightness temperature, True => All-sky)
         lsfcEmis   ! Control for surface emissivity calculation (true => compute surface emissivity,
                    ! provided that the field "month" is available)

#include "rttov_read_coefs.interface"
#include "rttov_read_scattcoeffs.interface"
#include "rttov_user_options_checkinput.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_dealloc_scattcoeffs.interface"
#include "rttov_setup_emis_atlas.interface"
#include "rttov_deallocate_emis_atlas.interface"
#include "rttov_print_opts.interface"
#include "rttov_print_profile.interface"
#include "rttov_boundaryconditions.interface"

    ! Initialize some things
    totalice   = 0._wp
    Tbs(:,:,:) = 0._wp
!    Tb(:,:)    = 0._wp
    error      = ''

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Setup for call to RTTOV
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! First, check to see if we are doing an all-sky or clear-sky calculation brightness
    ! temperature
    lallSky = .false.
    if (present(tca) .and. present(clw) .and. present(ciw) .and. present(rain)           &
        .and. present(snow)) lallSky=.true.

    ! Check to see if we need to compute the surface emissivity (defualt is to compute
    ! surface emissivity using the atlas tables)
    lsfcEmis = .true.
    if (present(surfem)) lsfcEmis = .false.
    
    ! We also need the month for the emissivity atlas, so check...
    if (.not. present(month)) lsfcEmis = .false.

    if (lsfcEmis .eqv. .false. .and. .not. present(surfem)) then
       error = 'ERROR (rttov_column): User did not provide surface emissivity and did not '//&
               'request the surface emissivity to be calculated!!!'
       return
    endif

    ! Convert specific humidity to ppmv
    sh       = ( q   / ( q   + eps * ( 1._wp - q ) ) ) * 1e6   
    sh_surf  = ( q2m / ( q2m + eps * ( 1._wp - q2m ) ) ) * 1e6

    ! Settings unique to all-sky call.
    use_totalice = .false.
    mmr_snowrain = .true.
    cfrac        = .true.
    opts_scatt%lusercfrac = cfrac

    nloop  =  npoints / maxlim
    rmod   =  mod( npoints, maxlim )
    if( rmod .ne. 0 ) then
       nloop = nloop + 1
    endif

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Initialize emissivity atlas data for chosen sensor.
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    IF (coef_rttov%coef%id_sensor == sensor_id_mw .OR. &
        coef_rttov%coef%id_sensor == sensor_id_po) THEN
      atlas_type = atlas_type_mw ! mw atlas
    ELSE
      atlas_type = atlas_type_ir ! ir atlas
    ENDIF
    
    DO j = 1,12
      CALL rttov_setup_emis_atlas(errorstatus, &
                                  opts,        &
                                  j,           &
                                  atlas_type,  & ! selects mw (1) or ir (2)
                                  emis_atlas(j),  & 
                                  path = TRIM(rttovdir) // 'emis_data', &
                                  coefs = coef_rttov)
    ENDDO

    IF (errorstatus /= errorstatus_success) THEN
      error = 'ERROR (rttov_column): Error reading emis atlas data!'
      RETURN
    ENDIF
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Some quality control prior to RTTOV call
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Ensure the options and coefficients are consistent
    if(opts_scatt%config%do_checkinput) then
       call rttov_user_options_checkinput(errorstatus, opts, coef_rttov)
       if (errorstatus /= errorstatus_success) then
          error =  'ERROR (rttov_column): Error when checking input data!'
          return
       endif
    endif
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Call to RTTOV
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Looping over maxlim number of profiles
    DO il = 1, nloop
       istart  =  (il - 1) * maxlim + 1
       istop   =  min(il * maxlim, npoints)     
       if( ( il .eq. nloop ) .and. ( rmod .ne. 0 ) ) then
          npts_it = rmod
       else
          npts_it = maxlim
       endif
       ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Clear-sky brightness temperature
       ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
       IF (.NOT. lallSky) THEN
         CALL rttov_multprof(nChannels,iChannel,surfem,npts_it,nLevels,platform,         &
                             satellite,sensor,opts,coef_rttov,zenang,                    &
                             p(istart:istop,:)/100._wp,t(istart:istop,:),                &
                             sh(istart:istop,:),(mdry/mo3)*o3(istart:istop,:)*1e6,       &
                             (mdry/mco2)*co2*1e6,(mdry/mch4)*ch4*1e6,(mdry/mn2o)*n2o*1e6,&
                             (mdry/mco)*co*1e6,h_surf(istart:istop),u_surf(istart:istop),&
                             v_surf(istart:istop),t_skin(istart:istop),                  &
                             p_surf(istart:istop)/100.,t2m(istart:istop),                &
                             sh_surf(istart:istop),lsmask(istart:istop),                 &
                             seaice(istart:istop),lat(istart:istop),lon(istart:istop),   &
                             Tbs(1,istart:istop,:),                                      &
                             month = month)  
       ELSE
         DO subcol = 1, nSubCols
           CALL rttov_multprof(nChannels,iChannel,surfem,npts_it,nLevels,platform,         &
                                  satellite,sensor,opts,coef_rttov,zenang,                    &
                                  p(istart:istop,:)/100._wp,t(istart:istop,:),                &
                                  sh(istart:istop,:),(mdry/mo3)*o3(istart:istop,:)*1e6,       &
                                  (mdry/mco2)*co2*1e6,(mdry/mch4)*ch4*1e6,(mdry/mn2o)*n2o*1e6,&
                                  (mdry/mco)*co*1e6,h_surf(istart:istop),u_surf(istart:istop),&
                                  v_surf(istart:istop),t_skin(istart:istop),                  &
                                  p_surf(istart:istop)/100.,t2m(istart:istop),                &
                                  sh_surf(istart:istop),lsmask(istart:istop),                 &
                                  seaice(istart:istop),lat(istart:istop),lon(istart:istop),   &
                                  Tbs(subcol,istart:istop,:),                                      &
                                  opts_scatt = opts_scatt, &
                                  coef_scatt = coef_scatt, &
                                  ph_in = ph(istart:istop,:)/100._wp,                                 &
!                                 tca = tca(istart:istop, :), &
                                  ciw_in = ciw(istart:istop,subcol,:),     &
                                  clw_in = clw(istart:istop,subcol,:),                          &
                                  rain_in = rain(istart:istop,subcol,:), &
                                  sp_in = snow(istart:istop,subcol,:),& 
                                  totalice_in = totalice(istart:istop,:), &
                                  use_totalice = use_totalice,  &
                                  mmr_snowrain = mmr_snowrain,  &
                                  cfrac = cfrac, &
                                  month = month)
         ENDDO
       ENDIF
    !    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !    ! All-sky brightness temperature
    !    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !    if (lallSky) then
    !       ! Loop over all subcolumns
    !       do subcol = 1, nSubCols	
    !          ! Call RTTOV
    !          call cosp_rttov_mwscatt(nChannels,iChannel,surfem,nPoints,nlevels,platform,  &
    !                                  satellite,sensor,opts,opts_scatt,coef_rttov,         &
    !                                  coef_scatt,zenang,p(istart:istop,:)/100._wp,         &
    !                                  ph(istart:istop,:)/100._wp,t(istart:istop, :),       &
    !                                  sh(istart:istop, :),                                 &
    !                                  (mdry/mo3)*o3(istart:istop,:)*1e6,                   &
    !                                  clw(istart:istop,subcol,:),                          &
    !                                  ciw(istart:istop,subcol,:),tca(istart:istop, :),     &
    !                                  totalice(istart:istop,:),snow(istart:istop,subcol,:),& 
    !                                  rain(istart:istop,subcol,:),(mdry/mco2)*co2*1e6,     &
    !                                  (mdry/mch4)*ch4*1e6,(mdry/mn2o)*n2o*1e6,             &
    !                                  (mdry/mco)*co*1e6,h_surf(istart:istop),              &
    !                                  u_surf(istart:istop),v_surf(istart:istop),           &
    !                                  t_skin(istart:istop), p_surf(istart:istop)/100.,     &
    !                                  t2m(istart:istop),sh_surf(istart:istop),             &
    !                                  lsmask(istart:istop),seaice(istart:istop),           &
    !                                  lat(istart:istop),lon(istart:istop), use_totalice,   &
    !                                  mmr_snowrain,cfrac,Tbs(subcol,istart:istop,:))
    !       enddo 
    !    endif 
     ENDDO

    ! For all-sky calculation we need to average together all of the cloudy subcolumns.
    if (lallSky) then
       do subcol = 1, nSubCols
          Tb = Tb + tbs(subcol,:,:)
       enddo
       Tb = Tb/nSubCols
     ELSE
       Tb = Tbs(1,:,:)
    endif
    
    ! Free up space
    if (lCleanup) then
       call rttov_dealloc_coefs(errorstatus,coef_rttov)
       DO j = 1,12
         CALL rttov_deallocate_emis_atlas(emis_atlas(j))
       ENDDO
       if (lallSky) call rttov_dealloc_scattcoeffs(coef_scatt)    
    endif
  end subroutine rttov_column

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE rttov_multprof
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine rttov_multprof( &
       nch_in,     & ! number of channels
       ichan_in,   & ! channel indices
       surfem_in,  & ! surface emissivity values
       prf_num_in, & ! number of profiles to simulate
       nlevels_in, & ! number of pressure levels
       plat_in,    & ! platform number
       sat_in,     & ! satellite number
       sens_in,    & ! instrument number
       opts,       &
       coef_rttov, &
       zenang_in,  & ! zenith angle
       p_in,       & ! pressure [hpa]
       t_in,       & ! temperature [ k ]
       q_in,       & ! specific humidity [ ppmv ]
       o3_in,      & ! ozone vmr [ ppmv ]
       co2_in,     & ! co2 vmr [ ppmv ] *this is a single value*
       ch4_in,     & ! ch4 vmr [ ppmv ] *this is a single value*
       n2o_in,     & ! n2o vmr [ ppmv ] *this is a single value*
       co_in,      & ! co vmr [ ppmv ]  *this is a single value*
       h_surf,     & ! surface height [ m ]
       u_surf,     & ! u wind at 10 m  [ m/s ]
       v_surf,     & ! v wind at 10 m [ m/s ]
       t_skin,     & ! skin temperatre [ k ]
       p_surf,     & ! surface pressure
       t_surf,     & ! 1.5 m temperature [ k ]
       q_surf,     & ! 1.5 m specific humidity [ ppmv ]
       lsmask,     & ! land sea mask
       seaice,     & ! seaice fraction
       latitude,   & ! latitude [ deg north ]
       longitude,  & ! longitude [ deg east ]
       tbs,        & ! brightness temperature [ k ] (output)
       month,      &
       opts_scatt, &
       coef_scatt, &
       ph_in,      & ! pressure on half levels [hpa]
       clw_in,     & ! cloud water [0-1]
       ciw_in,     & ! cloud ice [0-1]
       cc_in,      & ! effective cloud fraction [0-1]
       totalice_in,& ! total ice, except snow [kg/kg] or [kg/m2/s]
       sp_in,      & ! solid precip with snow [kg/kg] or [kg/m2/s]
       rain_in,    & ! total liquid water [kg/kg] or [kg/m2/s]
       use_totalice, & ! 
       mmr_snowrain,& ! set units for snow and rain: if true units are kg/kg (the default)
       cfrac)         ! opts_scatt%lusercfrac=true., supply the effective cloud fraction


    !------ input arguments. no rttov kinds should be used here -----------------
    integer, intent(in)  :: nch_in             ! number of channels to be computed
    integer, intent(in)  :: ichan_in(nch_in)   ! indices of selected channels
    real(wp), intent(in)     :: surfem_in(nch_in)  ! surface emissivities for the channels
    integer, intent(in)  :: prf_num_in
    integer, intent(in)  :: nlevels_in
    integer, intent(in)  :: plat_in  ! satellite platform
    integer, intent(in)  :: sat_in   ! satellite number
    integer, intent(in)  :: sens_in  ! satellite sensor
    real(wp), intent(in)     :: zenang_in          ! satellite zenith angle

    TYPE(rttov_options), intent(in)  :: opts
    TYPE(rttov_coefs),intent(in)    :: coef_rttov

    real(wp), intent(in)     :: p_in(prf_num_in, nlevels_in)  ! pressure profiles
    real(wp), intent(in)     :: t_in(prf_num_in, nlevels_in)  ! temperature profiles
    real(wp), intent(in)     :: q_in(prf_num_in, nlevels_in)  ! humidity profiles
    real(wp), intent(in)     :: o3_in(prf_num_in, nlevels_in) ! ozone profiles

    ! the following trace gases contain constant values
    real(wp), intent(in) ::  co2_in ! carbon dioxide
    real(wp), intent(in) ::  ch4_in ! methane
    real(wp), intent(in) ::  n2o_in ! n2o
    real(wp), intent(in) ::  co_in  ! carbon monoxide
    real(wp), intent(in) ::  h_surf(prf_num_in)         ! surface height
    real(wp), intent(in) ::  u_surf(prf_num_in)         ! u component of surface wind
    real(wp), intent(in) ::  v_surf(prf_num_in)         ! v component of surface wind
    real(wp), intent(in) ::  t_skin(prf_num_in)         ! surface skin temperature
    real(wp), intent(in) ::  p_surf(prf_num_in)                  ! surface pressure
    real(wp), intent(in) ::  t_surf(prf_num_in)         ! 1.5 m temperature
    real(wp), intent(in) ::  q_surf(prf_num_in)         ! 1.5 m specific humidity
    real(wp), intent(in) ::  lsmask(prf_num_in)         ! land-sea mask
    real(wp), intent(in) ::  seaice(prf_num_in)         ! sea-ice fraction
    real(wp), intent(in) ::  latitude(prf_num_in)       ! latitude
    real(wp), intent(in) ::  longitude(prf_num_in)      ! longitude

    REAL(wp), INTENT(inout) :: tbs(prf_num_in, nch_in)  ! tbs (in the right format)

    TYPE(rttov_options_scatt), OPTIONAL  :: opts_scatt
    TYPE(rttov_scatt_coef), OPTIONAL     :: coef_scatt

    INTEGER, INTENT(in), OPTIONAL :: month

    REAL(wp), INTENT(in), OPTIONAL :: clw_in(prf_num_in, nlevels_in)
    REAL(wp), INTENT(in), OPTIONAL :: ciw_in(prf_num_in, nlevels_in)
    REAL(wp), INTENT(in), OPTIONAL :: cc_in(prf_num_in, nlevels_in)
    REAL(wp), INTENT(in), OPTIONAL :: totalice_in(prf_num_in, nlevels_in)
    REAL(wp), INTENT(in), OPTIONAL :: sp_in(prf_num_in, nlevels_in)
    REAL(wp), INTENT(in), OPTIONAL :: rain_in(prf_num_in, nlevels_in)
    REAL(wp), INTENT(in), OPTIONAL :: ph_in(prf_num_in, nlevels_in+1)
    LOGICAL, INTENT(in), OPTIONAL  :: use_totalice, cfrac, mmr_snowrain
    LOGICAL :: lallsky

    TYPE(rttov_profile_cloud), ALLOCATABLE :: cld_profiles(:)
    integer       , allocatable :: frequencies (:)

!#include "rttov_scatt_setupindex.interface"
include "rttov_scatt.interface"
include "rttov_alloc_rad.interface"
include "rttov_alloc_prof.interface"
include "rttov_alloc_scatt_prof.interface"
include "rttov_get_emis.interface"
include "rttov_boundaryconditions.interface"


    !------ local variables. use only rttov kinds or derived types.
    !       logical variables are declared with the same kind
    !       as integers, as they are affected inthe same way by flags like -qintsize=8

    TYPE(rttov_chanprof),    POINTER :: chanprof(:)    ! input channel/profile list
    TYPE(rttov_profile),     POINTER :: profiles(:)    ! input profiles
    LOGICAL,                 POINTER :: calcemis(:)    ! flag to indicate calculation of emissivity within rttov
    TYPE(rttov_emissivity),  POINTER :: emissivity(:)  ! input/output surface emissivity
    type(rttov_transmission)         :: transmission   ! output transmittances
    type(rttov_radiance)             :: radiance       ! output radiances

    INTEGER, ALLOCATABLE :: instrument(:,:) ! instrument id (3 x n_instruments)
    INTEGER, ALLOCATABLE :: nchan(:) ! number of channels per instrument
    INTEGER, ALLOCATABLE :: ichan(:,:)   ! channel list per instrument

    integer :: asw
    integer :: mxchn
    integer :: nrttovid     ! maximum number of instruments
    integer :: no_id        ! instrument loop index
    integer :: i, j, jch
    integer :: nprof  ! number of calls to rttov
    integer :: nch ! intermediate variable
    integer :: errorstatus
    integer :: ich, ich_temp, nchanprof, nchannels, chan
    integer :: alloc_status(60)

    real(wp),    allocatable :: input_emissivity (:)

    character (len=14)  :: nameofroutine = 'rttov_multprof'

    logical :: refrac, solrad, laerosl, lclouds, lsun, all_channels

    ! local variables for input arguments that need type casting to avoid type-mismatch with
    ! rttov kinds. this happens with some compiler flags (-qintsize=8).
    integer  :: prof_num
    integer  :: nlevels

    ! --------------------------------------------------------------------------
    ! 0. initialise cosp-specific things
    ! --------------------------------------------------------------------------

    ! type-casting of input arguments that need to be passed to rttov
    prof_num = prf_num_in
    nlevels  = nlevels_in
    nprof = prof_num

    ! currently we plan to calculate only 1 instrument per call
    nrttovid  =  1
    mxchn  =  nch_in

    errorstatus     = 0
    alloc_status(:) = 0

    !     allocate(coefs(nrttovid), stat = alloc_status(1))

    !     allocate(instrument(3, nrttovid), stat = alloc_status(4))

    !maximum number of channels allowed for one instrument is mxchn
    !    allocate(surfem(nch_in, nrttovid), stat = alloc_status(11))
    allocate(ichan(nch_in, nrttovid), stat = alloc_status(12))
    call rttov_error('ichan mem allocation error for profile array' , lalloc = .true.)


    do no_id = 1, nrttovid
       ichan(:, no_id)   = ichan_in
    enddo
    no_id = 1
    ! --------------------------------------------------------------------------
    ! 3. build the list of profile/channel indices in chanprof
    ! --------------------------------------------------------------------------
    
    ALLOCATE(nchan(nprof))     ! number of channels per profile
    nchan(:) = SIZE(ichan(:,no_id))  ! = nch_in
    
    ! size of chanprof array is total number of channels over all profiles
    ! square in this case - here same channels done for all profiles
    nchanprof = SUM(nchan(:))    

    ! allocate structures for rttov_direct
    CALL rttov_alloc_direct(   &
      errorstatus,             &
      1_wp,                    &  ! 1 => allocate
      nprof,                   &  ! nprofiles
      nchanprof,               &
      nlevels,                 &
      chanprof,                &
      opts,                    &
      profiles,                &
      coef_rttov,              &
      transmission,            &
      radiance,                &
      calcemis = calcemis,     &
      emissivity = emissivity, &
      init = .TRUE.)

    CALL rttov_error('allocation error for direct' , lalloc = .TRUE.)
!    ! pack channels and input emissivity arrays
!    chanprof(:)%chan = 0

    nch = 0
    DO j = 1, nprof
      DO jch = 1, nchan(j)
        nch = nch + 1
        chanprof(nch)%prof = j
        IF(ichan(jch, no_id) < 1) THEN
          errorstatus = errorstatus_fatal
          CALL rttov_error('Sensor channel number must be 1 or greater' , lalloc = .TRUE.)
        ELSE
          chanprof(nch)%chan = ichan(jch, no_id)
        ENDIF
      ENDDO
    ENDDO

    IF(lallsky) THEN
      ALLOCATE (cld_profiles(prf_num_in))
      CALL rttov_alloc_scatt_prof(errorstatus, &
        prf_num_in, cld_profiles, nlevels_in, use_totalice, asw = 1, init = .TRUE.)   

    ! setup indices
    call rttov_scatt_setupindex ( &
                                 prf_num_in,      			& ! in
                                 nch_in,          			& ! in
                                 coef_rttov, 	                        & ! in
                                 nchanprof,       			& ! in
                                 chanprof,        			& ! out
                                 frequencies)       		      ! out

    ENDIF
    ! end of 3.

    ! --------------------------------------------------------------------------
    ! 5. store profile data in profile type
    ! --------------------------------------------------------------------------
    DO i = 1, nprof
       profiles(i)%p(:) =  p_in(i, :)
       profiles(i)%t(:) =  t_in(i, :)
       profiles(i)%q(:) =  q_in(i, :)

       where(profiles(i)%q(:) < 1e-4)
          profiles(i)%q(:) = 1e-4
       end where

       profiles(i)%cfraction  =  0.
       profiles(i)%ctp        =  500.

       ! 2m parameters
       profiles(i)%s2m%p  =  p_surf(i)
       profiles(i)%s2m%t  =  t_surf(i)
       profiles(i)%s2m%q  =  q_surf(i)
       profiles(i)%s2m%u  =  u_surf(i) 
       profiles(i)%s2m%v  =  v_surf(i) 
       profiles(i)%s2m%wfetc  =  10000. 

       ! skin variables for emissivity calculations
       profiles(i)%skin%t          =  t_skin(i)

       ! fastem coefficients - for mw calculations
       profiles(i)%skin%fastem(1)  =  3.0
       profiles(i)%skin%fastem(2)  =  5.0
       profiles(i)%skin%fastem(3)  =  15.0
       profiles(i)%skin%fastem(4)  =  0.1
       profiles(i)%skin%fastem(5)  =  0.3

       profiles(i)%zenangle      = zenang_in ! pass in from cosp

       profiles(i)%azangle       = 0. ! hard-coded in rttov9 int

       profiles(i)%latitude      = latitude(i)
       profiles(i)%longitude     = longitude(i)
       profiles(i)%elevation     = h_surf(i)

       profiles(i)%sunzenangle   = 0. 
       profiles(i)%sunazangle    = 0. 
       
       ! surface type
       ! land-sea mask indicates proportion of land in grid
       if (lsmask(i) < 0.5) then
          profiles(i)%skin%surftype  = surftype_sea
       else
          profiles(i)%skin%surftype  = surftype_land
       endif
       ! sea-ice fraction
       if (seaice(i) >= 0.5) then
          profiles(i)%skin%surftype  = surftype_seaice
       endif

       ! dar: hard-coded to 1 (=ocean water) in rttov 9 int
       profiles(i)%skin%watertype = 1
       profiles(i) % idg        = 0
       profiles(i) % ice_scheme = 0

       IF(lallsky) THEN
         cld_profiles(i)%ph(:)   = ph_in(j,:)
         cld_profiles(i)%cc(:)   = cc_in(j,:)
         cld_profiles(i)%clw(:)  = clw_in(j,:)
         cld_profiles(i)%ciw(:)  = ciw_in(j,:)
         cld_profiles(i)%rain(:) = rain_in(j,:)
         cld_profiles(i)%sp(:)   = sp_in(j,:)
         profiles(i)%s2m%p       = cld_profiles(i)%ph(nlevels_in+1)
       ENDIF
       
     ENDDO
     ! end of 5.

    ich_temp = 1
    nchannels = nch_in

    CALL rttov_get_emis(       &
      errorstatus, &
      opts,        &
      chanprof,    &
      profiles,    &
      coef_rttov,  &
      emis_atlas(month),  &
      emissivity=emissivity(:)%emis_in)

    ! calculate emissivity for missing and ocean location (fastem)
    calcemis(:) = .FALSE.
    WHERE (emissivity(:)%emis_in <= 0.0)
      calcemis(:) = .TRUE.
    endwhere

    IF(lallsky) THEN
      CALL rttov_scatt ( &
        errorstatus,         &! out   error flag
        opts_scatt,          &! in    RTTOV-SCATT options structure
        nlevels,             &! in    number of profile levels
        chanprof,            &! in    channel and profile index structure
        frequencies,         &! in    channel indexes for Mietable lookup
        profiles,            &! in    profile array
        cld_profiles,        &! in    cloud/hydrometeor profile array
        coef_rttov,               &! in    coefficients structure
        coef_scatt,         &! in    Mietable structure
        calcemis,            &! in    flag for internal emissivity calcs
        emissivity,          &! inout input/output emissivities per channel
        radiance)             ! inout computed radiances
      CALL rttov_error('rttov_scatt error', lalloc = .TRUE.)
    ELSE
      CALL rttov_direct(         &
                        errorstatus,               &! out
                        chanprof,                  &
                        opts,                      &
                        profiles,                  &! in
                        coef_rttov,                &! in
                        transmission,              &! out
                        radiance,                  &
                        calcemis = calcemis,       &! in
                        emissivity = emissivity)    ! inout
      CALL rttov_error('rttov_direct error', lalloc = .TRUE.)
    ENDIF

    WRITE(999,*) radiance%bt(1:nchanprof)
    WRITE(998,*) prof_num, ich_temp, SIZE(ichan(:,no_id))

    tbs(1:prof_num , ich_temp:ich_temp + SIZE(ichan(:,no_id)) - 1) = &
      TRANSPOSE(RESHAPE(radiance%bt(1:nchanprof), (/ SIZE(ichan(:,no_id)), prof_num/) ))

    WRITE(997,*) tbs

    ich_temp = ich_temp + SIZE(ichan(:,no_id))

    ! --------------------------------------------------------------------------
    ! 8. deallocate all rttov arrays and structures
    ! --------------------------------------------------------------------------

    ! allocate structures for rttov_direct
    CALL rttov_alloc_direct(   &
      errorstatus,             &
      0_wp,                    &  ! 0 => deallocate
      nprof,                   &  ! nprofiles
      nchanprof,               &
      nlevels,                 &
      chanprof,                &
      opts,                    &
      profiles,                &
      coef_rttov,              &
      transmission,            &
      radiance,                &
      calcemis = calcemis,     &
      emissivity = emissivity)

  CONTAINS

    subroutine rttov_error(msg, lalloc)
      character(*) :: msg
      logical  :: lalloc

      if(lalloc) then
         if (any(alloc_status /= 0)) then
            write(*,*) msg
            errorstatus = 1
            call rttov_exit(errorstatus)
         endif
      else
         if (errorstatus /= errorstatus_success) then
            write(*,*) msg
            call rttov_exit(errorstatus)
         endif
      endif
    end subroutine rttov_error

  end subroutine rttov_multprof
  function construct_rttov_coeffilename(platform,satellite,instrument)
    ! Inputs
    integer,intent(in) :: platform,satellite,instrument
    ! Outputs
    character(len=256) :: construct_rttov_coeffilename
    ! Local variables
    character(len=256) :: coef_file
    integer :: error

    ! Initialize
    error = 0
    
    ! Platform
    if (platform .eq. 1)  coef_file = 'rtcoef_noaa_'
    if (platform .eq. 10) coef_file = 'rtcoef_metop_'
    if (platform .eq. 11) coef_file = 'rtcoef_envisat_'
    if (platform .ne. 1 .and. platform .ne. 10 .and. platform .ne. 11) then
       error=error+1
       write ( *,* ) 'Unsupported platform ID ',platform
       return
    endif

    ! Satellite
    if (satellite .lt. 10) then
       coef_file = trim(coef_file) // char(satellite+48)
    else if (satellite .lt. 100) then
       coef_file = trim(coef_file) // char(int(satellite/10)+48)
       coef_file = trim(coef_file) // char(satellite-int(satellite/10)*10+48)
    else
       error=error+1
       write ( *,* ) 'Unsupported satellite number ',satellite
       return
    endif

    ! Sensor
    if (sensor .eq. 3)  coef_file = trim(coef_file) // '_amsua.dat'
    if (sensor .eq. 5)  coef_file = trim(coef_file) // '_avhrr.dat'
    if (sensor .eq. 49) coef_file = trim(coef_file) // '_mwr.dat'
    if (sensor .ne. 3 .and. sensor .ne. 5 .and. sensor .ne. 49) then
       error=error+1
       write ( *,* ) 'Unsupported sensor number ', sensor
       return
    endif

    if (error .eq. 0) construct_rttov_coeffilename=coef_file
    
  end function construct_rttov_coeffilename
  function construct_rttov_scatfilename(platform,satellite,instrument)
    ! Inputs
    integer,intent(in) :: platform,satellite,instrument
    ! Outputs
    character(len=256) :: construct_rttov_scatfilename
    ! Local variables
    character(len=256) :: coef_file
    integer :: error

    ! Initialize
    error = 0
    
    ! Platform
    if (platform .eq. 1)  coef_file = 'sccldcoef_noaa_'
    if (platform .eq. 10) coef_file = 'sccldcoef_metop_'
    if (platform .eq. 11) coef_file = 'sccldcoef_envisat_'
    if (platform .ne. 1 .and. platform .ne. 10 .and. platform .ne. 11) then
       error=error+1
       write ( *,* ) 'Unsupported platform ID ',platform
       return
    endif

    ! Satellite
    if (satellite .lt. 10) then
       coef_file = trim(coef_file) // char(satellite+48)
    else if (satellite .lt. 100) then
       coef_file = trim(coef_file) // char(int(satellite/10)+48)
       coef_file = trim(coef_file) // char(satellite-int(satellite/10)*10+48)
    else
       error=error+1
       write ( *,* ) 'Unsupported satellite number ',satellite
       return
    endif

    ! Sensor
    if (sensor .eq. 3)  coef_file = trim(coef_file) // '_amsua.dat'
    if (sensor .eq. 5)  coef_file = trim(coef_file) // '_avhrr.dat'
    if (sensor .eq. 49) coef_file = trim(coef_file) // '_mwr.dat'
    if (sensor .ne. 3 .and. sensor .ne. 5 .and. sensor .ne. 49) then
       error=error+1
       write ( *,* ) 'Unsupported sensor number ', sensor
       return
    endif

    if (error .eq. 0) construct_rttov_scatfilename=coef_file
    
  end function construct_rttov_scatfilename
  
end module mod_cosp_rttov
