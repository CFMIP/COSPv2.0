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
! JKS fill this in when working :)

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module mod_cosp_rttov
  use rttov_const,         only : errorstatus_success, errorstatus_fatal,                &
                                  platform_name,inst_name
! Use RTTOV v13 types here, more may be needed
  use rttov_types,         only : rttov_options,rttov_coefs,rttov_profile,                &
                                  rttov_transmission,rttov_radiance,rttov_chanprof,       &
                                  rttov_emissivity,rttov_reflectance,rttov_opt_param
!  use rttov_types,         only : rttov_options,rttov_coefs,profile_type,                &
!                                  transmission_type,radiance_type,rttov_chanprof,        &
!                                  rttov_emissivity,profile_cloud_type,rttov_scatt_coef,  &
!                                  rttov_options_scatt
  use rttov_const,         only : surftype_sea, surftype_land, surftype_seaice,           &
                                  errorstatus_success,errorstatus_fatal,                  &
                                  platform_name,inst_name,                                &
                                  sensor_id_mw,sensor_id_po
  use rttov_unix_env,      only : rttov_exit
  use cosp_kinds,          only : wp
  use mod_cosp_config,     only : RTTOV_MAX_CHANNELS,N_HYDRO,rttovDir
  use cosp_phys_constants, only : mdry=>amd,mO3=>amO3,mco2=>amCO2,mCH4=>amCH4,           &
                                  mn2o=>amN2O,mco=>amCO

  ! The rttov_emis_atlas_data type must be imported separately
  use mod_rttov_emis_atlas, ONLY : &
        rttov_emis_atlas_data, &
        atlas_type_ir, atlas_type_mw

  ! The rttov_brdf_atlas_data type must be imported separately
  use mod_rttov_brdf_atlas, ONLY : rttov_brdf_atlas_data

  ! jpim, jprb and jplm are the RTTOV integer, real and logical KINDs
  use parkind1, ONLY : jpim, jprb, jplm

  use rttov_unix_env, ONLY : rttov_exit
  
  implicit none


! Old includes
!#include "rttov_direct.interface"
!#include "rttov_alloc_prof.interface"
!#include "rttov_alloc_rad.interface"
!#include "rttov_alloc_transmission.interface"
!#include "rttov_dealloc_coefs.interface"
!#include "rttov_user_options_checkinput.interface"
!#include "rttov_read_coefs.interface"
!#include "rttov_get_emis.interface"
!#include "rttov_boundaryconditions.interface"

! New includes for v13 (will need to clean up others)
#include "rttov_direct.interface"
#include "rttov_parallel_direct.interface"
#include "rttov_read_coefs.interface"
#include "rttov_read_scattcoeffs.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_alloc_direct.interface"
#include "rttov_init_emis_refl.interface"
#include "rttov_user_options_checkinput.interface"
#include "rttov_print_opts.interface"
#include "rttov_print_profile.interface"
#include "rttov_skipcommentline.interface"

! Includes when directly inputting cloud optical parameters
#include "rttov_init_opt_param.interface"
#include "rttov_bpr_init.interface"
#include "rttov_bpr_calc.interface"
#include "rttov_bpr_dealloc.interface"
#include "rttov_legcoef_calc.interface"


  ! Module parameters
  integer, parameter :: maxlim =  10000
  real(wp),parameter :: eps    =  0.622

  ! Initialization parameters
  integer :: &
       platform,   & ! RTTOV platform
       sensor,     & ! RTTOV instrument
       satellite,  & ! RTTOV satellite
       nChannels     ! Number of channels
  integer,dimension(RTTOV_MAX_CHANNELS) :: &
       iChannel      ! RTTOV channel numbers

  ! Scattering coefficients (read in once during initialization)
  type(rttov_coefs) :: &
       coef_rttov
! JKS - KISS
!  type(rttov_scatt_coef) :: &
!       coef_scatt
  ! RTTOV setup and options (set during initialization)
  type(rttov_options) :: &
       opts     ! defaults to everything optional switched off
  logical             :: &
      do_rttov_cld,        & ! Include clouds in RTTOV calculations
      do_rttov_aer,        & ! Include aerosols in RTTOV calculations
      do_rttov_rad,        & ! Return radiances in RTTOV calculations
      rttov_cld_optparam,  & ! Use user-supplied optical cloud parameters
      rttov_aer_optparam     ! Use user-supplied optical aerosol parameters
       
  ! JKS should this be module-wide?
!  type(rttov_IN)      :: &
!       rttovIN

! JKS - KISS
!  type(rttov_options_scatt) :: &
!       opts_scatt

  ! module-wides variables for input
  !====================
  integer(kind=jpim) :: nthreads
  integer(kind=jpim) :: dosolar
  integer(kind=jpim) :: nchanprof ! JKS - jpim is RTTOV integer type
!  integer(kind=jpim), allocatable :: channel_list(:) ! JKS this needs to be specified
  
  TYPE(rttov_chanprof),    POINTER :: chanprof(:)    => NULL() ! Input channel/profile list
  LOGICAL(KIND=jplm),      POINTER :: calcemis(:)    => NULL() ! Flag to indicate calculation of emissivity within RTTOV
  TYPE(rttov_emissivity),  POINTER :: emissivity(:)  => NULL() ! Input/output surface emissivity
  LOGICAL(KIND=jplm),      POINTER :: calcrefl(:)    => NULL() ! Flag to indicate calculation of BRDF within RTTOV
  TYPE(rttov_reflectance), POINTER :: reflectance(:) => NULL() ! Input/output surface BRDF
  TYPE(rttov_profile),     POINTER :: profiles(:)    => NULL() ! Input profiles
  TYPE(rttov_transmission)         :: transmission             ! Output transmittances
  TYPE(rttov_radiance)             :: radiance                 ! Output radiances
  
  INTEGER(KIND=jpim)               :: errorstatus              ! Return error status of RTTOV subroutine calls

  INTEGER(KIND=jpim) :: alloc_status(60)

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! TYPE rttov_in
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  ! JKS - add additional COSP inputs here.
  type rttov_IN
     integer(kind=jpim),pointer :: & ! JKS trying this
          nPoints,      & ! Number of profiles to simulate
          nLevels,      & ! Number of levels
          nSubCols,     & ! Number of subcolumns
          nChannels,    & ! Number of channels to simulate
          month           ! Month (needed for surface emissivity calculation)
     real(wp),pointer :: &
          zenang,       & ! Satellite zenith angle
          co2,          & ! Carbon dioxide 
          ch4,          & ! Methane 
          n2o,          & ! n2o 
          co              ! Carbon monoxide
     real(wp),dimension(:),pointer :: &
          surfem           ! Surface emissivities for the channels
!          refl,         & ! Surface reflectances for the channels
     integer(kind=jpim),dimension(:),pointer :: &
          channels        ! Surface reflectances for the channels
     real(wp),dimension(:),pointer :: &
          h_surf,       & ! Surface height
          u_surf,       & ! U component of surface wind
          v_surf,       & ! V component of surface wind
          t_skin,       & ! Surface skin temperature
          p_surf,       & ! Surface pressure
          t2m,          & ! 2 m Temperature
          q2m,          & ! 2 m Specific humidity
          lsmask,       & ! land-sea mask
          latitude,     & ! Latitude
          longitude,    & ! Longitude
          seaice          ! Sea-ice? 
     real(wp),dimension(:,:),pointer :: &
          p,            & ! Pressure @ model levels
          ph,           & ! Pressure @ model half levels
          t,            & ! Temperature 
          q,            & ! Specific humidity
          o3              ! Ozone
     
     ! These fields below are needed ONLY for the RTTOV all-sky brightness temperature
     real(wp),dimension(:,:),pointer :: &
          tca,          & ! Cloud fraction
          cldIce,       & ! Cloud ice
          cldLiq,       & ! Cloud liquid
          fl_rain,      & ! Precipitation flux (startiform+convective rain) (kg/m2/s)
          fl_snow         ! Precipitation flux (stratiform+convective snow)
  end type rttov_IN

contains

  ! Wrapper function for exiting RTTOV and reporting the error
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


  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE rttov_allocate - JKS
  ! ------------------------------------------------------
  ! 3. Allocate RTTOV input and output structures
  ! 4. Build the list of profile/channel indices in chanprof
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine cosp_rttov_allocate(rttovIN &
                           )
                           
    type(rttov_in),intent(in) :: & ! What is the best way to do this? Should rttovIN be a module-wide DDT? Yes.
        rttovIN
        
    ! Loop variables
    integer(kind=jpim) :: j, jch, nch
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 3. Allocate RTTOV input and output structures
    ! ------------------------------------------------------
    ! Largely from RTTOV documentation.
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ! Determine the total number of radiances to simulate (nchanprof).
    ! We aren't doing subcolumn sampling (RTTOV already does this and it would be slow)
!      nchanprof = nchannels * nprof

    ! Allocate and fill in channel_list
!    allocate(channel_list(rttovIN%nChannels)) 
!    channel_list(1:nChannels) = rttovIN%channels(1:nChannels)
    
!    print*,'rttovIN%channels: ',rttovIN%channels
!    print*,'channel_list:  ',channel_list
    
    nchanprof = rttovIN%nChannels * rttovIN%nPoints

    ! Allocate structures for rttov_direct
    call rttov_alloc_direct( &
        errorstatus,             &
        1_jpim,                  &  ! 1 => allocate
!       nprof,                   &
        rttovIN%nPoints,         &
        nchanprof,               &
!       nlevels,                 &
        rttovIN%nLevels,         &
        chanprof,                &
        opts,                    &
        profiles,                &
        coef_rttov,              &
!       coefs,                   &
        transmission,            &
        radiance,                &
        calcemis=calcemis,       &
        emissivity=emissivity,   &
        calcrefl=calcrefl,       &
        reflectance=reflectance, &
        init=.TRUE._jplm)
    call rttov_error('allocation error for rttov_direct structures' , lalloc = .false.)

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 4. Build the list of profile/channel indices in chanprof
    ! ------------------------------------------------------
    ! Largely from RTTOV documentation.
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    nch = 0_jpim
    do j = 1, rttovIN%nPoints
      do jch = 1, rttovIN%nChannels
        nch = nch + 1_jpim
        chanprof(nch)%prof = j
        chanprof(nch)%chan = rttovIN%channels(jch) ! Example code used channel_list
      end do
    end do
        
  end subroutine cosp_rttov_allocate
  
  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! 5. rttov_construct_profiles: 5. Read profile data
  ! ------------------------------------------------------
  ! Largely from cosp_rttov_v11.F90 file.
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine cosp_rttov_construct_profiles(rttovIN &
                           )

    type(rttov_in),intent(in) :: & ! What is the best way to do this? Should rttovIN be a module-wide DDT? Yes.
        rttovIN
        
    ! Loop variables
    integer(kind=jpim) :: i, j ! Use i to iterate over profile, j for levels.
    
    logical                             :: &
        rttov_simulate_cld,                &
        rttov_simulate_aer
        
    ! Store profile data from rttovIN in profile type.
    ! See RTTOV user guide pg 163 for description of "profiles" type
    
    ! "The rttov_profile structure is composed of the atmospheric part 
    ! and two other structures for 2 meter air and skin surface. 
    ! If you are not able to provide ozone, CO2, etc profiles the flags
    ! ozone_data, co2_data and so on in the options structure should be 
    ! set to false."

    profiles(:)%gas_units  =  1 ! kg/kg over moist air (default)
    
    do i = 1, rttovIN%nPoints
        
      ! Initialize trace gas concentrations
      profiles(i)%co2(:)        =  rttovIN%co2
      profiles(i)%n2o(:)        =  rttovIN%n2o
      profiles(i)%co(:)         =  rttovIN%n2o
      profiles(i)%ch4(:)        =  rttovIN%ch4
      
      ! Initialize column pressure, temperature, and humidity
      profiles(i)%p(:) =  rttovIN%p(i, :) * 1e-2 ! convert Pa to hPa
      profiles(i)%t(:) =  rttovIN%t(i, :)
      profiles(i)%q(:) =  rttovIN%q(i, :)

      ! q coefficient limit is 09=.1e-10
      where(profiles(i)%q(:) < 0.1e-10)
         profiles(i)%q(:) = 0.11e-10
      end where
      
!      if (any(profiles%q < 0.1e-10)) then
!        write(*,*) 'q profile less than RTTOV minimum. Rewritten to 1e-11'
!        errorstatus = 1
!        call rttov_exit(errorstatus)
!      endif

      ! Gas profiles
      profiles(i)%o3         =  rttovIN%o3(i, :)
!      profiles(i)%so2        =  ! Sulfate not in COSP input files
       
      ! 2m parameters
      profiles(i)%s2m%p      =  rttovIN%p_surf(i) * 1e-2 ! convert Pa to hPa
      profiles(i)%s2m%t      =  rttovIN%t2m(i) ! JKS or rttovIN%t_skin
      profiles(i)%s2m%q      =  rttovIN%q2m(i) ! Should be the same as gas units (kg/kg)
      profiles(i)%s2m%u      =  rttovIN%u_surf(i)
      profiles(i)%s2m%v      =  rttovIN%v_surf(i)
      profiles(i)%s2m%wfetc  =  10000. ! only used by sea surface solar BRDF model.

      ! skin variables for emissivity calculations
      profiles(i)%skin%t          =  rttovIN%t_skin(i)

      ! fastem coefficients - for mw calculations
      profiles(i)%skin%fastem(1)  =  3.0
      profiles(i)%skin%fastem(2)  =  5.0
      profiles(i)%skin%fastem(3)  =  15.0
      profiles(i)%skin%fastem(4)  =  0.1
      profiles(i)%skin%fastem(5)  =  0.3

      ! Viewing angles
      profiles(i)%zenangle      = rttovIN%zenang ! pass in from cosp
      profiles(i)%azangle       = 0. ! hard-coded in rttov9 int JKS-?

      profiles(i)%latitude      = rttovIN%latitude(i)
      profiles(i)%longitude     = rttovIN%longitude(i)
      profiles(i)%elevation     = rttovIN%h_surf(i) * 1e-3 ! Convert m to km

      ! Solar angles. JKS - get this from COSP? Doesn't seem to be passed in.
      profiles(i)%sunzenangle   = 0. ! hard-coded in rttov9 int
      profiles(i)%sunazangle    = 0. ! hard-coded in rttov9 int

      ! surface type
      ! land-sea mask (lsmask) indicates proportion of land in grid
      if (rttovIN%lsmask(i) < 0.5) then
        profiles(i)%skin%surftype  = surftype_sea
      else
        profiles(i)%skin%surftype  = surftype_land
      endif
      ! sea-ice fraction
      if (rttovIN%seaice(i) >= 0.5) then
        profiles(i)%skin%surftype  = surftype_seaice
      endif

      ! dar: hard-coded to 1 (=ocean water) in rttov 9 int
      profiles(i)%skin%watertype = 1
      !profiles(i) %idg         = 0. ! Depreciated?
      !profiles(i) %ish         = 0. ! Depreciated?
    end do
        
    ! JKS - nothing to check here, this will never trigger.
    call rttov_error('error in profile initialization' , lalloc = .false.)
    
!    print*,'profiles(90)%p(:):  ',profiles(90)%p(:)
!    print*,'profiles(90)%q(:):  ',profiles(90)%q(:)
!    print*,'profiles(90)%t(:):  ',profiles(90)%t(:)
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Only add the cloud fields if simulating cloud.
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    rttov_simulate_cld = .false.
!    if (rttov_simulate_cld) then
    
    if (do_rttov_cld) then

      ! Set cloud mass mixing ratio units
      profiles(:)%mmr_cldaer =  .true. ! kg/kg for cloud and aerosol (default)

      profiles(:)%clw_scheme   = 2 ! Deff scheme avoids cloud types
    !    profiles%clwde_scheme = 1. ! Not implemented?
      profiles(:)%ice_scheme   = 1 !1:Baum 2:Baran(2014) 3:Baran(2018)
      profiles(:)%icede_param  = 2 ! 2:Wyser(recommended). Only used if ice effective diameter not input
        
      do i = 1, rttovIN%nPoints        
        ! Cloud scheme stuff
        profiles(i)%cfrac(:)   = rttovIN%tca(i,:)         ! Cloud fraction for each layer       
        profiles(i)%cloud(1,:) = rttovIN%cldLiq(i,:) ! Cloud water mixing ratio (all in the first type for Deff)
        profiles(i)%cloud(6,:) = rttovIN%cldIce(i,:) ! Cloud ice mixing ratio (1 type). See pg 74.

        ! Example UKMO input has effective radii for multiple cloud types, making identification of a single
        ! liquid droplet or ice crystal effective diameter difficult.
        ! I opt to let RTTOV decide on the effective radius values, but more complex implementation
        ! could do a more thorough conversion between UKMO output and RTTOV input
    !    profiles(i)%clwde = ! Cloud water effective diameter
    !    profiles(i)%icede = ! Cloud ice effective diameter

        ! Old code for simple cloud schemes only
    !    profiles(i)%cfraction  =  0.
    !    profiles(i)%ctp        =  500.

        ! Other options not implemented
!        profiles(i)%clw        = ! Cloud liquid water (kg/kg) – MW only,
      end do
    endif
    
    ! JKS - nothing to check here, this will never trigger.
    call rttov_error('error in cloud profile initialization' , lalloc = .false.)
    
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Only add the aerosol fields if simulating aerosol.
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    rttov_simulate_aer = .false.
!    if (rttov_simulate_aer) then
    
    if (do_rttov_aer) then
    
      ! Set aerosol mass mixing ratio units
      profiles%mmr_cldaer =  .true. ! kg/kg for cloud and aerosol (default)
        
      ! Read in aerosol profiles
!      do i = 1, rttovIN%nPoints
!      profiles(i)%aerosols(naertyp,nlayers) = ! Aerosols in different modes (see User Guide pg 80)
!
!      end do
    endif
    
    ! JKS - nothing to check here, this will never trigger.
    call rttov_error('error in aerosol profile initialization' , lalloc = .true.)
    
    ! JKS To-do: set up scattering profiles (MW only) (rttov_profile_cloud)

  end subroutine cosp_rttov_construct_profiles
  

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! 6. rttov_setup_emissivity_reflectance - Specify surface emissivity and reflectance
  ! ------------------------------------------------------
  ! From RTTOV example files. Will need to be expanded on to pass in values.
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine cosp_rttov_setup_emissivity_reflectance()
    
    ! In this example we have no values for input emissivities or reflectances
    ! so we initialise all inputs to zero
    call rttov_init_emis_refl(emissivity, reflectance)
    call rttov_error('error for emissivity/reflectance initialization' , lalloc = .true.)

    ! Calculate emissivity within RTTOV where the input emissivity value is
    ! zero or less (all channels in this case)
    calcemis(:) = (emissivity(:) % emis_in <= 0._jprb)
    
    ! Calculate reflectances within RTTOV where the input BRDF value is zero or
    ! less (all channels in this case)
    calcrefl(:) = (reflectance(:) % refl_in <= 0._jprb)  
  
  end subroutine cosp_rttov_setup_emissivity_reflectance
  
  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! 7. rttov_call_direct - Call RTTOV forward model (Woohoo!)
  ! ------------------------------------------------------
  ! From RTTOV example files.
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  subroutine cosp_rttov_call_direct(nthreads)
    
    integer(kind=jpim) :: nthreads

    nthreads = 1 ! Default not parallel for now. Can be optimized later. - JKS
    print*,'Calling rttov_direct'
    
    if (nthreads <= 1) then
      call rttov_direct(                &
              errorstatus,              &! out   error flag
              chanprof,                 &! in    channel and profile index structure
              opts,                     &! in    options structure
              profiles,                 &! in    profile array
              coef_rttov,               &! in    coefficients structure
              transmission,             &! inout computed transmittances
              radiance,                 &! inout computed radiances
              calcemis    = calcemis,   &! in    flag for internal emissivity calcs
              emissivity  = emissivity, &! inout input/output emissivities per channel
              calcrefl    = calcrefl,   &! in    flag for internal BRDF calcs
              reflectance = reflectance) ! inout input/output BRDFs per channel
    else
      call rttov_parallel_direct(       &
              errorstatus,              &! out   error flag
              chanprof,                 &! in    channel and profile index structure
              opts,                     &! in    options structure
              profiles,                 &! in    profile array
              coef_rttov,               &! in    coefficients structure
              transmission,             &! inout computed transmittances
              radiance,                 &! inout computed radiances
              calcemis    = calcemis,   &! in    flag for internal emissivity calcs
              emissivity  = emissivity, &! inout input/output emissivities per channel
              calcrefl    = calcrefl,   &! in    flag for internal BRDF calcs
              reflectance = reflectance,&! inout input/output BRDFs per channel
              nthreads    = nthreads)    ! in    number of threads to use
    endif
    call rttov_error('rttov_direct error', lalloc = .true.)
  
  end subroutine cosp_rttov_call_direct
  
  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! rttov_save_and_deallocate - 8. Save output data, 9. Deallocate all RTTOV arrays and structures
  ! ------------------------------------------------------
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine cosp_rttov_save_and_deallocate_profiles(rttovIN,Tb)

    type(rttov_in),intent(in) :: &
        rttovIN
    real(wp),dimension(rttovIN%nPoints,rttovIN%nChannels),intent(inout) :: & ! Can I do this? I guess so!
        Tb        ! RTTOV brightness temperature.
        
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 8. Save output data
    ! ------------------------------------------------------
    ! JKS - Need to allow options for Tb and radiance for clear- and cloudy-skies
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Tb(1:rttovIN%nPoints, 1:rttovIN%nChannels) = &
        transpose(reshape(radiance%bt(1:nchanprof), (/ rttovIN%nChannels, rttovIN%nPoints/) ))    
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 9. Deallocate all RTTOV arrays and structures
    ! ------------------------------------------------------
    ! From RTTOV example files.
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ! JKS no longer using channel_list
    !deallocate (channel_list, stat=alloc_status(1))
    !if (alloc_status(1) /= 0) then
    !  write(*,*) 'mem dellocation error'
    !endif

    ! Deallocate structures for rttov_direct
    call rttov_alloc_direct(     &
        errorstatus,             &
        0_jpim,                  &  ! 0 => deallocate
!       nprof,                   &
        rttovIN%nPoints,         &
        nchanprof,               &
!       nlevels,                 &
        rttovIN%nLevels,         &
        chanprof,                &
        opts,                    &
        profiles,                &
        coef_rttov,              &
        transmission,            &
        radiance,                &
        calcemis=calcemis,       &
        emissivity=emissivity,   &
        calcrefl=calcrefl,       &
        reflectance=reflectance)
    call rttov_error('deallocation error for rttov_direct structures', lalloc = .true.)
        
    call rttov_dealloc_coefs(errorstatus, coef_rttov)
    call rttov_error('coefs deallocation error', lalloc = .true.)

  end subroutine cosp_rttov_save_and_deallocate_profiles
  
  subroutine cosp_rttov_deallocate_coefs()
  
    call rttov_dealloc_coefs(errorstatus, coef_rttov)
    if (errorstatus /= errorstatus_success) then
      write(*,*) 'coefs deallocation error'
    endif

  end subroutine cosp_rttov_deallocate_coefs


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
