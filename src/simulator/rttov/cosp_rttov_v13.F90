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
                                  rttov_emissivity,rttov_reflectance,rttov_opt_param,     &
                                  rttov_pccomp
  use rttov_const,         only : surftype_sea, surftype_land, surftype_seaice,           &
                                  errorstatus_success,errorstatus_fatal,                  &
                                  platform_name,inst_name,                                &
                                  sensor_id_mw,sensor_id_po
  use rttov_unix_env,      only : rttov_exit
  use cosp_kinds,          only : wp
  use mod_cosp_config,     only : N_HYDRO
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


! New includes for v13 (will need to clean up others)
#include "rttov_direct.interface"
#include "rttov_parallel_direct.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_alloc_direct.interface"
#include "rttov_init_emis_refl.interface"
#include "rttov_user_options_checkinput.interface"
#include "rttov_print_opts.interface"
#include "rttov_print_profile.interface"

! checking inputs
#include "rttov_dealloc_coef_scatt.interface"
#include "rttov_dealloc_coef.interface"
#include "rttov_dealloc_coef_pccomp.interface"


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
  integer(kind=jpim) :: &
       nChannels     ! Number of channels
  integer(kind=jpim) ::  & ! Parallelization is default off
       rttov_direct_nthreads = 1_jpim
  integer(kind=jpim),allocatable :: &
       iChannel(:),      &  ! RTTOV channel indices
       emisChannel(:),   &  ! RTTOV channel emissivity
       reflChannel(:)       ! RTTOV channel reflectivity

  ! Scattering coefficients (read in once during initialization)
! JKS - KISS
!  type(rttov_scatt_coef) :: &
!       coef_scatt
  ! RTTOV setup and options (set during initialization)
  logical             :: &
      do_rttov_cld,        & ! Include clouds in RTTOV calculations
      do_rttov_aer,        & ! Include aerosols in RTTOV calculations
      do_rttov_bt,         & ! Return brightness temps in RTTOV calculations
      do_rttov_rad,        & ! Return radiances in RTTOV calculations
      do_rttov_refl,       & ! Return reflectances in RTTOV calculations
      do_rttov_pcrttov,    & ! Do computations using PC-RTTOV
      rttov_cld_optparam,  & ! Use user-supplied optical cloud parameters
      rttov_aer_optparam     ! Use user-supplied optical aerosol parameters
      
  ! --- Well-mixed trace gas mixing ratios from user via RTTOV namelist
  real(wp)  :: so2 = 0._wp
  real(wp)  :: ch4 = 0._wp
  real(wp)  :: co  = 0._wp
  real(wp)  :: co2 = 0._wp
  real(wp)  :: n2o = 0._wp
  real(wp)  :: zenang = 0._wp ! NADIR default
       
  character(len=256) :: &
      rttovDir                    ! Directory for the RTTOV source code
      
  character(len=256) :: PC_coef_filepath = ''
      
! JKS - KISS
!  type(rttov_options_scatt) :: &
!       opts_scatt

  ! module-wides variables for input
  !====================
  integer(kind=jpim) :: dosolar
  integer(kind=jpim) :: nchanprof ! JKS - jpim is RTTOV integer type

  TYPE(rttov_options)              :: opts                     ! defaults to everything optional switched off
  TYPE(rttov_coefs)                :: coef_rttov               ! Coefficients structure
  TYPE(rttov_chanprof),    POINTER :: chanprof(:)    => NULL() ! Input channel/profile list
  LOGICAL(KIND=jplm),      POINTER :: calcemis(:)    => NULL() ! Flag to indicate calculation of emissivity within RTTOV
  TYPE(rttov_emissivity),  POINTER :: emissivity(:)  => NULL() ! Input/output surface emissivity
  LOGICAL(KIND=jplm),      POINTER :: calcrefl(:)    => NULL() ! Flag to indicate calculation of BRDF within RTTOV
  TYPE(rttov_reflectance), POINTER :: reflectance(:) => NULL() ! Input/output surface BRDF
  TYPE(rttov_profile),     POINTER :: profiles(:)    => NULL() ! Input profiles
  TYPE(rttov_transmission)         :: transmission             ! Output transmittances
  TYPE(rttov_radiance)             :: radiance                 ! Output radiances
  TYPE(rttov_pccomp)               :: pccomp                    ! Output PC structure
  INTEGER(KIND=jpim),      POINTER :: channels_rec(:) => NULL() ! Reconstructed radiance channel list
  
  INTEGER(KIND=jpim)               :: errorstatus              ! Return error status of RTTOV subroutine calls

  INTEGER(KIND=jpim) :: alloc_status(60)
  
  ! JKS additional variables used in PC-RTTOV
  INTEGER(KIND=jpim), POINTER :: predictindex(:)
  INTEGER(KIND=jpim) :: nchannels_rec, nchannels_comp, npcscores
  INTEGER(KIND=jpim) :: lo, hi

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! TYPE rttov_in
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  ! JKS - add additional COSP inputs here.
  type rttov_IN
     integer(kind=jpim),pointer :: & ! JKS trying this
          nPoints,      & ! Number of profiles to simulate
          nLevels,      & ! Number of levels
          nSubCols,     & ! Number of subcolumns
          nChannels,    & ! Number of channels to simulate ! JKS
          month           ! Month (needed for surface emissivity calculation)
     real(wp),pointer :: & ! Could change the dimensionality of these in the future
          co2,          & ! Carbon dioxide 
          ch4,          & ! Methane 
          n2o,          & ! n2o 
          co              ! Carbon monoxide
!     real(wp),dimension(:),pointer :: &
!          surfem           ! Surface emissivities for the channels
!          refl,         & ! Surface reflectances for the channels
!     integer(kind=jpim),dimension(:),pointer :: &
!          channels        ! 
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
                           
    type(rttov_in),intent(in) :: &
        rttovIN
        
    ! Loop variables
    integer(kind=jpim) :: j, jch, nch
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 3. Allocate RTTOV input and output structures
    ! ------------------------------------------------------
    ! Largely from RTTOV documentation.
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ! Determine the total number of radiances to simulate (nchanprof).    
    nchanprof = rttovIN%nChannels * rttovIN%nPoints

    ! Allocate structures for rttov_direct
    call rttov_alloc_direct( &
        errorstatus,             &
        1_jpim,                  &  ! 1 => allocate
        rttovIN%nPoints,         &
        nchanprof,               &
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
      do jch = 1, rttovIN%nChannels ! nChannels
        nch = nch + 1_jpim
        chanprof(nch)%prof = j
        chanprof(nch)%chan = iChannel(jch) ! Example code used channel_list
      end do
    end do
        
  end subroutine cosp_rttov_allocate
  

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE pc_rttov_allocate - Subroutine for running PC-RTTOV
  ! ------------------------------------------------------
  ! 3. Allocate RTTOV input and output structures
  ! 4. Build the list of profile/channel indices in chanprof
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine cosp_pc_rttov_allocate(rttovIN &
                                   )
                           
    type(rttov_in),intent(in) :: &
        rttovIN
        
    ! Loop variables
    integer(kind=jpim) :: j, jch, nch
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 3. Allocate RTTOV input and output structures
    ! ------------------------------------------------------
    ! Largely from RTTOV documentation.
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
!    IF (errorstatus /= errorstatus_success) THEN
!      WRITE(*,*) 'rttov_get_pc_predictindex fatal error'
!      CALL rttov_exit(errorstatus)
!    ENDIF

! JKS clean this up when done.
!    nchannels_rec = rttovIN%nChannels             ! Number of channels to reconstruct indicated in namelist.
    nchannels_rec  = SIZE(predictindex)           ! Number of channels to reconstruct indicated in namelist.
    nchannels_comp = SIZE(predictindex)           ! Number of channels to compute is determined by ipcreg in init.
    nchanprof = nchannels_comp * rttovIN%nPoints  ! Size of chanprof array is total number of channels over all profiles.

    ! Determine the number of reconstructed radiances per profile (nchannels_rec)
    if (opts % rt_ir % pc % addradrec) then
      if (nchannels_rec < 0) then
        ! If the number of channels is negative, don't reconstruct radiances at all
        opts % rt_ir % pc % addradrec = .FALSE.
      else if (nchannels_rec == 0) then
        ! If the number of channels is set to 0 then reconstruct all instrument channels
        nchannels_rec = coef_rttov % coef % fmv_chn
        allocate(channels_rec(nchannels_rec))
        channels_rec(:) = (/ (j, j = 1, nchannels_rec) /)
      else
        ! Otherwise read the channel list from the file
        allocate(channels_rec(nchannels_rec))
        channels_rec(:) = iChannel ! channels_rec is just the index of the desired channels
      endif
    endif

    ! Ensure we don't have unassociated pointers below when addradrec is FALSE
    if (nchannels_rec <= 0) allocate(channels_rec(0))

    ! Allocate structures for rttov_direct
    CALL rttov_alloc_direct(                             &
          errorstatus,                                   &
          1_jpim,                                        &  ! 1 => allocate
          rttovIN%nPoints,                               &
          nchanprof,                                     &
          rttovIN%nLevels,                               &
          chanprof,                                      &
          opts,                                          &
          profiles,                                      &
          coef_rttov,                                    &
          transmission,                                  &
          radiance,                                      &
          calcemis=calcemis,                             &
          emissivity=emissivity,                         &
          npcscores=npcscores * rttovIN%nPoints,         &
          nchannels_rec=nchannels_rec * rttovIN%nPoints, &
          pccomp=pccomp,                                 &
          init=.TRUE._jplm)
    call rttov_error('allocation error for rttov_direct structures (PC-RTTOV)' , lalloc = .true.)
          
!    IF (errorstatus /= errorstatus_success) THEN
!      WRITE(*,*) 'allocation error for rttov_direct structures'
!      CALL rttov_exit(errorstatus)
!    ENDIF

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 4. Build the list of profile/channel indices in chanprof
    ! ------------------------------------------------------
    ! Largely from RTTOV documentation.
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Populate chanprof using the channel list obtained above in predictindex(:)
    do j = 1, rttovIN%nPoints
      lo = (j - 1) * nchannels_comp + 1
      hi = lo + nchannels - 1
      chanprof(lo:hi)%prof = j
      chanprof(lo:hi)%chan = predictindex(:)
    end do
        
  end subroutine cosp_pc_rttov_allocate
  

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
        
    ! Store profile data from rttovIN in profile type.
    ! See RTTOV user guide pg 163 for description of "profiles" type
    
    ! "The rttov_profile structure is composed of the atmospheric part 
    ! and two other structures for 2 meter air and skin surface. 
    ! If you are not able to provide ozone, CO2, etc profiles the flags
    ! ozone_data, co2_data and so on in the options structure should be 
    ! set to false."

!    print*,'co2: ',co2
!    print*,'n2o: ',n2o
!    print*,'co:  ',co
!    print*,'ch4: ',ch4
!    print*,'so2: ',so2
!    print*,'zenang: ',zenang

    profiles(:)%gas_units  =  1 ! kg/kg over moist air (default)
    
    do i = 1, rttovIN%nPoints
        
      ! Initialize trace gas concentrations from user input
      ! These gases are not in COSP input files but might be in the futre
      profiles(i)%co2(:)        = co2
      profiles(i)%n2o(:)        = n2o
      profiles(i)%co(:)         = co
      profiles(i)%ch4(:)        = ch4
      profiles(i)%so2           = so2
      
!      profiles(i)%co2(:)        =  rttovIN%co2
!      profiles(i)%n2o(:)        =  rttovIN%n2o
!      profiles(i)%co(:)         =  rttovIN%co
!      profiles(i)%ch4(:)        =  rttovIN%ch4
      
      ! Initialize column pressure, temperature, and humidity
      profiles(i)%p(:) =  rttovIN%p(i, :) * 1e-2 ! convert Pa to hPa
      profiles(i)%t(:) =  rttovIN%t(i, :)
      profiles(i)%q(:) =  rttovIN%q(i, :)

      ! q coefficient limit is 0.1e-10
      where(profiles(i)%q(:) < 0.1e-10)
         profiles(i)%q(:) = 0.11e-10
      end where
      
      ! Gas profiles
      profiles(i)%o3         =  rttovIN%o3(i, :)
       
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
      profiles(i)%zenangle      = zenang ! pass in from cosp
      profiles(i)%azangle       = 0. ! hard-coded in rttov9 int JKS-?

      profiles(i)%latitude      = rttovIN%latitude(i)
      profiles(i)%longitude     = rttovIN%longitude(i)
      profiles(i)%elevation     = rttovIN%h_surf(i) * 1e-3 ! Convert m to km

      ! Solar angles. JKS - get this from COSP/CESM? Doesn't seem to be passed in.
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
    
!    print*,'profiles(1)%p(:):  ',profiles(1)%p(:)
!    print*,'profiles(1)%q(:):  ',profiles(1)%q(:)
!    print*,'profiles(1)%t(:):  ',profiles(1)%t(:)
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Only add the cloud fields if simulating cloud.
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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


  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! 6. rttov_setup_emissivity_reflectance - Specify surface emissivity for PC-RTTOV
  ! ------------------------------------------------------
  ! From RTTOV example files.
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine cosp_pc_rttov_setup_emissivity()
    
    ! PC-RTTOV requires using RTTOV to calculate the surface emissivities.
    ! Reflectances are never calculated for hyper-spectral IR sounders
    call rttov_init_emis_refl(emissivity)
    calcemis(:) = .TRUE.
    call rttov_error('error for emissivity initialization (PC-RTTOV)' , lalloc = .true.)
  
  end subroutine cosp_pc_rttov_setup_emissivity

  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! 7. rttov_call_direct - Call RTTOV forward model (Woohoo!)
  ! ------------------------------------------------------
  ! From RTTOV example files.
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  subroutine cosp_rttov_call_direct()
    
    if (rttov_direct_nthreads <= 1) then
      print*,'Calling rttov_direct'
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
      print*,'Calling rttov_parallel_direct'
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
              nthreads    = rttov_direct_nthreads)    ! in    number of threads to use
    endif
    call rttov_error('rttov_direct error', lalloc = .true.)
  
  end subroutine cosp_rttov_call_direct
  

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! 7. rttov_call_direct - Call PC-RTTOV forward model (Woohoo!)
  ! ------------------------------------------------------
  ! From RTTOV example files.
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  subroutine cosp_pc_rttov_call_direct()
    
    if (rttov_direct_nthreads <= 1) then
      print*,'Calling rttov_direct (PC-RTTOV)'
      call rttov_direct(                 &
              errorstatus,               &! out   error flag
              chanprof,                  &! in    channel and profile index structure
              opts,                      &! in    options structure
              profiles,                  &! in    profile array
              coef_rttov,                &! in    coefficients structure
              transmission,              &! inout computed transmittances
              radiance,                  &! inout computed radiances
              calcemis     = calcemis,   &! in    flag for internal emissivity calcs
              emissivity   = emissivity, &! inout input/output emissivities per channel
              pccomp       = pccomp,     &! inout computed PC scores
              channels_rec = channels_rec) ! in    reconstructed channel list
    else
      print*,'Calling rttov_parallel_direct (PC-RTTOV)'
      call rttov_parallel_direct(         &
              errorstatus,                &! out   error flag
              chanprof,                   &! in    channel and profile index structure
              opts,                       &! in    options structure
              profiles,                   &! in    profile array
              coef_rttov,                 &! in    coefficients structure
              transmission,               &! inout computed transmittances
              radiance,                   &! inout computed radiances
              calcemis     = calcemis,    &! in    flag for internal emissivity calcs
              emissivity   = emissivity,  &! inout input/output emissivities per channel
              pccomp       = pccomp,      &! inout computed PC scores
              channels_rec = channels_rec,&! in    reconstructed channel list
              nthreads     = rttov_direct_nthreads)     ! in    number of threads to use
    endif
    call rttov_error('rttov_direct error (PC-RTTOV)', lalloc = .true.)
  
  end subroutine cosp_pc_rttov_call_direct


  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! 8. Save output data
  ! ------------------------------------------------------
  ! JKS - Need to allow options for Tb and radiance for clear- and cloudy-skies
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine cosp_rttov_save_output(rttovIN,                        &
                                    bt_total,bt_clear,              &
                                    rad_total,rad_clear,rad_cloudy, &
                                    refl_total,refl_clear)
    type(rttov_in),intent(in) :: &
        rttovIN
    real(wp),dimension(rttovIN%nPoints,rttovIN%nChannels),intent(inout) :: & ! Can I do this? I guess so!
        bt_total,       &
        bt_clear,       &
        rad_total,      &
        rad_clear,      &
        rad_cloudy,     &
        refl_total,     &
        refl_clear

    ! Documentation for RTTOV radiance structure in RTTOV User Guide pg 166
    
    ! Only save output if appropriate
    if (do_rttov_bt) then
        print*,'trying'
        bt_total(1:rttovIN%nPoints, 1:rttovIN%nChannels) = &
            transpose(reshape(radiance%bt(1:nchanprof), (/ rttovIN%nChannels, rttovIN%nPoints/) ))
    endif
    if (do_rttov_bt .and. (do_rttov_cld .or. do_rttov_aer)) then
        bt_clear(1:rttovIN%nPoints, 1:rttovIN%nChannels) = &
            transpose(reshape(radiance%bt_clear(1:nchanprof), (/ rttovIN%nChannels, rttovIN%nPoints/) ))
    endif
    
    if (do_rttov_rad) then
        rad_total(1:rttovIN%nPoints, 1:rttovIN%nChannels) = &
            transpose(reshape(radiance%total(1:nchanprof), (/ rttovIN%nChannels, rttovIN%nPoints/) ))
    endif
    if (do_rttov_rad .and. (do_rttov_cld .or. do_rttov_aer)) then
        rad_clear(1:rttovIN%nPoints, 1:rttovIN%nChannels) = &
            transpose(reshape(radiance%clear(1:nchanprof), (/ rttovIN%nChannels, rttovIN%nPoints/) )) 
        rad_cloudy(1:rttovIN%nPoints, 1:rttovIN%nChannels) = &
            transpose(reshape(radiance%cloudy(1:nchanprof), (/ rttovIN%nChannels, rttovIN%nPoints/) ))   
    endif

    if (do_rttov_refl) then
        refl_total(1:rttovIN%nPoints, 1:rttovIN%nChannels) = &
            transpose(reshape(radiance%refl(1:nchanprof), (/ rttovIN%nChannels, rttovIN%nPoints/) ))
    endif
    if (do_rttov_refl .and. (do_rttov_cld .or. do_rttov_aer)) then
        bt_clear(1:rttovIN%nPoints, 1:rttovIN%nChannels) = &
            transpose(reshape(radiance%refl_clear(1:nchanprof), (/ rttovIN%nChannels, rttovIN%nPoints/) ))
    endif
    print*,'done.'
          
  end subroutine cosp_rttov_save_output
  

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! 8. Save output data (PC-RTTOV)
  ! ------------------------------------------------------
  ! JKS - Need to allow options for Tb and radiance for clear- and cloudy-skies
  ! PC-RTTOV only does clear-sky IR calculations (can handle aerosols, but I'll ignore that for now.
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  subroutine cosp_pc_rttov_save_output(rttovIN,                        &
!                                       bt_total,bt_clear,              &
!                                       rad_total,rad_clear,rad_cloudy, &
!                                       refl_total,refl_clear)

  subroutine cosp_pc_rttov_save_output(rttovIN,                        &
                                       bt_total,rad_total)

    type(rttov_in),intent(in) :: &
        rttovIN
    real(wp),dimension(rttovIN%nPoints,rttovIN%nChannels),intent(inout) :: & ! Can I do this? I guess so!
        bt_total,       &
        rad_total

    ! Documentation for RTTOV radiance structure in RTTOV User Guide pg 166
        
    print*,'do_rttov_bt:   ',do_rttov_bt
    print*,'do_rttov_rad:  ',do_rttov_rad    

    ! Only save output if appropriate
    if (do_rttov_bt) then
        bt_total(1:rttovIN%nPoints, 1:rttovIN%nChannels) = &
            transpose(reshape(pccomp%bt_pccomp(1:nchanprof), (/ rttovIN%nChannels, rttovIN%nPoints/) ))
    endif
    
    if (do_rttov_rad) then
        rad_total(1:rttovIN%nPoints, 1:rttovIN%nChannels) = &
            transpose(reshape(pccomp%total_pccomp(1:nchanprof), (/ rttovIN%nChannels, rttovIN%nPoints/) ))
    endif
          
  end subroutine cosp_pc_rttov_save_output
  

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! 9. Deallocate all RTTOV arrays and structures
  ! ------------------------------------------------------
  ! From RTTOV example files.
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine cosp_rttov_deallocate_profiles(rttovIN)
          
    type(rttov_in),intent(in) :: &
        rttovIN
    
    ! JKS no longer using channel_list, rttovIN%channels instead
    !deallocate (channel_list, stat=alloc_status(1))
    !if (alloc_status(1) /= 0) then
    !  write(*,*) 'mem dellocation error'
    !endif

    ! Deallocate structures for rttov_direct
    call rttov_alloc_direct(     &
        errorstatus,             &
        0_jpim,                  &  ! 0 => deallocate
        rttovIN%nPoints,         &
        nchanprof,               &
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

  end subroutine cosp_rttov_deallocate_profiles


  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! 9. Deallocate all RTTOV arrays and structures (PC-RTTOV)
  ! ------------------------------------------------------
  ! From RTTOV example files.
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine cosp_pc_rttov_deallocate_profiles(rttovIN)
          
    type(rttov_in),intent(in) :: &
        rttovIN
    
    if (ASSOCIATED(predictindex)) deallocate (predictindex, stat=alloc_status(10))
    call rttov_error('mem dellocation error for "predictindex"', lalloc = .true.)

    if (ASSOCIATED(channels_rec)) deallocate (channels_rec, stat=alloc_status(11))
    call rttov_error('mem dellocation error for "channels_rec"', lalloc = .true.)
    
    ! Deallocate structures for rttov_direct
    call rttov_alloc_direct(     &
        errorstatus,             &
        0_jpim,                  &  ! 0 => deallocate
        rttovIN%nPoints,         &
        nchanprof,               &
        rttovIN%nLevels,         &
        chanprof,                &
        opts,                    &
        profiles,                &
        coef_rttov,              &
        transmission,            &
        radiance,                &
        calcemis=calcemis,       &
        emissivity=emissivity,   &
        npcscores=npcscores * rttovIN%nPoints,         &
        nchannels_rec=nchannels_rec * rttovIN%nPoints, &
        pccomp=pccomp)
    call rttov_error('deallocation error for rttov_direct structures (PC-RTTOV)', lalloc = .true.)

  end subroutine cosp_pc_rttov_deallocate_profiles
  

  subroutine cosp_rttov_deallocate_coefs()

    call rttov_dealloc_coefs(errorstatus, coef_rttov)
    if (errorstatus /= errorstatus_success) then
      write(*,*) 'coefs deallocation error'
    endif

  end subroutine cosp_rttov_deallocate_coefs
  
  
! Other subroutines used for bug chasing
  subroutine cosp_rttov_deallocate_coefs2()
  
!      print*,'dealloc2 test(1)'
      CALL rttov_dealloc_coef_scatt(errorstatus, coef_rttov%coef_scatt)
!      THROW(err.NE.0)

!      print*,'dealloc2 test(2)'
      IF (ASSOCIATED(coef_rttov%coef_pccomp%pcreg)) THEN
        print*,'running: "rttov_dealloc_coef_pccomp"'
        CALL rttov_dealloc_coef_pccomp(errorstatus, coef_rttov%coef_pccomp)
!        THROW(err.NE.0)
      ENDIF

!      IF (ASSOCIATED(coefs%coef_mfasis_cld%lut)) THEN
!        CALL rttov_dealloc_coef_mfasis(err, coefs%coef_mfasis_cld)
!        THROW(err.NE.0)
!      ENDIF

!      IF (ASSOCIATED(coefs%coef_mfasis_aer%lut)) THEN
!        CALL rttov_dealloc_coef_mfasis(err, coefs%coef_mfasis_aer)
!        THROW(err.NE.0)
!      ENDIF

!      IF (ASSOCIATED(coefs%coef_mfasis_nn%nn)) THEN
!        CALL rttov_dealloc_coef_mfasis_nn(err, coefs%coef_mfasis_nn)
!        THROW(err.NE.0)
!      ENDIF

!      IF (ASSOCIATED(coefs%coef_htfrtc%p)) THEN
!        CALL rttov_dealloc_coef_htfrtc(err, coefs%coef_htfrtc)
!        THROW(err.NE.0)
!      ENDIF

      print*,'dealloc2 test(3)'
      CALL rttov_dealloc_coef2(errorstatus, coef_rttov%coef)
!      CALL rttov_dealloc_coef(errorstatus, coef_rttov%coef)
!      THROW(err.NE.0)

      print*,'dealloc2 test(4)'
      coef_rttov%initialised = .FALSE.
      print*,'dealloc2 test(5)'
      
  end subroutine cosp_rttov_deallocate_coefs2
  

    SUBROUTINE rttov_dealloc_coef2(err, coef)
    !INTF_OFF
#include "throw.h"
    !INTF_ON
      USE rttov_types, ONLY : rttov_coef
      USE parkind1, ONLY : jpim
      IMPLICIT NONE

      INTEGER(KIND=jpim), INTENT(OUT)   :: err
      TYPE(rttov_coef),   INTENT(INOUT) :: coef
    !INTF_END
#include "rttov_nullify_coef.interface"
#include "rttov_errorreport.interface"
    !- End of header --------------------------------------------------------

      TRY

      print*,'dealloc_coef test(0)'

      IF (ASSOCIATED(coef%fmv_gas_id)) DEALLOCATE (coef%fmv_gas_id, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%fmv_gas_pos)) DEALLOCATE (coef%fmv_gas_pos, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%fmv_var)) DEALLOCATE (coef%fmv_var, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%fmv_coe)) DEALLOCATE (coef%fmv_coe, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%fmv_ncorr)) DEALLOCATE (coef%fmv_ncorr, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%fmv_lvl)) DEALLOCATE (coef%fmv_lvl, STAT = err)
      THROW(err.NE.0)

      print*,'dealloc_coef test(1)'

      IF (ASSOCIATED(coef%ff_ori_chn)) DEALLOCATE (coef%ff_ori_chn, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%ff_val_chn)) DEALLOCATE (coef%ff_val_chn, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%ff_cwn)) DEALLOCATE (coef%ff_cwn, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%ff_bco)) DEALLOCATE (coef%ff_bco, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%ff_bcs)) DEALLOCATE (coef%ff_bcs, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%ff_gam)) DEALLOCATE (coef%ff_gam, STAT = err)
      THROW(err.NE.0)

      print*,'dealloc_coef test(2)'

      IF (ASSOCIATED(coef%fastem_polar)) DEALLOCATE (coef%fastem_polar, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%pol_phi)) DEALLOCATE (coef%pol_phi, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%pol_fac_v)) DEALLOCATE (coef%pol_fac_v, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%pol_fac_h)) DEALLOCATE (coef%pol_fac_h, STAT = err)
      THROW(err.NE.0)

      print*,'dealloc_coef test(3)'

      IF (ASSOCIATED(coef%ssirem_a0)) DEALLOCATE (coef%ssirem_a0, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%ssirem_a1)) DEALLOCATE (coef%ssirem_a1, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%ssirem_a2)) DEALLOCATE (coef%ssirem_a2, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%ssirem_xzn1)) DEALLOCATE (coef%ssirem_xzn1, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%ssirem_xzn2)) DEALLOCATE (coef%ssirem_xzn2, STAT = err)
      THROW(err.NE.0)

      print*,'dealloc_coef test(4)'

      IF (ASSOCIATED(coef%iremis_coef)) DEALLOCATE (coef%iremis_coef, STAT = err)
      THROW(err.NE.0)


      IF (ASSOCIATED(coef%ref_prfl_p)) DEALLOCATE (coef%ref_prfl_p, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%ref_prfl_t)) DEALLOCATE (coef%ref_prfl_t, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%ref_prfl_mr)) DEALLOCATE (coef%ref_prfl_mr, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%bkg_prfl_mr)) DEALLOCATE (coef%bkg_prfl_mr, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%lim_prfl_p)) DEALLOCATE (coef%lim_prfl_p, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%env_prfl_tmax)) DEALLOCATE (coef%env_prfl_tmax, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%env_prfl_tmin)) DEALLOCATE (coef%env_prfl_tmin, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%env_prfl_gmin)) DEALLOCATE (coef%env_prfl_gmin, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%env_prfl_gmax)) DEALLOCATE (coef%env_prfl_gmax, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%lim_prfl_tmax)) DEALLOCATE (coef%lim_prfl_tmax, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%lim_prfl_tmin)) DEALLOCATE (coef%lim_prfl_tmin, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%lim_prfl_gmin)) DEALLOCATE (coef%lim_prfl_gmin, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%lim_prfl_gmax)) DEALLOCATE (coef%lim_prfl_gmax, STAT = err)
      THROW(err.NE.0)

      print*,'coef%fmv_gas:    ',coef%fmv_gas

      print*,'dealloc_coef test(5)'
!      print*,'coef%thermal:   ',coef%thermal
      CALL dealloc_fast_coefs(err, coef%thermal)
      THROW(err.NE.0)
      print*,'dealloc_coef test(5a)'

      CALL dealloc_fast_coefs(err, coef%thermal_corr)
      THROW(err.NE.0)
      print*,'dealloc_coef test(5b)'

      IF (coef%solarcoef) THEN
        print*,'dealloc_coef test(5c)'
        CALL dealloc_fast_coefs(err, coef%solar)
        THROW(err.NE.0)

        print*,'dealloc_coef test(5d)'
        CALL dealloc_fast_coefs(err, coef%solar_corr)
        THROW(err.NE.0)
      ENDIF
      
      print*,'dealloc_coef test(6)'

      IF (ASSOCIATED(coef%bounds)) DEALLOCATE (coef%bounds, STAT = err)
      THROW(err.NE.0)


      IF (coef%nltecoef) THEN
        IF (ASSOCIATED(coef%nlte_coef%coef)) DEALLOCATE(coef%nlte_coef%coef, STAT = err)
        THROW(err.NE.0)

        IF (ASSOCIATED(coef%nlte_coef%sol_zen_angle)) DEALLOCATE(coef%nlte_coef%sol_zen_angle, STAT = err)
        THROW(err.NE.0)

        IF (ASSOCIATED(coef%nlte_coef%sat_zen_angle)) DEALLOCATE(coef%nlte_coef%sat_zen_angle, STAT = err)
        THROW(err.NE.0)

        IF (ASSOCIATED(coef%nlte_coef%cos_sol)) DEALLOCATE(coef%nlte_coef%cos_sol, STAT = err)
        THROW(err.NE.0)

        IF (ASSOCIATED(coef%nlte_coef%sec_sat)) DEALLOCATE(coef%nlte_coef%sec_sat, STAT = err)
        THROW(err.NE.0)

        IF (ASSOCIATED(coef%nlte_coef)) DEALLOCATE(coef%nlte_coef, STAT = err)
        THROW(err.NE.0)
      ENDIF

      print*,'dealloc_coef test(7)'

      IF (coef%pmc_shift) THEN
        IF (ASSOCIATED(coef%pmc_ppmc)) DEALLOCATE(coef%pmc_ppmc, STAT = err)
        THROW(err.NE.0)

        IF (ASSOCIATED(coef%pmc_coef)) DEALLOCATE(coef%pmc_coef, STAT = err)
        THROW(err.NE.0)

        IF (ASSOCIATED(coef%pmc_pnominal)) DEALLOCATE(coef%pmc_pnominal, STAT = err)
        THROW(err.NE.0)
      ENDIF

      print*,'dealloc_coef test(8)'

      IF (ASSOCIATED(coef%planck1)) DEALLOCATE (coef%planck1, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%planck2)) DEALLOCATE (coef%planck2, STAT = err)
      THROW(err.NE.0)


      IF (ASSOCIATED(coef%frequency_ghz)) DEALLOCATE (coef%frequency_ghz, STAT = err)
      THROW(err.NE.0)


      IF (ASSOCIATED(coef%dp)) DEALLOCATE (coef%dp, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%dpp)) DEALLOCATE (coef%dpp, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%tstar)) &
        DEALLOCATE (coef%tstar, coef%tstar_r, coef%tstar_wsum_r, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%tstarmod_wsum_r)) &
        DEALLOCATE (coef%tstarmod_wsum_r, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%tstar_uwsum_r)) &
        DEALLOCATE (coef%tstar_uwsum_r, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%to3star)) &
        DEALLOCATE (coef%to3star, coef%to3star_r, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%wstar)) &
        DEALLOCATE (coef%wstar, coef%wstar_r, coef%wstar_wsum_r, coef%wtstar_wsum_r, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%ostar)) &
        DEALLOCATE (coef%ostar, coef%ostar_r, coef%ostar_wsum_r, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%co2star)) &
        DEALLOCATE (coef%co2star, coef%co2star_r, coef%co2star_wsum_r, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%n2ostar)) &
        DEALLOCATE (coef%n2ostar, coef%n2ostar_r, coef%n2ostar_wsum_r, coef%n2otstar_wsum_r, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%costar)) &
        DEALLOCATE (coef%costar, coef%costar_r, coef%costar_wsum_r, coef%cotstar_wsum_r, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%ch4star)) &
        DEALLOCATE (coef%ch4star, coef%ch4star_r, coef%ch4star_wsum_r, coef%ch4tstar_wsum_r, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%so2star)) &
        DEALLOCATE (coef%so2star, coef%so2star_r, coef%so2star_wsum_r, coef%so2tstar_wsum_r, STAT = err)
      THROW(err.NE.0)

      print*,'dealloc_coef test(9)'

      IF (ASSOCIATED(coef%tt_val_chn)) DEALLOCATE (coef%tt_val_chn, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%tt_a0)) DEALLOCATE (coef%tt_a0, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%tt_a1)) DEALLOCATE (coef%tt_a1, STAT = err)
      THROW(err.NE.0)

      print*,'dealloc_coef test(10)'

      IF (ASSOCIATED(coef%pw_val_chn)) DEALLOCATE (coef%pw_val_chn, STAT = err)
      THROW(err.NE.0)


      IF (ASSOCIATED(coef%ss_val_chn)) DEALLOCATE (coef%ss_val_chn, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%ss_solar_spectrum)) DEALLOCATE (coef%ss_solar_spectrum, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%ss_rayleigh_ext)) DEALLOCATE (coef%ss_rayleigh_ext, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%rayleigh_depol_gamma)) DEALLOCATE (coef%rayleigh_depol_gamma, STAT = err)
      THROW(err.NE.0)

      print*,'dealloc_coef test(11)'

      IF (ASSOCIATED(coef%refl_visnir_ow)) DEALLOCATE (coef%refl_visnir_ow, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%refl_visnir_fw)) DEALLOCATE (coef%refl_visnir_fw, STAT = err)
      THROW(err.NE.0)


      IF (ASSOCIATED(coef%woc_waopc_ow)) DEALLOCATE (coef%woc_waopc_ow, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%woc_waopc_fw)) DEALLOCATE (coef%woc_waopc_fw, STAT = err)
      THROW(err.NE.0)

      print*,'dealloc_coef test(12)'

      IF (ASSOCIATED(coef%ws_k_omega)) DEALLOCATE (coef%ws_k_omega, STAT = err)
      THROW(err.NE.0)

      IF (ASSOCIATED(coef%ws_npoint)) DEALLOCATE (coef%ws_npoint, STAT = err)
      THROW(err.NE.0)

      print*,'dealloc_coef test(13)'

      CALL rttov_nullify_coef(coef)

      print*,'dealloc_coef test(14)'
      
      CATCH

    CONTAINS

      SUBROUTINE dealloc_fast_coefs(err, fast_coefs)
        USE rttov_types, ONLY : rttov_fast_coef
        INTEGER(jpim),                  INTENT(OUT)   :: err
        TYPE(rttov_fast_coef), POINTER, INTENT(INOUT) :: fast_coefs(:)
        INTEGER(jpim) :: ichan, igas
        TRY
        print*,'dealloc_fast_coef test(0)'
        IF (ASSOCIATED(fast_coefs)) THEN ! If this part of the coef DDT is allocated
          print*,'dealloc_fast_coef test(1)'
          print*,'SIZE(fast_coefs):   ',SIZE(fast_coefs)
          DO ichan = 1, SIZE(fast_coefs) ! Iterate over each instrument channel (2378)
!            print*,ichan
!            print*,(ichan .eq. 52)
            IF (ASSOCIATED(fast_coefs(ichan)%gasarray)) THEN ! If this part of the coef DDT is allocated
!              if (ichan .eq. 52) print*,'igas:  ',igas
!              print*,'coef%fmv_gas:   ',coef%fmv_gas
!              print*,'SIZE(fast_coefs(ichan)%gasarray):   ',SIZE(fast_coefs(ichan)%gasarray)
              DO igas = 1, coef%fmv_gas ! Iterate over each gas type
                if (ichan .eq. 52) then
                    print*,'igas:  ',igas
                    print*,'igas, SIZE(fast_coefs(ichan)%gasarray(igas)%coef):  ',igas,'  ',SIZE(fast_coefs(ichan)%gasarray(igas)%coef)
                    print*,'ASSOCIATED(fast_coefs(ichan)%gasarray(igas)%coef):  ',ASSOCIATED(fast_coefs(ichan)%gasarray(igas)%coef)
                endif
!                if (ichan .eq. 52) print*,'igas:  ',igas
!                if (ichan .eq. 52) print*,'SIZE(fast_coefs(ichan)%gasarray(igas)%coef):  ',SIZE(fast_coefs(ichan)%gasarray(igas)%coef)
                IF (ASSOCIATED(fast_coefs(ichan)%gasarray(igas)%coef)) &
                  DEALLOCATE(fast_coefs(ichan)%gasarray(igas)%coef, STAT = err)
                  THROW(err.NE.0)
              ENDDO
              DEALLOCATE(fast_coefs(ichan)%gasarray, STAT = err)
              THROW(err.NE.0)
            ENDIF
          ENDDO
          DEALLOCATE(fast_coefs, STAT = err)
          THROW(err.NE.0)
        ENDIF
        CATCH
      END SUBROUTINE dealloc_fast_coefs

    END SUBROUTINE rttov_dealloc_coef2


  
end module mod_cosp_rttov
