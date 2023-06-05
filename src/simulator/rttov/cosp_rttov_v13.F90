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
#include "rttov_get_pc_predictindex.interface"

! checking inputs
!#include "rttov_dealloc_coef_scatt.interface"
!#include "rttov_dealloc_coef.interface"
!#include "rttov_dealloc_coef_pccomp.interface"

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

  ! Number of channels to reconstruct (-1: None, 1:All, otherwise the actual number)
  INTEGER(kind=jpim) :: nchannels_rec 
       
  integer(kind=jpim) ::  & ! Parallelization is default off
       rttov_direct_nthreads = 1_jpim
  integer(kind=jpim),allocatable :: &
       iChannel(:),      &  ! RTTOV channel indices
       iChannel_out(:)      ! Passing out the channel indices
  real(kind=jplm),allocatable    :: &
       emisChannel(:),   &  ! RTTOV channel emissivity
       reflChannel(:)       ! RTTOV channel reflectivity
       
  ! Scattering coefficients (read in once during initialization)
! JKS - KISS
!  type(rttov_scatt_coef) :: &
!       coef_scatt      

  ! module-wides variables for input
  !====================
  integer(kind=jpim) :: dosolar

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
  INTEGER(KIND=jpim) :: nchannels_comp, npcscores, npred_pc ! npred to go here

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! TYPE rttov_in
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  ! JKS - add additional COSP inputs here.
  type rttov_IN
     integer(kind=jpim),pointer :: & ! JKS trying this
          nPoints,      & ! Number of profiles to simulate
          nLevels,      & ! Number of levels
          nSubCols,     & ! Number of subcolumns
          month           ! Month (needed for surface emissivity calculation)
     real(wp),pointer :: & ! Could change the dimensionality of these in the future
          co2,          & ! Carbon dioxide 
          ch4,          & ! Methane 
          n2o,          & ! n2o 
          co              ! Carbon monoxide
!     real(wp),dimension(:),pointer :: &
!          surfem           ! Surface emissivities for the channels
!          refl,         & ! Surface reflectances for the channels
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

  subroutine cosp_rttov_allocate(rttovIN,inst_nChannels_rec,inst_opts,      &
                                 inst_coefs,inst_iChannel,inst_nchanprof)
                           
    type(rttov_in),intent(in)      :: &
        rttovIN
    integer(kind=jpim),intent(in)  :: &
        inst_nChannels_rec
    type(rttov_options),intent(in) :: &
        inst_opts
    type(rttov_coefs),intent(in)   :: &
        inst_coefs
    integer(kind=jpim),dimension(inst_nChannels_rec),intent(in) :: &
        inst_iChannel
    integer(kind=jpim),intent(inout) :: &
        inst_nchanprof    
    ! Loop variables
    integer(kind=jpim) :: j, jch, nch
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 3. Allocate RTTOV input and output structures
    ! ------------------------------------------------------
    ! Largely from RTTOV documentation.
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ! Determine the total number of radiances to simulate (nchanprof).    
    inst_nchanprof = inst_nChannels_rec * rttovIN%nPoints ! RTTOV (non-PC) needs nchan_out? JKS potential bug
    
    ! Allocate structures for rttov_direct
    call rttov_alloc_direct( &
        errorstatus,             &
        1_jpim,                  &  ! 1 => allocate
        rttovIN%nPoints,         &
        inst_nchanprof,          &
        rttovIN%nLevels,         &
        chanprof,                &
        inst_opts,               &
        profiles,                &
        inst_coefs,              &
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
      do jch = 1, inst_nChannels_rec ! nChannels
!      do jch = 1, rttovIN%nChannels ! nChannels
        nch = nch + 1_jpim
        chanprof(nch)%prof = j
        chanprof(nch)%chan = inst_iChannel(jch) ! Example code used channel_list
      end do
    end do
        
  end subroutine cosp_rttov_allocate

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE pc_rttov_allocate - Subroutine for running PC-RTTOV
  ! ------------------------------------------------------
  ! 3. Allocate RTTOV input and output structures
  ! 4. Build the list of profile/channel indices in chanprof
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine cosp_pc_rttov_allocate(rttovIN,inst_PC_coef_filepath, & !inst_npcscores,    &
                                    inst_coefs,inst_opts,                            &
                                    inst_nchannels_rec,inst_iChannel_in,             &
                                    inst_nchanprof, inst_iChannel_out) !,inst_predictindex)
                           
    type(rttov_in),intent(in) :: &
        rttovIN
    character(256),intent(in) :: &
        inst_PC_coef_filepath
    type(rttov_coefs),intent(in)   :: &
        inst_coefs        
    type(rttov_options),intent(inout) :: &
        inst_opts
    integer(kind=jpim),intent(inout) :: &
        inst_nchannels_rec
    integer(kind=jpim),intent(in),dimension(inst_nchannels_rec)     :: &
        inst_iChannel_in ! Channel indices the user initially requests.        
    integer(kind=jpim),intent(inout) :: &
        inst_nchanprof    
    integer(kind=jpim),intent(inout),allocatable  :: &
        inst_iChannel_out(:)      ! Passing out the channel indices

    ! Loop variables
    integer(kind=jpim) :: j, jch, nch
    integer(kind=jpim) :: lo, hi
    
    ! Local variables
    integer(kind=jpim) :: inst_npred_pc
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 3. Allocate RTTOV input and output structures
    ! ------------------------------------------------------
    ! Largely from RTTOV documentation.
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nullify(predictindex)
    call rttov_get_pc_predictindex(errorstatus, inst_opts, predictindex, file_pccoef=inst_PC_coef_filepath)
    call rttov_error('rttov_get_pc_predictindex fatal error' , lalloc = .false.)        

    ! npred_pc is only used in the pc_rttov_allocate step so I can remove the global definition later
    inst_npred_pc  = SIZE(predictindex)
    inst_nchanprof = inst_npred_pc * rttovIN%nPoints  ! Size of chanprof array is total number of predictors over all profiles
      
    ! Determine the number of reconstructed radiances per profile (nchannels_rec)
    if (inst_opts % rt_ir % pc % addradrec) then
      if (inst_nchannels_rec < 0) then
        ! If the number of channels is negative, don't reconstruct radiances at all
        print*,'radrec 1.'
        inst_opts % rt_ir % pc % addradrec = .FALSE.
      else if (inst_nchannels_rec == 0) then
        ! If the number of channels is set to 0 then reconstruct all instrument channels
        print*,'radrec 2.'
        inst_nchannels_rec = inst_coefs % coef % fmv_chn
        allocate(inst_iChannel_out(inst_nchannels_rec))
        inst_iChannel_out    = (/ (j, j = 1, inst_nchannels_rec) /)
      else
        ! Otherwise read the channel list from the file
        print*,'radrec 3.'
        print*,'inst_nchannels_rec:   ',inst_nchannels_rec
        allocate(inst_iChannel_out(inst_nchannels_rec))
        inst_iChannel_out    = inst_iChannel_in
      endif
    endif

    ! Ensure we don't have unassociated pointers below when addradrec is FALSE
    if (inst_nchannels_rec <= 0) allocate(inst_iChannel_out(0))

    ! Allocate structures for rttov_direct
    CALL rttov_alloc_direct(                             &
          errorstatus,                                   &
          1_jpim,                                        &  ! 1 => allocate
          rttovIN%nPoints,                               &
          inst_nchanprof,                                &
          rttovIN%nLevels,                               &
          chanprof,                                      & ! Make this instrument-specific? The rttov_config DDT would then be assigned to this value. Allocation difficulties?
          inst_opts,                                     &
          profiles,                                      &
          inst_coefs,                                    &
          transmission,                                  &
          radiance,                                      &
          calcemis=calcemis,                             &
          emissivity=emissivity,                         &
          npcscores=inst_opts%rt_ir%pc%npcscores * rttovIN%nPoints,         &
          nchannels_rec=inst_nchannels_rec * rttovIN%nPoints, &
          pccomp=pccomp,                                      &
          init=.TRUE._jplm)
    call rttov_error('allocation error for rttov_direct structures (PC-RTTOV)' , lalloc = .true.)

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 4. Build the list of profile/channel indices in chanprof
    ! ------------------------------------------------------
    ! Largely from RTTOV documentation.
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Populate chanprof using the channel list obtained above in predictindex(:)
    do j = 1, rttovIN%nPoints
      lo = (j - 1) * inst_npred_pc + 1
      hi = lo + inst_npred_pc - 1
      chanprof(lo:hi)%prof = j
      chanprof(lo:hi)%chan = predictindex(:)
    end do
        
  end subroutine cosp_pc_rttov_allocate


  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! 5. rttov_construct_profiles: 5. Read profile data
  ! ------------------------------------------------------
  ! Largely from cosp_rttov_v11.F90 file.
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine cosp_rttov_construct_profiles(rttovIN,        &
                                           Lrttov_cld,     &
                                           Lrttov_aer,     &
                                           inst_co2_mr,    &
                                           inst_ch4_mr,    &
                                           inst_co_mr,     &
                                           inst_n2o_mr,    &
                                           inst_so2_mr,    &
                                           inst_zenang)

    type(rttov_in),intent(in) :: & ! What is the best way to do this? Should rttovIN be a module-wide DDT? Yes.
        rttovIN
    logical,intent(in)        :: &
        Lrttov_cld,       &
        Lrttov_aer
    real(wp),intent(in)       :: &
        inst_co2_mr,      &
        inst_ch4_mr,      &
        inst_co_mr,       &
        inst_n2o_mr,      &
        inst_so2_mr,      &
        inst_zenang
    
    ! Loop variables
    integer(kind=jpim) :: i, j ! Use i to iterate over profile, j for levels.
        
    ! Store profile data from rttovIN in profile type.
    ! See RTTOV user guide pg 163 for description of "profiles" type
    
    ! "The rttov_profile structure is composed of the atmospheric part 
    ! and two other structures for 2 meter air and skin surface. 
    ! If you are not able to provide ozone, CO2, etc profiles the flags
    ! ozone_data, co2_data and so on in the options structure should be 
    ! set to false."

    profiles(:)%gas_units  =  1 ! kg/kg over moist air (default)
    
    do i = 1, rttovIN%nPoints
        
      ! Initialize trace gas concentrations from user input
      ! These gases are not in COSP input files but might be in the futre

      profiles(i)%co2(:)        = inst_co2_mr
      profiles(i)%n2o(:)        = inst_n2o_mr
      profiles(i)%co(:)         = inst_co_mr
      profiles(i)%ch4(:)        = inst_ch4_mr
      profiles(i)%so2           = inst_so2_mr ! syntax slightly different?
      
! For when trace gas columns are supplied   
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
      profiles(i)%zenangle      = inst_zenang ! pass in from cosp
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
    
    if (Lrttov_cld) then

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
    
    if (Lrttov_aer) then
    
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
  
  subroutine cosp_rttov_call_direct(inst_nthreads,   &
                                    inst_opts,       &
                                    inst_coefs)
  
    integer(KIND=jpim),intent(in)   :: &
        inst_nthreads
    type(rttov_options),intent(in) :: &
        inst_opts
    type(rttov_coefs),intent(in)   :: &
        inst_coefs
        
    if (inst_nthreads <= 1) then
      print*,'Calling rttov_direct'
      call rttov_direct(                &
              errorstatus,              &! out   error flag
              chanprof,                 &! in    channel and profile index structure
              inst_opts,                &! in    options structure
              profiles,                 &! in    profile array
              inst_coefs,               &! in    coefficients structure
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
              inst_opts,                &! in    options structure
              profiles,                 &! in    profile array
              inst_coefs,               &! in    coefficients structure
              transmission,             &! inout computed transmittances
              radiance,                 &! inout computed radiances
              calcemis    = calcemis,   &! in    flag for internal emissivity calcs
              emissivity  = emissivity, &! inout input/output emissivities per channel
              calcrefl    = calcrefl,   &! in    flag for internal BRDF calcs
              reflectance = reflectance,&! inout input/output BRDFs per channel
              nthreads    = inst_nthreads)    ! in    number of threads to use
    endif
    call rttov_error('rttov_direct error', lalloc = .true.)
  
  end subroutine cosp_rttov_call_direct
  

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! 7. rttov_call_direct - Call PC-RTTOV forward model (Woohoo!)
  ! ------------------------------------------------------
  ! From RTTOV example files.
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  subroutine cosp_pc_rttov_call_direct(inst_nthreads,                &
                                       inst_opts,                    &
                                       inst_coefs,                   &
                                       inst_nchannels_rec,           &
                                       inst_channels_rec)
  
    integer(KIND=jpim),intent(in)   :: &
        inst_nthreads
    type(rttov_options),intent(in) :: &
        inst_opts
    type(rttov_coefs),intent(in)   :: &
        inst_coefs
    integer(jpim),intent(in) :: &
        inst_nchannels_rec
    integer(jpim),dimension(inst_nchannels_rec),intent(in)  :: &
        inst_channels_rec
    
    if (inst_nthreads <= 1) then
      print*,'Calling rttov_direct (PC-RTTOV)'
      call rttov_direct(                 &
              errorstatus,               &! out   error flag
              chanprof,                  &! in    channel and profile index structure
              inst_opts,                 &! in    options structure
              profiles,                  &! in    profile array
              inst_coefs,                &! in    coefficients structure
              transmission,              &! inout computed transmittances
              radiance,                  &! inout computed radiances
              calcemis     = calcemis,   &! in    flag for internal emissivity calcs
              emissivity   = emissivity, &! inout input/output emissivities per channel
              pccomp       = pccomp,     &! inout computed PC scores
              channels_rec = inst_channels_rec) ! in    reconstructed channel list
    else
      print*,'Calling rttov_parallel_direct (PC-RTTOV)'
      call rttov_parallel_direct(         &
              errorstatus,                &! out   error flag
              chanprof,                   &! in    channel and profile index structure
              inst_opts,                  &! in    options structure
              profiles,                   &! in    profile array
              inst_coefs,                 &! in    coefficients structure
              transmission,               &! inout computed transmittances
              radiance,                   &! inout computed radiances
              calcemis     = calcemis,    &! in    flag for internal emissivity calcs
              emissivity   = emissivity,  &! inout input/output emissivities per channel
              pccomp       = pccomp,      &! inout computed PC scores
              channels_rec = inst_channels_rec,&! in    reconstructed channel list
              nthreads     = inst_nthreads)     ! in    number of threads to use
    endif
    call rttov_error('rttov_direct error (PC-RTTOV)', lalloc = .true.)
  
  end subroutine cosp_pc_rttov_call_direct

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! 8. Save output data
  ! ------------------------------------------------------
  ! JKS - Need to allow options for Tb and radiance for clear- and cloudy-skies
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  
  subroutine cosp_rttov_save_output(rttovIN,inst_nchan_out,         & ! Inputs
                                    Lrttov_bt,Lrttov_rad,Lrttov_refl,  &
                                    Lrttov_cld,Lrttov_aer,             &
                                    bt_total,bt_clear,                 &
                                    rad_total,rad_clear,rad_cloudy,    &
                                    refl_total,refl_clear)
    type(rttov_in),intent(in) :: &
        rttovIN
    integer,intent(in)        :: &
        inst_nchan_out
    logical,intent(in)        :: &
        Lrttov_bt,       &
        Lrttov_rad,      &
        Lrttov_refl,     &
        Lrttov_cld,      &
        Lrttov_aer
        
    real(wp),dimension(rttovIN%nPoints,inst_nchan_out),intent(inout) :: & ! Can I do this? I guess so!
        bt_total,       &
        bt_clear,       &
        rad_total,      &
        rad_clear,      &
        rad_cloudy,     &
        refl_total,     &
        refl_clear

    ! Documentation for RTTOV radiance structure in RTTOV User Guide pg 166
    
    ! Only save output if appropriate
    if (Lrttov_bt) then
        bt_total(1:rttovIN%nPoints, 1:inst_nchan_out) = &
            transpose(reshape(radiance%bt(1:inst_nchan_out * rttovIN%nPoints), (/ inst_nchan_out, rttovIN%nPoints/) ))
    endif
    if (Lrttov_bt .and. (Lrttov_cld .or. Lrttov_aer)) then
        bt_clear(1:rttovIN%nPoints, 1:inst_nchan_out) = &
            transpose(reshape(radiance%bt_clear(1:inst_nchan_out * rttovIN%nPoints), (/ inst_nchan_out, rttovIN%nPoints/) ))
    endif
    
    if (Lrttov_rad) then
        rad_total(1:rttovIN%nPoints, 1:inst_nchan_out) = &
            transpose(reshape(radiance%total(1:inst_nchan_out * rttovIN%nPoints), (/ inst_nchan_out, rttovIN%nPoints/) ))
    endif
    if (Lrttov_rad .and. (Lrttov_cld .or. Lrttov_aer)) then
        rad_clear(1:rttovIN%nPoints, 1:inst_nchan_out) = &
            transpose(reshape(radiance%clear(1:inst_nchan_out * rttovIN%nPoints), (/ inst_nchan_out, rttovIN%nPoints/) )) 
        rad_cloudy(1:rttovIN%nPoints, 1:inst_nchan_out) = &
            transpose(reshape(radiance%cloudy(1:inst_nchan_out * rttovIN%nPoints), (/ inst_nchan_out, rttovIN%nPoints/) ))   
    endif

    if (Lrttov_refl) then
        refl_total(1:rttovIN%nPoints, 1:inst_nchan_out) = &
            transpose(reshape(radiance%refl(1:inst_nchan_out * rttovIN%nPoints), (/ inst_nchan_out, rttovIN%nPoints/) ))
    endif
    if (Lrttov_refl .and. (Lrttov_cld .or. Lrttov_aer)) then
        bt_clear(1:rttovIN%nPoints, 1:inst_nchan_out) = &
            transpose(reshape(radiance%refl_clear(1:inst_nchan_out * rttovIN%nPoints), (/ inst_nchan_out, rttovIN%nPoints/) ))
    endif
          
  end subroutine cosp_rttov_save_output
  

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! 8. Save output data (PC-RTTOV)
  ! ------------------------------------------------------
  ! JKS - Need to allow options for Tb and radiance for clear- and cloudy-skies
  ! PC-RTTOV only does clear-sky IR calculations (can handle aerosols, but I'll ignore that for now.
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine cosp_pc_rttov_save_output(nPoints,              &
                                       inst_nchannels_rec,   &
                                       Lrttov_bt,            &
                                       Lrttov_rad,           &
                                       bt_clear,             &
                                       rad_clear)

    integer,intent(in)        :: &
        nPoints,            &
        inst_nchannels_rec
    logical,intent(in)        :: &
        Lrttov_bt,      &
        Lrttov_rad
    real(wp),dimension(nPoints,inst_nchannels_rec),intent(inout) :: & ! Can I do this? I guess so!
        bt_clear,       &
        rad_clear

!    print*,'shape(bt_total):   ',shape(bt_total)
!    print*,'shape(rad_total):  ',shape(rad_total)
!    print*,'rttovIN%nPoints:   ',rttovIN%nPoints
!    print*,'rttovIN%nChannels: ',rttovIN%nChannels
!    print*,'nchanprof:    ',nchanprof ! This is the number of predictors so not the reconstructed channel dimension
!    print*,'size(pccomp%bt_pccomp):   ',size(pccomp%bt_pccomp)
!    print*,'size(pccomp%total_pccomp):   ',size(pccomp%total_pccomp)
!    print*,'nchannels_rec * rttovIN%nPoints:   ',nchannels_rec * rttovIN%nPoints

    ! JKS why not just pass in rttovIN%nPoints and use nchannels_rec here? TO-DO?
    ! Documentation for RTTOV radiance structure in RTTOV User Guide pg 166
        
    ! Only save output if appropriate
    if (Lrttov_bt) then
        bt_clear(1:nPoints, 1:inst_nchannels_rec) = &
            transpose(reshape(pccomp%bt_pccomp(1:(inst_nchannels_rec * nPoints)), (/ inst_nchannels_rec, nPoints/) ))
    endif
    
    if (Lrttov_rad) then
        rad_clear(1:nPoints, 1:inst_nchannels_rec) = &
            transpose(reshape(pccomp%total_pccomp(1:(inst_nchannels_rec * nPoints)), (/ inst_nchannels_rec, nPoints/) ))
    endif
          
  end subroutine cosp_pc_rttov_save_output

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! 9. Deallocate all RTTOV arrays and structures
  ! ------------------------------------------------------
  ! From RTTOV example files.
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine cosp_rttov_deallocate_profiles(rttovIN,              &
                                            inst_opts,            &
                                            inst_coefs,           &
                                            inst_nchanprof)

    type(rttov_in),intent(in) :: &
        rttovIN
!    integer(kind=jpim),intent(in)  :: &
!        inst_nChannels_rec
    type(rttov_options),intent(in) :: &
        inst_opts
    type(rttov_coefs),intent(in)   :: &
        inst_coefs
    integer(kind=jpim),intent(in) :: &
        inst_nchanprof    
        
    ! Deallocate structures for rttov_direct
    call rttov_alloc_direct(     &
        errorstatus,             &
        0_jpim,                  &  ! 0 => deallocate
        rttovIN%nPoints,         &
        inst_nchanprof,          & 
        rttovIN%nLevels,         &
        chanprof,                & ! JKS
        inst_opts,               &
        profiles,                &
        inst_coefs,              &
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

  subroutine cosp_pc_rttov_deallocate_profiles(rttovIN,               &
                                               inst_nChannels_rec,    &
                                               inst_opts,             &
                                               inst_coefs,            &
                                               inst_nchanprof)
                                                  
    type(rttov_in),intent(in)      :: &
        rttovIN
    integer(kind=jpim),intent(in)  :: &
        inst_nChannels_rec
    type(rttov_options),intent(in) :: &
        inst_opts
    type(rttov_coefs),intent(in)   :: &
        inst_coefs
    integer(kind=jpim),intent(in)  :: &
        inst_nchanprof    
        
    if (ASSOCIATED(predictindex)) deallocate (predictindex, stat=alloc_status(10))
    call rttov_error('mem dellocation error for "predictindex"', lalloc = .true.)

    if (ASSOCIATED(channels_rec)) deallocate (channels_rec, stat=alloc_status(11))
    call rttov_error('mem dellocation error for "channels_rec"', lalloc = .true.)
    
    ! Deallocate structures for rttov_direct
    call rttov_alloc_direct(     &
        errorstatus,             &
        0_jpim,                  &  ! 0 => deallocate
        rttovIN%nPoints,         &
        inst_nchanprof,          &
        rttovIN%nLevels,         &
        chanprof,                &
        inst_opts,               &
        profiles,                &
        inst_coefs,              &
        transmission,            &
        radiance,                &
        calcemis=calcemis,       &
        emissivity=emissivity,   &
        npcscores=inst_opts%rt_ir%pc%npcscores * rttovIN%nPoints,         &
        nchannels_rec=inst_nChannels_rec * rttovIN%nPoints, &
        pccomp=pccomp)
    call rttov_error('deallocation error for rttov_direct structures (PC-RTTOV)', lalloc = .true.)

  end subroutine cosp_pc_rttov_deallocate_profiles
  
  subroutine cosp_rttov_deallocate_coefs(inst_coefs)

    type(rttov_coefs),intent(inout)   :: &
        inst_coefs

    call rttov_dealloc_coefs(errorstatus, inst_coefs)
    if (errorstatus /= errorstatus_success) then
      write(*,*) 'coefs deallocation error'
    endif

  end subroutine cosp_rttov_deallocate_coefs

!##########################
! Module End
!##########################
end module mod_cosp_rttov