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
#include "rttov_calc_solar_angles.interface"

  ! Module parameters
  integer, parameter :: maxlim =  10000
  real(wp),parameter :: eps    =  0.622

  ! Initialization parameters
       
  ! Scattering coefficients (read in once during initialization)
! JKS - KISS
!  type(rttov_scatt_coef) :: &
!       coef_scatt      

  ! module-wides variables for input
  !====================
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
  
  type rttov_IN
     integer(kind=jpim),pointer :: & ! JKS trying this
          nPoints,      & ! Number of profiles to simulate
          nLevels,      & ! Number of levels
          nSubCols        ! Number of subcolumns
     real(kind=wp),pointer :: &
          emis_grey => null()
     integer(kind=jpim),dimension(:),pointer :: &
          month
!     real(wp),dimension(:),pointer :: &
!          surfem           ! Surface emissivities for the channels
!          refl,         & ! Surface reflectances for the channels
     real(wp),dimension(:),pointer :: &
          h_surf,        & ! Surface height
          u_surf,        & ! U component of surface wind
          v_surf,        & ! V component of surface wind
          t_skin,        & ! Surface skin temperature
          p_surf,        & ! Surface pressure
          t2m => null(), & ! 2 m Temperature
          q2m => null(), & ! 2 m Specific humidity
          sfcmask,       & ! sea-land-ice mask (0=sea, 1=land, 2=seaice)
          latitude,      & ! Latitude (degrees)
          longitude,     & ! Longitude (degrees)
          time_frac,     & ! Fractional UTC time [0-1]
          sza => null()    ! Solar zenith angle (deg)
     real(wp),dimension(:,:),pointer :: &
          p,            & ! Pressure @ model levels
          ph,           & ! Pressure @ model half levels
          t,            & ! Temperature 
          q,            & ! Specific humidity
          o3,           & ! Ozone
          co2,          & ! Carbon dioxide 
          ch4,          & ! Methane 
          n2o,          & ! n2o 
          co,           & ! Carbon monoxide
          so2,          & ! Sulfur dioxide
          rttov_date,   & ! Date of the profile as year (e.g. 2013), month (1-12), and day (1-31)
          rttov_time,   & ! Time of profile as hour, minute, second.          
     ! These fields below are needed ONLY for the RTTOV all-sky brightness temperature
          tca,          & ! Cloud fraction
          cldIce,       & ! Cloud ice
          cldLiq,       & ! Cloud liquid
          DeffLiq,      & ! Cloud liquid effective diameter
          DeffIce,      & ! Cloud ice effective diameter
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
                                 inst_coefs,inst_iChannel,rttov_Nlocaltime, &
                                 rttov_localtime,rttov_localtime_width,     &
                                 inst_nchanprof,inst_nprof,swath_mask,debug)
                           
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
    integer(KIND=jpim),intent(in)           :: &
        rttov_Nlocaltime    
    real(kind=jprb), dimension(rttov_Nlocaltime), intent(in)    :: &
        rttov_localtime,       &
        rttov_localtime_width
    integer(kind=jpim),intent(inout) :: &
        inst_nchanprof, &
        inst_nprof          ! Now accounting for orbits
    logical(jplm),dimension(rttovIN % nPoints),intent(inout)    :: &
        swath_mask
    logical,intent(in),optional :: &
        debug
        
    !---- Local variables ----!
    ! Loop variables
    integer(kind=jpim) :: j, jch, nch

    real(kind=jprb),parameter                     :: &
        pi = 4.D0*DATAN(1.D0),  &  ! yum
        radius = 6371.0            ! Earth's radius in km (mean volumetric)

    real(kind=jprb), dimension(rttovIN % nPoints,rttov_Nlocaltime) :: &
        sat_lon,        & ! Central longitude of the instrument.
        swath_mask_all, & ! Mask of reals over all local times
        dlon,           & ! distance to satellite longitude in degrees
        dx                ! distance to satellite longitude in km?
        
    logical :: verbose = .false.
    
    if (present(debug)) verbose = debug
        
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 3. Allocate RTTOV input and output structures
    ! ------------------------------------------------------
    ! Largely from RTTOV documentation.
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Handle swathing here. Initial code from Genevieve with minor changes.
    if (rttov_Nlocaltime .gt. 0) then
        ! Iterate over local times
        do j=1,rttov_Nlocaltime
            ! Calculate the central longitude for each gridcell and orbit
            sat_lon(:,j) = 15.0 * (rttov_localtime(j) - (rttovIN%time_frac * 24.0)) 
            ! Calculate distance (in degrees) from each grid cell to the satellite central long
            dlon(:,j) = mod((rttovIN%longitude - sat_lon(:,j) + 180.0), 360.0) - 180.0             
            ! calculate distance to satellite in km. Remember to convert to radians for cos/sine calls
            dx(:,j)   = dlon(:,j) * (pi/180.0) * COS(rttovIN%latitude * pi / 180) * radius
        end do
        
        ! inside swath = 1, outside swath = 0 for "swath_mask_all"
        swath_mask_all = 0
        do j=1,rttov_Nlocaltime
            where (abs(dx(:,j))<(rttov_localtime_width(j)*0.5))
                 swath_mask_all(:,j) = 1
            end where
        end do

        ! Collapse along the Nlocaltimes dimension and shift to logicals
        swath_mask(:) = .false. ! Initialize to false
        do j = 1,rttovIN % nPoints
            if ( ANY( swath_mask_all(j,:) .eq. 1) ) then
                swath_mask(j) = .true.
            end if
        end do
    else
        swath_mask(:)  = .true. ! Compute on all columns in no local times are passed.
    end if
    ! Determine the total number of radiances to simulate (nchanprof).
    inst_nprof     = count(swath_mask)
    inst_nchanprof = inst_nChannels_rec * inst_nprof
    
    if (verbose) then
        print*,'inst_nprof:          ',inst_nprof
        print*,'inst_nChannels_rec:  ',inst_nChannels_rec
        print*,'inst_nchanprof:      ',inst_nchanprof
    end if 
    
    ! Allocate structures for rttov_direct
    call rttov_alloc_direct( &
        errorstatus,             &
        1_jpim,                  &  ! 1 => allocate
        inst_nprof,              &
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
    do j = 1, inst_nprof
      do jch = 1, inst_nChannels_rec ! nChannels
        nch = nch + 1_jpim
        chanprof(nch)%prof = j
        chanprof(nch)%chan = inst_iChannel(jch) ! Example code used channel_list
      end do
    end do
    if (verbose) print*,'Done with "cosp_rttov_allocate"'
        
  end subroutine cosp_rttov_allocate

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE pc_rttov_allocate - Subroutine for running PC-RTTOV
  ! ------------------------------------------------------
  ! 3. Allocate RTTOV input and output structures
  ! 4. Build the list of profile/channel indices in chanprof
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine cosp_pc_rttov_allocate(rttovIN,inst_PC_coef_filepath,                          &
                                    inst_coefs,inst_opts,                                   &
                                    inst_nchannels_rec,inst_iChannel_in,                    &
                                    rttov_Nlocaltime,rttov_localtime,rttov_localtime_width, &
                                    inst_nchanprof,inst_nprof,inst_iChannel_out,swath_mask, &
                                    debug)
                           
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
    integer(KIND=jpim),intent(in)           :: &
        rttov_Nlocaltime    
    real(kind=jprb), dimension(rttov_Nlocaltime), intent(in)    :: &
        rttov_localtime,       &
        rttov_localtime_width        
    integer(kind=jpim),intent(inout) :: &
        inst_nchanprof,        &
        inst_nprof
    integer(kind=jpim),intent(inout),allocatable  :: &
        inst_iChannel_out(:)      ! Passing out the channel indices
    logical(jplm),dimension(rttovIN % nPoints),intent(inout)    :: &
        swath_mask
    logical,intent(in),optional :: &
        debug

    ! Loop variables
    integer(kind=jpim) :: j, jch, nch
    integer(kind=jpim) :: lo, hi
    
    ! Local variables
    integer(kind=jpim) :: inst_npred_pc

    real(kind=jprb),parameter                     :: &
        pi = 4.D0*DATAN(1.D0),  &  ! yum
        radius = 6371.0            ! Earth's radius in km (mean volumetric)

    real(kind=jprb), dimension(rttovIN % nPoints,rttov_Nlocaltime) :: &
        sat_lon,        & ! Central longitude of the instrument.
        swath_mask_all, & ! Mask of reals over all local times
        dlon,           & ! distance to satellite longitude in degrees
        dx                ! distance to satellite longitude in km?
        
    logical :: verbose = .false.
    
    if (present(debug)) verbose = debug

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 3. Allocate RTTOV input and output structures
    ! ------------------------------------------------------
    ! Largely from RTTOV documentation.
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nullify(predictindex)
    call rttov_get_pc_predictindex(errorstatus, inst_opts, predictindex, file_pccoef=inst_PC_coef_filepath)
    call rttov_error('rttov_get_pc_predictindex fatal error' , lalloc = .false.)        

    ! Handle swathing here. Initial code from Genevieve with minor changes.
    if (rttov_Nlocaltime > 0) then
        ! Iterate over local times
        do j=1,rttov_Nlocaltime
            ! Calculate the central longitude for each gridcell and orbit
            sat_lon(:,j) = 15.0 * (rttov_localtime(j) - (rttovIN%time_frac * 24.0)) 
            ! Calculate distance (in degrees) from each grid cell to the satellite central long
            dlon(:,j) = mod((rttovIN%longitude - sat_lon(:,j) + 180.0), 360.0) - 180.0             
            ! calculate distance to satellite in km. Remember to convert to radians for cos/sine calls
            dx(:,j)   = dlon(:,j) * (pi/180.0) * COS(rttovIN%latitude * pi / 180) * radius
        end do
        
        ! inside swath = 1, outside swath = 0 for "swath_mask_all"
        swath_mask_all = 0
        do j=1,rttov_Nlocaltime
            where (abs(dx(:,j))<(rttov_localtime_width(j)*0.5))
                 swath_mask_all(:,j) = 1
            end where
        end do

        ! Collapse along the Nlocaltimes dimension and shift to logicals
        swath_mask(:) = .false. ! Initialize to false
        do j = 1,rttovIN % nPoints
            if ( ANY( swath_mask_all(j,:) .eq. 1) ) then
                swath_mask(j) = .true.
            end if
        end do
        ! Determine the total number of radiances to simulate (nchanprof).
    else
        swath_mask(:)  = .true. ! Compute on all columns in no local times are passed.
    end if
    ! Determine the total number of radiances to simulate (nchanprof).
    inst_nprof     = count(swath_mask)

    ! npred_pc is only used in the pc_rttov_allocate step so I can remove the global definition later
    inst_npred_pc  = SIZE(predictindex)
    inst_nchanprof = inst_npred_pc * inst_nprof  ! Size of chanprof array is total number of predictors over all profiles
    
    if (verbose) then
        print*,'inst_nprof:          ',inst_nprof
        print*,'inst_nChannels_rec:  ',inst_nChannels_rec
        print*,'inst_nchanprof:      ',inst_nchanprof
    end if

    ! Determine the number of reconstructed radiances per profile (nchannels_rec)    
    if (allocated(inst_iChannel_out))             deallocate(inst_iChannel_out) ! Reset because this variable is internal and used by multiple instruments.
    if (inst_opts % rt_ir % pc % addradrec) then
      if (inst_nchannels_rec < 0) then
        ! If the number of channels is negative, don't reconstruct radiances at all
        if (verbose) print*,'radrec 1.'
        inst_opts % rt_ir % pc % addradrec = .FALSE.
      else if (inst_nchannels_rec == 0) then
        ! If the number of channels is set to 0 then reconstruct all instrument channels
        if (verbose) print*,'radrec 2.'
        inst_nchannels_rec = inst_coefs % coef % fmv_chn
        allocate(inst_iChannel_out(inst_nchannels_rec))
        inst_iChannel_out    = (/ (j, j = 1, inst_nchannels_rec) /)
      else
        ! Otherwise read the channel list from the file
        if (verbose) print*,'radrec 3.'
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
          inst_nprof,                                    &
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
          npcscores=inst_opts%rt_ir%pc%npcscores * inst_nprof,         &
          nchannels_rec=inst_nchannels_rec * inst_nprof, &
          pccomp=pccomp,                                 &
          init=.TRUE._jplm)
    call rttov_error('allocation error for rttov_direct structures (PC-RTTOV)' , lalloc = .true.)

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 4. Build the list of profile/channel indices in chanprof
    ! ------------------------------------------------------
    ! Largely from RTTOV documentation.
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Populate chanprof using the channel list obtained above in predictindex(:)
    do j = 1, inst_nprof
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
                                           Lrttov_solar,   &
                                           Luser_tracegas, &
                                           Ldo_co2,        &
                                           Ldo_ch4,        &
                                           Ldo_co,         &
                                           Ldo_n2o,        &
                                           Ldo_so2,        &
                                           Ldo_o3,         &
                                           inst_co2_mr,    &
                                           inst_ch4_mr,    &
                                           inst_co_mr,     &
                                           inst_n2o_mr,    &
                                           inst_so2_mr,    &
                                           inst_zenang,    &
                                           inst_nprof,     &
                                           inst_swath_mask,&
                                           debug)

    type(rttov_in),intent(in) :: & ! What is the best way to do this? Should rttovIN be a module-wide DDT? Yes.
        rttovIN
    logical,intent(in)        :: &
        Lrttov_cld,       &
        Lrttov_aer,       &
        Lrttov_solar,     &
        Luser_tracegas,   & ! Use user-supplied trace gas columns from instrument namelists. 
        Ldo_co2,          &
        Ldo_ch4,          &
        Ldo_co,           &
        Ldo_n2o,          &
        Ldo_so2,          &
        Ldo_o3
    real(wp),intent(in)       :: &
        inst_co2_mr,      &
        inst_ch4_mr,      &
        inst_co_mr,       &
        inst_n2o_mr,      &
        inst_so2_mr,      &
        inst_zenang
    integer(kind=jpim),intent(in) :: &
        inst_nprof
    logical(kind=jplm),dimension(rttovIN % nPoints),intent(in) :: &
        inst_swath_mask
    logical,intent(in),optional :: &
        debug
    
    ! Loop variables
    integer(kind=jpim) :: i, j ! Use i to iterate over profile, j for swath_mask.
    logical :: verbose = .false.
          
    if (present(debug)) verbose = debug
          
    ! Store profile data from rttovIN in profile type.
    ! See RTTOV user guide pg 163 for description of "profiles" type
    
    ! "The rttov_profile structure is composed of the atmospheric part 
    ! and two other structures for 2 meter air and skin surface. 
    ! If you are not able to provide ozone, CO2, etc profiles the flags
    ! ozone_data, co2_data and so on in the options structure should be 
    ! set to false."

    profiles(:)%gas_units  =  1 ! kg/kg over moist air (default)
    
    if (verbose) then
        print*,'shape(rttovIN%co2):    ',shape(rttovIN%co2)
        print*,'shape(rttovIN%n2o):    ',shape(rttovIN%n2o)
        print*,'rttovIN%co2(1,:):    ',rttovIN%co2(1,1:10)
        print*,'rttovIN%n2o(1,:):    ',rttovIN%n2o(1,1:10)
!        print*,'rttovIN%t_skin:   ',rttovIN%t_skin
    end if
    
    ! Iterate over all columns
    j = 0 ! Initialize input
    do i = 1, rttovIN%nPoints
        
      if (inst_swath_mask(i)) then ! only added masked columns to profiles
          j = j + 1 ! Increment first
      
          ! Initialize trace gas concentrations from user input.
          if (Luser_tracegas) then
              if (Ldo_co2) profiles(j)%co2(:)        = inst_co2_mr
              if (Ldo_n2o) profiles(j)%n2o(:)        = inst_n2o_mr
              if (Ldo_co)  profiles(j)%co(:)         = inst_co_mr
              if (Ldo_ch4) profiles(j)%ch4(:)        = inst_ch4_mr
              if (Ldo_so2) profiles(j)%so2(:)        = inst_so2_mr
              if (Ldo_o3)  profiles(j)%o3(:)         = rttovIN%o3(i, :) ! no O3 user input set up
          else
              ! For when trace gas columns are supplied by the model. Units must match (kg/kg over moist air)
              if (Ldo_co2) profiles(j)%co2(:)        = rttovIN%co2(i,:)
              if (Ldo_n2o) profiles(j)%n2o(:)        = rttovIN%n2o(i,:)
              if (Ldo_co)  profiles(j)%co(:)         = rttovIN%co(i,:)
              if (Ldo_ch4) profiles(j)%ch4(:)        = rttovIN%ch4(i,:)
              if (Ldo_so2) profiles(j)%so2(:)        = rttovIN%so2(i,:)
              if (Ldo_o3)  profiles(j)%o3(:)         = rttovIN%o3(i, :)
          end if
          
          ! Initialize column pressure, temperature, and humidity
          profiles(j)%p(:) =  rttovIN%p(i, :) * 1e-2 ! convert Pa to hPa
          profiles(j)%t(:) =  rttovIN%t(i, :)
          profiles(j)%q(:) =  rttovIN%q(i, :)

          ! q coefficient limit is 0.1e-10
          where(profiles(j)%q(:) < 0.1e-10)
             profiles(j)%q(:) = 0.11e-10
          end where

          ! Gas profiles
          profiles(j)%o3         =  rttovIN%o3(i, :)

          ! 2m parameters
          profiles(j)%s2m%p      =  rttovIN%p_surf(i) * 1e-2 ! convert Pa to hPa
          profiles(j)%s2m%t      =  rttovIN%t2m(i)
          profiles(j)%s2m%q      =  rttovIN%q2m(i) ! Should be the same as gas units (kg/kg)
          profiles(j)%s2m%u      =  rttovIN%u_surf(i)
          profiles(j)%s2m%v      =  rttovIN%v_surf(i)
          profiles(j)%s2m%wfetc  =  10000. ! only used by sea surface solar BRDF model.

          ! skin variables for emissivity calculations
          profiles(j)%skin%t          =  rttovIN%t_skin(i)
          
          ! fastem coefficients - for mw calculations
          profiles(j)%skin%fastem(1)  =  3.0
          profiles(j)%skin%fastem(2)  =  5.0
          profiles(j)%skin%fastem(3)  =  15.0
          profiles(j)%skin%fastem(4)  =  0.1
          profiles(j)%skin%fastem(5)  =  0.3

          ! Viewing angles
          profiles(j)%zenangle      = inst_zenang ! pass in from cosp
          profiles(j)%azangle       = 0. ! hard-coded in rttov9 int JKS-?

          profiles(j)%latitude      = rttovIN%latitude(i)
          profiles(j)%longitude     = rttovIN%longitude(i)
          profiles(j)%elevation     = rttovIN%h_surf(i) * 1e-3 ! Convert m to km

          ! Solar angles.
          if (associated(rttovIN%sza)) then
             profiles(j)%sunzenangle   = rttovIN%sza(i) ! SZA in degrees
          else
             print*,'No solar zenith angle passed. Setting to zero.'
             profiles(j)%sunzenangle   = 0.
          end if
          profiles(j)%sunazangle    = 0. ! hard-coded in like rttov9

          ! surface type. sfcmask is 0 for ocean, 1 for land, and 2 for sea ice
          if (rttovIN%sfcmask(i) .lt. 0.5) then
             profiles(j)%skin%surftype  = surftype_land
          else if (rttovIN%sfcmask(i) .lt. 1.5) then
             profiles(j)%skin%surftype  = surftype_sea
          else
             profiles(j)%skin%surftype  = surftype_seaice
          end if
          
          ! land-sea mask (lsmask) indicates proportion of land in grid (not in CESM implementation! just a binary mask there)
!          if (rttovIN%lsmask(i) < 0.5) then
!            profiles(j)%skin%surftype  = surftype_sea
!          else
!            profiles(j)%skin%surftype  = surftype_land
!          endif
          ! sea-ice fraction
!          if (rttovIN%icefrac(i) >= 0.5) then
!            profiles(j)%skin%surftype  = surftype_seaice
!          endif

          ! dar: hard-coded to 1 (=ocean water) in rttov 9 int
          profiles(j)%skin%watertype = 1
          !profiles(j) %idg         = 0. ! Depreciated?
          !profiles(j) %ish         = 0. ! Depreciated?
      end if 
    end do     
        
!    if (verbose) then
!        print*,'profiles(1)%p(:):     ',profiles(1)%p(:)
!        print*,'profiles(1)%t(:):     ',profiles(1)%t(:)
!        print*,'profiles(1)%q(:):     ',profiles(1)%q(:)
!        print*,'profiles(1)%co2(:):   ',profiles(1)%co2(:)
!        print*,'profiles(1)%skin%t:   ',profiles(1)%skin%t
!        print*,'profiles(1)%s2m%t:    ',profiles(1)%s2m%t
!    end if
        
    ! JKS - nothing to check here, this will never trigger.
    call rttov_error('error in profile initialization' , lalloc = .false.)
        
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Only add the cloud fields if simulating cloud.
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if (Lrttov_cld) then

      ! Set cloud mass mixing ratio units
      profiles(:)%mmr_cldaer =  .true. ! kg/kg for cloud and aerosol (default)

      profiles(:)%clw_scheme   = 2 ! Deff scheme avoids cloud types but requires an effective diameter value
    !    profiles(:)%clwde_scheme = 1. ! Scheme for cloud liquid water cotent to effective diameter. User guide says do not change.
      profiles(:)%ice_scheme   = 1 !1:Baum 2:Baran(2014) 3:Baran(2018)
      profiles(:)%icede_param  = 2 ! 2:Wyser(recommended). Only used if ice effective diameter not input
        
      j = 0 ! Initialize input
      do i = 1,rttovIN%nPoints
      
        if (inst_swath_mask(i)) then ! only added masked columns to profiles
            j = j + 1 ! Increment profile counter      
            
            ! Cloud scheme stuff
            profiles(j)%cfrac(:)   = rttovIN%tca(i,:)    ! Cloud fraction for each layer       
            profiles(j)%cloud(1,:) = rttovIN%cldLiq(i,:) ! Cloud water mixing ratio (all in the first type for Deff)
            profiles(j)%cloud(6,:) = rttovIN%cldIce(i,:) ! Cloud ice mixing ratio (1 type). See pg 74.

            profiles(j)%clwde = rttovIN%DeffLiq(i,:) ! Cloud water effective diameter
            profiles(j)%icede = rttovIN%DeffIce(i,:) ! Cloud ice effective diameter

        ! Example UKMO input has effective radii for multiple cloud types, making identification of a single
        ! liquid droplet or ice crystal effective diameter difficult.
        ! I opt to let RTTOV decide on the effective radius values, but more complex implementation
        ! could do a more thorough conversion between UKMO output and RTTOV input
    !    profiles(j)%clwde = ! Cloud water effective diameter
    !    profiles(j)%icede = ! Cloud ice effective diameter

        ! Old code for simple cloud schemes only
    !    profiles(j)%cfraction  =  0.
    !    profiles(j)%ctp        =  500.

        ! Other options not implemented
!        profiles(j)%clw        = ! Cloud liquid water (kg/kg) – MW only,
        end if
      end do
    end if
    
    ! JKS - nothing to check here, this will never trigger.
    call rttov_error('error in cloud profile initialization' , lalloc = .false.)
    
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Only add the aerosol fields if simulating aerosol.
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if (Lrttov_aer) then
    
      ! Set aerosol mass mixing ratio units
      profiles%mmr_cldaer =  .true. ! kg/kg for cloud and aerosol (default)       
      ! Read in aerosol profiles
!      j = 0 ! Initialize input
!      do i = 1,rttovIN%nPoints
!        if (inst_swath_mask(i)) then ! only added masked columns to profiles
!            j = j + 1 ! Increment profile counter   
!            profiles(j)%aerosols(naertyp,nlayers) = rttovIN%aerosols ! Aerosols in different modes (see User Guide pg 80)
!        end if
!      end do
    end if
        
    if (Lrttov_solar) then

      print*,'rttovIN%sza(:):   ',rttovIN%sza(:)

      ! Populate longitude, latitude, time, and date profile fields
      ! Read in aerosol profiles
      j = 0 ! Initialize input
      do i = 1,rttovIN%nPoints
        if (inst_swath_mask(i)) then ! only added masked columns to profiles
            j = j + 1 ! Increment profile counter   
            profiles(j)%date(:) = rttovIN%rttov_date(i,:)
            profiles(j)%time(:) = rttovIN%rttov_time(i,:)
        end if
      end do

      ! Call functions to calculate the appropriate solar zenith and azimuthal angles.
      call RTTOV_CALC_SOLAR_ANGLES(errorstatus, profiles)
      call rttov_error('Error when calling RTTOV_CALC_SOLAR_ANGLES', lalloc = .false.)
    
      !call RTTOV_CALC_GEO_SAT_ANGLES
      
      print*,'profiles(:))%sunzenangle:    ',profiles(:)%sunzenangle
      print*,'profiles(:))%sunazangle:     ',profiles(:)%sunazangle      
      
    end if
    
    ! JKS - nothing to check here, this will never trigger.
    call rttov_error('error in aerosol profile initialization' , lalloc = .true.)
    
    ! JKS To-do: set up scattering profiles (MW only) (rttov_profile_cloud)

  end subroutine cosp_rttov_construct_profiles

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! 6. rttov_setup_emissivity_reflectance - Specify surface emissivity and reflectance
  ! ------------------------------------------------------
  ! From RTTOV example files. Will need to be expanded on to pass in values.
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  subroutine cosp_rttov_setup_emissivity_reflectance(emis_in,refl_in)
        
!    real(kind=jprb),intent(in),dimension(Npoints,Nchan),optional :: &
!        emis_in,   & ! User input emissivities
!        refl_in      ! User input reflectivities
                
!    real(kind=jprb),intent(in),optional    :: emis_in(SIZE(emissivity % emis_in))
!    real(kind=jprb),intent(in),optional    :: refl_in(SIZE(reflectance % refl_in))    

    ! Set emissivities/reflectivites here.
!    if (present(emis_in)) emissivity(:) % emis_in  = emis_in
!    if (present(refl_in)) reflectance(:) % refl_in = refl_in

  subroutine cosp_rttov_setup_emissivity_reflectance(emis_grey)
 
    real(kind=wp),intent(in),pointer,optional :: emis_grey
 
    ! In this example we have no values for input emissivities or reflectances
    ! so we initialise all inputs to zero
    call rttov_init_emis_refl(emissivity, reflectance)
    call rttov_error('error for emissivity/reflectance initialization' , lalloc = .true.)

    if (present(emis_grey)) emissivity(:) % emis_in  = emis_grey ! Use greybody emissivity.

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
                                    inst_coefs,      &
                                    debug)
  
    integer(KIND=jpim),intent(in)   :: &
        inst_nthreads
    type(rttov_options),intent(in)  :: &
        inst_opts
    type(rttov_coefs),intent(in)    :: &
        inst_coefs
    logical,intent(in),optional     :: &
        debug
    logical :: verbose = .false.
    
    if (present(debug)) verbose = debug
    
    if (verbose) then
        print*,'shape(chanprof%prof):      ',shape(chanprof%prof)
        print*,'shape(chanprof%chan):      ',shape(chanprof%chan)
        print*,'shape(profiles):           ',shape(profiles)
!    print*,'shape(profiles(:)%q):         ',shape(profiles(:)%q)
!    print*,'shape(bt_total):   ',shape(bt_total)
    end if
    
    
    if (inst_nthreads <= 1) then
      if (verbose) print*,'Calling rttov_direct'
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
      if (verbose) print*,'Calling rttov_parallel_direct'
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
                                       inst_channels_rec,            &
                                       debug)
  
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
    logical,intent(in),optional :: &
        debug
    logical :: verbose = .false.
    
    if (present(debug)) verbose = debug
    
    if (inst_nthreads <= 1) then
      if (verbose) print*,'Calling rttov_direct (PC-RTTOV)'
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
      if (verbose) print*,'Calling rttov_parallel_direct (PC-RTTOV)'
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
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  
  subroutine cosp_rttov_save_output(nPoints,inst_nchan_out,inst_swath_mask, &
                                    Lrttov_bt,Lrttov_rad,Lrttov_refl,       &
                                    Lrttov_cld,Lrttov_aer,                  &
                                    bt_total,bt_clear,                      &
                                    rad_total,rad_clear,rad_cloudy,         &
                                    refl_total,refl_clear)
    integer,intent(in)        :: &
        nPoints,    &
        inst_nchan_out
    logical,dimension(nPoints),intent(in) :: &
        inst_swath_mask
    logical,intent(in)        :: &
        Lrttov_bt,       &
        Lrttov_rad,      &
        Lrttov_refl,     &
        Lrttov_cld,      &
        Lrttov_aer
    real(wp),dimension(nPoints,inst_nchan_out),intent(inout) :: &
        bt_total,       &
        bt_clear,       &
        rad_total,      &
        rad_clear,      &
        rad_cloudy,     &
        refl_total,     &
        refl_clear
        
    ! Local iterators. i is the gridcell index. j is the swath cells index.
    integer :: i, j

    ! Documentation for RTTOV radiance structure in RTTOV User Guide pg 166
    
    ! Only save output if appropriate
    if (count(inst_swath_mask) .eq. nPoints) then ! No swathing, save all output
        if (Lrttov_bt) then
            bt_total(1:nPoints, 1:inst_nchan_out) = &
                transpose(reshape(radiance%bt(1:inst_nchan_out * nPoints), (/ inst_nchan_out, nPoints/) ))
        end if
        if (Lrttov_bt .and. (Lrttov_cld .or. Lrttov_aer)) then
            bt_clear(1:nPoints, 1:inst_nchan_out) = &
                transpose(reshape(radiance%bt_clear(1:inst_nchan_out * nPoints), (/ inst_nchan_out, nPoints/) ))
        end if

        if (Lrttov_rad) then
            rad_total(1:nPoints, 1:inst_nchan_out) = &
                transpose(reshape(radiance%total(1:inst_nchan_out * nPoints), (/ inst_nchan_out, nPoints/) ))
        end if
        if (Lrttov_rad .and. (Lrttov_cld .or. Lrttov_aer)) then
            rad_clear(1:nPoints, 1:inst_nchan_out) = &
                transpose(reshape(radiance%clear(1:inst_nchan_out * nPoints), (/ inst_nchan_out, nPoints/) )) 
            rad_cloudy(1:nPoints, 1:inst_nchan_out) = &
                transpose(reshape(radiance%cloudy(1:inst_nchan_out * nPoints), (/ inst_nchan_out, nPoints/) ))   
        end if

        if (Lrttov_refl) then
            refl_total(1:nPoints, 1:inst_nchan_out) = &
                transpose(reshape(radiance%refl(1:inst_nchan_out * nPoints), (/ inst_nchan_out, nPoints/) ))
        end if
        if (Lrttov_refl .and. (Lrttov_cld .or. Lrttov_aer)) then
            refl_clear(1:nPoints, 1:inst_nchan_out) = &
                transpose(reshape(radiance%refl_clear(1:inst_nchan_out * nPoints), (/ inst_nchan_out, nPoints/) ))
        end if
    else ! If swathing is occurring, assign the outputs appropriately
        j = 0
        do i=1,nPoints
          if (inst_swath_mask(i)) then ! only added masked columns to profiles
            if (Lrttov_bt) then
              bt_total(i, 1:inst_nchan_out) = radiance%bt(1 + (j * inst_nchan_out):(j+1) * inst_nchan_out)
            end if
            if (Lrttov_bt .and. (Lrttov_cld .or. Lrttov_aer)) then
              bt_clear(i, 1:inst_nchan_out) = radiance%bt_clear(1 + (j * inst_nchan_out):(j+1) * inst_nchan_out)
            end if            
            if (Lrttov_rad) then
              rad_total(i, 1:inst_nchan_out) = radiance%total(1 + (j * inst_nchan_out):(j+1) * inst_nchan_out)
            end if
            if (Lrttov_rad .and. (Lrttov_cld .or. Lrttov_aer)) then
              rad_clear(i, 1:inst_nchan_out) = radiance%clear(1 + (j * inst_nchan_out):(j+1) * inst_nchan_out)
              rad_cloudy(i, 1:inst_nchan_out) = radiance%cloudy(1 + (j * inst_nchan_out):(j+1) * inst_nchan_out)
            end if
            if (Lrttov_refl) then
              refl_total(i, 1:inst_nchan_out) = radiance%refl(1 + (j * inst_nchan_out):(j+1) * inst_nchan_out)
            end if
            if (Lrttov_refl .and. (Lrttov_cld .or. Lrttov_aer)) then
              refl_clear(i, 1:inst_nchan_out) = radiance%refl_clear(1 + (j * inst_nchan_out):(j+1) * inst_nchan_out)
            end if            
            j = j + 1 ! Increment profile counter afterwards
          end if
        end do
    end if
          
  end subroutine cosp_rttov_save_output
  

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! 8. Save output data (PC-RTTOV)
  ! ------------------------------------------------------
  ! PC-RTTOV only does clear-sky IR calculations (can handle aerosols, but I'll ignore that for now.
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine cosp_pc_rttov_save_output(nPoints,              &
                                       inst_nchannels_rec,   &
                                       inst_swath_mask,      &
                                       Lrttov_bt,            &
                                       Lrttov_rad,           &
                                       bt_clear,             &
                                       rad_clear,            &
                                       debug)

    integer,intent(in)        :: &
        nPoints,            &
        inst_nchannels_rec
    logical,dimension(nPoints),intent(in) :: &
        inst_swath_mask        
    logical,intent(in)        :: &
        Lrttov_bt,      &
        Lrttov_rad
    real(wp),dimension(nPoints,inst_nchannels_rec),intent(inout) :: & ! Can I do this? I guess so!
        bt_clear,       &
        rad_clear
    logical,intent(in),optional :: &
        debug
    
    ! Local iterators. i is the gridcell index. j is the swath cells index.
    integer :: i, j
    logical :: verbose = .false.
    
    if (present(debug)) verbose = debug
    
    if (verbose) then
        print*,'shape(bt_total):   ',shape(bt_clear)
        print*,'shape(rad_total):  ',shape(rad_clear)
        print*,'nPoints:   ',nPoints
        print*,'inst_nchannels_rec: ',inst_nchannels_rec
        print*,'size(pccomp%bt_pccomp):   ',size(pccomp%bt_pccomp)
        print*,'size(pccomp%total_pccomp):   ',size(pccomp%total_pccomp)
        print*,'inst_nchannels_rec * nPoints:   ',inst_nchannels_rec * nPoints
    end if

    ! Documentation for RTTOV radiance structure in RTTOV User Guide pg 166
        
    ! Only save output if appropriate
    if (count(inst_swath_mask) .eq. nPoints) then ! No swathing, save all output  
      if (Lrttov_bt) then
          bt_clear(1:nPoints, 1:inst_nchannels_rec) = &
              transpose(reshape(pccomp%bt_pccomp(1:(inst_nchannels_rec * nPoints)), (/ inst_nchannels_rec, nPoints/) ))
      end if
      if (Lrttov_rad) then
          rad_clear(1:nPoints, 1:inst_nchannels_rec) = &
              transpose(reshape(pccomp%total_pccomp(1:(inst_nchannels_rec * nPoints)), (/ inst_nchannels_rec, nPoints/) ))
      end if
    else ! If swathing is occurring, assign the outputs appropriately
      j = 0
      do i=1,nPoints
        if (inst_swath_mask(i)) then ! only added masked columns to profiles
          if (Lrttov_bt) then
            bt_clear(i, 1:inst_nchannels_rec) = pccomp%bt_pccomp(1 + (j * inst_nchannels_rec):(j+1) * inst_nchannels_rec)
          end if          
          if (Lrttov_rad) then
            rad_clear(i, 1:inst_nchannels_rec) = pccomp%total_pccomp(1 + (j * inst_nchannels_rec):(j+1) * inst_nchannels_rec)
          end if       
          j = j + 1 ! Increment profile counter afterwards
        end if
      end do
    end if
          
  end subroutine cosp_pc_rttov_save_output

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! 9. Deallocate all RTTOV arrays and structures
  ! ------------------------------------------------------
  ! From RTTOV example files.
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine cosp_rttov_deallocate_profiles(inst_nprof,           &
                                            inst_nchanprof,       &
                                            nLevels,              &
                                            inst_opts,            &
                                            inst_coefs)

    integer(kind=jpim),intent(in) :: &
        inst_nprof,        &
        inst_nchanprof,    &
        nLevels
    type(rttov_options),intent(in) :: &
        inst_opts
    type(rttov_coefs),intent(in)   :: &
        inst_coefs

        
    ! Deallocate structures for rttov_direct
    call rttov_alloc_direct(     &
        errorstatus,             &
        0_jpim,                  &  ! 0 => deallocate
        inst_nprof,              &
        inst_nchanprof,          & 
        nLevels,                 &
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

  subroutine cosp_pc_rttov_deallocate_profiles(inst_nprof,            &
                                               inst_nchanprof,        &
                                               nlevels,               &
                                               inst_nChannels_rec,    &
                                               inst_opts,             &
                                               inst_coefs)
 
    integer(kind=jpim),intent(in)  :: &
        inst_nprof,      &
        inst_nchanprof,  &
        nlevels,         &
        inst_nChannels_rec
    type(rttov_options),intent(in) :: &
        inst_opts
    type(rttov_coefs),intent(in)   :: &
        inst_coefs
        
    if (ASSOCIATED(predictindex)) deallocate (predictindex, stat=alloc_status(10))
    call rttov_error('mem dellocation error for "predictindex"', lalloc = .true.)

    if (ASSOCIATED(channels_rec)) deallocate (channels_rec, stat=alloc_status(11))
    call rttov_error('mem dellocation error for "channels_rec"', lalloc = .true.)
    
    ! Deallocate structures for rttov_direct
    call rttov_alloc_direct(     &
        errorstatus,             &
        0_jpim,                  &  ! 0 => deallocate
        inst_nprof,              &
        inst_nchanprof,          &
        nLevels,                 &
        chanprof,                &
        inst_opts,               &
        profiles,                &
        inst_coefs,              &
        transmission,            &
        radiance,                &
        calcemis=calcemis,       &
        emissivity=emissivity,   &
        npcscores=inst_opts%rt_ir%pc%npcscores * inst_nprof,         &
        nchannels_rec=inst_nChannels_rec * inst_nprof, &
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