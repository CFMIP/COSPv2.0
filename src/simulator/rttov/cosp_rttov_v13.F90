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
  use cosp_phys_constants, only : mdry=>amd,mH2O=>amw,mO3=>amO3,mCO2=>amCO2,              &
                                  mCH4=>amCH4,mN2O=>amN2O,mCO=>amCO,mSO2=>amSO2


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
       
  ! Scattering coefficients (read in once during initialization)
! JKS - KISS
!  type(rttov_scatt_coef) :: &
!       coef_scatt      

  ! module-wides variables for input. Not sure if unsafe for threading.
  !====================
  
  INTEGER(KIND=jpim) :: errorstatus              ! Return error status of RTTOV subroutine calls
  INTEGER(KIND=jpim) :: alloc_status(60)

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
  ! SUBROUTINE cosp_rttov_swath - JKS
  ! ------------------------------------------------------
  ! Determine which gridcells should be swathed.
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine cosp_rttov_swath(rttovIN,rttov_Nlocaltime,               &
                              rttov_localtime,rttov_localtime_width,  &
                              inst_swath_mask, debug)

    type(rttov_in),intent(in)      :: &
        rttovIN
    integer(KIND=jpim),intent(in)           :: &
        rttov_Nlocaltime
    real(kind=jprb), dimension(rttov_Nlocaltime), intent(in)    :: &
        rttov_localtime,       &
        rttov_localtime_width                                
    logical(jplm),dimension(rttovIN % nPoints),intent(inout)    :: &
        inst_swath_mask
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
        dlon,           & ! distance to satellite longitude in degrees
        dx                ! distance to satellite longitude in km?

    logical(kind=jplm), dimension(rttovIN % nPoints,rttov_Nlocaltime) :: &
        swath_mask_all    ! Mask of logicals over all local times

    integer, dimension(rttovIN % nPoints) :: &
        rttov_DOY         ! Array of day of year values
    real(kind=jprb), dimension(rttovIN % nPoints) :: &
        localtime_offsets ! Offset values to avoid striping with hourly RT calls. [hours]
    logical :: verbose = .false.
    
    if (present(debug)) verbose = debug

    ! Compute the day of the year and determine the localtime offset
    do j=1,rttovIN%nPoints
        call get_DOY(int(rttovIN%rttov_date(j,2)), int(rttovIN%rttov_date(j,3)), rttov_DOY(j))
    end do
    localtime_offsets = (mod(rttov_DOY(:), 5) - 2) / 5.0  ! Need to cast to real

    ! Handle swathing here. Initial code from Genevieve with implementation changes.
    swath_mask_all = .false.
    if (rttov_Nlocaltime .gt. 0) then
        ! Iterate over local times
        do j=1,rttov_Nlocaltime
            ! Calculate the central longitude for each gridcell and orbit
            sat_lon(:,j) = 15.0 * (rttov_localtime(j) + localtime_offsets(:) - (rttovIN%rttov_time(:,1) + rttovIN%rttov_time(:,2) / 60))
            ! Calculate distance (in degrees) from each grid cell to the satellite central long
            dlon(:,j) = mod((rttovIN%longitude - sat_lon(:,j) + 180.0), 360.0) - 180.0             
            ! calculate distance to satellite in km. Remember to convert to radians for cos/sine calls
            dx(:,j)   = dlon(:,j) * (pi/180.0) * COS(rttovIN%latitude * pi / 180) * radius
            ! Determine if a gridcell falls in the swath width
            where (abs(dx(:,j))<(rttov_localtime_width(j)*0.5))
                swath_mask_all(:,j) = .true.
            end where                    
        end do

        inst_swath_mask = ANY( swath_mask_all(:,:),2)
    else
        inst_swath_mask(:)  = .true. ! Compute on all columns in no local times are passed.
    end if
    if (verbose) print*,'inst_swath_mask:  ',inst_swath_mask

  end subroutine cosp_rttov_swath

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE rttov_allocate - JKS
  ! ------------------------------------------------------
  ! 3. Allocate RTTOV input and output structures
  ! 4. Build the list of profile/channel indices in inst_chanprof
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine cosp_rttov_allocate(rttovIN,inst_nChannels_rec,inst_opts,inst_coefs, &
                                 inst_profiles, inst_iChannel, inst_chanprof,     &
                                 inst_nchanprof,inst_nprof,inst_swath_mask,       &
                                 inst_extend_atmos,                               &
                                 inst_transmission,inst_radiance,inst_calcemis,   &
                                 inst_emissivity,inst_calcrefl,inst_reflectance,  &
                                 debug)
                           
    type(rttov_in),intent(in)      :: &
        rttovIN
    integer(kind=jpim),intent(in)  :: &
        inst_nChannels_rec
    type(rttov_options),intent(in) :: &
        inst_opts
    type(rttov_coefs),intent(in)   :: &
        inst_coefs
    type(rttov_profile),pointer,intent(out) :: &
        inst_profiles(:)  
    integer(kind=jpim),dimension(inst_nChannels_rec),intent(in) :: &
        inst_iChannel
    type(rttov_chanprof),pointer,intent(inout) :: &
        inst_chanprof(:)    
    integer(kind=jpim),intent(inout) :: &
        inst_nchanprof
    integer(kind=jpim),intent(in) :: &
        inst_nprof          ! Now accounting for orbits
    logical(jplm),dimension(rttovIN % nPoints),intent(inout)    :: &
        inst_swath_mask
    integer(kind=jpim),intent(in)  :: &
        inst_extend_atmos
    type(rttov_transmission),intent(out) :: &
        inst_transmission
    type(rttov_radiance),intent(out) :: &
        inst_radiance
    logical(kind=jplm),pointer,intent(out) :: &
        inst_calcemis(:)
    type(rttov_emissivity),pointer,intent(out) :: &
        inst_emissivity(:)
    logical(kind=jplm),pointer,intent(out) :: &
        inst_calcrefl(:)
    type(rttov_reflectance),pointer,intent(out) :: &
        inst_reflectance(:)
    logical,intent(in),optional :: &
        debug
        
    !---- Local variables ----!
    ! Loop variables
    integer(kind=jpim) :: j, jch, nch, nlevels_rttov

    logical :: verbose = .false.
    
    if (present(debug)) verbose = debug

    if (inst_extend_atmos .eq. 0) nlevels_rttov = rttovIN%nLevels+1 ! Just use pressure levels that are supplied.
    if (inst_extend_atmos .eq. 1) nlevels_rttov = rttovIN%nLevels+2 ! Simplying extend the atmosphere with a single top layer. CAM6-like.
    ! To-do: implement a SARTA-like interpolation to a standard atmosphere.
        
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 3. Allocate RTTOV input and output structures
    ! ------------------------------------------------------
    ! Largely from RTTOV documentation.
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    ! Determine the total number of radiances to simulate (nchanprof).
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
        nlevels_rttov,           & ! "levels" means interfaces, not layers
        inst_chanprof,           &
        inst_opts,               &
        inst_profiles,           &
        inst_coefs,              &
        inst_transmission,       &
        inst_radiance,           &
        calcemis=inst_calcemis,       &
        emissivity=inst_emissivity,   &
        calcrefl=inst_calcrefl,       &
        reflectance=inst_reflectance, &
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
        inst_chanprof(nch)%prof = j
        inst_chanprof(nch)%chan = inst_iChannel(jch) ! Example code used channel_list
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

  subroutine cosp_pc_rttov_allocate(rttovIN,inst_PC_coef_filepath,                            &
                                    inst_coefs,inst_opts,inst_profiles,                       &
                                    inst_nchannels_rec,inst_iChannel_in,inst_chanprof,        &
                                    inst_nchanprof,inst_nprof,inst_iChannel_out,              &
                                    inst_swath_mask,inst_extend_atmos,inst_transmission,      &
                                    inst_radiance, inst_calcemis,inst_emissivity,inst_pccomp, &
                                    inst_predictindex,debug)
                           
    type(rttov_in),intent(in) :: &
        rttovIN
    character(256),intent(in) :: &
        inst_PC_coef_filepath
    type(rttov_coefs),intent(in)   :: &
        inst_coefs        
    type(rttov_options),intent(inout) :: &
        inst_opts
    type(rttov_profile),pointer,intent(out) :: &
        inst_profiles(:)        
    integer(kind=jpim),intent(inout) :: &
        inst_nchannels_rec
    integer(kind=jpim),intent(in),dimension(inst_nchannels_rec)     :: &
        inst_iChannel_in ! Channel indices the user initially requests.
    type(rttov_chanprof),pointer,intent(inout) :: &
        inst_chanprof(:)
    integer(kind=jpim),intent(inout) :: &
        inst_nchanprof
    integer(kind=jpim),intent(in) :: &
        inst_nprof
    integer(kind=jpim),intent(inout),allocatable  :: &
        inst_iChannel_out(:)      ! Passing out the channel indices
    logical(jplm),dimension(rttovIN % nPoints),intent(inout)    :: &
        inst_swath_mask
    integer(kind=jpim),intent(in)  :: &
        inst_extend_atmos        
    type(rttov_transmission),intent(out) :: &
        inst_transmission
    type(rttov_radiance),intent(out) :: &
        inst_radiance        
    logical(kind=jplm),pointer,intent(out) :: &
        inst_calcemis(:)
    type(rttov_emissivity),pointer,intent(out) :: &
        inst_emissivity(:)        
    type(rttov_pccomp),intent(inout) :: &
        inst_pccomp     ! Output PC structure
    integer(kind=jpim),pointer,intent(inout) :: &
        inst_predictindex(:)
    logical,intent(in),optional :: &
        debug

    ! Loop variables
    integer(kind=jpim) :: j, jch, nch, nlevels_rttov
    integer(kind=jpim) :: lo, hi
    
    ! Local variables
    integer(kind=jpim) :: inst_npred_pc
        
    logical :: verbose = .false.
    
    if (present(debug)) verbose = debug

    if (inst_extend_atmos .eq. 0) nlevels_rttov = rttovIN%nLevels+1 ! Just use pressure levels that are supplied.
    if (inst_extend_atmos .eq. 1) nlevels_rttov = rttovIN%nLevels+2 ! Simplying extend the atmosphere with a single top layer. CAM6-like.
    ! To-do: implement a SARTA-like interpolation to a standard atmosphere.    

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 3. Allocate RTTOV input and output structures
    ! ------------------------------------------------------
    ! Largely from RTTOV documentation.
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nullify(inst_predictindex)
    call rttov_get_pc_predictindex(errorstatus, inst_opts, inst_predictindex, file_pccoef=inst_PC_coef_filepath)
    call rttov_error('rttov_get_pc_predictindex fatal error' , lalloc = .false.)

    ! Determine the total number of radiances to simulate (nchanprof).
    inst_npred_pc  = SIZE(inst_predictindex)
    inst_nchanprof = inst_npred_pc * inst_nprof  ! Size of chanprof array is total number of predictors over all profiles
    
    ! if (verbose) then
    !     print*,'inst_nprof:          ',inst_nprof
    !     print*,'inst_nChannels_rec:  ',inst_nChannels_rec
    !     print*,'inst_nchanprof:      ',inst_nchanprof
    ! end if

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
          nlevels_rttov,                                 & ! "levels" means interfaces, not layers
          inst_chanprof,                                 &
          inst_opts,                                     &
          inst_profiles,                                 &
          inst_coefs,                                    &
          inst_transmission,                             &
          inst_radiance,                                 &
          calcemis=inst_calcemis,                        &
          emissivity=inst_emissivity,                    &
          npcscores=inst_opts%rt_ir%pc%npcscores * inst_nprof,         &
          nchannels_rec=inst_nchannels_rec * inst_nprof, &
          pccomp=inst_pccomp,                            &
          init=.TRUE._jplm)
    call rttov_error('allocation error for rttov_direct structures (PC-RTTOV)' , lalloc = .true.)

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 4. Build the list of profile/channel indices in chanprof
    ! ------------------------------------------------------
    ! Largely from RTTOV documentation.
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Populate chanprof using the channel list obtained above in inst_predictindex(:)
    do j = 1, inst_nprof
      lo = (j - 1) * inst_npred_pc + 1
      hi = lo + inst_npred_pc - 1
      inst_chanprof(lo:hi)%prof = j
      inst_chanprof(lo:hi)%chan = inst_predictindex(:)
    end do
        
  end subroutine cosp_pc_rttov_allocate


  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! 5. rttov_construct_profiles: 5. Read profile data
  ! ------------------------------------------------------
  ! Largely from cosp_rttov_v11.F90 file.
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine cosp_rttov_construct_profiles(rttovIN,            &
                                           inst_profiles,      &
                                           Lrttov_cld,         &
                                           Lrttov_aer,         &
                                           Lrttov_solar,       &
                                           Luser_tracegas,     &
                                           Ldo_co2,            &
                                           Ldo_ch4,            &
                                           Ldo_co,             &
                                           Ldo_n2o,            &
                                           Ldo_so2,            &
                                           Ldo_o3,             &
                                           inst_co2_mr,        &
                                           inst_ch4_mr,        &
                                           inst_co_mr,         &
                                           inst_n2o_mr,        &
                                           inst_so2_mr,        &
                                           inst_zenang,        &
                                           inst_nprof,         &
                                           inst_swath_mask,    &
                                           inst_gas_units,     &
                                           inst_clw_scheme,    &
                                           inst_ice_scheme,    &
                                           inst_icede_param,   &
                                           inst_extend_atmos,  &
                                           debug)

    type(rttov_in),intent(in) :: &
        rttovIN
    type(rttov_profile),pointer,intent(inout) :: &    
        inst_profiles(:)
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
        inst_nprof,       &
        inst_gas_units,   &
        inst_clw_scheme,  &
        inst_ice_scheme,  &
        inst_icede_param, &
        inst_extend_atmos
    logical(kind=jplm),dimension(rttovIN % nPoints),intent(in) :: &
        inst_swath_mask
    logical,intent(in),optional :: &
        debug
    
    ! Loop variables
    integer(kind=jpim) :: i, j, k, kk          ! Use i to iterate over profile, j for swath_mask, k for vertical interpolation, kk to move from profile vertical coord to rttovIN coord
    integer(kind=jpim) :: modeltop_index
    logical :: verbose = .false.

    if (present(debug)) verbose = debug
          
    ! Store profile data from rttovIN in profile type.
    ! See RTTOV user guide pg 163 for description of "profiles" type
    
    ! "The rttov_profile structure is composed of the atmospheric part 
    ! and two other structures for 2 meter air and skin surface. 
    ! If you are not able to provide ozone, CO2, etc profiles the flags
    ! ozone_data, co2_data and so on in the options structure should be 
    ! set to false."
        
    ! Iterate over all columns
    j = 0 ! Initialize input
    do i = 1, rttovIN%nPoints
      if (i .gt. rttovIN%nPoints) exit  
      if (inst_swath_mask(i)) then ! only added masked columns to profiles
          j = j + 1 ! Increment first

          ! To imitate CAM6-RRTMG, set a model top index. If inst_extend_atmos==1 then just increment that index by one and retroactively apply it to the top layer.
          ! This will be the same index that the vertical interpolation operates over so it will fix that too!

          ! top layer thickness
          if (inst_extend_atmos == 0) then
              inst_profiles(j)%p(:) =  rttovIN%ph(i, :) * 1e-2 ! convert Pa to hPa. Pressure on levels.
              if (inst_profiles(j)%p(1) .le. 0) inst_profiles(j)%p(1) = 2.25 ! If the model top is set to zero (like the COSPv2 driver and CAM6) make it 2.25mb. CAM-like correction
              modeltop_index = 1
          else if (inst_extend_atmos==1) then 
              modeltop_index = 2
              inst_profiles(j)%p(1) = 1e-4
              inst_profiles(j)%p(2) = 2.25
              inst_profiles(j)%p(3:inst_profiles(j)%nlevels) = rttovIN%ph(i, 2:rttovIN%nlevels) * 1e-2 ! JKS confirm correct indices.
          end if

          ! Handle the top and bottom levels separately.
          ! Top interface supplied by the input. 
          call interpolate_logp(100*inst_profiles(j)%p(modeltop_index),rttovIN%p(i,1),rttovIN%p(i,2),rttovIN%q(i,1),rttovIN%q(i,2),inst_profiles(j)%q(modeltop_index))
          call interpolate_logp(100*inst_profiles(j)%p(modeltop_index),rttovIN%p(i,1),rttovIN%p(i,2),rttovIN%t(i,1),rttovIN%t(i,2),inst_profiles(j)%t(modeltop_index))
          ! Bottom
          call interpolate_logp(100*inst_profiles(j)%p(inst_profiles(j)%nlevels),rttovIN%p(i,rttovIN%nlevels-1),rttovIN%p(i,rttovIN%nlevels),rttovIN%q(i,rttovIN%nlevels-1),rttovIN%q(i,rttovIN%nlevels),inst_profiles(j)%q(inst_profiles(j)%nlevels))
          call interpolate_logp(100*inst_profiles(j)%p(inst_profiles(j)%nlevels),rttovIN%p(i,rttovIN%nlevels-1),rttovIN%p(i,rttovIN%nlevels),rttovIN%t(i,rttovIN%nlevels-1),rttovIN%t(i,rttovIN%nlevels),inst_profiles(j)%t(inst_profiles(j)%nlevels))

          do k=modeltop_index+1,inst_profiles(j)%nlevels-1 ! Iterate from the layer directly below the first supplied model layer (top of model) to the layer above the last supplied model layer (surface)
            kk = k + 1 - modeltop_index
            call interpolate_logp(100*inst_profiles(j)%p(k),rttovIN%p(i,kk-1),rttovIN%p(i,kk),rttovIN%q(i,kk-1),rttovIN%q(i,kk),inst_profiles(j)%q(k))
            call interpolate_logp(100*inst_profiles(j)%p(k),rttovIN%p(i,kk-1),rttovIN%p(i,kk),rttovIN%t(i,kk-1),rttovIN%t(i,kk),inst_profiles(j)%t(k))
          end do

          ! Trace gas concentrations on levels (not layers!)
          ! Initialize trace gas concentrations from user input.
          if (Luser_tracegas) then
            if (Ldo_co2) inst_profiles(j)%co2(:)        = inst_co2_mr
            if (Ldo_n2o) inst_profiles(j)%n2o(:)        = inst_n2o_mr
            if (Ldo_co)  inst_profiles(j)%co(:)         = inst_co_mr
            if (Ldo_ch4) inst_profiles(j)%ch4(:)        = inst_ch4_mr
            if (Ldo_so2) inst_profiles(j)%so2(:)        = inst_so2_mr
            ! if (Ldo_o3)  inst_profiles(j)%o3(:)         = rttovIN%o3(i, :)
            if (Ldo_o3)  then ! no O3 user input set up
                call interpolate_logp(100*inst_profiles(j)%p(modeltop_index),rttovIN%p(i,1),rttovIN%p(i,2),rttovIN%o3(i,1),rttovIN%o3(i,2),inst_profiles(j)%o3(modeltop_index))
                call interpolate_logp(100*inst_profiles(j)%p(inst_profiles(j)%nlevels),rttovIN%p(i,rttovIN%nlevels-1),rttovIN%p(i,rttovIN%nlevels),rttovIN%o3(i,rttovIN%nlevels-1),rttovIN%o3(i,rttovIN%nlevels),inst_profiles(j)%o3(inst_profiles(j)%nlevels))
                do k=modeltop_index+1,inst_profiles(j)%nlevels-1 
                    kk = k + 1 - modeltop_index
                    call interpolate_logp(100*inst_profiles(j)%p(k),rttovIN%p(i,kk-1),rttovIN%p(i,kk),rttovIN%o3(i,kk-1),rttovIN%o3(i,kk),inst_profiles(j)%o3(k)) 
                end do
            end if
          else ! For when trace gas columns are supplied by the model. Units must match (kg/kg over moist air) and concentration must be supplied on model levels (not layers), requiring interpolation.
            if (Ldo_co2)  then ! CO2
                call interpolate_logp(100*inst_profiles(j)%p(modeltop_index),rttovIN%p(i,1),rttovIN%p(i,2),rttovIN%co2(i,1),rttovIN%co2(i,2),inst_profiles(j)%co2(modeltop_index))
                call interpolate_logp(100*inst_profiles(j)%p(inst_profiles(j)%nlevels),rttovIN%p(i,rttovIN%nlevels-1),rttovIN%p(i,rttovIN%nlevels),rttovIN%co2(i,rttovIN%nlevels-1),rttovIN%co2(i,rttovIN%nlevels),inst_profiles(j)%co2(inst_profiles(j)%nlevels))
                do k=modeltop_index+1,inst_profiles(j)%nlevels-1
                    kk = k + 1 - modeltop_index
                    call interpolate_logp(100*inst_profiles(j)%p(k),rttovIN%p(i,kk-1),rttovIN%p(i,kk),rttovIN%co2(i,kk-1),rttovIN%co2(i,kk),inst_profiles(j)%co2(k)) 
                end do
            end if
            if (Ldo_n2o)  then ! N2O
                call interpolate_logp(100*inst_profiles(j)%p(modeltop_index),rttovIN%p(i,1),rttovIN%p(i,2),rttovIN%n2o(i,1),rttovIN%n2o(i,2),inst_profiles(j)%n2o(modeltop_index))
                call interpolate_logp(100*inst_profiles(j)%p(inst_profiles(j)%nlevels),rttovIN%p(i,rttovIN%nlevels-1),rttovIN%p(i,rttovIN%nlevels),rttovIN%n2o(i,rttovIN%nlevels-1),rttovIN%n2o(i,rttovIN%nlevels),inst_profiles(j)%n2o(inst_profiles(j)%nlevels))
                do k=modeltop_index+1,inst_profiles(j)%nlevels-1
                    kk = k + 1 - modeltop_index
                    call interpolate_logp(100*inst_profiles(j)%p(k),rttovIN%p(i,kk-1),rttovIN%p(i,kk),rttovIN%n2o(i,kk-1),rttovIN%n2o(i,kk),inst_profiles(j)%n2o(k)) 
                end do
            end if
            if (Ldo_co)  then ! CO
                call interpolate_logp(100*inst_profiles(j)%p(modeltop_index),rttovIN%p(i,1),rttovIN%p(i,2),rttovIN%co(i,1),rttovIN%co(i,2),inst_profiles(j)%co(modeltop_index))
                call interpolate_logp(100*inst_profiles(j)%p(inst_profiles(j)%nlevels),rttovIN%p(i,rttovIN%nlevels-1),rttovIN%p(i,rttovIN%nlevels),rttovIN%co(i,rttovIN%nlevels-1),rttovIN%co(i,rttovIN%nlevels),inst_profiles(j)%co(inst_profiles(j)%nlevels))
                do k=modeltop_index+1,inst_profiles(j)%nlevels-1
                    kk = k + 1 - modeltop_index
                    call interpolate_logp(100*inst_profiles(j)%p(k),rttovIN%p(i,kk-1),rttovIN%p(i,kk),rttovIN%co(i,kk-1),rttovIN%co(i,kk),inst_profiles(j)%co(k)) 
                end do
            end if
            if (Ldo_ch4)  then ! CH4
                call interpolate_logp(100*inst_profiles(j)%p(modeltop_index),rttovIN%p(i,1),rttovIN%p(i,2),rttovIN%ch4(i,1),rttovIN%ch4(i,2),inst_profiles(j)%ch4(modeltop_index))
                call interpolate_logp(100*inst_profiles(j)%p(inst_profiles(j)%nlevels),rttovIN%p(i,rttovIN%nlevels-1),rttovIN%p(i,rttovIN%nlevels),rttovIN%ch4(i,rttovIN%nlevels-1),rttovIN%ch4(i,rttovIN%nlevels),inst_profiles(j)%ch4(inst_profiles(j)%nlevels))
                do k=modeltop_index+1,inst_profiles(j)%nlevels-1
                    kk = k + 1 - modeltop_index
                    call interpolate_logp(100*inst_profiles(j)%p(k),rttovIN%p(i,kk-1),rttovIN%p(i,kk),rttovIN%ch4(i,kk-1),rttovIN%ch4(i,kk),inst_profiles(j)%ch4(k)) 
                end do
            end if
            if (Ldo_so2)  then ! SO2
                call interpolate_logp(100*inst_profiles(j)%p(modeltop_index),rttovIN%p(i,1),rttovIN%p(i,2),rttovIN%so2(i,1),rttovIN%so2(i,2),inst_profiles(j)%so2(modeltop_index))
                call interpolate_logp(100*inst_profiles(j)%p(inst_profiles(j)%nlevels),rttovIN%p(i,rttovIN%nlevels-1),rttovIN%p(i,rttovIN%nlevels),rttovIN%so2(i,rttovIN%nlevels-1),rttovIN%so2(i,rttovIN%nlevels),inst_profiles(j)%so2(inst_profiles(j)%nlevels))
                do k=modeltop_index+1,inst_profiles(j)%nlevels-1
                    kk = k + 1 - modeltop_index
                    call interpolate_logp(100*inst_profiles(j)%p(k),rttovIN%p(i,kk-1),rttovIN%p(i,kk),rttovIN%so2(i,kk-1),rttovIN%so2(i,kk),inst_profiles(j)%so2(k)) 
                end do
            end if                                                
            if (Ldo_o3)  then ! Ozone
                call interpolate_logp(100*inst_profiles(j)%p(modeltop_index),rttovIN%p(i,1),rttovIN%p(i,2),rttovIN%o3(i,1),rttovIN%o3(i,2),inst_profiles(j)%o3(modeltop_index))
                call interpolate_logp(100*inst_profiles(j)%p(inst_profiles(j)%nlevels),rttovIN%p(i,rttovIN%nlevels-1),rttovIN%p(i,rttovIN%nlevels),rttovIN%o3(i,rttovIN%nlevels-1),rttovIN%o3(i,rttovIN%nlevels),inst_profiles(j)%o3(inst_profiles(j)%nlevels))
                do k=modeltop_index+1,inst_profiles(j)%nlevels-1
                    kk = k + 1 - modeltop_index
                    call interpolate_logp(100*inst_profiles(j)%p(k),rttovIN%p(i,kk-1),rttovIN%p(i,kk),rttovIN%o3(i,kk-1),rttovIN%o3(i,kk),inst_profiles(j)%o3(k)) 
                end do
            end if
          end if

          ! If adding CAM-like model top layer, copy temperature and trace gas values to the new top level. 
          if (inst_extend_atmos==1) then ! Replicate RRTMG-LW operations in CAM6
              inst_profiles(j)%t(1)    = inst_profiles(j)%t(modeltop_index)
              inst_profiles(j)%t(inst_profiles(j)%nlevels) = rttovIN%t_skin(i) ! RRTMG sets the lowest atmospheric temperature equal to the surface skin temperature
              inst_profiles(j)%q(1)    = inst_profiles(j)%q(modeltop_index)
              if (Ldo_co2) inst_profiles(j)%co2(1)  = inst_profiles(j)%co2(modeltop_index)
              if (Ldo_n2o) inst_profiles(j)%n2o(1)  = inst_profiles(j)%n2o(modeltop_index)
              if (Ldo_co)  inst_profiles(j)%co(1)   = inst_profiles(j)%co(modeltop_index)
              if (Ldo_ch4) inst_profiles(j)%ch4(1)  = inst_profiles(j)%ch4(modeltop_index)
              if (Ldo_so2) inst_profiles(j)%so2(1)  = inst_profiles(j)%so2(modeltop_index)
              if (Ldo_o3)  inst_profiles(j)%o3(1)   = 0.0
              ! CAM sets the top layer O3 to a linear interpolation between the highest layer and zero O3 at 50 Pa. So the top interface should just be zero then (see rrtmg_state.F90)
            !   if (Ldo_o3) call interpolate_logp(100*inst_profiles(j)%p(1),50._wp,rttovIN%p(i,1),0.0_wp,rttovIN%o3(i,1),inst_profiles(j)%o3(1))
          end if
          
          ! q and o3 coefficient limit is 0.1e-10 for MMRs over moist air
          if (inst_profiles(j)%gas_units .eq. 1) then
            where(inst_profiles(j)%q(:) < 0.1e-10)
                inst_profiles(j)%q(:) = 0.11e-10
            end where
            where(inst_profiles(j)%o3(:) < 0.1e-10)
                inst_profiles(j)%o3(:) = 0.11e-10
            end where
          end if

          ! 2m parameters
          inst_profiles(j)%s2m%p      =  rttovIN%p_surf(i) * 1e-2 ! convert Pa to hPa
          if (inst_extend_atmos==1) then
              inst_profiles(j)%s2m%t      =  rttovIN%t_skin(i)
          else
              inst_profiles(j)%s2m%t      =  rttovIN%t2m(i)
          end if
          inst_profiles(j)%s2m%q      =  rttovIN%q2m(i) ! Should be the same as gas units (kg/kg)
          inst_profiles(j)%s2m%u      =  rttovIN%u_surf(i)
          inst_profiles(j)%s2m%v      =  rttovIN%v_surf(i)
          inst_profiles(j)%s2m%wfetc  =  10000. ! only used by sea surface solar BRDF model.

          ! skin variables for emissivity calculations
          inst_profiles(j)%skin%t          =  rttovIN%t_skin(i)
          
          ! fastem coefficients - for mw calculations
          inst_profiles(j)%skin%fastem(1)  =  3.0
          inst_profiles(j)%skin%fastem(2)  =  5.0
          inst_profiles(j)%skin%fastem(3)  =  15.0
          inst_profiles(j)%skin%fastem(4)  =  0.1
          inst_profiles(j)%skin%fastem(5)  =  0.3

          ! Viewing angles
          inst_profiles(j)%zenangle      = inst_zenang ! pass in from cosp
          inst_profiles(j)%azangle       = 0. ! hard-coded in rttov9 int JKS-?

          inst_profiles(j)%latitude      = rttovIN%latitude(i)
          inst_profiles(j)%longitude     = rttovIN%longitude(i)
          inst_profiles(j)%elevation     = rttovIN%h_surf(i) * 1e-3 ! Convert m to km

          ! Solar angles. This is overwritten in the solar code now.
          if (associated(rttovIN%sza)) then
             inst_profiles(j)%sunzenangle   = rttovIN%sza(i) ! SZA in degrees
          else
             if (verbose) print*,'No solar zenith angle passed. Setting to zero.'
             inst_profiles(j)%sunzenangle   = 0.
          end if
          inst_profiles(j)%sunazangle    = 0. ! hard-coded in like rttov9

          ! surface type. sfcmask is 0 for ocean, 1 for land, and 2 for sea ice
          if (rttovIN%sfcmask(i) .lt. 0.5) then
             inst_profiles(j)%skin%surftype  = surftype_land
          else if (rttovIN%sfcmask(i) .lt. 1.5) then
             inst_profiles(j)%skin%surftype  = surftype_sea
          else
             inst_profiles(j)%skin%surftype  = surftype_seaice
          end if
          
          ! land-sea mask (lsmask) indicates proportion of land in grid (not in CESM implementation! just a binary mask there)
!          if (rttovIN%lsmask(i) < 0.5) then
!            inst_profiles(j)%skin%surftype  = surftype_sea
!          else
!            inst_profiles(j)%skin%surftype  = surftype_land
!          endif
          ! sea-ice fraction
!          if (rttovIN%icefrac(i) >= 0.5) then
!            inst_profiles(j)%skin%surftype  = surftype_seaice
!          endif

          ! dar: hard-coded to 1 (=ocean water) in rttov 9 int
          inst_profiles(j)%skin%watertype = 1
          !inst_profiles(j) %idg         = 0. ! Depreciated?
          !inst_profiles(j) %ish         = 0. ! Depreciated?

          ! Correct units if dry mass mixing ratios in ppmv were supplied
          !   Units for gas abundances: (must be the same for all profiles)
          !   2 => ppmv over moist air
          !   1=> kg/kg over moist air (default)
          !   0 (or less) => ppmv over dry air
          ! JKS added 3 => kg/kg over dry air, which requires conversion.
          if (inst_gas_units .eq. 3) then
            ! Convert to ppmv over dry air
            inst_profiles(j)%s2m%q     = (inst_profiles(j)%s2m%q / (1.0 - inst_profiles(j)%s2m%q)) * mdry / mH2O * 1.e6
            inst_profiles(j)%q(:)      = (inst_profiles(j)%q(:) / (1.0 - inst_profiles(j)%q(:))) * mdry / mH2O * 1.e6
            inst_profiles(j)%o3(:)     = inst_profiles(j)%o3(:) * mdry / mO3 * 1.e6
            inst_profiles(j)%co2(:)    = inst_profiles(j)%co2(:) * mdry / mCO2 * 1.e6
            inst_profiles(j)%ch4(:)    = inst_profiles(j)%ch4(:) * mdry / mCH4 * 1.e6
            inst_profiles(j)%n2o(:)    = inst_profiles(j)%n2o(:) * mdry / mN2O * 1.e6
            inst_profiles(j)%co(:)     = inst_profiles(j)%co(:) * mdry / mCO * 1.e6
            inst_profiles(j)%so2(:)    = inst_profiles(j)%so2(:) * mdry / mSO2 * 1.e6
            inst_profiles(j)%gas_units  =  0
            ! Alternately, convert kg/kg/ over dry air to kg/kg over moist air
            ! inst_profiles(j)%o3(:)     = inst_profiles(j)%o3(:) * (1.0 - inst_profiles(j)%q(:))
            ! inst_profiles(j)%co2(:)    = inst_profiles(j)%co2(:) * (1.0 - inst_profiles(j)%q(:))
            ! inst_profiles(j)%ch4(:)    = inst_profiles(j)%ch4(:) * (1.0 - inst_profiles(j)%q(:))
            ! inst_profiles(j)%n2o(:)    = inst_profiles(j)%n2o(:) * (1.0 - inst_profiles(j)%q(:))
            ! inst_profiles(j)%co(:)     = inst_profiles(j)%co(:) * (1.0 - inst_profiles(j)%q(:))
            ! inst_profiles(j)%so2(:)    = inst_profiles(j)%so2(:) * (1.0 - inst_profiles(j)%q(:))
            ! inst_profiles(j)%gas_units  =  1
          else
            inst_profiles(j)%gas_units  =  inst_gas_units
          end if

      end if 
    end do

    ! JKS - nothing to check here, this will never trigger.
    call rttov_error('error in profile initialization' , lalloc = .false.)
        
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Only add the cloud fields if simulating cloud.
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if (Lrttov_cld) then

      ! Set cloud mass mixing ratio units
      inst_profiles(:)%mmr_cldaer =  .true. ! kg/kg for cloud and aerosol (default)

      ! See RTTOV documentation page 75 for details.
      inst_profiles(:)%clw_scheme   = inst_clw_scheme ! 1: OPAC 2: Deff
      inst_profiles(:)%ice_scheme   = inst_ice_scheme ! 1:Baum 2:Baran(2014) 3:Baran(2018)
      if (inst_icede_param .eq. 0) then ! Need a filler even if icede is suppled
        inst_profiles(:)%icede_param = 2
      else ! Only use the model diameters if icede_param supplied on [1:4]
        inst_profiles(:)%icede_param = inst_icede_param ! 1: Ou and Liou, 2: Wyser(recommended), 3: Boudala, 4: McFarquhar. Only used if ice effective diameter not input
      end if
    !   inst_profiles(:)%clw_scheme   = 2 ! Deff scheme avoids cloud types but requires an effective diameter value
    ! !    inst_profiles(:)%clwde_scheme = 1. ! Scheme for cloud liquid water cotent to effective diameter. User guide says do not change.
    !   inst_profiles(:)%ice_scheme   = 1 !1:Baum 2:Baran(2014) 3:Baran(2018)
    !   inst_profiles(:)%icede_param  = 2 ! 2:Wyser(recommended). Only used if ice effective diameter not input
        
      j = 0 ! Initialize input
      do i = 1,rttovIN%nPoints
        if (i .gt. rttovIN%nPoints) exit
        if (inst_swath_mask(i)) then ! only added masked columns to profiles
            j = j + 1 ! Increment profile counter      
            
            ! Cloud scheme stuff. Values are on layers, not levels like the gas concentrations.
            inst_profiles(j)%cfrac(modeltop_index:inst_profiles(j)%nlayers)   = rttovIN%tca(i,:)    ! Cloud fraction for each layer       
            inst_profiles(j)%cloud(1,modeltop_index:inst_profiles(j)%nlayers) = rttovIN%cldLiq(i,:) ! Cloud water mixing ratio (all in the first type for Deff)
            inst_profiles(j)%cloud(6,modeltop_index:inst_profiles(j)%nlayers) = rttovIN%cldIce(i,:) ! Cloud ice mixing ratio (1 type). See pg 74.

            inst_profiles(j)%clwde(modeltop_index:inst_profiles(j)%nlayers)   = rttovIN%DeffLiq(i,:) ! Cloud water effective diameter (um)
            if (inst_icede_param .eq. 0) then ! Only use the model diameters if icede_param supplied
                inst_profiles(j)%icede(modeltop_index:inst_profiles(j)%nlayers)   = rttovIN%DeffIce(i,:) ! Cloud ice effective diameter (um)
            end if

        ! Example UKMO input has effective radii for multiple cloud types, making identification of a single
        ! liquid droplet or ice crystal effective diameter difficult.
        ! I opt to let RTTOV decide on the effective radius values, but more complex implementation
        ! could do a more thorough conversion between UKMO output and RTTOV input
    !    inst_profiles(j)%clwde = ! Cloud water effective diameter
    !    inst_profiles(j)%icede = ! Cloud ice effective diameter

        ! Old code for simple cloud schemes only
    !    inst_profiles(j)%cfraction  =  0.
    !    inst_profiles(j)%ctp        =  500.

        ! Other options not implemented
!        inst_profiles(j)%clw        = ! Cloud liquid water (kg/kg) – MW only,
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
      inst_profiles%mmr_cldaer =  .true. ! kg/kg for cloud and aerosol (default)       
      ! Read in aerosol profiles
!      j = 0 ! Initialize input
!      do i = 1,rttovIN%nPoints
!        if (i .gt. rttovIN%nPoints)
!        if (inst_swath_mask(i)) then ! only added masked columns to profiles
!            j = j + 1 ! Increment profile counter   
!            inst_profiles(j)%aerosols(naertyp,nlayers) = rttovIN%aerosols ! Aerosols in different modes (see User Guide pg 80)
!        end if
!      end do
    end if
        
    if (Lrttov_solar) then

      ! Populate longitude, latitude, time, and date profile fields
      ! Read in aerosol profiles
      j = 0 ! Initialize input
      do i = 1,rttovIN%nPoints
        if (i .gt. rttovIN%nPoints) exit
        if (inst_swath_mask(i)) then ! only added masked columns to profiles
            j = j + 1 ! Increment profile counter   
            inst_profiles(j)%date(:) = rttovIN%rttov_date(i,:)
            inst_profiles(j)%time(:) = rttovIN%rttov_time(i,:)
        end if
      end do

      ! Call functions to calculate the appropriate solar zenith and azimuthal angles.
      call RTTOV_CALC_SOLAR_ANGLES(errorstatus, inst_profiles)
      call rttov_error('Error when calling RTTOV_CALC_SOLAR_ANGLES', lalloc = .false.)
    
!      if (verbose) then
!        print*,'inst_profiles(:))%sunzenangle:    ',inst_profiles(:)%sunzenangle
!        print*,'inst_profiles(:))%sunazangle:     ',inst_profiles(:)%sunazangle      
!        print*,'inst_profiles(:))%zenangle:       ',inst_profiles(:)%zenangle      
!        print*,'inst_profiles(:))%azangle:        ',inst_profiles(:)%azangle      
!      end if
    end if

    ! if (verbose) then ! JKS remove at some point
    !     print*,"inst_profiles(1)%gas_units: ",inst_profiles(1)%gas_units
    !     print*,'inst_profiles(1)%skin%t:   ',inst_profiles(1)%skin%t
    !     print*,"inst_profiles(1)%s2m%t:    ",inst_profiles(1)%s2m%t
    !     print*,'inst_profiles(1)%p(:):     ',inst_profiles(1)%p(:)
    !     print*,'inst_profiles(1)%t(:):     ',inst_profiles(1)%t(:)
    !     print*,'inst_profiles(1)%q(:):     ',inst_profiles(1)%q(:)
    !     print*,'inst_profiles(1)%co2(:):   ',inst_profiles(1)%co2(:)
    !     print*,'inst_profiles(1)%ch4(:):     ',inst_profiles(1)%ch4(:)
    !     print*,'inst_profiles(1)%o3(:):     ',inst_profiles(1)%o3(:)
    !     print*,'inst_profiles(1)%n2o(:):    ',inst_profiles(1)%n2o(:)
    ! end if

    ! if (verbose) then
        ! print*,'inst_profiles(1)%nlevels:  ',inst_profiles(1)%nlevels
        ! print*,'inst_profiles(1)%nlayers:  ',inst_profiles(1)%nlayers        
    !     print*,'shape(rttovIN%t): ', shape(rttovIN%t)
    !     print*,'shape(rttovIN%p): ', shape(rttovIN%p)
    !     print*,'shape(rttovIN%ph): ', shape(rttovIN%ph)
    !     print*,'shape(inst_profiles(1)%p(:)):     ',shape(inst_profiles(1)%p(:))
    !     print*,'shape(inst_profiles(1)%t(:)):     ',shape(inst_profiles(1)%t(:))
    !     print*,'shape(inst_profiles(1)%q(:)):     ',shape(inst_profiles(1)%q(:))
    !     print*,'shape(inst_profiles(1)%cfrac(:)):     ',shape(inst_profiles(1)%cfrac(:))
    !     print*,'shape(rttovIN%tca(1,:)):  ',shape(rttovIN%tca(1,:))
        ! print*,'rttovIN%ph(1,:): ', rttovIN%ph(1,:)
        ! print*,'rttovIN%p(1,:): ', rttovIN%p(1,:)
        ! print*,'rttovIN%t(1,:): ', rttovIN%t(1,:)
        ! print*,'rttovIN%q(1,:): ', rttovIN%q(1,:)
        ! print*,'rttovIN%o3(1,:): ', rttovIN%o3(1,:)
        ! print*,'inst_profiles(1)%p(:):     ',inst_profiles(1)%p(:)
        ! print*,'inst_profiles(1)%t(:):     ',inst_profiles(1)%t(:)
        ! print*,'inst_profiles(1)%q(:):     ',inst_profiles(1)%q(:)
        ! print*,'inst_profiles(1)%o3(:):     ',inst_profiles(1)%o3(:)
        ! print*,'inst_profiles(1)%cfrac:    ',inst_profiles(1)%cfrac
    !     print*,'inst_profiles(1)%co2(:):   ',inst_profiles(1)%co2(:)
    !     print*,'inst_profiles(1)%skin%t:   ',inst_profiles(1)%skin%t
    !     print*,'inst_profiles(1)%s2m%t:    ',inst_profiles(1)%s2m%t
    !     print*,'inst_profiles(1)%s2m%p:    ',inst_profiles(1)%s2m%p
    ! end if    
    
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

  subroutine cosp_rttov_setup_emissivity_reflectance(inst_calcemis,      &
                                                     inst_emissivity,    &
                                                     inst_calcrefl,      &
                                                     inst_reflectance,   &
                                                     emis_grey)
 
    logical(kind=jplm),pointer,intent(inout) :: &
        inst_calcemis(:)
    type(rttov_emissivity),pointer,intent(inout) :: &
        inst_emissivity(:)
    logical(kind=jplm),pointer,intent(inout) :: &
        inst_calcrefl(:)
    type(rttov_reflectance),pointer,intent(inout) :: &
        inst_reflectance(:)
 
    real(kind=wp),intent(in),pointer,optional :: emis_grey
 
    ! In this example we have no values for input emissivities or reflectances
    ! so we initialise all inputs to zero
    call rttov_init_emis_refl(inst_emissivity, inst_reflectance)
    call rttov_error('error for emissivity/reflectance initialization' , lalloc = .true.)

    if (present(emis_grey)) inst_emissivity(:) % emis_in  = emis_grey ! Use greybody emissivity.

    ! Calculate emissivity within RTTOV where the input emissivity value is
    ! zero or less (all channels in this case)
    inst_calcemis(:) = (inst_emissivity(:) % emis_in <= 0._jprb)
    
    ! Calculate reflectances within RTTOV where the input BRDF value is zero or
    ! less (all channels in this case)
    inst_calcrefl(:) = (inst_reflectance(:) % refl_in <= 0._jprb)  
  
  end subroutine cosp_rttov_setup_emissivity_reflectance


  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! 6. rttov_setup_emissivity_reflectance - Specify surface emissivity for PC-RTTOV
  ! ------------------------------------------------------
  ! From RTTOV example files.
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine cosp_pc_rttov_setup_emissivity(inst_calcemis,    &
                                            inst_emissivity)

    logical(kind=jplm),pointer,intent(inout) :: &
        inst_calcemis(:)
    type(rttov_emissivity),pointer,intent(inout) :: &
        inst_emissivity(:)

    ! PC-RTTOV requires using RTTOV to calculate the surface emissivities.
    ! Reflectances are never calculated for hyper-spectral IR sounders
    call rttov_init_emis_refl(inst_emissivity)
    inst_calcemis(:) = .TRUE.
    call rttov_error('error for emissivity initialization (PC-RTTOV)' , lalloc = .true.)
  
  end subroutine cosp_pc_rttov_setup_emissivity

  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! 7. rttov_call_direct - Call RTTOV forward model (Woohoo!)
  ! ------------------------------------------------------
  ! From RTTOV example files.
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
  
  subroutine cosp_rttov_call_direct(inst_nthreads,      &
                                    inst_opts,          &
                                    inst_profiles,      &
                                    inst_coefs,         &
                                    inst_chanprof,      &
                                    inst_transmission,  &
                                    inst_radiance,      &
                                    inst_calcemis,      &
                                    inst_emissivity,    &
                                    inst_calcrefl,      &
                                    inst_reflectance,   &
                                    debug)
  
    integer(KIND=jpim),intent(in)   :: &
        inst_nthreads
    type(rttov_options),intent(in)  :: &
        inst_opts
    type(rttov_profile),pointer,intent(in) :: &
        inst_profiles(:)
    type(rttov_coefs),intent(in)    :: &
        inst_coefs
    type(rttov_chanprof),pointer,intent(in) :: &
        inst_chanprof(:)
    type(rttov_transmission),intent(inout) :: &
        inst_transmission
    type(rttov_radiance),intent(inout) :: &
        inst_radiance
    logical(kind=jplm),pointer,intent(in) :: &
        inst_calcemis(:)
    type(rttov_emissivity),pointer,intent(in) :: &
        inst_emissivity(:)
    logical(kind=jplm),pointer,intent(in) :: &
        inst_calcrefl(:)
    type(rttov_reflectance),pointer,intent(in) :: &
        inst_reflectance(:)        
    logical,intent(in),optional     :: &
        debug

    ! Local variables
    logical :: verbose = .false.
    
    if (present(debug)) verbose = debug
    
    if (verbose) then
        print*,'shape(inst_chanprof%prof):      ',shape(inst_chanprof%prof)
        print*,'shape(inst_chanprof%chan):      ',shape(inst_chanprof%chan)
        print*,'shape(inst_profiles):           ',shape(inst_profiles)
!    print*,'shape(inst_profiles(:)%q):         ',shape(inst_profiles(:)%q)
    end if

!    print*,'NTHRDS tests'
!    print*,'inst_profiles(:)%s2m%p:    ',inst_profiles(:)%s2m%p
!    print*,'inst_profiles(1)%p_surf:         ',inst_profiles(1)%p    
    
    if (inst_nthreads <= 1) then
      if (verbose) print*,'Calling rttov_direct'
      call rttov_direct(                &
              errorstatus,              &! out   error flag
              inst_chanprof,            &! in    channel and profile index structure
              inst_opts,                &! in    options structure
              inst_profiles,            &! in    profile array
              inst_coefs,               &! in    coefficients structure
              inst_transmission,        &! inout computed transmittances
              inst_radiance,            &! inout computed radiances
              calcemis    = inst_calcemis,   &! in    flag for internal emissivity calcs
              emissivity  = inst_emissivity, &! inout input/output emissivities per channel
              calcrefl    = inst_calcrefl,   &! in    flag for internal BRDF calcs
              reflectance = inst_reflectance) ! inout input/output BRDFs per channel
    else
      if (verbose) print*,'Calling rttov_parallel_direct'
      call rttov_parallel_direct(       &
              errorstatus,              &! out   error flag
              inst_chanprof,            &! in    channel and profile index structure
              inst_opts,                &! in    options structure
              inst_profiles,            &! in    profile array
              inst_coefs,               &! in    coefficients structure
              inst_transmission,        &! inout computed transmittances
              inst_radiance,            &! inout computed radiances
              calcemis    = inst_calcemis,   &! in    flag for internal emissivity calcs
              emissivity  = inst_emissivity, &! inout input/output emissivities per channel
              calcrefl    = inst_calcrefl,   &! in    flag for internal BRDF calcs
              reflectance = inst_reflectance,&! inout input/output BRDFs per channel
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
                                       inst_profiles,                &
                                       inst_coefs,                   &
                                       inst_chanprof,                &
                                       inst_transmission,            &
                                       inst_nchannels_rec,           &
                                       inst_channels_rec,            &
                                       inst_radiance,                &
                                       inst_calcemis,                &
                                       inst_emissivity,              &
                                       inst_pccomp,                  &
                                       debug)
  
    integer(KIND=jpim),intent(in)   :: &
        inst_nthreads
    type(rttov_options),intent(in) :: &
        inst_opts
    type(rttov_profile),pointer,intent(in) :: &
        inst_profiles(:)
    type(rttov_coefs),intent(in)   :: &
        inst_coefs
    type(rttov_chanprof),pointer,intent(in) :: &
        inst_chanprof(:)
    type(rttov_transmission),intent(inout) :: &
        inst_transmission        
    integer(jpim),intent(in) :: &
        inst_nchannels_rec
    integer(jpim),dimension(inst_nchannels_rec),intent(in)  :: &
        inst_channels_rec
    type(rttov_radiance),intent(inout) :: &
        inst_radiance
    logical(kind=jplm),pointer,intent(in) :: &
        inst_calcemis(:)        
    type(rttov_emissivity),pointer,intent(in) :: &
        inst_emissivity(:)       
    type(rttov_pccomp),intent(inout) :: &
        inst_pccomp     ! Output PC structure  
        
    logical,intent(in),optional :: &
        debug
    logical :: verbose = .false.
    
    if (present(debug)) verbose = debug
    
    if (inst_nthreads <= 1) then
      if (verbose) print*,'Calling rttov_direct (PC-RTTOV)'
      call rttov_direct(                 &
              errorstatus,               &! out   error flag
              inst_chanprof,             &! in    channel and profile index structure
              inst_opts,                 &! in    options structure
              inst_profiles,             &! in    profile array
              inst_coefs,                &! in    coefficients structure
              inst_transmission,         &! inout computed transmittances
              inst_radiance,             &! inout computed radiances
              calcemis     = inst_calcemis,   &! in    flag for internal emissivity calcs
              emissivity   = inst_emissivity, &! inout input/output emissivities per channel
              pccomp       = inst_pccomp,&! inout computed PC scores
              channels_rec = inst_channels_rec) ! in    reconstructed channel list
    else
      if (verbose) print*,'Calling rttov_parallel_direct (PC-RTTOV)'
      call rttov_parallel_direct(         &
              errorstatus,                &! out   error flag
              inst_chanprof,              &! in    channel and profile index structure
              inst_opts,                  &! in    options structure
              inst_profiles,              &! in    profile array
              inst_coefs,                 &! in    coefficients structure
              inst_transmission,          &! inout computed transmittances
              inst_radiance,              &! inout computed radiances
              calcemis     = inst_calcemis,    &! in    flag for internal emissivity calcs
              emissivity   = inst_emissivity,  &! inout input/output emissivities per channel
              pccomp       = inst_pccomp, &! inout computed PC scores
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
                                    Lrttov_cld,Lrttov_aer,inst_radiance,    &
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
    type(rttov_radiance),intent(in) :: &
        inst_radiance
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
                transpose(reshape(inst_radiance%bt(1:inst_nchan_out * nPoints), (/ inst_nchan_out, nPoints/) ))
        end if
        if (Lrttov_bt .and. (Lrttov_cld .or. Lrttov_aer)) then
            bt_clear(1:nPoints, 1:inst_nchan_out) = &
                transpose(reshape(inst_radiance%bt_clear(1:inst_nchan_out * nPoints), (/ inst_nchan_out, nPoints/) ))
        end if

        if (Lrttov_rad) then
            rad_total(1:nPoints, 1:inst_nchan_out) = &
                transpose(reshape(inst_radiance%total(1:inst_nchan_out * nPoints), (/ inst_nchan_out, nPoints/) ))
        end if
        if (Lrttov_rad .and. (Lrttov_cld .or. Lrttov_aer)) then
            rad_clear(1:nPoints, 1:inst_nchan_out) = &
                transpose(reshape(inst_radiance%clear(1:inst_nchan_out * nPoints), (/ inst_nchan_out, nPoints/) )) 
            rad_cloudy(1:nPoints, 1:inst_nchan_out) = &
                transpose(reshape(inst_radiance%cloudy(1:inst_nchan_out * nPoints), (/ inst_nchan_out, nPoints/) ))   
        end if

        if (Lrttov_refl) then
            refl_total(1:nPoints, 1:inst_nchan_out) = &
                transpose(reshape(inst_radiance%refl(1:inst_nchan_out * nPoints), (/ inst_nchan_out, nPoints/) ))
        end if
        if (Lrttov_refl .and. (Lrttov_cld .or. Lrttov_aer)) then
            refl_clear(1:nPoints, 1:inst_nchan_out) = &
                transpose(reshape(inst_radiance%refl_clear(1:inst_nchan_out * nPoints), (/ inst_nchan_out, nPoints/) ))
        end if
    else ! If swathing is occurring, assign the outputs appropriately
        j = 0
        do i=1,nPoints
        !   if (i .gt. nPoints) exit
          if (inst_swath_mask(i)) then ! only added masked columns to profiles
            if (Lrttov_bt) then
              bt_total(i, 1:inst_nchan_out) = inst_radiance%bt(1 + (j * inst_nchan_out):(j+1) * inst_nchan_out)
            end if
            if (Lrttov_bt .and. (Lrttov_cld .or. Lrttov_aer)) then
              bt_clear(i, 1:inst_nchan_out) = inst_radiance%bt_clear(1 + (j * inst_nchan_out):(j+1) * inst_nchan_out)
            end if            
            if (Lrttov_rad) then
              rad_total(i, 1:inst_nchan_out) = inst_radiance%total(1 + (j * inst_nchan_out):(j+1) * inst_nchan_out)
            end if
            if (Lrttov_rad .and. (Lrttov_cld .or. Lrttov_aer)) then
              rad_clear(i, 1:inst_nchan_out) = inst_radiance%clear(1 + (j * inst_nchan_out):(j+1) * inst_nchan_out)
              rad_cloudy(i, 1:inst_nchan_out) = inst_radiance%cloudy(1 + (j * inst_nchan_out):(j+1) * inst_nchan_out)
            end if
            if (Lrttov_refl) then
              refl_total(i, 1:inst_nchan_out) = inst_radiance%refl(1 + (j * inst_nchan_out):(j+1) * inst_nchan_out)
            end if
            if (Lrttov_refl .and. (Lrttov_cld .or. Lrttov_aer)) then
              refl_clear(i, 1:inst_nchan_out) = inst_radiance%refl_clear(1 + (j * inst_nchan_out):(j+1) * inst_nchan_out)
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
                                       inst_pccomp,          &
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
    type(rttov_pccomp),intent(in) :: &
        inst_pccomp     ! Output PC structure          
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
        print*,'size(inst_pccomp%bt_pccomp):   ',size(inst_pccomp%bt_pccomp)
        print*,'size(inst_pccomp%total_pccomp):   ',size(inst_pccomp%total_pccomp)
        print*,'inst_nchannels_rec * nPoints:   ',inst_nchannels_rec * nPoints
    end if

    ! Documentation for RTTOV radiance structure in RTTOV User Guide pg 166
        
    ! Only save output if appropriate
    if (count(inst_swath_mask) .eq. nPoints) then ! No swathing, save all output  
      if (Lrttov_bt) then
          bt_clear(1:nPoints, 1:inst_nchannels_rec) = &
              transpose(reshape(inst_pccomp%bt_pccomp(1:(inst_nchannels_rec * nPoints)), (/ inst_nchannels_rec, nPoints/) ))
      end if
      if (Lrttov_rad) then
          rad_clear(1:nPoints, 1:inst_nchannels_rec) = &
              transpose(reshape(inst_pccomp%total_pccomp(1:(inst_nchannels_rec * nPoints)), (/ inst_nchannels_rec, nPoints/) ))
      end if
    else ! If swathing is occurring, assign the outputs appropriately
      j = 0
      do i=1,nPoints
        ! if (i .gt. nPoints) exit
        if (inst_swath_mask(i)) then ! only added masked columns to profiles
          if (Lrttov_bt) then
            bt_clear(i, 1:inst_nchannels_rec) = inst_pccomp%bt_pccomp(1 + (j * inst_nchannels_rec):(j+1) * inst_nchannels_rec)
          end if          
          if (Lrttov_rad) then
            rad_clear(i, 1:inst_nchannels_rec) = inst_pccomp%total_pccomp(1 + (j * inst_nchannels_rec):(j+1) * inst_nchannels_rec)
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
                                            inst_profiles,        &
                                            inst_coefs,           &
                                            inst_chanprof,        &
                                            inst_transmission,    &
                                            inst_radiance,        &
                                            inst_calcemis,        &
                                            inst_emissivity,      &
                                            inst_calcrefl,        &
                                            inst_reflectance)

    integer(kind=jpim),intent(in) :: &
        inst_nprof,        &
        inst_nchanprof,    &
        nLevels
    type(rttov_options),intent(in) :: &
        inst_opts
    type(rttov_profile),pointer,intent(in) :: &
        inst_profiles(:)
    type(rttov_coefs),intent(in)   :: &
        inst_coefs
    type(rttov_chanprof),pointer,intent(inout) :: &
        inst_chanprof(:)
    type(rttov_transmission),intent(inout) :: &
        inst_transmission
    type(rttov_radiance),intent(inout) :: &
        inst_radiance
    logical(kind=jplm),pointer,intent(inout) :: &
        inst_calcemis(:)
    type(rttov_emissivity),pointer,intent(inout) :: &
        inst_emissivity(:)
    logical(kind=jplm),pointer,intent(inout) :: &
        inst_calcrefl(:)
    type(rttov_reflectance),pointer,intent(inout) :: &
        inst_reflectance(:)          

    ! Deallocate structures for rttov_direct
    call rttov_alloc_direct(     &
        errorstatus,             &
        0_jpim,                  &  ! 0 => deallocate
        inst_nprof,              &
        inst_nchanprof,          & 
        nLevels+1,               & ! "levels" means interfaces, not layers
        ! nLevels,                 &
        inst_chanprof,           &
        inst_opts,               &
        inst_profiles,           &
        inst_coefs,              &
        inst_transmission,       &
        inst_radiance,           &
        calcemis=inst_calcemis,       &
        emissivity=inst_emissivity,   &
        calcrefl=inst_calcrefl,       &
        reflectance=inst_reflectance)
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
                                               inst_profiles,         &
                                               inst_coefs,            &
                                               inst_chanprof,         &
                                               inst_transmission,     &
                                               inst_radiance,         &
                                               inst_calcemis,         &
                                               inst_emissivity,       &
                                               inst_pccomp,           &
                                               inst_predictindex)
 
    integer(kind=jpim),intent(in)  :: &
        inst_nprof,      &
        inst_nchanprof,  &
        nlevels,         &
        inst_nChannels_rec
    type(rttov_options),intent(in) :: &
        inst_opts
    type(rttov_profile),pointer,intent(in) :: &
        inst_profiles(:)
    type(rttov_coefs),intent(in)   :: &
        inst_coefs
    type(rttov_chanprof),pointer,intent(inout) :: &
        inst_chanprof(:)
    type(rttov_transmission),intent(inout) :: &
        inst_transmission        
    type(rttov_radiance),intent(inout) :: &
        inst_radiance
    logical(kind=jplm),pointer,intent(inout) :: &
        inst_calcemis(:)
    type(rttov_emissivity),pointer,intent(inout) :: &
        inst_emissivity(:)       
    type(rttov_pccomp),intent(inout) :: &
        inst_pccomp     ! Output PC structure  
    integer(kind=jpim),pointer,intent(inout) :: &
        inst_predictindex(:)        
        
    if (ASSOCIATED(inst_predictindex)) deallocate (inst_predictindex, stat=alloc_status(10))
    call rttov_error('mem dellocation error for "inst_predictindex"', lalloc = .true.)
    
    ! Deallocate structures for rttov_direct
    call rttov_alloc_direct(     &
        errorstatus,             &
        0_jpim,                  &  ! 0 => deallocate
        inst_nprof,              &
        inst_nchanprof,          &
        nLevels+1,               & ! "levels" means interfaces, not layers
        ! nLevels,                 &
        inst_chanprof,           &
        inst_opts,               &
        inst_profiles,           &
        inst_coefs,              &
        inst_transmission,       &
        inst_radiance,           &
        calcemis=inst_calcemis,       &
        emissivity=inst_emissivity,   &
        npcscores=inst_opts%rt_ir%pc%npcscores * inst_nprof,         &
        nchannels_rec=inst_nChannels_rec * inst_nprof, &
        pccomp=inst_pccomp)
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

  subroutine interpolate_logp(ptarget,p1,p2,z1,z2,ztarget)

    real(jprb),intent(in) :: ptarget
    real(WP),intent(in) :: p1,p2,z1,z2
    real(jprb),intent(out) :: ztarget ! variable interpolated to the target pressure level

    ! Normal procedure where ptarget falls within [p1,p2]
    if ((ptarget .gt. p1) .and. (ptarget .lt. p2)) then 
        ztarget = z1 + (z2-z1) * log(ptarget / p1) / log(p2 / p1)
    elseif (ptarget .lt. p1) then ! Top of model level. ptarget may be zero...
        ztarget = z1 ! Just set it to the layer value?? Not sure how to handle this if ptarget=0. I think that this is fine. We're basically out of range.
    elseif (ptarget .gt. p2) then ! surface level
        ztarget = z2 + (z2-z1) * log(ptarget / p2) / log(p2 / p1)
    end if
        
  end subroutine interpolate_logp

  subroutine get_DOY(month, day, DOY)

    integer,intent(in) :: &
        month,   &
        day
    integer,intent(out) :: &
        DOY

    ! This subroutine does not handle leap years because it is not relevant to the purpose.
    ! Simple look-up table for DOY.
    if (month .eq. 1)  DOY = day
    if (month .eq. 2)  DOY = 31 + day
    if (month .eq. 3)  DOY = 59 + day
    if (month .eq. 4)  DOY = 90 + day
    if (month .eq. 5)  DOY = 120 + day
    if (month .eq. 6)  DOY = 151 + day
    if (month .eq. 7)  DOY = 181 + day
    if (month .eq. 8)  DOY = 212 + day
    if (month .eq. 9)  DOY = 243 + day
    if (month .eq. 10) DOY = 273 + day
    if (month .eq. 11) DOY = 304 + day
    if (month .eq. 12) DOY = 334 + day

  end subroutine get_DOY

!##########################
! Module End
!##########################
end module mod_cosp_rttov