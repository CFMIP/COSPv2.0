! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copyright (c) 2015, Regents of the University of Colorado
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
! May 2015 - D. Swales - Original version
! Apr 2015 - D. Swales - Modified for RTTOVv11.3
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE MOD_COSP_RTTOV_INTERFACE
  USE COSP_KINDS,       ONLY: wp
  USE MOD_COSP_CONFIG,  ONLY: RTTOV_MAX_CHANNELS,rttovDir
  use mod_cosp_rttov,   only: platform,satellite,sensor,nChannels,iChannel,coef_rttov,   &
                              opts,construct_rttov_coeffilename,rttov_in,                &
                              construct_rttov_scatfilename
!                              coef_scatt,opts,opts_scatt,construct_rttov_coeffilename,   &
                              
  ! rttov_const contains useful RTTOV constants
  USE rttov_const, ONLY :     &
         errorstatus_success, &
         errorstatus_fatal,   &
         platform_name,       &
         inst_name,           &
         surftype_sea,        &
         surftype_land,       &
         surftype_seaice
  ! rttov_types contains definitions of all RTTOV data types
  USE rttov_types, ONLY :     &
         rttov_options,       &
         rttov_coefs,         &
         rttov_profile,       &
         rttov_transmission,  &
         rttov_radiance,      &
         rttov_chanprof,      &
         rttov_emissivity,    &
         rttov_reflectance

  ! jpim, jprb and jplm are the RTTOV integer, real and logical KINDs
  USE parkind1, ONLY : jpim, jprb, jplm
  
  USE rttov_unix_env, ONLY : rttov_exit
                              
  IMPLICIT NONE

#include "rttov_direct.interface"
#include "rttov_parallel_direct.interface"
#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_alloc_direct.interface"
#include "rttov_init_emis_refl.interface"
#include "rttov_user_options_checkinput.interface"
#include "rttov_print_opts.interface"
#include "rttov_print_profile.interface"
#include "rttov_skipcommentline.interface"

  !--------------------------
  !
  INTEGER(KIND=jpim), PARAMETER :: iup   = 20   ! unit for input profile file
  INTEGER(KIND=jpim), PARAMETER :: ioout = 21   ! unit for output

  ! RTTOV variables/structures
  !====================
! JKS these are already imported from "mod_cosp_rttov"
!  TYPE(rttov_options)              :: opts                     ! Options structure
!  TYPE(rttov_coefs)                :: coefs                    ! Coefficients structure

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
!  integer :: alloc_status_cosp(60)
  

! Old
!#include "rttov_read_coefs.interface"
!#include "rttov_read_scattcoeffs.interface"

CONTAINS

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_rttov_init
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_RTTOV_INIT(NchanIN,platformIN,satelliteIN,instrumentIN,channelsIN,nlevels)
    integer,intent(in) :: & 
         NchanIN,     & ! Number of channels
         platformIN,  & ! Satellite platform
         satelliteIN, & ! Satellite
         instrumentIN, & ! Instrument
         nlevels
    integer,intent(in),dimension(RTTOV_MAX_CHANNELS) :: &
         channelsIN     ! RTTOV channels

    ! Local variables
    character(len=256) :: &
        coef_file, &
        scat_file, &
        rttov_coefDir, &
        rttov_predDir, &
        rttov_cldaerDir, &
        OD_coef_filepath, &
        aer_coef_file, &
        cld_coef_file, &
        aer_coef_filepath, &
        cld_coef_filepath

    ! Initialize fields in module memory (cosp_rttovXX.F90)
    nChannels  = NchanIN
    platform   = platformIN 
    satellite  = satelliteIN 
    sensor     = instrumentIN 
    iChannel   = channelsIN

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 1. Initialise RTTOV options structure
    ! ------------------------------------------------------
    ! See page 157 of RTTOV v13 user guide for documentation
    ! Initializing all options to defaults for consistency
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ! General configuration options
    opts%config%do_checkinput    = .true.
    opts%config%apply_reg_limits = .false. ! True in v11
    opts%config%verbose          = .true.
    opts%config%opdep13_gas_clip = .true.
    
    ! General Radiative Transfer Options
    ! Gas profiles
    opts%rt_all%ozone_data        = .true.
    opts%rt_all%so2_data          = .true.
    
    ! Well-mixed gases
    opts%rt_all%co2_data          = .true.
    opts%rt_all%n2o_data          = .true.
    opts%rt_all%co_data           = .true.
    opts%rt_all%ch4_data          = .true.
    
    ! Other general RT options (initializing to defaults for completeness)
    opts%rt_all%do_lambertian       = .false.
    opts%rt_all%switchrad           = .false.
    opts%rt_all%rad_down_lin_tau    = .true.
    opts%rt_all%use_t2m_opdep       = .true.
    opts%rt_all%use_q2m             = .true.
    opts%rt_all%use_tskin_eff       = .false.
    opts%rt_all%addrefrac           = .true.
    opts%rt_all%plane_parallel      = .false.
    opts%rt_all%transmittances_only = .false.
    
    ! MW-only radiative transfer options:
    ! JKS make this optional?
    opts%rt_mw%clw_data             = .false.
    opts%rt_mw%clw_scheme           = 2       ! Default = 2/Rosenkranz
    opts%rt_mw%clw_cloud_top        = 322     ! Default is 322 hPa
    opts%rt_mw%fastem_version       = 6       ! Default FASTEM-6
    opts%rt_mw%supply_foam_fraction = .false.

    ! UV/visible/IR-only radiative transfer options
    opts%rt_ir%addsolar                = .false.
    opts%rt_ir%rayleigh_max_wavelength = 2._wp ! 2um 
    opts%rt_ir%rayleigh_min_pressure   = 0._wp ! 0 hPa
    opts%rt_ir%rayleigh_single_scatt   = .true.
    opts%rt_ir%rayleigh_depol          = .true. ! Default false, recommended true
    opts%rt_ir%do_nlte_correction      = .false.
    opts%rt_ir%solar_sea_brdf_model    = 2
    opts%rt_ir%ir_sea_emis_model       = 2
    opts%rt_ir%addaerosl               = .false.
    opts%rt_ir%addclouds               = .false.
    opts%rt_ir%user_aer_opt_param      = .false. ! User specifies the aerosol scattering optical parameters 
    opts%rt_ir%user_cld_opt_param      = .false. ! User specifies the cloud scattering optical parameters 
    opts%rt_ir%grid_box_avg_cloud      = .true.
    opts%rt_ir%cldcol_threshold        = -1._wp
    opts%rt_ir%cloud_overlap           = 1 ! Maximum-random overlap
    opts%rt_ir%cc_low_cloud_top        = 750_wp ! 750 hPa. Only applies when cloud_overlap=2.
    opts%rt_ir%ir_scatt_model          = 2
    opts%rt_ir%vis_scatt_model         = 1
    opts%rt_ir%dom_nstreams            = 8
    opts%rt_ir%dom_accuracy            = 0._wp ! only applies when addclouds or addaerosl is true and DOM is selected as a scattering solver.
    opts%rt_ir%dom_opdep_threshold     = 0._wp
    opts%rt_ir%dom_rayleigh            = .false.
    
    ! Principal Components-only radiative transfer options:
    opts%rt_ir%pc%addpc     = .false.
    opts%rt_ir%pc%npcscores = -1
    opts%rt_ir%pc%addradrec = .false.
    opts%rt_ir%pc%ipcbnd    = 1
    opts%rt_ir%pc%ipcreg    = 1 ! The index of the required set of PC predictors
    
    ! Options related to interpolation and the vertical grid:
    opts%interpolation%addinterp         = .true.
    opts%interpolation%interp_mode       = 1
!    opts%interpolation%reg_limit_extrap = .true. ! Depreciated
    opts%interpolation%lgradp            = .false.
!    opts%interpolation%spacetop         = .true. ! Depreciated

    ! Options related to HTFRTC:
    opts%htfrtc_opts%htfrtc       = .false.
    opts%htfrtc_opts%n_pc_in      = -1
    opts%htfrtc_opts%reconstruct  = .false.
    opts%htfrtc_opts%simple_cloud = .false.
    opts%htfrtc_opts%overcast     = .false.
    
    ! Developer options that may be useful:
    opts%dev%do_opdep_calc = .true.
    
    ! JKS To-do: include opts_scatt settings (user guide pg 161)


    ! Old configured options

    ! Options common to RTTOV clear-sky Tb calculation
!    opts%interpolation%addinterp  = .true.  ! allow interpolation of input profile
!    opts%rt_all%use_q2m           = .true.
!    opts%config%do_checkinput     = .false.
!    opts%config%verbose           = .false.
!    opts%rt_all%addrefrac         = .true.  ! include refraction in path calc
!    opts%interpolation%reg_limit_extrap = .true.
    
!    opts%rt_mw%clw_data                = .true. 
    ! Options common to RTTOV clear-sky Tb calculation
    
    ! These scattering options are depreciated in v13
!    opts_scatt%config%do_checkinput    = .false.
!    opts_scatt%config%verbose          = .false.
!    opts_scatt%config%apply_reg_limits = .true.
!    opts_scatt%interp_mode             = 1
!    opts_scatt%reg_limit_extrap        = .true.
!    opts_scatt%use_q2m                 = .true.


    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 2. Read coefficients (from RTTOV example files)
    ! ------------------------------------------------------
    ! Using the GUI to figure out files that work together could be helpful here.
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Construct optical depth and cloud coefficient files
    
    ! rttovDir should be "/glade/u/home/jonahshaw/w/RTTOV/" passed from namelist
    
    ! Hardcoding these other paths for now, they should be input later.
    rttov_coefDir   = "rtcoef_rttov13/" ! directory for coefficient in RTTOV v13
    rttov_predDir   = "rttov13pred54L/" ! example directory for predictors. This should be input
    rttov_cldaerDir = "cldaer_visir/" ! This should be input. Also "cldaer_ir".
    
    ! Optical depth file
    OD_coef_filepath = trim(rttovDir)//trim(rttov_coefDir)//trim(rttov_predDir)// &
                       trim(construct_rttov_coeffilename(platform,satellite,sensor))

    ! Example coefficient files (hardcoded)
    aer_coef_file = "scaercoef_eos_2_airs_cams_chou-only.H5"
    cld_coef_file = "sccldcoef_eos_2_airs_chou-only.H5"
    
    ! Coefficient files from the "construct_rttov_scatfilename" function. Not sure if working.
!    aer_coef_file = construct_rttov_scatfilename(platform,satellite,sensor)
!    cld_coef_file = construct_rttov_scatfilename(platform,satellite,sensor)
    
    ! Cloud and Aerosol scattering (and absorption?) file(s)
    aer_coef_filepath = trim(rttovDir)//trim(rttov_coefDir)//trim(rttov_cldaerDir)// &
                        trim(aer_coef_file)
    cld_coef_filepath = trim(rttovDir)//trim(rttov_coefDir)//trim(rttov_cldaerDir)// &
                        trim(cld_coef_file)
         
    ! Read optical depth and cloud coefficient files together
    call rttov_read_coefs(errorstatus, coef_rttov, opts,    &
                          file_coef=OD_coef_filepath,       &
                          file_scaer=aer_coef_filepath,     &
                          file_sccld=cld_coef_filepath)
                          
    ! We aren't checking an allocation steps so this seems more appropriate.
    call rttov_error('fatal error reading coefficients' , lalloc = .false.)
     ! Old error check
!    if (errorstatus /= errorstatus_success) then
!        write(*,*) 'fatal error reading coefficients'
!        call rttov_exit(errorstatus)
!    endif
    
    ! Ensure input number of channels is not higher than number stored in coefficient file
    if (nchannels > coef_rttov % coef % fmv_chn) then
        nchannels = coef_rttov % coef % fmv_chn
    endif

    ! Ensure the options and coefficients are consistent
    call rttov_user_options_checkinput(errorstatus, opts, coef_rttov)
    
    ! We aren't checking an allocation steps so this seems more appropriate.
    call rttov_error('error in rttov options' , lalloc = .false.) 
    ! Old error check
!    if (errorstatus /= errorstatus_success) then
!        write(*,*) 'error in rttov options'
!        call rttov_exit(errorstatus)
!    endif
    
    ! Old code
    !    coef_file = trim(rttovDir)//"rtcoef_rttov11/rttov7pred54L/"// &
    !         trim(construct_rttov_coeffilename(platform,satellite,sensor))
    
        ! Read in scattering (clouds+aerosol) coefficient file. *ONLY NEEDED IF DOING RTTOV ALL-SKY.*
    !scat_file = trim(rttovDir)//"rtcoef_rttov11/cldaer/"//&
    !     trim(construct_rttov_scatfilename(platform,satellite,sensor))
    ! Can't pass filename to rttov_read_scattcoeffs!!!!!
    !call rttov_read_scattcoeffs (errorstatus, coef_rttov%coef, coef_scatt,)

!!!! Old version 

    ! Read in scattering coefficient file.
!    coef_file = trim(rttovDir)//"rtcoef_rttov11/rttov7pred54L/"// &
!         trim(construct_rttov_coeffilename(platform,satellite,sensor))
!    call rttov_read_coefs(errorstatus,coef_rttov, opts, file_coef=trim(coef_file))

    ! Read in scattering (clouds+aerosol) coefficient file. *ONLY NEEDED IF DOING RTTOV ALL-SKY.*
    !scat_file = trim(rttovDir)//"rtcoef_rttov11/cldaer/"//&
    !     trim(construct_rttov_scatfilename(platform,satellite,sensor))
    ! Can't pass filename to rttov_read_scattcoeffs!!!!!
    !call rttov_read_scattcoeffs (errorstatus, coef_rttov%coef, coef_scatt,)
    
    
    ! subsub routines
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
 
  END SUBROUTINE COSP_RTTOV_INIT
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_rttov_simulate - JKS should this move to cosp_rttov_v13?
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_RTTOV_SIMULATE(rttovIN,         & ! Inputs
                                 Tb,error)                    ! Outputs
  
    type(rttov_in),intent(in) :: &
        rttovIN
    real(wp),dimension(rttovIN%nPoints,rttovIN%nChannels) :: & ! Can I do this? I guess so!
        Tb        ! RTTOV brightness temperature.
    character(len=128) :: &
        error     ! Error messages (only populated if error encountered)  
    logical                             :: &
        rttov_simulate_cld,                &
        rttov_simulate_aer


    ! variables for input
    !====================
    integer(kind=jpim) :: nthreads
    integer(kind=jpim) :: dosolar
    integer(kind=jpim) :: nchanprof ! JKS - jpim is RTTOV integer type
    integer(kind=jpim), allocatable :: channel_list(:) ! JKS this needs to be specified
      
    ! Loop variables
    integer(kind=jpim) :: i, j, jch
    integer(kind=jpim) :: nch
    
    integer :: errorstatus
    integer :: ich, ich_temp, chan
    
    
    ! Initialize some things
    Tb(:,:)    = 0._wp
    error      = ''
      
  ! How do I want the interface to function? How should it to be consistent with the rest of COSP?

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 3. Allocate RTTOV input and output structures
    ! ------------------------------------------------------
    ! Largely from RTTOV documentation.
    ! This only needs to be performed once if all chunks are the same size, right? <-- We don't know chunk sizes at init, so this must be in the main loop.
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
      ! Determine the total number of radiances to simulate (nchanprof).
      ! We aren't doing subcolumn sampling (RTTOV already does this and it would be slow)
!      nchanprof = nchannels * nprof
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
    
    ! COSP v11 code has an example for multiple instruments
  
    nch = 0_jpim
    do j = 1, rttovIN%nPoints
      do jch = 1, rttovIN%nChannels
        nch = nch + 1_jpim
        chanprof(nch)%prof = j
        chanprof(nch)%chan = channel_list(jch) ! channel_list should be specified
      end do
    end do
      
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 5. Read profile data
    ! ------------------------------------------------------
    ! Largely from cosp_rttov_v11.F90 file.
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ! Create a new subroutine in cosp_rttov_v13 and call it to return the profiles
    
    !call rttov_profile_setup(rttovIN,profiles)
    
    ! Store profile data from rttovIN in profile type.
    ! See RTTOV user guide pg 163 for description of "profiles" type
    
! "The rttov_profile structure is composed of the atmospheric part 
! and two other structures for 2 meter air and skin surface. 
! If you are not able to provide ozone, CO2, etc profiles the flags
! ozone_data, co2_data and so on in the options structure should be 
! set to false."

    ! Want to set these up as arrays of size (nLevels,nPoints) with
    ! constant values. Not sure how to do it efficiently. i.e. no loops
    
    ! Initialize trace gas column concentrations (well-mixed so constant in input)
    do j = 1, rttovIN%nlevels
      profiles(:)%co2(j)        =  rttovIN%co2
      profiles(:)%n2o(j)        =  rttovIN%n2o
      profiles(:)%co(j)         =  rttovIN%n2o
      profiles(:)%ch4(j)        =  rttovIN%ch4
    end do

    profiles%gas_units  =  1 ! kg/kg over moist air (default)
    
    do i = 1, rttovIN%nPoints
      profiles(i)%p(:) =  rttovIN%p(i, :)
      profiles(i)%t(:) =  rttovIN%t(i, :)
      profiles(i)%q(:) =  rttovIN%q(i, :)

      ! JKS not sure if I should keep this
      where(profiles(i)%q(:) < 1e-4)
        profiles(i)%q(:) = 1e-4
      end where


       ! Add handling of more complex cloud schemes

      ! Gas profiles
      profiles(i)%o3         =  rttovIN%o3(i, :)
!      profiles(i)%so2        =  ! Sulfate not in COSP input files
       

!      profiles(i)%aerosols(naertyp,nlayers) = ! Aerosols in different modes (see User Guide pg 80)



      ! 2m parameters
      profiles(i)%s2m%p      =  rttovIN%p_surf(i)
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
      profiles(i)%elevation     = rttovIN%h_surf(i)

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
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Only add the cloud fields if simulating cloud.
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rttov_simulate_cld = .false.
    if (rttov_simulate_cld) then

      ! Set cloud mass mixing ratio units
      profiles%mmr_cldaer =  .true. ! kg/kg for cloud and aerosol (default)

      profiles%clw_scheme   = 2. ! Deff scheme avoids cloud types
    !    profiles%clwde_scheme = 1. ! Not implemented?
      profiles%ice_scheme   = 1. !1:Baum 2:Baran(2014) 3:Baran(2018)
      profiles%icede_param  = 2. ! 2:Wyser(recommended). Only used if ice effective diameter not input
        
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
    rttov_simulate_aer = .false.
    if (rttov_simulate_aer) then
    
      ! Set aerosol mass mixing ratio units
      profiles%mmr_cldaer =  .true. ! kg/kg for cloud and aerosol (default)
        
      ! Read in aerosol profiles
!      do i = 1, rttovIN%nPoints
!
!      end do
    endif
    
    ! JKS - nothing to check here, this will never trigger.
    call rttov_error('error in aerosol profile initialization' , lalloc = .true.)
    
    ! JKS To-do: set up scattering profiles (MW only) (rttov_profile_cloud)
    
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 6. Specify surface emissivity and reflectance
    ! ------------------------------------------------------
    ! From RTTOV example files.
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
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
    
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 7. Call RTTOV forward model (Woohoo!)
    ! ------------------------------------------------------
    ! From RTTOV example files.
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    nthreads = 1 ! Default not parallel for now. Can be optimized later. - JKS
    
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
    
!    if (errorstatus /= errorstatus_success) then
!      write (*,*) 'rttov_direct error'
!      call rttov_exit(errorstatus)
!    endif
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 8. Save output data
    ! ------------------------------------------------------
    ! 
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Tb(1:rttovIN%nPoints, 1:rttovIN%nChannels) = &
        transpose(reshape(radiance%bt(1:nchanprof), (/ rttovIN%nChannels, rttovIN%nPoints/) ))
    
!    tbs(1:prof_num , ich_temp:ich_temp + size(ichan(:,no_id)) - 1) = &
!        transpose(reshape(radiance%bt(1:nchanprof), (/ size(ichan(:,no_id)), prof_num/) ))

!    ich_temp = ich_temp + size(ichan(:,no_id))
    
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 9. Deallocate all RTTOV arrays and structures
    ! ------------------------------------------------------
    ! From RTTOV example files.
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    deallocate (channel_list, stat=alloc_status(1))
    if (alloc_status(1) /= 0) then
      write(*,*) 'mem dellocation error'
    endif

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
    
  
    ! subsub routines
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

  END SUBROUTINE COSP_RTTOV_SIMULATE
  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! END MODULE
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE MOD_COSP_RTTOV_INTERFACE
