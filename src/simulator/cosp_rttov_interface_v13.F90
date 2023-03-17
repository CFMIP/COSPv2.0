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
                              construct_rttov_scatfilename,do_rttov_cld,do_rttov_aer,    &
                              do_rttov_rad,rttov_cld_optparam,rttov_aer_optparam
                              
  ! rttov_const contains useful RTTOV constants
  USE rttov_const, ONLY :     &
         errorstatus_success, &
         errorstatus_fatal

  ! rttov_types contains definitions of all RTTOV data types
  USE rttov_types, ONLY :     &
         rttov_options,       &
         rttov_coefs

  ! jpim, jprb and jplm are the RTTOV integer, real and logical KINDs
  USE parkind1, ONLY : jpim, jprb, jplm
  
  USE rttov_unix_env, ONLY : rttov_exit
                              
  IMPLICIT NONE

#include "rttov_read_coefs.interface"
#include "rttov_user_options_checkinput.interface"
#include "rttov_print_opts.interface"
  
  ! RTTOV variables/structures
  !====================
  LOGICAL(KIND=jplm),      POINTER :: calcemis(:)    => NULL() ! Flag to indicate calculation of emissivity within RTTOV
  LOGICAL(KIND=jplm),      POINTER :: calcrefl(:)    => NULL() ! Flag to indicate calculation of BRDF within RTTOV
  INTEGER(KIND=jpim)               :: errorstatus              ! Return error status of RTTOV subroutine calls

  INTEGER(KIND=jpim) :: alloc_status(60)
  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! TYPE rttov_init_IN (RTTOV init DDT to be passed to cosp_init)
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
! I may remove this because it will require an additional dependency between
! cosp2_test and the RTTOV interface.

! Integers: NchanIN, platformIN, satelliteIN, instrumentIN, channelsIN,                       &
! Logicals: Lrttov_cld, Lrttov_aer, Lrttov_rad, Lrttov_cldparam, Lrttov_aerparam
  
  type rttov_init_IN
     logical,pointer :: &
          Lrttov_cld,      &
          Lrttov_aer,      &
          Lrttov_rad,      &
          Lrttov_cldparam, &
          Lrttov_aerparam
     integer,pointer :: &
          NchanIN,         & ! Number of spectral channels to simulate
          platformIN,      & ! Index of the platform
          satelliteIN,     & ! Index of the satellite
          instrumentIN       ! Index of the instrument
     integer,dimension(RTTOV_MAX_CHANNELS) :: &
         channelsIN          ! Indices of spectral channels
  end type rttov_init_IN
  
CONTAINS

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_rttov_init
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_RTTOV_INIT(NchanIN,platformIN,satelliteIN,instrumentIN,channelsIN,   &
                             nlevels,Lrttov_cld,Lrttov_aer,Lrttov_rad,Lrttov_cldparam, &
                             Lrttov_aerparam,                                          &
                             rttov_input_namelist)
    integer,intent(in) :: & 
         NchanIN,      & ! Number of channels
         platformIN,   & ! Satellite platform
         satelliteIN,  & ! Satellite
         instrumentIN, & ! Instrument
         nlevels
    integer,intent(in),dimension(RTTOV_MAX_CHANNELS) :: &
         channelsIN     ! RTTOV channels
    logical,intent(in)   :: &
         Lrttov_cld,       &
         Lrttov_aer,       &
         Lrttov_rad,       &
         Lrttov_cldparam,  &
         Lrttov_aerparam
    
    ! JKS testing using a RTTOV input namelist here 
    ! (default cosp_rttov namelist is set in cosp.F90)
    character(len=256),intent(in) :: rttov_input_namelist
    
    ! Local variables
    character(len=256) :: &
        coef_file,         &
        scat_file,         &
        rttov_coefDir,     &
        rttov_predDir,     &
        rttov_cldaerDir,   &
        OD_coef_file,      &
        aer_coef_file,     &
        cld_coef_file,     &
        OD_coef_filepath,  &
        aer_coef_filepath, &
        cld_coef_filepath

    ! Declare RTTOV namelist fields
    logical :: so2_data,n2o_data,co_data,ch4_data,co2_data,ozone_data
    character(len=256) :: cosp_status

    ! Read RTTOV namelist fields
    namelist/RTTOV_INPUT/rttov_coefDir,    &
                         OD_coef_filepath,  &
                         aer_coef_filepath,cld_coef_filepath,so2_data,n2o_data,      &
                         co_data,ch4_data,co2_data,ozone_data
!        rttov_Nlocaltime, rttov_localtime, rttov_localtimewindow !JKS

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Read in namelists
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    open(10,file=rttov_input_namelist,status='unknown')
    read(10,nml=RTTOV_INPUT)
    close(10)
    
    !print*,'OD_coef_file:    ',OD_coef_file
    !print*,'aer_coef_file:    ',aer_coef_file
    !print*,'cld_coef_file:    ',cld_coef_file
    !print*,'OD_coef_filepath:    ',OD_coef_filepath

    ! Initialize fields in module memory (cosp_rttovXX.F90)
    nChannels  = NchanIN
    platform   = platformIN 
    satellite  = satelliteIN 
    sensor     = instrumentIN 
    iChannel   = channelsIN
    
    ! Set logicals for RTTOV options
    do_rttov_cld       = Lrttov_cld
    do_rttov_aer       = Lrttov_aer
    do_rttov_rad       = Lrttov_rad       ! to be used in output
    rttov_cld_optparam = Lrttov_cldparam
    rttov_aer_optparam = Lrttov_aerparam

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
    ! Gas profile logicals
    opts%rt_all%ozone_data        = ozone_data
    opts%rt_all%so2_data          = so2_data    
    opts%rt_all%co2_data          = co2_data
    opts%rt_all%n2o_data          = n2o_data
    opts%rt_all%co_data           = co_data
    opts%rt_all%ch4_data          = ch4_data
    
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
    
    ! User options - JKS
    opts%rt_ir%addaerosl               = do_rttov_aer
    opts%rt_ir%addclouds               = do_rttov_cld
    opts%rt_ir%user_aer_opt_param      = rttov_aer_optparam ! User specifies the aerosol scattering optical parameters 
    opts%rt_ir%user_cld_opt_param      = rttov_cld_optparam ! User specifies the cloud scattering optical parameters 
    
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
    
    ! rttovDir should be "/glade/u/home/jonahshaw/w/RTTOV/" and is defined in cosp_config.F90    
    OD_coef_filepath = trim(rttovDir)//trim(rttov_coefDir)//trim(OD_coef_filepath)
    aer_coef_filepath = trim(rttovDir)//trim(rttov_coefDir)//trim(aer_coef_filepath)
    cld_coef_filepath = trim(rttovDir)//trim(rttov_coefDir)//trim(cld_coef_filepath)
         
    print*,'OD_coef_filepath:    ',OD_coef_filepath
    print*,'aer_coef_filepath:   ',aer_coef_filepath
    print*,'cld_coef_filepath:   ',cld_coef_filepath
         
    ! Read optical depth and cloud coefficient files together
    call rttov_read_coefs(errorstatus, coef_rttov, opts,    &
                          file_coef=OD_coef_filepath,       &
                          file_scaer=aer_coef_filepath,     &
                          file_sccld=cld_coef_filepath)
                          
    ! We aren't checking an allocation steps so this seems more appropriate.
    call rttov_error('fatal error reading coefficients' , lalloc = .false.)

    
    ! Ensure input number of channels is not higher than number stored in coefficient file
    if (nchannels > coef_rttov % coef % fmv_chn) then
        nchannels = coef_rttov % coef % fmv_chn
    endif

    ! Ensure the options and coefficients are consistent
    call rttov_user_options_checkinput(errorstatus, opts, coef_rttov)
    
    ! We aren't checking an allocation steps so this seems more appropriate.
    call rttov_error('error in rttov options' , lalloc = .false.) 
    
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
  ! SUBROUTINE cosp_rttov_simulate - Call subroutines in mod_cosp_rttov to run RTTOV
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_RTTOV_SIMULATE(rttovIN,lCleanup,          & ! Inputs
                                 Tb,error)                    ! Outputs
  
    use mod_cosp_rttov,             only:        &
        cosp_rttov_allocate,                     &
        cosp_rttov_construct_profiles,           &
        cosp_rttov_setup_emissivity_reflectance, &
        cosp_rttov_call_direct,                  &
        cosp_rttov_save_and_deallocate_profiles, &
        cosp_rttov_deallocate_coefs
  
    type(rttov_in),intent(in) :: &
        rttovIN
    logical,intent(in) :: &
         lCleanup   ! Flag to determine whether to deallocate RTTOV types
    real(wp),intent(inout),dimension(rttovIN%nPoints,rttovIN%nChannels) :: & ! Can I do this? I guess so!
        Tb        ! RTTOV brightness temperature.
    character(len=128) :: &
        error     ! Error messages (only populated if error encountered)  
    logical                             :: &
        rttov_simulate_cld,                &
        rttov_simulate_aer
    integer(kind=jpim) :: nthreads ! Parallelization, should become an input
    real(wp),dimension(10) :: driver_time

    ! Run each step for running RTTOV from mod_cosp_rttov (and time them)
    call cpu_time(driver_time(1))
    call cosp_rttov_allocate(rttovIN)
    call cpu_time(driver_time(2))
    call cosp_rttov_construct_profiles(rttovIN)
    call cpu_time(driver_time(3))
    call cosp_rttov_setup_emissivity_reflectance()
    call cpu_time(driver_time(4))
    call cosp_rttov_call_direct(nthreads)
    call cpu_time(driver_time(5))
    call cosp_rttov_save_and_deallocate_profiles(rttovIN,Tb)
    call cpu_time(driver_time(6))
    
    print*,'Time to run "cosp_rttov_allocate":     ',                    driver_time(2)-driver_time(1)
    print*,'Time to run "cosp_rttov_construct_profiles":     ',          driver_time(3)-driver_time(2)
    print*,'Time to run "cosp_rttov_setup_emissivity_reflectance":     ',driver_time(4)-driver_time(3)
    print*,'Time to run "cosp_rttov_call_direct":     ',                 driver_time(5)-driver_time(4)
    print*,'Time to run "cosp_rttov_save_and_deallocate_profiles":     ',driver_time(6)-driver_time(5)
    
    ! Deallocate the coefficient files if directed
    if (lCleanup) then
        call cpu_time(driver_time(7))
        call cosp_rttov_deallocate_coefs()
        call cpu_time(driver_time(8))
        print*,'Time to run "cosp_rttov_deallocate_coefs":     ',driver_time(8)-driver_time(7)
    endif

    print*,'Total RTTOV run time:     ',driver_time(8)-driver_time(1)

  END SUBROUTINE COSP_RTTOV_SIMULATE
  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! END MODULE
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE MOD_COSP_RTTOV_INTERFACE
