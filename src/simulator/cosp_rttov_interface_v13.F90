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

  use mod_cosp_rttov,   only: rttov_in                              
                              
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
#include "rttov_get_pc_predictindex.interface"
  
  ! RTTOV variables/structures
  !====================
  LOGICAL(KIND=jplm),      POINTER :: calcemis(:)    => NULL() ! Flag to indicate calculation of emissivity within RTTOV
  LOGICAL(KIND=jplm),      POINTER :: calcrefl(:)    => NULL() ! Flag to indicate calculation of BRDF within RTTOV
  INTEGER(KIND=jpim)               :: errorstatus              ! Return error status of RTTOV subroutine calls

  INTEGER(KIND=jpim) :: alloc_status(60)

  ! DDT for each instrument being simulated. Values to be assigned during the cosp_rttov_init subroutine
  type rttov_cfg
      logical(KIND=jplm)           :: &
          Lrttov_bt,           &
          Lrttov_rad,          &
          Lrttov_refl,         &
          Lrttov_cld,          &
          Lrttov_aer,          &
          Lrttov_pc
      character(len=256)           :: &
          rttov_srcDir,        &
          rttov_coefDir,       &
          OD_coef_filepath,    &
          aer_coef_filepath,   &
          cld_coef_filepath,   &
          PC_coef_filepath
      integer(KIND=jpim)           :: &
          nchanprof,           &
          rttov_direct_nthreads
      integer(KIND=jpim)           :: &
          nchan_out,           &
          nchannels_rec         
      real(wp)                     :: &
          CO2_mr,              &
          CH4_mr,              &
          CO_mr,               &
          N2O_mr,              &
          SO2_mr,              &
          ZenAng
      integer(kind=jpim), allocatable  :: &
          iChannel(:),      &  ! Requested channel indices
          iChannel_out(:)      ! Passing out the channel indices (actual output channels)
      real(kind=jprb),allocatable    :: &
          emisChannel(:),   &                ! RTTOV channel emissivity
          reflChannel(:)                     ! RTTOV channel reflectivity
      type(rttov_options)          :: &
          opts                               ! RTTOV options structure
      type(rttov_coefs)            :: &
          coefs                              ! RTTOV coefficients structure
  end type rttov_cfg
  
  type rttov_output
      integer             :: &
          nchan_out
      integer,pointer     :: &
          channel_indices(:)
      real(wp),pointer    :: &
          bt_total(:,:),    &
          bt_clear(:,:),    &
          rad_total(:,:),   &
          rad_clear(:,:),   &
          rad_cloudy(:,:),  &
          refl_total(:,:),  &
          refl_clear(:,:),  &
          bt_total_pc(:,:), &
          rad_total_pc(:,:)
  end type rttov_output    

  
CONTAINS

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_rttov_init
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_RTTOV_INIT(Lrttov, Nlevels,Ninstruments,instrument_namelists,       &
                             rttov_configs)

      logical,intent(in) :: &
          Lrttov
      integer,intent(in) :: &
          Nlevels,   &
          Ninstruments
      type(character(len=128)), dimension(Ninstruments)     :: & 
          instrument_namelists   ! Array of paths to RTTOV instrument namelists      
      type(rttov_cfg), dimension(:), allocatable :: & ! intent(out)?
          rttov_configs
      ! Local variables
      integer            :: &
          inst_idx ! iterator
        
      allocate(rttov_configs(Ninstruments))
        
      ! Create config objects for each instrument to be simulated by RTTOV. Return to the main subroutine.
      do inst_idx=1,Ninstruments
          print*,'inst_idx:    ',inst_idx
          print*,'instrument_namelists(inst_idx):    ',instrument_namelists(inst_idx)
!          print*,'kind(rttov_configs(inst_idx)):     ',kind(rttov_configs(inst_idx))
!          print*,'rttov_configs(inst_idx):           ',rttov_configs(inst_idx)
          call cosp_rttov_init_s(Nlevels,instrument_namelists(inst_idx),rttov_configs(inst_idx))
      end do
         
       
  END SUBROUTINE COSP_RTTOV_INIT

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_rttov_init
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_RTTOV_INIT_S(Nlevels,namelist_filepath,       &
                               rttov_config)
  
     integer,intent(in)                :: &
         Nlevels
     character(len=128),intent(in)     :: & 
         namelist_filepath   ! Array of paths to RTTOV instrument namelists      
     type(rttov_cfg),intent(out)       :: & ! intent(out)?
         rttov_config 
             
    ! Local variables
    character(len=256),target :: &
        channel_filepath,  &
        rttov_srcDir,      &
        rttov_coefDir,     &
        OD_coef_filepath,  &
        aer_coef_filepath, &
        cld_coef_filepath, &
        PC_coef_filepath
        
    real(wp), target  :: &
        CO2_mr,      &
        CH4_mr,      &
        CO_mr,       &
        N2O_mr,      &
        SO2_mr,      &
        rttov_ZenAng

    ! Declare RTTOV namelist fields
    logical         :: Lrttov_bt = .false.
    logical         :: Lrttov_rad = .false.
    logical         :: Lrttov_refl = .false.
    logical         :: Lrttov_cld = .false.
    logical         :: Lrttov_aer = .false.
    logical         :: Lrttov_cldparam = .false.
    logical         :: Lrttov_aerparam = .false.
    logical         :: Lrttov_pc = .false.

    logical         :: Lchannel_filepath = .false.
    
    logical         :: SO2_data = .false. 
    logical         :: N2O_data = .false. 
    logical         :: CO_data = .false.
    logical         :: CO2_data = .false.
    logical         :: CH4_data = .false. 
    logical         :: ozone_data = .false.
    logical         :: clw_data = .false.
    
    character(len=256) :: cosp_status
    integer ::             &
        i,                 &
        nchannels_rec
        
    integer, target ::     &
        rttov_nthreads
                    
    integer(kind=jpim)   :: &
        ipcbnd,        &
        ipcreg,        &
        npcscores
                
    ! Read RTTOV namelist fields
    namelist/RTTOV_INPUT/Lrttov_bt,Lrttov_rad,Lrttov_refl,Lrttov_cld,            & ! Logicals for RTTOV configuration
                         Lrttov_aer,Lrttov_cldparam,Lrttov_aerparam,             & ! 
                         Lrttov_pc,nchannels_rec,Lchannel_filepath,            &
                         channel_filepath,rttov_srcDir,rttov_coefDir,            &
                         OD_coef_filepath,aer_coef_filepath,cld_coef_filepath,   &
                         PC_coef_filepath,                                       &
                         CO2_data,CH4_data,CO_data,N2O_data,SO2_data,ozone_data, & ! User-supplied trace gas concentrations
                         clw_data,                                               & ! MW option
                         CO2_mr,CH4_mr,CO_mr,N2O_mr,SO2_mr,                      & ! Mixing ratios
                         ipcbnd,ipcreg,npcscores,                                & ! PC-RTTOV config values
                         rttov_nthreads,rttov_ZenAng
                                                      
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Read in namelists
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    open(10,file=namelist_filepath,status='unknown')
    read(10,nml=RTTOV_INPUT)
    close(10)
    
    print*,'Lrttov_pc:   ',Lrttov_pc

    ! Set logicals for RTTOV config
    rttov_config%Lrttov_bt         = Lrttov_bt
    rttov_config%Lrttov_rad        = Lrttov_rad
    rttov_config%Lrttov_refl       = Lrttov_refl
    rttov_config%Lrttov_pc         = Lrttov_pc
        
    ! Set paths for RTTOV config
    rttov_config%rttov_srcDir      = rttov_srcDir
    rttov_config%rttov_coefDir     = rttov_coefDir

    ! Construct optical depth and cloud coefficient files    
    rttov_config%OD_coef_filepath  = trim(rttov_config%rttov_srcDir)//trim(rttov_config%rttov_coefDir)//trim(OD_coef_filepath)
    rttov_config%aer_coef_filepath = trim(rttov_config%rttov_srcDir)//trim(rttov_config%rttov_coefDir)//trim(aer_coef_filepath)
    rttov_config%cld_coef_filepath = trim(rttov_config%rttov_srcDir)//trim(rttov_config%rttov_coefDir)//trim(cld_coef_filepath)
    rttov_config%PC_coef_filepath  = trim(rttov_config%rttov_srcDir)//trim(rttov_config%rttov_coefDir)//trim(PC_coef_filepath)
        
    ! Set other RTTOV config variables
    rttov_config%rttov_direct_nthreads = rttov_nthreads
    
    rttov_config%SO2_mr       = SO2_mr
    rttov_config%N2O_mr       = N2O_mr
    rttov_config%CO_mr        = CO_mr
    rttov_config%CO2_mr       = CO2_mr
    rttov_config%CH4_mr       = CH4_mr
    rttov_config%ZenAng       = rttov_ZenAng

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 1. Initialise RTTOV options structure
    ! ------------------------------------------------------
    ! See page 157 of RTTOV v13 user guide for documentation
    ! Initializing all options to defaults for consistency
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ! General configuration options
    rttov_config % opts % config % do_checkinput    = .true.
    rttov_config % opts % config % apply_reg_limits = .false. ! True in v11
    rttov_config % opts % config % verbose          = .false.  ! JKS suppress for now
    rttov_config % opts % config % opdep13_gas_clip = .true.

    ! Declare RTTOV namelist fields
    rttov_config % opts % rt_all % SO2_data   = SO2_data
    rttov_config % opts % rt_all % N2O_data   = N2O_data
    rttov_config % opts % rt_all % CO_data    = CO_data
    rttov_config % opts % rt_all % CO2_data   = CO2_data
    rttov_config % opts % rt_all % CH4_data   = CH4_data
    rttov_config % opts % rt_all % ozone_data = ozone_data
    rttov_config % opts % rt_mw  % clw_data   = clw_data
    
    ! Other general RT options (initializing to defaults for completeness)
    rttov_config % opts % rt_all % do_lambertian       = .false.
    rttov_config % opts % rt_all % switchrad           = .false.
    rttov_config % opts % rt_all % rad_down_lin_tau    = .true.
    rttov_config % opts % rt_all % use_t2m_opdep       = .true.
    rttov_config % opts % rt_all % use_q2m             = .true.
    rttov_config % opts % rt_all % use_tskin_eff       = .false.
    rttov_config % opts % rt_all % addrefrac           = .true.
    rttov_config % opts % rt_all % plane_parallel      = .false.
    rttov_config % opts % rt_all % transmittances_only = .false.
    
    ! MW-only radiative transfer options:
    ! JKS make this optional?
    rttov_config % opts % rt_mw % clw_data             = .false.
    rttov_config % opts % rt_mw % clw_scheme           = 2       ! Default = 2/Rosenkranz
    rttov_config % opts % rt_mw % clw_cloud_top        = 322     ! Default is 322 hPa
    rttov_config % opts % rt_mw % fastem_version       = 6       ! Default FASTEM-6
    rttov_config % opts % rt_mw % supply_foam_fraction = .false.

    ! UV/visible/IR-only radiative transfer options
    rttov_config % opts % rt_ir % addsolar                = .false.
    rttov_config % opts % rt_ir % rayleigh_max_wavelength = 2._wp ! 2um 
    rttov_config % opts % rt_ir % rayleigh_min_pressure   = 0._wp ! 0 hPa
    rttov_config % opts % rt_ir % rayleigh_single_scatt   = .true.
    rttov_config % opts % rt_ir % rayleigh_depol          = .true. ! Default false, recommended true
    rttov_config % opts % rt_ir % do_nlte_correction      = .false.
    rttov_config % opts % rt_ir % solar_sea_brdf_model    = 2
    rttov_config % opts % rt_ir % ir_sea_emis_model       = 2

    ! User options - JKS
    ! Duplicate for STUB functionality
    rttov_config % Lrttov_aer                             = Lrttov_aer
    rttov_config % Lrttov_cld                             = Lrttov_cld
    rttov_config % opts % rt_ir % addaerosl               = Lrttov_aer
    rttov_config % opts % rt_ir % addclouds               = Lrttov_cld
    rttov_config % opts % rt_ir % user_aer_opt_param      = Lrttov_aerparam ! User specifies the aerosol scattering optical parameters 
    rttov_config % opts % rt_ir % user_cld_opt_param      = Lrttov_cldparam ! User specifies the cloud scattering optical parameters 
    
    rttov_config % opts % rt_ir % grid_box_avg_cloud      = .true.
    rttov_config % opts % rt_ir % cldcol_threshold        = -1._wp
    rttov_config % opts % rt_ir % cloud_overlap           = 1 ! Maximum-random overlap
    rttov_config % opts % rt_ir % cc_low_cloud_top        = 750_wp ! 750 hPa. Only applies when cloud_overlap=2.
    rttov_config % opts % rt_ir % ir_scatt_model          = 2
    rttov_config % opts % rt_ir % vis_scatt_model         = 1
    rttov_config % opts % rt_ir % dom_nstreams            = 8
    rttov_config % opts % rt_ir % dom_accuracy            = 0._wp ! only applies when addclouds or addaerosl is true and DOM is selected as a scattering solver.
    rttov_config % opts % rt_ir % dom_opdep_threshold     = 0._wp
    rttov_config % opts % rt_ir % dom_rayleigh            = .false.
    
    ! Principal Components-only radiative transfer options:
    ! Default off
    rttov_config % opts % rt_ir % pc % addpc     = .false.
    rttov_config % opts % rt_ir % pc % npcscores = -1
    rttov_config % opts % rt_ir % pc % addradrec = .false.
    rttov_config % opts % rt_ir % pc % ipcbnd    = 1
    rttov_config % opts % rt_ir % pc % ipcreg    = 1 ! The index of the required set of PC predictors
    
    ! Options related to interpolation and the vertical grid:
    rttov_config % opts % interpolation % addinterp         = .true.
    rttov_config % opts % interpolation % interp_mode       = 1
!    rttov_config % opts % interpolation % reg_limit_extrap = .true. ! Depreciated
    rttov_config % opts % interpolation % lgradp            = .false.
!    rttov_config % opts % interpolation % spacetop         = .true. ! Depreciated

    ! Options related to HTFRTC:
    rttov_config % opts % htfrtc_opts % htfrtc       = .false.
    rttov_config % opts % htfrtc_opts % n_pc_in      = -1
    rttov_config % opts % htfrtc_opts % reconstruct  = .false.
    rttov_config % opts % htfrtc_opts % simple_cloud = .false.
    rttov_config % opts % htfrtc_opts % overcast     = .false.
    
    ! Developer options that may be useful:
    rttov_config % opts % dev % do_opdep_calc = .true.
    
    ! If using PC-RTTOV, some settings must be a certain way. This isn't always true though...
    if (Lrttov_pc) then
      rttov_config % opts % rt_ir % pc % addpc          = .true.
      rttov_config % opts % rt_ir % pc % ipcbnd         = ipcbnd
      rttov_config % opts % rt_ir % pc % ipcreg         = ipcreg
      rttov_config % opts % rt_ir % pc % npcscores      = npcscores

      rttov_config % opts % rt_ir % pc % addradrec      = .true. ! Alway reconstruct radiances

      rttov_config % opts % interpolation % addinterp   = .true.  ! Allow interpolation of input profile
      rttov_config % opts % interpolation % interp_mode = 1       ! Set interpolation method
      rttov_config % opts % rt_all % addrefrac          = .true.  ! Include refraction in path calc (always for PC)
      rttov_config % opts % rt_ir % addclouds           = .false. ! Don't include cloud effects     (always for PC?)
      rttov_config % opts % rt_ir % addaerosl           = .false. ! Don't include aerosol effects   (not always for PC)
      rttov_config % opts % rt_ir % addsolar            = .false. ! Do not include solar radiation  (always for PC?)

    endif
    
    ! JKS To-do: include opts_scatt settings (user guide pg 161)

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 2. Read coefficients (from RTTOV example files)
    ! ------------------------------------------------------
    ! Using the GUI to figure out files that work together could be helpful here.
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    ! Do I need logicals here to direct how to read the coefficients?
    if (rttov_config % Lrttov_pc) then ! PC-RTTOV cannot handle cloud, some aerosols
        call rttov_read_coefs(errorstatus, rttov_config % coefs, &
                              rttov_config % opts,               &
                              file_coef=rttov_config % OD_coef_filepath,        &
!                              file_scaer=rttov_config % aer_coef_filepath,     & ! Needs to be PC-RTTOV compatible
                              file_pccoef=rttov_config % PC_coef_filepath) 
    else ! Read optical depth and cloud coefficient files together
        call rttov_read_coefs(errorstatus, rttov_config % coefs, & 
                              rttov_config % opts,               &
                              file_coef=rttov_config % OD_coef_filepath,        &
                              file_scaer=rttov_config % aer_coef_filepath,      &
                              file_sccld=rttov_config % cld_coef_filepath)

        ! Ensure input number of channels is not higher than number stored in coefficient file
        if (nchannels_rec > rttov_config % coefs % coef % fmv_chn) then
            nchannels_rec = rttov_config % coefs % coef % fmv_chn
            print*,'nchannels_rec cap hit'
        endif            
    endif
                          
    ! We aren't checking an allocation steps so this seems more appropriate.
    call rttov_error('fatal error reading coefficients' , lalloc = .false.)

    ! Ensure the options and coefficients are consistent
    call rttov_user_options_checkinput(errorstatus, rttov_config % opts, rttov_config % coefs)
    
    ! We aren't checking an allocation steps so this seems more appropriate.
    call rttov_error('error in rttov options' , lalloc = .false.)


    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Figure out how many channels we actually want to reconstruct
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Handle different radiance reconstruction options
    ! JKS - I need to set these two separate variables here? Can one be set later?
    if (nchannels_rec < 0) then
        print*,'The namelist variable "nchannels_rec" is negative, rttov_direct call will fail. Exiting.'
        errorstatus = errorstatus_fatal
        call rttov_exit(errorstatus)
        ! If the number of channels is negative, don't reconstruct radiances at all
        rttov_config % nchan_out = 0
        rttov_config % nchannels_rec = 0 ! Avoid nchanprof set to a negative value
    else if (nchannels_rec == 0) then
        ! If the number of channels is set to 0 then reconstruct all instrument channels
        rttov_config % nchan_out     = rttov_config % coefs % coef % fmv_chn
        rttov_config % nchannels_rec = rttov_config % coefs % coef % fmv_chn ! Avoid nchanprof set to 0
    else
        ! Otherwise read the channel list from the file
        rttov_config % nchan_out     = nchannels_rec
        rttov_config % nchannels_rec = nchannels_rec
    endif


    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Read in channel indices, emissivities, and reflectivities from .csv if file is passed
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    if (Lchannel_filepath) then
        allocate(rttov_config % iChannel(rttov_config % nchan_out))
        allocate(rttov_config % emisChannel(rttov_config % nchan_out))
        allocate(rttov_config % reflChannel(rttov_config % nchan_out))
    
        open(18,file=channel_filepath,access='sequential',form="formatted")
        do i = 1, rttov_config % nchan_out
            read(18,*) rttov_config % iChannel(i), rttov_config % emisChannel(i), rttov_config % reflChannel(i)
        end do
        close(18)
    else ! If nothing is passed, compute the first "nchan_out" channels. Ignore emissivity and reflectivity for now.
        allocate(rttov_config % iChannel(rttov_config % nchan_out))
        rttov_config % iChannel(:) = (/ (i, i = 1, rttov_config % nchan_out) /)
    endif
        
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
 
  END SUBROUTINE COSP_RTTOV_INIT_S

  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_rttov_simulate - Call subroutines in mod_cosp_rttov to run RTTOV
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_RTTOV_SIMULATE(rttovIN,rttovConfig,lCleanup,                  & ! Inputs
                                 bt_total,bt_clear,                                & ! Brightness Temp Outputs
                                 rad_total,rad_clear,rad_cloudy,                   & ! Radiance Outputs
                                 refl_total,refl_clear,                            & ! Reflectance Outputs
                                 error)      

    type(rttov_in),intent(in) :: &
        rttovIN
    type(rttov_cfg),intent(inout) :: &
        rttovConfig
    logical,intent(in) :: &
        lCleanup   ! Flag to determine whether to deallocate RTTOV types
    real(wp),intent(inout),dimension(rttovIN%nPoints,rttovConfig%nchan_out) :: & ! Can I do this? I guess so! 
        bt_total,                          &        ! All-sky
        bt_clear,                          &        ! Clear-sky
        rad_total,                         &        ! All-sky
        rad_clear,                         &        ! Clear-sky
        rad_cloudy,                        &        ! Cloudy-sky
        refl_total,                        &        ! All-sky
        refl_clear                                  ! Clear-sky
    character(len=128) :: &
        error     ! Error messages (only populated if error encountered) 

    ! Check options to determine if the principal component approach should be run
    if (rttovConfig % opts % rt_ir % pc % addpc) then
        print*,'Running COSP_PC_RTTOV_SIMULATE in multi-inst. set-up'
        call COSP_PC_RTTOV_SIMULATE(rttovIN,rttovConfig,lCleanup,                     &
                                    bt_clear,rad_clear,                               &
                                    error)                                

        print*,'end of COSP_PC_RTTOV_SIMULATE'
    else
        print*,'Running COSP_REG_RTTOV_SIMULATE in multi-inst. set-up'
        call COSP_REG_RTTOV_SIMULATE(rttovIN,rttovConfig,lCleanup,                                 &
                                     bt_total,bt_clear,                                &
                                     rad_total,rad_clear,rad_cloudy,                   &
                                     refl_total,refl_clear,                            &
                                     error)
        
    endif

  END SUBROUTINE COSP_RTTOV_SIMULATE

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_rttov_simulate - Call regular subroutines in mod_cosp_rttov to run RTTOV
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_REG_RTTOV_SIMULATE(rttovIN,rttovConfig,lCleanup,                  & ! Inputs
                                     bt_total,bt_clear,                                & ! Brightness Temp Outputs
                                     rad_total,rad_clear,rad_cloudy,                   & ! Radiance Outputs
                                     refl_total,refl_clear,                            & ! Reflectance Outputs
                                     error)                                              
  
    use mod_cosp_rttov,             only:        &
        cosp_rttov_allocate,                     &
        cosp_rttov_construct_profiles,           &
        cosp_rttov_setup_emissivity_reflectance, &
        cosp_rttov_call_direct,                  &
        cosp_rttov_save_output,                  &
        cosp_rttov_deallocate_profiles,          &
        cosp_rttov_deallocate_coefs
  
    type(rttov_in),intent(in) :: &
        rttovIN
    type(rttov_cfg),intent(inout) :: &
        rttovConfig
    logical,intent(in) :: &
         lCleanup   ! Flag to determine whether to deallocate RTTOV types
    real(wp),intent(inout),dimension(rttovIN%nPoints,rttovConfig%nchan_out) :: & 
        bt_total,                          &        ! All-sky
        bt_clear,                          &        ! Clear-sky
        rad_total,                         &        ! All-sky
        rad_clear,                         &        ! Clear-sky
        rad_cloudy,                        &        ! Cloudy-sky
        refl_total,                        &        ! All-sky
        refl_clear                                  ! Clear-sky
    character(len=128) :: &
        error     ! Error messages (only populated if error encountered)  

    real(wp),dimension(10) :: driver_time

    ! Run each step for running RTTOV from mod_cosp_rttov (and time them)
    print*,'cosp_rttov_allocate begin' ! jks
    call cpu_time(driver_time(1))
    call cosp_rttov_allocate(rttovIN,                       &
                             rttovConfig % nChannels_rec,   &
                             rttovConfig % opts,            &
                             rttovConfig % coefs,           &
                             rttovConfig % iChannel,        &
                             rttovConfig % nchanprof)    
        
    print*,'cosp_rttov_allocate successful' ! jks
    call cpu_time(driver_time(2))
    call cosp_rttov_construct_profiles(rttovIN, &
                                       rttovConfig % Lrttov_cld,               &
                                       rttovConfig % Lrttov_aer,               &
                                       rttovConfig % CO2_mr,                   &
                                       rttovConfig % CH4_mr,                   &
                                       rttovConfig % CO_mr,                    &
                                       rttovConfig % N2O_mr,                   &
                                       rttovConfig % SO2_mr,                   &
                                       rttovConfig % ZenAng)
          
    print*,'cosp_rttov_construct_profiles successful' ! jks
    call cpu_time(driver_time(3))
    call cosp_rttov_setup_emissivity_reflectance() ! Config agnostic after allocate step.
    print*,'cosp_rttov_setup_emissivity_reflectance successful' ! jks
    call cpu_time(driver_time(4))
    call cosp_rttov_call_direct(rttovConfig % rttov_direct_nthreads,  &
                                rttovConfig % opts,                   &
                                rttovConfig % coefs)    
    
    print*,'cosp_rttov_call_direct successful' ! jks
    call cpu_time(driver_time(5))
    
    call cosp_rttov_save_output(rttovIN,                        &
                                rttovConfig % nchan_out,                &
                                rttovConfig % Lrttov_bt,                &
                                rttovConfig % Lrttov_rad,               &
                                rttovConfig % Lrttov_refl,              &
                                rttovConfig % Lrttov_cld,               &
                                rttovConfig % Lrttov_aer,               &                                
                                bt_total,bt_clear,                      &
                                rad_total,rad_clear,rad_cloudy,         &
                                refl_total,refl_clear)

    print*,'cosp_rttov_save_output successful' ! jks
    call cpu_time(driver_time(6))
    call cosp_rttov_deallocate_profiles(rttovIN,                       &
                                        rttovConfig % opts,            &
                                        rttovConfig % coefs,           &
                                        rttovConfig % nchanprof)    
    call cpu_time(driver_time(7))
    
    print*,'Time to run "cosp_rttov_allocate":     ',                    driver_time(2)-driver_time(1)
    print*,'Time to run "cosp_rttov_construct_profiles":     ',          driver_time(3)-driver_time(2)
    print*,'Time to run "cosp_rttov_setup_emissivity_reflectance":     ',driver_time(4)-driver_time(3)
    print*,'Time to run "cosp_rttov_call_direct":     ',                 driver_time(5)-driver_time(4)
    print*,'Time to run "cosp_rttov_save_output":     ',                 driver_time(6)-driver_time(5)
    print*,'Time to run "cosp_rttov_deallocate_profiles":     ',         driver_time(7)-driver_time(6)
    
    ! Deallocate the coefficient files if directed
    if (lCleanup) then
        call cpu_time(driver_time(8))
        call cosp_rttov_deallocate_coefs(rttovConfig % coefs)
        call cpu_time(driver_time(9))
        print*,'Time to run "cosp_rttov_deallocate_coefs":     ',driver_time(9)-driver_time(8)
    endif

!    print*,'Total RTTOV run time:     ',driver_time(8)-driver_time(1)

  END SUBROUTINE COSP_REG_RTTOV_SIMULATE


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_pc_rttov_simulate - Call subroutines in mod_cosp_rttov to run RTTOV
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_PC_RTTOV_SIMULATE(rttovIN,rttovConfig,lCleanup,                  & ! Inputs
                                    bt_clear,rad_clear,                               & ! Outputs
                                    error)                                              
  
    use mod_cosp_rttov,             only:   &
        cosp_pc_rttov_allocate,             &
        cosp_rttov_construct_profiles,      &
        cosp_pc_rttov_setup_emissivity,     &
        cosp_pc_rttov_call_direct,          &
        cosp_pc_rttov_save_output,          &
        cosp_pc_rttov_deallocate_profiles,  &
        cosp_rttov_deallocate_coefs
  
    type(rttov_in),intent(in) :: &
        rttovIN
    type(rttov_cfg),intent(inout) :: &
        rttovConfig        
    logical,intent(in) :: &
         lCleanup   ! Flag to determine whether to deallocate RTTOV types
    real(wp),intent(inout),dimension(rttovIN%nPoints,rttovConfig%nchan_out) :: & ! Can I do this? I guess so!
        bt_clear,                          &        ! All-sky
        rad_clear                                   ! All-sky
    character(len=128) :: &
        error     ! Error messages (only populated if error encountered)  

    real(wp),dimension(10) :: driver_time

    ! Run each step for running RTTOV from mod_cosp_rttov (and time them)
!    print*,'cosp_rttov_allocate begin' ! jks
    call cpu_time(driver_time(1))
    call cosp_pc_rttov_allocate(rttovIN, &
                                rttovConfig % PC_coef_filepath,              &
                                rttovConfig % coefs,                         &
                                rttovConfig % opts,                          &
                                rttovConfig % nchannels_rec,                 &
                                rttovConfig % iChannel,                      &
                                rttovConfig % nchanprof,                     &
                                rttovConfig % iChannel_out)
    print*,'cosp_pc_rttov_allocate successful' ! jks
    call cpu_time(driver_time(2))
    call cosp_rttov_construct_profiles(rttovIN, &
                                       rttovConfig % Lrttov_cld,               &
                                       rttovConfig % Lrttov_aer,               &
                                       rttovConfig % CO2_mr,                   &
                                       rttovConfig % CH4_mr,                   &
                                       rttovConfig % CO_mr,                    &
                                       rttovConfig % N2O_mr,                   &
                                       rttovConfig % SO2_mr,                   &
                                       rttovConfig % ZenAng)
    print*,'cosp_rttov_construct_profiles successful' ! jks
    call cpu_time(driver_time(3))
    call cosp_pc_rttov_setup_emissivity()
    print*,'cosp_pc_rttov_setup_emissivity successful' ! jks
    call cpu_time(driver_time(4))
    call cosp_pc_rttov_call_direct(rttovConfig % rttov_direct_nthreads,  &
                                   rttovConfig % opts,                   &
                                   rttovConfig % coefs,                  &
                                   rttovConfig % nchannels_rec,          &
                                   rttovConfig % iChannel_out) ! iChannel_out should have been updated

    print*,'cosp_pc_rttov_call_direct successful' ! jks
    call cpu_time(driver_time(5))
    call cosp_pc_rttov_save_output(rttovIN%nPoints,                       &
                                   rttovConfig%nchannels_rec,             &
                                   rttovConfig%Lrttov_bt,                 &
                                   rttovConfig%Lrttov_rad,                &
                                   bt_clear,                              &
                                   rad_clear)

    print*,'cosp_pc_rttov_save_output successful' ! jks
    call cpu_time(driver_time(6))
    call cosp_pc_rttov_deallocate_profiles(rttovIN,                       &
                                           rttovConfig % nChannels_rec,   &
                                           rttovConfig % opts,            &
                                           rttovConfig % coefs,           &
                                           rttovConfig % nchanprof)    
    call cpu_time(driver_time(7))
    
    print*,'Time to run "cosp_pc_rttov_allocate":     ',                    driver_time(2)-driver_time(1)
    print*,'Time to run "cosp_rttov_construct_profiles":     ',             driver_time(3)-driver_time(2)
    print*,'Time to run "cosp_pc_rttov_setup_emissivity":     ',            driver_time(4)-driver_time(3)
    print*,'Time to run "cosp_pc_rttov_call_direct":     ',                 driver_time(5)-driver_time(4)
    print*,'Time to run "cosp_pc_rttov_save_output":     ',                 driver_time(6)-driver_time(5)
    print*,'Time to run "cosp_pc_rttov_deallocate_profiles":     ',         driver_time(7)-driver_time(6)
    
    ! Deallocate the coefficient files if directed
    if (lCleanup) then
        call cpu_time(driver_time(8))
!        call cosp_rttov_deallocate_coefs(rttovConfig % coefs) ! JKS this again with PC-RTTOV...
        call cpu_time(driver_time(9))
        print*,'Time to run "cosp_rttov_deallocate_coefs":     ',driver_time(9)-driver_time(8)
    endif

  END SUBROUTINE COSP_PC_RTTOV_SIMULATE

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! END MODULE
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE MOD_COSP_RTTOV_INTERFACE
