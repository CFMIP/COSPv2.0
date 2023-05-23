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
  use mod_cosp_rttov,   only: nchannels_rec,iChannel,emisChannel,reflChannel,            &
                              coef_rttov,opts,rttov_in,                                  &
                              do_rttov_bt,do_rttov_rad,do_rttov_refl,                    &
                              do_rttov_cld,do_rttov_aer,do_rttov_pcrttov,                &
                              rttov_cld_optparam,rttov_aer_optparam,                     &
                              rttov_direct_nthreads,rttovDir,PC_coef_filepath,           &
                              so2,ch4,co,co2,n2o,zenang,npcscores,predictindex,          &
                              iChannel_out,Lchannel_filepath
                              
                              
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
          Lrttov_pc
      character(len=256)           :: &
          rttov_srcDir,        &
          rttov_coefDir,       &
          OD_coef_filepath,    &
          aer_coef_filepath,   &
          cld_coef_filepath,   &
          PC_coef_filepath
      integer(KIND=jpim)           :: &
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
          rttov_ZenAng
      integer(kind=jpim), allocatable  :: &
          iChannel(:),      &  ! Requested channel indices
          iChannel_out(:)      ! Passing out the channel indices (actual output channels)
      real(kind=jplm),allocatable    :: &
          emisChannel(:),   &                ! RTTOV channel emissivity
          reflChannel(:)                     ! RTTOV channel reflectivity
      type(rttov_options)          :: &
          opts                               ! RTTOV options structure
      type(rttov_coefs)            :: &
          coefs
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
  ! SUBROUTINE cosp_rttov_ini2
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_RTTOV_INI2(Nlevels,Ninstruments,instrument_namelists,       &
                             rttov_configs)

      integer,intent(in) :: &
          Nlevels,   &
          Ninstruments
      type(character(len=256)), dimension(Ninstruments)     :: & 
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
          print*,'Nlevels:     ',Nlevels
          print*,'instrument_namelists(inst_idx):    ',instrument_namelists(inst_idx)
!          print*,'kind(rttov_configs(inst_idx)):     ',kind(rttov_configs(inst_idx))
!          print*,'rttov_configs(inst_idx):           ',rttov_configs(inst_idx)
          call cosp_rttov_init_s(Nlevels,instrument_namelists(inst_idx),rttov_configs(inst_idx))
      end do
         
       
  END SUBROUTINE COSP_RTTOV_INI2

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_rttov_init
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_RTTOV_INIT_S(Nlevels,namelist_filepath,       &
                               rttov_config)
  
     integer,intent(in)                :: &
         Nlevels
     character(len=256),intent(in)     :: & 
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

    logical         :: Lchannel_filepath2 = .false.
    
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
        nchannels_rec2
        
    integer, target ::     &
        rttov_nthreads
                    
    integer(kind=jpim)   :: &
        ipcbnd,        &
        ipcreg,        &
        npcscores
                
    ! Read RTTOV namelist fields
    namelist/RTTOV_INPUT/Lrttov_bt,Lrttov_rad,Lrttov_refl,Lrttov_cld,            & ! Logicals for RTTOV configuration
                         Lrttov_aer,Lrttov_cldparam,Lrttov_aerparam,             & ! 
                         Lrttov_pc,nchannels_rec2,Lchannel_filepath2,            &
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
    rttov_config%rttov_ZenAng = rttov_ZenAng

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
        if (nchannels_rec2 > rttov_config % coefs % coef % fmv_chn) then
            nchannels_rec2 = rttov_config % coefs % coef % fmv_chn
            print*,'nchannels_rec2 cap hit'
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
    if (nchannels_rec2 < 0) then
        print*,'The namelist variable "nchannels_rec2" is negative, rttov_direct call will fail. Exiting.'
        errorstatus = errorstatus_fatal
        call rttov_exit(errorstatus)
        ! If the number of channels is negative, don't reconstruct radiances at all
        rttov_config % nchan_out = 0
        rttov_config % nchannels_rec = 0 ! Avoid nchanprof set to a negative value
    else if (nchannels_rec2 == 0) then
        ! If the number of channels is set to 0 then reconstruct all instrument channels
        rttov_config % nchan_out     = rttov_config % coefs % coef % fmv_chn
        rttov_config % nchannels_rec = rttov_config % coefs % coef % fmv_chn ! Avoid nchanprof set to 0
    else
        ! Otherwise read the channel list from the file
        rttov_config % nchan_out     = nchannels_rec2
        rttov_config % nchannels_rec = nchannels_rec2
    endif


    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Read in channel indices, emissivities, and reflectivities from .csv if file is passed
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    if (Lchannel_filepath2) then
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
  ! SUBROUTINE cosp_rttov_init
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_RTTOV_INIT(nchan_out,nlevels,Lrttov_bt,Lrttov_rad,Lrttov_refl,       &
                             Lrttov_cld,Lrttov_aer,Lrttov_cldparam,Lrttov_aerparam,    &
                             Lrttov_pc,rttov_input_namelist)
    integer,intent(inout) :: &
         nchan_out    ! JKS make this out only soon
    integer,intent(in) :: &
         nlevels
    logical,intent(in)   :: &
         Lrttov_bt,        &
         Lrttov_rad,       &
         Lrttov_refl,      &
         Lrttov_cld,       &
         Lrttov_aer,       &
         Lrttov_cldparam,  &
         Lrttov_aerparam,  &
         Lrttov_pc
    
    ! JKS testing using a RTTOV input namelist here 
    ! (default cosp_rttov namelist is set in cosp.F90)
    character(len=256),intent(in) :: rttov_input_namelist
    
    ! Local variables
    character(len=256) :: &
        rttov_srcDir,      &
        rttov_coefDir,     &
        OD_coef_file,      &
        aer_coef_file,     &
        cld_coef_file,     &
        OD_coef_filepath,  &
        aer_coef_filepath, &
        cld_coef_filepath, &
        channel_filepath
        
    real(wp)  :: &
        SO2_mr,      &
        N2O_mr,      &
        CO_mr,       &
        CH4_mr,      &
        CO2_mr,      &
        rttov_ZenAng

    ! Declare RTTOV namelist fields
    logical :: SO2_data = .false. 
    logical :: N2O_data = .false. 
    logical :: CO_data = .false.
    logical :: CO2_data = .false.
    logical :: CH4_data = .false. 
    logical :: ozone_data = .false. 
        
    character(len=256) :: cosp_status
    integer ::             &
        i,                 &
        rttov_nthreads
                    
    integer(kind=jpim) :: ipcbnd, ipcreg
            
    ! Read RTTOV namelist fields
    namelist/RTTOV_INPUT/Lchannel_filepath,rttov_srcDir,rttov_coefDir,           &
                         OD_coef_filepath,aer_coef_filepath,cld_coef_filepath,   &
                         SO2_mr,N2O_mr,CO_mr,CH4_mr,CO2_mr,rttov_ZenAng,         & ! Mixing ratios
                         SO2_data,N2O_data,CO_data,CH4_data,CO2_data,ozone_data, &
                         rttov_nthreads,nchannels_rec
                             
    ! Only read channel indices and emissivities if prompted
    if (Lchannel_filepath) then
        namelist/RTTOV_INPUT/channel_filepath
    endif
              
    ! Only read some namelist fields if PC-RTTOV will run  
    if (Lrttov_pc) then
        namelist/RTTOV_INPUT/PC_coef_filepath,ipcbnd,ipcreg,npcscores
    endif
    
    !! JKS - Hardcode in some options that will eventually be moved to the namelist
!    ipcbnd = 1 ! This should always be one per the User Guide
!    ipcreg = 2 ! 300 predictors (channels). See RTTOV user guide Table 31.
!    npcscores = 100 ! 100 principal component scores

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Read in namelists
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    open(10,file=rttov_input_namelist,status='unknown')
    read(10,nml=RTTOV_INPUT)
    close(10)
    

    ! Initialize fields in module memory (cosp_rttovXX.F90)
    rttovDir   = rttov_srcDir
    
    ! Set logicals for RTTOV options
    do_rttov_bt        = Lrttov_bt
    do_rttov_rad       = Lrttov_rad
    do_rttov_refl      = Lrttov_refl
    do_rttov_cld       = Lrttov_cld
    do_rttov_aer       = Lrttov_aer
    do_rttov_pcrttov   = Lrttov_pc
    rttov_cld_optparam = Lrttov_cldparam
    rttov_aer_optparam = Lrttov_aerparam
    
    so2 = SO2_mr
    n2o = N2O_mr
    co  = CO_mr
    co2 = CO2_mr
    ch4 = CH4_mr
    zenang = rttov_ZenAng
    
    ! Set degree of parallelization
    rttov_direct_nthreads = rttov_nthreads

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 1. Initialise RTTOV options structure
    ! ------------------------------------------------------
    ! See page 157 of RTTOV v13 user guide for documentation
    ! Initializing all options to defaults for consistency
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ! General configuration options
    opts%config%do_checkinput    = .true.
    opts%config%apply_reg_limits = .false. ! True in v11
    opts%config%verbose          = .false.  ! JKS suppress for now
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
    ! Default off
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
    
    ! If using PC-RTTOV, some settings must be a certain way. This isn't always true though...
    if (do_rttov_pcrttov) then
      opts % rt_ir % pc % addpc          = .true.
      opts % rt_ir % pc % ipcbnd         = ipcbnd
      opts % rt_ir % pc % ipcreg         = ipcreg
      opts % rt_ir % pc % npcscores      = npcscores

      ! In this example we reconstruct radiances if there is an input file
      ! containing a channel list
!      INQUIRE(file=radrec_filename, exist=exists) ! exists is a logical, but I don't think we need it
      opts % rt_ir % pc % addradrec      = .true. ! exists, I could add a logical in the namelist to determine this

      opts % interpolation % addinterp   = .true.  ! Allow interpolation of input profile
      opts % interpolation % interp_mode = 1       ! Set interpolation method
      opts % rt_all % addrefrac          = .true.  ! Include refraction in path calc (always for PC)
      opts % rt_ir % addclouds           = .false. ! Don't include cloud effects     (always for PC?)
      opts % rt_ir % addaerosl           = .false. ! Don't include aerosol effects   (not always for PC)
      opts % rt_ir % addsolar            = .false. ! Do not include solar radiation  (always for PC?)

    endif
    
    ! JKS To-do: include opts_scatt settings (user guide pg 161)

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 2. Read coefficients (from RTTOV example files)
    ! ------------------------------------------------------
    ! Using the GUI to figure out files that work together could be helpful here.
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Construct optical depth and cloud coefficient files
    OD_coef_filepath  = trim(rttovDir)//trim(rttov_coefDir)//trim(OD_coef_filepath)
    aer_coef_filepath = trim(rttovDir)//trim(rttov_coefDir)//trim(aer_coef_filepath)
    cld_coef_filepath = trim(rttovDir)//trim(rttov_coefDir)//trim(cld_coef_filepath)
                  
    ! Do I need logicals here to direct how to read the coefficients?
    
    
    if (do_rttov_pcrttov) then ! PC-RTTOV cannot handle cloud, some aerosols
        PC_coef_filepath  = trim(rttovDir)//trim(rttov_coefDir)//trim(PC_coef_filepath)
        
        call rttov_read_coefs(errorstatus, coef_rttov, opts,    &
                              file_coef=OD_coef_filepath,       &
!                              file_scaer=aer_coef_filepath,     & ! Needs to be PC-RTTOV compatible
                              file_pccoef=PC_coef_filepath)    
    else ! Read optical depth and cloud coefficient files together
        call rttov_read_coefs(errorstatus, coef_rttov, opts,    &
                              file_coef=OD_coef_filepath,       &
                              file_scaer=aer_coef_filepath,     &
                              file_sccld=cld_coef_filepath)
        ! Ensure input number of channels is not higher than number stored in coefficient file
        if (nchannels_rec > coef_rttov % coef % fmv_chn) then
            nchannels_rec = coef_rttov % coef % fmv_chn
            print*,'nchannels_rec cap hit'
        endif            
    endif
                          
    ! We aren't checking an allocation steps so this seems more appropriate.
    call rttov_error('fatal error reading coefficients' , lalloc = .false.)

    ! Ensure the options and coefficients are consistent
    call rttov_user_options_checkinput(errorstatus, opts, coef_rttov)
    
    ! We aren't checking an allocation steps so this seems more appropriate.
    call rttov_error('error in rttov options' , lalloc = .false.)


    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Figure out how many channels we actually want to reconstruct
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Handle different radiance reconstruction options
    if (nchannels_rec < 0) then
        print*,'The namelist variable "nchannels_rec" is negative, rttov_direct call will fail. Exiting.'
        errorstatus = errorstatus_fatal
        call rttov_exit(errorstatus)
        ! If the number of channels is negative, don't reconstruct radiances at all
        nchan_out = 0
        nchannels_rec = 0 ! Avoid nchanprof set to a negative value
    else if (nchannels_rec == 0) then
        ! If the number of channels is set to 0 then reconstruct all instrument channels
        nchan_out = coef_rttov % coef % fmv_chn
        nchannels_rec = coef_rttov % coef % fmv_chn ! Avoid nchanprof set to 0
    else
        ! Otherwise read the channel list from the file
        nchan_out = nchannels_rec
    endif


    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Read in channel indices, emissivities, and reflectivities from .csv if file is passed
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    print*,'nchan_out:   ',nchan_out
    
    if (Lchannel_filepath) then
        allocate(iChannel(nchan_out))
        allocate(emisChannel(nchan_out))
        allocate(reflChannel(nchan_out))
    
        open(18,file=channel_filepath,access='sequential',form="formatted")
        do i = 1, nchan_out
            read(18,*) iChannel(i), emisChannel(i), reflChannel(i)
        end do
        close(18)
    else ! If nothing is passed, compute the first "nchan_out" channels. Ignore emissivity and reflectivity for now.
        allocate(iChannel(nchan_out))
        iChannel(:) = (/ (i, i = 1, nchan_out) /)
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
 
  END SUBROUTINE COSP_RTTOV_INIT
    
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_rttov_simulate - Call subroutines in mod_cosp_rttov to run RTTOV
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_RTTOV_SIMULATE_MI(rttovIN,rttovConfig,lCleanup,                  & ! Inputs
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
    real(wp),intent(inout),dimension(rttovIN%nPoints,rttovIN%nChannels) :: & ! Can I do this? I guess so!
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
        call COSP_PC_RTTOV_SIMULATE_MI(rttovIN,rttovConfig,lCleanup,                                 &
                                    bt_total,rad_total,                               &
                                    error)                                

        print*,'end of COSP_PC_RTTOV_SIMULATE_MI'
    else
        print*,'Running COSP_REG_RTTOV_SIMULATE in multi-inst. set-up'
        call COSP_REG_RTTOV_SIMULATE_MI(rttovIN,rttovConfig,lCleanup,                                 &
                                     bt_total,bt_clear,                                &
                                     rad_total,rad_clear,rad_cloudy,                   &
                                     refl_total,refl_clear,                            &
                                     error)
        
    endif

  END SUBROUTINE COSP_RTTOV_SIMULATE_MI

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_rttov_simulate - Call regular subroutines in mod_cosp_rttov to run RTTOV
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_REG_RTTOV_SIMULATE_MI(rttovIN,rttovConfig,lCleanup,                  & ! Inputs
                                     bt_total,bt_clear,                                & ! Brightness Temp Outputs
                                     rad_total,rad_clear,rad_cloudy,                   & ! Radiance Outputs
                                     refl_total,refl_clear,                            & ! Reflectance Outputs
                                     error)                                              
  
    use mod_cosp_rttov,             only:        &
        cosp_rttov_allocate_mi,                  &
        cosp_rttov_construct_profiles_mi,        &
        cosp_rttov_setup_emissivity_reflectance, &
        cosp_rttov_call_direct_mi,               &
        cosp_rttov_save_output_mi,               &
        cosp_rttov_deallocate_profiles_mi,       &
        cosp_rttov_deallocate_coefs_mi
  
    type(rttov_in),intent(in) :: &
        rttovIN
    type(rttov_cfg),intent(inout) :: &
        rttovConfig
    logical,intent(in) :: &
         lCleanup   ! Flag to determine whether to deallocate RTTOV types
    real(wp),intent(inout),dimension(rttovIN%nPoints,rttovIN%nChannels) :: & ! Can I do this? I guess so!
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
    call cosp_rttov_allocate_mi(rttovIN,                       &
                                rttovConfig % nChannels_rec,   &
                                rttovConfig % opts,            &
                                rttovConfig % coefs,           &
                                rttovConfig % iChannel)
        
    print*,'cosp_rttov_allocate successful' ! jks
    call cpu_time(driver_time(2))
    call cosp_rttov_construct_profiles_mi(rttovIN, &
                                          rttovConfig % opts % rt_ir % addclouds, &
                                          rttovConfig % opts % rt_ir % addaerosl)
    print*,'cosp_rttov_construct_profiles_mi successful' ! jks
    call cpu_time(driver_time(3))
    call cosp_rttov_setup_emissivity_reflectance() ! Config agnostic after allocate step.
    print*,'cosp_rttov_setup_emissivity_reflectance successful' ! jks
    call cpu_time(driver_time(4))
    call cosp_rttov_call_direct_mi(rttovConfig % rttov_direct_nthreads,  &
                                   rttovConfig % opts,                   &
                                   rttovConfig % coefs)    
    
    print*,'cosp_rttov_call_direct successful' ! jks
    call cpu_time(driver_time(5))
    
    call cosp_rttov_save_output_mi(rttovIN,                        &
                                rttovConfig % nchan_out,                &
                                rttovConfig % Lrttov_bt,                &
                                rttovConfig % Lrttov_rad,               &
                                rttovConfig % Lrttov_refl,              &
                                rttovConfig % opts % rt_ir % addclouds, &
                                rttovConfig % opts % rt_ir % addaerosl, &
                                bt_total,bt_clear,                      &
                                rad_total,rad_clear,rad_cloudy,         &
                                refl_total,refl_clear)

    print*,'cosp_rttov_save_output successful' ! jks
    call cpu_time(driver_time(6))
    call cosp_rttov_deallocate_profiles_mi(rttovIN,                       &
                                           rttovConfig % nChannels_rec,   &
                                           rttovConfig % opts,            &
                                           rttovConfig % coefs)    
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
        call cosp_rttov_deallocate_coefs_mi(rttovConfig % coefs)
        call cpu_time(driver_time(9))
        print*,'Time to run "cosp_rttov_deallocate_coefs":     ',driver_time(9)-driver_time(8)
    endif

!    print*,'Total RTTOV run time:     ',driver_time(8)-driver_time(1)

  END SUBROUTINE COSP_REG_RTTOV_SIMULATE_MI


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_pc_rttov_simulate - Call subroutines in mod_cosp_rttov to run RTTOV
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_PC_RTTOV_SIMULATE_MI(rttovIN,rttovConfig,lCleanup,                  & ! Inputs
                                    bt_total,rad_total,                               & ! Outputs
                                    error)                                              
  
    use mod_cosp_rttov,             only:   &
        cosp_pc_rttov_allocate_mi,          &
        cosp_rttov_construct_profiles_mi,   &
        cosp_pc_rttov_setup_emissivity,     &
        cosp_pc_rttov_call_direct_mi,       &
        cosp_pc_rttov_save_output_mi,       &
        cosp_pc_rttov_deallocate_profiles_mi,  &
        cosp_rttov_deallocate_coefs_mi
  
    type(rttov_in),intent(in) :: &
        rttovIN
    type(rttov_cfg),intent(inout) :: &
        rttovConfig        
    logical,intent(in) :: &
         lCleanup   ! Flag to determine whether to deallocate RTTOV types
    real(wp),intent(inout),dimension(rttovIN%nPoints,rttovIN%nChannels) :: & ! Can I do this? I guess so!
        bt_total,                          &        ! All-sky
        rad_total                                   ! All-sky
    character(len=128) :: &
        error     ! Error messages (only populated if error encountered)  

    real(wp),dimension(10) :: driver_time

    ! Run each step for running RTTOV from mod_cosp_rttov (and time them)
!    print*,'cosp_rttov_allocate begin' ! jks
    call cpu_time(driver_time(1))
    call cosp_pc_rttov_allocate_mi(rttovIN, &
                                   rttovConfig%PC_coef_filepath,              &
                                   rttovConfig%coefs,                         &
                                   rttovConfig%opts,                          &
                                   rttovConfig%nchannels_rec,                 &
                                   rttovConfig%iChannel,                      &
                                   rttovConfig%iChannel_out)    
    print*,'cosp_pc_rttov_allocate_mi successful' ! jks
    call cpu_time(driver_time(2))
    call cosp_rttov_construct_profiles_mi(rttovIN, &
                                          rttovConfig % opts % rt_ir % addclouds, &
                                          rttovConfig % opts % rt_ir % addaerosl)    
    print*,'cosp_rttov_construct_profiles_mi successful' ! jks
    call cpu_time(driver_time(3))
    call cosp_pc_rttov_setup_emissivity()
    print*,'cosp_pc_rttov_setup_emissivity successful' ! jks
    call cpu_time(driver_time(4))
    call cosp_pc_rttov_call_direct_mi(rttovConfig % rttov_direct_nthreads,  &
                                      rttovConfig % opts,                   &
                                      rttovConfig % coefs,                  &
                                      rttovConfig % nchannels_rec,          &
                                      rttovConfig % iChannel_out) ! iChannel_out should have been updated

    print*,'cosp_pc_rttov_call_direct successful' ! jks
    call cpu_time(driver_time(5))
    call cosp_pc_rttov_save_output_mi(rttovIN%nPoints,                       &
                                      rttovConfig%nchannels_rec,             &
                                      rttovConfig%Lrttov_bt,                 &
                                      rttovConfig%Lrttov_rad,                &
                                      bt_total,                              &
                                      rad_total)

    print*,'cosp_pc_rttov_save_output successful' ! jks
    call cpu_time(driver_time(6))
    call cosp_pc_rttov_deallocate_profiles_mi(rttovIN,                       &
                                              rttovConfig % nChannels_rec,   &
                                              rttovConfig % opts,            &
                                              rttovConfig % coefs)
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
!        call cosp_rttov_deallocate_coefs_mi(rttovConfig % coefs) ! JKS this again with PC-RTTOV...
        call cpu_time(driver_time(9))
        print*,'Time to run "cosp_rttov_deallocate_coefs":     ',driver_time(9)-driver_time(8)
    endif

  END SUBROUTINE COSP_PC_RTTOV_SIMULATE_MI

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! END MODULE
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE MOD_COSP_RTTOV_INTERFACE
