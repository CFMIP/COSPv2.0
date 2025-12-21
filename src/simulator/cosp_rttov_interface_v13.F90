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
! Jun 2025 - J.K. Shaw - Modified for RTTOV v13
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE MOD_COSP_RTTOV_INTERFACE
  USE COSP_KINDS,          ONLY: wp
  use MOD_COSP_RTTOV,      ONLY: rttov_in                              
  USE MOD_COSP_RTTOV_UTIL, ONLY: rttov_cfg,             rttov_output
  use MOD_COSP_ERROR,      ONLY: errorMessage

  ! rttov_const contains useful RTTOV constants
  USE rttov_const, ONLY :     &
         errorstatus_success, &
         errorstatus_fatal

  ! rttov_types contains definitions of all RTTOV data types
  USE rttov_types, ONLY :     &
         rttov_options,       &
         rttov_options_scatt, &
         rttov_coefs,         &
         rttov_scatt_coef,    &
         rttov_pccomp,        &
         rttov_radiance,      &
         rttov_transmission,  &
         rttov_profile,       &
         rttov_emissivity,    &
         rttov_reflectance,   &
         rttov_chanprof

  ! jpim, jprb and jplm are the RTTOV integer, real and logical KINDs
  USE parkind1, ONLY : jpim, jprb, jplm
  
  USE rttov_unix_env, ONLY : rttov_exit
                              
  IMPLICIT NONE

#include "rttov_read_coefs.interface"
#include "rttov_user_options_checkinput.interface"
#include "rttov_print_opts.interface"


CONTAINS

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_rttov_init
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_RTTOV_INIT(Lrttov, Nlevels,Ninstruments,instrument_namelists,       &
                             rttov_configs, unitn, debug)

      logical,intent(in) :: &
          Lrttov
      integer,intent(in) :: &
          Nlevels,   &
          Ninstruments
      type(character(len=256)), dimension(Ninstruments)     :: & 
          instrument_namelists   ! Array of paths to RTTOV instrument namelists      
      type(rttov_cfg), dimension(:), intent(out), allocatable :: &
          rttov_configs
      integer,intent(in),Optional :: unitn
      logical,intent(in),Optional :: debug
          
      ! Local variables
      integer            :: &
          inst_idx ! iterator
      logical :: verbose
      if (present(debug)) verbose = debug
        
      allocate(rttov_configs(Ninstruments))
        
      ! Create config objects for each instrument to be simulated by RTTOV. Return to the main subroutine.
      do inst_idx=1,Ninstruments
          if (present(unitn)) then
              call cosp_rttov_init_s(Nlevels,instrument_namelists(inst_idx),rttov_configs(inst_idx),unitn=unitn,debug=verbose)
          else
              call cosp_rttov_init_s(Nlevels,instrument_namelists(inst_idx),rttov_configs(inst_idx),debug=verbose)
          endif
      end do
         
       
  END SUBROUTINE COSP_RTTOV_INIT

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_rttov_init
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_RTTOV_INIT_S(Nlevels,namelist_filepath,       &
                               rttov_config, unitn, debug)
  
     integer,intent(in)                :: &
         Nlevels
    character(len=256),intent(in)     :: & 
         namelist_filepath   ! Array of paths to RTTOV instrument namelists      
    type(rttov_cfg),intent(out)       :: &
         rttov_config 
    
    integer,intent(in),Optional :: unitn ! Used for io limits
    logical,intent(in),Optional :: debug
    
    ! Local variables
    character(len=256),target :: &
        channel_filepath,  &
        wavenum_filepath,  &
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
    logical         :: Lrttov_bt
    logical         :: Lrttov_rad
    logical         :: Lrttov_refl
    logical         :: Lrttov_cld
    logical         :: Lrttov_aer
    logical         :: Lrttov_cldparam
    logical         :: Lrttov_aerparam
    logical         :: Lrttov_gridbox_cldmmr
    logical         :: Ldo_nlte_correction
    logical         :: Lrttov_pc
    logical         :: Lrttov_solar
    logical         :: Lchannel_filepath
    logical         :: Lwavenum_filepath
    logical         :: user_tracegas_input 
    logical         :: SO2_data
    logical         :: N2O_data
    logical         :: CO_data
    logical         :: CO2_data
    logical         :: CH4_data
    logical         :: ozone_data
    logical         :: clw_data
    
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
        
    ! For orbital swathing
    integer(kind=jpim)     ::  &   
        rttov_Nlocaltime,          & ! Number of orbits
        rttov_gas_units,           & ! RTTOV units for trace gases: 0 ppmv over dry air, 1 kg/kg over moist air, 2 ppmv over moist air, 3 kg/kg over dry air
        rttov_clw_scheme,          & ! Scheme for determining cloud water optical properties.
        rttov_ice_scheme,          & ! Scheme for determining cloud ice optical properties.
        rttov_icede_param,         & ! Scheme for how cloud ice deff is parameterized
        rttov_extendatmos
    real(wp),dimension(20) ::  &         ! Reasonable but arbitrary limit at 10 local time orbits
        rttov_localtime,           & ! RTTOV subsetting by local time in hours [0,24]
        rttov_localtime_width        ! Width of satellite swath (km).
                
    ! For checking errors in filenames.
    character(len=256) :: imsg  !<-- some suitable length, say XX=256      
    integer            :: erro
    integer(kind=jplm) :: errorstatus              ! Return error status of RTTOV subroutine calls
    integer(kind=jpim) :: alloc_status(60)
               
    logical :: verbose = .false.
    character(len=32) :: flag
    character(len=80) :: msg
    if (present(debug)) verbose = debug
               
    ! Init. variables to false.
    rttov_Nlocaltime      = 0
    Lrttov_bt             = .false.
    Lrttov_rad            = .false.
    Lrttov_refl           = .false.
    Lrttov_cld            = .false.
    Lrttov_aer            = .false.
    Lrttov_cldparam       = .false.
    Lrttov_aerparam       = .false.
    Lrttov_gridbox_cldmmr = .true. ! Assume gridbox average MMRs. Most common for GCMs.
    Ldo_nlte_correction   = .false. ! Correct for non-local thermal equilibrium effects? Default false.
    Lrttov_pc             = .false.
    Lrttov_solar          = .false.
    Lchannel_filepath     = .false.
    Lwavenum_filepath     = .false.
    SO2_data              = .false. 
    N2O_data              = .false. 
    CO_data               = .false.
    CO2_data              = .false.
    CH4_data              = .false. 
    ozone_data            = .false.
    clw_data              = .false.    
    user_tracegas_input   = .false.
    rttov_Nlocaltime      = 0       ! Default: zero swath masking
    rttov_gas_units       = 1       ! Default: kg/kg over moist air (should be updated by user!)
    rttov_clw_scheme      = 2       ! 1: OPAC, 2: Deff
    rttov_ice_scheme      = 1       ! 1: Baum, 2: Baran (2014), 3: Baran (2018)
    rttov_icede_param     = 0       ! 0: Indicates that Deff is supplied but is rejected by RTTOV. 1: Ou and Liou, 2: Wyser(recommended), 3: Boudala, 4: McFarquhar.
    rttov_extendatmos     = 0       ! 0: do not extend above supplied pressure levels. 1: Simply top layer. 2: Not yet implemented.
    
    ! Read RTTOV namelist fields
    namelist/RTTOV_INPUT/Lrttov_bt,Lrttov_rad,Lrttov_refl,Lrttov_cld,            & ! Logicals for RTTOV configuration
                         Lrttov_aer,Lrttov_cldparam,Lrttov_aerparam,             & ! 
                         Lrttov_gridbox_cldmmr,Ldo_nlte_correction,              & ! Assume cloud water mixing ratios are gridbox average instead of in-cloud
                         Lrttov_pc,Lrttov_solar,nchannels_rec,Lchannel_filepath, &
                         channel_filepath,Lwavenum_filepath,wavenum_filepath,    &
                         rttov_srcDir,rttov_coefDir,            &
                         OD_coef_filepath,aer_coef_filepath,cld_coef_filepath,   &
                         PC_coef_filepath,                                       &
                         CO2_data,CH4_data,CO_data,N2O_data,SO2_data,ozone_data, & ! Use trace gases for radiative transfer?
                         clw_data,                                               & ! MW option
                         user_tracegas_input,                                    & ! User-supplied trace gas concentrations
                         CO2_mr,CH4_mr,CO_mr,N2O_mr,SO2_mr,                      & ! Mixing ratios
                         ipcbnd,ipcreg,npcscores,                                & ! PC-RTTOV config values
                         rttov_nthreads,rttov_ZenAng,rttov_Nlocaltime,           &
                         rttov_localtime,rttov_localtime_width,                  &
                         rttov_gas_units,rttov_clw_scheme,rttov_ice_scheme,      &
                         rttov_icede_param,rttov_extendatmos

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Read in namelists
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Handle indices of files already opened (for CESM integration)    
    if (present(unitn)) then
        open(unitn,file=namelist_filepath,status='unknown',iostat=erro,iomsg=imsg)
        if (erro > 0) then
            call errorMessage('Error reading in "namelist_filepath" in COSP_RTTOV_INIT_S')
            call errorMessage('erro:              ' // CHAR(erro))
            call errorMessage('imsg:              ' // imsg)
            call errorMessage('namelist_filepath: ' // namelist_filepath)
            errorstatus = 1
            call rttov_exit(errorstatus)                
        end if
        read(unitn,nml=RTTOV_INPUT)
        close(unitn)    
    else
        open(10,file=namelist_filepath,status='unknown',iostat=erro,iomsg=imsg)
        if (erro > 0) then
            call errorMessage('Error reading in "namelist_filepath" in COSP_RTTOV_INIT_S')
            call errorMessage('erro:              ' // CHAR(erro))
            call errorMessage('imsg:              ' // imsg)
            call errorMessage('namelist_filepath: ' // namelist_filepath)
            errorstatus = 1
            call rttov_exit(errorstatus)                
        end if      
        read(10,nml=RTTOV_INPUT)
        close(10)
    endif
    
    ! Set swath arrays correctly.
    allocate(rttov_config%rttov_localtime(rttov_Nlocaltime),rttov_config%rttov_localtime_width(rttov_Nlocaltime))
    rttov_config%rttov_Nlocaltime         = rttov_Nlocaltime
    rttov_config%rttov_localtime(:)       = rttov_localtime(1:rttov_Nlocaltime)
    rttov_config%rttov_localtime_width(:) = rttov_localtime_width(1:rttov_Nlocaltime)

    ! Extend atmosphere setting. If user-supplied values end too low, channels sounding the upper atmosphere will be bad.
    rttov_config%rttov_extendatmos = rttov_extendatmos

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
    
    ! Set to false in namelist if model supplies trace profiles
    rttov_config%user_tracegas_input = user_tracegas_input
    rttov_config%gas_units = rttov_gas_units

    ! Parametrization of cloud optical properties.
    rttov_config%clw_scheme   = rttov_clw_scheme
    rttov_config%ice_scheme   = rttov_ice_scheme
    rttov_config%icede_param  = rttov_icede_param

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
    rttov_config % opts % config % verbose          = .false.
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
    rttov_config % opts % rt_mw % clw_data             = .false.
    rttov_config % opts % rt_mw % clw_scheme           = 2       ! Default = 2/Rosenkranz
    rttov_config % opts % rt_mw % clw_cloud_top        = 322     ! Default is 322 hPa
    rttov_config % opts % rt_mw % fastem_version       = 6       ! Default FASTEM-6
    rttov_config % opts % rt_mw % supply_foam_fraction = .false.

    ! UV/visible/IR-only radiative transfer options
    rttov_config % opts % rt_ir % addsolar                = Lrttov_solar
    rttov_config % Lrttov_solar                           = Lrttov_solar
    rttov_config % opts % rt_ir % rayleigh_max_wavelength = 2._wp ! 2um 
    rttov_config % opts % rt_ir % rayleigh_min_pressure   = 0._wp ! 0 hPa
    rttov_config % opts % rt_ir % rayleigh_single_scatt   = .true.
    rttov_config % opts % rt_ir % rayleigh_depol          = .true. ! Default false, recommended true
    rttov_config % opts % rt_ir % do_nlte_correction      = Ldo_nlte_correction
    rttov_config % opts % rt_ir % solar_sea_brdf_model    = 2
    rttov_config % opts % rt_ir % ir_sea_emis_model       = 2

    ! User options
    ! Duplicate for STUB functionality
    rttov_config % Lrttov_aer                             = Lrttov_aer
    rttov_config % Lrttov_cld                             = Lrttov_cld
    rttov_config % opts % rt_ir % addaerosl               = Lrttov_aer
    rttov_config % opts % rt_ir % addclouds               = Lrttov_cld
    rttov_config % opts % rt_ir % user_aer_opt_param      = Lrttov_aerparam ! User specifies the aerosol scattering optical parameters 
    rttov_config % opts % rt_ir % user_cld_opt_param      = Lrttov_cldparam ! User specifies the cloud scattering optical parameters 
    
    rttov_config % opts % rt_ir % grid_box_avg_cloud      = Lrttov_gridbox_cldmmr
    rttov_config % opts % rt_ir % cldcol_threshold        = -1._wp
    rttov_config % opts % rt_ir % cloud_overlap           = 1 ! Maximum-random overlap
    rttov_config % opts % rt_ir % cc_low_cloud_top        = 750_wp ! 750 hPa. Only applies when cloud_overlap=2.
    rttov_config % opts % rt_ir % ir_scatt_model          = 2
    rttov_config % opts % rt_ir % vis_scatt_model         = 1 ! Scattering model to use for solar source term: 1 => DOM; 2 => single-scattering; 3 => MFASIS-LUT; 4 => MFASIS-NN (default = 1); only applies when addclouds or addaerosl is true and addsolar is true. JKS note, DOM is the most expensive, MFASIS-NN might be a better option
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
    
    if (rttov_config % Lrttov_mwscatt) then
        rttov_config % opts_scatt % config % do_checkinput = .true.
        rttov_config % opts_scatt % config % apply_reg_limits = .false.
        rttov_config % opts_scatt % config % verbose = .true.
        rttov_config % opts_scatt % ozone_data = .false.  ! Default
        rttov_config % opts_scatt % use_t2m_opdep = .true.
        rttov_config % opts_scatt % use_q2m = .true.
        rttov_config % opts_scatt % use_tskin_eff = .false.
        rttov_config % opts_scatt % addrefrac = .true.
        rttov_config % opts_scatt % rad_down_lin_tau = .false. ! Recommended
        rttov_config % opts_scatt % interp_mode = 1 ! Default
        rttov_config % opts_scatt % lgradp = .false.
        rttov_config % opts_scatt % fastem_version = 6 ! Default
        rttov_config % opts_scatt % supply_foam_fraction = .false.
        rttov_config % opts_scatt % lusercfrac = .false. ! User supplied cloud fraction in rttov_profile_cloud. Maybe set to true? pg 164
        rttov_config % opts_scatt % cc_threshold = 0.001_wp ! Default
!        rttov_config % opts_scatt % pol_mode = 
!        rttov_config % opts_scatt % ice_polarisation =
        rttov_config % opts_scatt % hydro_cfrac_tlad = .true. ! Default
        rttov_config % opts_scatt % zero_hydro_tlad = .false. ! Default
    end if

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
            if (verbose) call errorMessage('nchannels_rec cap hit')
        end if            
    end if   
                          
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
    if (nchannels_rec < 0) then
        call errorMessage('The namelist variable "nchannels_rec" is negative, rttov_direct call will fail. Exiting.')
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
    
    allocate(rttov_config % iChannel(rttov_config % nchan_out)) ! There is a need for these variables to be separate somewhere...
    allocate(rttov_config % iChannel_out(rttov_config % nchan_out))
    if (Lwavenum_filepath) then
        allocate(rttov_config % wavenumChannel(rttov_config % nchan_out))
        open(18,file=wavenum_filepath,access='sequential',form="formatted")
        do i = 1, rttov_config % nchan_out
            read(18,*) rttov_config % wavenumChannel(i)
        end do
        close(18)
    end if
    if (Lchannel_filepath) then
        allocate(rttov_config % emisChannel(rttov_config % nchan_out))
        allocate(rttov_config % reflChannel(rttov_config % nchan_out))
    
        open(18,file=channel_filepath,access='sequential',form="formatted")
        do i = 1, rttov_config % nchan_out
            read(18,*) rttov_config % iChannel(i), rttov_config % emisChannel(i), rttov_config % reflChannel(i)
        end do
        close(18)
        rttov_config % iChannel_out = rttov_config % iChannel
    else ! If nothing is passed, compute the first "nchan_out" channels. Ignore emissivity and reflectivity for now.
        rttov_config % iChannel(:) = (/ (i, i = 1, rttov_config % nchan_out) /)
        rttov_config % iChannel_out = rttov_config % iChannel
    endif
    
    if (verbose) then
        call errorMessage('RTTOV configuration:')
        call errorMessage('rttov_config % nchan_out:         ' // trim(adjustl(CHAR(rttov_config % nchan_out))))
        write(flag, '(L1)') rttov_config%Lrttov_bt
        call errorMessage('rttov_config % Lrttov_bt:         ' // trim(flag))
        write(flag, '(L1)') rttov_config%Lrttov_rad
        call errorMessage('rttov_config % Lrttov_rad:        ' // trim(flag))
        write(flag, '(L1)') rttov_config%Lrttov_refl
        call errorMessage('rttov_config % Lrttov_refl:       ' // trim(flag))
        write(flag, '(L1)') rttov_config%Lrttov_cld
        call errorMessage('rttov_config % Lrttov_cld:        ' // trim(flag))
        write(flag, '(L1)') rttov_config%Lrttov_aer
        call errorMessage('rttov_config % Lrttov_aer:        ' // trim(flag))
        write(flag, '(L1)') rttov_config%Lrttov_pc
        call errorMessage('rttov_config % Lrttov_pc:         ' // trim(flag))
        write(flag, '(L1)') rttov_config%Lrttov_solar
        call errorMessage('rttov_config % Lrttov_solar:      ' // trim(flag))
        write(flag, '(L1)') rttov_config % opts % rt_ir % grid_box_avg_cloud
        call errorMessage('rttov_config % opts % rt_ir % grid_box_avg_cloud:   ' // trim(flag))
        write(flag, '(L1)') rttov_config % opts % rt_ir % do_nlte_correction
        call errorMessage('rttov_config % opts % rt_ir % do_nlte_correction:     ' // trim(flag))
        call errorMessage('rttov_config % rttov_Nlocaltime:        ' // trim(adjustl(CHAR(rttov_config % rttov_Nlocaltime))))
        write(msg, '("shape(inst_chanprof%prof):", *(1x, i0))') rttov_config % rttov_localtime
        call errorMessage('rttov_config % rttov_localtime:         ' // trim(msg))
        write(msg, '("shape(inst_chanprof%prof):", *(1x, i0))') rttov_config % rttov_localtime_width
        call errorMessage('rttov_config % rttov_localtime_width:   ' // trim(msg))
        call errorMessage('rttov_config % rttov_extendatmos:       ' // trim(CHAR(rttov_config % rttov_extendatmos)))
        call errorMessage('rttov_config % gas_units:               ' // trim(CHAR(rttov_config % gas_units)))
        call errorMessage('rttov_config % clw_scheme:              ' // trim(CHAR(rttov_config % clw_scheme)))
        call errorMessage('rttov_config % ice_scheme:              ' // trim(CHAR(rttov_config % ice_scheme)))
        call errorMessage('rttov_config % icede_param:             ' // trim(CHAR(rttov_config % icede_param)))
        call rttov_print_opts(rttov_config % opts)
    end if

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
  
  
  SUBROUTINE DESTROY_RTTOV_CONFIG(rttovConfig)
      use mod_cosp_rttov,             only:   cosp_rttov_deallocate_coefs
  
      type(rttov_cfg),intent(inout) :: &
          rttovConfig
          
      if (allocated(rttovConfig % iChannel))              deallocate(rttovConfig % iChannel)
      if (allocated(rttovConfig % iChannel_out))          deallocate(rttovConfig % iChannel_out)
      if (allocated(rttovConfig % emisChannel))           deallocate(rttovConfig % emisChannel)
      if (allocated(rttovConfig % reflChannel))           deallocate(rttovConfig % reflChannel)
      if (allocated(rttovConfig % rttov_localtime))       deallocate(rttovConfig % rttov_localtime)
      if (allocated(rttovConfig % rttov_localtime_width)) deallocate(rttovConfig % rttov_localtime_width)
      if (allocated(rttovConfig % swath_mask))            deallocate(rttovConfig % swath_mask)
  
      call cosp_rttov_deallocate_coefs(rttovConfig % coefs)

  END SUBROUTINE DESTROY_RTTOV_CONFIG
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_rttov_simulate - Call subroutines in mod_cosp_rttov to run RTTOV
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_RTTOV_SIMULATE(rttovIN,rttovConfig,error,                        & ! Inputs
                                 bt_total,bt_clear,                                & ! Brightness Temp Outputs
                                 rad_total,rad_clear,rad_cloudy,                   & ! Radiance Outputs
                                 refl_total,refl_clear,                            & ! Reflectance Outputs
                                 debug)

    use mod_cosp_rttov,             only: cosp_rttov_swath                                 

    type(rttov_in),intent(in) :: &
        rttovIN
    type(rttov_cfg),intent(inout) :: &
        rttovConfig
    character(len=128),intent(inout) :: &
        error     ! Error messages (only populated if error encountered)         
    real(wp),intent(inout),dimension(rttovIN%nPoints,rttovConfig%nchan_out),optional :: & ! Can I do this? I guess so! 
        bt_total,                          &        ! All-sky
        bt_clear,                          &        ! Clear-sky
        rad_total,                         &        ! All-sky
        rad_clear,                         &        ! Clear-sky
        rad_cloudy,                        &        ! Cloudy-sky
        refl_total,                        &        ! All-sky
        refl_clear                                  ! Clear-sky
    logical,intent(in),optional :: &
        debug
        
    logical :: verbose = .false.
    if (present(debug)) verbose = debug

    if (allocated(rttovConfig % swath_mask)) deallocate(rttovConfig % swath_mask)
    allocate(rttovConfig % swath_mask(rttovIN % nPoints))

    call cosp_rttov_swath(rttovIN,                                &
                          rttovConfig % rttov_Nlocaltime,         &
                          rttovConfig % rttov_localtime,          &
                          rttovConfig % rttov_localtime_width,    &
                          rttovConfig % swath_mask,               &
                          debug)
    rttovConfig % nprof = count(rttovConfig % swath_mask)    

    if (rttovConfig % nprof .gt. 0) then ! Skip calculations if all values are swathed out
        ! Check options to determine if the principal component approach should be run
        if (rttovConfig % opts % rt_ir % pc % addpc) then
            call COSP_PC_RTTOV_SIMULATE(rttovIN,rttovConfig,                              &
                                        bt_clear,rad_clear,                               &
                                        error,verbose)                                
        else
            call COSP_REG_RTTOV_SIMULATE(rttovIN,rttovConfig,                             &
                                        bt_total,bt_clear,                               &
                                        rad_total,rad_clear,rad_cloudy,                  &
                                        refl_total,refl_clear,                           &
                                        error,verbose)
        endif
    else
        if (verbose) call errorMessage('empty chunk')
    endif

  END SUBROUTINE COSP_RTTOV_SIMULATE

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_reg_rttov_simulate - Call regular subroutines in mod_cosp_rttov to run RTTOV
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_REG_RTTOV_SIMULATE(rttovIN,rttovConfig,                              & ! Inputs
                                     bt_total,bt_clear,                                & ! Brightness Temp Outputs
                                     rad_total,rad_clear,rad_cloudy,                   & ! Radiance Outputs
                                     refl_total,refl_clear,                            & ! Reflectance Outputs
                                     error,verbose)                                              
  
    use mod_cosp_rttov,             only:        &
        cosp_rttov_allocate,                     &
        cosp_rttov_construct_profiles,           &
        cosp_rttov_setup_emissivity_reflectance, &
        cosp_rttov_call_direct,                  &
        cosp_rttov_save_output,                  &
        cosp_rttov_deallocate_profiles
  
    type(rttov_in),intent(in) :: &
        rttovIN
    type(rttov_cfg),intent(inout) :: &
        rttovConfig
    real(wp),intent(inout),dimension(rttovIN%nPoints,rttovConfig%nchan_out) :: & 
        bt_total,                          &        ! All-sky
        bt_clear,                          &        ! Clear-sky
        rad_total,                         &        ! All-sky
        rad_clear,                         &        ! Clear-sky
        rad_cloudy,                        &        ! Cloudy-sky
        refl_total,                        &        ! All-sky
        refl_clear                                  ! Clear-sky
    character(len=128),intent(inout) :: &
        error     ! Error messages (only populated if error encountered)
        
    ! Local variables
    type(rttov_radiance)         :: radiance
    type(rttov_transmission)     :: transmission
    type(rttov_profile),     pointer :: profiles(:)    => NULL() ! Input profiles
    logical(kind=jplm),      pointer :: calcemis(:)    => NULL() ! Flag to indicate calculation of emissivity within RTTOV
    type(rttov_emissivity),  pointer :: emissivity(:)  => NULL() ! Input/output surface emissivity
    logical(kind=jplm),      pointer :: calcrefl(:)    => NULL() ! Flag to indicate calculation of BRDF within RTTOV
    type(rttov_reflectance), pointer :: reflectance(:) => NULL() ! Input/output surface BRDF
    type(rttov_chanprof),    pointer :: chanprof(:)    => NULL() ! Input channel/profile list

    logical,intent(in) :: verbose

    real(wp),dimension(10) :: driver_time

    ! Run each step for running RTTOV from mod_cosp_rttov (and time them)
    call cpu_time(driver_time(1))
    call cosp_rttov_allocate(rttovIN,                             &
                            rttovConfig % nChannels_rec,         &
                            rttovConfig % opts,                  &
                            rttovConfig % coefs,                 &
                            profiles,                            &
                            rttovConfig % iChannel,              &
                            chanprof,                            &
                            rttovConfig % nchanprof,             &
                            rttovConfig % nprof,                 &
                            rttovConfig % swath_mask,            &
                            rttovConfig % rttov_extendatmos,     &
                            transmission,                        &
                            radiance,                            &
                            calcemis,                            &
                            emissivity,                          &
                            calcrefl,                            &
                            reflectance,                         &
                            verbose)

    call cpu_time(driver_time(2))
    call cosp_rttov_construct_profiles(rttovIN,                                &
                                    profiles,                               &
                                    rttovConfig % Lrttov_cld,               &
                                    rttovConfig % Lrttov_aer,               &
                                    rttovConfig % Lrttov_solar,             &
                                    rttovConfig % user_tracegas_input,      &
                                    rttovConfig % opts % rt_all % CO2_data, &
                                    rttovConfig % opts % rt_all % CH4_data, &
                                    rttovConfig % opts % rt_all % CO_data,  &
                                    rttovConfig % opts % rt_all % N2O_data, &
                                    rttovConfig % opts % rt_all % SO2_data, &
                                    rttovConfig % opts % rt_all % ozone_data, &
                                    rttovConfig % CO2_mr,                   &
                                    rttovConfig % CH4_mr,                   &
                                    rttovConfig % CO_mr,                    &
                                    rttovConfig % N2O_mr,                   &
                                    rttovConfig % SO2_mr,                   &
                                    rttovConfig % ZenAng,                   &
                                    rttovConfig % nprof,                    &
                                    rttovConfig % swath_mask,               &
                                    rttovConfig % gas_units,                &
                                    rttovConfig % clw_scheme,               &
                                    rttovConfig % ice_scheme,               &
                                    rttovConfig % icede_param,              &
                                    rttovConfig % rttov_extendatmos,        &
                                    verbose)
                                    
    call cpu_time(driver_time(3))
    
    if (associated(rttovIN % emis_grey)) then
        call cosp_rttov_setup_emissivity_reflectance(calcemis,    &
                                                    emissivity,  &
                                                    calcrefl,    &
                                                    reflectance, &
                                                    emis_grey = rttovIN % emis_grey) ! Config agnostic after allocate step.
    else
        call cosp_rttov_setup_emissivity_reflectance(calcemis,    &
                                                    emissivity,  &
                                                    calcrefl,    &
                                                    reflectance)
    end if
    call cpu_time(driver_time(4))
    
    call cosp_rttov_call_direct(rttovConfig % rttov_direct_nthreads,  &
                                rttovConfig % opts,                   &
                                profiles,                             &
                                rttovConfig % coefs,                  &
                                chanprof,                             &
                                transmission,                         &
                                radiance,                             &
                                calcemis,                             &
                                emissivity,                           &
                                calcrefl,                             &
                                reflectance,                          &
                                verbose)                                    
    
    call cpu_time(driver_time(5))
    
    call cosp_rttov_save_output(rttovIN % nPoints,                      &
                                rttovConfig % nchan_out,                &
                                rttovConfig % swath_mask,               &
                                rttovConfig % Lrttov_bt,                &
                                rttovConfig % Lrttov_rad,               &
                                rttovConfig % Lrttov_refl,              &
                                rttovConfig % Lrttov_cld,               &
                                rttovConfig % Lrttov_aer,               &
                                radiance,                               &
                                bt_total,bt_clear,                      &
                                rad_total,rad_clear,rad_cloudy,         &
                                refl_total,refl_clear)

    call cpu_time(driver_time(6))                                        
    call cosp_rttov_deallocate_profiles(rttovConfig % nprof,           &
                                        rttovConfig % nchanprof,       &
                                        rttovIN % nLevels,             &
                                        rttovConfig % opts,            &
                                        profiles,                      &
                                        rttovConfig % coefs,           &
                                        chanprof,                      &
                                        transmission,                  &
                                        radiance,                      &
                                        calcemis,                      &
                                        emissivity,                    &
                                        calcrefl,                      &
                                        reflectance)                                         
    call cpu_time(driver_time(7))

  END SUBROUTINE COSP_REG_RTTOV_SIMULATE


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_pc_rttov_simulate - Call subroutines in mod_cosp_rttov to run RTTOV
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_PC_RTTOV_SIMULATE(rttovIN,rttovConfig,                           & ! Inputs
                                    bt_clear,rad_clear,                            & ! Outputs
                                    error,verbose)                                              
  
    use mod_cosp_rttov,             only:   &
        cosp_pc_rttov_allocate,             &
        cosp_rttov_construct_profiles,      &
        cosp_pc_rttov_setup_emissivity,     &
        cosp_pc_rttov_call_direct,          &
        cosp_pc_rttov_save_output,          &
        cosp_pc_rttov_deallocate_profiles
  
    type(rttov_in),intent(in) :: &
        rttovIN
    type(rttov_cfg),intent(inout) :: &
        rttovConfig        
    real(wp),intent(inout),dimension(rttovIN%nPoints,rttovConfig%nchan_out) :: & ! Can I do this? I guess so!
        bt_clear,                          &        ! All-sky
        rad_clear                                   ! All-sky
    character(len=128),intent(inout) :: &
        error     ! Error messages (only populated if error encountered)

    ! Local variables
    type(rttov_radiance) :: radiance
    type(rttov_transmission) :: transmission
    type(rttov_profile),     pointer :: profiles(:)    => NULL() ! Input profiles
    logical(kind=jplm),      pointer :: calcemis(:)    => NULL() ! Flag to indicate calculation of emissivity within RTTOV
    type(rttov_emissivity),  pointer :: emissivity(:)  => NULL() ! Input/output surface emissivity
    logical(kind=jplm),      pointer :: calcrefl(:)    => NULL() ! Flag to indicate calculation of BRDF within RTTOV
    type(rttov_reflectance), pointer :: reflectance(:) => NULL() ! Input/output surface BRDF
    type(rttov_chanprof),    pointer :: chanprof(:)    => NULL() ! Input channel/profile list
    integer(KIND=jpim),      pointer :: predictindex(:)
    
    logical,intent(in) :: verbose
    real(wp),dimension(10) :: driver_time

    ! Run each step for running RTTOV from mod_cosp_rttov (and time them)
    call cpu_time(driver_time(1))
    call cosp_pc_rttov_allocate(rttovIN, &
                                rttovConfig % PC_coef_filepath,              &
                                rttovConfig % coefs,                         &
                                rttovConfig % opts,                          &
                                profiles,                                    &
                                rttovConfig % nchannels_rec,                 &
                                rttovConfig % iChannel,                      &
                                chanprof,                                    &
                                rttovConfig % nchanprof,                     &
                                rttovConfig % nprof,                         &
                                rttovConfig % iChannel_out,                  &
                                rttovConfig % swath_mask,                    &
                                rttovConfig % rttov_extendatmos,             &
                                transmission,                                &
                                radiance,                                    &
                                calcemis,                                    &
                                emissivity,                                  &
                                rttovConfig % pccomp,                        &
                                predictindex)
    call cpu_time(driver_time(2))

    call cosp_rttov_construct_profiles(rttovIN,                                &
                                    profiles,                               &
                                    rttovConfig % Lrttov_cld,               &
                                    rttovConfig % Lrttov_aer,               &
                                    rttovConfig % Lrttov_solar,             &
                                    rttovConfig % user_tracegas_input,      &
                                    rttovConfig % opts % rt_all % CO2_data, &
                                    rttovConfig % opts % rt_all % CH4_data, &
                                    rttovConfig % opts % rt_all % CO_data,  &
                                    rttovConfig % opts % rt_all % N2O_data, &
                                    rttovConfig % opts % rt_all % SO2_data, &
                                    rttovConfig % opts % rt_all % ozone_data, &
                                    rttovConfig % CO2_mr,                   &
                                    rttovConfig % CH4_mr,                   &
                                    rttovConfig % CO_mr,                    &
                                    rttovConfig % N2O_mr,                   &
                                    rttovConfig % SO2_mr,                   &
                                    rttovConfig % ZenAng,                   &
                                    rttovConfig % nprof,                    &
                                    rttovConfig % swath_mask,               &
                                    rttovConfig % gas_units,                &
                                    rttovConfig % clw_scheme,               &
                                    rttovConfig % ice_scheme,               &
                                    rttovConfig % icede_param,              &
                                    rttovConfig % rttov_extendatmos,        &
                                    verbose)
    call cpu_time(driver_time(3))
    call cosp_pc_rttov_setup_emissivity(calcemis,   &
                                        emissivity)
    call cpu_time(driver_time(4))
    call cosp_pc_rttov_call_direct(rttovConfig % rttov_direct_nthreads,  &
                                rttovConfig % opts,                   &
                                profiles,                             &
                                rttovConfig % coefs,                  &
                                chanprof,                             &
                                transmission,                         &
                                rttovConfig % nchannels_rec,          &
                                rttovConfig % iChannel_out,           &
                                radiance,                             &
                                calcemis,                             &
                                emissivity,                           &
                                rttovConfig % pccomp)

    call cpu_time(driver_time(5))
    call cosp_pc_rttov_save_output(rttovIN % nPoints,                       &
                                rttovConfig % nchannels_rec,             &
                                rttovConfig % swath_mask,                &
                                rttovConfig % pccomp,                    &
                                rttovConfig % Lrttov_bt,                 &
                                rttovConfig % Lrttov_rad,                &
                                bt_clear,                                &
                                rad_clear)

    call cpu_time(driver_time(6))                                           
    call cosp_pc_rttov_deallocate_profiles(rttovConfig % nprof,           &
                                        rttovConfig % nchanprof,       &
                                        rttovIN % nlevels,             &
                                        rttovConfig % nChannels_rec,   &
                                        rttovConfig % opts,            &
                                        profiles,                      &
                                        rttovConfig % coefs,           &
                                        chanprof,                      &
                                        transmission,                  &
                                        radiance,                      &
                                        calcemis,                      &
                                        emissivity,                    &
                                        rttovConfig % pccomp,          &
                                        predictindex)
                                        
    call cpu_time(driver_time(7))

  END SUBROUTINE COSP_PC_RTTOV_SIMULATE

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! END MODULE
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE MOD_COSP_RTTOV_INTERFACE
