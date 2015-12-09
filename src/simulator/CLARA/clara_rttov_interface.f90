! ########################################################################################
! This module contains subroutines used for calling rttov routines
! Code is strongly based on Ronald Scheirers and Thomas Landelius
! adaptions of the example_fwd.F90 from the rttov test routines
!
! $Id: rttov11_wrapper.F90 62 2015-10-19 19:27:06Z seliasson $
! Salomon Eliasson
! ########################################################################################

MODULE clara_rttov_interface

  use parkind1,       only : jpim,jprb,jplm
  use rttov_types,    only : rttov_options,rttov_coefs,profile_Type,transmission_Type,   &
                             radiance_Type,rttov_coef_scatt_ir,rttov_optpar_ir,          &
                             rttov_chanprof,rttov_emissivity,rttov_reflectance
  use rttov_const,    only : errorstatus_success,errorstatus_fatal,platform_name,inst_name
  use rttov_unix_env, only : rttov_exit
  use cosp_kinds,     only : wp
  IMPLICIT NONE
  
  ! Public subroutines
 ! public :: get_Tb_rttov_clear,get_Tb_rttov_black_clouds,get_Tb_rttov_clouds

  ! RTTOV11.3 specific parameters
  real(wp), parameter :: &
     tmin_baran = 193.1571_wp  ! Minimum temperature for baran parameterization in RTTOV  

  ! ######################################################################################
  ! Include files
  ! ######################################################################################
#include "rttov_alloc_rad.interface"
#include "rttov_alloc_transmission.interface"
#include "rttov_alloc_prof.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_direct.interface"
#include "rttov_print_opts.interface"
#include "rttov_print_profile.interface"
#include "rttov_read_coefs.interface"
#include "rttov_user_options_checkinput.interface"
#include "rttov_user_profile_checkinput.interface"

  ! ######################################################################################
  ! TYPE rttov_id
  ! ######################################################################################
  type rttov_id
     ! Satellite and instument id's:  name
     ! platform      : 1 = NOAA
     ! satellite     : 14 = number 14 in platform series
     ! sensor        : 5 = AVHRR
     ! nchannels     : 3 = number of channels 
     ! channels      : 1-3 = the IR channels on AVHRR
     ! coef_filename : name of the sensor coefficient file
     integer(Kind=jpim) :: &
        platform,  &
        satellite, &
        sensor,    &
        nchannels
     integer(Kind=jpim),allocatable,dimension(:) :: &
        channels
     character(len=256) :: &
        coef_filename,  & 
        sccld_filename, &
        scaer_filename, &
        rttov_dir
  end type rttov_id
  
  ! ######################################################################################
  ! TYPE rttov_other
  ! ######################################################################################
  type rttov_other
     logical(kind=jplm),allocatable,dimension(:) :: &
        calcemis, & ! Whether or not to provide your own emissivity
        calcrefl    ! Whether or not to provide your own reflectivity
     integer :: &
        dbg ! debug level 
  end type rttov_other
  
  ! ######################################################################################
  ! TYPE rttov_input
  ! ######################################################################################
  ! Input data to RTTOV
  type rttov_input
     ! Coefficients structure for cloud free, and cloudy conditions
     ! (simple or advanced cloud options are possible)
     type(rttov_coefs)                   :: coefs  
     ! profile from input model
     type(profile_type), allocatable     :: profiles(:)
     ! Input/output surface emissivity
     type(rttov_emissivity), allocatable :: emissivity(:)
     ! Input channel/profile list
     type(rttov_chanprof),   allocatable :: chanprof(:)    
     ! Flag to indicate calculation of emissivity (not the same length as in other_opts)
     logical(kind=jplm),     allocatable :: calcemis(:)       
     ! Flag to indicate calculation of BRDF
     logical(kind=jplm),     allocatable :: calcrefl(:)    
     ! some internal constants
     integer(kind=jpim)                  :: nchanprof ! nchannels x
     ! nprofiles
     integer(kind=jpim),     allocatable :: nchan(:)  ! Number of
     ! channels per
     ! profile
  end type rttov_input
  
  ! ######################################################################################
  ! TYPE rttov_output
  ! ######################################################################################  
  ! Output from RTTOV
  type rttov_output
     ! structure that contains the computed brightness temperature and
     type(radiance_type)                  :: radiance       
     type(transmission_type)              :: transmission  
     type(rttov_reflectance), allocatable :: reflectance(:)
  end type rttov_output

  ! ######################################################################################
  ! TYPE rttov
  ! ######################################################################################
  ! My main RTTOV structure 
  type rttovALL
     type(rttov_id)       :: id
     type(rttov_options)  :: opts   
     type(rttov_input)    :: input
     type(rttov_output)   :: output
     type(rttov_other)    :: other_opts
  end type rttovALL

  ! RTTOV setup for CLARA brightness temperature retrieval. Setup during initialization
  type(rttovALL),save :: R                        

contains

  ! ######################################################################################
  ! SUBROUTINE get_Tb_rttov_clear
  ! ######################################################################################
  subroutine get_Tb_rttov_clear(nprofIN,nlvlIN,nchanIN,pIN,tIN,qIN,t2mIN,p2mIN,sktIN,  &
                                waterTypeIN,zenangleIN,emisIN,surfTypeIN,TbOUT)
     ! Inputs
     integer,intent(in) :: &
        nprofIN,    & ! Number of points
        nlvlIN,     & ! Number of vertical levels
        nchanIN       ! Number of channels                                   
     real(wp),dimension(nprofIN,nlvlIN),intent(in) :: &
        pIN,        & ! Pressure (hPa)
        tIN,        & ! Temperature
        qIN           ! Specific humidity                           
     real(wp),dimension(nprofIN),intent(in) :: &
        t2mIN,      & ! 2-meter temperature (K)
        p2mIN,      & ! 2-meter pressure (hPa)
        sktIN,      & ! Skin Temperature (K)
        zenangleIN    ! Satellite zenith angle (degrees)
     real(wp),dimension(nprofIN*nchanIN),intent(in) :: &
        emisIN        ! Surface emissivity for all channels   
     integer,dimension(nprofIN),intent(in) :: &
        surfTypeIN, & ! Surface type
        waterTypeIN   ! Water type
    ! Outputs
    real(wp),dimension(nprofIN,nchanIN),intent(out) :: &
       TbOUT          ! RTTOV retrieved clear-sky brightness temperature

    ! Local variables
    real(kind=jprb),dimension(nprofIN,nlvlIN)  :: P,T,Q  
    real(kind=jprb),dimension(nprofIN)         :: T2M,P2M,SKT,zenangle
    real(kind=jprb),dimension(nprofIN*nchanIN) :: emis       
    integer(kind=jpim),dimension(nprofIN)      :: surftype,watertype         
    integer(kind=jpim)                         :: i,j,k,nprof,nlvl,nchan
    !real(kind=jprb),dimension(nprof,nchan)     :: Tb

    ! Convert inputs to RTTOV types
    nchan               = nchanIN
    nlvl                = nlvlIN
    nprof               = nprofIN
    emis(1:nprof*nchan) = emisIN(1:nprof*nchan)    
    watertype(1:nprof)  = watertypeIN(1:nprof)
    surftype(1:nprof)   = surftypeIN(1:nprof)
    t2m(1:nprof)        = t2mIN(1:nprof)
    p2m(1:nprof)        = p2mIN(1:nprof)
    skt(1:nprof)        = sktIN(1:nprof)
    zenangle(1:nprof)   = zenangleIN(1:nprof)
    P(1:nprof,1:nlvl)   = pIN(1:nprof,1:nlvl)

    ! Convert specific humidity to ppmv
    Q(1:nprof,1:nlvl)   = qIN(1:nprof,1:nlvl)*1.60771704e+6
    
    ! Set temperatures below some minimum to minimum.
    T(1:nprof,1:nlvl) = merge(tmin_baran,tIN(1:nProf,1:nlvl),tIN(1:nprof,1:nlvl) .lt. tmin_baran)         

    ! set options for clear sky
    R%opts%rt_ir%addclouds = .FALSE.

    ! Allocate space
    call allocate_rttov(nprof,nlvl)
  
    ! Initialize
    call initialise_rttov_mandatory(nprof,nlvl,nchan,P=P,T=T,Q=Q,T2M=T2M,P2M=P2M,      &
                                    SKT=SKT,emis=emis,surfType=surfType,                 &
                                    watertype=waterType,zenangle=zenangle)
    ! Call RTTOV
    call run_rttov() 

    ! DS: This can be done without looping. You can use the "reshape" intrinsic function.
    ! Get the brightness temperatures.
    k = 0 
    do i = 1, nprof
       do j = 1, nchan
          k = k+1
          TbOUT(i,j) = R % output % radiance % bt(k)
       end do
    end do
    
    ! Free up space
    call deallocate_rttov(nprof, nlvl)

  END SUBROUTINE get_Tb_rttov_clear
  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE rttov_simple
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  SUBROUTINE rttov_simple(nlevIN,nchanIN,tauIN,P_int,T_int,pIN,tIN,qIN,T2Min,P2Min,      &
                          SKTin,emisIN,satzenIN,minTau,waterTypeIN,sfcTypeIN,tauDepth,   &
                          TbOUT)

    ! Inputs
    integer,intent(in) :: &
       nlevIN,    & ! Number of vertical levels
       nchanIN,   & ! Number of RTTOV channels
       sfcTypeIN, & ! Surface type
       waterTypeIN  ! Water type
    real(wp),intent(in),dimension(nchanIN) :: &
       emisIN       ! Surface emissivity for RTTOV channels
    real(wp),intent(in),dimension(nlevIN) :: &
       tauIN,     & ! Optical-depth
       pIN,       & ! Pressure profile (hPa)
       tIN,       & ! Temperature profile (K)
       qIN          ! Specific-humidity profile (ppmv)
    real(wp),intent(in),dimension(nlevIN+1) :: &
       P_int,     & ! Pressure at model half-levels (hPa)
       T_int        ! Temperature at model half levels (K)
    real(wp),intent(in) :: &
       T2Min,     & ! Temperature at 2-meters (K)
       P2Min,     & ! Pressure at 2-meters (hPa)
       SKTin,     & ! Skin temperature (K)
       satzenIN,  & ! Satellite zenith angle
       minTau,    & ! Minimum detectable optical thickness
       tauDepth     ! Optical depth visible into the cloud
       
    ! Outputs
    real(wp),dimension(nchanIN),intent(inout) :: &
       TbOUT        ! Brightness temperature (K)

    ! Local variables
    integer(kind=jpim)                 :: nlev,nchan,sfcType(1),waterType(1)
    real(kind=jprb),dimension(nchanIN) :: emis,Tb
    real(kind=jprb),dimension(nlevIN)  :: P,T,Q
    real(kind=jprb)                    :: T2M(1),P2M(1),SKT(1),satzen(1)

    ! Parameters    
    REAL(kind=jprb), PARAMETER :: &
       ctpmin = 50.0_jprb,  & !
       ctpmax = 1100.0_jprb   !
    REAL(kind=jprb) :: tmp_ctp(1), cfraction(1)
    REAL :: ctp(1), ctt(1)

    ! Copy fields to RTTOV types
    nchan  = nchanIN
    nlev   = nlevIN
    P      = pIN
    T      = tIN
    q      = qIN*1.60771704e+6
    T2M    = T2Min
    P2M    = P2Min
    SKT    = SKTin
    emis   = emisIN
    satzen = satzenIN
    sfcType = sfcTypeIN
    waterType = waterTypeIN

    ! Initialize
    cfraction = 0._jprb
    ctp       = 0._jprb
    tmp_ctp   = 0._jprb

    ! First get the ctp
    CALL simple_ctth(nlevIN,tau=tauIN,P=P_int,T=T_int,tau_min=minTau,tau_eqCT = tauDepth,&
                     ctp=ctp,ctt=ctt)
    
    ! Bound the cloud-top pressure
    tmp_ctp = ctp
    tmp_ctp = maxval([tmp_ctp,ctpmin])
    tmp_ctp = minval([tmp_ctp,ctpmax])


    ! Cloud emissivity is assumed to be equivalent to cloud fraction. note: tau_IR approx tau_VIS/2
    cfraction = 1 - exp(-tauIN(nlev)/2)

    CALL get_Tb_rttov_black_clouds(Tb,1,nlev,nchan,P=P,T=T,Q=Q,T2M=T2M,P2M=P2M,SKT=SKT,  &
                                   emis=emis,surfType = sfcType,waterType = waterType,   &
                                   zenangle=satzen,ctp = tmp_ctp,cfraction = cfraction)
    TbOUT = Tb
    
  END SUBROUTINE rttov_simple
  
  ! ######################################################################################
  ! SUBROUTINE simple_ctth
  ! Find the cloud top heights per grid box by finding the level where e.g., 
  ! tau=1 (Holz et. al., 2006) integrated from the TOP DOWN (troposphere -> surface).
  !
  ! Originally coded by Salomon Eliasson (Salomon.Eliasson@smhi.se)
  ! Modified by Dustin Swales for use in COSPv2.0.0 (dustin.swales@noaa.gov)
  ! ######################################################################################    
  subroutine simple_ctth(nlev, tau, P, T, tau_min, tau_eqCT, ctp, ctt)

    ! Inputs
    integer,intent(in) :: &
       nlev ! Number of vertical levels
    real(wp),intent(in) :: &
       tau_min, & ! Minimum detectable optical-depth
       tau_eqCT   ! Threshold optical depth into cloud 
    real(wp),intent(in),dimension(nlev) :: &
       tau        ! Optical depth profile at 0.67 microns
    real(wp),intent(in),dimension(nlev+1) :: &
       P,       & ! Pressure profile (hPa) 
       T          ! Temperature profile (K)

    ! Outputs
    real(wp), intent(out) :: &
       ctp(1), & ! Cloud-top pressure
       ctt(1)    ! Cloud-top temperature

    ! Local variables
    integer :: &
       l2,          & ! Model level interface above the model layer where the total 
                      ! integrated optical depth is greater than tau_eqCT
       l1,          & ! Below
       inl,         & !
       inl2           !
    real(wp) :: &
       tau2,        & !
       tau1,        & !
       tau_ct_ratio   ! Where in the layer the integrated optical-depth=tau_eqCT
    real(wp),dimension(nlev) :: tau_int
    real(wp) :: ctp_, ctt_

    ! Initialize 
    ctp_         = -999._wp
    ctt_         = -999._wp
    ctp          = -999._wp
    ctt          = -999._wp
    tau2         = 0._wp
    tau1         = 0._wp
    tau_int      = 0._wp
    tau_CT_ratio = 0._wp
    L2           = 0
    L1           = 0

    ! For each model layer tau is centred around the middle between two layer
    ! interfaces. L2 and L1 are at the layer interfaces
    ! Loop from the top of the atmosphere down.
    do inl = 1, nlev

       ! gather tau from the top down. copy it over

       if (tau(inl) .le. 0) cycle 
       tau_int(inl) = tau(inl)

       ! Check if the sum of the tau_int (so far) is higher than
       ! the boundary we want. If it is, find the model layers
       ! between which the integrated tau = tau_eqCT

       ! this level is the upper bounds of tau (tau(L2)<tau1<tau(L1))
       tau1 = sum( tau_int(1:inl) )

       if (tau1 .GE. tau_eqCT) then

          ! L1 = the layer boundary below where tau_eqCT is exceeded
          L1 = inl +1 

          ! Look BACK UP the profile to find the closest level for
          !  which the integrated optical depth is less than
          !  tau_eqCT. Once this is found, add +1 to inl so that we
          !  pick <this> model layers' LOWER INTERFACE

          do inl2 = L1-2,1,-1
             if (tau_int(inl2)>0) then
                !doing this once
                tau2 = tau1 - tau_int(L1-1)
                L2 = inl2+1
                exit
             end if
          end do ! end looking back

          ! If the first level from the top down has tau >tau_eq
          ! then set L2 to L1 -1 (enough cloud is contained in
          ! this layer alone)
          if (L2 .eq. 0) then
             tau2 = 0.
             L2=L1-1
          end if

          ! Find the CTH between these layers using linear
          ! interpolation of log(pressure) according to 
          ! where tau_eqCT falls compared to the layer boundaries

          if ( (tau2 .gt. tau_eqCT) .or. (tau2 .gt. tau1) ) then
             stop "This has to be false: (tau2 .gt. tau_eq)&
                  & .or. (tau2 .gt. tau1)"
          end if

          tau_CT_ratio = (tau_eqCT-tau2)/(tau1-tau2)

          ! The pressures and Temperatures are at the layer interfaces
          ctp_ = EXP( LOG(P(L2)) + tau_CT_ratio*( LOG( P(L1)/P(L2) ) ) )
          ctt_ = T(L2) + tau_CT_ratio*( T(L1)-T(L2) )

          ctp = ctp_
          ctt = ctt_

       end if ! tau_int > tau_eqCT

       if (L2 .gt. 0) then
          exit
       end if
    end do ! loop levels

    if ( (sum(tau_int) .ge. tau_min) .and. (sum(tau_int) .lt. tau_eqCT) ) then
       ! Then there is a cloud that is semi-transparent and
       ! therefore we'll retrieve as far down as there are
       ! clouds and call that the equivalent cloud top.

       do inl = nlev,1,-1
          ! find the bottom of the cloud by looping up from
          !  the surface

          if (tau_int(inl)>0) then
             ctp = P(inl+1)
             ctt = T(inl+1)
             tau1=sum(tau_int)
             tau2=0.0
             tau_CT_ratio = 1
             exit
          end if
       end do

    end if ! if cloud is transparent
  end subroutine simple_ctth
  
  ! ######################################################################################
  ! SUBROUTINE get_Tb_rttov_clouds
  ! Subroutine to calculate radiances taking into account cloud particle scattering. Using
  ! the the most advanced scheme available in RTTOV.  
  ! ######################################################################################
  SUBROUTINE get_Tb_rttov_clouds(nprofIN,nlvlIN,nchanIN,pIN,tIN,qIN,T2Min,P2Min,SKTin,emisIN, &
                                 zenangleIN,waterTypeIN,sfcTypeIN,cloudIN,cfracIN,TbOUT)
                                 
    ! Inputs
    integer,intent(in) :: &
       nprofIN,   & ! Number of points
       nlvlIN,    & ! Number of vertical levels
       nchanIN,   & ! Number of RTTOV channels
       sfcTypeIN, & ! Surface type
       waterTypeIN  ! Water type
    real(wp),intent(in),dimension(nchanIN) :: &
       emisIN       ! Surface emissivity for RTTOV channels
    real(wp),intent(in),dimension(nlvlIN) :: &
       pIN,       & ! Pressure profile (hPa)
       tIN,       & ! Temperature profile (K)
       qIN          ! Specific-humidity profile (ppmv)
    real(wp),intent(in) :: &
       T2Min,     & ! Temperature at 2-meters (K)
       P2Min,     & ! Pressure at 2-meters (hPa)
       SKTin,     & ! Skin temperature (K)
       zenangleIN   ! Satellite zenith angle
    real(wp),intent(in),dimension(nprofIN,6,nlvlIN-1) :: &
       cloudIN      ! Cloud type (required by RTTOV)
    real(wp),intent(in),dimension(nprofIN,nlvlIN-1) :: &
       cfracIN      ! Cloud fraction  
    
    ! Outputs
    real(wp),intent(out),dimension(nprofIN,nchanIN) :: TbOUT                                 

    ! Local variables
    integer(kind=jpim)                            :: nprof, nlvl, nchan,i,j,k
    real(kind=jprb),dimension(nlvlIN)     :: P,T,Q   
    real(kind=jprb),dimension(1)            :: T2M,P2M,SKT,zenangle
    real(kind=jprb),dimension(nprofIN*nchanIN)    :: emis     
    integer(kind=jpim),dimension(1)         :: surftype,watertype
    real(kind=jprb),dimension(nprofIN,nchanIN)    :: Tb                              
    real(kind=jprb),dimension(nprofIN,6,nlvlIN-1) :: cloud
    real(kind=jprb),dimension(nprofIN,nlvlIN-1)   :: cfrac    

    ! Copy fields to RTTOV types
    nprof     = nprofIN
    nchan     = nchanIN
    nlvl      = nlvlIN
    P         = pIN
    T(1:nlvl) = merge(tmin_baran,tIN(1:nlvl),tIN(1:nlvl) .lt. tmin_baran)     
    q         = qIN*1.60771704e+6
    T2M       = T2Min
    P2M       = P2Min
    SKT       = SKTin
    emis      = emisIN
    zenangle  = zenangleIN
    surftype  = sfcTypeIN
    waterType = waterTypeIN
    cloud     = cloudIN
    cfrac     = cfracIN
	
    ! Initialize
    i=0; j=0; k=0; 
    
    ! Set options for semi-transparent clouds
    R % opts % rt_ir % addclouds = .TRUE.    
    
    ! ALLOCATE
    call allocate_rttov(nprof,nlvl)
    R % input % profiles % ish = 3 

    ! Initialize
    call initialise_rttov_mandatory(nprof,nlvl,nchan,P=P,T=T,Q=Q,T2M=T2M,P2M=P2M,      &
                                    SKT=SKT,emis=emis,surfType=surfType,                 &
                                    watertype=waterType,zenangle=zenangle)

    R % input % profiles(1) % cloud(:,:) = cloud(1,:,:)
    R % input % profiles(1) % cfrac(:) = cfrac(1,:)

    call run_rttov() 

    ! Get the brightness temperatures I need
    k = 0 
    do i = 1, nprof
       do j = 1, nchan
          k = k+1
          Tb(i,j) = R % output % radiance % bt(k)
       end do
    end do
    TbOut(1:nprof,1:nchan)=Tb(1:nprof,1:nchan)

    ! Deallocate
    call deallocate_rttov(nprof, nlvl)

  END SUBROUTINE get_Tb_rttov_clouds

  ! ######################################################################################
  ! SUBROUTINE get_TB_rttov_black_clouds
  ! Subroutine to calculate radiances treating cloud as a black body  
  ! ######################################################################################
  SUBROUTINE get_Tb_rttov_black_clouds(Tb,nprof,nlvl,nchan,P,T,Q,T2M,P2M,SKT,emis,     &
                                       surfType,watertype,zenangle, ctp, cfraction)

    integer(kind=jpim), intent(in) :: nprof, nlvl, nchan
    real(kind=jprb), intent(out) :: Tb (nprof, nchan)
    integer(kind=jpim) :: i,j,k,errorstatus

    real(kind=jprb), intent(in)    :: P        (nprof,nlvl)
    real(kind=jprb), intent(in)    :: T        (nprof,nlvl)
    real(kind=jprb), intent(in)    :: Q        (nprof,nlvl)
    real(kind=jprb), intent(in)    :: T2M      (nprof)
    real(kind=jprb), intent(in)    :: P2M      (nprof)
    real(kind=jprb), intent(in)    :: SKT      (nprof)
    real(kind=jprb), intent(in)    :: emis     (nprof * nchan)
    integer(kind=jpim), intent(in) :: surftype (nprof)
    integer(kind=jpim), intent(in) :: watertype(nprof)
    real(kind=jprb), intent(in)    :: zenangle (nprof)

    ! clouds
    real(kind=jprb), intent(in) :: ctp        (nprof)
    real(kind=jprb), intent(in) :: cfraction  (nprof)

    ! local
    k=0; i=0; j=0 

    ! set options for opaque clouds    
    R % opts % rt_ir % addclouds = .FALSE.

    ! ALLOCATE
    call allocate_rttov(nprof,nlvl)

    ! INITIALISE
    call initialise_rttov_mandatory(nprof, nlvl, nchan, &
            P=P,T=T,Q=Q,T2M=T2M,P2M=P2M,SKT=SKT,&
            emis=emis,surfType=surfType,watertype=waterType,&
            zenangle=zenangle)
    !do i=1,nlvl
    !   print*,i,q(1,i),(q(1,i)/(q(1,i)+0.622*(1-q(1,i))))*1e6
    !enddo
    !stop

    R % input % profiles(1) % ctp = ctp(1)
    R % input % profiles(1) % cfraction = cfraction(1)

    errorstatus = 0_jpim
    call run_rttov() 

    ! Get the brightness temperatures I need
    k = 0 
    do i = 1, nprof
       do j = 1, nchan
          k = k+1
          Tb(i,j) = R % output % radiance % bt(k)
       end do
    end do

    call deallocate_rttov(nprof, nlvl)
    
  END SUBROUTINE get_Tb_rttov_black_clouds

  ! ######################################################################################
  ! SUBROUTINE get_rttov_coeffs
  ! Read the coefficient table for AVHRR and set flags
  ! ######################################################################################
  SUBROUTINE get_rttov_coeffs (error)

    integer(kind=jpim), intent(out) :: error
    integer(kind=jpim)              :: errorstatus    ! Return error
    ! status of
    ! RTTOV
    ! subroutine
    ! calls
    character(len=256)              :: toAppend
    character(len=256)              :: rttov_dir
    logical                         :: file_exists
    logical                         :: clouds, aerosols

    ! ------- End header ----

    errorstatus     = 0_jpim
    clouds          = R % opts % rt_ir % addclouds
    aerosols        = R % opts % rt_ir % addaerosl

    ! ------------------
    ! PUT THE FILENAME TOGETHER
    ! 
    ! e.g. $SOME_PATH/rttov_11.1/rtcoef_rttov11/rttov7pred54L/rtcoef_noaa_14_avhrr.dat

    ! Fixme: It could be that RTTOV actually has a function that does
    ! what this subroutine is doing
    call getenv("RTTOV_DIR",rttov_dir)
    rttov_dir = R%id%rttov_dir
    
    ! Platform
    if (R % id % platform .eq. 1) then
       toAppend = 'noaa_'
    else if (R % id % platform .eq. 10) then
       toAppend = 'metop_'
    else if (R % id % platform .eq. 11) then
       toAppend = 'envisat_'
    else
       write ( *,* ) 'Unsupported platform ID ', R % id % platform
       error = 1
       stop
    endif

    ! Satellite
    if (R % id % satellite .lt. 10) then
       toAppend = &
            trim(toAppend) // char(R % id % satellite+48)
    else if (R % id % satellite .lt. 100) then
       toAppend = &
            trim(toAppend) // char(int(R % id % satellite/10)+48)
       toAppend = &
            trim(toAppend) // char(R % id % satellite-int(R % id % satellite/10)*10+48)
    else
       write ( *,* ) 'Unsupported satellite number ', R % id % satellite
       error = 1
       return
    endif

    ! Sensor
    if (R % id % sensor .eq. 3) then
       toAppend = trim(toAppend) // '_amsua.dat'
    elseif (R % id % sensor .eq. 5) then
       toAppend = trim(toAppend) // '_avhrr.dat'
    elseif (R % id % sensor .eq. 49) then
       toAppend = trim(toAppend) // '_mwr.dat'
    else
       write ( *,* ) 'Unsupported sensor number ', R % id % sensor
       error = 1
       return
    endif

    ! coef filename
    write(R % id % coef_filename,'(A,"/rtcoef_rttov11/rttov7pred54L/rtcoef_",A)') &
         trim(rttov_dir),trim(toAppend)
    inquire(FILE=trim(R % id % coef_filename), EXIST=file_exists)
    if ( .not. file_exists ) then
       write (*,*) "File does not exist", trim(R % id % coef_filename)
       error = 1
    else
       write (*,*) "RTTOV coefficients from: ", trim(toAppend)
    end if

    if (clouds) then
       ! cloud scattering coefs. Allways read this!
       write(R % id % sccld_filename,'(A,"/rtcoef_rttov11/cldaer/sccldcoef_",A)') &
            trim(rttov_dir),trim(toAppend)
       inquire(FILE=trim(R % id % sccld_filename), EXIST=file_exists)
       if ( .not. file_exists ) then
          print '("File:",1x,A,", ","does not exist!")',&
               trim(R % id % sccld_filename)
          stop
       end if
    end if

    if (aerosols) then
       write(R % id % scaer_filename,'(A,"/rtcoef_rttov11/cldaer/scaercoef_",A)') &
            trim(rttov_dir),trim(toAppend)
       inquire(FILE=trim(R % id % scaer_filename), EXIST=file_exists)
       if ( .not. file_exists ) then
          print '("File:",1x,A,", ","does not exist!")',&
               trim(R % id % scaer_filename)
          stop
       end if
    end if
    

    ! no clouds or scattering
    if (.not. aerosols .and. .not. clouds) then
       call rttov_read_coefs(errorstatus, &
            coefs=R % input % coefs, opts= R % opts, &
            form_coef='formatted', &
            file_coef=R % id % coef_filename)

    elseif (.not. aerosols .and. clouds ) then
       ! cloud scattering but no aerosols
       call rttov_read_coefs(errorstatus, &
            coefs=R % input % coefs, opts= R % opts, &
            form_coef='formatted', form_sccld = 'formatted', &
            file_coef=R % id % coef_filename, &
            file_sccld= R % id % sccld_filename)

    elseif (aerosols .and. .not. clouds ) then
       ! no cloud scattering but aerosols
       call rttov_read_coefs(errorstatus, &
            coefs=R % input % coefs, opts= R % opts, &
            form_coef='formatted', form_scaer='formatted', &
            file_coef=R % id % coef_filename, &
            file_scaer= R % id % scaer_filename)

    elseif (aerosols .and. clouds) then
       ! cloud scattering and cloud aeorosols
       call rttov_read_coefs(errorstatus, &
            coefs=R % input % coefs, opts= R % opts, &
            form_coef='formatted', form_sccld = 'formatted',form_scaer='formatted', &
            file_coef=R % id % coef_filename, &
            file_sccld= R % id % sccld_filename, &
            file_scaer= R % id % scaer_filename)
    end if

    if (errorstatus /= errorstatus_success) then
       write(*,*) 'fatal error reading coefficients'
       call rttov_exit(errorstatus)
    endif

    ! Ensure the options and coefficients are consistent
    !  call rttov_user_options_checkinput(errorstatus, opts, coefs)
    !  if (errorstatus /= errorstatus_success) then
    !    write(*,*) 'error in rttov options'
    !    call rttov_exit(errorstatus)
    !  endif

  END SUBROUTINE get_rttov_coeffs
  
  ! ######################################################################################
  ! SUBROUTINE allocate_rttov
  ! Allocate profiles, radiance, transmittance, reflectance and emissivity structures
  ! ######################################################################################
  SUBROUTINE allocate_rttov(nprof,nlevels)

    integer(kind=jpim), intent(in)  :: nprof
    integer(kind=jpim), intent(in)  :: nlevels
    integer(kind=jpim)              :: nlayers
    integer(kind=jpim)              :: errorstatus
    integer(kind=jpim)              :: alloc_status(10)

    errorstatus     = 0_jpim
    alloc_status(:) = 0_jpim
    nlayers = nlevels-1

    ! PROFILES
    allocate(R % input % profiles(nprof), stat=alloc_status(1))
    if (ANY(alloc_status /= 0)) then
       write(*,*) 'mem allocation error for profile array'
       call rttov_exit(errorstatus_fatal)
    endif

    call rttov_alloc_prof(                 &
         & errorstatus,                    &
         & nprof,                          &
         & profiles = R % input % profiles,&
         & nlevels = nlevels,              &
         & opts = R % opts,                &
         & asw = 1_jpim,                   & ! =allocate
         & coefs=R % input % coefs,        &
         & init=.TRUE._jplm                 )
    if (errorstatus /= errorstatus_success) then
       write(*,*) 'allocation error for profile arrays'
       call rttov_exit(errorstatus)
    endif

    ! RADIANCES
    allocate(R % input % nchan(nprof), stat=alloc_status(1))
    if (ANY(alloc_status /= 0)) then
       write(*,*) 'mem allocation error for nchan'
       call rttov_exit(errorstatus_fatal)
    endif

    R % input % nchan(:) = R % id % nchannels
    R % input % nchanprof = SUM(R % input % nchan(:))  ! Size of chanprof array
    ! is total number of
    ! channels over all
    ! profiles

    call rttov_alloc_rad(        &
         & errorstatus,          &
         & R % input % nchanprof,&
         & R % output % radiance,&
         & nlayers,              &
         & asw = 1_jpim,         &
         init=.TRUE._jplm)
    if (errorstatus /= errorstatus_success) then
       write(*,*) 'allocation error for radiance arrays'
       call rttov_exit(errorstatus)
    endif

    ! TRANSMITANCE    
    call rttov_alloc_transmission(            &
         & errorstatus,                       &
         & R % output % transmission,         &
         & nlayers,                           &
         & R % input % nchanprof,             &
         & asw = 1_jpim,                      &
         & init=.TRUE._jplm                    )

    if (errorstatus /= errorstatus_success) then
       write(*,*) 'allocation error for transmission arrays'
       call rttov_exit(errorstatus)
    endif

    ! EMISSIVITY
    ALLOCATE(R % input % calcemis(R % input % nchanprof), stat=alloc_status(1))
    ALLOCATE(R % input % emissivity(R % input % nchanprof), stat=alloc_status(2))
    IF (ANY(alloc_status /= 0)) THEN
       WRITE(*,*) 'mem allocation error for emissivity arrays'
       CALL rttov_exit(errorstatus_fatal)
    ENDIF

    ! REFLECTANCE
    allocate(R % input % calcrefl(R % input % nchanprof), stat=alloc_status(1))
    allocate(R % output % reflectance(R % input % nchanprof), stat=alloc_status(2))
    if (ANY(alloc_status /= 0)) then
       write(*,*) 'memallocation error for reflectance arrays'
       call rttov_exit(errorstatus_fatal)
    endif


    ! CHANPROF
    allocate ( R % input % chanprof(R % input % nchanprof), stat=alloc_status(1) )
    if (ANY(alloc_status /= 0)) then
       write(*,*) 'memallocation error for chanprof'
       call rttov_exit(errorstatus_fatal)
    endif

  END SUBROUTINE allocate_rttov

  ! ######################################################################################
  ! SUBROUTINE initialise_rttov_mandatory
  ! Initialize memory for RTTOV
  ! ######################################################################################
  SUBROUTINE initialise_rttov_mandatory(nprof,nlvl,nchan,P,T,Q,T2M,P2M,SKT,emis,       &
                                        surfType,watertype,zenangle)

    integer(kind=jpim), intent(in) :: nprof, nlvl, nchan
    integer(kind=jpim) :: jch, inp, nch
    real(kind=jprb), intent(in) :: P                    (nprof,nlvl)
    real(kind=jprb), intent(in) :: T                    (nprof,nlvl)
    real(kind=jprb), intent(in) :: Q                    (nprof,nlvl)
    real(kind=jprb), intent(in) :: T2M                  (nprof)
    real(kind=jprb), intent(in) :: P2M                  (nprof)
    real(kind=jprb), intent(in) :: SKT                  (nprof)
    real(kind=jprb), intent(in) :: emis                 (nprof * nchan)
    integer(kind=jpim), intent(in) :: surftype          (nprof)
    integer(kind=jpim), intent(in) :: watertype         (nprof)
    real(kind=jprb), intent(in) :: zenangle             (nprof)

! debug
    integer :: ii

    nch = 0_jpim
    do inp = 1,nprof
       do jch = 1, nchan
          ! CHANPROF, EMISSIVITY and REFLECTIVITY
          nch = nch + 1_jpim
          R % input % calcemis(nch) = R % other_opts % calcemis(jch)
          R % input % calcrefl(nch) = R % other_opts % calcrefl(jch)
          R % input % chanprof(nch) % prof = inp
          R % input % chanprof(nch) % chan = R % id % channels(jch)
       end do

       ! ---------
       !  Feed profiles with values
       ! ---------

       R % input % profiles(inp) % p   (:) = P(inp,:)
       R % input % profiles(inp) % t   (:) = T(inp,:)
       R % input % profiles(inp) % q   (:) = Q(inp,:)

       !2m variables
       R % input % profiles(inp) % s2m % t = T2M(inp)
       R % input % profiles(inp) % s2m % p = P2M(inp)
       R % input % profiles(inp) % skin % t         = SKT(inp)
       R % input % profiles(inp) % skin % surftype  = surfType(inp)
       R % input % profiles(inp) % skin % watertype = watertype(inp)
       R % input % profiles(inp) % zenangle = zenangle(inp)

       ! This has to always be in place
       R % input % profiles(inp) % ctp       = 500._jprb
       R % input % profiles(inp) % cfraction = 0._jprb
    end do ! nprof

    ! Emissivity
    do ii = 1, nprof * nchan
       IF (.not. R % input % calcemis(ii) ) then
          R % input % emissivity(ii) % emis_in = emis(ii) ! length = nchanprof
       ELSE
          R % input % emissivity(ii) % emis_in = 0._jprb
       END IF
    end do

  END SUBROUTINE initialise_rttov_mandatory

  ! ######################################################################################
  ! SUBROUTINE run_rttov
  ! ######################################################################################
  SUBROUTINE run_rttov ()

    ! internal
    integer(kind=jpim)  :: errorstatus
    type(rttov_input)   :: in
    type(rttov_options) :: opts
    type(rttov_output)  :: out
    errorstatus = 0_jpim
    in          = R % input
    out         = R % output
    opts        = R % opts
    ! debug


    if (R % other_opts % dbg > 0) then
       ! Ensure the options and coefficients are consistent

       call rttov_print_opts(opts)
       
       call rttov_user_options_checkinput(errorstatus, opts, in % coefs)
       
       if (errorstatus /= errorstatus_success) then
          write(*,*) 'error in rttov options'
          !call rttov_exit(errorstatus)
       endif
       
       call rttov_print_profile(in % profiles(1))
       ! check profile
       call rttov_user_profile_checkinput(errorstatus, opts, in % coefs, &
            in % profiles(1)) ! need (1) otherwise segfault
       
       if (errorstatus /= errorstatus_success) then
          write(*,*) 'error in rttov profile'
          !call rttov_exit(errorstatus)
       endif
       
    end if

    call rttov_direct(                        &
         & errorstatus,                       &! out   error flag
         & in  % chanprof,                    &! in    channel and profile index structure
         & R   % opts,                        &! in    options structure
         & in  % profiles,                    &! in    profile array
         & in  % coefs,                       &! in    coefficients strucutre
         & out % transmission,                &! inout computed transmittances
         & out % radiance,                    &! inout computed radiances
         & calcemis     = in  % calcemis,     &! in    flag for internal emissivity calcs
         & emissivity   = in  % emissivity,   &! inout input/output emissivities per channel
         & calcrefl     = in  % calcrefl,     &! in    flag for internal BRDF calcs
         & reflectance  = out % reflectance)   ! inout input/output BRDFs per channel      

  IF (errorstatus /= errorstatus_success) THEN
    WRITE (*,*) 'rttov_direct error'
    CALL rttov_exit(errorstatus)
  ENDIF

    ! Merge it in to my RTTOV structure
    R % output = out

    ! The radiance is found in               rad = R % output % radiance % total
    ! The brightness temperature is found in Tb  = R % output % radiance % bt
  END SUBROUTINE run_rttov
  
  ! ######################################################################################
  ! SUBROUTINE deallocate_rttov
  ! ######################################################################################
  SUBROUTINE deallocate_rttov(nprof,nlevels)

    integer(kind=jpim), intent(in)        :: nprof, nlevels
    integer(kind=jpim)                    :: nlayers
    integer(kind=jpim) :: alloc_status(6)
    integer(kind=jpim) :: errorstatus
    alloc_status(:) = 0_jpim
    errorstatus = 0_jpim

    nlayers = nlevels -1

    deallocate (R % input % nchan,        stat=alloc_status(1))
    deallocate (R % input % chanprof,     stat=alloc_status(2))
    deallocate (R % input % emissivity,   stat=alloc_status(3))
    deallocate (R % input % calcemis,     stat=alloc_status(4))
    deallocate (R % input % calcrefl,     stat=alloc_status(5))
    deallocate (R % output % reflectance, stat=alloc_status(6))
    !    deallocate (R % id % channels,        stat=alloc_status(7))
    if (any(alloc_status /= 0)) then
       write(*,*) 'mem dellocation error'
    endif

    ! deallocate radiance arrays
    call rttov_alloc_rad(errorstatus, &
         R % id % nchannels,          &
         R % output % radiance,       &
         nlevels-1_jpim,              &
         asw=0_jpim)
    if(errorstatus /= errorstatus_success) then
       write(*,*) 'radiance deallocation error'
    endif

    ! deallocate transmission arrays
    call rttov_alloc_transmission(errorstatus, &
         R % output % transmission,            &
         nlayers,                              &
         R % id % nchannels,                   &
         asw=0_jpim)
    if (errorstatus /= errorstatus_success) then
       write(*,*) 'transmission deallocation error'
    endif

    ! deallocate profile arrays
    ! first the innards
    call rttov_alloc_prof(errorstatus, &
         nprof,                        &
         R % input % profiles,         &
         nlevels,                      &
         R % opts,                     &
         asw=0_jpim)
    if (errorstatus /= errorstatus_success .or. alloc_status(1) /= 0) then
       write(*,*) 'profile deallocation error'
    endif
    ! then the structure
    deallocate(R % input % profiles, stat=alloc_status(1))
    if (alloc_status(1) /= 0) then
       write(*,*) 'mem deallocation error for profile array'
    endif

  END SUBROUTINE deallocate_rttov
  ! ######################################################################################
  ! END MODULE clara_rttov_interface
  ! ######################################################################################
END MODULE clara_rttov_interface
