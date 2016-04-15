module mod_cosp_io
  use cosp_kinds, only: wp
  use netcdf
  use mod_cosp,   only: cosp_outputs
  USE cmor_users_functions
  USE MOD_COSP_CONFIG, ONLY: vgrid_zu,vgrid_zl,vgrid_z,&!mgrid_zu,mgrid_zl,mgrid_z,&
                             R_UNDEF,DBZE_BINS,CFAD_ZE_MIN,CFAD_ZE_WIDTH,COSP_VERSION,   &
                             LIDAR_NTEMP,LIDAR_PHASE_TEMP,LIDAR_PHASE_TEMP_BNDS,misr_histHgtCenters,&
                             misr_histHgtEdges,numMISRHgtBins,PARASOL_NREFL,SR_BINS,   &
                             isccp_histTauEdges,isccp_histPresEdges,PARASOL_SZA,          &
                             isccp_histPresCenters,isccp_histTauCenters
  USE MOD_COSP_RTTOV, ONLY: iChannel

  implicit none
  
  ! Types to be used as arrays of pointers
  TYPE var1d
     character(len=16) :: name
     character(len=16) :: units
     integer :: dimsid(3)
     integer :: dimssz(2)
     integer :: vid
     logical :: lout
     real(wp),pointer,dimension(:) :: pntr
  END TYPE var1d
  TYPE var2d
     character(len=16) :: name
     character(len=16) :: units
     integer :: dimsid(4)
     integer :: dimssz(3)
     integer :: vid
     logical :: lout
     real(wp),pointer,dimension(:,:) :: pntr
  END TYPE var2d
  TYPE var3d
     character(len=16) :: name
     character(len=16) :: units
     integer :: dimsid(5)
     integer :: dimssz(4)
     integer :: vid
     logical :: lout
     real(wp),pointer,dimension(:,:,:) :: pntr
  END TYPE var3d
contains
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE construct_cospOutList
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine construct_cospOutList(Lpctisccp,Lclisccp,Lboxptopisccp,Lboxtauisccp,        &
                                   Ltauisccp,Lcltisccp,Lmeantbisccp,Lmeantbclrisccp,     &
                                   Lalbisccp,LclMISR,Lcltmodis,Lclwmodis,Lclimodis,      &
                                   Lclhmodis,Lclmmodis,Lcllmodis,Ltautmodis,Ltauwmodis,  &
                                   Ltauimodis,Ltautlogmodis,Ltauwlogmodis,Ltauilogmodis, &
                                   Lreffclwmodis,Lreffclimodis,Lpctmodis,Llwpmodis,      &
                                   Liwpmodis,Lclmodis,Latb532,LlidarBetaMol532,          &
                                   LcfadLidarsr532,Lclcalipso2,Lclcalipso,Lclhcalipso,   &
                                   Lcllcalipso,Lclmcalipso,Lcltcalipso,Lcltlidarradar,   &
                                   Lclcalipsoliq,Lclcalipsoice,Lclcalipsoun,Lclcalipsotmp,&
                                   Lclcalipsotmpliq,Lclcalipsotmpice,Lclcalipsotmpun,    &
                                   Lcltcalipsoliq,Lcltcalipsoice,Lcltcalipsoun,          &
                                   Lclhcalipsoliq,Lclhcalipsoice,Lclhcalipsoun,          &
                                   Lclmcalipsoliq,Lclmcalipsoice,Lclmcalipsoun,          &
                                   Lcllcalipsoliq,Lcllcalipsoice,Lcllcalipsoun,          & 
                                   LcfadDbze94,Ldbze94,Lparasolrefl,Ltbrttov,N_out_list, &
                                   out_list)
    ! Inputs
    logical,intent(in) :: &
         Lpctisccp,        & ! ISCCP mean cloud top pressure
         Lclisccp,         & ! ISCCP cloud area fraction
         Lboxptopisccp,    & ! ISCCP CTP in each column
         Lboxtauisccp,     & ! ISCCP optical epth in each column
         Ltauisccp,        & ! ISCCP mean optical depth
         Lcltisccp,        & ! ISCCP total cloud fraction
         Lmeantbisccp,     & ! ISCCP mean all-sky 10.5micron brightness temperature
         Lmeantbclrisccp,  & ! ISCCP mean clear-sky 10.5micron brightness temperature
         Lalbisccp,        & ! ISCCP mean cloud albedo         
         LclMISR,          & ! MISR cloud fraction
         Lcltmodis,        & ! MODIS total cloud fraction
         Lclwmodis,        & ! MODIS liquid cloud fraction
         Lclimodis,        & ! MODIS ice cloud fraction
         Lclhmodis,        & ! MODIS high-level cloud fraction
         Lclmmodis,        & ! MODIS mid-level cloud fraction
         Lcllmodis,        & ! MODIS low-level cloud fraction
         Ltautmodis,       & ! MODIS total cloud optical thicknes
         Ltauwmodis,       & ! MODIS liquid optical thickness
         Ltauimodis,       & ! MODIS ice optical thickness
         Ltautlogmodis,    & ! MODIS total cloud optical thickness (log10 mean)
         Ltauwlogmodis,    & ! MODIS liquid optical thickness (log10 mean)
         Ltauilogmodis,    & ! MODIS ice optical thickness (log10 mean)
         Lreffclwmodis,    & ! MODIS liquid cloud particle size
         Lreffclimodis,    & ! MODIS ice particle size
         Lpctmodis,        & ! MODIS cloud top pressure
         Llwpmodis,        & ! MODIS cloud liquid water path
         Liwpmodis,        & ! MODIS cloud ice water path
         Lclmodis,         & ! MODIS cloud area fraction
         Latb532,          & ! CALIPSO attenuated total backscatter (532nm)
         LlidarBetaMol532, & ! CALIPSO molecular backscatter (532nm)         
         LcfadLidarsr532,  & ! CALIPSO scattering ratio CFAD
         Lclcalipso2,      & ! CALIPSO cloud fraction undetected by cloudsat
         Lclcalipso,       & ! CALIPSO cloud area fraction
         Lclhcalipso,      & ! CALIPSO high-level cloud fraction
         Lcllcalipso,      & ! CALIPSO low-level cloud fraction
         Lclmcalipso,      & ! CALIPSO mid-level cloud fraction
         Lcltcalipso,      & ! CALIPSO total cloud fraction
         Lcltlidarradar,   & ! CALIPSO-CLOUDSAT total cloud fraction
         Lclcalipsoliq,    & ! CALIPSO liquid cloud area fraction
         Lclcalipsoice,    & ! CALIPSO ice cloud area fraction 
         Lclcalipsoun,     & ! CALIPSO undetected cloud area fraction
         Lclcalipsotmp,    & ! CALIPSO undetected cloud area fraction
         Lclcalipsotmpliq, & ! CALIPSO liquid cloud area fraction
         Lclcalipsotmpice, & ! CALIPSO ice cloud area fraction
         Lclcalipsotmpun,  & ! CALIPSO undetected cloud area fraction
         Lcltcalipsoliq,   & ! CALIPSO liquid total cloud fraction
         Lcltcalipsoice,   & ! CALIPSO ice total cloud fraction
         Lcltcalipsoun,    & ! CALIPSO undetected total cloud fraction
         Lclhcalipsoliq,   & ! CALIPSO high-level liquid cloud fraction
         Lclhcalipsoice,   & ! CALIPSO high-level ice cloud fraction
         Lclhcalipsoun,    & ! CALIPSO high-level undetected cloud fraction
         Lclmcalipsoliq,   & ! CALIPSO mid-level liquid cloud fraction
         Lclmcalipsoice,   & ! CALIPSO mid-level ice cloud fraction
         Lclmcalipsoun,    & ! CALIPSO mid-level undetected cloud fraction
         Lcllcalipsoliq,   & ! CALIPSO low-level liquid cloud fraction
         Lcllcalipsoice,   & ! CALIPSO low-level ice cloud fraction
         Lcllcalipsoun,    & ! CALIPSO low-level undetected cloud fraction
         LcfadDbze94,      & ! CLOUDSAT radar reflectivity CFAD
         Ldbze94,          & ! CLOUDSAT radar reflectivity
         LparasolRefl,     & ! PARASOL reflectance
         Ltbrttov            ! RTTOV mean clear-sky brightness temperature                                    
    integer, intent(in) :: n_out_list

    ! Outputs
    character(len=32),dimension(n_out_list),intent(out) :: out_list

    if (Lalbisccp)         out_list(1)  = 'albisccp'
    if (Latb532)           out_list(2)  = 'atb532'
    if (Lboxptopisccp)     out_list(3)  = 'boxptopisccp'
    if (Lboxtauisccp)      out_list(4)  = 'boxtauisccp'
    if (LcfadDbze94)       out_list(5)  = 'cfadDbze94'
    if (LcfadLidarsr532)   out_list(6)  = 'cfadLidarsr532'
    if (Lclcalipso2)       out_list(7)  = 'clcalipso2'
    if (Lclcalipso)        out_list(8)  = 'clcalipso'
    if (Lclhcalipso)       out_list(9)  = 'clhcalipso'
    if (Lclisccp)          out_list(10) = 'clisccp'
    if (Lcllcalipso)       out_list(11) = 'cllcalipso'
    if (Lclmcalipso)       out_list(12) = 'clmcalipso'
    if (Lcltcalipso)       out_list(13) = 'cltcalipso'
    if (Lcllcalipsoice)    out_list(14) = 'cllcalipsoice'
    if (Lclmcalipsoice)    out_list(15) = 'clmcalipsoice'
    if (Lclhcalipsoice)    out_list(16) = 'clhcalipsoice'
    if (Lcltcalipsoice)    out_list(17) = 'cltcalipsoice'
    if (Lcllcalipsoliq)    out_list(18) = 'cllcalipsoliq'
    if (Lclmcalipsoliq)    out_list(19) = 'clmcalipsoliq'
    if (Lclhcalipsoliq)    out_list(20) = 'clhcalipsoliq'
    if (Lcltcalipsoliq)    out_list(21) = 'cltcalipsoliq'
    if (Lcllcalipsoun)     out_list(22) = 'cllcalipsoun'
    if (Lclmcalipsoun)     out_list(23) = 'clmcalipsoun'
    if (Lclhcalipsoun)     out_list(24) = 'clhcalipsoun'
    if (Lcltcalipsoun)     out_list(25) = 'cltcalipsoun'
    if (Lclcalipsoice)     out_list(26) = 'clcalipsoice'
    if (Lclcalipsoliq)     out_list(27) = 'clcalipsoliq'
    if (Lclcalipsoun)      out_list(28) = 'clcalipsoun'
    if (Lclcalipsotmp)     out_list(29) = 'clcalipsotmp'
    if (Lclcalipsotmpice)  out_list(30) = 'clcalipsotmpice'
    if (Lclcalipsotmpliq)  out_list(31) = 'clcalipsotmpliq'
    if (Lclcalipsotmpun)   out_list(32) = 'clcalipsotmpun'
    if (Lcltlidarradar)    out_list(33) = 'cltlidarradar'
    if (Lpctisccp)         out_list(34) = 'pctisccp'
    if (Ldbze94)           out_list(35) = 'dbze94'
    if (Ltauisccp)         out_list(36) = 'tauisccp'
    if (Lcltisccp)         out_list(37) = 'cltisccp'
    if (LparasolRefl)      out_list(39) = 'parasolRefl'
    if (LclMISR)           out_list(40) = 'clMISR'
    if (Lmeantbisccp)      out_list(41) = 'meantbisccp'
    if (Lmeantbclrisccp)   out_list(42) = 'meantbclrisccp'
    if (LlidarBetaMol532)  out_list(44) = 'lidarBetaMol532'
    if (Ltbrttov)          out_list(45) = 'tbrttov'
    if (Lcltmodis)         out_list(46) = 'cltmodis'
    if (Lclwmodis)         out_list(47) = 'clwmodis'
    if (Lclimodis)         out_list(48) = 'climodis'
    if (Lclhmodis)         out_list(49) = 'clhmodis'
    if (Lclmmodis)         out_list(50) = 'clmmodis'
    if (Lcllmodis)         out_list(51) = 'cllmodis'
    if (Ltautmodis)        out_list(52) = 'tautmodis'
    if (Ltauwmodis)        out_list(53) = 'tauwmodis'
    if (Ltauimodis)        out_list(54) = 'tauimodis'
    if (Ltautlogmodis)     out_list(55) = 'tautlogmodis'
    if (Ltauwlogmodis)     out_list(56) = 'tauwlogmodis'
    if (Ltauilogmodis)     out_list(57) = 'tauilogmodis'
    if (Lreffclwmodis)     out_list(58) = 'reffclwmodis'
    if (Lreffclimodis)     out_list(59) = 'reffclimodis'
    if (Lpctmodis)         out_list(60) = 'pctmodis'
    if (Llwpmodis)         out_list(61) = 'lwpmodis'
    if (Liwpmodis)         out_list(62) = 'iwpmodis'
    if (Lclmodis)          out_list(63) = 'clmodis'

  end subroutine construct_cospOutList

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE nc_read_input_file
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE NC_READ_INPUT_FILE(fname,Npnts,Nl,Nhydro,lon,lat,p,ph,z,zh,T,qv,rh,tca,cca, &
                                mr_lsliq,mr_lsice,mr_ccliq,mr_ccice,fl_lsrain,fl_lssnow, &
                                fl_lsgrpl,fl_ccrain,fl_ccsnow,Reff,dtau_s,dtau_c,dem_s,  &
                                dem_c,skt,landmask,mr_ozone,u_wind,v_wind,sunlit,        &
                                emsfc_lw,mode,Nlon,Nlat,verbosity)
     
    ! Arguments
    character(len=512),intent(in) :: fname ! File name
    integer,intent(in) :: Npnts,Nl,Nhydro
    real(wp),dimension(Npnts),intent(out) :: lon,lat
    real(wp),dimension(Npnts,Nl),target,intent(out) :: p,ph,z,zh,T,qv,rh,tca,cca, &
         mr_lsliq,mr_lsice,mr_ccliq,mr_ccice,fl_lsrain,fl_lssnow,fl_lsgrpl, &
         fl_ccrain,fl_ccsnow,dtau_s,dtau_c,dem_s,dem_c,mr_ozone
    real(wp),dimension(Npnts,Nl,Nhydro),intent(out) :: Reff
    real(wp),dimension(Npnts),intent(out) :: skt,landmask,u_wind,v_wind,sunlit
    real(wp),intent(out) :: emsfc_lw
    integer,intent(out) :: mode,Nlon,Nlat
    integer,optional :: verbosity
    
    ! Local variables
    integer,parameter :: NMAX_DIM=5
    integer :: Npoints,Nlevels,i,j,k,vrank,vdimid(NMAX_DIM),ncid,vid,ndims,nvars,ngatts, &
               recdim,dimsize(NMAX_DIM),errst,Na,Nb,Nc,Nd,Ne
    integer,dimension(:),allocatable :: plon,plat
    logical :: Llat,Llon,Lpoint
    real(wp),dimension(Npnts) :: ll
    real(wp),allocatable :: x1(:),x2(:,:),x3(:,:,:),x4(:,:,:,:),x5(:,:,:,:,:) ! Temporary arrays
    character(len=128) :: vname
    character(len=256) :: dimname(NMAX_DIM) ! 256 hardcoded, instead of MAXNCNAM. This works for NetCDF 3 and 4.
    character(len=64) :: routine_name='NC_READ_INPUT_FILE'
    character(len=128) :: errmsg,straux
    
    mode = 0
    Nlon = 0
    Nlat = 0
    
    Npoints = Npnts
    Nlevels = Nl
    
    ! Open file
    errst = nf90_open(fname, nf90_nowrite, ncid)
    if (errst /= 0) then
       errmsg="Couldn't open "//trim(fname)
       call cosp_error(routine_name,errmsg)
    endif
    
    ! Get information about dimensions. Curtain mode or lat/lon mode?
    Llat  =.false.
    Llon  =.false.
    Lpoint=.false.
    errst = nf90_inquire(ncid, ndims, nvars, ngatts, recdim)
    if (errst /= 0) then
       errmsg="Error in  nf90_inquire"
       call cosp_error(routine_name,errmsg,errcode=errst)
    endif
    do i = 1,ndims
       errst = nf90_Inquire_Dimension(ncid,i,name=dimname(i),len=dimsize(i))
       if (errst /= 0) then
          write(straux, *)  i
          errmsg="Error in nf90_Inquire_Dimension, i: "//trim(straux)
          call cosp_error(routine_name,errmsg)
       endif
       if ((trim(dimname(i)).eq.'level').and.(Nlevels > dimsize(i))) then
          errmsg='Number of levels selected is greater than in input file '//trim(fname)
          call cosp_error(routine_name,errmsg)
       endif
       if (trim(dimname(i)).eq.'point') then
          Lpoint = .true.
          if (Npnts > dimsize(i)) then
             errmsg='Number of points selected is greater than in input file '//trim(fname)
             call cosp_error(routine_name,errmsg)
          endif
       endif
       if (trim(dimname(i)).eq.'lon') then
          Llon = .true.
          Nlon = dimsize(i)
       endif
       if (trim(dimname(i)).eq.'lat') then
          Llat = .true.
          Nlat = dimsize(i)
       endif
    enddo
    
    ! Get lon and lat
    if (Llon.and.Llat) then ! 2D mode
       if ((Npnts) > Nlon*Nlat) Npoints=Nlon*Nlat
       lon = -1.0E30
       lat = -1.0E30
       mode = 2 ! Don't know yet if (lon,lat) or (lat,lon) at this point
    else if (Lpoint) then ! 1D mode
       Nlon = Npoints
       Nlat = Npoints
       mode = 1
    else
       errmsg= trim(fname)//' file contains wrong dimensions'
       call cosp_error(routine_name,errmsg)
    endif
    errst = nf90_inq_varid(ncid, 'lon', vid)
    if (errst /= 0) then
       errmsg="Error in nf90_inq_varid, var: lon"
       call cosp_error(routine_name,errmsg,errcode=errst)
    endif
    errst = nf90_get_var(ncid, vid, lon, start = (/1/), count = (/Nlon/))
    if (errst /= 0) then
       errmsg="Error in nf90_get_var, var: lon"
       call cosp_error(routine_name,errmsg,errcode=errst)
    endif
    errst = nf90_inq_varid(ncid, 'lat', vid)
    if (errst /= 0) then
       errmsg="Error in nf90_inq_varid, var: lat"
       call cosp_error(routine_name,errmsg,errcode=errst)
    endif
    errst = nf90_get_var(ncid, vid, lat, start = (/1/), count = (/Nlat/))
    if (errst /= 0) then
       errmsg="Error in nf90_get_var, var: lat"
       call cosp_error(routine_name,errmsg,errcode=errst)
    endif
    
    ! Get all variables
    do vid = 1,nvars
       vdimid=0
       errst = nf90_Inquire_Variable(ncid, vid, name=vname, ndims=vrank, dimids=vdimid)
       if (errst /= 0) then
          write(straux, *)  vid
          errmsg='Error in nf90_Inquire_Variable, vid '//trim(straux)
          call cosp_error(routine_name,errmsg,errcode=errst)
       endif
       ! Read in into temporary array of correct shape
       !if (present(verbosity).and.(verbosity == 1)) print *, 'Reading '//trim(vname)//' ...'
       if (vrank == 1) then
          Na = dimsize(vdimid(1))
          allocate(x1(Na))
          errst = nf90_get_var(ncid, vid, x1, start=(/1/), count=(/Na/))
       endif
       if (vrank == 2) then
          Na = dimsize(vdimid(1))
          Nb = dimsize(vdimid(2))
          allocate(x2(Na,Nb))
          errst = nf90_get_var(ncid, vid, x2, start=(/1,1/), count=(/Na,Nb/))
       endif
       if (vrank == 3) then
          Na = dimsize(vdimid(1))
          Nb = dimsize(vdimid(2))
          Nc = dimsize(vdimid(3))
          allocate(x3(Na,Nb,Nc))
          errst = nf90_get_var(ncid, vid, x3, start=(/1,1,1/), count=(/Na,Nb,Nc/))
          if ((mode == 2).or.(mode == 3)) then
             if ((Na == Nlon).and.(Nb == Nlat)) then
                mode = 2
             else if ((Na == Nlat).and.(Nb == Nlon)) then
                mode = 3
             else
                errmsg='Wrong mode for variable '//trim(vname)
                call cosp_error(routine_name,errmsg)
             endif
          endif
       endif
       if (vrank == 4) then
          Na = dimsize(vdimid(1))
          Nb = dimsize(vdimid(2))
          Nc = dimsize(vdimid(3))
          Nd = dimsize(vdimid(4))
          allocate(x4(Na,Nb,Nc,Nd))
          errst = nf90_get_var(ncid, vid, x4, start=(/1,1,1,1/), count=(/Na,Nb,Nc,Nd/))
       endif
       if (vrank == 5) then
          Na = dimsize(vdimid(1))
          Nb = dimsize(vdimid(2))
          Nc = dimsize(vdimid(3))
          Nd = dimsize(vdimid(4))
          Ne = dimsize(vdimid(5))
          allocate(x5(Na,Nb,Nc,Nd,Ne))
          errst = nf90_get_var(ncid, vid, x5, start=(/1,1,1,1,1/), count=(/Na,Nb,Nc,Nd,Ne/))
       endif
       if (errst /= 0) then
          write(straux, *)  vid
          errmsg='Error in nf90_get_var, vid '//trim(straux)
          call cosp_error(routine_name,errmsg,errcode=errst)
       endif
       ! Map to the right input argument
       select case (trim(vname))
       case ('pfull')
          if (Lpoint) then
             p(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
          else
             call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=p)
          endif
       case ('phalf')
          if (Lpoint) then
             ph(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
          else
             call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=ph)
          endif
       case ('height')
          if (Lpoint) then
             z(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
          else
             call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=z)
          endif
       case ('height_half')
          if (Lpoint) then
             zh(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
          else
             call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=zh)
          endif
       case ('T_abs')
          if (Lpoint) then
             T(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
          else
             call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=T)
          endif
       case ('qv')
          if (Lpoint) then
             qv(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
          else
             call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=qv)
          endif
       case ('rh')
          if (Lpoint) then
             rh(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
          else
             call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=rh)
          endif
       case ('tca')
          if (Lpoint) then
             tca(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
          else
             call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=tca)
          endif
          tca = tca
       case ('cca')
          if (Lpoint) then
             cca(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
          else
             call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=cca)
          endif
          cca = cca
       case ('mr_lsliq')
          if (Lpoint) then
             mr_lsliq(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
          else
             call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=mr_lsliq)
          endif
       case ('mr_lsice')
          if (Lpoint) then
             mr_lsice(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
          else
             call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=mr_lsice)
          endif
       case ('mr_ccliq')
          if (Lpoint) then
             mr_ccliq(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
          else
             call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=mr_ccliq)
          endif
       case ('mr_ccice')
          if (Lpoint) then
             mr_ccice(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
          else
             call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=mr_ccice)
          endif
       case ('fl_lsrain')
          if (Lpoint) then
             fl_lsrain(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
          else
             call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=fl_lsrain)
          endif
       case ('fl_lssnow')
          if (Lpoint) then
             fl_lssnow(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
          else
             call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=fl_lssnow)
          endif
       case ('fl_lsgrpl')
          if (Lpoint) then
             fl_lsgrpl(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
          else
             call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=fl_lsgrpl)
          endif
       case ('fl_ccrain')
          if (Lpoint) then
             fl_ccrain(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
          else
             call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=fl_ccrain)
          endif
       case ('fl_ccsnow')
          if (Lpoint) then
             fl_ccsnow(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
          else
             call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=fl_ccsnow)
          endif
       case ('dtau_s')
          if (Lpoint) then
             dtau_s(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
          else
             call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=dtau_s)
          endif
       case ('dtau_c')
          if (Lpoint) then
             dtau_c(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
          else
             call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=dtau_c)
          endif
       case ('dem_s')
          if (Lpoint) then
             dem_s(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
          else
             call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=dem_s)
          endif
       case ('dem_c')
          if (Lpoint) then
             dem_c(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
          else
             call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=dem_c)
          endif
       case ('Reff')
          if (Lpoint) then
             Reff(1:Npoints,:,:) = x3(1:Npoints,1:Nlevels,:)
          else
             call map_ll_to_point(Na,Nb,Npoints,x4=x4,y3=Reff)
          endif
       case ('skt')
          if (Lpoint) then
             skt(1:Npoints) = x1(1:Npoints)
          else
             call map_ll_to_point(Na,Nb,Npoints,x2=x2,y1=skt)
          endif
       case ('landmask')
          if (Lpoint) then
             landmask(1:Npoints) = x1(1:Npoints)
          else
             call map_ll_to_point(Na,Nb,Npoints,x2=x2,y1=landmask)
          endif
       case ('mr_ozone')
          if (Lpoint) then
             mr_ozone(1:Npoints,:) = x2(1:Npoints,1:Nlevels)
          else
             call map_ll_to_point(Na,Nb,Npoints,x3=x3,y2=mr_ozone)
          endif
       case ('u_wind')
          if (Lpoint) then
             u_wind(1:Npoints) = x1(1:Npoints)
          else
             call map_ll_to_point(Na,Nb,Npoints,x2=x2,y1=u_wind)
          endif
       case ('v_wind')
          if (Lpoint) then
             v_wind(1:Npoints) = x1(1:Npoints)
          else
             call map_ll_to_point(Na,Nb,Npoints,x2=x2,y1=v_wind)
          endif
       case ('sunlit')
          if (Lpoint) then
             sunlit(1:Npoints) = x1(1:Npoints)
          else
             call map_ll_to_point(Na,Nb,Npoints,x2=x2,y1=sunlit)
          endif
       end select
       ! Free memory
       if (vrank == 1) deallocate(x1)
       if (vrank == 2) deallocate(x2)
       if (vrank == 3) deallocate(x3)
       if (vrank == 4) deallocate(x4)
       if (vrank == 5) deallocate(x5)
    enddo
    
    ! SFC emissivity
    errst = nf90_inq_varid(ncid, 'emsfc_lw', vid)
    if (errst /= 0) then
       if (errst == nf90_enotvar) then ! Does not exist, use 1.0
          emsfc_lw = 1.0
          print *, ' ********* COSP Warning:  emsfc_lw does not exist in input file. Set to 1.0.'
       else  ! Other error, stop
          errmsg='Error in nf90_inq_varid, var: emsfc_lw'
          call cosp_error(routine_name,errmsg,errcode=errst)
       endif
    else
       errst = nf90_get_var(ncid, vid, emsfc_lw)
       if (errst /= 0) then
          errmsg='Error in nf90_get_var, var: emsfc_lw'
          call cosp_error(routine_name,errmsg,errcode=errst)
       endif
    endif
    
    
    ! Fill in the lat/lon vectors with the right values for 2D modes
    ! This might be helpful if the inputs are 2D (gridded) and 
    ! you want outputs in 1D mode
    allocate(plon(Npoints),plat(Npoints))
    if (mode == 2) then !(lon,lat)
       ll = lat
       do j=1,Nb
          do i=1,Na
             k = (j-1)*Na + i
             plon(k) = i  
             plat(k) = j
          enddo
       enddo
       lon(1:Npoints) = lon(plon(1:Npoints))
       lat(1:Npoints) = ll(plat(1:Npoints))
    else if (mode == 3) then !(lat,lon)
       ll = lon
       do j=1,Nb
          do i=1,Na
             k = (j-1)*Na + i
             lon(k) = ll(j)
             lat(k) = lat(i)
          enddo
       enddo
       lon(1:Npoints) = ll(plon(1:Npoints))
       lat(1:Npoints) = lat(plat(1:Npoints))
    endif
    deallocate(plon,plat)
    
    ! Close file
    errst = nf90_close(ncid)
    if (errst /= 0) then
       errmsg='Error in nf90_close'
       call cosp_error(routine_name,errmsg,errcode=errst)
    endif 
  END SUBROUTINE NC_READ_INPUT_FILE

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE map_point_to_ll
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE MAP_POINT_TO_LL(Nx,Ny,geomode,x1,x2,x3,x4,y2,y3,y4,y5)
    ! Input arguments
    integer,intent(in) :: Nx,Ny,geomode
    real(wp),intent(in),optional :: x1(:),x2(:,:),x3(:,:,:), &
         x4(:,:,:,:)
    real(wp),intent(out),optional :: y2(:,:),y3(:,:,:), &
         y4(:,:,:,:),y5(:,:,:,:,:)
    ! Local variables
    integer :: Npoints
    integer :: px(Nx*Ny),py(Nx*Ny)
    integer :: i,j,k,l,m
    integer :: Ni,Nj,Nk,Nl
    integer :: Mi,Mj,Mk,Ml,Mm
     character(len=128) :: proname='MAP_POINT_TO_LL'
     
     Npoints = Nx*Ny
     
     px=0
     py=0
     ! Obtain pointers to do the mapping
     if (geomode == 2) then ! (lon,lat) mode
        do j=1,Ny
           do i=1,Nx
              k = (j-1)*Nx+i
              px(k) = i  
              py(k) = j  
           enddo
        enddo
     else if (geomode == 3) then ! (lon,lat) mode
        do j=1,Nx
           do i=1,Ny
              k = (j-1)*Ny+i
              px(k) = j
              py(k) = i  
           enddo
        enddo
     else
        print *, ' -- '//trim(proname)//': geomode not supported, ',geomode
        stop
     endif
     
     if (present(x1).and.present(y2)) then
        Ni = size(x1,1)
        Mi = size(y2,1)
        Mj = size(y2,2)
        if (Mi*Mj /= Ni) then
           print *, ' -- '//trim(proname)//': Nlon*Nlat /= Npoints (opt 1)'
           stop
        endif
        do i=1,Npoints
           y2(px(i),py(i)) = x1(i)
        enddo
     else if (present(x2).and.present(y3)) then
        Ni = size(x2,1)
        Nj = size(x2,2)
        Mi = size(y3,1)
        Mj = size(y3,2)
        Mk = size(y3,3)
        if (Mi*Mj /= Ni) then
           print *, ' -- '//trim(proname)//': Nlon*Nlat /= Npoints (opt 2)'
           stop
        endif
        if (Nj /= Mk) then
           print *, ' -- '//trim(proname)//': Nj /= Mk (opt 2)'
           print *, Ni,Nj,Mi,Mj,Mk
           stop
        endif
        do k=1,Mk
           do i=1,Npoints
              y3(px(i),py(i),k) = x2(i,k)
           enddo
        enddo
     else if (present(x3).and.present(y4)) then
        Ni = size(x3,1)
        Nj = size(x3,2)
        Nk = size(x3,3)
        Mi = size(y4,1)
        Mj = size(y4,2)
        Mk = size(y4,3)
        Ml = size(y4,4)
        if (Mi*Mj /= Ni) then
           print *, ' -- '//trim(proname)//': Nlon*Nlat /= Npoints (opt 3)'
           stop
        endif
        if (Nj /= Mk) then
           print *, ' -- '//trim(proname)//': Nj /= Mk (opt 3)'
           stop
        endif
        if (Nk /= Ml) then
           print *, ' -- '//trim(proname)//': Nk /= Ml (opt 3)'
           stop
        endif
        do l=1,Ml
           do k=1,Mk
              do i=1,Npoints
                 y4(px(i),py(i),k,l) = x3(i,k,l)
              enddo
           enddo
        enddo
     else if (present(x4).and.present(y5)) then
        Ni = size(x4,1)
        Nj = size(x4,2)
        Nk = size(x4,3)
        Nl = size(x4,4)
        Mi = size(y5,1)
        Mj = size(y5,2)
        Mk = size(y5,3)
        Ml = size(y5,4)
        Mm = size(y5,5)
        if (Mi*Mj /= Ni) then
           print *, ' -- '//trim(proname)//': Nlon*Nlat /= Npoints (opt 4)'
           stop
        endif
        if (Nj /= Mk) then
           print *, ' -- '//trim(proname)//': Nj /= Mk (opt 4)'
           stop
        endif
        if (Nk /= Ml) then
           print *, ' -- '//trim(proname)//': Nk /= Ml (opt 4)'
           stop
        endif
        if (Nl /= Mm) then
           print *, ' -- '//trim(proname)//': Nl /= Mm (opt 4)'
           stop
        endif
        do m=1,Mm
           do l=1,Ml
              do k=1,Mk
                 do i=1,Npoints
                    y5(px(i),py(i),k,l,m) = x4(i,k,l,m)
                 enddo
              enddo
           enddo
        enddo
     else
        print *, ' -- '//trim(proname)//': wrong option'
        stop
     endif
     
     
   END SUBROUTINE MAP_POINT_TO_LL

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE map_ll_to_point
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE MAP_LL_TO_POINT(Nx,Ny,Np,x2,x3,x4,x5,y1,y2,y3,y4)
    ! Input arguments
    integer,intent(in) :: Nx,Ny,Np
    real(wp),intent(in),optional :: x2(:,:),x3(:,:,:), &
         x4(:,:,:,:),x5(:,:,:,:,:)
    real(wp),intent(out),optional :: y1(:),y2(:,:),y3(:,:,:), &
         y4(:,:,:,:)
    ! Local variables
    integer :: px(Nx*Ny),py(Nx*Ny)
    integer :: i,j,k,l,m
    integer :: Ni,Nj,Nk,Nl,Nm
    integer :: Mi,Mj,Mk,Ml
    character(len=128) :: proname='MAP_LL_TO_POINT'
    
    px=0
    py=0
    if (Nx*Ny < Np) then
       print *, ' -- '//trim(proname)//': Nx*Ny < Np'
       stop
    endif
    do j=1,Ny
       do i=1,Nx
          k = (j-1)*Nx+i
          px(k) = i  
          py(k) = j  
       enddo
    enddo
    
    if (present(x2).and.present(y1)) then
       Ni = size(x2,1)
       Nj = size(x2,2)
       Mi = size(y1,1)
       if (Ni*Nj < Mi) then
          print *, ' -- '//trim(proname)//': Nlon*Nlat < Npoints (opt 1)'
          stop
       endif
       do j=1,Np
          y1(j) = x2(px(j),py(j))
       enddo
    else if (present(x3).and.present(y2)) then
       Ni = size(x3,1)
       Nj = size(x3,2)
       Nk = size(x3,3)
       Mi = size(y2,1)
       Mj = size(y2,2)
       if (Ni*Nj < Mi) then
          print *, ' -- '//trim(proname)//': Nlon*Nlat < Npoints (opt 2)'
          stop
       endif
       if (Nk /= Mj) then
          print *, ' -- '//trim(proname)//': Nk /= Mj (opt 2)'
          stop
       endif
       do k=1,Nk
          do j=1,Np
             y2(j,k) = x3(px(j),py(j),k)
          enddo
       enddo
    else if (present(x4).and.present(y3)) then
       Ni = size(x4,1)
       Nj = size(x4,2)
       Nk = size(x4,3)
       Nl = size(x4,4)
       Mi = size(y3,1)
       Mj = size(y3,2)
       Mk = size(y3,3)
       if (Ni*Nj < Mi) then
          print *, ' -- '//trim(proname)//': Nlon*Nlat < Npoints (opt 3)'
          stop
       endif
       if (Nk /= Mj) then
          print *, ' -- '//trim(proname)//': Nk /= Mj (opt 3)'
          stop
       endif
       if (Nl /= Mk) then
          print *, ' -- '//trim(proname)//': Nl /= Mk (opt 3)'
          stop
       endif
       do l=1,Nl
          do k=1,Nk
             do j=1,Np
                y3(j,k,l) = x4(px(j),py(j),k,l)
             enddo
          enddo
       enddo
    else if (present(x5).and.present(y4)) then
       Ni = size(x5,1)
       Nj = size(x5,2)
       Nk = size(x5,3)
       Nl = size(x5,4)
       Nm = size(x5,5)
       Mi = size(y4,1)
       Mj = size(y4,2)
       Mk = size(y4,3)
       Ml = size(y4,4)
       if (Ni*Nj < Mi) then
          print *, ' -- '//trim(proname)//': Nlon*Nlat < Npoints (opt 4)'
          stop
       endif
       if (Nk /= Mj) then
          print *, ' -- '//trim(proname)//': Nk /= Mj (opt 4)'
          stop
       endif
       if (Nl /= Mk) then
          print *, ' -- '//trim(proname)//': Nl /= Mk (opt 4)'
          stop
       endif
       if (Nm /= Ml) then
          print *, ' -- '//trim(proname)//': Nm /= Ml (opt 4)'
          stop
       endif
       do m=1,Nm
          do l=1,Nl
             do k=1,Nk
                do j=1,Np
                   y4(j,k,l,m) = x5(px(j),py(j),k,l,m)
                enddo
             enddo
          enddo
       enddo
    else
       print *, ' -- '//trim(proname)//': wrong option'
       stop
    endif 
  END SUBROUTINE MAP_LL_TO_POINT
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE nc_cmor_init
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE NC_CMOR_INIT(cmor_nl,wmode,nPoints,nColumns,nLevels,nChannels,nLvgrid,      &
                          longitude,latitude,mgrid_zl,mgrid_zu,mgrid_z,cospOUT,geomode,  &
                          Nlon,Nlat,N1,N2,N3,N_OUT_LIST,out_list,lon_axid,lat_axid,      &
                          time_axid,height_axid,height_mlev_axid,grid_id,lonvar_id,      &
                          latvar_id,column_axid,sza_axid,temp_axid,channel_axid,         &
                          dbze_axid,sratio_axid,MISR_CTH_axid,tau_axid,pressure2_axid,   &
                          v1d,v2d,v3d)
    ! Input arguments
    character(len=*),  intent(in) :: cmor_nl,wmode ! Writing mode 'replace' or 'append'
    character(len=32),dimension(n_out_list),intent(in) :: out_list
    type(cosp_outputs),intent(in) :: cospOUT
    integer,intent(in) :: geomode,Nlon,Nlat,N1,N2,N3,N_OUT_LIST,nPoints,nColumns,nLevels,&
                          nChannels,nLvgrid
    real(wp),dimension(nPoints),intent(in) :: longitude,latitude
    real(wp),dimension(nLevels),intent(in) :: mgrid_zl,mgrid_zu,mgrid_z
    
    ! Outputs
    integer,intent(out) :: grid_id,latvar_id,lonvar_id,column_axid,height_axid,dbze_axid,&
                            height_mlev_axid,sratio_axid,tau_axid,pressure2_axid,        &
                            lon_axid,lat_axid,time_axid,sza_axid,MISR_CTH_axid,          &
                            channel_axid,temp_axid
    type(var1d),intent(inout) :: v1d(N1)
    type(var2d),intent(inout) :: v2d(N2)
    type(var3d),intent(inout) :: v3d(N3)

    ! Local variables
    integer ::  profile_axid
    integer :: error_flag,i,j,maxtsteps,Dmax
    logical :: lfound
    real(wp) :: lon_ax(Nlon),lat_ax(Nlat)
    character(len=512) :: inpath,outpath,start_date,model_id,experiment_id,institution,  &
         institute_id,source,calendar,contact,history,comment,table,parent_experiment_id,&
         parent_experiment_rip,forcing
    character(len=2056) :: references
    integer :: realization, nc_action,initialization_method,physics_version
    double precision :: branch_time
    namelist/CMOR/inpath,outpath,start_date,model_id,experiment_id,branch_time,          &
         parent_experiment_id,parent_experiment_rip,forcing,institution,institute_id,    &
         source,calendar,realization,initialization_method,physics_version,contact,      &
         history,comment,references,table,maxtsteps
    real(wp),dimension(:),allocatable :: profile_ax,column_ax,dbze_ax,channel_ax
    real(wp),dimension(:,:),allocatable :: dbze_bounds,vgrid_bounds,sratio_bounds, &
         lon_bounds,lat_bounds,mgrid_bounds
    integer :: d2(2),d3(3),d4(4),d5(5)
    character(len=64) :: pro_name = 'NC_CMOR_INIT'

    
    nc_action = CMOR_APPEND_3
    if (trim(wmode) == 'replace') nc_action = CMOR_REPLACE_3
    

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Allocate memory and compute axes and bounds
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    allocate(column_ax(Ncolumns),dbze_ax(DBZE_BINS),channel_ax(Nchannels),               &
             dbze_bounds(2,DBZE_BINS),vgrid_bounds(2,Nlvgrid),mgrid_bounds(2,Nlevels),   &
             sratio_bounds(2,SR_BINS),lon_bounds(2,Nlon),lat_bounds(2,Nlat))
    
    ! Profile
    if (geomode == 1) then
       allocate(profile_ax(Npoints))
       do i=1,Npoints
          profile_ax(i) = i
       enddo
    endif
    ! Column
    do i=1,Ncolumns
       column_ax(i) = i
    enddo
    ! Channels
    channel_ax = float(iChannel(1:Nchannels))
    ! Radar Ze
    do i=1,DBZE_BINS
       dbze_ax(i) = CFAD_ZE_MIN + CFAD_ZE_WIDTH*(i - 0.5)
    enddo
    do i=1,DBZE_BINS
       dbze_bounds(1,i) = CFAD_ZE_MIN + CFAD_ZE_WIDTH*(i - 1)
       dbze_bounds(2,i) = CFAD_ZE_MIN + CFAD_ZE_WIDTH*i
    enddo
    ! Height of model levels
    do i=1,Nlevels
       mgrid_bounds(1,i) = mgrid_zl(i)
       mgrid_bounds(2,i) = mgrid_zu(i)
    enddo
    ! Height of std grid
    do i=1,Nlvgrid
       vgrid_bounds(1,i) = vgrid_zl(i)
       vgrid_bounds(2,i) = vgrid_zu(i)
    enddo
    ! Lidar scattering ratio bounds (They are output by cosp_cfad_sr->diag_lidar in lmd_ipsl_stats.f90)
    sratio_bounds(2,:)         = cospOUT%calipso_srbval(:) ! srbval contains the upper limits from lmd_ipsl_stats.f90
    sratio_bounds(1,2:SR_BINS) = cospOUT%calipso_srbval(1:SR_BINS-1)
    sratio_bounds(1,1)         = 0.0
    sratio_bounds(2,SR_BINS)   = 1.e5 ! This matches with Chepfer et al., JGR, 2009. However, it is not consistent 
    ! with the upper limit in lmd_ipsl_stats.f90, which is LIDAR_UNDEF-1=998.999
    ! Lat lon axes
    if (geomode == 2) then
       lon_ax = longitude(1:Nlon)
       lat_ax = latitude(1:Npoints:Nlon)
    else if (geomode == 3) then
       lon_ax = longitude(1:Npoints:Nlat)
       lat_ax = latitude(1:Nlat)
    else if (geomode == 4) then
       lon_ax = longitude(1:Nlon)
       lat_ax = latitude(1:Nlat)
    endif
    if (geomode > 1) then
       lon_bounds(1,2:Nlon) = (lon_ax(1:Nlon-1) + lon_ax(2:Nlon))/2.0
       lon_bounds(1,1) = lon_ax(1) - (lon_bounds(1,2) - lon_ax(1))
       lon_bounds(2,1:Nlon-1) = lon_bounds(1,2:Nlon)
       lon_bounds(2,Nlon) = lon_ax(Nlon) + (lon_ax(Nlon) - lon_bounds(2,Nlon-1))
       lat_bounds(1,2:Nlat) = (lat_ax(1:Nlat-1) + lat_ax(2:Nlat))/2.0
       lat_bounds(1,1) = lat_ax(1) - (lat_bounds(1,2) - lat_ax(1))
       lat_bounds(2,1:Nlat-1) = lat_bounds(1,2:Nlat)
       lat_bounds(2,Nlat) = lat_ax(Nlat) + (lat_ax(Nlat) - lat_bounds(2,Nlat-1))
    endif
    
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Read namelist with information for CMOR output file
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    open(10,file=cmor_nl,status='old')
    read(10,nml=cmor)
    close(10)
    
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Specify path for tables and set up other CMOR options
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    error_flag = cmor_setup(inpath=trim(inpath),netcdf_file_action=nc_action,create_subdirectories=0)
    
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Define dataset as output from COSP, and other model details
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    error_flag = cmor_dataset(outpath=trim(outpath),experiment_id=trim(experiment_id),institution=trim(institution), &
         source=trim(source),calendar=trim(calendar),realization=realization,contact=trim(contact), &
         history=trim(history),comment=trim(comment),references=trim(references),model_id=trim(model_id), &
         branch_time=branch_time,parent_experiment_id=trim(parent_experiment_id),forcing=trim(forcing), &
         institute_id=trim(institute_id),parent_experiment_rip=trim(parent_experiment_rip), &
         initialization_method=initialization_method,physics_version=physics_version)
    error_flag = cmor_set_cur_dataset_attribute('cosp_version',trim(COSP_VERSION))
    
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Define axis
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (geomode == 1) then
       profile_axid = cmor_axis(table=table, table_entry='location', units='1', length=Npoints, coord_vals=profile_ax)
    else
       lon_axid = cmor_axis(table=table, table_entry='longitude', units='degrees_east', length=Nlon, coord_vals=lon_ax, &
            cell_bounds = lon_bounds)
       lat_axid = cmor_axis(table=table, table_entry='latitude', units='degrees_north', length=Nlat, coord_vals=lat_ax, &
            cell_bounds = lat_bounds)
    endif
    column_axid  = cmor_axis(table=table, table_entry='column', units='1', length=Ncolumns, coord_vals=column_ax)
    channel_axid = cmor_axis(table=table, table_entry='channel', units='1', length=Nchannels, coord_vals=channel_ax)
    height_axid  = cmor_axis(table=table, table_entry='alt40', units='m', length=Nlvgrid, &
         coord_vals=vgrid_z,cell_bounds=vgrid_bounds)
    temp_axid    = cmor_axis(table=table, table_entry='temp', units='C', length=LIDAR_NTEMP, &
         coord_vals=LIDAR_PHASE_TEMP,cell_bounds=LIDAR_PHASE_TEMP_BNDS)
    dbze_axid    = cmor_axis(table=table, table_entry='dbze', units='dBZ', length=DBZE_BINS, &
         coord_vals=dbze_ax,cell_bounds=dbze_bounds)
    height_mlev_axid  = cmor_axis(table=table, table_entry='alevel', units='1', length=Nlevels, &
         coord_vals=mgrid_z,cell_bounds=mgrid_bounds)
    sratio_axid  = cmor_axis(table=table, table_entry='scatratio', units='1', length=SR_BINS, &
         coord_vals=(sratio_bounds(1,:)+sratio_bounds(2,:))/2.0,cell_bounds=sratio_bounds)
    tau_axid     = cmor_axis(table=table, table_entry='tau', units='1', length=7, &
         coord_vals=isccp_histTauCenters,cell_bounds=isccp_histTauEdges)
    pressure2_axid = cmor_axis(table=table, table_entry='plev7', units='Pa', length=7, &
         coord_vals=isccp_histPresCenters,cell_bounds=isccp_histPresEdges)
    sza_axid   = cmor_axis(table=table, table_entry='sza5', units='degree', length=PARASOL_NREFL, coord_vals=PARASOL_SZA)
    MISR_CTH_axid = cmor_axis(table=table, table_entry='cth16', units='m', length=numMISRHgtBins, &
         coord_vals=misr_histHgtCenters,cell_bounds=misr_histHgtEdges)
    time_axid  = cmor_axis(table=table, table_entry='time1', units='days since '//trim(start_date), length=maxtsteps)
    
    
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Define grid
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (geomode == 1) then
       grid_id = cmor_grid((/profile_axid, time_axid/))
       latvar_id = cmor_time_varying_grid_coordinate(grid_id,'latitude','degrees_north',R_UNDEF)
       lonvar_id = cmor_time_varying_grid_coordinate(grid_id,'longitude','degrees_east' ,R_UNDEF)
       if (grid_id > 0) then
          print *,  '*********** Error, grid_id: ', grid_id
          stop
       endif
    endif
    
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Associate table of variables. Needed here to fill in the table with names.
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (geomode == 1) then
       call nc_cmor_associate_1d(grid_id,time_axid,height_axid,height_mlev_axid,        &
                                 column_axid,sza_axid,temp_axid,channel_axid,dbze_axid, &
                                 sratio_axid,MISR_CTH_axid,tau_axid,pressure2_axid,Nlon,&
                                 Nlat,nPoints,nColumns,nLevels,nChannels,nLvgrid,       &
                                 cospOUT,N1,N2,N3,v1d,v2d,v3d)


    else
       call nc_cmor_associate_2d(lon_axid,lat_axid,time_axid,height_axid,               &
                                 height_mlev_axid,column_axid,sza_axid,temp_axid,       &
                                 channel_axid,dbze_axid,sratio_axid,MISR_CTH_axid,      &
                                 tau_axid,pressure2_axid,Nlon,Nlat,nPoints,nColumns,    &
                                 nLevels,nChannels,nLvgrid,cospOUT,N1,N2,N3,v1d,v2d,v3d)
    endif
    v1d(:)%lout = .false.
    v2d(:)%lout = .false.
    v3d(:)%lout = .false.
    
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Find list of outputs to be written
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    do i=1,N_OUT_LIST
       lfound = .false.
       if (trim(out_list(i)) /= '') then
          do j=1,N1
             if (trim(v1d(j)%name) == trim(out_list(i))) then
                v1d(j)%lout = .true.
                lfound = .true.
                exit
             endif
          enddo
          if (.not.lfound) then
             do j=1,N2
                if (trim(v2d(j)%name) == trim(out_list(i))) then
                   v2d(j)%lout = .true.
                   lfound = .true.
                   exit
                endif
             enddo
          endif
          if (.not.lfound) then
             do j=1,N3
                if (trim(v3d(j)%name) == trim(out_list(i))) then
                   v3d(j)%lout = .true.
                   lfound = .true.
                   exit
                endif
             enddo
          endif
       endif
    enddo
    
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Define variables. Fill in dimensions table first if needed
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! 1D variables
    Dmax = 3
    if (geomode == 1) Dmax=1
    do i=1,N1
       if (v1d(i)%lout) v1d(i)%vid = cmor_variable(table=table, table_entry=v1d(i)%name, units=v1d(i)%units, &
            axis_ids=v1d(i)%dimsid(1:Dmax), missing_value=R_UNDEF)
    enddo
    ! 2D variables
    Dmax = Dmax + 1
    do i=1,N2
       if (v2d(i)%lout) v2d(i)%vid = cmor_variable(table=table, table_entry=v2d(i)%name, units=v2d(i)%units, &
            axis_ids=v2d(i)%dimsid(1:Dmax), missing_value=R_UNDEF)
    enddo
    ! 3D variables
    Dmax = Dmax + 1
    do i=1,N3
       if (v3d(i)%lout) v3d(i)%vid = cmor_variable(table=table, table_entry=v3d(i)%name, units=v3d(i)%units, &
            axis_ids=v3d(i)%dimsid(1:Dmax), missing_value=R_UNDEF)
    enddo
    
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Deallocate memory
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (geomode == 1) deallocate(profile_ax)
    deallocate(column_ax,dbze_ax,channel_ax,dbze_bounds,sratio_bounds, &
         vgrid_bounds,mgrid_bounds,lon_bounds,lat_bounds)
    
  END SUBROUTINE NC_CMOR_INIT

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE nc_cmor_associate_1D
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE NC_CMOR_ASSOCIATE_1D(grid_id,time_axid,height_axid,height_mlev_axid,        &
                                  column_axid,sza_axid,temp_axid,channel_axid,dbze_axid, &
                                  sratio_axid,MISR_CTH_axid,tau_axid,pressure2_axid,Nlon,&
                                  Nlat,nPoints,nColumns,nLevels,nChannels,nLvgrid,       &
                                  cospOUT,N1D,N2D,N3D,v1d,v2d,v3d)
     
    ! Inputs
    integer,intent(in) :: grid_id,time_axid,height_axid,height_mlev_axid,column_axid,    &
                          sza_axid,temp_axid,channel_axid,dbze_axid,sratio_axid,         &
                          MISR_CTH_axid,tau_axid,pressure2_axid,Nlon,Nlat,N1D,N2D,N3D,   &
                          nPoints,nColumns,nLevels,nChannels,Nlvgrid
    type(cosp_outputs),intent(in) :: cospOUT
    type(var1d),intent(inout) :: v1d(N1D+1)
    type(var2d),intent(inout) :: v2d(N2D)
    type(var3d),intent(inout) :: v3d(N3D)

    ! Local variables
    integer :: d2(2),d3(3),d4(4),d5(5)
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Fill in variable info and associate pointers
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 1D variables
    d3 = (/grid_id,0,0/)
    d2 = (/Npoints,0/)
    call construct_var1d('cllcalipso',     d3, d2, cospOUT%calipso_cldlayer(:,1),                v1d(1), units='%')
    call construct_var1d('clmcalipso',     d3, d2, cospOUT%calipso_cldlayer(:,2),                v1d(2), units='%')
    call construct_var1d('clhcalipso',     d3, d2, cospOUT%calipso_cldlayer(:,3),                v1d(3), units='%')
    call construct_var1d('cltcalipso',     d3, d2, cospOUT%calipso_cldlayer(:,4),                v1d(4), units='%')
    call construct_var1d('cltlidarradar',  d3, d2, cospOUT%radar_lidar_tcc,                      v1d(5), units='%')
    call construct_var1d('cltisccp',       d3, d2, cospOUT%isccp_totalcldarea,                   v1d(6), units='%')
    call construct_var1d('pctisccp',       d3, d2, cospOUT%isccp_meanptop,                       v1d(7), units='Pa')
    call construct_var1d('tauisccp',       d3, d2, cospOUT%isccp_meantaucld,                     v1d(8), units='1')
    call construct_var1d('albisccp',       d3, d2, cospOUT%isccp_meanalbedocld,                  v1d(9), units='1')
    call construct_var1d('meantbisccp',    d3, d2, cospOUT%isccp_meantb,                         v1d(10),units='K')
    call construct_var1d('meantbclrisccp', d3, d2, cospOUT%isccp_meantbclr,                      v1d(11),units='K')
    call construct_var1d('cltmodis',       d3, d2, cospOUT%modis_Cloud_Fraction_Total_Mean,      v1d(12),units='%')
    call construct_var1d('clwmodis',       d3, d2, cospOUT%modis_Cloud_Fraction_Water_Mean,      v1d(13),units='%')
    call construct_var1d('climodis',       d3, d2, cospOUT%modis_Cloud_Fraction_Ice_Mean,        v1d(14),units='%')
    call construct_var1d('clhmodis',       d3, d2, cospOUT%modis_Cloud_Fraction_High_Mean,       v1d(15),units='%')
    call construct_var1d('clmmodis',       d3, d2, cospOUT%modis_Cloud_Fraction_Mid_Mean,        v1d(16),units='%')
    call construct_var1d('cllmodis',       d3, d2, cospOUT%modis_Cloud_Fraction_Low_Mean,        v1d(17),units='%')
    call construct_var1d('tautmodis',      d3, d2, cospOUT%modis_Optical_Thickness_Total_Mean,   v1d(18),units='1')
    call construct_var1d('tauwmodis',      d3, d2, cospOUT%modis_Optical_Thickness_Water_Mean,   v1d(19),units='1')
    call construct_var1d('tauimodis',      d3, d2, cospOUT%modis_Optical_Thickness_Ice_Mean,     v1d(20),units='1')
    call construct_var1d('tautlogmodis',   d3, d2, cospOUT%modis_Optical_Thickness_Total_LogMean,v1d(21),units='1')
    call construct_var1d('tauwlogmodis',   d3, d2, cospOUT%modis_Optical_Thickness_Water_LogMean,v1d(22),units='1')
    call construct_var1d('tauilogmodis',   d3, d2, cospOUT%modis_Optical_Thickness_Ice_LogMean,  v1d(23),units='1')
    call construct_var1d('reffclwmodis',   d3, d2, cospOUT%modis_Cloud_Particle_Size_Water_Mean, v1d(24),units='m')
    call construct_var1d('reffclimodis',   d3, d2, cospOUT%modis_Cloud_Particle_Size_Ice_Mean,   v1d(25),units='m')
    call construct_var1d('pctmodis',       d3, d2, cospOUT%modis_Cloud_Top_Pressure_Total_Mean,  v1d(26),units='Pa')
    call construct_var1d('lwpmodis',       d3, d2, cospOUT%modis_Liquid_Water_Path_Mean,         v1d(27),units='kg m-2')
    call construct_var1d('iwpmodis',       d3, d2, cospOUT%modis_Ice_Water_Path_Mean,            v1d(28),units='kg m-2')
    call construct_var1d('cllcalipsoice',  d3, d2, cospOUT%calipso_cldlayerphase(:,1,1),         v1d(30),units='%')
    call construct_var1d('clmcalipsoice',  d3, d2, cospOUT%calipso_cldlayerphase(:,2,1),         v1d(31),units='%')
    call construct_var1d('clhcalipsoice',  d3, d2, cospOUT%calipso_cldlayerphase(:,3,1),         v1d(32),units='%')
    call construct_var1d('cltcalipsoice',  d3, d2, cospOUT%calipso_cldlayerphase(:,4,1),         v1d(33),units='%')
    call construct_var1d('cllcalipsoliq',  d3, d2, cospOUT%calipso_cldlayerphase(:,1,2),         v1d(34),units='%')
    call construct_var1d('clmcalipsoliq',  d3, d2, cospOUT%calipso_cldlayerphase(:,2,2),         v1d(35),units='%')
    call construct_var1d('clhcalipsoliq',  d3, d2, cospOUT%calipso_cldlayerphase(:,3,2),         v1d(36),units='%')
    call construct_var1d('cltcalipsoliq',  d3, d2, cospOUT%calipso_cldlayerphase(:,4,2),         v1d(37),units='%')
    call construct_var1d('cllcalipsoun',   d3, d2, cospOUT%calipso_cldlayerphase(:,1,3),         v1d(38),units='%')
    call construct_var1d('clmcalipsoun',   d3, d2, cospOUT%calipso_cldlayerphase(:,2,3),         v1d(39),units='%')
    call construct_var1d('clhcalipsoun',   d3, d2, cospOUT%calipso_cldlayerphase(:,3,3),         v1d(40),units='%')
    call construct_var1d('cltcalipsoun',   d3, d2, cospOUT%calipso_cldlayerphase(:,4,3),         v1d(41),units='%')
    ! 2D variables
    d4 = (/grid_id,height_axid,0,0/)
    d3 = (/Npoints,Nlvgrid,0/)
    call construct_var2d('clcalipso',      d4, d3, cospOUT%calipso_lidarcld,                     v2d(1),units='%')
    call construct_var2d('clcalipso2',     d4, d3, cospOUT%lidar_only_freq_cloud,                v2d(2),units='%')
    d4 = (/grid_id,height_mlev_axid,0,0/)
    d3 = (/Npoints,Nlevels,0/)
    call construct_var2d('lidarBetaMol532',d4, d3, cospOUT%calipso_beta_mol,                     v2d(3),units='m-1 sr-1')
    d4 = (/grid_id,column_axid,0,0/)
    d3 = (/Npoints,Ncolumns,0/)
    call construct_var2d('boxtauisccp',    d4, d3, cospOUT%isccp_boxtau,                         v2d(4),units='1')
    call construct_var2d('boxptopisccp'  , d4, d3, cospOUT%isccp_boxptop,                        v2d(5),units='Pa')
    d4 = (/grid_id,sza_axid,0,0/)
    d3 = (/Npoints,PARASOL_NREFL,0/)
    call construct_var2d('parasolRefl',    d4, d3, cospOUT%parasolGrid_refl,                     v2d(6),units='1')
    d4 = (/grid_id,height_axid,0,0/)
    d3 = (/Npoints,Nlvgrid,0/)
    call construct_var2d('clcalipsoice',   d4, d3, cospOUT%calipso_lidarcldphase(:,:,1),         v2d(8),units='%')
    call construct_var2d('clcalipsoliq',   d4, d3, cospOUT%calipso_lidarcldphase(:,:,2),         v2d(7),units='%')
    call construct_var2d('clcalipsoun',    d4, d3, cospOUT%calipso_lidarcldphase(:,:,3),         v2d(9),units='%')
    d3 = (/Npoints,LIDAR_NTEMP,0/)
    d4 = (/grid_id,temp_axid,0,0/)
    call construct_var2d('clcalipsotmp',   d4, d3, cospOUT%calipso_lidarcldtmp(:,:,1),           v2d(10),units='%')
    call construct_var2d('clcalipsotmpice',d4, d3, cospOUT%calipso_lidarcldtmp(:,:,2),           v2d(11),units='%')
    call construct_var2d('clcalipsotmpliq',d4, d3, cospOUT%calipso_lidarcldtmp(:,:,3),           v2d(12),units='%')
    call construct_var2d('clcalipsotmpun', d4, d3, cospOUT%calipso_lidarcldtmp(:,:,4),           v2d(13),units='%')
    d4 = (/grid_id,channel_axid,0,0/) 
    d3 = (/Npoints,Nchannels,0/) 
    call construct_var2d('tbrttov',        d4, d3, cospOUT%rttov_tbs,                            v2d(14),units='K') 
    
    ! 3D variables
    d5 = (/grid_id,column_axid,height_mlev_axid,0,0/)
    d4 = (/Npoints,Ncolumns,Nlevels,0/)
    call construct_var3d('dbze94',         d5, d4, cospOUT%cloudsat_Ze_tot,                      v3d(1),units='1')
    call construct_var3d('atb532',         d5, d4, cospOUT%calipso_beta_tot,                     v3d(2),units='m-1 sr-1')
    d5 = (/grid_id,dbze_axid,height_axid,0,0/)
    d4 = (/Npoints,DBZE_BINS,Nlvgrid,0/)
    call construct_var3d('cfadDbze94',     d5, d4, cospOUT%cloudsat_cfad_ze,                      v3d(4),units='1')
    d5 = (/grid_id,sratio_axid,height_axid,0,0/)
    d4 = (/Npoints,SR_BINS,Nlvgrid,0/)
    call construct_var3d('cfadLidarsr532', d5, d4, cospOUT%calipso_cfad_sr,                       v3d(5),units='1')
    d5 = (/grid_id,tau_axid,pressure2_axid,0,0/)
    d4 = (/Npoints,7,7,0/)
    call construct_var3d('clisccp',        d5, d4, cospOUT%isccp_fq,                              v3d(6),units='%')
    call construct_var3d('clmodis',        d5, d4, cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure, v3d(7), units='%')
    d5 = (/grid_id,tau_axid,MISR_CTH_axid,0,0/)
    d4 = (/Npoints,7,numMISRHgtBins,0/)
    call construct_var3d('clMISR',         d5, d4, cospOUT%misr_fq,v3d(8),units='%')
    
  END SUBROUTINE NC_CMOR_ASSOCIATE_1D

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE nc_cmor_associate_2D
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE NC_CMOR_ASSOCIATE_2D(lon_axid,lat_axid,time_axid,height_axid,height_mlev_axid,&
                                  column_axid,sza_axid,temp_axid,channel_axid,dbze_axid,   &
                                  sratio_axid,MISR_CTH_axid,tau_axid,pressure2_axid,Nlon,  &
                                  Nlat,nPoints,nColumns,nLevels,nChannels,nLvgrid,cospOUT, &
                                  N1D,N2D,N3D,v1d,v2d,v3d)
     
    ! Arguments
    integer,intent(in) :: lon_axid,lat_axid,time_axid,height_axid,height_mlev_axid,         &
                          column_axid,sza_axid,temp_axid,channel_axid,dbze_axid,sratio_axid,&
                          MISR_CTH_axid,tau_axid,pressure2_axid,nPoints,nColumns,nLevels,   &
                          nChannels,nLvgrid
    integer,intent(in) :: Nlon,Nlat,N1D,N2D,N3D
    type(cosp_outputs),intent(in)   :: cospOUT
    type(var1d),intent(inout) :: v1d(N1D)
    type(var2d),intent(inout) :: v2d(N2D)
    type(var3d),intent(inout) :: v3d(N3D)

    ! Local variables
    integer :: d2(2),d3(3),d4(4),d5(5)
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Fill in variable info and associate pointers
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! 1D variables
     d3 = (/lon_axid,lat_axid,time_axid/)
     d2 = (/Nlon,Nlat/)
     call construct_var1d('cllcalipso',     d3, d2, cospOUT%calipso_cldlayer(:,1),                v1d(1), units='%')
     call construct_var1d('clmcalipso',     d3, d2, cospOUT%calipso_cldlayer(:,2),                v1d(2), units='%')
     call construct_var1d('clhcalipso',     d3, d2, cospOUT%calipso_cldlayer(:,3),                v1d(3), units='%')
     call construct_var1d('cltcalipso',     d3, d2, cospOUT%calipso_cldlayer(:,4),                v1d(4), units='%')
     call construct_var1d('cltlidarradar',  d3, d2, cospOUT%radar_lidar_tcc,                      v1d(5), units='%')
     call construct_var1d('cltisccp',       d3, d2, cospOUT%isccp_totalcldarea,                   v1d(6), units='%')
     call construct_var1d('pctisccp',       d3, d2, cospOUT%isccp_meanptop,                       v1d(7), units='Pa')
     call construct_var1d('tauisccp',       d3, d2, cospOUT%isccp_meantaucld,                     v1d(8), units='1')
     call construct_var1d('albisccp',       d3, d2, cospOUT%isccp_meanalbedocld,                  v1d(9), units='1')
     call construct_var1d('meantbisccp',    d3, d2, cospOUT%isccp_meantb,                         v1d(10),units='K')
     call construct_var1d('meantbclrisccp', d3, d2, cospOUT%isccp_meantbclr,                      v1d(11),units='K')
     call construct_var1d('cltmodis',       d3, d2, cospOUT%modis_Cloud_Fraction_Total_Mean,      v1d(12),units='%')
     call construct_var1d('clwmodis',       d3, d2, cospOUT%modis_Cloud_Fraction_Water_Mean,      v1d(13),units='%')
     call construct_var1d('climodis',       d3, d2, cospOUT%modis_Cloud_Fraction_Ice_Mean,        v1d(14),units='%')
     call construct_var1d('clhmodis',       d3, d2, cospOUT%modis_Cloud_Fraction_High_Mean,       v1d(15),units='%')
     call construct_var1d('clmmodis',       d3, d2, cospOUT%modis_Cloud_Fraction_Mid_Mean,        v1d(16),units='%')
     call construct_var1d('cllmodis',       d3, d2, cospOUT%modis_Cloud_Fraction_Low_Mean,        v1d(17),units='%')
     call construct_var1d('tautmodis',      d3, d2, cospOUT%modis_Optical_Thickness_Total_Mean,   v1d(18),units='1')
     call construct_var1d('tauwmodis',      d3, d2, cospOUT%modis_Optical_Thickness_Water_Mean,   v1d(19),units='1')
     call construct_var1d('tauimodis',      d3, d2, cospOUT%modis_Optical_Thickness_Ice_Mean,     v1d(20),units='1')
     call construct_var1d('tautlogmodis',   d3, d2, cospOUT%modis_Optical_Thickness_Total_LogMean,v1d(21),units='1')
     call construct_var1d('tauwlogmodis',   d3, d2, cospOUT%modis_Optical_Thickness_Water_LogMean,v1d(22),units='1')
     call construct_var1d('tauilogmodis',   d3, d2, cospOUT%modis_Optical_Thickness_Ice_LogMean,  v1d(23),units='1')
     call construct_var1d('reffclwmodis',   d3, d2, cospOUT%modis_Cloud_Particle_Size_Water_Mean, v1d(24),units='m')
     call construct_var1d('reffclimodis',   d3, d2, cospOUT%modis_Cloud_Particle_Size_Ice_Mean,   v1d(25),units='m')
     call construct_var1d('pctmodis',       d3, d2, cospOUT%modis_Cloud_Top_Pressure_Total_Mean,  v1d(26),units='Pa')
     call construct_var1d('lwpmodis',       d3, d2, cospOUT%modis_Liquid_Water_Path_Mean,         v1d(27),units='kg m-2')
     call construct_var1d('iwpmodis',       d3, d2, cospOUT%modis_Ice_Water_Path_Mean,            v1d(28),units='kg m-2')
     call construct_var1d('cllcalipsoice',  d3, d2, cospOUT%calipso_cldlayerphase(:,1,1),         v1d(29),units='%')
     call construct_var1d('clmcalipsoice',  d3, d2, cospOUT%calipso_cldlayerphase(:,2,1),         v1d(30),units='%')
     call construct_var1d('clhcalipsoice',  d3, d2, cospOUT%calipso_cldlayerphase(:,3,1),         v1d(31),units='%')
     call construct_var1d('cltcalipsoice',  d3, d2, cospOUT%calipso_cldlayerphase(:,4,1),         v1d(32),units='%')
     call construct_var1d('cllcalipsoliq',  d3, d2, cospOUT%calipso_cldlayerphase(:,1,2),         v1d(33),units='%')
     call construct_var1d('clmcalipsoliq',  d3, d2, cospOUT%calipso_cldlayerphase(:,2,2),         v1d(34),units='%')
     call construct_var1d('clhcalipsoliq',  d3, d2, cospOUT%calipso_cldlayerphase(:,3,2),         v1d(35),units='%')
     call construct_var1d('cltcalipsoliq',  d3, d2, cospOUT%calipso_cldlayerphase(:,4,2),         v1d(36),units='%')
     call construct_var1d('cllcalipsoun',   d3, d2, cospOUT%calipso_cldlayerphase(:,1,3),         v1d(37),units='%')
     call construct_var1d('clmcalipsoun',   d3, d2, cospOUT%calipso_cldlayerphase(:,2,3),         v1d(38),units='%')
     call construct_var1d('clhcalipsoun',   d3, d2, cospOUT%calipso_cldlayerphase(:,3,3),         v1d(39),units='%')
     call construct_var1d('cltcalipsoun',   d3, d2, cospOUT%calipso_cldlayerphase(:,4,3),         v1d(40),units='%')
     ! 2D variables
     d4 = (/lon_axid,lat_axid,height_axid,time_axid/)
     d3 = (/Nlon,Nlat,Nlvgrid/)
     call construct_var2d('clcalipso',      d4, d3, cospOUT%calipso_lidarcld,                     v2d(1), units='%')
     call construct_var2d('clcalipso2',     d4, d3, cospOUT%lidar_only_freq_cloud,                v2d(2), units='%')
     d4 = (/lon_axid,lat_axid,height_mlev_axid,time_axid/)
     d3 = (/Nlon,Nlat,Nlevels/)
     call construct_var2d('lidarBetaMol532',d4, d3, cospOUT%calipso_beta_mol,                     v2d(3), units='m-1 sr-1')
     d4 = (/lon_axid,lat_axid,column_axid,time_axid/)
     d3 = (/Nlon,Nlat,Ncolumns/)
     call construct_var2d('boxtauisccp',    d4, d3, cospOUT%isccp_boxtau,                         v2d(4), units='1')
     call construct_var2d('boxptopisccp',   d4, d3, cospOUT%isccp_boxptop,                        v2d(5), units='Pa')
     d4 = (/lon_axid,lat_axid,sza_axid,time_axid/)
     d3 = (/Nlon,Nlat,PARASOL_NREFL/)
     call construct_var2d('parasolRefl',    d4, d3, cospOUT%parasolGrid_refl,                     v2d(6), units='1')
     d4 = (/lon_axid,lat_axid,height_axid,time_axid/)
     d3 = (/Nlon,Nlat,Nlvgrid/)
     call construct_var2d('clcalipsoice',   d4, d3, cospOUT%calipso_lidarcldphase(:,:,1),         v2d(8), units='%')
     call construct_var2d('clcalipsoliq',   d4, d3, cospOUT%calipso_lidarcldphase(:,:,2),         v2d(7), units='%')
     call construct_var2d('clcalipsoun',    d4, d3, cospOUT%calipso_lidarcldphase(:,:,3),         v2d(9), units='%')
     d3 = (/Nlon,Nlat,LIDAR_NTEMP/)
     d4 = (/lon_axid,lat_axid,temp_axid,time_axid/)
     call construct_var2d('clcalipsotmp',   d4, d3, cospOUT%calipso_lidarcldtmp(:,:,1),           v2d(10),units='%')
     call construct_var2d('clcalipsotmpice',d4, d3, cospOUT%calipso_lidarcldtmp(:,:,2),           v2d(11),units='%')
     call construct_var2d('clcalipsotmpliq',d4, d3, cospOUT%calipso_lidarcldtmp(:,:,3),           v2d(12),units='%')
     call construct_var2d('clcalipsotmpun', d4, d3, cospOUT%calipso_lidarcldtmp(:,:,4),           v2d(13),units='%')
     d4 = (/lon_axid,lat_axid,channel_axid,time_axid/)
     d3 = (/Nlon,Nlat,Nchannels/)
     call construct_var2d('tbrttov',        d4, d3, cospOUT%rttov_tbs,                            v2d(14),units='K') 
     
     ! 3D variables
     d5 = (/lon_axid,lat_axid,column_axid,height_mlev_axid,time_axid/)
     d4 = (/Nlon,Nlat,Ncolumns,Nlevels/)
     call construct_var3d('dbze94',         d5, d4, cospOUT%cloudsat_Ze_tot,                      v3d(1), units='1')
     call construct_var3d('atb532',         d5, d4, cospOUT%calipso_beta_tot,                     v3d(2), units='m-1 sr-1')
     d5 = (/lon_axid,lat_axid,dbze_axid,height_axid,time_axid/)
     d4 = (/Nlon,Nlat,DBZE_BINS,Nlvgrid/)
     call construct_var3d('cfadDbze94',     d5, d4, cospOUT%cloudsat_cfad_ze,                     v3d(4), units='1')
     d5 = (/lon_axid,lat_axid,sratio_axid,height_axid,time_axid/)
     d4 = (/Nlon,Nlat,SR_BINS,Nlvgrid/)
     call construct_var3d('cfadLidarsr532', d5, d4, cospOUT%calipso_cfad_sr,                      v3d(5), units='1')
     d5 = (/lon_axid,lat_axid,tau_axid,pressure2_axid,time_axid/)
     d4 = (/Nlon,Nlat,7,7/)
     call construct_var3d('clisccp',        d5, d4, cospOUT%isccp_fq,                             v3d(6), units='%')
     call construct_var3d('clmodis',        d5, d4, cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure, v3d(7), units='%')
     d5 = (/lon_axid,lat_axid,tau_axid,MISR_CTH_axid,time_axid/)
     d4 = (/Nlon,Nlat,7,numMISRHgtBins/)
     call construct_var3d('clMISR',         d5, d4, cospOUT%misr_fq,                              v3d(8), units='%')
   END SUBROUTINE NC_CMOR_ASSOCIATE_2D

  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE construct_var1d
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_VAR1D(name,dimsid,dimssz,pntr,y,units)
    ! Input arguments
    character(len=*),intent(in) :: name
    integer,intent(in) :: dimsid(3)
    integer,intent(in) :: dimssz(2)
    real(wp),dimension(:),target,intent(in) :: pntr
    type(var1d),intent(out) :: y
    character(len=*),optional,intent(in) :: units
    
    y%name =  name
    if (present(units)) y%units   =  units
    y%dimsid =  dimsid
    y%dimssz =  dimssz
    y%pntr => pntr
    
  END SUBROUTINE CONSTRUCT_VAR1D
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE construct_var2d
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_VAR2D(name,dimsid,dimssz,pntr,y,units)
    ! Input arguments
    character(len=*),intent(in) :: name
    integer,intent(in) :: dimsid(4)
    integer,intent(in) :: dimssz(3)
    real(wp),dimension(:,:),target,intent(in) :: pntr
    type(var2d),intent(out) :: y
    character(len=*),optional,intent(in) :: units
    
    y%name =  name
    if (present(units)) y%units   =  units
    y%dimsid =  dimsid
    y%dimssz =  dimssz
    y%pntr => pntr
    
  END SUBROUTINE CONSTRUCT_VAR2D
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE construct_var3d
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE CONSTRUCT_VAR3D(name,dimsid,dimssz,pntr,y,units)
    ! Input arguments
    character(len=*),intent(in) :: name
    integer,intent(in) :: dimsid(5)
    integer,intent(in) :: dimssz(4)
    real(wp),dimension(:,:,:),target,intent(in) :: pntr
    type(var3d),intent(out) :: y
    character(len=*),optional,intent(in) :: units
    
    y%name =  name
    if (present(units)) y%units   =  units
    y%dimsid =  dimsid
    y%dimssz =  dimssz
    y%pntr => pntr
    
  END SUBROUTINE CONSTRUCT_VAR3D
  
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! SUBROUTINE nc_cmor_write_1D
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE NC_CMOR_WRITE_1D(nPoints,longitude,latitude,tbnds,lonvar_id,latvar_id,N1,N2,N3,v1d,v2d,v3d)
     ! Inputs
     double precision,intent(in) :: tbnds(2,1)
     integer,intent(in) :: lonvar_id,latvar_id,N1,N2,N3,nPoints
     real(wp),dimension(nPoints),intent(in) :: longitude,latitude
     type(var1d),intent(inout) :: v1d(N1)
     type(var2d),intent(inout) :: v2d(N1)
     type(var3d),intent(inout) :: v3d(N1)
     ! Local variables
     integer :: error_flag,i
     real(wp),allocatable :: y2(:,:),y3(:,:,:),y4(:,:,:,:)
     character(len=64) :: pro_name = 'NC_WRITE_COSP_1D'

     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! Write to file
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! 1D variables
     do i=1,N1
        
        if (v1d(i)%lout) then
           error_flag = cmor_write(var_id=v1d(i)%vid, data=v1d(i)%pntr, &
                ntimes_passed=1,time_vals=(/1D0/),time_bnds=tbnds)
           if (error_flag < 0) then
              print *,  trim(pro_name)//': Error writing '//trim(v1d(i)%name)
              stop
           endif
           error_flag = cmor_write(var_id=lonvar_id, data=longitude,store_with=v1d(i)%vid, &
                ntimes_passed=1,time_vals=(/1D0/))
           if (error_flag < 0) then
              print *,  trim(pro_name)//': Error writing lon for '//trim(v1d(i)%name)
              stop
           endif
           error_flag = cmor_write(var_id=latvar_id, data=latitude,store_with=v1d(i)%vid, &
                ntimes_passed=1,time_vals=(/1D0/))
           if (error_flag < 0) then
              print *,  trim(pro_name)//': Error writing lat for '//trim(v1d(i)%name)
              stop
           endif
        endif
     enddo
     ! 2D variables
     do i=1,N2
        if (v2d(i)%lout) then
           error_flag = cmor_write(var_id=v2d(i)%vid, data=v2d(i)%pntr, &
                ntimes_passed=1,time_vals=(/1D0/),time_bnds=tbnds)
           if (error_flag < 0) then
              print *,  trim(pro_name)//': Error writing '//trim(v2d(i)%name)
              stop
           endif
           error_flag = cmor_write(var_id=lonvar_id, data=longitude,store_with=v2d(i)%vid, &
                ntimes_passed=1,time_vals=(/1D0/))
           if (error_flag < 0) then
              print *,  trim(pro_name)//': Error writing lon for '//trim(v2d(i)%name)
              stop
           endif
           error_flag = cmor_write(var_id=latvar_id, data=latitude,store_with=v2d(i)%vid, &
                ntimes_passed=1,time_vals=(/1D0/))
           if (error_flag < 0) then
              print *,  trim(pro_name)//': Error writing lat for '//trim(v2d(i)%name)
              stop
           endif
        endif
     enddo
     ! 3D variables
     do i=1,N3
        if (v3d(i)%lout) then
           error_flag = cmor_write(var_id=v3d(i)%vid, data=v3d(i)%pntr, &
                ntimes_passed=1,time_vals=(/1D0/),time_bnds=tbnds)
           if (error_flag < 0) then
              print *,  trim(pro_name)//': Error writing '//trim(v3d(i)%name)
              stop
           endif
           error_flag = cmor_write(var_id=lonvar_id, data=longitude,store_with=v3d(i)%vid, &
                ntimes_passed=1,time_vals=(/1D0/))
           if (error_flag < 0) then
              print *,  trim(pro_name)//': Error writing lon for '//trim(v3d(i)%name)
              stop
           endif
           error_flag = cmor_write(var_id=latvar_id, data=latitude,store_with=v3d(i)%vid, &
                ntimes_passed=1,time_vals=(/1D0/))
           if (error_flag < 0) then
              print *,  trim(pro_name)//': Error writing lat for '//trim(v3d(i)%name)
              stop
           endif
        endif
     enddo
   END SUBROUTINE NC_CMOR_WRITE_1D
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! SUBROUTINE nc_cmor_write_2D
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE NC_CMOR_WRITE_2D(tbnds,geomode,Nlon,Nlat,N1,N2,N3,v1d,v2d,v3d)
     ! Input arguments
     double precision,intent(in) :: tbnds(2,1)
     integer,intent(in) :: geomode,Nlon,Nlat,N1,N2,N3
     type(var1d),intent(inout) :: v1d(N1)
     type(var2d),intent(inout) :: v2d(N2)
     type(var3d),intent(inout) :: v3d(N3)
     !--- Local variables ---
     integer :: error_flag,i
     real(wp),allocatable :: y2(:,:),y3(:,:,:),y4(:,:,:,:)
     character(len=64) :: pro_name = 'NC_WRITE_COSP_2D'
     
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! Write variables to file
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! 1D variables (2D output)
     do i=1,N1
        if (v1d(i)%lout) then
           allocate(y2(v1d(i)%dimssz(1),v1d(i)%dimssz(2)))
           call map_point_to_ll(Nlon,Nlat,geomode,x1=v1d(i)%pntr,y2=y2) ! Regridding
           error_flag = cmor_write(var_id=v1d(i)%vid, data=y2, &
                ntimes_passed=1,time_vals=(/1D0/),time_bnds=tbnds)
           if (error_flag < 0) then
              print *,  trim(pro_name)//': Error writing '//trim(v1d(i)%name)
              stop
           endif
           deallocate(y2)
        endif
     enddo
     ! 2D variables (3D output)
     do i=1,N2
        if (v2d(i)%lout) then
           allocate(y3(v2d(i)%dimssz(1),v2d(i)%dimssz(2),v2d(i)%dimssz(3)))
           call map_point_to_ll(Nlon,Nlat,geomode,x2=v2d(i)%pntr,y3=y3) ! Regridding
           error_flag = cmor_write(var_id=v2d(i)%vid, data=y3, &
                ntimes_passed=1,time_vals=(/1D0/),time_bnds=tbnds)
           if (error_flag < 0) then
              print *,  trim(pro_name)//': Error writing '//trim(v2d(i)%name)
              stop
           endif
           deallocate(y3)
        endif
     enddo
     ! 3D variables (4D output)
     do i=1,N3
        if (v3d(i)%lout) then
           allocate(y4(v3d(i)%dimssz(1),v3d(i)%dimssz(2),v3d(i)%dimssz(3),v3d(i)%dimssz(4)))
           call map_point_to_ll(Nlon,Nlat,geomode,x3=v3d(i)%pntr,y4=y4) ! Regridding
           error_flag = cmor_write(var_id=v3d(i)%vid, data=y4, &
                ntimes_passed=1,time_vals=(/1D0/),time_bnds=tbnds)
           if (error_flag < 0) then
              print *,  trim(pro_name)//': Error writing '//trim(v3d(i)%name)
              stop
           endif
           deallocate(y4)
        endif
     enddo
   END SUBROUTINE NC_CMOR_WRITE_2D
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Subrotuine cosp_error
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_ERROR(routine_name,message,errcode) 
    character(len = *), intent(in) :: routine_name
    character(len = *), intent(in) :: message
    integer,optional :: errcode
    
    write(6, *) " ********** Failure in ", trim(routine_name)
    write(6, *) " ********** ", trim(message)
    if (present(errcode)) write(6, *) " ********** errcode: ", errcode
    flush(6)
    stop
  END SUBROUTINE COSP_ERROR
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE nc_cmor_close
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE NC_CMOR_CLOSE()
    integer :: error_flag
    error_flag = cmor_close()
  END SUBROUTINE NC_CMOR_CLOSE
 
end module mod_cosp_io
