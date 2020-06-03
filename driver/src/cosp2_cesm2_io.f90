module mod_cosp2_cesm2_io
  use cosp_kinds, only: wp
  use mod_cosp,   only: cosp_outputs
  use netcdf
  USE MOD_COSP_CONFIG, ONLY:  Nlvgrid, LIDAR_NCAT, SR_BINS, PARASOL_NREFL, cloudsat_DBZE_BINS, &
       numMODISReffIceBins, numMODISReffLiqBins, ntau, tau_binBounds, tau_binCenters, &
       tau_binEdges,npres, pres_binBounds, pres_binCenters, pres_binEdges, nhgt,      &
       hgt_binBounds, hgt_binCenters, hgt_binEdges, reffLIQ_binCenters,vgrid_z,       &
       reffICE_binCenters, reffLIQ_binCenters, cloudsat_binCenters, PARASOL_SZA,      &
       calipso_binCenters, grLidar532_binCenters, atlid_binCenters,                   &
       CFODD_NDBZE,  CFODD_HISTDBZE, CFODD_HISTDBZEcenters,                           &
       CFODD_NICOD,  CFODD_HISTICOD, CFODD_HISTICODcenters
  implicit none

contains

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE write_cosp2_output
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine write_cosp2_cesm2_output(nLon, nLat, Ncolumns, Nlevels, lev, lon, lat, cospOUT, outFileName)
    integer,intent(in) :: nLon, nLat, Ncolumns, Nlevels
    real(wp),dimension(nLon*nLat),intent(in) :: lon,lat
    real(wp),dimension(Nlevels),intent(in) :: lev
    type(cosp_outputs),intent(in) :: cospOUT
    character(len=256),intent(in) :: outFileName

    integer :: fileID,status,ij,nPoints
    integer,dimension(20)  :: dimID
    integer,dimension(150) :: varID
    integer,dimension(nLon*nLat) :: loc
    integer,dimension(Ncolumns) :: cosp_scol
    integer,dimension(2) :: bnds
    
    nPoints = nLon*nLat
    loc=(/(ij,ij=1,Npoints)/)
    cosp_scol=(/(ij,ij=1,Ncolumns)/)
    bnds=(/(ij,ij=1,2)/)

    ! ---------------------------------------------------------------------------------------
    ! Create output file.
    ! ---------------------------------------------------------------------------------------
    status = nf90_create(path=trim(outFileName),cmode = nf90_clobber,ncid=fileID)
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))

    ! ---------------------------------------------------------------------------------------
    ! Define GLOBAL attributes.
    ! ---------------------------------------------------------------------------------------
    status = nf90_put_att(fileID,NF90_GLOBAL,"Conventions","CF-1.6")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))

    ! ---------------------------------------------------------------------------------------
    ! Define dimensions.
    ! ---------------------------------------------------------------------------------------
    status = nf90_def_dim(fileID,"lon",nLon,dimID(19))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_def_dim(fileID,"lat",nLat,dimID(20))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))    
    status = nf90_def_dim(fileID,"cosp_scol",Ncolumns,dimID(2))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_def_dim(fileID,"lev",Nlevels,dimID(3))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_def_dim(fileID,"levStat",Nlvgrid,dimID(4))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_def_dim(fileID,"tau7",ntau,dimID(5))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_def_dim(fileID,"bnds",2,dimID(6))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_def_dim(fileID,"pres7",npres,dimID(7))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_def_dim(fileID,"hgt16",nhgt,dimID(8))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_def_dim(fileID,"SR_BINS",SR_BINS,dimID(12))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_def_dim(fileID,"PARASOL_NREFL",PARASOL_NREFL,dimID(13))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_def_dim(fileID,"cloudsat_DBZE_BINS",cloudsat_DBZE_BINS,dimID(14))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_def_dim(fileID,"RELIQ_MODIS",numMODISReffLiqBins,dimID(15))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_def_dim(fileID,"REICE_MODIS",numMODISReffIceBins,dimID(16))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_def_dim(fileID,"CFODD_NDBZE",CFODD_NDBZE,dimID(17))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_def_dim(fileID,"CFODD_NICOD",CFODD_NICOD,dimID(18))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))

    ! ---------------------------------------------------------------------------------------
    ! Define varaibles
    ! ---------------------------------------------------------------------------------------
    ! Longitude
    status = nf90_def_var(fileID,"longitude",  nf90_float, (/dimID(19),dimID(20)/),varID(1))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(1),"long_name","longitude")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(1),"units",        "degrees_east")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(1),"standard_name", "longitude")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    ! Latitude
    status = nf90_def_var(fileID,"latitude",   nf90_float, (/dimID(19),dimID(20)/),varID(2))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(2),"long_name","latitude")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(2),"units",        "degrees_north")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(2),"standard_name", "latitude")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    ! Joint-histogram axis
    ! Tau
    status = nf90_def_var(fileID,"tau7",        nf90_float, (/dimID(5)/),varID(3))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(3),"long_name","cloud_optical_depth_bin_centers")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(3),"units",        "1")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(3),"standard_name", "atmosphere_optical_thickness_due_to_cloud")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    ! Tau edges
    status = nf90_def_var(fileID,"tau7_bnds",   nf90_float, (/dimID(6),dimID(5)/),varID(4))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(4),"long_name","cloud_optical_depth_bin_edges")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(4),"units",        "1")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(4),"standard_name", "atmosphere_optical_thickness_due_to_cloud")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    ! Pressure
    status = nf90_def_var(fileID,"pres7",      nf90_float, (/dimID(7)/),varID(5))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(5),"long_name","air_pressure_bin_centers")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(5),"units",        "Pa")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(5),"standard_name", "air_pressure")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    ! Pressure Edges
    status = nf90_def_var(fileID,"pres7_bnds", nf90_float, (/dimID(6),dimID(7)/),varID(6))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(6),"long_name","air_pressure_bin_edges")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(6),"units",        "Pa")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(6),"standard_name", "air_pressure")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    ! Height
    status = nf90_def_var(fileID,"hgt16",       nf90_float, (/dimID(8)/),varID(7))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(7),"long_name","altitude_bin_centers")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(7),"units",        "m")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(7),"standard_name", "altitude")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    ! Height Edges
    status = nf90_def_var(fileID,"hgt16_bnds",  nf90_float, (/dimID(6),dimID(8)/),varID(8))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(8),"long_name","altitude_bin_edges")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(8),"units",        "m")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(8),"standard_name", "altitude")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))    
    ! Levels
    status = nf90_def_var(fileID,"lev",  nf90_float, (/dimID(3)/),varID(84))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(84),"long_name","level indices")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(84),"units",        "1")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    ! Levels for statistical diagnostics (lidar and radar)
    status = nf90_def_var(fileID,"levStat",  nf90_float, (/dimID(4)/),varID(85))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(85),"long_name","level indices")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(85),"units",        "1")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    
    ! Subcolumms
    status = nf90_def_var(fileID,"cosp_scol",  nf90_float, (/dimID(2)/),varID(86))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(86),"long_name","subcolumn indices")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(86),"units",        "1")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    ! Bnds
    status = nf90_def_var(fileID,"bnds",  nf90_float, (/dimID(6)/),varID(82))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(82),"long_name","bounds")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(82),"units",        "1")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    ! loc
    status = nf90_def_var(fileID,"loc",  nf90_float, (/dimID(19),dimID(20)/),varID(83))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(83),"long_name","loc")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(83),"units",        "1")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    
    ! CALIPSO simulator output
    if (associated(cospOUT%calipso_betaperp_tot)) then
       status = nf90_def_var(fileID,"atb532_perp",nf90_float, (/dimID(19),dimID(20),dimID(2),dimID(3)/),varID(9))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(9),"long_name","CALIPSO Attenuated Total Perpendicular Backscatter (532nm)")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(9),"units",        "m-1 sr-1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))      
       status = nf90_put_att(fileID,varID(9),"standard_name", "volume_attenuated_backwards_scattering_function_in_air")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%calipso_beta_tot)) then
       status = nf90_def_var(fileID,"atb532",nf90_float, (/dimID(19),dimID(20),dimID(2),dimID(3)/),varID(10))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(10),"long_name","CALIPSO Attenuated Total Backscatter (532nm)")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(10),"units",        "m-1 sr-1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))      
       status = nf90_put_att(fileID,varID(10),"standard_name", "volume_attenuated_backwards_scattering_function_in_air")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%calipso_tau_tot)) then
       status = nf90_def_var(fileID,"calipso_tau",nf90_float, (/dimID(19),dimID(20),dimID(2),dimID(3)/),varID(11))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(11),"long_name","CALIPSO optical-thickness @ 0.67microns")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(11),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(11),"standard_name", "atmosphere_optical_thickness_due_to_cloud")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))       
    endif
    if (associated(cospOUT%calipso_lidarcldphase)) then
       ! Ice
       status = nf90_def_var(fileID,"clcalipsoice",nf90_float, (/dimID(19),dimID(20),dimID(4)/),varID(58))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(58),"long_name","CALIPSO Cloud Fraction Undetected by CloudSat")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(58),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(58),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       ! Liquid
       status = nf90_def_var(fileID,"clcalipsoliq",nf90_float, (/dimID(19),dimID(20),dimID(4)/),varID(59))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(59),"long_name","CALIPSO Cloud Fraction Undetected by CloudSat")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(59),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(59),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       ! Undefined
       status = nf90_def_var(fileID,"clcalipsoun",nf90_float, (/dimID(19),dimID(20),dimID(4)/),varID(60))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(60),"long_name","CALIPSO Cloud Fraction Undetected by CloudSat")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(60),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(60),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))       
    endif
    if (associated(cospOUT%calipso_cldlayerphase)) then
       ! Ice
       status = nf90_def_var(fileID,"cllcalipsoice",nf90_float, (/dimID(19),dimID(20)/),varID(61))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(61),"long_name","CALIPSO Ice Low Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(61),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(61),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))     
       status = nf90_def_var(fileID,"clmcalipsoice",nf90_float, (/dimID(19),dimID(20)/),varID(62))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(62),"long_name","CALIPSO Ice Mid Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(62),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(62),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))     
       status = nf90_def_var(fileID,"clhcalipsoice",nf90_float, (/dimID(19),dimID(20)/),varID(63))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(63),"long_name","CALIPSO Ice High Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(63),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(63),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))     
       status = nf90_def_var(fileID,"cltcalipsoice",nf90_float, (/dimID(19),dimID(20)/),varID(64))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(64),"long_name","CALIPSO Ice Total Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(64),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(64),"standard_name", "cloud_area_fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       ! Liquid
       status = nf90_def_var(fileID,"cllcalipsoliq",nf90_float, (/dimID(19),dimID(20)/),varID(65))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(65),"long_name","CALIPSO Liquid Low Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(65),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(65),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))     
       status = nf90_def_var(fileID,"clmcalipsoliq",nf90_float, (/dimID(19),dimID(20)/),varID(66))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(66),"long_name","CALIPSO Liquid Mid Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(66),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(66),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))     
       status = nf90_def_var(fileID,"clhcalipsoliq",nf90_float, (/dimID(19),dimID(20)/),varID(67))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(67),"long_name","CALIPSO Liquid High Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(67),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(67),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))     
       status = nf90_def_var(fileID,"cltcalipsoliq",nf90_float, (/dimID(19),dimID(20)/),varID(68))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(68),"long_name","CALIPSO Liquid Total Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(68),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(68),"standard_name", "cloud_area_fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))            
       ! Undetermined
       status = nf90_def_var(fileID,"cllcalipsoun",nf90_float, (/dimID(19),dimID(20)/),varID(69))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(69),"long_name","CALIPSO Undefined-Phase Low Level Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(69),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(69),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))     
       status = nf90_def_var(fileID,"clmcalipsoun",nf90_float, (/dimID(19),dimID(20)/),varID(70))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(70),"long_name","CALIPSO Undefined-Phase Mid Level Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(70),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(70),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))     
       status = nf90_def_var(fileID,"clhcalipsoun",nf90_float, (/dimID(19),dimID(20)/),varID(71))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(71),"long_name","CALIPSO Undefined-Phase High Level Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(71),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(71),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))     
       status = nf90_def_var(fileID,"cltcalipsoun",nf90_float, (/dimID(19),dimID(20)/),varID(72))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(72),"long_name","CALIPSO Undefined-Phase Total Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(72),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(72),"standard_name", "cloud_area_fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))            
    endif
    if (associated(cospOUT%calipso_lidarcldtmp)) then
       status = nf90_def_var(fileID,"clcalipsotmp",nf90_float, (/dimID(19),dimID(20),dimID(4)/),varID(77))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(77),"long_name","CALIPSO Cloud Fraction Undetected by CloudSat")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(77),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(77),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))     
       status = nf90_def_var(fileID,"clcalipsotmpice",nf90_float, (/dimID(19),dimID(20),dimID(4)/),varID(78))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(78),"long_name","CALIPSO Cloud Fraction Undetected by CloudSat")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(78),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(78),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))     
       status = nf90_def_var(fileID,"clcalipsotmpliq",nf90_float, (/dimID(19),dimID(20),dimID(4)/),varID(79))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(79),"long_name","CALIPSO Cloud Fraction Undetected by CloudSat")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(79),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(79),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))            
       status = nf90_def_var(fileID,"clcalipsotmpun",nf90_float, (/dimID(19),dimID(20),dimID(4)/),varID(80))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(80),"long_name","CALIPSO Cloud Fraction Undetected by CloudSat")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(80),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(80),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))            
    endif
    if (associated(cospOUT%calipso_cfad_sr)) then
       status = nf90_def_var(fileID,"cfadLidarsr532",nf90_float, (/dimID(19),dimID(20),dimID(12),dimID(4)/),varID(15))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(15),"long_name","CALIPSO Scattering Ratio CFAD")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(15),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(15),"standard_name", "histogram_of_backscattering_ratio_over_height_above_reference_ellipsoid")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))     
    endif
    if (associated(cospOUT%calipso_lidarcld)) then
       status = nf90_def_var(fileID,"clcalipso",nf90_float, (/dimID(19),dimID(20),dimID(4)/),varID(16))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(16),"long_name","CALIPSO Cloud Area Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(16),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(16),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))     
    endif
    if (associated(cospOUT%calipso_cldlayer)) then
       ! Low-level
       status = nf90_def_var(fileID,"cllcalipso",nf90_float, (/dimID(19),dimID(20)/),varID(73))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(73),"long_name","CALIPSO Low Level Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(73),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(73),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))           
       ! Mid-level
       status = nf90_def_var(fileID,"clmcalipso",nf90_float, (/dimID(19),dimID(20)/),varID(74))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(74),"long_name","CALIPSO Mid Level Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(74),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(74),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))    
       ! High-level
       status = nf90_def_var(fileID,"clhcalipso",nf90_float, (/dimID(19),dimID(20)/),varID(75))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(75),"long_name","CALIPSO High Level Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(75),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(75),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))    
       ! Total
       status = nf90_def_var(fileID,"cltcalipso",nf90_float, (/dimID(19),dimID(20)/),varID(76))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(76),"long_name","CALIPSO Total Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(76),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(76),"standard_name", "cloud_area_fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))           
    endif
    if (associated(cospOUT%calipso_beta_mol)) then
       status = nf90_def_var(fileID,"lidarBetaMol532",nf90_float, (/dimID(19),dimID(20),dimID(3)/),varID(18))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(18),"long_name","CALIPSO Molecular Backscatter Coefficient (532nm)")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(18),"units",        "m-1 sr-1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(18),"standard_name", "volume_attenuated_backwards_scattering_function_in_air_assuming_no_aerosol_or_cloud")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))           
    endif
    if (associated(cospOUT%calipso_srbval) .or. associated(cospOUT%calipso_cfad_sr)) then    
       status = nf90_def_var(fileID,"SR_BINS",nf90_float, (/dimID(12)/),varID(81))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(81),"long_name","CALIPSO Backscattering Ratio (SR) Bin Centers")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(81),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(81),"standard_name", "backscattering_ratio")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))    
       status = nf90_def_var(fileID,"SR_EDGES",nf90_float, (/dimID(6),dimID(12)/),varID(19))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(19),"long_name","CALIPSO Backscattering Ratio (SR) Bin Bounds")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(19),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(19),"standard_name", "backscattering_ratio")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))    
    endif

    !OPAQ diagnostics
    if (associated(cospOUT%calipso_cldtype)) then
       ! Opaque cloud cover
       status = nf90_def_var(fileID,"clopaquecalipso",nf90_float, (/dimID(19),dimID(20)/),varID(91))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(91),"long_name","CALIPSO Opaque Cloud Cover")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(91),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(91),"standard_name", "opaque_cloud_cover")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))           
       ! Thin cloud cover
       status = nf90_def_var(fileID,"clthincalipso",nf90_float, (/dimID(19),dimID(20)/),varID(92))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(92),"long_name","CALIPSO Thin Cloud Cover")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(92),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(92),"standard_name", "thin_cloud_cover")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))    
       ! z_opaque altitude
       status = nf90_def_var(fileID,"clzopaquecalipso",nf90_float, (/dimID(19),dimID(20)/),varID(93))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(93),"long_name","CALIPSO z_opaque Altitude")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(93),"units",        "m")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(93),"standard_name", "z_opaque")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))    
    endif
    !3D cloud fractions
    if (associated(cospOUT%calipso_lidarcldtype)) then
       ! Opaque profiles cloud fraction
       status = nf90_def_var(fileID,"clcalipsoopaque",nf90_float, (/dimID(19),dimID(20),dimID(4)/),varID(94))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(94),"long_name","CALIPSO Opaque Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(94),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(94),"standard_name", "opaque_cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       ! Non-Opaque profiles cloud fraction
       status = nf90_def_var(fileID,"clcalipsothin",nf90_float, (/dimID(19),dimID(20),dimID(4)/),varID(95))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(95),"long_name","CALIPSO Thin Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(95),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(95),"standard_name", "thin_cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       ! z_opaque fraction
       status = nf90_def_var(fileID,"clcalipsozopaque",nf90_float, (/dimID(19),dimID(20),dimID(4)/),varID(96))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(96),"long_name","CALIPSO z_opaque Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(96),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(96),"standard_name", "z_opaque_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))       
       ! Lidar opacity fraction
       status = nf90_def_var(fileID,"clcalipsoopacity",nf90_float, (/dimID(19),dimID(20),dimID(4)/),varID(97))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(97),"long_name","CALIPSO opacity Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(97),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(97),"standard_name", "opacity_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))       
    endif
    if (associated(cospOUT%calipso_cldtypetemp)) then
       ! Opaque cloud temperature
       status = nf90_def_var(fileID,"clopaquetemp",nf90_float, (/dimID(19),dimID(20)/),varID(98))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(98),"long_name","CALIPSO Opaque Cloud Temperature")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(98),"units",        "K")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(98),"standard_name", "opaque_cloud_temperature")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))           
       ! Thin cloud temperature
       status = nf90_def_var(fileID,"clthintemp",nf90_float, (/dimID(19),dimID(20)/),varID(99))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(99),"long_name","CALIPSO Thin Cloud Temperature")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(99),"units",        "K")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(99),"standard_name", "thin_cloud_temperature")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))    
       ! z_opaque temperature
       status = nf90_def_var(fileID,"clzopaquetemp",nf90_float, (/dimID(19),dimID(20)/),varID(100))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(100),"long_name","CALIPSO z_opaque Temperature")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(100),"units",        "K")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(100),"standard_name", "z_opaque_temperature")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))    
    endif
    if (associated(cospOUT%calipso_cldtypemeanz)) then
       ! Opaque cloud altitude
       status = nf90_def_var(fileID,"clopaquemeanz",nf90_float, (/dimID(19),dimID(20)/),varID(101))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(101),"long_name","CALIPSO Opaque Cloud Altitude")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(101),"units",        "m")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(101),"standard_name", "opaque_cloud_altitude")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))           
       ! Thin cloud altitude
       status = nf90_def_var(fileID,"clthinmeanz",nf90_float, (/dimID(19),dimID(20)/),varID(102))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(102),"long_name","CALIPSO Thin Cloud Altitude")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(102),"units",        "m")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(102),"standard_name", "thin_cloud_altitude")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))    
    endif
    if (associated(cospOUT%calipso_cldthinemis)) then
       ! Thin cloud emissivity
       status = nf90_def_var(fileID,"clthinemis",nf90_float, (/dimID(19),dimID(20)/),varID(103))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(103),"long_name","CALIPSO Thin Cloud Emissivity")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(103),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(103),"standard_name", "thin_cloud_emissivity")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))    
    endif
    if (associated(cospOUT%calipso_cldtypemeanzse)) then
       ! Opaque cloud altitude with respect to Surface Elevation
       status = nf90_def_var(fileID,"clopaquemeanzse",nf90_float, (/dimID(19),dimID(20)/),varID(104))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(104),"long_name","CALIPSO Opaque Cloud Altitude with respect to SE")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(104),"units",        "m")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(104),"standard_name", "opaque_cloud_altitude_se")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))           
       ! Thin cloud altitude with respect to Surface Elevation
       status = nf90_def_var(fileID,"clthinmeanzse",nf90_float, (/dimID(19),dimID(20)/),varID(105))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(105),"long_name","CALIPSO Thin Cloud Altitude with respect to SE")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(105),"units",        "m")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(105),"standard_name", "thin_cloud_altitude_se")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))    
       ! z_opaque altitude with respect to Surface Elevation
       status = nf90_def_var(fileID,"clzopaquecalipsose",nf90_float, (/dimID(19),dimID(20)/),varID(106))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(106),"long_name","CALIPSO z_opaque Altitude with respect to SE")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(106),"units",        "m")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(106),"standard_name", "z_opaque_se")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))    
    endif

    !GROUND LIDAR simulator output
    if (associated(cospOUT%grLidar532_cldlayer)) then
       ! Low-level cloud cover
       status = nf90_def_var(fileID,"cllgrLidar532",nf90_float, (/dimID(19),dimID(20)/),varID(107))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(107),"long_name","GROUND LIDAR Low Level Cloud Cover")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(107),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(107),"standard_name", "grLidar532_low_cloud_cover")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))           
       ! Mid-level cloud cover
       status = nf90_def_var(fileID,"clmgrLidar532",nf90_float, (/dimID(19),dimID(20)/),varID(108))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(108),"long_name","GROUND LIDAR Mid Level Cloud Cover")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(108),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(108),"standard_name", "grLidar532_mid_cloud_cover")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))    
       ! High-level cloud cover
       status = nf90_def_var(fileID,"clhgrLidar532",nf90_float, (/dimID(19),dimID(20)/),varID(109))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(109),"long_name","GROUND LIDAR High Level Cloud Cover")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(109),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(109),"standard_name", "grLidar532_high_cloud_cover")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))    
       ! Total cloud cover
       status = nf90_def_var(fileID,"cltgrLidar532",nf90_float, (/dimID(19),dimID(20)/),varID(110))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(110),"long_name","GROUND LIDAR Total Cloud Cover")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(110),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(110),"standard_name", "grLidar532_total_cloud_cover")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))           
    endif
    !3D cloud fraction
    if (associated(cospOUT%grLidar532_lidarcld)) then
       status = nf90_def_var(fileID,"clgrLidar532",nf90_float, (/dimID(19),dimID(20),dimID(4)/),varID(111))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(111),"long_name","GROUND LIDAR Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(111),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(111),"standard_name", "grLidar532_cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    !Molecular backscatter
    if (associated(cospOUT%grLidar532_beta_mol)) then
       status = nf90_def_var(fileID,"lidarBetaMol532gr",nf90_float, (/dimID(19),dimID(20),dimID(3)/),varID(112))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(112),"long_name","GROUND LIDAR  Molecular Backscatter Coefficient (532nm)")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(112),"units",        "m-1 sr-1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(112),"standard_name", "grLidar532_volume_attenuated_backwards_scattering_function_in_air_assuming_no_aerosol_or_cloud")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))           
    endif
    !Height-Intensity histogram (SR)
    if (associated(cospOUT%grLidar532_cfad_sr)) then
       status = nf90_def_var(fileID,"cfadLidarsr532gr",nf90_float, (/dimID(19),dimID(20),dimID(12),dimID(4)/),varID(113))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(113),"long_name","GROUND LIDAR Scattering Ratio CFAD")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(113),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(113),"standard_name", "grLidar532_histogram_of_backscattering_ratio_over_height_above_reference_ellipsoid")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))     
    endif
    if (associated(cospOUT%grLidar532_beta_tot)) then
       status = nf90_def_var(fileID,"atb532gr",nf90_float, (/dimID(19),dimID(20),dimID(2),dimID(3)/),varID(114))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(114),"long_name","GROUND LIDAR Attenuated Total Backscatter (532nm)")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(114),"units",        "m-1 sr-1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))      
       status = nf90_put_att(fileID,varID(114),"standard_name", "volume_attenuated_backwards_scattering_function_in_air")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif

    if (associated(cospOUT%grLidar532_srbval) .or. associated(cospOUT%grLidar532_cfad_sr)) then 
       status = nf90_def_var(fileID,"SR_BINS_GR",nf90_float, (/dimID(12)/),varID(115))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(115),"long_name","GROUND LIDAR Backscattering Ratio (SR) Bin Centers")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(115),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(115),"standard_name", "backscattering_ratio")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))    
       status = nf90_def_var(fileID,"SR_EDGES_GR",nf90_float, (/dimID(6),dimID(12)/),varID(116))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(116),"long_name","GROUND LIDAR Backscattering Ratio (SR) Bin Bounds")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(116),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(116),"standard_name", "backscattering_ratio")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))    
    endif

    !ATLID simulator output
    if (associated(cospOUT%atlid_cldlayer)) then
       ! Low-level cloud cover
       status = nf90_def_var(fileID,"cllatlid",nf90_float, (/dimID(19),dimID(20)/),varID(117))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(117),"long_name","ATLID Low Level Cloud Cover")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(117),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(117),"standard_name", "atlid_low_cloud_cover")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))           
       ! Mid-level cloud cover
       status = nf90_def_var(fileID,"clmatlid",nf90_float, (/dimID(19),dimID(20)/),varID(118))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(118),"long_name","ATLID Mid Level Cloud Cover")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(118),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(118),"standard_name", "atlid_mid_cloud_cover")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))    
       ! High-level cloud cover
       status = nf90_def_var(fileID,"clhatlid",nf90_float, (/dimID(19),dimID(20)/),varID(119))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(119),"long_name","ATLID High Level Cloud Cover")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(119),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(119),"standard_name", "atlid_high_cloud_cover")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))    
       ! Total cloud cover
       status = nf90_def_var(fileID,"cltatlid",nf90_float, (/dimID(19),dimID(20)/),varID(120))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(120),"long_name","ATLID Total Cloud Cover")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(120),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(120),"standard_name", "atlid_total_cloud_cover")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))           
    endif
    !3D cloud fraction
    if (associated(cospOUT%atlid_lidarcld)) then
       status = nf90_def_var(fileID,"clatlid",nf90_float, (/dimID(19),dimID(20),dimID(4)/),varID(121))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(121),"long_name","ATLID Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(121),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(121),"standard_name", "atlid_cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    !Molecular backscatter
    if (associated(cospOUT%atlid_beta_mol)) then
       status = nf90_def_var(fileID,"lidarBetaMol355",nf90_float, (/dimID(19),dimID(20),dimID(3)/),varID(122))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(122),"long_name","ATLID Molecular Backscatter Coefficient (355nm)")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(122),"units",        "m-1 sr-1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(122),"standard_name", "atlid_volume_attenuated_backwards_scattering_function_in_air_assuming_no_aerosol_or_cloud")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))           
    endif
    !Height-Intensity histogram (SR)
    if (associated(cospOUT%atlid_cfad_sr)) then
       status = nf90_def_var(fileID,"cfadLidarsr355",nf90_float, (/dimID(19),dimID(20),dimID(12),dimID(4)/),varID(123))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(123),"long_name","ATLID Scattering Ratio CFAD")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(123),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(123),"standard_name", "atlid_histogram_of_backscattering_ratio_over_height_above_reference_ellipsoid")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))     
    endif
    if (associated(cospOUT%atlid_beta_tot)) then
       status = nf90_def_var(fileID,"atb355",nf90_float, (/dimID(19),dimID(20),dimID(2),dimID(3)/),varID(124))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(124),"long_name","ATLID Attenuated Total Backscatter (355nm)")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(124),"units",        "m-1 sr-1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))      
       status = nf90_put_att(fileID,varID(124),"standard_name", "volume_attenuated_backwards_scattering_function_in_air")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif

    if (associated(cospOUT%atlid_srbval) .or. associated(cospOUT%atlid_cfad_sr)) then 
       status = nf90_def_var(fileID,"SR_BINS_ATLID",nf90_float, (/dimID(12)/),varID(125))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(125),"long_name","ATLID Backscattering Ratio (SR) Bin Centers")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(125),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(125),"standard_name", "backscattering_ratio")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))    
       status = nf90_def_var(fileID,"SR_EDGES_ATLID",nf90_float, (/dimID(6),dimID(12)/),varID(126))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(126),"long_name","ATLID Backscattering Ratio (SR) Bin Bounds")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(126),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(126),"standard_name", "backscattering_ratio")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))    
    endif

    ! PARASOL simulator output
    if (associated(cospOUT%parasolPix_refl)) then
       status = nf90_def_var(fileID,"parasolPix_refl",nf90_float, (/dimID(19),dimID(20),dimID(2),dimID(13)/),varID(20))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(20),"long_name","PARASOL Subcolumn Reflectance")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(20),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(20),"standard_name", "toa_bidirectional_reflectance")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))          
    endif
    if (associated(cospOUT%parasolGrid_refl)) then
       status = nf90_def_var(fileID,"parasolGrid_refl",nf90_float, (/dimID(19),dimID(20),dimID(13)/),varID(21))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(21),"long_name","PARASOL Reflectance")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(21),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(21),"standard_name", "toa_bidirectional_reflectance")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))       
    endif
    if (associated(cospOUT%parasolPix_refl) .or.associated(cospOUT%parasolGrid_refl)) then
       status = nf90_def_var(fileID,"PARASOL_NREFL",nf90_float, (/dimID(13)/),varID(87))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(87),"long_name","PARASOL Solar Zenith Angle")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(87),"units",        "degree")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(87),"standard_name", "solar_zenith_angle")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    
    ! Cloudsat simulator output
    if (associated(cospOUT%cloudsat_Ze_tot)) then
       status = nf90_def_var(fileID,"dbze94",nf90_float, (/dimID(19),dimID(20),dimID(2),dimID(3)/),varID(22))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(22),"long_name","CloudSat Radar Reflectivity")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(22),"units",        "dBZ")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(22),"standard_name", "equivalent_reflectivity_factor")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))       
    endif
    if (associated(cospOUT%cloudsat_cfad_ze)) then
       status = nf90_def_var(fileID,"cfadDbze94",nf90_float, (/dimID(19),dimID(20),dimID(14),dimID(4)/),varID(23))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(23),"long_name","CloudSat Radar reflectivity CFAD")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(23),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(23),"standard_name", "histogram_of_equivalent_reflectivity_factor_over_height_above_reference_ellipsoid")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status)) 
       status = nf90_def_var(fileID,"cloudsat_DBZE_BINS",nf90_float, (/dimID(14)/),varID(88))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(88),"long_name","CloudSat simulator equivalent radar reflectivity factor")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(88),"units",        "dBZ")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(88),"standard_name", "equivalent_reflectivity_factor")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))              
    endif
    if (associated(cospOUT%cloudsat_precip_cover)) then
       status = nf90_def_var(fileID,"ptcloudsatflag0",nf90_float, (/dimID(19),dimID(20)/),varID(127))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(127),"long_name","Cloudsat precipitation cover for flag0")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(127),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_def_var(fileID,"ptcloudsatflag1",nf90_float, (/dimID(19),dimID(20)/),varID(128))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(128),"long_name","Cloudsat precipitation cover for flag1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(128),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_def_var(fileID,"ptcloudsatflag2",nf90_float, (/dimID(19),dimID(20)/),varID(129))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(129),"long_name","Cloudsat precipitation cover for flag2")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(129),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_def_var(fileID,"ptcloudsatflag3",nf90_float, (/dimID(19),dimID(20)/),varID(130))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(130),"long_name","Cloudsat precipitation cover for flag3")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(130),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_def_var(fileID,"ptcloudsatflag4",nf90_float, (/dimID(19),dimID(20)/),varID(131))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(131),"long_name","Cloudsat precipitation cover for flag4")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(131),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_def_var(fileID,"ptcloudsatflag5",nf90_float, (/dimID(19),dimID(20)/),varID(132))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(132),"long_name","Cloudsat precipitation cover for flag5")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(132),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_def_var(fileID,"ptcloudsatflag6",nf90_float, (/dimID(19),dimID(20)/),varID(133))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(133),"long_name","Cloudsat precipitation cover for flag6")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(133),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_def_var(fileID,"ptcloudsatflag7",nf90_float, (/dimID(19),dimID(20)/),varID(134))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(134),"long_name","Cloudsat precipitation cover for flag7")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(134),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_def_var(fileID,"ptcloudsatflag8",nf90_float, (/dimID(19),dimID(20)/),varID(135))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(135),"long_name","Cloudsat precipitation cover for flag8")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(135),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_def_var(fileID,"ptcloudsatflag9",nf90_float, (/dimID(19),dimID(20)/),varID(136))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(136),"long_name","Cloudsat precipitation cover for flag9")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(136),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
    endif
    if (associated(cospOUT%cloudsat_pia)) then
       status = nf90_def_var(fileID,"cloudsatpia",nf90_float, (/dimID(19),dimID(20)/),varID(137))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(137),"long_name","Cloudsat path integrated attenuation")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(137),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))      
    endif

    ! ISCCP simulator outputs
    if (associated(cospOUT%isccp_totalcldarea)) then
       status = nf90_def_var(fileID,"cltisccp",nf90_float, (/dimID(19),dimID(20)/),varID(24))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(24),"long_name","ISCCP Total Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(24),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(24),"standard_name", "cloud_area_fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))              
    endif
    if (associated(cospOUT%isccp_meantb)) then
       status = nf90_def_var(fileID,"meantbisccp",nf90_float, (/dimID(19),dimID(20)/),varID(25))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(25),"long_name","ISCCP all-sky 10.5 micron brightness temperature")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(25),"units",        "K")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(25),"standard_name", "toa_brightness_temperature")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))   
    endif
    if (associated(cospOUT%isccp_meantbclr)) then
       status = nf90_def_var(fileID,"meantbclrisccp",nf90_float, (/dimID(19),dimID(20)/),varID(26))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(26),"long_name","ISCCP clear-sky 10.5 micron brightness temperature")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(26),"units",        "K")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(26),"standard_name", "toa_brightness_temperature_assuming_clear_sky")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))   
    endif
    if (associated(cospOUT%isccp_meanptop)) then
       status = nf90_def_var(fileID,"pctisccp",nf90_float, (/dimID(19),dimID(20)/),varID(27))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(27),"long_name","ISCCP Mean Cloud Top Pressure")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(27),"units",        "hPa")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(27),"standard_name", "air_pressure_at_cloud_top")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))         
    endif
    if (associated(cospOUT%isccp_meantaucld)) then
       status = nf90_def_var(fileID,"tauisccp",nf90_float, (/dimID(19),dimID(20)/),varID(28))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(28),"long_name","ISCCP Mean Optical Depth")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(28),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(28),"standard_name", "atmosphere_optical_thickness_due_to_cloud")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))         
    endif
    if (associated(cospOUT%isccp_meanalbedocld)) then
       status = nf90_def_var(fileID,"albisccp",nf90_float, (/dimID(19),dimID(20)/),varID(29))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(29),"long_name","ISCCP Mean Cloud Albedo")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(29),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(29),"standard_name", "cloud_albedo")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
    endif
    if (associated(cospOUT%isccp_boxtau)) then
       status = nf90_def_var(fileID,"boxtauisccp",nf90_float, (/dimID(19),dimID(20),dimID(2)/),varID(30))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(30),"long_name","ISCCP Subcolumn Optical Depth")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(30),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(30),"standard_name", "atmosphere_optical_thickness_due_to_cloud")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))         
    endif
    if (associated(cospOUT%isccp_boxptop)) then
       status = nf90_def_var(fileID,"boxptopisccp",nf90_float, (/dimID(19),dimID(20),dimID(2)/),varID(31))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(31),"long_name","ISCCP Subcolumn Cloud Top Pressure")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(31),"units",        "Pa")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(31),"standard_name", "air_pressure_at_cloud_top")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))         
    endif
    if (associated(cospOUT%isccp_fq)) then	    
       status = nf90_def_var(fileID,"clisccp",nf90_float, (/dimID(19),dimID(20),dimID(5),dimID(7)/),varID(32))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(32),"long_name","ISCCP joint-PDF of cloud top pressure and optical depth")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(32),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(32),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))         
    endif
    
    ! MISR simulator output
    if (associated(cospOUT%misr_fq)) then
       status = nf90_def_var(fileID,"clMISR",nf90_float, (/dimID(19),dimID(20),dimID(5),dimID(8)/),varID(33))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(33),"long_name","MISR joint-PDF of cloud top pressure and optical depth")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(33),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(33),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
    endif
    if (associated(cospOUT%misr_meanztop)) then
       status = nf90_def_var(fileID,"misr_meanztop",nf90_float, (/dimID(19),dimID(20)/),varID(34))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(34),"long_name","MISR Mean Cloud Top Height")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(34),"units",        "m")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(34),"standard_name", "cloud_top_altitude")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
    endif
    if (associated(cospOUT%misr_cldarea)) then	    
       status = nf90_def_var(fileID,"misr_cldarea",nf90_float, (/dimID(19),dimID(20)/),varID(35))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(35),"long_name","MISR cloud cover")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(35),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(35),"standard_name", "cloud_area_fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))        
    endif

    ! MODIS simulator output
    if (associated(cospOUT%modis_Cloud_Fraction_Total_Mean)) then
       status = nf90_def_var(fileID,"cltmodis",nf90_float, (/dimID(19),dimID(20)/),varID(36))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(36),"long_name","MODIS Total Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(36),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(36),"standard_name", "cloud_area_fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%modis_Cloud_Fraction_Water_Mean)) then
       status = nf90_def_var(fileID,"clwmodis",nf90_float, (/dimID(19),dimID(20)/),varID(37))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(37),"long_name","MODIS Liquid Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(37),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(37),"standard_name", "cloud_area_fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%modis_Cloud_Fraction_Ice_Mean)) then
       status = nf90_def_var(fileID,"climodis",nf90_float, (/dimID(19),dimID(20)/),varID(38))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(38),"long_name","MODIS Ice Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(38),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(38),"standard_name", "cloud_area_fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%modis_Cloud_Fraction_High_Mean)) then
       status = nf90_def_var(fileID,"clhmodis",nf90_float, (/dimID(19),dimID(20)/),varID(39))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(39),"long_name","MODIS High Level Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(39),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(39),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%modis_Cloud_Fraction_Mid_Mean)) then
       status = nf90_def_var(fileID,"clmmodis",nf90_float, (/dimID(19),dimID(20)/),varID(40))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(40),"long_name","MODIS Mid Level Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(40),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(40),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
     endif
    if (associated(cospOUT%modis_Cloud_Fraction_Low_Mean)) then
       status = nf90_def_var(fileID,"cllmodis",nf90_float, (/dimID(19),dimID(20)/),varID(41))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(41),"long_name","MODIS Low Level Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(41),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(41),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%modis_Optical_Thickness_Total_Mean)) then
       status = nf90_def_var(fileID,"tautmodis",nf90_float, (/dimID(19),dimID(20)/),varID(42))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(42),"long_name","MODIS Total Cloud Optical Thickness")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(42),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(42),"standard_name", "atmosphere_optical_thickness_due_to_cloud")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%modis_Optical_Thickness_Water_Mean)) then
       status = nf90_def_var(fileID,"tauwmodis",nf90_float, (/dimID(19),dimID(20)/),varID(43))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(43),"long_name","MODIS Liquid Cloud Optical Thickness")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(43),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(43),"standard_name", "atmosphere_optical_thickness_due_to_cloud")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%modis_Optical_Thickness_Ice_Mean)) then
       status = nf90_def_var(fileID,"tauimodis",nf90_float, (/dimID(19),dimID(20)/),varID(44))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(44),"long_name","MODIS Ice Cloud Optical Thickness")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(44),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(44),"standard_name", "atmosphere_optical_thickness_due_to_cloud")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status)) 
    endif
    if (associated(cospOUT%modis_Optical_Thickness_Total_logMean)) then
       status = nf90_def_var(fileID,"tautlogmodis",nf90_float, (/dimID(19),dimID(20)/),varID(45))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(45),"long_name","MODIS Total Cloud Optical Thickness (Log10 Mean)")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(45),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(45),"standard_name", "atmosphere_optical_thickness_due_to_cloud")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))        
    endif
    if (associated(cospOUT%modis_Optical_Thickness_Water_logMean)) then
       status = nf90_def_var(fileID,"tauwlogmodis",nf90_float, (/dimID(19),dimID(20)/),varID(46))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(46),"long_name","MODIS Liquid Cloud Optical Thickness (Log10 Mean)")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(46),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(46),"standard_name", "atmosphere_optical_thickness_due_to_cloud")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status)) 
    endif
    if (associated(cospOUT%modis_Optical_Thickness_Ice_logMean)) then
       status = nf90_def_var(fileID,"tauilogmodis",nf90_float, (/dimID(19),dimID(20)/),varID(47))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(47),"long_name","MODIS Ice Cloud Optical Thickness (Log10 Mean)")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(47),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(47),"standard_name", "atmosphere_optical_thickness_due_to_cloud")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))        
    endif
    if (associated(cospOUT%modis_Cloud_Particle_Size_Water_Mean)) then
       status = nf90_def_var(fileID,"reffclwmodis",nf90_float, (/dimID(19),dimID(20)/),varID(48))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(48),"long_name","MODIS Liquid Cloud Particle Size")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(48),"units",        "m")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(48),"standard_name", "effective_radius_of_cloud_liquid_water_particle")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status)) 
    endif
    if (associated(cospOUT%modis_Cloud_Particle_Size_Ice_Mean)) then
       status = nf90_def_var(fileID,"reffclimodis",nf90_float, (/dimID(19),dimID(20)/),varID(49))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(49),"long_name","MODIS Ice Cloud Particle Size")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(49),"units",        "m")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status)) 
       status = nf90_put_att(fileID,varID(49),"standard_name", "effective_radius_of_cloud_liquid_water_particle")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status)) 
    endif
    if (associated(cospOUT%modis_Cloud_Top_Pressure_Total_Mean)) then
       status = nf90_def_var(fileID,"pctmodis",nf90_float, (/dimID(19),dimID(20)/),varID(50))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(50),"long_name","MODIS Cloud Top Pressure")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(50),"units",        "hPa")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status)) 
       status = nf90_put_att(fileID,varID(50),"standard_name", "air_pressure_at_cloud_top")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%modis_Liquid_Water_Path_Mean)) then
       status = nf90_def_var(fileID,"lwpmodis",nf90_float, (/dimID(19),dimID(20)/),varID(51))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(51),"long_name","MODIS Cloud Liquid Water Path")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(51),"units",        "kg m-2")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status)) 
       status = nf90_put_att(fileID,varID(51),"standard_name", "atmosphere_cloud_liquid_water_content")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%modis_Ice_Water_Path_Mean)) then
       status = nf90_def_var(fileID,"iwpmodis",nf90_float, (/dimID(19),dimID(20)/),varID(52))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(52),"long_name","MODIS Cloud Ice Water Path")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(52),"units",        "kg m-2")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status)) 
       status = nf90_put_att(fileID,varID(52),"standard_name", "atmosphere_mass_content_of_cloud_ice")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure)) then
       status = nf90_def_var(fileID,"clmodis",nf90_float, (/dimID(19),dimID(20),dimID(5),dimID(7)/),varID(53))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(53),"long_name","MODIS joint-PDF of cloud top pressure and optical depth")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(53),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status)) 
       status = nf90_put_att(fileID,varID(53),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status)) 
    endif
    if (associated(cospOUT%modis_Optical_Thickness_vs_ReffICE)) then
       status = nf90_def_var(fileID,"modis_Optical_Thickness_vs_ReffICE",nf90_float, (/dimID(19),dimID(20),dimID(5),dimID(16)/),varID(54))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(54),"long_name","MODIS Joint-PDF of optical-depth and ice particle size")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(54),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_def_var(fileID,"REICE_MODIS",nf90_float, (/dimID(15)/),varID(89))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(89),"long_name","MODIS Joint-PDF ice particle size bin centers")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(89),"units",        "meters")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))      
    endif
    if (associated(cospOUT%modis_Optical_Thickness_vs_ReffLIQ)) then
       status = nf90_def_var(fileID,"modis_Optical_Thickness_vs_ReffLIQ",nf90_float, (/dimID(19),dimID(20),dimID(5),dimID(15)/),varID(55))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(55),"long_name","MODIS Joint-PDF of optical-depth and liquid particle size")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(55),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_def_var(fileID,"RELIQ_MODIS",nf90_float, (/dimID(15)/),varID(90))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(90),"long_name","MODIS Joint-PDF liquid particle size bin centers")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(90),"units",        "meters")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status)) 
    endif

    ! Joint simulator products.
    if (associated(cospOUT%lidar_only_freq_cloud)) then
       status = nf90_def_var(fileID,"clcalipso2",nf90_float, (/dimID(19),dimID(20),dimID(4)/),varID(56))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(56),"long_name","CALIPSO Cloud Fraction Undetected by CloudSat")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(56),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(56),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status)) 
    endif
    if (associated(cospOUT%radar_lidar_tcc)) then
       status = nf90_def_var(fileID,"cltlidarradar",nf90_float, (/dimID(19),dimID(20)/),varID(57))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(57),"long_name","CALIPSO and CloudSat Total Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(57),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status)) 
       status = nf90_put_att(fileID,varID(57),"standard_name", "cloud_area_fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status)) 
    endif
    if (associated(cospOUT%cloudsat_tcc)) then
       status = nf90_def_var(fileID,"cloudsat_tcc",nf90_float, (/dimID(19),dimID(20)/),varID(138))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(138),"long_name","CloudSat Total Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(138),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status)) 
       status = nf90_put_att(fileID,varID(138),"standard_name", "cloud_area_fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status)) 
    endif
    if (associated(cospOUT%cloudsat_tcc2)) then
       status = nf90_def_var(fileID,"cloudsat_tcc2",nf90_float, (/dimID(19),dimID(20)/),varID(139))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(139),"long_name","CloudSat Total Cloud Fraction (no 1km)")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(139),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status)) 
       status = nf90_put_att(fileID,varID(139),"standard_name", "cloud_area_fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status)) 
    endif    
    ! warm-rain occurrence frequency diagnostics
    if (associated(cospOUT%wr_occfreq_ntotal)) then
       status = nf90_def_var(fileID,"npdfcld",nf90_float, (/dimID(19),dimID(20)/),varID(140))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(140),"long_name","# of Non-Precipitating Clouds")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(140),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(140),"standard_name", "number_of_slwc_nonprecip")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_def_var(fileID,"npdfdrz",nf90_float, (/dimID(19),dimID(20)/),varID(141))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(141),"long_name","# of Drizzling Clouds")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(141),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(141),"standard_name", "number_of_slwc_drizzle")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_def_var(fileID,"npdfrain",nf90_float, (/dimID(19),dimID(20)/),varID(142))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(142),"long_name","# of Precipitating Clouds")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(142),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(142),"standard_name", "number_of_slwc_precip")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    ! Contoured Frequency by Optical Depth Diagram (CFODD)
    if (associated(cospOUT%cfodd_ntotal)) then
       status = nf90_def_var(fileID,"ncfodd1",nf90_float, (/dimID(19),dimID(20),dimID(17),dimID(18)/),varID(143))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(143),"long_name","# of CFODD (05 < Reff < 12 micron)")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(143),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(143),"standard_name", "cfodd_reff_small")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_def_var(fileID,"ncfodd2",nf90_float, (/dimID(19),dimID(20),dimID(17),dimID(18)/),varID(144))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(144),"long_name","# of CFODD (12 < Reff < 18 micron)")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(144),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(144),"standard_name", "cfodd_reff_medium")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_def_var(fileID,"ncfodd3",nf90_float, (/dimID(19),dimID(20),dimID(17),dimID(18)/),varID(145))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(145),"long_name","# of CFODD (18 < Reff < 35 micron)")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(145),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(145),"standard_name", "cfodd_reff_large")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       !! axes for CFODD
       status = nf90_def_var(fileID,"CFODD_NDBZE",nf90_float,(/dimID(17)/),varID(146))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(146),"long_name","CloudSat+MODIS dBZe vs ICOD joint PDF X-axis")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(146),"units",        "dBZ")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(146),"standard_name", "cloudsat_quivalent_reflectivity_factor")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_def_var(fileID,"CFODD_NICOD",nf90_float,(/dimID(18)/),varID(147))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(147),"long_name","CloudSat+MODIS dBZe vs ICOD joint PDF Y-axis")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(147),"units",        "none")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(147),"standard_name", "modis_in-cloud_optical_depth")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif 
    
    ! ---------------------------------------------------------------------------------------
    ! Exit define mode
    ! ---------------------------------------------------------------------------------------
    status = nf90_enddef(fileID)
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))

    ! ---------------------------------------------------------------------------------------
    ! Populate outputs
    ! ---------------------------------------------------------------------------------------
    ! Geo
    status = nf90_put_var(fileID,varID(1),reshape(lon,(/nLon,nLat/)))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_var(fileID,varID(2),reshape(lat,(/nLon,nLat/)))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    ! Joint-histogram axis variables
    status = nf90_put_var(fileID,varID(3),tau_binCenters)
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_var(fileID,varID(4),tau_binEdges)
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_var(fileID,varID(5),pres_binCenters)
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_var(fileID,varID(6),pres_binEdges)
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_var(fileID,varID(7),hgt_binCenters)
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_var(fileID,varID(8),hgt_binEdges)
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_var(fileID,varID(84),lev)
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_var(fileID,varID(85),vgrid_z)
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_var(fileID,varID(86),cosp_scol)
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_var(fileID,varID(82),bnds)
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_var(fileID,varID(83),loc)
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    
    ! CALIPSO simulator output
    if (associated(cospOUT%calipso_betaperp_tot)) then
       status = nf90_put_var(fileID,varID(9),reshape(cospOUT%calipso_betaperp_tot,(/nLon,nLat,nColumns/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%calipso_beta_tot)) then
       status = nf90_put_var(fileID,varID(10),reshape(cospOUT%calipso_beta_tot,(/nLon,nLat,nColumns,nLevels/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%calipso_tau_tot)) then
       status = nf90_put_var(fileID,varID(11),reshape(cospOUT%calipso_tau_tot,(/nLon,nLat,nColumns/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%calipso_lidarcldphase)) then
       status = nf90_put_var(fileID,varID(58),reshape(cospOUT%calipso_lidarcldphase(:,:,1),(/nLon,nLat,nLVgrid/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(59),reshape(cospOUT%calipso_lidarcldphase(:,:,2),(/nLon,nLat,nLVgrid/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(60),reshape(cospOUT%calipso_lidarcldphase(:,:,3),(/nLon,nLat,nLVgrid/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%calipso_cldlayerphase)) then
       ! Ice
       status = nf90_put_var(fileID,varID(61),reshape(cospOUT%calipso_cldlayerphase(:,1,1),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(62),reshape(cospOUT%calipso_cldlayerphase(:,2,1),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(63),reshape(cospOUT%calipso_cldlayerphase(:,3,1),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(64),reshape(cospOUT%calipso_cldlayerphase(:,4,1),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       ! Liquid
       status = nf90_put_var(fileID,varID(65),reshape(cospOUT%calipso_cldlayerphase(:,1,2),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(66),reshape(cospOUT%calipso_cldlayerphase(:,2,2),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(67),reshape(cospOUT%calipso_cldlayerphase(:,3,2),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(68),reshape(cospOUT%calipso_cldlayerphase(:,4,2),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       ! Undetermined
       status = nf90_put_var(fileID,varID(69),reshape(cospOUT%calipso_cldlayerphase(:,1,3),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(70),reshape(cospOUT%calipso_cldlayerphase(:,2,3),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(71),reshape(cospOUT%calipso_cldlayerphase(:,3,3),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(72),reshape(cospOUT%calipso_cldlayerphase(:,4,3),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%calipso_lidarcldtmp)) then
       status = nf90_put_var(fileID,varID(77),reshape(cospOUT%calipso_lidarcldtmp(:,:,1),(/nLon,nLat,nLVgrid/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(78),reshape(cospOUT%calipso_lidarcldtmp(:,:,2),(/nLon,nLat,nLVgrid/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(79),reshape(cospOUT%calipso_lidarcldtmp(:,:,3),(/nLon,nLat,nLVgrid/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(80),reshape(cospOUT%calipso_lidarcldtmp(:,:,4),(/nLon,nLat,nLVgrid/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%calipso_cfad_sr)) then
       status = nf90_put_var(fileID,varID(15),reshape(cospOUT%calipso_cfad_sr,(/nLon,nLat,SR_BINS,nLVgrid/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%calipso_lidarcld)) then
       status = nf90_put_var(fileID,varID(16),reshape(cospOUT%calipso_lidarcld,(/nLon,nLat,nLVgrid/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%calipso_cldlayer)) then
       status = nf90_put_var(fileID,varID(73),reshape(cospOUT%calipso_cldlayer(:,1),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(74),reshape(cospOUT%calipso_cldlayer(:,2),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(75),reshape(cospOUT%calipso_cldlayer(:,3),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(76),reshape(cospOUT%calipso_cldlayer(:,4),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%calipso_beta_mol)) then
       status = nf90_put_var(fileID,varID(18),reshape(cospOUT%calipso_beta_mol,(/nLon,nLat,nLevels/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%calipso_srbval)) then
       status = nf90_put_var(fileID,varID(19),reshape([cospOUT%calipso_srbval(1:SR_BINS),cospOUT%calipso_srbval(2:SR_BINS+1)],(/2,SR_BINS/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%calipso_srbval) .or. associated(cospOUT%calipso_cfad_sr)) then
       status = nf90_put_var(fileID,varID(81),calipso_binCenters)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif

       !OPAQ diagnostics
    if (associated(cospOUT%calipso_cldtype)) then
       status = nf90_put_var(fileID,varID(91),reshape(cospOUT%calipso_cldtype(:,1),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(92),reshape(cospOUT%calipso_cldtype(:,2),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(93),reshape(cospOUT%calipso_cldtype(:,3),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%calipso_lidarcldtype)) then
       status = nf90_put_var(fileID,varID(94),reshape(cospOUT%calipso_lidarcldtype(:,:,1),(/nLon,nLat,nLVgrid/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(95),reshape(cospOUT%calipso_lidarcldtype(:,:,2),(/nLon,nLat,nLVgrid/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(96),reshape(cospOUT%calipso_lidarcldtype(:,:,3),(/nLon,nLat,nLVgrid/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(97),reshape(cospOUT%calipso_lidarcldtype(:,:,4),(/nLon,nLat,nLVgrid/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%calipso_cldtypetemp)) then
       status = nf90_put_var(fileID,varID(98),reshape(cospOUT%calipso_cldtypetemp(:,1),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(99),reshape(cospOUT%calipso_cldtypetemp(:,2),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(100),reshape(cospOUT%calipso_cldtypetemp(:,3),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%calipso_cldtypemeanz)) then
       status = nf90_put_var(fileID,varID(101),reshape(cospOUT%calipso_cldtypemeanz(:,1),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(102),reshape(cospOUT%calipso_cldtypemeanz(:,2),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%calipso_cldthinemis)) then
       status = nf90_put_var(fileID,varID(103),reshape(cospOUT%calipso_cldthinemis,(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%calipso_cldtypemeanzse)) then
       status = nf90_put_var(fileID,varID(104),reshape(cospOUT%calipso_cldtypemeanzse(:,1),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(105),reshape(cospOUT%calipso_cldtypemeanzse(:,2),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(106),reshape(cospOUT%calipso_cldtypemeanzse(:,3),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif

    ! GROUND LIDAR simulator output
    if (associated(cospOUT%grLidar532_cldlayer)) then
       status = nf90_put_var(fileID,varID(107),reshape(cospOUT%grLidar532_cldlayer(:,1),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(108),reshape(cospOUT%grLidar532_cldlayer(:,2),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(109),reshape(cospOUT%grLidar532_cldlayer(:,3),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(110),reshape(cospOUT%grLidar532_cldlayer(:,4),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%grLidar532_lidarcld)) then
       status = nf90_put_var(fileID,varID(111),reshape(cospOUT%grLidar532_lidarcld,(/nLon,nLat,nLVgrid/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%grLidar532_beta_mol)) then
       status = nf90_put_var(fileID,varID(112),reshape(cospOUT%grLidar532_beta_mol,(/nLon,nLat,nLevels/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%grLidar532_cfad_sr)) then
       status = nf90_put_var(fileID,varID(113),reshape(cospOUT%grLidar532_cfad_sr,(/nLon,nLat,SR_BINS,nLVgrid/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%grLidar532_beta_tot)) then
       status = nf90_put_var(fileID,varID(114),reshape(cospOUT%grLidar532_beta_tot,(/nLon,nLat,nColumns,nLevels/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif

    if (associated(cospOUT%grLidar532_srbval)) then
       status = nf90_put_var(fileID,varID(116),reshape([cospOUT%grLidar532_srbval(1:SR_BINS),cospOUT%grLidar532_srbval(2:SR_BINS+1)],(/2,SR_BINS/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%grLidar532_srbval) .or. associated(cospOUT%grLidar532_cfad_sr)) then
       status = nf90_put_var(fileID,varID(115),grLidar532_binCenters)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif

    ! ATLID simulator output
    if (associated(cospOUT%atlid_cldlayer)) then
       status = nf90_put_var(fileID,varID(117),reshape(cospOUT%atlid_cldlayer(:,1),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(118),reshape(cospOUT%atlid_cldlayer(:,2),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(119),reshape(cospOUT%atlid_cldlayer(:,3),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(120),reshape(cospOUT%atlid_cldlayer(:,4),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%atlid_lidarcld)) then
       status = nf90_put_var(fileID,varID(121),reshape(cospOUT%atlid_lidarcld,(/nLon,nLat,nLVgrid/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%atlid_beta_mol)) then
       status = nf90_put_var(fileID,varID(122),reshape(cospOUT%atlid_beta_mol,(/nLon,nLat,nLevels/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%atlid_cfad_sr)) then
       status = nf90_put_var(fileID,varID(123),reshape(cospOUT%atlid_cfad_sr,(/nLon,nLat,SR_BINS,nLVgrid/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%atlid_beta_tot)) then
       status = nf90_put_var(fileID,varID(124),reshape(cospOUT%atlid_beta_tot,(/nLon,nLat,nColumns,nLevels/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif

    if (associated(cospOUT%atlid_srbval)) then
       status = nf90_put_var(fileID,varID(126),reshape([cospOUT%atlid_srbval(1:SR_BINS),cospOUT%atlid_srbval(2:SR_BINS+1)],(/2,SR_BINS/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%atlid_srbval) .or. associated(cospOUT%atlid_cfad_sr)) then
       status = nf90_put_var(fileID,varID(125),atlid_binCenters)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif

    ! PARASOL simulator output
    if (associated(cospOUT%parasolPix_refl)) then
       status = nf90_put_var(fileID,varID(20),reshape(cospOUT%parasolPix_refl,(/nLon,nLat,nColumns,PARASOL_NREFL/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%parasolGrid_refl)) then
       status = nf90_put_var(fileID,varID(21),reshape(cospOUT%parasolGrid_refl,(/nLon,nLat,PARASOL_NREFL/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%parasolPix_refl) .or.associated(cospOUT%parasolGrid_refl)) then
       status = nf90_put_var(fileID,varID(87),PARASOL_SZA)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    
    ! Cloudsat simulator output
    if (associated(cospOUT%cloudsat_Ze_tot)) then
       status = nf90_put_var(fileID,varID(22),reshape(cospOUT%cloudsat_Ze_tot,(/nLon,nLat,nColumns,nLevels/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%cloudsat_cfad_ze)) then
       status = nf90_put_var(fileID,varID(23),reshape(cospOUT%cloudsat_cfad_ze,(/nLon,nLat,cloudsat_DBZE_BINS,nLvgrid/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(88),cloudsat_binCenters)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%cloudsat_precip_cover)) then
       status = nf90_put_var(fileID,varID(127),reshape(cospOUT%cloudsat_precip_cover(:,1),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(128),reshape(cospOUT%cloudsat_precip_cover(:,2),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(129),reshape(cospOUT%cloudsat_precip_cover(:,3),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(130),reshape(cospOUT%cloudsat_precip_cover(:,4),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(131),reshape(cospOUT%cloudsat_precip_cover(:,5),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(132),reshape(cospOUT%cloudsat_precip_cover(:,6),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(133),reshape(cospOUT%cloudsat_precip_cover(:,7),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(134),reshape(cospOUT%cloudsat_precip_cover(:,8),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(135),reshape(cospOUT%cloudsat_precip_cover(:,9),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(136),reshape(cospOUT%cloudsat_precip_cover(:,10),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%cloudsat_pia)) then
       status = nf90_put_var(fileID,varID(137),reshape(cospOUT%cloudsat_pia,(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif

    if (associated(cospOUT%isccp_totalcldarea)) then
       status = nf90_put_var(fileID,varID(24),reshape(cospOUT%isccp_totalcldarea,(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%isccp_meantb)) then
       status = nf90_put_var(fileID,varID(25),reshape(cospOUT%isccp_meantb,(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%isccp_meantbclr)) then
       status = nf90_put_var(fileID,varID(26),reshape(cospOUT%isccp_meantbclr,(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%isccp_meanptop)) then
       status = nf90_put_var(fileID,varID(27),reshape(cospOUT%isccp_meanptop,(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%isccp_meantaucld)) then
       status = nf90_put_var(fileID,varID(28),reshape(cospOUT%isccp_meantaucld,(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%isccp_meanalbedocld)) then
       status = nf90_put_var(fileID,varID(29),reshape(cospOUT%isccp_meanalbedocld,(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%isccp_boxtau)) then
       status = nf90_put_var(fileID,varID(30),reshape(cospOUT%isccp_boxtau,(/nLon,nLat,nColumns/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%isccp_boxptop)) then
       status = nf90_put_var(fileID,varID(31),reshape(cospOUT%isccp_boxptop,(/nLon,nLat,nColumns/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%isccp_fq)) then	    
       status = nf90_put_var(fileID,varID(32),reshape(cospOUT%isccp_fq,(/nLon,nLat,ntau,npres/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    ! MISR simulator output
    if (associated(cospOUT%misr_fq)) then
       status = nf90_put_var(fileID,varID(33),reshape(cospOUT%misr_fq,(/nLon,nLat,ntau,nhgt/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%misr_meanztop)) then
       status = nf90_put_var(fileID,varID(34),reshape(cospOUT%misr_meanztop,(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%misr_cldarea)) then	    
       status = nf90_put_var(fileID,varID(35),reshape(cospOUT%misr_cldarea,(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    ! MODIS simulator output
    if (associated(cospOUT%modis_Cloud_Fraction_Total_Mean)) then
       status = nf90_put_var(fileID,varID(36),reshape(cospOUT%modis_Cloud_Fraction_Total_Mean,(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%modis_Cloud_Fraction_Water_Mean)) then
       status = nf90_put_var(fileID,varID(37),reshape(cospOUT%modis_Cloud_Fraction_Water_Mean,(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%modis_Cloud_Fraction_Ice_Mean)) then
       status = nf90_put_var(fileID,varID(38),reshape(cospOUT%modis_Cloud_Fraction_Ice_Mean,(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%modis_Cloud_Fraction_High_Mean)) then
       status = nf90_put_var(fileID,varID(39),reshape(cospOUT%modis_Cloud_Fraction_High_Mean,(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%modis_Cloud_Fraction_Mid_Mean)) then
       status = nf90_put_var(fileID,varID(40),reshape(cospOUT%modis_Cloud_Fraction_Mid_Mean,(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%modis_Cloud_Fraction_Low_Mean)) then
       status = nf90_put_var(fileID,varID(41),reshape(cospOUT%modis_Cloud_Fraction_Low_Mean,(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%modis_Optical_Thickness_Total_Mean)) then
       status = nf90_put_var(fileID,varID(42),reshape(cospOUT%modis_Optical_Thickness_Total_Mean,(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%modis_Optical_Thickness_Water_Mean)) then
       status = nf90_put_var(fileID,varID(43),reshape(cospOUT%modis_Optical_Thickness_Water_Mean,(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%modis_Optical_Thickness_Ice_Mean)) then 
       status = nf90_put_var(fileID,varID(44),reshape(cospOUT%modis_Optical_Thickness_Ice_Mean,(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))      
    endif
    if (associated(cospOUT%modis_Optical_Thickness_Total_LogMean)) then
       status = nf90_put_var(fileID,varID(45),reshape(cospOUT%modis_Optical_Thickness_Total_LogMean,(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%modis_Optical_Thickness_Water_LogMean)) then
       status = nf90_put_var(fileID,varID(46),reshape(cospOUT%modis_Optical_Thickness_Water_LogMean,(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (associated(cospOUT%modis_Optical_Thickness_Ice_LogMean)) then 
       status = nf90_put_var(fileID,varID(47),reshape(cospOUT%modis_Optical_Thickness_Ice_LogMean,(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))      
    endif
    if (associated(cospOUT%modis_Cloud_Particle_Size_Water_Mean)) then
       status = nf90_put_var(fileID,varID(48),reshape(cospOUT%modis_Cloud_Particle_Size_Water_Mean,(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
    endif
    if (associated(cospOUT%modis_Cloud_Particle_Size_Ice_Mean)) then
       status = nf90_put_var(fileID,varID(49),reshape(cospOUT%modis_Cloud_Particle_Size_Ice_Mean,(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
    endif
    if (associated(cospOUT%modis_Cloud_Top_Pressure_Total_Mean)) then
       status = nf90_put_var(fileID,varID(50),reshape(cospOUT%modis_Cloud_Top_Pressure_Total_Mean,(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
    endif
    if (associated(cospOUT%modis_Liquid_Water_Path_Mean)) then
       status = nf90_put_var(fileID,varID(51),reshape(cospOUT%modis_Liquid_Water_Path_Mean,(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
    endif
    if (associated(cospOUT%modis_Ice_Water_Path_Mean)) then
       status = nf90_put_var(fileID,varID(52),reshape(cospOUT%modis_Ice_Water_Path_Mean,(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
    endif
    if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure)) then          			   
       status = nf90_put_var(fileID,varID(53),reshape(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure,(/nLon,nLat,ntau,npres/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
    endif
    if (associated(cospOUT%modis_Optical_Thickness_vs_ReffICE)) then
       status = nf90_put_var(fileID,varID(54),reshape(cospOUT%modis_Optical_Thickness_vs_ReffICE,(/nLon,nLat,ntau,numMODISReffIceBins/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(89),reffICE_binCenters)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
    endif
    if (associated(cospOUT%modis_Optical_Thickness_vs_ReffLIQ)) then
       status = nf90_put_var(fileID,varID(55),reshape(cospOUT%modis_Optical_Thickness_vs_ReffLIQ,(/nLon,nLat,ntau,numMODISReffLiqBins/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_var(fileID,varID(90),reffLIQ_binCenters)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif

    if (associated(cospOUT%lidar_only_freq_cloud)) then
       status = nf90_put_var(fileID,varID(56),reshape(cospOUT%lidar_only_freq_cloud,(/nLon,nLat,nLvgrid/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
    endif
    if (associated(cospOUT%radar_lidar_tcc)) then
       status = nf90_put_var(fileID,varID(57),reshape(cospOUT%radar_lidar_tcc,(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
    endif
    if (associated(cospOUT%cloudsat_tcc)) then
       status = nf90_put_var(fileID,varID(138),reshape(cospOUT%cloudsat_tcc,(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
    endif
    if (associated(cospOUT%cloudsat_tcc2)) then
       status = nf90_put_var(fileID,varID(139),reshape(cospOUT%cloudsat_tcc2,(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
    endif

    ! Cloudsat+MODIS Joint simulators output
      !! warm-rain occurrence frequency diagnostics
    if (associated(cospOUT%wr_occfreq_ntotal)) then
       status = nf90_put_var(fileID,varID(140),reshape(cospOUT%wr_occfreq_ntotal(:,1),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(141),reshape(cospOUT%wr_occfreq_ntotal(:,2),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(142),reshape(cospOUT%wr_occfreq_ntotal(:,3),(/nLon,nLat/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
      !! Contoured Frequency by Optical Depth Diagram (CFODD)
    if (associated(cospOUT%cfodd_ntotal)) then
       status = nf90_put_var(fileID,varID(143),reshape(cospOUT%cfodd_ntotal(:,:,:,1),(/nLon,nLat,CFODD_NDBZE,CFODD_NICOD/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(144),reshape(cospOUT%cfodd_ntotal(:,:,:,2),(/nLon,nLat,CFODD_NDBZE,CFODD_NICOD/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(145),reshape(cospOUT%cfodd_ntotal(:,:,:,3),(/nLon,nLat,CFODD_NDBZE,CFODD_NICOD/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(146),CFODD_HISTDBZEcenters)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(147),CFODD_HISTICODcenters)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif

    ! Close file
    status = nf90_close(fileID)
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))

  end subroutine write_cosp2_cesm2_output

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Subroutine read_cesm2_input
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine read_cesm2_input(fname, nLon, nLat, nPts, nLev, skt, sunlit,        &
      landmask, lon, lat, p, ph, z, zh, T, qv, tca, cca, mr_lsliq, mr_lsice, mr_ccliq,   &
      mr_ccice, fl_lsrain, fl_lssnow, fl_lsgrpl, fl_ccrain, fl_ccsnow, dtau_s, dtau_c,   &
      dem_s, dem_c, dem_s_snow, dtau_s_snow, mr_ozone, Reff)
      
     ! Inputs
     character(len=512),intent(in) :: fname ! File name

	 ! Outputs
	 integer, intent(out) :: nLon, nLat, nPts, nLev
	 real(wp),dimension(:), allocatable, intent(out) :: skt, sunlit, landmask, lon, lat
     real(wp),dimension(:,:), allocatable, intent(out) :: p, ph, z, zh, T, qv, tca, cca, &
         mr_lsliq, mr_lsice, mr_ccliq, mr_ccice, fl_lsrain, fl_lssnow, fl_lsgrpl,        &
         fl_ccrain, fl_ccsnow, dtau_s, dtau_c, dem_s, dem_c, mr_ozone, dem_s_snow,       &
         dtau_s_snow
    real(wp),dimension(:,:,:), allocatable, intent(out) :: Reff
    
    ! Local vars
	integer :: err,ncid,varID
	integer, dimension(3) :: dimID
    real(wp),dimension(:,:,:),allocatable :: var3dtemp
    real(wp),dimension(:,:),allocatable :: var2dtemp
    real(wp),dimension(:),allocatable :: lon1d,lat1d
    character(len=128) :: errmsg
    character(len=24) :: routine_name = "read_cesm2_input"
    
    !
    ! Open file
    !
    err = nf90_open(fname, nf90_nowrite, ncid)
    if (err /= 0) then
       errmsg="Couldn't open "//trim(fname)
       call cosp_error(routine_name,errmsg)
    endif

	!
    ! Get information about dimensions.
	! 
    err = nf90_inq_dimid(ncid,'lon',dimID(1))
    err = nf90_Inquire_Dimension(ncid,dimID(1),len=nLon)    
	!
    err = nf90_inq_dimid(ncid,'lat',dimID(2))
    err = nf90_Inquire_Dimension(ncid,dimID(2),len=nLat)    
	!
    err = nf90_inq_dimid(ncid,'lev',dimID(3))
    err = nf90_Inquire_Dimension(ncid,dimID(3),len=nLev)    
    
    ! Problem size
    nPts = nLon*nLat
    
    ! Allocate space
    ! 1D
    allocate(skt(nPts),sunlit(nPts),landmask(nPts),lon(nPts),lat(nPts))
    ! 2D
	allocate(p(nPts,nLev), ph(nPts,nLev), z(nPts,nLev), zh(nPts,nLev), T(nPts,nLev),     &
	   qv(nPts,nLev), tca(nPts,nLev), cca(nPts,nLev), mr_lsliq(nPts,nLev),               &
	   mr_lsice(nPts,nLev), mr_ccliq(nPts,nLev), mr_ccice(nPts,nLev),                    &
	   fl_lsrain(nPts,nLev), fl_lssnow(nPts,nLev), fl_lsgrpl(nPts,nLev),                 &
       fl_ccrain(nPts,nLev), fl_ccsnow(nPts,nLev), dtau_s(nPts,nLev),                    &
       dtau_c(nPts,nLev), dem_s(nPts,nLev), dem_c(nPts,nLev), dtau_s_snow(nPts,nLev),    &
       dem_s_snow(nPts,nLev), mr_ozone(nPts,nLev))
    ! 3D
    allocate(Reff(nPts,nLev,9))
	
	! Read in data
	! 1D
	allocate(lon1D(nLon),lat1D(nLat))
	! lon
	err = nf90_inq_varid(ncid, 'lon', varID)
    err = nf90_get_var(ncid, varID, lon1D)
    lon = pack(spread(lon1D,1,nLat),.true.)
    ! lat
	err = nf90_inq_varid(ncid, 'lat', varID)
    err = nf90_get_var(ncid, varID, lat1D)    
	lat = pack(spread(lat1D,1,nLon),.true.)
	! 2D
	allocate(var2dtemp(nLon,nLat))
	! sunlit
	err = nf90_inq_varid(ncid, 'lsmask', varID)
    err = nf90_get_var(ncid, varID, var2dtemp)
    sunlit = reshape(var2dtemp, (/ nPts/))	
	! TS_COSP
	err = nf90_inq_varid(ncid, 'TS_COSP', varID)
    err = nf90_get_var(ncid, varID, var2dtemp)
    skt = reshape(var2dtemp, (/ nPts/))		
    ! LANDFRAC
	err = nf90_inq_varid(ncid, 'LANDFRAC', varID)
    err = nf90_get_var(ncid, varID, var2dtemp)
    landmask = nint(reshape(var2dtemp, (/ nPts/)))	
        
	! 3D
	allocate(var3dtemp(nLon,nLat,nLev))	
	! P_COSP
	err = nf90_inq_varid(ncid, 'P_COSP', varID)
    err = nf90_get_var(ncid, varID, var3dtemp)
    p   = reshape(var3dtemp, (/ nPts, nLev /))
	! PH_COSP	
	err = nf90_inq_varid(ncid, 'PH_COSP', varID)
    err = nf90_get_var(ncid, varID, var3dtemp)
    ph   = reshape(var3dtemp, (/ nPts, nLev /))
	! ZLEV_COSP
	err = nf90_inq_varid(ncid, 'ZLEV_COSP', varID)
    err = nf90_get_var(ncid, varID, var3dtemp)
    z   = reshape(var3dtemp, (/ nPts, nLev /))
	! ZLEV_HALF_COSP	
	err = nf90_inq_varid(ncid, 'ZLEV_HALF_COSP', varID)
    err = nf90_get_var(ncid, varID, var3dtemp)
    zh   = reshape(var3dtemp, (/ nPts, nLev /))
	! T_COSP
	err = nf90_inq_varid(ncid, 'T_COSP', varID)
    err = nf90_get_var(ncid, varID, var3dtemp)
    T   = reshape(var3dtemp, (/ nPts, nLev /))        
	! Q_COSP
	err = nf90_inq_varid(ncid, 'Q_COSP', varID)
    err = nf90_get_var(ncid, varID, var3dtemp)
    qv  = reshape(var3dtemp, (/ nPts, nLev /))
	! CLD_COSP
	err = nf90_inq_varid(ncid, 'CLD_COSP', varID)
    err = nf90_get_var(ncid, varID, var3dtemp)
    tca = reshape(var3dtemp, (/ nPts, nLev /))        
	! CONCLD_COSP
	err = nf90_inq_varid(ncid, 'CONCLD_COSP', varID)
    err = nf90_get_var(ncid, varID, var3dtemp)
    cca = reshape(var3dtemp, (/ nPts, nLev /)) 
    ! LSCLIQ_COSP
	err = nf90_inq_varid(ncid, 'LSCLIQ_COSP', varID)
    err = nf90_get_var(ncid, varID, var3dtemp)
    mr_lsliq = reshape(var3dtemp, (/ nPts, nLev /))     
    ! LSCICE_COSP
	err = nf90_inq_varid(ncid, 'LSCICE_COSP', varID)
    err = nf90_get_var(ncid, varID, var3dtemp)
    mr_lsice = reshape(var3dtemp, (/ nPts, nLev /))  
    ! CVCLIQ_COSP
	err = nf90_inq_varid(ncid, 'CVCLIQ_COSP', varID)
    err = nf90_get_var(ncid, varID, var3dtemp)
    mr_ccliq = reshape(var3dtemp, (/ nPts, nLev /))     
    ! CVCICE_COSP
	err = nf90_inq_varid(ncid, 'CVCICE_COSP', varID)
    err = nf90_get_var(ncid, varID, var3dtemp)
    mr_ccice = reshape(var3dtemp, (/ nPts, nLev /))   
    ! RAIN_LS_COSP
	err = nf90_inq_varid(ncid, 'RAIN_LS_COSP', varID)
    err = nf90_get_var(ncid, varID, var3dtemp)
    fl_lsrain = reshape(var3dtemp, (/ nPts, nLev /))       
    ! SNOW_LS_COSP
	err = nf90_inq_varid(ncid, 'SNOW_LS_COSP', varID)
    err = nf90_get_var(ncid, varID, var3dtemp)
    fl_lssnow = reshape(var3dtemp, (/ nPts, nLev /))      
    ! GFPL_LS_COSP
	err = nf90_inq_varid(ncid, 'GRPL_LS_COSP', varID)
    err = nf90_get_var(ncid, varID, var3dtemp)
    fl_lsgrpl = reshape(var3dtemp, (/ nPts, nLev /))      
    ! RAIN_CV_COSP
	err = nf90_inq_varid(ncid, 'RAIN_CV_COSP', varID)
    err = nf90_get_var(ncid, varID, var3dtemp)
    fl_ccrain = reshape(var3dtemp, (/ nPts, nLev /))       
    ! SNOW_CV_COSP
	err = nf90_inq_varid(ncid, 'SNOW_CV_COSP', varID)
    err = nf90_get_var(ncid, varID, var3dtemp)
    fl_ccsnow = reshape(var3dtemp, (/ nPts, nLev /))       
    ! DTAU_S_COSP
	err = nf90_inq_varid(ncid, 'DTAU_S_COSP', varID)
    err = nf90_get_var(ncid, varID, var3dtemp)
    dtau_s = reshape(var3dtemp, (/ nPts, nLev /))     
    ! DTAU_C_COSP
	err = nf90_inq_varid(ncid, 'DTAU_C_COSP', varID)
    err = nf90_get_var(ncid, varID, var3dtemp)
    dtau_c = reshape(var3dtemp, (/ nPts, nLev /)) 
    ! DEM_S_COSP
	err = nf90_inq_varid(ncid, 'DEM_S_COSP', varID)
    err = nf90_get_var(ncid, varID, var3dtemp)
    dem_s = reshape(var3dtemp, (/ nPts, nLev /))     
    ! DEM_C_COSP
	err = nf90_inq_varid(ncid, 'DEM_C_COSP', varID)
    err = nf90_get_var(ncid, varID, var3dtemp)
    dem_c = reshape(var3dtemp, (/ nPts, nLev /))     
	! O3_COSP
	err = nf90_inq_varid(ncid, 'O3_COSP', varID)
    err = nf90_get_var(ncid, varID, var3dtemp)
    mr_ozone = reshape(var3dtemp, (/ nPts, nLev /))  
	! 4D
	! REFF_COSP
	err = nf90_inq_varid(ncid, 'REFF_COSP_1', varID)
    err = nf90_get_var(ncid, varID, var3dtemp)
    Reff(:,:,1) = reshape(var3dtemp, (/ nPts, nLev /))          
	err = nf90_inq_varid(ncid, 'REFF_COSP_2', varID)
    err = nf90_get_var(ncid, varID, var3dtemp)
    Reff(:,:,2) = reshape(var3dtemp, (/ nPts, nLev /)) 
	err = nf90_inq_varid(ncid, 'REFF_COSP_3', varID)
    err = nf90_get_var(ncid, varID, var3dtemp)
    Reff(:,:,3) = reshape(var3dtemp, (/ nPts, nLev /)) 
	err = nf90_inq_varid(ncid, 'REFF_COSP_4', varID)
    err = nf90_get_var(ncid, varID, var3dtemp)
    Reff(:,:,4) = reshape(var3dtemp, (/ nPts, nLev /)) 
	err = nf90_inq_varid(ncid, 'REFF_COSP_5', varID)
    err = nf90_get_var(ncid, varID, var3dtemp)
    Reff(:,:,5) = reshape(var3dtemp, (/ nPts, nLev /)) 
	err = nf90_inq_varid(ncid, 'REFF_COSP_6', varID)
    err = nf90_get_var(ncid, varID, var3dtemp)
    Reff(:,:,6) = reshape(var3dtemp, (/ nPts, nLev /))     
	err = nf90_inq_varid(ncid, 'REFF_COSP_7', varID)
    err = nf90_get_var(ncid, varID, var3dtemp)
    Reff(:,:,7) = reshape(var3dtemp, (/ nPts, nLev /)) 
	err = nf90_inq_varid(ncid, 'REFF_COSP_8', varID)
    err = nf90_get_var(ncid, varID, var3dtemp)
    Reff(:,:,8) = reshape(var3dtemp, (/ nPts, nLev /)) 
	err = nf90_inq_varid(ncid, 'REFF_COSP_9', varID)
    err = nf90_get_var(ncid, varID, var3dtemp)
    Reff(:,:,9) = reshape(var3dtemp, (/ nPts, nLev /))             
    ! DTAU_S_COSP_SNOW
	err = nf90_inq_varid(ncid, 'DTAU_S_COSP_SNOW', varID)
    err = nf90_get_var(ncid, varID, var3dtemp)
    dtau_s_snow = reshape(var3dtemp, (/ nPts, nLev /))    
    ! DEM_S_COSP_SNOW
	err = nf90_inq_varid(ncid, 'DEM_S_COSP_SNOW', varID)
    err = nf90_get_var(ncid, varID, var3dtemp)
    dem_s_snow = reshape(var3dtemp, (/ nPts, nLev /))      
  end subroutine read_cesm2_input

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
  end module mod_cosp2_cesm2_io
