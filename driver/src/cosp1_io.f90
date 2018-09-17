module mod_cosp1_io
  use cosp_kinds, only: wp
  use mod_cosp,   only: cosp_outputs
  use netcdf
  USE MOD_COSP_CONFIG, ONLY:  Nlvgrid, LIDAR_NCAT, SR_BINS, PARASOL_NREFL, cloudsat_DBZE_BINS, &
       numMODISReffIceBins, numMODISReffLiqBins, ntau, tau_binBounds, tau_binCenters, &
       tau_binEdges,npres, pres_binBounds, pres_binCenters, pres_binEdges, nhgt,      &
       hgt_binBounds, hgt_binCenters, hgt_binEdges, reffLIQ_binCenters,vgrid_z,       &
       reffICE_binCenters, reffLIQ_binCenters, cloudsat_binCenters, PARASOL_SZA,      &
       calipso_binCenters
  USE MOD_COSP_INTERFACE_v1p4, ONLY: cosp_gridbox, cosp_config, cosp_subgrid,         &
       cosp_sglidar, cosp_lidarstats, cosp_isccp, cosp_misr, cosp_rttov, cosp_sgradar,&
       cosp_radarstats, cosp_modis, cosp_vgrid

  implicit none

contains

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE write_cosp1_output
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine write_cosp1_output(Npoints, Ncolumns, Nlevels, lev, lon, lat, cfg, vgrid,   &
       gbx, sgx, sgradar, sglidar, isccp, misr, modis, rttov, stradar, stlidar,  outFileName)
    
    integer,intent(in) :: Npoints, Ncolumns, Nlevels
    real(wp),dimension(Npoints),intent(in) :: lon,lat
    real(wp),dimension(Nlevels),intent(in) :: lev
    type(cosp_config),     intent(in) :: cfg      ! Configuration options
    type(cosp_vgrid),      intent(in) :: vgrid    ! Information on vertical grid of stats
    type(cosp_subgrid),    intent(in) :: sgx      ! Subgrid info
    type(cosp_sgradar),    intent(in) :: sgradar  ! Output from radar simulator (pixel)
    type(cosp_sglidar),    intent(in) :: sglidar  ! Output from lidar simulator (pixel)
    type(cosp_isccp),      intent(in) :: isccp    ! Output from ISCCP simulator
    type(cosp_misr),       intent(in) :: misr     ! Output from MISR simulator
    type(cosp_modis),      intent(in) :: modis    ! Output from MODIS simulator
    type(cosp_rttov),      intent(in) :: rttov    ! Output from RTTOV
    type(cosp_radarstats), intent(in) :: stradar  ! Summary statistics from cloudsatsimulator (gridbox)
    type(cosp_lidarstats), intent(in) :: stlidar  ! Output from LIDAR simulator (gridbox)
    type(cosp_gridbox),    intent(in) :: gbx      ! COSP gridbox type from v1.4


    character(len=256),intent(in) :: outFileName

    integer :: fileID,status,ij
    integer,dimension(20)  :: dimID
    integer,dimension(100) :: varID
    integer,dimension(Npoints) :: loc
    integer,dimension(Ncolumns) :: cosp_scol
    integer,dimension(2) :: bnds
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
    status = nf90_def_dim(fileID,"loc",Npoints,dimID(1))
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

    ! ---------------------------------------------------------------------------------------
    ! Define varaibles
    ! ---------------------------------------------------------------------------------------
    ! Longitude
    status = nf90_def_var(fileID,"longitude",  nf90_float, (/dimID(1)/),varID(1))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(1),"long_name","longitude")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(1),"units",        "degrees_east")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(1),"standard_name", "longitude")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    ! Latitude
    status = nf90_def_var(fileID,"latitude",   nf90_float, (/dimID(1)/),varID(2))
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
    status = nf90_def_var(fileID,"lev",  nf90_float, (/dimID(3)/),varID(86))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(86),"long_name","level indices")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(86),"units",        "1")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    ! Levels for statistical diagnostics (lidar and radar)
    status = nf90_def_var(fileID,"levStat",  nf90_float, (/dimID(4)/),varID(87))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(87),"long_name","level indices")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(87),"units",        "1")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    
    ! Subcolumms
    status = nf90_def_var(fileID,"cosp_scol",  nf90_float, (/dimID(2)/),varID(88))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(88),"long_name","subcolumn indices")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(88),"units",        "1")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    ! Bnds
    status = nf90_def_var(fileID,"bnds",  nf90_float, (/dimID(6)/),varID(89))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(89),"long_name","bounds")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(89),"units",        "1")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    ! loc
    status = nf90_def_var(fileID,"loc",  nf90_float, (/dimID(1)/),varID(90))
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(90),"long_name","loc")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_att(fileID,varID(90),"units",        "1")
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    
    ! CALIPSO simulator output
    if (cfg%Latb532) then
       status = nf90_def_var(fileID,"atb532",nf90_float, (/dimID(1),dimID(2),dimID(3)/),varID(10))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(10),"long_name","CALIPSO Attenuated Total Backscatter (532nm)")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(10),"units",        "m-1 sr-1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))      
       status = nf90_put_att(fileID,varID(10),"standard_name", "volume_attenuated_backwards_scattering_function_in_air")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))

       status = nf90_def_var(fileID,"atb532_perp",nf90_float, (/dimID(1),dimID(2),dimID(3)/),varID(9))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(9),"long_name","CALIPSO Attenuated Total Perpendicular Backscatter (532nm)")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(9),"units",        "m-1 sr-1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))      
       status = nf90_put_att(fileID,varID(9),"standard_name", "volume_attenuated_backwards_scattering_function_in_air")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lclcalipsoice) then
       status = nf90_def_var(fileID,"clcalipsoice",nf90_float, (/dimID(1),dimID(4)/),varID(58))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(58),"long_name","CALIPSO Cloud Fraction Undetected by CloudSat")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(58),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(58),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lclcalipsoliq) then
       status = nf90_def_var(fileID,"clcalipsoliq",nf90_float, (/dimID(1),dimID(4)/),varID(59))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(59),"long_name","CALIPSO Cloud Fraction Undetected by CloudSat")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(59),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(59),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lclcalipsoun) then
       status = nf90_def_var(fileID,"clcalipsoun",nf90_float, (/dimID(1),dimID(4)/),varID(60))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(60),"long_name","CALIPSO Cloud Fraction Undetected by CloudSat")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(60),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(60),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))       
    endif
    if (cfg%Lcllcalipsoice) then
       status = nf90_def_var(fileID,"cllcalipsoice",nf90_float, (/dimID(1)/),varID(61))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(61),"long_name","CALIPSO Ice Low Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(61),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(61),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))     
    endif
    if (cfg%Lclmcalipsoice) then
       status = nf90_def_var(fileID,"clmcalipsoice",nf90_float, (/dimID(1)/),varID(62))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(62),"long_name","CALIPSO Ice Mid Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(62),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(62),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lclhcalipsoice) then
       status = nf90_def_var(fileID,"clhcalipsoice",nf90_float, (/dimID(1)/),varID(63))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(63),"long_name","CALIPSO Ice High Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(63),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(63),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lcltcalipsoice) then
       status = nf90_def_var(fileID,"cltcalipsoice",nf90_float, (/dimID(1)/),varID(64))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(64),"long_name","CALIPSO Ice Total Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(64),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(64),"standard_name", "cloud_area_fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lcllcalipsoliq) then
       status = nf90_def_var(fileID,"cllcalipsoliq",nf90_float, (/dimID(1)/),varID(65))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(65),"long_name","CALIPSO Liquid Low Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(65),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(65),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))     
    endif
    if (cfg%Lclmcalipsoliq) then
       status = nf90_def_var(fileID,"clmcalipsoliq",nf90_float, (/dimID(1)/),varID(66))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(66),"long_name","CALIPSO Liquid Mid Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(66),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(66),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))     
    endif
    if (cfg%Lclhcalipsoliq) then
       status = nf90_def_var(fileID,"clhcalipsoliq",nf90_float, (/dimID(1)/),varID(67))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(67),"long_name","CALIPSO Liquid High Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(67),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(67),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))     
    endif
    if (cfg%Lcltcalipsoliq) then
       status = nf90_def_var(fileID,"cltcalipsoliq",nf90_float, (/dimID(1)/),varID(68))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(68),"long_name","CALIPSO Liquid Total Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(68),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(68),"standard_name", "cloud_area_fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))            
    endif
    if (cfg%Lcllcalipsoun) then
       status = nf90_def_var(fileID,"cllcalipsoun",nf90_float, (/dimID(1)/),varID(69))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(69),"long_name","CALIPSO Undefined-Phase Low Level Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(69),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(69),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))     
    endif
    if (cfg%Lclmcalipsoun) then
       status = nf90_def_var(fileID,"clmcalipsoun",nf90_float, (/dimID(1)/),varID(70))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(70),"long_name","CALIPSO Undefined-Phase Mid Level Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(70),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(70),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))     
    endif
    if (cfg%Lclhcalipsoun) then
       status = nf90_def_var(fileID,"clhcalipsoun",nf90_float, (/dimID(1)/),varID(71))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(71),"long_name","CALIPSO Undefined-Phase High Level Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(71),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(71),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))     
    endif
    if (cfg%Lcltcalipsoun) then
       status = nf90_def_var(fileID,"cltcalipsoun",nf90_float, (/dimID(1)/),varID(72))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(72),"long_name","CALIPSO Undefined-Phase Total Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(72),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(72),"standard_name", "cloud_area_fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))            
    endif
    if (cfg%Lclcalipsotmp) then
       status = nf90_def_var(fileID,"clcalipsotmp",nf90_float, (/dimID(1),dimID(4)/),varID(77))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(77),"long_name","CALIPSO Cloud Fraction Undetected by CloudSat")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(77),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(77),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))     
    endif
    if (cfg%Lclcalipsotmpice) then
       status = nf90_def_var(fileID,"clcalipsotmpice",nf90_float, (/dimID(1),dimID(4)/),varID(78))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(78),"long_name","CALIPSO Cloud Fraction Undetected by CloudSat")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(78),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(78),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))     
    endif
    if (cfg%Lclcalipsotmpliq) then
       status = nf90_def_var(fileID,"clcalipsotmpliq",nf90_float, (/dimID(1),dimID(4)/),varID(79))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(79),"long_name","CALIPSO Cloud Fraction Undetected by CloudSat")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(79),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(79),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))            
    endif
    if (cfg%Lclcalipsotmpun) then
       status = nf90_def_var(fileID,"clcalipsotmpun",nf90_float, (/dimID(1),dimID(4)/),varID(80))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(80),"long_name","CALIPSO Cloud Fraction Undetected by CloudSat")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(80),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(80),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))            
    endif
    if (cfg%LcfadLidarsr532) then
       status = nf90_def_var(fileID,"cfadLidarsr532",nf90_float, (/dimID(1),dimID(12),dimID(4)/),varID(15))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(15),"long_name","CALIPSO Scattering Ratio CFAD")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(15),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(15),"standard_name", "histogram_of_backscattering_ratio_over_height_above_reference_ellipsoid")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))

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
    if (cfg%Lclcalipso) then
       status = nf90_def_var(fileID,"clcalipso",nf90_float, (/dimID(1),dimID(4)/),varID(16))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(16),"long_name","CALIPSO Cloud Area Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(16),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(16),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))     
    endif
    if (cfg%Lcllcalipso) then
       ! Low-level
       status = nf90_def_var(fileID,"cllcalipso",nf90_float, (/dimID(1)/),varID(73))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(73),"long_name","CALIPSO Low Level Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(73),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(73),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))           
    endif
    if (cfg%Lclmcalipso) then
       ! Mid-level
       status = nf90_def_var(fileID,"clmcalipso",nf90_float, (/dimID(1)/),varID(74))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(74),"long_name","CALIPSO Mid Level Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(74),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(74),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))    
    endif
    if (cfg%Lclhcalipso) then
       ! High-level
       status = nf90_def_var(fileID,"clhcalipso",nf90_float, (/dimID(1)/),varID(75))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(75),"long_name","CALIPSO High Level Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(75),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(75),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lcltcalipso) then
       ! Total
       status = nf90_def_var(fileID,"cltcalipso",nf90_float, (/dimID(1)/),varID(76))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(76),"long_name","CALIPSO Total Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(76),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(76),"standard_name", "cloud_area_fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))           
    endif
    if (cfg%LlidarBetaMol532) then
       status = nf90_def_var(fileID,"lidarBetaMol532",nf90_float, (/dimID(1),dimID(3)/),varID(18))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(18),"long_name","CALIPSO Molecular Backscatter Coefficient (532nm)")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(18),"units",        "m-1 sr-1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(18),"standard_name", "volume_attenuated_backwards_scattering_function_in_air_assuming_no_aerosol_or_cloud")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))           
    endif

    ! PARASOL simulator output
    if (cfg%LparasolRefl) then
       status = nf90_def_var(fileID,"parasolPix_refl",nf90_float, (/dimID(1),dimID(2),dimID(13)/),varID(20))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(20),"long_name","PARASOL Subcolumn Reflectance")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(20),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(20),"standard_name", "toa_bidirectional_reflectance")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))          

       status = nf90_def_var(fileID,"parasolGrid_refl",nf90_float, (/dimID(1),dimID(13)/),varID(21))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(21),"long_name","PARASOL Reflectance")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(21),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(21),"standard_name", "toa_bidirectional_reflectance")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))       

       status = nf90_def_var(fileID,"PARASOL_NREFL",nf90_float, (/dimID(13)/),varID(82))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(82),"long_name","PARASOL Solar Zenith Angle")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(82),"units",        "degree")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(82),"standard_name", "solar_zenith_angle")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    
    ! Cloudsat simulator output
    if (cfg%Ldbze94) then
       status = nf90_def_var(fileID,"dbze94",nf90_float, (/dimID(1),dimID(2),dimID(3)/),varID(22))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(22),"long_name","CloudSat Radar Reflectivity")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(22),"units",        "dBZ")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(22),"standard_name", "equivalent_reflectivity_factor")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))       
    endif
    if (cfg%LcfadDbze94) then
       status = nf90_def_var(fileID,"cfadDbze94",nf90_float, (/dimID(1),dimID(14),dimID(4)/),varID(23))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(23),"long_name","CloudSat Radar reflectivity CFAD")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(23),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(23),"standard_name", "histogram_of_equivalent_reflectivity_factor_over_height_above_reference_ellipsoid")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status)) 
       status = nf90_def_var(fileID,"cloudsat_DBZE_BINS",nf90_float, (/dimID(14)/),varID(83))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(83),"long_name","CloudSat simulator equivalent radar reflectivity factor")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(83),"units",        "dBZ")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(83),"standard_name", "equivalent_reflectivity_factor")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))              
    endif
    
    ! ISCCP simulator outputs
    if (cfg%Lcltisccp) then
       status = nf90_def_var(fileID,"cltisccp",nf90_float, (/dimID(1)/),varID(24))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(24),"long_name","ISCCP Total Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(24),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(24),"standard_name", "cloud_area_fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))              
    endif
    if (cfg%Lmeantbisccp) then
       status = nf90_def_var(fileID,"meantbisccp",nf90_float, (/dimID(1)/),varID(25))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(25),"long_name","ISCCP all-sky 10.5 micron brightness temperature")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(25),"units",        "K")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(25),"standard_name", "toa_brightness_temperature")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))   
    endif
    if (cfg%Lmeantbclrisccp) then
       status = nf90_def_var(fileID,"meantbclrisccp",nf90_float, (/dimID(1)/),varID(26))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(26),"long_name","ISCCP clear-sky 10.5 micron brightness temperature")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(26),"units",        "K")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(26),"standard_name", "toa_brightness_temperature_assuming_clear_sky")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))   
    endif
    if (cfg%Lpctisccp) then
       status = nf90_def_var(fileID,"pctisccp",nf90_float, (/dimID(1)/),varID(27))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(27),"long_name","ISCCP Mean Cloud Top Pressure")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(27),"units",        "hPa")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(27),"standard_name", "air_pressure_at_cloud_top")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))         
    endif
    if (cfg%Ltauisccp) then
       status = nf90_def_var(fileID,"tauisccp",nf90_float, (/dimID(1)/),varID(28))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(28),"long_name","ISCCP Mean Optical Depth")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(28),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(28),"standard_name", "atmosphere_optical_thickness_due_to_cloud")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))         
    endif
    if (cfg%Lalbisccp) then
       status = nf90_def_var(fileID,"albisccp",nf90_float, (/dimID(1)/),varID(29))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(29),"long_name","ISCCP Mean Cloud Albedo")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(29),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(29),"standard_name", "cloud_albedo")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
    endif
    if (cfg%Lboxtauisccp) then
       status = nf90_def_var(fileID,"boxtauisccp",nf90_float, (/dimID(1),dimID(2)/),varID(30))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(30),"long_name","ISCCP Subcolumn Optical Depth")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(30),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(30),"standard_name", "atmosphere_optical_thickness_due_to_cloud")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))         
    endif
    if (cfg%Lboxptopisccp) then
       status = nf90_def_var(fileID,"boxptopisccp",nf90_float, (/dimID(1),dimID(2)/),varID(31))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(31),"long_name","ISCCP Subcolumn Cloud Top Pressure")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(31),"units",        "Pa")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(31),"standard_name", "air_pressure_at_cloud_top")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))         
    endif
    if (cfg%Lclisccp) then
       status = nf90_def_var(fileID,"clisccp",nf90_float, (/dimID(1),dimID(5),dimID(7)/),varID(32))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(32),"long_name","ISCCP joint-PDF of cloud top pressure and optical depth")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(32),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(32),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))         
    endif
    
    ! MISR simulator output
    if (cfg%LclMISR) then
       status = nf90_def_var(fileID,"clMISR",nf90_float, (/dimID(1),dimID(5),dimID(8)/),varID(33))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(33),"long_name","MISR joint-PDF of cloud top pressure and optical depth")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(33),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(33),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  

       status = nf90_def_var(fileID,"misr_meanztop",nf90_float, (/dimID(1)/),varID(34))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(34),"long_name","MISR Mean Cloud Top Height")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(34),"units",        "m")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(34),"standard_name", "cloud_top_altitude")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
	    
       status = nf90_def_var(fileID,"misr_cldarea",nf90_float, (/dimID(1)/),varID(35))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(35),"long_name","MISR cloud cover")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(35),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(35),"standard_name", "cloud_area_fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))        
    endif

    ! MODIS simulator output
    if (cfg%Lcltmodis) then
       status = nf90_def_var(fileID,"cltmodis",nf90_float, (/dimID(1)/),varID(36))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(36),"long_name","MODIS Total Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(36),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(36),"standard_name", "cloud_area_fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lclwmodis) then
       status = nf90_def_var(fileID,"clwmodis",nf90_float, (/dimID(1)/),varID(37))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(37),"long_name","MODIS Liquid Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(37),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(37),"standard_name", "cloud_area_fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lclimodis) then
       status = nf90_def_var(fileID,"climodis",nf90_float, (/dimID(1)/),varID(38))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(38),"long_name","MODIS Ice Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(38),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(38),"standard_name", "cloud_area_fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lclhmodis) then
       status = nf90_def_var(fileID,"clhmodis",nf90_float, (/dimID(1)/),varID(39))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(39),"long_name","MODIS High Level Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(39),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(39),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lclmmodis) then
       status = nf90_def_var(fileID,"clmmodis",nf90_float, (/dimID(1)/),varID(40))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(40),"long_name","MODIS Mid Level Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(40),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(40),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lcllmodis) then
       status = nf90_def_var(fileID,"cllmodis",nf90_float, (/dimID(1)/),varID(41))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(41),"long_name","MODIS Low Level Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(41),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(41),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Ltautmodis) then
       status = nf90_def_var(fileID,"tautmodis",nf90_float, (/dimID(1)/),varID(42))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(42),"long_name","MODIS Total Cloud Optical Thickness")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(42),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(42),"standard_name", "atmosphere_optical_thickness_due_to_cloud")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Ltauwmodis) then
       status = nf90_def_var(fileID,"tauwmodis",nf90_float, (/dimID(1)/),varID(43))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(43),"long_name","MODIS Liquid Cloud Optical Thickness")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(43),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(43),"standard_name", "atmosphere_optical_thickness_due_to_cloud")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Ltauimodis) then
       status = nf90_def_var(fileID,"tauimodis",nf90_float, (/dimID(1)/),varID(44))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(44),"long_name","MODIS Ice Cloud Optical Thickness")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(44),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(44),"standard_name", "atmosphere_optical_thickness_due_to_cloud")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status)) 
    endif
    if (cfg%Ltautlogmodis) then
       status = nf90_def_var(fileID,"tautlogmodis",nf90_float, (/dimID(1)/),varID(45))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(45),"long_name","MODIS Total Cloud Optical Thickness (Log10 Mean)")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(45),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(45),"standard_name", "atmosphere_optical_thickness_due_to_cloud")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))        
    endif
    if (cfg%Ltauwlogmodis) then
       status = nf90_def_var(fileID,"tauwlogmodis",nf90_float, (/dimID(1)/),varID(46))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(46),"long_name","MODIS Liquid Cloud Optical Thickness (Log10 Mean)")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(46),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(46),"standard_name", "atmosphere_optical_thickness_due_to_cloud")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status)) 
    endif
    if (cfg%Ltauilogmodis) then
       status = nf90_def_var(fileID,"tauilogmodis",nf90_float, (/dimID(1)/),varID(47))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(47),"long_name","MODIS Ice Cloud Optical Thickness (Log10 Mean)")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(47),"units",        "1")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(47),"standard_name", "atmosphere_optical_thickness_due_to_cloud")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))        
    endif
    if (cfg%Lreffclwmodis) then
       status = nf90_def_var(fileID,"reffclwmodis",nf90_float, (/dimID(1)/),varID(48))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(48),"long_name","MODIS Liquid Cloud Particle Size")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(48),"units",        "m")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
       status = nf90_put_att(fileID,varID(48),"standard_name", "effective_radius_of_cloud_liquid_water_particle")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status)) 
    endif
    if (cfg%Lreffclimodis) then
       status = nf90_def_var(fileID,"reffclimodis",nf90_float, (/dimID(1)/),varID(49))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(49),"long_name","MODIS Ice Cloud Particle Size")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(49),"units",        "m")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status)) 
       status = nf90_put_att(fileID,varID(49),"standard_name", "effective_radius_of_cloud_liquid_water_particle")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status)) 
    endif
    if (cfg%Lpctmodis) then
       status = nf90_def_var(fileID,"pctmodis",nf90_float, (/dimID(1)/),varID(50))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(50),"long_name","MODIS Cloud Top Pressure")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(50),"units",        "hPa")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status)) 
       status = nf90_put_att(fileID,varID(50),"standard_name", "air_pressure_at_cloud_top")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Llwpmodis) then
       status = nf90_def_var(fileID,"lwpmodis",nf90_float, (/dimID(1)/),varID(51))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(51),"long_name","MODIS Cloud Liquid Water Path")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(51),"units",        "kg m-2")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status)) 
       status = nf90_put_att(fileID,varID(51),"standard_name", "atmosphere_cloud_liquid_water_content")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Liwpmodis) then
       status = nf90_def_var(fileID,"iwpmodis",nf90_float, (/dimID(1)/),varID(52))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(52),"long_name","MODIS Cloud Ice Water Path")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(52),"units",        "kg m-2")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status)) 
       status = nf90_put_att(fileID,varID(52),"standard_name", "atmosphere_mass_content_of_cloud_ice")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lclmodis) then
       status = nf90_def_var(fileID,"clmodis",nf90_float, (/dimID(1),dimID(5),dimID(7)/),varID(53))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(53),"long_name","MODIS joint-PDF of cloud top pressure and optical depth")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(53),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status)) 
       status = nf90_put_att(fileID,varID(53),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status)) 
    endif
    ! Joint simulator products.
    if (cfg%Lclcalipso2) then
       status = nf90_def_var(fileID,"clcalipso2",nf90_float, (/dimID(1),dimID(4)/),varID(56))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(56),"long_name","CALIPSO Cloud Fraction Undetected by CloudSat")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(56),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(56),"standard_name", "cloud_area_fraction_in_atmosphere_layer")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status)) 
    endif
    if (cfg%Lcltlidarradar) then
       status = nf90_def_var(fileID,"cltlidarradar",nf90_float, (/dimID(1)/),varID(57))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(57),"long_name","CALIPSO and CloudSat Total Cloud Fraction")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_att(fileID,varID(57),"units",        "%")
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status)) 
       status = nf90_put_att(fileID,varID(57),"standard_name", "cloud_area_fraction")
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
    status = nf90_put_var(fileID,varID(1),lon)
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_var(fileID,varID(2),lat)
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
    status = nf90_put_var(fileID,varID(86),lev)
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_var(fileID,varID(87),vgrid_z)
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_var(fileID,varID(88),cosp_scol)
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_var(fileID,varID(89),bnds)
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    status = nf90_put_var(fileID,varID(90),loc)
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    
    ! CALIPSO simulator output
    if (cfg%Latb532) then
       status = nf90_put_var(fileID,varID(9),sglidar%betaperp_tot)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))

       status = nf90_put_var(fileID,varID(10),sglidar%beta_tot)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lclcalipsoice) then
       status = nf90_put_var(fileID,varID(58), stlidar%lidarcldphase(:,:,1))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lclcalipsoliq) then
       status = nf90_put_var(fileID,varID(59), stlidar%lidarcldphase(:,:,2))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lclcalipsoun) then    
       status = nf90_put_var(fileID,varID(60), stlidar%lidarcldphase(:,:,3))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lcllcalipsoice) then
       status = nf90_put_var(fileID,varID(61), stlidar%cldlayerphase(:,1,1))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lclmcalipsoice) then
       status = nf90_put_var(fileID,varID(62), stlidar%cldlayerphase(:,2,1))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lclhcalipsoice) then
       status = nf90_put_var(fileID,varID(63), stlidar%cldlayerphase(:,3,1))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lcltcalipsoice) then
       status = nf90_put_var(fileID,varID(64), stlidar%cldlayerphase(:,4,1))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lcllcalipsoliq) then       
       status = nf90_put_var(fileID,varID(65), stlidar%cldlayerphase(:,1,2))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lclmcalipsoliq) then
       status = nf90_put_var(fileID,varID(66), stlidar%cldlayerphase(:,2,2))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lclhcalipsoliq) then
       status = nf90_put_var(fileID,varID(67), stlidar%cldlayerphase(:,3,2))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lcltcalipsoliq) then
       status = nf90_put_var(fileID,varID(68), stlidar%cldlayerphase(:,4,2))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lcllcalipsoun) then       
       status = nf90_put_var(fileID,varID(69), stlidar%cldlayerphase(:,1,3))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lclmcalipsoun) then
       status = nf90_put_var(fileID,varID(70), stlidar%cldlayerphase(:,2,3))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lclhcalipsoun) then
       status = nf90_put_var(fileID,varID(71), stlidar%cldlayerphase(:,3,3))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lcltcalipsoun) then
       status = nf90_put_var(fileID,varID(72),stlidar%cldlayerphase(:,4,3))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lclcalipsotmp) then
       status = nf90_put_var(fileID,varID(77), stlidar%lidarcldtmp(:,:,1))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lclcalipsotmpice) then
       status = nf90_put_var(fileID,varID(78), stlidar%lidarcldtmp(:,:,2))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lclcalipsotmpliq) then
       status = nf90_put_var(fileID,varID(79), stlidar%lidarcldtmp(:,:,3))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lclcalipsotmpun) then
       status = nf90_put_var(fileID,varID(80), stlidar%lidarcldtmp(:,:,4))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%LcfadLidarsr532) then
       status = nf90_put_var(fileID,varID(15),stlidar%cfad_sr)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       
       status = nf90_put_var(fileID,varID(19),reshape([stlidar%srbval(1:SR_BINS), stlidar%srbval(2:SR_BINS+1)],(/2,SR_BINS/)))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       
       status = nf90_put_var(fileID,varID(81),calipso_binCenters)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lclcalipso) then
       status = nf90_put_var(fileID,varID(16),stlidar%lidarcld)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lcllcalipso) then
       status = nf90_put_var(fileID,varID(73),stlidar%cldlayer(:,1))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lclmcalipso) then
       status = nf90_put_var(fileID,varID(74),stlidar%cldlayer(:,2))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lclhcalipso) then
       status = nf90_put_var(fileID,varID(75),stlidar%cldlayer(:,3))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lcltcalipso) then
       status = nf90_put_var(fileID,varID(76),stlidar%cldlayer(:,4))
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%LlidarBetaMol532) then
       status = nf90_put_var(fileID,varID(18),sglidar%beta_mol)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    
    ! PARASOL simulator output
    if (cfg%LparasolRefl) then
       status = nf90_put_var(fileID,varID(20),sglidar%refl)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))

       status = nf90_put_var(fileID,varID(21),stlidar%parasolrefl)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))

       status = nf90_put_var(fileID,varID(82),PARASOL_SZA)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    
    ! Cloudsat simulator output
    if (cfg%Ldbze94) then
       status = nf90_put_var(fileID,varID(22),sgradar%Ze_tot)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%LcfadDbze94) then
       status = nf90_put_var(fileID,varID(23),stradar%cfad_ze)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
       status = nf90_put_var(fileID,varID(83),cloudsat_binCenters)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif

    if (cfg%Lcltisccp) then
       status = nf90_put_var(fileID,varID(24),isccp%totalcldarea)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lmeantbisccp) then
       status = nf90_put_var(fileID,varID(25),isccp%meantb)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lmeantbclrisccp) then
       status = nf90_put_var(fileID,varID(26),isccp%meantbclr)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lpctisccp) then
       status = nf90_put_var(fileID,varID(27),isccp%meanptop)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Ltauisccp) then
       status = nf90_put_var(fileID,varID(28),isccp%meantaucld)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lalbisccp) then
       status = nf90_put_var(fileID,varID(29),isccp%meanalbedocld)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lboxtauisccp) then
       status = nf90_put_var(fileID,varID(30),isccp%boxtau)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lboxptopisccp) then
       status = nf90_put_var(fileID,varID(31),isccp%boxptop)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lclisccp) then
       status = nf90_put_var(fileID,varID(32),isccp%fq_isccp)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    
    ! MISR simulator output
    if (cfg%LclMISR) then
       status = nf90_put_var(fileID,varID(33),misr%fq_MISR)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))

       status = nf90_put_var(fileID,varID(34),misr%MISR_meanztop)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    
       status = nf90_put_var(fileID,varID(35),misr%MISR_cldarea)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    
    ! MODIS simulator output
    if (cfg%Lcltmodis) then
       status = nf90_put_var(fileID,varID(36),modis%Cloud_Fraction_Total_Mean)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lclwmodis) then
       status = nf90_put_var(fileID,varID(37),modis%Cloud_Fraction_Water_Mean)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lclimodis) then
       status = nf90_put_var(fileID,varID(38),modis%Cloud_Fraction_Ice_Mean)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lclhmodis) then
       status = nf90_put_var(fileID,varID(39),modis%Cloud_Fraction_High_Mean)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lclmmodis) then
       status = nf90_put_var(fileID,varID(40),modis%Cloud_Fraction_Mid_Mean)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Lcllmodis) then
       status = nf90_put_var(fileID,varID(41),modis%Cloud_Fraction_Low_Mean)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Ltautmodis) then
       status = nf90_put_var(fileID,varID(42),modis%Optical_Thickness_Total_Mean)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Ltauwmodis) then
       status = nf90_put_var(fileID,varID(43),modis%Optical_Thickness_Water_Mean)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Ltauimodis) then
       status = nf90_put_var(fileID,varID(44),modis%Optical_Thickness_Ice_Mean)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))      
    endif
    if (cfg%Ltautlogmodis) then
       status = nf90_put_var(fileID,varID(45),modis%Optical_Thickness_Total_LogMean)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Ltauwlogmodis) then
       status = nf90_put_var(fileID,varID(46),modis%Optical_Thickness_Water_LogMean)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))
    endif
    if (cfg%Ltauilogmodis) then
       status = nf90_put_var(fileID,varID(47),modis%Optical_Thickness_Ice_LogMean)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))      
    endif
    if (cfg%Lreffclwmodis) then
       status = nf90_put_var(fileID,varID(48),modis%Cloud_Particle_Size_Water_Mean)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
    endif
    if (cfg%Lreffclimodis) then
       status = nf90_put_var(fileID,varID(49),modis%Cloud_Particle_Size_Ice_Mean)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
    endif
    if (cfg%Lpctmodis) then
       status = nf90_put_var(fileID,varID(50),modis%Cloud_Top_Pressure_Total_Mean)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
    endif
    if (cfg%Llwpmodis) then
       status = nf90_put_var(fileID,varID(51),modis%Liquid_Water_Path_Mean)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
    endif
    if (cfg%Liwpmodis) then
       status = nf90_put_var(fileID,varID(52),modis%Ice_Water_Path_Mean)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
    endif
    if (cfg%Lclmodis) then
       status = nf90_put_var(fileID,varID(53),modis%Optical_Thickness_vs_Cloud_Top_Pressure)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
    endif

    ! Joint simulator products
    if (cfg%Lclcalipso2) then
       status = nf90_put_var(fileID,varID(56),stradar%lidar_only_freq_cloud)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
    endif
    if (cfg%Lcltlidarradar) then
       status = nf90_put_var(fileID,varID(57),stradar%radar_lidar_tcc)
       if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))  
    endif

    ! Close file
    status = nf90_close(fileID)
    if (status .ne. nf90_NoERR) print*,trim(nf90_strerror(status))

  end subroutine write_cosp1_output

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE nc_read_input_file
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE NC_READ_INPUT_FILE(fname,Npnts,Nl,Nhydro,lon,lat,p,ph,z,zh,T,qv,rh,tca,cca, &
                                mr_lsliq,mr_lsice,mr_ccliq,mr_ccice,fl_lsrain,fl_lssnow, &
                                fl_lsgrpl,fl_ccrain,fl_ccsnow,Reff,dtau_s,dtau_c,dem_s,  &
                                dem_c,skt,landmask,mr_ozone,u_wind,v_wind,sunlit,        &
                                emsfc_lw,mode,Nlon,Nlat)
     
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

SUBROUTINE READ_COSP_OUTPUT_NL(cosp_nl,N_OUT_LIST,cfg)
     character(len=*),intent(in) :: cosp_nl
     type(cosp_config),intent(out) :: cfg
     integer,intent(in) :: N_OUT_LIST
     ! Local variables
     integer :: i
     logical :: Lradar_sim,Llidar_sim,Lisccp_sim,Lmodis_sim,Lmisr_sim,Lrttov_sim,Lparasol_sim, &
          Lalbisccp,Latb532,Lboxptopisccp,Lboxtauisccp,LcfadDbze94, &
          LcfadLidarsr532,Lclcalipso2,Lclcalipso,Lclhcalipso,Lclisccp,Lcllcalipso, &
          Lclmcalipso,Lcltcalipso,Lcltlidarradar,Lpctisccp,Ldbze94,Ltauisccp,Lcltisccp, &
          Lclcalipsoliq,Lclcalipsoice,Lclcalipsoun, &
          Lclcalipsotmp,Lclcalipsotmpliq,Lclcalipsotmpice,Lclcalipsotmpun, &
          Lcltcalipsoliq,Lcltcalipsoice,Lcltcalipsoun, &
          Lclhcalipsoliq,Lclhcalipsoice,Lclhcalipsoun, &
          Lclmcalipsoliq,Lclmcalipsoice,Lclmcalipsoun, &
          Lcllcalipsoliq,Lcllcalipsoice,Lcllcalipsoun, &
          LparasolRefl,LclMISR,Lmeantbisccp,Lmeantbclrisccp, &
          Lfracout,LlidarBetaMol532,Ltbrttov, &
          Lcltmodis,Lclwmodis,Lclimodis,Lclhmodis,Lclmmodis,Lcllmodis,Ltautmodis,Ltauwmodis,Ltauimodis,Ltautlogmodis, &
          Ltauwlogmodis,Ltauilogmodis,Lreffclwmodis,Lreffclimodis,Lpctmodis,Llwpmodis, &
          Liwpmodis,Lclmodis
     
     namelist/COSP_OUTPUT/Lradar_sim,Llidar_sim,Lisccp_sim,Lmodis_sim,Lmisr_sim,Lrttov_sim,Lparasol_sim, &
          Lalbisccp,Latb532,Lboxptopisccp,Lboxtauisccp,LcfadDbze94, &
          LcfadLidarsr532,Lclcalipso2,Lclcalipso,Lclhcalipso,Lclisccp, &
          Lcllcalipso,Lclmcalipso,Lcltcalipso,Lcltlidarradar,Lpctisccp,Ldbze94,Ltauisccp, &
          Lclcalipsoliq,Lclcalipsoice,Lclcalipsoun, &
          Lclcalipsotmp,Lclcalipsotmpliq,Lclcalipsotmpice,Lclcalipsotmpun, &
          Lcltcalipsoliq,Lcltcalipsoice,Lcltcalipsoun, &
          Lclhcalipsoliq,Lclhcalipsoice,Lclhcalipsoun, &
          Lclmcalipsoliq,Lclmcalipsoice,Lclmcalipsoun, &
          Lcllcalipsoliq,Lcllcalipsoice,Lcllcalipsoun, &
          Lcltisccp,LparasolRefl,LclMISR,Lmeantbisccp,Lmeantbclrisccp, &
          Lfracout,LlidarBetaMol532,Ltbrttov, &
          Lcltmodis,Lclwmodis,Lclimodis,Lclhmodis,Lclmmodis,Lcllmodis,Ltautmodis,Ltauwmodis,Ltauimodis,Ltautlogmodis, &
          Ltauwlogmodis,Ltauilogmodis,Lreffclwmodis,Lreffclimodis,Lpctmodis,Llwpmodis, &
          Liwpmodis,Lclmodis

     allocate(cfg%out_list(N_OUT_LIST))
     do i=1,N_OUT_LIST
        cfg%out_list(i)=''
     enddo
     open(10,file=cosp_nl,status='old')
     read(10,nml=cosp_output)
     close(10)
     
     ! Deal with dependencies
     if (.not.Lradar_sim) then
        LcfadDbze94   = .false.
        Lclcalipso2    = .false.
        Lcltlidarradar = .false. ! Needs radar & lidar
        Ldbze94        = .false.
        Lclcalipso2    = .false. ! Needs radar & lidar
     endif
     if (.not.Llidar_sim) then
        Latb532 = .false.
        LcfadLidarsr532 = .false.
        Lclcalipso2      = .false. ! Needs radar & lidar
        Lclcalipso       = .false.
        Lclhcalipso      = .false.
        Lcllcalipso      = .false.
        Lclmcalipso      = .false.
        Lcltcalipso      = .false.
        Lcltlidarradar   = .false.
        LparasolRefl    = .false.
        LlidarBetaMol532 = .false.
        Lcltlidarradar = .false. ! Needs radar & lidar
        Lclcalipsoliq    = .false.
        Lclcalipsoice    = .false.
        Lclcalipsoun    = .false.
        Lclcalipsotmp    = .false.
        Lclcalipsotmpun    = .false.
        Lclcalipsotmpliq    = .false.
        Lclcalipsotmpice    = .false.
        Lclhcalipsoliq      = .false.
        Lcllcalipsoliq      = .false.
        Lclmcalipsoliq     = .false.
        Lcltcalipsoliq      = .false.
        Lclhcalipsoice      = .false.
        Lcllcalipsoice      = .false.
        Lclmcalipsoice      = .false.
        Lcltcalipsoice      = .false.
        Lclhcalipsoun      = .false.
        Lcllcalipsoun      = .false.
        Lclmcalipsoun      = .false.
        Lcltcalipsoun      = .false.
     endif
     if (.not.Lisccp_sim) then
        Lalbisccp       = .false.
        Lboxptopisccp   = .false.
        Lboxtauisccp    = .false.
        Lclisccp        = .false.
        Lpctisccp       = .false.
        Ltauisccp       = .false.
        Lcltisccp       = .false.
        Lmeantbisccp    = .false.
        Lmeantbclrisccp = .false.
     endif
     if (.not.Lmisr_sim) then
        LclMISR = .false.
     endif
     if (.not.Lrttov_sim) then
        Ltbrttov = .false.
     endif
     if ((.not.Lradar_sim).and.(.not.Llidar_sim).and. &
          (.not.Lisccp_sim).and.(.not.Lmisr_sim)) then
        Lfracout = .false.
     endif
     if (.not.Lmodis_sim) then
        Lcltmodis=.false.
        Lclwmodis=.false.
        Lclimodis=.false.
        Lclhmodis=.false.
        Lclmmodis=.false.
        Lcllmodis=.false.
        Ltautmodis=.false.
        Ltauwmodis=.false.
        Ltauimodis=.false.
        Ltautlogmodis=.false.
        Ltauwlogmodis=.false.
        Ltauilogmodis=.false.
        Lreffclwmodis=.false.
        Lreffclimodis=.false.
        Lpctmodis=.false.
        Llwpmodis=.false.
        Liwpmodis=.false.
        Lclmodis=.false.
     endif
     if (Lmodis_sim) Lisccp_sim = .true.
     
     cfg%Lstats = .false.
     if ((Lradar_sim).or.(Llidar_sim).or.(Lisccp_sim)) cfg%Lstats = .true.
     
     ! Copy instrument flags to cfg structure
     cfg%Lradar_sim = Lradar_sim
     cfg%Llidar_sim = Llidar_sim
     cfg%Lisccp_sim = Lisccp_sim
     cfg%Lmodis_sim = Lmodis_sim
     cfg%Lmisr_sim  = Lmisr_sim
     cfg%Lrttov_sim = Lrttov_sim
     cfg%Lparasol_sim = Lparasol_sim
     
     ! Flag to control output to file
     cfg%Lwrite_output = .false.
     if (cfg%Lstats.or.cfg%Lmisr_sim.or.cfg%Lrttov_sim) then
        cfg%Lwrite_output = .true.
     endif
     
     ! Output diagnostics
     i = 1
     if (Lalbisccp)        cfg%out_list(i) = 'albisccp'
     i = i+1
     if (Latb532)          cfg%out_list(i) = 'atb532'
     i = i+1
     if (Lboxptopisccp)    cfg%out_list(i) = 'boxptopisccp'
     i = i+1
     if (Lboxtauisccp)     cfg%out_list(i) = 'boxtauisccp'
     i = i+1
     if (LcfadDbze94)      cfg%out_list(i) = 'cfadDbze94'
     i = i+1
     if (LcfadLidarsr532)  cfg%out_list(i) = 'cfadLidarsr532'
     i = i+1
     if (Lclcalipso2)      cfg%out_list(i) = 'clcalipso2'
     i = i+1
     if (Lclcalipso)       cfg%out_list(i) = 'clcalipso'
     i = i+1
     if (Lclhcalipso)      cfg%out_list(i) = 'clhcalipso'
     i = i+1
     if (Lclisccp)         cfg%out_list(i) = 'clisccp'
     i = i+1
     if (Lcllcalipso)      cfg%out_list(i) = 'cllcalipso'
     i = i+1
     if (Lclmcalipso)      cfg%out_list(i) = 'clmcalipso'
     i = i+1
     if (Lcltcalipso)      cfg%out_list(i) = 'cltcalipso'
     i = i+1
     
     if (Lcllcalipsoice)      cfg%out_list(i) = 'cllcalipsoice'
     i = i+1
     if (Lclmcalipsoice)      cfg%out_list(i) = 'clmcalipsoice'
     i = i+1
     if (Lclhcalipsoice)      cfg%out_list(i) = 'clhcalipsoice'
     i = i+1
     if (Lcltcalipsoice)      cfg%out_list(i) = 'cltcalipsoice'
     i = i+1
     if (Lcllcalipsoliq)      cfg%out_list(i) = 'cllcalipsoliq'
     i = i+1
     if (Lclmcalipsoliq)      cfg%out_list(i) = 'clmcalipsoliq'
     i = i+1
     if (Lclhcalipsoliq)      cfg%out_list(i) = 'clhcalipsoliq'
     i = i+1
     if (Lcltcalipsoliq)      cfg%out_list(i) = 'cltcalipsoliq'
     i = i+1
     if (Lcllcalipsoun)      cfg%out_list(i) = 'cllcalipsoun'
     i = i+1
     if (Lclmcalipsoun)      cfg%out_list(i) = 'clmcalipsoun'
     i = i+1
     if (Lclhcalipsoun)      cfg%out_list(i) = 'clhcalipsoun'
     i = i+1
     if (Lcltcalipsoun)      cfg%out_list(i) = 'cltcalipsoun'
     i = i+1
     
     if (Lclcalipsoice)       cfg%out_list(i) = 'clcalipsoice'
     i = i+1
     if (Lclcalipsoliq)       cfg%out_list(i) = 'clcalipsoliq'
     i = i+1
     if (Lclcalipsoun)       cfg%out_list(i) = 'clcalipsoun'
     i = i+1
     
     if (Lclcalipsotmp)       cfg%out_list(i) = 'clcalipsotmp'
     i = i+1
     if (Lclcalipsotmpice)       cfg%out_list(i) = 'clcalipsotmpice'
     i = i+1
     if (Lclcalipsotmpliq)       cfg%out_list(i) = 'clcalipsotmpliq'
     i = i+1
     if (Lclcalipsotmpun)       cfg%out_list(i) = 'clcalipsotmpun'
     i = i+1
     if (Lcltlidarradar)   cfg%out_list(i) = 'cltlidarradar'
     i = i+1
     if (Lpctisccp)        cfg%out_list(i) = 'pctisccp'
     i = i+1
     if (Ldbze94)          cfg%out_list(i) = 'dbze94'
     i = i+1
     if (Ltauisccp)        cfg%out_list(i) = 'tauisccp'
     i = i+1
     if (Lcltisccp)        cfg%out_list(i) = 'cltisccp'
     i = i+1
     !if (Ltoffset)         cfg%out_list(i) = 'toffset'
     i = i+1
     if (LparasolRefl)     cfg%out_list(i) = 'parasolRefl'
     i = i+1
     if (LclMISR)          cfg%out_list(i) = 'clMISR'
     i = i+1
     if (Lmeantbisccp)     cfg%out_list(i) = 'meantbisccp'
     i = i+1
     if (Lmeantbclrisccp)  cfg%out_list(i) = 'meantbclrisccp'
     i = i+1
     if (Lfracout)         cfg%out_list(i) = 'fracout'
     i = i+1
     if (LlidarBetaMol532) cfg%out_list(i) = 'lidarBetaMol532'
     i = i+1
     if (Ltbrttov)         cfg%out_list(i) = 'tbrttov'
     i = i+1
     if (Lcltmodis)        cfg%out_list(i) = 'cltmodis'
     i = i+1
     if (Lclwmodis)        cfg%out_list(i) = 'clwmodis'
     i = i+1
     if (Lclimodis)        cfg%out_list(i) = 'climodis'
     i = i+1
     if (Lclhmodis)        cfg%out_list(i) = 'clhmodis'
     i = i+1
     if (Lclmmodis)        cfg%out_list(i) = 'clmmodis'
     i = i+1
     if (Lcllmodis)        cfg%out_list(i) = 'cllmodis'
     i = i+1
     if (Ltautmodis)       cfg%out_list(i) = 'tautmodis'
     i = i+1
     if (Ltauwmodis)       cfg%out_list(i) = 'tauwmodis'
     i = i+1
     if (Ltauimodis)       cfg%out_list(i) = 'tauimodis'
     i = i+1
     if (Ltautlogmodis)    cfg%out_list(i) = 'tautlogmodis'
     i = i+1
     if (Ltauwlogmodis)    cfg%out_list(i) = 'tauwlogmodis'
     i = i+1
     if (Ltauilogmodis)    cfg%out_list(i) = 'tauilogmodis'
     i = i+1
     if (Lreffclwmodis)    cfg%out_list(i) = 'reffclwmodis'
     i = i+1
     if (Lreffclimodis)    cfg%out_list(i) = 'reffclimodis'
     i = i+1
     if (Lpctmodis)        cfg%out_list(i) = 'pctmodis'
     i = i+1
     if (Llwpmodis)        cfg%out_list(i) = 'lwpmodis'
     i = i+1
     if (Liwpmodis)        cfg%out_list(i) = 'iwpmodis'
     i = i+1
     if (Lclmodis)         cfg%out_list(i) = 'clmodis'
     
     if (i /= N_OUT_LIST) then
        print *, 'COSP_IO: wrong number of output diagnostics'
        print *, i,N_OUT_LIST
        stop
     endif
     
     ! Copy diagnostic flags to cfg structure
     ! ISCCP simulator  
     cfg%Lalbisccp = Lalbisccp
     cfg%Latb532 = Latb532
     cfg%Lboxptopisccp = Lboxptopisccp
     cfg%Lboxtauisccp = Lboxtauisccp
     cfg%Lmeantbisccp = Lmeantbisccp
     cfg%Lmeantbclrisccp = Lmeantbclrisccp
     cfg%Lclisccp = Lclisccp
     cfg%Lpctisccp = Lpctisccp
     cfg%Ltauisccp = Ltauisccp
     cfg%Lcltisccp = Lcltisccp
     ! CloudSat simulator  
     cfg%Ldbze94 = Ldbze94
     cfg%LcfadDbze94 = LcfadDbze94
     ! CALIPSO/PARASOL simulator  
     cfg%LcfadLidarsr532 = LcfadLidarsr532
     cfg%Lclcalipso2 = Lclcalipso2
     cfg%Lclcalipso = Lclcalipso
     cfg%Lclhcalipso = Lclhcalipso
     cfg%Lcllcalipso = Lcllcalipso
     cfg%Lclmcalipso = Lclmcalipso
     cfg%Lcltcalipso = Lcltcalipso
     cfg%Lclhcalipsoice = Lclhcalipsoice
     cfg%Lcllcalipsoice = Lcllcalipsoice
     cfg%Lclmcalipsoice = Lclmcalipsoice
     cfg%Lcltcalipsoice = Lcltcalipsoice
     cfg%Lclhcalipsoliq = Lclhcalipsoliq
     cfg%Lcllcalipsoliq = Lcllcalipsoliq
     cfg%Lclmcalipsoliq = Lclmcalipsoliq
     cfg%Lcltcalipsoliq = Lcltcalipsoliq
     cfg%Lclhcalipsoun = Lclhcalipsoun
     cfg%Lcllcalipsoun = Lcllcalipsoun
     cfg%Lclmcalipsoun = Lclmcalipsoun
     cfg%Lcltcalipsoun = Lcltcalipsoun
     cfg%Lclcalipsoice = Lclcalipsoice
     cfg%Lclcalipsoliq = Lclcalipsoliq
     cfg%Lclcalipsoun = Lclcalipsoun
     cfg%Lclcalipsotmp = Lclcalipsotmp
     cfg%Lclcalipsotmpice = Lclcalipsotmpice
     cfg%Lclcalipsotmpliq = Lclcalipsotmpliq
     cfg%Lclcalipsotmpun = Lclcalipsotmpun
     cfg%Lcltlidarradar = Lcltlidarradar
     cfg%LparasolRefl = LparasolRefl
     ! MISR simulator  
     cfg%LclMISR = LclMISR
     ! Other
     cfg%Ltoffset = .false.!Ltoffset
     cfg%Lfracout = Lfracout
     cfg%LlidarBetaMol532 = LlidarBetaMol532
     ! RTTOV
     cfg%Ltbrttov = Ltbrttov
     ! MODIS simulator  
     cfg%Lcltmodis=Lcltmodis
     cfg%Lclwmodis=Lclwmodis
     cfg%Lclimodis=Lclimodis
     cfg%Lclhmodis=Lclhmodis
     cfg%Lclmmodis=Lclmmodis
     cfg%Lcllmodis=Lcllmodis
     cfg%Ltautmodis=Ltautmodis
     cfg%Ltauwmodis=Ltauwmodis
     cfg%Ltauimodis=Ltauimodis
     cfg%Ltautlogmodis=Ltautlogmodis
     cfg%Ltauwlogmodis=Ltauwlogmodis
     cfg%Ltauilogmodis=Ltauilogmodis
     cfg%Lreffclwmodis=Lreffclwmodis
     cfg%Lreffclimodis=Lreffclimodis
     cfg%Lpctmodis=Lpctmodis
     cfg%Llwpmodis=Llwpmodis
     cfg%Liwpmodis=Liwpmodis
     cfg%Lclmodis=Lclmodis
   END SUBROUTINE READ_COSP_OUTPUT_NL

  
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
end module mod_cosp1_io
