; ###################################################################
;
; The purpose of this program is to read in the 1D sample input data.
;
; ###################################################################
pro read_data,file,data
  
fileID = ncdf_open(file)
ncdf_varget,fileID,ncdf_varid(fileID,'lat'),lat
ncdf_varget,fileID,ncdf_varid(fileID,'lon'),lon
ncdf_varget,fileID,ncdf_varid(fileID,'landmask'),lsmask
ncdf_varget,fileID,ncdf_varid(fileID,'orography'),orography
ncdf_varget,fileID,ncdf_varid(fileID,'psfc'),sfcP
ncdf_varget,fileID,ncdf_varid(fileID,'height'),hgt
ncdf_varget,fileID,ncdf_varid(fileID,'height_half'),hgt_half
ncdf_varget,fileID,ncdf_varid(fileID,'T_abs'),temp
ncdf_varget,fileID,ncdf_varid(fileID,'qv'),qv
ncdf_varget,fileID,ncdf_varid(fileID,'rh'),rh
ncdf_varget,fileID,ncdf_varid(fileID,'pfull'),pres
ncdf_varget,fileID,ncdf_varid(fileID,'phalf'),pres_half
ncdf_varget,fileID,ncdf_varid(fileID,'mr_lsliq'),mr_lsliq
ncdf_varget,fileID,ncdf_varid(fileID,'mr_lsice'),mr_lsice
ncdf_varget,fileID,ncdf_varid(fileID,'mr_ccliq'),mr_ccliq
ncdf_varget,fileID,ncdf_varid(fileID,'mr_ccice'),mr_ccice
ncdf_varget,fileID,ncdf_varid(fileID,'fl_lsrain'),fl_lsrain
ncdf_varget,fileID,ncdf_varid(fileID,'fl_lssnow'),fl_lssnow
ncdf_varget,fileID,ncdf_varid(fileID,'fl_lsgrpl'),fl_lsgrpl
ncdf_varget,fileID,ncdf_varid(fileID,'fl_ccrain'),fl_ccrain
ncdf_varget,fileID,ncdf_varid(fileID,'fl_ccsnow'),fl_ccsnow
ncdf_varget,fileID,ncdf_varid(fileID,'tca'),tca
ncdf_varget,fileID,ncdf_varid(fileID,'cca'),cca
ncdf_varget,fileID,ncdf_varid(fileID,'Reff'),Reff
ncdf_varget,fileID,ncdf_varid(fileID,'dtau_c'),dtau_c
ncdf_varget,fileID,ncdf_varid(fileID,'dtau_s'),dtau_s
ncdf_varget,fileID,ncdf_varid(fileID,'dem_c'),dem_c
ncdf_varget,fileID,ncdf_varid(fileID,'dem_s'),dem_s
ncdf_varget,fileID,ncdf_varid(fileID,'skt'),skt
ncdf_varget,fileID,ncdf_varid(fileID,'sunlit'),sunlit
ncdf_varget,fileID,ncdf_varid(fileID,'u_wind'),uwind
ncdf_varget,fileID,ncdf_varid(fileID,'v_wind'),vwind
ncdf_varget,fileID,ncdf_varid(fileID,'mr_ozone'),o3
ncdf_varget,fileID,ncdf_varid(fileID,'emsfc_lw'),emsfcLW
ncdf_close,fileID

data = {lon:lon,lat:lat,lsmask:lsmask,orography:orography,sfcP:sfcP,hgt:hgt,$
        hgt_half:hgt_half,temp:temp,qv:qv,rh:rh,pres:pres,pres_half:pres_half,$
        mr_lsliq:mr_lsliq,mr_lsice:mr_lsice,mr_ccliq:mr_ccliq,fl_lsrain:fl_lsrain,$
        fl_lssnow:fl_lssnow,fl_lsgrpl:fl_lsgrpl,fl_ccrain:fl_ccrain,fl_ccsnow:fl_ccsnow,$
        tca:tca,cca:cca,Reff:Reff,dtau_c:dtau_c,dtau_s:dtau_s,dem_c:dem_c,dem_s:dem_s,$
        skt:skt,sunlit:sunlit,uwind:uwind,vwind:vwind,o3:o3,emsfcLW:emsfcLW}




; ###################################################################
; END PROGRAM
; ###################################################################
end
