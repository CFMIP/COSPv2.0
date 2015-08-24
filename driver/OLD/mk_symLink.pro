; ########################################################################################
; The purpose of this program is to create symbolic links to the
; coefficient files used by RTTOV
; ########################################################################################

dir = '/Projects/Clouds/dswales/RTTOV/rttov93/rtcoef_rttov7/'

filesLONG  = file_search(dir,'rtcoef*.dat')
nfiles     = n_elements(filesLONG)
filesSHORT = strarr(nfiles)
for ij=0,nfiles-1 do begin
	filesSHORT(ij)=strmid(filesLONG(ij),strlen(dir),strlen(filesLONG(ij))-strlen(dir))
	spawn,'ln -s '+filesLONG(ij)+' '+filesSHORT(ij)
end






; END PROGRAM
end