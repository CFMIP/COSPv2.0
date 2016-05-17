; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
; Copyright (c) 2015, Regents of the University of Colorado
; All rights reserved.
;
; Redistribution and use in source and binary forms, with or without modification, are 
; permitted provided that the following conditions are met:
;
; 1. Redistributions of source code must retain the above copyright notice, this list of 
;    conditions and the following disclaimer.
;
; 2. Redistributions in binary form must reproduce the above copyright notice, this list
;    of conditions and the following disclaimer in the documentation and/or other 
;    materials provided with the distribution.
;
; 3. Neither the name of the copyright holder nor the names of its contributors may be 
;    used to endorse or promote products derived from this software without specific prior
;    written permission.
;
; THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
; EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
; MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
; THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
; SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT 
; OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
; INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
; LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
; OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
;
; History
; May 2015 - D. Swales - Original version
;
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;
; The purpose of this program is to compare COSP simulator data. This
; code is the IDL equivalent to the PYTHON code cosp_tool1.py.
;
; Coded by Dustin Swales 2015 CIRES/NOAA-ESRL-PSD
;
; VARIABLE NAME     TYPE        DESCRIPTION 
; dir               string      location of COSP simulator output data
; dirREF            string      location of COSP reference data
; dataset           string      either 1D or 2D
; thresh            float       error threshold(optional)
;
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

; Dataset dimensionality (1D or 2D)
dim='1D'

; Threshold for error reporting
thresh = 1e-5

; Reference data location
dirREF = '/home/dswales/Projects/COSP/COSPv2.0/driver/data/output/ref1D/'

; Data location
dir = '/home/dswales/Projects/COSP/COSPv2.0/driver/data/output/'+dim+'/'
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;                               No changes needed below
; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

; Files to compare
files = file_search(dir,'*'+dim+'*.nc')

; Names of variables to compare
vars = strarr(n_elements(files))
for ij=0,n_elements(files)-1 do vars(ij)=strmid(files(ij),strlen(dir),strlen(files(ij))-(55+strlen(dir)))

; Compare files (if available)
count=0
minV      = fltarr(n_elements(files))-999.
maxV      = fltarr(n_elements(files))-999.
occurence = fltarr(n_elements(files))-999.
for ij=0,n_elements(files)-1 do begin
   fileRef = file_search(dirREF,vars(ij)+'_*'+dim+'*.nc')
   if (strcmp(fileRef(0),'') eq 0) then begin

      ; Read in data
      fileID = ncdf_open(fileREF)
      ncdf_varget,fileID,ncdf_varid(fileID,vars(ij)),varREF
      ncdf_close,fileID
      fileID = ncdf_open(files(ij))
      ncdf_varget,fileID,ncdf_varid(fileID,vars(ij)),var
      ncdf_close,fileID
      
      ; Compare data
      diff  = varREF-var        ; Absolute difference
      diff2 = var/varREF-1      ; Relative difference
      fix1 = where(finite(diff2,/infinity) eq 1) 
      fix2 = where(finite(diff2,/Nan) eq 1)
      fix3 = where(diff2 eq -1)
      fix4 = where(diff eq 0)
      if (fix1(0) ne -1) then diff2(fix1) = 1
      if (fix2(0) ne -1) then diff2(fix2) = 0
      if (fix3(0) ne -1) then diff2(fix3) = 0
      if (fix4(0) ne -1) then diff2(fix4) = 0
      
      ; Compute frequency of error
      bad  = where(diff ne 0 and abs(diff2) gt thresh)
      if (bad(0) ne -1) then begin
         top           = n_elements(bad)
         bot           = n_elements(diff)
         occurence(ij) = 100*float(top)/float(bot)
      endif
      minV(ij) = min(diff2(where(diff2 ne 0)))
      maxV(ij) = max(diff2(where(diff2 ne 0)))
   endif
end

; Print results
count=0
print,'#######################################################################################'
print,'Comparing...'
print,'COSP output:    ',dir
print,'COSP reference: ',dirREF
print,'#######################################################################################'
for ij=0,n_elements(files)-1 do begin
   ;if (minV(ij) ne -999 and abs(minV(ij)) gt thresh or abs(maxV(ij)) gt thresh) then begin
   if (occurence(ij) ne -999) then begin
      if (count eq 0) then print,'DIFFERENCES EXIST FOR THE FOLLOWING:'
      slen1 = strlen(vars(ij))
      fstring = '(a'+string(slen1,format='(i2.2)')+',a'+string(17-slen1,format='(i2.2)')+',a7,f9.4,a24,e11.3,a4,e11.3)'
      print,vars(ij),'','differ ',occurence(ij),'% of the time.     From ',minV(ij),' to ',maxV(ij),format=fstring
      count=count+1
   endif
end
if (count eq 0) then print,'ALL FILES MATCH'
print,'#######################################################################################'

; ################################################################
; END PROGRAM
; ################################################################
end
