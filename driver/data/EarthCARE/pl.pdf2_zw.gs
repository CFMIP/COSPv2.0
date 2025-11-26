'reinit' ; 'cc'

'open pdf2_zw.ctl'

RGN.1=TR ; RGN.2=NM ; RGN.3=NP ; RGN.4=SM ; RGN.5=SP
TYP.1=WRM ; TYP.2=UPP
TYP.3=CNV ; TYP.4=STR ; TYP.5=OTH
TYP.6=ALL ; TYP.7=CLD

'c'

'set grads off'
'set mproj off'

'set xlopts 1 6 0.25'
'set ylopts 1 6 0.25'

r=1
c=4
xwid=2.
ywid=0.2

'set x 1 35'
'set y 1 80'
'set e 'r

if (c<=5)
'set z 'c
'nrm=sum(sum(cfed,x=1,x=35),y=1,y=80)'
'var=cfed/nrm/'xwid'/'ywid

else
if (c=6)
'set z 1'
'nrm=sum(sum(sum(cfed,z=1,z=5),x=1,x=35),y=1,y=80)'
'var=sum(cfed,z=1,z=5)/nrm/'xwid'/'ywid

else
if (c=7)
'set z 1'
'nrm=sum(sum(sum(cfed,z=2,z=4),x=1,x=35),y=1,y=80)'
'var=sum(cfed,z=2,z=4)/nrm/'xwid'/'ywid

endif
endif
endif

'set lon -40 30'
'set xlint 10'
'set lat -6 6'
'set ylint 1'

'set parea 1.0 10.0 1.2 7.0'

'color 0 3.0 0.05 -kind precip -gxout grfill'
'd var*100'
'xcbar 2.0 9.0 0.45 0.6 -edge triangle -fs 10 -ft 8 -fh 0.18 -fw 0.16'

'draw string 5.5 7.5 PDF Ze-Wa 'RGN.r' 'TYP.c

'gxprint ./figure/EC_pdf2_ZeWa_'RGN.r'_'TYP.c'.eps'
