'reinit' ; 'cc'

'open cfed_Ze.ctl'

RGN.1=TR ; RGN.2=NM ; RGN.3=NP ; RGN.4=SM ; RGN.5=SP
TYP.1=WRM ; TYP.2=UPP
TYP.3=CNV ; TYP.4=STR ; TYP.5=OTH
TYP.6=ALL ; TYP.7=CLD

r=1

c=1
while(c<=7)

'c'

'set grads off'
'set mproj off'

'set xlopts 1 6 0.25'
'set ylopts 1 6 0.25'

'set x 1 35'
'set y 1 57'
'set e 'r

wid=2.

if (c<=5)
'set z 'c
'nrm=sum(cfed,x=1,x=35)'
'var=cfed/nrm/'wid

else
if (c=6)
'set z 1'
'nrm=sum(sum(cfed,z=1,z=4),x=1,x=35)'
'var=sum(cfed,z=1,z=4)/nrm/'wid

else
if (c=7)
'set z 1'
'nrm=sum(sum(cfed,z=2,z=4),x=1,x=35)'
'var=sum(cfed,z=2,z=4)/nrm/'wid

endif
endif
endif

'set lon -40 30'
'set xlint 10'
'set lat -80 20'
'set ylint 10'
'set yflip on'

'set parea 1.0 10.0 1.2 7.0'

'color 0 0.1 0.002 -kind precip -gxout grfill'
'd var'
'xcbar 2.0 9.0 0.45 0.6 -edge triangle -fs 10 -ft 8 -fh 0.18 -fw 0.16'

'draw string 5.5 7.5 CFED Ze 'RGN.r' 'TYP.c

'gxprint ./figure/EC_CFED_Ze_'RGN.r'_'TYP.c'.eps'

c=c+1
endwhile
