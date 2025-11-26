'reinit' ; 'cc'

'open cfed_Rs.ctl'

RGN.1=SP ; RGN.2=TR ; RGN.3=NP
TYP.1=UPP ; TYP.2=WRM ; TYP.3=STR
TYP.4=5Ze ; TYP.5=MZe ; TYP.6=PIA

'c'

'set grads off'
'set mproj off'

'set xlopts 1 6 0.25'
'set ylopts 1 6 0.25'

r=2
c=4
wid=0.2

'set x 1 60'
'set y 1 57'
'set z 'r
'set e 'c

'nrm=sum(cfed,x=1,x=60)'
'var=cfed/nrm/'wid

'set lon -6 6'
'set xlint 1'
'set lat -80 20'
'set ylint 10'
'set yflip on'

'set parea 1.0 10.0 1.2 7.0'

'color 0 1 0.02 -kind precip -gxout grfill'
'd var'
'xcbar 2.0 9.0 0.45 0.6 -edge triangle -fs 10 -ft 8 -fh 0.18 -fw 0.16'

'draw string 5.5 7.5 CFED Res 'RGN.r' 'TYP.c

'gxprint ./figure/EC_CFED_Rs_'RGN.r'_'TYP.c'.eps'
