#grads -a 1.4142

'reinit' ; 'cc'

'open cfed_ptyp.ctl'

RGN.1=SP ; RGN.2=TR ; RGN.3=NP
TYP.1=UPP ; TYP.2=WRM ; TYP.3=STR
TYP.4=5Ze ; TYP.5=MZe ; TYP.6=PIA

'set grads off'
'set mproj off'

'set xlopts 1 6 0.25'
'set ylopts 1 6 0.25'

r=2
c=3
wid=2.

'set x 1 9'
'set y 1 57'
'set z 'r
'set e 'c

'nrm=sum(cfed,x=1,x=9)'

'set lon 0'
'other=cfed'

'set lon 1'
'warm=cfed'

'set lon 2'
'cool=cfed'

'set x 1'
'set lat -80 20'
'set gxout bar'
'set xyrev on'
'set yflip on'
'set ccolor 3'
'set vrange 0 1'

'set parea 1.0 10.0 1.2 7.0'

'd (warm+cool)/nrm'

'gxprint ./figure/particle_type.eps'
