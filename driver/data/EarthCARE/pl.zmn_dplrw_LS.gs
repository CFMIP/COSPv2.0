#grads -a 1.4142

'reinit' ; 'cc'

fname='zmn_dplrw_LS'

'open vd_zonal_mean.ctl'

'set grads off'

'set xlopts 1 6 0.25'
'set ylopts 1 6 0.25'

'set lat -90 90'
'set xlint 30'
'set lev 0 15'
'set ylint 1'

'set t 4'

'total=vdm'
'smpl=vdn'

'set parea 1.0 10.0 1.2 7.0'

'color -4 2 0.1 -kind (26,51,51)->(35,85,130)->(61,144,199)->(120,197,204)->(231,255,232)-(0)->(1,64,38)->(71,108,25)->(154,125,66)->(206,185,156)->(242,239,246)->(236,196,225)->(204,124,186)->(164,65,138)->(101,2,75)'
'd maskout(total/smpl,smpl)'
'xcbar 2.0 9.0 0.45 0.6 -edge triangle -fs 20 -ft 8 -fh 0.18 -fw 0.16'

'set gxout contour'
'set grid off' ; 'set frame off' ; 'set xlab off' ; 'set ylab off'

'Temp=tem/ten'

min=-100 ; max=30
int=10 ; bint=40

'set lwid 13 4.5'
wt=12  ; bt=5
wt2=13 ; bt2=12

'set cmin 'min
'set cmax 'max
'set cint 'int
'set ccolor 0'
'set cthick 'wt
'set cstyle 1'
'set clab off'
'd Temp-273.15'

'set cmin 'min
'set cmax 'max
'set cint 'int
'set ccolor 1'
'set cthick 'bt
'set cstyle 1'
'set clab off'
'd Temp-273.15'

'set cmin 'min
'set cmax 'max
'set cint 'bint
'set ccolor 0'
'set cthick 'wt2
'set cstyle 1'
'set clab off'
'd Temp-273.15'

'set cmin 'min
'set cmax 'max
'set cint 'bint
'set ccolor 1'
'set cthick 'bt2
'set cstyle 1'
'set clab on'
'set clopts 1 6 0.15'
'd Temp-273.15'


'gxprint ./figure/'fname'.eps'
