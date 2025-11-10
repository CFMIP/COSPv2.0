#!/bin/bash

dir="ctrl"
tmax=7

gfortran -fbacktrace -fbounds-check gt2input.f90 -L/opt/homebrew/lib -I/opt/homebrew/include -lnetcdff -lnetcdf -lhdf5_fortran -lhdf5_hl -lhdf5 -o exe.gt2input.out

cat <<EOF > setup.nml
&setup
 dir='$dir',
 tmax=$tmax,
&end

EOF

rm -f ../MIROC_inputs/{$dir}_MIROCinput.nc

./exe.gt2input.out
