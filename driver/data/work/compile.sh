#!/bin/bash

gfortran -fbacktrace -fbounds-check gt2input.f90 -L/opt/homebrew/lib -I/opt/homebrew/include -lnetcdff -lnetcdf -lhdf5_fortran -lhdf5_hl -lhdf5 -o exe.gt2input.out
