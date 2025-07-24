#!/bin/bash
initial_dir=${PWD}
cd data/outputs/UKMO/
files=$( ls cosp2_output.um_global.gfortran.kgo.v004.nc.*.png )

for file_a in ${files} ;
do
  echo ${file_a}
  file_b='cosp2_output.um_global.'$( echo ${file_a} | cut -d. -f6- )
  echo ${file_b}
  file_out='combined.'$( echo ${file_a} | cut -d. -f7- )
  echo ${file_out}
  convert -append ${file_a} ${file_b} ${file_out}
done
convert -page a5 combined.*.png combined.pdf
cd ${initial_dir}
