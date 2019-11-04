#!/bin/bash

flist=($(ls -1 /data/jain1/c/sshu3/SBGC/results/SBGC_regional/Global_1DSBGC.bgc-yearly-2d*.nc))
nf=${#flist[@]}
for((k=0; k<$nf; k++)) do
#for id in ${!flist[*]}; do
   i=${flist[k]}
   #echo ${i}
   Rscript get_psoc_tot.r ${i}
done
