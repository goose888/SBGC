#!/bin/bash

flist=($(ls -1 Global_1DSBGC.bgc-yearly-2d*.nc))
flist2=($(ls -1 ./fw_annual_mean/rcp85/fw_frac_mean_*.nc))
nf=${#flist[@]}
for((k=0; k<$nf; k++)) do
#for id in ${!flist[*]}; do
   i=${flist[k]}
   j=${flist2[k]}
   #echo ${i}
   #echo ${j}
   Rscript get_us_tot.r ${i}
   #Rscript get_us_sink.r ${i} ${j}
   #Rscript get_us_prod.r ${i} ${j}
done
