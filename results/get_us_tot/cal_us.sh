#!/bin/bash

flist=($(ls -1 /data/jain1/a/sshu3/CH4_CONUS/hist/Global_1DSBGC.bgc-yearly-2d*.nc))
flist2=($(ls -1 ./fw_annual_mean/hist/fw_frac_mean_*.nc))
nf=${#flist[@]}
for((k=0; k<$nf; k++)) do
#for id in ${!flist[*]}; do
   i=${flist[k]}
   j=${flist2[k]}
   #echo ${i}
   #echo ${j}
   Rscript get_us_tot.r ${i} ${j}
   #Rscript get_us_prod.r ${i} ${j}
   #Rscript get_us_oxid.r ${i} ${j}
   #Rscript get_us_aere.r ${i} ${j}
   #Rscript get_us_ebul.r ${i} ${j}
   #Rscript get_us_diff.r ${i} ${j}
done
