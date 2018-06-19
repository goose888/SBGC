#!/bin/bash

flist=`ls -1 FW_*.nc`
#for i in ${flist};
#do
#   Rscript mask_fw.r ${i}
#done
#for i in ${flist};
#do
#   mv ${i} fw_frac_${i}
#done
rename _FW_ _ fw_frac_*.nc
