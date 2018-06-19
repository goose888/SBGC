#!/bin/bash

## Please load the NCO before running this script
yr=("1992" "1993" "1994" "1995" "1996" "1997" "1998" "1999" "2000" "2001" "2002" "2003" "2004" "2005" "2006" "2007" "2008" "2009" "2010" "2011" "2012")
nyr=${#yr[@]}

for((k=0; k<$nyr; k++)) do
   # Merge fw into a year
   # ncecat -u time fw_frac_${yr[k]}_1.nc fw_frac_${yr[k]}_2.nc fw_frac_${yr[k]}_3.nc fw_frac_${yr[k]}_4.nc fw_frac_${yr[k]}_5.nc fw_frac_${yr[k]}_6.nc fw_frac_${yr[k]}_7.nc fw_frac_${yr[k]}_8.nc fw_frac_${yr[k]}_9.nc fw_frac_${yr[k]}_10.nc fw_frac_${yr[k]}_11.nc fw_frac_${yr[k]}_12.nc fw_frac_${yr[k]}_13.nc fw_frac_${yr[k]}_14.nc fw_frac_${yr[k]}_15.nc fw_frac_${yr[k]}_16.nc fw_frac_${yr[k]}_17.nc fw_frac_${yr[k]}_18.nc fw_frac_${yr[k]}_19.nc fw_frac_${yr[k]}_20.nc fw_frac_${yr[k]}_21.nc fw_frac_${yr[k]}_22.nc fw_frac_${yr[k]}_23.nc fw_frac_${yr[k]}_24.nc fw_frac_${yr[k]}_25.nc fw_frac_${yr[k]}_26.nc fw_frac_${yr[k]}_27.nc fw_frac_${yr[k]}_28.nc fw_frac_${yr[k]}_29.nc fw_frac_${yr[k]}_30.nc fw_frac_${yr[k]}_31.nc fw_frac_${yr[k]}_32.nc fw_frac_${yr[k]}_33.nc fw_frac_${yr[k]}_34.nc fw_frac_${yr[k]}_35.nc fw_frac_${yr[k]}_36.nc fw_frac_${yr[k]}.nc
   # Annual mean fw
   ncea fw_frac_${yr[k]}_1.nc fw_frac_${yr[k]}_2.nc fw_frac_${yr[k]}_3.nc fw_frac_${yr[k]}_4.nc fw_frac_${yr[k]}_5.nc fw_frac_${yr[k]}_6.nc fw_frac_${yr[k]}_7.nc fw_frac_${yr[k]}_8.nc fw_frac_${yr[k]}_9.nc fw_frac_${yr[k]}_10.nc fw_frac_${yr[k]}_11.nc fw_frac_${yr[k]}_12.nc fw_frac_${yr[k]}_13.nc fw_frac_${yr[k]}_14.nc fw_frac_${yr[k]}_15.nc fw_frac_${yr[k]}_16.nc fw_frac_${yr[k]}_17.nc fw_frac_${yr[k]}_18.nc fw_frac_${yr[k]}_19.nc fw_frac_${yr[k]}_20.nc fw_frac_${yr[k]}_21.nc fw_frac_${yr[k]}_22.nc fw_frac_${yr[k]}_23.nc fw_frac_${yr[k]}_24.nc fw_frac_${yr[k]}_25.nc fw_frac_${yr[k]}_26.nc fw_frac_${yr[k]}_27.nc fw_frac_${yr[k]}_28.nc fw_frac_${yr[k]}_29.nc fw_frac_${yr[k]}_30.nc fw_frac_${yr[k]}_31.nc fw_frac_${yr[k]}_32.nc fw_frac_${yr[k]}_33.nc fw_frac_${yr[k]}_34.nc fw_frac_${yr[k]}_35.nc fw_frac_${yr[k]}_36.nc fw_frac_mean_${yr[k]}.nc
done


