#!/bin/csh

## Please load the NCO before running this script
yr=("2091" "2092" "2093" "2094" "2095" "2096" "2097" "2098" "2099" "2100")
nyr=${#yr[@]}

for((k=0; k<$nyr; k++)) do
   # Mean fw for each 10-days
   ncecat -u time FW_${yr[k]}_1.nc FW_${yr[k]}_2.nc FW_${yr[k]}_3.nc FW_${yr[k]}_4.nc FW_${yr[k]}_5.nc FW_${yr[k]}_6.nc FW_${yr[k]}_7.nc FW_${yr[k]}_8.nc FW_${yr[k]}_9.nc FW_${yr[k]}_10.nc FW_${yr[k]}_11.nc FW_${yr[k]}_12.nc FW_${yr[k]}_13.nc FW_${yr[k]}_14.nc FW_${yr[k]}_15.nc FW_${yr[k]}_16.nc FW_${yr[k]}_17.nc FW_${yr[k]}_18.nc FW_${yr[k]}_19.nc FW_${yr[k]}_20.nc FW_${yr[k]}_21.nc FW_${yr[k]}_22.nc FW_${yr[k]}_23.nc FW_${yr[k]}_24.nc FW_${yr[k]}_25.nc FW_${yr[k]}_26.nc FW_${yr[k]}_27.nc FW_${yr[k]}_28.nc FW_${yr[k]}_29.nc FW_${yr[k]}_30.nc FW_${yr[k]}_31.nc FW_${yr[k]}_32.nc FW_${yr[k]}_33.nc FW_${yr[k]}_34.nc FW_${yr[k]}_35.nc FW_${yr[k]}_36.nc FW_${yr[k]}.nc
done

ncrcat FW_2091.nc FW_2092.nc FW_2093.nc FW_2094.nc FW_2095.nc FW_2096.nc FW_2097.nc FW_2098.nc FW_2099.nc FW_2100.nc FW_9100.nc

# Get max
ncra -y max FW_9100.nc FW_max.nc

# Get min
ncra -y min FW_9100.nc FW_min.nc
