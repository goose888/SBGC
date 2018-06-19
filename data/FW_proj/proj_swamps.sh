#!/bin/bash
datarcp45="./BC_Hist_plus_RCP45"
datarcp85="./BC_Hist_plus_RCP85"
precrcp45="./Precip_rcp45"
precrcp85="./Precip_rcp85"
fw="./FW"
molist=("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12")
nmon=${#molist[@]}

 # Extract precipitation datasets
 # RCP45
# cd ${datarcp45}
# flist=`ls -1 199[23456789]*.nc`
# cd ..
# for i in ${flist};
# do
#    ncks -v PRECTmms ${datarcp45}/${i} ${precrcp45}/${i}
# done
#
# cd ${datarcp45}
# flist=`ls -1 200[6789]*.nc`
# cd ..
# for i in ${flist};
# do
#    ncks -v PRECTmms ${datarcp45}/${i} ${precrcp45}/${i}
# done
#
# cd ${datarcp45}
# flist=`ls -1 201[012345]*.nc`
# cd ..
# for i in ${flist};
# do
#    ncks -v PRECTmms ${datarcp45}/${i} ${precrcp45}/${i}
# done
#
# cd ${datarcp45}
# flist=`ls -1 2???_*.nc`
# cd ..
# for i in ${flist};
# do
#    ncks -v PRECIP ${datarcp45}/${i} ${precrcp45}/${i}
#    #ncap2 -O -s 'PRECTmms=PRECIP*3600*6' ${precrcp45}/${i} ${precrcp45}/${i}
#    ncap2 -O -s 'PRECTmms=PRECIP' ${precrcp45}/${i} ${precrcp45}/${i}
# done

 # RCP85
# cd ${datarcp85}
# flist=`ls -1 199[23456789]*.nc`
# cd ..
# for i in ${flist};
# do
#    ncks -v PRECTmms ${datarcp85}/${i} ${precrcp85}/${i}
# done
#
# cd ${datarcp85}
# flist=`ls -1 200[12345]*.nc`
# cd ..
# for i in ${flist};
# do
#    ncks -v PRECTmms ${datarcp85}/${i} ${precrcp85}/${i}
# done
# cd ${datarcp85}
# flist=`ls -1 2???_*.nc`
# cd ..
# for i in ${flist};
# do
#    ncks -v PRECIP ${datarcp85}/${i} ${precrcp85}/${i}
#    #ncap2 -O -s 'PRECTmms=PRECIP*3600*6' ${precrcp85}/${i} ${precrcp85}/${i}
#    ncap2 -O -s 'PRECTmms=PRECIP' ${precrcp85}/${i} ${precrcp85}/${i}
# done

#  cd ${datarcp85}
#  flist=`ls -1 209[123456789]_*.nc`
#  cd ..
#  for i in ${flist};
#  do
#     ncks -v PRECIP ${datarcp85}/${i} ${precrcp85}/${i}
#     #ncap2 -O -s 'PRECTmms=PRECIP*3600*6' ${precrcp85}/${i} ${precrcp85}/${i}
#     ncap2 -O -s 'PRECTmms=PRECIP' ${precrcp85}/${i} ${precrcp85}/${i}
#  done
# 
#  cd ${datarcp85}
#  flist=`ls -1 2100_*.nc`
#  cd ..
#  for i in ${flist};
#  do
#     ncks -v PRECIP ${datarcp85}/${i} ${precrcp85}/${i}
#     #ncap2 -O -s 'PRECTmms=PRECIP*3600*6' ${precrcp85}/${i} ${precrcp85}/${i}
#     ncap2 -O -s 'PRECTmms=PRECIP' ${precrcp85}/${i} ${precrcp85}/${i}
#  done

#  # Accumulate the precipitation for specific time period
#  # For Jan, Mar, May, Jul, Aug, Oct, Dec: 1 - 10, 11 - 20, 21 - 31
#  # For Apr, Jun, Sep, Nov: 1 - 10, 11 - 20, 21 - 30
#  # For Feb: 1 - 10, 11 - 20, 21 - 28
#  # RCP45
#  cd ${precrcp45}
#  # Historical
#  yrlist=( $(ls -1 *-01.nc | cut -c 1-4) )
#  nyr=${#yrlist[@]}
#  for((p=0; p<$nyr; p++)) do
#     for((k=0; k<$nmon; k++)) do
#        # First 10 days
#        let "m = k * 3 + 1"
#        ncks -d Time,0,39 -v PRECTmms ${yrlist[p]}-${molist[k]}.nc PR_${yrlist[p]}_${m}.nc
#        ncks -O --mk_rec_dmn Time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        ncra -O -d Time,0,39 -y ttl PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        ncwa -O -a Time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        # Second 10 days
#        let "m += 1"
#        ncks -d Time,40,79 -v PRECTmms ${yrlist[p]}-${molist[k]}.nc PR_${yrlist[p]}_${m}.nc
#        ncks -O --mk_rec_dmn Time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        ncra -O -d Time,0,39 -y ttl PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        ncwa -O -a Time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        # The rest of the days
#        let "m += 1"
#        if (( (k==0)||(k==2)||(k==4)||(k==6)||(k==7)||(k==9)||(k==11) )); then
#           ncks -d Time,80,123 -v PRECTmms ${yrlist[p]}-${molist[k]}.nc PR_${yrlist[p]}_${m}.nc
#           ncks -O --mk_rec_dmn Time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#           ncra -O -d Time,0,43 -y ttl PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#           ncwa -O -a Time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        fi
#        if (( (k==3)||(k==5)||(k==8)||(k==10) )); then
#           ncks -d Time,80,119 -v PRECTmms ${yrlist[p]}-${molist[k]}.nc PR_${yrlist[p]}_${m}.nc
#           ncks -O --mk_rec_dmn Time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#           ncra -O -d Time,0,39 -y ttl PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#           ncwa -O -a Time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        fi
#        if (( (k==1) )); then
#           ncks -d Time,80,111 -v PRECTmms ${yrlist[p]}-${molist[k]}.nc PR_${yrlist[p]}_${m}.nc
#           ncks -O --mk_rec_dmn Time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#           ncra -O -d Time,0,31 -y ttl PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#           ncwa -O -a Time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        fi
#        echo "Month"
#        echo ${molist[k]}
#     done
#     echo "Year"
#     echo ${yrlist[p]}
#  done
#  #Projection
#  yrlist=( $(ls -1 *_01.nc | cut -c 1-4) )
#  nyr=${#yrlist[@]}
#  for((p=0; p<$nyr; p++)) do
#     for((k=0; k<$nmon; k++)) do
#        # First 10 days
#        let "m = k * 3 + 1"
#        ncks -d time,0,39 -v PRECIP ${yrlist[p]}_${molist[k]}.nc PR_${yrlist[p]}_${m}.nc
#        ncrename -v PRECIP,PRECTmms PR_${yrlist[p]}_${m}.nc
#        ncks -O --mk_rec_dmn time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        ncra -O -d time,0,39 -y ttl PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        ncwa -O -a time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        # Second 10 days
#        let "m += 1"
#        ncks -d time,40,79 -v PRECIP ${yrlist[p]}_${molist[k]}.nc PR_${yrlist[p]}_${m}.nc
#        ncrename -v PRECIP,PRECTmms PR_${yrlist[p]}_${m}.nc
#        ncks -O --mk_rec_dmn time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        ncra -O -d time,0,39 -y ttl PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        ncwa -O -a time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        # The rest of the days
#        let "m += 1"
#        if (( (k==0)||(k==2)||(k==4)||(k==6)||(k==7)||(k==9)||(k==11) )); then
#           ncks -d time,80,123 -v PRECIP ${yrlist[p]}_${molist[k]}.nc PR_${yrlist[p]}_${m}.nc
#           ncrename -v PRECIP,PRECTmms PR_${yrlist[p]}_${m}.nc
#           ncks -O --mk_rec_dmn time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#           ncra -O -d time,0,43 -y ttl PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#           ncwa -O -a time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        fi
#        if (( (k==3)||(k==5)||(k==8)||(k==10) )); then
#           ncks -d time,80,119 -v PRECIP ${yrlist[p]}_${molist[k]}.nc PR_${yrlist[p]}_${m}.nc
#           ncrename -v PRECIP,PRECTmms PR_${yrlist[p]}_${m}.nc
#           ncks -O --mk_rec_dmn time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#           ncra -O -d time,0,39 -y ttl PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#           ncwa -O -a time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        fi
#        if (( (k==1) )); then
#           ncks -d time,80,111 -v PRECIP ${yrlist[p]}_${molist[k]}.nc PR_${yrlist[p]}_${m}.nc
#           ncrename -v PRECIP,PRECTmms PR_${yrlist[p]}_${m}.nc
#           ncks -O --mk_rec_dmn time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#           ncra -O -d time,0,31 -y ttl PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#           ncwa -O -a time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        fi
#        echo "Month"
#        echo ${molist[k]}
#     done
#     echo "Year"
#     echo ${yrlist[p]}
#  done
#  cd ..
#  # RCP85
#  cd ${precrcp85}
#  # Historical
#  yrlist=( $(ls -1 *-01.nc | cut -c 1-4) )
#  nyr=${#yrlist[@]}
#  for((p=0; p<$nyr; p++)) do
#     for((k=0; k<$nmon; k++)) do
#        # First 10 days
#        let "m = k * 3 + 1"
#        ncks -d Time,0,39 -v PRECTmms ${yrlist[p]}-${molist[k]}.nc PR_${yrlist[p]}_${m}.nc
#        ncks -O --mk_rec_dmn Time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        ncra -O -d Time,0,39 -y ttl PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        ncwa -O -a Time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        # Second 10 days
#        let "m += 1"
#        ncks -d Time,40,79 -v PRECTmms ${yrlist[p]}-${molist[k]}.nc PR_${yrlist[p]}_${m}.nc
#        ncks -O --mk_rec_dmn Time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        ncra -O -d Time,0,39 -y ttl PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        ncwa -O -a Time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        # The rest of the days
#        let "m += 1"
#        if (( (k==0)||(k==2)||(k==4)||(k==6)||(k==7)||(k==9)||(k==11) )); then
#           ncks -d Time,80,123 -v PRECTmms ${yrlist[p]}-${molist[k]}.nc PR_${yrlist[p]}_${m}.nc
#           ncks -O --mk_rec_dmn Time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#           ncra -O -d Time,0,43 -y ttl PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#           ncwa -O -a Time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        fi
#        if (( (k==3)||(k==5)||(k==8)||(k==10) )); then
#           ncks -d Time,80,119 -v PRECTmms ${yrlist[p]}-${molist[k]}.nc PR_${yrlist[p]}_${m}.nc
#           ncks -O --mk_rec_dmn Time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#           ncra -O -d Time,0,39 -y ttl PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#           ncwa -O -a Time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        fi
#        if (( (k==1) )); then
#           ncks -d Time,80,111 -v PRECTmms ${yrlist[p]}-${molist[k]}.nc PR_${yrlist[p]}_${m}.nc
#           ncks -O --mk_rec_dmn Time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#           ncra -O -d Time,0,31 -y ttl PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#           ncwa -O -a Time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        fi
#        echo "Month"
#        echo ${molist[k]}
#     done
#     echo "Year"
#     echo ${yrlist[p]}
#  done
#  #Projection
#  yrlist=( $(ls -1 *_01.nc | cut -c 1-4) )
#  nyr=${#yrlist[@]}
#  for((p=0; p<$nyr; p++)) do
#     for((k=0; k<$nmon; k++)) do
#        # First 10 days
#        let "m = k * 3 + 1"
#        ncks -d time,0,39 -v PRECIP ${yrlist[p]}_${molist[k]}.nc PR_${yrlist[p]}_${m}.nc
#        ncrename -v PRECIP,PRECTmms PR_${yrlist[p]}_${m}.nc
#        ncks -O --mk_rec_dmn time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        ncra -O -d time,0,39 -y ttl PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        ncwa -O -a time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        # Second 10 days
#        let "m += 1"
#        ncks -d time,40,79 -v PRECIP ${yrlist[p]}_${molist[k]}.nc PR_${yrlist[p]}_${m}.nc
#        ncrename -v PRECIP,PRECTmms PR_${yrlist[p]}_${m}.nc
#        ncks -O --mk_rec_dmn time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        ncra -O -d time,0,39 -y ttl PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        ncwa -O -a time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        # The rest of the days
#        let "m += 1"
#        if (( (k==0)||(k==2)||(k==4)||(k==6)||(k==7)||(k==9)||(k==11) )); then
#           ncks -d time,80,123 -v PRECIP ${yrlist[p]}_${molist[k]}.nc PR_${yrlist[p]}_${m}.nc
#           ncrename -v PRECIP,PRECTmms PR_${yrlist[p]}_${m}.nc
#           ncks -O --mk_rec_dmn time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#           ncra -O -d time,0,43 -y ttl PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#           ncwa -O -a time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        fi
#        if (( (k==3)||(k==5)||(k==8)||(k==10) )); then
#           ncks -d time,80,119 -v PRECIP ${yrlist[p]}_${molist[k]}.nc PR_${yrlist[p]}_${m}.nc
#           ncrename -v PRECIP,PRECTmms PR_${yrlist[p]}_${m}.nc
#           ncks -O --mk_rec_dmn time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#           ncra -O -d time,0,39 -y ttl PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#           ncwa -O -a time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        fi
#        if (( (k==1) )); then
#           ncks -d time,80,111 -v PRECIP ${yrlist[p]}_${molist[k]}.nc PR_${yrlist[p]}_${m}.nc
#           ncrename -v PRECIP,PRECTmms PR_${yrlist[p]}_${m}.nc
#           ncks -O --mk_rec_dmn time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#           ncra -O -d time,0,31 -y ttl PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#           ncwa -O -a time PR_${yrlist[p]}_${m}.nc PR_${yrlist[p]}_${m}.nc
#        fi
#        echo "Month"
#        echo ${molist[k]}
#     done
#     echo "Year"
#     echo ${yrlist[p]}
#  done
#  cd ..

 # Select the data from year 1992 - 2005 as baseline and calculate the mean variability factor v
 # RCP45
 cd ${precrcp45}
 tlist1=( $(ls -1 PR_2006_?.nc | cut -c 9-9) )
 tlist2=( $(ls -1 PR_2006_??.nc | cut -c 9-10) )
 tlist=("${tlist1[@]}" "${tlist2[@]}")
 nt=${#tlist[@]}
# for((p=0; p<$nt; p++)) do
#    flist=`ls PR_199[23456789]_${tlist[p]}.nc ; ls PR_200[12345]_${tlist[p]}.nc`
#    ncea -v PRECTmms ${flist} PR_${tlist[p]}.nc
# done
 for((p=0; p<$nt; p++)) do
    flist=`ls PR_200[6789]_${tlist[p]}.nc ; ls PR_201[012345]_${tlist[p]}.nc`
    ncea -v PRECTmms ${flist} PR_${tlist[p]}.nc
 done
 cd ..
 # RCP85
# cd ${precrcp85}
# for((p=0; p<$nt; p++)) do
#    flist=`ls PR_199[23456789]_${tlist[p]}.nc ; ls PR_200[12345]_${tlist[p]}.nc`
#    ncea -v PRECTmms ${flist} PR_${tlist[p]}.nc
# done
# cd ..

# # FW
# cd ${fw}
# for ((p=0; p<$nt; p++)) do
#     flist=`ls fw_frac_199[23456789]_${tlist[p]}.nc ; ls fw_frac_200[12345]_${tlist[p]}.nc`
#     ncea ${flist} fw_frac_${tlist[p]}.nc
# done
# cd ..
# # Calculate the v factor
# # RCP45
# cd ${precrcp45}
#  cp ../${fw}/fw_frac_?.nc .
#  cp ../${fw}/fw_frac_??.nc .
# for ((p=0; p<$nt; p++)) do
#     ncrename -d Lon,longitude -d Lat,latitude fw_frac_${tlist[p]}.nc
#     ncks -A PR_${tlist[p]}.nc fw_frac_${tlist[p]}.nc
#     ncap2 -O -s 'v=FW/(PRECTmms+0.001)' fw_frac_${tlist[p]}.nc temp_${tlist[p]}.nc
#  #   ncks -v v temp_${tlist[p]}.nc v_${tlist[p]}.nc
#  #   rm temp_${tlist[p]}.nc
#     mv temp_${tlist[p]}.nc v_${tlist[p]}.nc
#     ncrename -v PRECTmms,PRECTstd v_${tlist[p]}.nc
# done
# cd ..
# # RCP85
# cd ${precrcp85}
#  cp ../${fw}/fw_frac_?.nc .
#  cp ../${fw}/fw_frac_??.nc .
# for ((p=0; p<$nt; p++)) do
#     ncrename -d Lon,longitude -d Lat,latitude fw_frac_${tlist[p]}.nc
#     ncks -A PR_${tlist[p]}.nc fw_frac_${tlist[p]}.nc
#     ncap2 -O -s 'v=FW/(PRECTmms+0.001)' fw_frac_${tlist[p]}.nc temp_${tlist[p]}.nc
#  #   ncks -v v temp_${tlist[p]}.nc v_${tlist[p]}.nc
#  #   rm temp_${tlist[p]}.nc
#     mv temp_${tlist[p]}.nc v_${tlist[p]}.nc
#     ncrename -v PRECTmms,PRECTstd v_${tlist[p]}.nc
# done
# cd ..

# Get the projection of the inundated map
# RCP45
cd ${precrcp45}
yrlist=( $(ls -1 PR_*_1.nc | cut -c 4-7) )
nyr=${#yrlist[@]}
for((p=0; p<$nyr; p++)) do
   for((k=0; k<$nt; k++)) do
      ncks -A v_${tlist[k]}.nc PR_${yrlist[p]}_${tlist[k]}.nc
   #   ncap2 -s 'FW=(PRECTmms+0.001)*v' PR_${yrlist[p]}_${tlist[k]}.nc temp_${yrlist[p]}_${tlist[k]}.nc
      ncap2 -s 'FW=2*FW/(1+exp(-((PRECTmms+0.001)/(PRECTstd+0.001))+1))' PR_${yrlist[p]}_${tlist[k]}.nc temp_${yrlist[p]}_${tlist[k]}.nc
      ncks -v FW temp_${yrlist[p]}_${tlist[k]}.nc FW_${yrlist[p]}_${tlist[k]}.nc
      rm temp_${yrlist[p]}_${tlist[k]}.nc
   done
done
cd ..
# # RCP85
# cd ${precrcp85}
# for((p=0; p<$nyr; p++)) do
#    for((k=0; k<$nt; k++)) do
#       ncks -A v_${tlist[k]}.nc PR_${yrlist[p]}_${tlist[k]}.nc
#   #    ncap2 -s 'FW=(PRECTmms+0.001)*v' PR_${yrlist[p]}_${tlist[k]}.nc temp_${yrlist[p]}_${tlist[k]}.nc
#       ncap2 -s 'FW=2*FW/(1+exp(-((PRECTmms+0.001)/(PRECTstd+0.001))+1))' PR_${yrlist[p]}_${tlist[k]}.nc temp_${yrlist[p]}_${tlist[k]}.nc
#       ncks -v FW temp_${yrlist[p]}_${tlist[k]}.nc FW_${yrlist[p]}_${tlist[k]}.nc
#       rm temp_${yrlist[p]}_${tlist[k]}.nc
#    done
# done
# cd ..
