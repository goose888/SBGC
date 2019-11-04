file_s21=ISAM_S2.1_nbp_monthly.nc
file_s22=ISAM_S2.2.0_nbp_monthly.nc
file_s23=ISAM_S2.3.0_nbp_monthly.nc
cdo selyear,2018 $file_s21 $file_s21.tmp

ncdiff $file_s21.tmp $file_s22 diff_$file_s22
ncdiff $file_s21.tmp $file_s23 diff_$file_s23

rm *.tmp
