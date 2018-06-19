#!/bin/bash
# Gather the isotopic carbon profile for each site

# Clear old files
if [ -e dc14.dat ]
then
   rm -f dc14.dat
fi
if [ -e site_4_eval.dat ]
then
   rm -f site_4_eval.dat
fi
if [ -e socprof.dat ]
then
   rm -f socprof.dat
fi

# Collect the model simulated D14C
for dir in `ls -d */`
do
   if [ -e ${dir}deltac14.txt ]
   then
     echo ${dir} >> site_4_eval.dat
#    tail -n 2 ${dir}/isam.log
#    head -n 50 ${dir}/isam.log
     tail -n 1 ${dir}deltac14.txt >> dc14.dat
     tail -n 1 ${dir}socprof.txt >> socprof.dat
   else
     echo ${dir} >> site_abandoned.dat
   fi
done

# Processing texts
sed 's/\///' site_4_eval.dat > site_w_prof.dat
sed 's/\///' site_abandoned.dat > site_out_mask.dat
paste site_w_prof.dat socprof.dat > isam_soc.dat
paste site_w_prof.dat dc14.dat > isam_dc14.dat
rm -f dc14.dat
rm -f socprof.dat
rm -f site_4_eval.dat
rm -f site_w_prof.dat
rm -f site_abandoned.dat

# Collect the D14C data from He et al., 2016

