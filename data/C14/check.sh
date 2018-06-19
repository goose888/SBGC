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

# Collect the model simulated D14C
for dir in `ls -d */`
do
   if [ -e ${dir}deltac14.txt ]
   then
     echo ${dir} >> site_4_eval.dat
#    tail -n 2 ${dir}/isam.log
#    head -n 50 ${dir}/isam.log
     tail -n 1 ${dir}deltac14.txt >> dc14.dat
   fi
done

# Collect the D14C data from He et al., 2016
/home/zsj/work/C14/C14processing/retrieve_d14c.py
