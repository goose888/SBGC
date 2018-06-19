#!/bin/tcsh

#----------------------------------------------------------------------------------
# Author: Shijie Shu
# Created on: Nov 2 2015
# A CSV file containing site ID, Lat, Lon, PFT, MOSSPFT and # of Moss layer 
# information must be presented
#----------------------------------------------------------------------------------
# User defined options and variables
#----------------------------------------------------------------------------------
# Paths
set workdir=`pwd`
set datarepo='/data/jain1/b/team/datasets4'
set isamdir='/data/keeling/a/sshu3/ISAM_master/ISAM'
# I/O
set siteinfo='case_info.csv'
## set outfile='isamtexture.csv'
set proffile='socprofeq.dat'
set difffile='diffurate.dat'
set aldfile='ald.dat'
set c14file='deltac14.dat'
# Options
set setup=false
set genmet=false
set runcases=false
set postprocessing=true
set clean=false
# Sub options for runcases
set runtype=historical
# Sub options for postprocess
set FTRENDY='ISAM_S3_cSoil.nc'
set D14_hist=0


#----------------------------------------------------------------------------------
# USE CAUTIONS WHEN EDITING BELOW
#----------------------------------------------------------------------------------

echo 'Retreive list of profile IDs:'
awk -F , 'NR>1 { print $1}' $workdir/$siteinfo > list_prof

if($setup == 'true') then

   echo '################ ISAM Cases Setup ###############'
 
   echo '(1) Building new cases ...'
   #----------------------------------------------------------------------------------
   # get profile ID from file
   #----------------------------------------------------------------------------------
   #set num=`wc -l < list_prof`
   #----------------------------------------------------------------------------------
   # Loop for each profile
   #----------------------------------------------------------------------------------
   set lines=`cat list_prof`
   set j=1
   while ( $j <= $#lines)
      echo 'Soil sample ID:' $lines[$j]
      if ( -d $lines[$j]) then
        echo 'Case directory already exists!'
      else
        echo 'Create new directory for soil profile!'
        mkdir $lines[$j]
      endif
      cd $lines[$j]
      #----------------------------------------------------------------------------------
      # Set specific options in namelist
      #----------------------------------------------------------------------------------
      set CASENAME=$lines[$j]
      echo $CASENAME
      set lat=`awk -F , -v casename=$CASENAME '$1==casename { print $3}' $workdir/$siteinfo`
      echo $lat
      set lon=`awk -F , -v casename=$CASENAME '$1==casename { print $2}' $workdir/$siteinfo`
      echo $lon
      #----------------------------------------------------------------------------------
      # Calculate Lat/Lon idx based on lat/lon from site info file
      #----------------------------------------------------------------------------------
      set lonid = `echo ${lon} | awk '{lon = $1; if (lon < 0.) {lonid = (720. + lon*2. + 0.5); res = int(lonid); if(lonid-res > 0.5) lonid= res+1; else lonid = res;} else {lonid = (lon*2 + 0.5); res=int(lonid); if(lonid-res > 0.5) lonid= res+1; else lonid = res;} printf "%i", lonid;}'`
      set latid = `echo ${lat} | awk '{lat = $1; latid = (lat*2. + 0.5 + 180.); res=int(latid); if(latid-res > 0.5) latid= res+1; else latid = res; printf "%i", latid;}'`
      #----------------------------------------------------------------------------------
      # Create ISAM NAMELIST 
      #----------------------------------------------------------------------------------
      set starttime=1901
      set runtime=110
      /bin/tcsh ${workdir}/generate_atms.csh ${CASENAME} ${starttime} ${runtime} ${lonid} ${latid}
      set starttime=1901
      set runtime=20
      /bin/tcsh ${workdir}/generate_spinup.csh ${CASENAME} ${starttime} ${runtime} ${lonid} ${latid}
      set starttime=1901
      set runtime=110
      /bin/tcsh ${workdir}/generate_hist.csh ${CASENAME} ${starttime} ${runtime} ${lonid} ${latid}
      # set starttime=2011
      # set runtime=90
      # /bin/tcsh ./generate_proj.csh ${CASENAME} ${starttime} ${runtime} ${lonid} ${latid}
      ### /bin/tcsh ./generate_fixedproj.csh ${CASENAME} ${starttime} ${runtime} ${lonid} ${latid}
      # /bin/tcsh ./generate_tproj.csh ${CASENAME} ${starttime} ${runtime} ${lonid} ${latid}
      # /bin/tcsh ./generate_wproj.csh ${CASENAME} ${starttime} ${runtime} ${lonid} ${latid}

      #----------------------------------------------------------------------------------
      # Create directory
      #----------------------------------------------------------------------------------
      if ( -d data ) then
        echo 'Data diresctory already exists!'
      else
        echo 'Create data directory for case:' ${CASENAME}
        mkdir data
        #----------------------------------------------------------------------------------
        # Link data from data repository and prepare input data
        #----------------------------------------------------------------------------------
        ln -s ${datarepo}/* data/
        # ln -s /data/jain1/b/team/datasets4/del_c14_annual_50000bp_2015.nc ./data/
        rm data/initial_bgp
        rm data/initial_bgc
        mkdir data/initial_bgp
        mkdir data/initial_bgc
        mkdir data/SITE_MET_OUTPUT_FROM_GLOBAL
      endif
 
      if ( -d output ) then
        echo 'Output directory already exists!'
      else
        echo 'Create output directory for case:' ${CASENAME}
        mkdir output
      endif

      #----------------------------------------------------------------------------------
      # Copy ISAM executable program
      #----------------------------------------------------------------------------------
      cp ${isamdir}/isam .

      #----------------------------------------------------------------------------------
      # Create Job script and submit the case
      #----------------------------------------------------------------------------------
      /bin/tcsh ${workdir}/new_jobscript.csh ${CASENAME} genmet
      /bin/tcsh ${workdir}/new_jobscript.csh ${CASENAME} spinup
      /bin/tcsh ${workdir}/new_jobscript.csh ${CASENAME} hist
      # /bin/tcsh ./new_jobscript.csh ${CASENAME} proj
      # /bin/tcsh ./new_jobscript.csh ${CASENAME} tproj
      # /bin/tcsh ./new_jobscript.csh ${CASENAME} wproj

      #----------------------------------------------------------------------------------
      # Return to the main work directory
      #----------------------------------------------------------------------------------
      echo 'Done ' ${CASENAME} '!'
      cd $workdir
      @ j = $j + 1

   end  # End of Loop inside case directory
   
   echo '################ ISAM Cases Setup Done ###############'

else

   echo '################ Skip ISAM Cases Setup ###############'

endif

#----------------------------------------------------------------------------------
# Submit the job script to generate met forcing from CRUNCEP
#----------------------------------------------------------------------------------
if($genmet == 'true') then
   if( -f list_prof ) then
      set cases=`cat list_prof`
      while ( $j <= $#cases)
        echo 'Now purge case:' $cases[$j]
        cd $workdir/$cases[$j]
 
        sbatch genmet_$cases[$j].sh
        mv ./output/*forcing-CRU_NCEP.nc ./data/SITE_MET_OUTPUT_FROM_GLOBAL/
        cd $workdir
        @ j = $j + 1
      end
   else
      echo "Case list not found. Check if file exists!"
   endif
endif

#----------------------------------------------------------------------------------
# Submit the job script to run simulation
#----------------------------------------------------------------------------------
if($runcases == 'true') then

   if( -f list_prof ) then
      set cases=`cat list_prof`
      set pft=`cat pftlist`
      set j=1
      while ( $j <= $#cases)
         echo 'Now run case:' $cases[$j]
         if ( -d $cases[$j]) then
           echo 'Found, submit the spinup job.'
           cd $workdir/$cases[$j]
           cp $isamdir/isam .

           echo 'Warning! Pruge all residual files from previous run!'
           echo 'If you want to keep these files please uncomment here scripts!'
           rm *.txt
           rm isam*.log

           # Set the latest initial file as model input
           cd data/initial_bgp
           mv `ls -t1 | head -n 1` bgp.isam_initial.nc
           cd ../initial_bgc
           mv `ls -t1 | head -n 1` bgc.isam_initial.nc
           # Set softlink for any other missing data
           ### ln -s ${workdir}/projected_climate/* data/SITE_MET_OUTPUT_FROM_GLOBAL/
           ### ln -s ~/team/datasets4/projected_co2/ data/
           ### ln -s ~/team/datasets4/projected_ndep/ data/
           ### ln -s ~/team/datasets4/RTM ./data/

           # Back to the case directory
           cd $workdir/$cases[$j]

           # Change settings in the namelist
           # <User can specify the change in the namelist as they like>
           set pft=`awk -F , -v casename=$cases[$j] '$1==casename { print $4}' $workdir/$siteinfo`
           set pmoss=`awk -F , -v casename=$cases[$j] '$1==casename { print $5}' $workdir/$siteinfo`
           set nmoss=`awk -F , -v casename=$cases[$j] '$1==casename { print $6}' $workdir/$siteinfo`
           /bin/tcsh ${workdir}/changenm.csh single_pft_num -9999 $pft namelist.hist
           /bin/tcsh ${workdir}/changenm.csh mosspft -9999 $pmoss namelist.hist
           /bin/tcsh ${workdir}/changenm.csh mosslayers -9999 $nmoss namelist.hist

           # Submit the run
           sbatch hist_$cases[$j].sh

           # Arrange the outputs
           mkdir spinup
           mv *.txt spinup/
           mv *.dat spinup/
           mv isam*.log spinup/
           # mv *.txt forced_proj_t
           # mv free_spin free_proj

           # Return
           cd $workdir

         else

           echo 'Failed to find the case directory. Check if setup was successful.'

         endif
         @ j = $j + 1
      end
   else

      echo 'Case list not found. Check if file exists!'
      exit 1

   endif

endif

#----------------------------------------------------------------------------------
# Postprocess for model data comparison
#----------------------------------------------------------------------------------
if($postprocessing == 'true') then

   # Clear old files
   if( -e dc14.dat ) rm -f dc14.dat
   if( -e site_4_eval.dat ) rm -f site_4_eval.dat
   if( -e socprof.dat ) rm -f socprof.dat
   if( -e totsoc.dat ) rm -f totsoc.dat
   if( -e totsoc_trendy.dat ) rm -f totsoc_trendy.dat

   # Counter for the line number
   set lin=1
   # Loop over every directory (i.e., site)
   foreach dir ( `ls -d */` )
      if( -e ${dir}deltac14.txt ) then
         # Collect the IDs of selected sites
         echo ${dir} >> site_4_eval.dat
         # Collect the model simulated D14C
         if( $D14_hist == 1 ) then
            # For historical case, get the year collecting D14C data
            set num=`sed -n "${lin}p" < ext_line`
            @ lin = $lin + 1
            sed -n "${num}p" < ${dir}deltac14.txt >> dc14.dat
         else
            # For spin-up case, get the last line of simulated D14C
            tail -n 1 ${dir}deltac14.txt >> dc14.dat
         endif
         # Collect simulated SOC profile (kgC/m2/lev)
         tail -n 1 ${dir}socprof.txt >> socprof.dat
         # Collect simulated Total SOC (kgC/m2)
         tail -n 1 ${dir}som.txt >> totsoc.dat
   
         # Extract the SOC from corresponding point from the trendy output
         cd ${dir}
         # Get Lon ID from the namelist
         set lonid=`sed -n -e 's/^.*single_x = //p' < namelist.spinup`
         @ lonid = $lonid - 1
         # Get Lat ID from the namelist
         set latid=`sed -n -e 's/^.*single_y = //p' < namelist.spinup`
         @ latid = $latid - 1
         cd ${workdir}
         set var='cSoil'
         set timeid=0
         ncks -v ${var} -d time,${timeid} -d lat,${latid} -d lon,${lonid} ${FTRENDY} > temp.dat
         sed -n -e 's/^.*cSoil\[/\[/p' < temp.dat | awk -F'[ =]' '{ print $2 }'  >> totsoc_trendy.dat
         rm temp.dat
      else
         echo "D14C model output of site ${dir} is missing."
      endif
   end   
   
      
   # Processing texts
   sed 's/\///' site_4_eval.dat > site_w_prof.dat
   paste site_w_prof.dat socprof.dat > isam_soc.dat
   paste site_w_prof.dat dc14.dat > isam_dc14.dat
   paste site_w_prof.dat totsoc.dat > isam_totsoc.dat
   paste site_w_prof.dat totsoc_trendy.dat > trendy_totsoc.dat
   rm -f dc14.dat
   rm -f socprof.dat
   rm -f totsoc.dat
   rm -f totsoc_trendy.dat
   rm -f site_4_eval.dat
   rm -f site_w_prof.dat
   rm -f temp.dat

endif

#----------------------------------------------------------------------------------
# CLEAN OLD RESULTS
#----------------------------------------------------------------------------------
if($clean == 'true') then
   if( -f list_prof ) then
      set cases=`cat list_prof`
      while ( $j <= $#cases )
         echo 'Now purge case:' $cases[$j]
         cd $cases[$j]
 
         rm output/*
         rm -r spinup
         rm -r historical
         rm -r projection
         rm *.txt
         rm *.dat
         rm isam*.log
         cd $workdir
         @ j = $j + 1
       end
   else
      echo "Case list not found. Please setup again."
   endif
endif

# Set the namelist for each variable
###(BASH) while read name
###do
###   echo "$name"
###done < list_histel


