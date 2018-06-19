#!/bin/tcsh

#----------------------------------------------------------------------------------
# Author: Shijie Shu
# Created on: Nov 2 2015
#----------------------------------------------------------------------------------
# User defined options and variables
set workdir=`pwd`
set datarepo='/data/jain1/b/team/datasets4'
set isamdir='/data/keeling/a/sshu3/ISAM_master/ISAM'
set siteinfo='site_info.csv'
set outfile='isamtexture.csv'
set proffile='socprofeq.dat'
set difffile='diffurate.dat'
set aldfile='ald.dat'
set c14file='deltac14.dat'
set setup=false
set runcases=true
set postprocessing=false

#----------------------------------------------------------------------------------
# DO NOT EDIT BELOW
#----------------------------------------------------------------------------------
if($setup == 'true') then

   echo '################ ISAM Cases Setup ###############'
 
   echo '(1) Check if directories missing and create new direcotories:'
   #----------------------------------------------------------------------------------
   # get profile ID from file
   #----------------------------------------------------------------------------------
   echo 'Extract profile lists:'
   awk -F , 'NR>1 { print $1}' $workdir/$siteinfo > list_prof
   #set num=`wc -l < list_prof`
   #----------------------------------------------------------------------------------
   # Loop for each profile
   #----------------------------------------------------------------------------------
   set lines=`cat list_prof`
   set j=1
   while ( $j <= $#lines)
      echo 'Soil sample ID:' $lines[$j]
      if ( -d $lines[$j]) then
        echo 'Found, no need to create directory.'
      else
        echo 'Create new directory for soil profile.'
        mkdir $lines[$j]
      endif
      cd $lines[$j]
      #----------------------------------------------------------------------------------
      # Set specific options in namelist
      #----------------------------------------------------------------------------------
      set CASENAME=$lines[$j]
      set starttime=1901
      set runtime=112
     ## set starttime=1901
     ## set runtime=20
      echo $CASENAME
      set lat=`awk -F , -v casename=$CASENAME '$1==casename { print $7}' $workdir/$siteinfo`
      echo $lat
      set lon=`awk -F , -v casename=$CASENAME '$1==casename { print $6}' $workdir/$siteinfo`
      echo $lon
      #----------------------------------------------------------------------------------
      # Calculate Lat/Lon idx based on lat/lon from site info file
      #----------------------------------------------------------------------------------
      set lonid = `echo ${lon} | awk '{lon = $1; if (lon < 0.) {lonid = (720. + lon*2. + 0.5); res = int(lonid); if(lonid-res > 0.5) lonid= res+1; else lonid = res;} else {lonid = (lon*2 + 0.5); res=int(lonid); if(lonid-res > 0.5) lonid= res+1; else lonid = res;} printf "%i", lonid;}'`
      set latid = `echo ${lat} | awk '{lat = $1; latid = (lat*2. + 0.5 + 180.); res=int(latid); if(latid-res > 0.5) latid= res+1; else latid = res; printf "%i", latid;}'`
      #----------------------------------------------------------------------------------
      # Create ISAM NAMELIST 
      #----------------------------------------------------------------------------------
      echo 'Generate namelist for case ' ${CASENAME}
         cat > namelist.genmet << EOF1
&isamcfg

        runname = '${CASENAME}'

        bgp_bgc_mode = 1
        restart = .false.
        run_type = 0
        CN_RATIO = 0
        initial_bgp_in_datadir = 'initial_bgp/bgp.isam_initial.nc'
        initial_bgc_in_datadir = 'initial_bgc/bgc.isam_initial.nc'

        repeat_cycles = 1
        start_time = ${starttime}
        run_time = ${runtime}

        offline_clim_data = 'CRU_NCEP'
        luc_data = 1

        clim_disturb= .true.
        co2_disturb= .false.
        luc_disturb= .false.
        n_disturb= .false.
        fna_func = .false.

        region_mask = 1,1,1,1,1,1,1,1,1,0,0

        datadir   = 'data'
        outputdir = 'output'

        single_point = .true.
        single_x = ${lonid}
        single_y = ${latid}
        generate_site_met_nc_from_global = .true.

        single_pft = .false.
        single_pft_num = -999
        use_site_met = .false.
        sitename   = '${CASENAME}'
        site_ref_ht = -9999

        fixed_harvest   = .false.
        fixed_planttime = .false.
        crop_mode       = 'generic'

        hist_freq_yr = 100
        save_hist_yr = .true.
        save_hist_mon = .false.
        save_hist_bgp_to_bgc = .false.

        restart_freq_yr = 100

        isam_deltim = 3600
/
EOF1
 
         echo 'Generate namelist for spinup case ' ${CASENAME}
         cat > namelist.spinup << EOF1
&isamcfg

        runname = '${CASENAME}'

        bgp_bgc_mode = 3
        restart = .false.
        run_type = 0
        CN_RATIO = 0

        initial_bgp_in_datadir = 'initial_bgp/bgp.isam_initial.nc'
        initial_bgc_in_datadir = 'initial_bgc/bgc.isam_initial.nc'

        repeat_cycles = 1000
        start_time = ${starttime}
        run_time = ${runtime}

        offline_clim_data = 'CRU_NCEP'
        luc_data = 1

        clim_disturb= .true.
        clim_random= .false.
        co2_disturb= .false.
        c14_disturb= .false.
        luc_disturb= .false.
        n_disturb= .false.
        fna_func = .true.

        region_mask = 1,1,1,1,1,1,1,1,1,0,0

        datadir   = 'data'
        outputdir = 'output'

        single_point = .true.
        single_x = ${lonid}
        single_y = ${latid}
        generate_site_met_nc_from_global = .false.

        single_pft = .true.
        single_pft_num = -9999
        use_site_met = .true.
        sitename   = '${CASENAME}'
        site_ref_ht = 30
        dp_drain = 1     ! Topographic impedance for lateral drainage
        tao = 8      ! Impedance for calculating cyroturbation rate
        site_met_option = 'CRU_NCEP'
        site_met_timezone = 'local'

        dyn_veg         = .true.
        dyn_soc         = .false.   ! Switch of the coupled SOC content - thermal properties dynamics
        coupledsoc      = .false.   ! Switch of the coupled SOC content - thermal properties dynamics
        projection      = .false.   ! Apply the RCP8.5 projected atm CO2 concentration
        tsep            = .false.   ! (Under projection) Only consider SOC's thermal impact
        wsep            = .false.   ! (Under projection) Only consider SOC's hydraulic impact
        soc_saved       = .false.   ! Read saved soc from file
        fixed_harvest   = .false.
        fixed_planttime = .false.
        crop_mode       = 'generic'

        param_from_ascii = .false.
        soil_param_file = 'soil_param.dat'  ! Read soil texture from file
        soil_data = 'Combined'    ! Combined / GSDE

        hist_freq_yr = 10000
        save_hist_yr = .true.
        save_hist_mon = .false.
        save_hist_bgp_to_bgc = .false.

        restart_freq_yr = 10000
        site_restart = .true.

        dt_ch4      = 3600
        isam_deltim = 3600
/
EOF1

# Simulation stage
 
         echo 'Generate namelist for historical case ' ${CASENAME}
         cat > namelist.hist << EOF1
&isamcfg

        runname = '${CASENAME}'

        bgp_bgc_mode = 3
        restart = .false.
        run_type = 2
        CN_RATIO = 0

        initial_bgp_in_datadir = 'initial_bgp/bgp.isam_initial.nc'
        initial_bgc_in_datadir = 'initial_bgc/bgc.isam_initial.nc'

        repeat_cycles = 1
        start_time = ${starttime}
        run_time = ${runtime}

        offline_clim_data = 'CRU_NCEP'
        luc_data = 1

        clim_disturb= .true.
        clim_random= .false.
        co2_disturb= .true.
        c14_disturb= .true.
        luc_disturb= .false.
        n_disturb= .true.
        fna_func = .true.

        region_mask = 1,1,1,1,1,1,1,1,1,0,0

        datadir   = 'data'
        outputdir = 'output'

        single_point = .true.
        single_x = ${lonid}
        single_y = ${latid}
        generate_site_met_nc_from_global = .false.

        single_pft = .true.
        single_pft_num = -9999
        use_site_met = .false.
        sitename   = '${CASENAME}'
        site_ref_ht = 30
        dp_drain = 1     ! Topographic impedance for lateral drainage
        tao = 8      ! Impedance for calculating cyroturbation rate
        site_met_option = 'CRU_NCEP'
        site_met_timezone = 'local'

        dyn_veg         = .true.
        dyn_soc         = .false.   ! Switch of the coupled SOC content - thermal properties dynamics
        coupledsoc      = .false.   ! Switch of the coupled SOC content - thermal properties dynamics
        projection      = .false.   ! Apply the RCP8.5 projected atm CO2 concentration
        tsep            = .false.   ! (Under projection) Only consider SOC's thermal impact
        wsep            = .false.   ! (Under projection) Only consider SOC's hydraulic impact
        soc_saved       = .false.   ! Read saved soc from file
        fixed_harvest   = .false.
        fixed_planttime = .false.
        crop_mode       = 'generic'

        param_from_ascii = .false.
        soil_param_file = 'soil_param.dat'  ! Read soil texture from file
        soil_data = 'Combined'    ! Combined / GSDE

        hist_freq_yr = 60
        save_hist_yr = .true.
        save_hist_mon = .false.
        save_hist_bgp_to_bgc = .false.

        restart_freq_yr = 60
        site_restart = .true.

        dt_ch4      = 3600
        isam_deltim = 3600
/
EOF1

# Projection stage
###         cat > namelist.proj << EOF1
###&isamcfg
###
###        runname = '${CASENAME}'
###
###        bgp_bgc_mode = 3
###        restart = .false.
###        run_type = 1
###        CN_RATIO = 0
###        initial_bgp_in_datadir = 'initial_bgp/bgp.isam_initial_free.nc'
###        initial_bgc_in_datadir = 'initial_bgc/bgc.isam_initial_free.nc'
###
###        repeat_cycles = 10
###        start_time = ${starttime}
###        run_time = ${runtime}
###
###        offline_clim_data = 'CESM_PROJ'
###        luc_data = 1
###
###        clim_disturb= .true.
###        clim_random= .false.
###        co2_disturb= .false.
###        luc_disturb= .false.
###        n_disturb= .false.
###        fna_func = .true.
###
###        region_mask = 1,1,1,1,1,1,1,1,1,0,0
###
###        datadir   = 'data'
###        outputdir = 'output'
###
###        single_point = .true.
###        single_x = ${lonid}
###        single_y = ${latid}
###        generate_site_met_nc_from_global = .false.
###
###        single_pft = .true.
###        single_pft_num = -9999
###        use_site_met = .true.
###        sitename   = '${CASENAME}'
###        site_ref_ht = 30
###        site_met_option = 'CESM_PROJ'
###        site_met_timezone = 'local'
###
###        dyn_veg         = .true.
###        dyn_soc         = .true.
###        forced_soc      = .false.
###        projection      = .true.
###        fixed_harvest   = .false.
###        fixed_planttime = .false.
###        crop_mode       = 'generic'
###
###        hist_freq_yr = 6000
###        save_hist_yr = .true.
###        save_hist_mon = .false.
###        save_hist_bgp_to_bgc = .false.
###
###        restart_freq_yr = 6000
###
###        isam_deltim = 3600
###/
###EOF1
### 
###         cat > namelist.fixedproj << EOF1
###&isamcfg
###
###        runname = '${CASENAME}'
###
###        bgp_bgc_mode = 3
###        restart = .false.
###        run_type = 1
###        CN_RATIO = 0
###        initial_bgp_in_datadir = 'initial_bgp/bgp.isam_initial_free.nc'
###        initial_bgc_in_datadir = 'initial_bgc/bgc.isam_initial_free.nc'
###
###        repeat_cycles = 10
###        start_time = ${starttime}
###        run_time = ${runtime}
###
###        offline_clim_data = 'CESM_PROJ'
###        luc_data = 1
###
###        clim_disturb= .true.
###        clim_random= .false.
###        co2_disturb= .false.
###        luc_disturb= .false.
###        n_disturb= .false.
###        fna_func = .true.
###
###        region_mask = 1,1,1,1,1,1,1,1,1,0,0
###
###        datadir   = 'data'
###        outputdir = 'output'
###
###        single_point = .true.
###        single_x = ${lonid}
###        single_y = ${latid}
###        generate_site_met_nc_from_global = .false.
###
###        single_pft = .true.
###        single_pft_num = -9999
###        use_site_met = .true.
###        sitename   = '${CASENAME}'
###        site_ref_ht = 30
###        site_met_option = 'CESM_PROJ'
###        site_met_timezone = 'local'
###
###        coupledsoc      = .false.
###        dyn_veg         = .true.
###        dyn_soc         = .true.
###        forced_soc      = .false.
###        projection      = .true.
###        fixed_harvest   = .false.
###        fixed_planttime = .false.
###        crop_mode       = 'generic'
###
###        hist_freq_yr = 6000
###        save_hist_yr = .true.
###        save_hist_mon = .false.
###        save_hist_bgp_to_bgc = .false.
###
###        restart_freq_yr = 6000
###
###        isam_deltim = 3600
###/
###EOF1
 
# T-W Seperation test Projection stage
###         cat > namelist.tproj << EOF1
###&isamcfg
###
###        runname = '${CASENAME}'
###
###        bgp_bgc_mode = 3
###        restart = .false.
###        run_type = 1
###        CN_RATIO = 0
###        initial_bgp_in_datadir = 'initial_bgp/bgp.isam_initial_free.nc'
###        initial_bgc_in_datadir = 'initial_bgc/bgc.isam_initial_free.nc'
###
###        repeat_cycles = 10
###        start_time = ${starttime}
###        run_time = ${runtime}
###
###        offline_clim_data = 'CESM_PROJ'
###        luc_data = 1
###
###        clim_disturb= .true.
###        clim_random= .false.
###        co2_disturb= .false.
###        luc_disturb= .false.
###        n_disturb= .false.
###        fna_func = .true.
###
###        region_mask = 1,1,1,1,1,1,1,1,1,0,0
###
###        datadir   = 'data'
###        outputdir = 'output'
###
###        single_point = .true.
###        single_x = ${lonid}
###        single_y = ${latid}
###        generate_site_met_nc_from_global = .false.
###
###        single_pft = .true.
###        single_pft_num = -9999
###        use_site_met = .true.
###        sitename   = '${CASENAME}'
###        site_ref_ht = 30
###        site_met_option = 'CESM_PROJ'
###        site_met_timezone = 'local'
###
###        wsep            = .false.
###        dyn_veg         = .true.
###        dyn_soc         = .true.
###        forced_soc      = .false.
###        projection      = .true.
###        fixed_harvest   = .false.
###        fixed_planttime = .false.
###        crop_mode       = 'generic'
###
###        hist_freq_yr = 6000
###        save_hist_yr = .true.
###        save_hist_mon = .false.
###        save_hist_bgp_to_bgc = .false.
###
###        restart_freq_yr = 6000
###
###        isam_deltim = 3600
###/
###EOF1
### 
###         cat > namelist.tfixedproj << EOF1
###&isamcfg
###
###        runname = '${CASENAME}'
###
###        bgp_bgc_mode = 3
###        restart = .false.
###        run_type = 1
###        CN_RATIO = 0
###        initial_bgp_in_datadir = 'initial_bgp/bgp.isam_initial_free.nc'
###        initial_bgc_in_datadir = 'initial_bgc/bgc.isam_initial_free.nc'
###
###        repeat_cycles = 10
###        start_time = ${starttime}
###        run_time = ${runtime}
###
###        offline_clim_data = 'CESM_PROJ'
###        luc_data = 1
###
###        clim_disturb= .true.
###        clim_random= .false.
###        co2_disturb= .false.
###        luc_disturb= .false.
###        n_disturb= .false.
###        fna_func = .true.
###
###        region_mask = 1,1,1,1,1,1,1,1,1,0,0
###
###        datadir   = 'data'
###        outputdir = 'output'
###
###        single_point = .true.
###        single_x = ${lonid}
###        single_y = ${latid}
###        generate_site_met_nc_from_global = .false.
###
###        single_pft = .true.
###        single_pft_num = -9999
###        use_site_met = .true.
###        sitename   = '${CASENAME}'
###        site_ref_ht = 30
###        site_met_option = 'CESM_PROJ'
###        site_met_timezone = 'local'
###
###        coupledsoc      = .false.
###        wsep            = .false.
###        dyn_veg         = .true.
###        dyn_soc         = .true.
###        forced_soc      = .false.
###        projection      = .true.
###        fixed_harvest   = .false.
###        fixed_planttime = .false.
###        crop_mode       = 'generic'
###
###        hist_freq_yr = 6000
###        save_hist_yr = .true.
###        save_hist_mon = .false.
###        save_hist_bgp_to_bgc = .false.
###
###        restart_freq_yr = 6000
###
###        isam_deltim = 3600
###/
###EOF1

      #----------------------------------------------------------------------------------
      # Create directory
      #----------------------------------------------------------------------------------
      if ( -d data) then
        echo 'Found, no need to create directory: data.'
      else
        echo 'Create data directory for case:' ${CASENAME}
        mkdir data
        #----------------------------------------------------------------------------------
        # Link data from data repository and prepare input data
        #----------------------------------------------------------------------------------
        ln -s ${datarepo}/* data/
        rm data/initial_bgp
        rm data/initial_bgc
        mkdir data/initial_bgp
        mkdir data/initial_bgc
        mkdir data/SITE_MET_OUTPUT_FROM_GLOBAL
      endif
 
      if ( -d output) then
        echo 'Found, no need to create directory: output.'
      else
        echo 'Create output directory for case:' ${CASENAME}
        mkdir output
      endif
      #----------------------------------------------------------------------------------
      # Copy ISAM executable program
      #----------------------------------------------------------------------------------
      cp ${isamdir}/isam .
      #----------------------------------------------------------------------------------
      # Create Job script and submit the case to let it run
      #----------------------------------------------------------------------------------
      echo 'Generate batch script for case ' ${CASENAME}
      cat > genmet_${CASENAME}.sh << EOF1
#!/bin/sh

ulimit -s unlimited
ulimit -s
ulimit -c 10000
hostname

cd `pwd`

##------------------------------------
## NORMAL RUN
##------------------------------------
./isam < namelist.genmet > isam.log

##------------------------------------
## DEBUG RUN
##------------------------------------
##valgrind --tool=memcheck --leak-check=no --num-callers=20 --undef-value-errors=yes --track-origins=yes --read-var-info=yes --smc-check=all ./isam < namelist 2>&1 >isam.log

EOF1

         cat > spinup.sh << EOF1
#!/bin/sh 
#    -S /bin/sh 
#SBATCH -n 12 
#SBATCH --mem-per-cpu=5gb
#SBATCH --time=96:00:00 
#SBATCH --mail-type=FAIL 
#SBATCH --mail-type=END 
#SBATCH --mail-user=sshu3@illinois.edu 

cd `pwd`

##------------------------------------
## NORMAL RUN
##------------------------------------
./isam < namelist.spinup > isam.log

##------------------------------------
## DEBUG RUN
##------------------------------------
##valgrind --tool=memcheck --leak-check=no --num-callers=20 --undef-value-errors=yes --track-origins=yes --read-var-info=yes --smc-check=all ./isam < namelist 2>&1 >isam.log

EOF1

      #----------------------------------------------------------------------------------
      # Submit the job script to generate met forcing from CRUNCEP
      #----------------------------------------------------------------------------------
      # sbatch genmet_${CASENAME}.sh
      mv ./output/*forcing-CRU_NCEP.nc ./data/SITE_MET_OUTPUT_FROM_GLOBAL/
 
      #----------------------------------------------------------------------------------
      # Create a list of cases for further purpose
      #----------------------------------------------------------------------------------
##       cat >> ${workdir}/caselist << EOF1
##`pwd`
##EOF1

      #----------------------------------------------------------------------------------
      # Go back to upper level
      #----------------------------------------------------------------------------------
      echo 'Done ' ${CASENAME} '!'
      cd $workdir
      @ j = $j + 1

   end  # End of Loop inside case directory
   
   echo '################ ISAM Cases Setup Done ###############'

else

   echo '################ Skip ISAM Cases Setup ###############'

endif

if($runcases == 'true') then

   if( -f caselist ) then
      set cases=`cat caselist`
      set pft=`cat pftlist`
      set j=1
      while ( $j <= $#cases)
         echo 'Now run case:' $cases[$j]
         if ( -d $cases[$j]) then
           echo 'Found, submit the spinup job.'
           cd $cases[$j]
         cat > spinup.sh << EOF1
#!/bin/sh 
#    -S /bin/sh 
#SBATCH -n 1
#  #SBATCH --mem-per-cpu=5gb
#SBATCH --time=96:00:00 
#SBATCH --mail-type=FAIL 
#SBATCH --mail-type=END 
#SBATCH --mail-user=sshu3@illinois.edu 
#SBATCH -p g,f,d

cd `pwd`

##------------------------------------
## NORMAL RUN
##------------------------------------
./isam < namelist.spinup > isam.log

##------------------------------------
## DEBUG RUN
##------------------------------------
##valgrind --tool=memcheck --leak-check=no --num-callers=20 --undef-value-errors=yes --track-origins=yes --read-var-info=yes --smc-check=all ./isam < namelist 2>&1 >isam.log

EOF1

         cat > historical.sh << EOF1
#!/bin/sh 
#    -S /bin/sh 
#SBATCH -n 1
#  #SBATCH --mem-per-cpu=5gb
#SBATCH --time=96:00:00 
#SBATCH --mail-type=FAIL 
#SBATCH --mail-type=END 
#SBATCH --mail-user=sshu3@illinois.edu 
#SBATCH -p g,f,d

cd `pwd`

##------------------------------------
## NORMAL RUN
##------------------------------------
./isam < namelist.hist > isam.log

##------------------------------------
## DEBUG RUN
##------------------------------------
##valgrind --tool=memcheck --leak-check=no --num-callers=20 --undef-value-errors=yes --track-origins=yes --read-var-info=yes --smc-check=all ./isam < namelist 2>&1 >isam.log

EOF1

         cat > wproj.sh << EOF1
#!/bin/sh 
#    -S /bin/sh 
#SBATCH -n 1
#  #SBATCH --mem-per-cpu=5gb
#SBATCH --time=96:00:00 
#SBATCH --mail-type=FAIL 
#SBATCH --mail-type=END 
#SBATCH --mail-user=sshu3@illinois.edu 
#SBATCH -p g,f,d

cd `pwd`

##------------------------------------
## NORMAL RUN
##------------------------------------
./isam < namelist.wproj > isam.log

##------------------------------------
## DEBUG RUN
##------------------------------------
##valgrind --tool=memcheck --leak-check=no --num-callers=20 --undef-value-errors=yes --track-origins=yes --read-var-info=yes --smc-check=all ./isam < namelist 2>&1 >isam.log

EOF1

         cat > wfixedproj.sh << EOF1
#!/bin/sh 
#    -S /bin/sh 
#SBATCH -n 1
#  #SBATCH --mem-per-cpu=5gb
#SBATCH --time=96:00:00 
#SBATCH --mail-type=FAIL 
#SBATCH --mail-type=END 
#SBATCH --mail-user=sshu3@illinois.edu 
#SBATCH -p g,f,d

cd `pwd`

##------------------------------------
## NORMAL RUN
##------------------------------------
./isam < namelist.wfixedproj > isam.log

##------------------------------------
## DEBUG RUN
##------------------------------------
##valgrind --tool=memcheck --leak-check=no --num-callers=20 --undef-value-errors=yes --track-origins=yes --read-var-info=yes --smc-check=all ./isam < namelist 2>&1 >isam.log

EOF1

           cp $isamdir/isam .

           ## If you have records from previous run
           ## Please uncomment here to store these records
           # mv output/*forcing-CRU_NCEP.nc data/SITE_MET_OUTPUT_FROM_GLOBAL
           # echo 'Warning! Clean all txt files!'
            mkdir spinup
            mv *.txt spinup/
            mv isam.log spinup/
           # mv 1st/ old_case/
           # mv *.txt forced_proj_t
           # mv free_spin free_proj

           # Rename the latest initial file
            cd data/initial_bgp
            mv `ls -t1 | head -n 1` bgp.isam_initial.nc
           # mv bgp.isam_initial.nc bgp.isam_initial_forced.nc
            cd ../initial_bgc
            mv `ls -t1 | head -n 1` bgc.isam_initial.nc
           # mv bgc.isam_initial.nc bgc.isam_initial_forced.nc
           # ln -s ${workdir}/projected_climate/* data/SITE_MET_OUTPUT_FROM_GLOBAL/
           # ln -s ~/team/datasets4/projected_co2/ data/
           # ln -s ~/team/datasets4/projected_ndep/ data/

           # Back to the case directory
           cd ../..

           # If want to change namelist options please uncomment here
           # sed -i "s/hist_freq_yr\ =\ 12000/hist_freq_yr\ =\ 6000/" namelist.spinup
           # sed -i 's/restart_freq_yr\ =\ 12000/restart_freq_yr\ =\ 6000/' namelist.spinup
           sed -i "s/single_pft_num\ =\ -9999/single_pft_num\ =\ $pft[$j]/" namelist.hist

           sbatch historical.sh
           cd $workdir

         else
           echo 'Failed to find the directory. Check if setup was successful.'
         endif
         @ j = $j + 1
      end
   else
      exit 1
   endif

else

   echo '############# Do not run cases this time ############'

endif

if($postprocessing == 'true') then

   if( -f caselist ) then
      set cases=`cat caselist`
      set j=1
      while ( $j <= $#cases)
         echo 'Now processing case:' $cases[$j]
         if ( -d $cases[$j]) then
           echo 'Found, start post-processing.'
           cd $cases[$j]
           if ( -e socprof.txt) then
              tail -n 1 socprof.txt >> $workdir/$proffile
              tail -n 1 diffurate.txt >> $workdir/$difffile
              tail -n 1 ald.txt >> $workdir/$aldfile
              tail -n 1 deltac14.txt >> $workdir/$c14file
           else
              echo "-9999." >> $workdir/$proffile
              echo "-9999." >> $workdir/$difffile
              echo "-9999." >> $workdir/$aldfile
              echo "-9999." >> $workdir/$c14file
           endif
           # Purge the txts
           # rm *.txt
           cd $workdir
           ##awk '{ print $1","$2","$3","$4","$5}' texture_prof.txt >> $workdir/$outfile
         else
           echo 'Failed to find the directory. Check if setup was successful.'
         endif
         @ j = $j + 1
      end
   else
      exit 1
   endif

else

   echo '############# No post-processing will be executed ############'

endif


# Set the namelist for each variable

###(BASH) while read name
###do
###   echo "$name"
###done < list_histel
