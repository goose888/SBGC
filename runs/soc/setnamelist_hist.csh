#!/bin/csh

cat > namelist.hist << EOF1
&isamcfg

        runname = '${1}'

        bgp_bgc_mode = 3
        restart = .false.
        run_type = 2
        CN_RATIO = 0

        initial_bgp_in_datadir = 'initial_bgp/bgp.isam_initial.nc'
        initial_bgc_in_datadir = 'initial_bgc/bgc.isam_initial.nc'

        repeat_cycles = 1
        start_time = ${2}
        run_time = ${3}

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
        single_x = ${4}
        single_y = ${5}
        generate_site_met_nc_from_global = .false.

        single_pft = .true.
        single_pft_num = -9999
        use_site_met = .false.
        sitename   = '${1}'
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
