#!/bin/csh

cat > namelist.tproj << EOF1
&isamcfg

        runname = '${1}'

        bgp_bgc_mode = 3
        restart = .false.
        run_type = 1
        CN_RATIO = 0
        initial_bgp_in_datadir = 'initial_bgp/bgp.isam_initial_free.nc'
        initial_bgc_in_datadir = 'initial_bgc/bgc.isam_initial_free.nc'

        repeat_cycles = 10
        start_time = ${2}
        run_time = ${3}

        offline_clim_data = 'CESM_PROJ'
        luc_data = 1

        clim_disturb= .true.
        clim_random= .false.
        co2_disturb= .false.
        luc_disturb= .false.
        n_disturb= .false.
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
        use_site_met = .true.
        sitename   = '${1}'
        site_ref_ht = 30
        site_met_option = 'CESM_PROJ'
        site_met_timezone = 'local'

        coupledsoc      = .true.
        tsep            = .true.
        wsep            = .false.
        dyn_veg         = .true.
        dyn_soc         = .true.
        forced_soc      = .false.
        projection      = .true.
        fixed_harvest   = .false.
        fixed_planttime = .false.
        crop_mode       = 'generic'

        hist_freq_yr = 6000
        save_hist_yr = .true.
        save_hist_mon = .false.
        save_hist_bgp_to_bgc = .false.

        restart_freq_yr = 6000

        isam_deltim = 3600
/
EOF1
