#!/bin/csh

cat > namelist.genmet << EOF1
&isamcfg

        runname = '${1}'

        bgp_bgc_mode = 1
        restart = .false.
        run_type = 0
        CN_RATIO = 0
        initial_bgp_in_datadir = 'initial_bgp/bgp.isam_initial.nc'
        initial_bgc_in_datadir = 'initial_bgc/bgc.isam_initial.nc'

        repeat_cycles = 1
        start_time = ${2}
        run_time = ${3}

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
        single_x = ${4}
        single_y = ${5}
        generate_site_met_nc_from_global = .true.

        single_pft = .false.
        single_pft_num = -999
        use_site_met = .false.
        sitename   = '${1}'
        site_ref_ht = -9999

        fixed_harvest   = .false.
        fixed_planttime = .false.
        crop_mode       = 'generic'

        hist_freq_yr = 200
        save_hist_yr = .true.
        save_hist_mon = .false.
        save_hist_bgp_to_bgc = .false.

        restart_freq_yr = 200

        isam_deltim = 3600
/
EOF1
