#!/bin/bash
grasscr Global Permafrost
v.in.ascii --overwrite input=allsites.asc output=N_circum_sites fs='space'
v.db.addtable map=N_circum_sites
v.db.addcolumn map=N_circum_sites column='PF_status INTEGER'
v.what.rast map=N_circum_sites raster=pf_extent_int column=PF_status
v.out.ascii --overwrite input=N_circum_sites output=N_circum_sites_pf separator='space' columns='PF_status'

