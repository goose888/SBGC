#!/bin/bash
grasscr Global global_lc
v.in.ascii --overwrite input=allsites.asc output=N_circum_sites fs='space'
v.db.addtable map=N_circum_sites
v.db.addcolumn map=N_circum_sites column='Veg_from_Foleymap INTEGER'
v.what.rast map=N_circum_sites raster=foley_potveg column=Veg_from_Foleymap
v.out.ascii --overwrite input=N_circum_sites output=N_circum_sites_veg columns='Veg_from_Foleymap' separator='comma'
