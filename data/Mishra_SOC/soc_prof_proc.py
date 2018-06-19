#!/bin/python

import numpy as np
import pandas as pd
import soc_analysis_lib as soca
import isamcalc_lib as isam
import socplot_lib as socplt
import auxiliary_lib as au

# 1. Open and read in original data
filename='N_circum_qaed.xlsx'
data = pd.read_excel(filename,index_col='Profile ID', skiprows=[1])
all_profid = data.index.unique()
lons = soca.getvarxls(data,'Long',all_profid,0)
lats = soca.getvarxls(data,'Lat N',all_profid,0)

## 2. Plot profile locations on the map
## Sort the lons to avoid the bug from matplot
#lon_s = lons[lons.argsort()]
#lat_s = lats[lons.argsort()]
#tit = 'Location of all SOC samples'
#path = 'Location_map.png'
#status = socplt.plot_profmap(-180, -60, 180, 80, lon_s, lat_s, tit, path)

## 3. Use Grass to get the permafrost status of each sites
## Store coordinates into ascii file, needed by GRASS
#fname = './allsites.asc'
#pp = np.column_stack((lons,lats))
#np.savetxt(fname, pp, fmt='%8.8f', delimiter=' ', newline='\n', header='', footer='')
#status = au.exec_shell('bash', './get_pf_status.sh', [])

# 4. Get the permafrost status of all sites and select cont/discont profiles only
pfst_table = np.loadtxt(fname='N_circum_sites_pf', dtype='float', delimiter=' ')
pfst = pfst_table[:,3].astype('int')
# Only select sites with continuous and discontinuous permafrost
# 1, 2, 5, 6, 9, 10, 13, 14, 17, 18.
lgc = np.ones((len(pfst)), dtype=bool)
for i in range(0, len(pfst)):
    lgc[i] = pfst[i] in [1, 2, 5, 6, 9, 10, 13, 14, 17, 18]
sel_profid = pfst_table[:,2].astype('int')[lgc]
data_sel = soca.subset_by_ids(data, sel_profid)
# Replace space with underscore
data_sel.columns = au.rep_sp_by_sym(data_sel, sym='_')

# 5. Extract vegetation information based on Foley pot veg map
sel_profid = data_sel.index.unique()
pfst_veg = np.loadtxt(fname='N_circum_sites_veg', dtype='float', delimiter=',')
site_veg = pfst_veg[:,3].astype('int')
# Transfer into pd dataframe
prof_veg = pd.DataFrame(np.stack((all_profid, site_veg), axis=1), index=all_profid, columns=['all_profid', 'veg_id'])
# Veg code based on foley pot veg map
sel_veg = soca.subset_by_ids(prof_veg, sel_profid)

## 6. Get the interpolated SOC profiles
## Extract data fields for mass preserving interpolation
## We only want Lat, Lon, basal depth, layer depth, soil order, 
## suborder, Horizon type, bulk density and C density 
#data_out = data_sel[['Veg_Class', 'Lat_N', 'Long', 'Basal_depth', 'Layer_depth', 'Soil_Order', \
#                                'Suborder', 'Horizon_type', 'bulk_density', 'C_Density']]
## Write out CSV file, needed by R script
#data_out.to_csv('./sel_sites_carbon.csv')
## Then call R script to run the mass-preseving spline interpolation
#status = au.exec_shell('Rscript', 'SOCinterp.r', [])
# Read in the interpolated SOC profile
filename = 'SOCprofile.csv'
soc_interped = pd.read_csv(filename,encoding='iso-8859-1',index_col=0)
pid = soc_interped.index
soc_interped = soc_interped * 1000.   # gcm-3 to kgm-3
soc_interped.index = sel_profid

## 7. Plot individual profile
#for i in range(len(sel_profid)):
#    X1 = soc_interped.loc[pid[i],:]   # SOC profile kgCm-3
#    Y = np.arange(1,201)     # 1cm to 200cm
#    if (data_sel.Profile_Name[sel_profid[i]].__class__.__name__ == 'Series'):
#        tit = data_sel.Profile_Name[sel_profid[i]].unique().astype('string')[0]
#    else:
#        tit = data_sel.Profile_Name[sel_profid[i]].encode('ascii','ignore')
#    path = './Figs/'+str(sel_profid[i])+'_'+tit+'.png'
#    status = socplt.plot_soilprof(X1, Y, tit, path)

# 8. Separate into different soil order + suborder
# Orthel
data_orthel, orthel_profid = soca.subset_by_cols(data_sel, 'Suborder', 'Orthel', True)
soc_orthel = soca.subset_by_ids(soc_interped, orthel_profid)
# Turbel
data_turbel, turbel_profid = soca.subset_by_cols(data_sel, 'Suborder', 'Turbel', True)
soc_turbel = soca.subset_by_ids(soc_interped, turbel_profid)
# Histel
data_histel, histel_profid = soca.subset_by_cols(data_sel, 'Suborder', 'Histel', True)
soc_histel = soca.subset_by_ids(soc_interped, histel_profid)
# Others
data_others, others_profid = soca.subset_by_cols(data_sel, 'Soil_Order', 'Gelisol', False)
soc_others = soca.subset_by_ids(soc_interped, others_profid)
# Derive mean and std.
m_orthel = soc_orthel.mean(axis=0)
s_orthel = soc_orthel.std(axis=0)
m_turbel = soc_turbel.mean(axis=0)
s_turbel = soc_turbel.std(axis=0)
m_histel = soc_histel.mean(axis=0)
s_histel = soc_histel.std(axis=0)
m_others = soc_others.mean(axis=0)
s_others = soc_others.std(axis=0)
## Make figure for each different soil order / suborder
#tit = 'Orthel'
#path = './Figs/'+tit+'.png'
#status = socplt.plot_prof_with_errbar(m_orthel[0:199], np.arange(1,200), s_orthel[0:199], tit, path)
#tit = 'Turbel'
#path = './Figs/'+tit+'.png'
#status = socplt.plot_prof_with_errbar(m_turbel[0:199], np.arange(1,200), s_turbel[0:199], tit, path)
#tit = 'Histel'
#path = './Figs/'+tit+'.png'
#status = socplt.plot_prof_with_errbar(m_histel[0:199], np.arange(1,200), s_histel[0:199], tit, path)
#tit = 'Others'
#path = './Figs/'+tit+'.png'
#status = socplt.plot_prof_with_errbar(m_others[0:199], np.arange(1,200), s_others[0:199], tit, path)

# 9. Read in model output and select sites overlapping with model estimation
# Read in the model output
socm = pd.read_table('isam_soc.dat', header=None, delimiter=r"\s+")
socm.columns = ['ID', 'Layer1', 'Layer2', 'Layer3', 'Layer4', 'Layer5', 'Layer6', 'Layer7', 'Layer8', 'Layer9', 'Layer10']
socm = socm.set_index('ID')
mod_profid = socm.index
# Unit conversion from kgCm-2 to kgCm-3
z, dz, zsoih = isam.get_isam_soildp(10)
socm = socm / dz
# Get the observed profiles by using IDs from model
soco_m = np.zeros(shape=(len(mod_profid),200) , dtype=float)
for i in range(len(mod_profid)):
    soco_m[i,:] = soc_interped[sel_profid == mod_profid[i]].astype(float)
# Transfer into pandas dataframe
soco = pd.DataFrame(soco_m, index=mod_profid, columns=soc_interped.columns)
## Generate figures
#for i in range(len(mod_profid)):
#    Xobs = soco.iloc[i,:]   # SOC profile kgCm-3
#    Yobs = np.arange(1,201)     # 1cm to 200cm
#    Xmod = socm.iloc[i,:]   # SOC profile kgCm-3
#    Ymod = z*100.     # 1cm to 200cm
#    if (data_sel.Profile_Name[mod_profid[i]].__class__.__name__ == 'Series'):
#        tit = data_sel.Profile_Name[mod_profid[i]].unique().astype('string')[0]
#    else:
#        tit = data_sel.Profile_Name[mod_profid[i]].encode('ascii','ignore')
#    path = './Figs_obsvsmod_calibrate/'+str(mod_profid[i])+'_'+tit+'.png'
#    status = socplt.plot_obsvsmod(Xobs, Yobs, Xmod, Ymod, tit, path)

# 10. Separate into different vegetation cover
# First, we can only examine the dataset with vegetation types being described
# Open the veg types file, veg types are coded in ISAM convention
filename = 'halfdeg_veg_gapfilled.csv'
vegcode_mod = pd.read_csv(filename,encoding='iso-8859-1',index_col=0)
# Only consider the samples with vegetation described in the original data
vegcode, veg_pid = soca.subset_by_cols(vegcode_mod, 'FROMSITE', 1, True)
# Grassland - 7 in ISAM
grass_site, grassid = soca.subset_by_cols(vegcode, 'VEGCODE', 7, True)
# Shrubland - 8 in ISAM
shrub_site, shrubid = soca.subset_by_cols(vegcode, 'VEGCODE', 8, True)
# Evergreen boreal forest - 5 in ISAM
everfor_site, everforid = soca.subset_by_cols(vegcode, 'VEGCODE', 5, True)
# Deciduous boreal forest - 20 in ISAM
decidfor_site, decidforid = soca.subset_by_cols(vegcode, 'VEGCODE', 20, True)
# Tundra - 9 in ISAM
tundra_site, tundraid = soca.subset_by_cols(vegcode, 'VEGCODE', 9, True)

soco_grass = soca.subset_by_ids(soco, grassid)
socm_grass = soca.subset_by_ids(socm, grassid)
soco_shrub = soca.subset_by_ids(soco, shrubid)
socm_shrub = soca.subset_by_ids(socm, shrubid)
soco_everfor = soca.subset_by_ids(soco, everforid)
socm_everfor = soca.subset_by_ids(socm, everforid)
soco_decidfor = soca.subset_by_ids(soco, decidforid)
socm_decidfor = soca.subset_by_ids(socm, decidforid)
soco_tundra = soca.subset_by_ids(soco, tundraid)
socm_tundra = soca.subset_by_ids(socm, tundraid)

m_grass_obs = soco_grass.mean(axis=0)
s_grass_obs = soco_grass.std(axis=0)
m_grass_mod = socm_grass.mean(axis=0)
s_grass_mod = socm_grass.std(axis=0)
m_shrub_obs = soco_shrub.mean(axis=0)
s_shrub_obs = soco_shrub.std(axis=0)
m_shrub_mod = socm_shrub.mean(axis=0)
s_shrub_mod = socm_shrub.std(axis=0)
m_everfor_obs = soco_everfor.mean(axis=0)
s_everfor_obs = soco_everfor.std(axis=0)
m_everfor_mod = socm_everfor.mean(axis=0)
s_everfor_mod = socm_everfor.std(axis=0)
m_decidfor_obs = soco_decidfor.mean(axis=0)
s_decidfor_obs = soco_decidfor.std(axis=0)
m_decidfor_mod = socm_decidfor.mean(axis=0)
s_decidfor_mod = socm_decidfor.std(axis=0)
m_tundra_obs = soco_tundra.mean(axis=0)
s_tundra_obs = soco_tundra.std(axis=0)
m_tundra_mod = socm_tundra.mean(axis=0)
s_tundra_mod = socm_tundra.std(axis=0)

## Create figures
#tit = 'Grassland'
#path = './Figs_obsvsmod/'+tit+'.png'
#status = socplt.plot_obsvsmod_with_errbar(m_grass_obs[0:199], np.arange(1,200), m_grass_mod, z*100, s_grass_obs[0:199], s_grass_mod, tit, path)
#tit = 'Shrubland'
#path = './Figs_obsvsmod/'+tit+'.png'
#status = socplt.plot_obsvsmod_with_errbar(m_shrub_obs[0:199], np.arange(1,200), m_shrub_mod, z*100, s_shrub_obs[0:199], s_shrub_mod, tit, path)
#tit = 'Evergreen_Forest'
#path = './Figs_obsvsmod/'+tit+'.png'
#status = socplt.plot_obsvsmod_with_errbar(m_everfor_obs[0:199], np.arange(1,200), m_everfor_mod, z*100, s_everfor_obs[0:199], s_everfor_mod, tit, path)
#tit = 'Deciduous_Forest'
#path = './Figs_obsvsmod/'+tit+'.png'
#status = socplt.plot_obsvsmod_with_errbar(m_decidfor_obs[0:199], np.arange(1,200), m_decidfor_mod, z*100, s_decidfor_obs[0:199], s_decidfor_mod, tit, path)
#tit = 'Tundra'
#path = './Figs_obsvsmod/'+tit+'.png'
#status = socplt.plot_obsvsmod_with_errbar(m_tundra_obs[0:199], np.arange(1,200), m_tundra_mod, z*100, s_tundra_obs[0:199], s_tundra_mod, tit, path)

# 11. Check the soc prediction under each veg type for Histel and Non-histel separately
#soco_tun_nonhis = soca.subset_by_ids(soco_tundra, orthel_profid)
soco_tun_hist = soca.subset_by_ids(soco_tundra, histel_profid)
soco_tun_orth = soca.subset_by_ids(soco_tundra, orthel_profid)
soco_tun_turb = soca.subset_by_ids(soco_tundra, turbel_profid)
socm_tun_hist = soca.subset_by_ids(socm_tundra, histel_profid)
socm_tun_orth = soca.subset_by_ids(socm_tundra, orthel_profid)
socm_tun_turb = soca.subset_by_ids(socm_tundra, turbel_profid)

soco_gra_hist = soca.subset_by_ids(soco_grass, histel_profid)
soco_gra_orth = soca.subset_by_ids(soco_grass, orthel_profid)
soco_gra_turb = soca.subset_by_ids(soco_grass, turbel_profid)
socm_gra_hist = soca.subset_by_ids(socm_grass, histel_profid)
socm_gra_orth = soca.subset_by_ids(socm_grass, orthel_profid)
socm_gra_turb = soca.subset_by_ids(socm_grass, turbel_profid)

soco_srb_hist = soca.subset_by_ids(soco_shrub, histel_profid)
soco_srb_orth = soca.subset_by_ids(soco_shrub, orthel_profid)
soco_srb_turb = soca.subset_by_ids(soco_shrub, turbel_profid)
socm_srb_hist = soca.subset_by_ids(socm_shrub, histel_profid)
socm_srb_orth = soca.subset_by_ids(socm_shrub, orthel_profid)
socm_srb_turb = soca.subset_by_ids(socm_shrub, turbel_profid)

soco_evf_hist = soca.subset_by_ids(soco_everfor, histel_profid)
soco_evf_orth = soca.subset_by_ids(soco_everfor, orthel_profid)
soco_evf_turb = soca.subset_by_ids(soco_everfor, turbel_profid)
socm_evf_hist = soca.subset_by_ids(socm_everfor, histel_profid)
socm_evf_orth = soca.subset_by_ids(socm_everfor, orthel_profid)
socm_evf_turb = soca.subset_by_ids(socm_everfor, turbel_profid)

soco_def_hist = soca.subset_by_ids(soco_decidfor, histel_profid)
soco_def_orth = soca.subset_by_ids(soco_decidfor, orthel_profid)
soco_def_turb = soca.subset_by_ids(soco_decidfor, turbel_profid)
socm_def_hist = soca.subset_by_ids(socm_decidfor, histel_profid)
socm_def_orth = soca.subset_by_ids(socm_decidfor, orthel_profid)
socm_def_turb = soca.subset_by_ids(socm_decidfor, turbel_profid)

# Get mean and std
m_gra_hist_obs = soco_gra_hist.mean(axis=0)
s_gra_hist_obs = soco_gra_hist.std(axis=0)
m_gra_hist_mod = socm_gra_hist.mean(axis=0)
s_gra_hist_mod = socm_gra_hist.std(axis=0)
m_gra_turb_obs = soco_gra_turb.mean(axis=0)
s_gra_turb_obs = soco_gra_turb.std(axis=0)
m_gra_turb_mod = socm_gra_turb.mean(axis=0)
s_gra_turb_mod = socm_gra_turb.std(axis=0)
m_gra_orth_obs = soco_gra_orth.mean(axis=0)
s_gra_orth_obs = soco_gra_orth.std(axis=0)
m_gra_orth_mod = socm_gra_orth.mean(axis=0)
s_gra_orth_mod = socm_gra_orth.std(axis=0)

m_srb_hist_obs = soco_srb_hist.mean(axis=0)
s_srb_hist_obs = soco_srb_hist.std(axis=0)
m_srb_hist_mod = socm_srb_hist.mean(axis=0)
s_srb_hist_mod = socm_srb_hist.std(axis=0)
m_srb_turb_obs = soco_srb_turb.mean(axis=0)
s_srb_turb_obs = soco_srb_turb.std(axis=0)
m_srb_turb_mod = socm_srb_turb.mean(axis=0)
s_srb_turb_mod = socm_srb_turb.std(axis=0)
m_srb_orth_obs = soco_srb_orth.mean(axis=0)
s_srb_orth_obs = soco_srb_orth.std(axis=0)
m_srb_orth_mod = socm_srb_orth.mean(axis=0)
s_srb_orth_mod = socm_srb_orth.std(axis=0)

m_tun_hist_obs = soco_tun_hist.mean(axis=0)
s_tun_hist_obs = soco_tun_hist.std(axis=0)
m_tun_hist_mod = socm_tun_hist.mean(axis=0)
s_tun_hist_mod = socm_tun_hist.std(axis=0)
m_tun_turb_obs = soco_tun_turb.mean(axis=0)
s_tun_turb_obs = soco_tun_turb.std(axis=0)
m_tun_turb_mod = socm_tun_turb.mean(axis=0)
s_tun_turb_mod = socm_tun_turb.std(axis=0)
m_tun_orth_obs = soco_tun_orth.mean(axis=0)
s_tun_orth_obs = soco_tun_orth.std(axis=0)
m_tun_orth_mod = socm_tun_orth.mean(axis=0)
s_tun_orth_mod = socm_tun_orth.std(axis=0)

m_evf_hist_obs = soco_evf_hist.mean(axis=0)
s_evf_hist_obs = soco_evf_hist.std(axis=0)
m_evf_hist_mod = socm_evf_hist.mean(axis=0)
s_evf_hist_mod = socm_evf_hist.std(axis=0)
m_evf_turb_obs = soco_evf_turb.mean(axis=0)
s_evf_turb_obs = soco_evf_turb.std(axis=0)
m_evf_turb_mod = socm_evf_turb.mean(axis=0)
s_evf_turb_mod = socm_evf_turb.std(axis=0)
m_evf_orth_obs = soco_evf_orth.mean(axis=0)
s_evf_orth_obs = soco_evf_orth.std(axis=0)
m_evf_orth_mod = socm_evf_orth.mean(axis=0)
s_evf_orth_mod = socm_evf_orth.std(axis=0)

m_def_hist_obs = soco_def_hist.mean(axis=0)
s_def_hist_obs = soco_def_hist.std(axis=0)
m_def_hist_mod = socm_def_hist.mean(axis=0)
s_def_hist_mod = socm_def_hist.std(axis=0)
m_def_turb_obs = soco_def_turb.mean(axis=0)
s_def_turb_obs = soco_def_turb.std(axis=0)
m_def_turb_mod = socm_def_turb.mean(axis=0)
s_def_turb_mod = socm_def_turb.std(axis=0)
m_def_orth_obs = soco_def_orth.mean(axis=0)
s_def_orth_obs = soco_def_orth.std(axis=0)
m_def_orth_mod = socm_def_orth.mean(axis=0)
s_def_orth_mod = socm_def_orth.std(axis=0)

## Create figures
#tit = 'Grassland_Orthel'
#path = './Figs_obsvsmod/'+tit+'.png'
#status = socplt.plot_obsvsmod_with_errbar(m_gra_orth_obs[0:199], np.arange(1,200), m_gra_orth_mod, z*100, s_gra_orth_obs[0:199], s_gra_orth_mod, tit, path)
#tit = 'Grassland_Turbel'
#path = './Figs_obsvsmod/'+tit+'.png'
#status = socplt.plot_obsvsmod_with_errbar(m_gra_turb_obs[0:199], np.arange(1,200), m_gra_turb_mod, z*100, s_gra_turb_obs[0:199], s_gra_turb_mod, tit, path)
#tit = 'Grassland_Histel'
#path = './Figs_obsvsmod/'+tit+'.png'
#status = socplt.plot_obsvsmod_with_errbar(m_gra_hist_obs[0:199], np.arange(1,200), m_gra_hist_mod, z*100, s_gra_hist_obs[0:199], s_gra_hist_mod, tit, path)
#
#tit = 'Shrubland_Orthel'
#path = './Figs_obsvsmod/'+tit+'.png'
#status = socplt.plot_obsvsmod_with_errbar(m_srb_orth_obs[0:199], np.arange(1,200), m_srb_orth_mod, z*100, s_srb_orth_obs[0:199], s_srb_orth_mod, tit, path)
#tit = 'Shrubland_Turbel'
#path = './Figs_obsvsmod/'+tit+'.png'
#status = socplt.plot_obsvsmod_with_errbar(m_srb_turb_obs[0:199], np.arange(1,200), m_srb_turb_mod, z*100, s_srb_turb_obs[0:199], s_srb_turb_mod, tit, path)
#tit = 'Shrubland_Histel'
#path = './Figs_obsvsmod/'+tit+'.png'
#status = socplt.plot_obsvsmod_with_errbar(m_srb_hist_obs[0:199], np.arange(1,200), m_srb_hist_mod, z*100, s_srb_hist_obs[0:199], s_srb_hist_mod, tit, path)
#
#tit = 'Evergreen_Forest_Orthel'
#path = './Figs_obsvsmod/'+tit+'.png'
#status = socplt.plot_obsvsmod_with_errbar(m_evf_orth_obs[0:199], np.arange(1,200), m_evf_orth_mod, z*100, s_evf_orth_obs[0:199], s_evf_orth_mod, tit, path)
#tit = 'Evergreen_Forest_Turbel'
#path = './Figs_obsvsmod/'+tit+'.png'
#status = socplt.plot_obsvsmod_with_errbar(m_evf_turb_obs[0:199], np.arange(1,200), m_evf_turb_mod, z*100, s_evf_turb_obs[0:199], s_evf_turb_mod, tit, path)
#tit = 'Evergreen_Forest_Histel'
#path = './Figs_obsvsmod/'+tit+'.png'
#status = socplt.plot_obsvsmod_with_errbar(m_evf_hist_obs[0:199], np.arange(1,200), m_evf_hist_mod, z*100, s_evf_hist_obs[0:199], s_evf_hist_mod, tit, path)
#
#tit = 'Deciduous_Forest_Orthel'
#path = './Figs_obsvsmod/'+tit+'.png'
#status = socplt.plot_obsvsmod_with_errbar(m_def_orth_obs[0:199], np.arange(1,200), m_def_orth_mod, z*100, s_def_orth_obs[0:199], s_def_orth_mod, tit, path)
#tit = 'Deciduous_Forest_Turbel'
#path = './Figs_obsvsmod/'+tit+'.png'
#status = socplt.plot_obsvsmod_with_errbar(m_def_turb_obs[0:199], np.arange(1,200), m_def_turb_mod, z*100, s_def_turb_obs[0:199], s_def_turb_mod, tit, path)
#tit = 'Deciduous_Forest_Histel'
#path = './Figs_obsvsmod/'+tit+'.png'
#status = socplt.plot_obsvsmod_with_errbar(m_def_hist_obs[0:199], np.arange(1,200), m_def_hist_mod, z*100, s_def_hist_obs[0:199], s_def_hist_mod, tit, path)
#
#tit = 'Tundra_Orthel'
#path = './Figs_obsvsmod/'+tit+'.png'
#status = socplt.plot_obsvsmod_with_errbar(m_tun_orth_obs[0:199], np.arange(1,200), m_tun_orth_mod, z*100, s_tun_orth_obs[0:199], s_tun_orth_mod, tit, path)
#tit = 'Tundra_Turbel'
#path = './Figs_obsvsmod/'+tit+'.png'
#status = socplt.plot_obsvsmod_with_errbar(m_tun_turb_obs[0:199], np.arange(1,200), m_tun_turb_mod, z*100, s_tun_turb_obs[0:199], s_tun_turb_mod, tit, path)
#tit = 'Tundra_Histel'
#path = './Figs_obsvsmod/'+tit+'.png'
#status = socplt.plot_obsvsmod_with_errbar(m_tun_hist_obs[0:199], np.arange(1,200), m_tun_hist_mod, z*100, s_tun_hist_obs[0:199], s_tun_hist_mod, tit, path)

# 12. Select profile ID for model calibration
## Find the profile that matches one profile the best
score, bid_gra_orth = soca.choose_best_prof(soco_gra_orth, m_gra_orth_obs)
score, bid_gra_turb = soca.choose_best_prof(soco_gra_turb, m_gra_turb_obs)
score, bid_gra_hist = soca.choose_best_prof(soco_gra_hist, m_gra_hist_obs)

score, bid_srb_orth = soca.choose_best_prof(soco_srb_orth, m_srb_orth_obs)
score, bid_srb_turb = soca.choose_best_prof(soco_srb_turb, m_srb_turb_obs)
score, bid_srb_hist = soca.choose_best_prof(soco_srb_hist, m_srb_hist_obs)

score, bid_evf_orth = soca.choose_best_prof(soco_evf_orth, m_evf_orth_obs)
score, bid_evf_turb = soca.choose_best_prof(soco_evf_turb, m_evf_turb_obs)
score, bid_evf_hist = soca.choose_best_prof(soco_evf_hist, m_evf_hist_obs)

score, bid_def_orth = soca.choose_best_prof(soco_def_orth, m_def_orth_obs)
score, bid_def_turb = soca.choose_best_prof(soco_def_turb, m_def_turb_obs)
score, bid_def_hist = soca.choose_best_prof(soco_def_hist, m_def_hist_obs)

score, bid_tun_orth = soca.choose_best_prof(soco_tun_orth, m_tun_orth_obs)
score, bid_tun_turb = soca.choose_best_prof(soco_tun_turb, m_tun_turb_obs)
score, bid_tun_hist = soca.choose_best_prof(soco_tun_hist, m_tun_hist_obs)

#13. Check all the citations and the profiles being related to these citations
data_selected_profs = soca.subset_by_ids(data, soco.index)
citations = data_selected_profs.citation.unique()
profs_sel_citation = list()
# Retreive the profile ID that from the same publication
for i in range(0,len(citations)):
# [0, 10, 225, 235, 240, 290, 295, 300, 325, 330, 340, 580, 1000, 1200, 1400, 1600, 1750, 1754, 1760, 1930, 1990, 2200, 2270, 2300, 2320]
    if(citations.iloc[i]==citations.iloc[0]):
        profs_sel_citation.append(int(citations.index[i]))
profs_sel_citation = np.unique(np.asarray(profs_sel_citation))
o=0
t=0
h=0
for i in range(0,len(profs_sel_citation)):
    if(data_selected_profs.loc[profs_sel_citation[i],:].Suborder.iloc[0].encode('utf8') == 'Orthel'):
        o = o + 1
    if(data_selected_profs.loc[profs_sel_citation[i],:].Suborder.iloc[0].encode('utf8') == 'Turbel'):
        t = t + 1
    if(data_selected_profs.loc[profs_sel_citation[i],:].Suborder.iloc[0].encode('utf8') == 'Histel'):
        h = h + 1

#lons = soca.getvarxls(data_sel,'Long',sel_profid,0)
#lats = soca.getvarxls(data_sel,'Lat_N',sel_profid,0)
#data_sel_meta=np.column_stack((sel_profid,lons,lats))
#fname='site_latlon.asc'
#np.savetxt(fname, data_sel_meta, fmt='%8.8f', delimiter=' ', newline='\n', header='', footer='')

