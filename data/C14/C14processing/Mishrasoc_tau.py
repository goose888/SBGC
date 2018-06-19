# -*- Python source file -*-
"""
File Name : Mishrasoc_tau.py

Description : Calculate the SOC stock and turnover time and create related plots for Mishra's SOC profiles

Created on : Mon Jun 18 23:34:49 2018

Last Modified : Mon Jun 18 23:41:08 2018

Author : Shijie Shu
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import C14preproc as prep
import pylab
import C14utils
import matplotlib
import socplot_lib as socplt
import isamcalc_lib as isam
import SOCtools as soc

# =====================================================================
#  Plot the simulated SOC profiles for Umakant's samples
#  Compare with the observation also
# =====================================================================
fobs = "Mishra_soc_interp.csv"
fobs_pid = "sel_sites_carbon.csv"
fmod_pid = "caselist"
fmod_biome = "pftlist"
fmod = "Mishra_soc_mod.dat"
fhistel_filter = "nohistel"

# Read in obs
# 354 samples in total
obs = pd.read_csv(fobs, encoding='iso-8859-1', index_col=0)
obs_pid = pd.read_csv(fobs_pid, encoding='iso-8859-1', index_col=0)
obs_profid = obs_pid.index.unique()
obs.index = obs_profid
# Read in the model output
# 255 samples in total
mod = pd.read_csv(fmod, delim_whitespace=True, header=None)
mod_pid = pd.read_csv(fmod_pid, header=None, index_col=0)
mod_biome = pd.read_csv(fmod_biome, header=None, index_col=0)
mod.index = mod_pid.index

# Read in the histel filter
histel_f = pd.read_csv(fhistel_filter, header=None, index_col=0)

sel_profid = mod_pid.index
obs_sel = obs.loc[sel_profid]

# Filter out the histel before the comparison
obs_sel.loc[histel_f.index] = -9999.
mod.loc[histel_f.index] = -9999.

# Convert the interpolated SOC into ISAM depth
z, dz, zsoih = isam.get_isam_soildp(10)
obs_as_isam = isam.mean_by_depth(10, zsoih, sel_profid.size, obs_sel.as_matrix())
obs_as_isam[obs_as_isam<=0]=float('nan')
# Unit convert: g/cm3 -> kg/m3
obs_isam = obs_as_isam * 1000.
obs_isam[obs_isam<0.] = float("nan")
# Aggregated SOC density to SOC stock profile
obs_agg_prof = soc.aggre_profden_to_profstock(10, dz, obs_isam)
obs_agg = soc.aggre_profden_to_stock(7, dz, obs_isam)

mod_isam = mod.as_matrix()/dz
mod_isam[mod_isam<0.] = float("nan")
mod_agg_prof = soc.aggre_profden_to_profstock(10, dz, mod_isam)
mod_agg = soc.aggre_profden_to_stock(7, dz, mod_isam)

# Filter out the observation if model did not calculate
obs_agg[np.isnan(mod_agg)] = float("nan")
obs_agg_prof[np.isnan(mod_agg_prof)] = float("nan")

# Prepare figures
# 1) The comparison of the total SOC stock
# do not separate into different biomes (all in one)
obs_agg_avg = np.nanmean(obs_agg)
mod_agg_avg = np.nanmean(mod_agg)
# Boxplot of the total SOC in the first 1m
# Shijie: Shall filter out peat before plotting!
data_for_bplot = [obs_agg[~np.isnan(obs_agg)], mod_agg[~np.isnan(mod_agg)]]

# multiple box plots on one figure
plt.figure()
plt.boxplot(data_for_bplot)
plt.xticks([1, 2], ['Obs', 'ISAM_1DSBGC'])
plt.savefig('SOC_stock.png')

# 2) The comparison of the accumulated SOC profile
obs_aggprof_avg = np.nanmean(obs_agg_prof, axis=0)
mod_aggprof_avg = np.nanmean(mod_agg_prof, axis=0)
obs_aggprof_std = np.nanstd(obs_isam*dz, axis=0)
mod_aggprof_std = np.nanstd(mod_isam*dz, axis=0)
obs_aggprof_avg[9] = float("nan")
obs_aggprof_std[9] = float("nan")
mod_aggprof_std[:] = 0.
obs_aggprof_std = obs_aggprof_std #/np.sqrt(146.)  # 146 is the sample size
mod_aggprof_std = mod_aggprof_std #/np.sqrt(146.)  # 146 is the sample size

soildp = z[0:8]*100
#tit = "Accumulated SOC profile"
tit = ""
path = "./SOC_profile.png"

status = socplt.plot_obsvsmod_with_errbar(obs_aggprof_avg[0:8], soildp, mod_aggprof_avg[0:8], soildp, obs_aggprof_std[0:8], mod_aggprof_std[0:8], tit, path)

# 3) The comparison of the accumulated SOC profile by different biomes
# 5 biomes in total:
# 7 - grassland
# 8 - shrub
# 9 - tundra
# 5 - boreal evergreen
# 20 - boreal deciduous
obs_agg_prof_pft = obs_agg_prof[mod_biome.index == 7]
mod_agg_prof_pft = 1.7*mod_agg_prof[mod_biome.index == 7]
obs_aggprof_avg = np.nanmean(obs_agg_prof_pft, axis=0)
mod_aggprof_avg = np.nanmean(mod_agg_prof_pft, axis=0)
obs_aggprof_std = np.nanstd(obs_isam[mod_biome.index == 7]*dz, axis=0)
mod_aggprof_std = np.nanstd(mod_isam[mod_biome.index == 7]*dz, axis=0)
obs_aggprof_avg[9] = float("nan")
obs_aggprof_std[9] = float("nan")
mod_aggprof_std[:] = 0.
obs_aggprof_std = obs_aggprof_std #/np.sqrt(146.)  # 146 is the sample size
mod_aggprof_std = mod_aggprof_std #/np.sqrt(146.)  # 146 is the sample size

soildp = z[0:8]*100
#tit = "Accumulated SOC profile"
tit = ""
path = "./SOC_grassland_profile.png"

status = socplt.plot_obsvsmod_with_errbar(obs_aggprof_avg[0:8], soildp, mod_aggprof_avg[0:8], soildp, obs_aggprof_std[0:8], mod_aggprof_std[0:8], tit, path)

# =====================================================================
#  Plot the simulated 1-m integrated D14C profiles for Umakant's samples
#  Compare with the others?
# =====================================================================
# ================================
## Recalculate the Tau
# ================================
fmod_pid = "caselist"
fmod_biome = "pftlist"
fmod = "Mishra_d14c_mod.dat"
fmod_soc = "Mishra_soc_mod.dat"
fhistel_filter = "nohistel"

# Read in D14C from observation
data = pd.read_csv(fmod, delim_whitespace=True, header=None)
pid = pd.read_csv(fmod_pid, header=None, index_col=0)
data.index = pid.index

# Read in the SOC profile
dsoc = pd.read_csv(fmod_soc, delim_whitespace=True, header=None)
dsoc.index = pid.index
soc = dsoc.as_matrix()
frac = soc * np.float("nan")
for i in np.arange(0,7):
    frac[:,i] = soc[:,i]/np.sum(soc[:,0:6], axis=1)

# Filter the Histel data
histel_f = pd.read_csv(fhistel_filter, header=None, index_col=0)
data.loc[histel_f.index] = -9999.

# Assign Nan value for the all the data with -9999.
data[data<-1200.] = float("nan")
# Get the C-weighted D14
d14c = np.nansum(frac * data, axis=1)
d14c[d14c==0.]=np.float("nan")

# Calculate tau
z, dz, zsoih = isam.get_isam_soildp(10)
# First get the C-weighted mean D14C for each profile
 
sampleyr = 2015 * np.ones(len(d14c))

tau, cost = C14utils.cal_tau(d14c, sampleyr, 1, 0)
tau[tau==2.00000000e+03] = np.float("nan")
data['tau'] = pd.Series(tau[:,0], index=data.index)
data.nodedepth = z*100
ttt=tau.reshape(255)
data_bplot = ttt[~np.isnan(ttt)]

# Plot boxplot for tau
plt.figure()
plt.boxplot(data_bplot)
plt.xticks([1], ['ISAM turnover'])
plt.savefig('tau_ISAM.png')

# Save the calculated tau into CSV table
data.to_csv("Mishra_tau_from_ISAM.csv")

# ===============================================
# Now check the turnover by different biome!
# ===============================================
fmod_biome = "pftlist"
biome = pd.read_csv(fmod_biome, header=None, index_col=0)
tt1 = ttt[biome.index == 5]
data_bplot = tt1[~np.isnan(tt1)]

plt.figure()
plt.boxplot(data_bplot)
plt.xticks([1], ['Boreal turnover'])
plt.savefig('tau_ISAM.png')

# =========================================================================
#  Calculate the tau using the simulated D14C of Umakant's samples
#  Check under different biome how are the simulated turnover look like?
# =========================================================================
filename = 'Non_peat_data_synthesis.csv'
data = pd.read_csv(filename,encoding='iso-8859-1',index_col='ProfileID')
biome = {1:'Boreal Forest',2:'Temperate Forest',3:'Tropical Forest',4:'Grassland', \
         5:'Cropland',6:'Shrublands',7:'Peatland',8:'Savannas',9:'Tundra',10:'Desert'}
var = ['D14C_BulkLayer','BulkDensity','pct_C','totalSOCgcm3', 'tau']
varlabel = [r"$\Delta14C$ ("+ u"\u2030)",
            r"Bulk Density $(g\ cm^{-3})$",
            r"percent C (%)",
            r"SOC content $(gC\ cm^{-3})$",
            r"turnover time (yr)"]
# customize
pltbiome = biome.keys()[8:10]
var1 = 1 # plot on the top x axis
var2 = 2 # plot on the bottom x axis

dum = 0
fig, axes = plt.subplots(nrows=1, ncols=len(set(pltbiome)), figsize=(12,8))
for i in set(pltbiome): # loop over biomes
    biomedataid = data[data.VegTypeCode_Local==i].index
    ax1 = fig.axes[dum] # bottom
    ax2 = ax1.twiny() # top
    cm = plt.get_cmap('Set1')
    numcolr = len(set(biomedataid)) # no repeat in color
    ax1.set_color_cycle([cm(1.*jj/numcolr) for jj in range(numcolr)])
    ax2.set_color_cycle([cm(1.*jj/numcolr) for jj in range(numcolr)])
    for kk in set(biomedataid): # loop over profiles in current biome
        Y = (data[data.index==kk]['Layer_top_norm'].values.astype(float)+\
            data[data.index==kk]['Layer_bottom_norm'].values.astype(float))/2.0
        X1 = data[data.index==kk][var[var1]].astype(float)
        if var2 == 3:
            X2 = data[data.index==kk]['BulkDensity'].astype(float) * \
                 data[data.index==kk]['pct_C'].astype(float)/100.0 # total SOC g/cm3
        else:
            X2 = data[data.index==kk][var[var2]].astype(float)
        h1 = ax1.plot(X2,Y,'-.',lw=3,label='%s,%s,%s' % (data[data.index==kk].Site.values[0],\
            data[data.index==kk].Country.values[0],data[data.index==kk].State_City.values[0]))
        ax1.set_xlabel(varlabel[var2]+'(dashed line)',color='g')
        ax1.set_ylabel(r"Depth $(cm)$")
        for tl in ax1.get_xticklabels():
                tl.set_color('g')
        h2 = ax2.plot(X1,Y)
        ax2.set_xlabel(varlabel[var1]+'(solid line)')
    plt.gca().invert_yaxis()
    pylab.text(0.8, 0.1,biome[i]+"\n"+"N = "+str(len(set(biomedataid))),
               horizontalalignment='center',verticalalignment='center',
               transform = ax1.transAxes,fontsize=16)
    dum = dum + 1
fig.tight_layout() # Or equivalently,  "plt.tight_layout()"
matplotlib.rcParams.update({'font.size': 10})
fig.savefig('./BD_SOC_%s_%s.png'%(biome[pltbiome[0]],biome[pltbiome[1]]))

