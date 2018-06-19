# -*- Python source file -*-
"""
File Name : SOC.py

Description : Codes processing SOC profile from He et al., 2016

Created on : Tue Jun 19 00:11:12 2018

Last Modified : Tue Jun 19 00:12:26 2018

Author : Shijie Shu
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import C14preproc as prep
from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset
from scipy.interpolate import interp1d
import auxiliary_lib as au
import isamcalc_lib as isam
import socplot_lib as socplt
import time

#=========================================================================
## Plot simulated SOC results and compare with the original observation
#=========================================================================
cutdep = 100.
filename = 'Non_peat_data_permafrost.csv'
data = pd.read_csv(filename,encoding='iso-8859-1',index_col='ProfileID', skiprows=[1])
all_profid = data.index.unique()
lons = prep.getvarxls(data,'Lon',all_profid[0:16],0)
lats = prep.getvarxls(data,'Lat',all_profid[0:16],0)

# Original SOC profile
soc_obs = data.C_Density

# Read in model outputs
soccm = pd.read_table('isam_soc_cali.dat', header=None, delimiter=r"\s+")
soccm.columns = ['ID', 'Layer1', 'Layer2', 'Layer3', 'Layer4', 'Layer5', 'Layer6', 'Layer7', 'Layer8', 'Layer9', 'Layer10']
soccm = soccm.set_index('ID')
mod_profid = soccm.index

# Get the corresponding SOC comparable against ISAM
z, dz, zsoih = isam.get_isam_soildp(10)
zsoih = zsoih * 100

# Extract corresponding measurements for each site and make figure
data.nodedepth = data.Layer_top + (data.Layer_bottom - data.Layer_top)/2.
# myarray = soccm.index.unique()
myarray = [43, 110, 143, 146, 197]
# myarray = [143] 

# Omit site 155 since there's only one value for that profile
for i in myarray:
    # Observation
    dn = data.nodedepth.loc[i].as_matrix()
    dp = (data.Layer_bottom.loc[i] - data.Layer_top.loc[i]).as_matrix()
    num_obs = soc_obs[i].__len__()
    accsoc = soc_obs[i].as_matrix() * float('nan')
    sitesoc = soc_obs.loc[i].as_matrix()    
    accsoc[0] = sitesoc[0] * 1000. * dp[0] / 100.
    accsoc[-1] = float('nan')
    for j in range(1,num_obs):
        accsoc[j] = accsoc[j-1] + sitesoc[j] * 1000. * dp[j] / 100.
    # Model results
    modsoc = soccm.loc[i].as_matrix()
    for j in range(1,10):
        modsoc[j] = modsoc[j-1] + modsoc[j]
 
    # Make plot
    Xobs = accsoc   # SOC profile kgCm-3
    Yobs = dn       # 1cm to 200cm
    Xmod = modsoc   # SOC profile kgCm-3
    Ymod = z*100.     # 1cm to 200cm
    # Xmod2 = modsoc   # SOC profile kgCm-3
    # Ymod2 = z*100.     # 1cm to 200cm
    if (data.Site[i].__class__.__name__ == 'Series'):
        tit = data.Site[i].unique().astype('string')[0]
    else:
        tit = data.Site[i].encode('ascii','ignore')
    path = './Figs_obsvsmod_original_soc/'+str(i)+'_'+tit+'_soc.png'
    xticks = (0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
    #status = socplt.plot_obsvsmod(Xobs, Yobs, Xmod, Ymod, tit, path, xticks)
    status = socplt.plot_obsvsmod(Xobs, Yobs, Xmod, Ymod, tit, path, None, None, xticks, 24, 30, 'upper right')

#===================================================================
## Plot simulated SOC results and compare with the interpolated SOC
#===================================================================

filename = 'Non_peat_data_permafrost.csv'
data = pd.read_csv(filename,encoding='iso-8859-1',index_col='ProfileID', skiprows=[1])
all_profid = data.index.unique()
lons = prep.getvarxls(data,'Lon',all_profid[0:16],0)
lats = prep.getvarxls(data,'Lat',all_profid[0:16],0)

# Read in model outputs
soccm = pd.read_table('isam_soc.dat', header=None, delimiter=r"\s+")
soccm.columns = ['ID', 'Layer1', 'Layer2', 'Layer3', 'Layer4', 'Layer5', 'Layer6', 'Layer7', 'Layer8', 'Layer9', 'Layer10']
soccm = soccm.set_index('ID')
mod_profid = soccm.index
# Read in another case from the ISAM model output if required
# soccm2 = pd.read_table('isam_soc_wt_ohorizon.dat', header=None, delimiter=r"\s+")
# soccm2.columns = ['ID', 'Layer1', 'Layer2', 'Layer3', 'Layer4', 'Layer5', 'Layer6', 'Layer7', 'Layer8', 'Layer9', 'Layer10']
# soccm2 = soccm2.set_index('ID')
# mod2_profid = soccm2.index

# Get the model depth
z, dz, zsoih = isam.get_isam_soildp(10)

## Shall use the mass perserving slpine (Bishop et al., 1999, Geoderma) to interpolate the bulk density data
## Okay, seems like only R has the offical package that can perform this interpolation
## so we shall use the R code to get this done instead.
#data_out = data[['Veg_Class', 'Lat', 'Lon', 'Basal_depth', 'Layer_depth', 'Soil_Order', \
#                                'Suborder', 'Horizon_type', 'bulk_density', 'C_Density']]
#data_out.to_csv('./sel_sites_soc.csv')
## Then call R script to run the mass-preseving spline interpolation
#status = au.exec_shell('Rscript', 'SOCinterp.r', [])

# Read in the interpolated SOC profile
sel_profid = data.index.unique()
filename = 'SOCprofile.csv'
soc_interped = pd.read_csv(filename,encoding='iso-8859-1',index_col=0)
pid = soc_interped.index
soc_interped.index = sel_profid

# Get the corresponding SOC comparable against ISAM
z, dz, zsoih = isam.get_isam_soildp(10)
zsoih = zsoih * 100

# Extract corresponding measurements for each site and make figure
data.nodedepth = data.Layer_top + (data.Layer_bottom - data.Layer_top)/2.
myarray = soccm.index.unique()
for i in myarray:
    accsoc = dz * float('nan')
    sitesoc = soc_interped.loc[i].as_matrix()
    modsoc = soccm.loc[i].as_matrix()
    accsoc[0] = 1000 * np.nanmean(sitesoc[0]) * dz[0]
    accsoc[9] = float('nan')
    for j in range(1,7):
        if (np.isnan(sitesoc[np.ceil(zsoih[j+1])])):
            accsoc[j] = float('nan')
        else:
            lb = np.floor(zsoih[j])
            ub = np.ceil(zsoih[j+1])
            accsoc[j] = accsoc[j-1] + 1000 * np.nanmean(sitesoc[lb:ub]) * dz[j]
    for j in range(1,10):
        modsoc[j] = modsoc[j-1] + modsoc[j]
 
    Xobs = accsoc   # SOC profile kgCm-3
    Yobs = z*100   # 1cm to 200cm
    Xmod = modsoc   # SOC profile kgCm-3
    Ymod = z*100.     # 1cm to 200cm
    # Xmod2 = modsoc   # SOC profile kgCm-3
    # Ymod2 = z*100.     # 1cm to 200cm
    if (data.Site[i].__class__.__name__ == 'Series'):
        tit = data.Site[i].unique().astype('string')[0]
    else:
        tit = data.Site[i].encode('ascii','ignore')
    path = './Figs_obsvsmod_calibrate_soc/'+str(i)+'_'+tit+'_soc.png'
    xticks = (0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
    #status = socplt.plot_obsvsmod(Xobs, Yobs, Xmod, Ymod, tit, path, xticks)
    status = socplt.plot_obsvsmod(Xobs, Yobs, Xmod, Ymod, tit, path, None, None, xticks, 24, 30, 'upper right')


#===================================================================
## Plot all simulated SOC profile
#===================================================================

filename = 'socprofeq.dat'
data = pd.read_csv(filename,encoding='iso-8859-1',index_col='ProfileID', skiprows=[1])
all_profid = data.index.unique()
lons = prep.getvarxls(data,'Lon',all_profid[0:16],0)
lats = prep.getvarxls(data,'Lat',all_profid[0:16],0)

# Read in model outputs
soccm = pd.read_table('isam_soc.dat', header=None, delimiter=r"\s+")
soccm.columns = ['ID', 'Layer1', 'Layer2', 'Layer3', 'Layer4', 'Layer5', 'Layer6', 'Layer7', 'Layer8', 'Layer9', 'Layer10']
soccm = soccm.set_index('ID')
mod_profid = soccm.index
# Read in another case from the ISAM model output if required
# soccm2 = pd.read_table('isam_soc_wt_ohorizon.dat', header=None, delimiter=r"\s+")
# soccm2.columns = ['ID', 'Layer1', 'Layer2', 'Layer3', 'Layer4', 'Layer5', 'Layer6', 'Layer7', 'Layer8', 'Layer9', 'Layer10']
# soccm2 = soccm2.set_index('ID')
# mod2_profid = soccm2.index

# Get the model depth
z, dz, zsoih = isam.get_isam_soildp(10)

## Shall use the mass perserving slpine (Bishop et al., 1999, Geoderma) to interpolate the bulk density data
## Okay, seems like only R has the offical package that can perform this interpolation
## so we shall use the R code to get this done instead.
#data_out = data[['Veg_Class', 'Lat', 'Lon', 'Basal_depth', 'Layer_depth', 'Soil_Order', \
#                                'Suborder', 'Horizon_type', 'bulk_density', 'C_Density']]
#data_out.to_csv('./sel_sites_soc.csv')
## Then call R script to run the mass-preseving spline interpolation
#status = au.exec_shell('Rscript', 'SOCinterp.r', [])

# Read in the interpolated SOC profile
sel_profid = data.index.unique()
filename = 'SOCprofile.csv'
soc_interped = pd.read_csv(filename,encoding='iso-8859-1',index_col=0)
pid = soc_interped.index
soc_interped.index = sel_profid

# Get the corresponding SOC comparable against ISAM
z, dz, zsoih = isam.get_isam_soildp(10)
zsoih = zsoih * 100

# Extract corresponding measurements for each site and make figure
data.nodedepth = data.Layer_top + (data.Layer_bottom - data.Layer_top)/2.
myarray = soccm.index.unique()
for i in myarray:
    accsoc = dz * float('nan')
    sitesoc = soc_interped.loc[i].as_matrix()
    modsoc = soccm.loc[i].as_matrix()
    accsoc[0] = 1000 * np.nanmean(sitesoc[0]) * dz[0]
    accsoc[9] = float('nan')
    for j in range(1,7):
        if (np.isnan(sitesoc[np.ceil(zsoih[j+1])])):
            accsoc[j] = float('nan')
        else:
            lb = np.floor(zsoih[j])
            ub = np.ceil(zsoih[j+1])
            accsoc[j] = accsoc[j-1] + 1000 * np.nanmean(sitesoc[lb:ub]) * dz[j]
    for j in range(1,10):
        modsoc[j] = modsoc[j-1] + modsoc[j]
 
    Xobs = accsoc   # SOC profile kgCm-3
    Yobs = z*100   # 1cm to 200cm
    Xmod = modsoc   # SOC profile kgCm-3
    Ymod = z*100.     # 1cm to 200cm
    # Xmod2 = modsoc   # SOC profile kgCm-3
    # Ymod2 = z*100.     # 1cm to 200cm
    if (data.Site[i].__class__.__name__ == 'Series'):
        tit = data.Site[i].unique().astype('string')[0]
    else:
        tit = data.Site[i].encode('ascii','ignore')
    path = './Figs_obsvsmod_calibrate_soc/'+str(i)+'_'+tit+'_soc.png'
    xticks = (0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
    #status = socplt.plot_obsvsmod(Xobs, Yobs, Xmod, Ymod, tit, path, xticks)
    status = socplt.plot_obsvsmod(Xobs, Yobs, Xmod, Ymod, tit, path, None, None, xticks, 24, 30, 'upper right')

#===================================================================
## Plot all simulated SOC profile from Umakant
#===================================================================

filename = 'socprofeq.dat'
data = pd.read_csv(filename,encoding='iso-8859-1', skiprows=[1])
prof = np.zeros([231,10])
l = 0
with open(filename, 'r') as fobj:
    for line in fobj:
        numbers = [num for num in line.split()]
        # do something with this line of numbers before moving on to the next.
        for i in np.arange(0,len(numbers)):
            prof[l, i] = float(numbers[i])
        l=l+1
