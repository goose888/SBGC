# -*- coding: utf-8 -*-
"""
Module for plotting box chart

Created on Thu Feb 12 20:28:02 2018

@author: Shijie Shu
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

#======================================
## Plot simulated D14C results
#======================================
filename = 'Non_peat_data_permafrost.csv'
data = pd.read_csv(filename,encoding='iso-8859-1',index_col='ProfileID', skiprows=[1])
all_profid = data.index.unique()
lons = prep.getvarxls(data,'Lon',all_profid[0:16],0)
lats = prep.getvarxls(data,'Lat',all_profid[0:16],0)

# Read in model output
d14cm = pd.read_table('isam_dc14.dat', header=None, delimiter=r"\s+")
# Read in another case from the ISAM model output
# d14cm2 = pd.read_table('isam_dc14.dat', header=None, delimiter=r"\s+")
d14cm.columns = ['ID', 'Layer1', 'Layer2', 'Layer3', 'Layer4', 'Layer5', 'Layer6', 'Layer7', 'Layer8', 'Layer9', 'Layer10']
d14cm = d14cm.set_index('ID')
mod_profid = d14cm.index
# d14cm2.columns = ['ID', 'Layer1', 'Layer2', 'Layer3', 'Layer4', 'Layer5', 'Layer6', 'Layer7', 'Layer8', 'Layer9', 'Layer10']
# d14cm2 = d14cm2.set_index('ID')
# mod2_profid = d14cm2.index
z, dz, zsoih = isam.get_isam_soildp(10)

# Extract corresponding measurements for each site and make figure
data.nodedepth = data.Layer_top + (data.Layer_bottom - data.Layer_top)/2.
myarray = d14cm.index.unique()
for i in myarray:
    d14co = data.D14C_BulkLayer.loc[i]
    nd = data.nodedepth.loc[i]
    if(d14co.__class__.__name__ == 'float64'):
        Xobs = d14co   # SOC profile kgCm-3
        Yobs = nd      # 1cm to 200cm
    else:
        Xobs = d14co.as_matrix()   # SOC profile kgCm-3
        Yobs = nd.as_matrix()     # 1cm to 200cm
    Xmod = d14cm.loc[i].as_matrix()   # SOC profile kgCm-3
    Ymod = z*100.     # 1cm to 200cm
    # Xmod2 = d14cm2.loc[i].as_matrix()   # SOC profile kgCm-3
    # Ymod2 = z*100.     # 1cm to 200cm
    if (data.Site[i].__class__.__name__ == 'Series'):
        tit = data.Site[i].unique().astype('string')[0]
    else:
        tit = data.Site[i].encode('ascii','ignore')
    path = './Figs_obsvsmod_calibrate/'+str(i)+'_'+tit+'.png'
    xticks = (-1000, -800, -600, -400, -200, 0, 200)
    #status = socplt.plot_obsvsmod(Xobs, Yobs, Xmod, Ymod, tit, path, xticks)
    status = socplt.plot_obsvsmod(Xobs, Yobs, Xmod, Ymod, tit, path, None, None, xticks)

#========================================================
## Plot the results of sensitivity analysis for D14C
#========================================================
# User options:
siteid = 146

filename = 'Non_peat_data_permafrost.csv'
data = pd.read_csv(filename,encoding='iso-8859-1',index_col='ProfileID', skiprows=[1])
all_profid = data.index.unique()
lons = prep.getvarxls(data,'Lon',all_profid[0:16],0)
lats = prep.getvarxls(data,'Lat',all_profid[0:16],0)

# Read in model output from the sensitivity tests cases
fopen = 'isam_'+str(siteid)+'_sensitivity.dat'
d14cm = pd.read_table(fopen, header=None, delimiter=r"\s+")
# d14cm2 = pd.read_table('isam_143_soc.dat', header=None, delimiter=r"\s+")
d14cm.columns = ['ID', 'Layer1', 'Layer2', 'Layer3', 'Layer4', 'Layer5', 'Layer6', 'Layer7', 'Layer8', 'Layer9', 'Layer10']
d14cm = d14cm.set_index('ID')
mod_profid = d14cm.index
# d14cm2.columns = ['ID', 'Layer1', 'Layer2', 'Layer3', 'Layer4', 'Layer5', 'Layer6', 'Layer7', 'Layer8', 'Layer9', 'Layer10']
# d14cm2 = d14cm2.set_index('ID')
# mod2_profid = d14cm2.index
z, dz, zsoih = isam.get_isam_soildp(10)

# Extract corresponding measurements for each site and make figure
data.nodedepth = data.Layer_top + (data.Layer_bottom - data.Layer_top)/2.
myarray = d14cm.index.unique()

# Check BD perturbation
d14co = data.D14C_BulkLayer.loc[siteid]
nd = data.nodedepth.loc[siteid]
if(d14co.__class__.__name__ == 'float64'):
    Xobs = d14co   # SOC profile kgCm-3
    Yobs = nd      # 1cm to 200cm
else:
    Xobs = d14co.as_matrix()   # SOC profile kgCm-3
    Yobs = nd.as_matrix()     # 1cm to 200cm
Xmod = d14cm.iloc[0].as_matrix()   # SOC profile kgCm-3
Ymod = z*100.     # 1cm to 200cm
Xmod2 = d14cm.iloc[1].as_matrix()   # SOC profile kgCm-3
Ymod2 = z*100.     # 1cm to 200cm
if (data.Site[siteid].__class__.__name__ == 'Series'):
    tit = data.Site[siteid].unique().astype('string')[0]
else:
    tit = data.Site[siteid].encode('ascii','ignore')
path = './Figs_obsvsmod_calibrate_sensitivity/'+str(siteid)+'_'+tit+'_bd.png'
xticks = (-1000, -800, -600, -400, -200, 0, 200)
# status = socplt.plot_sensitivity(Xobs, Yobs, Xmod, Ymod, tit, path, Xmod2, Ymod2, xticks)
status = socplt.plot_sensitivity(Xobs, Yobs, Xmod, Ymod, tit, path, Xmod2, Ymod2, xticks)

# Check Diffusivity perturbation
d14co = data.D14C_BulkLayer.loc[siteid]
nd = data.nodedepth.loc[siteid]
if(d14co.__class__.__name__ == 'float64'):
    Xobs = d14co   # SOC profile kgCm-3
    Yobs = nd      # 1cm to 200cm
else:
    Xobs = d14co.as_matrix()   # SOC profile kgCm-3
    Yobs = nd.as_matrix()     # 1cm to 200cm
Xmod = d14cm.iloc[3].as_matrix()   # SOC profile kgCm-3
Ymod = z*100.     # 1cm to 200cm
Xmod2 = d14cm.iloc[4].as_matrix()   # SOC profile kgCm-3
Ymod2 = z*100.     # 1cm to 200cm
if (data.Site[siteid].__class__.__name__ == 'Series'):
    tit = data.Site[siteid].unique().astype('string')[0]
else:
    tit = data.Site[siteid].encode('ascii','ignore')
path = './Figs_obsvsmod_calibrate_sensitivity/'+str(siteid)+'_'+tit+'_d.png'
xticks = (-1000, -800, -600, -400, -200, 0, 200)
#status = socplt.plot_obsvsmod(Xobs, Yobs, Xmod, Ymod, tit, path, xticks)
status = socplt.plot_sensitivity(Xobs, Yobs, Xmod, Ymod, tit, path, Xmod2, Ymod2, xticks)

# Check Cryoturbation depth perturbation
d14co = data.D14C_BulkLayer.loc[siteid]
nd = data.nodedepth.loc[siteid]
if(d14co.__class__.__name__ == 'float64'):
    Xobs = d14co   # SOC profile kgCm-3
    Yobs = nd      # 1cm to 200cm
else:
    Xobs = d14co.as_matrix()   # SOC profile kgCm-3
    Yobs = nd.as_matrix()     # 1cm to 200cm
Xmod = d14cm.iloc[5].as_matrix()   # SOC profile kgCm-3
Ymod = z*100.     # 1cm to 200cm
Xmod2 = d14cm.iloc[6].as_matrix()   # SOC profile kgCm-3
Ymod2 = z*100.     # 1cm to 200cm
if (data.Site[siteid].__class__.__name__ == 'Series'):
    tit = data.Site[siteid].unique().astype('string')[0]
else:
    tit = data.Site[siteid].encode('ascii','ignore')
path = './Figs_obsvsmod_calibrate_sensitivity/'+str(siteid)+'_'+tit+'_cdepth.png'
xticks = (-1000, -800, -600, -400, -200, 0, 200)
#status = socplt.plot_obsvsmod(Xobs, Yobs, Xmod, Ymod, tit, path, xticks)
status = socplt.plot_sensitivity(Xobs, Yobs, Xmod, Ymod, tit, path, Xmod2, Ymod2, xticks)

# Check Q10 perturbation
d14co = data.D14C_BulkLayer.loc[siteid]
nd = data.nodedepth.loc[siteid]
if(d14co.__class__.__name__ == 'float64'):
    Xobs = d14co   # SOC profile kgCm-3
    Yobs = nd      # 1cm to 200cm
else:
    Xobs = d14co.as_matrix()   # SOC profile kgCm-3
    Yobs = nd.as_matrix()     # 1cm to 200cm
Xmod = d14cm.iloc[7].as_matrix()   # SOC profile kgCm-3
Ymod = z*100.     # 1cm to 200cm
Xmod2 = d14cm.iloc[8].as_matrix()   # SOC profile kgCm-3
Ymod2 = z*100.     # 1cm to 200cm
if (data.Site[siteid].__class__.__name__ == 'Series'):
    tit = data.Site[siteid].unique().astype('string')[0]
else:
    tit = data.Site[siteid].encode('ascii','ignore')
path = './Figs_obsvsmod_calibrate_sensitivity/'+str(siteid)+'_'+tit+'_q10.png'
xticks = (-1000, -800, -600, -400, -200, 0, 200)
#status = socplt.plot_obsvsmod(Xobs, Yobs, Xmod, Ymod, tit, path, xticks)
status = socplt.plot_sensitivity(Xobs, Yobs, Xmod, Ymod, tit, path, Xmod2, Ymod2, xticks)

# Check S perturbation
d14co = data.D14C_BulkLayer.loc[siteid]
nd = data.nodedepth.loc[siteid]
if(d14co.__class__.__name__ == 'float64'):
    Xobs = d14co   # SOC profile kgCm-3
    Yobs = nd      # 1cm to 200cm
else:
    Xobs = d14co.as_matrix()   # SOC profile kgCm-3
    Yobs = nd.as_matrix()     # 1cm to 200cm
Xmod = d14cm.iloc[11].as_matrix()   # SOC profile kgCm-3
Ymod = z*100.     # 1cm to 200cm
Xmod2 = d14cm.iloc[12].as_matrix()   # SOC profile kgCm-3
Ymod2 = z*100.     # 1cm to 200cm
if (data.Site[siteid].__class__.__name__ == 'Series'):
    tit = data.Site[siteid].unique().astype('string')[0]
else:
    tit = data.Site[siteid].encode('ascii','ignore')
path = './Figs_obsvsmod_calibrate_sensitivity/'+str(siteid)+'_'+tit+'_fz.png'
xticks = (-1000, -800, -600, -400, -200, 0, 200)
#status = socplt.plot_obsvsmod(Xobs, Yobs, Xmod, Ymod, tit, path, xticks)
status = socplt.plot_sensitivity(Xobs, Yobs, Xmod, Ymod, tit, path, Xmod2, Ymod2, xticks)

# Check RD perturbation
d14co = data.D14C_BulkLayer.loc[siteid]
nd = data.nodedepth.loc[siteid]
if(d14co.__class__.__name__ == 'float64'):
    Xobs = d14co   # SOC profile kgCm-3
    Yobs = nd      # 1cm to 200cm
else:
    Xobs = d14co.as_matrix()   # SOC profile kgCm-3
    Yobs = nd.as_matrix()     # 1cm to 200cm
Xmod = d14cm.iloc[9].as_matrix()   # SOC profile kgCm-3
Ymod = z*100.     # 1cm to 200cm
Xmod2 = d14cm.iloc[10].as_matrix()   # SOC profile kgCm-3
Ymod2 = z*100.     # 1cm to 200cm
if (data.Site[siteid].__class__.__name__ == 'Series'):
    tit = data.Site[siteid].unique().astype('string')[0]
else:
    tit = data.Site[siteid].encode('ascii','ignore')
path = './Figs_obsvsmod_calibrate_sensitivity/'+str(siteid)+'_'+tit+'_rd.png'
xticks = (-1000, -800, -600, -400, -200, 0, 200)
#status = socplt.plot_obsvsmod(Xobs, Yobs, Xmod, Ymod, tit, path, xticks)
status = socplt.plot_sensitivity(Xobs, Yobs, Xmod, Ymod, tit, path, Xmod2, Ymod2, xticks)

# Check Tao perturbation
d14co = data.D14C_BulkLayer.loc[siteid]
nd = data.nodedepth.loc[siteid]
if(d14co.__class__.__name__ == 'float64'):
    Xobs = d14co   # SOC profile kgCm-3
    Yobs = nd      # 1cm to 200cm
else:
    Xobs = d14co.as_matrix()   # SOC profile kgCm-3
    Yobs = nd.as_matrix()     # 1cm to 200cm
Xmod = d14cm.iloc[13].as_matrix()   # SOC profile kgCm-3
Ymod = z*100.     # 1cm to 200cm
Xmod2 = d14cm.iloc[14].as_matrix()   # SOC profile kgCm-3
Ymod2 = z*100.     # 1cm to 200cm
if (data.Site[siteid].__class__.__name__ == 'Series'):
    tit = data.Site[siteid].unique().astype('string')[0]
else:
    tit = data.Site[siteid].encode('ascii','ignore')
path = './Figs_obsvsmod_calibrate_sensitivity/'+str(siteid)+'_'+tit+'_tao.png'
xticks = (-1000, -800, -600, -400, -200, 0, 200)
#status = socplt.plot_obsvsmod(Xobs, Yobs, Xmod, Ymod, tit, path, xticks)
status = socplt.plot_sensitivity(Xobs, Yobs, Xmod, Ymod, tit, path, Xmod2, Ymod2, xticks)

#========================================================
## Calculate teh refined Willmott's index
## Here we must use the interpolated D14C profile
#========================================================
## Read in D14C csv file
filename = 'Non_peat_data_permafrost.csv'
data = pd.read_csv(filename,encoding='iso-8859-1',index_col='ProfileID', skiprows=[1])
all_profid = data.index.unique()
lons = prep.getvarxls(data,'Lon',all_profid[0:16],0)
lats = prep.getvarxls(data,'Lat',all_profid[0:16],0)

z, dz, zsoih = isam.get_isam_soildp(10)
z = z*100
zsoih = zsoih*100
# First test extracting site 143

# Read in model output
d14cm = pd.read_table('isam_dc14.dat', header=None, delimiter=r"\s+")
# Read in another case from the ISAM model output
# d14cm2 = pd.read_table('isam_dc14.dat', header=None, delimiter=r"\s+")
d14cm.columns = ['ID', 'Layer1', 'Layer2', 'Layer3', 'Layer4', 'Layer5', 'Layer6', 'Layer7', 'Layer8', 'Layer9', 'Layer10']
d14cm = d14cm.set_index('ID')
mod_profid = d14cm.index
# d14cm2.columns = ['ID', 'Layer1', 'Layer2', 'Layer3', 'Layer4', 'Layer5', 'Layer6', 'Layer7', 'Layer8', 'Layer9', 'Layer10']
# d14cm2 = d14cm2.set_index('ID')
# mod2_profid = d14cm2.index

# Extract corresponding measurements for each site and calculate the refined Willmott's index
data.nodedepth = data.Layer_top + (data.Layer_bottom - data.Layer_top)/2.
myarray = d14cm.index.unique()
# The refined willmott's index is poor for sample 197
# myarray=[197]
for i in myarray:
    # interpolate d14C for missing values if any
    data.nodedepth = data.Layer_top + (data.Layer_bottom - data.Layer_top)/2.
    d14C = data.D14C_BulkLayer.loc[i]
    nd = data.nodedepth.loc[i]

    notNANs = ~np.isnan(d14C)
    f_i = interp1d(nd[notNANs], d14C[notNANs])
    f_x = prep.extrap1d(f_i)

    # Get the linear interpolated D14C profile at the same depth as model soil layers
    d14Cinterp = f_x(z)

    # Only truncate the interpolated values!
    for j in range(0,10):
        if(z[j]>nd.values[-1]):
            d14Cinterp[j]=float("nan") 

    d14co = data.D14C_BulkLayer.loc[i]
    nd = data.nodedepth.loc[i]
    Xobs = d14Cinterp   # D14C profile
    Yobs = nd
    Xmod = d14cm.loc[i].as_matrix()   # D14C profile   
    Ymod = z
    absdiff_obs = np.nansum(np.abs(Xobs - np.nanmean(Xobs)))
    absdiff_mod = np.nansum(np.abs(Xmod - Xobs))
    if(absdiff_mod <= 2*absdiff_obs):
        wmidx = 1 - absdiff_mod/(2*absdiff_obs)
        #print('tag1')
    else:
        wmidx = 2*absdiff_obs/absdiff_mod - 1
        #print('tag2')
    print(wmidx)

#========================================================
## Linear interpolation of all D14C profile
## These interpolated D14C are only used for 
## model calibration
#========================================================
## Read in D14C csv file
filename = 'Non_peat_data_permafrost.csv'
data = pd.read_csv(filename,encoding='iso-8859-1',index_col='ProfileID', skiprows=[1])
all_profid = data.index.unique()
lons = prep.getvarxls(data,'Lon',all_profid[0:16],0)
lats = prep.getvarxls(data,'Lat',all_profid[0:16],0)

# interpolate d14C for missing values if any
z, dz, zsoih = isam.get_isam_soildp(10)
z = z*100
zsoih = zsoih*100
# First test extracting site 143

myarray = [43, 44, 45, 47, 110, 143, 144, 145, 146, 154, 197, 198, 199]
#myarray=[143]
for i in myarray:
    data.nodedepth = data.Layer_top + (data.Layer_bottom - data.Layer_top)/2.
    d14C = data.D14C_BulkLayer.loc[i]
    nd = data.nodedepth.loc[i]

    notNANs = ~np.isnan(d14C)
    f_i = interp1d(nd[notNANs], d14C[notNANs])
    f_x = prep.extrap1d(f_i)

    # Get the linear interpolated D14C profile at the same depth as model soil layers
    d14Cinterp = f_x(z)

    # Only truncate the interpolated values!
    for j in range(0,10):
        if(z[j]>nd.values[-1]):
            d14Cinterp[j]=float("nan") 

    # Write out as observations for the usage of FFSQP
    fout = open("delc14_obs." + str(i), "w")
    for j in range(0,10):
        if(zsoih[j+1]<nd.values[-1]):
            fout.write("%d  %.8f\n" % \
                      (j+1, d14Cinterp[j]))

