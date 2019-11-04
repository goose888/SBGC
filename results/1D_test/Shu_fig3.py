# -*- Python source file -*-
"""
File Name : Shu_fig3.py

Description : Code for creating the Figure 3 in Shu et al., 2018 SBGC paper

Created on : Mon Jun 18 23:59:02 2018

Last Modified : Tue Jun 19 00:02:08 2018

Author : Shijie
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

#================================================================================
## Plot the results of ccalibration + sensitivity analysis for D14C in 1 figure
#================================================================================
# obs
obs_fname = 'Non_peat_data_permafrost.csv'
data = pd.read_csv(obs_fname,encoding='iso-8859-1',index_col='ProfileID', skiprows=[1])
all_profid = data.index.unique()
lons = prep.getvarxls(data,'Lon',all_profid[0:16],0)
lats = prep.getvarxls(data,'Lat',all_profid[0:16],0)

# Calibration
isam_cal_d14c_fname = 'isam_dc14_cali.dat'
isam_cal_soc_fname = 'isam_soc_cali.dat'

# Sensitivity test
#isam_sen43_d14c_fname = 'isam_dc14_43sen.dat'
#isam_sen43_soc_fname = 'isam_soc_43sen.dat'
#isam_sen110_d14c_fname = 'isam_dc14_110sen.dat'
#isam_sen110_soc_fname = 'isam_soc_110sen.dat'
#isam_sen143_d14c_fname = 'isam_dc14_143sen.dat'
#isam_sen143_soc_fname = 'isam_soc_143sen.dat'
#isam_sen146_d14c_fname = 'isam_dc14_146sen.dat'
#isam_sen146_soc_fname = 'isam_soc_146sen.dat'
#isam_sen197_d14c_fname = 'isam_dc14_197sen.dat'
#isam_sen197_soc_fname = 'isam_soc_197sen.dat'

# Um's results
isam_um_d14c_fname = 'isam_um_dc14.dat'
isam_um_soc_fname = 'isam_um_soc.dat'

# Create site dict
site = {
    "43" : 0,
    "110" : 1,
    "143" : 2,
    "197" : 3,
    "146" : 4,
    
}
#site = {
#    "197" : 4
#}

## DO NOT EDIT BELOW
## D14C
# Read in the calibrated model output from the optimized site cases
d14cm = pd.read_table(isam_cal_d14c_fname, header=None, delimiter=r"\s+")
d14cm.columns = ['ID', 'Layer1', 'Layer2', 'Layer3', 'Layer4', 'Layer5', 'Layer6', 'Layer7', 'Layer8', 'Layer9', 'Layer10']
d14cm = d14cm.set_index('ID')
mod_profid = d14cm.index
z, dz, zsoih = isam.get_isam_soildp(10)

#initialize
d14cm_test = [None] * 5

# Read in model output from the sensitivity test cases
for siteid, ind in site.items():
    fname = 'isam_dc14_'+str(siteid)+'sen.dat'
    d14cm_test[ind] = pd.read_table(fname, header=None, delimiter=r"\s+")
    # d14cm2 = pd.read_table('isam_143_soc.dat', header=None, delimiter=r"\s+")
    d14cm_test[ind].columns = ['ID', 'Layer1', 'Layer2', 'Layer3', 'Layer4', 'Layer5', 'Layer6', 'Layer7', 'Layer8', 'Layer9', 'Layer10']
    d14cm_test[ind] = d14cm_test[ind].set_index('ID')
    # mod_profid = d14cm.index

# Figure with the added shaded area representing the uncertainty for D14C
data.nodedepth = data.Layer_top + (data.Layer_bottom - data.Layer_top)/2.
myarray = [43, 110, 143, 146, 197]
for siteid, ind in site.items():
    d14co = data.D14C_BulkLayer.loc[int(float(siteid))]
    nd = data.nodedepth.loc[int(float(siteid))]
    if(d14co.__class__.__name__ == 'float64'):
        Xobs = d14co   # SOC profile kgCm-3
        Yobs = nd      # 1cm to 200cm
    else:
        Xobs = d14co.as_matrix()   # SOC profile kgCm-3
        Yobs = nd.as_matrix()     # 1cm to 200cm
    Xmod = d14cm.loc[int(float(siteid))].as_matrix()   # SOC profile kgCm-3
    Ymod = z*100.     # 1cm to 200cm
    Xmod[Xmod<-999.] = float("nan")
    Xobs[np.isnan(Xmod)] = float("nan")
    # Sensitivity tests
    # Mask the impact from bd
    d14cm_test[ind].iloc[0] = float("nan")
    d14cm_test[ind].iloc[1] = float("nan")
    # Mask the impact from max_cryo
    #d14cm_test[ind].iloc[5] = float("nan")
    #d14cm_test[ind].iloc[6] = float("nan")
    # Mask the impact from rd
    d14cm_test[ind].iloc[9] = float("nan")
    d14cm_test[ind].iloc[10] = float("nan")
    sen_lower = d14cm_test[ind].min().as_matrix()   # SOC profile kgCm-3
    sen_upper = d14cm_test[ind].max().as_matrix()
    Ysen = z*100.     # 1cm to 200cm
    if (data.Site[int(float(siteid))].__class__.__name__ == 'Series'):
        tit = data.Site[int(float(siteid))].unique().astype('string')[0]
    else:
        tit = data.Site[int(float(siteid))].encode('ascii','ignore')
    path = './Fig3_Shu/'+siteid+'_'+tit+'.png'
    xticks = (-1000, -800, -600, -400, -200, 0, 200)
    
    # status = socplt.plot_obsvsmod(Xobs, Yobs, Xmod, Ymod, tit, path, None, None, xticks)
    status = socplt.plot_obsvsmodshades(Xobs, Yobs, Xmod, Ymod, sen_lower, sen_upper, Ysen, tit, path, None, None, xticks)
    plt.savefig(path)
    #plt.close("all")

# SOC
# Original SOC profile
soc_obs = data.C_Density

# Read in model outputs
soccm = pd.read_table(isam_cal_soc_fname, header=None, delimiter=r"\s+")
soccm.columns = ['ID', 'Layer1', 'Layer2', 'Layer3', 'Layer4', 'Layer5', 'Layer6', 'Layer7', 'Layer8', 'Layer9', 'Layer10']
soccm = soccm.set_index('ID')
mod_profid = soccm.index

# Get the corresponding SOC comparable against ISAM
z, dz, zsoih = isam.get_isam_soildp(10)
zsoih = zsoih * 100

#initialize
soccm_test = [None] * 5

# Read in model output from the sensitivity test cases
for siteid, ind in site.items():
    fname = 'isam_soc_'+str(siteid)+'sen.dat'
    soccm_test[ind] = pd.read_table(fname, header=None, delimiter=r"\s+")
    # socm2 = pd.read_table('isam_143_soc.dat', header=None, delimiter=r"\s+")
    soccm_test[ind].columns = ['ID', 'Layer1', 'Layer2', 'Layer3', 'Layer4', 'Layer5', 'Layer6', 'Layer7', 'Layer8', 'Layer9', 'Layer10']
    soccm_test[ind] = soccm_test[ind].set_index('ID')
    # mod_profid = soccm.index
    # Mask the impact from bd
    soccm_test[ind].iloc[0] = float("nan")
    soccm_test[ind].iloc[1] = float("nan")
    # Mask the impact from max_cryo
    #soccm_test[ind].iloc[5] = float("nan")
    #soccm_test[ind].iloc[6] = float("nan")
    # Mask the impact from rd
    soccm_test[ind].iloc[9] = float("nan")
    soccm_test[ind].iloc[10] = float("nan")
    # soccm_test[3].iloc[11] = float("nan")

# Need to mask the soc profile result for Site 146 with high s (stronger depth modifier)
soccm_test[3].iloc[11] = float("nan")

# Extract corresponding measurements for each site and make figure
data.nodedepth = data.Layer_top + (data.Layer_bottom - data.Layer_top)/2.
# myarray = soccm.index.unique()
myarray = [43, 110, 143, 146, 197]
# myarray = [143] 

# Figure with the added shaded area representing the uncertainty for SOC
for siteid, ind in site.items():
    # Observation
    dn = data.nodedepth.loc[int(float(siteid))].as_matrix()
    dp = (data.Layer_bottom.loc[int(float(siteid))] - data.Layer_top.loc[int(float(siteid))]).as_matrix()
    num_obs = soc_obs[int(float(siteid))].__len__()
    accsoc = soc_obs[int(float(siteid))].as_matrix() * float('nan')
    sitesoc = soc_obs.loc[int(float(siteid))].as_matrix()    
    accsoc[0] = sitesoc[0] * 1000. * dp[0] / 100.
    accsoc[-1] = float('nan')
    for j in range(1,num_obs):
        accsoc[j] = accsoc[j-1] + sitesoc[j] * 1000. * dp[j] / 100.
    # Model results
    modsoc = soccm.loc[int(float(siteid))].as_matrix()
    acc_test = soccm_test[ind].as_matrix()
    for j in range(1,10):
        modsoc[j] = modsoc[j-1] + modsoc[j]
    # Sensitivity test results
        acc_test[:,j] = acc_test[:,j-1] + acc_test[:,j]
        
    # Make plot
    Xobs = accsoc   # SOC profile kgCm-3
    Yobs = dn       # 1cm to 200cm
    Xmod = modsoc   # SOC profile kgCm-3
    Ymod = z*100.     # 1cm to 200cm
    # Sensitivity tests
    sen_lower = np.nanmin(acc_test, axis=0)   # SOC profile kgCm-3
    sen_upper = np.nanmax(acc_test, axis=0)
    sen_upper[sen_upper>140] = float("nan")
    
    if (data.Site[int(float(siteid))].__class__.__name__ == 'Series'):
        tit = data.Site[int(float(siteid))].unique().astype('string')[0]
    else:
        tit = data.Site[int(float(siteid))].encode('ascii','ignore')
    path = './Figs_obsvsmod_original_soc/'+siteid+'_'+tit+'_soc.png'
    #xticks = (0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
    xticks = (0, 20, 40, 60, 80, 100, 120, 140)
    #status = socplt.plot_obsvsmod(Xobs, Yobs, Xmod, Ymod, tit, path, xticks)
    status = socplt.plot_obsvsmodshades(Xobs, Yobs, Xmod, Ymod, sen_lower, sen_upper, Ysen, tit, path, None, None, xticks)
    plt.savefig(path)
    #plt.show()
    #plt.close("all")

#=============================================================
# Plot the results of sensitivity test for each parameters
#=============================================================
## D14C
# Read in the calibrated model output from the optimized site cases
d14cm = pd.read_table(isam_cal_d14c_fname, header=None, delimiter=r"\s+")
d14cm.columns = ['ID', 'Layer1', 'Layer2', 'Layer3', 'Layer4', 'Layer5', 'Layer6', 'Layer7', 'Layer8', 'Layer9', 'Layer10']
d14cm = d14cm.set_index('ID')
mod_profid = d14cm.index
z, dz, zsoih = isam.get_isam_soildp(10)

#initialize
d14cm_test = [None] * 5

# Read in model output from the sensitivity test cases
for siteid, ind in site.items():
    fname = 'isam_dc14_'+str(siteid)+'sen.dat'
    d14cm_test[ind] = pd.read_table(fname, header=None, delimiter=r"\s+")
    # d14cm2 = pd.read_table('isam_143_soc.dat', header=None, delimiter=r"\s+")
    d14cm_test[ind].columns = ['ID', 'Layer1', 'Layer2', 'Layer3', 'Layer4', 'Layer5', 'Layer6', 'Layer7', 'Layer8', 'Layer9', 'Layer10']
    d14cm_test[ind] = d14cm_test[ind].set_index('ID')
    # mod_profid = d14cm.index

# Figure with the added shaded area representing the uncertainty for D14C
data.nodedepth = data.Layer_top + (data.Layer_bottom - data.Layer_top)/2.
myarray = [43, 110, 143, 146, 197]
for siteid, ind in site.items():
    d14co = data.D14C_BulkLayer.loc[int(float(siteid))]
    nd = data.nodedepth.loc[int(float(siteid))]
    if(d14co.__class__.__name__ == 'float64'):
        Xobs = d14co   # SOC profile kgCm-3
        Yobs = nd      # 1cm to 200cm
    else:
        Xobs = d14co.as_matrix()   # SOC profile kgCm-3
        Yobs = nd.as_matrix()     # 1cm to 200cm
    Xmod = d14cm.loc[int(float(siteid))].as_matrix()   # SOC profile kgCm-3
    Ymod = z*100.     # 1cm to 200cm
    Xmod[Xmod<-999.] = float("nan")
    Xobs[np.isnan(Xmod)] = float("nan")
    # Sensitivity tests
    # Mask the impact from bd
    #d14cm_test[ind].iloc[0] = float("nan")
    #d14cm_test[ind].iloc[1] = float("nan")
    # Mask the impact from max_cryo
    #d14cm_test[ind].iloc[5] = float("nan")
    #d14cm_test[ind].iloc[6] = float("nan")
    # Mask the impact from rd
    #d14cm_test[ind].iloc[9] = float("nan")
    #d14cm_test[ind].iloc[10] = float("nan")
    sen_lower = d14cm_test[ind].iloc[0].as_matrix()   # SOC profile kgCm-3
    sen_upper = d14cm_test[ind].iloc[1].as_matrix()
    Ysen = z*100.     # 1cm to 200cm
    if (data.Site[int(float(siteid))].__class__.__name__ == 'Series'):
        tit = data.Site[int(float(siteid))].unique().astype('string')[0]
    else:
        tit = data.Site[int(float(siteid))].encode('ascii','ignore')
    path = './Figs_obsvsmod_calibrate/'+siteid+'_'+tit+'bd.png'
    xticks = (-1000, -600, -200, 200)
    
    # status = socplt.plot_obsvsmod(Xobs, Yobs, Xmod, Ymod, tit, path, None, None, xticks)
    status = socplt.plot_obsvsmodshades_color(Xobs, Yobs, Xmod, Ymod, sen_lower, sen_upper, Ysen, tit, path, None, None, xticks)
    plt.savefig(path)
    
    sen_lower = d14cm_test[ind].iloc[3].as_matrix()   # SOC profile kgCm-3
    sen_upper = d14cm_test[ind].iloc[4].as_matrix()
    path = './Figs_obsvsmod_calibrate/'+siteid+'_'+tit+'d.png'
    status = socplt.plot_obsvsmodshades_color(Xobs, Yobs, Xmod, Ymod, sen_lower, sen_upper, Ysen, tit, path, None, None, xticks)
    plt.savefig(path)
    
    sen_lower = d14cm_test[ind].iloc[5].as_matrix()   # SOC profile kgCm-3
    sen_upper = d14cm_test[ind].iloc[6].as_matrix()
    path = './Figs_obsvsmod_calibrate/'+siteid+'_'+tit+'depth.png'
    status = socplt.plot_obsvsmodshades_color(Xobs, Yobs, Xmod, Ymod, sen_lower, sen_upper, Ysen, tit, path, None, None, xticks)
    plt.savefig(path)

    sen_lower = d14cm_test[ind].iloc[7].as_matrix()   # SOC profile kgCm-3
    sen_upper = d14cm_test[ind].iloc[8].as_matrix()
    path = './Figs_obsvsmod_calibrate/'+siteid+'_'+tit+'q.png'
    status = socplt.plot_obsvsmodshades_color(Xobs, Yobs, Xmod, Ymod, sen_lower, sen_upper, Ysen, tit, path, None, None, xticks)
    plt.savefig(path)
    
    sen_lower = d14cm_test[ind].iloc[9].as_matrix()   # SOC profile kgCm-3
    sen_upper = d14cm_test[ind].iloc[10].as_matrix()
    path = './Figs_obsvsmod_calibrate/'+siteid+'_'+tit+'rd.png'
    status = socplt.plot_obsvsmodshades_color(Xobs, Yobs, Xmod, Ymod, sen_lower, sen_upper, Ysen, tit, path, None, None, xticks)
    plt.savefig(path)
    
    sen_lower = d14cm_test[ind].iloc[11].as_matrix()   # SOC profile kgCm-3
    sen_upper = d14cm_test[ind].iloc[12].as_matrix()
    path = './Figs_obsvsmod_calibrate/'+siteid+'_'+tit+'s.png'
    status = socplt.plot_obsvsmodshades_color(Xobs, Yobs, Xmod, Ymod, sen_lower, sen_upper, Ysen, tit, path, None, None, xticks)
    plt.savefig(path)
    
    sen_lower = d14cm_test[ind].iloc[13].as_matrix()   # SOC profile kgCm-3
    sen_upper = d14cm_test[ind].iloc[14].as_matrix()
    path = './Figs_obsvsmod_calibrate/'+siteid+'_'+tit+'t.png'
    status = socplt.plot_obsvsmodshades_color(Xobs, Yobs, Xmod, Ymod, sen_lower, sen_upper, Ysen, tit, path, None, None, xticks)
    plt.savefig(path)

# SOC
# Original SOC profile
soc_obs = data.C_Density

# Read in model outputs
soccm = pd.read_table(isam_cal_soc_fname, header=None, delimiter=r"\s+")
soccm.columns = ['ID', 'Layer1', 'Layer2', 'Layer3', 'Layer4', 'Layer5', 'Layer6', 'Layer7', 'Layer8', 'Layer9', 'Layer10']
soccm = soccm.set_index('ID')
mod_profid = soccm.index

# Get the corresponding SOC comparable against ISAM
z, dz, zsoih = isam.get_isam_soildp(10)
zsoih = zsoih * 100

#initialize
soccm_test = [None] * 5

# Read in model output from the sensitivity test cases
for siteid, ind in site.items():
    fname = 'isam_soc_'+str(siteid)+'sen.dat'
    soccm_test[ind] = pd.read_table(fname, header=None, delimiter=r"\s+")
    # socm2 = pd.read_table('isam_143_soc.dat', header=None, delimiter=r"\s+")
    soccm_test[ind].columns = ['ID', 'Layer1', 'Layer2', 'Layer3', 'Layer4', 'Layer5', 'Layer6', 'Layer7', 'Layer8', 'Layer9', 'Layer10']
    soccm_test[ind] = soccm_test[ind].set_index('ID')
    # mod_profid = soccm.index
    # Mask the impact from bd
    #soccm_test[ind].iloc[0] = float("nan")
    #soccm_test[ind].iloc[1] = float("nan")
    # Mask the impact from max_cryo
    #d14cm_test[ind].iloc[5] = float("nan")
    #d14cm_test[ind].iloc[6] = float("nan")
    # Mask the impact from rd
    #soccm_test[ind].iloc[9] = float("nan")
    #soccm_test[ind].iloc[10] = float("nan")
    # soccm_test[3].iloc[11] = float("nan")

# Need to mask the soc profile result for Site 146 with high s (stronger depth modifier)
soccm_test[3].iloc[11] = float("nan")

# Extract corresponding measurements for each site and make figure
data.nodedepth = data.Layer_top + (data.Layer_bottom - data.Layer_top)/2.
# myarray = soccm.index.unique()
myarray = [43, 110, 143, 146, 197]
# myarray = [143] 

# Figure with the added shaded area representing the uncertainty for SOC
for siteid, ind in site.items():
    # Observation
    dn = data.nodedepth.loc[int(float(siteid))].as_matrix()
    dp = (data.Layer_bottom.loc[int(float(siteid))] - data.Layer_top.loc[int(float(siteid))]).as_matrix()
    num_obs = soc_obs[int(float(siteid))].__len__()
    accsoc = soc_obs[int(float(siteid))].as_matrix() * float('nan')
    sitesoc = soc_obs.loc[int(float(siteid))].as_matrix()    
    accsoc[0] = sitesoc[0] * 1000. * dp[0] / 100.
    accsoc[-1] = float('nan')
    for j in range(1,num_obs):
        accsoc[j] = accsoc[j-1] + sitesoc[j] * 1000. * dp[j] / 100.
    # Model results
    modsoc = soccm.loc[int(float(siteid))].as_matrix()
    acc_test = soccm_test[ind].as_matrix()
    for j in range(1,10):
        modsoc[j] = modsoc[j-1] + modsoc[j]
    # Sensitivity test results
        acc_test[:,j] = acc_test[:,j-1] + acc_test[:,j]
        
    # Make plot
    Xobs = accsoc   # SOC profile kgCm-3
    Yobs = dn       # 1cm to 200cm
    Xmod = modsoc   # SOC profile kgCm-3
    Ymod = z*100.     # 1cm to 200cm
    # Sensitivity tests
    sen_lower = acc_test[0,:]   # SOC profile kgCm-3
    sen_upper = acc_test[1,:]
    
    if (data.Site[int(float(siteid))].__class__.__name__ == 'Series'):
        tit = data.Site[int(float(siteid))].unique().astype('string')[0]
    else:
        tit = data.Site[int(float(siteid))].encode('ascii','ignore')
    path = './Figs_obsvsmod_original_soc/'+siteid+'_'+tit+'_soc_bd.png'
    #xticks = (0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
    xticks = (0, 20, 40, 60, 80, 100, 120, 140)
    #status = socplt.plot_obsvsmod(Xobs, Yobs, Xmod, Ymod, tit, path, xticks)
    status = socplt.plot_obsvsmodshades_color(Xobs, Yobs, Xmod, Ymod, sen_lower, sen_upper, Ysen, tit, path, None, None, xticks)
   
    sen_lower = acc_test[3,:]    # SOC profile kgCm-3
    sen_upper = acc_test[4,:] 
    path = './Figs_obsvsmod_calibrate/'+siteid+'_'+tit+'d.png'
    status = socplt.plot_obsvsmodshades_color(Xobs, Yobs, Xmod, Ymod, sen_lower, sen_upper, Ysen, tit, path, None, None, xticks)
    
    sen_lower = acc_test[5,:]    # SOC profile kgCm-3
    sen_upper = acc_test[6,:] 
    path = './Figs_obsvsmod_calibrate/'+siteid+'_'+tit+'depth.png'
    status = socplt.plot_obsvsmodshades_color(Xobs, Yobs, Xmod, Ymod, sen_lower, sen_upper, Ysen, tit, path, None, None, xticks)

    sen_lower = acc_test[7,:]    # SOC profile kgCm-3
    sen_upper = acc_test[8,:] 
    path = './Figs_obsvsmod_calibrate/'+siteid+'_'+tit+'q.png'
    status = socplt.plot_obsvsmodshades_color(Xobs, Yobs, Xmod, Ymod, sen_lower, sen_upper, Ysen, tit, path, None, None, xticks)
    
    sen_lower = acc_test[9,:]    # SOC profile kgCm-3
    sen_upper = acc_test[10,:] 
    path = './Figs_obsvsmod_calibrate/'+siteid+'_'+tit+'rd.png'
    status = socplt.plot_obsvsmodshades_color(Xobs, Yobs, Xmod, Ymod, sen_lower, sen_upper, Ysen, tit, path, None, None, xticks)    
    
    sen_lower = acc_test[11,:]    # SOC profile kgCm-3
    sen_upper = acc_test[12,:] 
    path = './Figs_obsvsmod_calibrate/'+siteid+'_'+tit+'s.png'
    status = socplt.plot_obsvsmodshades_color(Xobs, Yobs, Xmod, Ymod, sen_lower, sen_upper, Ysen, tit, path, None, None, xticks)    
    
    sen_lower = acc_test[13,:]    # SOC profile kgCm-3
    sen_upper = acc_test[14,:] 
    path = './Figs_obsvsmod_calibrate/'+siteid+'_'+tit+'t.png'
    status = socplt.plot_obsvsmodshades_color(Xobs, Yobs, Xmod, Ymod, sen_lower, sen_upper, Ysen, tit, path, None, None, xticks)
    
#=============================================================
# Make figure and add all sensitivity experiments together
#=============================================================
cases = {
    "bd_high" : 0,
    "bd_low" : 1,
    "d_high" : 3,
    "d_low" : 4,
    "dep_high" : 5,
    "dep_low" : 6,
    "q_high" : 7,
    "q_low" : 8,
    "rd_high" : 9,
    "rd_low" : 10,
    "s_high" : 11,
    "s_low" : 12,
    "tao_high" : 13,
    "tao_low" : 14
}

palette = {
    "chocolate" : 0,
    "darkorange" : 1,
    "tan" : 3,
    "oldlace" : 4,
    "gold" : 5,
    "khaki" : 6,
    "olive" : 7,
    "yellowgreen" : 8,
    "limegreen" : 9,
    "aquamarine" : 10,
    "teal" : 11,
    "royalblue" : 12,
    "blueviolet" : 13,
    "magenta" : 14
}

data.nodedepth = data.Layer_top + (data.Layer_bottom - data.Layer_top)/2.
myarray = [43, 110, 143, 146, 197]
for siteid, ind in site.items():
    d14co = data.D14C_BulkLayer.loc[int(float(siteid))]
    nd = data.nodedepth.loc[int(float(siteid))]
    if(d14co.__class__.__name__ == 'float64'):
        Xobs = d14co   # SOC profile kgCm-3
        Yobs = nd      # 1cm to 200cm
    else:
        Xobs = d14co.as_matrix()   # SOC profile kgCm-3
        Yobs = nd.as_matrix()     # 1cm to 200cm
    Xmod = d14cm.loc[int(float(siteid))].as_matrix()   # SOC profile kgCm-3
    Ymod = z*100.     # 1cm to 200cm
    # Sensitivity tests
    Xsen = d14cm_test[ind].as_matrix()   # SOC profile kgCm-3
    Ysen = z*100.     # 1cm to 200cm
    if (data.Site[int(float(siteid))].__class__.__name__ == 'Series'):
        tit = data.Site[int(float(siteid))].unique().astype('string')[0]
    else:
        tit = data.Site[int(float(siteid))].encode('ascii','ignore')
    path = './Figs_obsvsmod_calibrate/'+siteid+'_'+tit+'.png'
    xticks = (-1000, -800, -600, -400, -200, 0, 200)
    # status = socplt.plot_obsvsmod(Xobs, Yobs, Xmod, Ymod, tit, path, None, None, xticks)
    status = socplt.plot_obsvsmodplussensi(Xobs, Yobs, Xmod, Ymod, Xsen, Ysen, cases, palette, tit, path, None, None, xticks)
    

