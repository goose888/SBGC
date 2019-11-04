# -*- Python source file -*-
"""
File Name : Get_BD.py

Description : Get the bulk density profile from He et al., 2016 data and write into netcdf file for ISAM simulations.

Created on : Mon Jun 18 23:23:25 2018

Last Modified : Mon Jun 18 23:25:00 2018

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
import time

cutdep = 100.
filename = 'Non_peat_data_permafrost.csv'
data = pd.read_csv(filename,encoding='iso-8859-1',index_col='ProfileID', skiprows=[1])
all_profid = data.index.unique()
lons = prep.getvarxls(data,'Lon',all_profid[0:16],0)
lats = prep.getvarxls(data,'Lat',all_profid[0:16],0)

## Shall use the mass perserving slpine (Bishop et al., 1999, Geoderma) to interpolate the bulk density data
## Okay, seems like only R has the offical package that can perform this interpolation
## so we shall use the R code to get this done instead.
#data_out = data[['Veg_Class', 'Lat', 'Lon', 'Basal_depth', 'Layer_depth', 'Soil_Order', \
#                                'Suborder', 'Horizon_type', 'bulk_density', 'C_Density']]
#data_out.to_csv('./sel_sites_bd.csv')
## Then call R script to run the mass-preseving spline interpolation
#status = au.exec_shell('Rscript', 'SOCinterp.r', [])

# Read in the interpolated SOC profile
sel_profid = data.index.unique()
filename = 'BDprofile.csv'
bd_interped = pd.read_csv(filename,encoding='iso-8859-1',index_col=0)
pid = bd_interped.index
# bd_interped = bd_interped * 1000.   # gcm-3 to kgm-3
bd_interped.index = sel_profid

# Average the interpolated SOC profile onto ISAM model soil layer depth
z, dz, zsoih = isam.get_isam_soildp(10)
bd_isam = isam.mean_by_depth(10, zsoih, sel_profid.size, bd_interped.as_matrix())
bd_isam[bd_isam<=0]=float('nan')

# Store BD into CSV file
bd_isam.to_csv('./sel_sites_bd.csv')

# Or store BD into a new NC file
nc = Dataset('BD_sites.nc', 'w', format ='NETCDF4_CLASSIC')
level = nc.createDimension('level', 10) 
site = nc.createDimension('site', sel_profid.size)
# Coordinate variables
level = nc.createVariable('level', np.float64, ('level',)) 
site = nc.createVariable('site', np.int32, ('site',)) 
# Actual variable
siteid = nc.createVariable('siteid', np.int32, ('site')) 
BD = nc.createVariable('BD', np.float64, ('site','level')) 
# Global Attributes 
nc.description = 'Bulk density of the sampled gelisol sites'
nc.history = 'Created ' + time.ctime(time.time())  
nc.source = 'Gelisol sites from He et al., 2016'
# Variable Attributes  
level.units = 'm'
site.units = 'N/A'
siteid.units = 'N/A'
BD.units = 'g/cm3'

# Assign values
level[:] = z
site[:] = np.arange(1,sel_profid.size+1)
siteid[:] = np.unique(data.index.values)
BD[:,:] = bd_isam

# Finilize
nc.close()

# np.savetxt('./isam_all_bd.csv', bd_isam, delimiter=',')   # X is an array('./isam_all_bd.csv')
