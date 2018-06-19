from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp
from netCDF4 import Dataset as NetCDFFile
import numpy as N
from calendar import isleap
import matplotlib.pyplot as plt
import glob
import numpy.ma as ma
import scipy.stats
from matplotlib.dates import DateFormatter
import datetime
import pandas as pd
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
import csv




data=open('wtd_hr_leap.csv','r')
sh=csv.reader(data)
#print sh
tw=[]
ch4=[]
fch4=[]
next(sh, None)  # skip the headers

for row in sh:
    tw1=row[1]
    ch41=row[2]
    fch41=row[3]
    tw.append(tw1)
    ch4.append(ch41)
    fch4.append(fch41)
#print tw
tw = [float(x) for x in tw]
ch4 = [float(x) for x in ch4]
fch4 = [float(x) for x in fch4]
#tw=[ x if x != -9999 else 'nan' for x in tw ]
#print tw
ch4=N.asarray(ch4)
ch4=ma.masked_where(ch4==-9999,ch4)
fch4=N.asarray(fch4)
fch4=ma.masked_where(fch4==-9999,fch4)
tw=N.asarray(tw)
for x in range (0,87648):
	if tw[x]==-9999.0:
		#print tw[x],tw[x-1]
		tw[x]=tw[x-1]
		#print 'after',tw[x]
print tw
tw=ma.masked_where(tw==-9999,tw)
tw=ma.filled(tw, fill_value=-0.001)
tw1=-tw
tw1=ma.masked_where(tw1<=0,tw1)
tw1=ma.filled(tw1, fill_value=0.001)
print tw1


#writing nc file
ncfile=NetCDFFile('US-Twt_riceleap_flood.nc','w',format='NETCDF3_64BIT_OFFSET')

# dimensions
ncfile.createDimension('lat', 1)
ncfile.createDimension('lon', 1)
ncfile.createDimension('time',87648)
#
latitudes = ncfile.createVariable('lat', 'f8', ('lat',))
longitudes = ncfile.createVariable('lon', 'f8', ('lon',))
times=ncfile.createVariable('time', 'f8', ('time',))
#
maize4= ncfile.createVariable('FCH4', 'f8', ('time','lat','lon'),fill_value=-9999.)
maize5= ncfile.createVariable('CH4', 'f8', ('time','lat','lon'),fill_value=-9999.)
maize6= ncfile.createVariable('WT', 'f8', ('time','lat','lon'),fill_value=-9999.)
latitudes[:] = 38.0498
longitudes[:] = -121.7651
maize4[:]=fch4
maize5[:]=ch4
maize6[:]=0.001
#
#
latitudes.units = 'degrees_north'
longitudes.units = 'degrees_east'
maize4.units='nmol CH4 m-2 s-1'
maize5.units='nmol CH4 mol-1'
maize6.units='m'

maize4.long_name = 'Methane fluxes'
maize5.long_name = 'Soil methane concentration'
maize6.long_name = 'Water table depth(+below soil)'
#
times.long_name = 'local time -8 half-hour start at 2010-01-01 00'
latitudes.long_name = 'latitude'
longitudes.long_name = 'longitude'
#

ncfile.close()

