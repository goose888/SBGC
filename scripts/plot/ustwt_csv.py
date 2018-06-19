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

#data strating from 2010 to 2015 local time every half-hourly -8 wetland from shijie
data=NetCDFFile('/data/jain1/b/team/datasets4/SITE_WT/US-Twt_WT.nc','r')
tws=data.variables['WT'][:,0,0]
data1=NetCDFFile('US-Twt_riceleap.nc','r')
twsrice=data1.variables['WT'][:,0,0]




#2010-2014 excluding 2/29/2012 us-twt rice from fluxnet
#prescribed wetland water table depth
isam=open('/data/jain1/d/tslin2/ISAM/samcrop/edison/rice/rice_isam/methane/site_daily_10_tgas_wetland.txt')
#prescribed rice water table depth
isam1=open('/data/jain1/d/tslin2/ISAM/samcrop/edison/rice/rice_isam/methane/site_daily_10_tgas.txt')
#dynamic water table depth
isam2=open('/data/jain1/d/tslin2/ISAM/samcrop/edison/rice/rice_isam/methane/fort.204')

isam2hy=open('/data/jain1/d/tslin2/ISAM/samcrop/edison/rice/rice_isam/methane/site_daily_5_hydro.txt')

ch4a = []
ch4a1= []
ch4a2= []
wta2=[]
for line in isam.readlines():
    y = [value for value in line.split()]
    ch4a.append( y[0] )
ch4a = [float(x) for x in ch4a]


for line in isam1.readlines():
    y = [value for value in line.split()]
    ch4a1.append( y[0] )
#print ch4a1
ch4a1 = [float(x) for x in ch4a1]

for line in isam2.readlines():
    y = [value for value in line.split()]
    ch4a2.append( y[15] )
ch4a2 = [float(x) for x in ch4a2]
a=1
for line in isam2hy.readlines():
	y = [value for value in line.split()]
	a=a+1
	if a%2==0:
    		wta2.append( y[9] )
wta2 = [float(x) for x in wta2]

data=open('wtd_hr.csv','r')
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
tw=N.asarray(tw)
tw=ma.masked_where(tw==-9999,tw)
ch4=N.asarray(ch4)
ch4=ma.masked_where(ch4==-9999,ch4)
fch4=N.asarray(fch4)
fch4=ma.masked_where(fch4==-9999,fch4)


a=1
all_twsrice=0
all_tws=0
all_tw=0
all_ch4=0
all_fch4=0
daily_twsrice=N.zeros(1825)
daily_tws=N.zeros(1825)
daily_tw=N.zeros(1825)
daily_ch4=N.zeros(1825)
daily_fch4=N.zeros(1825)

x=0
for num in range(0,87600):

	if a>=1 & a<=48:
                all_twsrice=twsrice[num]+all_twsrice
		all_tws=tws[num]+all_tws
		all_tw=tw[num]+all_tw
	        all_ch4=ch4[num]+all_ch4
	        all_fch4=fch4[num]+all_fch4

		a=a+1
	if a==49:
                daily_twsrice[x]=all_twsrice/(a-1)
		daily_tws[x]=all_tws/(a-1)
		daily_tw[x]=all_tw/(a-1)
		daily_ch4[x]=all_ch4/(a-1)
                daily_fch4[x]=all_fch4/(a-1)
		all_twsrice=0
		all_tws=0
		all_tw=0
		all_ch4=0
		all_fch4=0
		x=x+1
		a=1

#print daily_tw,x

b=4
ch4ay=ch4a[0+1825*b:1825+1825*b]
#print ch4ay
ch4ay=N.asarray(ch4ay)

ch4a1y=ch4a1[0+1825*b:1825+1825*b]
ch4a1y=N.asarray(ch4a1y)
ch4a2y=ch4a2[0+1825*b:1825+1825*b]
ch4a2y=N.asarray(ch4a2y)

wta2y=wta2[0+1825*b:1825+1825*b]
wta2y=N.asarray(wta2y)


fig = plt.figure(figsize=(20,10))
ax = fig.add_subplot(311)

#xdates = [datetime.datetime.strptime(str(int(date)),'%Y') for date in xx]
#xdates = [datetime.date(2010, 1, 1), datetime.date(2015, 12, 27)] 
months    = MonthLocator(range(1,13), bymonthday=1, interval=4)
monthsFmt = DateFormatter("%b '%y")

xdates = pd.date_range('2010-01-01', periods=1826, freq='D')
#print xdates
#xdates_new=pd.DatetimeIndex(data=(t for t in xdates if not isleap(t.year)), freq="D")
fig.autofmt_xdate()

leap = []
for each in xdates:
    if each.month==2 and each.day ==29:
        leap.append(each)

xdates = xdates.drop(leap)
#print xdates 
plt.ylim(-100,800)

ax.plot_date(xdates,daily_fch4*12*86400/1000000,"ko",label="OBS")
#ax.plot_date(xdates,ch4ay*1000,"r-",label="MOD-Wetland")
ax.plot_date(xdates,ch4a1y*1000,"r-",label="MOD")
#ax.plot_date(xdates,ch4a2y*1000,"g-",label="MOD-Flood")

x1, y1 = ['2010-04-16', '2010-04-16'], [-90, 790]
x2, y2 = ['2010-10-28','2010-10-28'],[-90, 790]
x3, y3 = ['2011-04-22', '2011-04-22'], [-90, 790]
x4, y4 = ['2011-10-13','2011-10-13'],[-90, 790]
x5, y5 = ['2012-05-17', '2012-05-17'], [-90, 790]
x6, y6 = ['2012-11-14','2012-11-14'],[-90, 790]
x7, y7 = ['2013-04-02', '2013-04-02'], [-90, 790]
x8, y8 = ['2013-09-20','2013-09-20'],[-90, 790]
x9, y9 = ['2014-04-28', '2014-04-28'], [-90, 790]
x10, y10 = ['2014-10-02','2014-10-02'],[-90, 790]


plt.plot(x1,y1,marker='>',color='b',linestyle='--',linewidth=1.5)
plt.plot(x2,y2,marker='<',color='b',linestyle=':',linewidth=1.5)
plt.plot(x3,y3,marker='>',color='b',linestyle='--',linewidth=1.5)
plt.plot(x4,y4,marker='<',color='b',linestyle=':',linewidth=1.5)
plt.plot(x5,y5,marker='>',color='b',linestyle='--',linewidth=1.5)
plt.plot(x6,y6,marker='<',color='b',linestyle=':',linewidth=1.5)
plt.plot(x7,y7,marker='>',color='b',linestyle='--',linewidth=1.5)
plt.plot(x8,y8,marker='<',color='b',linestyle=':',linewidth=1.5)
plt.plot(x9,y9,marker='>',color='b',linestyle='--',linewidth=1.5)
plt.plot(x10,y10,marker='<',color='b',linestyle=':',linewidth=1.5)



#fig.autofmt_xdate()
ax.xaxis.set_major_locator(months)
ax.xaxis.set_major_formatter(DateFormatter('%b %y'))

plt.ylabel("FCH4 (mgC/m2/d)",fontsize=18)
plt.title("US-Twt",fontsize=18)
plt.tick_params(axis='both',labelsize=18)

leg = plt.legend(loc=2,fancybox=True, fontsize=16)
leg.get_frame().set_alpha(0.5)



ax = fig.add_subplot(312)
plt.tick_params(axis='both',labelsize=18)
ax.tick_params('y', colors='k')

ax.plot_date(xdates,ch4a2y,"r-",label="MOD")

ax.xaxis.set_major_locator(months)
ax.xaxis.set_major_formatter(DateFormatter('%b %y'))

plt.ylabel("NPP (gC/$\mathregular{m^{2}}$/d)",fontsize=18)
plt.tick_params(axis='both',labelsize=18)
x1, y1 = ['2010-04-16', '2010-04-16'], [-1.6, 15.6]
x2, y2 = ['2010-10-28','2010-10-28'],[-1.6, 15.6]
x3, y3 = ['2011-04-22', '2011-04-22'], [-1.6, 15.6]
x4, y4 = ['2011-10-13','2011-10-13'],[-1.6, 15.6]
x5, y5 = ['2012-05-17', '2012-05-17'], [-1.6, 15.6]
x6, y6 = ['2012-11-14','2012-11-14'],[-1.6, 15.6]
x7, y7 = ['2013-04-02', '2013-04-02'], [-1.6, 15.6]
x8, y8 = ['2013-09-20','2013-09-20'],[-1.6, 15.6]
x9, y9 = ['2014-04-28', '2014-04-28'], [-1.6, 15.6]
x10, y10 = ['2014-10-02','2014-10-02'],[-1.6, 15.6]


plt.plot(x1,y1,marker='>',color='b',linestyle='--',linewidth=1.5)
plt.plot(x2,y2,marker='<',color='b',linestyle=':',linewidth=1.5)
plt.plot(x3,y3,marker='>',color='b',linestyle='--',linewidth=1.5)
plt.plot(x4,y4,marker='<',color='b',linestyle=':',linewidth=1.5)
plt.plot(x5,y5,marker='>',color='b',linestyle='--',linewidth=1.5)
plt.plot(x6,y6,marker='<',color='b',linestyle=':',linewidth=1.5)
plt.plot(x7,y7,marker='>',color='b',linestyle='--',linewidth=1.5)
plt.plot(x8,y8,marker='<',color='b',linestyle=':',linewidth=1.5)
plt.plot(x9,y9,marker='>',color='b',linestyle='--',linewidth=1.5)
plt.plot(x10,y10,marker='<',color='b',linestyle=':',linewidth=1.5)


ax = fig.add_subplot(313)
plt.ylim(-100,40)
plt.tick_params(axis='both',labelsize=18)
ax.tick_params('y', colors='k')

ax.plot_date(xdates,daily_tw*100,"ko",label="OBS")
ax.plot_date(xdates,-daily_twsrice*100,"r-",label="MOD-Rice")
ax.plot_date(xdates,-wta2y*100,"g-",label="MOD-Flood")

#ax2=ax.twinx()
#ax2.tick_params('y', colors='r')
#ax2.plot_date(xdates,-daily_tws*100,"r-",label="MOD-Wetland")
#plt.ylim(-1,1)
#plt.tick_params(axis='both',labelsize=18)


x1, y1 = ['2010-04-16', '2010-04-16'], [-98, 38]
x2, y2 = ['2010-10-28','2010-10-28'],[-98, 38]
x3, y3 = ['2011-04-22', '2011-04-22'], [-98, 38]
x4, y4 = ['2011-10-13','2011-10-13'],[-98, 38]
x5, y5 = ['2012-05-17', '2012-05-17'], [-98, 38]
x6, y6 = ['2012-11-14','2012-11-14'],[-98, 38]
x7, y7 = ['2013-04-02', '2013-04-02'], [-98, 38]
x8, y8 = ['2013-09-20','2013-09-20'],[-98, 38]
x9, y9 = ['2014-04-28', '2014-04-28'], [-98, 38]
x10, y10 = ['2014-10-02','2014-10-02'],[-98, 38]


plt.plot(x1,y1,marker='>',color='b',linestyle='--',linewidth=1.5)
plt.plot(x2,y2,marker='<',color='b',linestyle=':',linewidth=1.5)
plt.plot(x3,y3,marker='>',color='b',linestyle='--',linewidth=1.5)
plt.plot(x4,y4,marker='<',color='b',linestyle=':',linewidth=1.5)
plt.plot(x5,y5,marker='>',color='b',linestyle='--',linewidth=1.5)
plt.plot(x6,y6,marker='<',color='b',linestyle=':',linewidth=1.5)
plt.plot(x7,y7,marker='>',color='b',linestyle='--',linewidth=1.5)
plt.plot(x8,y8,marker='<',color='b',linestyle=':',linewidth=1.5)
plt.plot(x9,y9,marker='>',color='b',linestyle='--',linewidth=1.5)
plt.plot(x10,y10,marker='<',color='b',linestyle=':',linewidth=1.5)


fig.autofmt_xdate()
ax.xaxis.set_major_locator(months)
ax.xaxis.set_major_formatter(DateFormatter('%b %y'))

ax.set_ylabel("Water table depth (cm)",fontsize=18)

#leg = plt.legend(loc=2,fancybox=True, fontsize=16)
#leg.get_frame().set_alpha(0.5)


#ax1 = fig.add_subplot(212)
#
#plt.ylim(-200,800)
#
#xdates1 = pd.date_range('2010-01-01', periods=87600, freq='30min')
##print xdates1.shape
#ax1.plot_date(xdates1,fch4,"ro",label="OBS")
#fig.autofmt_xdate()
#ax1.xaxis.set_major_locator(months)
#ax1.xaxis.set_major_formatter(DateFormatter('%b %y'))
#
#plt.ylabel("FCH4 (nmol/m2/s)",fontsize=18)
#plt.tick_params(axis='both',labelsize=18)
#
#leg = plt.legend(loc=2,fancybox=True, fontsize=16)
#leg.get_frame().set_alpha(0.5)

plt.savefig('ustwt_fixed.jpg',bbox_inches='tight')

plt.show()


#
##writing nc file
#ncfile=NetCDFFile('US-Twt_rice.nc','w',format='NETCDF3_64BIT_OFFSET')
#
## dimensions
#ncfile.createDimension('lat', 1)
#ncfile.createDimension('lon', 1)
#ncfile.createDimension('time',87600)
#
#latitudes = ncfile.createVariable('lat', 'f8', ('lat',))
#longitudes = ncfile.createVariable('lon', 'f8', ('lon',))
#times=ncfile.createVariable('time', 'f8', ('time',))
#
#maize4= ncfile.createVariable('FCH4', 'f8', ('time','lat','lon'),fill_value=-9999.)
#maize5= ncfile.createVariable('CH4', 'f8', ('time','lat','lon'),fill_value=-9999.)
#maize6= ncfile.createVariable('WT', 'f8', ('time','lat','lon'),fill_value=-9999.)
#latitudes[:] = 38.0498
#longitudes[:] = -121.7651
#maize4[:]=fch4
#maize5[:]=ch4
#maize6[:]=-tw
#
#
#latitudes.units = 'degrees_north'
#longitudes.units = 'degrees_east'
#maize4.units='nmol CH4 m-2 s-1'
#maize5.units='nmol CH4 mol-1'
#maize6.units='m'

#maize4.long_name = 'Methane fluxes'
#maize5.long_name = 'Soil methane concentration'
#maize6.long_name = 'Water table depth(+below soil)'
#
#times.long_name = 'local time -8 half-hour start at 2010-01-01 00'
#latitudes.long_name = 'latitude'
#longitudes.long_name = 'longitude'
#

#ncfile.close()

