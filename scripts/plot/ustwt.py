from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import glob
import numpy.ma as ma
import scipy.stats
from matplotlib.dates import DateFormatter
import datetime
import pandas as pd
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter

#data strating from 2010 to 2015 local time every half-hourly -8
data=NetCDFFile('/data/jain1/b/team/datasets4/SITE_WT/US-Twt_WT.nc','r')
tw=data.variables['WT'][:,0,0]
print tw.shape
ch4=data.variables['CH4'][:,0,0]
fch4=data.variables['FCH4'][:,0,0]
fig = plt.figure(figsize=(30,10))

a=1
all_tw=0
all_ch4=0
all_fch4=0
daily_tw=N.zeros(2191)
daily_ch4=N.zeros(2191)
daily_fch4=N.zeros(2191)

x=0
for num in range(0,105168):

	if a>=1 & a<=48:
		all_tw=tw[num]+all_tw
	        all_ch4=ch4[num]+all_ch4
	        all_fch4=fch4[num]+all_fch4

		a=a+1
	if a==49:
		daily_tw[x]=all_tw/(a-1)
		daily_ch4[x]=all_ch4/(a-1)
                daily_fch4[x]=all_fch4/(a-1)
		all_tw=0
		all_ch4=0
		all_fch4=0
		x=x+1
		a=1

print daily_tw,x
ax = fig.add_subplot(311)


#xdates = [datetime.datetime.strptime(str(int(date)),'%Y') for date in xx]
#xdates = [datetime.date(2010, 1, 1), datetime.date(2015, 12, 31)] 
months    = MonthLocator(range(1,13), bymonthday=1, interval=6)
monthsFmt = DateFormatter("%b '%y")


xdates = pd.date_range('2010-01-01', periods=2191, freq='D')
ax.plot_date(xdates,daily_tw,"r-o",label="OBS")
fig.autofmt_xdate()
ax.xaxis.set_major_locator(months)
ax.xaxis.set_major_formatter(DateFormatter('%b %y'))

plt.ylabel("Water table depth (m)",fontsize=18)
plt.title("US-Twt",fontsize=18)
plt.tick_params(axis='both',labelsize=18)

leg = plt.legend(loc=2,fancybox=True, fontsize=16)
leg.get_frame().set_alpha(0.5)


#ax = fig.add_subplot(312)




#xdates1 = pd.date_range('2010-01-01', periods=105168, freq='30min')
#print xdates1.shape
#ax.plot_date(xdates1,tw,"r-o",label="OBS")
#fig.autofmt_xdate()
##ax.xaxis.set_major_locator(months)
##ax.xaxis.set_major_formatter(DateFormatter('%b %y'))
#
#plt.ylabel("Water table depth (m)",fontsize=18)
#plt.title("US-Twt",fontsize=18)
#plt.tick_params(axis='both',labelsize=18)

#leg = plt.legend(loc=2,fancybox=True, fontsize=16)
#leg.get_frame().set_alpha(0.5)


plt.show()

