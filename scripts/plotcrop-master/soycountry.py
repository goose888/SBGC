from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import glob
import numpy.ma as ma
from scipy.interpolate import griddata

iizumi=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/iizumi/iizumi.2013JAN29.soybean.1982-2006.30min.nc4','r')
#print iizumi
iyield = iizumi.variables['yield50'][18,:,:]
iarea =iizumi.variables['area'][18,:,:]
la=iizumi.variables['lat'][:]
lo=iizumi.variables['lon'][:]
iyield=N.flipud(iyield)
iarea=N.flipud(iarea)
#print iyield
country=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/Ctry_halfdeg.nc','r')
#print iizumi
coun = country.variables['MASK_Country'][:,:]

country3=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/Ctry_M3.nc','r')
#print iizumi
counm3 = country3.variables['MASK_Country'][:,:]

region1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/HistoricalGLM_crop_150901.nc','r')
maitrop = region1.variables['soy_trop'][99,:,:]
maitemp = region1.variables['soy_temp'][99,:,:]
maizeto = maitrop+maitemp


clm=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/soytrop_historical_co2_rf_nofert_0.5x0.5.nc','r')
clmtropf = clm.variables['yield'][99,:,:]
#clmtropf=N.average(clmtrop,axis=0)

clm1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/soytemp_historical_co2_rf_nofert_0.5x0.5.nc','r')
clmtempf = clm1.variables['yield'][99,:,:]
#clmtempf=N.average(clmtemp,axis=0)

#print clmtropf
#print clmtropf.shape
clmtropf=N.flipud(clmtropf)
clmtempf=N.flipud(clmtempf)


clmtropf= ma.masked_where(maitrop<=0,clmtropf)
clmtempf= ma.masked_where(maitemp<=0,clmtempf)
clmtropf=ma.filled(clmtropf, fill_value=0.)
clmtempf=ma.filled(clmtempf, fill_value=0.)




area=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/gridareahalf.nc','r')
gridarea = area.variables['cell_area'][:,:]
gridlon = area.variables['lon'][:]

gridarea,gridlon = shiftgrid(180.5,gridarea,gridlon,start=False)

area1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/maizegridarea.nc','r')
gridaream3 = area1.variables['cell_area'][:,:]
gridlonm3 = area1.variables['lon'][:]
gridlatm3 = area1.variables['lat'][:]

#print gridarea

nclu=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/soybean_AreaYieldProduction.nc','r')
#print(nclu)
ncvar_maize = nclu.variables['soybeanData'][:]
latnc = nclu.variables['latitude'][:]
znc = nclu.variables['level'][:]
lonnc = nclu.variables['longitude'][:]
timenc = nclu.variables['time'][:]
#lon,lat = N.meshgrid(lonnc,latnc)
lat_new=N.flipud(latnc)
lon,lat = N.meshgrid(lonnc,lat_new)
ncvar_maizea = nclu.variables['soybeanData'][0,1,:,:]
ncvar_maizearea = nclu.variables['soybeanData'][0,4,:,:]

ncvar_maizea=N.flipud(ncvar_maizea)

#print lat_new
maize_new = N.flipud(ncvar_maize)
ncvar_maizearea = N.flipud(ncvar_maizearea)



yieldf= N.zeros((1, 360, 720))
yieldf2= N.zeros((1, 360, 720))
years = range(2000, 2001)
for i, year in enumerate(years):
    base = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_his/cruhis/output/cruhis.bgp-yearly_crop_{0}.nc".format(year), mode='r')
    base2 = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_his/cesmhis/output/cesmhis.bgp-yearly_crop_{0}.nc".format(year), mode='r')
    lona1 = base.variables["lon"][:]
    lata1 = base.variables["lat"][:]
   
    yield1 = base.variables["yield"][1,:,:]
    yieldf[i, :, :] = yield1
    yield2 = base2.variables["yield"][1,:,:]
    yieldf2[i, :, :] = yield2

yielda=N.average(yieldf,axis=0)
yielda2=N.average(yieldf2,axis=0)

yield_new,lona11 = shiftgrid(180.5,yielda,lona1,start=False)
yield_new2,lona11 = shiftgrid(180.5,yielda2,lona1,start=False)

#print lona11
lon2,lat2 = N.meshgrid(lona11,lata1)
#print lon2.shape

#yield_smooth = interp(maize_new, lonnc,lat_new, lons_sub,lats_sub , order=1)
m3coun = interp(coun,lona11,lata1,lon,lat,order=1)

#http://matplotlib.org/basemap/users/mapsetup.html
iarea = ma.masked_where(iarea<=0,iarea)
iarea=ma.filled(iarea, fill_value=0.)

iyield = ma.masked_where(iyield<=0,iyield)
iyield = ma.masked_where(iarea<=0,iyield)
iyield=ma.filled(iyield, fill_value=0.)


yield_new[N.isnan(yield_new)] = -9999
yield_new = ma.masked_where(yield_new<=0,yield_new)
yield_new = ma.masked_where(iyield<=0,yield_new)
yield_new=ma.filled(yield_new, fill_value=0.)

yield_new2[N.isnan(yield_new2)] = -9999
yield_new2 = ma.masked_where(yield_new2<=0,yield_new2)
yield_new2 = ma.masked_where(iyield<=0,yield_new2)
yield_new2=ma.filled(yield_new2, fill_value=0.)


yield_clmtf=clmtropf+clmtempf
yield_clmtf = ma.masked_where(yield_clmtf<=0,yield_clmtf)
yield_clmtf = ma.masked_where(iyield<=0,yield_clmtf)
yield_clmtf=ma.filled(yield_clmtf, fill_value=0.)


ncvar_maizea[N.isnan(ncvar_maizea)] = -9999
ncvar_maizearea[N.isnan(ncvar_maizearea)] = -9999
ncvar_maizea = ma.masked_where(ncvar_maizea<=0,ncvar_maizea)
ncvar_maizearea = ma.masked_where(ncvar_maizearea<=0,ncvar_maizearea)
ncvar_maizearea=ma.filled(ncvar_maizearea, fill_value=0.)
ncvar_maizea=ma.filled(ncvar_maizea, fill_value=0.)

def annualyield(couna,counb):

	yieldagf=0.
	yieldg=0.
	harea=0.
	a=0
	yieldagfi=0.
	yieldgi=0.

	yieldagfc=0.
	yieldgc=0.
	yieldgi1=0.
	yieldagfi1=0.

#USA 11501~11550
	for xx in range(0,360):
	    for yy in range(0,720):
	        if coun[xx,yy] >=couna and  coun[xx,yy] <=counb:
	            yieldg=iyield[xx,yy]*iarea[xx,yy]/100*gridarea[xx,yy]/10000+yieldg
	            harea=iarea[xx,yy]/100*gridarea[xx,yy]/10000+ harea
	            yieldgi=yield_new[xx,yy]*iarea[xx,yy]/100*gridarea[xx,yy]/10000+yieldgi
	            yieldgc=yield_clmtf[xx,yy]*iarea[xx,yy]/100*gridarea[xx,yy]/10000+yieldgc
                    yieldgi1=yield_new2[xx,yy]*iarea[xx,yy]/100*gridarea[xx,yy]/10000+yieldgi1
        	    a=a+1

	yieldagf=yieldg/harea
	yieldagfi=yieldgi/harea
	yieldagfc=yieldgc/harea
	yieldagfi1=yieldgi1/harea

	a=0
	m3harea=0.
	m3yieldg=0.
	m3yieldag=0.
	for xx in range(0,2160):
	    for yy in range(0,4320):
	        if m3coun[xx,yy] >=couna and  m3coun[xx,yy] <=counb:
	            m3yieldg=ncvar_maizea[xx,yy]*ncvar_maizearea[xx,yy]+m3yieldg
	            m3harea=ncvar_maizearea[xx,yy]+ m3harea
                     
	            a=a+1

	m3yieldag=m3yieldg/m3harea
	return "Soybean","iizumi harvested area",harea, "ha","iizumi yield",yieldagf,"t/ha","production",yieldg,"tonnes","ISAM-CRU yield",yieldagfi,"t/ha","production",yieldgi,"tonnes","CLM yield",yieldagfc,"t/ha","production",yieldgc,"tonnes","ISAM-CESM yield",yieldagfi1,"t/ha","production",yieldgi1,"tonnes","M3 harvested area",m3harea, "ha","M3 yield",m3yieldag,"t/ha","production",m3yieldg,"tonnes"

name=["globe","us","brazil","argentina","china","india","paraguay","canada","uruguay","italy","bolivia"]
range1=[10100,11501,20301,20101,41501,42901,20900,10201,21300,43400,20200]
range2=[50700,11550,20327,20124,41529,42929,20900,10212,21300,43400,20200]
a=0
x=11
toarea= N.zeros(x)
zumiy= N.zeros(x)
clmy= N.zeros(x)
isamy= N.zeros(x)
zumip= N.zeros(x)
clmp= N.zeros(x)
isamp= N.zeros(x)
isamy1= N.zeros(x)
isamp1= N.zeros(x)
m3y=N.zeros(x)
m3p=N.zeros(x)
for a, name1 in enumerate(name):

	reu=annualyield(range1[a],range2[a])
        toarea[a]=reu[2]
        zumiy[a]=reu[5]
        isamy[a]=reu[11]
        clmy[a]=reu[17]
        zumip[a]=reu[8]
        isamp[a]=reu[14]
        clmp[a]=reu[20]
        isamy1[a]=reu[23]
        isamp1[a]=reu[26]
        m3y[a]=reu[32]
        m3p[a]=reu[35]
        

# data to plot
n_groups = 11
#means_frank = (yieldagf, yieldagf1, yieldagf2, yieldagf3,yieldagf4)
#means_guido = (yieldagfi, yieldagfi1, yieldagfi2, yieldagfi3,yieldagfi4)
#means_final = (yieldagfc, yieldagfc1, yieldagfc2, yieldagfc3,yieldagfc4)
#means_fao = (85910*0.0001, 45974*0.0001, 27447*0.0001, 54329*0.0001,24620*0.0001)
#means_m3 = (m3yieldag, m3yieldag1, m3yieldag2, m3yieldag3,m3yieldag4)

means_frank = zumiy
means_guido = isamy
means_final = clmy
means_fao = (21689*0.0001,25613*0.0001, 23999*0.0001,23312*0.0001, 16559*0.0001,8222*0.0001,25331*0.0001,25483*0.0001,7640*0.0001,35761*0.0001,19406*0.0001)
means_m3 = m3y
means_cesm=isamy1
#means_frank1 = (yieldg, yieldg1, yieldg2, yieldg3,yieldg4)
#means_guido1 = (yieldgi, yieldgi1, yieldgi2, yieldgi3,yieldgi4)
#means_final1 = (yieldgc, yieldgc1, yieldgc2, yieldgc3,yieldgc4)
#pro_fao = (251852210, 106000000, 31879392, 16780650,17556900)
#pro_m3 = (m3yieldg, m3yieldg1, m3yieldg2, m3yieldg3,m3yieldg4)


means_frank1 = zumip
means_guido1 = isamp
means_final1 = clmp
pro_fao = (161297357,75053800, 32734958,20135800 ,15411200 ,5275800,2980060,2703000,6800,903490,1197251)
pro_m3 = m3p
pro_cesm=isamp1

# create plot
#fig, ax = plt.subplots(211)
fig = plt.figure(figsize=(18,18))


ax = fig.add_subplot(211)
index = N.arange(n_groups)
bar_width = 0.15
opacity = 0.8
rects0 = plt.bar(index, means_fao, bar_width,
                 alpha=opacity,
                 color='k',
                 label='FAO')
rects4 = plt.bar(index+bar_width, means_m3, bar_width,
                 alpha=opacity,
                 color='c',
                 label='M3')
rects1 = plt.bar(index+bar_width*2, means_frank, bar_width,
                 alpha=opacity,
                 color='r',
                 label='Iizumi')
 
rects2 = plt.bar(index + bar_width*3, means_guido, bar_width,
                 alpha=opacity,
                 color='g',
                 label='ISAM-NCEP')
rects0 = plt.bar(index + bar_width*4, means_cesm, bar_width,
                 alpha=opacity,
                 color='y',
                 label='ISAM-CESM')
rects3 = plt.bar(index + bar_width*5, means_final, bar_width,
                 alpha=opacity,
                 color='b',
                 label='CLM')
 
plt.xlabel('Country',fontsize=18)
plt.ylabel('Soybean yield (t/ha)',fontsize=18)
plt.title('2000',fontsize=18)
plt.xticks(index + bar_width+0.2, ('Global','USA','Brazil','Argentina','China','India','Paraguay','Canada','Uruguay','Italy','Bolivia'))
leg=plt.legend()
leg.get_frame().set_alpha(0.5)
plt.tick_params(axis='both',labelsize=15)

plt.tight_layout()


ax = fig.add_subplot(212)
index = N.arange(n_groups)
bar_width = 0.15
opacity = 0.8
rects0 = plt.bar(index, pro_fao, bar_width,
                 alpha=opacity,
                 color='k',
                 label='FAO')
rects4 = plt.bar(index+bar_width, pro_m3, bar_width,
                 alpha=opacity,
                 color='c',
                 label='M3')
rects1 = plt.bar(index+ bar_width*2, means_frank1, bar_width,
                 alpha=opacity,
                 color='r',
                 label='Iizumi')
 
rects2 = plt.bar(index + bar_width*3, means_guido1, bar_width,
                 alpha=opacity,
                 color='g',
                 label='ISAM-NCEP')
rects0 = plt.bar(index + bar_width*4, pro_cesm, bar_width,
                 alpha=opacity,
                 color='y',
                 label='ISAM-CESM')

rects3 = plt.bar(index + bar_width*5, means_final1, bar_width,
                 alpha=opacity,
                 color='b',
                 label='CLM')
 
plt.xlabel('Country',fontsize=18)
plt.ylabel('Soybean production (tonnes)',fontsize=18)
plt.title('2000',fontsize=18)
plt.xticks(index + bar_width+0.2, ('Global','USA','Brazil','Argentina','China','India','Paraguay','Canada','Uruguay','Italy','Bolivia'))
leg=plt.legend()
leg.get_frame().set_alpha(0.5)
plt.tick_params(axis='both',labelsize=15)

plt.tight_layout()
plt.savefig('soy2000edison.png')
#plt.show()


