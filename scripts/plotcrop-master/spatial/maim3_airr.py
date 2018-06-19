from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp,maskoceans
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import glob
import numpy.ma as ma
from scipy.interpolate import griddata
import matplotlib.colors as colors


iizumi=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/iizumi/iizumi.2013JAN29.maize.1982-2006.30min.nc4','r')
#print iizumi
iyield = iizumi.variables['yield50'][18,:,:]
iarea =iizumi.variables['area'][18,:,:]
la=iizumi.variables['lat'][:]
lo=iizumi.variables['lon'][:]
iyield=N.flipud(iyield)
iarea=N.flipud(iarea)

region1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/HistoricalGLM_crop_150901.nc','r')
maitrop = region1.variables['maize_trop'][99,:,:]
maitemp = region1.variables['maize_temp'][99,:,:]
maitropi=region1.variables['maize_trop_irrig'][99,:,:]
maitempi=region1.variables['maize_temp_irrig'][99,:,:]
maitrop=ma.masked_where(maitrop<=0,maitrop)
maitrop=ma.filled(maitrop, fill_value=0.)
maitemp=ma.masked_where(maitemp<=0,maitemp)
maitemp=ma.filled(maitemp, fill_value=0.)
maitropi=ma.masked_where(maitropi<=0,maitropi)
maitropi=ma.filled(maitropi, fill_value=0.)
maitempi=ma.masked_where(maitempi<=0,maitempi)
maitempi=ma.filled(maitempi, fill_value=0.)
maizetor=maitrop+maitemp
maizetoi=maitropi+maitempi
maizeto = maitrop+maitemp+maitropi+maitempi


clm=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/maizetrop_historical_co2_rf_fert_0.5x0.5.nc','r')
clmtrop = clm.variables['yield'][96:103,:,:]
clmtropf=N.average(clmtrop,axis=0)

clm1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/maizetemp_historical_co2_rf_fert_0.5x0.5.nc','r')
clmtemp = clm1.variables['yield'][96:103,:,:]
clmtempf=N.average(clmtemp,axis=0)

clm2=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/maizetrop_historical_co2_irrig_fert_0.5x0.5.nc','r')
clmtropfi1 = clm2.variables['yield'][96:103,:,:]
clmtropfi=N.average(clmtropfi1,axis=0)

clm3=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/maizetemp_historical_co2_irrig_fert_0.5x0.5.nc','r')
clmtempfi2 = clm3.variables['yield'][96:103,:,:]
clmtempfi=N.average(clmtempfi2,axis=0)

clmtropf=N.flipud(clmtropf)
clmtempf=N.flipud(clmtempf)
clmtropfi=N.flipud(clmtropfi)
clmtempfi=N.flipud(clmtempfi)

clmtropf= ma.masked_where(maitrop<=0,clmtropf)
clmtempf= ma.masked_where(maitemp<=0,clmtempf)
clmtropf=ma.filled(clmtropf, fill_value=0.)
clmtempf=ma.filled(clmtempf, fill_value=0.)

clmtropfi= ma.masked_where(maitropi<=0,clmtropfi)
clmtempfi= ma.masked_where(maitempi<=0,clmtempfi)
clmtropfi=ma.filled(clmtropfi, fill_value=0.)
clmtempfi=ma.filled(clmtempfi, fill_value=0.)

yield_clmtf=clmtropf+clmtempf
yield_clmtf = ma.masked_where(yield_clmtf<=0,yield_clmtf)
#yield_clmtf  = ma.masked_where(maizetor<=0,yield_clmtf )
yield_clmtf=ma.filled(yield_clmtf, fill_value=0.)


yield_clmtfi=clmtropfi+clmtempfi
yield_clmtfi = ma.masked_where(yield_clmtfi<=0,yield_clmtfi)
#yield_clmtfi = ma.masked_where(maizetoi<=0,yield_clmtfi)
yield_clmtfi=ma.filled(yield_clmtfi, fill_value=0.)



area=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/maizegridarea.nc','r')
gridarea = area.variables['cell_area'][:,:]

nclu=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/maize_AreaYieldProduction.nc','r')
ncvar_maize = nclu.variables['maizeData'][:]
latnc = nclu.variables['latitude'][:]
znc = nclu.variables['level'][:]
lonnc = nclu.variables['longitude'][:]
timenc = nclu.variables['time'][:]
#lon,lat = N.meshgrid(lonnc,latnc)
lon,lat = N.meshgrid(lonnc,latnc)

ncvar_maizef= N.zeros((2160, 4320))
ncvar_maizef=ncvar_maize[0,1,:,:]
ncvar_maize1=ncvar_maize[0,4,:,:]
ncvar_mask= N.zeros((2160, 4320))
ncvar_mask=ncvar_maize[0,0,:,:]


ncvar_maizef[N.isnan(ncvar_maizef)] = -9999
ncvar_mask[N.isnan(ncvar_mask)] = -9999
ncvar_maize1[N.isnan(ncvar_maize1)] = -9999

ncvar_maizef = ma.masked_where(ncvar_maizef<=0,ncvar_maizef)
#ncvar_maizef= ma.masked_where(ncvar_mask<0.01,ncvar_maizef)
ncvar_maizef = ma.masked_where(ncvar_maize1<=0,ncvar_maizef)
#ncvar_maizef = ma.masked_where(mask_clm<=0,ncvar_maizef)


yieldf= N.zeros((7, 360, 720))
yieldf2= N.zeros((7, 360, 720))
years = range(1997, 2004)
for i, year in enumerate(years):
    base = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_his/cesmhis/output/cesmhis.bgp-yearly_crop_{0}.nc".format(year), mode='r')
    lona1 = base.variables["lon"][:]
    lata1 = base.variables["lat"][:]
    basei = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_his/cesmhis_irr/output/cesmhisirr.bgp-yearly_crop_{0}.nc".format(year), mode='r')
    yield1 = base.variables["yield"][0,:,:]
    yieldf[i, :, :] = yield1
    yield2 = basei.variables["yield"][0,:,:]
    yieldf2[i, :, :] = yield2

yielda=N.average(yieldf,axis=0)
yielda2=N.average(yieldf2,axis=0)

yield_new,lona11 = shiftgrid(180.5,yielda,lona1,start=False)
yield_new2,lona11 = shiftgrid(180.5,yielda2,lona1,start=False)

lon2,lat2 = N.meshgrid(lona11,lata1)
#print lon2.shape
#ncvar_maize2 = interp(ncvar_maizea,lonnc,lat_new,lon,lat,order=1)

#yield_smooth = interp(maize_new, lonnc,lat_new, lons_sub,lats_sub , order=1)
yield_fine = interp(yield_new, lona11,lata1,lon,lat  , order=1)
yield_fine2 = interp(yield_new2, lona11,lata1,lon,lat  , order=1)

yield_clmtf1 = interp(yield_clmtf, lona11,lata1,lon,lat  , order=1)
yield_clmtfi1 = interp(yield_clmtfi, lona11,lata1,lon,lat  , order=1)

mask_clm = interp(maizeto, lona11,lata1,lon,lat  , order=1)
mask_clmra=interp(maizetor, lona11,lata1,lon,lat  , order=1)
mask_clmii=interp(maizetoi, lona11,lata1,lon,lat  , order=1)



fig, ax = plt.subplots(figsize=(8,6))


ax.set_title("M3 Maize Yield (t/ha)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()



lone,late = N.meshgrid(lonnc,latnc) #Returns coordinate matrices from coordinate vectors
x,y = map(lone,late)
ncvar_maizef=maskoceans(x,y,ncvar_maizef)
ncvar_maizef=ma.filled(ncvar_maizef, fill_value=0.)
ncvar_maizef = ma.masked_where(ncvar_maizef<=0,ncvar_maizef)
ncvar_maizef = ma.masked_where(mask_clm<=0,ncvar_maizef)

#print maskus_new
m3y=ncvar_maizef
cs = map.pcolormesh(x,y,m3y,cmap=plt.cm.YlGn,vmin=0,vmax=16)
cbar = map.colorbar(cs,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15) 
plt.axis('off')
plt.savefig('maiirrm3.jpg',dpi=300,bbox_inches='tight')


fig = plt.figure(figsize=(20,15))

ax2 = fig.add_subplot(321)
ax2.set_title("ISAM-CESM Maize Yield (t/ha)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(projection='robin',lon_0=0,resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
yield_fine=maskoceans(x,y,yield_fine)
yield_fine=ma.filled(yield_fine, fill_value=0.)
yield_fine = ma.masked_where(yield_fine<=0,yield_fine)
yield_fine = ma.masked_where(ncvar_maizef<=0,yield_fine)
yield_fine2=maskoceans(x,y,yield_fine2)
yield_fine2=ma.filled(yield_fine2, fill_value=0.)
yield_fine2 = ma.masked_where(yield_fine2<=0,yield_fine2)
yield_fine2 = ma.masked_where(ncvar_maizef<=0,yield_fine2)

isamy=((yield_fine*mask_clmra*gridarea)+(yield_fine2*mask_clmii*gridarea))/((mask_clmra*gridarea)+(mask_clmii*gridarea))


cs1 = map.pcolormesh(x,y,isamy,cmap=plt.cm.YlGn,vmin=0,vmax=16)
#cs1 = map.pcolormesh(x,y,yield_fine,cmap=plt.cm.jet,vmin=0.0,vmax=15)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15) 
plt.axis('off')


ax2 = fig.add_subplot(322)
ax2.set_title("CLM Maize Yield (t/ha)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(projection='robin',lon_0=0,resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()

yield_clmtf1=maskoceans(x,y,yield_clmtf1)
yield_clmtf1 = ma.masked_where(mask_clm<=0,yield_clmtf1)
yield_clmtfi1=maskoceans(x,y,yield_clmtfi1)
yield_clmtfi1 = ma.masked_where(mask_clm<=0,yield_clmtfi1)
clmy=((yield_clmtf1*mask_clmra*gridarea)+(yield_clmtfi1*mask_clmii*gridarea))/((mask_clmra*gridarea)+(mask_clmii*gridarea))
clmy = ma.masked_where(ncvar_maizef<=0,clmy)


cs1 = map.pcolormesh(x,y,clmy,cmap=plt.cm.YlGn,vmin=0,vmax=16)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15) 
plt.axis('off')



ax2 = fig.add_subplot(323)
ax2.set_title("ISAM-M3 Maize Yield (t/ha)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,isamy-m3y,cmap=plt.cm.bwr,vmin=-5,vmax=5)
#cs1 = map.pcolormesh(x,y,(yield_fine)-(ncvar_maizef),cmap=plt.cm.bwr,vmin=-0.5,vmax=0.5)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15) 
plt.axis('off')


ax2 = fig.add_subplot(324)
ax2.set_title("CLM-M3 Maize Yield (t/ha)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(projection='robin',lon_0=0,resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,clmy-m3y,cmap=plt.cm.bwr,vmin=-5,vmax=5)
#cs1 = map.pcolormesh(x,y,(yield_fine)-(ncvar_maizef),cmap=plt.cm.bwr,vmin=-0.5,vmax=0.5)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15) 
plt.axis('off')


ax2 = fig.add_subplot(325)
ax2.set_title("ISAM-M3 Maize Yield gridcell (%)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(projection='robin',lon_0=0,resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,(isamy-m3y)/m3y*100,cmap=plt.cm.bwr,vmin=-100,vmax=100)
#cs1 = map.pcolormesh(x,y,(yield_fine)-(ncvar_maizef),cmap=plt.cm.bwr,vmin=-0.5,vmax=0.5)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15)
plt.axis('off')


ax2 = fig.add_subplot(326)
ax2.set_title("CLM-M3 Maize Yield gridcell (%)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(projection='robin',lon_0=0,resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,(clmy-m3y)/m3y*100,cmap=plt.cm.bwr,vmin=-100,vmax=100)
#cs1 = map.pcolormesh(x,y,(yield_fine)-(ncvar_maizef),cmap=plt.cm.bwr,vmin=-0.5,vmax=0.5)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15)
plt.axis('off')



plt.savefig('maim3irr1_cesm.jpg',dpi=300,bbox_inches='tight')
#plt.show()



fig = plt.figure(figsize=(8,12))
colors = (0,0,1)
colorsr = (1,0,0)
ax = fig.add_subplot(211)
#ax.plot([0,700],[0,700], 'k--',label='1:1')
#ax.scatter(m3y, isamy, c=colors,alpha=1,label='ISAM')
ax.scatter(m3y, clmy, c=colorsr,alpha=0.5,label='CLM')
plt.xlim(0, 16)
plt.ylim(0, 16)
ax.plot([0,16],[0,16], 'k--',label='1:1')

#ax.set_title('Maize yield over gridcell',fontsize=18)
#ax.legend()
plt.tick_params(axis='both',labelsize=15)
plt.xlabel('M3-Crops (t/ha)',fontsize=18)
plt.ylabel('Model (t/ha)',fontsize=18)

ax = fig.add_subplot(212)
#ax.plot([0,700],[0,700], 'k--',label='1:1')
ax.scatter(m3y, isamy, c=colors,alpha=0.5,label='ISAM')
#ax.scatter(m3y, clmy, c=colorsr,alpha=1,label='CLM')
plt.xlim(0, 16)
plt.ylim(0, 16)
#ax.set_title('Soybean yield over gridcell',fontsize=18)
ax.plot([0,16],[0,16], 'k--',label='1:1')

#ax.legend()
#ax.plot([0,1000],[0,1000], 'k--',label='1:1')
plt.xticks([])
plt.tick_params(axis='both',labelsize=15)
#plt.xlabel('M3-Crops (g/$\mathregular{m^2}$)',fontsize=18)
plt.ylabel('Model (t/ha)',fontsize=18)


plt.savefig('scatter_maim3irr1_cesm.png',bbox_inches='tight')
#plt.show()



