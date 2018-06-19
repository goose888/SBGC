from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp,maskoceans
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import glob
import numpy.ma as ma
from scipy.interpolate import griddata
import matplotlib.colors as colors


area=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/gridareahalf.nc','r')
gridarea = area.variables['cell_area'][:,:]
gridlon = area.variables['lon'][:]

gridarea,gridlon = shiftgrid(180.5,gridarea,gridlon,start=False)



iizumi=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/iizumi/iizumi.2013JAN29.maize.1982-2006.30min.nc4','r')
#print iizumi
iyield = iizumi.variables['yield50'][18,:,:]
iarea =iizumi.variables['area'][18,:,:]
la=iizumi.variables['lat'][:]
lo=iizumi.variables['lon'][:]
iyield=N.flipud(iyield)
iarea=N.flipud(iarea)
la=N.flipud(la)
lon2,lat2 = N.meshgrid(lo,la)
region1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/HistoricalGLM_crop_150901.nc','r')
maitrop = region1.variables['maize_trop'][99,:,:]
maitemp = region1.variables['maize_temp'][99,:,:]
maitropi=region1.variables['maize_trop_irrig'][99,:,:]
maitempi=region1.variables['maize_temp_irrig'][99,:,:]
gridarea = region1.variables['area'][:,:]
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
clmtropf = clm.variables['yield'][99,:,:]
#clmtropf=N.average(clmtrop,axis=0)

clm1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/maizetemp_historical_co2_rf_fert_0.5x0.5.nc','r')
clmtempf = clm1.variables['yield'][99,:,:]

clm2=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/maizetrop_historical_co2_irrig_fert_0.5x0.5.nc','r')
clmtropfi = clm2.variables['yield'][99,:,:]

clm3=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/maizetemp_historical_co2_irrig_fert_0.5x0.5.nc','r')
clmtempfi = clm3.variables['yield'][99,:,:]

clmtropf=N.flipud(clmtropf)
clmtempf=N.flipud(clmtempf)
clmtropfi=N.flipud(clmtropfi)
clmtempfi=N.flipud(clmtempfi)

#clmtropf= ma.masked_where(maitrop<=0,clmtropf)
#clmtempf= ma.masked_where(maitemp<=0,clmtempf)
clmtropf=ma.filled(clmtropf, fill_value=0.)
clmtempf=ma.filled(clmtempf, fill_value=0.)

#clmtropfi= ma.masked_where(maitropi<=0,clmtropfi)
#clmtempfi= ma.masked_where(maitempi<=0,clmtempfi)
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


area=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/gridareahalf.nc','r')
gridarea = area.variables['cell_area'][:,:]
gridlon = area.variables['lon'][:]

gridarea,gridlon = shiftgrid(180.5,gridarea,gridlon,start=False)

nclu=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/maize_AreaYieldProduction.nc','r')
ncvar_maize = nclu.variables['maizeData'][:]
latnc = nclu.variables['latitude'][:]
znc = nclu.variables['level'][:]
lonnc = nclu.variables['longitude'][:]
timenc = nclu.variables['time'][:]




fig = plt.figure(figsize=(20,15))


ax2 = fig.add_subplot(321)
ax2.set_title("CLM Maize Yield rainfed (t/ha)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
x,y = map(lon2,lat2)
iyield = ma.masked_where(iyield<=0,iyield)
iarea = ma.masked_where(iarea<=0,iarea)
#iyield = ma.masked_where(iarea<=0,iyield)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
#yield_new = ma.masked_where(yield_new<=0,yield_new)
yield_clmtf=maskoceans(x,y,yield_clmtf)

cs1 = map.pcolormesh(x,y,yield_clmtf,cmap=plt.cm.jet,vmin=0,vmax=25)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15) 
plt.axis('off')


ax2 = fig.add_subplot(322)
ax2.set_title("CLM Maize Yield irrigated (t/ha)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(projection='robin',lon_0=0,resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()

#yield_clmtf = ma.masked_where(yield_clmtf<=0,yield_clmtf)

#yield_clmtfi = ma.masked_where(yield_clmtfi<=0,yield_clmtfi)

yield_clmtf=maskoceans(x,y,yield_clmtf)
#yield_clmtf = ma.masked_where(maizeto<=0,yield_clmtf)

yield_clmtfi=maskoceans(x,y,yield_clmtfi)
#yield_clmtfi = ma.masked_where(maizeto<=0,yield_clmtfi)


cs1 = map.pcolormesh(x,y,yield_clmtfi,cmap=plt.cm.jet,vmin=0,vmax=25)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15) 
plt.axis('off')


ax2 = fig.add_subplot(323)
ax2.set_title("CLM Maize Yield rainfed production (t)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
x,y = map(lon2,lat2)
iyield = ma.masked_where(iyield<=0,iyield)
iarea = ma.masked_where(iarea<=0,iarea)
#iyield = ma.masked_where(iarea<=0,iyield)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
#yield_new = ma.masked_where(yield_new<=0,yield_new)
yield_clmtf=maskoceans(x,y,yield_clmtf)
#yield_clmtf = ma.masked_where(yield_clmtf<=0,yield_clmtf)
#yield_clmtfi = ma.masked_where(yield_clmtfi<=0,yield_clmtfi)

cc=yield_clmtf*maizetor*gridarea
cc = ma.masked_where(cc<=0,cc)

cs1 = map.pcolormesh(x,y,cc,cmap=plt.cm.jet,vmin=0,vmax=10000000000)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15)
plt.axis('off')


ax2 = fig.add_subplot(324)
ax2.set_title("CLM Maize Yield irrigated production (t)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(projection='robin',lon_0=0,resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()

cc1=yield_clmtfi*maizetoi*gridarea
cc1 = ma.masked_where(cc1<=0,cc1)

cs1 = map.pcolormesh(x,y,cc1,cmap=plt.cm.jet,vmin=0,vmax=10000000000)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15)
plt.axis('off')

ax2 = fig.add_subplot(325)
ax2.set_title("CLM Maize Yield production (t)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(projection='robin',lon_0=0,resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()

aa=(yield_clmtf*maizetor*gridarea)+(yield_clmtfi*maizetoi*gridarea)
aa = ma.masked_where(aa<=0,aa)
cs1 = map.pcolormesh(x,y,aa,cmap=plt.cm.jet,vmin=0,vmax=10000000000)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15)
plt.axis('off')


plt.savefig('maiexample.jpg',dpi=300,bbox_inches='tight')
plt.show()



