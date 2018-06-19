from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp,maskoceans
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import glob
import numpy.ma as ma
from scipy.interpolate import griddata
import matplotlib.colors as colors


iizumi=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/iizumi/iizumi.2013JAN29.soybean.1982-2006.30min.nc4','r')
#print iizumi
iyield = iizumi.variables['yield50'][18,:,:]
iarea =iizumi.variables['area'][18,:,:]
la=iizumi.variables['lat'][:]
lo=iizumi.variables['lon'][:]
iyield=N.flipud(iyield)
iarea=N.flipud(iarea)

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

nclu=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/soybean_AreaYieldProduction.nc','r')
ncvar_maize = nclu.variables['soybeanData'][:]
latnc = nclu.variables['latitude'][:]
znc = nclu.variables['level'][:]
lonnc = nclu.variables['longitude'][:]
timenc = nclu.variables['time'][:]



yieldf= N.zeros((1, 360, 720))
yieldf2= N.zeros((1, 360, 720))
years = range(2000, 2001)
for i, year in enumerate(years):
    base = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_his/cruhis/output/cruhis.bgp-yearly_crop_{0}.nc".format(year), mode='r')
    lona1 = base.variables["lon"][:]
    lata1 = base.variables["lat"][:]
   
    yield1 = base.variables["yield"][1,:,:]
    yieldf[i, :, :] = yield1
    #yield2 = base2.variables["totalyield"][:]
    #yieldf2[i, :, :] = yield2

yielda=N.average(yieldf,axis=0)
#yielda2=N.average(yieldf2,axis=0)

yield_new,lona11 = shiftgrid(180.5,yielda,lona1,start=False)
#yield_new2,lona11 = shiftgrid(180.,yielda2,lona1,start=False)

lon2,lat2 = N.meshgrid(lona11,lata1)
#print lon2.shape
#ncvar_maize2 = interp(ncvar_maizea,lonnc,lat_new,lon,lat,order=1)



fig, ax = plt.subplots(figsize=(8,6))

ax.set_title("Iizumi Soybean Yield gridecll (g/$\mathregular{m^2}$)",fontsize=20)

map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
x,y = map(lon2,lat2)
iyield = ma.masked_where(iyield<=0,iyield)
iarea = ma.masked_where(iarea<=0,iarea)
iyield = ma.masked_where(iarea<=0,iyield)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs = map.pcolormesh(x,y,iyield*iarea/100/10*1000,cmap=plt.cm.YlGn,norm=colors.PowerNorm(gamma=1./2.),vmin=0,vmax=50)

cbar = map.colorbar(cs,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15)
plt.axis('off')
plt.savefig('soyiiz.jpg',dpi=300,bbox_inches='tight')

fig = plt.figure(figsize=(20,15))


ax2 = fig.add_subplot(321)
ax2.set_title("ISAM-NCEP Soybean Yield gridcell (g/$\mathregular{m^2}$)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
x,y = map(lon2,lat2)
iyield = ma.masked_where(iyield<=0,iyield)
iarea = ma.masked_where(iarea<=0,iarea)
iyield = ma.masked_where(iarea<=0,iyield)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
#yield_new = ma.masked_where(yield_new<=0,yield_new)
yield_new=maskoceans(x,y,yield_new)
#yield_new[N.isnan(yield_new)] = -9999
yield_new=ma.filled(yield_new, fill_value=0.)

#yield_new[N.isnan(yield_new)] = -9999
yield_new = ma.masked_where(yield_new<=0,yield_new)
yield_new = ma.masked_where(iyield<=0,yield_new)
#cs1 = map.pcolormesh(x,y,yield_new*iarea/100/10,cmap=plt.cm.YlGn,vmin=0.01,vmax=0.15)
cs1 = map.pcolormesh(x,y,yield_new*iarea/100/10*1000,cmap=plt.cm.YlGn,norm=colors.PowerNorm(gamma=1./2.),vmin=0,vmax=50)
cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15) 
plt.axis('off')
#print N.max(yield_fine*ncvar_maize1*1000/gridarea)


ax2 = fig.add_subplot(322)
ax2.set_title("CLM Soybean Yield gridcell (g/$\mathregular{m^2}$)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(projection='robin',lon_0=0,resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()

yield_clmtf=clmtropf+clmtempf
yield_clmtf = ma.masked_where(yield_clmtf<=0,yield_clmtf)
yield_clmtf = ma.masked_where(iyield<=0,yield_clmtf)

#cs1 = map.pcolormesh(x,y,yield_clmtf*mask_clm/10,cmap=plt.cm.YlGn,vmin=0.0,vmax=0.35)

#cs1 = map.pcolormesh(x,y,yield_clmtf*iarea/100/10,cmap=plt.cm.YlGn,vmin=0.01,vmax=0.15)
cs1 = map.pcolormesh(x,y,yield_clmtf*iarea/100/10*1000,cmap=plt.cm.YlGn,norm=colors.PowerNorm(gamma=1./2.),vmin=0,vmax=50)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15) 
plt.axis('off')
#print N.max(yield_fine*ncvar_maize1*1000/gridarea)



ax2 = fig.add_subplot(323)
ax2.set_title("ISAM-Iizumi Soybean Yield gridcell (g/$\mathregular{m^2}$)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(projection='robin',lon_0=0,resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,(yield_new*iarea/100/10*1000)-(iyield*iarea/100/10*1000),cmap=plt.cm.bwr,vmin=-50,vmax=50)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15) 
plt.axis('off')


ax2 = fig.add_subplot(324)
ax2.set_title("CLM-Iizumi Soybean Yield gridcell (g/$\mathregular{m^2}$)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(projection='robin',lon_0=0,resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,(yield_clmtf*iarea/100/10*1000)-(iyield*iarea/100/10*1000),cmap=plt.cm.bwr,vmin=-50,vmax=50)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15) 
plt.axis('off')



ax2 = fig.add_subplot(325)
ax2.set_title("ISAM-Iizumi Soybean Yield gridcell (%)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(projection='robin',lon_0=0,resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,((yield_new*iarea/100/10*1000)-(iyield*iarea/100/10*1000))/(iyield*iarea/100/10*1000)*100,cmap=plt.cm.bwr,vmin=-100,vmax=100)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15)
plt.axis('off')


ax2 = fig.add_subplot(326)
ax2.set_title("CLM-Iizumi Soybean Yield gridcell (%)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(projection='robin',lon_0=0,resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
diff=(yield_clmtf*iarea/100/10*1000)-(iyield*iarea/100/10*1000)
cs1 = map.pcolormesh(x,y,diff/(iyield*iarea/100/10*1000)*100,cmap=plt.cm.bwr,vmin=-100,vmax=100)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15)
plt.axis('off')



plt.savefig('soyiiz_ncep.jpg',dpi=300,bbox_inches='tight')
#plt.show()


fig, ax = plt.subplots(figsize=(6,6))
colors = (0,0,1)
colorsr = (1,0,0)

ax.plot([0,250],[0,250], 'k--',label='1:1')

ax.scatter(iyield*iarea/100/10*1000, yield_new*iarea/100/10*1000, c=colors,alpha=1,label='ISAM')
ax.scatter(iyield*iarea/100/10*1000, yield_clmtf*iarea/100/10*1000, c=colorsr,alpha=0.5,label='CLM')
plt.xlim(0, 250)
plt.ylim(0, 250)
ax.set_title('Soybean yield over gridcell',fontsize=18)
ax.legend()
plt.tick_params(axis='both',labelsize=15)

plt.xlabel('Iizumi-Crops (g/$\mathregular{m^2}$)',fontsize=18)
plt.ylabel('Model (g/$\mathregular{m^2}$)',fontsize=18)
plt.savefig('scatter_soyiiz_ncep.png',bbox_inches='tight')
#plt.show()


