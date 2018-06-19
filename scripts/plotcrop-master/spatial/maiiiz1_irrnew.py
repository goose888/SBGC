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



yieldfi= N.zeros((1, 360, 720))
yieldf2i= N.zeros((1, 360, 720))

yieldf= N.zeros((1, 360, 720))
yieldf2= N.zeros((1, 360, 720))
years = range(2000, 2001)
for i, year in enumerate(years):
    base = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_his/cruhis/output/cruhis.bgp-yearly_crop_{0}.nc".format(year), mode='r')
    lona1 = base.variables["lon"][:]
    lata1 = base.variables["lat"][:]
   
    yield1 = base.variables["yield"][0,:,:]
    yieldf[i, :, :] = yield1

    basei = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_his/cruhis_irr/output/cruhisirr.bgp-yearly_crop_{0}.nc".format(year), mode='r')
    yield2 = basei.variables["yield"][0,:,:]
    yieldf2[i, :, :] = yield2

yielda=N.average(yieldf,axis=0)
yielda2=N.average(yieldf2,axis=0)

yield_new,lona11 = shiftgrid(180.5,yielda,lona1,start=False)
yield_new2,lona11 = shiftgrid(180.5,yielda2,lona1,start=False)

lon2,lat2 = N.meshgrid(lona11,lata1)



fig, ax = plt.subplots(figsize=(8,6))

ax.set_title("Iizumi Maize Yield (t/ha)",fontsize=20)

map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
x,y = map(lon2,lat2)
iyield = ma.masked_where(iyield<=0,iyield)
iarea = ma.masked_where(iarea<=0,iarea)
#iyield = ma.masked_where(iarea<=0,iyield)
iyield = ma.masked_where(maizeto<=0,iyield)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
iizumy=iyield
cs = map.pcolormesh(x,y,iizumy,cmap=plt.cm.YlGn,vmin=0,vmax=16)

cbar = map.colorbar(cs,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15)
plt.axis('off')
plt.savefig('maiiizirr.jpg',dpi=300,bbox_inches='tight')

fig = plt.figure(figsize=(20,15))


ax2 = fig.add_subplot(321)
ax2.set_title("ISAM-NCEP Maize Yield (t/ha)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
x,y = map(lon2,lat2)
iyield = ma.masked_where(iyield<=0,iyield)
iarea = ma.masked_where(iarea<=0,iarea)
#iyield = ma.masked_where(iarea<=0,iyield)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
#yield_new = ma.masked_where(yield_new<=0,yield_new)
yield_new=maskoceans(x,y,yield_new)
yield_new=ma.filled(yield_new, fill_value=0.)
yield_new = ma.masked_where(yield_new<=0,yield_new)
yield_new = ma.masked_where(maizeto<=0,yield_new)

yield_new2=maskoceans(x,y,yield_new2)
yield_new2=ma.filled(yield_new2, fill_value=0.)
yield_new2 = ma.masked_where(yield_new2<=0,yield_new2)
yield_new2 = ma.masked_where(maizeto<=0,yield_new2)

isamy=((yield_new*maizetor*gridarea)+(yield_new2*maizetoi*gridarea))/((maizetoi*gridarea)+(maizetor*gridarea))
isamy = ma.masked_where(iizumy<=0,isamy)

cs1 = map.pcolormesh(x,y,isamy,cmap=plt.cm.YlGn,vmin=0,vmax=16)
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

#yield_clmtf = ma.masked_where(yield_clmtf<=0,yield_clmtf)

#yield_clmtfi = ma.masked_where(yield_clmtfi<=0,yield_clmtfi)

yield_clmtf=maskoceans(x,y,yield_clmtf)
#yield_clmtf=ma.filled(yield_clmtf, fill_value=0.)
#yield_clmtf = ma.masked_where(yield_clmtf<=0,yield_clmtf)
yield_clmtf = ma.masked_where(maizeto<=0,yield_clmtf)

yield_clmtfi=maskoceans(x,y,yield_clmtfi)
#yield_clmtfi=ma.filled(yield_clmtfi, fill_value=0.)
#yield_clmtfi = ma.masked_where(yield_clmtfi<=0,yield_clmtfi)
yield_clmtfi = ma.masked_where(maizeto<=0,yield_clmtfi)

clmy=((yield_clmtf*maizetor*gridarea)+(yield_clmtfi*maizetoi*gridarea))/((maizetoi*gridarea)+(maizetor*gridarea))
clmy = ma.masked_where(iizumy<=0,clmy)

cs1 = map.pcolormesh(x,y,clmy,cmap=plt.cm.YlGn,vmin=0,vmax=16)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15) 
plt.axis('off')
#print N.max(yield_fine*ncvar_maize1*1000/gridarea)



ax2 = fig.add_subplot(323)
ax2.set_title("ISAM-Iizumi Maize Yield (t/ha)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(projection='robin',lon_0=0,resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cd=isamy-iizumy
cd = ma.masked_where(cd==0,cd)

cs1 = map.pcolormesh(x,y,cd,cmap=plt.cm.bwr,vmin=-5,vmax=5)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15) 
plt.axis('off')


ax2 = fig.add_subplot(324)
ax2.set_title("CLM-Iizumi Maize Yield (t/ha)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(projection='robin',lon_0=0,resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,clmy-iizumy,cmap=plt.cm.bwr,vmin=-5,vmax=5)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15) 
plt.axis('off')



ax2 = fig.add_subplot(325)
ax2.set_title("ISAM-Iizumi Maize Yield gridcell (%)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(projection='robin',lon_0=0,resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs1 = map.pcolormesh(x,y,(isamy-iizumy)/iizumy*100,cmap=plt.cm.bwr,vmin=-100,vmax=100)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15)
plt.axis('off')


ax2 = fig.add_subplot(326)
ax2.set_title("CLM-Iizumi Maize Yield gridcell (%)",fontsize=20)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
#map = Basemap(projection='robin',lon_0=0,resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
#diff=(yield_clmtf*iarea/100/10*1000)-(iyield*iarea/100/10*1000)
cs1 = map.pcolormesh(x,y,(clmy-iizumy)/iizumy*100,cmap=plt.cm.bwr,vmin=-100,vmax=100)

cbar = map.colorbar(cs1,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15)
plt.axis('off')



plt.savefig('maiiizirr1_ncep.jpg',dpi=300,bbox_inches='tight')
plt.show()


#fig, ax = plt.subplots(figsize=(6,6))
colors = (0,0,1)
colorsr = (1,0,0)
fig = plt.figure(figsize=(8,12))
ax = fig.add_subplot(211)

#ax.plot([0,1000],[0,1000], 'k--',label='1:1')

#ax.scatter(iizumy, isamy, c=colors,alpha=1,label='ISAM')
ax.scatter(iizumy, clmy, c=colorsr,alpha=0.5,label='CLM')
plt.xlim(0, 16)
plt.ylim(0, 16)
ax.plot([0,16],[0,16], 'k--',label='1:1')

#ax.set_title('Maize yield over gridcell',fontsize=18)
#ax.legend()

plt.tick_params(axis='both',labelsize=15)
plt.xlabel('Iizumi-Crops (t/ha)',fontsize=18)
plt.ylabel('Model (t/ha)',fontsize=18)


ax = fig.add_subplot(212)

#ax.plot([0,1000],[0,1000], 'k--',label='1:1')

ax.scatter(iizumy, isamy, c=colors,alpha=0.5,label='ISAM')
#ax.scatter(iizumy, clmy, c=colorsr,alpha=0.5,label='CLM')
ax.plot([0,16],[0,16], 'k--',label='1:1')

plt.xlim(0, 16)
plt.ylim(0, 16)
#ax.set_title('Maize yield over gridcell',fontsize=18)
#ax.legend()
plt.xticks([])
plt.tick_params(axis='both',labelsize=15)
#plt.xlabel('Iizumi-Crops (g/$\mathregular{m^2}$)',fontsize=18)
plt.ylabel('Model (t/ha)',fontsize=18)

plt.savefig('scatter_maiiizirr1_ncep.png',bbox_inches='tight')
#plt.show()


