from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import glob
import numpy.ma as ma
from scipy.interpolate import griddata
from scipy.stats import ttest_ind

region1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/RCP45_crop_150901.nc','r')
maitrop = region1.variables['soy_trop'][4,:,:]
maitemp = region1.variables['soy_temp'][4,:,:]
maitropi = region1.variables['soy_trop_irrig'][4,:,:]
maitempi = region1.variables['soy_temp_irrig'][4,:,:]
maitrop= ma.masked_where(maitrop<=0,maitrop)
maitropi= ma.masked_where(maitropi<=0,maitropi)
maitemp= ma.masked_where(maitemp<=0,maitemp)
maitempi= ma.masked_where(maitempi<=0,maitempi)
maitrop=ma.filled(maitrop, fill_value=0.)
maitropi=ma.filled(maitropi,fill_value=0.)
maitemp=ma.filled(maitemp, fill_value=0.)
maitempi=ma.filled(maitempi, fill_value=0.)
maizeto = maitrop+maitemp
maizetoi = maitropi+maitempi
maitoatemp=maitemp+maitempi
maitoatrop=maitrop+maitropi
maizetotal = maizeto+maizetoi
clm=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp45/soytrop_rcp45_constco2_rf_nofert_0.5x0.5.nc','r')
#clm=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp45/maizetrop_rcp45_co2_rf_nofert_0.5x0.5.nc','r')
#clm=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp45/maizetrop_rcp45_co2_rf_fert_0.5x0.5.nc','r')
#clm=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp45/maizetrop_rcp45_co2_irrig_fert_0.5x0.5.nc','r')


clmtrop = clm.variables['yield'][4:13,:,:]#2010-2019
clmtropf=N.average(clmtrop,axis=0)
clmtropa = clm.variables['yield'][44:53,:,:]#2050-2059
clmtropfa=N.average(clmtropa,axis=0)


clm1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp45/soytemp_rcp45_constco2_rf_nofert_0.5x0.5.nc','r')
#clm1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp45/maizetemp_rcp45_co2_rf_nofert_0.5x0.5.nc','r')
#clm1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp45/maizetemp_rcp45_co2_rf_fert_0.5x0.5.nc','r')
#clm1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp45/maizetemp_rcp45_co2_irrig_fert_0.5x0.5.nc','r')



clmtemp = clm1.variables['yield'][4:13,:,:]
clmtempf=N.average(clmtemp,axis=0)
clmtempa = clm1.variables['yield'][44:53,:,:]
clmtempfa=N.average(clmtempa,axis=0)

clmtropf=N.flipud(clmtropf)
clmtempf=N.flipud(clmtempf)
clmtropfa=N.flipud(clmtropfa)
clmtempfa=N.flipud(clmtempfa)

clmtropf= ma.masked_where(maitoatrop<=0,clmtropf)
clmtempf= ma.masked_where(maitoatemp<=0,clmtempf)
clmtropf=ma.filled(clmtropf, fill_value=0.)
clmtempf=ma.filled(clmtempf, fill_value=0.)

clmtropfa= ma.masked_where(maitoatrop<=0,clmtropfa)
clmtempfa= ma.masked_where(maitoatemp<=0,clmtempfa)
clmtropfa=ma.filled(clmtropfa, fill_value=0.)
clmtempfa=ma.filled(clmtempfa, fill_value=0.)

clmhis=clmtropf+clmtempf
clmfuture=clmtropfa+clmtempfa
clmhis= ma.masked_where(clmhis[:,:]<=0,clmhis)
clmfuture= ma.masked_where(clmfuture[:,:]<=0,clmfuture)


clmhist=clmtrop+clmtemp
clmfutt=clmtropa+clmtempa

tc, pTc = ttest_ind(clmhist,clmfutt, axis = 0, equal_var = False)

tc=N.flipud(tc)
pTc=N.flipud(pTc)

yieldclm=clmfuture-clmhis
yieldclm= ma.masked_where(yieldclm==0.,yieldclm)

yieldclm1= ma.masked_where( pTc[:,:]>0.1,yieldclm)


yieldf= N.zeros((10, 360, 720))
yieldf2= N.zeros((10, 360, 720))
yieldfa= N.zeros((10, 360, 720))
yieldf2a= N.zeros((10, 360, 720))
yieldfb= N.zeros((10, 360, 720))
yieldf2b= N.zeros((10, 360, 720))
yieldfc= N.zeros((10, 360, 720))
yieldf2c= N.zeros((10, 360, 720))



years = range(2010, 2020)
years2 = range(2050,2060)

for i, year in enumerate(years):
   
    base = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_future/rcp45as/output/rcp45as.bgp-yearly_crop_{0}.nc".format(year), mode='r')
    basea = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_future/rcp45bs/output/rcp45bs.bgp-yearly_crop_{0}.nc".format(year), mode='r')
    baseb = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_future/rcp45cs/output/rcp45cs.bgp-yearly_crop_{0}.nc".format(year), mode='r')
    basec = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_future/rcp45ds/output/rcp45ds.bgp-yearly_crop_{0}.nc".format(year), mode='r')
    lona1 = base.variables["lon"][:]
    lata1 = base.variables["lat"][:]
   
    yield1 = base.variables["g_Temp"][7,1,:,:]
    yieldf[i, :, :] = yield1
    yield1a = basea.variables["g_Temp"][7,1,:,:]
    yieldfa[i, :, :] = yield1a
    yield1b = baseb.variables["g_Temp"][7,1,:,:]
    yieldfb[i, :, :] = yield1b
    yield1c = basec.variables["g_Temp"][7,1,:,:]
    yieldfc[i, :, :] = yield1c

    
for i, year1 in enumerate(years2):
    base2 = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_future/rcp45as/output/rcp45as.bgp-yearly_crop_{0}.nc".format(year1), mode='r')
    base2a = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_future/rcp45bs/output/rcp45bs.bgp-yearly_crop_{0}.nc".format(year1), mode='r')
    base2b = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_future/rcp45cs/output/rcp45cs.bgp-yearly_crop_{0}.nc".format(year1), mode='r')
    base2c = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_future/rcp45ds/output/rcp45ds.bgp-yearly_crop_{0}.nc".format(year1), mode='r')    
    yield2 = base2.variables["g_Temp"][7,1,:,:]
    yieldf2[i, :, :] = yield2
    yield2a = base2a.variables["g_Temp"][7,1,:,:]
    yieldf2a[i, :, :] = yield2a
    yield2b = base2b.variables["g_Temp"][7,1,:,:]
    yieldf2b[i, :, :] = yield2b
    yield2c = base2c.variables["g_Temp"][7,1,:,:]
    yieldf2c[i, :, :] = yield2c



yielda=N.average(yieldf,axis=0)
yielda2=N.average(yieldf2,axis=0)
yieldab=N.average(yieldfa,axis=0)
yielda2b=N.average(yieldf2a,axis=0)
yieldac=N.average(yieldfb,axis=0)
yielda2c=N.average(yieldf2b,axis=0)
yieldad=N.average(yieldfc,axis=0)
yielda2d=N.average(yieldf2c,axis=0)


yield_new,lona11 = shiftgrid(180.5,yielda,lona1,start=False)
yield_new2,lona11 = shiftgrid(180.5,yielda2,lona1,start=False)
yield_newb,lona11 = shiftgrid(180.5,yieldab,lona1,start=False)
yield_new2b,lona11 = shiftgrid(180.5,yielda2b,lona1,start=False)
yield_newc,lona11 = shiftgrid(180.5,yieldac,lona1,start=False)
yield_new2c,lona11 = shiftgrid(180.5,yielda2c,lona1,start=False)
yield_newd,lona11 = shiftgrid(180.5,yieldad,lona1,start=False)
yield_new2d,lona11 = shiftgrid(180.5,yielda2d,lona1,start=False)







lon2,lat2 = N.meshgrid(lona11,lata1)

yield_new= ma.masked_where( clmhis[:,:]<=0,yield_new)
yield_new2=ma.masked_where( clmfuture[:,:]<=0,yield_new2)

yield_newb= ma.masked_where( clmhis[:,:]<=0,yield_newb)
yield_new2b=ma.masked_where( clmfuture[:,:]<=0,yield_new2b)
yield_newc= ma.masked_where( clmhis[:,:]<=0,yield_newc)
yield_new2c=ma.masked_where( clmfuture[:,:]<=0,yield_new2c)
yield_newd= ma.masked_where( clmhis[:,:]<=0,yield_newd)
yield_new2d=ma.masked_where( clmfuture[:,:]<=0,yield_new2d)



t, pT = ttest_ind(yieldf, yieldf2, axis = 0, equal_var = False)
t1,lona11 = shiftgrid(180.5,t,lona1,start=False)
pT1,lona11 = shiftgrid(180.5,pT,lona1,start=False)

yieldisam=yield_new2-yield_new
yieldisam1= ma.masked_where( pT1[:,:]>0.1,yieldisam)
yieldisam= ma.masked_where(yieldisam==0.,yieldisam)

fig = plt.figure(figsize=(12,6))


ax1 = fig.add_subplot(221)
ax1.set_title("CO2 (K/yr)",fontsize=18)

map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
x,y = map(lon2,lat2)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs = map.pcolormesh(x,y,yield_new2-yield_new,cmap=plt.cm.seismic,vmin=-8,vmax=8)
cbar = map.colorbar(cs,location='bottom',size="4%",pad="2%")
cbar.ax.tick_params(labelsize=12) 
plt.axis('off')

ax2 = fig.add_subplot(222)
ax2.set_title("CO2+CLIMATE (K/yr)",fontsize=18)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs = map.pcolormesh(x,y,yield_new2b-yield_newb,cmap=plt.cm.seismic,vmin=-8,vmax=8)
cbar = map.colorbar(cs,location='bottom',size="4%",pad="2%")
cbar.ax.tick_params(labelsize=12)
plt.axis('off')

ax5 = fig.add_subplot(223)
ax5.set_title("CO2+CLIMATE+N (K/yr)",fontsize=18)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs = map.pcolormesh(x,y,yield_new2c-yield_newc,cmap=plt.cm.seismic,vmin=-8,vmax=8)
cbar = map.colorbar(cs,location='bottom',size="4%",pad="2%")
cbar.ax.tick_params(labelsize=12)
plt.axis('off')



ax5 = fig.add_subplot(224)
ax5.set_title("CO2+CLIMATE+N+I (K/yr)",fontsize=18)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs = map.pcolormesh(x,y,yield_new2d-yield_newd,cmap=plt.cm.seismic,vmin=-8,vmax=8)
cbar = map.colorbar(cs,location='bottom',size="4%",pad="2%")
cbar.ax.tick_params(labelsize=12)
plt.axis('off')

plt.savefig('scp45soyt.jpg',dpi=300,bbox_inches='tight')
plt.show()

