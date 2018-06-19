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
clma=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp45/soytrop_rcp45_co2_rf_nofert_0.5x0.5.nc','r')
#clm=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp45/soytrop_rcp45_co2_rf_fert_0.5x0.5.nc','r')
#clm=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp45/soytrop_rcp45_co2_irrig_fert_0.5x0.5.nc','r')


clmtrop = clm.variables['yield'][4:13,:,:]#2010-2019
clmtropf=N.average(clmtrop,axis=0)
clmtropa = clm.variables['yield'][44:53,:,:]#2050-2059
clmtropfa=N.average(clmtropa,axis=0)

clmtrop1 = clma.variables['yield'][4:13,:,:]#2010-2019
clmtropf1=N.average(clmtrop1,axis=0)
clmtropa1 = clma.variables['yield'][44:53,:,:]#2050-2059
clmtropfa1=N.average(clmtropa1,axis=0)

clm1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp45/soytemp_rcp45_constco2_rf_nofert_0.5x0.5.nc','r')
clm1a=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp45/soytemp_rcp45_co2_rf_nofert_0.5x0.5.nc','r')
#clm1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp45/soytemp_rcp45_co2_rf_fert_0.5x0.5.nc','r')
#clm1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp45/soytemp_rcp45_co2_irrig_fert_0.5x0.5.nc','r')



clmtemp = clm1.variables['yield'][4:13,:,:]
clmtempf=N.average(clmtemp,axis=0)
clmtempa = clm1.variables['yield'][44:53,:,:]
clmtempfa=N.average(clmtempa,axis=0)

clmtemp1 = clm1a.variables['yield'][4:13,:,:]
clmtempf1=N.average(clmtemp1,axis=0)
clmtempa1 = clm1a.variables['yield'][44:53,:,:]
clmtempfa1=N.average(clmtempa1,axis=0)

clmtropf=N.flipud(clmtropf)
clmtempf=N.flipud(clmtempf)
clmtropfa=N.flipud(clmtropfa)
clmtempfa=N.flipud(clmtempfa)

clmtropf1=N.flipud(clmtropf1)
clmtempf1=N.flipud(clmtempf1)
clmtropfa1=N.flipud(clmtropfa1)
clmtempfa1=N.flipud(clmtempfa1)





clmtropf= ma.masked_where(maitoatrop<=0,clmtropf)
clmtempf= ma.masked_where(maitoatemp<=0,clmtempf)
clmtropf=ma.filled(clmtropf, fill_value=0.)
clmtempf=ma.filled(clmtempf, fill_value=0.)

clmtropfa= ma.masked_where(maitoatrop<=0,clmtropfa)
clmtempfa= ma.masked_where(maitoatemp<=0,clmtempfa)
clmtropfa=ma.filled(clmtropfa, fill_value=0.)
clmtempfa=ma.filled(clmtempfa, fill_value=0.)


clmtropf1= ma.masked_where(maitoatrop<=0,clmtropf1)
clmtempf1= ma.masked_where(maitoatemp<=0,clmtempf1)
clmtropf1=ma.filled(clmtropf1, fill_value=0.)
clmtempf1=ma.filled(clmtempf1, fill_value=0.)

clmtropfa1= ma.masked_where(maitoatrop<=0,clmtropfa1)
clmtempfa1= ma.masked_where(maitoatemp<=0,clmtempfa1)
clmtropfa1=ma.filled(clmtropfa1, fill_value=0.)
clmtempfa1=ma.filled(clmtempfa1, fill_value=0.)




clmhis=clmtropf+clmtempf
clmfuture=clmtropfa+clmtempfa
clmhis= ma.masked_where(clmhis[:,:]<=0,clmhis)
clmfuture= ma.masked_where(clmfuture[:,:]<=0,clmfuture)


clmhis1=clmtropf1+clmtempf1
clmfuture1=clmtropfa1+clmtempfa1
clmhis1= ma.masked_where(clmhis1[:,:]<=0,clmhis1)
clmfuture1= ma.masked_where(clmfuture1[:,:]<=0,clmfuture1)


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
yieldf1= N.zeros((10, 360, 720))
yieldf21= N.zeros((10, 360, 720))
yieldfet= N.zeros((10, 360, 720))
yieldf2et= N.zeros((10, 360, 720))
yieldf1et= N.zeros((10, 360, 720))
yieldf21et= N.zeros((10, 360, 720))

years = range(2010, 2020)
years2 = range(2050,2060)

for i, year in enumerate(years):
   
    base = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_future/rcp45as/output/rcp45as.bgp-yearly_crop_{0}.nc".format(year), mode='r')
    basea = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_future/rcp45bs/output/rcp45bs.bgp-yearly_crop_{0}.nc".format(year), mode='r')
#    base = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_future/rcp45cs/output/rcp45cs.bgp-yearly_crop_{0}.nc".format(year), mode='r')
#    base = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_future/rcp45ds/output/rcp45ds.bgp-yearly_crop_{0}.nc".format(year), mode='r')
    lona1 = base.variables["lon"][:]
    lata1 = base.variables["lat"][:]
   
    yield1 = base.variables["yield"][1,:,:]
    yieldf[i, :, :] = yield1
    yield11 = basea.variables["yield"][1,:,:]
    yieldf1[i, :, :] = yield11

    yield1e = base.variables["g_ET"][1,:,:]
    yieldfet[i, :, :] = yield1e
    yield11e = basea.variables["g_ET"][1,:,:]
    yieldf1et[i, :, :] = yield11e

    
for i, year1 in enumerate(years2):
    base2 = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_future/rcp45as/output/rcp45as.bgp-yearly_crop_{0}.nc".format(year1), mode='r')
    base2a = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_future/rcp45bs/output/rcp45bs.bgp-yearly_crop_{0}.nc".format(year1), mode='r')
#    base2 = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_future/rcp45cs/output/rcp45cs.bgp-yearly_crop_{0}.nc".format(year1), mode='r')
#    base2 = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_future/rcp45ds/output/rcp45ds.bgp-yearly_crop_{0}.nc".format(year1), mode='r')    
    yield2 = base2.variables["yield"][1,:,:]
    yieldf2[i, :, :] = yield2

    yield21 = base2a.variables["yield"][1,:,:]
    yieldf21[i, :, :] = yield21

    yield2e = base2.variables["g_ET"][1,:,:]
    yieldf2et[i, :, :] = yield2e

    yield21e = base2a.variables["g_ET"][1,:,:]
    yieldf21et[i, :, :] = yield21e



yielda=N.average(yieldf,axis=0)
yielda2=N.average(yieldf2,axis=0)

yieldaet=N.average(yieldfet,axis=0)
yielda2et=N.average(yieldf2et,axis=0)


yield_new,lona11 = shiftgrid(180.5,yielda,lona1,start=False)
yield_new2,lona11 = shiftgrid(180.5,yielda2,lona1,start=False)

yield_newet,lona11 = shiftgrid(180.5,yieldaet,lona1,start=False)
yield_new2et,lona11 = shiftgrid(180.5,yielda2et,lona1,start=False)

yielda1=N.average(yieldf1,axis=0)
yielda21=N.average(yieldf21,axis=0)

yielda1et=N.average(yieldf1et,axis=0)
yielda21et=N.average(yieldf21et,axis=0)



yield_new1,lona11 = shiftgrid(180.5,yielda1,lona1,start=False)
yield_new21,lona11 = shiftgrid(180.5,yielda21,lona1,start=False)

yield_new1et,lona11 = shiftgrid(180.5,yielda1et,lona1,start=False)
yield_new21et,lona11 = shiftgrid(180.5,yielda21et,lona1,start=False)

lon2,lat2 = N.meshgrid(lona11,lata1)

yield_new= ma.masked_where( clmhis[:,:]<=0,yield_new)
yield_new2=ma.masked_where( clmfuture[:,:]<=0,yield_new2)

yield_new1= ma.masked_where( clmhis1[:,:]<=0,yield_new1)
yield_new21=ma.masked_where( clmfuture1[:,:]<=0,yield_new21)

yield_newet= ma.masked_where( clmhis[:,:]<=0,yield_newet)
yield_new2et=ma.masked_where( clmfuture[:,:]<=0,yield_new2et)

yield_new1et= ma.masked_where( clmhis1[:,:]<=0,yield_new1et)
yield_new21et=ma.masked_where( clmfuture1[:,:]<=0,yield_new21et)


t, pT = ttest_ind(yieldf, yieldf2, axis = 0, equal_var = False)
t1,lona11 = shiftgrid(180.5,t,lona1,start=False)
pT1,lona11 = shiftgrid(180.5,pT,lona1,start=False)

yieldisam=yield_new2-yield_new
yieldisam1= ma.masked_where( pT1[:,:]>0.1,yieldisam)
yieldisam= ma.masked_where(yieldisam==0.,yieldisam)

fig = plt.figure(figsize=(12,6))


ax1 = fig.add_subplot(221)
ax1.set_title("CLM yield difference 2010s (t/ha)",fontsize=18)

map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
x,y = map(lon2,lat2)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
clmhis= ma.masked_where(clmhis<=0.,clmhis)
clmhis1= ma.masked_where(clmhis1<=0.,clmhis1)
cc=clmhis1-clmhis
cc= ma.masked_where(cc==0.,cc)

cs = map.pcolormesh(x,y,cc,cmap=plt.cm.bwr,vmin=-1,vmax=1)
cbar = map.colorbar(cs,location='bottom',size="4%",pad="2%")
cbar.ax.tick_params(labelsize=12) 
plt.axis('off')

ax2 = fig.add_subplot(222)
ax2.set_title("ISAM CWP relative change 2010s (%)",fontsize=18)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
yield_new= ma.masked_where(yield_new<=0.,yield_new)
yield_new1= ma.masked_where(yield_new1<=0.,yield_new1)
cisa=(yield_new1/yield_new1et)-(yield_new/yield_newet)
cisa= ma.masked_where(cisa==0.,cisa)


cs = map.pcolormesh(x,y,cisa/(yield_new/yield_newet)*100,cmap=plt.cm.YlGnBu,vmin=0,vmax=40)
cbar = map.colorbar(cs,location='bottom',size="4%",pad="2%")
cbar.ax.tick_params(labelsize=12)
plt.axis('off')

ax5 = fig.add_subplot(223)
ax5.set_title("CLM yield difference 2050s (t/ha)",fontsize=18)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
clmfuture= ma.masked_where(clmfuture<=0.,clmfuture)
clmfuture1= ma.masked_where(clmfuture1<=0.,clmfuture1)
cas=clmfuture1-clmfuture
cas= ma.masked_where(cas==0.,cas)

cs = map.pcolormesh(x,y,cas,cmap=plt.cm.bwr,vmin=-1,vmax=1)
cbar = map.colorbar(cs,location='bottom',size="4%",pad="2%")
cbar.ax.tick_params(labelsize=12)
plt.axis('off')



ax5 = fig.add_subplot(224)
ax5.set_title("ISAM CWP relative change 2050s (%)",fontsize=18)
map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
yield_new2= ma.masked_where(yield_new2<=0.,yield_new2)
yield_new21= ma.masked_where(yield_new21<=0.,yield_new21)
cas=(yield_new21/yield_new21et)-(yield_new2/yield_new2et)
cas= ma.masked_where(cas==0.,cas)
cs = map.pcolormesh(x,y,cas/(yield_new2/yield_new2et)*100,cmap=plt.cm.YlGnBu,vmin=0,vmax=40)
cbar = map.colorbar(cs,location='bottom',size="4%",pad="2%")
cbar.ax.tick_params(labelsize=12)
plt.axis('off')

plt.savefig('scp45soyybdiffcwp.jpg',dpi=300,bbox_inches='tight')
plt.show()

