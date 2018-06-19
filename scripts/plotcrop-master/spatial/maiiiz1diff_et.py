from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp,maskoceans
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import glob
import numpy.ma as ma
from scipy.interpolate import griddata
import matplotlib.colors as colors
from scipy.stats import ttest_ind


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



yieldfi= N.zeros((100, 360, 720))
yieldf2i= N.zeros((100, 360, 720))

yieldf= N.zeros((100, 360, 720))
yieldf2= N.zeros((100, 360, 720))
years = range(1901, 2001)
for i, year in enumerate(years):
    base = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_his/cruhis/output/cruhis.bgp-yearly_crop_{0}.nc".format(year), mode='r')
    lona1 = base.variables["lon"][:]
    lata1 = base.variables["lat"][:]
   
    yield1 = base.variables["g_ET"][0,:,:]
    yieldf[i, :, :] = yield1

    basei = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_his/cruhis_irr/output/cruhisirr.bgp-yearly_crop_{0}.nc".format(year), mode='r')
    yield2 = basei.variables["g_ET"][0,:,:]
    yieldf2[i, :, :] = yield2

    base1 = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_his/cesmhis/output/cesmhis.bgp-yearly_crop_{0}.nc".format(year), mode='r')

    yield12 = base1.variables["g_ET"][0,:,:]
    yieldfi[i, :, :] = yield12

    basei1 = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_his/cesmhis_irr/output/cesmhisirr.bgp-yearly_crop_{0}.nc".format(year), mode='r')
    yield23 = basei1.variables["g_ET"][0,:,:]
    yieldf2i[i, :, :] = yield23

yielda=N.average(yieldf,axis=0)
yielda2=N.average(yieldf2,axis=0)

yield_new,lona11 = shiftgrid(180.5,yielda,lona1,start=False)
yield_new2,lona11 = shiftgrid(180.5,yielda2,lona1,start=False)

yieldai=N.average(yieldfi,axis=0)
yielda2i=N.average(yieldf2i,axis=0)

yield_newf,lona11 = shiftgrid(180.5,yieldai,lona1,start=False)
yield_new2f,lona11 = shiftgrid(180.5,yielda2i,lona1,start=False)

lon2,lat2 = N.meshgrid(lona11,lata1)



fig, ax = plt.subplots(figsize=(12,12))

ax.set_title("ISAM CESM-NCEP Maize ET (mm/yr)",fontsize=18)

map = Basemap(projection ='cyl', llcrnrlat=-65, urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
x,y = map(lon2,lat2)
iyield = ma.masked_where(iyield<=0,iyield)
iarea = ma.masked_where(iarea<=0,iarea)
#iyield = ma.masked_where(iarea<=0,iyield)
iyield = ma.masked_where(maizeto<=0,iyield)
iizumy=iyield*maizeto*100


yield_new=maskoceans(x,y,yield_new)
yield_new=ma.filled(yield_new, fill_value=0.)
yield_new = ma.masked_where(yield_new<=0,yield_new)
yield_new = ma.masked_where(maizeto<=0,yield_new)

yield_new2=maskoceans(x,y,yield_new2)
yield_new2=ma.filled(yield_new2, fill_value=0.)
yield_new2 = ma.masked_where(yield_new2<=0,yield_new2)
yield_new2 = ma.masked_where(maizeto<=0,yield_new2)

isamy=((yield_new*maizetor*gridarea)+(yield_new2*maizetoi*gridarea))*100/gridarea
isamy = ma.masked_where(iizumy<=0,isamy)

yield_newf=maskoceans(x,y,yield_newf)
yield_newf=ma.filled(yield_newf, fill_value=0.)
yield_newf = ma.masked_where(yield_newf<=0,yield_newf)
yield_newf = ma.masked_where(maizeto<=0,yield_newf)

yield_new2f=maskoceans(x,y,yield_new2f)
yield_new2f=ma.filled(yield_new2f, fill_value=0.)
yield_new2f = ma.masked_where(yield_new2f<=0,yield_new2f)
yield_new2f = ma.masked_where(maizeto<=0,yield_new2f)

isamyf=((yield_newf*maizetor*gridarea)+(yield_new2f*maizetoi*gridarea))*100/gridarea
isamyf = ma.masked_where(iizumy<=0,isamyf)
cc=yield_newf-yield_new
cc = ma.masked_where(iizumy<=0,cc)
tc, pTc = ttest_ind(yieldf,yieldfi, axis = 0, equal_var = False)
pTc,lona11 = shiftgrid(180.5,pTc,lona1,start=False)



cc= ma.masked_where( pTc[:,:]>0.1,cc)

map.drawcoastlines()
map.drawcountries()
map.drawmapboundary()
cs = map.pcolormesh(x,y,cc,cmap=plt.cm.bwr,vmin=-100,vmax=100)
cbar = map.colorbar(cs,location='bottom',size="5%",pad="2%")
cbar.ax.tick_params(labelsize=15)
plt.axis('off')
#plt.savefig('maiiizisam.jpg')
plt.savefig('maiiizisam_et100.jpg',dpi=300,bbox_inches='tight')

plt.show()
#fig = plt.figure(figsize=(20,15))

