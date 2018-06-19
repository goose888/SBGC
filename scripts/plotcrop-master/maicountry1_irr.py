from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import glob
import numpy.ma as ma
from scipy.interpolate import griddata
import scipy.stats

iizumi=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/iizumi/iizumi.2013JAN29.maize.1982-2006.30min.nc4','r')
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
clmtropf = clm.variables['yield'][99,:,:]
#clmtropf=N.average(clmtrop,axis=0)

clm1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/maizetemp_historical_co2_rf_fert_0.5x0.5.nc','r')
clmtempf = clm1.variables['yield'][99,:,:]
#clmtempf=N.average(clmtemp,axis=0)

clm2=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/maizetrop_historical_co2_irrig_fert_0.5x0.5.nc','r')
clmtropfi = clm2.variables['yield'][99,:,:]


clm3=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/maizetemp_historical_co2_irrig_fert_0.5x0.5.nc','r')
clmtempfi = clm3.variables['yield'][99,:,:]


#print clmtropf
#print clmtropf.shape
clmtropf=N.flipud(clmtropf)
clmtempf=N.flipud(clmtempf)

clmtropfi=N.flipud(clmtropfi)
clmtempfi=N.flipud(clmtempfi)

clmtropf= ma.masked_where(maitrop<=0,clmtropf)
clmtempf= ma.masked_where(maitemp<=0,clmtempf)
clmtropf=ma.filled(clmtropf, fill_value=0.)
clmtempf=ma.filled(clmtempf, fill_value=0.)

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
yield_clmtf  = ma.masked_where(maizetor<=0,yield_clmtf )
yield_clmtf=ma.filled(yield_clmtf, fill_value=0.)

yield_clmtfi=clmtropfi+clmtempfi
yield_clmtfi = ma.masked_where(yield_clmtfi<=0,yield_clmtfi)
yield_clmtfi = ma.masked_where(maizetoi<=0,yield_clmtfi)
yield_clmtfi=ma.filled(yield_clmtfi, fill_value=0.)






area=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/gridareahalf.nc','r')
gridarea = area.variables['cell_area'][:,:]
gridlon = area.variables['lon'][:]

gridarea,gridlon = shiftgrid(180.5,gridarea,gridlon,start=False)

area1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/maizegridarea.nc','r')
gridaream3 = area1.variables['cell_area'][:,:]
gridlonm3 = area1.variables['lon'][:]
gridlatm3 = area1.variables['lat'][:]

#print gridarea

nclu=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/maize_AreaYieldProduction.nc','r')
#print(nclu)
ncvar_maize = nclu.variables['maizeData'][:]
latnc = nclu.variables['latitude'][:]
znc = nclu.variables['level'][:]
lonnc = nclu.variables['longitude'][:]
timenc = nclu.variables['time'][:]
#lon,lat = N.meshgrid(lonnc,latnc)
lat_new=N.flipud(latnc)
lon,lat = N.meshgrid(lonnc,lat_new)
ncvar_maizea = nclu.variables['maizeData'][0,1,:,:]
ncvar_maizearea = nclu.variables['maizeData'][0,4,:,:]

ncvar_maizea=N.flipud(ncvar_maizea)

#print lat_new
maize_new = N.flipud(ncvar_maize)
ncvar_maizearea = N.flipud(ncvar_maizearea)


yieldfi= N.zeros((1, 360, 720))
yieldfai= N.zeros((1, 360, 720))
yieldf= N.zeros((1, 360, 720))
yieldf2= N.zeros((1, 360, 720))
years = range(2000, 2001)
for i, year in enumerate(years):
    base = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_his/cruhis/output/cruhis.bgp-yearly_crop_{0}.nc".format(year), mode='r')
    base2 = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_his/cesmhis/output/cesmhis.bgp-yearly_crop_{0}.nc".format(year), mode='r')
    lona1 = base.variables["lon"][:]
    lata1 = base.variables["lat"][:]
   
    yield1 = base.variables["yield"][0,:,:]
    yieldf[i, :, :] = yield1
    yield2 = base2.variables["yield"][0,:,:]
    yieldf2[i, :, :] = yield2
    basei = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_his/cruhis_irr/output/cruhisirr.bgp-yearly_crop_{0}.nc".format(year), mode='r')
    base2i = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_his/cesmhis_irr/output/cesmhisirr.bgp-yearly_crop_{0}.nc".format(year), mode='r')

    yield1i = basei.variables["yield"][0,:,:]
    yieldfi[i, :, :] = yield1i
    yield2i = base2i.variables["yield"][0,:,:]
    yieldfai[i, :, :] = yield2i


yielda=N.average(yieldf,axis=0)
yielda2=N.average(yieldf2,axis=0)

yield_new,lona11 = shiftgrid(180.5,yielda,lona1,start=False)
yield_new1,lona11 = shiftgrid(180.5,yielda2,lona1,start=False)

yieldai=N.average(yieldfi,axis=0)
yield_newi,lona11 = shiftgrid(180.5,yieldai,lona1,start=False)
yieldbi=N.average(yieldfai,axis=0)
yield_new1i,lona11 = shiftgrid(180.5,yieldbi,lona1,start=False)


yield_new[N.isnan(yield_new)] = -9999
yield_new = ma.masked_where(yield_new<=0,yield_new)
yield_new = ma.masked_where(maizetor<=0,yield_new)
yield_new=ma.filled(yield_new, fill_value=0.)

yield_new1[N.isnan(yield_new1)] = -9999
yield_new1 = ma.masked_where(yield_new1<=0,yield_new1)
yield_new1 = ma.masked_where(maizetor<=0,yield_new1)
yield_new1=ma.filled(yield_new1, fill_value=0.)

yield_newi[N.isnan(yield_newi)] = -9999
yield_newi = ma.masked_where(yield_newi<=0,yield_newi)
yield_newi = ma.masked_where(maizetoi<=0,yield_newi)
yield_newi=ma.filled(yield_newi, fill_value=0.)


yield_new1i[N.isnan(yield_new1i)] = -9999
yield_new1i = ma.masked_where(yield_new1i<=0,yield_new1i)
yield_new1i = ma.masked_where(maizetoi<=0,yield_new1i)
yield_new1i=ma.filled(yield_new1i, fill_value=0.)




lon2,lat2 = N.meshgrid(lona11,lata1)

#yield_smooth = interp(maize_new, lonnc,lat_new, lons_sub,lats_sub , order=1)
m3coun = interp(coun,lona11,lata1,lon,lat,order=1)
#yield_newm3= interp(yield_new,lona11,lata1,lon,lat,order=1)
#yield_new2m3= interp(yield_new2,lona11,lata1,lon,lat,order=1)

#http://matplotlib.org/basemap/users/mapsetup.html
iarea = ma.masked_where(iarea<=0,iarea)
iarea=ma.filled(iarea, fill_value=0.)

iyield = ma.masked_where(iyield<=0,iyield)
#iyield = ma.masked_where(iarea<=0,iyield)
iyield=ma.filled(iyield, fill_value=0.)



#yield_clmtf=clmtropf+clmtempf
#yield_clmtfm3= interp(yield_clmtf,lona11,lata1,lon,lat,order=1)
#yield_clmtf = ma.masked_where(yield_clmtf<=0,yield_clmtf)
#yield_clmtf = ma.masked_where(iyield<=0,yield_clmtf)
#yield_clmtf=ma.filled(yield_clmtf, fill_value=0.)


ncvar_maizea[N.isnan(ncvar_maizea)] = -9999
ncvar_maizearea[N.isnan(ncvar_maizearea)] = -9999
ncvar_maizea = ma.masked_where(ncvar_maizea<=0,ncvar_maizea)
ncvar_maizearea = ma.masked_where(ncvar_maizearea<=0,ncvar_maizearea)
ncvar_maizearea=ma.filled(ncvar_maizearea, fill_value=0.)
ncvar_maizea=ma.filled(ncvar_maizea, fill_value=0.)

m3yy = interp(ncvar_maizea,lonnc,lat_new,lon2,lat2,order=1)


def annualyield(couna,counb):

	yieldagf=0.
	yieldg=0.
	harea=0.
	a=0
	yieldagfi=0.
	yieldgi=0.
        m3yieldg=0.
        m3yieldag=0.

	yieldagfc=0.
	yieldgc=0.
	yieldgi1=0.
	yieldagfi1=0.

#USA 11501~11550
	for xx in range(0,360):
	    for yy in range(0,720):
	        if coun[xx,yy] >=couna and  coun[xx,yy] <=counb:
	            yieldg=iyield[xx,yy]*maizeto[xx,yy]*gridarea[xx,yy]+yieldg
                    m3yieldg=m3yy[xx,yy]*maizeto[xx,yy]*gridarea[xx,yy]+m3yieldg
	            harea=maizeto[xx,yy]*gridarea[xx,yy]+ harea
                    yieldgi=(yield_new[xx,yy]*gridarea[xx,yy]*maizetor[xx,yy])+(yield_newi[xx,yy]*gridarea[xx,yy]*maizetoi[xx,yy])+yieldgi
                    yieldgc=(yield_clmtf[xx,yy]*gridarea[xx,yy]*maizetor[xx,yy])+(yield_clmtfi[xx,yy]*gridarea[xx,yy]*maizetoi[xx,yy])+yieldgc
                    yieldgi1=(yield_new1[xx,yy]*gridarea[xx,yy]*maizetor[xx,yy])+(yield_new1i[xx,yy]*gridarea[xx,yy]*maizetoi[xx,yy])+yieldgi1
                
	

	yieldagf=yieldg/harea
	yieldagfi=yieldgi/harea
	yieldagfc=yieldgc/harea
	yieldagfi1=yieldgi1/harea
        m3yieldag=m3yieldg/harea
	return "Maize","iizumi harvested area",harea, "ha","iizumi yield",yieldagf,"t/ha","production",yieldg,"tonnes","ISAM-CRU yield",yieldagfi,"t/ha","production",yieldgi,"tonnes","CLM yield",yieldagfc,"t/ha","production",yieldgc,"tonnes","ISAM-CESM yield",yieldagfi1,"t/ha","production",yieldgi1,"tonnes","M3 harvested area","no", "ha","M3 yield",m3yieldag,"t/ha","production",m3yieldg,"tonnes"

name=["globe","italy","chile","spain","germany","france","us","canada","argentina","china","hungary","thailand","ukraine","southafrica","indonesia","vietnam","brazil","mexico","india","philippines","romania","nigeria"]
range1=[10100,43400,20400,46800,42500,42200,11501,10201,20101,41501,42700,47500,47800,33900,50300,48200,20301,11101,42901,50700,46000,33400]
range2=[50700,43400,20400,46800,42500,42200,11550,10212,20124,41529,42700,47500,47800,33900,50300,48200,20327,11132,42929,50700,46000,33400]

a=0
x=22
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
m3area=N.zeros(x)
for a, name1 in enumerate(name):

	reu=annualyield(range1[a],range2[a])
	print reu
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
#        m3area[a]=reu[29]

# data to plot
n_groups = 22

means_frank = zumiy
means_guido = isamy
means_final = clmy
means_fao = (43244*0.0001,95277*0.0001,94120*0.0001,92157*0.0001,92119*0.0001,90768*0.0001,85910*0.0001,62844*0.0001,54329*0.0001 ,45992*0.0001,41790*0.0001,36715*0.0001,30091*0.0001,28492*0.0001,27649*0.0001,27471*0.0001 ,27447*0.0001,24620*0.0001,18216*0.0001,17970*0.0001,16061*0.0001,13001*0.0001)

means_m3 = m3y
means_cesm=isamy1
means_all=(means_fao+zumiy+m3y)/3


means_frank1 = zumip
means_guido1 = isamp
means_final1 = clmp
pro_m3 = m3p
pro_cesm=isamp1

# create plot
#fig, ax = plt.subplots(211)
fig = plt.figure(figsize=(33,10))


ax = fig.add_subplot(111)
index = N.arange(n_groups)
bar_width = 0.21
opacity = 0.8
#rects0 = plt.bar(index, means_fao, bar_width,
#                 alpha=opacity,
#                 color='k',
#                 label='FAO')
#rects4 = plt.bar(index+bar_width, means_m3, bar_width,
#                 alpha=opacity,
#                 color='c',
#                 label='M3')
#rects1 = plt.bar(index+bar_width*2, means_frank, bar_width,
#                 alpha=opacity,
#                 color='r',
#                 label='Iizumi')
rects1 = plt.bar(index, means_fao, bar_width,
                 alpha=opacity,
                 color='r',
                 label='FAO')
 
rects2 = plt.bar(index + bar_width*1, means_guido, bar_width,
                 alpha=opacity,
                 color='g',
                 label='ISAM-NCEP')
rects0 = plt.bar(index + bar_width*2, means_cesm, bar_width,
                 alpha=opacity,
                 color='y',
                 label='ISAM-CESM')
rects3 = plt.bar(index + bar_width*3, means_final, bar_width,
                 alpha=opacity,
                 color='b',
                 label='CLM')
#plt.ylim(0,4)
 
plt.xlabel('Country',fontsize=40)
plt.ylabel('Maize yield (t/ha)',fontsize=40)
#plt.title('2000',fontsize=18)
plt.xticks(index + bar_width+0.2, ("Globe","Italy","Chile","Spain","Germany","France","USA","Canada","Argentina","China","Hungary","Thailand","Ukraine","South Africa","Indonesia","Vietnam","Brazil","Mexico","India","Philippines","Romania","Nigeria"),rotation='vertical')
#plt.xticks(index + bar_width+0.2, ("Nigeria","Thailand","Japan","South Korea","France","South Africa","Vietnam","Russia","Zimbabwe","Iran","Uganda"))
leg=plt.legend(fontsize=40)
leg.get_frame().set_alpha(0.5)
plt.tick_params(axis='both',labelsize=40)


plt.tight_layout()
plt.savefig('mai2000edison_irr1.png')
#plt.show()

fig, ax = plt.subplots(figsize=(6,6))
colors = (1,1,0)
colorsr = (0,0,1)
color=(0,1,0)
ax.plot([0,12],[0,12], 'k--',label='1:1')
ax.scatter(means_fao, means_cesm,s=124, c=colors,alpha=0.5,label='ISAM-CESM')
ax.scatter(means_fao, means_final,s=124, c=colorsr,alpha=0.5,label='CLM')
ax.scatter(means_fao, means_guido,s=124, c=color,alpha=0.5,label='ISAM-NCEP')
ccp1=scipy.stats.pearsonr(means_cesm,means_fao)
ccp2=scipy.stats.pearsonr(means_final,means_fao)
ccp3=scipy.stats.pearsonr(means_guido,means_fao)
leg=ax.legend(['1:1', 'ISAM-CESM {:04.2f}'.format(ccp1[0]),'CLM {:04.2f}'.format(ccp2[0]),'ISAM-NCEP {:04.2f}'.format(ccp3[0])],fontsize=12)


plt.xlim(0, 12)
plt.ylim(0, 12)
#leg=ax.legend(fontsize=15)
leg.get_frame().set_alpha(0.5)

plt.tick_params(axis='both',labelsize=16)
plt.xlabel('FAO (t/ha)',fontsize=18)
plt.ylabel('Model (t/ha)',fontsize=18)
plt.savefig('scatter_mai_count.png',bbox_inches='tight')



