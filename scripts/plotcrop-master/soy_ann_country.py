from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib.pyplot as plt
import glob
import numpy.ma as ma
from scipy.interpolate import griddata
import scipy.stats
from matplotlib.dates import DateFormatter
import datetime


country=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/Ctry_halfdeg.nc','r')
#print iizumi
coun = country.variables['MASK_Country'][:,:]
area=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/gridareahalf.nc','r')
gridarea = area.variables['cell_area'][:,:]
gridlon = area.variables['lon'][:]

gridarea,gridlon = shiftgrid(180.5,gridarea,gridlon,start=False)



def annualyield(year,couna,counb):
    aa=year-1982
    bb=year-1901
    iizumi=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/iizumi/iizumi.2013JAN29.soybean.1982-2006.30min.nc4','r')

    iyield = iizumi.variables['yield50'][aa,:,:]
    iarea =iizumi.variables['area'][aa,:,:]
    la=iizumi.variables['lat'][:]
    lo=iizumi.variables['lon'][:]
    iyield=N.flipud(iyield)
    iarea=N.flipud(iarea)
    
    region1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/HistoricalGLM_crop_150901.nc','r')
    maitrop = region1.variables['soy_trop'][bb,:,:]
    maitemp = region1.variables['soy_temp'][bb,:,:]
    maizeto = maitrop+maitemp


    clm=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/soytrop_historical_co2_rf_nofert_0.5x0.5.nc','r')
    clmtropf = clm.variables['yield'][bb,:,:]


    clm1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/soytemp_historical_co2_rf_nofert_0.5x0.5.nc','r')
    clmtempf = clm1.variables['yield'][bb,:]

    clmtropf=N.flipud(clmtropf)
    clmtempf=N.flipud(clmtempf)

    clmtropf= ma.masked_where(maitrop<=0,clmtropf)
    clmtempf= ma.masked_where(maitemp<=0,clmtempf)
    clmtropf=ma.filled(clmtropf, fill_value=0.)
    clmtempf=ma.filled(clmtempf, fill_value=0.)
    yieldfa= N.zeros((1, 360, 720))
    yieldf= N.zeros((1, 360, 720))
    year1=year+1
    years = range(year,year1)
    for i, year in enumerate(years):
        base = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_his/cruhis/output/cruhis.bgp-yearly_crop_{0}.nc".format(year), mode='r')
        base2 = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_his/cesmhis/output/cesmhis.bgp-yearly_crop_{0}.nc".format(year), mode='r')

        lona1 = base.variables["lon"][:]
        lata1 = base.variables["lat"][:]
        yield1 = base.variables["yield"][1,:,:]
        yieldf[i, :, :] = yield1
        yield2 = base2.variables["yield"][1,:,:]
        yieldfa[i, :, :] = yield2

    yielda=N.average(yieldf,axis=0)
    yield_new,lona11 = shiftgrid(180.5,yielda,lona1,start=False)
    yieldb=N.average(yieldfa,axis=0)
    yield_new1,lona11 = shiftgrid(180.5,yieldb,lona1,start=False)


    
    iarea = ma.masked_where(iarea<=0,iarea)
    iarea=ma.filled(iarea, fill_value=0.)

    iyield = ma.masked_where(iyield<=0,iyield)
    iyield = ma.masked_where(iarea<=0,iyield)
    iyield=ma.filled(iyield, fill_value=0.)


    yield_new[N.isnan(yield_new)] = -9999
    yield_new = ma.masked_where(yield_new<=0,yield_new)
    yield_new = ma.masked_where(iyield<=0,yield_new)
    yield_new=ma.filled(yield_new, fill_value=0.)

    yield_new1[N.isnan(yield_new1)] = -9999
    yield_new1 = ma.masked_where(yield_new1<=0,yield_new1)
    yield_new1 = ma.masked_where(iyield<=0,yield_new1)
    yield_new1=ma.filled(yield_new1, fill_value=0.)



    yield_clmtf=clmtropf+clmtempf
    yield_clmtf = ma.masked_where(yield_clmtf<=0,yield_clmtf)
    yield_clmtf = ma.masked_where(iyield<=0,yield_clmtf)
    yield_clmtf=ma.filled(yield_clmtf, fill_value=0.)
    
    yieldagf=0.
    yieldg=0.
    harea=0.
    a=0
    yieldagfi=0.
    yieldgi=0.
    yieldgi1=0.
    yieldagfi1=0.

    yieldagfc=0.
    yieldgc=0.

    #USA 11501~11550
    for xx in range(0,360):
        for yy in range(0,720):
            if coun[xx,yy] >=couna and  coun[xx,yy] <=counb:
                yieldg=iyield[xx,yy]*iarea[xx,yy]/100*gridarea[xx,yy]/10000+yieldg
                harea=iarea[xx,yy]/100*gridarea[xx,yy]/10000+ harea
                yieldgi=yield_new[xx,yy]*iarea[xx,yy]/100*gridarea[xx,yy]/10000+yieldgi
                yieldgc=yield_clmtf[xx,yy]*iarea[xx,yy]/100*gridarea[xx,yy]/10000+yieldgc
                yieldgi1=yield_new1[xx,yy]*iarea[xx,yy]/100*gridarea[xx,yy]/10000+yieldgi1

                a=a+1

    yieldagf=yieldg/harea
    yieldagfi=yieldgi/harea
    yieldagfc=yieldgc/harea
    yieldagfi1=yieldgi1/harea

    return "Soybean","iizumi harvested area",harea, "ha","iizumi yield",yieldagf,"t/ha","production",yieldg,"tonnes","ISAM-CRU yield",yieldagfi,"t/ha","production",yieldgi,"tonnes","CLM yield",yieldagfc,"t/ha","production",yieldgc,"tonnes","ISAM-CESM yield",yieldagfi1,"t/ha","production",yieldgi1,"tonnes"
    #return harea

def runmean(x,input):
        import pandas as pd
        #mean_zumiy1 = pd.rolling_mean(zumiy, window=5).shift(-2)
        meanout=pd.rolling_mean(input, window=5, center=True)
        mean_zumiy1=pd.rolling_mean(input, window=3, center=True)
        #print mean_zumiy1
        #print meanout
        meanout[1]=mean_zumiy1[1]
        meanout[x-2]=mean_zumiy1[x-2]
        meanout1=input-meanout
        return meanout1

#illzmui only 1983~2005
a1=1983
a2=2005
x=a2-a1
toarea= N.zeros(x)
zumiy= N.zeros(x)
clmy= N.zeros(x)
isamy= N.zeros(x)
zumip= N.zeros(x)
clmp= N.zeros(x)
isamp= N.zeros(x)
isamy1= N.zeros(x)
isamp1= N.zeros(x)

#print toarea

name=["Global","USA","Brazil","Argentina","China","India","Paraguay","Canada","Uruguay","Italy","Bolivia"]
range1=[10100,11501,20301,20101,41501,42901,20900,10201,21300,43400,20200]
range2=[50700,11550,20327,20124,41529,42929,20900,10212,21300,43400,20200]

for i, name1 in enumerate(name):
	a=0
	toarea= N.zeros(x)
	zumiy= N.zeros(x)
	clmy= N.zeros(x)
	isamy= N.zeros(x)
	zumip= N.zeros(x)
	clmp= N.zeros(x)
	isamp= N.zeros(x)
	isamy1= N.zeros(x)
	isamp1= N.zeros(x)

	for num in range(a1,a2):
    
#    reu=annualyield(num,11501,11550)#usa
#    reu=annualyield(num,10100,50700)#globe
	    reu=annualyield(num,range1[i],range2[i])

#    print reu
	    toarea[a]=reu[2]
	    zumiy[a]=reu[5]
	    isamy[a]=reu[11]
	    clmy[a]=reu[17]
	    zumip[a]=reu[8]
	    isamp[a]=reu[14]
	    clmp[a]=reu[20]
	    isamy1[a]=reu[23]
	    isamp1[a]=reu[26]
	    a=a+1

	azumiy=zumiy-N.average(zumiy)
	aisamy=isamy-N.average(isamy)
	aclmy=clmy-N.average(clmy)
	azumip=zumip-N.average(zumip)
	aisamp=isamp-N.average(isamp)
	aclmp=clmp-N.average(clmp)
	aisamy1=isamy1-N.average(isamy1)
	aisamp1=isamp1-N.average(isamp1)


	mean_zumiy=runmean(x,zumiy)
	mean_isamy=runmean(x,isamy)
	mean_clmy=runmean(x,clmy)
	mean_zumip=runmean(x,zumip)
	mean_isamp=runmean(x,isamp)
	mean_clmp=runmean(x,clmp)
	mean_isamy1=runmean(x,isamy1)
	mean_isamp1=runmean(x,isamp1)
#print mean_zumiy

	fig = plt.figure(figsize=(26,20))


	ax = fig.add_subplot(321)
	xx=range(a1,a2)

#plt.ylim((0,12))

	xdates = [datetime.datetime.strptime(str(int(date)),'%Y') for date in xx]
#plt.xticks(xdates, xdates)

	ax.plot_date(xdates,zumiy,"ro-",label="Iizumi",linewidth=2)
	ax.plot_date(xdates,isamy,"go-",label="ISAM-NCEP")
	ax.plot_date(xdates,isamy1,"yo-",label="ISAM-CESM")
	ax.plot_date(xdates,clmy,"bo-",label="CLM")
	ax.xaxis.set_major_formatter(DateFormatter('%Y'))
	ccp=scipy.stats.pearsonr(isamy,zumiy)
	ccp1=scipy.stats.pearsonr(isamy1,zumiy)
	ccp2=scipy.stats.pearsonr(clmy,zumiy)

	leg=plt.legend(['Iizumi', 'ISAM-NCEP {:05.3f}'.format(ccp[0]),'ISAM-CESM {:05.3f}'.format(ccp1[0]),'CLM {:05.3f}'.format(ccp2[0])])
#leg = plt.legend(loc=2,fancybox=True, fontsize=14)
	leg.get_frame().set_alpha(0.5)

	plt.xlabel("Year (AD)",fontsize=18)
	plt.ylabel("Soybean yield (t/ha)",fontsize=18)
	plt.title("Yield {0}".format(name1),fontsize=18)
	plt.tick_params(axis='both',labelsize=15)


	ax = fig.add_subplot(322)
	ax.plot_date(xdates,zumip,"ro-",label="Iizumi",linewidth=2)
	ax.plot_date(xdates,isamp,"go-",label="ISAM-NCEP")
	ax.plot_date(xdates,isamp1,"yo-",label="ISAM-CESM")
	ax.plot_date(xdates,clmp,"bo-",label="CLM")
	ax.xaxis.set_major_formatter(DateFormatter('%Y'))

	ccp=scipy.stats.pearsonr(isamp,zumip)
	ccp1=scipy.stats.pearsonr(isamp1,zumip)
	ccp2=scipy.stats.pearsonr(clmp,zumip)

	leg=plt.legend(['Iizumi', 'ISAM-NCEP {:05.3f}'.format(ccp[0]),'ISAM-CESM {:05.3f}'.format(ccp1[0]),'CLM {:05.3f}'.format(ccp2[0])])
#leg = plt.legend(loc=2,fancybox=True, fontsize=14)
	leg.get_frame().set_alpha(0.5)

	plt.title("Production {0}".format(name1),fontsize=18)
	plt.xlabel("Year (AD)",fontsize=18)
	plt.ylabel("Soybean production (tonnes)",fontsize=18)
	plt.tick_params(axis='both',labelsize=15)



	ax = fig.add_subplot(323)
	ax.plot_date(xdates,azumiy,"ro-",label="Iizumi",linewidth=2)
	ax.plot_date(xdates,aisamy,"go-",label="ISAM-NCEP")
	ax.plot_date(xdates,aisamy1,"yo-",label="ISAM-CESM")
	ax.plot_date(xdates,aclmy,"bo-",label="CLM")
	ax.xaxis.set_major_formatter(DateFormatter('%Y'))

	ccp=scipy.stats.pearsonr(aisamy,azumiy)
	ccp1=scipy.stats.pearsonr(aisamy1,azumiy)
	ccp2=scipy.stats.pearsonr(aclmy,azumiy)

	leg=plt.legend(['Iizumi', 'ISAM-NCEP {:05.3f}'.format(ccp[0]),'ISAM-CESM {:05.3f}'.format(ccp1[0]),'CLM {:05.3f}'.format(ccp2[0])])
#leg = plt.legend(loc=2,fancybox=True, fontsize=14)
	leg.get_frame().set_alpha(0.5)


	plt.title("Anormaly by subtracting average {0}".format(name1),fontsize=18)
	plt.xlabel("Year (AD)",fontsize=18)
	plt.ylabel("Soybean yield (t/ha)",fontsize=18)
	plt.tick_params(axis='both',labelsize=15)



	ax = fig.add_subplot(324)
	ax.plot_date(xdates,azumip,"ro-",label="Iizumi",linewidth=2)
	ax.plot_date(xdates,aisamp,"go-",label="ISAM-NCEP")
	ax.plot_date(xdates,aisamp1,"yo-",label="ISAM-CESM")
	ax.plot_date(xdates,aclmp,"bo-",label="CLM")
	ax.xaxis.set_major_formatter(DateFormatter('%Y'))

	ccp=scipy.stats.pearsonr(aisamp,azumip)
	ccp1=scipy.stats.pearsonr(aisamp1,azumip)
	ccp2=scipy.stats.pearsonr(aclmp,azumip)

	leg=plt.legend(['Iizumi', 'ISAM-NCEP {:05.3f}'.format(ccp[0]),'ISAM-CESM {:05.3f}'.format(ccp1[0]),'CLM {:05.3f}'.format(ccp2[0])])
#leg = plt.legend(loc=2,fancybox=True, fontsize=14)
	leg.get_frame().set_alpha(0.5)


	plt.title("Anormaly by subtracting average {0}".format(name1),fontsize=18)
	plt.xlabel("Year (AD)",fontsize=18)
	plt.ylabel("Soybean production (tonnes)",fontsize=18)
	plt.tick_params(axis='both',labelsize=15)



	ax = fig.add_subplot(325)
	ax.plot_date(xdates,mean_zumiy,"ro-",label="Iizumi",linewidth=2)
	ax.plot_date(xdates,mean_isamy,"go-",label="ISAM-NCEP")
	ax.plot_date(xdates,mean_isamy1,"yo-",label="ISAM-CESM")
	ax.plot_date(xdates,mean_clmy,"bo-",label="CLM")
	ax.xaxis.set_major_formatter(DateFormatter('%Y'))

	mean1_isamy=N.zeros([x-2])
	mean1_zumiy=N.zeros([x-2])
	mean1_isamy1=N.zeros([x-2])
	mean1_clmy=N.zeros([x-2])


	for i in range(1,x-2):
		mean1_isamy[i]=mean_isamy[i]
		mean1_zumiy[i]=mean_zumiy[i]
		mean1_isamy1[i]=mean_isamy1[i]
		mean1_clmy[i]=mean_clmy[i]

	ccp=scipy.stats.pearsonr(mean1_isamy,mean1_zumiy)
	ccp1=scipy.stats.pearsonr(mean1_isamy1,mean1_zumiy)
	ccp2=scipy.stats.pearsonr(mean1_clmy,mean1_zumiy)

	leg=plt.legend(['Iizumi', 'ISAM-NCEP {:05.3f}'.format(ccp[0]),'ISAM-CESM {:05.3f}'.format(ccp1[0]),'CLM {:05.3f}'.format(ccp2[0])])
#leg = plt.legend(loc=2,fancybox=True, fontsize=14)
	leg.get_frame().set_alpha(0.5)


	plt.title("Anormaly by subtracting a moving mean average {0}".format(name1),fontsize=18)
	plt.xlabel("Year (AD)",fontsize=18)
	plt.ylabel("Soybean yield (t/ha)",fontsize=18)
	plt.tick_params(axis='both',labelsize=15)



	ax = fig.add_subplot(326)
	ax.plot_date(xdates,mean_zumip,"ro-",label="Iizumi",linewidth=2)
	ax.plot_date(xdates,mean_isamp,"go-",label="ISAM-NCEP")
	ax.plot_date(xdates,mean_isamp1,"yo-",label="ISAM-CESM")
	ax.plot_date(xdates,mean_clmp,"bo-",label="CLM")
	ax.xaxis.set_major_formatter(DateFormatter('%Y'))

	mean1_isamp=N.zeros([x-2])
	mean1_zumip=N.zeros([x-2])
	mean1_isamp1=N.zeros([x-2])
	mean1_clmp=N.zeros([x-2])

	for i in range(1,x-2):
		mean1_isamp[i]=mean_isamp[i]
		mean1_zumip[i]=mean_zumip[i]
		mean1_isamp1[i]=mean_isamp1[i]
		mean1_clmp[i]=mean_clmp[i]

	ccp=scipy.stats.pearsonr(mean1_isamp,mean1_zumip)
	ccp1=scipy.stats.pearsonr(mean1_isamp1,mean1_zumip)
	ccp2=scipy.stats.pearsonr(mean1_clmp,mean1_zumip)

	leg=plt.legend(['Iizumi', 'ISAM-NCRP {:05.3f}'.format(ccp[0]),'ISAM-CESM {:05.3f}'.format(ccp1[0]),'CLM {:05.3f}'.format(ccp2[0])])
#leg = plt.legend(loc=2,fancybox=True, fontsize=14)
	leg.get_frame().set_alpha(0.5)


	plt.title("Anormaly by subtracting a moving mean average {0}".format(name1),fontsize=18)
	plt.xlabel("Year (AD)",fontsize=18)
	plt.ylabel("Soybean production (tonnes)",fontsize=18)
	plt.tick_params(axis='both',labelsize=15)


	plt.savefig('soy_ag_{0}.png'.format(name1),dpi=300,bbox_inches='tight')
#	plt.show()


