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
#area=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/gridareahalf.nc','r')
#gridarea = area.variables['cell_area'][:,:]
#gridlon = area.variables['lon'][:]

#gridarea,gridlon = shiftgrid(180.5,gridarea,gridlon,start=False)



def annualyield(year,couna,counb):
    bb=year-1901
    
    region1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/HistoricalGLM_crop_150901.nc','r')
    maitrop = region1.variables['soy_trop'][bb,:,:]
    maitemp = region1.variables['soy_temp'][bb,:,:]
    maitropi=region1.variables['soy_trop_irrig'][bb,:,:]
    maitempi=region1.variables['soy_temp_irrig'][bb,:,:]
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


    clm=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/soytrop_historical_co2_rf_fert_0.5x0.5.nc','r')
    clmtropf = clm.variables['yield'][bb,:,:]


    clm1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/soytemp_historical_co2_rf_fert_0.5x0.5.nc','r')
    clmtempf = clm1.variables['yield'][bb,:,:]

    clm2=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/soytrop_historical_co2_irrig_fert_0.5x0.5.nc','r')
    clmtropfi = clm2.variables['yield'][bb,:,:]


    clm3=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/soytemp_historical_co2_irrig_fert_0.5x0.5.nc','r')
    clmtempfi = clm3.variables['yield'][bb,:,:]


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
    yield_clmtf  = ma.masked_where(maizetor<=0,yield_clmtf )
    yield_clmtf=ma.filled(yield_clmtf, fill_value=0.)

    yield_clmtfi=clmtropfi+clmtempfi
    yield_clmtfi = ma.masked_where(yield_clmtfi<=0,yield_clmtfi)
    yield_clmtfi = ma.masked_where(maizetoi<=0,yield_clmtfi)
    yield_clmtfi=ma.filled(yield_clmtfi, fill_value=0.)


    clmn=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/soytrop_historical_co2_rf_nofert_0.5x0.5.nc','r')
    clmtropfn = clmn.variables['yield'][bb,:,:]


    clm1n=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/soytemp_historical_co2_rf_nofert_0.5x0.5.nc','r')
    clmtempfn = clm1n.variables['yield'][bb,:,:]

    clm2n=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/soytrop_historical_co2_irrig_nofert_0.5x0.5.nc','r')
    clmtropfin = clm2n.variables['yield'][bb,:,:]


    clm3n=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45historical/soytemp_historical_co2_irrig_nofert_0.5x0.5.nc','r')
    clmtempfin = clm3n.variables['yield'][bb,:,:]


    clmtropfn=N.flipud(clmtropfn)
    clmtempfn=N.flipud(clmtempfn)
    clmtropfin=N.flipud(clmtropfin)
    clmtempfin=N.flipud(clmtempfin)

    clmtropfn= ma.masked_where(maitrop<=0,clmtropfn)
    clmtempfn= ma.masked_where(maitemp<=0,clmtempfn)
    clmtropfn=ma.filled(clmtropfn, fill_value=0.)
    clmtempfn=ma.filled(clmtempfn, fill_value=0.)

    clmtropfin= ma.masked_where(maitropi<=0,clmtropfin)
    clmtempfin= ma.masked_where(maitempi<=0,clmtempfin)
    clmtropfin=ma.filled(clmtropfin, fill_value=0.)
    clmtempfin=ma.filled(clmtempfin, fill_value=0.)

    yield_clmtfn=clmtropfn+clmtempfn
    yield_clmtfn = ma.masked_where(yield_clmtfn<=0,yield_clmtfn)
    yield_clmtfn  = ma.masked_where(maizetor<=0,yield_clmtfn )
    yield_clmtfn=ma.filled(yield_clmtfn, fill_value=0.)

    yield_clmtfin=clmtropfin+clmtempfin
    yield_clmtfin = ma.masked_where(yield_clmtfin<=0,yield_clmtfin)
    yield_clmtfin = ma.masked_where(maizetoi<=0,yield_clmtfin)
    yield_clmtfin=ma.filled(yield_clmtfin, fill_value=0.)


    yieldfa= N.zeros((1, 360, 720))
    yieldf= N.zeros((1, 360, 720))
    yieldfi= N.zeros((1, 360, 720))
    yieldfai= N.zeros((1, 360, 720))
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

        basei = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_his/cruhis_irr/output/cruhisirr.bgp-yearly_crop_{0}.nc".format(year), mode='r')
        base2i = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_his/cesmhis_irr/output/cesmhisirr.bgp-yearly_crop_{0}.nc".format(year), mode='r')

        yield1i = basei.variables["yield"][1,:,:]
        yieldfi[i, :, :] = yield1i
        yield2i = base2i.variables["yield"][1,:,:]
        yieldfai[i, :, :] = yield2i



    yielda=N.average(yieldf,axis=0)
    yield_new,lona11 = shiftgrid(180.5,yielda,lona1,start=False)
    yieldb=N.average(yieldfa,axis=0)
    yield_new1,lona11 = shiftgrid(180.5,yieldb,lona1,start=False)

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
    yieldagfcn=0.
    yieldgcn=0.
    #USA 11501~11550
    for xx in range(0,360):
        for yy in range(0,720):
            if coun[xx,yy] >=couna and  coun[xx,yy] <=counb:
                yieldg=0+yieldg
                harea=maizeto[xx,yy]*gridarea[xx,yy]+ harea
                yieldgi=(yield_new[xx,yy]*gridarea[xx,yy]*maizetor[xx,yy])+(yield_newi[xx,yy]*gridarea[xx,yy]*maizetoi[xx,yy])+yieldgi
                yieldgc=(yield_clmtf[xx,yy]*gridarea[xx,yy]*maizetor[xx,yy])+(yield_clmtfi[xx,yy]*gridarea[xx,yy]*maizetoi[xx,yy])+yieldgc
                yieldgi1=(yield_new1[xx,yy]*gridarea[xx,yy]*maizetor[xx,yy])+(yield_new1i[xx,yy]*gridarea[xx,yy]*maizetoi[xx,yy])+yieldgi1
                yieldgcn=(yield_clmtfn[xx,yy]*gridarea[xx,yy]*maizetor[xx,yy])+(yield_clmtfin[xx,yy]*gridarea[xx,yy]*maizetoi[xx,yy])+yieldgcn

                a=a+1

    yieldagf=yieldg/harea
    yieldagfi=yieldgi/harea
    yieldagfc=yieldgc/harea
    yieldagfcn=yieldgcn/harea
    yieldagfi1=yieldgi1/harea

    return "Soy ","iizumi harvested area",harea, "ha","iizumi yield",yieldagf,"t/ha","production",yieldg,"tonnes","ISAM-CRU yield",yieldagfi,"t/ha","production",yieldgi,"tonnes","CLM-N yield",yieldagfc,"t/ha","production",yieldgc,"tonnes","ISAM-CESM yield",yieldagfi1,"t/ha","production",yieldgi1,"tonnes","CLM yield",yieldagfcn,"production",yieldgcn
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
a1=1981
a2=2006
x=a2-a1
toarea= N.zeros(x)
zumiy= N.zeros(x)
clmy= N.zeros(x)
clmyn=N.zeros(x)
clmpn=N.zeros(x)
isamy= N.zeros(x)
zumip= N.zeros(x)
clmp= N.zeros(x)
isamp= N.zeros(x)
isamy1= N.zeros(x)
isamp1= N.zeros(x)
faoy=N.zeros(x)
faop=N.zeros(x)
faoy1=[17538,17586,16201,17143,19063,18199,19052,17047,18288,18958,18791,20379,19354,21833,20310,21313,21567,22561,21899,21689,23209,23009,22795,22437,23179]
faop1=[88525040,92121684,79467061,90752915,101156845,94446350,100102463,93521958,107254261,108456443,103322537,114467678,115146195,136447144,126949296,130202572,144356733,160135008,157776694,161297357,178242605,181677687,190650566,205522972,214559072]
c=0.0001
faoy=N.multiply(faoy1,c)
faop=N.multiply(faop1,1)
#name=["Global","USA","China","Brazil","Argentina","Mexico","India","Ukraine","Indonesia","France","Southafrica"]
#range1=[10100,11501,41501,20301,20101,11101,42901,47800,50300,42200,33900]
#range2=[50700,11550,41529,20327,20124,11132,42929,47800,50300,42200,33900]

#name=["Italy","Canada","Vietnam","Hungary","Romania","Philippines","Thailand","Chile","Spain","Nigeria","Germany"]
#range1=[43400,10201,48200,42700,46000,50700,47500,20400,46800,33400,42500]
#range2=[43400,10212,48200,42700,46000,50700,47500,20400,46800,33400,42500]

name=["Global"]
range1=[10100]
range2=[50700]




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
        clmyn=N.zeros(x)
        clmpn=N.zeros(x)
	for num in range(a1,a2):
    
	    reu=annualyield(num,range1[i],range2[i])

	    toarea[a]=reu[2]
	    zumiy[a]=reu[5]
	    isamy[a]=reu[11]
	    clmy[a]=reu[17]
	    zumip[a]=reu[8]
	    isamp[a]=reu[14]
	    clmp[a]=reu[20]
	    isamy1[a]=reu[23]
	    isamp1[a]=reu[26]
            clmyn[a]=reu[29]
            clmpn[a]=reu[31]

	    a=a+1

	azumiy=zumiy-N.average(zumiy)
	aisamy=isamy-N.average(isamy)
	aclmy=clmy-N.average(clmy)
	azumip=zumip-N.average(zumip)
	aisamp=isamp-N.average(isamp)
	aclmp=clmp-N.average(clmp)
	aisamy1=isamy1-N.average(isamy1)
	aisamp1=isamp1-N.average(isamp1)
        afaoy=faoy-N.average(faoy)
        afaop=faop-N.average(faop)
        aclmyn=clmyn-N.average(clmyn)
        aclmpn=clmpn-N.average(clmpn)


	mean_zumiy=runmean(x,zumiy)
	mean_isamy=runmean(x,isamy)
	mean_clmy=runmean(x,clmy)
	mean_zumip=runmean(x,zumip)
	mean_isamp=runmean(x,isamp)
	mean_clmp=runmean(x,clmp)
	mean_isamy1=runmean(x,isamy1)
	mean_isamp1=runmean(x,isamp1)
        mean_faoy=runmean(x,faoy)
        mean_faop=runmean(x,faop)
        mean_clmyn=runmean(x,clmyn)
        mean_clmpn=runmean(x,clmpn)


	fig = plt.figure(figsize=(26,20))


	ax = fig.add_subplot(321)
	xx=range(a1,a2)

#plt.ylim((0,12))

	xdates = [datetime.datetime.strptime(str(int(date)),'%Y') for date in xx]
#plt.xticks(xdates, xdates)

	ax.plot_date(xdates,faoy,"ro-",label="FAO",linewidth=2)
	ax.plot_date(xdates,isamy,"go-",label="ISAM-NCEP")
	ax.plot_date(xdates,isamy1,"yo-",label="ISAM-CESM")
	ax.plot_date(xdates,clmy,"bo-",label="CLM-N")
        ax.plot_date(xdates,clmyn,"ko-",label="CLM")

	ax.xaxis.set_major_formatter(DateFormatter('%Y'))
	ccp=scipy.stats.pearsonr(isamy,faoy)
	ccp1=scipy.stats.pearsonr(isamy1,faoy)
	ccp2=scipy.stats.pearsonr(clmy,faoy)
        ccp3=scipy.stats.pearsonr(clmyn,faoy)

#	leg=plt.legend(['FAO', 'ISAM-NCEP {:05.3f}'.format(ccp[0]),'ISAM-CESM {:05.3f}'.format(ccp1[0]),'CLM-N {:05.3f}'.format(ccp2[0]),'CLM {:05.3f}'.format(ccp3[0])])
        leg=plt.legend(['FAO', 'ISAM-NCEP {:04.2f}'.format(ccp[0]),'ISAM-CESM {:04.2f}'.format(ccp1[0]),'CLM w/ N {:04.2f}'.format(ccp2[0]),'CLM w/o N {:04.2f}'.format(ccp3[0])],loc=4,fontsize=18)
#leg = plt.legend(loc=2,fancybox=True, fontsize=14)
	leg.get_frame().set_alpha(0.5)

	plt.xlabel("Year",fontsize=18)
	plt.ylabel("Soybean yield (t/ha)",fontsize=18)
	plt.title("Yield {0}".format(name1),fontsize=18)
	plt.tick_params(axis='both',labelsize=18)


	ax = fig.add_subplot(322)
	ax.plot_date(xdates,faop,"ro-",label="FAO",linewidth=2)
	ax.plot_date(xdates,isamp,"go-",label="ISAM-NCEP")
	ax.plot_date(xdates,isamp1,"yo-",label="ISAM-CESM")
	ax.plot_date(xdates,clmp,"bo-",label="CLM-N")
        ax.plot_date(xdates,clmpn,"ko-",label="CLM")

	ax.xaxis.set_major_formatter(DateFormatter('%Y'))

	ccp=scipy.stats.pearsonr(isamp,faop)
	ccp1=scipy.stats.pearsonr(isamp1,faop)
	ccp2=scipy.stats.pearsonr(clmp,faop)
        ccp3=scipy.stats.pearsonr(clmpn,faop)
#	leg=plt.legend(['FAO', 'ISAM-NCEP {:05.3f}'.format(ccp[0]),'ISAM-CESM {:05.3f}'.format(ccp1[0]),'CLM-N {:05.3f}'.format(ccp2[0]),'CLM {:05.3f}'.format(ccp3[0])])
        leg=plt.legend(['FAO', 'ISAM-NCEP {:04.2f}'.format(ccp[0]),'ISAM-CESM {:04.2f}'.format(ccp1[0]),'CLM w/ N {:04.2f}'.format(ccp2[0]),'CLM w/o N {:04.2f}'.format(ccp3[0])],loc=4,fontsize=18)
#leg = plt.legend(loc=2,fancybox=True, fontsize=14)
	leg.get_frame().set_alpha(0.5)

	plt.title("Production {0}".format(name1),fontsize=18)
	plt.xlabel("Year",fontsize=18)
	plt.ylabel("Soybean production (tonnes)",fontsize=18)
	plt.tick_params(axis='both',labelsize=18)



	ax = fig.add_subplot(323)
	ax.plot_date(xdates,afaoy,"ro-",label="FAO",linewidth=2)
	ax.plot_date(xdates,aisamy,"go-",label="ISAM-NCEP")
	ax.plot_date(xdates,aisamy1,"yo-",label="ISAM-CESM")
	ax.plot_date(xdates,aclmy,"bo-",label="CLM-N")
        ax.plot_date(xdates,aclmyn,"ko-",label="CLM")
	ax.xaxis.set_major_formatter(DateFormatter('%Y'))

	ccp=scipy.stats.pearsonr(aisamy,afaoy)
	ccp1=scipy.stats.pearsonr(aisamy1,afaoy)
	ccp2=scipy.stats.pearsonr(aclmy,afaoy)
        ccp3=scipy.stats.pearsonr(aclmyn,afaoy)
#	leg=plt.legend(['FAO', 'ISAM-NCEP {:05.3f}'.format(ccp[0]),'ISAM-CESM {:05.3f}'.format(ccp1[0]),'CLM-N {:05.3f}'.format(ccp2[0]),'CLM {:05.3f}'.format(ccp3[0])])
        leg=plt.legend(['FAO', 'ISAM-NCEP {:04.2f}'.format(ccp[0]),'ISAM-CESM {:04.2f}'.format(ccp1[0]),'CLM w/ N {:04.2f}'.format(ccp2[0]),'CLM w/o N {:04.2f}'.format(ccp3[0])],loc=4,fontsize=18)
#leg = plt.legend(loc=2,fancybox=True, fontsize=14)
	leg.get_frame().set_alpha(0.5)


	plt.title("Anormaly by subtracting average {0}".format(name1),fontsize=18)
	plt.xlabel("Year",fontsize=18)
	plt.ylabel("Soybean yield (t/ha)",fontsize=18)
	plt.tick_params(axis='both',labelsize=18)



	ax = fig.add_subplot(324)
	ax.plot_date(xdates,afaop,"ro-",label="FAO",linewidth=2)
	ax.plot_date(xdates,aisamp,"go-",label="ISAM-NCEP")
	ax.plot_date(xdates,aisamp1,"yo-",label="ISAM-CESM")
	ax.plot_date(xdates,aclmp,"bo-",label="CLM-N")
        ax.plot_date(xdates,aclmpn,"ko-",label="CLM")
	ax.xaxis.set_major_formatter(DateFormatter('%Y'))

	ccp=scipy.stats.pearsonr(aisamp,afaop)
	ccp1=scipy.stats.pearsonr(aisamp1,afaop)
	ccp2=scipy.stats.pearsonr(aclmp,afaop)
        ccp3=scipy.stats.pearsonr(aclmpn,afaop)
#	leg=plt.legend(['FAO', 'ISAM-NCEP {:05.3f}'.format(ccp[0]),'ISAM-CESM {:05.3f}'.format(ccp1[0]),'CLM-N {:05.3f}'.format(ccp2[0]),'CLM {:05.3f}'.format(ccp3[0])])
        leg=plt.legend(['FAO', 'ISAM-NCEP {:04.2f}'.format(ccp[0]),'ISAM-CESM {:04.2f}'.format(ccp1[0]),'CLM w/ N {:04.2f}'.format(ccp2[0]),'CLM w/o N {:04.2f}'.format(ccp3[0])],loc=4,fontsize=18)
#leg = plt.legend(loc=2,fancybox=True, fontsize=14)
	leg.get_frame().set_alpha(0.5)


	plt.title("Anormaly by subtracting average {0}".format(name1),fontsize=18)
	plt.xlabel("Year",fontsize=18)
	plt.ylabel("Soybean production (tonnes)",fontsize=18)
	plt.tick_params(axis='both',labelsize=18)



	ax = fig.add_subplot(325)
	ax.plot_date(xdates,mean_faoy,"ro-",label="FAO",linewidth=2)
	ax.plot_date(xdates,mean_isamy,"go-",label="ISAM-NCEP")
	ax.plot_date(xdates,mean_isamy1,"yo-",label="ISAM-CESM")
	ax.plot_date(xdates,mean_clmy,"bo-",label="CLM-N")
        ax.plot_date(xdates,mean_clmyn,"ko-",label="CLM")
	ax.xaxis.set_major_formatter(DateFormatter('%Y'))

	mean1_isamy=N.zeros([x-2])
	mean1_faoy=N.zeros([x-2])
	mean1_isamy1=N.zeros([x-2])
	mean1_clmy=N.zeros([x-2])
        mean1_clmyn=N.zeros([x-2])

	for i in range(1,x-2):
		mean1_isamy[i]=mean_isamy[i]
		mean1_faoy[i]=mean_faoy[i]
		mean1_isamy1[i]=mean_isamy1[i]
		mean1_clmy[i]=mean_clmy[i]
                mean1_clmyn[i]=mean_clmyn[i]
	ccp=scipy.stats.pearsonr(mean1_isamy,mean1_faoy)
	ccp1=scipy.stats.pearsonr(mean1_isamy1,mean1_faoy)
	ccp2=scipy.stats.pearsonr(mean1_clmy,mean1_faoy)
        ccp3=scipy.stats.pearsonr(mean1_clmyn,mean1_faoy)
#	leg=plt.legend(['FAO', 'ISAM-NCEP {:05.3f}'.format(ccp[0]),'ISAM-CESM {:05.3f}'.format(ccp1[0]),'CLM-N {:05.3f}'.format(ccp2[0]),'CLM {:05.3f}'.format(ccp3[0])])
        leg=plt.legend(['FAO', 'ISAM-NCEP {:04.2f}'.format(ccp[0]),'ISAM-CESM {:04.2f}'.format(ccp1[0]),'CLM w/ N {:04.2f}'.format(ccp2[0]),'CLM w/o N {:04.2f}'.format(ccp3[0])],loc=4,fontsize=18)
#leg = plt.legend(loc=2,fancybox=True, fontsize=14)
	leg.get_frame().set_alpha(0.5)


	plt.title("Anormaly by subtracting a moving mean average {0}".format(name1),fontsize=18)
	plt.xlabel("Year",fontsize=18)
	plt.ylabel("Soybean yield (t/ha)",fontsize=18)
	plt.tick_params(axis='both',labelsize=18)



	ax = fig.add_subplot(326)
	ax.plot_date(xdates,mean_faop,"ro-",label="FAO",linewidth=2)
	ax.plot_date(xdates,mean_isamp,"go-",label="ISAM-NCEP")
	ax.plot_date(xdates,mean_isamp1,"yo-",label="ISAM-CESM")
	ax.plot_date(xdates,mean_clmp,"bo-",label="CLM-N")
        ax.plot_date(xdates,mean_clmpn,"ko-",label="CLM")

	ax.xaxis.set_major_formatter(DateFormatter('%Y'))

	mean1_isamp=N.zeros([x-2])
	mean1_faop=N.zeros([x-2])
	mean1_isamp1=N.zeros([x-2])
	mean1_clmp=N.zeros([x-2])
        mean1_clmpn=N.zeros([x-2])

	for i in range(1,x-2):
		mean1_isamp[i]=mean_isamp[i]
		mean1_faop[i]=mean_faop[i]
		mean1_isamp1[i]=mean_isamp1[i]
		mean1_clmp[i]=mean_clmp[i]
                mean1_clmpn[i]=mean_clmpn[i]

	ccp=scipy.stats.pearsonr(mean1_isamp,mean1_faop)
	ccp1=scipy.stats.pearsonr(mean1_isamp1,mean1_faop)
	ccp2=scipy.stats.pearsonr(mean1_clmp,mean1_faop)
        ccp3=scipy.stats.pearsonr(mean1_clmpn,mean1_faop)

#	leg=plt.legend(['FAO', 'ISAM-NCRP {:05.3f}'.format(ccp[0]),'ISAM-CESM {:05.3f}'.format(ccp1[0]),'CLM-N {:05.3f}'.format(ccp2[0]),'CLM {:05.3f}'.format(ccp3[0])])
        leg=plt.legend(['FAO', 'ISAM-NCEP {:04.2f}'.format(ccp[0]),'ISAM-CESM {:04.2f}'.format(ccp1[0]),'CLM w/ N {:04.2f}'.format(ccp2[0]),'CLM w/o N {:04.2f}'.format(ccp3[0])],loc=4,fontsize=18)
#leg = plt.legend(loc=2,fancybox=True, fontsize=14)
	leg.get_frame().set_alpha(0.5)


	plt.title("Anormaly by subtracting a moving mean average {0}".format(name1),fontsize=18)
	plt.xlabel("Year",fontsize=18)
	plt.ylabel("Soybean production (tonnes)",fontsize=18)
	plt.tick_params(axis='both',labelsize=18)


	plt.savefig('soy_fao1_{0}.png'.format(name1),dpi=300,bbox_inches='tight')
#plt.show()


