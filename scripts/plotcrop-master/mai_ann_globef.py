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
    bb=year-2006
    
    region1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/RCP85_crop_150901.nc','r')
    maitrop = region1.variables['maize_trop'][bb,:,:]
    maitemp = region1.variables['maize_temp'][bb,:,:]
    maitropi=region1.variables['maize_trop_irrig'][bb,:,:]
    maitempi=region1.variables['maize_temp_irrig'][bb,:,:]
    gridarea = region1.variables['area'][:,:]

    maitrop=ma.masked_where(maitrop<=0,maitrop)
    maitrop=ma.filled(maitrop, fill_value=0.)
    maitemp=ma.masked_where(maitemp<=0,maitemp)
    maitemp=ma.filled(maitemp, fill_value=0.)

    maitropi=ma.masked_where(maitropi<=0,maitropi)
    maitropi=ma.filled(maitropi, fill_value=0.)
    maitempi=ma.masked_where(maitempi<=0,maitempi)
    maitempi=ma.filled(maitempi, fill_value=0.)

    maizetro=maitrop+maitropi
    maizetem=maitemp+maitempi
    maizetor=maitrop+maitemp
    maizetoi=maitropi+maitempi
    maizeto = maitrop+maitemp+maitropi+maitempi


    clm=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp85/maizetrop_rcp85_constco2_rf_nofert_0.5x0.5.nc','r')
    clmtropf = clm.variables['yield'][bb,:,:]


    clm1=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp85/maizetemp_rcp85_constco2_rf_nofert_0.5x0.5.nc','r')
    clmtempf = clm1.variables['yield'][bb,:,:]

    clm2=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp85/maizetrop_rcp85_co2_rf_nofert_0.5x0.5.nc','r')
    clmtropfi = clm2.variables['yield'][bb,:,:]


    clm3=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp85/maizetemp_rcp85_co2_rf_nofert_0.5x0.5.nc','r')
    clmtempfi = clm3.variables['yield'][bb,:,:]


    clmtropf=N.flipud(clmtropf)
    clmtempf=N.flipud(clmtempf)
    clmtropfi=N.flipud(clmtropfi)
    clmtempfi=N.flipud(clmtempfi)

    clmtropf= ma.masked_where(maizetro<=0,clmtropf)
    clmtempf= ma.masked_where(maizetem<=0,clmtempf)
    clmtropf=ma.filled(clmtropf, fill_value=0.)
    clmtempf=ma.filled(clmtempf, fill_value=0.)

    clmtropfi= ma.masked_where(maizetro<=0,clmtropfi)
    clmtempfi= ma.masked_where(maizetem<=0,clmtempfi)
    clmtropfi=ma.filled(clmtropfi, fill_value=0.)
    clmtempfi=ma.filled(clmtempfi, fill_value=0.)

    yield_clmtf=clmtropf+clmtempf
    yield_clmtf = ma.masked_where(yield_clmtf<=0,yield_clmtf)
    yield_clmtf  = ma.masked_where(maizeto<=0,yield_clmtf )
    yield_clmtf=ma.filled(yield_clmtf, fill_value=0.)

    yield_clmtfi=clmtropfi+clmtempfi
    yield_clmtfi = ma.masked_where(yield_clmtfi<=0,yield_clmtfi)
    yield_clmtfi = ma.masked_where(maizeto<=0,yield_clmtfi)
    yield_clmtfi=ma.filled(yield_clmtfi, fill_value=0.)


    clmn=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp85/maizetrop_rcp85_co2_rf_fert_0.5x0.5.nc','r')
    clmtropfn = clmn.variables['yield'][bb,:,:]


    clm1n=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp85/maizetemp_rcp85_co2_rf_fert_0.5x0.5.nc','r')
    clmtempfn = clm1n.variables['yield'][bb,:,:]

    clm2n=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp85/maizetrop_rcp85_co2_irrig_fert_0.5x0.5.nc','r')
    clmtropfin = clm2n.variables['yield'][bb,:,:]


    clm3n=NetCDFFile('/scratch2/scratchdirs/tslin2/plot/globalcrop/data/clm/clm45rcp85/maizetemp_rcp85_co2_irrig_fert_0.5x0.5.nc','r')
    clmtempfin = clm3n.variables['yield'][bb,:,:]


    clmtropfn=N.flipud(clmtropfn)
    clmtempfn=N.flipud(clmtempfn)
    clmtropfin=N.flipud(clmtropfin)
    clmtempfin=N.flipud(clmtempfin)

    clmtropfn= ma.masked_where(maizetro<=0,clmtropfn)
    clmtempfn= ma.masked_where(maizetem<=0,clmtempfn)
    clmtropfn=ma.filled(clmtropfn, fill_value=0.)
    clmtempfn=ma.filled(clmtempfn, fill_value=0.)

    clmtropfin= ma.masked_where(maizetro<=0,clmtropfin)
    clmtempfin= ma.masked_where(maizetem<=0,clmtempfin)
    clmtropfin=ma.filled(clmtropfin, fill_value=0.)
    clmtempfin=ma.filled(clmtempfin, fill_value=0.)

    yield_clmtfn=clmtropfn+clmtempfn
    yield_clmtfn = ma.masked_where(yield_clmtfn<=0,yield_clmtfn)
    yield_clmtfn  = ma.masked_where(maizeto<=0,yield_clmtfn )
    yield_clmtfn=ma.filled(yield_clmtfn, fill_value=0.)

    yield_clmtfin=clmtropfin+clmtempfin
    yield_clmtfin = ma.masked_where(yield_clmtfin<=0,yield_clmtfin)
    yield_clmtfin = ma.masked_where(maizeto<=0,yield_clmtfin)
    yield_clmtfin=ma.filled(yield_clmtfin, fill_value=0.)


    yieldfa= N.zeros((1, 360, 720))
    yieldf= N.zeros((1, 360, 720))
    yieldfi= N.zeros((1, 360, 720))
    yieldfai= N.zeros((1, 360, 720))
    year1=year+1
    years = range(year,year1)
    for i, year in enumerate(years):
        base = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_future/rcp85am/output/rcp85am.bgp-yearly_crop_{0}.nc".format(year), mode='r')
        base2 = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_future/rcp85bm/output/rcp85bm.bgp-yearly_crop_{0}.nc".format(year), mode='r')

        lona1 = base.variables["lon"][:]
        lata1 = base.variables["lat"][:]
        yield1 = base.variables["yield"][0,:,:]
        yieldf[i, :, :] = yield1
        yield2 = base2.variables["yield"][0,:,:]
        yieldfa[i, :, :] = yield2

        basei = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_future/rcp85cm/output/rcp85cm.bgp-yearly_crop_{0}.nc".format(year), mode='r')
        base2i = NetCDFFile ("/scratch2/scratchdirs/tslin2/isam/model/global/isam_future/rcp85dm/output/rcp85dm.bgp-yearly_crop_{0}.nc".format(year), mode='r')

        yield1i = basei.variables["yield"][0,:,:]
        yieldfi[i, :, :] = yield1i
        yield2i = base2i.variables["yield"][0,:,:]
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
    yield_new = ma.masked_where(maizeto<=0,yield_new)
    yield_new=ma.filled(yield_new, fill_value=0.)

    yield_new1[N.isnan(yield_new1)] = -9999
    yield_new1 = ma.masked_where(yield_new1<=0,yield_new1)
    yield_new1 = ma.masked_where(maizeto<=0,yield_new1)
    yield_new1=ma.filled(yield_new1, fill_value=0.)

    yield_newi[N.isnan(yield_newi)] = -9999
    yield_newi = ma.masked_where(yield_newi<=0,yield_newi)
    yield_newi = ma.masked_where(maizeto<=0,yield_newi)
    yield_newi=ma.filled(yield_newi, fill_value=0.)

    yield_new1i[N.isnan(yield_new1i)] = -9999
    yield_new1i = ma.masked_where(yield_new1i<=0,yield_new1i)
    yield_new1i = ma.masked_where(maizeto<=0,yield_new1i)
    yield_new1i=ma.filled(yield_new1i, fill_value=0.)



    
    yieldagf=0.
    yieldg=0.
    harea=0.
    a=0
    yieldgid=0.
    yieldgia=0.
    yieldgib=0.
    yieldgic=0.

    yieldgca=0.
    yieldgcb=0.
    yieldgcc=0.
    yieldgcd=0.

    yieldagfa=0.
    yieldagfb=0.
    yieldagfc=0.
    yieldagfd=0.
    yieldagfa1=0.
    yieldagfb1=0.
    yieldagfc1=0.
    yieldagfd1=0.



    #USA 11501~11550
    for xx in range(0,360):
        for yy in range(0,720):
            if coun[xx,yy] >=couna and  coun[xx,yy] <=counb:
                yieldg=0+yieldg
                harea=maizeto[xx,yy]*gridarea[xx,yy]+ harea

                yieldgia=(yield_new[xx,yy]*gridarea[xx,yy]*maizeto[xx,yy])+yieldgia
                yieldgca=(yield_clmtf[xx,yy]*gridarea[xx,yy]*maizeto[xx,yy])+yieldgca
                yieldgib=(yield_new1[xx,yy]*gridarea[xx,yy]*maizeto[xx,yy])+yieldgib
                yieldgcc=(yield_clmtfn[xx,yy]*gridarea[xx,yy]*maizeto[xx,yy])+yieldgcc

                yieldgic=(yield_newi[xx,yy]*gridarea[xx,yy]*maizeto[xx,yy])+yieldgic
                yieldgcb=(yield_clmtfi[xx,yy]*gridarea[xx,yy]*maizeto[xx,yy])+yieldgcb
                yieldgid=(yield_new1i[xx,yy]*gridarea[xx,yy]*maizeto[xx,yy])+yieldgid
                yieldgcd=(yield_clmtfin[xx,yy]*gridarea[xx,yy]*maizeto[xx,yy])+yieldgcd

                a=a+1

    yieldagf=yieldg/harea

    yieldagfa=yieldgia/harea
    yieldagfb=yieldgib/harea
    yieldagfc=yieldgic/harea
    yieldagfd=yieldgid/harea
   
    yieldagfa1=yieldgca/harea
    yieldagfb1=yieldgcb/harea
    yieldagfc1=yieldgcc/harea
    yieldagfd1=yieldgcd/harea

    return "Maize ","iizumi harvested area",harea, "ha","iizumi yield",yieldagf,"t/ha","production",yieldg,"tonnes","ISAM-a yield",yieldagfa,"t/ha","CLM-a yield",yieldagfa1,"t/ha","ISAM-b yield",yieldagfb,"t/ha","clm-a yield",yieldagfb1,"tonnes","ISAM-c yield",yieldagfc,"t/ha","clm-c yield",yieldagfc1,"t/ha","ISAM-d yield",yieldagfd,"clm-a yield",yieldagfd1
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

a1=2010
a2=2060
x=a2-a1
toarea= N.zeros(x)
zumiy= N.zeros(x)
clmya= N.zeros(x)
clmyb=N.zeros(x)
clmyc=N.zeros(x)
isamyd= N.zeros(x)
zumip= N.zeros(x)
clmyd= N.zeros(x)
isamya= N.zeros(x)
isamyb= N.zeros(x)
isamyc= N.zeros(x)
#faoy1=[34933.0,36090.0,29452.0,35256.0,37202.0,36279.,34863.,31001.,36186.,36907.,36967.,39009.,36295.,41236.,38092.,42223.,41489.,44362.,44247.,43244.,44756.,43855.,44582.,49437.,48193.]
#faop1=[446772517,448932280,347082034,450449992,485527301,478176622,453115794,403050234,476874503,483620724,494476131,533599651,476783169,569024180,517293135,589483113,585521420,615783169,607190213,592467104,615534635,603480524,645096760,728898458,713663038]
#c=0.0001
#faoy=N.multiply(faoy1,c)
#faop=N.multiply(faop1,1)
#name=["Global","USA","China","Brazil","Argentina","Mexico","India","Ukraine","Indonesia","France","Southafrica"]
#range1=[10100,11501,41501,20301,20101,11101,42901,47800,50300,42200,33900]
#range2=[50700,11550,41529,20327,20124,11132,42929,47800,50300,42200,33900]

#name=["Italy","Canada","Vietnam","Hungary","Romania","Philippines","Thailand","Chile","Spain","Nigeria","Germany"]
#range1=[43400,10201,48200,42700,46000,50700,47500,20400,46800,33400,42500]
#range2=[43400,10212,48200,42700,46000,50700,47500,20400,46800,33400,42500]

name=["Global"]
range1=[10100]
range2=[50700]

#name=["USA"]
#range1=[11501]
#range2=[11550]





for i, name1 in enumerate(name):
        a=0
	toarea= N.zeros(x)
	zumiy= N.zeros(x)
	clmya= N.zeros(x)
	isamya= N.zeros(x)
	zumip= N.zeros(x)
	clmyb= N.zeros(x)
	isamyb= N.zeros(x)
	isamyc= N.zeros(x)
	isamyd= N.zeros(x)
        clmyc=N.zeros(x)
        clmyd=N.zeros(x)
	for num in range(a1,a2):
    
	    reu=annualyield(num,range1[i],range2[i])

	    toarea[a]=reu[2]
	    zumiy[a]=reu[5]
	    isamya[a]=reu[11]
	    isamyb[a]=reu[17]
	    zumip[a]=reu[8]
	    clmya[a]=reu[14]
	    clmyb[a]=reu[20]
	    isamyc[a]=reu[23]
	    clmyc[a]=reu[26]
            isamyd[a]=reu[29]
            clmyd[a]=reu[31]

	    a=a+1

	azumiy=zumiy-N.average(zumiy)
	azumip=zumip-N.average(zumip)
	aisamya=isamya-N.average(isamya)
	aisamyb=isamyb-N.average(isamyb)
	aisamyc=isamyc-N.average(isamyc)
	aisamyd=isamyd-N.average(isamyd)

        aclmya=clmya-N.average(clmya)
        aclmyb=clmyb-N.average(clmyb)
        aclmyc=clmyc-N.average(clmyc)
        aclmyd=clmyd-N.average(clmyd)


	mean_zumiy=runmean(x,zumiy)
	mean_isamya=runmean(x,isamya)
	mean_clmya=runmean(x,clmya)
	mean_zumip=runmean(x,zumip)
	mean_isamyb=runmean(x,isamyb)
	mean_clmyb=runmean(x,clmyb)
        mean_isamyc=runmean(x,isamyc)
        mean_clmyc=runmean(x,clmyc)
        mean_isamyd=runmean(x,isamyd)
        mean_clmyd=runmean(x,clmyd)


	fig = plt.figure(figsize=(26,20))


	ax = fig.add_subplot(321)
	xx=range(a1,a2)

#plt.ylim((0,12))

	xdates = [datetime.datetime.strptime(str(int(date)),'%Y') for date in xx]
#plt.xticks(xdates, xdates)

	ax.plot_date(xdates,isamya,"go-",label="ISAM-A")
	ax.plot_date(xdates,isamyb,"yo-",label="ISAM-B")
	ax.plot_date(xdates,isamyc,"bo-",label="ISAM-C")
        ax.plot_date(xdates,isamyd,"ko-",label="ISAM-D")

	ax.xaxis.set_major_formatter(DateFormatter('%Y'))

	plt.xlabel("Year",fontsize=18)
	plt.ylabel("Maize yield (t/ha)",fontsize=18)
	plt.title("Yield",fontsize=18)
	plt.tick_params(axis='both',labelsize=18)

        leg = plt.legend(loc=2,fancybox=True, fontsize=16)
        leg.get_frame().set_alpha(0.5)

	ax = fig.add_subplot(322)
        ax.plot_date(xdates,clmya,"g*-",label="CLM-A")
        ax.plot_date(xdates,clmyb,"y*-",label="CLM-B")
        ax.plot_date(xdates,clmyc,"b*-",label="CLM-C")
        ax.plot_date(xdates,clmyd,"k*-",label="CLM-D")

	ax.xaxis.set_major_formatter(DateFormatter('%Y'))

	plt.title("Maize yield",fontsize=18)
	plt.xlabel("Year",fontsize=18)
	plt.ylabel("Maize Yield (t/ha)",fontsize=18)
	plt.tick_params(axis='both',labelsize=18)

        leg = plt.legend(loc=2,fancybox=True, fontsize=16)
        leg.get_frame().set_alpha(0.5)


	ax = fig.add_subplot(323)


        ax.plot_date(xdates,aisamya,"go-",label="ISAM-A")
        ax.plot_date(xdates,aisamyb,"yo-",label="ISAM-B")
        ax.plot_date(xdates,aisamyc,"bo-",label="ISAM-C")
        ax.plot_date(xdates,aisamyd,"ko-",label="ISAM-D")


	ax.xaxis.set_major_formatter(DateFormatter('%Y'))


	plt.title("Anormaly by subtracting average",fontsize=18)
	plt.xlabel("Year",fontsize=18)
	plt.ylabel("Maize yield (t/ha)",fontsize=18)
	plt.tick_params(axis='both',labelsize=18)

        leg = plt.legend(loc=2,fancybox=True, fontsize=16)
        leg.get_frame().set_alpha(0.5)


	ax = fig.add_subplot(324)

        ax.plot_date(xdates,aclmya,"g*-",label="CLM-A")
        ax.plot_date(xdates,aclmyb,"y*-",label="CLM-B")
        ax.plot_date(xdates,aclmyc,"b*-",label="CLM-C")
        ax.plot_date(xdates,aclmyd,"k*-",label="CLM-D")

	plt.title("Anormaly by subtracting average",fontsize=18)
	plt.xlabel("Year",fontsize=18)
	plt.ylabel("Maize yield (t/ha)",fontsize=18)
	plt.tick_params(axis='both',labelsize=18)

        leg = plt.legend(loc=2,fancybox=True, fontsize=16)
        leg.get_frame().set_alpha(0.5)


	ax = fig.add_subplot(325)

        ax.plot_date(xdates,mean_isamya,"go-",label="ISAM-A")
        ax.plot_date(xdates,mean_isamyb,"yo-",label="ISAM-B")
        ax.plot_date(xdates,mean_isamyc,"bo-",label="ISAM-C")
        ax.plot_date(xdates,mean_isamyd,"ko-",label="ISAM-D")

	ax.xaxis.set_major_formatter(DateFormatter('%Y'))


	plt.title("Anormaly by subtracting a moving mean average",fontsize=18)
	plt.xlabel("Year",fontsize=18)
	plt.ylabel("Maize yield (t/ha)",fontsize=18)
	plt.tick_params(axis='both',labelsize=18)

        leg = plt.legend(loc=2,fancybox=True, fontsize=16)
        leg.get_frame().set_alpha(0.5)


	ax = fig.add_subplot(326)
        ax.plot_date(xdates,mean_clmya,"g*-",label="CLM-A")
        ax.plot_date(xdates,mean_clmyb,"y*-",label="CLM-B")
        ax.plot_date(xdates,mean_clmyc,"b*-",label="CLM-C")
        ax.plot_date(xdates,mean_clmyd,"k*-",label="CLM-D")


	ax.xaxis.set_major_formatter(DateFormatter('%Y'))



	plt.title("Anormaly by subtracting a moving mean average",fontsize=18)
	plt.xlabel("Year ",fontsize=18)
	plt.ylabel("Maize yield (t/ha)",fontsize=18)
	plt.tick_params(axis='both',labelsize=18)

        leg = plt.legend(loc=2,fancybox=True, fontsize=16)
        leg.get_frame().set_alpha(0.5)

	plt.savefig('maize_rcp85_{0}.png'.format(name1),dpi=300,bbox_inches='tight')
plt.show()


