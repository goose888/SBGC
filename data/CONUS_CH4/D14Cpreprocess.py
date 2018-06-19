# -*- coding: utf-8 -*-
"""
Module to preprocess D14C synthesis data

Created on Thu Feb 12 20:28:02 2015

@author: Yujie
"""

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy.ma as ma
import math 
import matplotlib
from scipy.interpolate import interp1d
from scipy import array
from scipy.optimize import curve_fit

def extrap1d(interpolator):
    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
            if np.all(~np.isnan([xs[0], xs[1], ys[0], ys[1]])):
                return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
            else:
                return ys[0]
        elif x > xs[-1]:
            if np.all(~np.isnan([ys[-1],ys[-2],xs[-1],xs[-2]])):
                return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
            else:
                return ys[-1]
        else:
            return interpolator(x)

    def ufunclike(xs):
        return array(map(pointwise, array(xs)))

    return ufunclike

def extrap1d_constby(interpolator):
    ''' extrapolation same as 'extrap1d' but use constant beyond x[0] and x[-1]
    '''
    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0]
        elif x > xs[-1]:
            return ys[-1]
        else:
            return interpolator(x)

    def ufunclike(xs):
        return array(map(pointwise, array(xs)))

    return ufunclike
    
def getCweightedD14C(filename,cutdep=None):
    """ Get SOC averaged D14C for each profile. 
    arguments: 1 - csv filename
    return: 6 col nparray, 
            1 - layer top    
            2 - layer bottom
            3 - Profile ID
            4 - SOC-averaged D14C
            5 - totsoc for the whole profile in kgC/m2
    """
    data = pd.read_csv(filename,encoding='iso-8859-1')  
    pid = data[data['Start_of_Profile']==1].index   # index of profile start
    out = np.zeros((np.nansum(data['Start_of_Profile']),6)) 
    for i in range(0,pid.shape[0]-1):
        #print 'i is :',i
        for j in range(pid[i],pid[i+1]):
            if np.all(~np.isnan((
                    data.loc[j,['D14C_BulkLayer','BulkDensity','pct_C']].astype(float)))):   
                if cutdep is not None:
                    if data.loc[j,['Layer_bottom']].values.astype(float) <= cutdep + 30.:
                        cdi = j
                else:
                    cdi = j
        if cdi >= pid[i]:
                out[i,1] = data.loc[pid[i],['Layer_top']]
                out[i,2] = data.loc[cdi,['Layer_bottom']]
                out[i,3] = data.loc[pid[i],['ProfileID']]
                totsoc = np.array(data.loc[pid[i]:cdi,['BulkDensity']]).astype(float) * \
                        (np.array(data.loc[pid[i]:cdi,['Layer_bottom']].astype(float)) - \
                        np.array(data.loc[pid[i]:cdi,['Layer_top']].astype(float))) * \
                        np.array(data.loc[pid[i]:cdi,['pct_C']].astype(float))/100.0
                out[i,4] = np.nansum(np.array(
                        data.loc[pid[i]:cdi,['D14C_BulkLayer']].astype(float))*totsoc,0)/np.nansum(totsoc)
                out[i,5] = np.nansum(totsoc)*10   #gC/cm2 converted to kgC/m2 for the given depth
    out[np.all(out==0,1),:] = np.nan
    return out

def getCweightedD14C2(filename,cutdep=None):
    """ Get SOC averaged D14C for each profile. with linear interpolation to cutdep.
    for profiles that are shallower than cutdep, the whole profile is used
    params:
        data: {DataFrame} df
    return: {DataFrame}. ProfileID is index
            1 - layer top    
            2 - layer bottom
            3 - SOC-averaged D14C
            4 - totsoc for the whole profile (if no cutdep) in kgC/m2 or up to cutdep
            5 - 10: Lon, Lat, MAT, MAP, VegCode, SoilOrder
    """
    data = pd.read_csv(filename,encoding='iso-8859-1', index_col='ProfileID',skiprows=[1])  
    pid = data.index.unique()   # index of profile start
    out = pd.DataFrame(columns=['Layer_top','Layer_bottom','SOCave_14C','totSOC',
                                'Lon','Lat','MAT','MAP','VegISAM',
                                'SoilOrder_LEN_USDA'],index=pid)
    marker = myplt._markeriter
    for n,i in enumerate(pid):
        print 'profile is :', i
        # number of layer = 1, calculate whole profile
        if data.loc[i:i,'Layer_top'].shape[0] == 1:
            out.iloc[n,0] = data.loc[i:i,'Layer_top'].values[0]
            out.iloc[n,1] = data.loc[i:i,'Layer_bottom'].values[0]
            totsoc = np.array(data.loc[i:i,'BulkDensity']).astype(float) * \
                     (np.array(data.loc[i:i,'Layer_bottom']).astype(float) - \
                     np.array(data.loc[i:i,'Layer_top']).astype(float)) * \
                     np.array(data.loc[i:i,'pct_C']).astype(float)/100.0
            out.iloc[n,2] = np.nansum(np.array(
                       data.loc[i:i,'D14C_BulkLayer']).astype(float)*totsoc,0)/np.nansum(totsoc)
            out.iloc[n,3] = np.nansum(totsoc)*10   #gC/cm2 converted to kgC/m2 for the given depth
            out.iloc[n,4:10] = data.loc[i:i,['Lon','Lat','MAT','MAP','VegISAM',
                                        'SoilOrder_LEN_USDA']].values[0]
        else: # multiple layers
            # calculate to cutdep 
            if cutdep is not None and float(data.loc[i:i,'Layer_bottom'].values[-1]) >= cutdep:
                out.iloc[n,0] = data.loc[i:i,'Layer_top'].values[0]
                out.iloc[n,1] = cutdep
                cutdep = float(cutdep) 
                layerbotori = np.array(data.loc[i:i,'Layer_bottom']).astype(float)
                layertop = np.array(data.loc[i:i,'Layer_top']).astype(float)
                layerbot = np.mean(np.c_[layertop, layerbotori], axis=1)
                bd = np.array(data.loc[i:i,'BulkDensity']).astype(float)
                d14C = np.array(data.loc[i:i,'D14C_BulkLayer']).astype(float)
                pctC = np.array(data.loc[i:i,'pct_C']).astype(float)  
                idx = layerbot > cutdep
                if cutdep in set(layerbot):
                    layerbotinterp = layerbot
                else:
                    layerbotinterp = np.hstack((layerbot[idx==False], cutdep, layerbot[idx==True]))
                
                notNANs = ~np.isnan(d14C)
                if sum(notNANs) < 2: # no value! skip this profile
                    continue
                else:
                    f_i = interp1d(layerbot[notNANs], d14C[notNANs]); f_x = extrap1d(f_i); 
                    d14Cinterp = f_x(layerbotinterp)

                notNANs = ~np.isnan(bd)
                if sum(notNANs) < 2: # no value!
                    continue
                else:
                    f_i = interp1d(layerbot[notNANs], bd[notNANs]); f_x = extrap1d(f_i); 
                    bdinterp = f_x(layerbotinterp)
                
                notNANs = ~np.isnan(pctC)
                if sum(notNANs) < 2: # no value!
                    continue
                else:
                    f_i = interp1d(layerbot[notNANs], pctC[notNANs]); f_x = extrap1d(f_i); 
                    pctCinterp = f_x(layerbotinterp)
                totsoc = bdinterp[layerbotinterp <= cutdep] * \
                         (layerbotinterp[layerbotinterp <= cutdep] - \
                          np.hstack((layertop[0],layerbotinterp[layerbotinterp < cutdep]))) * \
                         pctCinterp[layerbotinterp <= cutdep]/100.0
                out.iloc[n,2] = np.nansum(d14Cinterp[layerbotinterp <= cutdep]*totsoc,0)/np.nansum(totsoc)
                out.iloc[n,3] = np.nansum(totsoc)*10   #gC/cm2 converted to kgC/m2 for the given depth
                out.iloc[n,4:10] = data.loc[i:i,['Lon','Lat','MAT','MAP','VegISAM',
                                            'SoilOrder_LEN_USDA']].values[0]
                print 'interpolate to cutdep...'
#                fig = plt.figure()
#                ax = fig.add_axes([0.05,0.05,0.9,0.9])
#                ax.scatter(d14C, layerbot, marker=marker.next())
#                ax.scatter(d14Cinterp, layerbotinterp, marker=marker.next())
#                plt.gca().invert_yaxis()
#                raw_input('press Enter to continue...') 
#                plt.close()

            # calculate the whole profile
            else: 
                out.iloc[n,0] = data.loc[i:i,'Layer_top'].values[0]
                out.iloc[n,1] = data.loc[i:i,'Layer_bottom'].values[-1]
                # no missing value
                if np.all(~np.isnan(
                          data.loc[i:i,['D14C_BulkLayer','BulkDensity','pct_C']].astype(float))):
                    totsoc = np.array(data.loc[i:i,'BulkDensity'].values).astype(float) * \
                             (np.array(data.loc[i:i,'Layer_bottom'].values).astype(float) - \
                             np.array(data.loc[i:i,'Layer_top'].values).astype(float)) * \
                             np.array(data.loc[i:i,'pct_C'].values).astype(float)/100.0
                    out.iloc[n,2] = np.nansum(np.array(
                               data.loc[i:i,'D14C_BulkLayer'].values).astype(float)*totsoc,0)/np.nansum(totsoc)
                    out.iloc[n,3] = np.nansum(totsoc)*10   #gC/cm2 converted to kgC/m2 for the whole profile
                    out.iloc[n,4:10] = data.loc[i:i,['Lon','Lat','MAT','MAP','VegISAM',
                                                'SoilOrder_LEN_USDA']].values[0]
                # linear interpolate missing values at given layer
                else:   
                    layerbotori = np.array(data.loc[i:i,'Layer_bottom']).astype(float)
                    layertop = np.array(data.loc[i:i,'Layer_top']).astype(float)
                    layerbot = np.mean(np.c_[layertop, layerbotori], axis=1)                    
                    bd = data.loc[i:i,'BulkDensity'].values.astype(float)
                    d14C = data.loc[i:i,'D14C_BulkLayer'].values.astype(float)
                    pctC = data.loc[i:i,'pct_C'].values.astype(float)  
                    
                    notNANs = ~np.isnan(d14C)                    
                    if sum(notNANs) < 2: # no value! skip this profile
                        continue                        
                    else:
                        f_i = interp1d(layerbot[notNANs], d14C[notNANs]); f_x = extrap1d(f_i); 
                        d14Cinterp = f_x(layerbot)

                    notNANs = ~np.isnan(bd)
                    if sum(notNANs) < 2: # no value!
                        continue                        
                    else:
                        f_i = interp1d(layerbot[notNANs], bd[notNANs]); f_x = extrap1d(f_i); 
                        bdinterp = f_x(layerbot)

                    notNANs = ~np.isnan(pctC)
                    if sum(notNANs) < 2: # no value!
                        continue
                    else:
                        f_i = interp1d(layerbot[notNANs], pctC[notNANs]); f_x = extrap1d(f_i); 
                        pctCinterp = f_x(layerbot)

                    totsoc = bdinterp * (layerbotori - layertop) * pctCinterp/100.0
                    out.iloc[n,2] = np.nansum(d14Cinterp[layerbot <= layerbot[-1]]*totsoc,0)/np.nansum(totsoc)
                    out.iloc[n,3] = np.nansum(totsoc)*10   #gC/cm2 converted to kgC/m2 for the whole profile
                    out.iloc[n,4:10] = data.loc[i:i,['Lon','Lat','MAT','MAP','VegISAM',
                                                'SoilOrder_LEN_USDA']].values[0] 
                    print 'interpolate for missing value...'
#                    fig = plt.figure()
#                    ax = fig.add_axes([0.05,0.05,0.9,0.9])
#                    ax.scatter(d14C, layerbot, marker=marker.next())
#                    ax.scatter(d14Cinterp, layerbot, marker=marker.next())
#                    raw_input('press Enter to continue...') 
#                    plt.close()
    out.iloc[np.all(out.iloc[:,[2,3]]==0 ,1),[2,3]] = np.nan
    # use np.ix_ to broadcast index
    # or equivalently: out[:,[1,2,4,5]][np.all(out[:,[1,2,4,5]]==0,1),:] = np.nan
    return out
    
def getSOCbyLayer(filename):
    """ Get SOC (g/cm2) at each layer. 
    arguments: csv file name
    return: 2 col nparray, 1st col is profile ID, 2nd col is soc @ each layer        
    """
    data = pd.read_csv(filename,encoding='iso-8859-1')  
    pid = data[data['Start_of_Profile']==1].index   # index of profile start
    out = np.zeros((data.shape[0],2)) 
    for i in range(0,pid.shape[0]-1):
        layerdepth = range(pid[i],pid[i+1])        
        out[layerdepth,0] = np.repeat(
            data.loc[pid[i],['ProfileID']].values.astype(float),np.shape(layerdepth)[0])
        totsoc = data.loc[layerdepth,['BulkDensity']].values.astype(float)* \
            (data.loc[layerdepth,['Layer_bottom']].values.astype(float) - \
            data.loc[layerdepth,['Layer_top']].values.astype(float)) * \
            data.loc[layerdepth,['pct_C']].values.astype(float)/100.0
        out[layerdepth,1] = totsoc.flatten()*10.  # gC/cm2 converted to kgC/m2 for the given depth
    return out
    
def sweepdata(filename):
    """  sweep whole .csv file to check if all entries are in utf-8
    arguments: csv file name
    """
    data = pd.read_csv(filename,encoding='iso-8859-1',index_col='ProfileID')      
    for i in data.columns:
        row = 0
        for j in data[i]:  
            row = row + 1    
            print "row = %d" % row
            try:
                str(j).decode('utf-8')
            except UnicodeError:
                print "string is not UTF-8, row = %d,  col = %s" % (row,i)
                print "value is ", j
                raw_input("Press Enter to continue...")
                
def getprofile4modeling(filename,latdim,londim,cutdep=100.0,outfilename='sitegridid.txt'):
    """ from .csv get profile information of those deeper than cut-off depth (cutdep).
    arguments: csv filename, cut-off depth, latdim of model, londim of model
    return: write to '\t' delimited sitegridid.txt file with
            1-D array with [lon,
                            lat,
                            profile ID,
                            grid ID (according to sub2ind),
                            year of measurement,
                            soc averaged D14C (upto cutdep if specified),
                            totsoc (upto cutdep if specified)]
    """
    data = pd.read_csv(filename,encoding='iso-8859-1',index_col='ProfileID')      
    Cave14C = getCweightedD14C2(filename,cutdep=cutdep)
    selprof = Cave14C[Cave14C[:,2]==cutdep,3] # selected profiles ID
    pflonlat = data.loc[selprof][['Lon','Lat','SampleYear']]
    outf = open(outfilename,"w")
    # calculate gridid
    latstep = 180.0 / latdim; lonstep = 360.0 / londim;
    lonmax = 180.0 - lonstep/2.0;
    lonmin = -180.0 + lonstep/2.0;
    latmax = 90.0 - latstep/2.0;
    latmin = -90.0 + latstep/2.0;
    lat = np.arange(latmax,latmin-0.1*latstep,-latstep)
    lon = np.arange(lonmin,lonmax+0.1*lonstep,lonstep)
    for i in pflonlat.index.unique():
        print "i is:",i
        if ~np.any(np.isnan(pflonlat.loc[i:i].iloc[0,:].values.astype(float))):
            ii = np.argmin((float(pflonlat.loc[i:i,'Lat'].iloc[0]) - lat)**2) 
            jj = np.argmin((float(pflonlat.loc[i:i,'Lon'].iloc[0]) - lon)**2)
            gridid = np.ravel_multi_index((ii,jj),dims=(latdim,londim),order='F') + 1
                            # +1 because python starts from 0, this is to be feed to matlab            
            outf.write("%.2f, %.2f, %d, %d, %d, %.2f, %.2f\n" % \
                (float(pflonlat.loc[i:i].iloc[0,0]), \
                float(pflonlat.loc[i:i].iloc[0,1]), i, gridid, int(pflonlat.loc[i:i].iloc[0,2]), \
                Cave14C[Cave14C[:,3]==i,4], Cave14C[Cave14C[:,3]==i,5]))
    outf.close()
    
def getCRUannualdata(var='tmp',
                     CRUpath='C:\\download\\work\\!manuscripts\\' + \
                             'C14_synthesis\\CRUdata\\'):
    """ get annual cru data from 1921 to 2013. 
    output saved to variable.npy. output is top-down 90S-90N
    output:
        temp: scalar 1/10.
        precip: scalar *12/10
    """
    crulist = os.listdir(CRUpath + var)  
    cruyrdata = []
    for i, fn in enumerate(crulist):
        print "current file: ",i
        dum = np.loadtxt(CRUpath + var + '\\' + fn)
        dum[dum <= -990] = np.nan
        step = 360.
        totyr = dum.shape[0]/step/12.
        for yr in range(int(totyr)):
            dumyr = np.nanmean(np.reshape(dum[12*step*yr:step*12*(yr+1),:],(12,step,720),order='C'),axis=0)
            cruyrdata.append(dumyr)
            print "current yr: ",yr
    np.save(CRUpath+'cruannual'+var+'.npy',cruyrdata)
            
def findCRUdataforgrid(lon,lat,year,toprint=False,
                       CRUpath='C:\\download\\work\\!manuscripts\\' + \
                       'C14_synthesis\\AncillaryData\\CRUdata\\'):
    """ find the closest CRU grid to given lon, lat, year
    print the CRU ten-year average var
    params:
        lon, lat: array-like
    return:
        [tmp, pre]: same dimensin as lon,lat
    """
    crudata_tmp = np.load(CRUpath + 'cruannual' + 'tmp' + '.npy')
    crudata_pre = np.load(CRUpath + 'cruannual' + 'pre' + '.npy')
    outtmp = []; outpre = []
    latstep = 0.5; lonstep = 0.5;
    lonmax = 180.0 - lonstep/2.0;
    lonmin = -180.0 + lonstep/2.0;
    latmax = 90.0 - latstep/2.0;
    latmin = -90.0 + latstep/2.0;
    lats = np.arange(latmax,latmin-0.1*latstep,-latstep)
    lons = np.arange(lonmin,lonmax+0.1*lonstep,lonstep)
    
    for lo,la,yr in zip(lon,lat,year):
        if np.any(np.isnan([lo,la,yr])):
            outtmp += np.nan,
            outpre += np.nan,
            continue
        yr = int(yr)
        # get tmp
        yridx = range(yr-1921-10-1+1,yr-1921+1)
        tenyrmean = np.flipud(np.nanmean(crudata_tmp[yridx,:,:],axis=0))
        cols = np.argmin((lo - lons)**2)
        rows = np.argmin((la - lats)**2)
        if toprint:
            print '%s for grid lon=%.5f, lat=%.5f is: %.2f C' % ('temperature',lon,lat,tenyrmean[rows,cols]/10.)
        outtmp += tenyrmean[rows,cols]/10.,        
        # get pre
        yridx = range(yr-1921-10-1+1,yr-1921+1)
        tenyrmean = np.flipud(np.nanmean(crudata_pre[yridx,:,:],axis=0))
        cols = np.argmin((lo - lons)**2)
        rows = np.argmin((la - lats)**2)
        if toprint:
            print '%s for grid lon=%.5f, lat=%.5f is: %.2f mm' % ('precip',lon,lat,tenyrmean[rows,cols]*12./10.)
        outpre += tenyrmean[rows,cols]*12./10.,
    return outtmp, outpre

def fillmissingMATMAP_printonly():
    """ fill missing MATMAP using CRU in synthesis data. using print method
    """
    filename = 'Non_peat_data_synthesis.csv'
    data = pd.read_csv(filename,encoding='iso-8859-1',index_col='ProfileID')
    for i in data[['MAT','MAP']].index.unique():
        if ~np.isnan(i):
            print 'profile ID is: ',i
            if np.any(np.isnan(data.loc[i][['MAT','MAP']].astype(float))):
                try:
                    data.loc[i]['Lon'].unique()                
                    findCRUdataforgrid(data.loc[i]['Lon'].unique().astype(float),
                                       data.loc[i]['Lat'].unique().astype(float),
                                       data.loc[i]['SampleYear'].unique().astype(int),
                                       toprint=True)
                except AttributeError:
                    findCRUdataforgrid(float(data.loc[i]['Lon']),
                                       float(data.loc[i]['Lat']),
                                       int(data.loc[i]['SampleYear']),
                                       toprint=True)
                raw_input('Press enter to continue....')
        else:
            continue

def findSoilOrderforgrid(lon,lat,toprint=False,
                       SOpath='C:\\download\\work\\!manuscripts\\' + \
                       'C14_synthesis\\AncillaryData\\Glb_SoilOrder\\'):
    """ find the closest SoilOrder grid to given lon, lat
    params:
        lon, lat: array-like, same dimension
    return:
        soilorder: {list of str}same dimension as lon/lat
    """
    dic = {159:'Alf', 153:'And', 156:'Ard', 161:'Ent', 160:'Ept', 155:'Ert',
           150:'Gel', 151:'His', 158:'Oll', 154:'Ox',  152:'Spo', 157:'Ult'}
    glbso = np.load(SOpath + 'Glb_soilorder.npy') # [360,720], 90N-90S
    latstep = 0.5; lonstep = 0.5;
    lonmax = 180.0 - lonstep/2.0;
    lonmin = -180.0 + lonstep/2.0;
    latmax = 90.0 - latstep/2.0;
    latmin = -90.0 + latstep/2.0;
    lats = np.arange(latmax,latmin-0.1*latstep,-latstep)
    lons = np.arange(lonmin,lonmax+0.1*lonstep,lonstep)
    so = []    
    for loc_lon, loc_lat in zip(lon,lat):
        cols = np.argmin((loc_lon - lons)**2)
        rows = np.argmin((loc_lat - lats)**2)
        so += glbso[rows, cols],
    return map(lambda x: dic.get(x,np.nan), so)

def findLCforgrid(lon,lat,toprint=False,
                       LCpath='C:\\download\\work\\!manuscripts\\' + \
                       'C14_synthesis\\AncillaryData\\landcoverdata\\'):
    """ find the closest MODIS LC grid to given lon, lat
    params:
        lon, lat: array-like, same dimension
    return:
        soilorder: {list of str}same dimension as lon/lat
    """
    glblc = np.load(LCpath + '14CsynthesisLC.npy') # [360,720], 90N-90S
    latstep = 0.5; lonstep = 0.5;
    lonmax = 180.0 - lonstep/2.0;
    lonmin = -180.0 + lonstep/2.0;
    latmax = 90.0 - latstep/2.0;
    latmin = -90.0 + latstep/2.0;
    lats = np.arange(latmax,latmin-0.1*latstep,-latstep)
    lons = np.arange(lonmin,lonmax+0.1*lonstep,lonstep)
    lc = []    
    for loc_lon, loc_lat in zip(lon,lat):
        cols = np.argmin((loc_lon - lons)**2)
        rows = np.argmin((loc_lat - lats)**2)
        lc += glblc[rows, cols],
    return lc
    
def getHWSD(filename,lon,lat):
    """ get HWSD var for given lon,lat
    lon, lat can be 1-D ndarray, filename is the file of var of interest
    arguments: filename, lon, lat
    return: 1-D ndarray of var of interest
    """
    def findvar(var,lon,lat):        
        cols = np.argmin((lon - lons)**2)
        rows = np.argmin((lat - lats)**2)
        print 'cols:',cols
        print 'rows:',rows
        if ~var.mask[rows,cols]:
            return  var.data[rows,cols]
        else:
            return np.nan
    lon = np.array(lon.astype(float)); lat = np.array(lat.astype(float))
    ncfid = Dataset(filename, 'r')
    nc_attrs, nc_dims, nc_vars = mync.ncdump(ncfid)
    lats = ncfid.variables['lat'][:]
    lons = ncfid.variables['lon'][:]
    var = ncfid.variables[nc_vars[-1]][:]
    try:
        outlen = lon.shape[0]
        out = np.zeros(outlen)
        for n, i in enumerate(zip(lon, lat)):
            print 'n is ',n
            print 'i is ',i
            out[n] = findvar(var,*i)
    except AttributeError:
        print 'singleton dimension encountered, return scalar...'
        out = findvar(var, lon, lat)
    return out
    
def getvarxls(data, var, profid, loc):
    """ get variable from synthesis xlssheet corresponds to proifle ID
    params: 
        data: {DataFrame} whole data sheet read into DF
        var: {string} variable of interest
        profid: {array} 1-D array of profile ID
        loc: 0 for top value, -1 for last value, or ':' for all values
    return: 
        1-D array (same length of profid) of var
    """
    if isinstance(profid, int):
        l = data.loc[profid][var]
    else:
        l = [data.loc[i:i][var] for i in profid if ~np.isnan(i)]
    mat = []
    for i in l:
        try: # if has multi-layers
            if isinstance(loc, str):
                dum = eval('i.iloc[' + loc + '].values')                
                try: 
                    mat.append(list(dum.astype(float)))
                except ValueError: # is unicode
                    mat.append(list(dum))
            else:
                mat.append(i.iloc[loc])
        except AttributeError:
            if isinstance(loc, str):
                try:
                    mat.append([float(i)])
                except ValueError:
                    mat.append([i])
            else:
                mat.append(i)       
    if isinstance(loc, str):
        return np.array(sum(mat, []))
    else:        
        return np.array(mat)    
    
def getslopeFM_lnC(filename, plot=False):
    """
    do linear fit to FM ~ ln(C%), 
    return profile ID and slope in 2-D array.
    """
    data = pd.read_csv(filename,encoding='iso-8859-1',index_col='ProfileID')  
    pid = np.unique(data.index.values.astype(float))   # index of profile start
    pid = pid[np.logical_not(np.isnan(pid))]
    slp = np.zeros((len(pid),2))
    dum = 0
    if plot == True:
        nrow = 4; ncol = 5
        nfig = math.ceil(len(pid)/(nrow*ncol*1.))
    for i in pid:
        print "profile %d ...."%i
        d14C = np.array(data.loc[i:i,'D14C_BulkLayer'].values).astype(float)
        fm = d14C/1000. + 1
        lnC = np.log(np.array(data.loc[i:i,'pct_C'].values).astype(float))
        slp[dum, 0] = i
        notNANs = ~np.any(np.isnan([lnC,fm]),0)
        if len(lnC[notNANs]) != 0:
            slp[dum, 1] = np.polyfit(lnC[notNANs], fm[notNANs], 1)[0]
            lsfit = np.poly1d(np.polyfit(lnC[notNANs], fm[notNANs], 1))
        else:
            slp[dum, 1] = np.polyfit(lnC, fm, 1)[0]
            lsfit = np.poly1d(np.polyfit(lnC, fm, 1))
        yfit = lsfit(lnC)
        # plot line
        if plot == True:
            axesn = int(math.fmod(dum, nrow*ncol))
            if axesn == 0:
                fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=(12,8))
                curfign = 0
            ax = fig.axes[axesn]
            ax.scatter(lnC, fm)
            ax.plot(lnC, yfit, '--')
            if isinstance(data.loc[i,'Site'], unicode):
                ax.set_title(str(int(i)) + ": " + data.loc[i,'Site'])
            else:
                ax.set_title(str(int(i)) + ": " + data.loc[i,'Site'].values[0])
            ax.set_xlabel(r"ln C(%)")
            ax.set_ylabel(r"Fraction Modern")
            curfign = curfign + 1
            if curfign == nrow*ncol or i == len(pid):
                matplotlib.rcParams.update({'font.size': 6}) 
                plt.tight_layout() 
                fig.savefig('./figures/FmvslnC_' + str(i) + '.png')   
                plt.close()
        dum = dum + 1
    return slp
    
def getslopeD14C_cumC(filename, plot=False):
    """
    do linear fit to D14C ~ cumulative C content, 
    return profile ID and slope in 2-D array.
    """
    data = pd.read_csv(filename,encoding='iso-8859-1',index_col='ProfileID')  
    pid = np.unique(data.index.values.astype(float))   # index of profile start
    pid = pid[np.logical_not(np.isnan(pid))]
    slp = np.zeros((len(pid),2))
    dum = 0
    if plot == True:
        nrow = 4; ncol = 5
        nfig = math.ceil(len(pid)/(nrow*ncol)*1.)
    for i in pid:
        print "profile %d ...."%i
        d14C = np.array(data.loc[i:i,'D14C_BulkLayer'].values).astype(float)
        Ccnt = (np.array(data.loc[i:i,'pct_C'].values).astype(float))/100. * \
               (np.array(data.loc[i:i,'BulkDensity'].values).astype(float)) * \
               (np.array(data.loc[i:i,'Layer_bottom'].values).astype(float) - \
                   np.array(data.loc[i:i,'Layer_top'].values).astype(float))
        cumC = np.cumsum(Ccnt)*10. # g/cm2 to kg/m2
        slp[dum, 0] = i
        notNANs = ~np.any(np.isnan([cumC,d14C]),0)
        if len(cumC[notNANs]) != 0:
            slp[dum, 1] = np.polyfit(cumC[notNANs], d14C[notNANs], 1)[0]
            lsfit = np.poly1d(np.polyfit(cumC[notNANs], d14C[notNANs], 1))
        else:
            slp[dum, 1] = np.polyfit(cumC, d14C, 1)[0]
            lsfit = np.poly1d(np.polyfit(cumC, d14C, 1))
        yfit = lsfit(cumC)
        # plot line
        if plot == True:
            axesn = int(math.fmod(dum, nrow*ncol))
            if axesn == 0:
                fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=(12,8))
                curfign = 0
            ax = fig.axes[axesn]
            ax.scatter(cumC, d14C)
            ax.plot(cumC, yfit, '--')
            if isinstance(data.loc[i,'Site'], unicode):
                ax.set_title(str(int(i)) + ": " + data.loc[i,'Site'])
            else:
                ax.set_title(str(int(i)) + ": " + data.loc[i,'Site'].values[0])
            ax.set_xlabel(r"cumulative C$(kg/m^{2})$")
            ax.set_ylabel(r"$\Delta 14C$ ("+ u"\u2030)")
            curfign = curfign + 1
            if curfign == nrow*ncol or i == len(pid):
                matplotlib.rcParams.update({'font.size': 6}) 
                plt.tight_layout() 
                fig.savefig('./figures/D14CvscumC_' + str(i) + '.png')     
                plt.close()
        dum = dum + 1
#        raw_input('press Enter to continue...')        
    return slp    

def plotFM_lnC(filename, pid, leglabel1='Site', leglabel2=None, legloc=3, markersize=25):
    """
    scatter plot profiles in pid in one figure, with linear fit
    arguments: leglabel1 and 2 are the column names you want in the legend
               legloc is location of legend, 3 is lowerleft
    """
    data = pd.read_csv(filename,encoding='iso-8859-1',index_col='ProfileID')  
    fig, ax = plt.subplots(1,1, figsize=(8,6))
    cm = plt.get_cmap('Paired')
    numcolr = len(pid) # no repeat in color
    ax.set_color_cycle([cm(1.*jj/numcolr) for jj in range(numcolr)])
    marker = myplt._markeriter
    for jj,i in enumerate(pid):
        print "profile %d ...."%i
        d14C = np.array(data.loc[i:i,'D14C_BulkLayer'].values).astype(float)
        fm = d14C/1000. + 1
        lnC = np.log(np.array(data.loc[i:i,'pct_C'].values).astype(float))
        notNANs = ~np.any(np.isnan([lnC,fm]),0)
        if len(lnC[notNANs]) != 0:
            lsfit = np.poly1d(np.polyfit(lnC[notNANs], fm[notNANs], 1))
        else:
            lsfit = np.poly1d(np.polyfit(lnC, fm, 1))
        yfit = lsfit(lnC)
        # plot line
        if leglabel2 is not None:
            if not isinstance(data.loc[i:i,leglabel2].values[0], unicode):
                dumlab = ''
            else:
                dumlab = data.loc[i:i,leglabel2].values[0]
            lab = data.loc[i:i,leglabel1].values[0] + '(' + \
                  dumlab + ') (' + \
                  data.loc[i:i,'Layer_top'].values[0] + '-' + \
                  data.loc[i:i,'Layer_bottom'].values[-1] + 'cm)'
        else:
            lab = data.loc[i:i,leglabel1].values[0] + '(' + \
                  data.loc[i:i,'Layer_top'].values[0] + '-' + \
                  data.loc[i:i,'Layer_bottom'].values[-1] + 'cm)'          
        ax.scatter(lnC, fm, marker=marker.next(), color=cm(1.*jj/numcolr), \
                   label=lab.encode("ascii"), s=markersize)
        ax.plot(lnC, yfit, '--')
        ax.set_xlabel('ln (C(%))')
        ax.set_ylabel('Fraction Modern')
    plt.legend(scatterpoints=1,loc=legloc,fontsize=10)
    matplotlib.rcParams.update({'font.size': 12}) 

def plotD14C_cumC(filename, pid, leglabel1='Site', leglabel2=None, legloc=3, markersize=25):
    """
    scatter plot profiles in pid in one figure, with linear fit
    arguments: leglabel1 and 2 are the column names you want in the legend
               legloc is location of legend, 3 is lowerleft
    """
    data = pd.read_csv(filename,encoding='iso-8859-1',index_col='ProfileID')  
    fig, ax = plt.subplots(1,1, figsize=(8,6))
    cm = plt.get_cmap('Paired')
    numcolr = len(pid) # no repeat in color
    ax.set_color_cycle([cm(1.*jj/numcolr) for jj in range(numcolr)])
    marker = myplt._markeriter
    for jj,i in enumerate(pid):
        print "profile %d ...."%i
        d14C = np.array(data.loc[i:i,'D14C_BulkLayer'].values).astype(float)
        Ccnt = (np.array(data.loc[i:i,'pct_C'].values).astype(float))/100. * \
               (np.array(data.loc[i:i,'BulkDensity'].values).astype(float)) * \
               (np.array(data.loc[i:i,'Layer_bottom'].values).astype(float) - \
                   np.array(data.loc[i:i,'Layer_top'].values).astype(float))
        cumC = np.cumsum(Ccnt)*10. # g/cm2 to kg/m2
        notNANs = ~np.any(np.isnan([cumC,d14C]),0)
        if len(cumC[notNANs]) != 0:
            lsfit = np.poly1d(np.polyfit(cumC[notNANs], d14C[notNANs], 1))
        else:
            lsfit = np.poly1d(np.polyfit(cumC, d14C, 1))
        yfit = lsfit(cumC)
        # plot line
        if leglabel2 is not None:
            if not isinstance(data.loc[i:i,leglabel2].values[0], unicode): # is nan
                dumlab = ''
            else:
                dumlab = data.loc[i:i,leglabel2].values[0]
            lab = data.loc[i:i,leglabel1].values[0] + '(' + \
                  dumlab + ') (' + \
                  data.loc[i:i,'Layer_top'].values[0] + '-' + \
                  data.loc[i:i,'Layer_bottom'].values[-1] + 'cm)'
        else:
            lab = data.loc[i:i,leglabel1].values[0] + '(' + \
                  data.loc[i:i,'Layer_top'].values[0] + '-' + \
                  data.loc[i:i,'Layer_bottom'].values[-1] + 'cm)'          
        ax.scatter(cumC, d14C, marker=marker.next(), color=cm(1.*jj/numcolr), \
                   label=lab.encode("ascii"), s=markersize)
        ax.plot(cumC, yfit, '--')
        ax.set_xlabel(r"cumulative C$(kg/m^{2})$")
        ax.set_ylabel(r"$\Delta 14C$ ("+ u"\u2030)")
    plt.legend(scatterpoints=1,loc=legloc,fontsize=10)
    matplotlib.rcParams.update({'font.size': 12}) 
    
def plotDepth_cumC(filename, pid, leglabel1='Site', leglabel2=None, legloc=3, markersize=25):
    """
    plot depth vs. cumulative C content
    """
    data = pd.read_csv(filename,encoding='iso-8859-1',index_col='ProfileID')  
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_axes([0.1,0.1,0.6,0.85])
    cm = plt.get_cmap('Paired')
    numcolr = len(pid) # no repeat in color
    ax.set_color_cycle([cm(1.*jj/numcolr) for jj in range(numcolr)])
    marker = myplt._markeriter
    for jj,i in enumerate(pid):
        print "profile %d ...."%i
        dep = np.array(data.loc[i:i,'Layer_bottom'].values).astype(float)
        Ccnt = (np.array(data.loc[i:i,'pct_C'].values).astype(float))/100. * \
               (np.array(data.loc[i:i,'BulkDensity'].values).astype(float)) * \
               (np.array(data.loc[i:i,'Layer_bottom'].values).astype(float) - \
                   np.array(data.loc[i:i,'Layer_top'].values).astype(float))
        cumC = np.cumsum(Ccnt)*10. # g/cm2 to kg/m2
        # plot line
        if leglabel2 is not None:
            if not isinstance(data.loc[i:i,leglabel2].values[0], unicode):
                dumlab = ''
            else:
                dumlab = data.loc[i:i,leglabel2].values[0]
            lab = data.loc[i:i,leglabel1].values[0] + '(' + \
                  dumlab + ') (' + \
                  data.loc[i:i,'Layer_top'].values[0] + '-' + \
                  data.loc[i:i,'Layer_bottom'].values[-1] + 'cm)'
        else:
            lab = data.loc[i:i,leglabel1].values[0] + '(' + \
                  data.loc[i:i,'Layer_top'].values[0] + '-' + \
                  data.loc[i:i,'Layer_bottom'].values[-1] + 'cm)'          
        ax.plot(cumC, dep, marker=marker.next(), color=cm(1.*jj/numcolr), \
                label=lab.encode("ascii"))
    ax.set_xlabel(r'cumulative C $(kg/m^{2})$')
    ax.set_ylabel('Depth (cm)')
    plt.gca().invert_yaxis()
    plt.legend(scatterpoints=1, loc=legloc, bbox_to_anchor=(1., 0., .5, 1), \
               mode="expand", borderaxespad=0., fontsize=10)
    plt.tight_layout() 
    matplotlib.rcParams.update({'font.size': 12}) 
    
def getprofSOCD14C_interp(filename,profid,cutdep=None):
    """ Get profile of SOC and D14C for one profile id with linear interpolation to cutdep.
    for profiles that are shallower than cutdep, the whole profile is used
    arguments: 1 - csv filename
    return: 3 col nparray, 
            0 - layer depth, bottom
            1 - D14C at each layer
            2 - SOC at each layer for the whole profile (if no cutdep) in kgC/m2 up to cutdep
    """
    data = pd.read_csv(filename,encoding='iso-8859-1',index_col='ProfileID',skiprows=[1])  
    i = profid
    out = []
    # number of layer = 1, calculate whole profile
    if data.loc[i:i,'Layer_top'].shape[0] == 1:
        out.append(np.array(data.loc[i:i,'Layer_bottom']).astype(float))
        totsoc = np.array(data.loc[i:i,'BulkDensity']).astype(float) * \
                 (np.array(data.loc[i:i,'Layer_bottom']).astype(float) - \
                 np.array(data.loc[i:i,'Layer_top']).astype(float)) * \
                 np.array(data.loc[i:i,'pct_C']).astype(float)/100.0
                 #gC/cm2 converted to kgC/m2 for the given depth
        out.append(np.array(data.loc[i:i,'D14C_BulkLayer']).astype(float))
        out.append(totsoc*10.)
    else: # multiple layers
        # needs interpolation
        if cutdep is not None and float(data.loc[i:i,'Layer_bottom'].values[-1]) >= cutdep:
            cutdep = cutdep * 1.  
            layerbot = np.array(data.loc[i:i,'Layer_bottom']).astype(float)
            layertop = np.array(data.loc[i:i,'Layer_top']).astype(float)
            bd = np.array(data.loc[i:i,'BulkDensity']).astype(float)
            d14C = np.array(data.loc[i:i,'D14C_BulkLayer']).astype(float)
            pctC = np.array(data.loc[i:i,'pct_C']).astype(float)  
            idx = layerbot > cutdep
            if cutdep in set(layerbot):
                layerbotinterp = layerbot
            else:
                layerbotinterp = np.hstack((layerbot[idx==False], cutdep, layerbot[idx==True]))
            notNANs = ~np.isnan(d14C)
            f_i = interp1d(layerbot[notNANs], d14C[notNANs]); f_x = extrap1d(f_i); 
            d14Cinterp = f_x(layerbotinterp)
            notNANs = ~np.isnan(bd); 
            if sum(notNANs) == 0:
                bdinterp = bd
            else:
                f_i = interp1d(layerbot[notNANs], bd[notNANs]); f_x = extrap1d(f_i); 
                bdinterp = f_x(layerbotinterp)
            notNANs = ~np.isnan(pctC)
            if sum(notNANs) == 0:
                pctCinterp = pctC
            else:
                f_i = interp1d(layerbot[notNANs], pctC[notNANs]); f_x = extrap1d(f_i); 
                pctCinterp = f_x(layerbotinterp)
            out.append(layerbotinterp[layerbotinterp <= cutdep])
            totsoc = bdinterp[layerbotinterp <= cutdep] * \
                     (layerbotinterp[layerbotinterp <= cutdep] - \
                      np.hstack((layertop[0],layerbotinterp[layerbotinterp < cutdep]))) * \
                     pctCinterp[layerbotinterp <= cutdep]/100.0
            out.append(d14Cinterp[layerbotinterp <= cutdep])
            out.append(totsoc*10.) #gC/cm2 converted to kgC/m2 for the given depth
            print 'interpolate to cutdep...'
#                fig = plt.figure()
#                ax = fig.add_axes([0.05,0.05,0.9,0.9])
#                ax.scatter(d14C, layerbot, marker=marker.next())
#                ax.scatter(d14Cinterp, layerbotinterp, marker=marker.next())
#                plt.gca().invert_yaxis()
#                raw_input('press Enter to continue...') 
#                plt.close()
        # calculate the whole profile
        else: 
            # no missing value
            if np.all(~np.isnan(
                      data.loc[i:i,['D14C_BulkLayer','BulkDensity','pct_C']].astype(float))):
                out.append(np.array(data.loc[i:i,'Layer_bottom'].values).astype(float))
                totsoc = np.array(data.loc[i:i,'BulkDensity'].values).astype(float) * \
                         (np.array(data.loc[i:i,'Layer_bottom'].values).astype(float) - \
                         np.array(data.loc[i:i,'Layer_top'].values).astype(float)) * \
                         np.array(data.loc[i:i,'pct_C'].values).astype(float)/100.0
                out.append(np.array(data.loc[i:i,'D14C_BulkLayer'].values).astype(float))
                out.append(totsoc*10.) #gC/cm2 converted to kgC/m2 for the whole profile
            # linear interpolate missing values at given layer
            else:   
                layerbot = np.array(data.loc[i:i,'Layer_bottom']).astype(float)
                layertop = np.array(data.loc[i:i,'Layer_top']).astype(float)
                bd = np.array(data.loc[i:i,'BulkDensity']).astype(float)
                d14C = np.array(data.loc[i:i,'D14C_BulkLayer']).astype(float)
                pctC = np.array(data.loc[i:i,'pct_C']).astype(float)  
                notNANs = ~np.isnan(d14C)                    
                f_i = interp1d(layerbot[notNANs], d14C[notNANs]); f_x = extrap1d(f_i); 
                d14Cinterp = f_x(layerbot)
                notNANs = ~np.isnan(bd)
                if sum(notNANs) == 0:
                    bdinterp = bd
                else:
                    f_i = interp1d(layerbot[notNANs], bd[notNANs]); f_x = extrap1d(f_i); 
                    bdinterp = f_x(layerbot)
                notNANs = ~np.isnan(pctC)
                if sum(notNANs) == 0:
                    pctCinterp = pctC
                else:
                    f_i = interp1d(layerbot[notNANs], pctC[notNANs]); f_x = extrap1d(f_i); 
                    pctCinterp = f_x(layerbot)
                out.append(layerbot)
                out.append(d14Cinterp[layerbot <= layerbot[-1]])
                totsoc = bdinterp * (layerbot - layerbot) * pctCinterp/100.0
                out.append(totsoc*10.) #gC/cm2 converted to kgC/m2 for the whole profile
                print 'interpolate for missing value...'
#                    fig = plt.figure()
#                    ax = fig.add_axes([0.05,0.05,0.9,0.9])
#                    ax.scatter(d14C, layerbot, marker=marker.next())
#                    ax.scatter(d14Cinterp, layerbot, marker=marker.next())
#                    raw_input('press Enter to continue...') 
#                    plt.close()
    return np.array(out).T

def getDDelta14C(sampleyr, d14C):
    '''calculate deltaDelta14C, defined as D14Catm - D14Csoil
    param: 
        {np array, float | int} input scalar/array of sampleyr and the same size array of d14Csoil
    outputs: scalar/array of dd14C
    '''
    import scipy.io
    atm14C = scipy.io.loadmat('atmD14C_50kBP-2012.mat')
    atmD14C = atm14C['atmD14C']
    if sampleyr.ndim == 0:
        return atmD14C[-(2013-sampleyr),1] - d14C
    dd14C = np.zeros((d14C.shape[0]))
    for i in range(dd14C.shape[0]):
        if ~np.isnan(sampleyr[i]):
            if sampleyr[i] < 2013:
                dd14C[i] = atmD14C[-int(2013-sampleyr[i]),1] - d14C[i]
            else:
                dd14C[i] = atmD14C[-1,1] - d14C[i]
        else:
            dd14C[i] = np.nan
    return dd14C

def to_number(s):
    '''
    convert series that are numerical to floats
    '''
    try:
        s1 = s.astype(float)
        return s1
    except ValueError:
        return s
        
def get_biome_ave(biome, interpdf, cols):
    '''
    Get biome averaged value using data matrix that is interpolated to fixed increments
    parameters:
        biome    : {int} biome index (1-10)
        interpdf : {dataframe} df that is interperlated to fixed increments
        cols     : {str | list of str} cols names that you want to average on
    return:
        outave : list of averaged profile for each biome
        outstd : list of std for each biome
    '''
    grouped = interpdf.groupby('VegTypeCode_Local').get_group(biome)
    tmp = grouped.groupby('Layer_depth_incre')
    return {'Layer_depth_incre':tmp[cols].mean().index, 
            'mean':tmp[cols].mean().values,
            'std':tmp[cols].std().values,
            'sem':tmp[cols].sem().values}

def get_climzone_ave(kind, clim_code, interpdf, cols):
    '''
    Get climate-zone averaged value using data matrix that is interpolated to fixed increments
    parameters:
        kind     : {str}, 'whittaker' or 'self-defined'
        biome    : {int} clim_code index (1-9)
        interpdf : {dataframe} df that is interperlated to fixed increments
        cols     : {str | list of str} cols names that you want to average on
    return:
        outave : list of averaged profile for each biome
        outstd : list of std for each biome
    '''
    if kind == 'whittaker':
        grouped = interpdf.groupby('clim_code_whittaker').get_group(clim_code)
    elif kind == 'self-defined':
        grouped = interpdf.groupby('clim_code').get_group(clim_code)        
    tmp = grouped.groupby('Layer_depth_incre')
    return {'Layer_depth_incre':tmp[cols].mean().index, 
            'mean':tmp[cols].mean().values,
            'std':tmp[cols].std().values,
            'sem':tmp[cols].sem().values}


def get_soilorder_ave(soilorder, interpdf, cols):
    '''
    Get soil order averaged value using data matrix that is interpolated to fixed increments
    parameters:
        soilorder: {int|str} soil order index (1-12 or the string)
        interpdf : {dataframe} df that is interperlated to fixed increments
        cols     : {str | list of str} cols names that you want to average on. 
                   col should be numerical variable
    return:
        outave : list of averaged profile for each soil order
        outstd : list of std for each soil order
    '''
    interpdf[cols] = interpdf[cols].astype(float)
    grouped = interpdf.groupby('SoilOrder_LEN_USDA').get_group(soilorder)
    tmp = grouped.groupby('Layer_depth_incre')
    return {'Layer_depth_incre':tmp[cols].mean().index, 
            'mean':tmp[cols].mean().values,
            'std':tmp[cols].std().values}
            
            
def getD14C_cumC(filename):
    ''' 
    Add cumC content (kgC/m2) to the dataframe, including only columns that have complete 
    pctC information
    '''
    data = pd.read_csv(filename,encoding='iso-8859-1',index_col='ProfileID',skiprows=[1])
    pid = np.unique(data.index.values)
    data['cumC'] = np.repeat(np.nan, data.shape[0])
    for i in pid:
        d14C = np.array(data.loc[i:i,'D14C_BulkLayer'].values).astype(float)
        Ccnt = (np.array(data.loc[i:i,'pct_C'].values).astype(float))/100. * \
               (np.array(data.loc[i:i,'BulkDensity'].values).astype(float)) * \
               (np.array(data.loc[i:i,'Layer_bottom'].values).astype(float) - \
                   np.array(data.loc[i:i,'Layer_top'].values).astype(float))
        if Ccnt.shape[0] <=3:
            continue
        if Ccnt.shape[0] < data.loc[i:i,'Layer_bottom'].shape[0]:
            continue # skip missing value profile for now
        cumC = np.cumsum(Ccnt)*10. # g/cm2 to kg/m2
        data.loc[i:i, 'cumC'] = cumC
    return data
        
def prep_increment(newdf, colname, inc=1):
    '''
        pass in csv data, interpolate to increment (inc=1cm) by assigning all 
        increments within a layer the same values from that layer. 
    param: 
        newdf   : {dataframe}, cols from csv, makes sure you pass in depth 
                  (normalized depth, top, bot)
        colname : {list}, col names of the newdf passed in
        inc     : self-defined increment to interpolate to
    return:
        outdf   : {dataframe} interpolated newdf, with the continuous depth appended
                     col 0  : continuous depth, 'Layer_depth_incre'
                     col 1  : profid
                     col 2-3: depth
                     col 4+ : data
    '''
    uniqprofid = newdf.index.unique()
    out = []
    idx = []
    for idd in uniqprofid:
        lay_bot = np.r_[np.ceil(newdf.loc[idd:idd,'Layer_top_norm'].values[0]), 
                        np.ceil(newdf.loc[idd:idd,'Layer_bottom_norm'].values)]
        cont_lay_bot = np.arange(lay_bot[0], lay_bot[-1]+1)        
        tmpdata = newdf.loc[idd:idd, colname].values
        diff = np.diff(lay_bot)
        repn = diff/inc; repn = repn.astype(int)
        newdat = np.repeat(tmpdata, repn, axis=0)
        out += list(np.c_[cont_lay_bot,np.vstack((newdat,newdat[-1]))])
        idx += list(np.repeat(idd, len(cont_lay_bot))) # profid
    out = np.asarray(out)
    outdf = pd.DataFrame(data=out, index=idx, columns=['Layer_depth_incre']+colname) 
    converted = outdf.apply(lambda col : to_number(col) , axis = 1)  
    return converted

def prep_lnrinterp(data, min_layer=3):
    '''
    Given a profie id, if BD and pctC data are more than 3 for each, linear interp at 1cm
    interval of whole profile. return new dataframe.
    In contrast to 'prep_increment' which assumes homogeneous for each layer.
    params: 
        {dataframe} data
    return:
        {dataframe} new df with BD, pctC, D14C and depth linear interpolated at 1cm increment
    '''   
    pid = data.index.unique()   # index of profile start
    out = pd.DataFrame(columns=['Layer_depth_incre','D14C_BulkLayer','BulkDensity',
                                'pct_C','Cmass','Lon','Lat','MAT','MAP','VegTypeCode_Local',
                                'SoilOrder_LEN_USDA'])
    statprop = pd.DataFrame(columns=['D14C_R2','D14C_RMSE','D14C_pctERR',
                                     'BD_R2','BD_RMSE','BD_pctERR',
                                     'pctC_R2','pctC_RMSE','pctC_pctERR',
                                     'Cmass_R2','Cmass_RMSE','Cmass_pctERR'],index=pid)

    marker = myplt._markeriter
    for n,i in enumerate(pid):
        print 'profile is :',i
        if data.loc[i:i,'Layer_top'].shape[0] < min_layer: # number of layer < 3, skip
            print 'valid layer number less than %d, skip profile...'%min_layer
            continue            
        else: # multiple layers
            print 'interpolate ...'
            incre = np.arange(np.round(data.loc[i:i,'Layer_top_norm'].values[0]),
                              np.round(data.loc[i:i,'Layer_bottom_norm'].values[-1])+1, 1)
            tmpdf = pd.DataFrame(columns=['Layer_depth_incre','D14C_BulkLayer','BulkDensity',
                                'pct_C','Lon','Lat','MAT','MAP','VegTypeCode_Local',
                                'SoilOrder_LEN_USDA'],index=np.repeat(i,incre.shape[0]))   
            tmpdf.loc[i:i,'Layer_depth_incre'] = incre
            layerbotori = np.array(data.loc[i:i,'Layer_bottom_norm']).astype(float)
            layertop = np.array(data.loc[i:i,'Layer_top_norm']).astype(float)
            layerbot = np.mean(np.c_[layertop, layerbotori], axis=1)
            bd = np.array(data.loc[i:i,'BulkDensity']).astype(float)
            d14C = np.array(data.loc[i:i,'D14C_BulkLayer']).astype(float)
            pctC = np.array(data.loc[i:i,'pct_C']).astype(float)
            Cmass = np.array(data.loc[i:i,'Cmass']).astype(float)
            
            notNANs = ~np.isnan(d14C)
            if sum(notNANs) >= 3: # otherwise not enough value! skip this variable
                f_i = interp1d(layerbot[notNANs], d14C[notNANs]); f_x = extrap1d_constby(f_i); 
                tmpdf.loc[i:i,'D14C_BulkLayer'] = f_x(incre)
                yhat = f_x(layerbot)
                statprop.loc[i:i,'D14C_R2'] = mysm.cal_R2(d14C[notNANs], yhat[notNANs])
                statprop.loc[i:i,'D14C_RMSE'] = mysm.cal_RMSE(d14C[notNANs], yhat[notNANs])
                statprop.loc[i:i,'D14C_pctERR'] = mysm.cal_pctERR(d14C[notNANs], yhat[notNANs])
                
            notNANs = ~np.isnan(bd)
            if sum(notNANs) >= 3: # otherwise not enough value! skip this variable
                f_i = interp1d(layerbot[notNANs], bd[notNANs]); f_x = extrap1d_constby(f_i); 
                tmpdf.loc[i:i,'BulkDensity'] = f_x(incre)
                yhat = f_x(layerbot)
                statprop.loc[i:i,'BD_R2'] = mysm.cal_R2(bd[notNANs], yhat[notNANs])
                statprop.loc[i:i,'BD_RMSE'] = mysm.cal_RMSE(bd[notNANs], yhat[notNANs])
                statprop.loc[i:i,'BD_pctERR'] = mysm.cal_pctERR(bd[notNANs], yhat[notNANs])
            
            notNANs = ~np.isnan(pctC)
            if sum(notNANs) >= 3: # otherwise not enough value! skip this variable
                f_i = interp1d(layerbot[notNANs], pctC[notNANs]); f_x = extrap1d_constby(f_i); 
                tmpdf.loc[i:i,'pct_C'] = f_x(incre)
                yhat = f_x(layerbot)
                statprop.loc[i:i,'pctC_R2'] = mysm.cal_R2(pctC[notNANs], yhat[notNANs])
                statprop.loc[i:i,'pctC_RMSE'] = mysm.cal_RMSE(pctC[notNANs], yhat[notNANs])
                statprop.loc[i:i,'pctC_pctERR'] = mysm.cal_pctERR(pctC[notNANs], yhat[notNANs])

            notNANs = ~np.isnan(Cmass)
            if sum(notNANs) >= 3: # otherwise not enough value! skip this variable
                f_i = interp1d(layerbot[notNANs], Cmass[notNANs]); f_x = extrap1d_constby(f_i); 
                tmpdf.loc[i:i,'Cmass'] = f_x(incre)
                yhat = f_x(layerbot)
                statprop.loc[i:i,'Cmass_R2'] = mysm.cal_R2(Cmass[notNANs], yhat[notNANs])
                statprop.loc[i:i,'Cmass_RMSE'] = mysm.cal_RMSE(Cmass[notNANs], yhat[notNANs])
                statprop.loc[i:i,'Cmass_pctERR'] = mysm.cal_pctERR(Cmass[notNANs], yhat[notNANs])
                
            rep = incre.shape[0]
            for name in ['Lon','Lat','MAT','MAP','VegTypeCode_Local','SoilOrder_LEN_USDA']:
                tmpdf.loc[i:i, name] =  np.repeat(data.loc[i:i, name].iloc[0], rep)
            out = pd.concat([out, tmpdf],axis=0) 
#            fig, ax = plt.figure()
#            ax = fig.add_axes([0.05,0.05,0.9,0.9])
#            ax.scatter(d14C, layerbot, marker=marker.next())
#            ax.scatter(d14Cinterp, layerbotinterp, marker=marker.next())
#            plt.gca().invert_yaxis()
#            raw_input('press Enter to continue...') 
            plt.close()
    out = out.apply(lambda col : to_number(col) , axis = 1)  
    return out, statprop

def prep_expinterp(data, min_layer=3):
    '''
    Given a profie id, if BD and pctC data are more than 3 for each, exp interp at 1cm
    interval of whole profile. return new dataframe.
    In contrast to 'prep_increment' which assumes homogeneous for each layer.
    params: 
        {dataframe} data
    return:
        {dataframe} new df with BD, pctC, D14C and depth exp interpolated at 1cm increment
    '''   
    def expfunc(x, K, I):
        '''
        Z* function. not forcing through Csurf. 
        @params: K, I
            z*      = -1/K
            csurf   = I
        '''    
        return I*np.exp(K*x)    

    pid = data.index.unique()   # index of profile start
    out = pd.DataFrame(columns=['Layer_depth_incre','D14C_BulkLayer','BulkDensity',
                                'pct_C','Lon','Lat','MAT','MAP','VegTypeCode_Local',
                                'SoilOrder_LEN_USDA'])
    statprop = pd.DataFrame(columns=['D14C_R2','D14C_RMSE','D14C_pctERR','D14C_zstar','D14C_surf',
                                     'BD_R2','BD_RMSE','BD_pctERR','BD_zstar','BD_surf',
                                     'pctC_R2','pctC_RMSE','pctC_pctERR','pctC_zstar','pctC_surf',
                                     ''],index=pid)
    marker = myplt._markeriter
    for n,i in enumerate(pid):
        print 'profile is :',i
        if data.loc[i:i,'Layer_top'].shape[0] < min_layer: # number of layer < 3, skip
            print 'valid layer number less than 3, skip profile...'
            continue            
        else: # multiple layers
            print 'interpolate ...'
            incre = np.arange(data.loc[i:i,'Layer_top_norm'].values[0],
                              data.loc[i:i,'Layer_bottom_norm'].values[-1]+1, 1)
            tmpdf = pd.DataFrame(columns=['Layer_depth_incre','D14C_BulkLayer','BulkDensity',
                                'pct_C','Lon','Lat','MAT','MAP','VegTypeCode_Local',
                                'SoilOrder_LEN_USDA'],index=np.repeat(i,incre.shape[0])) 
            tmpdf.loc[i:i,'Layer_depth_incre'] = incre
            layerbotori = np.array(data.loc[i:i,'Layer_bottom_norm']).astype(float)
            layertop = np.array(data.loc[i:i,'Layer_top_norm']).astype(float)
            layerbot = np.mean(np.c_[layertop, layerbotori], axis=1) # use mid-point depth
            bd = np.array(data.loc[i:i,'BulkDensity']).astype(float)
            d14C = np.array(data.loc[i:i,'D14C_BulkLayer']).astype(float)
            pctC = np.array(data.loc[i:i,'pct_C']).astype(float)
            
            notNANs = ~np.isnan(d14C)
            if sum(notNANs) >= 3: # otherwise not enough value! skip this variable
                popt, pcov = curve_fit(expfunc, layerbot[notNANs], d14C[notNANs], maxfev=10000,
                                       p0=(-0.01,d14C[notNANs][0])) 
                tmpdf.loc[i:i,'D14C_BulkLayer'] = expfunc(incre, *popt)
                yhat = expfunc(layerbot[notNANs], *popt)
                statprop.loc[i:i,'D14C_R2'] = mysm.cal_R2(d14C[notNANs], yhat)
                statprop.loc[i:i,'D14C_RMSE'] = mysm.cal_RMSE(d14C[notNANs], yhat)
                statprop.loc[i:i,'D14C_pctERR'] = mysm.cal_pctERR(d14C[notNANs], yhat)
                statprop.loc[i:i,'D14C_zstar'] = -1./popt[1]
                statprop.loc[i:i,'D14C_surf'] = np.exp(popt[0])
                
            notNANs = ~np.isnan(bd)
            if sum(notNANs) >= 3: # otherwise not enough value! skip this variable
                popt, pcov = curve_fit(expfunc, layerbot[notNANs], bd[notNANs], maxfev=10000,
                                       p0=(-0.01,bd[notNANs][0])) 
                tmpdf.loc[i:i,'BulkDensity'] = expfunc(incre, *popt)
                yhat = expfunc(layerbot[notNANs], *popt)
                statprop.loc[i:i,'BD_R2'] = mysm.cal_R2(bd[notNANs], yhat)
                statprop.loc[i:i,'BD_RMSE'] = mysm.cal_RMSE(bd[notNANs], yhat)
                statprop.loc[i:i,'BD_pctERR'] = mysm.cal_pctERR(bd[notNANs], yhat)
                statprop.loc[i:i,'BD_zstar'] = -1./popt[1]
                statprop.loc[i:i,'BD_surf'] = np.exp(popt[0])
 
            notNANs = ~np.isnan(pctC)
            if sum(notNANs) >= 3: # otherwise not enough value! skip this variable
                popt, pcov = curve_fit(expfunc, layerbot[notNANs], pctC[notNANs], maxfev=10000,
                                       p0=(-0.01,pctC[notNANs][0])) 
                tmpdf.loc[i:i,'pct_C'] = expfunc(incre, *popt)
                yhat = expfunc(layerbot[notNANs], *popt)
                statprop.loc[i:i,'pctC_R2'] = mysm.cal_R2(pctC[notNANs], yhat)
                statprop.loc[i:i,'pctC_RMSE'] = mysm.cal_RMSE(pctC[notNANs], yhat)
                statprop.loc[i:i,'pctC_pctERR'] = mysm.cal_pctERR(pctC[notNANs], yhat)
                statprop.loc[i:i,'pctC_zstar'] = -1./popt[1]
                statprop.loc[i:i,'pctC_surf'] = np.exp(popt[0])
                
            rep = incre.shape[0]
            for name in ['Lon','Lat','MAT','MAP','VegTypeCode_Local','SoilOrder_LEN_USDA']:
                tmpdf.loc[i:i, name] =  np.repeat(data.loc[i:i, name].iloc[0], rep)
            out = pd.concat([out, tmpdf],axis=0) 
#            fig = plt.figure()
#            ax = fig.add_axes([0.05,0.05,0.9,0.9])
#            ax.scatter(d14C, layerbot, marker=marker.next())
#            ax.scatter(d14C[notNANs], layerbot[notNANs], marker=marker.next())
#            plt.gca().invert_yaxis()
#            raw_input('press Enter to continue...') 
            plt.close()
    out = out.apply(lambda col : to_number(col) , axis = 1)  
    return out, statprop 

def prep_polyinterp(data, deg, min_layer=3):
    '''
    Given a profie id, if BD and pctC data are more than 3 for each, do 2nd order polynomial
    interp at 1cm interval of whole profile. return new dataframe.
    In contrast to 'prep_increment' which assumes homogeneous for each layer.
    params: 
        data        : {dataf frame}
        deg         : {int} degree, 2 or 3
        min_layer   : {int} minimum number of layers required. > min_layer
    return:
        out         : {dataframe} new df with BD, pctC, D14C and depth exp 
                      interpolated at 1cm increment
        statprop    : {dataframe} fitting properties/statistics
    '''   
    def polyfunc(x, p, deg):
        '''
        @params: p, p(x) = p[0] * x**deg + ... + p[deg]; deg = 2 or 3
        '''    
        n = x.shape[0]
        if deg == 2:
            return np.sum(np.c_[x**2, x, np.tile(1.,(n,))]*p,axis=1)
        elif deg == 3:
            return np.sum(np.c_[x**3, x**2, x, np.tile(1.,(n,))]*p,axis=1)
        
    pid = data.index.unique()   # index of profile start
    out = pd.DataFrame(columns=['Layer_depth_incre','D14C_BulkLayer','BulkDensity',
                                'pct_C','Lon','Lat','MAT','MAP','VegTypeCode_Local',
                                'SoilOrder_LEN_USDA'])
    if deg == 2:
        statprop = pd.DataFrame(columns=['D14C_R2','D14C_RMSE','D14C_pctERR','D14C_p0','D14C_p1','D14C_p2',
                                         'BD_R2','BD_RMSE','BD_pctERR','BD_p0','BD_p1','BD_p2',
                                         'pctC_R2','pctC_RMSE','pctC_pctERR','pctC_p0','pctC_p1','pctC_p2'],
                                         index=pid)
    elif deg == 3:
        statprop = pd.DataFrame(columns=['D14C_R2','D14C_RMSE','D14C_pctERR','D14C_p0','D14C_p1','D14C_p2','D14C_p3',
                                         'BD_R2','BD_RMSE','BD_pctERR','BD_p0','BD_p1','BD_p2','BD_p3',
                                         'pctC_R2','pctC_RMSE','pctC_pctERR','pctC_p0','pctC_p1','pctC_p2','pctC_p3'],
                                         index=pid)
        
    marker = myplt._markeriter
    for n,i in enumerate(pid):
        print 'profile is :',i
        if data.loc[i:i,'Layer_top'].shape[0] < min_layer: # number of layer < 3, skip
            print 'valid layer number less than %d, skip profile...'%min_layer
            continue
        else: # multiple layers
            print 'interpolate ...'
            incre = np.arange(np.round(data.loc[i:i,'Layer_top_norm'].values[0]),
                              np.round(data.loc[i:i,'Layer_bottom_norm'].values[-1]+1), 1)
            tmpdf = pd.DataFrame(columns=['Layer_depth_incre','D14C_BulkLayer','BulkDensity',
                                'pct_C','Lon','Lat','MAT','MAP','VegTypeCode_Local',
                                'SoilOrder_LEN_USDA'],index=np.repeat(i,incre.shape[0])) 
            tmpdf.loc[i:i,'Layer_depth_incre'] = incre
            layerbotori = np.array(data.loc[i:i,'Layer_bottom_norm']).astype(float)
            layertop = np.array(data.loc[i:i,'Layer_top_norm']).astype(float)
            layerbot = np.mean(np.c_[layertop, layerbotori], axis=1) # use mid-point depth
            bd = np.array(data.loc[i:i,'BulkDensity']).astype(float)
            d14C = np.array(data.loc[i:i,'D14C_BulkLayer']).astype(float)
            pctC = np.array(data.loc[i:i,'pct_C']).astype(float)
            
            notNANs = ~np.isnan(d14C)
            if sum(notNANs) >= 3: # otherwise not enough value! skip this variable
                p = np.polyfit(layerbot[notNANs], d14C[notNANs], deg) 
                tmpdf.loc[i:i,'D14C_BulkLayer'] = polyfunc(incre, p, deg)
                yhat = polyfunc(layerbot[notNANs], p, deg)
                statprop.loc[i:i,'D14C_R2'] = mysm.cal_R2(d14C[notNANs], yhat)
                statprop.loc[i:i,'D14C_RMSE'] = mysm.cal_RMSE(d14C[notNANs], yhat)
                statprop.loc[i:i,'D14C_pctERR'] = mysm.cal_pctERR(d14C[notNANs], yhat)
                if deg == 2:
                    statprop.loc[i:i,['D14C_p0','D14C_p1','D14C_p2']] = p
                elif deg == 3:
                    statprop.loc[i:i,['D14C_p0','D14C_p1','D14C_p2','D14C_p3']] = p
                    
            notNANs = ~np.isnan(bd)
            if sum(notNANs) >= 3: # otherwise not enough value! skip this variable
                p = np.polyfit(layerbot[notNANs], bd[notNANs], deg) 
                tmpdf.loc[i:i,'BulkDensity'] = polyfunc(incre, p, deg)
                yhat = polyfunc(layerbot[notNANs], p, deg)
                statprop.loc[i:i,'BD_R2'] = mysm.cal_R2(bd[notNANs], yhat)
                statprop.loc[i:i,'BD_RMSE'] = mysm.cal_RMSE(bd[notNANs], yhat)
                statprop.loc[i:i,'BD_pctERR'] = mysm.cal_pctERR(bd[notNANs], yhat)
                if deg == 2:
                    statprop.loc[i:i,['BD_p0','BD_p1','BD_p2']] = p
                elif deg == 3:
                    statprop.loc[i:i,['BD_p0','BD_p1','BD_p2','BD_p3']] = p
                    
            notNANs = ~np.isnan(pctC)
            if sum(notNANs) >= 3: # otherwise not enough value! skip this variable
                p = np.polyfit(layerbot[notNANs], pctC[notNANs], deg) 
                tmpdf.loc[i:i,'pct_C'] = polyfunc(incre, p, deg)
                yhat = polyfunc(layerbot[notNANs], p, deg)
                statprop.loc[i:i,'pctC_R2'] = mysm.cal_R2(pctC[notNANs], yhat)
                statprop.loc[i:i,'pctC_RMSE'] = mysm.cal_RMSE(pctC[notNANs], yhat)
                statprop.loc[i:i,'pctC_pctERR'] = mysm.cal_pctERR(pctC[notNANs], yhat)
                if deg == 2:
                    statprop.loc[i:i,['pctC_p0','pctC_p1','pctC_p2']] = p
                elif deg == 3:
                    statprop.loc[i:i,['pctC_p0','pctC_p1','pctC_p2','pctC_p3']] = p
                    
            rep = incre.shape[0]
            for name in ['Lon','Lat','MAT','MAP','VegTypeCode_Local','SoilOrder_LEN_USDA']:
                tmpdf.loc[i:i, name] = np.repeat(data.loc[i:i, name].iloc[0], rep)
            out = pd.concat([out, tmpdf],axis=0) 
#            fig = plt.figure()
#            ax = fig.add_axes([0.05,0.05,0.9,0.9])
#            ax.scatter(d14C, layerbot, marker=marker.next())
#            ax.scatter(d14C[notNANs], layerbot[notNANs], marker=marker.next())
#            plt.gca().invert_yaxis()
#            raw_input('press Enter to continue...') 
            plt.close()
    out = out.apply(lambda col : to_number(col) , axis = 1)  
    return out, statprop 
    
def get_originalLayerdepth(df):
    ''' pass in normalized depth, get back original depth (start from 0)
    '''
    for idd in df.index.unique():
        print idd
        if df.loc[idd:idd, 'Layer_top_norm'].values[0] < 0:
            offset = np.abs(df.loc[idd:idd, 'Layer_top_norm'].values[0])
            df.loc[idd:idd, 'Layer_top'] = df.loc[idd:idd, 'Layer_top_norm'] + offset
            df.loc[idd:idd, 'Layer_bottom'] = df.loc[idd:idd, 'Layer_bottom_norm'] + offset
    return df

def cal_pfNumber(top, bottom, norm=False):
    ''' calculate number of profiles that spans the given
        top and bottom depth
        params:
            top/bottom: {float} top/bottom you want to constrain
            norm: {bool} whether to consider normed depth
        return:
            newdf: with the selected profiles
    '''
    filename = 'Non_peat_data_synthesis.csv'
    data = pd.read_csv(filename,encoding='iso-8859-1',index_col='ProfileID', skiprows=[1]) 
    pid = data.index.unique() # index of profile start    
    Is_complete = []
    for pf in pid:
        if norm is True:
            if isinstance(data.loc[pf:pf, 'Layer_top_norm'], float):
                Is_complete.append((data.loc[pf:pf,'Layer_top_norm'].values<=top) \
                                & (data.loc[pf:pf, 'Layer_bottom_norm'].values>=bottom))
            else:
                Is_complete.append((data.loc[pf:pf,'Layer_top_norm'].values[0]<=top) \
                                & (data.loc[pf:pf, 'Layer_bottom_norm'].values[-1]>=bottom))
        else:
            if isinstance(data.loc[pf:pf, 'Layer_top'], float):
                Is_complete.append((data.loc[pf:pf,'Layer_top'].values<=top) \
                                & (data.loc[pf:pf, 'Layer_bottom'].values>=bottom))
            else:
                Is_complete.append((data.loc[pf:pf,'Layer_top'].values[0]<=top) \
                                & (data.loc[pf:pf, 'Layer_bottom'].values[-1]>=bottom))            
            print pf
    print 'total number of selected profiles: ', sum(Is_complete)    
    return data[data.index.isin(pid[np.array(Is_complete)])]

def fill_LC_SoilOrder(df):
    '''
    Fill the MISSING LC/SoilOrder using gridded dataset. 
    params:
        df: {DataFrame} the data with .csv directly read in
    return:
        df: {DataFrame} with MISSING LC/SoilOrder filled
    '''
    profid = df.index.unique() # index of profile start
    lon = getvarxls(df,'Lon',profid,':')
    lat = getvarxls(df,'Lat',profid,':')   
    so = findSoilOrderforgrid(lon, lat)
    lc = findLCforgrid(lon, lat)
    newlc = []
    for n,i in enumerate(df.VegTypeCode_Local):
        if np.isnan(i):
            newlc += lc[n],
        else:
            newlc += i,
    df.VegTypeCode_Local = newlc
    newso = []
    for n,i in enumerate(df.SoilOrder_LEN_USDA):
        if isinstance(i, float) and np.isnan(i):
            newso += unicode(so[n]),
        else:
            newso += i,
    df.SoilOrder_LEN_USDA = newso
    df[['VegTypeCode_Local','SoilOrder_LEN_USDA']].to_csv('lcso.csv')

def fill_MATMAP(df):
    '''
    Fill the MISSING MAT/MAP using gridded dataset. 
    params:
        df: {DataFrame} the data with .csv directly read in
    return:
        df: {DataFrame} with MISSING MAT/MAP filled and save to csv
    '''
    profid = df.index.unique() # index of profile start
    lon = getvarxls(df,'Lon',profid,':')
    lat = getvarxls(df,'Lat',profid,':')   
    # fill MAT/MAP, slow, only done once
    yr = getvarxls(df, 'SampleYear', profid, ':')    
    matt, mapp = findCRUdataforgrid(lon,lat,yr)        
    newt = []; newp = []
    for n,(i,j) in enumerate(zip(df.MAT, df.MAP)):
        if isinstance(i, float) and np.isnan(i):
            newt += matt[n],
        else:
            newt += i,
        if isinstance(j, float) and np.isnan(j):
            newp += mapp[n],
        else:
            newp += j,
    df.MAT = newt
    df.MAP = newp
    df[['MAT','MAP']].to_csv('matmap.csv')
    
if __name__ == "__main__":
    print "nothing to do"
