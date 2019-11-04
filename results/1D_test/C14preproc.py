# -*- coding: utf-8 -*-
"""
Module to preprocess D14C data

Created on Thu Feb 12 20:28:02 2015

@Original author: Yujie
@Original source: D14Cpreprocess.py
@Modified by: Shijie Shu
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
