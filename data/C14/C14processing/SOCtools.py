"""

* File Name : SOCtools.py

* Purpose : Lib for processing SOC data

* Creation Date : 18-07-2017

* Last Modified : Tue 01 Aug 2017 11:15:32 PM EDT

* Created By : Shijie Shu

"""

import pandas as pd
import numpy as np
import math

def aggre_profden_to_profstock(nlev, dz, profile):

    """ Aggregate the SOC density profile into a stock profile. 
    Input:
        nlev --- number of soil layers, scalar
        dz --- depth of each soil layer, nlev*1
        profile --- SOC profile in kgC/m3, nlev*nprof
    Output:
        val --- the SOC profile after aggregation, here is kgC/m2, nlev*nprof
    """

    nprof = len(profile)
    val = np.zeros((nprof, nlev))
    for i in range(0, nprof):
        if(np.isnan(profile[i, 0])):
            val[i, :] = float("nan")
        else:
            val[i, 0] = profile[i, 0]*dz[0]
        for j in range(1, nlev):
            if(np.isnan(profile[i, j])):
                break
            else:
                val[i, j] = val[i, j-1] + profile[i,j]*dz[j]
    return val

def aggre_profden_to_stock(endlev, dz, profile):

    """ Aggregate the SOC density profile into a total stock value
    Input:
        endlev --- number of soil layers to be aggregated
        dz --- depth of each soil layer (m), nlev*1
        profile --- SOC profile in kgC/m3, nlev*nprof
    Output:
        val --- The aggregated SOC stock, 1*nprof
    """

    nprof = len(profile)
    val = np.zeros(nprof)
    for i in range(0, nprof):
        val[i] = profile[i, 0]*dz[0]
        for j in range(1, endlev):
            if(np.isnan(profile[i, j])):
                break
            else:
                val[i] = val[i] + profile[i,j]*dz[j]
    return val

def aggre_profden_to_stock_1m(zsoih, dz, profile):

    """ Aggregate the SOC density profile into a total stock value
        till 1m specifically for ISAM output
    Input:
        zsoih --- depth of the interface of each soil layer, (nlev+1) * 1
        dz --- depth of each soil layer (m), nlev*1
        profile --- SOC profile in kgC/m3, nlev*nprof
    Output:
        val --- The aggregated SOC stock, 1*nprof
    """

    nprof = len(profile)
    val = np.zeros(nprof)
    for i in range(0, nprof):
        val[i] = profile[i, 0]*dz[0]
        for j in range(1, 7):
            if(np.isnan(profile[i, j])):
                break
            else:
                val[i] = val[i] + profile[i,j]*dz[j]
        if(~np.isnan(profile[i, 7])):
            val[i] = val[i] + profile[i, 7] * (1-zsoih[7])
    return val

