"""

* File Name : SOCtools.py

* Purpose : Lib for processing SOC data

* Creation Date : 18-07-2017

* Last Modified : Thu May 30 05:43:23 2019

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

def aggre_prof_den_to_stock(depth, zsoih, dz, profile):
    """ Aggregate the SOC density profile into a total stock value till the preset 
        depth. Only accepts one profile as input at one time
    Input:
        depth --- the depth where SOC aggregate to (m)
        zsoih --- depth of the interface of each soil layer, nlev+1
        dz --- thickness of each soil layer (m), nlev
        profile --- SOC profile in kgC/m3, nlev
    Output:
        val --- The aggregated SOC stock (kgC/m2), scalar. Will be nan if dz is not 
                deep enough.
    """

    nlev = len(profile)
    val = 0.
    pt = 0
    # Point to the layer right below 30cm
    for j in range(0, nlev):
        if(zsoih[j] >= depth):
            pt = j
            break
    if(zsoih[nlev] < depth):
        val = float('nan')
    else:
        # Get the carbon amount
        for j in range(0, pt):
            if(zsoih[j+1] < depth):
                val = val + dz[j] * profile[j]
            else:
                val = val + (depth - zsoih[j]) * profile[j]

    return val


