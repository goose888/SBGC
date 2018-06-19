"""

* File Name : isamcalc_lib.py

* Purpose : Lib for processing ISAM output

* Creation Date : 18-07-2017

* Last Modified : Tue 01 Aug 2017 11:15:32 PM EDT

* Created By : Shijie Shu

"""

import pandas as pd
import numpy as np
import math

def get_isam_soildp(nlev):

    """ Get the node depth, layer thickness and layer interface used
    in ISAM model
    Input:
        nlev --- number of soil layers
    Output:
        z --- node depth (m), nlev*1
        dz --- layer thickness (m), nlev*1
        zsoih --- layer interface (m), (nlev+1)*1
    """

    z = np.zeros(nlev)
    dz = np.zeros(nlev)
    zsoih = np.zeros(nlev+1)
    # Node depth
    for j in range(0, nlev):  
        z[j] = 0.025*(np.exp(0.5*(j+0.5))-1.)  # node depths
    # Layer thickness 
    dz[0] = 0.5*(z[0]+z[1])   #=zsoih(1)
    for j in range(1,(nlev-1)):
        dz[j]= 0.5*(z[j+1]-z[j-1])    # thickness b/n two interfaces
    dz[nlev-1]= z[nlev-1]-z[nlev-2]
    # Layer interface
    zsoih[0] = 0. 
    for j in range(0,(nlev-1)):
        zsoih[j+1]= 0.5*(z[j]+z[j+1])   # interface depths
    zsoih[nlev] = z[nlev-1] + 0.5*dz[nlev-1]

    return z, dz, zsoih

def get_depth_modifier(nlev, z):

    """ Calculate the logistic depth modifier at different ISAM layer
    Input:
        nlev --- number of soil layers
        z --- node depth (m), nlev*1
    Output:
        z_mod --- The dimensionless depth modifier (m), nlev*1
    """

    z_mod = np.zeros(nlev)
    for k in range(0,nlev):
        z_mod[k] = (-1./(1. + np.exp(6.0*((z[k]-z[0])-z[0])))) / (-1./(1.+np.exp(6.0*(-z[0])) ) )

    return z_mod

def mean_by_depth(nlev, zsoih, nprof, profile):

    """ Get the mean value of a profile by averaging the values within each layer defined in ISAM model
    Input:
        nlev --- number of soil layers
        zsoih --- layer interface (m), (nlev+1)*1
        nprof --- Number of profiles , (nlev+1)*1
        profile --- Profile at vertical resolution of 1cm and for 2m (unit depends), 200*nprof
    Output:
        val --- Averaged values (unit depends), nlev*nprof
    """

    val = np.zeros((nprof, nlev))
    for i in range(0, nprof-1):
        val[i, 0] = profile[i, 0]
        for j in range(1, nlev-1):  
            if(zsoih[j] >= 2.0):
                val[i, j] = val[i, j-1]
            else:
                ub = np.floor(zsoih[j] * 100.0)
                lb = np.floor(zsoih[j+1] * 100.0)
                val[i, j] = np.nanmean(profile[i, ub:lb])
                # if(np.isnan(val[i, j])):
                #     val[i, j] = val[i, j-1]
                # else:
                #     val[i, j] = val[i, j-1]

    return val

def mean_by_depth_sep8(nlev, zsoih, nprof, profile):

    """ Get the mean value of a profile by averaging the values within each layer defined in ISAM model
    Input:
        nlev --- number of soil layers
        zsoih --- layer interface (m), (nlev+1)*1
        nprof --- Number of profiles , (nlev+1)*1
        profile --- Profile at vertical resolution of 1cm and for 2m (unit depends), 200*nprof
    Output:
        val --- Averaged values (unit depends), nlev*nprof
    """

    val = np.zeros((nprof, nlev))
    for i in range(0, nprof-1):
        val[i, 0] = profile[i, 0]
        for j in range(1, nlev-1):  
            if(zsoih[j] >= 2.0):
                val[i, j] = val[i, j-1]
            else:
                ub = np.floor(zsoih[j] * 100.0)
                lb = np.floor(zsoih[j+1] * 100.0)
                val[i, j] = np.nanmean(profile[i, ub:lb])
                # if(np.isnan(val[i, j])):
                #     val[i, j] = val[i, j-1]
                # else:
                #     val[i, j] = val[i, j-1]
            # Only consider the mean of the SOC density above 1m for level 8
            if(j == 7):
                ub = np.floor(zsoih[j] * 100.0)
                lb = np.floor(1 * 100.0)
                val[i, j] = np.nanmean(profile[i, ub:lb])
    return val

def latlon_2_idx(lat, lon):

    """ Get the mean value of a profile by averaging the values within each layer defined in ISAM model
    Input:
        lat --- latitude
        lon --- longitude
    Output:
        loc --- lat and lon index used by ISAM
    """

    latid = int(round(lat * 2 + 0.5) + 180)
    if (lon>=0):
        lonid = int(round(lon * 2 + 0.5))
    else:
        lonid = int(360 + round((360 + lon) * 2 + 0.5))
    loc = [latid, lonid]

    return loc
