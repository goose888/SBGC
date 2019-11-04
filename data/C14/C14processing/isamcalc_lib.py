"""

* File Name : isamcalc_lib.py

* Purpose : Lib for calculations related to ISAM or ISAM output

* Creation Date : 18-07-2017

* Last Modified : Tue Sep 24 05:36:35 2019

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
                ub = np.int(np.floor(zsoih[j] * 100.0))
                lb = np.int(np.floor(zsoih[j+1] * 100.0))
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
                ub = np.int(np.floor(zsoih[j] * 100.0))
                lb = np.int(np.floor(zsoih[j+1] * 100.0))
                val[i, j] = np.nanmean(profile[i, ub:lb])
                # if(np.isnan(val[i, j])):
                #     val[i, j] = val[i, j-1]
                # else:
                #     val[i, j] = val[i, j-1]
            # Only consider the mean of the SOC density above 1m for level 8
            if(j == 7):
                ub = np.int(np.floor(zsoih[j] * 100.0))
                lb = np.int(np.floor(1 * 100.0))
                val[i, j] = np.nanmean(profile[i, ub:lb])
    return val

def latlon_2_idx(lat, lon):

    """ Get the mean value of a profile by averaging the values within each layer defined in ISAM model
    Input:
        lat --- latitude (-90 ~ 90)
        lon --- longitude (-180 ~ 180)
    Output:
        loc --- lat and lon index used by ISAM
    """

    latid = int(round(lat * 2 + 0.5) + 180)
    if (lon>=0):
        lonid = int(round(lon * 2 + 0.5))
    else:
        lonid = int(360 + round((180 + lon) * 2 + 0.5))
    loc = [latid, lonid]

    return loc

def agg_1m_soil(prof):
    """ Aggregate the value of the 1m soil (0-100cm) by adding the first 7 layers and 30% of the 8th 
        layer of ISAM soil output
    Input:
        prof --- soil profile (needs at least 8 layers)
    Output:
        val --- scalar
    """

    val = np.sum(prof[0:7]) + 0.3*prof[7]

    return val

def avg_wt_1m_soil(weight, prof):

    """ Get the mean value of the 1m soil (0-100cm) by averaging the first 7 layers and 30% of the 8th
        layer of ISAM soil output using the predefined weight
    Input:
        weight --- the values used as weight. No need to standardize before calling this subroutine
                   but the length must be consistent to the prof
        prof --- soil profile (needs at least 5 layers)
    Output:
        val --- scalar
    """

    weight_tot = np.sum(weight[0:7]) + 0.3*weight[7]
    wt = np.zeros((8))
    wt[0:7] = weight[0:7] / weight_tot
    wt[7] = 0.3*weight[7] / weight_tot
    val = np.nansum(wt*prof[0:8])

    return val

def agg_topsoil(prof):
    """ Aggregate the value of the top soil (0-30cm) by adding the first 5 layer of ISAM soil output
    Input:
        prof --- soil profile (needs at least 5 layers)
    Output:
        val --- scalar
    """

    val = np.sum(prof[0:5])

    return val

def avg_wt_topsoil(weight, prof):

    """ Get the mean value of the top soil (0-30cm) by averaging the first 5 layer of ISAM soil output
        using the predefined weight
    Input:
        weight --- the values used as weight. No need to standardize before calling this subroutine
                   but the length must be consistent to the prof
        prof --- soil profile (needs at least 5 layers)
    Output:
        val --- scalar
    """

    weight_tot = np.sum(weight[0:5])
    wt = weight[0:5] / weight_tot
    val = np.nansum(wt*prof[0:5])

    return val

def agg_subsoil(prof):
    """ Aggregate the value of the subsoil (30-100cm) by adding the 6-8 layer of ISAM soil output
    Input:
        prof --- soil profile (needs at least 8 layers), no NA is acceptable
    Output:
        val --- scalar
    """

    val = prof[5]+prof[6]+0.3*prof[7]

    return val


def avg_wt_subsoil(weight, prof):

    """ Get the mean value of the subsoil (30-100cm) by averaging the first 5 layer of ISAM soil output
        using the predefined weight
    Input:
        weight --- the values used as weight. No need to standardize before calling this subroutine
                   but the length must be consistent to the prof
        prof --- soil profile (needs at least 8 layers)
    Output:
        val --- scalar
    """

    weight_tot = weight[5]+weight[6]+0.3*weight[7]
    wt = [0, 0, 0]
    for i in np.arange(0,3):
        if (i<=1):
            wt[i] =  weight[(i+5)] / weight_tot
        else:
            wt[i] = 0.3*weight[(i+5)] / weight_tot
    val = np.nansum(wt*prof[5:8])

    return val

def q10_1dag(coef, temp):

    """ Get the Aboveground Q10 temperature factor based on ISAM-1D calculation.
    Input:
        coef --- the Q10 coefficient
        temp --- mean soil temperature
    Output:
        val --- scalar, q10 factor
    """

    val = 1.0 * coef**((temp-25.0)/10.0)

    return val

def q10_1dbg(coef, temp):

    """ Get the Belowground Q10 temperature factor based on ISAM-1D calculation.
    Input:
        coef --- the Q10 coefficient
        temp --- mean soil temperature
    Output:
        val --- scalar, q10 factor
    """

    val = 1.3 * coef**((temp-10.0)/10.0) - 1.3*2.0**(-2.0)

    return val

def tsen_0dag(temp):

    """ Get the Aboveground Q10 temperature factor based on ISAM-1D calculation.
    Input:
        coef --- the Q10 coefficient
        temp --- mean soil temperature
    Output:
        val --- scalar, q10 factor
    """

    tmod = 0.73
    val = tmod*0.09*np.exp(0.095*temp)

    return val

def tsen_0dbg(temp):

    """ Get the Belowground Q10 temperature factor based on ISAM-1D calculation.
    Input:
        coef --- the Q10 coefficient
        temp --- mean soil temperature
    Output:
        val --- scalar, q10 factor
    """

    temp_in = temp
    if(temp < -10.0):
        temp_in = -10.0

    val = (47.91 / (1.0 + np.exp(106.06/(temp_in+18.27))))

    return val


