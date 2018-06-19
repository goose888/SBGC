"""

* File Name : Fluxtools.py

* Purpose : Lib for processing site levle flux data

* Creation Date : 23-03-2018

* Last Modified : Tue 01 Aug 2017 11:15:32 PM EDT

* Created By : Shijie Shu

"""

import pandas as pd
import numpy as np
import math

def mdv(srs, window=5, nstep=48):

    """ Mean diurnal variation method for filling data gaps
    in ISAM model
    Input:
        srs --- time series, missing values are represented by -9999 or -6999
        window --- Size of the moving window for MDV method, unit in day
        nstep --- Timesteps per day
    Output:
        srs_filled --- gap-filled time series
    """
    srs_filled = srs
    for i in np.arange(0, len(srs)):
        if (srs[i] < -5000):
            # Get all the available data within the window
            if (i < window*48):
                subarr = srs[i:(i+window*nstep)]
            else:
                if (i > (len(srs)-window*nstep)):
                    subarr = srs[(i-window*nstep):i]
                else:
                    subarr = srs[(i-window*nstep):(i+window*nstep)]
            ex_subarr = subarr[np.arange(0,len(subarr),nstep)]
            ex_subarr[ex_subarr < -5000] = np.float("nan")
            srs_filled[i] = np.nanmean(ex_subarr)
    srs_filled[np.isnan(srs_filled)] = -9999.

    return srs_filled

def hr2daily(srs, method="mean", nstep=48):

    """ Transfer half hourly observation to daily mean or daily total
    Input:
        srs --- time series, missing values are represented by nan, must be a multiple of nstep
        method --- way to get the daily value. mean - daily mean. agg - daily total estimated by daily mean
        nstep --- Timesteps per day, 48 stands for half hourly and 24 stands for hourly 
    Output:
        daily_srs --- daily time series
    """
    totday = np.floor(len(srs)/nstep).astype("int")
    daily_srs = float("nan") * np.ones(totday)
    for i in np.arange(0,totday):
        daily_srs[i] = np.nanmean(srs[i*nstep:(i+1)*nstep])
        if(method == "agg"):
            daily_srs[i] = daily_srs[i] * 24.

    return daily_srs

