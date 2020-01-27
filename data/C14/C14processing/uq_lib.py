"""

* File Name : uq_lib.py

* Purpose : Lib for examine uncertainties of Data or ISAM model output

* Creation Date : 18-07-2017

* Last Modified : Fri Nov 29 18:07:25 2019

* Created By : Shijie Shu

"""

# This package requires "isamcalc_lib" package, please make sure you obtain it from the author

import pandas as pd
import numpy as np
import math
import isamcalc_lib as isam

def get_value_SOC(idx, soc):

    """ Get the topsoil (0-30cm) and subsoil (30-100cm) SOC
    Input:
        idx_u --- index of the line corresponds to the SOC profile simulated using high end parameter value
        idx_l --- index of the line corresponds to the SOC profile simulated using low end parameter value
        soc --- SOC profiles under sensitivity test
    Output:
        rg_topsoil --- SOC sensitivity range, computed as the half of the full range
        rg_subsoil --- SOC sensitivity range, computed as the half of the full range
    """
    val_top = isam.agg_topsoil(soc[idx,:])    # SOC profile kgCm-3

    val_sub = isam.agg_subsoil(soc[idx,:])

    return val_top, val_sub

def get_sensitivity_range_SOC(idx_u, idx_l, soc):

    """ Get the uncertainty range of the topsoil (0-30cm) and subsoil (30-100cm) 
        This is simply acheived by subtarcting the lower bounds by using upper bounds of all sensitivity simulations.
    Input:
        idx_u --- index of the line corresponds to the SOC profile simulated using high end parameter value
        idx_l --- index of the line corresponds to the SOC profile simulated using low end parameter value
        soc --- SOC profiles under sensitivity test
    Output:
        rg_topsoil --- SOC sensitivity range, computed as the half of the full range
        rg_subsoil --- SOC sensitivity range, computed as the half of the full range
    """
    tmp_u = isam.agg_topsoil(soc[idx_u,:])    # SOC profile kgCm-3
    tmp_l = isam.agg_topsoil(soc[idx_l,:])
    rg_topsoil = (tmp_u-tmp_l)/2.

    tmp_u = isam.agg_subsoil(soc[idx_u,:])
    tmp_l = isam.agg_subsoil(soc[idx_l,:])
    rg_subsoil = (tmp_u-tmp_l)/2.

    return rg_topsoil, rg_subsoil

def get_sensitivity_range_D14C(idx_u, idx_l, soc, d14c):

    """ Get the D14C uncertainty range of the topsoil (0-30cm) and subsoil (30-100cm) 
        This is simply acheived by subtarcting the lower bounds by using upper bounds of all sensitivity simulations.
        Is it correct to use SOC as weight to directly scale D14C?
        D14C is linearly transformed from 14C mass, so I believe this method is applicable.
    Input:
        idx_u --- index of the line corresponds to the SOC profile simulated using high end parameter value
        idx_l --- index of the line corresponds to the SOC profile simulated using low end parameter value
        soc --- SOC profiles under sensitivity test
        d14c --- D14C profiles under sensitivity test
    Output:
        rg_topsoil --- SOC sensitivity range, computed as the half of the full range
        rg_subsoil --- SOC sensitivity range, computed as the half of the full range
    """

    tmp_u = isam.avg_wt_topsoil(soc[idx_u,:], d14c[idx_u,:])
    tmp_l = isam.avg_wt_topsoil(soc[idx_l,:], d14c[idx_l,:])
    rg_topsoil = (tmp_u-tmp_l)/2.

    tmp_u = isam.avg_wt_subsoil(soc[idx_u,:], d14c[idx_u,:])
    tmp_l = isam.avg_wt_subsoil(soc[idx_l,:], d14c[idx_l,:])
    rg_subsoil = (tmp_u-tmp_l)/2.

    return rg_topsoil, rg_subsoil
