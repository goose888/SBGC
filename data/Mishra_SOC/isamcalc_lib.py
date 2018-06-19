"""

* File Name : isamcalc_lib.py

* Purpose : Lib for processing ISAM output

* Creation Date : 18-07-2017

* Last Modified : Tue 01 Aug 2017 11:15:32 PM EDT

* Created By : Shijie Shu

"""

import pandas as pd
import numpy as np

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
        z_mod[k] = (-1./(1. + np.exp(9.5*((z[k]-z[0])-z[0])))) / (-1./(1.+np.exp(9.5*(-z[0])) ) )

    return z_mod
