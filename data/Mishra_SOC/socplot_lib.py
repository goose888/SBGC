"""

* File Name : socplot_lib.py

* Purpose : Lib for making plots for soc data

* Creation Date : 18-07-2017

* Last Modified : Wed 19 Jul 2017 10:55:37 PM EDT

* Created By : Shijie Shu

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm

def plot_soilprof(socprof, Y, tit, path):
    """ Make plot for a single SOC profile
    Input:
        socprof --- Soil profile (kgC m-3)
        Y --- Corresponding soil depth (cm)
        tit --- Title of the figure
        path --- Name and path of the output figure (PNG)
    Output:
        status --- Return value if the figure is closed successfully
    """

    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(16,10))
    ax1 = fig.axes[0]
    plt.gca().invert_yaxis()
    ax2 = ax1.twiny()
    # Y = np.arange(1,len(socprof))     # 1cm to 200cm
    # total SOC
    h1 = ax1.plot(socprof.astype(float),Y.astype(float),color='g')
    ax1.set_xlabel(r"SOC (kgC m-3)",color='g')
    start, end = ax1.get_xlim()
    ax1.xaxis.set_ticks(np.arange(start, end,(end-start)/6))
    for tl in ax1.get_xticklabels():
        tl.set_color('g')
    fig.suptitle(tit, fontsize=20)
    fig.tight_layout() # Or equivalently,  "plt.tight_layout()"
    fig.savefig(path)
    status = plt.close("all")

    return status

def plot_obsvsmod(obs, Yobs, mod, Ymod, tit, path):
    """ Make plot for a single SOC profile with both observation and modeled
    output
    Input:
        obs --- Observed SOC profile (kgC m-3)
        Yobs --- Corresponding soil depth for observation (cm)
        mod --- Modeled SOC profile (kgC m-3)
        Ymod --- Corresponding soil depth for model (cm)
        tit --- Title of the figure
        path --- Name and path of the output figure (PNG)
    Output:
        status --- Return value if the figure is closed successfully
    """

    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(16,10))
    ax1 = fig.axes[0]
    plt.gca().invert_yaxis()
    ax2 = ax1.twiny()
    # Y = np.arange(1,len(socprof))     # 1cm to 200cm
    # total SOC
    h1 = ax1.plot(obs.astype(float),Yobs.astype(float),color='g')
    h2 = ax1.plot(mod.astype(float),Ymod.astype(float),color='b')
    ax1.legend(['obs', 'modeled'], loc='upper left')
    ax1.set_xlabel(r"SOC (kgC m-3)",color='g')
    start, end = ax1.get_xlim()
    ax1.xaxis.set_ticks(np.arange(start, end,(end-start)/6))
    for tl in ax1.get_xticklabels():
        tl.set_color('g')
    fig.suptitle(tit, fontsize=20)
    fig.tight_layout() # Or equivalently,  "plt.tight_layout()"
    fig.savefig(path)
    status = plt.close("all")

    return status

def plot_profmap(west, south, east, north, lon_s, lat_s, tit, path):
    """ Make map to show the location of soil samples
    Input:
        west --- Western bound of the map (deg)
        south --- Southern bound of the map (deg)
        east --- Eastern bound of the map (deg)
        north --- Northern bound of the map (deg)
        lon_s --- Vector contains the longitude of all samples
        lat_s --- Vector contains the latitude of all samples, shall has
                  the same length as lon_s
        tit --- Title of the figure
        path --- Name and path of the output figure (PNG)
    Output:
        status --- Return value if the figure is closed successfully
    """

    fig = plt.figure()
    ax = fig.add_axes([0.05,0.05,0.9,0.9])
    m = Basemap(llcrnrlon=west, llcrnrlat=south, urcrnrlon=east, urcrnrlat=north, projection='mill', lon_0=0, lat_0=0)
    #lon, lat = np.meshgrid(lons, lats)
    x, y = m(lon_s,lat_s)
    m.drawcoastlines(linewidth=0.25)
    m.drawcountries(linewidth=0.25)
    m.drawmapboundary(fill_color='#99ffff')
    m.fillcontinents(color='grey',lake_color='#99ffff',zorder=0)
    m.scatter(lon_s,lat_s,15,marker='^',color='r',alpha=0.7,latlon=True)
    # draw parallels.
    parallels = np.arange(-90.,90.,30.)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    # draw meridians
    meridians = np.arange(0.,360.,45.)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
    ax.set_title(tit)
    fig.savefig(path)
    status = plt.close("all")

    return status

def plot_prof_with_errbar(X1, Y, Xstd, tit, path):
    """ plot SOC profile with error bar
    Input:
        X1 --- mean SOC profile (kgCm-3)
        Y --- Soil depth (cm)
        Xstd --- std error (kgCm-3)
        tit --- Title of the figure
        path --- Name and path of the output figure (PNG)
    Output:
        status --- Return value if the figure is closed successfully
    """

    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(16,10))
    ax1 = fig.axes[0]
    plt.gca().invert_yaxis()
    ax2 = ax1.twiny()
    font = {'family' : 'Times New Rome',
            'weight' : 'bold',
            'size'   : 22}
    plt.rc('font', **font)
    h1=ax1.errorbar(X1.astype(float),Y.astype(float),color='g',xerr=Xstd.astype(float))
    ax1.set_xlabel(r"SOC (kgC m-3)",color='g')
    ax1.set_ylabel(r"Depth $(cm)$")
    start, end = ax1.get_xlim()
    ax1.xaxis.set_ticks(np.arange(start, end,(end-start)/6))
    for tl in ax1.get_xticklabels():
        tl.set_color('g')
    fig.suptitle(tit, fontsize=20)
    fig.tight_layout() # Or equivalently,  "plt.tight_layout()"
    fig.savefig(path)
    status = plt.close("all")

    return status

def plot_obsvsmod_with_errbar(Xobs, Yobs, Xmod, Ymod, Xobs_std, Xmod_std, tit, path):
    """ plot SOC profile with error bar
    Input:
        Xobs --- observed mean SOC profile (kgCm-3)
        Yobs --- Soil depth (cm)
        Xmod --- modeled mean SOC profile (kgCm-3)
        Ymod --- Soil depth (cm)
        Xobs_std --- observed std error (kgCm-3)
        Xmod_std --- modeled std error (kgCm-3)
        tit --- Title of the figure
        path --- Name and path of the output figure (PNG)
    Output:
        status --- Return value if the figure is closed successfully
    """

    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(16,10))
    ax1 = fig.axes[0]
    plt.gca().invert_yaxis()
    ax2 = ax1.twiny()
    font = {'family' : 'Times New Rome',
            'weight' : 'bold',
            'size'   : 22}
    plt.rc('font', **font)
    h1=ax1.errorbar(Xobs.astype(float),Yobs.astype(float),color='g',xerr=Xobs_std.astype(float))
    h1=ax1.errorbar(Xmod.astype(float),Ymod.astype(float),color='b',xerr=Xmod_std.astype(float))
    ax1.legend(['obs', 'modeled'], loc='upper left')
    ax1.set_xlabel(r"SOC (kgC m-3)",color='g')
    ax1.set_ylabel(r"Depth $(cm)$")
    start, end = ax1.get_xlim()
    ax1.xaxis.set_ticks(np.arange(start, end,(end-start)/6))
    for tl in ax1.get_xticklabels():
        tl.set_color('g')
    fig.suptitle(tit, fontsize=20)
    fig.tight_layout() # Or equivalently,  "plt.tight_layout()"
    fig.savefig(path)
    status = plt.close("all")

    return status


