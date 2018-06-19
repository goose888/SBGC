"""

* File Name : socplot_lib.py

* Purpose : Lib for making plots for soc data

* Creation Date : 18-07-2017

* Last Modified : Wed 19 Jul 2017 10:55:37 PM EDT

* Created By : Shijie Shu

"""

import pandas as pd
import numpy as np
import matplotlib
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

def plot_obsvsmod(obs, Yobs, mod, Ymod, tit, path, mod2=None, Ymod2=None, xticks=None, fontsize_ax=24, fontsize_lg=30, legendps='upper left'):
    """ Make plot for a single D14C profile with both observation and modeled
    outputs
    This updated version included more settings on axes and title 
    Input:
        obs --- Observed SOC profile (kgC m-3)
        Yobs --- Corresponding soil depth for observation (cm)
        mod --- Modeled SOC profile (kgC m-3)
        Ymod --- Corresponding soil depth for model (cm)
        mod2 --- Modeled SOC profile for another ISAM case(kgC m-3)
        Ymod2 --- Corresponding soil depth for model (cm)
        tit --- Title of the figure
        path --- Name and path of the output figure (PNG)
        xticks --- The ticks for x axis (if set) 
        fontsize_ax --- Font size of the axis
        fontsize_lg --- Font size of the legend
        legendps --- Position of the legend
    Output:
        status --- Return value if the figure is closed successfully
    """

    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(16,10))
    ax1 = fig.axes[0]
    plt.gca().invert_yaxis()
    # ax2 = ax1.twiny()
    # Define resources, will be added later
    #axfont = {'family' : 'normal',
    #          'weight' : 'bold',
    #          'size'   : 22}
    #matplotlib.rc('font', **font)
    # Y = np.arange(1,len(socprof))     # 1cm to 200cm
    # total SOC
    h1 = ax1.plot(obs.astype(float),Yobs.astype(float),marker='o',color='g',linewidth=2.2)
    h2 = ax1.plot(mod.astype(float),Ymod.astype(float),color='b',linewidth=2.2)
    if mod2 is not None:
        h3 = ax1.plot(mod2.astype(float),Ymod2.astype(float),color='b',linewidth=2.2)
    ax1.grid(color='gray', which='major', axis='both', alpha=0.3)
    if mod2 is not None:
        ax1.legend(['obs', 'ISAM', 'ISAM_org'], loc=legendps, fontsize=fontsize_lg)
    else:
        ax1.legend(['obs', 'modeled'], loc=legendps, fontsize=fontsize_lg)
    if xticks is not None:
        ax1.set_xticks(xticks, minor=False)
    else:
        start, end = ax1.get_xlim()
        ax1.xaxis.set_ticks(np.arange(start, end,(end-start)/6))
    # Set font size
    for item in ([ax1.title, ax1.xaxis.label, ax1.yaxis.label] +
             ax1.get_xticklabels() + ax1.get_yticklabels()):
        item.set_fontsize(20)
    # ax1.xaxis.label.set_fontsize(fontsize_ax)
    # ax1.yaxis.label.set_fontsize(fontsize_ax)
    # ax1.get_xticklabels().set_fontsize(fontsize_ax)
    # ax1.get_yticklabels().set_fontsize(fontsize_ax)
    for tl in ax1.get_xticklabels():
        tl.set_color('g')
    fig.suptitle(tit, fontsize=30)
    fig.tight_layout() # Or equivalently,  "plt.tight_layout()"
    # fig.savefig(path)
    plt.show()
    status = plt.close("all")

    return status

def plot_sensitivity(obs, Yobs, mod, Ymod, tit, path, mod2=None, Ymod2=None, xticks=None, fontsize_ax=24, fontsize_lg=30):
    """ Make plot for a single SOC profile with both observation and modeled
    output
    This updated version included more settings on axes and title 
    Input:
        obs --- Observed SOC profile (kgC m-3)
        Yobs --- Corresponding soil depth for observation (cm)
        mod --- Modeled SOC profile for high value case (kgC m-3)
        Ymod --- Corresponding soil depth for model (cm)
        mod2 --- Modeled SOC profile for low value case (kgC m-3)
        Ymod2 --- Corresponding soil depth for model (cm)
        tit --- Title of the figure
        path --- Name and path of the output figure (PNG)
        xticks --- The ticks for x axis (if set) 
        fontsize_ax --- Font size of the axis
        fontsize_lg --- Font size of the legend
    Output:
        status --- Return value if the figure is closed successfully
    """

    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(16,10))
    ax1 = fig.axes[0]
    plt.gca().invert_yaxis()
    # ax2 = ax1.twiny()
    # Define resources, will be added later
    #axfont = {'family' : 'normal',
    #          'weight' : 'bold',
    #          'size'   : 22}
    #matplotlib.rc('font', **font)
    # Y = np.arange(1,len(socprof))     # 1cm to 200cm
    # total SOC
    h1 = ax1.plot(obs.astype(float),Yobs.astype(float),marker='o',color='g',linewidth=2.2)
    h2 = ax1.plot(mod.astype(float),Ymod.astype(float),color='b',linewidth=2.2)
    if mod2 is not None:
        h3 = ax1.plot(mod2.astype(float),Ymod2.astype(float),color='r',linewidth=2.2)
    ax1.grid(color='gray', which='major', axis='both', alpha=0.3)
    if mod2 is not None:
        ax1.legend(['obs', 'high end', 'low end'], loc='upper left', fontsize=fontsize_lg)
    else:
        ax1.legend(['obs', 'modeled'], loc='upper left', fontsize=fontsize_lg)
    ax1.set_xlabel(r"SOC (kgC m-3)",color='g')
    if xticks is not None:
        ax1.set_xticks(xticks, minor=False)
    else:
        start, end = ax1.get_xlim()
        ax1.xaxis.set_ticks(np.arange(start, end,(end-start)/6))
    # Set font size
    for item in ([ax1.title, ax1.xaxis.label, ax1.yaxis.label] +
             ax1.get_xticklabels() + ax1.get_yticklabels()):
        item.set_fontsize(20)
    # ax1.xaxis.label.set_fontsize(fontsize_ax)
    # ax1.yaxis.label.set_fontsize(fontsize_ax)
    # ax1.get_xticklabels().set_fontsize(fontsize_ax)
    # ax1.get_yticklabels().set_fontsize(fontsize_ax)
    for tl in ax1.get_xticklabels():
        tl.set_color('g')
    fig.suptitle(tit, fontsize=30)
    fig.tight_layout() # Or equivalently,  "plt.tight_layout()"
    fig.savefig(path)
    status = plt.close("all")

    return status

def plot_obsvsmodplussensi(obs, Yobs, mod, Ymod, sen, Ysen, cases, palette, tit, path, xticks=None, fontsize_ax=24, fontsize_lg=30, legendps='upper left'):
    """ Make plot for a single D14C profile with both observation and modeled
    outputs
    This updated version included more settings on axes and title 
    Input:
        obs --- Observed SOC profile (kgC m-3)
        Yobs --- Corresponding soil depth for observation (cm)
        mod --- Modeled SOC profile (kgC m-3)
        Ymod --- Corresponding soil depth for model (cm)
        sen --- Modeled SOC profile with puertubed parameters
        Ysen --- Corresponding soil depth for model (cm)
        cases --- Dictionary recording the case name and the corresponding case output
        palette --- Dictionary recording the color for the corresponding case output
        tit --- Title of the figure
        path --- Name and path of the output figure (PNG)
        xticks --- The ticks for x axis (if set) 
        fontsize_ax --- Font size of the axis
        fontsize_lg --- Font size of the legend
        legendps --- Position of the legend
    Output:
        status --- Return value if the figure is closed successfully
    """

    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(16,10))
    ax1 = fig.axes[0]
    plt.gca().invert_yaxis()
    # ax2 = ax1.twiny()
    # Define resources, will be added later
    #axfont = {'family' : 'normal',
    #          'weight' : 'bold',
    #          'size'   : 22}
    #matplotlib.rc('font', **font)
    # Y = np.arange(1,len(socprof))     # 1cm to 200cm
    h1 = ax1.plot(obs.astype(float),Yobs.astype(float),marker='o',color='g',linewidth=2.2)
    h2 = ax1.plot(mod.astype(float),Ymod.astype(float),color='b',linewidth=2.2)
    # Here Ysen shall have totally seven sets of experiments
    ax1.grid(color='gray', which='major', axis='both', alpha=0.3)
    for cname, ind in palette.items():
        ax1.plot(sen[ind].astype(float), Ysen.astype(float),color=cname,linewidth=2.2)
    if xticks is not None:
        ax1.set_xticks(xticks, minor=False)
    else:
        start, end = ax1.get_xlim()
        ax1.xaxis.set_ticks(np.arange(start, end,(end-start)/6))
    # Set font size
    for item in ([ax1.title, ax1.xaxis.label, ax1.yaxis.label] +
             ax1.get_xticklabels() + ax1.get_yticklabels()):
        item.set_fontsize(20)
    # ax1.xaxis.label.set_fontsize(fontsize_ax)
    # ax1.yaxis.label.set_fontsize(fontsize_ax)
    # ax1.get_xticklabels().set_fontsize(fontsize_ax)
    # ax1.get_yticklabels().set_fontsize(fontsize_ax)
    for tl in ax1.get_xticklabels():
        tl.set_color('g')
    fig.suptitle(tit, fontsize=30)
    fig.tight_layout() # Or equivalently,  "plt.tight_layout()"
    # fig.savefig(path)
    plt.show()
    status = plt.close("all")

    return status

def plot_obsvsmodshades(obs, Yobs, mod, Ymod, sen1, sen2, Ysen, cases, palette, tit, path, xticks=None, fontsize_ax=32, fontsize_lg=48, legendps='upper left'):
    """ Make plot for a single D14C profile with both observation and modeled
    outputs
    This updated version included more settings on axes and title 
    Input:
        obs --- Observed SOC profile (kgC m-3)
        Yobs --- Corresponding soil depth for observation (cm)
        mod --- Modeled SOC profile (kgC m-3)
        Ymod --- Corresponding soil depth for model (cm)
        sen --- Modeled SOC profile with puertubed parameters
        Ysen --- Corresponding soil depth for model (cm)
        cases --- Dictionary recording the case name and the corresponding case output
        palette --- Dictionary recording the color for the corresponding case output
        tit --- Title of the figure
        path --- Name and path of the output figure (PNG)
        xticks --- The ticks for x axis (if set) 
        fontsize_ax --- Font size of the axis
        fontsize_lg --- Font size of the legend
        legendps --- Position of the legend
    Output:
        status --- Return value if the figure is closed successfully
    """

    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(16,10))
    ax1 = fig.axes[0]
    plt.ylim((-2, 142))
    plt.gca().invert_yaxis()
    # ax2 = ax1.twiny()
    # Define resources, will be added later
    #axfont = {'family' : 'normal',
    #          'weight' : 'bold',
    #          'size'   : 22}
    #matplotlib.rc('font', **font)
    # Y = np.arange(1,len(socprof))     # 1cm to 200cm
    h1 = ax1.plot(obs.astype(float),Yobs.astype(float),marker='o',color='g',linewidth=2.2)
    h2 = ax1.plot(mod.astype(float),Ymod.astype(float),color='b',linewidth=2.2)
    # Here Ysen shall have totally seven sets of experiments
    ax1.grid(color='gray', which='major', axis='both', alpha=0.3)
   # yticks = (0, 20, 40, 60, 80, 100, 120, 140)
   # ax1.set_yticklabels(yticks, fontsize=32, minor=False)
    ax1.fill_betweenx(Ysen, sen1, sen2, alpha=.25, label='Uncertainty from parameters')
    ax1.plot(sen1, Ysen, alpha=.25, marker='o',color='b', markersize=10, linewidth=1.2)
    ax1.plot(sen2, Ysen, alpha=.25, marker='o',color='b', markersize=10, linewidth=1.2)
    # Set legend!
    # ax1.legend(['obs', 'ISAM'], loc=legendps, fontsize=fontsize_lg )
    if xticks is not None:
        ax1.set_xticks(xticks, minor=False)
   #     ax1.set_xticklabels(xticks, fontsize=32, minor=False)
    else:
        start, end = ax1.get_xlim()
        ax1.xaxis.set_ticks(np.arange(start, end,(end-start)/6))
    # Set font size
    for item in ([ax1.title, ax1.xaxis.label, ax1.yaxis.label] +
             ax1.get_xticklabels() + ax1.get_yticklabels()):
        item.set_fontsize(40)
    # ax1.xaxis.label.set_fontsize(fontsize_ax)
    # ax1.yaxis.label.set_fontsize(fontsize_ax)
    # ax1.get_xticklabels().set_fontsize(fontsize_ax)
    # ax1.get_yticklabels().set_fontsize(fontsize_ax)
    for tl in ax1.get_xticklabels():
        tl.set_color('g')
    # fig.suptitle(tit, fontsize=30)
    fig.tight_layout() # Or equivalently,  "plt.tight_layout()"
    #plt.savefig(path)
    #plt.show()
    status = 0
    #status = plt.close("all")
    #status = plt.close(fig)

    return status

def plot_obsvsmodshades_color(obs, Yobs, mod, Ymod, sen1, sen2, Ysen, cases, palette, tit, path, xticks=None, fontsize_ax=32, fontsize_lg=48, legendps='upper left'):
    """ Make plot for a single D14C profile with both observation and modeled
    outputs
    This updated version included more settings on axes and title 
    Input:
        obs --- Observed SOC profile (kgC m-3)
        Yobs --- Corresponding soil depth for observation (cm)
        mod --- Modeled SOC profile (kgC m-3)
        Ymod --- Corresponding soil depth for model (cm)
        sen --- Modeled SOC profile with puertubed parameters
        Ysen --- Corresponding soil depth for model (cm)
        cases --- Dictionary recording the case name and the corresponding case output
        palette --- Dictionary recording the color for the corresponding case output
        tit --- Title of the figure
        path --- Name and path of the output figure (PNG)
        xticks --- The ticks for x axis (if set) 
        fontsize_ax --- Font size of the axis
        fontsize_lg --- Font size of the legend
        legendps --- Position of the legend
    Output:
        status --- Return value if the figure is closed successfully
    """

    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(16,10))
    ax1 = fig.axes[0]
    plt.ylim((-5, 150))
    yticks = (0, 50, 100, 150)
    ax1.set_yticks(yticks, minor=False)
    plt.gca().invert_yaxis()
    plt.xlim((-1100, 300))
    # ax2 = ax1.twiny()
    # Define resources, will be added later
    #axfont = {'family' : 'normal',
    #          'weight' : 'bold',
    #          'size'   : 22}
    #matplotlib.rc('font', **font)
    # Y = np.arange(1,len(socprof))     # 1cm to 200cm
    h1 = ax1.plot(obs.astype(float),Yobs.astype(float),marker='o',markersize=18,color='g',linewidth=3.2)
    h2 = ax1.plot(mod.astype(float),Ymod.astype(float),color='b',linewidth=3.2)
    # Here Ysen shall have totally seven sets of experiments
    ax1.grid(color='gray', which='major', axis='both', alpha=0.3)
    ax1.fill_betweenx(Ysen, sen1, sen2, alpha=.25, label='Uncertainty from parameters')
    ax1.plot(sen1, Ysen, alpha=.25, marker='o',color='r', markersize=18, linewidth=2.2)
    ax1.plot(sen2, Ysen, alpha=.25, marker='o',color='b', markersize=18, linewidth=2.2)
    # Set legend!
    # ax1.legend(['obs', 'ISAM'], loc=legendps, fontsize=fontsize_lg )
    if xticks is not None:
        ax1.set_xticks(xticks, minor=False)
   #     ax1.set_xticklabels(xticks, fontsize=32, minor=False)
    else:
        start, end = ax1.get_xlim()
        ax1.xaxis.set_ticks(np.arange(start, end,(end-start)/6))
    # Set font size
    for item in ([ax1.title, ax1.xaxis.label, ax1.yaxis.label] +
             ax1.get_xticklabels() + ax1.get_yticklabels()):
        item.set_fontsize(40)
    # ax1.xaxis.label.set_fontsize(fontsize_ax)
    # ax1.yaxis.label.set_fontsize(fontsize_ax)
    # ax1.get_xticklabels().set_fontsize(fontsize_ax)
    # ax1.get_yticklabels().set_fontsize(fontsize_ax)
    for tl in ax1.get_xticklabels():
        tl.set_color('g')
    # fig.suptitle(tit, fontsize=30)
    fig.tight_layout() # Or equivalently,  "plt.tight_layout()"
    status = 0
    #fig.savefig(path)
    #plt.show()
    #status = plt.close("all")

    return status

def plot_tau(obs, Yobs, obs14, Yobs14, mod, Ymod, tit, path, xticks1=None, xticks2=None, yticks=None, fontsize_ax=24, fontsize_lg=30):
    """ Make turnover plot for a single SOC profile with both observation and modeled
    outputs
    The tau figure contains two x axes and plot both the observed D14C and tau in one figure
    Input:
        obs --- Tau profile of the sample (yr-1)
        Yobs --- Corresponding soil depth for obs (cm)
        obs14 --- Observed D14C profile
        Yobs14 --- Corresponding soil depth for obs14 (cm)
        mod --- Modeled tau profile (yr-1)
        Ymod --- Corresponding soil depth for mod (cm)
        tit --- Title of the figure
        path --- Name and path of the output figure (PNG)
        xticks1 --- The ticks for x axis (if set) 
        xticks2 --- The ticks for x axis (if set) 
        yticks --- The ticks for y axis (if set) 
        fontsize_ax --- Font size of the axis
        fontsize_lg --- Font size of the legend
    Output:
        status --- Return value if the figure is closed successfully
    """

    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(16,10))

    ax1 = fig.axes[0]
    ax2 = ax1.twiny()
    axes = plt.gca()

    h1 = ax1.plot(np.array(obs14.astype(float)),np.array(Yobs14.astype(float)),marker='o',color='g',linewidth=2.2)
    ax1.set_xlabel(r"$\Delta^{14}C$ $(permil)$",color='g')
    start, end = ax1.get_xlim()
    for tl in ax1.get_xticklabels():
        tl.set_color('g')
    h2 = ax2.plot(np.array(obs.astype(float)),np.array(Yobs.astype(float)),marker='o',color='b',linewidth=2.2)
    h3 = ax2.plot(np.array(mod.astype(float)),np.array(Ymod.astype(float)),color='r',linewidth=2.2)
    ax2.set_xlabel(r"\tau (yr)")
    ax2.set_xlim([0,30000])
    axes.set_ylim([0,100])
    ax1.grid(color='gray', which='major', axis='both', alpha=0.3)
    plt.gca().invert_yaxis()

    # Ticks
    if xticks1 is not None:
        ax1.xaxis.set_ticks(xticks1, minor=False)
    else:
        start, end = ax1.get_xlim()
        ax1.xaxis.set_ticks(np.arange(start, end,(end-start)/6), minor=False)

    if xticks2 is not None:
        ax2.xaxis.set_ticks(xticks2, minor=False)
    else:
        start, end = ax2.get_xlim()
        ax2.xaxis.set_ticks(np.arange(start, end,(end-start)/6), minor=False)

    if yticks is not None:
        ax1.set_yticks(yticks, minor=False)
    else:
        start, end = ax1.get_ylim()
        ax1.set_yticks(np.arange(start, end,(end-start)/5))

    # Legend
    axes.legend(['obs D14C', 'obs tau', 'model tau'], loc='lower right', fontsize=fontsize_lg)

    # Other font sizes
    for item in ([ax1.title, ax1.xaxis.label, ax1.yaxis.label] +
             ax1.get_xticklabels() + ax1.get_yticklabels()):
        item.set_fontsize(20)


    # fig.tight_layout() # Or equivalently,  "plt.tight_layout()"

    fig.suptitle(tit, fontsize=30)
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
    plt.ylim((0, 110))
    plt.gca().invert_yaxis()
   # ax2 = ax1.twiny()
    font = {'family' : 'Times New Rome',
            'weight' : 'bold',
            'size'   : 22}
    plt.rc('font', **font)
    h1=ax1.errorbar(Xobs.astype(float),Yobs.astype(float),color='g',xerr=Xobs_std.astype(float))
    h1=ax1.errorbar(Xmod.astype(float),Ymod.astype(float),color='b',xerr=Xmod_std.astype(float))
    ax1.legend(['obs', 'modeled'], loc='upper right')
    ax1.set_xlabel(r"SOC (kgC m-3)",color='g')
    ax1.set_ylabel(r"Depth $(cm)$")
    start, end = ax1.get_xlim()
    ax1.xaxis.set_ticks(np.arange(start, end,(end-start)/6))
    for tl in ax1.get_xticklabels():
        tl.set_color('g')

    xticks = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70]
    yticks = [0, 20, 40, 60, 80, 100]
    ax1.set_xticks(xticks, minor=False)
    ax1.set_yticks(yticks, minor=False)

    fig.suptitle(tit, fontsize=20)
    fig.tight_layout() # Or equivalently,  "plt.tight_layout()"
    fig.savefig(path)
    status = plt.close("all")

    return status


