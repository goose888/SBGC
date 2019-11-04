# -*- Python source file -*-
"""
File Name : mycolorbars.py

Description : Libs for creating and managing self-defined colorbars.

Created on : Mon Feb 25 04:17:48 2019

Last Modified : Tue Aug 27 20:51:54 2019

Author : Shijie Shu
"""

'''
NAME
    Custom Colormaps for Matplotlib
PURPOSE
    This program shows how to implement make_cmap which is a function that
    generates a colorbar.  If you want to look at different color schemes,
    check out https://kuler.adobe.com/create.
PROGRAMMER(S)
    Chris Slocum
REVISION HISTORY
    20130411 -- Initial version created
    20140313 -- Small changes made and code posted online
    20140320 -- Added the ability to set the position of each color
'''

def make_cmap(colors, position=None, bit=False):
    '''
    make_cmap takes a list of tuples which contain RGB values. The RGB
    values may either be in 8-bit [0 to 255] (in which bit must be set to
    True when called) or arithmetic [0 to 1] (default). make_cmap returns
    a cmap with equally spaced colors.
    Arrange your tuples so that the first color is the lowest value for the
    colorbar and the last is the highest.
    position contains values from 0 to 1 to dictate the location of each color.
    '''
    import matplotlib as mpl
    import numpy as np
    bit_rgb = np.linspace(0,1,256)
    if position == None:
        position = np.linspace(0,1,len(colors))
    else:
        if len(position) != len(colors):
            sys.exit("position length must be the same as colors")
        elif position[0] != 0 or position[-1] != 1:
            sys.exit("position must start with 0 and end with 1")
    if bit:
        for i in range(len(colors)):
            colors[i] = (bit_rgb[colors[i][0]],
                         bit_rgb[colors[i][1]],
                         bit_rgb[colors[i][2]])
    cdict = {'red':[], 'green':[], 'blue':[]}
    for pos, color in zip(position, colors):
        cdict['red'].append((pos, color[0], color[0]))
        cdict['green'].append((pos, color[1], color[1]))
        cdict['blue'].append((pos, color[2], color[2]))

    cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return cmap

def color_ch4diff():
    '''
    Simply generate the colorbar specifically for the CH4 diffusion spatial plot.
    '''
    from matplotlib import cm
    import numpy as np
    mycolors = []
    for i in np.arange(0,128):
        idx = i*2
        #if(i == 32 or i == 33 or i == 34):
        #if(i == 43 or i == 44 or i == 45):
        if(i == 63 or i == 64 or i == 65):
            mycolors.append((1.0,1.0,1.0))
            #mycolors.append(cm.jet(idx)[0:3])
        else:
            mycolors.append(cm.jet(idx)[0:3])
    mycmap = make_cmap(mycolors)
    # mycolors2 = []
    # for i in np.arange(0,8):
    #     idx = i*32
    #     if(i == 2):
    #         mycolors2.append((1.0,1.0,1.0))
    #         #mycolors.append(cm.jet(idx)[0:3])
    #     else:
    #         mycolors.append(cm.jet(idx)[0:3])
    #mycolors[60:66] = mycolors2
    # # Retrieve all the colors from the current modified cmap and 
    # # adjust the white color to match the actual zero values from 
    # # the map. Currently the white color falls on the positive side
    # # Set step to 2
    # mycolors = []
    # for i in np.arange(0,128):
    #     idx = i*2
    #     mycolors.append(mycmap(idx)[0:3])
    # mycolors[29] = (1.0,1.0,1.0,1.0)
    # mycolors[30] = (1.0,1.0,1.0,1.0)
    # mycolors[31] = (1.0,1.0,1.0,1.0)
    # mycolors[32] = (1.0,1.0,1.0,1.0)
    # mycolors[33] = (1.0,1.0,1.0,1.0)
    # mycolors[34] = (1.0,1.0,1.0,1.0)
    # # Update the cmap
    # mycmap = make_cmap(mycolors)

    return mycmap

def color_ch4emission():
    '''
    Simply generate the colorbar specifically for the net CH4 flux spatial plot.
    '''
    from matplotlib import cm
    import numpy as np
    mycolors = []
    for i in np.arange(0,128):
        idx = i*2
        if(i == 32 or i == 33 or i == 34):
        #if(i == 43 or i == 44 or i == 45):
        #if(i == 63 or i == 64 or i == 65):
            mycolors.append((1.0,1.0,1.0))
            #mycolors.append(cm.jet(idx)[0:3])
        else:
            mycolors.append(cm.jet(idx)[0:3])
    mycmap = make_cmap(mycolors)
    # mycolors2 = []
    # for i in np.arange(0,8):
    #     idx = i*32
    #     if(i == 2):
    #         mycolors2.append((1.0,1.0,1.0))
    #         #mycolors.append(cm.jet(idx)[0:3])
    #     else:
    #         mycolors.append(cm.jet(idx)[0:3])
    #mycolors[60:66] = mycolors2
    # # Retrieve all the colors from the current modified cmap and 
    # # adjust the white color to match the actual zero values from 
    # # the map. Currently the white color falls on the positive side
    # # Set step to 2
    # mycolors = []
    # for i in np.arange(0,128):
    #     idx = i*2
    #     mycolors.append(mycmap(idx)[0:3])
    # mycolors[29] = (1.0,1.0,1.0,1.0)
    # mycolors[30] = (1.0,1.0,1.0,1.0)
    # mycolors[31] = (1.0,1.0,1.0,1.0)
    # mycolors[32] = (1.0,1.0,1.0,1.0)
    # mycolors[33] = (1.0,1.0,1.0,1.0)
    # mycolors[34] = (1.0,1.0,1.0,1.0)
    # # Update the cmap
    # mycmap = make_cmap(mycolors)

    return mycmap

