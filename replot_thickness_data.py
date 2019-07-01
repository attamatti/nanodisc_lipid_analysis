#!/usr/bin/env python

import sys
import numpy as np
import glob
import matplotlib.pyplot as plt

try:
    gridx = np.fromfile('plotfiles/gridx.npy')
    print('FOUND: plotfiles/gridx.npy')
    gridy = np.fromfile('plotfiles/gridy.npy')
    print('FOUND: plotfiles/gridy.npy')
    drawMSPx = np.fromfile('plotfiles/msps_x.npy')
    print('FOUND: plotfiles/mspsx.npy')
    drawMSPy = np.fromfile('plotfiles/msps_y.npy')
    print('FOUND: plotfiles/mspsy.npy')
    num_extra_elements = len(glob.glob('plotfiles/elems_x_*.npy'))

    draw_elementsx,draw_elementsy = [],[]
    xee = glob.glob('plotfiles/elems_x_*.npy')
    xee.sort()
    for i in xee:
        vals= np.fromfile(i)
        draw_elementsx.append(vals)
        print('FOUND: {0}'.format(i))
    yee = glob.glob('plotfiles/elems_y_*.npy')
    yee.sort()
    for i in yee:
        vals= np.fromfile(i)
        draw_elementsy.append(vals)
        print('FOUND: {0}'.format(i))

    thickness_mean = np.fromfile('plotfiles/thickness.npy')
    print('FOUND: plotfiles/thickness.npy')

    thickness_std = np.fromfile('plotfiles/STD.npy')
    print('FOUND: pltofiles/STD.npy')

except:
    sys.exit('ERROR READING FILES!! - they should be in a directoty called plotfiles/ \nUSAGE: replot_thickness_data.py <thickness min> <thickness max> <thickness increment> <std min> <std max> <std increment> ')
try:
    thmin = int(sys.argv[1])
    thmax = int(sys.argv[2])
    thinc = int(sys.argv[3])
    stdmin = int(sys.argv[4])
    stdmax = int(sys.argv[5])
    stdinc = int(sys.argv[6])

except:
    sys.exit('''\nUSAGE: replot_thickness_data.py <thickness min> <thickness max> <thickness increment> <std min> <std max> <std increment>
             values analysis:
             Thickness min = {0}
             Thickness max = {1}
             STD min = {2}
             STD max = {3}'''.format(np.min(thickness_mean),np.max(thickness_mean),np.min(thickness_std),np.max(thickness_std)))

thickness_mean = thickness_mean.reshape(150,150)
thickness_std = thickness_std.reshape(150,150)

#### plot mean thickness
h = plt.contourf(gridx,gridy,thickness_mean,np.arange(thmin,thmax,thinc),vmin=thmin,vmax=thmax,cmap='coolwarm')
plt.xlim = (np.mean(drawMSPx)-75,np.mean(drawMSPx)+75)
plt.ylim = (np.mean(drawMSPy)-75,np.mean(drawMSPy)+75)
plt.gca().set_aspect('equal', adjustable='datalim')
plt.axis('off')
plt.scatter(drawMSPx,drawMSPy,c='K')
markers = ('v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd')
colors = ['k','w']
mcount,ccount = (0,0)
for i in range(num_extra_elements):
    plt.scatter(draw_elementsx[i],draw_elementsy[i],marker=markers[mcount],c=colors[ccount],edgecolors='k')
    mcount+=1
    ccount +=1
    if mcount > len(markers)-1:
        mcount=0
    if ccount >1:
        ccount=0
plt.colorbar(h)
plt.savefig('thickness.png')
plt.show()
plt.close()


#### plot the std
h = plt.contourf(gridx,gridy,thickness_std,np.arange(stdmin,stdmax,stdinc),vmin=stdmin,vmax=stdmax,cmap='YlOrRd')
plt.scatter(drawMSPx,drawMSPy,c='K')
markers = ('v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd')
colors = ['k','w']
mcount,ccount = (0,0)
plt.xlim = (np.mean(drawMSPx)-75,np.mean(drawMSPx)+75)
plt.ylim = (np.mean(drawMSPy)-75,np.mean(drawMSPy)+75)
plt.gca().set_aspect('equal', adjustable='datalim')
plt.axis('off')
for i in range(num_extra_elements):
    plt.scatter(draw_elementsx[i],draw_elementsy[i],marker=markers[mcount],c=colors[ccount],edgecolors='k')
    mcount+=1
    ccount +=1
    if mcount > len(markers)-1:
        mcount=0
    if ccount >1:
        ccount=0
plt.colorbar(h)
plt.savefig('STD.png')
plt.show()
plt.close()