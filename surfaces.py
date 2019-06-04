#!/usr/bin/env python

import sys
import numpy as np
from mayavi import mlab

data = open(sys.argv[1],'r').readlines()

strands = {}
n=1
for i in data[1:]:
    if i.split()[0] == 'strand':
        strands[int(i.split()[-1])]= data[n+1].split()[:-1]
    n+=1


X = np.array([float(x) for x in data[34].split()])
Y = np.array([float(x) for x in data[35].split()])
Z = np.array([float(x) for x in data[36].split()])
d = np.array([float(x) for x in data[37].split()])
 

Xt = np.array([float(x) for x in data[39].split()])
Yt = np.array([float(x) for x in data[40].split()])
Zt = np.array([float(x) for x in data[41].split()])
dt = np.array([float(x) for x in data[42].split()])


# Define the points in 3D space
# including color code based on Z coordinate.

pts = mlab.points3d(X, Y, Z, scale_factor=5)
pts.glyph.scale_mode = 'scale_by_vector'
pts.mlab_source.dataset.point_data.scalars = d

ptst = mlab.points3d(Xt, Yt, Zt, scale_factor=5)
ptst.glyph.scale_mode = 'scale_by_vector'
ptst.mlab_source.dataset.point_data.scalars = dt

# Triangulate based on X, Y with Delaunay 2D algorithm.
# Save resulting triangulation.

mesh = mlab.pipeline.delaunay2d(pts)
mesh2 = mlab.pipeline.delaunay2d(ptst)

# Remove the point representation from the plot
pts.remove()
ptst.remove()

def plotstrand(stranddata):
    xs,ys,zs = [],[],[]
    for i in stranddata:
        j = i.split(',')
        xs.append(float(j[0]))
        ys.append(float(j[1]))
        zs.append(float(j[2]))
    print(xs)
    print(ys)
    print(zs)
    strand = mlab.plot3d(xs, ys, zs, tube_radius=0.25)

for i in strands:
    print strands[i]
    plotstrand(strands[i])

# Draw a surface based on the triangulation
surf = mlab.pipeline.surface(mesh)
surf2 = mlab.pipeline.surface(mesh2)
for i in strands:
    plotstrand(strands[i])
    
# Simple plot.
mlab.axes(x_axis_visibility=False,y_axis_visibility=False,z_axis_visibility=False)
mlab.show()
