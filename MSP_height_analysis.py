#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt

def plane_from_points(p1,p2,p3):
    '''returns ax+by+cy+d=0'''
    v1 = p3 - p1
    v2 = p2 - p1
    cp = np.cross(v1, v2)
    a, b, c = cp
    d = np.dot(cp, p3)
    
    print('{0}x + {1}y + {2}z = {3}'.format(a, b, c, d))
    return([a,b,c,-d])

def proj_point_to_plane(plane,point):
    '''
    project a point on to a plane
    a(at+x)+b(bt+y)+c(ct+z)+d = 0
    '''
    a= plane[0]
    b= plane[1]
    c= plane[2]
    d= plane[3]
    x= point[0]
    y= point[1]
    z= point[2]
    aat = a**2
    ax = a*x
    bbt = b**2
    by = b*y
    cct = c**2
    cz = c*z
    ttot = aat+bbt+cct
    xyz = -1*(ax+by+cz+d)
    t= xyz/ttot
    xo = (a*t)+x
    yo = (b*t)+y
    zo = (c*t)+z
    return([xo,yo,zo])

def get_dist(p1,p2):
    dist = (p1[0]-p2[0])+(p1[1]-p2[1])+(p1[2]-p2[2])
    return(dist)

## read bild files
bildfile = open(sys.argv[1],'r').readlines()

## get CA and plane coordinates
CAs = []
for i in bildfile:
    line = i.split()
    if line[0] == '.sphere':
        CAs.append([float(x) for x in line[1:-1]])
    elif line[0] == '.polygon':
        pp1 = np.array([float(x) for x in line[1:4]])
        pp2 = np.array([float(x) for x in line[4:7]])
        pp3 = np.array([float(x) for x in line[7:10]])
# debug
for i in CAs:
    print(i)
print(pp1,pp2,pp3)

# get the plane equation
theplane = plane_from_points(pp1,pp2,pp3)
dists = []
ppointsx,ppointsy = [],[]

## get coords for ND CAs
for i in CAs:
    projpoint = proj_point_to_plane(theplane,i)
    ppointsx.append(projpoint[0])
    ppointsy.append(projpoint[1])
    diff = get_dist(i,projpoint)
    dists.append(diff)

## barrel and gate for plotting:

def read_pdb_get_coords(pdbfile,chainl,atoml):
    '''get the xy coordinates of atoms from pdb chain and atomtype comma separated
    returns {chain:[atom#atomtype:[x,y,z]],atom#atomtype:[x,y,z]]}'''
    
    filename = pdbfile.split('/')[-1].split('.')[0]
    filedata = open(pdbfile,'r').readlines()
    if ',' in atoml:
        atoms = atoml.split(',')
    else:
        atoms = [atoml]
    if ',' in chainl:
        chains = chainl.split(',')
    else:
        chains = [chainl]
    
    chaindic = {}            #{chain{atomtypeatomno:[x,y,z]}}
    selected_atoms = []
    for i in filedata:
        line= i.split()
        if len(i) > 25:
            if i[:4] == 'ATOM': 
                chain = i[21]
                atomtype = i[13:15].replace(' ','')
                atomno = int(i[23:26])
                if atomtype in atoms and chain in chains:
                    try:
                        chaindic[chain][str(atomno)+atomtype] =[i[31:38],i[38:47],i[47:55]]
                    except:
                        print('ERROR? Found multiple sets of coords for atom: {0}{1}{2}'.format(chain,atomtype,atomno))
                        chaindic[chain] = {str(atomno)+atomtype:([i[31:38],i[38:47],i[47:55]])}
    return(chaindic)

## get coords for barrel and lateral gate 
pdbcoords = read_pdb_get_coords(sys.argv[2],'A','CA')
barrel = range(735,746)+range(710,721)+range(628,640)+range(608,621)+range(590,601)+range(565,579)+range(523,538)+range(505,515)+range(485,494)+range(467,476)+range(457,463)
gate = range(767,779)+range(781,790)+range(442,446)+range(424,430)
barrelcoordsx,barrelcoordsy,gatecoordsx,gatecoordsy = [],[],[],[]
for i in pdbcoords['A']:
    if int(i.replace('CA','')) in barrel:
        coords = [float(x) for x in pdbcoords['A'][i]] 
        pcoords = proj_point_to_plane(theplane,coords)
        barrelcoordsx.append(pcoords[0])
        barrelcoordsy.append(pcoords[1])
    elif int(i.replace('CA','')) in gate:
        coords = [float(x) for x in pdbcoords['A'][i]] 
        pcoords = proj_point_to_plane(theplane,coords)
        gatecoordsx.append(pcoords[0])
        gatecoordsy.append(pcoords[1])   

## plot
plt.scatter(barrelcoordsx,barrelcoordsy,c='black', s=20, edgecolors='face')
plt.scatter(gatecoordsx,gatecoordsy,c='white',marker='^', s=30, edgecolors='black')
plt.scatter(ppointsx,ppointsy,c=dists,cmap='coolwarm', s=20, edgecolors='face',vmin=-10.0,vmax=10.0)
cbar = plt.colorbar()
cbar.solids.set_edgecolor("face")
filename= sys.argv[1].split('.')[0]
plt.savefig('mspplan_{0}.png'.format(filename))
