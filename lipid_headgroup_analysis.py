#!/usr/bin/env python

# to do:
# clean up screen barf
# claculate and graph angles between leaflet planes

import sys
import numpy as np
import scipy
import matplotlib.pyplot as plt
import os
import subprocess

def distance_to_plane(planeeq,point):
    '''using abs(ax+by+cz+d)/sqrt(a^2+b^2+c^2)'''
    
    a,b,c,d = planeeq[0],planeeq[1],planeeq[2],planeeq[3]
    ax,by,cz = point[0]*a,point[1]*b,point[2]*c
    dist = abs(ax+by+cz+d)/np.sqrt(a**2+b**2+c**2)
    return(dist)

def plane_from_points(p1,p2,p3):
    '''returns ax+by+cy+d=0'''
    v1 = p3 - p1
    v2 = p2 - p1
    cp = np.cross(v1, v2)
    a, b, c = cp
    d = np.dot(cp, p3)
    
    print('{0}x + {1}y + {2}z = {3}'.format(a, b, c, d))
    return(a,b,c,d)

def fit_plane(initial,XYZ):
    '''returns plane normal vector and two parallel vectors'''
    p0 = [0.506645455682, -0.185724560275, -1.43998120646, 1.37626378129]
    
    def f_min(X,p):
        plane_xyz = p[0:3]
        distance = (plane_xyz*X.T).sum(axis=1) + p[3]
        return distance / np.linalg.norm(plane_xyz)
    
    def residuals(params, signal, X):
        return f_min(X, params)
    
    from scipy.optimize import leastsq
    sol = leastsq(residuals, p0, args=(None, XYZ))[0]
    
    print("Solution: ", sol)
    print("Old Error: ", (f_min(XYZ, p0)**2).sum())
    print("New Error: ", (f_min(XYZ, sol)**2).sum())

    a,b,c,d = sol[0],sol[1],sol[2],sol[3]    
    print('{0}x + {1}y + {2}z = {3}'.format(a, b, c, d))
    
    ## calculate plane normal:
    PNx = a/(np.sqrt(2**2+b**2+c**2))
    PNy = b/(np.sqrt(2**2+b**2+c**2))
    PNz = c/(np.sqrt(2**2+b**2+c**2))
    
    print ('plane normal',PNx,PNy,PNz)

    #calculate plane vectors
    e1 = [c-b,a-c,b-a]
    e2 = [a*(b+c)-b**2-c**2,   b*(a+c)-a**2-c**2,   c*(a+b)-a**2-b**2]
    print ('e1',e1)
    print ('e2',e2)
    
    # get unit vector for plane vector
    scale = 100.0
    e1mag = np.sqrt(e1[0]**2+e1[1]**2+e1[2]**2)
    e1 = [(x/e1mag)*scale for x in e1]

    e2mag = np.sqrt(e2[0]**2+e2[1]**2+e2[2]**2)
    e2 = [(x/e2mag)*scale for x in e2]

    return([PNx,PNy,PNz],e1,e2,(a,b,c,d))
    
# DO IT!
def read_pdb_get_Ps(pdbfile,xmin,xmax,ymin,ymax):
    filedata = open(pdbfile,'r').readlines()
    bildout = open('bildfiles/HG_{0}.bild'.format(pdbfile.split('.')[0]),'w')
    coldic = {'PVPG':'red','PVCL':'yellow','PVPE':'blue'}
    HGcoords = []
    
    # plot the lipid headgroups write pdb and bild
    for i in filedata:
        line= i.split()
        if len(i) > 25:
            if i[17:19] =='PV' and i[13] == 'P':
                bildout.write('.color {0} \n.sphere {1} 1.0\n'.format(coldic[i[17:21]],i[31:56]))
                HGcoords.append([float(x) for x in i[31:56].split()])
    xs = [float(x[0]) for x in HGcoords]
    ys = [float(x[1]) for x in HGcoords]
    zs = [float(x[2]) for x in HGcoords]
    meanpoint = [np.mean(xs),np.mean(ys),np.mean(zs)]
    bildout.write('.color orange \n.sphere {0} {1} {2} 1.0\n'.format(meanpoint[0],meanpoint[1],meanpoint[2]))
    bildout.write('.color red \n.cylinder {0} {1} {2} {3} {4} {5} 0.5\n'.format(meanpoint[0],meanpoint[1],meanpoint[2],meanpoint[0]+20,meanpoint[1],meanpoint[2]))
    bildout.write('.color blue \n.cylinder {0} {1} {2} {3} {4} {5} 0.5\n'.format(meanpoint[0],meanpoint[1],meanpoint[2],meanpoint[0],meanpoint[1]+20,meanpoint[2]))
    bildout.write('.color yellow  \n.cylinder {0} {1} {2} {3} {4} {5} 0.5\n'.format(meanpoint[0],meanpoint[1],meanpoint[2],meanpoint[0],meanpoint[1],meanpoint[2]+20))
    
    

    # inital guess at the starting plane three points are true mean and true mean +20x +20y and true mean -20x and -20y
    
    startingplane = plane_from_points(np.array([meanpoint[0]+20,meanpoint[1]+20,meanpoint[2]]),np.array([meanpoint[0]-20,meanpoint[1]-20,meanpoint[2]]),np.array(meanpoint))

    # calculate and draw plane
    xyz = np.array([[x[0] for x in HGcoords],[x[1] for x in HGcoords],[x[2] for x in HGcoords]])
    PN,e1,e2,abcde = fit_plane(startingplane,xyz)
    pp1 = [meanpoint[0]+e1[0],meanpoint[1]+e1[1],meanpoint[2]+e1[2]]
    pp3 = [meanpoint[0]-e1[0],meanpoint[1]-e1[1],meanpoint[2]-e1[2]]
    pp2 = [meanpoint[0]+e2[0],meanpoint[1]+e2[1],meanpoint[2]+e2[2]]
    pp4 = [meanpoint[0]-e2[0],meanpoint[1]-e2[1],meanpoint[2]-e2[2]]
    bildout.write('.polygon {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11}\n'.format(pp1[0],pp1[1],pp1[2],pp2[0],pp2[1],pp2[2],pp3[0],pp3[1],pp3[2],pp4[0],pp4[1],pp4[2]))
    
    # find points above/below plane with dot product of point-plane point and plane normal vector
    top,bottom,ambig = [],[],[]
    for i in HGcoords:
        check = ((i[0]-pp1[0])*PN[0])+((i[1]-pp1[1])*PN[1])+((i[2]-pp1[2])*PN[2])
        if check > 0.2:
            top.append(i)
        elif check < -0.2:
            bottom.append(i)
        else:
            ambig.append(i)
            
    tbbild = open('bildfiles/TB_{0}.bild'.format(pdbfile.split('.')[0]),'w')
    tbbild.write('.color blue\n')
    for i in top:
        tbbild.write('.sphere {0} {1} {2} 1.0\n'.format(i[0],i[1],i[2]))
    tbbild.write('.color yellow\n')
    for i in bottom:
        tbbild.write('.sphere {0} {1} {2} 1.0\n'.format(i[0],i[1],i[2]))
    tbbild.write('.color red\n')
    for i in ambig:
        tbbild.write('.sphere {0} {1} {2} 1.0\n'.format(i[0],i[1],i[2]))
    
    # fit a plane to the bottom leaflet originall labelled top:
    txs = [float(x[0]) for x in top]
    tys = [float(x[1]) for x in top]
    tzs = [float(x[2]) for x in top]
    meanpoint = [np.mean(txs),np.mean(tys),np.mean(tzs)]
    startingplane = plane_from_points(np.array(meanpoint),np.array(top[0]),np.array(top[1]))
    xyz = np.array([[x[0] for x in top],[x[1] for x in top],[x[2] for x in top]])
    tPN,te1,te2,tabcde = fit_plane(startingplane,xyz)
    tpp1 = [meanpoint[0]+te1[0],meanpoint[1]+te1[1],meanpoint[2]+te1[2]]
    tpp3 = [meanpoint[0]-te1[0],meanpoint[1]-te1[1],meanpoint[2]-te1[2]]
    tpp2 = [meanpoint[0]+te2[0],meanpoint[1]+te2[1],meanpoint[2]+te2[2]]
    tpp4 = [meanpoint[0]-te2[0],meanpoint[1]-te2[1],meanpoint[2]-te2[2]]
    tbbild.write('.polygon {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11}\n'.format(tpp1[0],tpp1[1],tpp1[2],tpp2[0],tpp2[1],tpp2[2],tpp3[0],tpp3[1],tpp3[2],tpp4[0],tpp4[1],tpp4[2]))
    
    # make plot for bottom leaflet originall labelled top:
    pltxs,pltys,distances = [],[],[]
    for i in top:
        pltxs.append(i[0])
        pltys.append(i[1])
        dist = distance_to_plane(tabcde,i)
        check = ((i[0]-tpp1[0])*tPN[0])+((i[1]-tpp1[1])*tPN[1])+((i[2]-tpp1[2])*tPN[2])
        if check >= 0:
            distances.append(dist)
        if check < 0:
            distances.append(-dist)
    distnormval = max([abs(x) for x in distances])
    distances = [x/distnormval for x in distances]
    #get values for plot 
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    plt.title('bottom_{0}'.format(pdbfile.split('.')[0]))
    plt.scatter(pltxs,pltys,c=distances,cmap='coolwarm')
    plt.savefig('bottom/bottom_{0}.png'.format(pdbfile.split('.')[0]))
    plt.close()
    
# fit a plane to the top leaflet originally labelled bottom:
    bxs = [float(x[0]) for x in bottom]
    bys = [float(x[1]) for x in bottom]
    bzs = [float(x[2]) for x in bottom]
    meanpoint = [np.mean(bxs),np.mean(bys),np.mean(bzs)]
    startingplane = plane_from_points(np.array(meanpoint),np.array(bottom[0]),np.array(bottom[1]))
    xyz = np.array([[x[0] for x in bottom],[x[1] for x in bottom],[x[2] for x in bottom]])
    bPN,be1,be2,babcde = fit_plane(startingplane,xyz)
    bpp1 = [meanpoint[0]+be1[0],meanpoint[1]+be1[1],meanpoint[2]+be1[2]]
    bpp3 = [meanpoint[0]-be1[0],meanpoint[1]-be1[1],meanpoint[2]-be1[2]]
    bpp2 = [meanpoint[0]+be2[0],meanpoint[1]+be2[1],meanpoint[2]+be2[2]]
    bpp4 = [meanpoint[0]-be2[0],meanpoint[1]-be2[1],meanpoint[2]-be2[2]]
    tbbild.write('.polygon {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11}\n'.format(bpp1[0],bpp1[1],bpp1[2],bpp2[0],bpp2[1],bpp2[2],bpp3[0],bpp3[1],bpp3[2],bpp4[0],bpp4[1],bpp4[2]))

    # make plot for top leaflet originally labelled bottom:
    pltxs,pltys,distances = [],[],[]
    for i in bottom:
        pltxs.append(i[0])
        pltys.append(i[1])
        dist = distance_to_plane(babcde,i)
        check = ((i[0]-bpp1[0])*bPN[0])+((i[1]-bpp1[1])*bPN[1])+((i[2]-bpp1[2])*bPN[2])
        if check >= 0:
            distances.append(dist)
        if check < 0:
            distances.append(-dist)
    distnormval = max([abs(x) for x in distances])
    distances = [x/distnormval for x in distances]
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    plt.title('top_{0}'.format(pdbfile.split('.')[0]))
    plt.scatter(pltxs,pltys,c=distances,cmap='coolwarm')
    plt.savefig('top/top_{0}.png'.format(pdbfile.split('.')[0]))
    plt.close()

    bildout.close()
    

# make the necessary directories -- keep shit organised
if os.path.isdir('top') == False:
    subprocess.call(['mkdir','top'])

if os.path.isdir('bottom') == False:
    subprocess.call(['mkdir','bottom'])

if os.path.isdir('bildfiles') == False:
    subprocess.call(['mkdir','bildfiles'])

# get the x and y lims
limxs = []
limys = []
for i in sys.argv[1:]:
    filedata = open(i,'r').readlines()
    for i in filedata:
        line= i.split()
        if len(i) > 25:
            if i[17:19] =='PV' and i[13] == 'P':
                limxs.append(float(i[31:56].split()[0]))
                limys.append(float(i[31:56].split()[1]))
xmin,xmax = min(limxs)-20,max(limxs)+20
ymin,ymax = min(limys)-20,max(limys)+20

### DO IT!!!
for i in sys.argv[1:]:
    read_pdb_get_Ps(i,xmin,xmax,ymin,ymax)
