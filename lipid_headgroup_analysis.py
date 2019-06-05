#!/usr/bin/env python

# to do:
# clean up screen barf
# claculate and graph angles between leaflet planes

#### fix lazy indexing on xyz coordinates!!!!!!!
#### fix projection of spots on to leaflet plane rather than xy axis

import sys
import numpy as np
import scipy
import matplotlib.pyplot as plt
import os
import subprocess

def intersect_line_plabe(PNV,plane,point):
    '''PNV = [i,j,k], plane = [a,b,c,d], point = [x,y,z]
        ix = x - i(Ax+By+Cz+D)/Ai+Bj+Ck
        iy = y - j(Ax+By+Cz+D)/Ai+Bj+Ck'''
    print(PNV)
    print(plane)
    print(point)
    ix = point[0] - PNV[0]*(plane[0]*point[0]+plane[1]*point[1]+plane[2]*point[2]+plane[3])/plane[0]*PNV[0]+plane[1]*PNV[1]+plane[2]*PNV[2]
    iy = point[1] - PNV[1]*(plane[0]*point[0]+plane[1]*point[1]+plane[2]*point[2]+plane[3])/plane[0]*PNV[0]+plane[1]*PNV[1]+plane[2]*PNV[2]
    print(ix,iy)
    return((ix,iy))
    
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
def read_pdb_get_Ps(pdbfile):
    filedata = open(pdbfile,'r').readlines()
    bildout = open('bildfiles/HG_{0}.bild'.format(pdbfile.split('.')[0]),'w')
    coldic = {'PVPG':'red','PVCL':'yellow','PVPE':'blue'}
    HGcoords = []
    
    # plot the lipid headgroups write pdb and bild
    cas = {}            #{chain{atomno:[x,y,z]}}
    for i in filedata:
        line= i.split()
        if len(i) > 25:
            if i[17:19] =='PV' and i[13] == 'P':
                bildout.write('.color {0} \n.sphere {1} 1.0\n'.format(coldic[i[17:21]],i[31:56]))
                HGcoords.append([float(x) for x in i[31:56].split()])
            # get cAs for strands
            if i[:4] == 'ATOM': 
                chain = i[21]
                atomtype = i[13:15]
                atomno = int(i[23:26])
                if atomtype == 'CA':
                    try:
                        cas[chain][atomno] = (i[31:56].split())
                    except:
                        cas[chain] = {atomno:(i[31:56].split())}
    # plot the lipids in bildfiles
    xs = [float(x[0]) for x in HGcoords]
    ys = [float(x[1]) for x in HGcoords]
    zs = [float(x[2]) for x in HGcoords]
    meanpoint = [np.mean(xs),np.mean(ys),np.mean(zs)]
    bildout.write('.color orange \n.sphere {0} {1} {2} 1.0\n'.format(meanpoint[0],meanpoint[1],meanpoint[2]))
    bildout.write('.color red \n.cylinder {0} {1} {2} {3} {4} {5} 0.5\n'.format(meanpoint[0],meanpoint[1],meanpoint[2],meanpoint[0]+20,meanpoint[1],meanpoint[2]))
    bildout.write('.color blue \n.cylinder {0} {1} {2} {3} {4} {5} 0.5\n'.format(meanpoint[0],meanpoint[1],meanpoint[2],meanpoint[0],meanpoint[1]+20,meanpoint[2]))
    bildout.write('.color yellow  \n.cylinder {0} {1} {2} {3} {4} {5} 0.5\n'.format(meanpoint[0],meanpoint[1],meanpoint[2],meanpoint[0],meanpoint[1],meanpoint[2]+20))
    
    #write the cas to the 3d plotting file
    
    ##### define the strands to draw here ################
    strands = [             [range(1,171),'A']              ]
    ######################################################
    
    strdic = {}             #{strandno:[[x,y,z],[x,y,z], ..., [x,y,z]]}
    strcoords = {}
    for chain in cas:
        for atom in cas[chain]:
            strandn = 1    
            for strand in strands:
                if atom in strand[0] and chain == strand[1]:
                    print(strandn,atom,cas[chain][atom])
                    try:
                    	strcoords[strandn].append(cas[chain][atom])
                    except:
                        strcoords[strandn] = [cas[chain][atom]]
                strandn+=1
    print strcoords
    # write the strand data to the 3d plot file
    pltout = open('plotdata/{0}_plot.txt'.format(pdbfile.split('.')[0]),'w')
    for i in strcoords:
        pltout.write('\nstrand {0}\n'.format(i))
        for j in strcoords[i]:
            pltout.write('{0} '.format(','.join(j)))

    
    # inital guess at the starting plane three points are true mean and true mean +20x +20y and true mean -20x and -20y
    # this is only vaild for matt's aligned BAM ND structures with Bam roughly aligned with the zaxis perpendicular to the Nanodisc
    
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
    
    # write bottomdata for 3d plotting

    meanpoint = [np.mean(txs),np.mean(tys),np.mean(tzs)]
    startingplane = plane_from_points(np.array(meanpoint),np.array(top[0]),np.array(top[1]))
    xyz = np.array([[x[0] for x in top],[x[1] for x in top],[x[2] for x in top]])
    tPN,te1,te2,tabcd = fit_plane(startingplane,xyz)
    tpp1 = [meanpoint[0]+te1[0],meanpoint[1]+te1[1],meanpoint[2]+te1[2]]
    tpp3 = [meanpoint[0]-te1[0],meanpoint[1]-te1[1],meanpoint[2]-te1[2]]
    tpp2 = [meanpoint[0]+te2[0],meanpoint[1]+te2[1],meanpoint[2]+te2[2]]
    tpp4 = [meanpoint[0]-te2[0],meanpoint[1]-te2[1],meanpoint[2]-te2[2]]
    tbbild.write('.polygon {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11}\n'.format(tpp1[0],tpp1[1],tpp1[2],tpp2[0],tpp2[1],tpp2[2],tpp3[0],tpp3[1],tpp3[2],tpp4[0],tpp4[1],tpp4[2]))
    
    # make plot for bottom leaflet originally labelled top:

        
    tpltxs,tpltys,tdistances = [],[],[]
    # project points onto the plane and save for drawing
    for i in top:
        projxy = intersect_line_plabe(tPN,tabcd,i)
        tpltxs.append(projxy[0])
        tpltys.append(projxy[1])
        dist = distance_to_plane(tabcd,i)
        check = ((i[0]-tpp1[0])*tPN[0])+((i[1]-tpp1[1])*tPN[1])+((i[2]-tpp1[2])*tPN[2])
        if check >= 0:
            tdistances.append(dist)
        if check < 0:
            tdistances.append(-dist)
    #normalize distance to 0-1 scale
    distnormval = max([abs(x) for x in tdistances])
    tdistances = [x/distnormval for x in tdistances]
    
    # get values for 3D plot 
    pltout.write('\nbottom\n{0}\n{1}\n{2}\n{3}'.format(' '.join([str(x) for x in txs]),' '.join([str(x) for x in tys]),' '.join([str(x) for x in tzs]),' '.join([str(x) for x in tdistances])))

# fit a plane to the top leaflet originally labelled bottom:
    bxs = [float(x[0]) for x in bottom]
    bys = [float(x[1]) for x in bottom]
    bzs = [float(x[2]) for x in bottom]
    meanpoint = [np.mean(bxs),np.mean(bys),np.mean(bzs)]
    startingplane = plane_from_points(np.array(meanpoint),np.array(bottom[0]),np.array(bottom[1]))
    xyz = np.array([[x[0] for x in bottom],[x[1] for x in bottom],[x[2] for x in bottom]])
    bPN,be1,be2,babcd = fit_plane(startingplane,xyz)
    bpp1 = [meanpoint[0]+be1[0],meanpoint[1]+be1[1],meanpoint[2]+be1[2]]
    bpp3 = [meanpoint[0]-be1[0],meanpoint[1]-be1[1],meanpoint[2]-be1[2]]
    bpp2 = [meanpoint[0]+be2[0],meanpoint[1]+be2[1],meanpoint[2]+be2[2]]
    bpp4 = [meanpoint[0]-be2[0],meanpoint[1]-be2[1],meanpoint[2]-be2[2]]
    tbbild.write('.polygon {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11}\n'.format(bpp1[0],bpp1[1],bpp1[2],bpp2[0],bpp2[1],bpp2[2],bpp3[0],bpp3[1],bpp3[2],bpp4[0],bpp4[1],bpp4[2]))

    # make plot for top leaflet originally labelled bottom:
    bpltxs,bpltys,bdistances = [],[],[]
    for i in bottom:
        projxy = intersect_line_plabe(bPN,babcd,i)
        bpltxs.append(projxy[0])
        bpltys.append(projxy[1])
        dist = distance_to_plane(babcd,i)
        check = ((i[0]-bpp1[0])*bPN[0])+((i[1]-bpp1[1])*bPN[1])+((i[2]-bpp1[2])*bPN[2])
        if check >= 0:
            bdistances.append(dist)
        if check < 0:
            bdistances.append(-dist)
    distnormval = max([abs(x) for x in bdistances])
    bdistances = [x/distnormval for x in bdistances]
    
    
    # calculat the xy limits
    tcent = (np.mean(tpltxs),np.mean(tpltys))
    
    txmin,txmax = min(tpltxs),max(tpltxs)
    txrange = max([abs(tcent[0]-txmin),abs(tcent[0]-txmax)])
    txrmin,txrmax = tcent[0]-txrange-20,tcent[0]+txrange+20
    txr = txrmax - txrmin
    
    tymin,tymax = min(tpltys),max(tpltys)
    tyrange = max([abs(tcent[1]-tymin),abs(tcent[1]-tymax)])
    tyrmin,tyrmax = tcent[1]-tyrange-20,tcent[1]+tyrange+20
    tyr =  tyrmax - tyrmin
    
    bcent = (np.mean(bpltxs),np.mean(bpltys))
    
    bxmin,bxmax = min(bpltxs),max(bpltxs)
    bxrange = max([abs(bcent[0]-bxmin),abs(bcent[0]-bxmax)])
    bxrmin,bxrmax = bcent[0]-bxrange-20,bcent[0]+bxrange+20
    bxr = bxrmax - bxrmin
    
    bymin,bymax = min(bpltys),max(bpltys)
    byrange = max([abs(bcent[1]-bymin),abs(bcent[1]-bymax)])
    byrmin,byrmax = bcent[1]-byrange-20,bcent[1]+byrange+20
    byr =  byrmax - byrmin  
    
    xr = 0.5*max([txr,bxr])
    yr = 0.5*max([tyr,byr])
        
    # plot the top (called bottom)
    plt.xlim(bcent[0]-xr,bcent[0]+xr)
    plt.ylim(bcent[1]-yr,bcent[1]+yr)
    plt.title('top_{0}'.format(pdbfile.split('.')[0]))
    plt.scatter(bpltxs,bpltys,c=bdistances,cmap='coolwarm')
    plt.savefig('top/top_{0}.png'.format(pdbfile.split('.')[0]))
    plt.axis('off')
    plt.close()
    
    # plot the bottom (called top - top and bottom are arbitrary anyway!)
    plt.xlim(tcent[0]-xr,tcent[0]+xr)
    plt.ylim(tcent[1]-yr,tcent[1]+yr)
    plt.title('bottom_{0}'.format(pdbfile.split('.')[0]))
    plt.scatter(tpltxs,tpltys,c=tdistances,cmap='coolwarm')
    plt.savefig('bottom/bottom_{0}.png'.format(pdbfile.split('.')[0]))
    plt.axis('off')
    plt.close()

    bildout.close()
    
    # make a file to write data from 3D plotting
    pltout.write('\ntop\n{0}\n{1}\n{2}\n{3}'.format(' '.join([str(x) for x in bxs]),' '.join([str(x) for x in bys]),' '.join([str(x) for x in bzs]),' '.join([str(x) for x in bdistances])))

# make the necessary directories -- keep shit organised
if os.path.isdir('top') == False:
    subprocess.call(['mkdir','top'])

if os.path.isdir('bottom') == False:
    subprocess.call(['mkdir','bottom'])

if os.path.isdir('bildfiles') == False:
    subprocess.call(['mkdir','bildfiles'])

if os.path.isdir('plotdata') == False:
    subprocess.call(['mkdir','plotdata'])

### DO IT!!!
for i in sys.argv[1:]:
    read_pdb_get_Ps(i)
