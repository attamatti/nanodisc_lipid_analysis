#!/usr/bin/env python


import sys
import numpy as np
import scipy
import matplotlib.pyplot as plt


MSP1_chain = 'F'
MSP2_chain = 'G'
def iLP2(plane,point):
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
    #print (point,'point')
    #print(plane,'plane')
    #print (t,'t')
    xo = (a*t)+x
    yo = (b*t)+y
    zo = (c*t)+z
    return(xo,yo,zo)

def get_distance_to_plane(point, CPX):
    '''point = (x,y,z) CPX = (a,b,c,d)'''
    # project the point on to the centreplane
    cpp = iLP2(CPX,point)
    dist_to_plane = abs(point[0]-cpp[0])+abs(point[1]-cpp[1])+abs(point[2]-cpp[2])
    return([cpp[0],cpp[1]],dist_to_plane)

def plane_from_points(p1,p2,p3):
    '''returns ax+by+cy+d=0'''
    v1 = p3 - p1
    v2 = p2 - p1
    cp = np.cross(v1, v2)
    a, b, c = cp
    d = np.dot(cp, p3)
    e1 = [c-b,a-c,b-a]
    e2 = [a*(b+c)-b**2-c**2,   b*(a+c)-a**2-c**2,   c*(a+b)-a**2-b**2]
    print('{0}x + {1}y + {2}z = {3}'.format(a, b, c, d))
    return(a,b,c,-d)

def read_pdb_get_Ps(pdbfile):
    '''get the xy coordinates of teh two MSP proteins'''
    filename = pdbfile.split('/')[-1].split('.')[0]
    filedata = open(pdbfile,'r').readlines()
    
    cas = {}            #{chain{atomno:[x,y,z]}}
    for i in filedata:
        line= i.split()
        if len(i) > 25:
            if i[:4] == 'ATOM': 
                chain = i[21]
                atomtype = i[13:15]
                atomno = int(i[23:26])
                if atomtype == 'CA':
                    try:
                        cas[chain][atomno] =[i[31:38],i[38:47],i[47:55]]
                    except:
                        cas[chain] = {atomno:([i[31:38],i[38:47],i[47:55]])}
    ## get the 2 MSPs
    MSP1Cas,MSP2Cas = [],[]         # just the x and y coords of all MSP CAs
    for chain in cas:
        if chain == MSP1_chain:
            for aa in cas[chain]:
                #print(chain)
                #print(aa)
                MSP1Cas.append([float(x) for x in cas[chain][aa]])
        elif chain == MSP2_chain:
            for aa in cas[chain]:
                #print(chain)
                print(aa)
                MSP2Cas.append([float(x) for x in cas[chain][aa]])
    return(MSP1Cas,MSP2Cas)

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
    #
    #print("Solution: ", sol)
    #print("Old Error: ", (f_min(XYZ, p0)**2).sum())
    #print("New Error: ", (f_min(XYZ, sol)**2).sum())

    a,b,c,d = sol[0],sol[1],sol[2],sol[3]    
    print('{0}x + {1}y + {2}z + {3} = 0'.format(a, b, c, d))
    
    ## calculate plane normal:
    PNx = a
    PNy = b
    PNz = c
    
    print ('plane normal',PNx,PNy,PNz)

    #calculate plane vectors
    e1 = [c-b,a-c,b-a]
    e2 = [a*(b+c)-b**2-c**2,   b*(a+c)-a**2-c**2,   c*(a+b)-a**2-b**2]
    print ('e1',e1)
    print ('e2',e2)
    
    # get unit vector for plane vector
    scale = 100.0
    e1mag = np.sqrt((e1[0]**2)+(e1[1]**2)+(e1[2]**2))
    e1 = [(x/e1mag)*scale for x in e1]

    e2mag = np.sqrt((e2[0]**2)+(e2[1]**2)+(e2[2]**2))
    e2 = [(x/e2mag)*scale for x in e2]

    return([PNx,PNy,PNz],e1,e2,(a,b,c,d))

def fit_MSP_plane(MSP_coords,mspname):
    '''fit a plane to the MSP coordinates - starting with a plane perpendicular to the z axis as 1st guess'''
    bildout = open('{0}_diag.bild'.format(mspname),'w')
    for ca in MSP_coords:
        bildout.write('.sphere {0} {1} {2} 1.0\n'.format(ca[0],ca[1],ca[2]))
    meanpoint = [np.mean([x[0] for x in MSP_coords]),np.mean([x[1] for x in MSP_coords]),np.mean([x[2] for x in MSP_coords])]
    startingplane = plane_from_points(np.array([meanpoint[0]+20,meanpoint[1]+20,meanpoint[2]]),np.array([meanpoint[0]-20,meanpoint[1]-20,meanpoint[2]]),np.array([meanpoint[0]-10,meanpoint[1]-15,meanpoint[2]]))
    print(startingplane)
    xyz = np.array([[x[0] for x in MSP_coords],[x[1] for x in MSP_coords],[x[2] for x in MSP_coords]])
    PN,e1,e2,abcde = fit_plane(startingplane,xyz)
    ## get the distances to the centre plane
    tpp1 = [meanpoint[0]+e1[0],meanpoint[1]+e1[1],meanpoint[2]+e1[2]]
    tpp3 = [meanpoint[0]-e1[0],meanpoint[1]-e1[1],meanpoint[2]-e1[2]]
    tpp2 = [meanpoint[0]+e2[0],meanpoint[1]+e2[1],meanpoint[2]+e2[2]]
    tpp4 = [meanpoint[0]-e2[0],meanpoint[1]-e2[1],meanpoint[2]-e2[2]]
    bildout.write('.polygon {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11}\n'.format(tpp1[0],tpp1[1],tpp1[2],tpp2[0],tpp2[1],tpp2[2],tpp3[0],tpp3[1],tpp3[2],tpp4[0],tpp4[1],tpp4[2]))
    newCP = plane_from_points(np.array(tpp1),np.array(tpp2),np.array(tpp3)) 
    dists = [get_distance_to_plane(x,newCP) for x in MSP_coords]

    
    bildout.close()
    
    
    
    return(dists)

    
def do_it(file):
    MSP1,MSP2 = read_pdb_get_Ps(file)
    dists = fit_MSP_plane(MSP1,'MSP1')
    xs1,ys1,dist1 =[],[],[]
    for ca in dists:
        xs1.append(ca[0][0])
        ys1.append(ca[0][1])
        dist1.append(ca[1])
    #plt.scatter(xs,ys,c=dist)
    #plt.show()
    #plt.close()
    
    dists = fit_MSP_plane(MSP2,'MSP2')
    xs2,ys2,dist2 =[],[],[]
    for ca in dists:
        xs2.append(ca[0][0])
        ys2.append(ca[0][1])
        dist2.append(ca[1])
    #plt.scatter(xs,ys,c=dist)
    #plt.show()
    #plt.close()
    return(xs1,ys1,dist1,xs2,ys2,dist2)

fxs1,fx2s,fys1,fys2,fd1,fd2 = [],[],[],[],[],[]
for i in sys.argv[1:]:
    xs1,ys1,dist1,xs2,ys2,dist2 = do_it(i)
    fxs1.append(xs1)
    fx2s.append(xs2)
    fys1.append(ys1)
    fys2.append(ys2)
    fd1.append(dist1)
    fd2.append(dist2)
xzip1 = zip(*fxs1)
yzip1 = zip(*fys1)
dzip1 = zip(*fd1)


meansx1 = [np.mean(x) for x in xzip1]
meansy1 = [np.mean(x) for x in yzip1]
sumd1 = [sum(x) for x in dzip1]

plt.scatter(meansx1,meansy1,c=sumd1)
plt.show()