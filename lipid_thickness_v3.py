#!/usr/bin/env python

# to do:
# clean up screen barf
# make all plots have normalized color values - has a workaround in right now
# need to write out grid data to file then write script to parse all of the output files
# and redraw plots

# top and bottom leaflets are switched relative to BAM, but this designation is arbirtary anyways

#########  CAs to draw on the final map  ############
## for BAM
#MSP_chains = ['F','G']
#barreldraw = {'barrel':['A',range(436,794)],'lateralgate':['A',[423,424,425,426,427,428,429,430,421,432,433,434,435,436,437,794,795,796,797,798,799,800,801,802,803,804,805,806,807,808,809,810]]}            # {Name1:[Chain[AAs]],Name2:[Chain2,[AAs2]]}
#draworder = ['barrel', 'lateralgate']

###for tOmpA
#MSP_chains = ['B','C']
#barreldraw = {'barrel':['A',range(0,170)]}            # {Name1:[Chain[AAs]],Name2:[Chain2,[AAs2]]}
#draworder = ['barrel']


###for tOmpA simulationframes
MSP_chains = ['X']
barreldraw = {'barrel':['X',range(0,171)]}            # {Name1:[Chain[AAs]],Name2:[Chain2,[AAs2]]}
draworder = ['barrel']


#####################################################


import sys
import numpy as np
import scipy
import matplotlib.pyplot as plt
import os
import subprocess



def distance(xy1,xy2):      # xy1 and xy2 = (x,y) 
    '''distance between two points in 2D'''
    distance = abs(xy1[0]-xy2[0])+abs(xy1[1]-xy2[1])
    return (distance)

def find_nn(point,coords): # point = (x,y) coords = [[(x,y),dist],[(x,y),dist],[(x,y),dist]]
    distdic = {}
    for i in coords:
        try:
            distdic[distance(point,i[0])].append(i)
        except:
            distdic[distance(point,i[0])] = [i]
    distances = distdic.keys()
    distances.sort()
    nns = []    
    for i in distances:
        for j in distdic[i]:
            nns.append(j[1])
        if len(nns) == 3:
            break    
    mh = np.mean([x for x in nns])
    # print('point',point,mh)
    # print('nns',nns)
    # print('mean height', mh)
    return(mh,point[0],point[1])
    
def find_nn_height(point,coords):       # point = (x,y) coords = [(x,y,z),(x,y,z),(x,y,z)]
    '''given a point on a grid determine the height of the leaflet there as mean of three nearest neighbor points'''
    distdic = {}
    for i in coords:
        try:
            distdic[distance(point,i[0:2])].append(i)
        except:
            distdic[distance(point,i[0:2])] = [i]
        print(i,point,distance(point,i[0:2]))
    distances = distdic.keys()
    distances.sort()
    nns = []
    for i in distances:
        for j in distdic[i]:
            nns.append(j)
        if len(nns) == 3:
            break
    mh = np.mean([x[2] for x in nns])
    print('point',point,mh)
    print('nns',nns)
    print('mean height', mh)
    return(mh,point[0],point[1])
    
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
    print (point,'point')
    print(plane,'plane')
    print (t,'t')
    xo = (a*t)+x
    yo = (b*t)+y
    zo = (c*t)+z
    diag.write('.color green\n')
    diag.write('.sphere {0} {1} {2} 0.5\n'.format(x,y,z))
    diag.write('.color grey\n')
    diag.write('.sphere {0} {1} {2} 0.5\n'.format(xo,yo,zo))
    return(xo,yo,zo)
    
def get_distance_from_centreplane(point, CPX):
    '''point = (x,y,z) CPX = (a,b,c,d)'''
    # project the point on to the centreplane
    cpp = iLP2(CPX,point)
    diag.write('.arrow {0} {1} {2} {3} {4} {5} 0.2 0.6 \n'.format(point[0],point[1],point[2],cpp[0],cpp[1],cpp[2]))
    dist_to_plane = abs(point[0]-cpp[0])+abs(point[1]-cpp[1])+abs(point[2]-cpp[2])
    return([cpp[0],cpp[1]],dist_to_plane)

def make_leaflet_grid(points):
    '''points = ([x,y,dist],[x,y,dist])'''
    xs = [x[0] for x in points]
    ys = [x[1] for x in points]
    dists = [x[2] for x in points]

def plane_from_points(p1,p2,p3):
    '''returns ax+by+cy+d=0'''
    v1 = p3 - p1
    v2 = p2 - p1
    cp = np.cross(v1, v2)
    a, b, c = cp
    d = np.dot(cp, p3)
    
    print('{0}x + {1}y + {2}z = {3}'.format(a, b, c, d))
    return(a,b,c,-d)

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
    
# DO IT!
def read_pdb_get_Ps(pdbfile):
    filename = pdbfile.split('/')[-1].split('.')[0]
    filedata = open(pdbfile,'r').readlines()
    bildout = open('bildfiles/HG_{0}.bild'.format(filename),'w')
    coldic = {'PVPG':'red','PVCL':'yellow','PVPE':'blue'}
    HGcoords = []
    
    # plot the lipid headgroups write pdb and bild
    cas = {}            #{chain{atomno:[x,y,z]}}
    for i in filedata:
        line= i.split()
        if len(i) > 25:
            if i[17:19] =='PV' and i[13] == 'P':
                bildout.write('.color {0} \n.sphere {1} 1.0\n'.format(coldic[i[17:21]],i[31:56]))
                #HGcoords.append([float(x) for x in i[31:56].split()])
                HGcoords.append([float(i[31:38]),float(i[38:47]),float(i[47:55])])
            # get cAs for strands
            if i[:4] == 'ATOM': 
                chain = i[21]
                atomtype = i[13:15]
                atomno = int(i[23:26])
                if atomtype == 'CA':
                    try:
                        cas[chain][atomno] =[i[31:38],i[38:47],i[47:55]]
                    except:
                        cas[chain] = {atomno:([i[31:38],i[38:47],i[47:55]])}
    # plot the lipids in bildfiles
    xs = [float(x[0]) for x in HGcoords]
    ys = [float(x[1]) for x in HGcoords]
    zs = [float(x[2]) for x in HGcoords]
    meanpoint = [np.mean(xs),np.mean(ys),np.mean(zs)]
    bildout.write('.color orange \n.sphere {0} {1} {2} 1.0\n'.format(meanpoint[0],meanpoint[1],meanpoint[2]))
    bildout.write('.color red \n.cylinder {0} {1} {2} {3} {4} {5} 0.5\n'.format(meanpoint[0],meanpoint[1],meanpoint[2],meanpoint[0]+20,meanpoint[1],meanpoint[2]))
    bildout.write('.color blue \n.cylinder {0} {1} {2} {3} {4} {5} 0.5\n'.format(meanpoint[0],meanpoint[1],meanpoint[2],meanpoint[0],meanpoint[1]+20,meanpoint[2]))
    bildout.write('.color yellow  \n.cylinder {0} {1} {2} {3} {4} {5} 0.5\n'.format(meanpoint[0],meanpoint[1],meanpoint[2],meanpoint[0],meanpoint[1],meanpoint[2]+20))
    
    # get the MSP CAS for drawing:
    
    MSPCas = []         # just the x and y coords of all MSP CAs
    otherdraw = {}
    for chain in cas:
        #print(chain)
        if chain in MSP_chains:
            for aa in cas[chain]:
                print(chain)
                print(aa)
                MSPCas.append([float(x) for x in cas[chain][aa]])
        for i in barreldraw:
            if chain == barreldraw[i][0]:
                for aa in cas[chain]:
                    if aa in barreldraw[i][1]:
                        try:
                            otherdraw[i].append(cas[chain][aa])
                        except:
                            otherdraw[i]=[cas[chain][aa]]

       
    ### inital guess at the starting plane three points are true mean and true mean +20x +20y and true mean -20x and -20y
    ### optimized for matt's aligned BAM ND structures with BAM roughly aligned with the zaxis perpendicular to the Nanodisc
    
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
    
    print(bottom)
    
    
    # write the bild files       
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
    
    # # write bottomdata for plotting
    # 
    meanpoint = [np.mean(txs),np.mean(tys),np.mean(tzs)]
    startingplane = plane_from_points(np.array(meanpoint),np.array(top[0]),np.array(top[1]))
    xyz = np.array([[x[0] for x in top],[x[1] for x in top],[x[2] for x in top]])
    tPN,te1,te2,tabcd = fit_plane(startingplane,xyz)
    tpp1 = [meanpoint[0]+te1[0],meanpoint[1]+te1[1],meanpoint[2]+te1[2]]
    tpp3 = [meanpoint[0]-te1[0],meanpoint[1]-te1[1],meanpoint[2]-te1[2]]
    tpp2 = [meanpoint[0]+te2[0],meanpoint[1]+te2[1],meanpoint[2]+te2[2]]
    tpp4 = [meanpoint[0]-te2[0],meanpoint[1]-te2[1],meanpoint[2]-te2[2]]
    tbbild.write('.polygon {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11}\n'.format(tpp1[0],tpp1[1],tpp1[2],tpp2[0],tpp2[1],tpp2[2],tpp3[0],tpp3[1],tpp3[2],tpp4[0],tpp4[1],tpp4[2]))
    
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

    # get new center plane from average of the two leaflet planes
    newCP = [(tabcd[0]+babcd[0])/2,(tabcd[1]+babcd[1])/2,(tabcd[2]+babcd[2])/2,(tabcd[3]+babcd[3])/2]
    
    # illustrate the new centre plane
    a,b,c,d = (newCP[0],newCP[1],newCP[2],newCP[3])
    e1 = [c-b,a-c,b-a]
    e2 = [a*(b+c)-b**2-c**2,   b*(a+c)-a**2-c**2,   c*(a+b)-a**2-b**2]
    print ('e1',e1)
    print ('e2',e2)
    scale = 100.0
    
    # normalize magnitudes to make a square
    e1mag = np.sqrt(e1[0]**2+e1[1]**2+e1[2]**2)
    e1 = [(x/e1mag)*scale for x in e1]

    e2mag = np.sqrt(e2[0]**2+e2[1]**2+e2[2]**2)
    e2 = [(x/e2mag)*scale for x in e2]
    
    meanpoint = [np.mean(bxs+txs),np.mean(bys+tys),np.mean(bzs+tzs)]
    bpp1 = [meanpoint[0]+e1[0],meanpoint[1]+e1[1],meanpoint[2]+e1[2]]
    bpp3 = [meanpoint[0]-e1[0],meanpoint[1]-e1[1],meanpoint[2]-e1[2]]
    bpp2 = [meanpoint[0]+e2[0],meanpoint[1]+e2[1],meanpoint[2]+e2[2]]
    bpp4 = [meanpoint[0]-e2[0],meanpoint[1]-e2[1],meanpoint[2]-e2[2]]
    
    tbbild.write('.color yellow\n')
    tbbild.write('.polygon {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11}\n'.format(bpp1[0],bpp1[1],bpp1[2],bpp2[0],bpp2[1],bpp2[2],bpp3[0],bpp3[1],bpp3[2],bpp4[0],bpp4[1],bpp4[2]))
    tbbild.write('.color red\n')
    # get the distances to the  new centre plane for each headgroup
    
    newCP = plane_from_points(np.array(bpp1),np.array(bpp2),np.array(bpp3)) 
    pbottom = [get_distance_from_centreplane(x,newCP) for x in bottom]
    ptop = [get_distance_from_centreplane(x,newCP) for x in top] 
    
    # project the MSP CAs onto the center plane 
    MSPs = [iLP2(newCP,x) for x in MSPCas]
    MSPx =[x[0] for x in MSPs]
    MSPy =[x[1] for x in MSPs] 
    
    # do the same projectiosn for the barrel portions to be plotted
    bos = {}        #{name:[[x,y,z],[x,y,z]], name2:[[x,y,z],[x,y,z]]}
    for i in draworder:
        bos[i] = []
        for j in otherdraw[i]:
            bos[i].append(iLP2(newCP,[float(x) for x in j]))
    
    txs = [x[0][0] for x in ptop]
    tys = [x[0][1] for x in ptop]
    tds = [x[1] for x in ptop]
    bxs = [x[0][0] for x in pbottom]
    bys = [x[0][1] for x in pbottom]
    bds = [x[1] for x in pbottom]
    
    # make a square around the bilayer
    txrange = (abs(min(txs)-max(txs)))
    tyrange = (abs(min(tys)-max(tys)))
    meantx = np.mean([min(txs),max(txs)])
    meanty = np.mean([min(tys),max(tys)])
    bxrange = (abs(min(bxs)-max(bxs)))
    byrange = (abs(min(bys)-max(bys)))
    meanbx = np.mean([min(bxs),max(bxs)])
    meanby = np.mean([min(bys),max(bys)])
    txcent = np.mean(txs)
    tycent = np.mean(tys)
    bxcent = np.mean(bxs)
    bycent = np.mean(bys)
    #therange = 0.5*max([txrange,tyrange,bxrange,byrange])+10
    #switched the range toa fixed value so the different frames can be compared
    therange = 75
    # graph the top
    
    print ((txcent-therange)-(txcent+therange))
    print((tycent-therange)-(tycent+therange))
    plt.xlim = (txcent-therange,txcent+therange)
    plt.ylim = (tycent-therange,tycent+therange)
    plt.gca().set_aspect('equal', adjustable='datalim')
    plt.scatter(txs,tys,c=tds)
    #plt.show()
    plt.savefig('1.png')
    
    # graph thebottom
    
    print ((bxcent-therange)-(bxcent+therange))
    print((bycent-therange)-(bycent+therange))
    plt.xlim = (bxcent-therange,bxcent+therange)
    plt.ylim = (bycent-therange,bycent+therange)
    plt.gca().set_aspect('equal', adjustable='datalim')
    plt.scatter(bxs,bys,c=bds)
    #plt.show()
    plt.savefig('2.png')
    print(therange)
    # make the two surfaces
    
    # make the grid in the square - bottom (labelled top)
    tgridx = np.arange(txcent-therange,txcent+therange,1)
    tgridy = np.arange(tycent-therange,tycent+therange,1)     
    xx, yy = np.meshgrid(tgridx, tgridy)
    
    z = [find_nn((x[0],x[1]),ptop)[0] for x in zip(list(xx.flat),list(yy.flat))]
    z=np.array(z)
    z = np.reshape(z,xx.shape)
    h = plt.contourf(tgridx,tgridy,z)
    plt.title('bottom')
    plt.savefig('bottom_plot.png')    
    
    # then bottom
    bz = [find_nn((x[0],x[1]),pbottom)[0] for x in zip(list(xx.flat),list(yy.flat))]
    bz=np.array(bz)
    bz = np.reshape(bz,xx.shape)
    h = plt.contourf(tgridx,tgridy,bz)
    plt.title('top')
    plt.savefig('top_plot.png')   
    thickness = z+bz
    h = plt.contourf(tgridx,tgridy,thickness,vmin=20,vmax=75,cmap='coolwarm')
    cbar = plt.colorbar(h)
    
    ##plot the MSPs
    plt.scatter(MSPx,MSPy,c='k')
    
    ## plot the other barrel components
    markers = ('v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd')
    colors = ['k','w']
    mcount = 0
    ccount = 0
    boses = []
    for i in bos:
        bosx = [x[0] for x in bos[i]]
        bosy = [x[1] for x in bos[i]]
        plt.scatter(bosx,bosy,marker=markers[mcount],c=colors[ccount],edgecolors='k')
        mcount+=1
        ccount +=1
        if mcount > len(markers)-1:
            mcount=0
        if ccount >1:
            ccount=0
        boses.append((bosx,bosy))
    plt.title('thickness {0}'.format(filename))
    plt.savefig('results/{0}_thickness_plot.png'.format(filename))
    #plt.show()
    plt.close()

    print(newCP)
    print(plane_from_points(np.array(bpp1),np.array(bpp2),np.array(bpp3)))
    print(startingplane)
    print(tabcd)
    print(babcd)
    subprocess.call(['mv','diag.bild','bildfiles/{0}_diag.bild'.format(filename)])
    return(thickness,MSPx,MSPy,boses)
##### main ######
# make the necessary directories -- keep shit organised
if os.path.isdir('top') == False:
    subprocess.call(['mkdir','top'])

if os.path.isdir('bottom') == False:
    subprocess.call(['mkdir','bottom'])

if os.path.isdir('bildfiles') == False:
    subprocess.call(['mkdir','bildfiles'])

if os.path.isdir('results') == False:
    subprocess.call(['mkdir','results'])

### DO IT!!!
final_data = []
final_MSPx,final_MSPy = [],[]
final_drawx,final_drawy = [],[]
for i in sys.argv[1:]:
    diag = open('diag.bild','w')
    thick,MSPx,MSPy,boses =  read_pdb_get_Ps(i)
    final_data.append(thick)
    for i in boses:
        final_drawx.append(np.asarray(i[0]))
        final_drawy.append(np.asarray(i[1]))
    final_MSPx.append(np.asarray(MSPx))
    final_MSPy.append(np.asarray(MSPy))
    diag.close()

thickness_mean = np.mean(final_data, axis=0)
thickness_std =  np.std(final_data, axis=0)
drawMSPx = np.mean(final_MSPx,axis=0)
drawMSPy = np.mean(final_MSPy,axis=0)
num_extra_elements = len(final_drawx)/len(final_data)
print(num_extra_elements)


draw_elementsx,draw_elementsy = [],[]
for i in range(0,num_extra_elements):
    draw_elementsx.append(np.mean(final_drawx[i::num_extra_elements],axis=0))
    draw_elementsy.append(np.mean(final_drawy[i::num_extra_elements],axis=0))

## make the summary plots
xcent = np.mean(drawMSPx)
ycent = np.mean(drawMSPy)
gridx = np.arange(xcent-75,xcent+75,1)
gridy = np.arange(ycent-75,ycent+75,1)



# draw the thickness plot
h = plt.contourf(gridx,gridy,thickness_mean,vmin=np.min(thickness_mean),vmax=np.max(thickness_mean),cmap='coolwarm')
plt.colorbar(h)
plt.scatter(drawMSPx,drawMSPy,c='K')
np.savetxt('msps_x_{0}.txt'.format(i),drawMSPx,header='#o/k')
np.savetxt('msps_y_{0}.txt'.format(i),drawMSPy,header='#o/k')
markers = ('v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd')
colors = ['k','w']
mcount = 0
ccount = 0
print('****')
print (draw_elementsx)
for i in range(num_extra_elements):
    plt.scatter(draw_elementsx[i],draw_elementsy[i],marker=markers[mcount],c=colors[ccount],edgecolors='k')
    np.savetxt('elems_x_{0}.txt'.format(i),draw_elementsx[i],header='{0}/{1}'.format(markers[mcount],colors[ccount]))
    np.savetxt('elems_y_{0}.txt'.format(i),draw_elementsy[i],header='{0}/{1}'.format(markers[mcount],colors[ccount]))

    mcount+=1
    ccount +=1
    if mcount > len(markers)-1:
        mcount=0
    if ccount >1:
        ccount=0
    
plt.savefig('mean_thickness.png')
plt.show()
plt.close()

# draw the STD plot
h = plt.contourf(gridx,gridy,thickness_std,vmin=np.min(thickness_std),vmax=np.max(thickness_std),cmap='YlOrRd')
#### fake out line for color normalization
#h = plt.contourf(gridx,gridy,thickness_std,vmin=2.0,vmax=8.0,cmap='YlOrRd')
plt.colorbar(h)
plt.scatter(drawMSPx,drawMSPy,c='K')
markers = ('v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd')
colors = ['k','w']
mcount = 0
ccount = 0
for i in range(num_extra_elements):
    plt.scatter(draw_elementsx[i],draw_elementsy[i],marker=markers[mcount],c=colors[ccount],edgecolors='k')
    mcount+=1
    ccount +=1
    if mcount > len(markers)-1:
        mcount=0
    if ccount >1:
        ccount=0
plt.scatter(drawMSPx,drawMSPy,c='K')
plt.savefig('std_thickness.png')
plt.show()
plt.close()

