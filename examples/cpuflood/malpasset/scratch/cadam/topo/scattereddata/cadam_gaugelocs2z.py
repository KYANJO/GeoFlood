#!/usr/bin/python

"""
cadam_gaugelocs2z.py uses the original topo file for the cadam project. It finds the z value
of the 29 gauge locations.
"""
from numpy import *
from scipy import *  
import matplotlib.pyplot as pyplot
import matplotlib.mlab as mlab
from geotools import *
from datatools import * 
import string
import pdb

topofile='malpasset.geo.txt'

xlow=-4000.0
xhi= 18000.0
ylow=-3500.0
yhi=7000.0
npts=500

di=iotools.datafile2array
XYZ=di(topofile,skiplines=2)

xyzb=[[1.,-5000 ,-4000, 0.], \
      [2.,-5000 ,8000, 0.], \
      [3.,20000 ,-4000, 0.], \
      [3.,20000 ,8000, 0.]]

xyzb=array(xyzb)
XYZ=vstack((XYZ,xyzb))


#coordinates of points
pts=[ \
(4913.1,4244.0), \
(5159.7,4369.6), \
(5790.6,4177.7), \
(5886.5,4503.9),\
(6763.0,3429.6),\
(6929.9,3591.8),\
(7326.0,2948.7),\
(7451.0,3232.1),\
(8735.9,3264.6),\
(8628.6,3604.6),\
(9761.1,3480.3),\
(9832.9,2414.7),\
(10957.2,2651.9),\
(11115.7,3800.7),\
(11689.0,2592.3),\
(11626.0,3406.8),\
(12333.7,2269.7),\
(5550.,4400.),\
(11900.,3250.),\
(13000.,2700.),\
(4947.4,4289.7),\
(5717.30,4407.6),\
(6775.1,3869.2),\
(7128.2,3162.00),\
(8585.3,3443.1),\
(9675.,3085.9),\
(10939.1,3044.8),\
(11724.4,2810.4),\
(12723.7,2485.1)]

z=[]
if True:
    for i in range(len(pts)):

        x=pts[i][0]
        y=pts[i][1]

        npts=100
        delta=5.
        xgrid=linspace(x-delta,x+delta,npts)
        ygrid=linspace(y-delta,y+delta,npts)
        X,Y=meshgrid(xgrid,ygrid)

        Z = mlab.griddata(XYZ[:,1],XYZ[:,2],XYZ[:,3],X,Y)
        Y=flipud(Y)
        Z=flipud(Z)

        if ((x<X[0,0])|(x>X[0,-1])|(y<Y[-1,0])|(y>Y[0,0])):
            print('WARNING: point %i is outside data file: (x,y)= (%g,%g)' % (i+1,x,y))
            print('**file corners are (X0,Y0)= (%g,%g), and (X1,Y1) = (%g,%g)' % (X[0,0],Y[0,0],X[-1,-1],Y[-1,-1]))
            z.append(nan)
        else:

            #find indices of four corners
            #some corners might be the same, if x or y happen to intersect X or Y
            i0 = where(X[0,:]<=x)[0][-1]
            i1 = where(X[0,:]>=x)[0][0]

            j0 = where(Y[:,0]<=y)[0][0]
            j1 = where(Y[:,0]>=y)[0][-1]

            #find height of four corners
            Z00=Z[j0,i0]
            Z01=Z[j0,i1]
            Z10=Z[j1,i0]
            Z11=Z[j1,i1]

            X00=X[j0,i0]
            X01=X[j0,i1]
            X10=X[j1,i0]
            X11=X[j1,i1]

            Y00=Y[j0,i0]
            Y01=Y[j0,i1]
            Y10=Y[j1,i0]
            Y11=Y[j1,i1]

            #find slopes of opposing lines. 
            if i0==i1:
                dzdx0=0.0
                dzdx1=0.0
            else:
                dzdx0 = (Z01-Z00)/(X11-X00)
                dzdx1 = (Z11-Z10)/(X11-X00)

            #find height of points on lines
            zy0 = Z00 + (x-X00)*dzdx0
            zy1 = Z10 + (x-X10)*dzdx1

            if j0==j1:
                dzdy=0.0
            else:
                dzdy = (zy1-zy0)/(Y11-Y00)

            z.append( zy0 + (y-Y00)*dzdy)

    z=array(z)
    z=reshape(z,(len(z),1))
    print z
#    pdb.set_trace()

dout=iotools.array2datafile

dout(z,'zvaluesofgauges.txt')



