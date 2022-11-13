#!/usr/bin/python

"""
cadam_gaugelocs2z.py uses the original topo file for the cadam project. It finds the z value
of the 29 gauge locations.
"""
from numpy import *
from scipy import *  
import matplotlib.pyplot as pyplot
from geotools import *
from datatools import * 
import string

topofile='../goutal/malpasset.geo.txt'

xgrid=linspace(xlow,xhi,npts)
ygrid=linspace(ylow,yhi,npts)
X,Y=meshgrid(xgrid,ygrid)

di=iotools.datafile2array
XYZ=di(topofile,skiplines=2)

Z = mlab.griddata(XYZ[:,0],XYZ[:,1],XYZ[:,2],X,Y)
Y=flipud(Y)
Z=flipud(Z)

XYZ=hstack((X,Y,Z))
dout=iotools.array2datafile
dout(XYZ,'malpasset.geo_gridded.txt')

