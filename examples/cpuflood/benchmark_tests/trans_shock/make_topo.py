
"""
Module to create topo and qinit data files for this example.
"""

from __future__ import absolute_import
from clawpack.geoclaw.topotools import Topography
from numpy import *

def maketopo():
    """
    Output topography file for the entire domain
    """
    nxpoints = 550
    nypoints = 1050
    xlower = -40
    xupper = 1100
    yupper = 2100
    ylower = -40
    outfile= "bathy2.topotype2"     

    topography = Topography(topo_func=topo)
    topography.x = linspace(xlower,xupper,nxpoints)
    topography.y = linspace(ylower,yupper,nypoints)
    topography.write(outfile, topo_type=2, Z_format="%22.15e")

def makeqinit():
    """
    Create qinit data file
    """
    nxpoints = 70
    nypoints = 10
    xlower = 0
    xupper = 700
    yupper = 100
    ylower = 0
    outfile= "test2.xyz"     

    topography = Topography(topo_func=qinit)
    topography.x = linspace(xlower,xupper,nxpoints)
    topography.y = linspace(ylower,yupper,nypoints)
    topography.write(outfile, topo_type=1)

def topo(x,y):
    """
    rectangular domain
    """
    # value of z at origin:  Try zmin = 80 for shoreline or 250 for no shore
    zmin = 0.0
    z = (x**2 + y**2) - zmin
    return z*0.0


def qinit(x,y):
    """
    Constant water level = 9.7m
    """
    from numpy import where
    # mx = 70
    # my = 10
    # xlower = 0
    # ylower = 0
    # dx = 700./mx
    # dy = 100./my
    # for i in range(mx):
    #     for j in range(my):
    #         y = ylower + (j-0.5)*dy
    #         x = xlower + (i-0.5)*dx
    ze = (x**2 + y**2)
    z = where(ze<0, 9.7, 9.7)
            # z = 9.7
    return z

if __name__=='__main__':
    maketopo()
    makeqinit()