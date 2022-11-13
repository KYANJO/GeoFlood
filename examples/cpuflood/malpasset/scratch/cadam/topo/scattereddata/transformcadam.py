#! /usr/bin/python
"""
      The data given by the CADAM description (Goutal) is rotated and translated from the
	  topo data that I have tediously compiled. transformcadam.py transforms the CADAM (X,Y) data to the
	  set of axes I have been using.

          The transformation is as follows:
	  if the CADAM coordinates are denoted (X,Y,Z) and my coordinates are denoted (x,y,z) the 
	  transformation is as follows:

        |  0  1  0  |  |X| + Tx  = x
	    | -1  0  0  |  |Y| + Ty  = y
	    |  0  0  1  |  |Z| + Tz  = z

           or x = R*X + T, where

            Tx= 953595
            Ty= 1849222
	        Tz=0.

"""

from numpy import * 
import iotools

#======================================================================
def transformfile (infile,outfile):
    """
    transform data in inputfilename, with cadam coords, X,Y,Z
    to new coords (x,y,z) and output into outfile.
    """
    XYZ=iotools.datafile2array(infile)

    for row in xrange(len(XYZ)) :
        (XYZ[row,0],XYZ[row,1])=transformcoords(XYZ[row,0],XYZ[row,1])

    iotools.array2datafile(XYZ,outfile)


#========================================================================

#===========================================================================
def transformcoords (X,Y):
    """
    transform cadam coords (X,Y) to other coords (x,y).
    """
    Tx= 953595.
    Ty= 1849222.

    R=array([[0,1],[-1,0]])
    V=array([[X],[Y]])

    v=dot(R,V)+array([[Tx],[Ty]])

    x=v[0]
    y=v[1]

    return x,y
#===================================================================================



    

