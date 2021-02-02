# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from matplotlib import path
import numpy as np
from scipy.interpolate import griddata
import sys
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import numpy.ma as ma
from netCDF4 import Dataset
from optparse import OptionParser
from matplotlib.colors import LogNorm
import glob

# USER
#contour
contour = 'contour/AmundsenBasin_FG_NODES_copy.dat'
name = 'ASE'

#resolution of the regular grid
res = 1000.

# END USER

#get the contour
xy = np.loadtxt(contour,delimiter=' ')
xy = xy[:,0:2]
lx = xy[:,0]
ly = xy[:,1]
pa = path.Path(xy)

#create x & y
x = np.arange(lx.min(),lx.max(),res)
y = np.arange(ly.min(),ly.max(),res)

#create and build mask
mask = np.zeros((len(y),len(x)))

for i in range(0,len(x),1):
  print i+1, ' / ', len(x)
  for j in range(0,len(y),1):
    pt = [[x[i],y[j]]]
    if (pa.contains_points(pt)==True):
      mask[j,i] = 1

# creating ncfile
ncfile = Dataset('mask'+name+'.nc','w',format='NETCDF4')

ncfile.createDimension('nx',len(x))
ncfile.createDimension('ny',len(y))

xnc = ncfile.createVariable('x','f',('nx'))
ync = ncfile.createVariable('y','f',('ny'))
masknc = ncfile.createVariable('mask','f',('ny','nx'))

# writing variables
xnc[:] = x
ync[:] = y
masknc[:,:] = mask

ncfile.close()

