# -*- coding: utf-8 -*-
"""
2017/01/12
@authors: Lionel Favier & Nacho Merino
"""

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
import glob

# Check arguments
parser = OptionParser()
options, args = parser.parse_args()
if(len(args)<3):
  print "#######################################################################"
  print "Usage: python displayIdealParaviewResults.py 1 2 3 (4) ...             "
  print "you need at least 3 arguments                                          "
  print "1 is the name of the .pvtu file (no extension)                         "
  print "ex : RUN10_????.pvtu -> RUN10_ as the argument"
  print "2 is the name of the contour file in .csv (with extension)             "
  print "3 is the variable you want to display and record in a netcdf file      "
  print "optional 4 and more are for more data to display and record            "
  print ""
  print "example given here can be ran using :                                  "
  print "python fromVtuToNetcdf.py RUN10 AmundsenBasin_Limits_FG.csv zs groundedmask bedrock ssavelocity"
  print "#######################################################################"
  sys.exit(1)

#Do you want to plot ? 0 or 1
plot = 0

cmap = plt.get_cmap('RdYlBu')

filein = 'real_vtu/'+args[0]
fileout = 'real_vtu/'+args[0]
contour = 'real_vtu/'+args[1]

#resolution of the regular grid
res = 2000.

#get the contour
xy = np.loadtxt(contour,delimiter=',',skiprows=1)
xy = xy[:,0:2]
lx = xy[:,0]
ly = xy[:,1]
pa = path.Path(xy)

#get the number of steps
files = glob.glob(filein+'*pvtu')
files.sort()

nbfiles = len(files)

for ff in range(nbfiles):

  fich = files[ff]
  print "The pvtu file processed is", fich

  #get the vtu data
  reader = vtk.vtkXMLPUnstructuredGridReader()
  reader.SetFileName(fich)
  reader.Update()
  output=reader.GetOutput()
  PointData=output.GetPointData()#données du nodes
  Coords=vtk_to_numpy(output.GetPoints().GetData())

  #pour obtenir les indices de nodes a la surface
  cellData=output.GetCellData()
  pointData=output.GetPointData()
  GeometryIDS=vtk_to_numpy(cellData.GetArray(0))

  #la condition de borde 1 correspond a la surface
  indexGEO=np.where(GeometryIDS==1)

  listPoints=set()#pour garde les index
  for i in indexGEO[0]:

    celda1=output.GetCell(i)
    ids=celda1.GetPointIds()
    if ids.GetNumberOfIds()==3:#utiliser que les elementes que ont 3 ids == triangles (pas le lignes ou les nodes) 
      listPoints.add(ids.GetId(0))
      listPoints.add(ids.GetId(1))
      listPoints.add(ids.GetId(2))

  indexSurface=list(listPoints)#indices de nodes de surface

  numArrays=PointData.GetNumberOfArrays()#pour obtenir les indices des propietes (Zb,Zs,Vx,Vy....)

  #getting names of the variables
  names = []
  for ii in range(2,len(args)):
    names.append(args[ii])

  varindex = []

  for ii in range(0,len(names)):

    for jj in np.arange(numArrays):
      if PointData.GetArrayName(jj)==names[ii]:
        varindex.append(jj)
        break

  varis = []
  for ii in range(0,len(names)):
    varis.append(vtk_to_numpy(PointData.GetArray(varindex[ii])))

  results = []
  for ii in range(0,len(names)):
    if (len(np.shape(varis[ii][indexSurface])) < 2):
      results.append(varis[ii][indexSurface])
    else:
      print 'you have a vector accounted for, I assume this is the velocity vector and calculate the norm'
      tmp = np.sqrt(varis[ii][indexSurface,0]*varis[ii][indexSurface,0]+varis[ii][indexSurface,1]*varis[ii][indexSurface,1])
      results.append(tmp)

  #operations only needed at the initial step
  if (ff==0):

    #make the regular grid
    xold = Coords[indexSurface,0]
    yold = Coords[indexSurface,1]
    xg = np.arange(xold.min(),xold.max(),res)
    yg = np.arange(yold.min(),yold.max(),res)
    xx,yy = np.meshgrid(xg,yg)

    mask = xx*0.

    print "Creating the mask from contour file"
    for i in range(0,len(xg),1):
      print i+1, ' / ', len(xg)
      for j in range(0,len(yg),1):
        pt = [[xg[i],yg[j]]]
        if (pa.contains_points(pt)==True):
          mask[j,i] = 1

    #creating netcdffile
    resnc = []
    #saving in netcdf file
    test = Dataset(fileout+'.nc','w',format='NETCDF4')
    xnc = test.createDimension('x',len(xg))
    ync = test.createDimension('y',len(yg))
    tnc = test.createDimension('t',nbfiles)
    #create the variables
    xncs = test.createVariable('x','f4',('x'))
    yncs = test.createVariable('y','f4',('y'))
    tncs = test.createVariable('t','i4',('t'))
    for ii in range(len(names)):
      resnc.append(test.createVariable(names[ii],'f4',('t','y','x')))
    #create the attributes
    xncs.units = 'm'
    yncs.units = 'm'
    #writing data
    xncs[:] = xg
    yncs[:] = yg


  #interpolation on regular grid
  resinterp = []

  print "Step", ff
  print "Interpolation in regular grid"
  for ii in range(len(names)):
    print "Doing ", names[ii]
    resinterp.append(griddata((Coords[indexSurface,0],Coords[indexSurface,1]),results[ii][:],(xx,yy), method='linear'))
    resinterp[ii] = np.ma.array(resinterp[ii],mask=(mask==0))

  #ploting step
  if (plot==1):

    for ii in range(len(names)):
      plt.figure(ii)
      tmp = resinterp[ii][np.isnan(resinterp[ii])==0]
      #plt1 = plt.pcolormesh(xx,yy,resinterp[ii],cmap=cmap,vmin=resinterp[ii].min(),vmax=resinterp[ii].max())
      plt1 = plt.pcolormesh(xx,yy,resinterp[ii],cmap=cmap,vmin=tmp.min(),vmax=tmp.max())
      plt.colorbar(plt1)
      plt.title(names[ii])

    plt.tight_layout()
    plt.show()

  for ii in range(len(names)):
    resnc[ii][ff,:,:] = resinterp[ii][:,:]

#closing file
test.close()

