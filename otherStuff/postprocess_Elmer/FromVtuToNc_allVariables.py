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
import os

#USER
#which mask (to be done beforehand)
domain = 'ASE'
simu = 'RUN1'
pathsimu = 'result/'+simu

cmap = plt.get_cmap('RdBu_r')

rhoi = 917.
g = 9.81
rgamma = 100

# END USER

# if netcdf already created, juste postprocess the rest
if (not os.path.isfile(simu+'.nc')):
  #collecting results

  #get the vtu data
  filespvtu = glob.glob(pathsimu+'*pvtu')
  filespvtu.sort()

  print filespvtu

  filein = filespvtu[-1]

  #get the vtu data
  reader = vtk.vtkXMLPUnstructuredGridReader()
  reader.SetFileName(filein)
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

  #minx miny maxx maxy
  xold = Coords[indexSurface,0]
  yold = Coords[indexSurface,1]

  #collecting mask
  maskname = 'mask'+domain+'.nc'

  if (not os.path.isfile(maskname)):
    xg = np.arange(xold.min(),xold.max(),2000)
    yg = np.arange(yold.min(),yold.max(),2000)
  else:
    ncid = Dataset(maskname,'r')
    xg = ncid.variables['x'][:]
    yg = ncid.variables['y'][:]
    mask = ncid.variables['mask'][:,:]
    ncid.close()

  xx,yy = np.meshgrid(xg,yg)

  if (not os.path.isfile(maskname)):
    mask = np.ones((len(yg),len(xg)))

  varis = []
  names = []

  for ii in range(numArrays):
    names.append(PointData.GetArrayName(ii))
    varis.append(vtk_to_numpy(PointData.GetArray(ii)))

  cpt = 0

  #creating and opening nc file
  test = Dataset(simu+'.nc','w',format='NETCDF4')
  xnc = test.createDimension('x',len(xg))
  ync = test.createDimension('y',len(yg))
  #create the variables
  xncs = test.createVariable('x','f8',('x'))
  yncs = test.createVariable('y','f8',('y'))
  #create the attributes
  xncs.units = 'm'
  yncs.units = 'm'
  #writing data
  xncs[:] = xg
  yncs[:] = yg

  #interpolations
  for ii in range(numArrays):
    #interpolating in regular grid

    print 'Interpolating ', names[cpt]

    if (len(np.shape(varis[ii][indexSurface])) < 2):
      #calculating
      resinterp = griddata((xold,yold),varis[ii][indexSurface],(xx,yy),method='linear')
      resinterp = np.ma.array(resinterp,mask=(mask==0))
      #vmin vmax ...
      vmin = resinterp[mask==1].min()
      vmax = resinterp[mask==1].max()
      #saving to netcdf
      resnc = test.createVariable(names[cpt],'f8',('y','x'))    
      #writing data
      resnc[:,:] = resinterp[:,:]
      cpt = cpt+1
    else:
      #calculating
      resinterpx = griddata((xold,yold),varis[ii][indexSurface,0],(xx,yy),method='linear')
      resinterpx = np.ma.array(resinterpx,mask=(mask==0))
      resinterpy = griddata((xold,yold),varis[ii][indexSurface,1],(xx,yy),method='linear')
      resinterpy = np.ma.array(resinterpy,mask=(mask==0))
      #vmin vmax ...
      vminx = resinterpx[mask==1].min()
      vmaxx = resinterpx[mask==1].max()
      vminy = resinterpy[mask==1].min()
      vmaxy = resinterpy[mask==1].max()
      #saving to netcdf
      resncx = test.createVariable(names[cpt]+':x','f8',('y','x'))
      resncy = test.createVariable(names[cpt]+':y','f8',('y','x'))
      #writing data
      resncx[:,:] = resinterpx[:,:]
      resncy[:,:] = resinterpy[:,:]
      cpt = cpt+1

  #closing file
  test.close()

