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
import subprocess

# Check arguments
parser = OptionParser()
options, args = parser.parse_args()
if(len(args)<1):
  print "#######################################################################"
  print "Usage: python ImportResultsInversion.py family ..."
  print "you need at least 1 arguments:"
  print "The family name of your simulations category (ex: BRONDEX)"
  print "#######################################################################"
  sys.exit(1)

# input
Ga = <Ga>
Rcg = <rkcg>
Rdhdt = <rkdhdt>

pathbase = '<ELMER_PATH_INVERSE>'
nbparts = <nbpart>

#name of the inversion result file
#to be adapted
namefile = pathbase+'/RUNS_'+args[0]+'/run_OPTIM_Ga'+str(Ga)+'_Rcg'+str(Rcg)+'_Rdhdt'+str(Rdhdt)+'/mesh_'+str(nbparts)+'/'

filespvtu = glob.glob(namefile+'*pvtu')
filespvtu.sort()
filein = filespvtu[-1]

#print filein

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

varname1 = 'mu'
varname2 = 'var'

for i in np.arange(numArrays):
  if PointData.GetArrayName(i)==varname1:
    VarIndex1=i
    break
for i in np.arange(numArrays):
  if PointData.GetArrayName(i)==varname2:
    VarIndex2=i
    break

eta0 = vtk_to_numpy(PointData.GetArray(VarIndex1))
var = vtk_to_numpy(PointData.GetArray(VarIndex2))

alpha = var[indexSurface,0]
eta = var[indexSurface,1]

x = Coords[indexSurface,0]
y = Coords[indexSurface,1]

alpha_out = 10.**alpha
eta_out = (eta*eta)*eta0*0.01

alphadata = np.zeros((len(x),3))
etadata = np.zeros((len(x),3))

alphadata[:,0] = x
alphadata[:,1] = y
alphadata[:,2] = alpha_out

etadata[:,0] = x
etadata[:,1] = y
etadata[:,2] = eta_out

etafilename = 'Eta_OPTIM_Ga'+str(Ga)+'_Rcg'+str(Rcg)+'_Rdhdt'+str(Rdhdt)
alphafilename = 'Cwl_OPTIM_Ga'+str(Ga)+'_Rcg'+str(Rcg)+'_Rdhdt'+str(Rdhdt)

np.savetxt(alphafilename+'.dat',alphadata)
np.savetxt(etafilename+'.dat',etadata)

