# -*- coding: utf-8 -*-
import numpy as np
from scipy.interpolate import griddata
import sys
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import numpy.ma as ma
from netCDF4 import Dataset
import glob
import subprocess

# friction parameters
m = <M>
Cmax = <CMAX>  # Iken limit (taub/Neff<Cmax to be able to find the Schoof friction coefficient)
Clim = <CLIM>  # interpolation/extrapolation for all points where taub/Neff > Clim * Cmax

pathbase = '<ELMER_WORKDIRrel>'

pathbaseoutputs = '<ELMER_WORKDIRproSCHOOF>'
nameoutput = '<NAMEFILE_CSSCHOOF>'

nbparts = <nbpart>

namefile = pathbase+'/mesh_'+str(nbparts)+'/'+'RUN1'

filespvtu = glob.glob(namefile+'*pvtu')
filespvtu.sort()
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

varGm = 'groundedmask'
varCwl = 'cwl'
varEffPres = 'effective pressure'
varVel = 'ssavelocity'

for i in np.arange(numArrays):
  if PointData.GetArrayName(i)==varGm:
    VarIndexGm=i
    break
for i in np.arange(numArrays):
  if PointData.GetArrayName(i)==varCwl:
    VarIndexCwl=i
    break
for i in np.arange(numArrays):
  if PointData.GetArrayName(i)==varEffPres:
    VarIndexEffPres=i
    break
for i in np.arange(numArrays):
  if PointData.GetArrayName(i)==varVel:
    VarIndexVel=i
    break

gm = vtk_to_numpy(PointData.GetArray(VarIndexGm))
cwl = vtk_to_numpy(PointData.GetArray(VarIndexCwl))
Neff = vtk_to_numpy(PointData.GetArray(VarIndexEffPres))
vel = vtk_to_numpy(PointData.GetArray(VarIndexVel))

velx = vel[indexSurface,0]
vely = vel[indexSurface,1]

x = Coords[indexSurface,0]
y = Coords[indexSurface,1]

## Forcing Neff=0 to floating nodes
Neff[ gm==-1 ] = 0 #(TMP)

## Initial Taub with Linear Weertman friction law
taubInit = cwl*np.sqrt(velx**2+vely**2)

cond1 = ((taubInit/Neff) < (Clim*Cmax))
cond2 = (gm == -1)

x = x[ cond1 | cond2 ]
y = y[ cond1 | cond2 ]
gm = gm[ cond1 | cond2 ]
cwl = cwl[ cond1 | cond2 ]
Neff = Neff[ cond1 | cond2 ]
velx = velx[ cond1 | cond2 ]
vely = vely[ cond1 | cond2 ]
taubInit = taubInit[ cond1 | cond2 ]

vel = np.sqrt(velx**2+vely**2)

Cs = taubInit/(vel**m*(1.-(taubInit/(Cmax*Neff))**(1./m))**m)
Cs[ gm == -1 ] = 1.e-3

csdata = np.zeros((len(x),3))

csdata[:,0] = x
csdata[:,1] = y
csdata[:,2] = Cs

np.savetxt(pathbaseoutputs+'/'+nameoutput,csdata)


