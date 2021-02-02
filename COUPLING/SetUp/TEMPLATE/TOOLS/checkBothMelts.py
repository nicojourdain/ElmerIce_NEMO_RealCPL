import numpy as np
from netCDF4 import Dataset
import pyproj
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import glob
from optparse import OptionParser
from os import path as os_path

#Check arguments
parser = OptionParser()
options, args = parser.parse_args()

#
yeartosec = <YEARTOSEC>
rhoi = <RHOI_SI>
rhow = <RHOW_SI>

#
PATH = os_path.abspath(os_path.split(__file__)[0])
nameMask = PATH+"/<ELMER_MASK_IN_NEMO>"

maskElmer = np.loadtxt(nameMask,delimiter=' ')
maskElmer = maskElmer[:,2]

# computing total melt rates from Nemo on Elmer domain
meltRates = np.loadtxt(PATH+"/MELT_RATES/<MELT_RATES_NEMO_XY>",delimiter=' ')
meltRatesNemo = meltRates[:,2]
areasNemo = meltRates[:,3]

condMeltNemo = (maskElmer==1)

meltNemoMasked = (meltRatesNemo[condMeltNemo]*areasNemo[condMeltNemo]).sum()

# computing total melt rates from Elmer interpolated from Nemo
file1=str(args[0])
elmerdir=str(args[1])

filespvtu = glob.glob(elmerdir+'/mesh_<NBPART>/'+file1+'*pvtu')
filespvtu.sort()
file1 = filespvtu[-1]

#read the file
reader = vtk.vtkXMLPUnstructuredGridReader()
reader.SetFileName(file1)
reader.Update()
output = reader.GetOutput()
Coords = vtk_to_numpy(output.GetPoints().GetData())
PointData = output.GetPointData()

numArrays=PointData.GetNumberOfArrays()

for i in np.arange(numArrays):
  if PointData.GetArrayName(i)=='melt':
    VarIndex=i
    break
meltElmer=vtk_to_numpy(PointData.GetArray(VarIndex))

#get the indexes at the nodes' surface
cellData=output.GetCellData()
pointData=output.GetPointData()
GeometryIDS=vtk_to_numpy(cellData.GetArray(0))

#cond 1 corresponds to surface
indexGEO=np.where(GeometryIDS==1)

listPoints=set()#pour garde les index
for i in indexGEO[0]:

  celda1=output.GetCell(i)
  ids=celda1.GetPointIds()
  if ids.GetNumberOfIds()==3:#utiliser que les elementes que ont 3 ids == triangles (pas le lignes ou les nodes) 
    listPoints.add(ids.GetId(0))
    listPoints.add(ids.GetId(1))
    listPoints.add(ids.GetId(2))

indexSurface = list(listPoints)#indices de nodes de surface

#calculate the area
listPoints=set()
MeltTotal=0
for i in indexGEO[0]:

  #La primera celda de todas
  celda1 = output.GetCell(i)
  ids = celda1.GetPointIds()
  if ids.GetNumberOfIds()==3:
    mean = 0.0
    for j in range(3):
      mean = mean + meltElmer[ids.GetId(j)]/3.0
    MeltLocal = celda1.ComputeArea() * mean

  MeltTotal = MeltTotal + MeltLocal

print 'Melt rates differences between Elmer and Nemo'
print 'meltNemoMasked (Gt/a) = ', meltNemoMasked/1.0e9
print 'MeltTotal (Gt/a) = ', MeltTotal/1.0e9

print '(meltNemoMasked - MeltTotal) / MeltTotal * 100 = ', (meltNemoMasked - MeltTotal) / MeltTotal * 100.
print 'end....'


