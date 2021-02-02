import numpy as np
from scipy import interpolate  #this module needs to be before vtk
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from netCDF4 import Dataset
import pyproj
import glob
from optparse import OptionParser
from os import path as os_path
import subprocess

# Check arguments
# 1st-iith run Elmer
# 2nd-rel or pro
# 3rd-ElmerWorkDir

parser = OptionParser()
options, args = parser.parse_args()

#
rhoi = <RHOI_SI>
rhow = <RHOW_SI>

#
PATH = os_path.abspath(os_path.split(__file__)[0])

#
file1='RUN'+str(args[0])
elmerdir=str(args[2])

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

numArrays = PointData.GetNumberOfArrays()

#get the variables from ElmerIce results
for i in np.arange(numArrays):
  if PointData.GetArrayName(i)=='zb':
    VarIndex=i
    break
zb=vtk_to_numpy(PointData.GetArray(VarIndex))

for i in np.arange(numArrays):
  if PointData.GetArrayName(i)=='groundedmask':
    VarIndex=i
    break
gm=vtk_to_numpy(PointData.GetArray(VarIndex))

for i in np.arange(numArrays):
  if PointData.GetArrayName(i)=='bedrock':
    VarIndex=i
    break
bedrock=vtk_to_numpy(PointData.GetArray(VarIndex))

#change the sign because this is done the ocean manner
zb=-zb
zb[zb<0] = 0

bedrock=-bedrock
bedrock[bedrock<0] = 0

#get the coordinates, and do a new grid for displaying
xElmer = Coords[indexSurface,0]
yElmer = Coords[indexSurface,1]

#
nameMask = PATH+"/<ELMER_MASK_IN_NEMO>"

maskElmer = np.loadtxt(nameMask,delimiter=' ')
xNemo = maskElmer[:,0]
yNemo = maskElmer[:,1]
maskElmer = maskElmer[:,2]

# Get the init draft from Nemo
# this file should be taken from Nemo, then saved after every coupled period, then updated every time
bathyMeterFich = "<BATHY_METER_FICH>"

data = Dataset(bathyMeterFich,'r')
lat = data.variables['nav_lat'][:,:]
lon = data.variables['nav_lon'][:,:]
initDraft = data.variables['isf_draft'][:,:]
initbathy = data.variables['Bathymetry'][:,:]
initbathyisf = data.variables['Bathymetry_isf'][:,:]
data.close()

initDraftf = initDraft.flatten()
initbathyf = initbathy.flatten()
initbathyisff = initbathyisf.flatten()

## interpolate Elmer new draft onto Nemo grid
# zb
zbOnNemo = interpolate.griddata((xElmer,yElmer),zb,(xNemo,yNemo),method='linear')

newDraftf = np.copy(initDraftf)
newDraftf[maskElmer==1] = zbOnNemo[maskElmer==1]

#for ii in range(len(newDraftf)):
#  print xNemo[ii], yNemo[ii] ,newDraftf[ii]

sn = np.shape(initDraft)
newDraft = np.reshape(newDraftf,(sn[0],sn[1]))

#correct nans ????
#x, y = np.indices(newDraft.shape)
#
#interp = np.array(newDraft)
#interp[np.isnan(interp)] = interpolate.griddata((x[~np.isnan(newDraft)], y[~np.isnan(newDraft)]),newDraft[~np.isnan(newDraft)],(x[np.isnan(newDraft)], y[np.isnan(newDraft)]),method='linear')
#newDraft = np.array(interp)

#gm
gmonnemo = interpolate.griddata((xElmer,yElmer),gm,(xNemo,yNemo),method='nearest')
gmonnemo[maskElmer==0] = 3
gmonnemo[maskElmer==2] = 2

newgmf = np.array(gmonnemo)
newgm = np.reshape(newgmf,(sn[0],sn[1]))

# bathy
bathyOnNemo = interpolate.griddata((xElmer,yElmer),bedrock,(xNemo,yNemo),method='linear')
bathyOnNemo[maskElmer==0] = initbathyf[maskElmer==0]
bathyOnNemo[maskElmer==2] = 10

newBathyf = np.array(bathyOnNemo)
#newBathy = np.reshape(newBathyf,(sn[0],sn[1]))
newBathy = np.array(initbathy)

# bathy isf
#bathyOnNemo = interpolate.griddata((xElmer,yElmer),zb,(xNemo,yNemo),method='linear')
newBathyisff = np.copy(initbathyisff)
newBathyisff[maskElmer==1] = bathyOnNemo[maskElmer==1]
newBathyisf = np.reshape(newBathyisff,(sn[0],sn[1]))

# create the netcdf file for Nemo
if (str(args[1]) == 'rel'):
  nameFich = "isf_draft_meter_"+str(args[1])+".nc"
elif (str(args[1]) == 'pro'):
  nameFich = "isf_draft_meter_"+str(args[1])+str(args[0])+".nc"

data = Dataset(PATH+"/ISF_DRAFT/"+nameFich,'w',format='NETCDF4')
#dimensions
ync = data.createDimension('y',sn[0])
xnc = data.createDimension('x',sn[1])
#create variables
latnc = data.createVariable('nav_lat','f8',('y','x'))
lonnc = data.createVariable('nav_lon','f8',('y','x'))
isfdraftnc = data.createVariable('isf_draft','f8',('y','x'))
bathymeterisfnc = data.createVariable('Bathymetry_isf','f8',('y','x'))
bathymeternc = data.createVariable('Bathymetry','f8',('y','x'))
gmnc = data.createVariable('groundedmask','f8',('y','x'))
#writing data
latnc[:,:] = lat[:,:]
lonnc[:,:] = lon[:,:]
isfdraftnc[:,:] = newDraft[:,:]
bathymeterisfnc[:,:] = newBathyisf[:,:]
bathymeternc[:,:] = newBathy[:,:]
gmnc[:,:] = newgm[:,:]
#closing file
data.close()

subprocess.call(["ln","-sf",PATH+"/ISF_DRAFT/"+nameFich,PATH+"/ISF_DRAFT/"+"isf_draft_meter_current.nc"])

