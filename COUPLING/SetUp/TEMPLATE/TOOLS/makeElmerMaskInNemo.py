import numpy as np
from matplotlib import path
from netCDF4 import Dataset
import pyproj

#
yeartosec = <YEARTOSEC>
rhoi = <RHOI_SI>
rhow = <RHOW_SI>

# Get Nemo lat, lon, make it x, y
maskFich = "<NEMO_MESH_MASK_FILE>"

data = Dataset(maskFich,'r')
lat = data.variables['nav_lat'][:,:]
lon = data.variables['nav_lon'][:,:]
data.close()

projec = pyproj.Proj('+init=EPSG:3031')
xNemo,yNemo = projec(lon,lat)

sn = np.shape(xNemo)

xNemo = xNemo.flatten()
yNemo = yNemo.flatten()

# Get Elmer contour
elmerContourFich = "<ELMER_CONTOUR_FILE>"

contour = np.loadtxt(elmerContourFich,delimiter=' ')
cx = contour[:,0]
cy = contour[:,1]
pa = path.Path(contour)

maskElmer = np.zeros(len(xNemo))

for ii in range(0,len(xNemo),1):
  print ii, ' / ', len(xNemo)
  pt = [[xNemo[ii],yNemo[ii]]]
  if (pa.contains_points(pt)==True):
    maskElmer[ii] = 1

# to deal with overlapping parts between elmer and nemo
maskMat = np.reshape(maskElmer,(sn[0],sn[1]))
maskMatRef = np.array(maskMat)

maskMatTmp = np.array(maskMatRef)
maskMatTmp[1:,:] = maskMatRef[0:-1,:]
maskMat[abs(maskMatRef-maskMatTmp)==1] = 2

maskMatTmp = np.array(maskMatRef)
maskMatTmp[:,1:] = maskMatRef[:,0:-1]
maskMat[abs(maskMatRef-maskMatTmp)==1] = 2

maskMatTmp = np.array(maskMatRef)
maskMatTmp[0:-1,:] = maskMatRef[1:,:]
maskMat[abs(maskMatRef-maskMatTmp)==1] = 2

maskMatTmp = np.array(maskMatRef)
maskMatTmp[:,0:-1] = maskMatRef[:,1:]
maskMat[abs(maskMatRef-maskMatTmp)==1] = 2

###
maskElmer = maskMat.flatten()

maskAll = np.zeros((len(xNemo),3))
maskAll[:,0] = xNemo
maskAll[:,1] = yNemo
maskAll[:,2] = maskElmer

np.savetxt("<ELMER_MASK_IN_NEMO>",maskAll,delimiter=' ',fmt='%1.6e')



