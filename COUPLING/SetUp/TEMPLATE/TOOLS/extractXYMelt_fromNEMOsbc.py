import numpy as np
from netCDF4 import Dataset
import pyproj
from os import path as os_path
from optparse import OptionParser

# Check arguments
parser = OptionParser()
options, args = parser.parse_args()

#
yeartosec = <YEARTOSEC>
rhoi = <RHOI_SI>
rhow = <RHOW_SI>

#
PATH = os_path.abspath(os_path.split(__file__)[0])

#
nameMask = PATH+"/<ELMER_MASK_IN_NEMO>"

maskElmer = np.loadtxt(nameMask,delimiter=' ')
maskElmer = maskElmer[:,2]

#for lat, lon, fwfisf
namefich = PATH+"/MELT_RATES/<SBC_FILE>"

data = Dataset(namefich,'r')
lat = data.variables['nav_lat'][:,:]
lon = data.variables['nav_lon'][:,:]
melt = data.variables['fwfisf'][:,:,:]
data.close()

#for area of cells
namefich = "<NEMO_MESH_MASK_FILE>"

data = Dataset(namefich,'r')
e1t = data.variables['e1t'][:,:]
e2t = data.variables['e2t'][:,:]
data.close()

#
areasNemo = e1t*e2t
melt1 = np.mean(melt,0)*yeartosec*rhoi/rhow*0.001

projec = pyproj.Proj('+init=EPSG:3031')
xNemo,yNemo = projec(lon,lat)

xfin = xNemo.flatten()
yfin = yNemo.flatten()
meltfin = melt1.flatten()
areafin = areasNemo.flatten()

#
condMeltNemo = (maskElmer==1)
meltNemoMasked = (meltfin[condMeltNemo]*areafin[condMeltNemo]).sum()

data = np.zeros((len(xfin),4))
data[:,0] = xfin
data[:,1] = yfin
data[:,2] = meltfin
data[:,3] = areafin

datat = np.zeros((len(xfin),3))
datat[:,0] = xfin
datat[:,1] = yfin
datat[:,2] = np.ones(len(xfin))*meltNemoMasked

if ( args[0] == 'pro' ):
  np.savetxt(PATH+"/MELT_RATES/melt_rates_nemo_xy_"+str(args[0])+str(args[1]),data,delimiter=' ',fmt='%1.6e')
else:
  np.savetxt(PATH+"/MELT_RATES/melt_rates_nemo_xy_"+str(args[0]),data,delimiter=' ',fmt='%1.6e')

