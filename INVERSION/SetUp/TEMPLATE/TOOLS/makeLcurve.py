from scipy import interpolate
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d import Axes3D

base = 'run_OPTIM_Ga'
Ga = 100

imin = 1
imax = 20

#writing final with rms, cost Jreg Beta, lambda Jreg Beta, cost Jreg Eta, lambda Jreg Eta

final = np.zeros((imax-imin+1,6))

lbds = np.loadtxt('../INPUT/LREG_SansDHDt_Ga100.IN')
lbdBeta = lbds[:,0]
lbdEta = lbds[:,1]

for ii in range(imin,imax+1):

  print 'PostProcessing ', ii

  fileCostRms = base+str(Ga)+'_Rcg'+str(ii)+'_Rdhdt0/Cost_Ga'+str(Ga)+'_Rcg'+str(ii)+'_Rdhdt0.dat'
  fileCostBeta = base+str(Ga)+'_Rcg'+str(ii)+'_Rdhdt0/CostReg_Ga'+str(Ga)+'_Rcg'+str(ii)+'_Rdhdt0.dat'
  fileCostEta = base+str(Ga)+'_Rcg'+str(ii)+'_Rdhdt0/CostRegEta_Ga'+str(Ga)+'_Rcg'+str(ii)+'_Rdhdt0.dat'

  print 'Considering the following files'
  print 'fileCostRms = ', fileCostRms
  print 'fileCostBeta =', fileCostBeta 
  print 'fileCostEta = ', fileCostEta

  dataCostRms = np.loadtxt(fileCostRms)
  dataCostBeta = np.loadtxt(fileCostBeta)
  dataCostEta = np.loadtxt(fileCostEta)

  final[ii-1,:] = np.array([[dataCostRms[-1,2],dataCostBeta[-1,1],lbdBeta[ii-1],dataCostEta[-1,1],lbdEta[ii-1],dataCostRms[-1,0]]])

#remove simulation that did not go to the end
final = final[final[:,5]==600]
imax = len(final)

rms = final[:,0]
jregBeta = final[:,1]
lbdJregBeta = final[:,2]
jregEta = final[:,3]
lbdJregEta = final[:,4]
lastTimestep = final[:,5]

res = 2000.

x = np.arange(jregBeta.min(),jregBeta.max(),(jregBeta.max()-jregBeta.min())/res)
y = np.arange(jregEta.min(),jregEta.max(),(jregEta.max()-jregEta.min())/res)
xx,yy = np.meshgrid(x,y)

rmsInterp = griddata((jregBeta,jregEta),rms,(xx,yy),method='linear')
#rmsInterp = interpolate.interp2d(jregBeta,jregEta,rms,kind='linear')(x,y)

#plotting
plt.figure(1,figsize=(15,7))

ax = plt.subplot(121)
#ax = plt.subplot(111,adjustable='box',aspect=1.)
ax.set_xscale('log')
ax.set_yscale('log')
plt1 = ax.scatter(jregBeta,jregEta,100,rms)
plt.colorbar(plt1)

for ii in range(imin,imax+1):
  ax.text(final[ii-1,1],final[ii-1,3],ii)

ax = plt.subplot(122)
ax.set_xscale('log')
ax.set_yscale('log')
plt1 = ax.pcolormesh(xx,yy,rmsInterp)
cs1 = ax.contourf(xx,yy,rmsInterp,np.arange(rms.min(),rms.max(),2))
cs2 = ax.contour(xx,yy,rmsInterp,np.arange(rms.min(),rms.max(),2),colors='black',linewidths=0.5)
ax.clabel(cs2,inline=1,fontsize=10)

plt.colorbar(plt1)

plt.show()


 
  

  
