#!/usr/bin/env python2

import sys
#from pykrige.ok import OrdinaryKriging
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d
from scipy.interpolate import griddata as gd


base = 'run_OPTIM_Ga'
Ga = 100

Rcg = 3

imin = 1
imax = 40

#writing final with rms, cost Jreg Beta, lambda Jreg Beta, cost Jreg Eta, lambda Jreg Eta

final = np.zeros((imax-imin+1,5))

lbd = np.loadtxt('../INPUT/LREG_AvecDHDt_Ga100.IN')

for ii in range(imin,imax+1):

  print 'PostProcessing ', ii

  fileCostRms = base+str(Ga)+'_Rcg'+str(Rcg)+'_Rdhdt'+str(ii)+'/Cost_Ga'+str(Ga)+'_Rcg'+str(Rcg)+'_Rdhdt'+str(ii)+'.dat'
  fileCostDhdt = base+str(Ga)+'_Rcg'+str(Rcg)+'_Rdhdt'+str(ii)+'/Cost_dHdt_Ga'+str(Ga)+'_Rcg'+str(Rcg)+'_Rdhdt'+str(ii)+'.dat'

  print 'Considering the following files'
  print 'fileCostRms = ', fileCostRms
  print 'fileCostDhdt =', fileCostDhdt 

  dataCostRms = np.loadtxt(fileCostRms)
  dataCostDhdt = np.loadtxt(fileCostDhdt)

  # 1-Jdiv, 2-Jv, 3-rms, 4-lambdaDhdt, 5-iteration
  final[ii-1,:] = np.array([[dataCostDhdt[-1,1],dataCostRms[-1,1],dataCostRms[-1,2],lbd[ii-1],dataCostRms[-1,0]]])

  #print '1-Jdiv, 2-Jv, 3-rms, 4-lambdaDhdt, 5-iteration'
  #print final[ii-1,:]

#remove simulation that did not go to the end
final = final[final[:,4]==600]
imax = len(final)

jdiv = final[:,0]
jv = final[:,1]
rms = final[:,2]

#plotting
plt.figure(1,figsize=(10,7))

ax = plt.subplot(111)
#ax = plt.subplot(111,adjustable='box',aspect=1.)
#ax.set_xscale('log')
#ax.set_yscale('log')
plt1 = ax.scatter(jv,jdiv,10,rms,cmap=cm.coolwarm)
plt.colorbar(plt1)

ax.set_xlabel('Jv')
ax.set_ylabel('Jdiv')

for ii in range(imin,imax+1):
  ax.text(jv[ii-1],jdiv[ii-1],str(ii))

plt.show()


