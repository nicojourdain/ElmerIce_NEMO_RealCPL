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

  final[ii-1,:] = np.array([[dataCostBeta[-1,1],dataCostEta[-1,1],dataCostRms[-1,2],lbdBeta[ii-1],lbdEta[ii-1],dataCostRms[-1,0]]])

#remove simulation that did not go to the end
final = final[final[:,5]==600]
imax = len(final)

jregBeta = final[:,0]
jregEta = final[:,1]
rms = final[:,2]

lbdJregBeta = final[:,3]
lbdJregEta = final[:,4]
lastTimestep = final[:,5]

top = final[:,0:3]

xmin = jregBeta.min()
xmax = jregBeta.max()
ymin = jregEta.min()
ymax = jregEta.max()
zmin = rms.min()
zmax = rms.max()

resh = 5.
resv = 2.

interpolationmethod = 'cubic'
p = 2
extrapolation_interval = 20

def main():
    extrapolation_spots = get_plane(0, 100, 0, 100, extrapolation_interval)
    nearest_analysis(extrapolation_spots)
    plt.show()

def nearest_analysis(extrapolation_spots):
    top_extra = extrapolation(top, extrapolation_spots, method='nearest')
    gridx_top, gridy_top, gridz_top = interpolation(top_extra)
    plot(top, gridx_top, gridy_top, gridz_top, method='rotate', title='_top_nearest')

def nearest_neighbor_interpolation(data, x, y, p=0.5):
    """
    Nearest Neighbor Weighted Interpolation
    http://paulbourke.net/miscellaneous/interpolation/
    http://en.wikipedia.org/wiki/Inverse_distance_weighting

    :param data: numpy.ndarray
        [[float, float, float], ...]
    :param p: float=0.5
        importance of distant samples
    :return: interpolated data
    """
    n = len(data)
    vals = np.zeros((n, 2), dtype=np.float64)
    distance = lambda x1, x2, y1, y2: (x2 - x1)**2 + (y2 - y1)**2
    for i in range(n):
        vals[i, 0] = data[i, 2] / (distance(data[i, 0], x, data[i, 1], y))**p
        vals[i, 1] = 1          / (distance(data[i, 0], x, data[i, 1], y))**p
    z = np.sum(vals[:, 0]) / np.sum(vals[:, 1])
    return z

def get_plane(xl, xu, yl, yu, i):
    xx = np.arange(xl, xu, i)
    yy = np.arange(yl, yu, i)
    extrapolation_spots = np.zeros((len(xx) * len(yy), 2))
    count = 0
    for i in xx:
        for j in yy:
            extrapolation_spots[count, 0] = i
            extrapolation_spots[count, 1] = j
            count += 1
    return extrapolation_spots

def extrapolation(data, extrapolation_spots, method='nearest'):

    if method == 'nearest':
        new_points = np.zeros((len(extrapolation_spots), 3))
        new_points[:, 0] = extrapolation_spots[:, 0]
        new_points[:, 1] = extrapolation_spots[:, 1]
        for i in range(len(extrapolation_spots)):
            new_points[i, 2] = nearest_neighbor_interpolation(data,
                                    extrapolation_spots[i, 0], extrapolation_spots[i, 1], p=p)
        combined = np.concatenate((data, new_points))
        return combined


def interpolation(data):
    gridx, gridy = np.mgrid[xmin:xmax:resh, ymin:ymax:resh]
    gridz = gd(data[:, :2],data[:, 2], (gridx, gridy), method=interpolationmethod)
    return gridx, gridy, gridz

def plot(data, gridx, gridy, gridz, method='rotate', title='nearest'):
    def update(i):
        ax.view_init(azim=i)
        return ax,

    fig = plt.figure(figsize=(10,7))

    angles = [30]

    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection='3d')

    ax.plot_wireframe(gridx, gridy, gridz, cmap=cm.coolwarm, alpha=0.5)
    ax.scatter(data[:, 0], data[:, 1], data[:, 2], c='red')

    ax.view_init(azim=angles[0])

    ax.set_xlabel('Jreg, alpha')
    ax.set_ylabel('Jreg, gamma')
    ax.set_zlabel('Cost')

    animation.FuncAnimation(fig, update, np.arange(360 * 5), interval=1)

    for ii in range(imin,imax+1):
        ax.text(data[ii-1,0],data[ii-1,1],data[ii-1,2],str(ii))

    plt.savefig('Lsurface.png')
    plt.show()

if __name__ == '__main__':
#    sys.exit(main())
  main()
