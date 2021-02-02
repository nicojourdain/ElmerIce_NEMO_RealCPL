# -*- coding: utf-8 -*-
"""
plot grounded area against time for all 3 MISMIP+ experiments (7 branches)
plot the grounding lines at 0,  100, 200 years

Created on Tue Sep  8 09:29:09 2015

@author: s.l.cornford@bris.ac.uk
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import sys
from optparse import OptionParser

# for the arguments
#parser = OptionParser()
#options, args = parser.parse_args()
#if(len(args)<1):
#  print "Usage: python plot_example.py <param>"
#  print "param can be bg, pdc, pdcstar, lazer, picob2, picob5"
#  sys.exit(1)

dirinp='NMP_CALIBZb10/PARAMS/'
dirinpza='NMP_CALIBZb10/PARAMSZa/'
dirincpl='../CPL/'

dirout='NMP_CALIBZb10/figs_GAR_GL/'

#params
par1=['BG','PDC','PDCSTAR','LAZER1','LAZER2','LAZER3','LAZER4']
par2=['BGZa500','BGZa700','PDCZa500','PDCZa700','PDCSTARZa500','PDCSTARZa700',
      'PICO2','PICO5','PICO10','PICO2Za500','PICO5Za500','PICO10Za500']
#cpls
cplsdirs=['SIMUS_CPL_COM/','SIMUS_CPL_Utide/','SIMUS_CPL03_TKE_10m/','SIMUS_CPL02_TKE_1km/']
cplsres=['COM','UTIDE','CPL03_TKE_10m','CPL02_TKE_1km']

#exps
exps=['EXP10','EXP11','EXP12','EXP13','EXP20','EXP22']

#last steps for params and cpls
parstep = 10
cplstep = 1188

#functions definitions
def tscale(time):
    """
    scale time to sqrt(time) to emphasize earlier times
    """
    return np.sqrt(time)
def intscale(time):
    """
    inverse of tscale
    """
    return time**2

def garplot(ncfile, label, color, marker, typ):
    """
    add a plot of grounded area aggainst time to current axes
    """
    if (typ=='cpl'):
      rks = range(0,1188,120)
    else:
      rks = range(0,11,1)

    ncid = Dataset(ncfile, 'r')
    gar = ncid.variables["groundedArea"][rks]*1e-6*1e-3

    time = np.arange(len(rks))*10

    #if (typ=='cpl'):
    #  time = ncid.variables["time"][rks] / 120.
    #  print time
    #else:
    #  time = ncid.variables["time"][rks]

    #print gar

    plt.plot(tscale(time), gar, 'o-', mfc=color,
             color='black', label=label, marker=marker)
    ncid.close()
    return np.max(gar)

def glplot(ncfile, seq, colora, label, marker):
    """
    add a plot of grounding line points to current axes.
    makes use of the numpy.ma.MaskedArray when reading xGL,yGL
    """
    freq=5
    ncid = Dataset(ncfile, 'r')
    time = ncid.variables["time"][:]
    lxmax = 0.0
    lxmin = 800.0
    #for i in range(0, len(times)):
    #seq = (time == times[i])
    xGL = ncid.variables["xGL"][:, seq[0]]*1e-3
    lxmax = max(np.max(xGL), lxmax)
    lxmin = min(np.min(xGL), lxmin)
    yGL = ncid.variables["yGL"][:, seq[0]]*1e-3
    plt.plot(xGL[0:len(yGL):freq], yGL[0:len(yGL):freq], 's', ms=6, mfc=colora[0],
             mec='black', label=label,marker=marker)
    return lxmin, lxmax

for ii in range(len(par1)+len(par2)+len(cplsdirs)):

  ######### first figure ############"
  plt.figure(1,figsize=(15, 5))

  plt.subplot(111)

  if (ii < len(par1)):

    simu = par1[ii]
    print 'ii = ', ii
    print 'processing of ',dirinp, ' and ', simu

    #print initial grounding line
    xmin, xmax = glplot(dirinp+exps[0]+par1[ii]+'.nc', [0], ['black'], 'Initial GL','+')

    xmin, xmax = glplot(dirinp+exps[4]+par1[ii]+'.nc', [parstep], ['black'],exps[4]+' - '+par1[ii],'<')
    xmin, xmax = glplot(dirinp+exps[0]+par1[ii]+'.nc', [parstep], ['red'], exps[0]+' - '+par1[ii],'o')
    #plt.xlim([xmin-50.0, xmax+50.0])
    xmin, xmax = glplot(dirinp+exps[1]+par1[ii]+'.nc', [parstep], ['green'],exps[1]+' - '+par1[ii],'o')
    xmin, xmax = glplot(dirinp+exps[3]+par1[ii]+'.nc', [parstep], ['blue'],exps[3]+' - '+par1[ii],'o')

    xmin, xmax = glplot(dirinp+exps[5]+par1[ii]+'.nc', [parstep], ['orange'],exps[5]+' - '+par1[ii],'<')
    xmin, xmax = glplot(dirinp+exps[2]+par1[ii]+'.nc', [parstep], ['orange'],exps[2]+' - '+par1[ii],'o')

  elif (ii < len(par1)+len(par2)):

    rkii = ii-len(par1)

    simu = par2[rkii]
    print 'rkii = ', rkii
    print 'processing of ',dirinpza, ' and ', simu

    #print initial grounding line
    xmin, xmax = glplot(dirinpza+exps[0]+par2[rkii]+'.nc', [0], ['black'], 'Initial GL','+')

    xmin, xmax = glplot(dirinpza+exps[4]+par2[rkii]+'.nc', [parstep], ['black'],exps[4]+' - '+par2[rkii],'<')
    xmin, xmax = glplot(dirinpza+exps[0]+par2[rkii]+'.nc', [parstep], ['red'], exps[0]+' - '+par2[rkii],'o')
    #plt.xlim([xmin-50.0, xmax+50.0])
    xmin, xmax = glplot(dirinpza+exps[1]+par2[rkii]+'.nc', [parstep], ['green'],exps[1]+' - '+par2[rkii],'o')
    xmin, xmax = glplot(dirinpza+exps[3]+par2[rkii]+'.nc', [parstep], ['blue'],exps[3]+' - '+par2[rkii],'o')

    xmin, xmax = glplot(dirinpza+exps[5]+par2[rkii]+'.nc', [parstep], ['orange'],exps[5]+' - '+par2[rkii],'<')
    xmin, xmax = glplot(dirinpza+exps[2]+par2[rkii]+'.nc', [parstep], ['orange'],exps[2]+' - '+par2[rkii],'o')

  else:

    rkii = ii-len(par1)-len(par2)

    simu = cplsres[rkii]
    print 'rkii = ', rkii
    print 'processing of ',dirincpl, ' and ', simu

    #print initial grounding line
    xmin, xmax = glplot(dirincpl+cplsdirs[rkii]+cplsres[rkii]+'_'+exps[0]+'.nc', [0], ['black'], 'Initial GL','+')

    xmin, xmax = glplot(dirincpl+cplsdirs[rkii]+cplsres[rkii]+'_'+exps[4]+'.nc', [cplstep], ['black'],exps[4]+' - '+cplsres[rkii],'<')
    xmin, xmax = glplot(dirincpl+cplsdirs[rkii]+cplsres[rkii]+'_'+exps[0]+'.nc', [cplstep], ['red'], exps[0]+' - '+cplsres[rkii],'o')
    #plt.xlim([xmin-50.0, xmax+50.0])
    xmin, xmax = glplot(dirincpl+cplsdirs[rkii]+cplsres[rkii]+'_'+exps[1]+'.nc', [cplstep], ['green'],exps[1]+' - '+cplsres[rkii],'o')
    xmin, xmax = glplot(dirincpl+cplsdirs[rkii]+cplsres[rkii]+'_'+exps[3]+'.nc', [cplstep], ['blue'],exps[3]+' - '+cplsres[rkii],'o')

    xmin, xmax = glplot(dirincpl+cplsdirs[rkii]+cplsres[rkii]+'_'+exps[5]+'.nc', [cplstep], ['orange'],exps[5]+' - '+cplsres[rkii],'<')
    xmin, xmax = glplot(dirincpl+cplsdirs[rkii]+cplsres[rkii]+'_'+exps[2]+'.nc', [cplstep], ['orange'],exps[2]+' - '+cplsres[rkii],'o')

  plt.xlim([350,550])

  plt.legend(frameon=True, borderaxespad=0, fontsize='small', loc='right')
  plt.xlabel(r'$x$ (km)')
  plt.ylabel(r'$y$ (km)')
  #plt.savefig("example-gl.pdf")

  plt.title('GL_'+simu)

  plt.grid()

  plt.savefig(dirout+'plot_GL_'+simu+'.png')
  plt.close(1)

  ###### second figure ##########
  plt.figure(1,figsize=(7,5))
  plt.subplot(111)

  plt.plot(tscale([100, 100]), [0, 100], color="grey")
  #plt.plot(tscale([200, 200]), [0, 100], color="grey")
  plt.xlim(tscale([0, 100]))
  plt.ylim([32, 38])

  xtlocs = tscale([0, 10, 50, 100])
  #plt.xticks(xtlocs, intscale(xtlocs))
  plt.xticks(xtlocs, ['0', '10', '50', '100'])
  plt.xlabel('Time (yr)')
  plt.ylabel('Grounded Area (x 1000 km$^3$)')

  if (ii < len(par1)):

    maxa = garplot(dirinp+exps[4]+par1[ii]+'.nc', exps[4]+' - '+par1[ii], 'black', '<','par')
    maxa = garplot(dirinp+exps[0]+par1[ii]+'.nc', exps[0]+' - '+par1[ii], 'red', 'o','par')
    maxa = garplot(dirinp+exps[1]+par1[ii]+'.nc', exps[1]+' - '+par1[ii], 'green', 'o','par')
    maxa = garplot(dirinp+exps[3]+par1[ii]+'.nc', exps[3]+' - '+par1[ii], 'blue', 'o','par')

    maxa = garplot(dirinp+exps[5]+par1[ii]+'.nc', exps[5]+' - '+par1[ii], 'orange', '<','par')
    maxa = garplot(dirinp+exps[2]+par1[ii]+'.nc', exps[2]+' - '+par1[ii], 'orange', 'o','par')

  elif (ii < len(par1)+len(par2)):

    rkii = ii-len(par1)

    maxa = garplot(dirinpza+exps[4]+par2[rkii]+'.nc', exps[4]+' - '+par2[rkii], 'black', '<','par')
    maxa = garplot(dirinpza+exps[0]+par2[rkii]+'.nc', exps[0]+' - '+par2[rkii], 'red', 'o','par')
    maxa = garplot(dirinpza+exps[1]+par2[rkii]+'.nc', exps[1]+' - '+par2[rkii], 'green', 'o','par')
    maxa = garplot(dirinpza+exps[3]+par2[rkii]+'.nc', exps[3]+' - '+par2[rkii], 'blue', 'o','par')

    maxa = garplot(dirinpza+exps[5]+par2[rkii]+'.nc', exps[5]+' - '+par2[rkii], 'orange', '<','par')
    maxa = garplot(dirinpza+exps[2]+par2[rkii]+'.nc', exps[2]+' - '+par2[rkii], 'orange', 'o','par')

  else:

    rkii = ii-len(par1)-len(par2)

    maxa = garplot(dirincpl+cplsdirs[rkii]+cplsres[rkii]+'_'+exps[4]+'.nc', exps[4]+' - '+cplsres[rkii], 'black', '<','cpl')
    maxa = garplot(dirincpl+cplsdirs[rkii]+cplsres[rkii]+'_'+exps[0]+'.nc', exps[0]+' - '+cplsres[rkii], 'red', 'o','cpl')
    maxa = garplot(dirincpl+cplsdirs[rkii]+cplsres[rkii]+'_'+exps[1]+'.nc', exps[1]+' - '+cplsres[rkii], 'green', 'o','cpl')
    maxa = garplot(dirincpl+cplsdirs[rkii]+cplsres[rkii]+'_'+exps[3]+'.nc', exps[3]+' - '+cplsres[rkii], 'blue', 'o','cpl')

    maxa = garplot(dirincpl+cplsdirs[rkii]+cplsres[rkii]+'_'+exps[5]+'.nc', exps[5]+' - '+cplsres[rkii], 'orange', '<','cpl')
    maxa = garplot(dirincpl+cplsdirs[rkii]+cplsres[rkii]+'_'+exps[2]+'.nc', exps[2]+' - '+cplsres[rkii], 'orange', 'o','cpl')

  plt.legend(loc='lower left', ncol=1, frameon=True,
         borderaxespad=0, fontsize='small')

  plt.title('GAR_'+simu)

  plt.grid()

  plt.savefig(dirout+'plot_GAR_'+simu+'.png')
  plt.close(1)

