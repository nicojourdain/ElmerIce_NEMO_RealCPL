#!/usr/bin/python
import numpy
from optparse import OptionParser
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from netCDF4 import Dataset

'''
This script can be used to plot MISOMIP1 ice data in the standard
NetCDF format.  The script produces images for each time index of
each field, which can be used to produce movies if desired. Image 
files are skipped if already exist already.  (To regenerate images, 
the user must first delete the existing images.)
Usage: plotMISOIP1IceData.py <in_file.nc> <out_dir>
'''
def plotIce(fileName, outFolder):
  def plot(varName, label, cmap, scale=None):
    # the file name is the variable followed by the zero-padded time intex
    imageFileName = '%s/%s_%04i.png'%(outFolder, varName, timeIndex)
    if(os.path.exists(imageFileName)):
      # the image exists so we're going to save time and not replot it
      return
  
    # get the variable from the netCDF file
    var = ncFile.variables[varName]
  
    # the axes are 'xy', 'xz', or 'yz'
    axes = '%s%s'%(var.dimensions[2][1],var.dimensions[1][1])
  
    # get the extent based on the axes
    extent = extents[axes]
  
    # aspect ratio
    if(axes == 'xy'):
      # pixels are 1:1
      aspectRatio = None
    else:
      # stretch the axes to fill the plot area
      aspectRatio = 'auto'
  
    field = ncFile.variables[varName][timeIndex,:,:]
    
    # scale the variable if a scale is given
    if scale is not None:
      field *= scale
  
    # get the plotting limits from the dictionary of limits we created
    (lower, upper) = limits[varName]
  
    # make a figure
    fig = plt.figure(1, figsize=[8,3], dpi=100, facecolor='w')
  
    # activate the specified colorbar and set the background color
    cmap = plt.get_cmap(cmap)
    cmap.set_bad(backgroundColor)
  
    # clear the figure from the last plot
    plt.clf()
  
    # plot the data as an image
    plt.imshow(field, extent=extent, cmap=cmap, vmin=lower, vmax=upper, 
               aspect=aspectRatio, interpolation='nearest')
  
    #plt.colorbar()
    plt.title(label)
  
    if(axes == 'xy'):
      # y axis will be upside down in imshow, which we don't want for xy
      plt.gca().invert_yaxis()
      plt.xlabel('x (km)')
      plt.ylabel('y (km)')
    elif(axes == 'xz'):
      # upside-down y axis is okay
      plt.xlabel('x (km)')
      plt.ylabel('z (m)')
    else:
      # upside-down y axis is okay
      plt.xlabel('y (km)')
      plt.ylabel('z (m)')

    plt.xlim([400.,640.])
 
    # save the figure as an image
    plt.tight_layout()
    plt.draw()
    plt.savefig(imageFileName, dpi=100)
    plt.close()
  
  def makeFerretColormap():
    red = numpy.array([[0,0.6],
                       [0.15,1],
                       [0.35,1],
                       [0.65,0],
                       [0.8,0],
                       [1,0.75]])
    
    green = numpy.array([[0,0],
                         [0.1,0],
                         [0.35,1],
                         [1,0]])
    
    
    blue = numpy.array([[0,0],
                       [0.5,0],
                       [0.9,0.9],
                       [1,0.9]])
    
    colorCount = 21
    ferretColorList = numpy.ones((colorCount,4),float)
    ferretColorList[:,0] = numpy.interp(numpy.linspace(0,1,colorCount),red[:,0],red[:,1])
    ferretColorList[:,1] = numpy.interp(numpy.linspace(0,1,colorCount),green[:,0],green[:,1])
    ferretColorList[:,2] = numpy.interp(numpy.linspace(0,1,colorCount),blue[:,0],blue[:,1])
    ferretColorList = ferretColorList[::-1,:]
    
    cmap = plt.get_cmap('RdYlBu') #colors.LinearSegmentedColormap.from_list('ferret',ferretColorList,N=255)
    return cmap
  
  
  
  try:
    os.makedirs(outFolder)
  except OSError:
    pass
  
  # a rainbow colormap based on the default map in the application ferret
  cmap = makeFerretColormap()
  
  # get the filename without the path
  baseName = os.path.basename(fileName)
  # pull out the first part of the file name as the experiment name,
  # following the ISOMIP+ and MISOMIP1 filenaming requirements
  experiment = baseName.split('_')[0]
  
  # make sure this is a valid experiment (i.e. that the file has been
  # named correctly)
  #if(experiment not in ['IceOcean1','IceOcean2']):
  #  print "Unknown experiment", experiment
  #  exit(1)
  
  # open the netCDF file with the ISOMIP+ or MISOMIP1 ocean data
  ncFile = Dataset(fileName,'r')
  
  # convert x and y to km
  x = 1e-3*ncFile.variables['x'][:]
  y = 1e-3*ncFile.variables['y'][:]
  time = ncFile.variables['time'][:]
  print len(time)
  
  
  # the extents of the different plots for use in imshow
  extents = {}
  # the y extent is max then min because the y axis then gets flipped 
  # (imshow is weird that way)
  extents['xy'] = [numpy.amin(x),numpy.amax(x),numpy.amax(y),numpy.amin(y)]
  
  # set the limits for the colorbars for each field
  # NOTE: Users should feel free to modify these limits to suit their needs.
  limits = {}
  uLimits = [0., 1000.]
  vLimits = [-200., 200.]
  limits['iceThickness'] = [0, 1800.]
  limits['upperSurface'] = [0, 1600.]
  limits['lowerSurface'] = [-720., 400.]
  limits['basalMassBalance'] = [-50.,5.]
  limits['groundedMask'] = [0., 1.]
  limits['floatingMask'] = [0., 1.]
  limits['basalTractionMagnitude'] = [0., 80.]
  limits['uBase'] = uLimits
  limits['vBase'] = vLimits
  limits['uSurface'] = uLimits
  limits['vSurface'] = vLimits
  limits['uMean'] = uLimits
  limits['vMean'] = vLimits
  
  sPerYr = 365.*24.*60.*60.
  
  # light gray for use as an "invalid" background value wherever
  # data has been masked out in the NetCDF file
  backgroundColor = (0.9,0.9,0.9)
 
  for timeIndex in range(0,len(time)):
    yr=time[timeIndex]/sPerYr
    print timeIndex, 'year: ', yr
      
    # make plots for all the fields
    # meltRate needs to be scaled to m/a (as we are used to seeing)
    plot('iceThickness', 'ice thickness (m), time = '+str(yr), cmap) 
    plot('upperSurface', 'upper ice surface (m), time = '+str(yr), cmap) 
    plot('lowerSurface', 'lower ice surface (m), time = '+str(yr), cmap) 
    #plot('basalMassBalance', 'basal mass balance (m/a ice), time = '+str(yr), cmap, scale=sPerYr)
    plot('basalMassBalance', '', cmap, scale=sPerYr)
    plot('groundedMask', 'grounded mask, time = '+str(yr), cmap)
    plot('floatingMask', 'floating mask, time = '+str(yr), cmap)
    plot('uBase', 'basal x-velocity (m/a), time = '+str(yr), cmap, scale=sPerYr)
    plot('vBase', 'basal y-velocity (m/a), time = '+str(yr), cmap, scale=sPerYr)
    plot('uSurface', 'surface x-velocity (m/a), time = '+str(yr), cmap, scale=sPerYr)
    plot('vSurface', 'surface y-velocity (m/a), time = '+str(yr), cmap, scale=sPerYr)
    plot('uMean', 'mean x-velocity (m/a), time = '+str(yr), cmap, scale=sPerYr)
    plot('vMean', 'mean y-velocity (m/a), time = '+str(yr), cmap, scale=sPerYr)
               
  ncFile.close()

if __name__ == "__main__":
  # we could add some optional command-line argument here but so far none...
  parser = OptionParser()
  options, args = parser.parse_args()

  if(len(args) < 2):
    print "usage: plotMISOMIP1IceData.py <in_file.nc> <out_dir>"
    exit(1)

  # arguments are the COM file name and the directory for images
  fileName = args[0]
  outFolder = args[1]

  plotIce(fileName, outFolder)