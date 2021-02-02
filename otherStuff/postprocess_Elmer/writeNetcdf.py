# -*- coding: utf-8 -*-
"""
2017
Initial author for the mandatory fields: Nacho Merino
2018/03/09
Second author for the optional fields: Lionel Favier
"""

#still to be done
#- check the problem with seg fault when using linear interpolation in griddata
#- recompute GroundedMask and FloatingMask with the fraction of Grounded and Floating ice, respectively

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata #this module needs to be before vtk
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import netCDF4
import sys
#import os
import glob
#from optparse import OptionParser

# for the arguments
#parser = OptionParser()
#options, args = parser.parse_args()
#if(len(args)<0):
#  print "Usage: python write2dfields.py"
#  sys.exit(1)

cat = 'NMP_CALIBZb10'
cat2 = 'PARAMS'

#parameters values
rhow=1.028
rhoi=0.918
secpyear=1. #86400.0*365.2422
mpatopa=1.e6
#c=3.160e6 #(or Beta^2)power law friction coefficient (Pa m^-1/3 s^1/3)
c=1.e4 #Pa m^-1/3 a^1/3

#is calibration made below Zb=300 or Zb=10
if cat == 'CALIBZb300':
  #basis simulations
  #params=['BGZa500','BGZa700','PDCZa500','PDCZa700','PDCSTARZa500','PDCSTARZa700']
  #cs=['c3.4','c2.52','c28.97','c15.92','c33.8','c24.6']
  #params=['PICO2Za500','PICO5Za500','PICO10Za500','PICO2','PICO5','PICO10']
  #cs=['c5.6','c4.25','c4.65','c2.55','c2.25','c2.4']
  params=['BG','PDC','PDCSTAR','LAZER1','LAZER2','LAZER3','LAZER4']
  cs=['c3.7','c32.8','c64.52','c1.613','c2.53','c2.13','c2.7']
  exps=['EXP10','EXP11','EXP12','EXP13','EXP20','EXP22']

if cat == 'CALIBZb10':
  #basis simulations

  #params=['BGZa500','BGZa700','PDCZa500','PDCZa700','PDCSTARZa500','PDCSTARZa700']
  #cs=['c1.06','c0.77','c9.67','c5.13','c9.69','c5.13']

  #params=['PICO2Za500','PICO5Za500','PICO10Za500','PICO2','PICO5','PICO10']
  #cs=['c1.05','c1.12','c1.42','c0.6','c0.625','c0.72']

  #params=['BG','PDC','PDCSTAR','LAZER1','LAZER2','LAZER3','LAZER4']
  #cs=['c2.03','c26.51','c35.39','c0.32','c0.75','c0.53','c0.63']

  exps=['EXP10','EXP11','EXP12','EXP13','EXP20','EXP22']

if cat == 'NMP_CALIBZb10':
  #basis simulations

  #params=['BGZa500','BGZa700','PDCZa500','PDCZa700','PDCSTARZa500','PDCSTARZa700']
  #cs=['c1.06','c0.77','c9.67','c5.13','c9.69','c5.13']

  #params=['PICO2Za500','PICO5Za500','PICO10Za500','PICO2','PICO5','PICO10']
  #cs=['c1.05','c1.12','c1.42','c0.6','c0.625','c0.72']

  #params=['PICO2NEWZa500','PICO5NEWZa500','PICO10NEWZa500','PICO2NEW','PICO5NEW','PICO10NEW']
  #cs=['c1.05','c1.12','c1.42','c0.6','c0.625','c0.72']

  #params=['BG','PDC','PDCSTAR','LAZER1','LAZER2','LAZER3','LAZER4']
  #cs=['c2.03','c26.51','c35.48','c0.32','c0.75','c0.53','c0.63']

  #params=['LAZER5','LAZER6']
  #cs=['c0.34','c0.65']

  params=['PICO2']
  cs=['c0.6']
  exps=['EXP12']

  #exps=['EXP10','EXP11','EXP12','EXP13','EXP20','EXP22']

for ii in range(0,len(params)):

  #if (params[ii] != 'LAZER'):
  #  continue

  for jj in range(0,len(exps)):

    namefile=cat+'/'+cat2+'/'+exps[jj]+params[ii]+'.nc'
    print 'Creation of '+namefile

    #netcdf create
    path='./'
    pathNewFile= path+namefile
    ncfile = netCDF4.Dataset(pathNewFile,'w',format='NETCDF4')

    # mandatory variables for MISMIP+
    ncfile.createDimension('nPointGL',None)
    ncfile.createDimension('nTime',12)

    ncfile.createVariable('time','f',('nTime'))
    ncfile.createVariable('iceVolume','f',('nTime'))
    ncfile.createVariable('iceVAF','f',('nTime'))
    ncfile.createVariable('groundedArea','f',('nTime'))

    ncfile.createVariable('xGL','f',('nPointGL','nTime'))
    ncfile.createVariable('yGL','f',('nPointGL','nTime'))
    ncfile.createVariable('iceThicknessGL','f',('nPointGL','nTime'))
    ncfile.createVariable('uBaseGL','f',('nPointGL','nTime'))
    ncfile.createVariable('vBaseGL','f',('nPointGL','nTime'))
    ncfile.createVariable('vSurfaceGL','f',('nPointGL','nTime'))
    ncfile.createVariable('uSurfaceGL','f',('nPointGL','nTime'))
    ncfile.createVariable('uMeanGL','f',('nPointGL','nTime'))
    ncfile.createVariable('vMeanGL','f',('nPointGL','nTime'))

    # global attributes
    ncfile.description = "All units are SI"

    ncfile.close()

    ##################
    #reopen the netcdf
    ncfile = netCDF4.Dataset(namefile,'a')
    
    #dimensions
    timevar=ncfile.variables['time']

    #GL fields
    xGL= ncfile.variables['xGL']
    yGL= ncfile.variables['yGL']
    iceThicknessGL = ncfile.variables['iceThicknessGL']
    uBaseGL = ncfile.variables['uBaseGL']
    vBaseGL = ncfile.variables['vBaseGL']
    uSurfaceGL = ncfile.variables['uSurfaceGL']
    vSurfaceGL = ncfile.variables['vSurfaceGL']
    uMeanGL = ncfile.variables['uMeanGL']
    vMeanGL = ncfile.variables['vMeanGL']

    #get the files
    #basis1='/home/favierli/Calculs/2017MISMIP+/SSAStar_Schoof/'
    base='/home/favierli/Calculs/MISOMIP_params_melt/ELMERICE/'+cat2+'/'
    #base='/media/favierli/LaCie/CALCULS/PARAMS/'

    caseFirst='Results_NMP_'+exps[jj]+params[ii]+cs[ii]
    cases=['Results_NMP_'+exps[jj]+params[ii]+cs[ii]]
    runs=['Run0','Ice1r1','Ice1r2','Ice1r3','Ice1r4','Ice1r5','Ice1r6','Ice1r7','Ice1r8','Ice1r9','Ice1r10','Ice1r11']

    indexCases=0
    #Integrales
    Voldata=[]
    VolSubdata=[]
    AreaGdata=[]
    time=[0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100,110]
    indexTime=0
    for case in cases:
      for run in runs:
        indexFile=0
        if run==runs[0]: #First data comes from previous Run
          path=base+caseFirst+'/'+run+'/'
          filesIce=glob.glob(path+'*.pvtu')
          filesIce.sort()
          for fileTest in filesIce:
            file1=fileTest
        else:
          path=base+case+'/'+run+'/'
          filesIce=glob.glob(path+'*.pvtu')
          filesIce.sort()
          for fileTest in filesIce:
            file1=fileTest

        print file1

        reader = vtk.vtkXMLPUnstructuredGridReader()
        reader.SetFileName(file1)
        reader.Update()
        output=reader.GetOutput()
        Coords=vtk_to_numpy(output.GetPoints().GetData())
        PointData=output.GetPointData()

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

        numArrays=PointData.GetNumberOfArrays()   

        #getting the variables        
        for i in np.arange(numArrays):
          if PointData.GetArrayName(i)=='h':
            VarIndex=i
            break      
        h=vtk_to_numpy(PointData.GetArray(VarIndex))

        for i in np.arange(numArrays):
          if PointData.GetArrayName(i)=='zs':
            VarIndex=i
            break      
        zs=vtk_to_numpy(PointData.GetArray(VarIndex))

        for i in np.arange(numArrays):
          if PointData.GetArrayName(i)=='zb':
            VarIndex=i
            break      
        zb=vtk_to_numpy(PointData.GetArray(VarIndex))

        for i in np.arange(numArrays):
          if PointData.GetArrayName(i)=='melt':
            VarIndex=i
            break
        melt=vtk_to_numpy(PointData.GetArray(VarIndex))
	
        for i in np.arange(numArrays):
          if PointData.GetArrayName(i)=='groundedmask':
            VarIndex=i
            break
        gm=vtk_to_numpy(PointData.GetArray(VarIndex))
	
        for i in np.arange(numArrays):
          if PointData.GetArrayName(i)=='ssavelocity':
            VarIndex=i
            break
        veltmp=vtk_to_numpy(PointData.GetArray(VarIndex))
	
        #define the grounding line
        indexGL=np.where(gm==0)

        #units change MPa m a -> Pa m s
        #meltps=melt*(secpyear)**(-1.)

        #GL fields
        xGL[:,indexTime]=Coords[indexGL,0][0]
        yGL[:,indexTime]=Coords[indexGL,1][0]
        iceThicknessGL[:,indexTime]=h[indexGL]
        uBaseGL[:,indexTime]=veltmp[indexGL,0][0]
        vBaseGL[:,indexTime]=veltmp[indexGL,1][0]
        uMeanGL[:,indexTime]=veltmp[indexGL,0][0]
        vMeanGL[:,indexTime]=veltmp[indexGL,1][0]
        uSurfaceGL[:,indexTime]=veltmp[indexGL,0][0]
        vSurfaceGL[:,indexTime]=veltmp[indexGL,1][0]

        #integrals
        cellData=output.GetCellData()
        pointData=output.GetPointData()
        GeometryIDS=vtk_to_numpy(cellData.GetArray(0))
        indexGEO=np.where(GeometryIDS==1)

        #Volume
        listPoints=set()
        VolInteg=0
        for i in indexGEO[0]:
          #La primera celda de todas
          celda1=output.GetCell(i)
          ids=celda1.GetPointIds()
          if ids.GetNumberOfIds()==3:
            mean=0
            for j in np.arange(3):
              mean=mean+h[ids.GetId(j)]/3
            VolLocal=mean*celda1.ComputeArea()
          VolInteg=VolInteg+VolLocal
	
        Voldata.append(VolInteg)
	
        #Volume_Submerged
        listPoints=set()
        VolInteg=0
        for i in indexGEO[0]:
          #La primera celda de todas
          celda1=output.GetCell(i)
          ids=celda1.GetPointIds()
          if ids.GetNumberOfIds()==3:
            mean=0
            for j in np.arange(3):
              mean=mean+zb[ids.GetId(j)]/3
          VolLocal=mean*celda1.ComputeArea()
          if VolLocal>0:
            VolLocal=0
          else:
            VolLocal=-VolLocal
          VolInteg=VolInteg+VolLocal
	
        VolSubdata.append(VolInteg)

        #Grounded Area
        listPoints=set()
        AreaInteg=0
        for i in indexGEO[0]:
          #La primera celda de todas
          celda1=output.GetCell(i)
          ids=celda1.GetPointIds()
          if ids.GetNumberOfIds()==3:
            mean=0
            for j in np.arange(3):
              mean=mean+gm[ids.GetId(j)]/3
          if mean>0:
            AreaLocal=celda1.ComputeArea()
          else:
            AreaLocal=0.
          AreaInteg=AreaInteg+AreaLocal

        AreaGdata.append(AreaInteg)
	
        #incrementing
        indexFile=indexFile+1
        indexTime=indexTime+1

    #dimensions
    print np.shape(Voldata),indexTime

    timeVar = ncfile.variables['time']

    #integrals
    iceVolume=ncfile.variables['iceVolume']
    iceVAF=ncfile.variables['iceVAF']
    groundedArea=ncfile.variables['groundedArea']

    iceVolume[0:len(Voldata)]=np.array(Voldata)
    groundedArea[0:len(Voldata)]=np.array(AreaGdata)
    iceVAF[0:len(Voldata)]=np.array(Voldata)-np.array(VolSubdata)*(rhow/rhoi)
    timeVar[0:len(time)]=np.array(time)

    ncfile.close()
