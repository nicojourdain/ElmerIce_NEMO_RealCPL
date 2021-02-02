python makeMask.py # mask the mask from a contour file given in the python script

# To put all the variables from a vtu file to a netcdf file:
python FromVtuToNc_allVariables.py   # so far, no 3rd dimension is accounted for

# Same for one variable (but doing 3rd dimension):
python fromVtuToNetcdf.py


NB: ideally, merge the 2 scripts to do all variables and 3rd dimension.


# Previous MISOMIP scripts:
writeNetcdf.py
writeNetcdf2dfields.py
saveAllPlotMISOMIPIceData.sh
plotMISOMIPIceData.py
plot_example.py
