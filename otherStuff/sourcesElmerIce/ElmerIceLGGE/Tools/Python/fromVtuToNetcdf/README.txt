####

Would you like to export your Elmer's Paraview output (output.pvtu) as netcdf files ?
Here is a solution under the form of Python routines

They also plots the asked variables but works onnly in parallel
Also contains a script to generate a .csv file, which will be used for the contour, from a .shp file (under a node form, not a polygon)

Authors: Lionel Favier & Nacho Merino

####

Two options available
1- displayIdealParaviewResults.py
2- displayRealParaviewResults.py

####

Option 1

Usage: python displayIdealParaviewResults.py 1 2 3 (4) ...
you need at least 3 arguments
1 is the name of the .pvtu file (no extension)
2 is the name of the contour file in .csv (with extension)
3 is the variable you want to display and record in a netcdf file
optional 4 and more are for more data to display and record

example given here can be ran using
python displayRealParaviewResults.py exp1_0001 AmundsenBasin_Limits_FG.csv zs groundedmask bedrock ssavelocity

####

Option 2

Usage: python displayIdealParaviewResults.py 1 2 (3) ...
you need at least 2 arguments
1 is the name of the .pvtu file (no extension)
2 is the variable you want to display and record in a netcdf file
optional 3 and more are for more data to display and record

example given here can be ran using :
python displayIdealParaviewResults.py run00001 zs groundedmask bedrock ssavelocity

####
