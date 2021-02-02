## ElmerIce_NEMO_RealCPL

ElmerIce / NEMO coupling for  the real world
The current example given is for the Amundsen Sea Sector

*Contributors:* Nacho Merino, Lionel Favier, Nicolas Jourdain

### CONTOU
The contours of your real basins, i.e. ASCII file containing the list of (x y) values (in meter on the stereographic grid) that limit the area simulated by Elmer/Ice:
* One contour so far "AmundsenBasin\_FG\_NODES.dat"

### COUPLING
The framework for coupling simulations

### DATA\_SETS
For your data

### INVERSION
The framework for ElmerIce inversions

### PARAMETERS
A list of FILES.IN parameters files for your configuration including :
* Physical parameters (rhoi, yeartosec ...)
* The path to your ElmerIce mesh
* The path to your ElmerIce inversion results
* The path to your data

### otherStuff
Directories
* completion : Enables autocompletion based on the configuration files that are in the **PARAMETERS** directory
* notImportantStuff : Stuff that is not important :) so far a logo for the coupled model
* sourcesElmerIce : How to set up your ElmerIce environment in Occigen. The ElmerIce sources are divided into two kinds of origins :
  * The first is the official repository of ElmerIce from the CSC
  * The second is the local repository of ElmerIce (renater) -> SO FAR I COPIED PASTED THE SORCES NEEDED FROM THE REPOSITORY, BUT ACTUALLY THE Elmer/Ice-LGGE Git repository SHOULD BE GIT CLONED INSTEAD (https://groupes.renater.fr/wiki/elmerice/elmericegit)


