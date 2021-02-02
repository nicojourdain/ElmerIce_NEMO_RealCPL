## To add in the bashrc or bash_profile

#Using a version installed by Fab or Mondher
#from the official ElmerIce (from CSC)

#Version: 8.4 (Rev: 9f6699b, Compiled: 2019-03-07)
export ELMER_FEM_REVISION=9f6699b
export ELMER_HOME=${SHAREDSCRATCHDIR}/local/elmer_${ELMER_FEM_REVISION}
source ${SHAREDSCRATCHDIR}/local/env/elmervars.sh

## For using solvers that are in the local ElmerIce (from Renater)

-> https://groupes.renater.fr/wiki/elmerice/elmericegit

Or use the Sources that are in "ElmerIce_Real_TMP"

For coupling, so far
* SSAStar.F90
* EffectivePressure.F90
* USF_WaterPressure.F90


