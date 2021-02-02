#!/bin/bash
#SBATCH -J exportVTUtoNEMOdraft
#SBATCH --nodes=1
#SBATCH --constraint=HSW24
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
#SBATCH -e exportVTUtoNEMOdraft.e%j
#SBATCH -o exportVTUtoNEMOdraft.o%j
#=========================================================
# Used to put Elmer's ice draft (vtu format)
# onto NEMO's input format (netcdf).
#
# Arguments:
#   1st- family
#   2nd- ELMER_WORKDIR
#   3rd- 'rel1', 'rel2' or 'pro'
#=========================================================

#modules
module purge
module load intel
module load python
#end modules

export OMP_NUM_THREADS=1

ulimit -s unlimited

time
python exportVTUtoNEMOdraft.py $1 $3 $2
time

