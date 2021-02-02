#!/bin/bash
#SBATCH -J makeElmerMaskInNemo
#SBATCH --nodes=1
#SBATCH --constraint=HSW24
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
#SBATCH -e makeElmerMaskInNemo.e%j
#SBATCH -o makeElmerMaskInNemo.o%j
#==========================================================
# Used to compute the mask used to interpolate Elmer's
# ice draft onto NEMO's grid. The rest of NEMO's domain
# keeps the original ice draft (e.g. BedMachine). 
#
# Arguments:
#   1st- simulations family
#   2nd- EXCHANGE_DIR
#==========================================================

#modules
module purge
module load intel
module load python
#end modules

export OMP_NUM_THREADS=1

ulimit -s unlimited

echo "MAKING THE ELMER MASK WITHIN NEMO DOMAIN"
time
python makeElmerMaskInNemo.py
time

echo "END OF STAGEBLOCK 1"
