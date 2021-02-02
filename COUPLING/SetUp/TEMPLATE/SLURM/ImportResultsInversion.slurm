#!/bin/bash
#SBATCH -J ImportResultsInversion
#SBATCH --nodes=1
#SBATCH --constraint=HSW24
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
#SBATCH -e ImportResultsInversion.e%j
#SBATCH -o ImportResultsInversion.o%j
#==================================================================
# Used to import the result of an inversion through Elmer/Ice.
#
# Arguments:
#   1st- simulations family (e.g. TEST_COUPLING in PARAMETERS directory)
#   2nd- EXCHANGE_DIR
#==================================================================

#modules
module purge
module load intel
module load python
#end modules

export OMP_NUM_THREADS=1

ulimit -s unlimited

echo "IMPORTING RESULTS FROM ELMERICE INVERSION"

time
python ImportResultsInversion.py $1
time

# Launching following job
echo "RUN0.slurm $1 $2"
sbatch RUN0.slurm $1 $2



