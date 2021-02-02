#!/bin/bash
#SBATCH -J ImportResultsInversion

#SBATCH --nodes=1
#SBATCH --constraint=HSW24

#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

#SBATCH --time=00:05:00

#SBATCH --output ImportResultsInversion.%j.output
#SBATCH --error  ImportResultsInversion.%j.error

#modules
module purge
module load intel
module load python
#end modules

export OMP_NUM_THREADS=1

ulimit -s unlimited

#arguments received
#1 simulations family (e.g. TEST_COUPLING in PARAMETERS directory)
#2 EXCHANGE_DIR

echo "IMPORTING RESULTS FROM ELMERICE INVERSION"

time
python ImportResultsInversion.py $1
time

# Launching following job
sbatch RUN0.slurm



