#!/bin/bash
#SBATCH -J checkBothMelts
#SBATCH --nodes=1
#SBATCH --constraint=HSW24
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
#SBATCH -e checkBothMelts.e%j
#SBATCH -o checkBothMelts.o%j
#==========================================================
# Used to check that the total melt sent by NEMO is the 
# same as the total melt received by Elmer.
#
# Arguments:
#   1st - number of job
#   2nd - ELMER_WORKDIR
#==========================================================

#modules
module purge
module load intel
module load python
#end modules

export OMP_NUM_THREADS=1

ulimit -s unlimited

time
python checkBothMelts.py RUN$1 $2
time

