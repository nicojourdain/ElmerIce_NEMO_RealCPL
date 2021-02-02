#!/bin/bash
#SBATCH -J Calcul_Cs_LoiSchoof

#SBATCH --nodes=1
#SBATCH --constraint=HSW24

#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

#SBATCH --time=00:05:00

#SBATCH --output Calcul_Cs_LoiSchoof.%j.output
#SBATCH --error  Calcul_Cs_LoiSchoof.%j.error

#modules
module purge
module load intel
module load python
#end modules

export OMP_NUM_THREADS=1

ulimit -s unlimited

echo "COMPUTING THE CS FROM NON LINEAR SCHOOF FRICTION LAW"

time
python Calcul_Cs_LoiSchoof.py
time

# Launching following job
sbatch RUNPRO0.slurm 0



