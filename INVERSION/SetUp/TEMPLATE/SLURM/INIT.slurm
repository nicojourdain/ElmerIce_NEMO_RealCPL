#!/bin/bash
#SBATCH -J INIT_INVERSION

#SBATCH --nodes=1
#SBATCH --constraint=HSW24

#SBATCH --ntasks=<nbpart>
#SBATCH --ntasks-per-node=<nbpart>
##SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=1

#SBATCH --time=00:01:59

#SBATCH --output INIT_INVERSION.%j.output
#SBATCH --error  INIT_INVERSION.%j.error

export OMP_NUM_THREADS=1
#export MKL_NUM_THREADS=1

ulimit -s unlimited

export ELMER_FEM_REVISION=9f6699b
export ELMER_HOME=<ELMER_HOME>/elmer_$ELMER_FEM_REVISION
source <ELMER_HOME>/env/elmervars.sh
#For other sources : Renater + Personal
export LD_LIBRARY_PATH=<ELMERICELGGElibs>:$LD_LIBRARY_PATH

time
srun --mpi=pmi2 -K1 --resv-ports -n $SLURM_NTASKS ElmerSolver_mpi
time
