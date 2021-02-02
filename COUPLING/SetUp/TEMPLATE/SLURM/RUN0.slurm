#!/bin/bash
#SBATCH -J ELMER_CPL_RUN0
#SBATCH --nodes=<nbnode>
##SBATCH --mem=118000MB
#SBATCH --constraint=HSW24
#SBATCH --ntasks=<nbpart>
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
#SBATCH -e run_elmerice_0.e%j
#SBATCH -o run_elmerice_0.o%j
#=================================================================
# Used to run the first Elmer bit.
#
# Arguments:
#   1st- simulations family
#   2nd- EXCHANGE_DIR
#=================================================================

#modules
#Version: 8.4 (Rev: 9f6699b, Compiled: 2019-03-07)
export ELMER_FEM_REVISION=9f6699b
export ELMER_HOME=<ELMER_HOME>/elmer_$ELMER_FEM_REVISION
export NC_INC="-I `nc-config --includedir`"
export NC_LIB="`nc-config --libs` -lnetcdff"
source <ELMER_HOME>/env/elmervars.sh
#For other sources : Renater + Personal
export LD_LIBRARY_PATH=<ELMERICELGGElibs>:$LD_LIBRARY_PATH
module load python
module load nco
#end modules

export OMP_NUM_THREADS=1

ulimit -s unlimited

echo RUN0.sif > ELMERSOLVER_STARTINFO

echo "LAUNCHING ELMERICE RUN0 TO INITIALISE AFTER INVERSION"
time
srun --mpi=pmi2 -K1 --resv-ports -n $SLURM_NTASKS ElmerSolver_mpi
time

echo "GOING TO $2 DIRECTORY"
cd $2

# Launching following job
echo "sbatch makeElmerMaskInNemo.slurm $1 $2"
sbatch makeElmerMaskInNemo.slurm $1 $2


