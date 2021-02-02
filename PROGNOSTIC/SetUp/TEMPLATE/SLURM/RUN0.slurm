#!/bin/bash
#SBATCH -J RUN0

#SBATCH --nodes=1
#SBATCH --constraint=HSW24

#SBATCH --ntasks=<nbpart>
#SBATCH --ntasks-per-node=<nbpart>
##SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=1

#SBATCH --time=00:05:00

#SBATCH --output RUN0.%j.output
#SBATCH --error  RUN0.%j.error

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
module load nco/4.7.9-gcc-4.8.5-hdf5-1.8.18-openmpi-2.0.4
#end modules

export OMP_NUM_THREADS=1
#export MKL_NUM_THREADS=1
#export KMP_AFFINITY=granularity=fine,compact,1,0,verbose

ulimit -s unlimited

echo RUN0.sif > ELMERSOLVER_STARTINFO

time
srun --mpi=pmi2 -K1 --resv-ports -n $SLURM_NTASKS ElmerSolver_mpi
time


