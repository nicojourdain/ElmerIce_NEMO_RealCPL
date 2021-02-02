#!/bin/bash
#SBATCH -J ELMER_CPL_RUN<ID>_<SBATCH_ID>
#SBATCH --nodes=<nbnode>
#SBATCH --constraint=HSW24
#SBATCH --ntasks=<nbpart>
#SBATCH --ntasks-per-node=24
##SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:09:00
#SBATCH -e run_elmerice_<ID>_<SBATCH_ID>.e%j
#SBATCH -o run_elmerice_<ID>_<SBATCH_ID>.o%j
#==============================================================
# Used to run Elmer.
#
# Arguments:
#   1st- ID
#   2nd- EXT_NEMO
#==============================================================

#multi-threads
export OMP_NUM_THREADS=1
#export MKL_NUM_THREADS=4
#export KMP_AFFINITY=granularity=fine,compact,1,0,verbose

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

ulimit -s unlimited

EXCHANGE_DIR=<EXCHANGE_DIR>
ELMER_WORKDIR=<ELMER_WORKDIR>
NEMO_WORKDIRpro=<NEMO_WORKDIRpro>

echo RUN$1.sif > ELMERSOLVER_STARTINFO

time
srun --mpi=pmi2 -K1 --resv-ports -n $SLURM_NTASKS ElmerSolver_mpi
time

if [ '<SBATCH_ID>' == 'rel' ]; then

  next=$(($1+1))

  if [ -f RUN${next}.sif ] && [ -f RUN${next}.slurm ]; then
    echo "sbatch RUN${next}.slurm ${next}"
    sbatch RUN${next}.slurm ${next}
  else
    echo "GOING TO $EXCHANGE_DIR DIRECTORY"
    cd $EXCHANGE_DIR
    echo "sbatch exportVTUtoNEMOdraft.slurm $1 $ELMER_WORKDIR <SBATCH_ID>"
    sbatch exportVTUtoNEMOdraft.slurm $1 $ELMER_WORKDIR <SBATCH_ID>
    echo "sbatch checkBothMelts.slurm $1 $ELMER_WORKDIR"
    sbatch checkBothMelts.slurm $1 $ELMER_WORKDIR
  fi

elif [ '<SBATCH_ID>' == 'pro' ]; then

  echo "GOING TO $EXCHANGE_DIR DIRECTORY"
  cd $EXCHANGE_DIR
  echo "sbatch --parsable exportVTUtoNEMOdraft.slurm $1 $ELMER_WORKDIR <SBATCH_ID>"
  jobid=$(sbatch --parsable exportVTUtoNEMOdraft.slurm $1 $ELMER_WORKDIR <SBATCH_ID>)
  echo "sbatch checkBothMelts.slurm $1 $ELMER_WORKDIR"
  sbatch checkBothMelts.slurm $1 $ELMER_WORKDIR

  echo "GOING TO $NEMO_WORKDIRpro DIRECTORY"
  cd $NEMO_WORKDIRpro
  #this run_nemo.sh needs to wait for "exportVTUtoNEMOdraft"
  echo "sbatch --parsable --dependency=afterany:$jobid run_nemo.sh $2"
  newjobid=$(sbatch --parsable --dependency=afterany:$jobid run_nemo.sh $2)

fi

