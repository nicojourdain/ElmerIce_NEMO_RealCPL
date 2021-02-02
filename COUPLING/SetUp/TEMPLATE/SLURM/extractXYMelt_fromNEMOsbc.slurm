#!/bin/bash
#SBATCH -J extractXYMelt_fromNEMOsbc
#SBATCH --nodes=1
#SBATCH --constraint=HSW24
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
#SBATCH -e extractXYMelt_fromNEMOsbc.e%j
#SBATCH -o extractXYMelt_fromNEMOsbc.o%j
#============================================================
# Used to put NEMO's melt rate (netcdf) onto a format 
# readable by Elmer/Ice (ASCII).
#
# Arguments to be passed:
#    1st- rel1, rel2 or pro
#    2nd- NEMO_WORKDIRrel1, NEMO_WORKDIRrel2 or NEMO_WORKDIRpro
#    3th- step of nemo simulation, only for prognostic simulations
#============================================================

#modules
module purge
module load intel/17.0
module load openmpi/intel/2.0.1
module load hdf5/1.8.17
module load netcdf/4.4.0_fortran-4.4.2
module load python
module load nco
#end modules

export OMP_NUM_THREADS=1

ulimit -s unlimited

EXCHANGE_DIR=<EXCHANGE_DIR>
ELMER_WORKDIRpro=<ELMER_WORKDIRpro>
MELT_RATES_NEMO_XY=<MELT_RATES_NEMO_XY>
INIT_YEAR=<INIT_YEAR>
INIT_MONTH=<INIT_MONTH>
MONTHS_TO_SAVE_RELAX1=<MONTHS_TO_SAVE_RELAX1>
MONTHS_TO_SAVE_RELAX2=<MONTHS_TO_SAVE_RELAX2>

CONFIG=AMUXL025
CASE=GNJ002_BM01mv

DIR=$2/output_sbc

if [ $1 == 'rel1' ]; then
  MONTHS_TO_SAVE=$MONTHS_TO_SAVE_RELAX1
elif [ $1 == 'rel2' ]; then
  MONTHS_TO_SAVE=$MONTHS_TO_SAVE_RELAX2
elif [ $1 == 'pro' ]; then
  MONTHS_TO_SAVE=1
fi

rootfile=${CONFIG}-${CASE}_1m_fwfisf_

#which is the last fwfisf file
imax=0
for file in ${DIR}/*nc
do
  number=`echo $file | rev | cut -d"@" -f1 | rev | cut -d"." -f1`
  if [ $number -gt $imax ]; then
    imax=$number
  fi
done
imin=$(($imax-${MONTHS_TO_SAVE}+1))

if [ $1 == 'rel1' ] || [ $1 == 'rel2' ]; then

  #now make the mean of the last 36 months from $max-35
  echo "Processing the last $MONTHS_TO_SAVE months"
  cmd=ncrcat
  for ((ii=$imin ; ii<=$imax ; ii++))
  do
    cmd=$(echo $cmd)" ${DIR}/${rootfile}@$ii.nc"
  done
  cmd=$(echo $cmd)" ${EXCHANGE_DIR}/MELT_RATES/${rootfile}@steps${imin}to${imax}_$1.nc"

  $cmd

  ncra -O -F -d time_counter,1,${MONTHS_TO_SAVE} ${EXCHANGE_DIR}/MELT_RATES/${rootfile}@steps${imin}to${imax}_$1.nc ${EXCHANGE_DIR}/MELT_RATES/${rootfile}@steps${imin}to${imax}_$1_TMP.nc
  mv ${EXCHANGE_DIR}/MELT_RATES/${rootfile}@steps${imin}to${imax}_$1_TMP.nc ${EXCHANGE_DIR}/MELT_RATES/${rootfile}@steps${imin}to${imax}_$1.nc

  ln -sf ${EXCHANGE_DIR}/MELT_RATES/${rootfile}@steps${imin}to${imax}_$1.nc ${EXCHANGE_DIR}/MELT_RATES/melt_rates_current.nc

elif [ $1 == 'pro' ]; then

  cp ${DIR}/${rootfile}@${imax}.nc ${EXCHANGE_DIR}/MELT_RATES/${rootfile}@steps${imax}_$1.nc

  ln -sf ${EXCHANGE_DIR}/MELT_RATES/${rootfile}@steps${imax}_$1.nc ${EXCHANGE_DIR}/MELT_RATES/melt_rates_current.nc

fi


#check for NaNs
FILE_TEST=${EXCHANGE_DIR}/MELT_RATES/melt_rates_current.nc
NBNAN=`ncdump -v fwfisf ${FILE_TEST} |grep NaN |wc -l`
if [ $NBNAN -gt 0 ]; then
  echo '**************************************************'
  echo '      NaNs were found in the output files         '
  echo '             >>>>>>>>>   stopping here            '
  echo '**************************************************'
  rm ${EXCHANGE_DIR}/MELT_RATES/melt_rates_current.nc
  exit
fi

time
python extractXYMelt_fromNEMOsbc.py $1 $3
time

ln -sf ${EXCHANGE_DIR}/MELT_RATES/${MELT_RATES_NEMO_XY}_$1$3 ${EXCHANGE_DIR}/MELT_RATES/${MELT_RATES_NEMO_XY}

if [ $1 == 'pro' ]; then

  echo "GOING TO $ELMER_WORKDIRpro DIRECTORY"
  cd $ELMER_WORKDIRpro

  next=$(($3+1))
  echo "sbatch RUN${next}.slurm ${next} $1"
  sbatch RUN${next}.slurm ${next} $1

fi



