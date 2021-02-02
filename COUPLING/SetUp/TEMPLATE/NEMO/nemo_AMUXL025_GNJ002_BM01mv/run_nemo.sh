#!/bin/bash
#SBATCH -C HSW24
#SBATCH --nodes=4
#SBATCH --ntasks=96
#SBATCH --ntasks-per-node=24
#SBATCH --threads-per-core=1
#SBATCH -J NEMO_CPL_RUN_<SBATCH_ID>
#SBATCH -e run_nemo_<SBATCH_ID>.e%j
#SBATCH -o run_nemo_<SBATCH_ID>.o%j
#SBATCH --time=00:19:00
#=================================================================================
# Used to run NEMO.
#
# Argument:
#   1st- block name: rel1, rel2, or pro
#=================================================================================

module purge
module load intel/17.0
module load openmpi/intel/2.0.1
module load hdf5/1.8.17
module load netcdf/4.4.0_fortran-4.4.2
module load python
module load nco

date

set -x
ulimit -s unlimited

#=================================================================================
#=================================================================================
# 0- User's choices
#=================================================================================
#=================================================================================

CONFIG='AMUXL025'  ## FULL CONFIG NAME (e.g. "trop075" or "trop075_nest025")
                 ## NB: THIS NAME SHOULD NOT START WITH A NUMBER

CONFPAR=$CONFIG #- IF NO NEST SHOULD BE EQUAL TO $CONFIG
                #  IF NESTS, SHOULD BE THE ABSOLUTE PARENT CONFIG NAME
                #  (e.g. CONFPAR="trop075" when CONFIG="trop075_nest025")

CONFEXE='AMUcpldom'  # only for nemo.exe

CASE='GNJ002_BM01mv' # should not be too long [>15 char.] otherwise, NEMO file names are affected

YEAR0=<INIT_YEAR>      #- initial year of the long experiment  [ 4 digits ]
MONTH0=<INIT_MONTH>    #- initial month of the long experiment [ 2 digits ]
DAY0=01         #- initial day of the long experiment   [ 2 digits ]

YEAR_MAX=<FINAL_YEAR>   #- stop after $YEAR_MAX is completed

NRUN_MAX=500     #- stop after $NRUN_MAX re-submissions

NDAYS=35         #- Split year by slices of $NDAYS days 

BYMONTH=1     # =0 to stick to a duration of NDAYS (e.g. keep for initial spin up)
              # =1 to cut runs to entire months

WORKDIR="<NEMO_WORKDIR>"

RST_START=<NEMO_START_FROM_RST>   #only concerns the initial step, always restart for further steps
RST_DIR="<NEMO_RST_DIR>"   #initial restart

#STOCKDIR="$SHAREDELMER"  #- restart, output directory
STOCKDIR="<NEMO_WORKDIR>"

EXCHANGE_DIR="<EXCHANGE_DIR>"    # for exchanges between Nemo and ElmerIce

INPUTDIR="<BASE_DIR>/DATA_SETS/NEMO/input/nemo_${CONFIG}"

GLOBAL_SIM='GNJ002'  # ORCA simulation used for BDY, runoff and sss restoring
BDYDIR="${INPUTDIR}/BDY_${GLOBAL_SIM}"  #- input directory for BDYs
SSSDIR="${INPUTDIR}/SSS_${GLOBAL_SIM}"  #- input directory for SSS relaxation (if any)
RNFDIR="${INPUTDIR}/RNF_${GLOBAL_SIM}"  #- input directory for runoff relaxation (if any)

TOPO="BedMachineAntarctica-2019-05-24"  # to use bathy_meter_${CONFPAR}_${TOPO}.nc

#- Netcdf library for small fortran scripts (not for NEMO)
export NC_INC='-I/opt/software/occigen/libraries/netcdf/4.4.0_fortran-4.4.2/hdf5/1.8.17/intel/17.0/openmpi/intel/2.0.1/include'
export NC_LIB='-L/opt/software/occigen/libraries/netcdf/4.4.0_fortran-4.4.2/hdf5/1.8.17/intel/17.0/openmpi/intel/2.0.1/lib -lnetcdf -lnetcdff'

NEMOdir="/home/`whoami`/models/nemo_v3_6_STABLE_r6402/NEMOGCM" # NEMO model directory
XIOSdir="/home/`whoami`/models/xios-1.0" # XIOS directory

FORCINGdir='<BASE_DIR>/DATA_SETS/NEMO/input/FORCING_SETS/DFS5.2'

NZOOM=0  # nb of agrif nests (0 if no agrif nest)

NB_NPROC_XIOS_PER_NODE=1 # Number of core used per xios on each node (should typically be in the 1-3 range).

#=================================================================================
#=================================================================================
# 1- Initialization
#=================================================================================
#=================================================================================

PWDDIR=`pwd`
if [ ! `basename $PWDDIR` == nemo_${CONFIG}_${CASE} ]; then
 echo '~!@#%^&* CHECK CONFIG and CASE in run_nemo.sh >>>>>>>>>> Stop !!'
 exit
fi

export NB_NODES=`echo "${SLURM_NTASKS} / 24" |bc`
export NB_NPROC_IOS=$(( NB_NODES * NB_NPROC_XIOS_PER_NODE ))
export NB_NPROC=$(( SLURM_NTASKS - NB_NPROC_IOS ))

# { unset initiaux 
unset    OMPI_MCA_ess
#
unset    OMPI_MCA_pml
unset    OMPI_MCA_mtl
unset    OMPI_MCA_mtl_mxm_np 
unset    OMPI_MCA_pubsub  
# }

############################################################
##-- create links to executables :

rm -f nemo.exe
if [ $NZOOM -gt 0 ]; then
  ln -s ${NEMOdir}/CONFIG/${CONFEXE}_agrif/BLD/bin/nemo.exe
else
  ln -s ${NEMOdir}/CONFIG/${CONFEXE}/BLD/bin/nemo.exe
fi

rm -f xios_server.exe
ln -s ${XIOSdir}/bin/xios_server.exe

rm -f iodef.xml
if [ $BYMONTH -eq 1 ]; then
  ln -s -v iodef_monthly_daily.xml iodef.xml
else
  ln -s -v iodef_daily.xml iodef.xml
fi

##############################################
##-- define current year and nb of days

if [ $BYMONTH -eq 1 ] && [ $NDAYS -lt 28 ]; then
  echo "~!@#% You need NDAYS > 27 if BYMONTH = 1"
  echo "             >>>>> STOP"
  exit
fi

if [ -f prod_nemo.db ]; then
read NRUN YEAR MONTH DAY NITENDM1 NITENDM1ZOOM << EOF
`tail -1 prod_nemo.db`
EOF
else
echo "1 ${YEAR0} ${MONTH0} ${DAY0} 0" > prod_nemo.db
YEAR=${YEAR0}
MONTH=${MONTH0}
DAY=${DAY0}
NRUN=1
NITENDM1=0  ## last time step of previous run
## add last time step of previous runs on children domains 
## at the end of the line in prod_nemo.db :
for iZOOM in $(seq 1 ${NZOOM})
do
  sed -e "s/$/ 0/g" prod_nemo.db > tmp
  mv tmp prod_nemo.db
done
fi

if [ $YEAR -gt ${YEAR_MAX} ]; then
  echo " "

  #coupling purposes
  echo "GOING TO $EXCHANGE_DIR DIRECTORY"
  cd $EXCHANGE_DIR
  echo "EXTRACTING MELTING FROM NEMO OUTPUT FILES TO MAKE IT READABLE BY ELMERICE"
  sbatch extractXYMelt_fromNEMOsbc.slurm $1 $WORKDIR

  echo "Year greater then YEAR_MAX >>>>>>>> stop !!"

  exit
fi

if [ $NRUN -gt ${NRUN_MAX} ]; then
  echo " "

  #coupling purposes
  echo "GOING TO $EXCHANGE_DIR DIRECTORY"
  cd $EXCHANGE_DIR
  echo "EXTRACTING MELTING FROM NEMO OUTPUT FILES TO MAKE IT READABLE BY ELMERICE"
  sbatch extractXYMelt_fromNEMOsbc.slurm $1 $WORKDIR

  echo "Nb of re-submission greater then NRUN_MAX >>>>>>>> stop !!"
  exit
fi

#####
# adjust nb of days to finish at the end of the year
ISLEAP=`grep nn_leapy namelist_nemo_GENERIC_${CONFIG} | awk '{print $3}'`
if [ $BYMONTH -eq 0 ]; then
if [ ! -f calculate_end_date ]; then
  ifort -o calculate_end_date calculate_end_date.f90
fi
read YEARf MONTHf DAYf NDAYScorr << EOF
`./calculate_end_date $YEAR $MONTH $DAY $NDAYS $ISLEAP`
EOF
if [ $NDAYScorr -ne $NDAYS ]; then
echo "Adjusting run length to finish at the end of current year"
NDAYS=$NDAYScorr
read YEARf MONTHf DAYf NDAYScorr << EOF
`./calculate_end_date $YEAR $MONTH $DAY $NDAYS $ISLEAP`
EOF
fi
elif [ $BYMONTH -eq 1 ]; then
if [ ! -f calculate_end_date_month ]; then
  ifort -o calculate_end_date_month calculate_end_date_month.f90
fi
read YEARf MONTHf DAYf NDAYScorr << EOF
`./calculate_end_date_month $YEAR $MONTH $DAY $NDAYS $ISLEAP`
EOF
if [ $NDAYScorr -ne $NDAYS ]; then
echo "Adjusting run length to finish at the end of a month"
NDAYS=$NDAYScorr
read YEARf MONTHf DAYf NDAYScorr << EOF
`./calculate_end_date_month $YEAR $MONTH $DAY $NDAYS $ISLEAP`
EOF
fi
else
  echo "~!@@#$%^& WRONG VALUE FOR 'BYMONTH'  >>>> stop"
  exit
fi

if [ $NDAYS -gt 0 ]; then
  echo " -> Run duration = $NDAYS days"
  echo " "
else
  echo " -> Run duration < 1   --> STOP !@#%^&"
  echo " (check BYMONTH option)"
  exit
fi

##-- calculate corresponding number of time steps for NEMO:
RN_DT=`grep "rn_rdt " namelist_nemo_GENERIC_${CONFPAR} |grep 'dynamics' |grep 'step' |cut -d '=' -f2 | cut -d '!' -f1 | sed -e "s/ //g"`
NIT000=`echo "$NITENDM1 + 1" | bc`
NITEND=`echo "$NITENDM1 + ${NDAYS} * 86400 / ${RN_DT}" | bc`

echo "****************************************************"
echo "*          NEMO SIMULATION                          "
echo "*   config  $CONFIG                                 "
echo "*   case    $CASE                                   "
echo "*   from    ${DAY}/${MONTH}/${YEAR}                 "
echo "*   to      ${DAYf}/${MONTHf}/${YEARf}              "
echo "*   i.e. step $NIT000 to $NITEND (for mother grid)  "
echo "*                                                   "
echo "*   total number of tasks >>>>> ${SLURM_NTASKS}     "
echo "*   number of xios tasks  >>>>> ${NB_NPROC_IOS}     "
echo "*                                                   "
echo "****************************************************"
echo " "
date
echo " "

#####################################################################
##-- create executable and rebuild namelist to rebuild restart files

rm -f rebuild_nemo.exe
ln -s ${NEMOdir}/TOOLS/REBUILD_NEMO/BLD/bin/rebuild_nemo.exe

###############################################################
##-- edit NEMO's namelist

echo "Editing namelist..."

rm -f namelist_ref namelist_cfg

if [ $NRUN -gt 1 ]; then
  sed -e "s#<RESTNEM>#.true.#g"  namelist_nemo_GENERIC_${CONFPAR} > namelist_ref
  RST=1
else
  if [ $RST_START == 1 ]; then
    sed -e "s#<RESTNEM>#.true.#g"  namelist_nemo_GENERIC_${CONFPAR} > namelist_ref
    RST=1
  else
    sed -e "s#<RESTNEM>#.false.#g" namelist_nemo_GENERIC_${CONFPAR} > namelist_ref
    RST=0
  fi
fi

#if [ $NRUN -gt 1 ] || [ $CONFIG == "trop075" ]; then
#  sed -e "s#<RESTNEM>#.true.#g"  namelist_nemo_GENERIC_${CONFPAR} > namelist_ref
#  RST=1
#else
#  sed -e "s#<RESTNEM>#.false.#g" namelist_nemo_GENERIC_${CONFPAR} > namelist_ref
#  RST=0
#fi

if [ $NRUN -eq 1 ]; then
  if [ $RST_START == 1 ]; then
    sed -e "s#<AAAA>#  0 #g"  namelist_ref > tmp
    mv -f tmp namelist_ref
  else
    sed -e "s#<AAAA>#  0 #g" namelist_ref > tmp
    mv -f tmp namelist_ref
  fi
else
  sed -e "s#<AAAA>#  2 #g"  namelist_ref > tmp
  mv -f tmp namelist_ref
fi
##- Specific treatment for TROP075's restart/initial state:
#if [ $NRUN -eq 1 ] && [ $CONFIG == "trop075" ]; then
#  sed -e "s#<AAAA>#  0 #g" namelist_ref > tmp
#  mv -f tmp namelist_ref
#else
#  sed -e "s#<AAAA>#  2 #g"  namelist_ref > tmp
#  mv -f tmp namelist_ref
#fi

IS_ISCPL=`grep ln_iscpl namelist_nemo_GENERIC_${CONFIG} | awk '{print $3}' |sed -e "s/\.//g"`
if [ $IS_ISCPL == 'true' ]; then
  echo "WARNING : enabling the ice shelf geometry to move (ln_iscpl=true) !!!"
fi

##- create mesh_mask file for first iteration only (or if ice shelf coupling) :
if [ $NRUN -eq 1 ] || [ $IS_ISCPL == 'true' ]; then
  sed -e "s#<MSH>#  1  #g"  namelist_ref > tmp
  mv -f tmp namelist_ref
else
  sed -e "s#<MSH>#  0  #g"  namelist_ref > tmp
  mv -f tmp namelist_ref
fi

echo "s#<CCCC>#${CONFIG}#g ; s#<OOOO>#${CASE}#g ; s#<IIII>#${YEAR0}${MONTH0}${DAY0}#g ; s#<NIT000>#${NIT000}#g ; s#<NITEND>#${NITEND}#g"
sed -e "s#<CCCC>#${CONFIG}#g ; s#<OOOO>#${CASE}#g ; s#<IIII>#${YEAR0}${MONTH0}${DAY0}#g ; s#<NIT000>#${NIT000}#g ; s#<NITEND>#${NITEND}#g" namelist_ref > tmp
mv -f tmp namelist_ref

ln -s namelist_ref namelist_cfg

for iZOOM in $(seq 1 $NZOOM)
do
  echo "Editing ${iZOOM}_namelist..."
  rm -f ${iZOOM}_namelist_ref ${iZOOM}_namelist_cfg
  if [ $NRUN -gt 1 ] || [ $CONFIG == "trop075" ]; then
    sed -e "s#<RESTNEM>#.true.#g"  ${iZOOM}_namelist_nemo_GENERIC_${CONFPAR} > ${iZOOM}_namelist_ref
    RST=1
  else
    sed -e "s#<RESTNEM>#.false.#g" ${iZOOM}_namelist_nemo_GENERIC_${CONFPAR} > ${iZOOM}_namelist_ref
    RST=0
  fi
  ##- Specific treatment for TROP075's restart/initial state:
  if [ $NRUN -eq 1 ] && [ $CONFIG == "trop075" ]; then
    sed -e "s#<AAAA>#  0 #g" ${iZOOM}_namelist_ref > tmp
    mv -f tmp ${iZOOM}_namelist_ref
  else
    sed -e "s#<AAAA>#  2 #g"  ${iZOOM}_namelist_ref > tmp
    mv -f tmp ${iZOOM}_namelist_ref
  fi
  ##- calculate initial and last time step for the child domains :
  ##-- calculate corresponding number of time steps for NEMO:
  RN_DT_ZOOM=`grep "rn_rdt " ${iZOOM}_namelist_nemo_GENERIC_${CONFIG} |cut -d '=' -f2 | cut -d '!' -f1 | sed -e "s/ //g"`
  NIT000_ZOOM=`echo "( ${NITENDM1} * ${RN_DT} / ${RN_DT_ZOOM} ) + 1" | bc`
  NITEND_ZOOM=`echo "( ${NITENDM1} * ${RN_DT} / ${RN_DT_ZOOM} ) + ${NDAYS} * 86400 / ${RN_DT_ZOOM}" | bc`
  ##--
  sed -e "s#<CCCC>#${CONFIG}#g ; s#<OOOO>#${CASE}#g ; s#<IIII>#${YEAR0}0101#g ; s#<NIT000>#${NIT000_ZOOM}#g ; s#<NITEND>#${NITEND_ZOOM}#g" ${iZOOM}_namelist_ref > tmp
  mv -f tmp ${iZOOM}_namelist_ref
  ln -s ${iZOOM}_namelist_ref ${iZOOM}_namelist_cfg
done

rm -f namelist_ice_ref namelist_ice_cfg
cp -p namelist_ice_nemo_GENERIC_${CONFPAR} namelist_ice_ref
ln -s namelist_ice_ref namelist_ice_cfg

#############################################################
###-- prepare script that will be used to compress outputs :

STOCKDIR2=`echo $STOCKDIR |sed -e "s/\//\\\\\\\\\//g"`
WORKDIR2=`echo $WORKDIR  |sed -e "s/\//\\\\\\\\\//g"`

sed -e "s/CCCC/${CONFIG}/g ; s/cccc/${CONFPAR}/g ; s/OOOO/${CASE}/g ; s/SSSS/${STOCKDIR2}/g ; s/WWWW/${WORKDIR2}/g ; s/YYYY/${YEAR}/g ; s/NNNN/${NRUN}/g ; s/ZZZZ/${NZOOM}/g ; s/UUUU/${NITEND}/g" compress_nemo_GENERIC.sh > compress_nemo_${NRUN}.sh

chmod +x compress_nemo_${NRUN}.sh

#=======================================================================================
#=======================================================================================
# 2- Manage links to input files
#=======================================================================================
#=======================================================================================

echo " "
date
echo " "
echo " Linking input files from ${INPUTDIR}"

DATE0=`grep nn_date0 namelist_cfg | head -1 | awk '{print $3}'`
Y0=`echo $DATE0 | cut -c 1-4`
M0=`echo $DATE0 | cut -c 5-6`

YEARm1=`expr $YEAR - 1`
if [ $YEARm1 -lt 1000 ]; then
  YEARm1="0$YEARm1"
fi
# because no data before:
if [ $YEAR -eq $Y0 ]; then
  YEARm1=$Y0
fi
YEARp1=`expr $YEAR + 1`
if [ $YEARp1 -lt 1000 ]; then
  YEARp1="0$YEARp1"
fi

##########
##-- import files that are not time dependent if not already there

rm -f coordinates.nc
ln -s -v ${INPUTDIR}/coordinates_${CONFPAR}.nc coordinates.nc

if [ $CONFPAR == "trop075" ]; then
  rm -f ahmcoef.nc
  ln -s -v ${INPUTDIR}/ahmcoef_${CONFPAR}.nc ahmcoef.nc
fi

EXT_NEMO="<EXT_NEMO>"

rm -f bathy_meter.nc
rm -f isf_draft_meter.nc

if [ $EXT_NEMO == "rel1" ]; then
  ln -s -v ${INPUTDIR}/bathy_meter_${CONFPAR}_${TOPO}.nc bathy_meter.nc
  ln -s -v ${INPUTDIR}/bathy_meter_${CONFPAR}_${TOPO}.nc isf_draft_meter.nc
else
  ln -s -v <EXCHANGE_DIR>/ISF_DRAFT/isf_draft_meter_current.nc bathy_meter.nc
  ln -s -v <EXCHANGE_DIR>/ISF_DRAFT/isf_draft_meter_current.nc isf_draft_meter.nc
fi

rm -f chlorophyll.nc
ln -s -v ${INPUTDIR}/chlorophyll_${CONFPAR}.nc chlorophyll.nc

IS_RNF=`grep " ln_rnf " namelist_cfg | cut -d '=' -f2 | cut -d '!' -f1 | sed -e "s/ //g"`
if [ $IS_RNF == ".true." ]; then
  IS_CLIM=`grep " sn_rnf " namelist_cfg | cut -d ',' -f5 | sed -e "s/ //g"`
  if [ $IS_CLIM == ".true" ]; then
    rm -f runoff.nc runoff_y????.nc
    ln -s -v ${INPUTDIR}/runoff_${CONFPAR}.nc runoff.nc
  else
    for AN in $YEARm1 $YEAR $YEARp1
    do
      rm -f runoff.nc runoff_y${AN}.nc
      if [ -f ${RNFDIR}/runoff_y${AN}_${CONFPAR}.nc ]; then 
        ln -s -v ${RNFDIR}/runoff_y${AN}_${CONFPAR}.nc runoff_y${AN}.nc
      fi
    done
  fi
fi

IS_TRADMP=`grep ln_tradmp namelist_cfg | awk '{print $3}'`
if [ $IS_TRADMP == ".true." ]; then
  rm -f resto.nc
  ln -s -v ${INPUTDIR}/resto_${CONFPAR}.nc resto.nc
fi

# Same for AGRIF NESTS :

for iZOOM in $(seq 1 ${NZOOM})
do

  rm -f ${iZOOM}_coordinates.nc
  ln -s -v ${INPUTDIR}/${iZOOM}_coordinates_${CONFPAR}.nc ${iZOOM}_coordinates.nc

  rm -f ${iZOOM}_bathy_meter.nc
  ln -s -v ${INPUTDIR}/${iZOOM}_bathy_meter_${CONFPAR}.nc ${iZOOM}_bathy_meter.nc

  rm -f ${iZOOM}_isf_draft_meter.nc
  if [ -f ${INPUTDIR}/${iZOOM}_isf_draft_meter_${CONFPAR}.nc ]; then
    ln -s -v ${INPUTDIR}/${iZOOM}_isf_draft_meter_${CONFPAR}.nc ${iZOOM}_isf_draft_meter.nc
  else
    ln -s -v ${INPUTDIR}/${iZOOM}_bathy_meter_${CONFPAR}.nc ${iZOOM}_isf_draft_meter.nc
  fi

  rm -f ${iZOOM}_chlorophyll.nc
  ln -s -v ${INPUTDIR}/${iZOOM}_chlorophyll_${CONFPAR}.nc ${iZOOM}_chlorophyll.nc

  IS_RNF=`grep " ln_rnf " ${iZOOM}_namelist_cfg | cut -d '=' -f2 | cut -d '!' -f1 | sed -e "s/ //g"`
  if [ $IS_RNF == ".true." ]; then
    IS_CLIM=`grep " sn_rnf " ${iZOOM}_namelist_cfg | cut -d ',' -f5 | sed -e "s/ //g"`
    if [ $IS_CLIM == ".true" ]; then
      rm -f ${iZOOM}_runoff.nc
      ln -s -v ${RNFDIR}/${iZOOM}_runoff_${CONFPAR}.nc ${iZOOM}_runoff.nc
    else
      for AN in $YEARm1 $YEAR $YEARp1
      do
        rm -f ${iZOOM}_runoff_y${AN}.nc
        if [ -f ${RNFDIR}/RNF/${iZOOM}_runoff_y${AN}_${CONFPAR}.nc ]; then
          ln -s -v ${RNFDIR}/RNF/${iZOOM}_runoff_y${AN}_${CONFPAR}.nc ${iZOOM}_runoff_y${AN}.nc
        fi
      done
    fi
  fi

  IS_TRADMP=`grep ln_tradmp ${iZOOM}_namelist_cfg | awk '{print $3}'`
  if [ $IS_TRADMP == ".true." ]; then
    rm -f ${iZOOM}_resto.nc
    ln -s -v ${INPUTDIR}/${iZOOM}_resto_${CONFPAR}.nc ${iZOOM}_resto.nc
  fi

done

##########
##-- import Boundary Conditions (BDYs) 
##   (it is assumed that you use a coordinates_bdy.nc file)

rm -f coordinates_bdy.nc
ln -s -v ${INPUTDIR}/coordinates_bdy_${CONFIG}.nc coordinates_bdy.nc

for AN in $YEARm1 $YEAR $YEARp1
do
  for BDYNAM in bdyT_tra bdyU_u2d bdyU_u3d bdyV_u2d bdyV_u3d bdyT_ice bdyT_ssh
  do
    rm -f ${BDYNAM}_y${AN}.nc
    if [ -f ${BDYDIR}/${BDYNAM}_y${AN}_${CONFIG}.nc ]; then
      ln -s -v ${BDYDIR}/${BDYNAM}_y${AN}_${CONFIG}.nc ${BDYNAM}_y${AN}.nc
    fi
  done
done

## TIDES :
for HARM in `grep clname namelist_nemo_GENERIC_${CONFIG} | awk '{print $3}' |sed -e "s/'//g"`
do
  ln -s -v ${INPUTDIR}/BDY_TIDES/bdytide_${CONFIG}_${HARM}_grid_T.nc bdytide_${HARM}_grid_T.nc
  ln -s -v ${INPUTDIR}/BDY_TIDES/bdytide_${CONFIG}_${HARM}_grid_U.nc bdytide_${HARM}_grid_U.nc
  ln -s -v ${INPUTDIR}/BDY_TIDES/bdytide_${CONFIG}_${HARM}_grid_V.nc bdytide_${HARM}_grid_V.nc
done

### PATCH ##
#rm -f bdyT_ice.nc
#ln -s -v ${BDYDIR}/bdyT_ice_y${YEAR}_${CONFIG}.nc bdyT_ice.nc

##########
##-- atmospheric forcing

FORDTA=`basename $FORCINGdir`
rm -f w_bilin.nc w_bicubic.nc
ln -s -v ${INPUTDIR}/weights_bilin_${FORDTA}_${CONFIG}.nc    w_bilin.nc
ln -s -v ${INPUTDIR}/weights_bicubic_${FORDTA}_${CONFIG}.nc  w_bicubic.nc

# The following only work if in thefollowing order: namsbc_clio, namsbc_core, namsbc_mfs
IS_BLK_CORE=`grep " ln_blk_core " namelist_cfg | cut -d '=' -f2 | cut -d '!' -f1 |sed -e "s/ //g"`
if [ $IS_BLK_CORE == ".true." ]; then
  LINE_BLK_CORE=`grep -n namsbc_core namelist_cfg | tail -1 | cut -d ':' -f1`
  for NAMAT in sn_wndi sn_wndj sn_qsr sn_qlw sn_tair sn_humi sn_prec sn_snow sn_tdif
  do
    ATM_FILE=`awk "/${NAMAT}/ && NR >= ${LINE_BLK_CORE}" namelist_cfg | cut -d "'" -f2 | head -1`
    IS_CLIM=`awk "/${NAMAT}/ && NR >= ${LINE_BLK_CORE}" namelist_cfg | cut -d ',' -f5 | sed -e "s/ //g" | head -1`
    if [ $IS_CLIM == ".true." ]; then
      rm -f ${ATM_FILE}.nc
      if [ -f ${FORCINGdir}/${ATM_FILE}.nc ]; then
        ln -s -v ${FORCINGdir}/${ATM_FILE}.nc
      fi
    else
      for AN in  $YEARm1 $YEAR $YEARp1 
      do
        rm -f ${ATM_FILE}_y${AN}.nc
        if [ -f ${FORCINGdir}/${ATM_FILE}_y${AN}.nc ]; then
          ln -s -v ${FORCINGdir}/${ATM_FILE}_y${AN}.nc ${ATM_FILE}_y${AN}.nc
        fi
      done
    fi
  done
fi

for iZOOM in $(seq 1 ${NZOOM})
do

  rm -f ${iZOOM}_w_bilin.nc ${iZOOM}_w_bicubic.nc
  ln -s ${INPUTDIR}/${iZOOM}_weights_bilin_${FORDTA}_${CONFIG}.nc    ${iZOOM}_w_bilin.nc
  ln -s ${INPUTDIR}/${iZOOM}_weights_bicubic_${FORDTA}_${CONFIG}.nc  ${iZOOM}_w_bicubic.nc

  IS_BLK_CORE=`grep " ln_blk_core " ${iZOOM}_namelist_cfg | cut -d '=' -f2 | cut -d '!' -f1 |sed -e "s/ //g"`
  if [ $IS_BLK_CORE == ".true." ]; then
    LINE_BLK_CORE=`grep -n namsbc_core ${iZOOM}_namelist_cfg | tail -1 | cut -d ':' -f1`
    for NAMAT in sn_wndi sn_wndj sn_qsr sn_qlw sn_tair sn_humi sn_prec sn_snow sn_tdif
    do
      ATM_FILE=`awk "/${NAMAT}/ && NR >= ${LINE_BLK_CORE}" ${iZOOM}_namelist_cfg | cut -d "'" -f2 | head -1`
      IS_CLIM=`awk "/${NAMAT}/ && NR >= ${LINE_BLK_CORE}" ${iZOOM}_namelist_cfg | cut -d ',' -f5 | sed -e "s/ //g" | head -1`
      if [ $IS_CLIM == ".true." ]; then
        rm -f ${iZOOM}_${ATM_FILE}.nc
        if [ -f ${FORCINGdir}/${iZOOM}_${ATM_FILE}.nc ]; then
          ln -s -v ${FORCINGdir}/${iZOOM}_${ATM_FILE}.nc
        fi
      else
        for AN in  $YEARm1 $YEAR $YEARp1 
        do
          rm -f ${iZOOM}_${ATM_FILE}_y${AN}.nc
          if [ -f ${FORCINGdir}/${iZOOM}_${ATM_FILE}_y${AN}.nc ]; then
            ln -s -v ${FORCINGdir}/${iZOOM}_${ATM_FILE}_y${AN}.nc ${iZOOM}_${ATM_FILE}_y${AN}.nc
          fi
        done
      fi
    done
  fi

done

##########
##-- SSS restoring if any :

SSSREL=`grep " nn_sssr " namelist_cfg | awk '{print $3}'`

if [ $SSSREL -gt 0 ]; then
  echo "Run with SSS relaxation : nn_sssr = $SSSREL"
  SSS_FILE=`grep sn_sss namelist_cfg | cut -d "'" -f2 | head -1`
  IS_CLIM=`grep sn_sss namelist_cfg | cut -d ',' -f5 | sed -e "s/ //g" | head -1`
  if [ $IS_CLIM == ".true." ]; then
    rm -f ${SSS_FILE}.nc
    if [ -f ${INPUTDIRDIR}/sss.nc ]; then
      ln -s -v ${INPUTDIRDIR}/sss.nc ${SSS_FILE}.nc
    fi
  else
    for AN in  $YEARm1 $YEAR $YEARp1 
    do
      rm -f ${SSS_FILE}_y${AN}.nc
      if [ -f ${SSSDIR}/sss_y${AN}_${CONFIG}.nc ]; then
        ln -s -v ${SSSDIR}/sss_y${AN}_${CONFIG}.nc ${SSS_FILE}_y${AN}.nc
      fi
    done
  fi
else
  echo "Run without SSS relaxation ( nn_sssr = $SSSREL )"
fi

##########
##-- Initial state or Restart

rm -f restart.nc restart_ice.nc #restart.obc
rm -f dta_temp_y????m??.nc dta_sal_y????m??.nc dta_temp.nc dta_sal.nc
RSTN=`grep "from a restart file" namelist_cfg | awk '{print $3}' | sed -e "s/\.//g"`
NIT_RST=${NITENDM1}
if [ $RSTN == "true" ]; then
  if [ $NIT_RST -eq 0 ]; then
    # for restart
    NRST=`ls -1 $RST_DIR/restart_????????.nc | tail -1`
    ln -s -v $NRST restart.nc
    # for restart ice
    NRSTICE=`ls -1 $RST_DIR/restart_ice_????????.nc | tail -1`
    ln -s -v $NRSTICE restart_ice.nc
    #ln -s -v ${INPUTDIR}/${CONFPAR}_restart_00000000.nc restart.nc
    #ln -s -v ${INPUTDIR}/${CONFPAR}_restart_ice_00000000.nc restart_ice.nc
  else
    if [ ! -f restart_${NIT_RST}.nc ]; then
      echo "Copy ocean restart file from ${STOCKDIR}/restart/nemo_${CONFIG}-${CASE}"
      cp -p ${STOCKDIR}/restart/nemo_${CONFIG}_${CASE}/restart_${NIT_RST}.nc .
    fi
    ln -s -v restart_${NIT_RST}.nc   restart.nc
    if [ ! -f restart_ice_${NIT_RST}.nc ]; then
      echo "Copy ice restart file from ${STOCKDIR}/restart/nemo_${CONFIG}-${CASE}"
      cp -p ${STOCKDIR}/restart/nemo_${CONFIG}_${CASE}/restart_ice_${NIT_RST}.nc .
    fi
    ln -s -v restart_ice_${NIT_RST}.nc   restart_ice.nc
  fi
else
  echo "Not in restart mode -> import initial T,S state"
  ln -s -v ${INPUTDIR}/dta_temp_${CONFPAR}_${GLOBAL_SIM}_y${Y0}m${M0}.nc dta_temp.nc
  ln -s -v ${INPUTDIR}/dta_sal_${CONFPAR}_${GLOBAL_SIM}_y${Y0}m${M0}.nc  dta_sal.nc
fi

for iZOOM in $(seq 1 ${NZOOM})
do

  rm -f ${iZOOM}_restart.nc ${iZOOM}_restart_ice.nc
  rm -f ${iZOOM}_dta_temp_y????m??.nc ${iZOOM}_dta_sal_y????m??.nc ${iZOOM}_dta_temp.nc ${iZOOM}_dta_sal.nc
  RSTN=`grep "from a restart file" namelist_cfg | awk '{print $3}' | sed -e "s/\.//g"`
  NIT_RST=${NITENDM1}
  if [ $RSTN == "true" ]; then
    if [ $NIT_RST -eq 0 ]; then
      ln -s -v ${INPUTDIR}/${iZOOM}_${CONFPAR}_restart_00000000.nc ${iZOOM}_restart.nc
      ln -s -v ${INPUTDIR}/${iZOOM}_${CONFPAR}_restart_ice_00000000.nc ${iZOOM}_restart_ice.nc
    else
      if [ ! -f ${iZOOM}_restart_${NIT_RST}.nc ]; then
        echo "Copy zoom ocean restart file from ${STOCKDIR}/restart/nemo_${CONFIG}-${CASE}"
        cp -p ${STOCKDIR}/restart/nemo_${CONFIG}_${CASE}/${iZOOM}_restart_${NIT_RST}.nc .
      fi
      ln -s -v ${iZOOM}_restart_${NIT_RST}.nc   ${iZOOM}_restart.nc
      if [ ! -f ${iZOOM}_restart_ice_${NIT_RST}.nc ]; then
        echo "Copy zoom ice restart file from ${STOCKDIR}/restart/nemo_${CONFIG}-${CASE}"
        cp -p ${STOCKDIR}/restart/nemo_${CONFIG}_${CASE}/${iZOOM}_restart_ice_${NIT_RST}.nc .
      fi
      ln -s -v ${iZOOM}_restart_ice_${NIT_RST}.nc   ${iZOOM}_restart_ice.nc
    fi
  else
    echo "Not in restart mode -> import initial T,S state for the zoom"
    ln -s -v ${iZOOM}_dta_temp_${CONFPAR}_${GLOBAL_SIM}_y${Y0}m${M0}.nc ${iZOOM}_dta_temp_y${Y0}m${M0}.nc
    ln -s -v ${iZOOM}_dta_sal_${CONFPAR}_${GLOBAL_SIM}_y${Y0}m${M0}.nc  ${iZOOM}_dta_sal_y${Y0}m${M0}.nc
  fi

done

echo " "
echo "Import (links+copy) of input files completed."
echo " "
echo "Launching the long nemo simulation"
echo " "

#=======================================================================================
#=======================================================================================
# 3- Run script
#=======================================================================================
#=======================================================================================

rm -f app.conf
echo "0-$(( NB_NPROC_IOS - 1 )) xios_server.exe"          >  app.conf
echo "${NB_NPROC_IOS}-$(( SLURM_NTASKS - 1 )) nemo.exe "  >> app.conf

date
echo " "

srun --mpi=pmi2  -m cyclic \
    --cpu_bind=map_cpu:0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23\
    --multi-prog  ./app.conf

#srun --mpi=pmi2 --multi-prog  ./app.conf

echo " "
date
echo " "

##-- Rebuild mesh_mask, output.abbort and output.init :

for STUF in mesh_mask output.abort output.init
do
  if [ -f ${STUF}_0000.nc ]; then
    ./rebuild.sh ${STUF}
    for iZOOM in $(seq 1 ${NZOOM})
    do
      if [ -f ${iZOOM}_${STUF}_0000.nc ]; then
        ./rebuild.sh ${STUF} $iZOOM
      fi
    done
  fi
done

##-- export and compress output files:

if [ ! -d OUTPUT_${NRUN} ]; then
  mkdir OUTPUT_${NRUN}
fi

if [ $IS_ISCPL == 'true' ]; then
  if [ ! -f mesh_mask.nc ]; then
    echo '~!@#$%^&* ERROR: mesh_mask.nc not created >>> stop because compulsory in iscpl mode...'
    exit 
  else
    mv mesh_mask.nc OUTPUT_${NRUN}/mesh_mask_${YEAR}${MONTH}${DAY}.nc
  fi
fi

mv -f ${CONFIG}-${CASE}_[1-5][d-m]_*nc OUTPUT_${NRUN}/.
mv -f namelist_ref                     OUTPUT_${NRUN}/namelist.${NRUN}
mv -f namelist_ice_ref                 OUTPUT_${NRUN}/namelist_ice.${NRUN}
mv -f ocean.output                     OUTPUT_${NRUN}/ocean.output.${NRUN}
rm -f namelist_ice_cfg namelist_cfg

for iZOOM in $(seq 1 ${NZOOM})
do
  mv -f ${iZOOM}_${CONFPAR}-${CASE}_[1-5][d-m]_*nc  OUTPUT_${NRUN}/.
  mv -f ${iZOOM}_namelist_ref                       OUTPUT_${NRUN}/${iZOOM}_namelist.${NRUN}
  mv -f ${iZOOM}_namelist_ice_ref                   OUTPUT_${NRUN}/${iZOOM}_namelist_ice.${NRUN}
  mv -f ${iZOOM}_ocean.output                       OUTPUT_${NRUN}/${iZOOM}_ocean.output.${NRUN}
  rm -f ${iZOOM}_namelist_ice_cfg ${iZOOM}_namelist_cfg
done

##-- for coupling purposes, copy SBC files in the output_sbc directory

if [ ! -d output_sbc ]; then
  mkdir output_sbc
fi

echo "Copying SBC file, fwfisf variable in output_sbc directory"

ncks -O -v fwfisf OUTPUT_${NRUN}/${CONFIG}-${CASE}_1m_*_SBC.nc output_sbc/${CONFIG}-${CASE}_1m_fwfisf_@${NRUN}.nc
#ncks -O -v fwfisf OUTPUT_${NRUN}/${CONFIG}-${CASE}_1d_*_SBC.nc output_sbc/${CONFIG}-${CASE}_1d_fwfisf_${NRUN}.nc

#echo "Making the average over a month"
#ncra -O -F -d time_counter,1,40 output_sbc/${CONFIG}-${CASE}_1d_fwfisf_${NRUN}.nc output_sbc/${CONFIG}-${CASE}_1d_mean_fwfisf_@${NRUN}.nc
#rm output_sbc/${CONFIG}-${CASE}_1d_fwfisf_${NRUN}.nc

## used to know how many multiple output files are created (in xios mode "multiple_file")
echo "xxx $NB_NPROC_IOS xios_server.exe xxx" > OUTPUT_${NRUN}/app.copy

##########################################################
##-- prepare next run if every little thing went all right

NTEST_O=`ls -1 OUTPUT_${NRUN}/${CONFIG}-${CASE}_[1-5][d-m]_*nc |wc -l`
NTEST_R=`ls -1 ${CONFIG}-${CASE}_*_restart_*.nc |wc -l`

FILE_TEST=`ls -1 ${CONFIG}-${CASE}_*_SBC.nc | tail -1`
NBNAN=`ncdump -v fwfisf ${FILE_TEST} |grep NaN |wc -l`
if [ $NBNAN -gt 0 ]; then
  echo '**************************************************'
  echo '      NaNs were found in the output files         '
  echo '             >>>>>>>>>   stopping here            '
  echo '**************************************************'
  exit
fi 

if [ ${NTEST_O} -gt 0 ] && [ ${NTEST_R} -gt 0 ] && [ $NBNAN -eq 0 ]; then

  sbatch ./compress_nemo_${NRUN}.sh

  ##-- write last restart time step of mother grid in prod_nemo.db:
  LAST_RESTART_NIT=`ls -lrt ${CONFIG}-${CASE}_*_restart_*.nc |tail -1 | sed -e "s/${CONFIG}-${CASE}//g" | cut -d '_' -f2`
  echo " "
  echo "Last restart created at ocean time step ${LAST_RESTART_NIT}"
  echo "  ---> writting this date in prod_nemo.db"
  echo " "
  echo "$LAST_RESTART_NIT" > restart_nit.txt
  ##-- add last restart time step on chidren grids (at the end of last line in prod_nemo.db):
  for iZOOM in $(seq 1 ${NZOOM})
  do
    LAST_RESTART_NIT_ZOOM=`ls -lrt ${iZOOM}_${CONFIG}-${CASE}_*_restart_*.nc |tail -1 | sed -e "s/${iZOOM}_${CONFIG}-${CASE}//g" | cut -d '_' -f2`
    sed -e "`wc -l prod_nemo.db|cut -d ' ' -f1`s/$/ ${LAST_RESTART_NIT_ZOOM}/g" prod_nemo.db > tmp
    mv tmp prod_nemo.db
    echo " "
    echo "Last restart created for zoom nb $iZOOM at ocean time step ${LAST_RESTART_NIT_ZOOM}"
    echo "  ---> writting this date in prod_nemo.db"
    echo " "
    echo "$LAST_RESTART_NIT_ZOOM" > ${iZOOM}_restart_nit.txt
  done

  echo " "
  date
  echo " "

  ## rebuild restart files for mother grid :
  FILEBASE=`ls -1 ${CONFIG}-${CASE}_[0-9]???????_restart_0000.nc | sed -e "s/_0000.nc//g"`
  ./rebuild.sh $FILEBASE
  mv  ${FILEBASE}.nc restart_${LAST_RESTART_NIT}.nc
  FILEBASE=`ls -1 ${CONFIG}-${CASE}_[0-9]???????_restart_ice_0000.nc | sed -e "s/_0000.nc//g"`
  ./rebuild.sh $FILEBASE
  mv  ${FILEBASE}.nc restart_ice_${LAST_RESTART_NIT}.nc

  ## rebuild restart files for child domains :
  for iZOOM in $(seq 1 ${NZOOM})
  do
    LAST_RESTART_NIT_ZOOM=`cat ${iZOOM}_restart_nit.txt`
    FILEBASE=`ls -1 ${iZOOM}_${CONFIG}-${CASE}_[0-9]???????_restart_0000.nc | sed -e "s/_0000.nc//g"`
    ./rebuild.sh $FILEBASE $iZOOM
    mv  ${FILEBASE}.nc ${iZOOM}_restart_${LAST_RESTART_NIT_ZOOM}.nc
    FILEBASE=`ls -1 ${iZOOM}_${CONFIG}-${CASE}_[0-9]???????_restart_ice_0000.nc | sed -e "s/_0000.nc//g"`
    ./rebuild.sh $FILEBASE $iZOOM
    mv  ${FILEBASE}.nc ${iZOOM}_restart_ice_${LAST_RESTART_NIT_ZOOM}.nc
  done

  echo " "
  date
  echo " "

  # prepare initial state for following iteration:
  NRUNm1=$NRUN 
  NRUNm2=`expr $NRUN - 1`
  NRUN=`expr $NRUN + 1`
  TMPTMP="${LAST_RESTART_NIT}"
  for iZOOM in $(seq 1 ${NZOOM})
  do
    LAST_RESTART_NIT_ZOOM=`cat ${iZOOM}_restart_nit.txt`
    TMPTMP="${TMPTMP} ${LAST_RESTART_NIT_ZOOM}"
  done
  echo "${NRUN} ${YEARf} ${MONTHf} ${DAYf} ${TMPTMP}" >> prod_nemo.db    ## new line

  # clean WORKDIR:
  YEARm2=`expr $YEAR - 2`
  if [ $YEARm2 -lt 1000 ]; then
    YEARm2="0$YEARm2"
  fi
  if [ $IS_BLK_CORE == ".true." ]; then
    for NAMAT in sn_wndi sn_wndj sn_qsr sn_qlw sn_tair sn_humi sn_prec sn_snow sn_tdif
    do
      ATM_FILE=`grep $NAMAT namelist_cfg | cut -d "'" -f2 | head -1`
      IS_CLIM=`grep $NAMAT namelist_cfg | cut -d ',' -f5 | sed -e "s/ //g" | head -1`
      if [ $IS_CLIM == ".false." ]; then
        rm -f ${ATM_FILE}_y${YEARm2}.nc
        for iZOOM in $(seq 1 $NZOOM)
        do
           rm -f ${iZOOM}_${ATM_FILE}_y${YEARm2}.nc
        done
      fi
    done
  fi

  rm -f bdyT_tra_y${YEARm2}.nc
  rm -f bdyU_u?d_y${YEARm2}.nc
  rm -f bdyV_u?d_y${YEARm2}.nc
  rm -f bdyT_ice_y${YEARm2}.nc
  rm -f bdyT_ssh_y${YEARm2}.nc


elif [ $NBNAN -gt 0 ]; then

  echo '**************************************************'
  echo '*     NaNs were found in the output files        *'
  echo '*            >>>>>>>>>   stopping here           *'
  echo '**************************************************'
  exit

else

  echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  echo '!@#$%^&* BIG PROBLEM : no output or no restart files created for NEMO !!  >>>>>>> STOP '
  echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  exit

fi

##########################################################
##-- launch next year of the experiment

#Arguments
#1 EXT_NEMO (arguments because it is changed between experiments rel1, rel2, pro)

ELMER_WORKDIRpro=<ELMER_WORKDIRpro>

if [ $1 == 'rel1' ] || [ $1 == 'rel2' ]; then
 
  sbatch run_nemo.sh $1

elif [ $1 == 'pro' ]; then

  echo " "
  echo "GOING TO $EXCHANGE_DIR DIRECTORY"
  cd $EXCHANGE_DIR
  echo "EXTRACTING MELTING FROM NEMO OUTPUT FILE TO MAKE IT READABLE BY ELMERICE"
  sbatch extractXYMelt_fromNEMOsbc.slurm $1 $WORKDIR $NRUNm1

fi

echo " "
date

