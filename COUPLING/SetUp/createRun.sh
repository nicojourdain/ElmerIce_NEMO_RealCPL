#!/bin/bash

export LANG=C

#BASE_DIR=${PWD}/../../    #need to be similar to BASE_DIR in BASE_DIR/PARAMETERS/$1.IN
BASE_DIR=/scratch/shared/egige60/ElmerIce_NEMO_RealCPL

#This is for testing purposes
#If you want to use it, your need to name your COUPLING directory such as ../COUPLING@blabla/..
#Use the @ as a delimiter
playground=@`echo $PWD | cut -d'@' -f2 | cut -d'/' -f1`
if [ $playground == '@' ]; then
  playground=
fi

echo "############################################################"
echo "COUPLING"
echo "The family argument is mandatory (e.g. BEDMACHINE_FINE_MELTPLUS)"
echo "in any case you need a FAMILY in the PARAMETERS directory"
echo "#############################################################"

#test the number of arguments
echo $# arguments
if [ $# -ne 1 ];then
  echo "illegal number of arguments"
  echo "ABORT"
  exit
fi

#test if $1 file exists in '$BASE_DIR/PARAMETERS/'
fileIN=$BASE_DIR/PARAMETERS/$1
if [ ! -e $fileIN ]; then
  echo "the family file does not exist"
  echo "there should be a $1 file"
  echo "in your ../../PARAMETERS directory"
  echo "ABORT"
  exit
fi
FAMILY=`echo $1 | cut -d"." -f1`

###################################################
## User's choices :
###################################################

# DO IT BLOCK BY BLOCK
# Choose the ith block and everything behind will be done
# 0st block : Preliminary stuff
# 1st block : NEMOrel1
# 2nd block : ELMERICErel
# 3rd block : NEMOrel2
# 4th block : NEMO-ELMERICE pro
stageblock=4

#For NEMO, the period of the simulations is always 1 month
#For ELMERICE, 1 year for the relaxation and 1 month for coupling

INIT_YEAR=1995              #initial year of calculation used by Nemo relaxations and prognostic
INIT_MONTH=03

COUPLING_PERIOD=0.083333333 # coupling period in years, used in prognostics by Nemo and Elmerice
COUPLING_NBYEARS=5          # for prognostic run (4th block)

#ELMERICE
ELMER_RELAX_NBYEARS=20      #number of years in total for ELMERICE relaxation

ELMER_RELAX_DT=0.1          #time step (in years) for ELMERICE in stand alone jobs
ELMER_CPL_DT=0.0083333333   #time step (in years) for ELMERICE in coupled jobs (should divide COUPLING_PERIOD, e.g. COUPLING_PERIOD/10 to have 10 Elmer time step in each coupling window)

ELMER_PATH_INVERSE=$BASE_DIR/INVERSION$playground/Results                               #the inverse result that is listed in $1.IN

ELMER_WORKDIRrel=$BASE_DIR/COUPLING$playground/SetUp/RUNS_$FAMILY/RELAXATION/ELMERICE   #ElmerIce working directory for relaxation step
ELMER_WORKDIRpro=$BASE_DIR/COUPLING$playground/SetUp/RUNS_$FAMILY/PROGNOSTIC/ELMERICE   #ElmerIce working directory for prognostic step

ELMER_HOME=/scratch/cnt0021/egi6035/SHARED/local                          #ElmerIce main code from CSC
ELMERICELGGElibs=$BASE_DIR/otherStuff/sourcesElmerIce/ElmerIceLGGE/lib    #ElmerIce LGGE from Renater + personal sources

#NEMO
CONFIG=AMUXL025
CASE=GNJ002_BM01mv

NEMO_RELAX1_NBYEARS=1         #nemo first period of relaxation in years
MONTHS_TO_SAVE_RELAX1=12      #nemo how many month to be saved for next step in months

NEMO_RELAX2_NBYEARS=1        #nemo second period of relaxation in years
MONTHS_TO_SAVE_RELAX2=12     #nemo how many month to be saved for next step in months

NEMO_CONF='AMUcpl'    #the NEMO instal that is used

NEMO_WORKDIRrel1=$BASE_DIR/COUPLING$playground/SetUp/RUNS_$FAMILY/RELAXATION/NEMO_rel1/nemo_AMUXL025_GNJ002_BM01mv    #Nemo working directory for relaxation step
NEMO_WORKDIRrel2=$BASE_DIR/COUPLING$playground/SetUp/RUNS_$FAMILY/RELAXATION/NEMO_rel2/nemo_AMUXL025_GNJ002_BM01mv    #Nemo working directory for relaxation step
NEMO_WORKDIRpro=$BASE_DIR/COUPLING$playground/SetUp/RUNS_$FAMILY/PROGNOSTIC/nemo_AMUXL025_GNJ002_BM01mv    #Nemo working directory for prognostic step

TEMPLATE_DIR=$BASE_DIR/COUPLING$playground/SetUp/TEMPLATE
EXCHANGE_DIR=$BASE_DIR/COUPLING$playground/SetUp/RUNS_$FAMILY/EXCHANGE

#Files for NEMO
#BATHY_FILE=$BASE_DIR/DATA_SETS/NEMO/input/nemo_AMUXL025/-2019-05-24.nc    #initial bathyaetry & ice shelf draft
#ISF_DRAFT_FILE=$BATHY_FILE

###################################################
## End of User's choices
###################################################

echo "############################################################"
echo ""
echo "Creating coupled simulation $FAMILY"
echo "=> Coupling period of $COUPLING_PERIOD year is used"
echo "=> Be confident and it may work"
echo "=> Other stuff will follow..."
echo ""
echo "############################################################"

echo "CREATING DIRECTORIES"

#create directories
mkdir -vp RUNS_$FAMILY
mkdir -vp RUNS_$FAMILY/RELAXATION
mkdir -vp $ELMER_WORKDIRrel
mkdir -vp $NEMO_WORKDIRrel1
mkdir -vp $NEMO_WORKDIRrel2
mkdir -vp RUNS_$FAMILY/PROGNOSTIC
mkdir -vp $ELMER_WORKDIRpro
mkdir -vp $NEMO_WORKDIRpro
# for exchanging informations between ElmerIce and Nemo
mkdir -vp $EXCHANGE_DIR
mkdir -vp $EXCHANGE_DIR/MELT_RATES
mkdir -vp $EXCHANGE_DIR/ISF_DRAFT

cp -v createRun.sh RUNS_$FAMILY/createRun_$FAMILY.sh   #copy the main script

cat $BASE_DIR/PARAMETERS/$1 | sed -e "s#<BASE_DIR>#$BASE_DIR#g" > RUNS_$FAMILY/$1
source RUNS_$FAMILY/$1

echo "############################################################"
echo ""
echo "PREPARING RELAXATION STEP"
echo ""
echo "############################################################"

echo "PREPARING NEMO FIRST RELAXATION"
echo "GOING TO $NEMO_WORKDIRrel1 DIRECTORY"
cd $NEMO_WORKDIRrel1

NEMO_START_FROM_RST=0       #initially set to 0
NEMO_RST_DIR=    #should link the last restart file from nemo init simulations

EXT_NEMO=rel1
SBATCH_ID=$EXT_NEMO

cp -vr $TEMPLATE_DIR/NEMO/nemo_AMUXL025_GNJ002_BM01mv/* .    #all nemo scripts, already executables

NEMO_RELAX1_FINAL_YEAR=$(expr $INIT_YEAR+$NEMO_RELAX1_NBYEARS | bc)

cat $TEMPLATE_DIR/NEMO/nemo_AMUXL025_GNJ002_BM01mv/run_nemo.sh | sed -e "s#<SBATCH_ID>#$SBATCH_ID#g"\
  -e "s#<INIT_YEAR>#$INIT_YEAR#g"\
  -e "s#<INIT_MONTH>#$INIT_MONTH#g"\
  -e "s#<FINAL_YEAR>#$NEMO_RELAX1_FINAL_YEAR#g"\
  -e "s#<NEMO_WORKDIR>#$NEMO_WORKDIRrel1#g"\
  -e "s#<NEMO_START_FROM_RST>#$NEMO_START_FROM_RST#g"\
  -e "s#<NEMO_RST_DIR>#$NEMO_RST_DIR#g"\
  -e "s#<EXCHANGE_DIR>#$EXCHANGE_DIR#g"\
  -e "s#<BASE_DIR>#$BASE_DIR#g"\
  -e "s#<EXT_NEMO>#$EXT_NEMO#g"\
  -e "s#<ELMER_WORKDIRpro>#$ELMER_WORKDIRpro#g" > run_nemo.sh

####################################################
echo "PREPARING ELMERICE RELAXATION"
echo "GOING TO $ELMER_WORKDIRrel DIRECTORY"
cd $ELMER_WORKDIRrel

echo "SYNCING MESH USED AT THE INVERSION STEP"

rsync -av --exclude=*vtu --exclude=*result* $ELMER_PATH_INVERSE/RUNS_${FAMILY}/run_INIT_OPTIM/mesh_$nbpart .

echo "PREPARING INVERSION STUFF"

cat $TEMPLATE_DIR/TOOLS/ImportResultsInversion.py | sed -e "s#<Ga>#$Gastr#g"\
  -e "s#<rkcg>#$rkcg#g"\
  -e "s#<rkdhdt>#$rkdhdt#g"\
  -e "s#<ELMER_PATH_INVERSE>#$ELMER_PATH_INVERSE#g"\
  -e "s#<nbpart>#$nbpart#g" > ImportResultsInversion.py

cp -v $TEMPLATE_DIR/SLURM/ImportResultsInversion.slurm .

CWL_DATA=Cwl_OPTIM_$invstate".dat"
ETA_DATA=Eta_OPTIM_$invstate".dat"

echo "PREPARING RUN0"

ln -sf $THICKNESS_DATA
ln -sf $BEDROCK_DATA
ln -sf $SMB_CLIM_DATA
THICKNESS_DATAsif=`echo $THICKNESS_DATA | rev | cut -d"/" -f1 | rev`
BEDROCK_DATAsif=`echo $BEDROCK_DATA | rev | cut -d"/" -f1 | rev`
SMB_CLIM_DATAsif=`echo $SMB_CLIM_DATA | rev | cut -d"/" -f1 | rev`

cat $TEMPLATE_DIR/SIF/RUN0.sif | sed -e "s#<nbpart>#$nbpart#g"\
  -e "s#<RHOI_SI>#$rhoi_SI#g"\
  -e "s#<RHOW_SI>#$rhow_SI#g"\
  -e "s#<GRAVITY_SI>#$gravity_SI#g"\
  -e "s#<YEARTOSEC>#$yeartosec#g"\
  -e "s#<ZSL>#$zsl#g"\
  -e "s#<N>#$n#g"\
  -e "s#<THICKNESS_DATA>#$THICKNESS_DATAsif#g"\
  -e "s#<BEDROCK_DATA>#$BEDROCK_DATAsif#g"\
  -e "s#<nameThicknessVar>#$nameThicknessVar#g"\
  -e "s#<nameBedVar>#$nameBedVar#g"\
  -e "s#<SMB_CLIM_DATA>#$SMB_CLIM_DATAsif#g"\
  -e "s#<CWL_DATA>#$CWL_DATA#g"\
  -e "s#<ETA_DATA>#$ETA_DATA#g" > RUN0.sif

nbnode=$(expr $nbpart/24 | bc)    #nbpart thus should be a multiple of 24

cat $TEMPLATE_DIR/SLURM/RUN0.slurm | sed -e "s#<nbpart>#$nbpart#g"\
  -e "s#<ELMER_HOME>#$ELMER_HOME#g"\
  -e "s#<ELMERICELGGElibs>#$ELMERICELGGElibs#g"\
  -e "s#<nbnode>#$nbnode#g" > RUN0.slurm

cp -v $TEMPLATE_DIR/SIF/LINEAR_SOLVER.txt .

echo "PREPARING RUNii"
cd $ELMER_WORKDIRrel

SBATCH_ID=rel

imin=1
imax=$ELMER_RELAX_NBYEARS

for ((ii=$imin ; ii<=$imax ; ii++))
do 

  echo "Preparing ElmerIce $ii st simulation"

  ln -sf $SMB_ANOM_DATA
  ln -sf $EXCHANGE_DIR/MELT_RATES/$MELT_RATES_NEMO_XY
  SMB_ANOM_DATAsif=`echo $SMB_ANOM_DATA | rev | cut -d"/" -f1 | rev`
  MELT_RATES_NEMO_XYsif=`echo $MELT_RATES_NEMO_XY | rev | cut -d"/" -f1 | rev`

  cat $TEMPLATE_DIR/SIF/RUNii.sif | sed -e "s#<ID>#$ii#g"\
    -e "s#<ID-1>#$((ii-1))#g"\
    -e "s#<RHOI_SI>#$rhoi_SI#g"\
    -e "s#<RHOW_SI>#$rhow_SI#g"\
    -e "s#<GRAVITY_SI>#$gravity_SI#g"\
    -e "s#<YEARTOSEC>#$yeartosec#g"\
    -e "s#<ZSL>#$zsl#g"\
    -e "s#<N>#$n#g"\
    -e "s#<nbpart>#$nbpart#g"\
    -e "s#<ELMER_DT>#$ELMER_RELAX_DT#g"\
    -e "s#<SMB_ANOM_DATA>#$SMB_ANOM_DATAsif#g"\
    -e "s#<MELT_RATES_NEMO_XY>#$MELT_RATES_NEMO_XYsif#g" > RUN$ii.sif

  cat $TEMPLATE_DIR/SLURM/RUNii.slurm | sed -e "s#<SBATCH_ID>#$SBATCH_ID#g"\
    -e "s#<ELMER_HOME>#$ELMER_HOME#g"\
    -e "s#<ELMERICELGGElibs>#$ELMERICELGGElibs#g"\
    -e "s#<ID>#$ii#g"\
    -e "s#<nbnode>#$nbnode#g"\
    -e "s#<nbpart>#$nbpart#g"\
    -e "s#<EXCHANGE_DIR>#$EXCHANGE_DIR#g"\
    -e "s#<ELMER_WORKDIR>#$ELMER_WORKDIRrel#g"\
    -e "s#<NEMO_WORKDIRpro>#$NEMO_WORKDIRpro#g" > RUN$ii.slurm

done

####################################################
echo "PREPARING NEMO SECOND RELAXATION"
echo "GOING TO $NEMO_WORKDIRrel2 DIRECTORY"
cd $NEMO_WORKDIRrel2

NEMO_START_FROM_RST=0       #initially set to 0
NEMO_RST_DIR=    #should link the last restart file from nemo init simulations

EXT_NEMO=rel2
SBATCH_ID=$EXT_NEMO

cp -vr $TEMPLATE_DIR/NEMO/nemo_AMUXL025_GNJ002_BM01mv/* .    #all nemo scripts, already executables

NEMO_RELAX2_FINAL_YEAR=$(expr $INIT_YEAR+$NEMO_RELAX2_NBYEARS | bc)

cat $TEMPLATE_DIR/NEMO/nemo_AMUXL025_GNJ002_BM01mv/run_nemo.sh | sed -e "s#<SBATCH_ID>#$SBATCH_ID#g"\
  -e "s#<INIT_YEAR>#$INIT_YEAR#g"\
  -e "s#<INIT_MONTH>#$INIT_MONTH#g"\
  -e "s#<FINAL_YEAR>#$NEMO_RELAX2_FINAL_YEAR#g"\
  -e "s#<NEMO_WORKDIR>#$NEMO_WORKDIRrel2#g"\
  -e "s#<NEMO_START_FROM_RST>#$NEMO_START_FROM_RST#g"\
  -e "s#<NEMO_RST_DIR>#$NEMO_RST_DIR#g"\
  -e "s#<EXCHANGE_DIR>#$EXCHANGE_DIR#g"\
  -e "s#<BASE_DIR>#$BASE_DIR#g"\
  -e "s#<EXT_NEMO>#$EXT_NEMO#g"\
  -e "s#<ELMER_WORKDIRpro>#$ELMER_WORKDIRpro#g" > run_nemo.sh

echo "############################################################"
echo ""
echo "PREPARING PROGNOSTIC STEP"
echo ""
echo "############################################################"

echo "PREPARING ELMERICE PROGNOSTIC"
echo "GOING TO $ELMER_WORKDIRpro DIRECTORY"
cd $ELMER_WORKDIRpro

echo "SYNCING MESH USED AT THE INVERSION STEP"

rsync -av --exclude=*vtu --exclude=*result* $ELMER_PATH_INVERSE/RUNS_${FAMILY}/run_INIT_OPTIM/mesh_$nbpart .

echo "PREPARING RUNii"
echo "GOING TO $ELMER_WORKDIRpro DIRECTORY"
cd $ELMER_WORKDIRpro

SBATCH_ID=pro

imin=1
imax=$(expr 12*$COUPLING_NBYEARS | bc)

for ((ii=$imin ; ii<=$imax ; ii++))
do

  echo "Preparing ElmerIce $ii st simulation"

  ln -sf $SMB_ANOM_DATA
  ln -sf $EXCHANGE_DIR/MELT_RATES/$MELT_RATES_NEMO_XY
  SMB_ANOM_DATAsif=`echo $SMB_ANOM_DATA | rev | cut -d"/" -f1 | rev`
  MELT_RATES_NEMO_XYsif=`echo $MELT_RATES_NEMO_XY | rev | cut -d"/" -f1 | rev`

  cat $TEMPLATE_DIR/SIF/RUNii.sif | sed -e "s#<ID>#$ii#g"\
    -e "s#<ID-1>#$((ii-1))#g"\
    -e "s#<RHOI_SI>#$rhoi_SI#g"\
    -e "s#<RHOW_SI>#$rhow_SI#g"\
    -e "s#<GRAVITY_SI>#$gravity_SI#g"\
    -e "s#<YEARTOSEC>#$yeartosec#g"\
    -e "s#<ZSL>#$zsl#g"\
    -e "s#<N>#$n#g"\
    -e "s#<nbpart>#$nbpart#g"\
    -e "s#<ELMER_DT>#$ELMER_CPL_DT#g"\
    -e "s#<SMB_ANOM_DATA>#$SMB_ANOM_DATAsif#g"\
    -e "s#<MELT_RATES_NEMO_XY>#$MELT_RATES_NEMO_XYsif#g" > RUN$ii.sif

  cat $TEMPLATE_DIR/SLURM/RUNii.slurm | sed -e "s#<SBATCH_ID>#$SBATCH_ID#g"\
    -e "s#<ELMER_HOME>#$ELMER_HOME#g"\
    -e "s#<ELMERICELGGElibs>#$ELMERICELGGElibs#g"\
    -e "s#<ID>#$ii#g"\
    -e "s#<nbnode>#$nbnode#g"\
    -e "s#<nbpart>#$nbpart#g"\
    -e "s#<EXCHANGE_DIR>#$EXCHANGE_DIR#g"\
    -e "s#<ELMER_WORKDIR>#$ELMER_WORKDIRpro#g"\
    -e "s#<NEMO_WORKDIRpro>#$NEMO_WORKDIRpro#g" > RUN$ii.slurm

done

cp -v $TEMPLATE_DIR/SIF/LINEAR_SOLVER.txt .

####################################################
echo "PREPARING NEMO PROGNOSTIC"
echo "GOING TO $NEMO_WORKDIRpro DIRECTORY"
cd $NEMO_WORKDIRpro

NEMO_START_FROM_RST=1            #you start from a restart
NEMO_RST_DIR=$NEMO_WORKDIRrel2   #link to the fdirectory in which restart files are present
EXT_NEMO=pro

cp -vr $TEMPLATE_DIR/NEMO/nemo_AMUXL025_GNJ002_BM01mv/* .    #all nemo scripts, already executables

NEMO_CPL_FINAL_YEAR=$(expr $INIT_YEAR+$COUPLING_NBYEARS | bc)

cat $TEMPLATE_DIR/NEMO/nemo_AMUXL025_GNJ002_BM01mv/run_nemo.sh | sed -e "s#<SBATCH_ID>#$SBATCH_ID#g"\
  -e "s#<INIT_YEAR>#$INIT_YEAR#g"\
  -e "s#<INIT_MONTH>#$INIT_MONTH#g"\
  -e "s#<FINAL_YEAR>#$NEMO_CPL_FINAL_YEAR#g"\
  -e "s#<NEMO_WORKDIR>#$NEMO_WORKDIRpro#g"\
  -e "s#<NEMO_START_FROM_RST>#$NEMO_START_FROM_RST#g"\
  -e "s#<NEMO_RST_DIR>#$NEMO_RST_DIR#g"\
  -e "s#<EXCHANGE_DIR>#$EXCHANGE_DIR#g"\
  -e "s#<BASE_DIR>#$BASE_DIR#g"\
  -e "s#<EXT_NEMO>#$EXT_NEMO#g"\
  -e "s#<ELMER_WORKDIRpro>#$ELMER_WORKDIRpro#g" > run_nemo.sh

echo "############################################################"
echo ""
echo "PREPARING PYTHON STUFF FOR NEMO/ELMERICE EXCHANGE"
echo ""
echo "############################################################"

echo "GOING TO $EXCHANGE_DIR DIRECTORY"
cd $EXCHANGE_DIR

echo "PREPARING PYTHON ROUTINE TO MAKE NEMO MELT READABLE BY ELMERICE"

cat $TEMPLATE_DIR/TOOLS/extractXYMelt_fromNEMOsbc.py | sed -e "s#<SBC_FILE>#melt_rates_current.nc#g"\
  -e "s#<NEMO_MESH_MASK_FILE>#$NEMO_MESH_MASK_FILE#g"\
  -e "s#<ELMER_MASK_IN_NEMO>#$ELMER_MASK_IN_NEMO#g"\
  -e "s#<MELT_RATES_NEMO_XY>#$MELT_RATES_NEMO_XY#g"\
  -e "s#<YEARTOSEC>#$yeartosec#g"\
  -e "s#<RHOI_SI>#$rhoi_SI#g"\
  -e "s#<RHOW_SI>#$rhow_SI#g" > extractXYMelt_fromNEMOsbc.py

cat $TEMPLATE_DIR/SLURM/extractXYMelt_fromNEMOsbc.slurm | sed -e "s#<EXCHANGE_DIR>#$EXCHANGE_DIR#g"\
  -e "s#<MELT_RATES_NEMO_XY>#$MELT_RATES_NEMO_XY#g"\
  -e "s#<INIT_YEAR>#$INIT_YEAR#g"\
  -e "s#<INIT_MONTH>#$INIT_MONTH#g"\
  -e "s#<ELMER_WORKDIRpro>#$ELMER_WORKDIRpro#g"\
  -e "s#<MONTHS_TO_SAVE_RELAX1>#$MONTHS_TO_SAVE_RELAX1#g"\
  -e "s#<MONTHS_TO_SAVE_RELAX2>#$MONTHS_TO_SAVE_RELAX2#g" > extractXYMelt_fromNEMOsbc.slurm

echo "PREPARING PYTHON ROUTINE TO CREATE THE ELMERICE MASK WITHIN THE NEMO DOMAIN"

cat $TEMPLATE_DIR/TOOLS/makeElmerMaskInNemo.py | sed -e "s#<YEARTOSEC>#$yeartosec#g"\
  -e "s#<RHOI_SI>#$rhoi_SI#g"\
  -e "s#<RHOW_SI>#$rhow_SI#g"\
  -e "s#<NEMO_MESH_MASK_FILE>#$NEMO_MESH_MASK_FILE#g"\
  -e "s#<ELMER_CONTOUR_FILE>#$ELMER_CONTOUR_FILE#g"\
  -e "s#<ELMER_MASK_IN_NEMO>#$ELMER_MASK_IN_NEMO#g" > makeElmerMaskInNemo.py

cp -v $TEMPLATE_DIR/SLURM/makeElmerMaskInNemo.slurm .

echo "PREPARING PYTHON ROUTINE TO CHECK IF TOTAL MELT IS SIMILAR BETWEEN ELMERICE AND NEMO"

cat $TEMPLATE_DIR/TOOLS/checkBothMelts.py | sed -e "s#<YEARTOSEC>#$yeartosec#g"\
  -e "s#<RHOI_SI>#$rhoi_SI#g"\
  -e "s#<RHOW_SI>#$rhow_SI#g"\
  -e "s#<ELMER_MASK_IN_NEMO>#$ELMER_MASK_IN_NEMO#g"\
  -e "s#<MELT_RATES_NEMO_XY>#$MELT_RATES_NEMO_XY#g"\
  -e "s#<NBPART>#$nbpart#g" > checkBothMelts.py

cp -v $TEMPLATE_DIR/SLURM/checkBothMelts.slurm .

echo "PREPARING PYTHON ROUTINE TO CREATE ISF_DRAFT_METER NC FILE TO BE FURTHER USERD BY NEMO"

cat $TEMPLATE_DIR/TOOLS/exportVTUtoNEMOdraft.py | sed -e "s#<RHOI_SI>#$rhoi_SI#g"\
  -e "s#<RHOW_SI>#$rhow_SI#g"\
  -e "s#<NBPART>#$nbpart#g"\
  -e "s#<ELMER_MASK_IN_NEMO>#$ELMER_MASK_IN_NEMO#g"\
  -e "s#<BATHY_METER_FICH>#$BATHY_METER_FICH#g" > exportVTUtoNEMOdraft.py

cp -v $TEMPLATE_DIR/SLURM/exportVTUtoNEMOdraft.slurm .

echo "############################################################"
echo ""
echo "LAUNCHING RELAXATION STEP"
echo ""
echo "############################################################"

if [ $stageblock -eq 0 ]; then

  echo "DOING StageBlock00"
  echo "GOING TO $ELMER_WORKDIRrel DIRECTORY"
  cd $ELMER_WORKDIRrel

  sbatch ImportResultsInversion.slurm ${FAMILY} $EXCHANGE_DIR

fi

####################################################

if [ $stageblock -eq 1 ]; then

  echo "DOING StageBlock01"
  echo "LAUNCHING NEMO 1ST RELAXATION"
  echo "GOING TO $NEMO_WORKDIRrel1 DIRECTORY"
  cd $NEMO_WORKDIRrel1

  EXT_NEMO=rel1

  sbatch run_nemo.sh $EXT_NEMO

fi

####################################################

if [ $stageblock -eq 2 ]; then

  echo "DOING StageBlock02"
  echo "LAUNCHING ELMERICE RELAXATION"
  echo "GOING TO $ELMER_WORKDIRrel DIRECTORY"
  cd $ELMER_WORKDIRrel

  EXT_ELMER=rel

  # ElmerIce previous last step
  laststep=0

  imin=1
  imax=$ELMER_RELAX_NBYEARS

  #Check what ElmerIce step relaxation was previously done
  for ((ii=$imin ; ii<=$imax ; ii++))
  do
    tmpfile=$ELMER_WORKDIRrel/mesh_$nbpart/RUN${ii}_0002.pvtu
    if [ ! -f "$tmpfile" ]; then
      laststep=$((ii-1))
      break
    fi
  done

  echo "Last step previously done is $laststep"
  echo "Doing the rest if any"

  if [ $laststep -lt $imax ]; then

    #Redo the link to NEMOrel1 output files
    rootfile=${CONFIG}-${CASE}_1m_fwfisf_

    DIR=$NEMO_WORKDIRrel1/output_sbc
    imax=0
    for file in ${DIR}/*nc
    do
      number=`echo $file | rev | cut -d"@" -f1 | rev | cut -d"." -f1`
      if [ $number -gt $imax ]; then
        imax=$number
      fi
    done
    imin=$(($imax-${MONTHS_TO_SAVE_RELAX1}+1))

    ln -sf $EXCHANGE_DIR/MELT_RATES/${rootfile}@steps${imin}to${imax}_rel1.nc $EXCHANGE_DIR/MELT_RATES/melt_rates_current.nc
    ln -sf $EXCHANGE_DIR/MELT_RATES/melt_rates_nemo_xy_rel1 $EXCHANGE_DIR/MELT_RATES/melt_rates_nemo_xy

    echo "Starting relaxation at step $(($laststep+1))"
    sbatch RUN$(($laststep+1)).slurm $(($laststep+1))

  fi

fi

####################################################

if [ $stageblock -eq 3 ]; then

  echo "DOING StageBlock03"
  echo "LAUNCHING NEMO 2ND RELAXATION"
  echo "GOING TO $NEMO_WORKDIRrel2 DIRECTORY"
  cd $NEMO_WORKDIRrel2

  EXT_NEMO=rel2

  #Remove old rel2.nc files, they should be redone by run_nemo.sh
  rm $EXCHANGE_DIR/MELT_RATES/*rel2.nc

  sbatch run_nemo.sh $EXT_NEMO

fi


echo "############################################################"
echo ""
echo "LAUNCHING PROGNOSTIC STEP"
echo ""
echo "############################################################"

if [ $stageblock -eq 4 ]; then

  echo "DOING StageBlock04"
  echo "LAUNCHING 1-ELMERICE, 2-NEMO"

  EXT_NEMO=pro
  EXT_ELMER=pro
  rootfile=${CONFIG}-${CASE}_1m_fwfisf_

  #imin calculation
  iminpro=1

  echo "Check if already done"
  echo "Checking existence of prod_nemo.db file"
  file1=$NEMO_WORKDIRpro/prod_nemo.db

  if [ -f "$file1" ];then
    #check the last iteration
    iminpro=$(tail -n 1 $file1 | cut -d" " -f1)
    if [ $iminpro -eq 1 ]; then
      echo "Coupling job not started yet"
      echo "Will thus start at step 1"
      rm $file1
    else
      echo "Coupling already started"
      echo "Will restart at step $iminpro"
    fi
  fi

  if [ $iminpro -eq 1 ]; then
    
    # re-linking last result file from ELMERICE relaxation to prognostic directory, so ElmerIce prognostic run can restart from it
    iminrst=0
    imaxrst=$(expr $nbpart-1 | bc)
     for ((ii=$iminrst ; ii<=$imaxrst ; ii++))
    do
      ln -sf $ELMER_WORKDIRrel/mesh_$nbpart/RUN${ELMER_RELAX_NBYEARS}_.result.$ii $ELMER_WORKDIRpro/mesh_$nbpart/RUN0_.result.$ii
    done
      
    # relinking the rel2 file in Exchange
    # recompute the last fwfisf file in NEMOrel2
    imax=0
    DIR=$NEMO_WORKDIRrel2/output_sbc
    for file in ${DIR}/*nc
    do
      number=`echo $file | rev | cut -d"@" -f1 | rev | cut -d"." -f1`
      if [ $number -gt $imax ]; then
        imax=$number
      fi
    done
    imin=$(($imax-${MONTHS_TO_SAVE_RELAX2}+1))

    # The file should have been created at StageBlock03
    # Doing only the linking of melt files for ElmerIce
    ln -sf $EXCHANGE_DIR/MELT_RATES/${rootfile}@steps${imin}to${imax}_rel2.nc $EXCHANGE_DIR/MELT_RATES/melt_rates_current.nc
    ln -sf $EXCHANGE_DIR/MELT_RATES/melt_rates_nemo_xy_rel2 $EXCHANGE_DIR/MELT_RATES/melt_rates_nemo_xy
    
  else

    # relinking the pro file in Exchange 
    ln -sf $EXCHANGE_DIR/MELT_RATES/${rootfile}@steps$(($iminpro-1))_pro.nc $EXCHANGE_DIR/MELT_RATES/melt_rates_current.nc
    ln -sf $EXCHANGE_DIR/MELT_RATES/melt_rates_nemo_xy_pro$(($iminpro-1)) $EXCHANGE_DIR/MELT_RATES/melt_rates_nemo_xy

  fi

  # last coupling step for coupling simulations
  imaxpro=$(expr 12*$COUPLING_NBYEARS | bc)

  if [ $iminpro -lt $imaxpro ]; then

    echo "starting coupled prognostic at step $iminpro"
    echo "GOING TO $ELMER_WORKDIRpro DIRECTORY"
    cd $ELMER_WORKDIRpro
    sbatch RUN${iminpro}.slurm ${iminpro} $EXT_NEMO

  else

    echo "Coupled simulation already done, check your last time"

  fi

fi
