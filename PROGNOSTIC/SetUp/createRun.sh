#!/bin/bash

export LANG=C

#BASE_DIR=${PWD}/../../    #need to be similar to BASE_DIR in BASE_DIR/PARAMETERS/$1.IN
BASE_DIR=/scratch/shared/egige60/ElmerIce_NEMO_RealCPL

#This is for testing purposes
#If you want to use it, your need to name your PROGNOSTIC directory such as ../PROGNOSTIC@blabla/..
#Use the @ as a delimiter
playground=@`echo $PWD | cut -d'@' -f2 | cut -d'/' -f1`
if [ $playground == '@' ]; then
  playground=
fi

echo "############################################################"
echo "PROGNOSTIC"
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
# 1st block : Relaxation (Weertman linear)
# Either 2nd block : Prognostic (Weertman non linear)
# Or     3rd block : Prognostic (Schoof nonlinear)
stageblock=3

RELAX_NBYEARS=15    #multiple of 5, the duration of every Elmer run

PROGNOS_NBYEARS=5   #multiple of 5, the duration of every Elmer run

ELMER_PATH_INVERSE=$BASE_DIR/INVERSION$playground/Results                               #the inverse result that is listed in $1.IN

ELMER_WORKDIRrel=$BASE_DIR/PROGNOSTIC$playground/SetUp/RUNS_$FAMILY/RELAXATION              #ElmerIce working directory for relaxation step
ELMER_WORKDIRproWNL=$BASE_DIR/PROGNOSTIC$playground/SetUp/RUNS_$FAMILY/PROGNOSTIC_WNL       #ElmerIce working directory for prognostic step and non linear weertman law
ELMER_WORKDIRproSCHOOF=$BASE_DIR/PROGNOSTIC$playground/SetUp/RUNS_$FAMILY/PROGNOSTIC_SCHOOF #ElmerIce working directory for prognostic step and non linear schoof law

ELMER_HOME=/scratch/cnt0021/egi6035/SHARED/local                          #ElmerIce main code from CSC
ELMERICELGGElibs=$BASE_DIR/otherStuff/sourcesElmerIce/ElmerIceLGGE/lib    #ElmerIce LGGE from Renater + personal sources

TEMPLATE_DIR=$BASE_DIR/PROGNOSTIC$playground/SetUp/TEMPLATE

###################################################
## End of User's choices
###################################################

echo "############################################################"
echo ""
echo "Creating ElmerIce simulation $FAMILY"
echo "=> Be confident and it may work"
echo "=> Other stuff will follow..."
echo ""
echo "############################################################"

echo "CREATING DIRECTORIES"

#create directories
mkdir -vp RUNS_$FAMILY
mkdir -vp $ELMER_WORKDIRrel
mkdir -vp $ELMER_WORKDIRproWNL
mkdir -vp $ELMER_WORKDIRproSCHOOF

cp -v createRun.sh RUNS_$FAMILY/createRun_$FAMILY.sh   #copy the main script

cat $BASE_DIR/PARAMETERS/$1 | sed -e "s#<BASE_DIR>#$BASE_DIR#g" > RUNS_$FAMILY/$1
source RUNS_$FAMILY/$1

echo "############################################################"
echo ""
echo "PREPARING RELAXATION STEP"
echo ""
echo "############################################################"

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

echo "PREPARING VARIABLES FOR TYPE OF MELT"

if [ ${MELT_TYPE} == pdc_melt ]; then
  EXEC_MISMIP_MELT=Never
  EXEC_NC_MELT=Never
elif [ ${MELT_TYPE} == nc_melt ]; then
  EXEC_MISMIP_MELT=Never
  EXEC_NC_MELT=Always
elif [ ${MELT_TYPE} == mismip_melt ]; then
  EXEC_MISMIP_MELT=Always
  EXEC_NC_MELT=Never
fi

echo "PREPARING RUN0 TO INITIALISE RELAXATION"

ln -svf $THICKNESS_DATA
ln -svf $BEDROCK_DATA
ln -svf $SMB_CLIM_DATA
ln -svf $SMB_ANOM_DATA
ln -svf $BMB_DATA
ln -svf $ICESHELVES_MASK_DATA
ln -svf $BASINS
ln -svf $TF
ln -svf $DELTAT
ln -svf $MELT_RATES_NEMO_XY
ln -svf $MELT_RATES_NEMO_TOTAL
ln -svf $NC_MELT_FILE
THICKNESS_DATAsif=`echo $THICKNESS_DATA | rev | cut -d"/" -f1 | rev`
BEDROCK_DATAsif=`echo $BEDROCK_DATA | rev | cut -d"/" -f1 | rev`
SMB_CLIM_DATAsif=`echo $SMB_CLIM_DATA | rev | cut -d"/" -f1 | rev`
SMB_ANOM_DATAsif=`echo $SMB_ANOM_DATA | rev | cut -d"/" -f1 | rev`
BMB_DATAsif=`echo $BMB_DATA | rev | cut -d"/" -f1 | rev`
ICESHELVES_MASK_DATAsif=`echo $ICESHELVES_MASK_DATA | rev | cut -d"/" -f1 | rev`
BASINSsif=`echo $BASINS | rev | cut -d"/" -f1 | rev`
TFsif=`echo $TF | rev | cut -d"/" -f1 | rev`
DELTATsif=`echo $DELTAT | rev | cut -d"/" -f1 | rev`
MELT_RATES_NEMO_XYsif=`echo $MELT_RATES_NEMO_XY | rev | cut -d"/" -f1 | rev`
MELT_RATES_NEMO_TOTALsif=`echo $MELT_RATES_NEMO_TOTAL | rev | cut -d"/" -f1 | rev`
NC_MELT_FILEsif=`echo $NC_MELT_FILE | rev | cut -d"/" -f1 | rev`

cat $TEMPLATE_DIR/SIF/RUN0.sif | sed -e "s#<nbpart>#$nbpart#g"\
          -e "s#<nameINIT>#$nameINIT#g"\
          -e "s#<RHOI_SI>#$rhoi_SI#g"\
          -e "s#<RHOW_SI>#$rhow_SI#g"\
          -e "s#<GRAVITY_SI>#$gravity_SI#g"\
          -e "s#<YEARTOSEC>#$yeartosec#g"\
          -e "s#<LF_SI>#$lf_SI#g"\
          -e "s#<CW_SI>#$cw_SI#g"\
          -e "s#<EXEC_MISMIP_MELT>#$EXEC_MISMIP_MELT#g"\
          -e "s#<KT>#$kt#g"\
          -e "s#<ZSL>#$zsl#g"\
          -e "s#<N>#$n#g"\
          -e "s#<MISMIP_MELT_TYPE>#$MISMIP_MELT_TYPE#g"\
          -e "s#<GAMMA0>#$gamma0#g"\
          -e "s#<GLM>#$glm#g"\
          -e "s#<BASINNB>#$basinnb#g"\
          -e "s#<TIME_INIT>#$time_init#g"\
          -e "s#<EPSZ>#$epsz#g"\
          -e "s#<THICKNESS_DATA>#$THICKNESS_DATAsif#g"\
          -e "s#<BEDROCK_DATA>#$BEDROCK_DATAsif#g"\
          -e "s#<nameThicknessVar>#$nameThicknessVar#g"\
          -e "s#<nameBedVar>#$nameBedVar#g"\
          -e "s#<SMB_DATA>#$SMB_CLIM_DATAsif#g"\
          -e "s#<BMB_DATA>#$BMB_DATAsif#g"\
          -e "s#<CWL_DATA>#$CWL_DATA#g"\
          -e "s#<ETA_DATA>#$ETA_DATA#g"\
          -e "s#<ICESHELVES_MASK_DATA>#$ICESHELVES_MASK_DATAsif#g"\
          -e "s#<BASINS>#$BASINSsif#g"\
          -e "s#<TF>#$TFsif#g"\
          -e "s#<DELTAT>#$DELTATsif#g"\
          -e "s#<MELT_RATES_NEMO_XY>#$MELT_RATES_NEMO_XYsif#g"\
          -e "s#<MELT_RATES_NEMO_TOTAL>#$MELT_RATES_NEMO_TOTALsif#g"\
          -e "s#<BMB_NAME>#$BMB_NAME#g" > RUN0.sif

cat $TEMPLATE_DIR/SLURM/RUN0.slurm | sed -e "s#<nbpart>#$nbpart#g"\
  -e "s#<ELMER_HOME>#$ELMER_HOME#g"\
  -e "s#<ELMERICELGGElibs>#$ELMERICELGGElibs#g"\
  -e "s#<nbnode>#$nbnode#g" > RUN0.slurm

echo "PREPARING RUNii TO PERFORM RELAXATION"

cp -v $TEMPLATE_DIR/SIF/LINEAR_SOLVER.txt .

imin=1
imax=$(expr $RELAX_NBYEARS/5 | bc)

for ((ii=$imin ; ii<=$imax ; ii++))
do

  echo $ii

  cat $TEMPLATE_DIR/SIF/RUNii.sif | sed -e "s#<nbpart>#$nbpart#g"\
          -e "s#<ID>#$ii#g"\
          -e "s#<ID-1>#$(($ii-1))#g"\
          -e "s#<RHOI_SI>#$rhoi_SI#g"\
          -e "s#<RHOW_SI>#$rhow_SI#g"\
          -e "s#<GRAVITY_SI>#$gravity_SI#g"\
          -e "s#<YEARTOSEC>#$yeartosec#g"\
          -e "s#<LF_SI>#$lf_SI#g"\
          -e "s#<CW_SI>#$cw_SI#g"\
          -e "s#<EXEC_MISMIP_MELT>#$EXEC_MISMIP_MELT#g"\
          -e "s#<MELT_TYPE>#$MELT_TYPE#g"\
          -e "s#<EXEC_NC_MELT>#$EXEC_NC_MELT#g"\
          -e "s#<NC_MELT_FILE>#$NC_MELT_FILEsif#g"\
          -e "s#<KT>#$kt#g"\
          -e "s#<ZSL>#$zsl#g"\
          -e "s#<N>#$n#g"\
          -e "s#<MISMIP_MELT_TYPE>#$MISMIP_MELT_TYPE#g"\
          -e "s#<GAMMA0>#$gamma0#g"\
          -e "s#<GLM>#$glm#g"\
          -e "s#<BASINNB>#$basinnb#g"\
          -e "s#<TIME_INIT>#$time_init#g"\
          -e "s#<EPSZ>#$epsz#g"\
          -e "s#<TF>#$TFsif#g"\
          -e "s#<SMB_DATA>#$SMB_ANOM_DATAsif#g"\
          -e "s#<BMB_NAME>#$BMB_NAME#g" > RUN$ii.sif

  cat $TEMPLATE_DIR/SLURM/RUNii.slurm | sed -e "s#<nbpart>#$nbpart#g"\
          -e "s#<ID>#$ii#g"\
          -e "s#<ELMER_HOME>#$ELMER_HOME#g"\
          -e "s#<ELMERICELGGElibs>#$ELMERICELGGElibs#g"\
          -e "s#<nbpart>#$nbpart#g" > RUN$ii.slurm

done

echo "############################################################"
echo ""
echo "PREPARING PROGNOSTIC STEP"
echo "WEERTMAN NON LINEAR AND SCHOOF FRICTIONS LAWS"
echo ""
echo "############################################################"

echo "PREPARING ELMERICE PROGNOSTIC NON LINEAR WEERTMAN FRICTION LAW"
echo "GOING TO $ELMER_WORKDIRproWNL DIRECTORY"
cd $ELMER_WORKDIRproWNL

echo "SYNCING MESH USED AT THE INVERSION STEP"
rsync -av --exclude=*vtu --exclude=*result* $ELMER_PATH_INVERSE/RUNS_${FAMILY}/run_INIT_OPTIM/mesh_$nbpart .

echo "PREPARING RUNii TO PERFORM RELAXATION"

cp -v $TEMPLATE_DIR/SIF/LINEAR_SOLVER.txt .

ln -svf $SMB_ANOM_DATA
ln -svf $NC_MELT_FILE
SMB_ANOM_DATAsif=`echo $SMB_ANOM_DATA | rev | cut -d"/" -f1 | rev`
NC_MELT_FILEsif=`echo $NC_MELT_FILE | rev | cut -d"/" -f1 | rev`

cat $TEMPLATE_DIR/SIF/RUNPRO0_WNL.sif | sed -e "s#<nbpart>#$nbpart#g"\
        -e "s#<RHOI_SI>#$rhoi_SI#g"\
        -e "s#<RHOW_SI>#$rhow_SI#g"\
        -e "s#<GRAVITY_SI>#$gravity_SI#g"\
        -e "s#<YEARTOSEC>#$yeartosec#g"\
        -e "s#<LF_SI>#$lf_SI#g"\
        -e "s#<CW_SI>#$cw_SI#g"\
        -e "s#<KT>#$kt#g"\
        -e "s#<ZSL>#$zsl#g"\
        -e "s#<N>#$n#g"\
        -e "s#<GAMMA0>#$gamma0#g"\
        -e "s#<GLM>#$glm#g"\
        -e "s#<BASINNB>#$basinnb#g"\
        -e "s#<TIME_INIT>#$time_init#g"\
        -e "s#<EPSZ>#$epsz#g"\
        -e "s#<TF>#$TFsif#g" > RUNPRO0.sif

cat $TEMPLATE_DIR/SLURM/RUNPRO0_WNL.slurm | sed -e "s#<nbpart>#$nbpart#g"\
        -e "s#<ELMER_HOME>#$ELMER_HOME#g"\
        -e "s#<ELMERICELGGElibs>#$ELMERICELGGElibs#g" > RUNPRO0.slurm

imin=1
imax=$(expr $PROGNOS_NBYEARS/5 | bc)

for ((ii=$imin ; ii<=$imax ; ii++))
do

  echo $ii

  cat $TEMPLATE_DIR/SIF/RUNPROii_WNL.sif | sed -e "s#<nbpart>#$nbpart#g"\
          -e "s#<ID>#$ii#g"\
          -e "s#<ID-1>#$(($ii-1))#g"\
          -e "s#<RHOI_SI>#$rhoi_SI#g"\
          -e "s#<RHOW_SI>#$rhow_SI#g"\
          -e "s#<GRAVITY_SI>#$gravity_SI#g"\
          -e "s#<YEARTOSEC>#$yeartosec#g"\
          -e "s#<LF_SI>#$lf_SI#g"\
          -e "s#<CW_SI>#$cw_SI#g"\
          -e "s#<EXEC_MISMIP_MELT>#$EXEC_MISMIP_MELT#g"\
          -e "s#<MELT_TYPE>#$MELT_TYPE#g"\
          -e "s#<EXEC_NC_MELT>#$EXEC_NC_MELT#g"\
          -e "s#<NC_MELT_FILE>#$NC_MELT_FILEsif#g"\
          -e "s#<KT>#$kt#g"\
          -e "s#<ZSL>#$zsl#g"\
          -e "s#<N>#$n#g"\
          -e "s#<MISMIP_MELT_TYPE>#$MISMIP_MELT_TYPE#g"\
          -e "s#<GAMMA0>#$gamma0#g"\
          -e "s#<GLM>#$glm#g"\
          -e "s#<BASINNB>#$basinnb#g"\
          -e "s#<TIME_INIT>#$time_init#g"\
          -e "s#<EPSZ>#$epsz#g"\
          -e "s#<TF>#$TFsif#g"\
          -e "s#<SMB_DATA>#$SMB_ANOM_DATAsif#g"\
          -e "s#<BMB_NAME>#$BMB_NAME#g" > RUNPRO$ii.sif

  cat $TEMPLATE_DIR/SLURM/RUNPROii_WNL.slurm | sed -e "s#<nbpart>#$nbpart#g"\
          -e "s#<ID>#$ii#g"\
          -e "s#<ELMER_HOME>#$ELMER_HOME#g"\
          -e "s#<ELMERICELGGElibs>#$ELMERICELGGElibs#g" > RUNPRO$ii.slurm

done


echo "PREPARING ELMERICE PROGNOSTIC NON LINEAR SCHOOF FRICTION LAW"
echo "GOING TO $ELMER_WORKDIRproSCHOOF DIRECTORY"
cd $ELMER_WORKDIRproSCHOOF

echo "SYNCING MESH USED AT THE INVERSION STEP"
rsync -av --exclude=*vtu --exclude=*result* $ELMER_PATH_INVERSE/RUNS_${FAMILY}/run_INIT_OPTIM/mesh_$nbpart .

echo "PREPARING RUNii TO PERFORM RELAXATION"

echo "PREPARING PYTHON SCRIPTS TO COMPUTE CS FROM SCHOOF"

NAMEFILE_CSSCHOOF=csSchoof_Ga${Gastr}_Rcg${rkcg}_Rdhdt${rkdhdt}.dat

cat $TEMPLATE_DIR/TOOLS/Calcul_Cs_LoiSchoof.py | sed -e "s#<M>#$M#g"\
        -e "s#<CMAX>#$CMAX#g"\
        -e "s#<CLIM>#$CLIM#g"\
        -e "s#<ELMER_WORKDIRrel>#$ELMER_WORKDIRrel#g"\
        -e "s#<nbpart>#$nbpart#g"\
        -e "s#<ELMER_WORKDIRproSCHOOF>#$ELMER_WORKDIRproSCHOOF#g"\
        -e "s#<NAMEFILE_CSSCHOOF>#$NAMEFILE_CSSCHOOF#g" > Calcul_Cs_LoiSchoof.py

cp -v $TEMPLATE_DIR/SLURM/Calcul_Cs_LoiSchoof.slurm .


cp -v $TEMPLATE_DIR/SIF/LINEAR_SOLVER.txt .

ln -svf $SMB_ANOM_DATA
ln -svf $NC_MELT_FILE
SMB_ANOM_DATAsif=`echo $SMB_ANOM_DATA | rev | cut -d"/" -f1 | rev`
NC_MELT_FILEsif=`echo $NC_MELT_FILE | rev | cut -d"/" -f1 | rev`

cat $TEMPLATE_DIR/SIF/RUNPRO0_SCHOOF.sif | sed -e "s#<nbpart>#$nbpart#g"\
        -e "s#<RHOI_SI>#$rhoi_SI#g"\
        -e "s#<RHOW_SI>#$rhow_SI#g"\
        -e "s#<GRAVITY_SI>#$gravity_SI#g"\
        -e "s#<YEARTOSEC>#$yeartosec#g"\
        -e "s#<LF_SI>#$lf_SI#g"\
        -e "s#<CW_SI>#$cw_SI#g"\
        -e "s#<KT>#$kt#g"\
        -e "s#<ZSL>#$zsl#g"\
        -e "s#<N>#$n#g"\
        -e "s#<GAMMA0>#$gamma0#g"\
        -e "s#<GLM>#$glm#g"\
        -e "s#<BASINNB>#$basinnb#g"\
        -e "s#<TIME_INIT>#$time_init#g"\
        -e "s#<EPSZ>#$epsz#g"\
        -e "s#<NAMEFILE_CSSCHOOF>#$NAMEFILE_CSSCHOOF#g" > RUNPRO0.sif

cat $TEMPLATE_DIR/SLURM/RUNPRO0_SCHOOF.slurm | sed -e "s#<nbpart>#$nbpart#g"\
        -e "s#<ELMER_HOME>#$ELMER_HOME#g"\
        -e "s#<ELMERICELGGElibs>#$ELMERICELGGElibs#g" > RUNPRO0.slurm

imin=1
imax=$(expr $PROGNOS_NBYEARS/5 | bc)

for ((ii=$imin ; ii<=$imax ; ii++))
do

  echo $ii

  cat $TEMPLATE_DIR/SIF/RUNPROii_SCHOOF.sif | sed -e "s#<nbpart>#$nbpart#g"\
          -e "s#<ID>#$ii#g"\
          -e "s#<ID-1>#$(($ii-1))#g"\
          -e "s#<RHOI_SI>#$rhoi_SI#g"\
          -e "s#<RHOW_SI>#$rhow_SI#g"\
          -e "s#<GRAVITY_SI>#$gravity_SI#g"\
          -e "s#<YEARTOSEC>#$yeartosec#g"\
          -e "s#<LF_SI>#$lf_SI#g"\
          -e "s#<CW_SI>#$cw_SI#g"\
          -e "s#<EXEC_MISMIP_MELT>#$EXEC_MISMIP_MELT#g"\
          -e "s#<MELT_TYPE>#$MELT_TYPE#g"\
          -e "s#<EXEC_NC_MELT>#$EXEC_NC_MELT#g"\
          -e "s#<NC_MELT_FILE>#$NC_MELT_FILEsif#g"\
          -e "s#<KT>#$kt#g"\
          -e "s#<ZSL>#$zsl#g"\
          -e "s#<N>#$n#g"\
          -e "s#<MISMIP_MELT_TYPE>#$MISMIP_MELT_TYPE#g"\
          -e "s#<GAMMA0>#$gamma0#g"\
          -e "s#<GLM>#$glm#g"\
          -e "s#<BASINNB>#$basinnb#g"\
          -e "s#<TIME_INIT>#$time_init#g"\
          -e "s#<EPSZ>#$epsz#g"\
          -e "s#<TF>#$TFsif#g"\
          -e "s#<SMB_DATA>#$SMB_ANOM_DATAsif#g"\
          -e "s#<BMB_NAME>#$BMB_NAME#g" > RUNPRO$ii.sif

  cat $TEMPLATE_DIR/SLURM/RUNPROii_SCHOOF.slurm | sed -e "s#<nbpart>#$nbpart#g"\
          -e "s#<ID>#$ii#g"\
          -e "s#<ELMER_HOME>#$ELMER_HOME#g"\
          -e "s#<ELMERICELGGElibs>#$ELMERICELGGElibs#g" > RUNPRO$ii.slurm

done

echo "############################################################"
echo ""
echo "LAUNCHING RELAXATION STEP"
echo ""
echo "############################################################"

if [ $stageblock -eq 0 ]; then

  echo "DOING StageBlock00"
  echo "GOING TO $ELMER_WORKDIRrel DIRECTORY"
  cd $ELMER_WORKDIRrel

  sbatch ImportResultsInversion.slurm ${FAMILY}

fi

####################################################

if [ $stageblock -eq 1 ]; then

  echo "DOING StageBlock01"
  echo "LAUNCHING ELMERICE RELAXATION"
  echo "GOING TO $ELMER_WORKDIRrel DIRECTORY"
  cd $ELMER_WORKDIRrel

  # ElmerIce previous last step
  laststep=0

  imin=1
  imax=$(expr $RELAX_NBYEARS/5 | bc)

  for ((ii=$imin ; ii<=$imax ; ii++))
  do
    tmpfile=$ELMER_WORKDIRrel/mesh_$nbpart/RUN${ii}_0006.pvtu
    if [ ! -f "$tmpfile" ]; then
      laststep=$(($ii-1))
      break
    fi
  done

  echo "Last step previously done is $laststep"
  echo "Doing the rest if any"

  if [ $laststep -lt $imax ]; then

    echo "Starting relaxation at step $(($laststep+1))"
    sbatch RUN$(($laststep+1)).slurm $(($laststep+1))

  fi  

fi

echo "############################################################"
echo ""
echo "LAUNCHING PROGNOSTIC STEP"
echo "WITH NON LINEAR WEERTMAN FRICTION LAW"
echo ""
echo "############################################################"


if [ $stageblock -eq 2 ]; then

  echo "DOING StageBlock02"
  echo "LAUNCHING ELMERICE PROGNOSTIC WITH NON LINEAR WEERTMAN FRICTION LAW"
  echo "GOING TO $ELMER_WORKDIRproWNL DIRECTORY"
  cd $ELMER_WORKDIRproWNL

  # re-linking last result file from ELMERICE relaxation to prognostic directory, so ElmerIce prognostic run can restart from it
  iminrst=0
  imaxrst=$(expr $nbpart-1 | bc)

  #imaxprev=$(expr $RELAX_NBYEARS/5 | bc)
  imaxprev=1 # assuming that the linear relaxation over $RELAX_NBYEARS years has been done in one shot

  for ((ii=$iminrst ; ii<=$imaxrst ; ii++))
  do
    ln -svf $ELMER_WORKDIRrel/mesh_$nbpart/RUN${imaxprev}_.result.$ii $ELMER_WORKDIRproWNL/mesh_$nbpart/RUN0_.result.$ii
  done

  sbatch RUNPRO0.slurm 0

fi


echo "############################################################"
echo ""
echo "LAUNCHING PROGNOSTIC STEP"
echo "WITH NON LINEAR SCHOOF FRICTION LAW"
echo ""
echo "############################################################"


if [ $stageblock -eq 3 ]; then

  echo "DOING StageBlock03"
  echo "LAUNCHING ELMERICE PROGNOSTIC WITH NON LINEAR SCHOOF FRICTION LAW"
  echo "GOING TO $ELMER_WORKDIRproSCHOOF DIRECTORY"
  cd $ELMER_WORKDIRproSCHOOF

  # re-linking last result file from ELMERICE relaxation to prognostic directory, so ElmerIce prognostic run can restart from it
  iminrst=0
  imaxrst=$(expr $nbpart-1 | bc)

  #imaxprev=$(expr $RELAX_NBYEARS/5 | bc)
  imaxprev=1 # assuming that the linear relaxation over $RELAX_NBYEARS years has been done in one shot

  for ((ii=$iminrst ; ii<=$imaxrst ; ii++))
  do
    ln -svf $ELMER_WORKDIRrel/mesh_$nbpart/RUN${imaxprev}_.result.$ii $ELMER_WORKDIRproSCHOOF/mesh_$nbpart/RUN0_.result.$ii
  done

  sbatch Calcul_Cs_LoiSchoof.slurm

fi


