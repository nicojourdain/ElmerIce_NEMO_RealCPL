#!/bin/bash

BASE_DIR=/scratch/shared/egige60/ElmerIce_NEMO_RealCPL
ELMER_HOME=/scratch/cnt0021/egi6035/SHARED/local                          #ElmerIce main code from CSC
ELMERICELGGElibs=$BASE_DIR/otherStuff/sourcesElmerIce/ElmerIceLGGE/lib    #ElmerIce LGGE from Renater + personal sources

#This is for testing purposes
#If you want to use it, your need to name your COUPLING directory such as ../COUPLING@blabla/..
#Use the @ as a delimiter
playground=@`echo $PWD | cut -d'@' -f2 | cut -d'/' -f1`
if [ $playground == '@' ]; then
  playground=
fi

echo "############################################################"
echo "INVERSION"
echo "The family argument is mandatory (i.e. BRONDEX)"
echo "The coupled / non coupled argument is optional"
echo "type anything to have coupled simulation (e.g. cpl)"
echo "in any case you need a FAMILY.IN in the PARAMETERS directory"
echo "#############################################################"

#test the number of arguments
echo $# arguments
if [ $# -lt 1 ];then
  echo "illegal number of parameters"
  echo "ABORT"
  exit
fi

#test if $1 file exists in '$BASE_DIR/PARAMETERS/'
fileIN=$BASE_DIR/PARAMETERS/$1
if [ ! -e $fileIN ]; then
  echo "the family file does not exist"
  echo "there should be a "$1".IN file"
  echo "in your ../../PARAMETERS directory"
  echo "ABORT"
  exit
fi
FAMILY=`echo $1 | cut -d"." -f1`

echo "############################################################"
echo ""
echo "Creating inverse simulation $FAMILY"
echo ""
echo "############################################################"

echo "CREATING DIRECTORIES"

mkdir -p RUNS_$FAMILY

cp -v createRun.sh RUNS_$FAMILY/createRun_$FAMILY.sh
cat $BASE_DIR/PARAMETERS/$1 | sed -e "s#<BASE_DIR>#$BASE_DIR#g" > RUNS_$FAMILY/$1

source RUNS_$FAMILY/$1

echo "GOING TO RUNS_$FAMILY DIRECTORY"
cd RUNS_$FAMILY

#Name of init job
nameINIT=INIT

echo "##############################"
echo "INITIALISATION"
echo "##############################"

DIRECTORY_INIT=run_INIT_OPTIM

#check if init has already been done
tocheck=$DIRECTORY_INIT/mesh_$nbpart/"$nameINIT"0001.pvtu
if [ ! -e "$tocheck" ]; then

  mkdir -p $DIRECTORY_INIT
  cd $DIRECTORY_INIT

  cp -rf $mesh mesh_$nbpart
  ElmerGrid 2 2 mesh_$nbpart -autoclean -order 1.0 0.1 0.01 -metis $nbpart 4

  #create init sif
  cat ../../TEMPLATE/SIF/INIT_OPTIM.sif | sed -e "s#<nbpart>#$nbpart#g"\
          -e "s#<nameINIT>#$nameINIT#g"\
          -e "s#<RHOI_SI>#$rhoi_SI#g"\
          -e "s#<RHOW_SI>#$rhow_SI#g"\
          -e "s#<GRAVITY_SI>#$gravity_SI#g"\
          -e "s#<YEARTOSEC>#$yeartosec#g"\
          -e "s#<LF_SI>#$lf_SI#g"\
          -e "s#<CW_SI>#$cw_SI#g"\
          -e "s#<KT>#$kt#g"\
          -e "s#<ZSL>#$zsl#g"\
          -e "s#<N>#$n#g"\
          -e "s#<PARMLT>#$parmlt#g"\
          -e "s#<GAMMA0>#$gamma0#g"\
          -e "s#<GLM>#$glm#g"\
          -e "s#<BASINNB>#$basinnb#g"\
          -e "s#<TIME_INIT>#$time_init#g"\
          -e "s#<EPSZ>#$epsz#g"\
          -e "s#<THICKNESS_DATA>#$THICKNESS_DATA#g"\
          -e "s#<BEDROCK_DATA>#$BEDROCK_DATA#g"\
          -e "s#<nameThicknessVar>#$nameThicknessVar#g"\
          -e "s#<nameBedVar>#$nameBedVar#g"\
          -e "s#<SMB_DATA>#$SMB_CLIM_DATA#g"\
          -e "s#<VISCOSITY_DATA>#$VISCOSITY_DATA#g"\
          -e "s#<VELOCITY_DATA_X>#$VELX4INI#g"\
          -e "s#<VELOCITY_DATA_Y>#$VELY4INI#g"\
          -e "s#<nameVxVar>#$nameVxVar#g"\
          -e "s#<nameVyVar>#$nameVyVar#g" > INIT_OPTIM.sif

  cat ../../TEMPLATE/SLURM/INIT.slurm | sed -e "s#<nbpart>#$nbpart#g"\
                 -e "s#<ELMER_HOME>#$ELMER_HOME#g"\
                 -e "s#<ELMERICELGGElibs>#$ELMERICELGGElibs#g"  > INIT.slurm

  #do initialisation
  echo INIT_OPTIM.sif > ELMERSOLVER_STARTINFO

  jobid0=$(sbatch --parsable INIT.slurm)
  echo "INIT" "id" $jobid0

  cd ..

else
  echo "INITIALISATION ALREADY DONE, GOING TO INVERSION STEP"
fi

if [ ! -e "$tocheck" ]; then
  echo "INITIALISATION NOT DONE YET"
  echo "RESTART createRun.sh"
  echo "ABORT"
  exit
fi

echo "##############################"
echo "INVERSIONS"
echo "##############################"

#inversions
for ((c=$rkcmin ; c<=$rkcmax ; c++))
do
  for ((dhdt=$rkdhdtmin ; dhdt<=$rkdhdtmax ; dhdt++))
  do

    name="Ga${Gastr}_Rcg${c}_Rdhdt${dhdt}"

    DIRECTORY_INV=run_OPTIM_$name
    #if [ -d "$DIRECTORY_INV" ]; then
    #  echo "DIRECTORY ALREADY EXISTING: " $DIRECTORY_INV
    #  echo "ABORT"
    #  exit
    #fi

    #create the directory
    mkdir -p $DIRECTORY_INV
    cd $DIRECTORY_INV

    #copy init solution + mesh
    rsync -a --exclude=*vtu ../$DIRECTORY_INIT/mesh_$nbpart .

    #lambda1 and lambda2
    echo $(awk -v n=$c 'NR == n' ../../INPUT/$nameLregSans) > tmp
    lbd1=$(cut -d" " -f1 tmp)
    lbd2=$(cut -d" " -f2 tmp)
    rm tmp

    #lambda3
    if [ $dhdt -eq 0 ]; then
      lbd3=0
    else
      echo $(awk -v n=$dhdt 'NR == n' ../../INPUT/$nameLregAvec) > tmp
      lbd3=$(cut -d" " -f1 tmp)
      rm tmp
    fi

    #create inversion sif
    cat ../../TEMPLATE/SIF/OPTIM.sif | sed -e "s#<nbpart>#$nbpart#g"\
          -e "s#<name>#$name#g"\
          -e "s#<lambda1>#$lbd1#g"\
          -e "s#<lambda2>#$lbd2#g"\
          -e "s#<lambda3>#$lbd3#g"\
          -e "s#<Ga>#$Ga#g"\
          -e "s#<nameINIT>#$nameINIT#g"\
          -e "s#<RHOI_SI>#$rhoi_SI#g"\
          -e "s#<RHOW_SI>#$rhow_SI#g"\
          -e "s#<GRAVITY_SI>#$gravity_SI#g"\
          -e "s#<YEARTOSEC>#$yeartosec#g"\
          -e "s#<LF_SI>#$lf_SI#g"\
          -e "s#<CW_SI>#$cw_SI#g"\
          -e "s#<KT>#$kt#g"\
          -e "s#<ZSL>#$zsl#g"\
          -e "s#<N>#$n#g"\
          -e "s#<PARMLT>#$parmlt#g"\
          -e "s#<GAMMA0>#$gamma0#g"\
          -e "s#<GLM>#$glm#g"\
          -e "s#<BASINNB>#$basinnb#g"\
          -e "s#<TIME_INIT>#$time_init#g"\
          -e "s#<EPSZ>#$epsz#g"\
          -e "s#<VELOCITY_DATA>#$VEL4INV#g" > OPTIM_$name.sif

    cp ../../TEMPLATE/SIF/LINEAR_SOLVER.txt .

    cat ../../TEMPLATE/SLURM/RUN.slurm | sed -e "s#<nbpart>#$nbpart#g"\
                 -e "s#<ELMER_HOME>#$ELMER_HOME#g"\
                 -e "s#<ELMERICELGGElibs>#$ELMERICELGGElibs#g"  > RUN.slurm

    #launch simulations
    echo OPTIM_$name.sif > ELMERSOLVER_STARTINFO

    jobid=$(sbatch --parsable RUN.slurm)
    echo "INVERSION: gamma=", $Ga ", rk c=" $c ", rk dhdt=" $dhdt 

    cd ..

  done
done
