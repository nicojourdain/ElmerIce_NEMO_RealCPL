#!/bin/bash

#results not to be deleted
#the rest is deleted
good1=run_OPTIM_Ga100_Rcg3_Rdhdt0 # without flux divergence
good2=run_OPTIM_Ga100_Rcg3_Rdhdt6 # with    flux divergence

for ii in run_OPTIM* 
do
  if [ $ii == $good1 ] || [ $ii == $good2 ]; then
    echo 'Not deleting the $ii results'
  else
    rm -v $ii/mesh_24/*vtu
    rm -v $ii/mesh_24/*result*
  fi
done

#cp results in ../Results/...
folder=`pwd | rev | cut -d'/' -f1 | rev`
mkdir -p ../../Results
mkdir -p ../../Results/$folder

mkdir -p ../../Results/$folder/run_INIT_OPTIM
mkdir -p ../../Results/$folder/run_INIT_OPTIM/mesh_24

cp -v run_INIT_OPTIM/mesh_24/mesh.* ../../Results/$folder/run_INIT_OPTIM/mesh_24/
cp -rv run_INIT_OPTIM/mesh_24/partitioning* ../../Results/$folder/run_INIT_OPTIM/mesh_24/

for name in $good1 $good2
do

  mkdir -p ../../Results/$folder/$name
  mkdir -p ../../Results/$folder/$name/mesh_24

  cp -v $name/mesh_24/*7.vtu ../../Results/$folder/$name/mesh_24/
  cp -v $name/mesh_24/*7.pvtu ../../Results/$folder/$name/mesh_24/

done

