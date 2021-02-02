#!/bin/bash

for exps in 2Dfields/EXP10*nc
do
  echo $exps
  namedir=$(echo $exps | cut -f 1 -d '.')
  rm -fr $namedir
  echo $namedir
  python plotMISOMIPIceData.py $exps $namedir
done  
