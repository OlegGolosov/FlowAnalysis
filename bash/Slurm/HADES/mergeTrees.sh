#!/bin/bash

IN=/lustre/nyx/cbm/users/ogolosov/HADES_conv_protons_
OUT=/lustre/nyx/cbm/users/ogolosov/HADES_conv_protons
FIRSTRUN=12108160806
LASTRUN=12108235650

mkdir $OUT

for (( nRun=$FIRSTRUN; nRun<=$LASTRUN; nRun++ ));
do
file=${IN}/${nRun}_1_.root
	if [ -e "$file" ]
	then
    		hadd -f ${OUT}/${nRun}_.root $IN/${nRun}*_.root
    		hadd -f ${OUT}/${nRun}_QA.root $IN/${nRun}*_QA.root 
	fi
done
