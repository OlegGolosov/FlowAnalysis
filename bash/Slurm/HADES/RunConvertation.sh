#!/bin/bash

#Script to convert trees
user=$(whoami)   
currentDir=/lustre/nyx/cbm/users/ogolosov

jobscript=${currentDir}/FlowAnalysis/bash/Slurm/HADES/convert.sh
inputdir=/lustre/nyx/cbm/users/ogolosov/HADES_data/treeMaker/output/Feb_19_23_17/AuAu_1_23AGev_gen9_108.list
outputdir=/lustre/nyx/cbm/users/ogolosov/HADES_conv_protons
source=${currentDir}/FlowAnalysis
env=${currentDir}/env.sh
filelist=$outputdir/filelist.txt
logdir=$outputdir/log

if [ ! -d $outputdir ]
then
   echo "===> CREATE OUTPUTDIR : $outputdir"
   mkdir -p $outputdir
else
   echo "===> USE OUTPUTDIR : $outputdir"
fi


if [ ! -d $logdir ]
then
   echo "===> CREATE LOGDIR : $logdir"
   mkdir -p $logdir
else
   echo "===> USE LOGDIR : $logdir"
fi


if [ "$#" -eq "0" ] || [ "$1" -eq "0" ]
then
    echo "All runs from $filelist will be submitted"
    maxNumberOfRuns=1000000
else
    maxNumberOfRuns=${1}
fi


ls $inputdir/*.root>$filelist
echo "$filelist created"

numberOfRuns=0 

for file in $(cat $filelist | sort | sed -r 's/.root//')
do
    runnumber=$(echo $file | sort| sed -r 's/.*\/tree_//'| sed -r 's/.root//')
    log_err=${logdir}/${runnumber}_$2.err
    log_out=${logdir}/${runnumber}_$2.out
    sbatch --error=${log_err} --output=${log_out} ${jobscript} ${env} ${source} ${file} ${outputdir}/${runnumber}_${2}
    ((numberOfRuns+=1))
    if [ ${numberOfRuns} -eq ${maxNumberOfRuns} ]
    then
        echo " "
        echo "${numberOfRuns} jobs have been submitted"
        exit 0
    fi
done
