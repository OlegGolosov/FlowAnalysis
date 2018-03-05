#!/bin/bash

#Script to convert trees
user=$(whoami)
currentDir=/lustre/nyx/cbm/users/ogolosov

jobscript=${currentDir}/FlowAnalysis/bash/Slurm/HADES/getCorr.sh
inputdir=${currentDir}/HADES_conv_weight
outputdir=${currentDir}/HADES_corr_weight
env=${currentDir}/FlowAnalysis/bash/Slurm/HADES/env.sh
source=${currentDir}/FlowAnalysis
logdir=$outputdir/log
filelist=$outputdir/filelist.txt

if [ ! -d $outputdir ]
then
   echo "===> CREATE OUTPUTDIR : $outputdir"
   mkdir -p $outputdir
else
   echo "===> USE OUTPUTDIR : $outputdir"
fi


if [ "$#" -eq "0" ] || [ "$1" -eq "0" ]
then
    echo "All runs from $filelist will be submitted"
    maxNumberOfRuns=100000
else
    maxNumberOfRuns=${1}
fi


if [ ! -d $logdir ]
then
   echo "===> CREATE LOGDIR : $logdir"
   mkdir -p $logdir
else
   echo "===> USE LOGDIR : $logdir"
fi

cd $inputdir
ls *_.root>$filelist
echo "$filelist created"

numberOfRuns=0
while read F  ; do
    file=${F%.root}
    log_err=${logdir}/${file}_$2.err
    log_out=${logdir}/${file}_$2.out
    sbatch --error=${log_err} --output=${log_out} ${jobscript} ${env} ${source} ${outputdir}/${file}_$2 ${inputdir}/${file}
    ((numberOfRuns+=1))
    if [ ${numberOfRuns} -eq ${maxNumberOfRuns} ]
    then
        echo " "
        echo "${numberOfRuns} jobs have been submitted"
	break
    fi
done <$filelist

cd -
