#!/bin/bash

currentDir=/lustre/nyx/cbm/users/ogolosov
jobscript=$currentDir/FlowAnalysis/bash/Slurm/HADES/getFlow.sh
inputdir=$currentDir/HADES_corr_weight
outputdir=$currentDir/HADES_flow_weight
source=$currentDir/FlowAnalysis
env=$currentDir/env.sh
logdir=$outputdir/log

if [ $# -eq 1 ]
then
   ext=$1
else
   ext="_"
fi

echo ext = $ext

filename=HADES_$ext

log_err=${logdir}/$filename.err
log_out=${logdir}/$filename.out

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

rm -f "$outputdir/$filename"_corr.root
cd $inputdir
hadd -f "$outputdir/$filename"_corr.root *${ext}_corr.root
#rm -f *_Q_0.root
#rm -f *_Q_1.root
#rm -f *_Q_2.root
#rm -f *_sample.root

cd -

sbatch --error=${log_err} --output=${log_out} ${jobscript} ${env} ${source} $outputdir/$filename

