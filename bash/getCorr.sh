#!/bin/sh
#
#PBS -N getcorr
#PBS -q short
#PBS -l nodes=1:ppn=1,walltime=06:00:00

cd $SRC;

root.exe -b -l -q 'Flow.cxx+ ("correlations","'$OUT'/'$FILE'","'$NUIN'/'$FILE'")';
#root.exe -b -l -q 'GetFlow.cxx+ ("correlations","'$IN'/'$FILE'","'$NUIN'/'$FILE'","'$UIN'/'$FILE'")' > $OUT/Logs/"$FILE"_corr.log;
