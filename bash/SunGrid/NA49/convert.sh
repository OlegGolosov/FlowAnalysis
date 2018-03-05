#!/bin/sh
#
#PBS -N convert
#PBS -q short
#PBS -l nodes=1:ppn=1,walltime=06:00:00

cd $SRC;
root.exe -b -l -q 'RunConvertationNA49.C ("'$IN'/'$FILE'","'$OUT'/'$FILE'")';
#root.exe -b -l -q 'Convert.cxx+("'$IN'/'$FILE'","'$OUT'"/'$FILE'")' > $OUT/Logs/"$FILE"_conv.log;
