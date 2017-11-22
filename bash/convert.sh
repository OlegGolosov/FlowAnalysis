#!/bin/sh
#
#SBATCH -J convert
#SBATCH -o %3/%4_out.log
#SBATCH -e %3/%4_err.log
#SBATCH --time=00:30:00

SRC=%1
IN=%2
OUT=%3
FILE=%4

cd $SRC;
root.exe -b -l -q 'RunConvertation.C++ ("'$IN'/'$FILE'","'$OUT'/'$FILE'")';
#root.exe -b -l -q 'Convert.cxx+("'$IN'/'$FILE'","'$OUT'"/'$FILE'")' > $OUT/Logs/"$FILE"_conv.log;
