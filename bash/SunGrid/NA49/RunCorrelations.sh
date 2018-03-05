#!/bin/sh

SOURCE=/home/basov/ovgol/FlowAnalysis
INPUT_NONUNI=/mnt/pool/rhic/2/ovgol/NA49_conv
INPUT_UNI=/mnt/pool/rhic/2/ovgol/NA49_conv
OUTPUT=/mnt/pool/rhic/2/ovgol/NA49_flow

cd $INPUT_NONUNI
ls *$1.root>filelist.txt

cd $SOURCE/bash/SunGrid
while read F  ; do
    echo "${F%.root}";
        qsub -v FILE=${F%.root},SRC=$SOURCE,NUIN=$INPUT_NONUNI,UIN=$INPUT_UNI,OUT=$OUTPUT getCorr.sh
done <$INPUT_NONUNI/filelist.txt
