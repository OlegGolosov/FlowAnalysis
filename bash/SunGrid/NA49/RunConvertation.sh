#!/bin/sh

SOURCE=/home/basov/ovgol/FlowAnalysis
BASH=$SOURCE/bash/SunGrid
INPUT=/mnt/pool/rhic/2/ovgol/NA49_data
OUTPUT=/mnt/pool/rhic/2/ovgol/NA49_conv

cd $INPUT
ls *.root>$OUTPUT/filelist.txt

cd $BASH
while read F  ; do
    echo "${F%.root}";
    qsub -v FILE=${F%.root},SRC=$SOURCE,IN=$INPUT,OUT=$OUTPUT convert.sh
done <$OUTPUT/filelist.txt
