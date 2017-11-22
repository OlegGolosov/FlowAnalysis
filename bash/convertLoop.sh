#!/bin/sh

SOURCE=/u/ogolosov/FlowAnalysis/bash
INPUT=/lustre/nyx/cbm/users/ogolosov/NA49_data
OUTPUT=/lustre/nyx/cbm/users/ogolosov/NA49_conv

cd $INPUT
ls *.root>$OUTPUT/filelist.txt

cd $SOURCE

while read F  ; do
    echo "${F%.root}";
    #sbatch -v FILE=${F%.root},SRC=$SOURCE,IN=$INPUT,OUT=$OUTPUT convert.sh
    #sbatch convert.sh $SOURCE/.. $INPUT $OUTPUT ${F%.root} 
    sbatch convert.sh $SOURCE/.. $INPUT $OUTPUT 3154 
done <$OUTPUT/filelist.txt
