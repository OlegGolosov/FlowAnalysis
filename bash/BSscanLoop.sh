#!/bin/bash

NFILES=6

for (( i=1; i<=$NFILES; i++ ));
do
    echo $i
    qsub -v NFILE=$i BSscan.sh
    #sleep 1;
done
