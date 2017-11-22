#!/bin/bash

#1 nfiles
SOURCE=/home/basov/ovgol/FVC
OUTPUT=/home/basov/ovgol/Generated

for (( i=1; i<=2; i++ )); do
    qsub -q short -v NUM=$i,SRC=$SOURCE,OUT=$OUTPUT generate.sh;
    sleep 1;
done
