#!/bin/bash

IN=/mnt/pool/rhic/4/parfenovpeter/NA49/Data/RootFiles
OUT=/mnt/pool/rhic/2/ovgol/NA49_data
#DATASET=01d
#FIRSTRUN=3133
#LASTRUN=3166
DATASET=02c
FIRSTRUN=3003
LASTRUN=3060

hadd -f $OUT/$DATASET/$DATASET.root $IN/$DATASET/t49run*tree.root

for (( nRun=$FIRSTRUN; nRun<=$LASTRUN; nRun++ ));
do
    hadd -f $OUT/$DATASET/$nRun.root $IN/$DATASET/t49run$nRun.*tree.root
done
