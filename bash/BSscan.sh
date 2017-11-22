#!/bin/sh
#
#PBS -N BSscan
#PBS -q short
#PBS -l nodes=1:ppn=1,walltime=00:60:00

SRC=/home/basov/ovgol/FVC
OUT=/mnt/pool/rhic/2/ovgol

cd $SRC;
root.exe -b -l -q "BSscan.cxx+($NFILE)";
