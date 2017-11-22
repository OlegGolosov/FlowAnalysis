#!/bin/sh
#
#PBS -N convert
#PBS -q short
#PBS -l nodes=1:ppn=1,walltime=00:06:00

SRC=/home/basov/ovgol/FVC
IN=/mnt/pool/rhic/2/ovgol/NA49_raw
OUT=/mnt/pool/rhic/2/ovgol/NA49_y_16

cd $SRC;

#root.exe -b -l -q 'Convert.cxx+ ("'$IN'/'$FILE'","'$OUT'"/'$FILE'")';
root.exe -b -l -q 'Convert.cxx'
