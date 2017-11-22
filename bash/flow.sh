#!/bin/bash

SOURCE=/home/basov/ovgol/FVC
INPUT=/mnt/pool/rhic/2/ovgol/NA49_flow
OUTPUT=/mnt/pool/rhic/2/ovgol/NA49_flow
FILENAME=NA49

cd $INPUT
hadd -f "$OUTPUT/$FILENAME"_corr.root *_full_mh_corr.root

cd $SOURCE
root.exe -b -l -q 'Flow.cxx+ ("flow","'$OUTPUT/$FILENAME'")'
rm -f getcorr.o*
rm -f getcorr.e*

cd $INPUT
rm -f *_Q_0.root
rm -f *_Q_1.root
rm -f *_Q_2.root
rm -f *_sample.root
