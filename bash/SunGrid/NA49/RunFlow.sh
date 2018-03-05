#!/bin/bash

SOURCE=/home/basov/ovgol/FlowAnalysis
INPUT=/mnt/pool/rhic/2/ovgol/NA49_flow
OUTPUT=/mnt/pool/rhic/2/ovgol/NA49_flow
FILENAME=NA49_$1

rm -f getcorr.o*
rm -f getcorr.e*

cd $INPUT

rm -f "$OUTPUT/$FILENAME"_corr.root
hadd -f "$OUTPUT/$FILENAME"_corr.root *$1_corr.root

cd $SOURCE
root.exe -b -l -q 'RunFlowAnalysisNA49.C ("flow","'$OUTPUT/$FILENAME'")'

cd $INPUT
rm -f *_Q_0.root
rm -f *_Q_1.root
rm -f *_Q_2.root
rm -f *_sample.root
