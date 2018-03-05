#!/bin/bash

LOCATION=/mnt/pool/rhic/2/ovgol/NA49_conv

rm -f convert.o*
rm -f convert.e*

cd $LOCATION
#rm -f QA/NA49_$1_QA.root
hadd -f QA/NA49_$1_QA.root *_$1_QA.root
#rm -f *_$1_QA.root
