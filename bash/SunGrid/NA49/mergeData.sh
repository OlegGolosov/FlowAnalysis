#!/bin/bash
LOCATION=/mnt/pool/rhic/2/ovgol/NA49_data

cd $LOCATION
for i in `seq 3003 3166`;
do
    hadd $i.root t49run$i.*.root
done
