
#!/bin/bash

NFILES=8

for (( i=1; i<=$NFILES; i++ ));
do
    echo $i
    qsub -v NFILE=$i,SE=3 RapScan.sh
    sleep 1;
done
