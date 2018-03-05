#!/bin/bash
#SBATCH -J getFlow
#SBATCH --time=0:20:00
#SBATCH --partition=debug

. env.sh

cd /lustre/nyx/cbm/users/ogolosov/FlowAnalysis

root.exe << EOF
.x RunFlowAnalysisHADES.C("correlations")
EOF
status=$?

if [ $status -ne 0 ]
then
    echo "JOB: $JOB_ID CRASHED ON HOST: $host WITH OUTFILE $outfile_wo_path"
fi

format='+%Y/%m/%d-%H:%M:%S'

date $format
