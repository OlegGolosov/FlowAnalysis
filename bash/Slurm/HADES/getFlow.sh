#!/bin/bash
#SBATCH -J getFlow
#SBATCH --partition=debug
#SBATCH --time=0:20:00

echo "==> running enironment script ${1}"
bash $1

cd $2

echo "==> root"

root.exe << EOF
.x RunFlowAnalysisHADES.C("flow","$3")
EOF
status=$?

if [ $status -ne 0 ]
then
    echo "JOB: $JOB_ID CRASHED ON HOST: $host WITH OUTFILE $outfile_wo_path"
fi

format='+%Y/%m/%d-%H:%M:%S'

echo ""
echo "--------------------------------"
echo "Job with params "
echo "par1 = ${1}"
echo "par2 = ${2}"
echo "par3 = ${3}"
echo "finished!"
echo "--------------------------------"
echo ""
date $format
