#!/bin/bash
#SBATCH -J getCorr
#SBATCH --partition=main
#SBATCH --time=8:00:00


if [ $# -lt 1 ]
then
        echo "Usage : jobScript.sh par1 par2 par3 par4"
        exit
fi

echo "==> running enironment script ${1}"
. $1

which root

cd $2
echo "==> root"
root.exe << EOF
.x RunFlowAnalysisHADES.C("correlations","$3","$4")
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
echo "par4 = ${4}"
echo "finished!"
echo "--------------------------------"
echo ""
date $format
