#!/bin/bash
#SBATCH -J convert
#SBATCH --time=8:00:00
#SBATCH --partition=main

if [ $# -lt 1 ]
then
    echo "Usage : jobScript.sh par1 par2 par3"
    sleep 3
    exit
fi

echo "==> running enironment script ${1}"
bash $1

cd $2
echo "==> root"
root.exe << EOF
.x RunConvertationHADES.C ("$3","$4")
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
