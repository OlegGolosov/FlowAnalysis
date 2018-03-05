#!/bin/bash -x

#SRC=/home/basov/ovgol/FVC

cd $SRC;
root.exe -b -l -q 'Generate.cxx+ ("'$OUT'/generated_'$NUM'")' > $OUT/Logs/generated_$NUM.log;
#root.exe -b -l -q 'Generate.cxx+';


#$1 ecm
#$2 nev


#TMPDIR=$MAINPATH/TMP

#echo 1
#cp {inputfile,create_reduced_tree.C,reco.C,runMC.C,runqmd.bash,urqmd.x86_64} $TMPDIR

#cd $TMPDIR

#echo 2
#cat inputfile | sed "s/NEV/$2/g" >> tmp && mv tmp inputfile -f
#cat inputfile | sed "s/ECM/$1/g" >> tmp && mv tmp inputfile -f

#. runqmd.bash
#echo 3
#cat runMC.C | sed "s/NEV/$2/g" >> tmp && mv tmp runMC.C -f
#cat reco.C | sed "s/NEV/$2/g" >> tmp && mv tmp reco.C -f

#echo 4
#. $MAINPATH/Soft/MPDRoot/build/config.sh
#. /mnt/pool/1/svintsov/MPDRoot/build/config.sh

#SEED=$(od -An -N4 -i < /dev/urandom | awk '{print $1}'| cut -c2-7)


#root.exe -b -l -q runMC.C > $MAINPATH/mpd_data/"$1"gev/reco/log/mpddst_reduced_"$SEED"._runmc.o
#echo 6
#root.exe -b -l -q reco.C > $MAINPATH/mpd_data/"$1"gev/reco/log/mpddst_reduced_"$SEED"._reco.o
#echo 7
#root.exe -b -l -q 'create_reduced_tree.C("mpddst.root","mpddst_reduced.root")' > $MAINPATH/mpd_data/"$1"gev/reco/log/mpddst_reduced_"$SEED"._createreduced.o
#echo 8
#mv mpddst_reduced.root $MAINPATH/mpd_data/"$1"gev/reco/mpddst_reduced_"$SEED".root
#mv mpddst.root $MAINPATH/mpd_data/"$1"gev/reco/mpddst_"$SEED".root

#cd $INPUT
