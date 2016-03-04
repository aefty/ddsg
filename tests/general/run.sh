#!/bin/sh
echo "RUNNING TESTING"

METHOD='hdmr'
PROCESSPERGROUP=1
DIM=6
SGLEVEL=$2
GRIDTYPE=$1
HDMRmaxOrder=4
HDMRcutOff=0.0


echo 'mpirun -n' 4 './main' $METHOD $PROCESSPERGROUP $DIM $SGLEVEL $GRIDTYPE $HDMRmaxOrder $HDMRcutOff
#make clean
#make

mpirun -n 4 ./main $METHOD $PROCESSPERGROUP $DIM $SGLEVEL $GRIDTYPE $HDMRmaxOrder $HDMRcutOff