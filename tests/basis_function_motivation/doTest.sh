#!/bin/sh
GRID_TYPE=2

echo "Batch Test for Grid Type : "$GRID_TYPE

echo "Clear Surplus Files ...."
rm -rf surplus_hdmr/f1/*
rm -rf surplus_hdmr/f2/*
rm -rf surplus_hdmr/f3/*
rm -rf surplus_hdmr/f4/*
rm -rf surplus_hdmr/f5/*
rm -rf surplus_hdmr/f6/*

rm -rf surplus_sg/f1/*
rm -rf surplus_sg/f2/*
rm -rf surplus_sg/f3/*
rm -rf surplus_sg/f4/*
rm -rf surplus_sg/f5/*
rm -rf surplus_sg/f6/*

nodeCount=1

echo '## Simple Test'
echo '=============='
action=3

threads=$nodeCount
dim=2
pointCout=100

HDMRmaxOrder=1
HDMRsampleSize=0
HDMRcutOff=0.0001

SGmaxLevel=4 #
SGcutOff=0.0

SGgridType=$GRID_TYPE
folderName='results/unitTest'
mkdir $folderName
rm -rf $folderName/*

##                        <dim><pointCout><HDMRmaxOrder><HDMRsampleSize><HDMRcutOff><SGmaxLevel><SGcutOff><SGgridType><outputFolder>
time -p mpirun -n $threads ./main $dim $pointCout $HDMRmaxOrder $HDMRsampleSize $HDMRcutOff $SGmaxLevel $SGcutOff $SGgridType $folderName $action