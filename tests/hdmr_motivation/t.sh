#!/bin/sh

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
dim=10
pointCout=1

HDMRmaxOrder=4
HDMRsampleSize=0
HDMRcutOff=0.0001

SGmaxLevel=3 #
SGcutOff=0.0

SGgridType=1
folderName='results/test99'
mkdir $folderName
rm -rf $folderName/*

##                        <dim><pointCout><HDMRmaxOrder><HDMRsampleSize><HDMRcutOff><SGmaxLevel><SGcutOff><SGgridType><outputFolder>
time -p mpirun -n $threads ./main $dim $pointCout $HDMRmaxOrder $HDMRsampleSize $HDMRcutOff $SGmaxLevel $SGcutOff $SGgridType $folderName $action