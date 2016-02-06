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

echo '========='
echo '## Test 1'
echo '========='
action1=3

threads1=$nodeCount
dim1=10
pointCout1=1000

HDMRmaxOrder1=4
HDMRsampleSize1=0
HDMRcutOff1=0.0001

SGmaxLevel1=7 #
SGcutOff1=0.0001

SGgridType1=1
folderName1='results/test1'
mkdir $folderName1
rm -rf $folderName1/*

##                        <dim><pointCout><HDMRmaxOrder><HDMRsampleSize><HDMRcutOff><SGmaxLevel><SGcutOff><SGgridType><outputFolder>
mpirun -n $threads1 ./main $dim1 $pointCout1 $HDMRmaxOrder1 $HDMRsampleSize1 $HDMRcutOff1 $SGmaxLevel1 $SGcutOff1 $SGgridType1 $folderName1 $action1
clear

echo '========='
echo '## Test 2'
echo '========='
action2=3

threads2=$nodeCount
dim2=10
pointCout2=1000

HDMRmaxOrder2=4
HDMRsampleSize2=0
HDMRcutOff2=0.0001

SGmaxLevel2=6 #
SGcutOff2=0.0001

SGgridType2=1
folderName2='results/test2'
mkdir $folderName2
rm -rf $folderName2/*

mpirun -n $threads2 ./main $dim2 $pointCout2 $HDMRmaxOrder2 $HDMRsampleSize2 $HDMRcutOff2 $SGmaxLevel2 $SGcutOff2 $SGgridType2 $folderName2 $action2
clear

echo '========='
echo '## Test 3'
echo '========='
action3=3

threads3=$nodeCount
dim3=10
pointCout3=1000

HDMRmaxOrder3=4
HDMRsampleSize3=0
HDMRcutOff3=0.0001

SGmaxLevel3=5 #
SGcutOff3=0.0001

SGgridType3=1
folderName3='results/test3'
mkdir $folderName3
rm -rf $folderName3/*

mpirun -n $threads3 ./main $dim3 $pointCout3 $HDMRmaxOrder3 $HDMRsampleSize3 $HDMRcutOff3 $SGmaxLevel3 $SGcutOff3 $SGgridType3 $folderName3 $action3
clear

echo '========='
echo '## Test 4'
echo '========='
action4=3

threads4=$nodeCount
dim4=10
pointCout4=1000

HDMRmaxOrder4=4
HDMRsampleSize4=0
HDMRcutOff4=0.0001

SGmaxLevel4=4 #
SGcutOff4=0.0001

SGgridType4=1
folderName4='results/test4'
mkdir $folderName4
rm -rf $folderName4/*

mpirun -n $threads4 ./main $dim4 $pointCout4 $HDMRmaxOrder4 $HDMRsampleSize4 $HDMRcutOff4 $SGmaxLevel4 $SGcutOff4 $SGgridType4 $folderName4 $action4
clear

echo '========='
echo '## Test 5'
echo '========='
action5=3

threads5=$nodeCount
dim5=10
pointCout5=1000

HDMRmaxOrder5=4
HDMRsampleSize5=0
HDMRcutOff5=0.0001

SGmaxLevel5=3 #
SGcutOff5=0.0001

SGgridType5=1
folderName5='results/test5'
mkdir $folderName5
rm -rf $folderName5/*

##                        <dim><pointCout><HDMRmaxOrder><HDMRsampleSize><HDMRcutOff><SGmaxLevel><SGcutOff><SGgridType><outputFolder>
mpirun -n $threads5 ./main $dim5 $pointCout5 $HDMRmaxOrder5 $HDMRsampleSize5 $HDMRcutOff5 $SGmaxLevel5 $SGcutOff5 $SGgridType5 $folderName5 $action5
clear

echo '========='
echo '## Test 6'
echo '========='
action6=3

threads6=$nodeCount
dim6=10
pointCout6=1000

HDMRmaxOrder6=4
HDMRsampleSize6=0
HDMRcutOff6=0.0001

SGmaxLevel6=2 #
SGcutOff6=0.0001

SGgridType6=1
folderName6='results/test6'
mkdir $folderName6
rm -rf $folderName6/*

##                        <dim><pointCout><HDMRmaxOrder><HDMRsampleSize><HDMRcutOff><SGmaxLevel><SGcutOff><SGgridType><outputFolder>
mpirun -n $threads6 ./main $dim6 $pointCout6 $HDMRmaxOrder6 $HDMRsampleSize6 $HDMRcutOff6 $SGmaxLevel6 $SGcutOff6 $SGgridType6 $folderName6 $action6