#!/bin/sh


maxLevel=2;
dataPoints=100; # The count is zero based , so add +1

echo '## Testing Grid Type 1'
echo '==========================='

SGgridType=1
folderName='results/test1'
mkdir $folderName
rm -rf $folderName/*

##                        <SGgridType><maxLevel><folderName>
time -p mpirun -n 1 ./main $SGgridType $maxLevel $dataPoints $folderName


echo '## Testing Grid Type 2'
echo '==========================='

SGgridType=2
folderName='results/test2'
mkdir $folderName
rm -rf $folderName/*

##                        <SGgridType><maxLevel><folderName>
time -p mpirun -n 1 ./main $SGgridType $maxLevel $dataPoints $folderName



echo '## Testing Grid Type 3'
echo '==========================='

SGgridType=3
folderName='results/test3'
mkdir $folderName
rm -rf $folderName/*

##                        <SGgridType><maxLevel><folderName>
time -p mpirun -n 1 ./main $SGgridType $maxLevel $dataPoints $folderName



echo '## Testing Grid Type 4'
echo '==========================='

SGgridType=4
folderName='results/test4'
mkdir $folderName
rm -rf $folderName/*

##                        <SGgridType><maxLevel><folderName>
time -p mpirun -n 1 ./main $SGgridType $maxLevel $dataPoints $folderName
