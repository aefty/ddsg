#!/bin/sh

sleepTime=2

START=$(date +%s.%N)
git pull
sleep $sleepTime
END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
echo DIFF

source refersh.sh