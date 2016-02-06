#!/bin/sh

sleepTime=2

START=$(date +%s.%N)
git pull
sleep $sleepTime
END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
RATE= $(1/($ echo $DIFF  ) )

echo $RATE

source refersh.sh