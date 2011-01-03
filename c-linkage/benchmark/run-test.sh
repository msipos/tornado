#!/bin/bash

for i in {0..15}
do
echo "Test run $i"

TEST=n1000
echo "Test $TEST"
echo "c-linkage"
./run-clinkage.sh $TEST
echo "mothur"
./run-mothur.sh $TEST

TEST=n2000
echo "Test $TEST"
echo "c-linkage"
./run-clinkage.sh $TEST
echo "mothur"
./run-mothur.sh $TEST

TEST=n4000
echo "Test $TEST"
echo "c-linkage"
./run-clinkage.sh $TEST
echo "mothur"
./run-mothur.sh $TEST

TEST=n6000
echo "Test $TEST"
echo "c-linkage"
./run-clinkage.sh $TEST
echo "mothur"
./run-mothur.sh $TEST

TEST=n8000
echo "Test $TEST"
echo "c-linkage"
./run-clinkage.sh $TEST
echo "mothur"
./run-mothur.sh $TEST

done
