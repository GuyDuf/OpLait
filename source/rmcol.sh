#!/bin/bash

DIRECTORY=$1
cd $DIRECTORY
for f in ./*.csv; do
echo $f
echo $DIRECTORY
OUT="${f}_dropped.csv1"
cut -d"	" -f1,4-7,9-12 $f > $OUT
done
