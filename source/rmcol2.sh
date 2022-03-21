#!/bin/bash

DIRECTORY=$1
cd $DIRECTORY
for f in ./*.csv; do
echo $f
echo $DIRECTORY
OUT="${f}_dropped.csv2"
cut -d"	" -f1,4,6,7,9-13,48 $f > $OUT
done
