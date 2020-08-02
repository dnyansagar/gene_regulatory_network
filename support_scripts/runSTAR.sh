#!/bin/sh
set -eu

for F in $@ ; do
  R1=`echo $F | sed s/_2/_1/`
  R2=`echo $F | sed s/_1/_2/`
  echo "Processing $R1 $R2"
  ./star-wrapper.sh $R1 $R2
  echo "Done"
done
