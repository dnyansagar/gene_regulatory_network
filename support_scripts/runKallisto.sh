#!/bin/sh
set -eu
DBDIR=/proj/dnyansagar/brachyury/mm10/index/kallisto/mus.GRCm38.cds
for F in $@ ; do
  R1=`echo $F | sed s/_2/_1/`
  R2=`echo $F | sed s/_1/_2/`
  BASE=$(basename $R1 _1.fastq)
  echo "Processing $R1 $R2"
  kallisto quant -i $DBDIR -o $BASE -b 100 $R1 $R2 2>&1 | tee $BASE.log
  echo "Quantifying $BASE Done"
done
