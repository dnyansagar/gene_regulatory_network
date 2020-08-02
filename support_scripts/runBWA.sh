#!/bin/sh
set -eu
DBDIR=/proj/dnyansagar/brachyury/mm10/ensembl/Mus_musculus.GRCm38.dna.toplevel.fa
for F in $@ ; do
  R1=`echo $F | sed s/_2/_1/`
  R2=`echo $F | sed s/_1/_2/`
  BASE=$(basename $R1 _1.fastq).srt.bam
  echo "Processing $R1 $R2"
  echo "bwa mem -t 16 $DBDIR $R1 $R2 |samtools view -Sb - |samtools sort - > $BASE"
done
