#!/bin/sh
set -eu

## this is a wrapper for mapping single end and paired end data

STAR=`which STAR` # adapt to your needs
### lsalmonis genome for 76bp reads
DBDIR=/proj/dnyansagar/brachyury/mm10/index/star/Overhang74/
#licebase/genomedata/lsal76ribo

BASE=$(basename $1 _1.fastq)
DNAME=`dirname $1`
#BASE=VOD/$BASE
#ZCAT=
#ZCAT="--readFilesCommand zcat"
THREADS=16
OPTS_P1="--outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx  --outStd Log --outFileNamePrefix $DNAME/${BASE}"

if [ "$#" -ge 3  ]; then
   echo "illegal number of parameters"
   exit 1
fi

if [ "$#" -eq 2 ]; then
CMD="$STAR --runThreadN $THREADS --outBAMsortingThreadN 10 $OPTS_P1 --genomeDir $DBDIR --readFilesIn $1 $2 2>&1 | tee $BASE.log"
echo $CMD
#$CMD
fi

if [ "$#" -eq 1 ]; then
CMD="$STAR --runThreadN $THREADS --outBAMsortingThreadN 10 $OPTS_P1 --genomeDir $DBDIR --readFilesIn $1 2>&1 | tee $BASE.log"
echo $CMD
#$CMD
fi
