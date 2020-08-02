#!/usr/bin/python

import os
import re
import sys
import argparse
import pandas as pd
from Bio import SeqIO

parser=argparse.ArgumentParser(
    description=''' Extract peak sequence from genome based on bed file supplied
                        USAGE: python getPeakSequences.py -ge GENOME -be BED ''')
parser.add_argument('-ge', '--genome', required=True, help='Genome file (with path)')
parser.add_argument('-be', '--bed', required=True, help='Bed file (with path)')
parser.add_argument('--mode',
                                help="summit: get sequence 'around' the summit coords: get full peak sequence",
                                )
parser.add_argument('--around', type=int, required=False,default=100)
args=parser.parse_args()
bed = pd.read_csv(args.bed, sep='\t')
genome  = SeqIO.to_dict(SeqIO.parse(args.genome,'fasta'))

for index, row in bed.iterrows():
	sequence = genome[row['seqnames']]
	func = args.mode
	if func == 'summit':
		seq_strt, seq_end = [int(row['Summit'])-args.around, 
                                    int(row['Summit'])+args.around]
		sequence_subset = sequence.seq[seq_strt:seq_end]
		print '>'+sequence.id+'_'+str(seq_strt)+'-'+str(seq_end)+'\n'+\
        		sequence_subset
	elif func =='coords':
		seq_strt, seq_end = [int(row['start']), int(row['end'])]
		sequence_subset = sequence.seq[seq_strt:seq_end]
		print '>'+sequence.id+'_'+str(seq_strt)+'-'+str(seq_end)+'\n'+\
        	sequence_subset
