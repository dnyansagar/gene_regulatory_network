#!/usr/bin/python
import sys
import os
import argparse
from Bio import SeqIO
parser=argparse.ArgumentParser(
description=''' Filters fasta sequences based on length and creates two files 
				USAGE: lengthFilter.py -ge fasta -len INT ''')
parser.add_argument('-ge', '--fasta_file', required=True, help='Multi fasta file')
parser.add_argument('-len', '--length', required=True)
args=parser.parse_args()
basename  = os.path.splitext(args.fasta_file)[0]

smallfl = open(basename+ str(args.length)+ 'short' + '.fa','w')
largefl = open(basename+ str(args.length)+ 'long' + '.fa','w')

for record in SeqIO.parse(args.fasta_file,'fasta'):
    line  = '>'+ record.description+'\n'+ record.seq
    if len(record.seq) < int(args.length):
        print>>smallfl, line
    else:
        print>>largefl, line
