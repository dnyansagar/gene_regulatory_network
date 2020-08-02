#!/usr/bin/python
import re, sys, os 

def file_to_dict(filename, flip=False):
  ''' 
  Converts a file with two columns to dictionary with 
  first column as key and second column as value 
  IF there is a second column
  '''
  key_val = dict()
  with open(filename) as fh: 
    for line in fh: 
      spl = re.split('\s+', line)
      if len(spl)<2:continue
      else:
        if flip : key_val[spl[0]] = spl[1]
        else:key_val[spl[1]] = spl[0]
  return key_val

tra_gene =file_to_dict('/mnt/evo/Bra/processing/00_data/mmu_transcript_gene.txt', flip=False)
targetfl = 'Mmu_T_ChIP.tra_targets1.bed'
for i in open(targetfl).readlines()[1:]:
	spl = i.split('\t')
	t1 = spl[6].split('.')[0]
	t2 = spl[12].split('.')[0]
	if tra_gene[t1] ==tra_gene[t2]: continue
	else:	print tra_gene[t1], tra_gene[t2]

