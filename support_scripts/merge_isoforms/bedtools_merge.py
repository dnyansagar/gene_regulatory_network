#!/usr/bin/python

import re
import sys
fl = sys.argv[1] #'Mmu_Bra_peaks_common_peaks.bedtools.2targets.txt'
di = dict()
for i in open(fl):
	spl = re.split('\s+', i.strip())
	peak = '#'.join(spl[0:6])
	target = spl[6:]
	if peak in di: di[peak].append(target)
	else: di[peak]=[target]
for k, v in di.items():
	if len(v)>1:
		key = k.split('#'); value1 = v[0];	value2 = v[1]
		print '\t'.join(key + value1 + value2)
	else:
		key = k.split('#'); value1 = v[0]
		value2 = ['NA','NA','NA','NA','NA','NA']
		print '\t'.join(key + value1 + value2)

