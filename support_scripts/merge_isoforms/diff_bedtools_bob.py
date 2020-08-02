#!/usr/bin/python 
import re, sys, os

bedtools = 'Mmu_Bra_peaks_common_peaks.bedtools.2targets.merged.txt'
bpython = 'Mmu_Bra_peaks_common_peaks.bob_script.bothtargets.srt.txt'

di = dict()
for i in open(bedtools):
	spl = re.split('\s+', i.strip())
	peak = '#'.join(spl[0:3]) 
	target1 = spl[6];target2 = spl[12]
	di[peak] = [target1, target2]
di2 = dict()
for j in open(bpython):
	spl = re.split('\s+',j.strip())
	peak2 = '#'.join(spl[0:3])
	tar1 = spl[6] ; tar2 = spl[12]
	di2[peak2] = [tar1, tar2]
	if peak2 in di:
		if di[peak2][0] == tar1 and di[peak2][1]==tar2:
			continue
		else:
			print j,
			print di[peak2]
			#print peak.split('#'), di[peak2], [tar1, tar2]"""
