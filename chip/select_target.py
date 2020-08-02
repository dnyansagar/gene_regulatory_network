#!/usr/bin/python
import re
from datetime  import datetime
datestring = datetime.strftime(datetime.now(), '%Y-%m-%d')

def annotation_of(gene_description_file):
	xdict = dict()
	for i in open(gene_description_file):
		sl = re.split('\s+',i.strip())
		if len(sl) >1:
			xdict[sl[0]] = sl[1].lower()
		else: continue
	return xdict

def updown(fl_name):
	'''
	UP or DOWN regulation of gene based on Log2foldChange/'b' value of sleuth
	'''
	d = {}
	with open(fl_name) as f:
		for line in f:
			(key, val) = line.split()
			if float(val) >0: d[key] = 'Repressor'
			elif float(val) <0: d[key] = 'Activator'
	return d

def motif_dict(mot):
	'''
	ChIP peaks with Bra motifs
	Palindromic (full) and half palindromic (short)
	'''
	mot_dict = dict()
	short = ['TBX5','TBX4','TBX21',\
					'TBX2','TBX15','TBX1',\
					'TBR1' ,'MGA','EOMES']
	full = ['TBX19','TBX20', 'T']
	for i in open(mot):
		mot_spl = re.split('\s+', i)
		if mot_spl[1] in short or mot_spl[1] in full:
			mot_dict[mot_spl[0]]= mot_spl[1]
		else:continue
	return mot_dict

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

mmu_tra_gene_file     = '../00_data/mmu_transcript_gene.txt'
xtr_tra_gene_file 		= '../00_data/xtr_transcript_gene.txt'
mmu_tra_gene  = file_to_dict(mmu_tra_gene_file, flip=False)
xtr_tra_gene  = file_to_dict(xtr_tra_gene_file, flip=True)	

def parse_targets(chip_file_line):
	targ1 = ''; targ2 = '' 
	spl = re.split('\t', chip_file_line)
	if '_' in spl[6]:
		targ1 = xtr_tra_gene[spl[6]]; targ2 = xtr_tra_gene[spl[12]]
	elif 'ENSMUST'  in spl[6]:
		targ1 = mmu_tra_gene[spl[6]]; targ2 = mmu_tra_gene[spl[12]]
	else:
		targ1 = spl[6]; targ2 = spl[12]
	return(targ1, targ2)

def makeList(chip_file):
	'''
	Makes single list of ChIP targets including both targets called 
	using Bob's script associate_peaks.py
	/scratch/dnyansagar/brachyury/associate_genes/associate_genes.py
	/proj/rpz/Bra/strpu/01_strpu_peak_reassociate/associate_genes.py
	'''
	lst = set()
	for i in open(chip_file):
		lst.add(parse_targets(i)[0])
		lst.add(parse_targets(i)[1])
		# spl = re.split('\s+', i)
		"""if 'Supergene' in spl[6]:
			lst.add(spl[6].replace('Supergene','')) 
			lst.add(spl[12].replace('Supergene',''))
		elif '_' in spl[6]:
			lst.add(xtr_tra_gene[spl[6]])
			lst.add(xtr_tra_gene[spl[12]])
		else:
		lst.add(spl[6])
		lst.add(spl[12])"""
	return(list(sorted(lst)))

def binding_sites_per_gene(chip_file):
	'''
	Counts number of binding sites per gene in CHIP file
	'''
	bind_sites = dict()
	for i in open(chip_file):	
		spl  = re.split('\s+', i)
		peak = spl[0]+'_'+spl[1]+'-'+spl[2]
		# targ = spl[6].replace('Supergene','')
		targ = parse_targets(i)[0]
		if targ in bind_sites:
			bind_sites[targ].append(peak)
		else:
			bind_sites[targ] = [peak]
	return bind_sites

def print_binding_sites(binding_sites_dict, specie):
	'''
	Print Binding sites per Gene to a file
	'''
	bsite = open(specie+'bind_sites_per_gene_'+datestring+'.txt', 'w')
	for k, v in binding_sites_dict.items():
		print>> bsite, k,'\t',len(v)

"""def target_distance(chip_file):
	'''
	Makes dictionary of gene and its distance from closet peak
	'''
	dist_dict = dict()
	for  i in open(chip_file):
		spl = re.split('\s+', i)

		ta1 = spl[6] if '_' not in spl[6] else xtr_tra_gene[spl[6]]
		ta2 = spl[12] if '_' not in spl[12] else xtr_tra_gene[spl[12]]
		if 'ENSMUST' in 
		mu_tra_gene
		ta1 = spl[6].replace('Supergene','') if '_' not in spl[6] else xtr_tra_gene[spl[6]]
		ta2 = spl[12].replace('Supergene','') if '_' not in spl[12] else xtr_tra_gene[spl[12]]
		if spl[12] =='NA':
			dist_dict[ta1] = spl[11]
		else:
			dist_dict[ta1] = spl[11]
			dist_dict[ta2] = spl[17]
	return dist_dict
"""

def chip_lookup(vlist):
	'''
	Looks if orthologs of the ChIP target are also targets
	'''
	count = 0
	if (vlist[0]!='*' and 
			len(set([vlist[0]]).intersection(capsaspora))!=0):
		count += 1
	if (vlist[1]!='*' and 
			len(set([vlist[1]]).intersection(mmu_chip))!=0):
		count += 1
	if (vlist[2]!='*' and 
			len(set([vlist[2]]).intersection(nve_chip))!=0):
		count +=1 
	if (vlist[3]!='*' and 
			len(set([vlist[3]]).intersection(spu_chip))!=0):
		count +=1 
	if (vlist[4]!='*' and 
			len(set([vlist[4]]).intersection(xtr_chip))!=0):
		count +=1
	return count

def de_lookup(id_list,species):
	'''
	Looks if ChIP target is differentially expressed
	'''
	de_dict = dspecies[species]
	sel = []
	if id_list[1] in de_dict: sel.append(id_list[1])
	elif id_list[1] in de_dict and id_list[0] in de_dict:
		sel.extend([id_list[0], id_list[1]])
	else:sel.append(id_list[0])
	return sel

def mot_lookup(peak, species):
	'''
	Looks if ChIP targets has Bra motif
	'''
	mo_dic = mspecies[species]
	if peak in mo_dic:return True
	else:return False

tget1 =  set() ; tget2 = set()

def selection(list1, species):
	'''
	Selection of target based on score assigned
	based on how many of its orthologs are also
	ChIP targets
	'''
	if ((list1[1] > list1[3]) or 
			(list1[1] >0 and list1[3] == 0) or 
			(list1[1] >0 and list1[3] < 0) or 
			(list1[1]  == list1[3] == 0) or 
			(list1[1]  == list1[3] < 0)):
		tget1.add(list1[0])
		return [list1[0]]
	elif ((list1[1] == 0 and list1[3] >0) or 
				(list1[1] <0 and list1[3] >0)): 
		tget2.add(list1[2])
		return [list1[2]] 

def selection2(j, targ1_set, targ2_set):
	''' 
	if first target of current peak is a second target  of another gene 
	likewise if second target of current peak is a first target  of another gene

	'''
	if j[0] in targ2_set: return j[0]
	elif j[1] in targ1_set:	return j[1]
	elif (j[0] in targ2_set and j[1] in targ1_set):
		return [j[0], j[1]]
	else: None

def selection3(gene_pair, target_set, id_to_gene_dict = None):
	'''
	Selection of target based on whether it is 
	(transcription factor / mesodermal gene/ neural gene)
	'''
	if id_to_gene_dict !=None:
		if (gene_pair[0] in target_set or 
				id_to_gene_dict.get(gene_pair[0]) in target_set):
			return gene_pair[0]
		elif (gene_pair[1] in target_set or 
				id_to_gene_dict.get(gene_pair[1]) in target_set):
			return gene_pair[1]
		elif ((gene_pair[0] in target_set or 
					id_to_gene_dict.get(gene_pair[0]) in target_set) and
					(gene_pair[1] in target_set or 
					id_to_gene_dict.get(gene_pair[1]) in target_set)):
			return [gene_pair[0], gene_pair[1]]
		else:
			return None
	else:
		if gene_pair[0] in target_set:
			return gene_pair[0]
		elif gene_pair[1] in target_set:
			return gene_pair[1]
		elif ((gene_pair[0] in target_set) and (gene_pair[1] in target_set)):
			return [gene_pair[0], gene_pair[1]]
		else:
			return None
	
def undecided(list1, species):
	if list1[1] > 0 and list1[3] > 0: 
		return [list1[0],list1[2]] 
	else: return []

distance = dict()

def assign_score(chip_file_line, dictionary, length):
	spl = re.split('\t', chip_file_line.strip())
	peak = spl[0]+'_'+spl[1]+'-'+spl[2]
	tar1 = ''; score1 = 0 ; tar2 = ''; score2 = 0
	
	if len(spl) >int(length):
		tar1 = parse_targets(chip_file_line)[0]
		tar2 = parse_targets(chip_file_line)[1]
		#tar1 = spl[6] if '_' not in spl[6] else xtr_tra_gene[spl[6]]
		#tar2 = spl[12] if '_' not in spl[12] else xtr_tra_gene[spl[12]]
		#tar1 = spl[6].replace('Supergene','') if '_' not in spl[6] else xtr_tra_gene[spl[6]]
		#tar2 = spl[12].replace('Supergene','') if '_' not in spl[12] else xtr_tra_gene[spl[12]]
		#dist1 = spl[11]; dist2 = spl[17]
		#distance[tar1] = dist1; distance[tar2] = dist2
		for k,v in dictionary.iteritems():
			if tar1 in v:
				score1 = chip_lookup(v)-1
			elif tar2 in v:
				score2 = chip_lookup(v)-1
	else:
		tar1 = spl[6]  if '_' not in spl[6] else xtr_tra_gene[spl[6]]
		dist1 = spl[11]
		#distance[tar1] = dist1
		for k,v in dictionary.iteritems():
			if tar1 in v: score1 = chip_lookup(v)-1
	return [peak, tar1, score1, tar2, score2]

def primary_selection(chip_file, dictionary, length, species, de):
	'''
	Process ChIP file and assign score based 
	on how many of its orthologs are also targets
	'''
	listoflists = []; 
	selected = set(); moti = list(); und = dict()
	first = 0; second = 0;  both = 0; total = 0
	for line in open(chip_file):
		total +=1
		scored = assign_score(line, dictionary, length)
		if mot_lookup(scored[0], species):
			mot_selection1 = selection(scored[1:], species)
			#mot_selection2 = undecided(scored[1:], species)
			if mot_selection1:moti.extend(mot_selection1)
			#if mot_selection2:moti.extend(mot_selection2)
		out_selection = selection(scored[1:], species)
		if out_selection: 
			selected.add(out_selection[0])
			if out_selection[0] == scored[1]: first +=1
			elif out_selection[0] == scored[3] : second +=1
		und_selection = undecided(scored[1:], species)
		if und_selection != None and len(und_selection) >0:
			both +=1
			und[scored[0]] = und_selection
			#und.append(und_selection)
	return(und, list(selected), moti, [first, second, both, total])

def secondary_selections(undecided_dict, first_selection, motif, counts, targets, ids_to_gn_pro, 
					tget1, tget2,	specie):
	sec_selection = list(); thir_selection = list(); unselected = list()
	first = counts[0]; second = counts[1] ; 
	both = counts[2]; total = counts[3]
	#dis_dict = dispecies[specie]
	out = open(specie+'target'+datestring+'.txt', 'w')
	mot = open(specie+'motif'+datestring+'.txt', 'w')
	#dis = open(specie+'distance'+datestring+'.txt', 'w')
	undec = open(specie+'undecided'+datestring+'.txt', 'w')
	#temp = open(specie+'undecided_genes'+datestring+'.txt', 'w')
	for i, j in undecided_dict.items():
		total +=1
		second_round = selection2(j, tget1, tget2)
		#print j , '\t', second_round
		if second_round:
			if second_round == j[0]: 
				first +=1
				both -=1
				motif.append(second_round)
			elif second_round == j[1]: 
				second +=1
				both -=1
				motif.append(second_round)
			sec_selection.append(second_round)
			if mot_lookup(i, specie):motif.append(second_round)
			#if other_second in motif: motif.remove(other_third)
		else:
			unselected.append(j)
			sec_selection.append(j[0])
			first +=1
			"""third_round = selection3(i, targets, ids_to_gn_pro)
			if third_round : 
				if third_round == i[0]: first +=1
				elif third_round == i[1]: second +=1
				other_thir = i.remove(third_round)
				thir_selection.append(third_round)
				if other_thir in motif: motif.remove(other_thir)
			else: unselected.append(i)"""
	print specie, 'First selection', len(first_selection)
	print specie, 'Second selection ',len(sec_selection)
	print specie, 'Unselected', len(unselected),'\n'
	all_tar = first_selection + sec_selection + thir_selection
	#distance = [ '\t'.join([k, dis_dict[k]]) for k in all_tar]
	revised_counts = [first, second, both, total]
	for i in list(set(all_tar)): print>>out, i
	#for j in list(set(distance)): print>>dis, j
	for k in list(set(motif)):print>>mot, k
	for l in unselected: 
		index = dindex[specie]
		tmp1 = [t for t in di.values() if t[index] == l[0]]
		tmp2 = [t for t in di.values() if t[index] == l[1]]
		l0 = ''; l1 = ''
		if tmp1[0][1] !='*' and tmp1[0][4] !='*': 
			l0= mmu_gene_name[tmp1[0][1]] if tmp1[0][1] in mmu_gene_name else 'NA'
		elif tmp1[0][1] =='*' and tmp1[0][4] !='*':
			l0 = xtr_gene_name[tmp1[0][4]] if tmp1[0][4] in  xtr_gene_name else 'NA'
		elif tmp1[0][1] !='*' and tmp1[0][4] =='*':
			l0=  mmu_gene_name[tmp1[0][1]] if tmp1[0][1] in mmu_gene_name else 'NA'
		else: l0 =  'NA'
		if tmp2[0][1] !='*' and tmp2[0][4] !='*': 
			l1=  mmu_gene_name[tmp2[0][1]] if tmp2[0][1] in mmu_gene_name else 'NA' 
		elif tmp2[0][1] =='*' and tmp2[0][4] !='*':
			l1=  xtr_gene_name[tmp2[0][4]] if tmp2[0][4] in  xtr_gene_name else 'NA'
		elif tmp2[0][1] !='*' and tmp2[0][4] =='*':
			l1 =  mmu_gene_name[tmp2[0][1]] if tmp2[0][1] in mmu_gene_name else 'NA'
		else: l1= 'NA'
		print >>undec, l[0],"\t", l0,"\t", l[1],"\t",l1, "\t", tmp1[0],"\t",tmp2[0]
	return revised_counts

def ortholog_dict(oma_file):
	'''
	Create dictionary of OMA orthologs file
	'''
	d = dict()
	for line in open(oma_file):
		if line.startswith('#'):continue
		spl = re.split('\t', line.strip())
		nve_spl = filter(re.compile("^nvec").match, spl)
		cow_spl = filter(re.compile("^cowc").match, spl)
		mus_spl = filter(re.compile("^mmus").match, spl)
		xen_spl = filter(re.compile("^xtro").match, spl)
		spu_spl = filter(re.compile("^spur").match, spl)
		cow = '*' if len(cow_spl)==0 else cow_spl[0].split(':')[1].split('T')[0]
		nve = '*' if len(nve_spl)==0 else nve_spl[0].split(':')[1]
		spu = '*' if len(spu_spl)==0 else spu_spl[0].split(':')[1].split(' ')[0]
		#xen = '*' if len(xen_spl)==0 else xen_spl[0].split(':')[1]
		xen = '*' if len(xen_spl)==0 else xtr_tra_gene[xen_spl[0].split(':')[1]]
		#mus = '*' if len(mus_spl)==0 else mus_spl[0].split(':')[1]
		mus = '*' if len(mus_spl)==0 else mmu_tra_gene[mus_spl[0].split(':')[1]]
		d[spl[0]]= [cow,mus,nve,spu,xen]	
	return d

##################################################################
###################### INPUT FILES ###############################
##################################################################
### DE files (optional)
#mmu_de_filtered2018-12-26.csv nve_de_filtered2018-10-30.csv spu_de_filtered2018-10-30.csv xtr_de_filtered2018-12-26.csv
dm_de     = updown(glob.glob('../01_select_deseq/mmu_de_filtered*.csv')[0])
dn_de     = updown(glob.glob('../01_select_deseq/nve_de_filtered*.csv')[0])
ds_de     = updown(glob.glob('../01_select_deseq/spu_de_filtered*.csv')[0])
dx_de     = updown(glob.glob('../01_select_deseq/xtr_de_filtered*.csv')[0])
dspecies = {'mmu':dm_de, 'nve':dn_de, 'spu':ds_de, 'xtr':dx_de }
ofile  	= '../00_data/OMA/Output/OrthologousGroups.txt'
### Peak files 
chmmu = '../00_data/Mmu_T_ChIP.targets3.bed'
chnve = '../00_data/nematostella_Bra_peaks_score_cutoff_7_160108.targets1.txt'
chspu = '../00_data/sea_urchin_Bra_Ma_366peaks_score_cutoff_0.5_Spur2_170228.targets1.txt'
chxtr = '../00_data/Xenopus_Bra_peaks_score_cutoff_181121_fil.targets1.bed'
###
nb = print_binding_sites(binding_sites_per_gene(chnve), 'nve')
sb = print_binding_sites(binding_sites_per_gene(chspu), 'spu')
xb = print_binding_sites(binding_sites_per_gene(chxtr), 'xtr')
mb = print_binding_sites(binding_sites_per_gene(chmmu), 'mmu')
### ChIP lists
mmu_chip    = makeList(chmmu) 
nve_chip    = makeList(chnve) 
spu_chip    = makeList(chspu) 
xtr_chip    = makeList(chxtr) 
capsaspora  = set(l.strip() for l in open('../00_data/bra_ATAC.txt'))
### Motif files

nvemot = motif_dict('../../motif_analysis/Final_meme-chip-Nve/centrimo_out/Table_1e-20.csv')
xtrmot = motif_dict('../../motif_analysis/Final_meme-chip-Xtr_new_genome/centrimo_out/centrimo.csv')
mmumot = motif_dict('../../motif_analysis/Final_meme-chip-Mmu_new_paper/centrimo_out/centrimo2.csv')
spumot = motif_dict('../../motif_analysis/Final_meme-chip-Spu-All-337Peak/centrimo_out/Table_Motif.csv')

### transcription factors
human_tfs = set(l.split()[0].lower() for l in open('../00_data/human_tfs.csv'))
other_tfs = set(l.strip() for l in open('../00_data/other_tfs.txt'))
tfs 			= set(list(human_tfs)+ list(other_tfs))

### mesodermal and neural
mesodermal 	= set(l.strip().lower() for l in open('../00_data/MESODERMAL_GENES'))
neural 			= set(l.strip().lower() for l in open('../00_data/NEURAL_GENES'))
other_tar  	= set(list(mesodermal) +list(neural))
targets 		= set(list(tfs)+ list(other_tar))

### Motif  files dictionary
mspecies = {'mmu':mmumot, 'nve':nvemot, 'spu':spumot, 'xtr':xtrmot }
### Distance dictionary 
"""mmu_dist = target_distance(chmmu)
nve_dist = target_distance(chnve)
spu_dist = target_distance(chspu)
xtr_dist = target_distance(chxtr)
dispecies = {'mmu':mmu_dist, 'nve':nve_dist, 'spu':spu_dist, 'xtr':xtr_dist }"""
###
mmu_gene_name = annotation_of('../00_data/mouse_gene_names.txt')
xtr_gene_name = annotation_of('../00_data/xtr_transcript_annotation.txt')
ids_to_gn_pro = dict(mmu_gene_name, **xtr_gene_name)
###
di = ortholog_dict(ofile)
#for k, v in di.items():print k, '\t', v
dindex = {'cow':0, 'mmu':1, 'nve':2, 'spu':3, 'xtr':4}
##################################################################
##################################################################

n = primary_selection(chnve, di, 13, 'nve', de=False)
m = primary_selection(chmmu, di, 12, 'mmu', de=False)
x = primary_selection(chxtr, di, 12, 'xtr', de=False)
s = primary_selection(chspu, di, 12, 'spu', de=False)

n_select = secondary_selections(n[0], n[1], n[2], n[3], 
						targets, ids_to_gn_pro, tget1, tget2, 'nve')
m_select = secondary_selections(m[0], m[1], m[2], m[3], 
						targets, ids_to_gn_pro, tget1, tget2, 'mmu')
x_select = secondary_selections(x[0], x[1], x[2], x[3], 
						targets, ids_to_gn_pro, tget1, tget2, 'xtr')
s_select = secondary_selections(s[0], s[1], s[2], s[3], 
						targets, ids_to_gn_pro, tget1, tget2, 'spu')

print '%s\t%.2f\t%.2f\t%.2f'%('Nematostella',
				float(n_select[0])/n_select[3]*100,
				float(n_select[1])/n_select[3]*100, 
				float(n_select[2])/n_select[3]*100)
print '%s\t\t%.2f\t%.2f\t%.2f'%('Mouse',     
				float(m_select[0])/m_select[3]*100,
				float(m_select[1])/m_select[3]*100, 
				float(m_select[2])/m_select[3]*100)
print '%s\t\t%.2f\t%.2f\t%.2f'%('Xenopus',   
				float(x_select[0])/x_select[3]*100,
				float(x_select[1])/x_select[3]*100, 
				float(x_select[2])/x_select[3]*100)
print '%s\t%.2f\t%.2f\t%.2f'%('Sea urchin',  
				float(s_select[0])/s_select[3]*100,
				float(s_select[1])/s_select[3]*100, 
				float(s_select[2])/s_select[3]*100)
