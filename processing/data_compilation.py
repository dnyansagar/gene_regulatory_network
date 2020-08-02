#!/usr/bin/python 
import re
import glob
from itertools import *
from datetime  import datetime

datestring = datetime.strftime(datetime.now(), '%Y-%m-%d')

def take(n, iterable):
    "Return first n items of the iterable as a list"
    return list(islice(iterable, n))

def annotation_of(gene_description_file):
	'''
	Convert the files of annotation to dictionaries 
	'''
	xdict = dict()
	for i in open(gene_description_file):
		sl = re.split('\s+',i.strip())
		if len(sl) >1:
			desc = ' '.join(sl[1:-1])
			if '[' in desc:
				gdesc = desc.split('[')[0]
				xdict[sl[0]] = gdesc
			else:
				gdesc=desc
				xdict[sl[0]] = gdesc
		else: continue
	return xdict

def nve_annot(tsv_file, nve_dict):
    nved = dict() 
    for i in open(tsv_file):
        spl = i.split('\t')
        nv_id  =  spl[0]
        if nv_id in nve_dict: annot = nve_dict[nv_id]
        elif spl[7] != '':annot = spl[7] 
        elif spl[7] ==''  and spl[8] !='': annot = spl[8]
        elif spl[7] ==''  and spl[8] =='' and spl[10] != '': annot = spl[10]
        elif spl[7] == '' and spl[8] == '' and spl[10] == '' and spl[12] != '': annot = spl[12]
        else: annot = 'NA'
        nved[nv_id] = annot
    return nved

def file_to_dict(filename, flip=False):
  ''' 
  Converts a file with two columns to dictionary with 
  first column as key and second column as value 
  IF there is a second column
  '''
  key_val = dict()
  with open(filename) as fh: 
    for line in fh: 
      spl = re.split('\s+', line.strip())
      if len(spl)<2:continue
      else:
        if flip : key_val[spl[0]] = spl[1]
        else:key_val[spl[1]] = spl[0]
  return key_val

def file_to_dict_nve(filename, flip=False):
  ''' 
  Converts a file with two columns to dictionary with 
  first column as key and second column as value 
  IF there is a second column
  '''
  key_val = dict()
  with open(filename) as fh: 
    for line in fh: 
      spl = re.split('\t', line.strip())
      if len(spl)<2:continue
      else:
        if flip : key_val[spl[0]] = spl[1]
        else:key_val[spl[1]] = spl[0]
  return key_val



def updown(fl_name):
	'''
	Process differential expression file, 
	add 'Rep' to upregulated genes and 'Act' to Down regulated genes
	'''
	d = {}
	with open(fl_name) as f:
		for line in f:
			(key, val) = line.split()
			if float(val) >0: d[key] = 'Rep:'+str("%.2f"%float(val))
			elif float(val) <0: d[key] = 'Act:'+str("%.2f"%float(val))
	return d

def lookup_mouse(mtrans):
	'''
	Given the mouse id from the ortholog group (OMA file)
	this function checks if 
	1) it is a ChIP target 
	2) If that ChIP peak has a Bra motif
	3) If that gene is dfferentially expressed in mouse knockout
	'''
	#mgene = mmu_tra_gene[mtrans]
	mchip = 'Yes' if mtrans in mmu_chip else 'No'
	mmotif = 'Yes' if mtrans in mmu_motif else 'No'
	mde = dm_de[mtrans] if mtrans in dm_de else 'No'
	mmuname  = mmu_gene_name[mtrans] if mtrans in mmu_gene_name else 'NA'
	return [mchip, mmotif, mde, mmuname]

def lookup_xenopus(xtrans):
	'''
	Given the xenopus id from the ortholog group (OMA file)
	this function checks if 
	1) it is a ChIP target 
	2) If that ChIP peak has a Bra motif
	3) If that gene is dfferentially expressed in xenopus knockout
	'''
	# xgene = xtr_tra_gene[xtrans] #.split('|')[0]
	xchip = 'Yes' if xtrans in xtr_chip else 'No'
	xmotif = 'Yes' if xtrans in xtr_motif else 'No'
	xde = dx_de[xtrans] if xtrans in dx_de else 'No'
	xtrName = xtr_gene_name[xtrans] if xtrans in xtr_gene_name else 'NA'
	return [xchip, xmotif, xde, xtrName]

def lookup_seaurchin(urchingene):
	'''
	Given the sea urchin id from the ortholog group (OMA file)
	this function checks if 
	1) it is a ChIP target 
	2) If that ChIP peak has a Bra motif
	3) If that gene is dfferentially expressed in mouse knockdown
	'''
	spgn = urchingene.split(' ')[0]
	schip = 'Yes' if spgn in spu_chip else 'No'
	smotif ='Yes' if spgn in spu_motif else 'No'
	sde = ds_de[spgn] if spgn in ds_de else 'No'
	return [schip, smotif, sde]

def lookup_nematostella(ngene):
	'''
	Given the nematostella id from the ortholog group (OMA file)
	this function checks if 
	1) it is a ChIP target 
	2) If that ChIP peak has a Bra motif
	3) If that gene is dfferentially expressed in nematostella knockdown
	'''
	nchip = 'Yes' if ngene in nve_chip else 'No'
	nmotif = 'Yes' if ngene in nve_motif else 'No'
	nde = dn_de[ngene] if ngene in dn_de else 'No'
	nvename = nve_gene_name[ngene] if ngene in nve_gene_name else 'NA'
	return [nchip, nmotif, nde, nvename]

def lookup_capsaspora(cagn):
	'''
	Given the capsaspora id from the ortholog group (OMA file)
	this function checks if it is in ATAC-seq peaks targets with Bra motifs 
	
	'''
	cchip = cmotif ='Yes' if cagn.split('T')[0] in capsaspora else 'No'
	return [cchip, cmotif]

def other_lists(species, acc):
	chi = False; dif = False
	"""if species == 'mmus':
		if mmu_tra_gene[acc] in chipList[species]:chi = True
		#if mmu_tra_gene[acc] in deList[species]: dif =  True
		if acc in deList[species]: 
			dif =  True
			print acc
	else:"""		
	if acc in chipList[species]:chi = True
	if acc in deList[species]: dif = True
	return(chi,dif)

def ortholog_dict(oma_file):
	'''
	Create dictionary of OMA orthologs file wherein 
	the OMA_id is the key and the accessions as values
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
		xen = '*' if len(xen_spl)==0 else xtr_tra_gene[xen_spl[0].split(':')[1]]
		mus = '*' if len(mus_spl)==0 else mmu_tra_gene[mus_spl[0].split(':')[1]]
		d[spl[0]]= [cow,mus,nve,spu,xen]	
	return d
	
def create_table(gene_list, dictionary, species, sp_list, file_write ):
	'''
	The main function which calls other function and 
	saves their return value in a dictionary. 
	'''
	print 'Ortholog Table for: ', species
	fl_write = open(file_write, 'w')
	chfl = open(species+'chip'+datestring+'.txt', 'w')
	defl = open(species+'de'+datestring+'.txt', 'w')
	for i in open(gene_list):
		gnid  = i.strip()
		all_list = ['mchip', 'mmotif', 'mde', 
                'nchip', 'nmotif', 'nde', 
                'schip', 'smotif', 'sde', 
                'xchip', 'xmotif', 'xde', 
                'cchip', 'cmotif',
                'mname', 'nname', 'xname']
		# set default values 
		ortho = 'NA' ;found = False
		fix = dict.fromkeys(all_list, 'NA')
		for k, v in dictionary.iteritems():
			if gnid in v:
				found = True; ortho = k
				ch,de = other_lists(species, gnid)
				if ch :print>>chfl, ortho
				if de :print>>defl, ortho
				if v[0] == '*': 
					fix['cchip'] = fix['cmotif'] ='NA'
				else: 
					fix['cchip'], fix['cmotif'] = lookup_capsaspora(v[0]) 
				if v[1] =='*': 	
					fix['mchip'] = fix['mmotif'] = fix['mde'] =  fix['mname']= 'NA'
				else: 
					fix['mchip'], fix['mmotif'], fix['mde'], fix['mname'] = lookup_mouse(v[1])
				if v[2] == '*':	
					fix['nchip']= fix['nmotif']= fix['nde'] = fix['nname']= 'NA' 
				else: 
					fix['nchip'], fix['nmotif'], fix['nde'], fix['nname'] = lookup_nematostella(v[2])
				if v[3] == '*': 
					fix['schip']= fix['smotif']= fix['sde'] ='NA'
				else: 
					fix['schip'], fix['smotif'], fix['sde'] =	lookup_seaurchin(v[3])
				if v[4]== '*': 
					fix['xchip']= fix['xmotif']= fix['xde'] ='NA'
				else: 
					fix['xchip'], fix['xmotif'], fix['xde'], fix['xname'] = lookup_xenopus(v[4])
				break		
			if not found:
				ortho = 'NA'
				sp_dict = dict(zip(sp_list, functionList[species](gnid)))
				fix.update(sp_dict)
		if ('Yes' in fix.values() or 
				len(filter(re.compile("^Act:").match, fix.values())) !=0 or 
				len(filter(re.compile("^Rep:").match, fix.values())) !=0):
			print>>fl_write, '\t'.join([gnid, ortho, \
												fix['mchip'],fix['mmotif'],fix['mde'],\
												fix['xchip'],fix['xmotif'],fix['xde'],\
			 	  							fix['schip'],fix['smotif'],fix['sde'],\
												fix['nchip'],fix['nmotif'],fix['nde'],\
												fix['cchip'],fix['cmotif'],
   				  						fix['nname'], fix['mname'], fix['xname']])

###################################################################
####################### INPUT FILES ###############################
###################################################################
## Mouse ChIP, Motif, DE
mmu_chip	= set(l.strip() for l in open(glob.glob('../02_select_target/mmutarget*.txt')[0]))
mmu_motif	= set(l.strip() for l in open(glob.glob('../02_select_target/mmumotif*.txt')[0]))
dm_de     = updown(glob.glob('../01_select_deseq/mmu_de_filtered*.csv')[0]) 
## Nematostella 
nve_chip	= set(l.strip() for l in open(glob.glob('../02_select_target/nvetarget*.txt')[0]))
nve_motif	= set(l.strip() for l in open(glob.glob('../02_select_target/nvemotif*.txt')[0]))
dn_de     = updown(glob.glob('../01_select_deseq/nve_de_filtered*.csv')[0])
## Sea urchin 
spu_chip	= set(l.strip() for l in open(glob.glob('../02_select_target/sputarget*.txt')[0]))
spu_motif	= set(l.strip() for l in open(glob.glob('../02_select_target/spumotif*.txt')[0]))
ds_de     = updown(glob.glob('../01_select_deseq/spu_de_filtered*.csv')[0])
## Xenopus
xtr_chip	= set(l.strip() for l in open(glob.glob('../02_select_target/xtrtarget*.txt')[0]))
xtr_motif	= set(l.strip() for l in open(glob.glob('../02_select_target/xtrmotif*.txt')[0]))
dx_de    	= updown(glob.glob('../01_select_deseq/xtr_de_filtered*.csv')[0])
## Capsaspora 
capsaspora 	= set(l.strip() for l in open('../00_data/bra_ATAC.txt'))
## Other files 
mmu_tra_gene_file			= '../00_data/mmu_transcript_gene.txt'
xtr_tra_gene_file 		= '../00_data/xtr_transcript_gene.txt'
mmu_gene_description	= '../00_data/mouse_gene_names.txt'
xtr_gene_description  = '../00_data/xtr_gene_annotation_uniq.txt'
nve_annotation 	= '../00_data/nve_master_ortholog_table.csv'
ann 		= '../00_data/nve_master_ortholog_table.20190111.tsv'
nve_gri = '/mnt/evo/Bra/processing/00_data/nve_master_ortholog_grisha.csv'
nve_new = file_to_dict_nve(nve_gri, flip = True)


## Orthologs with OMA 
ofile   	    = '../00_data/OMA/Output/OrthologousGroups.txt'
## Annotation files 
mmu_gene_name = annotation_of(mmu_gene_description)
xtr_gene_name = annotation_of(xtr_gene_description)
#nve_gene_name = annotation_of(nve_annotation)
nve_gene_name = nve_annot(ann, nve_new)
## Mouse transcript to gene 
mmu_tra_gene  = file_to_dict(mmu_tra_gene_file, flip = False)
xtr_tra_gene 	= file_to_dict(xtr_tra_gene_file, flip = True)
''' 	
I made dictionaries of some of the input data so as to be 
easier to call the right function based on species. 
'''
functionList = {'capsa': lookup_capsaspora, 'nemve':lookup_nematostella,
	'spur':lookup_seaurchin, 'xtro':lookup_xenopus, 'mmus':lookup_mouse}

chipList = {'capsa':[x.split('T')[0] for x in capsaspora],  
						'xtro':xtr_chip,'nemve':nve_chip, 
						'spur':spu_chip,'mmus':mmu_chip}

deList = {'capsa':dict(),'nemve':dn_de, 'spur':ds_de, 'xtro':dx_de, 'mmus':dm_de}

##################################################################
##################################################################
di = ortholog_dict(ofile)
n_items = take(10, di.iteritems())
print n_items
### MAKE TABLES ###
cow_list = create_table('../00_data/cow_list', di, 'capsa', ['cchip', 'cmotif', 'cde'], 
	'cow_table'+datestring+'.txt')
nve_list = create_table('../00_data/nve_list', di, 'nemve', ['nchip', 'nmotif', 'nde', 'nname'], 
  'nve_table'+datestring+'.txt')
spu_list = create_table('../00_data/spu_list', di, 'spur', ['schip', 'smotif', 'sde'], 
	'spu_table'+datestring+'.txt')
xtr_list = create_table('../00_data/xtr_list', di, 'xtro', ['xchip', 'xmotif', 'xde', 'xname'], 
	'xtr_table'+datestring+'.txt')
mmu_list = create_table('../00_data/mmu_list', di, 'mmus', ['mchip', 'mmotif', 'mde', 'mname'], 
	'mmu_table'+datestring+'.txt')

##################################################################
##################################################################
