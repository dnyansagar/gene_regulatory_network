import re
import os
import sys
import scipy 
import pandas
import itertools
#import upsetplot
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from venn import venn
from collections import OrderedDict
from matplotlib_venn import venn2
#from upsetplot import from_contents
#from upsetplot import generate_counts, plot
from datetime  import datetime
datestring = datetime.strftime(datetime.now(), '%Y-%m-%d')

######### Plot Styles ############
#plt.xkcd()
plt.style.use('ggplot')
#plt.style.use('fivethirtyeight')
#plt.style.use('dark_background')
#plt.style.use('seaborn-dark')
#plt.style.use('seaborn-dark-palette')
#plt.style.use('bmh')
#plt.style.use('seaborn-bright')
#plt.style.use("seaborn-colorblind")
#################################

path = sys.argv[1] #  Takes meme-chip result directory as INPUT
plots = path+'plots/'
try:
    os.mkdir(plots)
except OSError:
    print ("Creation of the directory %s failed" % plots)
else:
    print ("Successfully created the directory %s " % plots)

####################################################

def save(path, ext='png', close=True, verbose=True, transparent=False):
    """
    jhamrick/savefig.py at (https://gist.github.com/jhamrick/5320734)
    Save a figure from pyplot.
    Parameters
    ----------
    path : string
        The path (and filename, without the extension) to save the
        figure to.

    ext : string (default='png')
        The file extension. This must be supported by the active
        matplotlib backend (see matplotlib.backends module).  Most
        backends support 'png', 'pdf', 'ps', 'eps', and 'svg'.

    close : boolean (default=True)
        Whether to close the figure after saving.  If you want to save
        the figure multiple times (e.g., to multiple formats), you
        should NOT close it in between saves or you will have to
        re-plot it.

    verbose : boolean (default=True)
        Whether to print information about when and where the image
        has been saved.
    """
    # Extract the directory and filename from the given path
    directory = os.path.split(path)[0]
    filename = "%s.%s" % (os.path.split(path)[1], ext)
    if directory == '': directory = '.'
    # If the directory does not exist, create it
    if not os.path.exists(directory): os.makedirs(directory)
    # The final path to save to
    savepath = os.path.join(directory, filename)
    if verbose: print("# Saving figure to '%s'..." % savepath),
    # Actually save the figure
    plt.savefig(savepath)
    # Close it
    if close: plt.close()
    if verbose: print("# Done")
"""
d = {'pax': ["PAX5","Pax2","PAX4","Pax6"],
        'hox' : ["HOXA5","Rhox11" ,"HOXA2","HOXB2",\
            "HOXB3","Hoxb5","HOXC13","Hoxd3"],
        'pou'     : ["POU6F1","POU2F1","POU3F2","POU3F3",\
            "POU3F4","POU6F2"],
        'fox'     : ["FOXF2","FOXD1","FOXC1","FOXL1","Foxq1",\
            "FOXI1","Foxa2","FOXA1","FOXO3","FOXH1",\
            "FOXG1","Foxj2","FOXB1","FOXC2","FOXD2",\
            "FOXO4","FOXO6","FOXP3","Foxk1"],
        'sox'     : ['Sox2','Sox17','Sox3','Sox6','Sox1','SOX4',\
            'SOX8','SOX9','Sox11','Sox1','SOX21', 'D'],
        'bra-half'   : ['TBX5', 'TBX4', 'TBX21', 'TBX2','TBX15', \
            'TBX1','TBR1', 'MGA', 'EOMES'],
        'bra-full'   : ['T', 'TBX19', 'TBX20']}
"""
def merge_family_motif_counts(dictionary, preDict): # Approach 1
    new_dict = dict()
    for k, v in dictionary.items():
        names = [name for name, seq in preDict.items() if k in seq]
        if names !=[]:
            if names[0] in new_dict: new_dict[names[0]] += int(v)
            else: new_dict[names[0]] = int(v)
        else:
            if "others" in new_dict: new_dict["others"] += int(v)
            else: new_dict["others"] = int(v)
    return new_dict

def file_to_dict(filename, flip=False): #Approach 1 & 2 

    ''' Converts a file with two columns to dictionary with 
    first column as key and second column as value IF there is a second column '''
    key_val = dict()
    with open(filename) as fh:
        for line in fh:
            spl = re.split('\s+', line)
            if len(spl)<2:continue
            else:
                if flip : key_val[spl[0]] = spl[1]
                else:key_val[spl[1]] = spl[0]
    return key_val
motif_names = file_to_dict('../JASPAR/jasparNames.txt', flip=True)

def file_column_to_list(fl_name, column_no): # Approach 2

    lst = []
    for line in open(fl_name):
        spl = re.split("\t", line.strip().replace('"',''))[column_no-1]
        lst.append(spl.replace(" ",""))
    return lst

family = {
        'bHLH'  : file_column_to_list('../JASPAR/JASPAR_bHLH.tsv',2),
        'hmg'   : file_column_to_list('../JASPAR/JASPAR_High_mobility.tsv',2),
        'pax'   : file_column_to_list('../JASPAR/JASPAR_Paired_box.tsv',2),
        'fox'   : file_column_to_list('../JASPAR/JASPAR_Forkhead_box.tsv',2),
        'homeo' : file_column_to_list('../JASPAR/JASPAR_Homeo_domain.tsv',2),
        'tbox1' : file_column_to_list('../JASPAR/JASPAR_T-Box1.tsv',2),
        'tbox2' : file_column_to_list('../JASPAR/JASPAR_T-Box2.tsv',2)
        }


def distance_from_center(m_start, m_end, seq_start, seq_end): # Approach 1 & 2

    '''Distance form middle of the sequence to the center of motif'''
    middle = round((int(seq_end)-int(seq_start))/2)
    location  = round((int(m_end) - int(m_start))/2) + int(m_start)
    distance = abs(middle-location)
    return distance

def merge_dicts(x, y): # Approach 1 & 2

    '''Given two dicts, merge them into a new dict'''
    z = x.copy()
    z.update(y)
    return z

def cm_names(lis): # Approach 1

    '''Give the list of jaspar ids get their names from centrimo'''
    new_lis = list()
    for i in lis:
        if i in cm:new_lis.append(cm[i])
        else:new_lis.append(i)
    return new_lis


def get_names(dictionary, centrimo, dtomtom, mtomtom): # Approach 1

    ''' Get names of motifs from centrimo, dreme_tomtom or meme_tomtom'''
    new_dictionary = dict()
    for k, v in dictionary.items():
        if k in centrimo: new_dictionary[centrimo[k][0]] =v
        elif k in dtomtom: new_dictionary[dtomtom[k]] =v
        elif k in mtomtom: new_dictionary[mtomtom[k]] =v
        else: new_dictionary[k] = v
    return new_dictionary

def overlap(ranges): # Approach 1 & 2

    ''' Find overlapping ranges from the list of ranges '''
    oList = []
    noList = []
    for i, j in itertools.combinations(ranges,2):
        start1, end1 = i
        start2, end2 = j
        if set(range(start1, end1)).intersection(range(start2, end2)):
            oList.append([[start1, end1], [start2, end2]])
        else:
            if not any(i in e for e in oList): noList.append(i)
            if not any(j in e for e in oList): noList.append(j)
    return (oList,noList)

def overlap_select(family_dict, lst_of_lst): # Approach  2
    ''' From the overlapping pair select the brachyury over other motifs OR
        brachyury palindromic over half palindrome motif OR
        if not Bra select based on q_value(FDR adjusted p_value)
    '''
    list_of_pairs = lst_of_lst[0]
    list_of_dist = lst_of_lst[1]
    list_of_qval = lst_of_lst[2]
    slist = [[],[],[]]
    for index,  p in enumerate(list_of_pairs):
        pair = p
        names = [[name for name, seq in family_dict.items() if i in seq] for i in pair]
        if names[0] != [] and names[1] != []:
            if names[0][0] == 'tbox1' and names[1][0] == 'tbox2': 
                slist[0].append(pair[0])
                slist[1].append(list_of_dist[index][0])
                slist[2].append(list_of_qval[index][0])
            elif names[0][0] == 'tbox2' and names[1][0] == 'tbox1': 
                slist[0].append(pair[1])
                slist[1].append(list_of_dist[index][1])
                slist[2].append(list_of_qval[index][1])
            else:
                minq = min(list_of_qval[index])
                index_minq = list_of_qval[index].index(minq)
                slist[0].append(pair[index_minq])
                slist[1].append(list_of_dist[index][index_minq])
                slist[2].append(list_of_qval[index][index_minq])
        elif names[0] == [] and names[1] != []:
            slist[0].append(pair[1])
            slist[1].append(list_of_dist[index][1])
            slist[2].append(list_of_qval[index][1])
        elif names[0] != [] and names[1] == []:
            slist[0].append(pair[0])
            slist[1].append(list_of_dist[index][0])
            slist[2].append(list_of_qval[index][0])
    return(slist)

def coOccurrence(preDict, fmo_dict, mtomtom, dtomtom, motif_names): # Approach 1
	'''Check the coccurneces of motifs''' # Superceded 
	other = merge_dicts(mtomtom, dtomtom)
	dic = dict()
	for k, v  in fmo_dict.items():		
		for i in v[0]:
			if i in motif_names:
				mo = motif_names[i]
				names = [name for name, seq in preDict.items() if mo in seq]
			elif i in other: 
				mo = motif_names[other[i][0]]
				names = [name for name, seq in preDict.items() if mo in seq]
			if names !=[]:
				
				if names[0] in dic: dic[names[0]].add(k)
				else: dic[names[0]] = set([k])
			else:continue
	return dic

def select(criteria, q25, lst_Of_lst): # Approach 1 
    ''' Based on the criteria provided get the best motif in peak''' # Superceded
    f_name, f_score, f_distance, f_qval, f_coords, f_cscore, f_pdf = lst_Of_lst
    if criteria == 'qval':
        ind = f_qval.index(np.array(f_qval, dtype=float).min())
    elif criteria == 'score':
        ind = f_score.index(max(f_score)) 
    elif criteria == 'distance':
        ind = f_distance.index(min(f_distance))
    elif criteria == 'cscore': 
        ind = f_cscore.index(min(f_cscore))
    elif criteria =='combined':
        cPval = [i*j for i, j in zip(f_qval, f_cscore)]
        ind = cPval.index(min(cPval))
    elif criteria == 'centrimo':
        ind = f_cscore.index(min(f_cscore))
    elif criteria == 'pdf':
        pdf_val = [i*j for i, j in zip(f_qval, f_pdf)]
        ind = pdf_val.index(min(pdf_val))
        # if float(f_score[ind]) > q25:     
    return (f_name[ind], f_score[ind], f_distance[ind], f_qval[ind], f_cscore[ind])

def tomtom_parse(path, folder, fil): # Approach 1 & 2 
    '''Get the database entry for DREME or MEME predicted motif'''
    tom = dict()
    for dirname, dirs, files in os.walk(path):
        if folder in dirname:
            for filename in files:
                if filename == fil:
                    for line in open(dirname+'/'+filename).readlines()[1:]:
                        spl = line.strip().split('\t')
                        criteria = float(spl[3])
                        if spl[0] in tom:
                            if criteria < float(tom[spl[0]][1]):
                                tom[spl[0]][0] = spl[1]
                            else:
                                continue
                        else: 
                            tom[spl[0]]= [spl[1], spl[3]]
    return tom

mt = tomtom_parse(path, 'meme_tomtom_out', 'tomtom.txt')
dt = tomtom_parse(path, 'dreme_tomtom_out', 'tomtom.txt')

def centrimo_parse(path): # Approach  1 & 2
    ''' Get centrimo motifs list along with their names'''
    cent = dict()
    for dirname, dirs, files in os.walk(path):
        if 'centrimo_out' in dirname:
            for filename in files:
                if filename == 'centrimo.txt':
                    for line in open(dirname+'/'+filename).readlines()[1:]:
                        if line.startswith('#'):continue
                        cspl = re.split('\s+', line)
                        if float(cspl[4]) < float(1e-1):
                            cent[cspl[2]] = [cspl[3],float(cspl[5])]
    return cent

cm = centrimo_parse(path)

fimo_manual = './Spur5/fimo.tsv'

def web_fimo(fimo_file, centrimo, centrimo_test =True): # Approach 1 & 2

        ''' Select the motifs based on the Fimo score and 
                calculate the distance from the middle of the sequence'''
        fimo_motif = dict(); fimo_score = list()
        for line in open(fimo_file).readlines():
                if line.startswith('#'):continue
                spl = line.split('\t')
                st = spl[2].split("_")[-1].split("-")[0]                        
                en = spl[2].split("_")[-1].split("-")[1]
                dist = distance_from_center(spl[3], spl[4], st, en)
                pdf = 1
                qval = float(spl[7])
                coords = [int(spl[3]), int(spl[4])]
                fimo_score.append(qval)
                
                if spl[2] in fimo_motif and float(spl[-2]) <1:
                        if centrimo_test and spl[0] in centrimo:
                                cscore = centrimo[spl[0]][1]
                                lst = [spl[0], spl[6], dist, qval, coords,cscore, pdf]
                                [fimo_motif[spl[2]][ix].append(lst[ix]) for ix, (a, b) in \
                                        enumerate(zip(fimo_motif[spl[2]], lst))]
                        elif not centrimo_test:  
                                lst = [spl[0], spl[6], dist, qval, coords, 1 , pdf]
                                [fimo_motif[spl[2]][ix].append(lst[ix]) for ix, \
                                        (a, b) in enumerate(zip(fimo_motif[spl[2]], lst))]
                elif float(spl[-2]) <1:
                        if  centrimo_test and spl[0] in centrimo:
                                cscore = centrimo[spl[0]][1]
                                fimo_motif[spl[2]] = [[spl[0]], [spl[6]], 
                                        [dist], [qval], [coords], [cscore], [pdf]]
                        elif not centrimo_test: 
                                fimo_motif[spl[2]] = [[spl[0]], [spl[6]], 
                                [dist], [qval], [coords], [1] , [pdf]]
                
        #print(fimo_motif)
        return (fimo_motif,fimo_score)

def fimo_parse(path, centrimo, centrimo_test =True): # Approach 1 & 2

    ''' Select the motifs based on the Fimo score and 
    calculate the distance from the middle of the sequence'''
    fimo_motif = dict(); fimo_score = list()
    for dirname, dirs, files in os.walk(path):
        if 'fimo_out' in dirname:
            for filename in files:
                if filename == 'fimo.txt':
                    for line in open(dirname+'/'+filename).readlines():
                        if line.startswith('#'):continue
                        spl = line.split('\t')
                        st,en = spl[1].split("_")[-1].split("-")
                        dist = distance_from_center(spl[2], spl[3], st, en)
                        pdf = float(scipy.stats.norm(0,100).pdf(dist))
                        motLength = int(spl[3])-int(spl[2]) # Superceded
                        normScore = float(spl[5])/motLength # Superceded
                        qval = float(spl[6])
                        coords = [int(spl[2]), int(spl[3])]
                        fimo_score.append(qval)
                        if spl[1] in fimo_motif:
                            if centrimo_test and spl[0] in centrimo:
                                cscore = centrimo[spl[0]][1]
                                lst = [spl[0], spl[5], dist, qval, coords,cscore, pdf]
                                [fimo_motif[spl[1]][ix].append(lst[ix]) for ix, (a, b) in \
                                        enumerate(zip(fimo_motif[spl[1]], lst))]
                            elif not centrimo_test:  
                                lst = [spl[0], spl[5], dist, qval, coords, 1 , pdf]
                                [fimo_motif[spl[1]][ix].append(lst[ix]) for ix, \
                                       (a, b) in enumerate(zip(fimo_motif[spl[1]], lst))]
                        else:
                            if  centrimo_test and spl[0] in centrimo:
                                cscore = centrimo[spl[0]][1]
                                fimo_motif[spl[1]] = [[spl[0]], [spl[5]], 
                                    [dist], [qval], [coords], [cscore], [pdf]]
                            elif not centrimo_test: 
                                fimo_motif[spl[1]] = [[spl[0]], [spl[5]], 
                                    [dist], [qval], [coords], [1] , [pdf]]
    return (fimo_motif,fimo_score)

# fimo_manual = '/export/home/dnyansagar/Bra/motif_analysis/Spur5/fimo.tsv'
#def web_fimo(fimo_file, centrimo, centrimo_test =True): # Approach
#print(cm)
fm = web_fimo(fimo_manual, cm, centrimo_test=True)[0]

#fm = fimo_parse(path, cm, centrimo_test=True)[0]
fms =  np.array(fimo_parse(path, cm, centrimo_test=True)[1], dtype=float)
q25 = np.quantile(fms, 0.25)

"""
plt.hist(fms, normed=True, bins=100)
plt.ylabel('FDR adjusted p-value')
print(q25)
q75 = np.quantile(fms, 0.75)
print(q75)
print(min(fms),max(fms))
"""
def get_counts(fimo, meme_tomtom, dreme_tomtom, centrimo): # Approach 1

    ''' Combine Fimo Centrimo results and select one motif per sequence '''
    count = dict()
    distance = dict()
    score = dict()
    oth = merge_dicts(meme_tomtom,dreme_tomtom)
    cooc = coOccurrence(d, fimo, meme_tomtom,  dreme_tomtom, motif_names)
    for k, v in fimo.items():
        f_name, f_score, f_distance, f_qval, f_coords, f_cscore , f_pdf = [
                v[0], v[1], v[2], v[3], v[4], v[5], v[6]]
        selection = select('combined', q25, 
                [f_name, f_score, f_distance, f_qval, f_coords, f_cscore, f_pdf]) 
        if selection != None: 
            score[k] = selection[1]
            distance[k] = selection[2]
            if selection[0].startswith('MA') and re.search(r'\d+$', selection[0]) != None:
                mot = motif_names[selection[0]]
                infostr1 = "\t".join(str(e) for e in \
                [k, mot, selection[1],  selection[2],  selection[3]])
                if mot in count: count[mot] +=1
                else: count[mot] =1
            else:
                if selection[0] in oth:
                    mot = motif_names[oth[selection[0]][0]]
                    infostr2 = "\t".join(str(e) for e in \
                    [k,  mot,  selection[1],  selection[2], selection[3]])
                    if mot in count: count[mot] +=1
                    else: count[mot] =1
                else:
                    stri =  'MEME' if selection[0].isdigit() else selection[0]
                    mot = stri+selection[0]
                    infostr3 = "\t".join( str(e) for e in \
                    [ k,  mot,  selection[1],  selection[2], selection[3]])
                    if mot in count: count[mot] +=1
                    else: count[mot] =1
    return (count,distance,score, cooc)

def count_distance(peak, dictionary, lst_of_lst): # Approach 2

    temp_dict1 = {'bHLH':0, 'hmg':0, 'pax':0, 'fox':0, 'homeo':0,'tbox1':0,'tbox2':0, 'other':0}
    temp_dict2 = {'bHLH':'NA', 'hmg':'NA', 'pax':'NA', 'fox':'NA', 'homeo':'NA','tbox1':'NA','tbox2':'NA'}
    dict_order = ('bHLH', 'hmg', 'pax', 'fox', 'homeo','tbox1','tbox2') # Fix order of the dictionary
    for i in lst_of_lst[0]:
        names = [name for name, seq in dictionary.items() if i in seq]
        if names !=[]:
            if temp_dict1[names[0]] ==1: 
                continue
            else: 
                temp_dict1[names[0]] +=1
                temp_dict2[names[0]] = lst_of_lst[1][lst_of_lst[0].index(i)]
        else:
            if temp_dict1['other'] ==1: 
                continue
            else: 
                temp_dict1['other'] +=1
    str_lst = [peak] + [str(temp_dict1[i]) for i in dict_order] + [str(temp_dict2[i]) for i in dict_order]
    #str_lst = [peak] + [str(i) for i in temp_dict1.values()] + [str(j) for j in temp_dict2.values()]
    ##print('\t'.join(str_lst))
    return('\t'.join(str_lst+['\n']))

def sp_case1(family_dict, list_of_list): # Approach 2 Specific
	''' This is the special case where we want to check 
	if the sox and T-box motifs have overlapping ranges'''
	mot_names = list_of_list[0]
	mot_dists = list_of_list[2]
	t1dist = "NA"; t2dist = "NA"; sox_dist = "NA"
	tonly = "NA"; sox_only = "NA"	
	names = [[name for name, seq in family_dict.items() if i in seq] for i in mot_names]
	for ix, i in enumerate(names):
		if i !=[]:
			if (i[0] =="tbox1" or i[0] =="tbox2" or i[0] =="hmg"):
				if i[0] == "tbox2" : t2dist = mot_dists[ix]
				elif i[0] == "tbox1": t1dist = mot_dists[ix]
				elif (i[0] =="hmg"): sox_dist = mot_dists[ix]
	if t1dist =="NA" and t2dist == "NA" and sox_dist !="NA":
		sox_only = sox_dist
	elif sox_dist  == "NA" and (t1dist !="NA" or t2dist != "NA"):
		tonly = t2dist if t2dist != "NA" else t1dist
	return (t1dist, t2dist, sox_dist, tonly, sox_only)

tandsox = [[],[]] ; ton = [[],[]] ;soxon = [[],[]]
def feed_scatter(t_sox_results):
	t_sox =  t_sox_results
	#[nList, f_score, f_distance, f_qval, f_coords, f_cscore , f_pdf])
	if t_sox[0] !='NA' and t_sox[1] !='NA' and t_sox[2] !='NA': # t and sox  != 0
		tandsox[0].append(min([t_sox[0], t_sox[1]]))
		tandsox[1].append(t_sox[2])
	elif t_sox[0] =='NA' and t_sox[1] =='NA' and t_sox[2] !='NA': # t = 0  sox !=0
		soxon[0].append(0)
		soxon[1].append(t_sox[2])
	elif (t_sox[0] !='NA' or t_sox[1] !='NA') and t_sox[2] =='NA': # t !=0 sox =0
		ton[0].append(t_sox[1] if t_sox[1] != 'NA' else t_sox[0])
		ton[1].append(0)
	elif (t_sox[0] !='NA' or t_sox[1] !='NA') and t_sox[2] !='NA': # t !=0 sox =0
		tandsox[0].append(t_sox[1] if t_sox[1] != 'NA' else t_sox[0])
		tandsox[1].append(t_sox[2])

thalf = set(); tfull = set()
def feed_bra_venn(family_dict, fin, peak):
	ful = [val for val in fin if val in family_dict['tbox1']]
	hal = [val for val in fin if val in family_dict['tbox2']]
	if len(ful) !=0 and len(hal) !=0:
		tfull.add(peak)
		thalf.add(peak)
	elif  len(ful) !=0: tfull.add(peak)
	elif len(hal)!=0: thalf.add(peak)


def get_counts2(fimo, meme_tomtom, dreme_tomtom, centrimo, motNameDict, family_dict): # Approach 2 
	oth = merge_dicts(meme_tomtom,dreme_tomtom)
	myfile = open(plots+"Upset.txt", "w")
	t1_d = []; t2_d = [] ; hmg = [] ; tonly_d = []; sox_only_d = []
	table = []
	for k, v in fimo.items(): # Get names of motifs from Jaspar Ids
		nList = []
		for i in v[0]:
			if i in oth: 
				#print(motNameDict[oth[i][0]])
				nList.append(motNameDict[oth[i][0]])
			elif i in motNameDict: nList.append(motNameDict[i])
			else: nList.append(i + "other")
		f_name, f_score, f_distance, f_qval, f_coords, f_cscore , f_pdf = [
																v[0], v[1], v[2], v[3], v[4], v[5], v[6]]
		t_sox =  sp_case1(family, [nList, f_score, f_distance, f_qval, f_coords, f_cscore , f_pdf])
		feed_scatter(t_sox) #For scatter plot
		feed_bra_venn(family_dict,nList, k)
		t1_d.append(t_sox[0]); t2_d.append(t_sox[1])
		hmg.append(t_sox[2]);  tonly_d.append(t_sox[3]); sox_only_d.append(t_sox[4])
		if overlap(f_coords)[0] != []: # Peaks with Overlapping motifs 
			olist = []; nolist = []
			odist = []; nodist = []
			oqval = []; noqval = []
			for i in overlap(f_coords)[0]:
				ovr_ind0 = f_coords.index(list(i[0]))
				ovr_ind1 = f_coords.index(list(i[1]))
				olist.append((nList[ovr_ind0], nList[ovr_ind1]))
				odist.append((f_distance[ovr_ind0], f_distance[ovr_ind1]))
				oqval.append((f_qval[ovr_ind0], f_qval[ovr_ind1]))
			for j in overlap(f_coords)[1]:
				noovr_ind = f_coords.index(list(j))
				nolist.append(nList[noovr_ind])
				nodist.append(f_distance[noovr_ind])
				noqval.append(f_qval[noovr_ind])
				selected_overlap = overlap_select(family, [olist, odist, oqval])
				#feed_bra_venn(family_dict,selected_overlap[0]+nolist, k)

			lst1 = [selected_overlap[0]+nolist,  selected_overlap[1]+nodist, 
                                                            selected_overlap[2]+noqval]
			print(count_distance(k, family, lst1))
			table.append(count_distance(k, family, lst1))
		else: # Peaks with no overlapping motifs
			lst2 = [nList, f_distance, f_qval]
			#feed_bra_venn(family_dict, nList, k)
			table.append(count_distance(k,family, lst2))
			print(count_distance(k,family, lst2))
	
    # Write counting to file for Upset plot
	for line in table:
		myfile.write(line)
	return(t1_d, t2_d,  hmg, tonly_d, sox_only_d)

######################
#### Print Counts ####
######################
count, distanceFromCenter, selectedScore, co_occurence =  get_counts(fm , mt, dt, cm)
dd = OrderedDict(sorted(count.items(), key=lambda x: x[1]))
merged_count = merge_family_motif_counts(count,d)
t1, t2, sox, tonly, sonly = get_counts2(fm , mt, dt, cm, motif_names, family)
ndict = merge_family_motif_counts(dd,d)

###############################
# Q-value distance density plot 
###############################

f0 = plt.figure(0)
plt.figure(figsize=(8,8))
dst2 = sns.distplot(fms,
            hist = False,
            kde = True,
            kde_kws = {'linewidth': 2},
            label = 'q-values')
plt.xlabel("Motif q-value")
save(path+'plots/fig0_qvalue_density_plot', ext='png', close=True, verbose=True)

########################
########################

#######################
###### Pie charts######
#######################

# Pie Chart 1
labels = dd.keys()
sizes = dd.values()
explode =[]
for i in dd.values():
    if i < 10: explode.append(0.2)
    else: explode.append(0)
#print(sum(dd.values()))
f1 = plt.figure(1)
plt.figure(figsize=(8,8))
plt.pie(sizes, labels=labels, 
        explode = explode, autopct='%1.1f%%',
        shadow=True, startangle=90)
plt.axis('equal')
plt.title("Motif pie chart")
save(path+'plots/fig1_piechart1'+datestring, ext='png', close=True, verbose=True)

# Pie Chart 2
f2 = plt.figure(2)
plt.figure(figsize=(8,8))
labels2 = ndict.keys()
sizes2 = ndict.values()
explode = list(itertools.repeat(0.05, len(labels2)))
plt.pie(sizes2, labels=labels2,
        autopct='%1.1f%%',
        pctdistance=0.8, 
        explode = explode,
        labeldistance=1.1,
        shadow=True, startangle=90)
centre_circle = plt.Circle((0,0),0.70,fc='white')
fig = plt.gcf()
fig.gca().add_artist(centre_circle)
plt.axis('equal')
plt.title("Motif pie chart")
save(path+'plots/fig2_piechart2'+datestring, ext='png', close=True, verbose=True)

#######################
###### Bar plot ######
#######################
val = np.array(list(sizes), dtype='float64')/sum(sizes)*100
f3 = plt.figure(3)
plt.figure(figsize=(8,8))
plt.barh(list(labels), list(val))
plt.xlabel("counts (%)")
plt.title("Motif bar chart")
save(path+'plots/fig3_barplot'+datestring, ext='png', close=True, verbose=True)

#######################
###### Histogram ######
#######################
# distanceFromCenter = get_counts(fm , mt, dt)[1].values()
f4 = plt.figure(4)
plt.figure(figsize=(8,8))
plt.hist(distanceFromCenter.values(), bins=20)
plt.xlabel("distance")
plt.title("Distance of motif from center of sequence")
save(path+'plots/fig4_histogram'+datestring, ext='png', close=True, verbose=True)

#######################
###### Box  plot ######
#######################
# selectedScore = get_counts(fm , mt, dt)[2].values()
ss = np.array(list(selectedScore.values()), dtype=float)
f5 = plt.figure(5)
plt.figure(figsize=(8,8))
bp0 = plt.boxplot(np.array(list(selectedScore.values())).astype(np.float),
        patch_artist=True, notch=True)
for element in ['boxes', 'whiskers', 'fliers',\
        'means', 'medians', 'caps']:
        plt.setp(bp0[element], color='red')
for patch in bp0['boxes']:
        patch.set(facecolor='grey') 
plt.xlabel("Fimo score")
plt.title("Score of the selected motif")
save(path+'plots/fig5_distance_boxplot'+datestring, ext='png', close=True, verbose=True)

#######################
# Distance density plot
#######################
def plotDensity(lists, names):
    plt.figure(figsize=(13,13))
    for i, j in zip(lists, names):
        dst = sns.distplot([int(x) for x in i if x !='NA'], 
            hist = False, 
            kde = True,
            kde_kws = {'linewidth': 7},
            label = j)#.set(xlim=(0))
        plt.title("Bra and Sox motif distribution",fontsize=40)
        plt.xlabel('Distance from summit,(bp)',fontsize=40)
        plt.tick_params(axis='both',labelsize=40)
        plt.legend(fontsize=40)
        #plt.suptitle("Bra and Sox motif distribution",fontsize=40)
        #dst.axes.set_title("Bra and Sox motif distribution",fontsize=40)
        #plt.set_xlabel("Distance from summit",fontsize=40)
        #plt.tick_params(labelsize=20)
        #plt.setp(dst.get_legend().get_texts(), fontsize=40)
    save(path+'plots/fig7_distanceDensityPlot'+datestring, ext='svg', close=False, verbose=True, transparent=True)
density_plot = plotDensity([t1,t2, sox, tonly, sonly], 
                ['Bra-half', 'Bra-Full', 'Sox', 'Bra only','Sox only'])

# new_dict = {'Bra-half': t1, 'Bra-Full':t2}
# venn(new_dict, )
#######################
##### venn diagram#####
#######################
def plotVenn(dictionary):
    dict2 = {k: v for k, v in dictionary.items() if len(v) > 10}
    if 'bra-half' in dict2:
        all = dict2['bra-half'].union(dict2['bra-full'])
        dict2['brachyury'] = all
        del dict2['bra-half']
        del dict2['bra-full']
        plt.figure(figsize=(8,8))
        venn(dict2, 
            #fmt="{percentage:.1f}%", 
            fontsize=16, 
            legend_loc="upper left")
        plt.title("Shared motifs in the peaks")
        save(path+'plots/fig6_intersection_venn'+datestring, ext='png', close=True, verbose=True)
    else:
        plt.figure(figsize=(8,8))
        venn(dict2, 
            #fmt="{percentage:.1f}%", 
            fontsize=16, 
            legend_loc="upper left")
        plt.title("Shared motifs in the peaks")
        save(path+'plots/fig6_intersection_venn'+datestring, ext='png', close=True, verbose=True)
#v = plotVenn(co_occurence)
f9 = plt.figure(0)
plt.figure(figsize=(10,10))
bra_venn = venn2([thalf, tfull], set_labels = ('Half palindrome', 'Full palindrome'))
plt.title('Brachyury motifs\n', fontsize=40)
for text in bra_venn.set_labels:
    text.set_fontsize(40)
for text in bra_venn.subset_labels:
    text.set_fontsize(40)
save(path+'plots/fig9_bra_venn', ext='svg', close=True, verbose=True)

data = (tandsox, ton, soxon)
colors = ("red", "green", "blue")
groups = ("Bra and Sox", "Bra only", "Sox Only")

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(1, 1, 1, facecolor="1.0")
for data, color, group in zip(data, colors, groups):
	x, y = data
	ax.scatter(x, y, alpha=0.8, c=color, edgecolors='none', s=30, label=group)
	plt.title('Matplot scatter plot')
	#plt.legend(loc=2)
save(path+'plots/fig8_distance_scatter'+datestring, ext='png', close=True, verbose=True)

