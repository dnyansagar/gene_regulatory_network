#! /usr/bin/env python

import collections
import itertools
import csv
import docopt
import IntervalTree
import merge_isoforms


__doc__ = '''
Usage:  associate_genes.py [options] PEAKFILE GENEFILE [DESEQ ...]

This script first determines genes in the given GENEFILE and then associates the
peaks with the genes, finally reporting all transcripts in the gene. Transcripts
are grouped into genes if they share at least one exon's coordinates exactly.
TSS are determined by the most common TSS among all the transcripts in the gene.
If there is a tie, the most upstream TSS is considered.

Options:
    --format=FORMAT             format of gene file (one of either BED or PSL) [default: PSL]
    --max-distance=DIST         maximum distance to associate a gene with [default: 50000]
    --effective-width=WIDTH     width to consider when deciding if exonic [default: 20]
    --better-descrs=FILE        specify a 2-column tab-separated file with manual annotations
    --gene-anchor=ANCHOR        define the anchor of a gene to measure the
                                distance from. Options are TSS (only the
                                beginning of the transcript), 
                                STARTEND (TSS and TTS)  [default: TSS]
    --start-only                only consider the 5' end of genes when considering intergenic peaks
    --ignore-utr                consider the boundaries of coding exons rather
                                than the transcription start site in the gene
                                annotations.
    --blast-file=FILE           a tabular blast file that has descriptions
                                associated with the gene IDs 
    --bed-out=FILE              output to bed file
'''

args = docopt.docopt(__doc__)
max_distance = int(args['--max-distance'])
eff_width = int(args['--effective-width'])

if args['--format'] not in ['BED','PSL']:
    raise SystemExit('--format must be one of either BED or PSL')

peakfile = args['PEAKFILE']

manual_descrs = {}
if args['--better-descrs'] is not None:
    manual_descrs = dict([l.strip().split('\t') for l in open(args['--better-descrs'])])

class Peak(object):
    def __init__(self,chrom,start,end,summit,score,enrichment,*args):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.summit = int(summit)
        self.score = score
        self.enrichment = enrichment
        
        self.rest = list(args)

        self.closest_gene = ''
        self.multiple = False
        self.exonic = False

if args['--blast-file'] is not None:
    bf = open(args['--blast-file'])
    bf.next()
    blast = {l.split('\t')[1]: l.split('\t')[4] for l in bf}
else:
    blast = {}

pf = open(args['PEAKFILE'])
#new_header = pf.next().strip() + '\tMultiple\tExonic\tGene\tBlast Hit'
pf.next()
'''
for deseq_group in args['DESEQ']:
    new_header += '\t' + '\t'.join([deseq_group+'_'+type for type in
                             ['log2FoldChange','padj']])
                             '''

#print new_header
peaks = collections.defaultdict(list)
for line in pf:
    l = line.split()
    peaks[l[0]].append(Peak(*l))

all_genes = collections.defaultdict(list)
intervals = collections.defaultdict(list)
tx_dict = {}
gf = open(args['GENEFILE'])
for l in gf:
    fields = l.split()
    if args['--format'] == 'BED':
        gene = merge_isoforms.Gene(*fields)
    else:
        gene = merge_isoforms.Gene.from_psl(*fields)
    all_genes[gene.chrom].append(gene)
    txid = fields[3] #fields[3].split('.')[0] if args.ignore_versions else 
    tx_dict[txid] = gene
    intervals[gene.chrom] += [IntervalTree.Interval(s,e) for s,e in gene.exons]

it = {c: IntervalTree.IntervalTree(i) for c,i in intervals.items()}

for chr,genes in all_genes.items():
    for i,geneA in enumerate(genes):
        if geneA.group is None:
            #if geneA.name.startswith('WHL22.628628'):
                #print 'checking %s: %s' % (geneA.name,','.join(g.name for g in
                                                               #genes))
            geneA.group = merge_isoforms.GeneGroup(geneA)
            for geneB in genes[i+1:]:
                if geneA.is_isoform(geneB):
                    geneA.group.add_gene(geneB)


DeseqResult = collections.namedtuple('DeseqResult',[
    'name','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj'])
deseqs = collections.defaultdict(dict)
for deseq_tab in args['DESEQ']:
    f = open(deseq_tab)
    r = csv.reader(f,delimiter='\t')
    r.next()
    for row in r:
        deseqs[deseq_tab][row[0]] = DeseqResult(*row)

def overlap_or_dist(coord, tss, tts, strand,name):
	#print name, coord, tss, tts, strand, strand == '+', strand =='-', coord < tts, tss - coord, tts - coord
	if tss < coord < tts: return 0
	if strand=='+':
		if coord < tss: return coord-tss
		elif tss < coord < tts: return 0
		elif coord > tts: return coord-tts
	elif strand=='-':
		if coord < tts: return tts-coord
		elif tts < coord < tss: return 0
		elif coord > tss: return tss-coord

#def dist_to_tss(gene,peak): return peak.summit-gene.consensus_tss
def dist_to_tss(gene,peak):
	#print gene.name, peak.summit, gene.start, gene.end, gene.strand, gene.consensus_tss
	if gene.strand =='+':
		return peak.summit-gene.consensus_tss		
	elif gene.strand == '-':
		return gene.consensus_tss-peak.summit
def dist_to_gene(gene,peak): 
	return overlap_or_dist(peak.summit,gene.consensus_tss, gene.consensus_tts, gene.strand, gene.name)

if args['--gene-anchor'] == 'TSS':
    gene_dist = dist_to_tss
elif args['--gene-anchor'] == 'STARTEND':
    gene_dist = dist_to_gene
else:
    raise ValueError, '''Invalid --gene-anchor value: use onen of "TSS" or "STARTEND"'''

all_groups = collections.defaultdict(list)
for group in merge_isoforms.GeneGroup.all_groups:
    all_groups[group.chrom].append(group)

print '\t'.join(['seqnames','start','end','summit','score','enrichment',
                 'GeneID','seqnames','start','end','strand','distance_to_tss','distance_to_gene',
                 'GeneID','seqnames','start','end','strand','distance_to_tss','distance_to_gene',
                 'exonic'])
bed_file = None
if args['--bed-out']:
    bed_file = open(args['--bed-out'],'w')

def sign(i): return i/abs(i)
for chrom,p in peaks.items():
    for peak in p:
        track = False
        multiple = False
        cand_genes = []
        for gene in all_groups[chrom]:
            metric_dist = gene_dist(gene,peak)
            dist = dist_to_gene(gene,peak)

            if dist <= max_distance:
                cand_genes.append((gene,metric_dist,dist))


        exonic = chrom in it and len(it[chrom].search(peak.summit-eff_width/2,
                                                      peak.summit+eff_width/2)) > 0

        cand_genes = list(sorted(cand_genes,cmp=lambda a,b:
                                 cmp(abs(a[1]),abs(b[1]))))

		


        # this code is used if you want to output a bedfile with the associated
        # gene names
        name = 'X'
        row = [chrom,peak.start,peak.end,peak.summit,peak.score,peak.enrichment]
        if len(cand_genes) > 0:
            gene = cand_genes[0][0]
            name = gene.representative.name
            row += [gene.representative.name,gene.chrom,
                    gene.representative.start,gene.representative.end,
                    gene.strand,cand_genes[0][1],cand_genes[0][2]]
        else:
            row += ['NA','NA','NA','NA','NA','NA','NA']

        if len(cand_genes) > 1:
            gene = cand_genes[1][0]
            name += ','+gene.representative.name
            row += [gene.representative.name,gene.chrom,
                    gene.representative.start,gene.representative.end,
                    gene.strand,cand_genes[1][1],cand_genes[1][2]]
        else:
            row += ['NA','NA','NA','NA','NA','NA','NA']


        for cand_gene in cand_genes:
            if not any(d == 0 for d in cand_gene[1:]) and sign(cand_gene[1]) != sign(cand_gene[2]):
                raise ValueError, 'problem with peak'
        '''
        for deseq_group in args['DESEQ']:
            if closest != nogene and closest.name in deseqs[deseq_group]:
                new_row += [deseqs[deseq_group][closest.name].log2FoldChange,
                            deseqs[deseq_group][closest.name].padj]
            else:
                new_row += ['','']
                '''

        print '\t'.join(str(x) for x in row)
        if bed_file:
            bed_file.write(
                '\t'.join(str(x) for x in [chrom,peak.start,peak.end,name]) 
                + '\n')
