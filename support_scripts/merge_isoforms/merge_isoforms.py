#! /usr/bin/env python

import argparse
import collections
import itertools

class GeneGroup(object):
    all_groups = collections.OrderedDict()
    def __init__(self,gene):
        self.genes = set([gene])
        self.chrom = gene.chrom
        self.strand = gene.strand
        self.name = gene.name
        self.dirty = True
        GeneGroup.all_groups[self] = gene
    def add_gene(self,gene):
        self.genes.add(gene)
        gene.group = self
        self.dirty = True
    def join(self,group):
        if self is group:
            return
        for gene in group.genes:
            self.add_gene(gene)
        if group in GeneGroup.all_groups:
            del GeneGroup.all_groups[group]
        self.dirty = True

    def _update(self):
        if self.dirty:
            self.make_supergene()
            self.find_representative()
            self.find_consensus_ends()
            self.find_upstream_tss()

    def find_consensus_ends(self):
        self._tss_counter = collections.Counter(
            gene.fivep for gene in self.genes)
        consensus,consensus_count = self._tss_counter.most_common(1)[0]
        self._tss_stat = '%s,%d,%d,%.2f' % (self.name,consensus_count,len(self.genes),
                                         float(consensus_count)/len(self.genes))
        cmpfun = min if self.strand == '+' else max
        self._consensus_tss = cmpfun(
            c[0] for c in self._tss_counter.most_common() 
            if c[1] == consensus_count)

        self._tts_counter = collections.Counter(
            gene.threep for gene in self.genes)
        consensus,consensus_count = self._tts_counter.most_common(1)[0]
        self._tts_stat = '%s,%d,%d,%.2f' % (self.name,consensus_count,len(self.genes),
                                         float(consensus_count)/len(self.genes))
        cmpfun = max if self.strand == '+' else min
        self._consensus_tts = cmpfun(
            c[0] for c in self._tts_counter.most_common() 
            if c[1] == consensus_count)



    def find_upstream_tss(self):
        if self.strand == '+':
            self._upstream_tss = min(gene.start for gene in self.genes)
        else:
            self._upstream_tss = max(gene.end for gene in self.genes)

    def make_supergene(self):
        '''
        algorithm to merge all segments: 
            create a stack with merged exons
            sort all exons by start point and add the first exon to it
            iterating through all exons,
                if the exon overlaps the current exon on the stack, merge the
                exons
                otherwise, add the exon to the stack
        '''
        all_exons = collections.deque(sorted(itertools.chain.from_iterable(
            [gene.exons for gene in self.genes]),
            cmp=lambda a,b: cmp(a[0],b[0])))
        merged = [all_exons.popleft()]
        start = merged[0][0]
        end = merged[0][0]
        while len(all_exons):
            exon = all_exons.popleft()
            end = max(exon[1],end)
            if exon[0] < merged[-1][1]:
                merged[-1] = (merged[-1][0],exon[1])
            else:
                merged.append(exon)

        self._supergene = Gene(
            self.chrom,start,end,
            self.name+'Supergene','0',self.strand,
            start,end,'0',len(merged),
            ','.join([str(exon[1]-exon[0]) for exon in merged]),
            ','.join([str(exon[0]-start) for exon in merged]))

    def find_representative(self):
        rep = None
        for gene in self.genes:
            if rep is None or gene.length > rep.length:
                rep = gene
            #if self.name == 'ENSMUSG00000079434':
                #print gene.name, gene.length, rep.name
        self._representative = rep


    @property
    def start(self):
        return self.supergene.start
        
    @property
    def end(self):
        return self.supergene.end

    @property
    def supergene(self):
        self._update()
        return self._supergene

    @property
    def representative(self):
        self._update()
        return self._representative
    @property
    def upstream_tss(self):
        self._update()
        return self._upstream_tss
    @property
    def consensus_tss(self):
        self._update()
        return self._consensus_tss
    def _tss_to_bed(self,tss):
        return '%s\t%d\t%d\t%s\t0\t%s' % (
            self.chrom,tss,tss,self.name,self.strand)
    @property
    def consensus_tts(self):
        self._update()
        return self._consensus_tts

    @property
    def upstream_tss_bed(self):
        self._update()
        return self._tss_to_bed(self._upstream_tss)
    @property
    def consensus_tss(self):
        self._update()
        return self._consensus_tss
    @property
    def consensus_tss_bed(self):
        self._update()
        return self._tss_to_bed(self._consensus_tss)
    @property
    def tss_statistics(self):
        self._update()
        return self._tss_stat

    def __str__(self):
        rep = self.representative
        return '%d %s: %s' % (len(rep.cds_exons), rep.name,
                          ','.join(gene.name for gene in self.genes))

class Gene(object):
    @classmethod
    def from_psl(cls,matches,misMatches,repMatches,nCount,qNumInsert,
                 qBaseInsert,tNumInsert,tBaseInsert,strand,
                 name,qSize,qStart,qEnd,chrom,tSize,start,end,block_count,
                 block_sizes,qStarts,block_starts):
        corr_block_starts = ','.join([str(int(s)-int(start)) 
                                      for s in block_starts.split(',') if len(s)])
        return cls(chrom,start,end,name,matches,strand,0,0,0,
                   block_count,block_sizes,corr_block_starts)

    def __init__(self,chrom,start,end,name,score,strand,cds_start,cds_end, 
                 item_rgb,block_count,block_sizes,block_starts,*args):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.name = name
        self.score = int(score)
        self.strand = strand
        self.cds_start = int(cds_start)
        self.cds_end = int(cds_end)
        self.item_rgb = int(item_rgb)
        self.block_count = int(block_count)
        self.block_sizes = [int(s) for s in block_sizes.split(',') if len(s)]
        self.length = sum(self.block_sizes)
        self.block_starts = [int(s) for s in block_starts.split(',') if len(s)]

        self.blast = ''
        self.group = None

        self.exons = []
        self.cds_exons = []

        for start,size in zip(self.block_starts,self.block_sizes):
            start += self.start
            end = start+size
            sc = sorted([start,end,self.cds_start,self.cds_end])
            self.exons.append((start,end))
            if (sc[0] == start and sc[1] == end) or (
                sc[0] == self.cds_start and sc[1] == self.cds_end):
                continue
            self.cds_exons.append((max(self.cds_start,start),min(self.cds_end,end)))

        #if name.startswith('WHL22.628628'):
            #print 'HI',name,self.block_starts,self.block_sizes,self.exons
        if strand == '+':
            self.fivep = self.start
            self.cds_fivep = self.cds_start
            self.threep = self.end
            self.cds_threep = self.cds_end
        else:
            self.fivep = self.end
            self.cds_fivep = self.cds_end
            self.threep = self.start
            self.cds_threep = self.cds_start

    def set_tss(self,tss):
        self.fivep = tss
        self.threep = tss

    def is_isoform(self,other,check=False):
        #optimization 
        #print self.name, other.name
        if self.start > other.end or other.start > self.end:
            return False
        if self.chrom != other.chrom:
            return False
        if self.strand != other.strand:
            return False
        #print self.exons
        #print other.exons
        for a,b in itertools.product(self.exons,other.exons):
            if a[0] == b[0] and a[1] == b[1]:
                return True
        return False
    def __str__(self):
        return '\t'.join(str(item) for item in 
                         [self.chrom,self.start,self.end,self.name,self.score,self.strand,
                          self.cds_start,self.cds_end,self.item_rgb,self.block_count,
                          ','.join(str(s) for s in self.block_sizes),
                          ','.join(str(s) for s in self.block_starts)])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='''Summarize genes in a bed file by picking the longest 
        isoform of any given gene provided the genes are on the same strand 
        and share at least one exon''')
    parser.add_argument('BED_FILE')
    parser.add_argument('--ignore-versions', action='store_true',
                        help='''Ignore version numbers (e.g.  ENSMUST00000086465.5 
                        is compared as ENSMUST00000086465)''')
    #XXX: if mapping file does not contain a transcript/gene, it will not be output
    parser.add_argument('--gene-mapping-file', metavar='TSV',
                        help='''A file (perhaps downloaded from biomart), whose
                        first column contains transcript IDs and second gene IDs,
                        to be used as the method of grouping genes instead of
                        the shared-exon method in the description here.''')
    parser.add_argument('--tss-file', metavar='TSV',
                        help='''A file (perhaps downloaded from biomart), whose
                        first column contains Transcript IDs and second
                        Transcription Start Sites.''')
    parser.add_argument('--supergene-bed-out', metavar='BED',
                        help='''output a "supergene" for each gene with all 
                        exons merged, so as to detect all possible exon 
                        overlaps''')
    parser.add_argument('--upstream-tss-bed-out', metavar='BED',
                        help='''output the 5' most upstream TSS among all a 
                        gene's transcripts for each gene to BED''')
    parser.add_argument('--consensus-tss-bed-out', metavar='BED',
                        help='''output the most common among all a gene's 
                        transcripts for each gene to BED''')
    parser.add_argument('--consensus-tss-stats-out', metavar='FILE', 
                        help='''output the statistics about the most common TSS 
                        for each gene to FILE''')
    parser.add_argument('--output-gene-mapping', metavar='OUT', help='''output 
                        gene to transcript id TSV''')
    args = parser.parse_args()

    all_genes = collections.defaultdict(list)
    tx_dict = {}
    for l in open(args.BED_FILE):
        fields = l.split()
        gene = Gene(*fields)
        all_genes[fields[0]].append(gene)
        txid = fields[3].split('.')[0] if args.ignore_versions else fields[3]
        tx_dict[txid] = gene

    if args.gene_mapping_file is not None:
        gene_groups = {}
        for line in open(args.gene_mapping_file):
            tx,gene = line.strip().split()[:2]
            if tx not in tx_dict:
                continue
            if gene not in gene_groups:
                gene_groups[gene] = GeneGroup(tx_dict[tx])
                gene_groups[gene].name = gene
            else:
                gene_groups[gene].add_gene(tx_dict[tx])
    else:
        for chr,genes in all_genes.items():
            for i,geneA in enumerate(genes):
                if geneA.group is None:
                    geneA.group = GeneGroup(geneA)
                    for geneB in genes[i+1:]:
                        if geneA.is_isoform(geneB):
                            geneA.group.add_gene(geneB)

    if args.tss_file is not None:
        for line in open(args.tss_file):
            tx,tss = line.strip().split()[:2]
            if not tx in tx_dict:
                continue
            tss = int(tss)
            tx_dict[tx].set_tss(tss)

    if args.output_gene_mapping is not None:
        with open(args.output_gene_mapping,'w') as out:
            for group in GeneGroup.all_groups.keys():
                for gene in group.genes:
                    out.write('\t'.join([group.name, gene.name]) + '\n')


    for group in GeneGroup.all_groups.keys():
        print group.representative

    def output_group_property(property,filename=None):
        if filename is not None:
            open(filename,'w').write('\n'.join(
                str(getattr(group,property)) 
                for group in GeneGroup.all_groups.keys()))

    output_group_property('supergene',args.supergene_bed_out)
    output_group_property('upstream_tss_bed',args.upstream_tss_bed_out)
    output_group_property('consensus_tss_bed',args.consensus_tss_bed_out)
    output_group_property('tss_statistics',args.consensus_tss_stats_out)
