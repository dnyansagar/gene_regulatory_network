
# Mouse genes

## Single exon sharing approach
We want to get a stable list of transcription start sites and a reliable Gene groupings. The first approach was to take genes from the UCSC Ensembl genes which share at least one exon boundary start and end. Preferably it would have been internal coding exons only but many genes had no coding sequence annotated. This was accomplished with an earlier version of the `merge_isoforms.py` script. 


## Ensembl dataset merging approach
We would prefer to use Ensembl data as it is well curated, but it was initially incompatible with Michi's scripts. In order to get around this problem, we merge the data with the UCSC BED file (which is compatible with Michi's scripts), and extract information about TSSs and exons using a newer version of the `merge_isoforms.py` script.

Downloaded a BED file from the UCSC table browser using mouse genome mm10 using the track "GENCODE VM11 (Ensembl 86)", table "Comprehensive (wgEncodeGencodeCompVM11)" outputting the whole gene for BED into the file `Mouse_mm10_ensembl_comprehensive.bed`.

Downloaded file `gene_to_tx.biomart.txt` from the BioMart website constraining all mm10 genes to type `protein_coding`. Included the features "Gene Stable ID", "Transcript Stable ID" and "Transcription Start Site". These were crossreferenced to the file Rohit gave me (UCSC Ensembl genes mm10) with `Mouse_mm10_ensembl_comprehensive.bed`:


```bash
txids = set(l.split()[0] for l in open('tx_to_tss.biomart.txt'))
w = open('mm10_ensembl_compr.protein_coding.bed','w')
for line in open('mm10_ensembl_compr.bed'):
    # note we forgo version numbers in our comparison
    if line.split()[3].split('.')[0] in txids: 
        w.write(line)
```

Some genes are lost in the process, specifically, there are 55269 IDs in the biomart set but when merged with the UCSC set it is reduced to 52266.

First sort the genes so they can be viewed in IGV and compared downstream with bedtools:


```bash
bedtools sort -i mm10_ensembl_compr.protein_coding.bed \
    > mm10_ensembl_compr.protein_coding.srt.bed
```

Now run `merge_isoforms.py`:


```bash
python /proj/rpz/src/merge_isoforms/merge_isoforms.py \
    --ignore-versions \
    --gene-mapping-file gene_to_tx.biomart.txt \
    --tss-file tx_to_tss.biomart.txt \
    --supergene-bed-out mm10_ensembl_compr.protein_coding.supergene.bed \
    --upstream-tss-bed-out mm10_ensembl_compr.protein_coding.upstream_tss.bed \
    --consensus-tss-bed-out mm10_ensembl_compr.protein_coding.consensus_tss.bed \
    --consensus-tss-stats-out mm10_ensembl_compr.protein_coding.consensus_tss.stats.txt \
    mm10_ensembl_compr.protein_coding.srt.bed \
    > mm10_ensembl_compr.protein_coding.representative.bed
```

The current version of `merge_isoforms.py` does not keep the output sorted. For now, leaving it. It creates a "supergene" BED file which contains all the exons of the genes (in this case, as defined by the biomart download of the genes from Ensembl.) The consensus TSS is defined as the most common TSS among the transcripts in the gene defaulting to the most upstream site in cases of a tie. The upstream approach and consensus approach differed in 1821 out of 21933 cases.  Often these were small 5' UTR exons in a single isoform which made the TSS jump far ahead of the consensus.

Now we will move ahead with the following:

1. filter out low scoring peaks (Michi's script)
2. filter exonic peaks using `mm10_ensembl_compr.protein_coding.supergene.bed` (or equivalently `mm10_ensembl_compr.protein_coding.srt.bed`). Would suggest trying the following for eliminating exonic peaks:
```
bedtools intersect -v -a PEAK_FILE -b mm10_ensembl_compr.protein_coding.supergene.bed
``` 
3. and `mm10_ensembl_compr.protein_coding.consensus_tss.bed` for detecting the closest gene (`bedtools closest`).

