# gene_regulatory_network
Study gene regulatory network with a combination of ChIP-seq, RNA-seq, and comparative genomics
This repository is a collection of scripts used for the study of the transcription factor brachyury. 
The project overall is divided into three parts
<br />
**1. chip (ChIP-seq)** 
Where scripts used to 
  i. identify TF binding targets 
  ii. Selection between two closest genes to the target based on the ortholog information obtained using OMA Standalone program
  iii. Motif analysis which consists of two parts, first part where centrally enriched motifs are identified in the dataset with Centrimo program from MEME-ChIP suite 
       In second aspect we look at the other motifs in the peak reagion to get insights into possible cofactors to the TF with FIMO program from the MEME-ChIP suite
       The parameters for the MEME-ChIP and FIMO programs are listed in motif_parameters.txt file.
       The output of these two programs was parsed with assign_centrimo_motifs.py and web_fimo.py respectively. 
<br />
**2. rnaseq (RNA-seq)**
  i. For RNA-seq a wrapper script around the two most widely used DEG analysis packages i.e DESeq2 and edgeR called SARTools is used. 
     See repository of SARTools (https://github.com/PF2-pasteur-fr/SARTools) 
  ii. Additional filtering criterion for the genes that have already passed the alpha significance filter.

<br />

**3. Data processing**
In this section the ortholog data along with chip target data and DE genes information is compiled into a single table


![alt text](https://github.com/dnyansagar/gene_regulatory_network/blob/master/support_scripts/projectLayout.png?raw=true)


