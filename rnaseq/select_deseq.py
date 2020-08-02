#!/usr/bin/python 
import re
import xlrd
import pandas
from datetime  import datetime
datestring = datetime.strftime(datetime.now(), '%Y-%m-%d')

xtr = '/home/rohit/Projects/2/Xtr/Differentially_expressed_mmc2_tab.xlsx'
spu = '/home/rohit/Projects/2/Spu/seaurchinGeni_diff_exp_all_with_ref_id.xlsx'
mmu = '/home/rohit/Projects/2/Mmu/Differentially_expressed_mouse_genes_default.xlsx'
nve = '/home/rohit/Projects/2/Nve/bra_travsctl.xlsx'

spudf = pandas.read_excel(spu)
spu_filt = spudf[(abs(spudf.Log2foldchange) >0.5) & (abs(spudf.pval) < 0.05)]
print '%s\t%d\t%s\t%d'%('Spu Original ',len(spudf[spudf.columns[0]]),\
												'Filtered', len(spu_filt[spu_filt.columns[0]]))
xtrdf = pandas.read_excel(xtr)
xtr_filt = xtrdf[(abs(xtrdf.log2foldchange) >0.5) & (abs(xtrdf.pvalue) < 0.05)]
xtr_filt.loc[:,'Ensembl_ID_new'] = xtr_filt['Ensembl_ID'].str.split('|').str[0].fillna('')
#print xtr_filt['Ensembl_ID_new'].str[0]
print '%s\t%d\t%s\t%d'%('Xtr Original ',len(xtrdf[xtrdf.columns[0]]),\
												'Filtered', len(xtr_filt[xtr_filt.columns[0]]))
mmudf = pandas.read_excel(mmu)
mmu_filt = mmudf[(abs(mmudf.Log2foldchange) >0.5) & (abs(mmudf.pval) < 0.05)]
print '%s\t%d\t%s\t%d'%('Mmu Original ',len(mmudf[mmudf.columns[0]]),\
												'Filtered',  len(mmu_filt[mmu_filt.columns[0]]))
nvedf = pandas.read_excel(nve)
nve_filt = nvedf[((nvedf.log2FoldChange >0.5) | (nvedf.log2FoldChange < -0.5)) & (nvedf.pvalue < 0.05)]
print '%s\t%d\t%s\t%d'%('Nve Original ',len(nvedf[nvedf.columns[0]]),\
												'Filtered',  len(nve_filt[nve_filt.columns[0]]))

nve_filt.to_csv('nve_de_filtered'+datestring+'.csv', header=False,\
								sep='\t', index=False, columns=('Id', 'log2FoldChange'))
mmu_filt.to_csv('mmu_de_filtered'+datestring+'.csv', header=False,\
								sep='\t', index=False, columns=('Gene_id', 'Log2foldchange'))
xtr_filt.to_csv('xtr_de_filtered'+datestring+'.csv', header=False,\
								sep='\t', index=False, columns=('Ensembl_ID_new', 'log2foldchange'))
spu_filt.to_csv('spu_de_filtered'+datestring+'.csv', header=False,\
								sep='\t', index=False, columns=('Tabella_input_DESeq_RC_cutoff_I_quartile_Gene', 'Log2foldchange'))

