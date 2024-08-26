import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import gseapy
from matplotlib_venn import venn2_unweighted, venn2_circles 


#################### Venn diagrams ####################

rank_list = pd.read_csv('./H_vs_N_DEG_table.csv') # Read DEG table for Hypoxia vs Normoxia conditions
rank_list = rank_list.sort_values('log2FoldChange', ascending=False)
rank_list.index = rank_list.iloc[:,0]
rank_list = rank_list.dropna()

FOSL2_module = ['VIM','NPC2','LAPTM5','CTSD','HMOX1','MCL1','CD163','CXCL8','PLIN2','FOSB','CXCR4','S100A10',
                'S100A6','ATP1B3','TFRC','IL1B','STAB1','ANXA1','SQSTM1','TGFBI','CAST','ARL4C','IQGAP1','ANXA2',
                'RGCC','SERPINE1','SLC16A10','VEGFA','CD44','LY6E','HK2','METRNL','ATP6V1G1','CD9','MDM2','FOSL2',
                'SOCS3','CD93','THBD','NRP1','NDRG1','ANKH','UPP1','LMNA','DHRS3','ATXN1','MYO5A','RPN1','MIR24-2',
                'VKORC1','ASPH','FAM20C','PDE8A','SLC6A6']

print(len(FOSL2_module))
print(len(rank_list[(rank_list.log2FoldChange >= .5) & (rank_list.pvalue <= .05)].index))
print(len(set(rank_list[(rank_list.log2FoldChange >= .5) & (rank_list.pvalue <= .05)].index).intersection(FOSL2_module)))
#54
#686
#14

font2 = {'family': 'Arial', 'size': 15} # use for labels
plt.rc('font', **font2) # sets the default font 

# depict venn diagram 
venn2_unweighted(subsets=(54-14, 686-14, 14),
      set_labels=('FOSL2 Module', 'Upregulated in\nHypoxia Only'),  
      set_colors=("orange", "blue", "red"), alpha=0.2,) 
  
# outline of circle line style and width 
venn2_circles(subsets=(1, 1, 1), linewidth=2) 

# title of the venn diagram 
plt.tight_layout()
plt.savefig('./venn_diag_FOSL2_H.pdf',    # Set path and filename
            dpi = 600,                     # Set dots per inch
            bbox_inches="tight",           # Remove extra whitespace around plot
            facecolor='white')             # Set background color to white

## Repeat above steps for normoxia co-culture, and hypoxia co-culture

#################### GSEA Analysis for FOSL2 Module ####################

gmt = dict()
gmt['FOSL2_module'] = FOSL2_module

pre_res = gseapy.prerank(rnk = rank_list,
                     gene_sets= gmt,
                     threads=4,
                     min_size=5,
                     max_size=500,
                     permutation_num=1000, 
                     outdir=None,
                     seed=6,
                     verbose=True, 
                    )

terms = pre_res.res2d.Term

from gseapy import gseaplot

gseaplot(rank_metric=pre_res.ranking, term=terms[0], **pre_res.results[terms[0]], ofname='./H_vs_N_GSEAplot.pdf')

## Repeat above steps for normoxia co-culture, and hypoxia co-culture

