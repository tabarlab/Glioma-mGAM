import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import scipy.stats as stats
import seaborn as sns

#################### CPM Heatmap ####################
# Read cpm table averaged by condition from EdgeR output for all bulk RNA seq samples
cpm_table = pd.read_csv('./cpm_table_bulkRNA_grouped.csv')
cpm_table = cpm_table.dropna()
cpm_table.index = cpm_table.hgnc_symbol
cpm_table = cpm_table.iloc[:, 2:10]

log_cpm = np.log2(cpm_table + 1)

# Calculate zscore per row 
cpm_row_zscore = stats.zscore(log_cpm, axis=1)

gene_list = ['IL1B', 'VEGFA', 'SERPINE1', 'EMP2', 'HK2', 'ANGPTL4', 'HMOX1', 'ADAM8', 'NDRG1', 'THBS1', 
             'CXCL8', 'CCL20', 'CD44', 'CXCL3', 'CXCL2', 'ADAM8', 'CXCR4', 'LY6E', 'ANXA1', 'TFRC', 'CD9', 
             'ADM', 'STAB1', 'MARCO', 'CD151', 'PCOLCE2', 'P4HA2', 'ACTN1', 'CAPN2', 'PLOD2', 'FN1', 'VIM', 'CTSD', 'FOSL2']

heatmap_df  = cpm_row_zscore.loc[[a for a in gene_list if a in cpm_row_zscore.index]]
heatmap_df.head()

#### Heatmap plot script ####

fig, axs = plt.subplots(6, 1, figsize=(3,14), height_ratios=[10, 8, 6, 9, 1, .5], constrained_layout=True)

g = sns.heatmap(heatmap_df.iloc[:10,:], cmap='bwr', ax=axs[0], cbar=False, linewidths=1, linecolor='#46494c')
axs[0].set_xlabel(None)
axs[0].set_ylabel(None)
axs[0].xaxis.set_tick_params(labeltop=True,      
                         labelbottom=False,  
                         bottom=False,      
                         labelsize=12)  
# Reformat y-axis tick labels
axs[0].yaxis.set_tick_params(labelsize=12,     
                            bottom=False, labelleft = False, labelright = True, rotation = 0)  
bottom, top = axs[0].get_ylim()
axs[0].set_ylim(bottom + 0.1, top)
left, right = axs[0].get_xlim()
axs[0].set_xlim(left, right + 0.1)

g2 = sns.heatmap(heatmap_df.iloc[10:18,:], cmap='bwr', ax=axs[1], cbar=False, linewidths=1, linecolor='#46494c')
axs[1].set_xlabel(None)
axs[1].set_ylabel(None)
axs[1].set_xticklabels([], ha='center')
axs[1].xaxis.set_tick_params(labeltop=False,     
                         labelbottom=False, 
                         bottom=False,      
                         labelsize=12)  

# Reformat y-axis tick labels
axs[1].yaxis.set_tick_params(labelsize=12,     
                            bottom=False, labelleft = False, labelright = True, rotation = 0)   
bottom, top = axs[1].get_ylim()
axs[1].set_ylim(bottom + 0.1, top)
left, right = axs[1].get_xlim()
axs[1].set_xlim(left, right + 0.1)

kws = dict(cbar_kws=dict(ticks=[-1.5, -1, 0, 1, 1.5], orientation='horizontal'))
g3 = sns.heatmap(heatmap_df.iloc[18:24,:], cmap='bwr', ax=axs[2], cbar=False, linewidths=1, linecolor='#46494c')
axs[2].set_xlabel(None)
axs[2].set_ylabel(None)

axs[2].xaxis.set_tick_params(labeltop=False,     
                         labelbottom=False, 
                         bottom=False,    
                         labelsize=12)  

# Reformat y-axis tick labels
axs[2].yaxis.set_tick_params(labelsize=12,     
                            bottom=False, labelleft = False, labelright = True, rotation = 0)   
bottom, top = axs[2].get_ylim()
axs[2].set_ylim(bottom + 0.1, top)
left, right = axs[2].get_xlim()
axs[2].set_xlim(left, right + 0.1)


g4 = sns.heatmap(heatmap_df.iloc[24:33,:], cmap='bwr', ax=axs[3], cbar=False, linewidths=1, linecolor='#46494c')
axs[3].set_xlabel(None)
axs[3].set_ylabel(None)
axs[3].xaxis.set_tick_params(labeltop=False,    
                         labelbottom=False, 
                         bottom=False,    
                         labelsize=12)  

# Reformat y-axis tick labels
axs[3].yaxis.set_tick_params(labelsize=12,     
                            bottom=False, labelleft = False, labelright = True, rotation = 0)  
bottom, top = axs[3].get_ylim()
axs[3].set_ylim(bottom + 0.1, top)
left, right = axs[3].get_xlim()
axs[3].set_xlim(left, right + 0.1)


g5 = sns.heatmap(heatmap_df.iloc[33:,:], cmap='bwr', ax=axs[4], cbar=False, linewidths=1, linecolor='#46494c')
axs[4].set_xlabel(None)
axs[4].set_ylabel(None)
axs[4].xaxis.set_tick_params(labeltop=False,   
                         labelbottom=False,  
                         bottom=False,     
                         labelsize=12)  

# Reformat y-axis tick labels
axs[4].yaxis.set_tick_params(labelsize=12,     
                            bottom=False, labelleft = False, labelright = True, rotation = 0)  
bottom, top = axs[4].get_ylim()
axs[4].set_ylim(bottom + 0.1, top)
left, right = axs[4].get_xlim()
axs[4].set_xlim(left, right + 0.1)

axs[0].set_ylabel('Angiogenesis &\nHypoxia Response', fontsize=14)
axs[1].set_ylabel('Inflammatory\nResponse & Cytokines', fontsize=14)
axs[2].set_ylabel('Receptor-Mediated\n Endocytosis', fontsize=14)
axs[3].set_ylabel('ECM\nOrganization', fontsize=14)
axs[4].set_ylabel('', fontsize=14)

fig.colorbar(axs[1].get_children()[0], cax=axs[5], orientation="horizontal")
axs[5].set_title('Row z-score')

plt.subplots_adjust(hspace=0)

plt.savefig('./heatmap_bulk_mono_mGAMs_FOSL2.pdf',   
           dpi = 600,                    
           bbox_inches="tight",  
           facecolor='white')  
