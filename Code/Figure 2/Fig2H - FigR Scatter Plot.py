#------------------------------------------------------------
# FigR Scatter Plot
#------------------------------------------------------------
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rc
from math import ceil, floor

rc('text', usetex=False)
plt.rcParams['figure.dpi'] = 100
col_list = ['#006BA2', '#DB444B', '#EBB434', '#379A8B', '#9A607F', '#3EBCD2']

FigR_results = pd.read_table('./FigR_res_all_tf_DORCS_separate_all.txt')

# List of differetially expressed genes for each of the archetypal cell groups
DE_list = pd.read_csv('./PURE_DE.csv')
DE_list = DE_list[(DE_list['p_val_adj'] <= 0.05)].iloc[:, 6:]
DE_list.columns = ['cluster_motif', 'Motif']

FigR_results = pd.merge(FigR_results, DE_list, on='Motif', how="outer")

DE_list.columns = ['cluster_dorc', 'DORC']
FigR_results = pd.merge(FigR_results, DE_list, on='DORC', how="outer")

match_cluster = []

for i, row in FigR_results.iterrows():
    if row['Corr.log10P'] < 0:
        if row['cluster_motif'] != row['cluster_dorc']:
            match_cluster.append(1)
        else:
            match_cluster.append(0)
    else:
        if row['cluster_motif'] == row['cluster_dorc']:
            match_cluster.append(1) 
        else:
            match_cluster.append(0)

# match_cluster column contains information on whether the TF and downstream target are differentially expressed in the same archteypal cell type.
FigR_results['match_cluster'] = match_cluster


filtered_table = FigR_results[FigR_results['match_cluster'] == 1]

fig, ax = plt.subplots(figsize=(9,9))

# Create grid 
ax.grid(which="major", axis='both', color='#758D99', alpha=0.4, zorder=1)

# Remove splines.
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)

xmin=floor(float(filtered_table['Corr.log10P'].min())) 
xmax=ceil(float(filtered_table['Corr.log10P'].max()))
ymin=floor(float(filtered_table['Enrichment.log10P'].min()))
ymax=ceil(float(filtered_table['Enrichment.log10P'].max()))


ax.scatter(filtered_table[[a == 'WT' for a in filtered_table.cluster_motif]]['Corr.log10P'], 
           filtered_table[[a == 'WT' for a in filtered_table.cluster_motif]]['Enrichment.log10P'], 
           label = 'WT',s=30, marker='o', zorder=4, c='#DB444B', alpha=0.4)
ax.scatter(filtered_table[[a == 'MUT' for a in filtered_table.cluster_motif]]['Corr.log10P'], 
           filtered_table[[a == 'MUT' for a in filtered_table.cluster_motif]]['Enrichment.log10P'], 
           label = 'MUT',s=30, marker='^', zorder=4, c='#006BA2', alpha=0.4)
ax.scatter(filtered_table[[a == 'NB' for a in filtered_table.cluster_motif]]['Corr.log10P'], 
           filtered_table[[a == 'NB' for a in filtered_table.cluster_motif]]['Enrichment.log10P'], 
           label = 'CTR',s=30, marker='s', zorder=3, c='#EBB434', alpha=0.4)

ax.vlines(np.array([0]), ymin-.5, ymax+.5, linestyles=(0,(5,5)), colors=['k', 'k'], lw=1.5, zorder=2)
ax.hlines(np.array([0]), -(np.array([abs(xmin), abs(xmax)]).max())-.5, 
          (np.array([abs(xmin), abs(xmax)]).max())+.5, linestyles=(0,(5,5)), colors=['k', 'k'], lw=1.5, zorder=2)

# Set xlim
ax.set_xlim(xmin-0.5, 5.2+0.5)

# Set ylim
ax.set_ylim(ymin-.5, ymax+.5)

# Reformat x-axis tick labels
ax.xaxis.set_tick_params(labeltop=False,   
                         labelbottom=True,
                         bottom=False,     
                         labelsize=15,      
                         pad=5)                

# Reformat y-axis tick labels
ax.yaxis.set_tick_params(pad=5,          
                         labelsize=15,    
                         bottom=False)     

ax.legend(loc=(.31,1.02), ncol=5, markerscale = 1.5, fontsize=13, 
          frameon=False, handletextpad=.02, handleheight=1)

ax.set_xlabel('Corr.log10p', fontsize=15)
ax.set_ylabel('Enrichment.log10P', fontsize=15)

ax.text(x=0.06, y=.935, s=(''), transform=fig.transFigure, ha='left', fontsize=20, weight='bold', alpha=.8)

plt.tight_layout()

plt.savefig('./FigR_results.pdf',   
            dpi = 300,               
            bbox_inches="tight",         
            facecolor='white')           


