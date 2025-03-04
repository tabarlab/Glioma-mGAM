import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
import seaborn as sns
from matplotlib.patches import Rectangle

plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.family'] = 'Arial'


def cluster_corr(corr_array, inplace=False):
    # From https://wil.yegelwel.com/cluster-correlation-matrix/
    """
    Rearranges the correlation matrix, corr_array, so that groups of highly 
    correlated variables are next to eachother 
    
    Parameters
    ----------
    corr_array : pandas.DataFrame or numpy.ndarray
        a NxN correlation matrix 
        
    Returns
    -------
    pandas.DataFrame or numpy.ndarray
        a NxN correlation matrix with the columns and rows rearranged
    """
    pairwise_distances = sch.distance.pdist(corr_array)
    linkage = sch.linkage(pairwise_distances, method='complete')
    cluster_distance_threshold = pairwise_distances.max()/2
    idx_to_cluster_array = sch.fcluster(linkage, cluster_distance_threshold, 
                                        criterion='distance')
    idx = np.argsort(idx_to_cluster_array)
    
    if not inplace:
        corr_array = corr_array.copy()
    
    if isinstance(corr_array, pd.DataFrame):
        return corr_array.iloc[idx, :].T.iloc[idx, :]
    return corr_array[idx, :][:, idx]

# Read in table with TPR, FPR and Precision values for each marker pair
df = pd.read_csv('./TPR_FPR_markers_table.csv', index_col=0)
df['pair'] = df['gene1'] + ":" + df['gene2']

# Calculate G-mean value
df['G-mean'] = np.sqrt(df['TPRs']*(1-df['FPRs']))

# Create a pairwise matrix and replace nan (same gene pairs) with 0
df_pivot = df.pivot(index='gene1', columns='gene2', values='G-mean')
df_new = df_pivot.fillna(0)

plt.figure(figsize=(16,14))
matrix = np.triu(cluster_corr(df_new))
ax = sns.heatmap(cluster_corr(df_new), annot=False, cmap='YlOrRd',cbar_kws={"shrink": 0.5}, vmin=0.18, vmax=.6, mask=matrix)
plt.tick_params(bottom=False, left=False)
plt.xticks(fontsize='13')
plt.yticks(fontsize='13')
plt.xlabel(None)
plt.ylabel(None)
# Uncomment the follwing line in order to add a rectangle highlighting any specific cell
#ax.add_patch(Rectangle((21, 22), 1, 1, fill=False, edgecolor='green', lw=3))

plt.savefig('./Gmean_heatmap_markers_mGAMs.pdf',    # Set path and filename
            dpi = 600,                     # Set dots per inch
            bbox_inches="tight",           # Remove extra whitespace around plot
            facecolor='white')   




