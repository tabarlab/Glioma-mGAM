import alluvial
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

################## Import input data ##################
input_data = pd.read_csv('./labels_ABD_gbmap_alluvial.csv') # From Seurat object metadata including the labels transfered from GBMap and Abdelfattah et al.

################## Color Maps ##################
small_color_map = c('#00B2CA','#FBDF88', '#DF5368')
cluster_colors = ["#ff595e", "#ffca3a", "#ff924c", "#c5ca30", "#36949d", "#8ac926", "#1982c4", "#565aa0"]

################## Figure 2B ##################
abd_prop_subset = input_data[['ABD_clusters', 'RNA_snn_res.0.3']]

ax = alluvial.plot(abd_prop_subset, color_side=1, alpha=.7, colors=cluster_colors)

fig = ax.get_figure()
fig.set_size_inches(10,20)

#plt.show()
plt.savefig('./alluvial_abd_to_prop.pdf',
            dpi = 600,                    
            bbox_inches="tight",          
            facecolor='white')  

prop_gbmap_subset = input_data[['RNA_snn_res.0.3', 'predicted.annotation_level_4']]

ax = alluvial.plot(prop_gbmap_subset, color_side=0, alpha=.7, colors=cluster_colors)

fig = ax.get_figure()
fig.set_size_inches(10,20)

#plt.show()
plt.savefig('./alluvial_prop_to_gbmap.pdf',
            dpi = 600,                    
            bbox_inches="tight",          
            facecolor='white')  


################## Figure 2C ##################

type_gbmap_subset = input_data[['type', 'predicted.annotation_level_4']]

ax = alluvial.plot(type_gbmap_subset, color_side=0, alpha=.7, colors=small_colors)

fig = ax.get_figure()
fig.set_size_inches(10,20)

#plt.show()
plt.savefig('./alluvial_type_to_gbmap.pdf',
            dpi = 600,                    
            bbox_inches="tight",          
            facecolor='white')  