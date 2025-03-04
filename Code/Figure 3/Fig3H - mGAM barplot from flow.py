import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.family'] = 'Arial'


path = './mGAM_flow_report_percentages.xlsx'
df = pd.read_excel(path)

y = np.arange(len(df.index))  # the label locations
width = 0.80  # the width of the bars
col_list = ['#ed7672', '#7cb0f4']

df = df.sort_values(by='Proportion', ascending=True)

# Setup plot size.
fig, (ax, ax2) = plt.subplots(2, 1, figsize=(15,8), 
                              gridspec_kw={'height_ratios': [5, .5]},  sharex='col')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(True)
ax.spines['bottom'].set_visible(False)

# Make left spine slightly thicker
ax.spines['left'].set_linewidth(1.5)

bottom = np.zeros(len(df))

ax.bar(y - width/0.9, df['Proportion'], width, color=col_list[0], 
       zorder=2, ecolor='black', capsize=4)
ax.bar(y - width/0.9, -df['DN%']/100, width, color=col_list[1], 
       zorder=2, ecolor='black', capsize=4)
ax.hlines([0], -5, y[-1]+10, colors = 'black', lw = 1.5, zorder = 4)

ax.set_xlim(-2, 33.5)  

# Reformat x-axis tick labels
ax.set_xticks(y- width/0.9)
ax.set_xticklabels([], ha = 'right')           
ax.xaxis.set_tick_params(labeltop=False,   
                         labelbottom=False,
                         bottom=False, 
                         labelsize=12, 
                         pad=1,
                         rotation=45)  


ax.set_ylim(-1.05, 1.05)

# Reformat y-axis tick labels
ax.set_yticks([-1, -.8, -.6, -.4, -.2, 0, .20, .40, .60, .80, 1])
ax.set_yticklabels([1, .8, .6, .4, .2, 0, .20, .40, .60, .80, 1],  
                   ha = 'right')         
ax.yaxis.set_tick_params(pad=3,          
                         labelsize=12,   
                         bottom=True) 

ax.set_ylabel("mGAM-" + " "*33 + "mGAM+", size=15, weight='bold')

## Bottom plot
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['bottom'].set_visible(False)

# Make left spine slightly thicker
ax2.spines['left'].set_linewidth(1.5)

group = df['Type']
cdict = {'Wildtype': '#DF5368', 'Mutant': '#00B2CA', 'PBMC': '#f8f9fa'}
cdict4 = {'Wildtype': '#DF5368', 'Mutant': '#00B2CA', 'PBMC': 'black'}

group2 = df['Grade']
cdict2 = {'4': '#212529', "3": "#6c757d", "2": "#dee2e6", "0": '#f8f9fa'}
cdict3 = {'4': '#212529', "3": "#6c757d", "2": "#dee2e6", "0": 'black'}

# Plot each bar, without edges
ax2.scatter(y-width/0.9, [0]*len(df), s = 300, c = [cdict[x] for x in list(df['Type'])], 
            edgecolors = [cdict4[x] for x in list(df['Type'])], marker='o', zorder=1)
ax2.scatter(y-width/0.9, [1]*len(df), s = 300, c = [cdict2[str(x)] for x in list(df['Grade'])], 
            edgecolors = [cdict3[str(x)] for x in list(df['Grade'])], marker='o', zorder=2)

# Reformat x-axis tick labels
ax2.set_xticks(y- width/0.9, minor=False)
ax2.set_xticklabels(df['MSK-ID'], ha='right',
                         rotation_mode='anchor')
ax2.xaxis.set_tick_params(labeltop=False, 
                         labelbottom=True,
                         bottom=False,  
                         labelsize=15,  
                         pad=1, 
                         rotation = 45)            

# Shrink y-lim to make plot a bit tighter
ax2.set_ylim(-.5, 1.5)

# Reformat y-axis tick labels
ax2.set_yticks([0, 1])
ax2.set_yticklabels(['IDH Status', 'Grade'])
ax2.yaxis.set_tick_params(pad=3,            
                         labelsize=12,      
                         left=True)    

# Set Legend
markers = [plt.Line2D([0,0],[0,0], color=color, marker='o', 
                      linestyle='') for color in cdict.values()] 
markers2 = [plt.Line2D([0,0],[0,0], color=color, marker='o', 
                       linestyle='') for color in cdict2.values()]
ax.legend(markers2, cdict2.keys(), loc=(1.05,.23), ncol=2, frameon=False, 
          markerscale=2, handletextpad=.1, handleheight=1, prop={'size':15}, 
          title = 'Grade', title_fontsize='15')
ax2.legend(markers, cdict.keys(), loc=(1.05,.95), ncol=1, frameon=False, 
           markerscale=2, handletextpad=.1, handleheight=1, prop={'size':15}, 
           title = 'IDH Status', title_fontsize='15')

fig.tight_layout()
