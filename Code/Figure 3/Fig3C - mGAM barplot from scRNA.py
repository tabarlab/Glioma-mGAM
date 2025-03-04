import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.family'] = 'Arial'


df1 = pd.read_csv('./FOSL2_high_per_patient.csv')
df = df1.pivot(index='Sample', columns='FOSL2_module_positive', 
               values='pct.high')
df = df.sort_values(by=1)

df2 = pd.read_csv('./Patient_metadata.csv')

x = np.arange(len(df.index))  # the label locations
width = 0.9 # the width of the bars

# Setup plot size.
fig, (ax, ax2) = plt.subplots(2, 1, figsize=(10,5), 
                              gridspec_kw={'height_ratios': [4, 1]},  sharex='col')

# Create grid 
ax.grid(which="major", axis='y', color='#758D99', alpha=0.6, zorder=1)

# Remove spines.
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)


# Make spines slightly thicker
ax.spines['bottom'].set_linewidth(1.2)
ax.spines['left'].set_linewidth(1.2)

# Initialize the bottom at zero for the first set of bars.
bottom = np.zeros(len(df))

# Plot each bar, without edges
ax.bar(x, df[1], width, color='#21B4FD', 
       zorder=2, label = 'FOSL2 High')
ax.bar(x, df[0], width, bottom=df[1], color='#082C45', 
       zorder=2, label = 'FOSL2 Low')

  
# Reformat x-axis tick labels
ax.set_xticks(range(0, len(df.index)))
ax.set_xticklabels([])
ax.xaxis.set_tick_params(labeltop=False,      
                         labelbottom=True,  
                         bottom=False,       
                         labelsize=15,       
                         pad=1, 
                         rotation = 45)     

# Reformat y-axis tick labels
ax.set_yticks([0, 20, 40, 60, 80, 100])
ax.set_yticklabels(['0%', '20%', '40%', '60%', '80%', '100%'])
ax.yaxis.set_tick_params(pad=3,
                         labelsize=12,
                         left=True)  

# Set Legend
ax.legend(loc=(1.02,0.4), ncol=1, frameon=False, 
          handletextpad=.2, handleheight=1, prop={'size':15})
ax.set_ylabel("Percent (%)", size=15, weight="bold")

## Bottom plot

ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['bottom'].set_visible(False)

ax2.spines['left'].set_linewidth(1.2)

cdict = {'WT': '#DF5368', 'MUT': '#00B2CA', 'NB': '#FBDF88'}
cdict2 = {'IV': '#212529', "III": "#6c757d", "II": "#dee2e6", 'NB': '#f8f9fa'}
cdict3 = {'IV': '#212529', "III": "#6c757d", "II": "#dee2e6", 'NB': 'black'}

ax2.scatter(df2['X'], [0]*len(df2), 
            c = [cdict[x] for x in list(df2['Type'][::-1])], s = 250)
ax2.scatter(df2['X'], [1]*len(df2), 
            c = [cdict2[x] for x in list(df2['Grade'][::-1])], 
            edgecolors = [cdict3[x] for x in list(df2['Grade'][::-1])], 
            s = 250)
  
# Reformat x-axis tick labels
ax2.set_xticks(range(0, len(df.index)), minor=False)
ax2.set_xticklabels(df.index, ha='right',
                         rotation_mode='anchor')
ax2.xaxis.set_tick_params(labeltop=False, 
                         labelbottom=True,
                         bottom=False,    
                         labelsize=15,    
                         pad=1, 
                         rotation = 45)   

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

legend1 = plt.legend(markers, cdict.keys(), loc=(1,-.95), ncol=1, frameon=False, 
                     markerscale=2, handletextpad=.1, handleheight=1, prop={'size':15}, 
                     title = 'IDH Status',title_fontsize='15')
legend2 = plt.legend(markers2, cdict2.keys(), loc=(1.15,-1.4), ncol=1, frameon=False, 
                     markerscale=2, handletextpad=.1, handleheight=1, prop={'size':15}, 
                     title = 'Grade',title_fontsize='15')

ax2.add_artist(legend2)
ax2.add_artist(legend1)

plt.setp(ax.get_xticklabels(), visible=False)

fig.tight_layout()
