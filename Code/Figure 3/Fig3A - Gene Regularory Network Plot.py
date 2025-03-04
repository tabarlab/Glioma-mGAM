import pandas as pd
import networkx
import numpy as np
from sklearn import preprocessing

from bokeh.io import show, save
from bokeh.io import export_svgs
import svglib.svglib as svglib
from reportlab.graphics import renderPDF
from bokeh.models import Range1d, Circle, MultiLine
from bokeh.plotting import figure
from bokeh.plotting import from_networkx
from bokeh.palettes import Blues8
from bokeh.transform import linear_cmap

df = pd.read_excel('./FigR_results_all.xlsx')

edge_df = df.loc[:,['Motif', 'DORC', 'Score', 'Corr']]
edge_df.Score = abs(edge_df.Score)
edge_df = edge_df[edge_df.Score >= 1.2]
edge_df.columns = ['from', 'to', 'weight', 'corr']
edge_df = edge_df.reset_index()

x = edge_df.loc[:,'weight'].to_numpy() #returns a numpy array
x = x.reshape(-1,1)
scaler = preprocessing.MinMaxScaler()
x_scaled = scaler.fit_transform(x)
edge_df['weight_scaled'] = abs(x_scaled)

G = networkx.from_pandas_edgelist(edge_df, 'from', 'to', edge_attr=True)

degrees = dict(networkx.degree(G))
networkx.set_node_attributes(G, name='degree', values=degrees)

# bin Node sizes 
bins = np.array([5, 25, 75, 125])
inds = np.digitize(list(degrees.values()), bins)
number_to_adjust_by = 5
adjusted_node_size = dict([(node, np.log2(degree)+number_to_adjust_by) for node, degree in networkx.degree(G)])
networkx.set_node_attributes(G, name='adjusted_node_size', values=adjusted_node_size)
networkx.set_node_attributes(G, name='id', values=dict([(node, node) for node, degree in networkx.degree(G)]))

# color node by condition it's DEG in
df1 = df.loc[:,['Motif', 'cluster_motif']]
df2 = df.loc[:,['DORC', 'cluster_dorc']]
df1.columns = ['id', 'Expression']
df2.columns = ['id', 'Expression']

node_df = pd.concat([df1, df2])
node_df = node_df.drop_duplicates()
node_df.index = node_df.id
node_df.pop('id')
node_df.head().to_dict()['Expression']

condition_group = node_df.to_dict()['Expression']
networkx.set_node_attributes(G, name='condition', values=condition_group)

#Choose attributes from G network to size and color by
size_by_this_attribute = 'adjusted_node_size'
color_by_this_attribute = 'adjusted_node_size'

#Pick a color palette
color_palette = Blues8

#Choose a title
title = 'GAM GRN'

#Establish which categories will appear when hovering over each node
HOVER_TOOLTIPS = [
       ("Character", "@id"),
        ("Degree", "@degree")
]

#Create a plot — set dimensions, toolbar, and title
plot = figure(tooltips = HOVER_TOOLTIPS,
              tools="pan,wheel_zoom,save,reset", active_scroll='wheel_zoom',
            x_range=Range1d(-10.1, 10.1), y_range=Range1d(-10.1, 10.1), title=title)

mapping = dict((n, i) for i, n in enumerate(G.nodes))
H = networkx.relabel_nodes(G, mapping)

#Create a network graph object
network_graph = from_networkx(H, networkx.kamada_kawai_layout, scale=10, center=(0, 0))

#Set node sizes and colors according to node degree (color as spectrum of color palette)
minimum_value_color = min(network_graph.node_renderer.data_source.data[color_by_this_attribute])
maximum_value_color = max(network_graph.node_renderer.data_source.data[color_by_this_attribute])
network_graph.node_renderer.glyph = Circle(size=size_by_this_attribute, fill_color=linear_cmap(color_by_this_attribute, color_palette, minimum_value_color, maximum_value_color))

#Set edge opacity and width
network_graph.edge_renderer.glyph = MultiLine(line_alpha=0.5, line_width=1)

plot.renderers.append(network_graph)

show(plot)

# Create empty dictionaries
cdict = {'WT': '#DF5368', 'MUT': '#00B2CA', 'CTR': '#FBDF88', 'Not DEG': '#379A8B'}
modularity_color = {}
#Loop through each community in the network
for gene, condition in networkx.get_node_attributes(G, 'condition').items():
    #For each member of the community, add their community number and a distinct color
    modularity_color[gene] = cdict[condition]

networkx.set_node_attributes(G, modularity_color, 'modularity_color')

relationship = networkx.get_edge_attributes(G,'corr')

edge_color = {}
for edge in relationship:
    if relationship[edge] > 0:
        edge_color[edge] = 'green'
    else: 
        edge_color[edge] = 'red'

networkx.set_edge_attributes(G, edge_color, 'relation_color')

#Choose attributes from G network to size and color by — setting manual size (e.g. 10) or color (e.g. 'skyblue') also allowed
size_by_this_attribute = 'adjusted_node_size'
color_by_this_attribute = 'modularity_color'
line_color_by_this_attribute = 'relation_color'

#Choose a title
title = 'Glioma-associated Macrophages Network'

#Establish which categories will appear when hovering over each node
HOVER_TOOLTIPS = [
       ("Character", "@id"),
        ("Degree", "@degree"),
         ("Condition", "@condition"),
        ("Color", "$color[swatch]:modularity_color"),
]

#Create a plot — set dimensions, toolbar, and title
plot = figure(tooltips = HOVER_TOOLTIPS,
              tools="pan,wheel_zoom,save,reset, tap", active_scroll='wheel_zoom',
            x_range=Range1d(-10.1, 10.1), y_range=Range1d(-10.1, 10.1), width=900, height=900, title=title)

mapping = dict((n, i) for i, n in enumerate(G.nodes))
H = networkx.relabel_nodes(G, mapping)

network_graph = from_networkx(H, networkx.kamada_kawai_layout, scale=10, center=(0, -2))

#Set node sizes and colors
network_graph.node_renderer.glyph = Circle(size=size_by_this_attribute, fill_color=color_by_this_attribute, line_color = color_by_this_attribute)

#Set edge opacity and width
network_graph.edge_renderer.glyph = MultiLine(line_alpha=.8, line_width='weight_scaled', line_color = line_color_by_this_attribute)

x, y = zip(*network_graph.layout_provider.graph_layout.values())

#plot.add_layout(labels)
plot.xgrid.grid_line_color = None  # make the x-axis grid lines invisible
plot.ygrid.grid_line_color = None  # make the x-axis grid lines invisible
plot.renderers.append(network_graph)

# step 1: bokeh save as svg
plot.output_backend = "svg"
export_svgs(plot, filename = './GRN_GAMs' + '.svg')

# step 2: read in svg
svg = svglib.svg2rlg('./GRN_GAMs'+".svg")

# step 3: save as pdf
renderPDF.drawToFile(svg, './GRN_GAMs'+".pdf")
save(plot, './GRN_GAMs'+".html")

