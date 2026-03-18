import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from matplotlib.gridspec import GridSpecFromSubplotSpec

def plot_graph3(ax):
    df = pd.read_csv('consistent_residue_pairs_interprotein interactions_gp120&gp41.csv')
    df['Group1_resid'] = df['Group1_resid'].astype('Int64').astype(str)
    df['Group2_resid'] = df['Group2_resid'].astype('Int64').astype(str)
    df['Group1_node'] = 'G1_' + df['Group1_resid']
    df['Group2_node'] = 'G2_' + df['Group2_resid']

    G = nx.DiGraph()
    G.add_edges_from(zip(df['Group1_node'], df['Group2_node']))

    top_nodes = sorted(df['Group1_node'].unique(), key=lambda x: int(x.split('_')[1]))
    bottom_nodes = sorted(df['Group2_node'].unique(), key=lambda x: int(x.split('_')[1]))

    def horizontal_pos(nodes, y):
        n = len(nodes)
        return {node: (i / (n - 1) if n > 1 else 0.5, y) for i, node in enumerate(nodes)}

    pos = {}
    pos.update(horizontal_pos(top_nodes, y=1))
    pos.update(horizontal_pos(bottom_nodes, y=0))

    colors = ['tomato' if node.startswith('G1_') else 'yellow' for node in G.nodes()]

    nx.draw(G, pos, ax=ax, node_color=colors, node_size=90,
            font_size=8, edge_color='gray', width=1.5, arrows=False)

    for node, (x, y) in pos.items():
        label = node.split('_')[1]
        ax.text(x, y, label, fontsize=9, ha='center', va='center')

    red_patch = mpatches.Patch(color='tomato', label='gp120')
    yellow_patch = mpatches.Patch(color='yellow', label='gp41')
    ax.legend(handles=[red_patch, yellow_patch], title="Key", loc='center left', fontsize=6, title_fontsize=7)

    ax.set_xlim(-0.009, 1.009)
    ax.set_ylim(-0.1, 1.1)
    ax.set_title('gp120-gp41 inter-protein interactions', size=9, fontweight='bold')
    ax.axis('off')

def plot_graph4(ax):
    df = pd.read_csv('consistent_residue_pairsinterprotein interactions_RT&I.csv')
    df['Group1_resid'] = df['Group1_resid'].astype('Int64').astype(str)
    df['Group2_resid'] = df['Group2_resid'].astype('Int64').astype(str)
    df['Group1_node'] = 'G1_' + df['Group1_resid']
    df['Group2_node'] = 'G2_' + df['Group2_resid']

    G = nx.DiGraph()
    G.add_edges_from(zip(df['Group1_node'], df['Group2_node']))

    top_nodes = sorted(df['Group1_node'].unique(), key=lambda x: int(x.split('_')[1]))
    bottom_nodes = sorted(df['Group2_node'].unique(), key=lambda x: int(x.split('_')[1]))

    def horizontal_pos(nodes, y):
        n = len(nodes)
        return {node: (i / (n - 1) if n > 1 else 0.5, y) for i, node in enumerate(nodes)}

    pos = {}
    pos.update(horizontal_pos(top_nodes, y=1))
    pos.update(horizontal_pos(bottom_nodes, y=0))

    colors = ['tomato' if node.startswith('G1_') else 'yellow' for node in G.nodes()]

    # Draw graph WITHOUT labels
    nx.draw(G, pos, ax=ax, node_color=colors, node_size=90,
            font_size=8, edge_color='gray', width=1.5, arrows=False)

    # Manually draw rotated labels
    for node, (x, y) in pos.items():
        label = node.split('_')[1]
        ax.text(x, y, label, fontsize=9, ha='center', va='center')

    red_patch = mpatches.Patch(color='tomato', label='RT')
    yellow_patch = mpatches.Patch(color='yellow', label='Integrase')
    ax.legend(handles=[red_patch, yellow_patch], title="Key", loc='center left', fontsize=6, title_fontsize=7)

    ax.set_xlim(-0.009, 1.009)
    ax.set_ylim(-0.1, 1.1)
    ax.set_title('RT and Integrase complex inter-protein interactions', size=9, fontweight='bold')
    ax.axis('off')

def plot_heatmap(ax, filepath, xlabel, ylabel, title):
    df = pd.read_csv(filepath)
    pivot_table = df.pivot(index="Group1_resid", columns="Group2_resid", values="Average_Distance")
    heatmap = ax.imshow(pivot_table, cmap='viridis', aspect='auto', origin='lower')
    cbar = plt.colorbar(heatmap, ax=ax)
    cbar.set_label('Distance in Å', fontsize=8)
    cbar.ax.tick_params(labelsize=8)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.tick_params(axis='both', labelsize=8)
    ax.set_title(title, fontsize=9, fontweight='bold')

row1 = 1.96
row2 = 3.0  
row3 = 1.96

fig = plt.figure(figsize=(15, 8))

gs = GridSpec(3, 3, figure=fig, height_ratios=[3.0, 3.0, 1.96])

ax3 = fig.add_subplot(gs[0, :])
ax4 = fig.add_subplot(gs[1, :])

inner_gs = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[2, :], wspace=0.3, width_ratios=[1,1])

ax5 = fig.add_subplot(inner_gs[0, 0]) 
ax6 = fig.add_subplot(inner_gs[0, 1])  

plot_graph3(ax3)
plot_graph4(ax4)
plot_heatmap(ax5, 'average_distance_interprotein_gp120&gp41.csv', 'gp41', 'gp120', 'Inter-protein interactions')
plot_heatmap(ax6, 'average_distance_interprotein_RT&I.csv', 'Integrase', 'RT', 'Inter-protein interactions')

ax3.text(-0.02, 1.02, 'a.', transform=ax3.transAxes, fontsize=10, fontweight='bold', va='bottom', ha='left')
ax4.text(-0.02, 1.02, 'b.', transform=ax4.transAxes, fontsize=10, fontweight='bold', va='bottom', ha='left')
ax5.text(-0.02, 1.02, 'c.', transform=ax5.transAxes, fontsize=10, fontweight='bold', va='bottom', ha='left')
ax6.text(-0.02, 1.02, 'd.', transform=ax6.transAxes, fontsize=10, fontweight='bold', va='bottom', ha='left')

plt.tight_layout()
plt.savefig('interprotein.png', dpi=300, bbox_inches='tight', pad_inches=0)
plt.show()