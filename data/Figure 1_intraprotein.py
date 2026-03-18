import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def plot_graph1(ax):
    df = pd.read_csv('consistent_residue_pairs_gp120.csv')
    x, y = df['Group1_resid'].values, df['Group2_resid'].values
    ax.scatter(x, y, s=5, color='black', alpha=0.6)
    ax.scatter(y, x, s=5, color='black', alpha=0.6)
    ax.set_xlabel('gp120 residues', size=8)
    ax.set_ylabel('gp120 residues', size=8)
    ax.set_title('gp120 intra-protein interactions', size=9, fontweight='bold')
    ax.set_aspect('equal', adjustable='box')
    ax.set_box_aspect(1)
    min_val, max_val = min(x.min(), y.min()) - 1, max(x.max(), y.max()) + 1
    ax.set_xlim(min_val, max_val)
    ax.set_ylim(min_val, max_val)
    ax.tick_params(axis='both', labelsize=8)
    ax.grid(True)

def plot_graph2(ax):
    df = pd.read_csv('consistent_residue_pairs_gp41.csv')
    x, y = df['Group1_resid'].values, df['Group2_resid'].values
    ax.scatter(x, y, s=5, color='black', alpha=0.6)
    ax.scatter(y, x, s=5, color='black', alpha=0.6)
    ax.set_xlabel('gp41 residues', size=8)
    ax.set_ylabel('gp41 residues', size=8)
    ax.set_title('gp41 intra-protein interactions', size=9, fontweight='bold')
    ax.set_aspect('equal', adjustable='box')
    ax.set_box_aspect(1)
    min_val, max_val = min(x.min(), y.min()) - 1, max(x.max(), y.max()) + 1
    ax.set_xlim(min_val, max_val)
    ax.set_ylim(min_val, max_val)
    ax.tick_params(axis='both', labelsize=8)
    ax.grid(True)

def plot_graph3(ax):
    df = pd.read_csv('consistent_residue_pairs_reverse transcriptase.csv')
    x, y = df['Group1_resid'].values, df['Group2_resid'].values
    ax.scatter(x, y, s=5, color='black', alpha=0.6)
    ax.scatter(y, x, s=5, color='black', alpha=0.6)
    ax.set_xlabel('RT residues', size=8)
    ax.set_ylabel('RT residues', size=8)
    ax.set_title('RT intra-protein interactions', size=9, fontweight='bold')
    ax.set_aspect('equal', adjustable='box')
    ax.set_box_aspect(1)
    min_val, max_val = min(x.min(), y.min()) - 1, max(x.max(), y.max()) + 1
    ax.set_xlim(min_val, max_val)
    ax.set_ylim(min_val, max_val)
    ax.tick_params(axis='both', labelsize=8)
    ax.grid(True)

def plot_graph4(ax):
    df = pd.read_csv('consistent_residue_pairs_Integrase.csv')
    x, y = df['Group1_resid'].values, df['Group2_resid'].values
    ax.scatter(x, y, s=5, color='black', alpha=0.6)
    ax.scatter(y, x, s=5, color='black', alpha=0.6)
    ax.set_xlabel('Integrase residues', size=8)
    ax.set_ylabel('Integrase residues', size=8)
    ax.set_title('Integrase intra-protein interactions', size=9, fontweight='bold')
    ax.set_aspect('equal', adjustable='box')
    ax.set_box_aspect(1)
    min_val, max_val = min(x.min(), y.min()) - 1, max(x.max(), y.max()) + 1
    ax.set_xlim(min_val, max_val)
    ax.set_ylim(min_val, max_val)
    ax.tick_params(axis='both', labelsize=8)
    ax.grid(True)

def plot_symmetric_heatmap(ax, filepath, xlabel, ylabel, title):
    df = pd.read_csv(filepath)
    pivot1 = df.pivot(index="Group1_resid", columns="Group2_resid", values="Average_Distance")
    pivot2 = df.pivot(index="Group2_resid", columns="Group1_resid", values="Average_Distance")
    full_index = sorted(set(pivot1.index).union(pivot1.columns).union(pivot2.index).union(pivot2.columns))
    pivot1 = pivot1.reindex(index=full_index, columns=full_index)
    pivot2 = pivot2.reindex(index=full_index, columns=full_index)
    combined = pd.concat([pivot1, pivot2.T]).groupby(level=0).mean()
    heatmap = ax.imshow(combined, cmap='viridis', aspect='equal', origin='lower')
    cbar = plt.colorbar(heatmap, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Distance in Å', fontsize=8)
    cbar.ax.tick_params(labelsize=8)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.tick_params(axis='both', labelsize=8)
    ax.set_title(title, fontsize=9, fontweight='bold')

fig = plt.figure(figsize=(15, 6))
gs = GridSpec(2, 4, figure=fig, height_ratios=[2, 2])

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[0, 2])
ax4 = fig.add_subplot(gs[0, 3])

ax5 = fig.add_subplot(gs[1, 0])
ax6 = fig.add_subplot(gs[1, 1])
ax7 = fig.add_subplot(gs[1, 2])
ax8 = fig.add_subplot(gs[1, 3])

plot_graph1(ax1)
plot_graph2(ax2)
plot_graph3(ax3)
plot_graph4(ax4)

plot_symmetric_heatmap(ax5,  'average_distance_gp120.csv', 'gp120 residues', 'gp120 residues', 'gp120 intra-protein interactions')
plot_symmetric_heatmap(ax6,  'average_distance_gp41.csv', 'gp41 residues', 'gp41 residues', 'gp41 intra-protein interactions')
plot_symmetric_heatmap(ax7,  'average_distance_reverse transcriptase.csv', 'RT residues', 'RT residues', 'RT intra-protein interactions')
plot_symmetric_heatmap(ax8, 'average_distance_Integrase.csv', 'Integrase residues', 'Integrase residues', 'Integrase intra-protein interactions')

labels = ['a.', 'b.', 'c.', 'd.', 'e.', 'f.', 'g.', 'h.']
axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]

for label, ax in zip(labels, axes):
    ax.text(-0.35, 1.05, label, transform=ax.transAxes, fontsize=9, fontweight='bold', va='bottom', ha='right')

plt.tight_layout()
plt.savefig('intra.png', dpi=300, bbox_inches='tight', pad_inches=0)
plt.show()
