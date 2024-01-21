#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 13:10:05 2024

@author: shame
"""
#%%
import os
os.chdir("/home/shame/Goldenboy_Project/scripts/")
#%%
#!mkdir data
#!wget http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz -O data/pbmc3k_filtered_gene_bc_matrices.tar.gz
#!cd data; tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz
#!mkdir write
#%%
import numpy as np
import pandas as pd
import scanpy as sc
#%%
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')
#%%
results_file = 'write/pbmc3k.h5ad'  # the file that will store the analysis results
#%%
adata = sc.read_10x_mtx(
    'data/filtered_gene_bc_matrices/hg19/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)    
                          # write a cache file for faster subsequent reading
#%%
adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`
#%%
adata
#%%
sc.pl.highest_expr_genes(adata, n_top=20)
#%%
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
#%%
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
#%%
adata
#Look at customizing plots section for correction of the axes
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
#%%
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]
#%%
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)
adata.raw = adata
#%%
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
#%%
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata, color='CST3')
sc.pl.pca_variance_ratio(adata, log=True)
#%%
adata.write(results_file)
adata
#%%
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=8)
#%%
sc.tl.umap(adata)
sc.pl.umap(adata, color=['CST3', 'NKG7', 'PPBP'])
sc.pl.umap(adata, color=['CST3', 'NKG7', 'PPBP'], use_raw=False)
#%%
sc.tl.leiden(adata)
sc.pl.umap(adata, color=['leiden', 'CST3', 'NKG7'])
adata.write(results_file)
#%%
sc.tl.rank_genes_groups(adata, 'leiden', method='logreg')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
#%%
marker_genes = ['IL7R', 'CD79A', 'MS4A1', 'CD8A', 'CD8B', 'LYZ', 'CD14',
                'LGALS3', 'S100A8', 'GNLY', 'NKG7', 'KLRB1',
                'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'PPBP']
#%%
pd.DataFrame(adata.uns['rank_genes_groups']['names'])
#%%
result = adata.uns['rank_genes_groups']
result
groups = result['names'].dtype.names
adata
groups
#%%
pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names']}).head(5)
#%%
#comparison between specific groups
sc.tl.rank_genes_groups(adata, 'leiden', groups=['9'], reference='1', method='wilcoxon')
sc.pl.rank_genes_groups(adata, groups=['9'], n_genes=20)
sc.pl.rank_genes_groups_violin(adata, groups=['9'], n_genes=8)
#%%
#Comparison between rest of groups
sc.tl.rank_genes_groups(adata, 'leiden', method='logreg')
sc.pl.rank_genes_groups_violin(adata, groups='0', n_genes=8)
#%%
#Comparison between genes and the rest of the groups
sc.pl.violin(adata, ['CST3', 'NKG7', 'PPBP'], groupby='leiden')
#%%
new_cluster_names = [
    'CD4 T', 'CD14 Monocytes',
    'B', 'CD8 T',
    'NK', 'FCGR3A Monocytes',
    'Dendritic', 'Megakaryocytes', "8", "9", "10", "11", "12", "13", "14"]
adata.rename_categories('leiden', new_cluster_names)
sc.pl.umap(adata, color='leiden', legend_loc='on data', title='', frameon=False)
#%%
sc.pl.dotplot(adata, marker_genes, groupby='leiden');
sc.pl.stacked_violin(adata, marker_genes, groupby='leiden', var_group_rotation=90);
#%%
adata.write(results_file, compression='gzip')  # `compression='gzip'` saves disk space, but slows down writing and subsequent reading
adata.raw.to_adata().write('./write/pbmc3k_withoutX.h5ad')
#%%
# Export single fields of the annotation of observations
# adata.obs[['n_counts', 'louvain_groups']].to_csv(
#     './write/pbmc3k_corrected_louvain_groups.csv')

# Export single columns of the multidimensional annotation
# adata.obsm.to_df()[['X_pca1', 'X_pca2']].to_csv(
#     './write/pbmc3k_corrected_X_pca.csv')

# Or export everything except the data using `.write_csvs`.
# Set `skip_data=False` if you also want to export the data.
# adata.write_csvs(results_file[:-5], )

#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
import scanpy as sc
#%%
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
results_file = './write/paul15.h5ad'
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(3, 3), facecolor='white')  # low dpi (dots per inch) yields small inline figures
#%%
adata = sc.datasets.paul15()
#%%
adata
adata.X = adata.X.astype('float64')  # this is not required and results will be comparable without it
#%%
sc.pp.recipe_zheng17(adata)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)
sc.tl.draw_graph(adata)
sc.pl.draw_graph(adata, color='paul15_clusters', legend_loc='on data')
#%%
# How to denoise the graph
sc.tl.diffmap(adata)
sc.pp.neighbors(adata, n_neighbors=10, use_rep='X_diffmap')
sc.tl.draw_graph(adata)
sc.pl.draw_graph(adata, color='paul15_clusters', legend_loc='on data')
#%%
sc.tl.leiden(adata, resolution=1.0)
sc.tl.paga(adata, groups='leiden')
sc.pl.paga(adata, color=['leiden', 'Hba-a2', 'Elane', 'Irf8'])
sc.pl.paga(adata, color=['leiden', 'Itga2b', 'Prss34', 'Cma1'])
#%%
adata.obs['leiden'].cat.categories
adata.obs['leiden_anno'] = adata.obs['leiden']
adata.obs['leiden_anno']
#%%
adata.obs['leiden_anno'].cat.categories = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10/Ery', '11', '12',
       '13', '14', '15', '16/Stem', '17', '18', '19/Neu', '20/Mk', '21', '22/Baso', '23', '24/Mo', '25', '26', '27', '28', '29', '30', '31', '32']
#%%
sc.tl.paga(adata, groups='leiden_anno')
sc.pl.paga(adata, threshold=0.03, show=False)
#%%
sc.tl.draw_graph(adata, init_pos='paga')
sc.pl.draw_graph(adata, color=['leiden_anno', 'Itga2b', 'Prss34', 'Cma1'], legend_loc='on data')
#%%
pl.figure(figsize=(8, 2))
for i in range(28):
    pl.scatter(i, 1, c=sc.pl.palettes.zeileis_28[i], s=200)
pl.show()
#%%
zeileis_colors = np.array(sc.pl.palettes.zeileis_28)
new_colors = np.array(adata.uns['leiden_anno_colors'])
new_colors[[16]] = zeileis_colors[[12]]  # Stem colors / green
new_colors[[10, 17, 5, 3, 15, 6, 18, 13, 7, 12]] = zeileis_colors[[5, 5, 5, 5, 11, 11, 10, 9, 21, 21]]  # Ery colors / red
new_colors[[20, 8]] = zeileis_colors[[17, 16]]  # Mk early Ery colors / yellow
new_colors[[4, 0]] = zeileis_colors[[2, 8]]  # lymph progenitors / grey
new_colors[[22]] = zeileis_colors[[18]]  # Baso / turquoise
new_colors[[19, 14, 2]] = zeileis_colors[[6, 6, 6]]  # Neu / light blue
new_colors[[24, 9, 1, 11]] = zeileis_colors[[0, 0, 0, 0]]  # Mo / dark blue
new_colors[[21, 23]] = zeileis_colors[[25, 25]]  # outliers / grey
sc.pl.paga_compare(
    adata, threshold=0.03, title='', right_margin=0.2, size=10, edge_width_scale=0.5,
    legend_fontsize=12, fontsize=12, frameon=False, edges=True, save=True)
#%%
adata.uns['iroot'] = np.flatnonzero(adata.obs['leiden_anno']  == '16')[0]
sc.tl.dpt(adata)
gene_names = ['Gata2', 'Gata1', 'Klf1', 'Epor', 'Hba-a2',  # erythroid
              'Elane', 'Cebpe', 'Gfi1',                    # neutrophil
              'Irf8', 'Csf1r', 'Ctsg']                     # monocyte
adata_raw = sc.datasets.paul15()
sc.pp.log1p(adata_raw)
sc.pp.scale(adata_raw)
adata.raw = adata_raw
sc.pl.draw_graph(adata, color=['leiden_anno', 'dpt_pseudotime'], legend_loc='on data')
#%%
paths = [('erythrocytes', [16, 12, 7, 13, 18, 6, 5, 10]),
         ('neutrophils', [16, 0, 4, 2, 14, 19]),
         ('monocytes', [16, 0, 4, 11, 1, 9, 24])]
adata.obs['distance'] = adata.obs['dpt_pseudotime']
adata.obs['clusters'] = adata.obs['leiden_anno']  # just a cosmetic change
adata.uns['clusters_colors'] = adata.uns['leiden_anno_colors']
#%%
_, axs = pl.subplots(ncols=3, figsize=(6, 2.5), gridspec_kw={'wspace': 0.05, 'left': 0.12})
pl.subplots_adjust(left=0.05, right=0.98, top=0.82, bottom=0.2)
for ipath, (descr, path) in enumerate(paths):
    _, data = sc.pl.paga_path(
        adata, path, gene_names,
        show_node_names=False,
        ax=axs[ipath],
        ytick_fontsize=12,
        left_margin=0.15,
        n_avg=50,
        annotations=['distance'],
        show_yticks=True if ipath==0 else False,
        show_colorbar=False,
        color_map='Greys',
        groups_key='clusters',
        color_maps_annotations={'distance': 'viridis'},
        title='{} path'.format(descr),
        return_data=True,
        show=False)
    data.to_csv('./write/paga_path_{}.csv'.format(descr))
pl.savefig('./figures/paga_path_paul15.pdf')
pl.show()
#%%
import scanpy as sc
import pandas as pd
from matplotlib.pyplot import rc_context
#%%
sc.set_figure_params(dpi=100, color_map = 'viridis_r')
sc.settings.verbosity = 1
sc.logging.print_header()
#%%
pbmc = sc.datasets.pbmc68k_reduced()
# inspect pbmc contents
pbmc
#%%
# rc_context is used for the figure size, in this case 4x4
with rc_context({'figure.figsize': (4, 4)}):
    sc.pl.umap(pbmc, color='CD79A')
#%%
with rc_context({'figure.figsize': (3, 3)}):
    sc.pl.umap(pbmc, color=['CD79A', 'MS4A1', 'IGJ', 'CD3D', 'FCER1A', 'FCGR3A', 'n_counts', 'bulk_labels'], s=50, frameon=False, ncols=4, vmax='p99')
#%%
# compute clusters using the leiden method and store the results with the name `clusters`
sc.tl.leiden(pbmc, key_added='clusters', resolution=0.5)
with rc_context({'figure.figsize': (5, 5)}):
    sc.pl.umap(pbmc, color='clusters', add_outline=True, legend_loc='on data',
               legend_fontsize=12, legend_fontoutline=2,frameon=False,
               title='clustering of cells', palette='Set1')
#%%
marker_genes_dict = {
    'B-cell': ['CD79A', 'MS4A1'],
    'Dendritic': ['FCER1A', 'CST3'],
    'Monocytes': ['FCGR3A'],
    'NK': ['GNLY', 'NKG7'],
    'Other': ['IGLL1'],
    'Plasma': ['IGJ'],
    'T-cell': ['CD3D'],
}
#%%
sc.pl.dotplot(pbmc, marker_genes_dict, 'clusters', dendrogram=True)
#%%
# create a dictionary to map cluster to annotation label
cluster2annotation = {
     '0': 'Monocytes',
     '1': 'Dendritic',
     '2': 'T-cell',
     '3': 'NK',
     '4': 'B-cell',
     '5': 'Dendritic',
     '6': 'Plasma',
     '7': 'Other',
     '8': 'Dendritic',
}

# add a new `.obs` column called `cell type` by mapping clusters to annotation using pandas `map` function
pbmc.obs['cell type'] = pbmc.obs['clusters'].map(cluster2annotation).astype('category')
sc.pl.dotplot(pbmc, marker_genes_dict, 'cell type', dendrogram=True)
#%%
sc.pl.umap(pbmc, color='cell type', legend_loc='on data',
           frameon=False, legend_fontsize=10, legend_fontoutline=2)
#%%
with rc_context({'figure.figsize': (4.5, 3)}):
    sc.pl.violin(pbmc, ['CD79A', 'MS4A1'], groupby='clusters' )
#%%
with rc_context({'figure.figsize': (4.5, 3)}):
    sc.pl.violin(pbmc, ['n_genes', 'percent_mito'], groupby='clusters', stripplot=False, inner='box')  # use stripplot=False to remove the internal dots, inner='box' adds a boxplot inside violins  
#%%
ax = sc.pl.stacked_violin(pbmc, marker_genes_dict, groupby='clusters', swap_axes=False, dendrogram=True)
#%%
sc.pl.matrixplot(pbmc, marker_genes_dict, 'clusters', dendrogram=True, cmap='Blues', standard_scale='var', colorbar_title='column scaled\nexpression')
#%%
# scale and store results in layer
pbmc.layers['scaled'] = sc.pp.scale(pbmc, copy=True).X
sc.pl.matrixplot(pbmc, marker_genes_dict, 'clusters', dendrogram=True,
                 colorbar_title='mean z-score', layer='scaled', vmin=-2, vmax=2, cmap='RdBu_r')
#%%
import matplotlib.pyplot as plt

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20,4), gridspec_kw={'wspace':0.9})

ax1_dict = sc.pl.dotplot(pbmc, marker_genes_dict, groupby='bulk_labels', ax=ax1, show=False)
ax2_dict = sc.pl.stacked_violin(pbmc, marker_genes_dict, groupby='bulk_labels', ax=ax2, show=False)
ax3_dict = sc.pl.matrixplot(pbmc, marker_genes_dict, groupby='bulk_labels', ax=ax3, show=False, cmap='viridis')
#%%
ax = sc.pl.heatmap(pbmc, marker_genes_dict, groupby='clusters', cmap='viridis', dendrogram=True)
#%%
ax = sc.pl.heatmap(pbmc, marker_genes_dict, groupby='clusters', layer='scaled', vmin=-2, vmax=2, cmap='RdBu_r', dendrogram=True, swap_axes=True, figsize=(11,4))
#%%
ax = sc.pl.tracksplot(pbmc, marker_genes_dict, groupby='clusters', dendrogram=True)
#%%
#This Line is going to be the most important for the following plots
sc.tl.rank_genes_groups(pbmc, groupby='clusters', method='wilcoxon')
sc.pl.rank_genes_groups_dotplot(pbmc, n_genes=4)
#%%
sc.pl.rank_genes_groups_dotplot(pbmc, n_genes=4, values_to_plot='logfoldchanges', min_logfoldchange=3, vmax=7, vmin=-7, cmap='bwr')
#%% 
# You can have the dotplto focus on particular groups you can also do this for violin, heatmap, and matrix plots
sc.pl.rank_genes_groups_dotplot(pbmc, n_genes=30, values_to_plot='logfoldchanges', min_logfoldchange=4, vmax=7, vmin=-7, cmap='bwr', groups=['1', '5'])
#%%
sc.pl.rank_genes_groups_matrixplot(pbmc, n_genes=3, use_raw=False, vmin=-3, vmax=3, cmap='bwr', layer='scaled')
#%%
sc.pl.rank_genes_groups_stacked_violin(pbmc, n_genes=3, cmap='viridis_r')
#%%
sc.pl.rank_genes_groups_heatmap(pbmc, n_genes=3, use_raw=False, swap_axes=True, vmin=-3, vmax=3, cmap='bwr', layer='scaled', figsize=(10,7), show=False);
#%%
sc.pl.rank_genes_groups_heatmap(pbmc, n_genes=10, use_raw=False, swap_axes=True, show_gene_labels=False,
                                vmin=-3, vmax=3, cmap='bwr')
#%%
sc.pl.rank_genes_groups_tracksplot(pbmc, n_genes=3)
#%%
with rc_context({'figure.figsize': (9, 1.5)}):
    sc.pl.rank_genes_groups_violin(pbmc, n_genes=20, jitter=False)
#%%
# compute hierarchical clustering using PCs (several distance metrics and linkage methods are available).
sc.tl.dendrogram(pbmc, 'bulk_labels')   
ax = sc.pl.dendrogram(pbmc, 'bulk_labels') 
#%%
ax = sc.pl.correlation_matrix(pbmc, 'bulk_labels', figsize=(5,3.5))
#%%
import scanpy as sc
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as mcolors

# Inital setting for plot size
from matplotlib import rcParams
FIGSIZE=(3,3)
rcParams['figure.figsize']=FIGSIZE

import warnings
warnings.filterwarnings('ignore')
#%%
adata = sc.datasets.pbmc68k_reduced()
#%%
# Examples of returned objects from the UMAP function

print('Categorical plots:')
axes=sc.pl.umap(adata,color=["bulk_labels"],show=False)
print('Axis from a single category plot:',axes)
plt.close()
#%%
axes=sc.pl.umap(adata,color=["bulk_labels",'S_score'],show=False)
print('Axes list from two categorical plots:',axes)
plt.close()
#%%
fig=sc.pl.umap(adata,color=["bulk_labels"],return_fig=True)
print('Axes list from a figure with one categorical plot:',fig.axes)
plt.close()
#%%
print('\nContinous plots:')
axes=sc.pl.umap(adata,color=["IGJ"],show=False)
print('Axes from one continuous plot:',axes)
plt.close()
#%%
fig=sc.pl.umap(adata,color=["IGJ"],return_fig=True)
print('Axes list from a figure of one continous plot:',fig.axes)
plt.close()
#%%
axes=sc.pl.dotplot(adata, ['CD79A', 'MS4A1'], 'bulk_labels', show=False)
print('Axes returned from dotplot object:',axes)
dp=sc.pl.dotplot(adata, ['CD79A', 'MS4A1'], 'bulk_labels', return_fig=True)
print('DotPlot object:',dp)
plt.close()
#%%
# Define matplotlib Axes
# Number of Axes & plot size
ncols=2
nrows=1
figsize=4
wspace=0.5
fig,axs = plt.subplots(nrows=nrows, ncols=ncols,
                       figsize=(ncols*figsize+figsize*wspace*(ncols-1),nrows*figsize))
plt.subplots_adjust(wspace=wspace)
# This produces two Axes objects in a single Figure
print('axes:',axs)

# We can use these Axes objects individually to plot on them
# We need to set show=False so that the Figure is not displayed before we
# finished plotting on all Axes and making all plot adjustments
sc.pl.umap(adata,color='louvain',ax=axs[0],show=False)
# Example zoom-in into a subset of louvain clusters
sc.pl.umap(adata[adata.obs.louvain.isin(['0','3','9']),:],color='S_score',ax=axs[1])
#%%
# In this example we want to show UMAPs of different cell type markers,
# with markers of a single cell type in one row
# and with a different number of markers per cell type (row)

# Marker genes
marker_genes= {
    'B-cell': ['CD79A', 'MS4A1'],
    'Dendritic': ['FCER1A', 'CST3'],
    'Monocytes': ['FCGR3A'],
    'NK': ['GNLY', 'NKG7'],
    'Other': ['IGLL1'],
    'Plasma': ['IGJ'],
    'T-cell': ['CD3D'],
}
# Make Axes
# Number of needed rows and columns (based on the row with the most columns)
nrow=len(marker_genes)
ncol=max([len(vs) for vs in marker_genes.values()])
fig,axs=plt.subplots(nrow,ncol,figsize=(2*ncol,2*nrow))
# Plot expression for every marker on the corresponding Axes object
for row_idx,(cell_type,markers) in enumerate(marker_genes.items()):
    col_idx=0
    for marker in markers:
        ax=axs[row_idx,col_idx]
        sc.pl.umap(adata,color=marker,ax=ax,show=False,frameon=False,s=20)
        # Add cell type as row label - here we simply add it as ylabel of
        # the first Axes object in the row
        if col_idx==0:
            # We disabled axis drawing in UMAP to have plots without background and border
            # so we need to re-enable axis to plot the ylabel
            ax.axis('on')
            ax.tick_params(
                top='off', bottom='off', left='off', right='off',
                labelleft='on', labelbottom='off')
            ax.set_ylabel(cell_type+'\n', rotation=90, fontsize=14)
            ax.set(frame_on=False)
        col_idx+=1
    # Remove unused column Axes in the current row
    while col_idx<ncol:
        axs[row_idx,col_idx].remove()
        col_idx+=1
# Alignment within the Figure
fig.tight_layout()
#%%
rcParams['figure.figsize']=(2,2)
sc.pl.umap(adata,color='bulk_labels')
# Set back to value selected above
rcParams['figure.figsize']=FIGSIZE
with plt.rc_context({'figure.figsize':(5,5)}):
    sc.pl.umap(adata,color='bulk_labels')
#%%
fig,ax=plt.subplots(figsize=(4,4))
sc.pl.umap(adata,color='bulk_labels',ax=ax)
#%%
ncol=2
nrow=1
figsize=3
wspace=1
# Adapt figure size based on number of rows and columns and added space between them
# (e.g. wspace between columns)
fig,axs=plt.subplots(nrow,ncol,figsize=(ncol*figsize+(ncol-1)*wspace*figsize,nrow*figsize))
plt.subplots_adjust(wspace=wspace)
sc.pl.umap(adata,color='louvain',ax=axs[0],show=False)
sc.pl.umap(adata,color='phase',ax=axs[1])
#%%
# Default, legend is overlapping
sc.pl.umap(adata,color=['bulk_labels','phase'])
#%%
# Increase gap size between plots
sc.pl.umap(adata,color=['bulk_labels','phase'],wspace=1)
#%%
# Set title with the title parameter
# Return Axes to further modify the plot
ax=sc.pl.umap(adata,color='bulk_labels', title='Cell type',show=False)
# Modify xlabel
_=ax.set_xlabel('umap1',fontsize=20)
#%%
# Make title italic
ax=sc.pl.umap(adata,color='IGJ', show=False)
_=ax.set_title('IGJ', style='italic')
#%%
# Transparent background and no borders/axis labels with frameon=False
sc.pl.umap(adata,color='bulk_labels',frameon=False)
#%%
dp=sc.pl.dotplot(adata, ['CD79A', 'MS4A1'], 'bulk_labels', show=False)
# All Axes used in dotplot
print('Dotplot axes:',dp)
# Select the Axes object that contains the subplot of interest
ax=dp['mainplot_ax']
# Loop through ticklabels and make them italic
for l in ax.get_xticklabels():
    l.set_style('italic')
    g=l.get_text()
    # Change settings (e.g. color) of certain ticklabels based on their text (here gene name)
    if g =='MS4A1':
        l.set_color('#A97F03')
#%%
# The default ordering of cell cycle phases is alphabetical
# To ensure that the ordering corresponds to cell cycle define order of categories;
# this should include all categories in the corresponding pandas table column
phases=['G1','S','G2M']
adata.obs['phase_ordered']=pd.Categorical(
    values=adata.obs.phase, categories=phases, ordered=True)
sc.pl.umap(adata,color=['phase','phase_ordered'], wspace=0.5)
# This just removes the newly added ordered column from adata as we do not need it below
adata.obs.drop('phase_ordered',axis=1,inplace=True)
#%%
fig=sc.pl.umap(adata,color=['bulk_labels'],return_fig=True)
ax=fig.axes[0]
ax.legend_.set_title('Cell type')
# Change Legend location
ax.legend_.set_bbox_to_anchor((-0.2,-0.7))
#%%
from matplotlib.lines import Line2D
fig=sc.pl.umap(adata,color=['bulk_labels'],return_fig=True)
ax=fig.axes[0]
# Remove original Legend
ax.legend_.remove()
# Make new Legend
l1=ax.legend(

    # Add Legend element for each color group
    handles=[
        # Instead of Line2D we can also use other matplotlib objects, such as Patch, etc.
        Line2D([0], [0], marker='x', color=c,lw=0,
               label=l, markerfacecolor=c, markersize=7)
        # Color groups in adata
        for l,c in zip(
            list(adata.obs.bulk_labels.cat.categories),
            adata.uns['bulk_labels_colors'])],

    # Customize Legend outline

    # Remove background
    frameon=False,
    # Make more Legend columns
    ncols=2,
    # Change location to not overlap with the plot
    bbox_to_anchor=(1,1),
    # Set title
    title='Cell type'
)
#%%
fig,ax=plt.subplots(figsize=(3,3))
sc.pl.umap(adata,color=['bulk_labels'],ax=ax,show=False)

# Encircle part of the plot

# Find location on the plot where circle should be added
location_cells=adata[adata.obs.bulk_labels=='CD56+ NK',:].obsm['X_umap']
x=location_cells[:,0].mean()
y=location_cells[:,1].mean()
size=1.5 # Set circle size
# Plot circle
circle = plt.Circle((x,y), size, color='k', clip_on=False,fill=False)
ax.add_patch(circle)

# Add annother Legend for the mark

# Save the original Legend
l1=ax.get_legend()
l1.set_title('Cell type')
# Make a new Legend for the mark
l2=ax.legend(handles=[Line2D([0],[0],marker='o', color='k',  markerfacecolor='none',
                        markersize=12,markeredgecolor='k',lw=0,label='selected')],
          frameon=False, bbox_to_anchor=(3,1),title='Annotation')
# Add back the original Legend which was overwritten by the new Legend
_=plt.gca().add_artist(l1)
#%%
# Package used for adding well aligned labels on the plot
from adjustText import adjust_text

with plt.rc_context({'figure.figsize':(5,5)}):
    x='means'
    y='dispersions'
    color='is_highly_variable'
    adata.var['is_highly_variable']=adata.var['highly_variable'].astype(bool).astype(str)
    ax=sc.pl.scatter(adata,x=x,y=y,color=color,show=False)
    print('Axes:',ax)
    # Move plot title from Axes to Legend
    ax.set_title('')
    ax.get_legend().set_title(color)

    # Labels

    # Select genes to be labeled
    texts = []
    genes= ['CD79A', 'MS4A1','FCER1A', 'CST3','FCGR3A','GNLY', 'NKG7','IGLL1','IGJ','CD3D']
    for gene in genes:
        # Position of object to be marked
        x_loc=adata.var.at[gene,x]
        y_loc=adata.var.at[gene,y]
        # Text color
        color_point='k'
        texts.append(ax.text(x_loc, y_loc, gene, color=color_point, fontsize=10))

    # Label selected genes on the plot
    _=adjust_text(texts, expand=(2, 2),
        arrowprops=dict(arrowstyle="->",  color='gray',  lw=1), ax=ax)
#%%
sc.pl.umap(adata,color='phase',s=20,
           palette={'S':'tab:cyan','G1':'tab:olive','G2M':'tab:red'})
#%%
# Center palette with vcenter

# Make mock column for plotting, here we use random values from normal distribution
loc=0
adata.obs['normal']=np.random.normal(loc=loc,size=adata.shape[0])

# Center at mean (loc) of the distribution with vcenter parameter
sc.pl.umap(adata,color='normal',cmap='coolwarm',s=20,vcenter=loc)
adata.obs.drop('normal',axis=1,inplace=True)
#%%
# Make symmetric palette with vmin and vmax

# Make mock column for plotting, here we use B cell score
sc.tl.score_genes(adata,['CD79A', 'MS4A1'],score_name='B_cell_score')

# To make a symmetric palette centerd around 0 we set vmax to maximal absolut value and vmin to
# the negative value of maxabs
maxabs=max(abs(adata.obs['B_cell_score']))
sc.pl.umap(adata,color='B_cell_score',cmap='coolwarm',s=20,vmin=-maxabs,vmax=maxabs)
adata.obs.drop('B_cell_score',axis=1,inplace=True)
#%%
# Log-scaled palette

# Make mock column with log-normally distirbuited values
adata.obs['lognormal']=np.random.lognormal(3, 1, adata.shape[0])

# Log scaling of the palette
norm=mcolors.LogNorm()
sc.pl.umap(adata,color='lognormal',s=20, norm=norm)

adata.obs.drop('lognormal',axis=1,inplace=True)
#%%
# Centered non-symmetric palette

# Make mock column for plotting, here we use B cell score
sc.tl.score_genes(adata,['CD79A', 'MS4A1'],score_name='B_cell_score')

# Palette normalization with centering and adapted dynamic range to correspond to
# the distance of vmin and vmax from the cenetr
# Adapted from https://stackoverflow.com/a/50003503
class MidpointNormalize(mcolors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=0, clip=False):
        self.midpoint = midpoint
        mcolors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        value=np.array(value).astype(float)
        normalized_min = max(0.0, 0.5 * (1.0 - abs((self.midpoint - self.vmin) / (self.midpoint - self.vmax))))
        normalized_max = min(1.0, 0.5 * (1.0 + abs((self.vmax - self.midpoint) / (self.midpoint - self.vmin))))
        normalized_mid = 0.5
        x, y = [self.vmin, self.midpoint, self.vmax], [normalized_min, normalized_mid, normalized_max]
        return np.ma.masked_array(np.interp(value, x, y))
# Add padding arround vmin and vmax as Colorbar sets value limits to round numbers below and
# above the vmin and vmax, respectively, which means that they can not be assigned the correct
# color with our nomalisation function that is limited to vmin and vmax
# However, this padding reduces the dynamic range as we set a broad padding and
# then later discard values that are not needed for the rounding up and down
# of the vmin and vmax on the Colorbar, respectively
vmin=adata.obs['B_cell_score'].min()
vmax=adata.obs['B_cell_score'].max()
vpadding=(vmax-vmin)*0.2
norm = MidpointNormalize(
    vmin=vmin-vpadding, vmax=vmax+vpadding,
    midpoint=0)
# Plot umap
fig=sc.pl.umap(adata,color='B_cell_score',cmap='coolwarm',s=20,
               norm=norm,return_fig=True,show=False)
# Adjust Colorbar ylim to be just outside of vmin,vmax and not far outside of this range
# as the padding we set initially may be too broad
cmap_yticklabels=np.array([t._y for t in fig.axes[1].get_yticklabels()])
fig.axes[1].set_ylim(max(cmap_yticklabels[cmap_yticklabels<vmin]),
                     min(cmap_yticklabels[cmap_yticklabels>vmax]))

adata.obs.drop('B_cell_score',axis=1,inplace=True)
#%%
ax=sc.pl.umap(adata,color=['bulk_labels'],groups=['Dendritic'], show=False)

# We can change the 'NA' in the legend that represents all cells outside of the
# specified groups
legend_texts=ax.get_legend().get_texts()
# Find legend object whose text is "NA" and change it
for legend_text in legend_texts:
    if legend_text.get_text()=="NA":
        legend_text.set_text('other cell types')
#%%
# Define dot size for all plot parts
dot_size=40
# Plot all cells as background
ax=sc.pl.umap(adata, show=False,s=dot_size)
# Plot ontop expression of a single cell group by subsetting adata
sc.pl.umap(adata[adata.obs.bulk_labels=='CD19+ B',:],color='IGJ',ax=ax,s=dot_size)        
#%%
# Make two batches in the adata object for the plot example
adata.obs['batch']=['a']*int(adata.shape[0]/2)+['b']*(adata.shape[0]-int(adata.shape[0]/2))

fig,axs=plt.subplots(1,2,figsize=(9,3))
plt.subplots_adjust(wspace=1)
sc.pl.umap(adata,color='batch',ax=axs[0],title='Default ordering',show=False)
# Randomly order cells by making a random index and subsetting AnnData based on it
# Set a random seed to ensure that the cell ordering will be reproducible
np.random.seed(0)
random_indices=np.random.permutation(list(range(adata.shape[0])))
sc.pl.umap(adata[random_indices,:],color='batch',ax=axs[1],title='Random re-ordering')
#%%
# Copy adata not to modify UMAP in the original adata object
adata_temp=adata.copy()
# Loop through different umap parameters, recomputting and replotting UMAP for each of them
for min_dist in [0.1,1,2]:
    for spread in [0.5,1,5]:
        param_str=' '.join(['min_dist =',str(min_dist),'and spread =',str(spread)])
        sc.tl.umap(adata_temp, min_dist=min_dist, spread=spread)
        # Display plot and then immediately close it to ensure that
        # we do not open too many plot windows at once
        g=sc.pl.umap(adata_temp,color=['louvain'],title=param_str,s=40,
                   show=False, return_fig=True)
        display(g)
        plt.close()
del adata_temp
#%%
sc.tl.paga(adata,groups='louvain')
# Distribution of PAGA connectivities for determining the cutting threshold
fig,axs=plt.subplots(1,2,figsize=(6,3))
paga_conn=adata.uns['paga']['connectivities'].toarray().ravel()
a=axs[0].hist(paga_conn,bins=30)
sns.violinplot(paga_conn,ax=axs[1], inner=None)
sns.swarmplot(paga_conn,ax=axs[1],color='k')
thr=0.5
_=axs[1].axhline(thr,c='r')
_=axs[0].axvline(thr,c='r')
#%%
# Compare PAGA with and without prunning
fig,axs=plt.subplots(1,2,figsize=(6,3))
sc.pl.paga(adata,ax=axs[0],title='PAGA',show=False)
sc.pl.paga(adata,ax=axs[1],title='PAGA - prunned',threshold=thr)    
#%%
# Compare UMAP and PAGA layouts
fig,axs=plt.subplots(1,2,figsize=(6,3))
sc.pl.umap(adata,color='louvain',ax=axs[0],show=False,title='UMAP',legend_loc='on data')
sc.pl.paga(adata,ax=axs[1],title='PAGA')
#%%
# Define PAGA positions based on the UMAP layout -
# for each cluster we use the mean of the UMAP positions from the cells in that cluster
pos=pd.DataFrame(adata.obsm['X_umap'],index=adata.obs_names)
pos['group']=adata.obs[adata.uns['paga']['groups']]
pos=pos.groupby('group').mean()

# Plot UMAP in the background
ax=sc.pl.umap(adata,show=False)
# Plot PAGA ontop of the UMAP
sc.pl.paga(adata, color='louvain',threshold=thr,
           node_size_scale=1, edge_width_scale=0.7,
           pos=pos.values,
           random_state=0, ax=ax)
#%%
import scanpy as sc
import pandas as pd
import seaborn as sns
#%%
sc.settings.verbosity = 1             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(3, 3), facecolor='white')
#%%
adata_ref = sc.datasets.pbmc3k_processed()  # this is an earlier version of the dataset from the pbmc3k tutorial
adata = sc.datasets.pbmc68k_reduced()
#%%
var_names = adata_ref.var_names.intersection(adata.var_names)
adata_ref = adata_ref[:, var_names]
adata = adata[:, var_names]
#%%
sc.pp.pca(adata_ref)
sc.pp.neighbors(adata_ref)
sc.tl.umap(adata_ref)
#%%
sc.pl.umap(adata_ref, color='louvain')
#%%
sc.tl.ingest(adata, adata_ref, obs='louvain')
adata.uns['louvain_colors'] = adata_ref.uns['louvain_colors']  # fix colors
sc.pl.umap(adata, color=['louvain', 'bulk_labels'], wspace=0.5)
#%%
adata_concat = adata_ref.concatenate(adata, batch_categories=['ref', 'new'])
adata_concat.obs.louvain = adata_concat.obs.louvain.astype('category')
adata_concat.obs.louvain.cat.reorder_categories(adata_ref.obs.louvain.cat.categories, inplace=True)  # fix category ordering
adata_concat.uns['louvain_colors'] = adata_ref.uns['louvain_colors']  # fix category colors
sc.pl.umap(adata_concat, color=['batch', 'louvain'])
#%%
sc.tl.pca(adata_concat)
!time
sc.external.pp.bbknn(adata_concat, batch_key='batch')  # running bbknn 1.3.6
sc.tl.umap(adata_concat)
sc.pl.umap(adata_concat, color=['batch', 'louvain'])
#%%
# note that this collection of batches is already intersected on the genes
adata_all = sc.read('data/pancreas.h5ad', backup_url='https://www.dropbox.com/s/qj1jlm9w10wmt0u/pancreas.h5ad?dl=1')
adata_all.shape
counts = adata_all.obs.celltype.value_counts()
counts
minority_classes = counts.index[-5:].tolist()        # get the minority classes
adata_all = adata_all[                               # actually subset
    ~adata_all.obs.celltype.isin(minority_classes)]
adata_all.obs.celltype.cat.reorder_categories(       # reorder according to abundance
    counts.index[:-5].tolist(), inplace=True)
sc.pp.pca(adata_all)
sc.pp.neighbors(adata_all)
sc.tl.umap(adata_all)
sc.pl.umap(adata_all, color=['batch', 'celltype'], palette=sc.pl.palettes.vega_20_scanpy)
#%%
!time
sc.external.pp.bbknn(adata_all, batch_key='batch')
#%%
sc.tl.umap(adata_all)
sc.pl.umap(adata_all, color=['batch', 'celltype'])
#%%
adata_ref = adata_all[adata_all.obs.batch == '0']
sc.pp.pca(adata_ref)
sc.pp.neighbors(adata_ref)
sc.tl.umap(adata_ref)
sc.pl.umap(adata_ref, color='celltype')
#%%
adatas = [adata_all[adata_all.obs.batch == i].copy() for i in ['1', '2', '3']]
sc.settings.verbosity = 2  # a bit more logging
for iadata, adata in enumerate(adatas):
    print(f'... integrating batch {iadata+1}')
    adata.obs['celltype_orig'] = adata.obs.celltype  # save the original cell type
    sc.tl.ingest(adata, adata_ref, obs='celltype')
#%%
adata_concat = adata_ref.concatenate(adatas)
adata_concat.obs.celltype = adata_concat.obs.celltype.astype('category')
adata_concat.obs.celltype.cat.reorder_categories(adata_ref.obs.celltype.cat.categories, inplace=True)  # fix category ordering
adata_concat.uns['celltype_colors'] = adata_ref.uns['celltype_colors']  # fix category coloring
sc.pl.umap(adata_concat, color=['batch', 'celltype'])
#%%
adata_query = adata_concat[adata_concat.obs.batch.isin(['1', '2', '3'])]
sc.pl.umap(
    adata_query, color=['batch', 'celltype', 'celltype_orig'], wspace=0.4)
#%%
obs_query = adata_query.obs
conserved_categories = obs_query.celltype.cat.categories.intersection(obs_query.celltype_orig.cat.categories)  # intersected categories
obs_query_conserved = obs_query.loc[obs_query.celltype.isin(conserved_categories) & obs_query.celltype_orig.isin(conserved_categories)]  # intersect categories
obs_query_conserved.celltype.cat.remove_unused_categories(inplace=True)  # remove unused categoriyes
obs_query_conserved.celltype_orig.cat.remove_unused_categories(inplace=True)  # remove unused categoriyes
obs_query_conserved.celltype_orig.cat.reorder_categories(obs_query_conserved.celltype.cat.categories, inplace=True)  # fix category ordering
pd.crosstab(obs_query_conserved.celltype, obs_query_conserved.celltype_orig)
#%%
pd.crosstab(adata_query.obs.celltype, adata_query.obs.celltype_orig)
#%%
sc.tl.embedding_density(adata_concat, groupby='batch')
sc.pl.embedding_density(adata_concat, groupby='batch')
#%%
for batch in ['1', '2', '3']:
    sc.pl.umap(adata_concat, color='batch', groups=[batch])
    
#%%
#From scanpy API section
#Also see [Data integration]. Note that a simple batch correction method is available via pp.regress_out(). Checkout scanpy.external for more.
pp.combat(adata[, key, covariates, inplace])

#ComBat function for batch effect correction [Johnson07] [Leek12] [Pedersen12].
