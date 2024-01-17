#!/usr/bin/env python3
# -*- coding: utf-8 -*-
@ Author Jake D. Lehle Ph.D.

#%%
# Loading screen chunk

#Imports from python
import os
os.chdir("/work/sdz852/goldenboy/SC")
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sb
import scipy as sp
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri
import anndata2ri
anndata2ri.activate()
rpy2.robjects.numpy2ri.activate()
anndata2ri.scipy2ri.activate()

# We need a way to import in R packages using rpy2 in a python env
import rpy2
from rpy2.robjects.packages import importr

# Now we use the importr() function to import in R packages and assign them a variable name in python
# The first teo packages you have to alyways bring in base and utility 
base = importr("base")
utils = importr("utils")
scran = importr("scran")
RColorBrewer = importr("RColorBrewer")
slingshot = importr("slingshot")
gam = importr("gam")
clusterExperiment = importr("clusterExperiment")
ggplot2 = importr("ggplot2")
plyr = importr("plyr")
MAST = importr("MAST")

from matplotlib import colors
from matplotlib import rcParams

#monocle = importr("monocle")

sc.logging.print_versions()


#%%
#This is just to set up how the figures will display
plt.rcParams['figure.figsize']=(8,8) #rescale figures
sc.settings.verbosity = 3
#sc.set_figure_params(dpi=200, dpi_save=300)
sc.logging.print_versions()
# Set up data loading
#Data files
# Different sample groups this comes from the naming convention of the files external of this python script.
sample_file_strings = ['E18_S4_L002_/', 'E18_S5_L002_/']
sample_py_strings = ['E18_E18R1', 'E18_E18R2']
# Last number in the sample id that is unique all other id numbers are the same. Wonder if this would work with more than 10 groups.
sample_id_strings = ['1', '2']
#Path to the files. Look how the file ends with the id beginning of the file and the only thing missing is the sample id value that is unique for each sample.
file_base = '/work/sdz852/goldenboy/SC/fastq/E18.5/'
#different file types tsv = files with feature and barcode sequences corresponding to row and column indices respectively. mtx
data_file_end = 'matrix.mtx'
barcode_file_end = 'barcodes.tsv'
gene_file_end = "features.tsv"
cc_genes_file = '../Macosko_cell_cycle_genes.txt'

# First data set load & annotation
#Parse Filenames
sample = sample_file_strings.pop(0)
sample_py = sample_py_strings.pop(0)
sample_id = sample_id_strings.pop(0)
data_file = file_base+sample+sample+'outs/filtered_feature_bc_matrix/'+data_file_end
barcode_file = file_base+sample+sample+'outs/filtered_feature_bc_matrix/'+barcode_file_end
gene_file = file_base+sample+sample+'outs/filtered_feature_bc_matrix/'+gene_file_end
#Load data
adata = sc.read(data_file, cache=True)
adata = adata.transpose()
adata.X = adata.X.toarray()
barcodes = pd.read_csv(barcode_file, header=None, sep='\t')
genes = pd.read_csv(gene_file, header=None, sep='\t')
#Annotate data
barcodes.rename(columns={0:'barcode'}, inplace=True)
barcodes.set_index('barcode', inplace=True)
adata.obs = barcodes
adata.obs['sample'] = [sample_py]*adata.n_obs
adata.obs['condition'] = [sample_py.split("_")[0]]*adata.n_obs
adata.obs['rep'] = [sample_py.split("_")[1]]*adata.n_obs
genes.rename(columns={0:'gene_id', 1:'gene_symbol'}, inplace=True)
genes.set_index('gene_symbol', inplace=True)
adata.var = genes
adata.var_names_make_unique()


# Loop to load rest of data sets
for i in range(len(sample_file_strings)):
    #Parse Filenames
    sample = sample_file_strings[i]
    sample_py = sample_py_strings[i]
    sample_id = sample_id_strings[i]
    data_file = file_base+sample+sample+'outs/filtered_feature_bc_matrix/'+data_file_end
    barcode_file = file_base+sample+sample+'outs/filtered_feature_bc_matrix/'+barcode_file_end
    gene_file = file_base+sample+sample+'outs/filtered_feature_bc_matrix/'+gene_file_end
    #Load data
    adata_tmp = sc.read(data_file, cache=True)
    adata_tmp = adata_tmp.transpose()
    adata_tmp.X = adata_tmp.X.toarray()
    barcodes_tmp = pd.read_csv(barcode_file, header=None, sep='\t')
    genes_tmp = pd.read_csv(gene_file, header=None, sep='\t')
    #Annotate data
    barcodes_tmp.rename(columns={0:'barcode'}, inplace=True)
    barcodes_tmp.set_index('barcode', inplace=True)
    adata_tmp.obs = barcodes_tmp
    adata_tmp.obs['sample'] = [sample_py]*adata_tmp.n_obs
    adata_tmp.obs['condition'] = [sample_py.split("_")[0]]*adata_tmp.n_obs
    adata_tmp.obs['rep'] = [sample_py.split("_")[1]]*adata_tmp.n_obs
    genes_tmp.rename(columns={0:'gene_id', 1:'gene_symbol'}, inplace=True)
    genes_tmp.set_index('gene_symbol', inplace=True)
    adata_tmp.var = genes_tmp
    adata_tmp.var_names_make_unique()
    # Concatenate to main adata object
    adata = adata.concatenate(adata_tmp, batch_key='sample_id')
    if 'gene_id-1' in adata.var.columns:
        adata.var['gene_id'] = adata.var['gene_id-1']
        adata.var.drop(columns = ['gene_id-1', 'gene_id-0'], inplace=True)
    adata.obs.drop(columns=['sample_id'], inplace=True)
    adata.obs_names = [c.split("-")[0] for c in adata.obs_names]
    adata.obs_names_make_unique(join='_')


print('out of loop')
#%%
#Assign variable names and gene id columns
#adata.var
#adata.var_names = [g.split("_")[1] for g in adata.var_names]
#adata.var['gene_id'] = [g.split("_")[1] for g in adata.var['gene_id']]

# Check data set annotations
print(adata.obs['sample'].value_counts())
print(adata.obs['condition'].value_counts())
print(adata.obs['rep'].value_counts())

# Checking the total size of the data set the first number is the number of single cells in the data set the second is the number of differentially expressed genes.
adata.shape

# Quality control - calculate QC covariates
adata.obs['n_counts'] = adata.X.sum(1) #first sum the counts
adata.obs['log_counts'] = np.log(adata.obs['n_counts']) #now take the log of the n_counts
adata.obs['n_genes'] = (adata.X > 0).sum(1) #find the number of expressed genes
mt_gene_mask = [gene.startswith('mt-') for gene in adata.var_names] #add the gene names
adata.obs['mt_frac'] = adata.X[:, mt_gene_mask].sum(1)/adata.obs['n_counts'] #filter the gene names using the expression > 0

#%%
# Quality control - plot QC metrics
#Sample quality plots
t1 = sc.pl.violin(adata, 'n_counts', groupby='sample', size=2, log=True, cut=0) #violin plot of the log caounts of reads for each sample. You can see that worst sample is the jej_M2 due to having the poorest quality of reads. However, there are still enough reads to process these cells. 
t2 = sc.pl.violin(adata, 'mt_frac', groupby='sample') #This shows the fraction of mitochondrial reads (MT frac) this should always be below 20-25%. Again you can see the jej_M2 is again the wost sample in the group but good enough to continue processing.

#Data quality summary plots
p1 = sc.pl.scatter(adata, 'n_counts', 'n_genes', color='mt_frac') #scatterplot of the number of genes vs the number of counts with MT fraction information. You can see in the lower portion of the graph there are cells with high read counts but not many genes the first throught would be to remove these cells as dying outliers but you can see that these cells are purple and have a low MT fraction. We will still probably filter out some of these in the future since they will be difficult to annotate (1000-4000 counts and <~500 genes) just keep this in mind right now.
# This looks good the heat map of the points indicate the sample prep went very wel and only the cells with the fewsest number of gene counts had mitochondrial genes present which is another indaication they are dying and need to be filtered away.
p2 = sc.pl.scatter(adata[adata.obs['n_counts']<10000], 'n_counts', 'n_genes', color='mt_frac')
#Thresholding decision: counts. Check to see for the updated use of distplot which will be discontinued in the future.
p3 = sb.distplot(adata.obs['n_counts'], kde=False)
plt.show()
p4 = sb.distplot(adata.obs['n_counts'][adata.obs['n_counts']<10000], kde=False, bins=60)
plt.show()
p5 = sb.distplot(adata.obs['n_counts'][adata.obs['n_counts']>10000], kde=False, bins=60)
plt.show()
#Thresholding decision: genes
p6 = sb.distplot(adata.obs['n_genes'], kde=False, bins=60)
plt.show()
p7 = sb.distplot(adata.obs['n_genes'][adata.obs['n_genes']<1000], kde=False, bins=60)
plt.show()
# Filter cells according to identified QC thresholds:
print('Total number of cells: {:d}'.format(adata.n_obs))
sc.pp.filter_cells(adata, min_counts = 1000)
print('Number of cells after min count filter: {:d}'.format(adata.n_obs))
sc.pp.filter_cells(adata, max_counts = 40000)
print('Number of cells after max count filter: {:d}'.format(adata.n_obs))
adata = adata[adata.obs['mt_frac'] < 0.2]
print('Number of cells after MT filter: {:d}'.format(adata.n_obs))
sc.pp.filter_cells(adata, min_genes = 700)
print('Number of cells after gene filter: {:d}'.format(adata.n_obs))
#Filter genes:
print('Total number of genes: {:d}'.format(adata.n_vars))
# Min 20 cells - filters out 0 count genes
sc.pp.filter_genes(adata, min_cells=20)
print('Number of genes after cell filter: {:d}'.format(adata.n_vars))
#Perform a clustering for scran normalization in clusters
adata_pp = adata.copy()
sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=1e6)
sc.pp.log1p(adata_pp)
sc.pp.pca(adata_pp, n_comps=15)
adata_pp.var_names
sc.pp.neighbors(adata_pp)
sc.tl.louvain(adata_pp, key_added='groups', resolution=0.5) 
#Preprocess variables for scran normalization in R. First set the counts file that has a matrix attached to it as a R object by saving it to the global environment so R fucntions can find it.

robjects.globalenv['data_mat'] = adata.X.T # This is our SingleCellExperiment object containing a count matrix.
adata_pp.obs['groups']
robjects.globalenv['input_groups'] = adata_pp.obs['groups'] # Here are the groups that match up with the count matrix.
print(robjects.r.ls())

size_factors = robjects.r(f'sizeFactors(computeSumFactors(SingleCellExperiment(list(counts=data_mat)), clusters = input_groups, min.mean = 0.1))')
#Pause! is the size_factors object in the R or python environment now??? remeber think about what rpy2 functions allow you to do in the 'python' environment...

#Delete adata_pp
del adata_pp
# Visualize the estimated size factors
adata.obs['size_factors'] = size_factors

sc.pl.scatter(adata, 'size_factors', 'n_counts')
sc.pl.scatter(adata, 'size_factors', 'n_genes')

sb.displot(size_factors, bins=50, kde=False)

#Keep the count data in a counts layer
adata.layers["counts"] = adata.X.copy()
#Normalize adata by first adding the size factors object to the adata.obs.
#adata.obs['size_factors'] = size_factors
adata.X /= adata.obs['size_factors'].values[:,None]
sc.pp.log1p(adata)
# Store the full data set in 'raw' as log-normalised data for statistical testing
adata.raw = adata
# ComBat batch correction
sc.pp.combat(adata, key='sample')
sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=4000)
print('\n','Number of highly variable genes: {:d}'.format(np.sum(adata.var['highly_variable'])))
sc.pl.highly_variable_genes(adata)
# Calculate the visualizations
sc.pp.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack')
sc.pp.neighbors(adata)
sc.tl.tsne(adata, n_jobs=50) #Note n_jobs works for MulticoreTSNE, but not regular implementation)
sc.tl.umap(adata)
sc.tl.diffmap(adata)
sc.tl.draw_graph(adata) #install fa2 $pip3 install fa2 

sc.pl.pca_scatter(adata, color='n_counts')
sc.pl.tsne(adata, color='n_counts')
sc.pl.umap(adata, color='n_counts')
sc.pl.diffmap(adata, color='n_counts', components=['1,2','1,3'])
sc.pl.draw_graph(adata, color='n_counts')
#%%
#Score cell cycle and visualize the effect:
cc_genes = pd.read_table(cc_genes_file, delimiter='\t') #The file can be found on the single-cell-tutorial github repository, or be taken from the supplementary material of the paper. Save it in the Documents folder or SC_RNA_SEQ file that has had ../ applied to it.
s_genes = cc_genes['S'].dropna()
g2m_genes = cc_genes['G2.M'].dropna()

s_genes_mm = [gene.lower().capitalize() for gene in s_genes]
g2m_genes_mm = [gene.lower().capitalize() for gene in g2m_genes]

s_genes_mm_ens = adata.var_names[np.in1d(adata.var_names, s_genes_mm)]
g2m_genes_mm_ens = adata.var_names[np.in1d(adata.var_names, g2m_genes_mm)]

sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes_mm_ens, g2m_genes=g2m_genes_mm_ens)
sc.pl.umap(adata, color=['S_score', 'G2M_score'], use_raw=False)
sc.pl.umap(adata, color='phase', use_raw=False)
#%%
# Perform clustering - using highly variable genes
sc.tl.louvain(adata, key_added='louvain_r1')

#%%
sc.tl.louvain(adata, resolution=0.75, key_added='louvain_r0.75', random_state=10) #play around with resolution to get the 9 total sub-groups
sc.tl.louvain(adata, resolution=0.95, key_added='louvain_r0.95', random_state=10) #play around with resolution to get the 9 total sub-groups
adata.obs['louvain_r0.95'].value_counts()
adata.obs['louvain_r0.75'].value_counts()
#Visualize the clustering and how this is reflected by different technical covariates
sc.pl.umap(adata, color=['louvain_r0.75', 'louvain_r1'], palette=sc.pl.palettes.vega_20) #module 'scanpy.plotting.palettes' has no attribute 'default_64'
sc.pl.umap(adata, color=['rep', 'n_counts'])
sc.pl.umap(adata, color=['log_counts', 'mt_frac'])
#Calculate marker genes this is a very important section the key added will be a dictionary that you can use to find raw relative gene expression values and plot figures with. The keys are stored in the adata.uns object
sc.tl.rank_genes_groups(adata, groupby='louvain_r0.75', key_added='rank_genes_r0.75')
sc.tl.rank_genes_groups(adata, groupby='louvain_r0.75', n_genes=adata.raw.shape[1])

type(adata.uns['rank_genes_r0.75'])

#Plot marker genes plotting
sc.pl.rank_genes_groups(adata, key='rank_genes_r0.75', groups=['4','10'], fontsize=12)
sc.pl.rank_genes_groups(adata, key='rank_genes_r0.75', groups=['5', '6', '7', '8', '9', '10'], fontsize=12)

#Known marker genes:
cell_death = pd.read_csv('/work/sdz852/goldenboy/SC/Cell_Death.csv')

dna_repair = pd.read_csv('/work/sdz852/goldenboy/SC/DNA_Repair.csv')
dna_repair

#Here I will make a new dictionary with known marker genes specifically for all the DNA and Cell death genes that came from Danny's project
marker_genes = dict()

marker_genes['Cell Death'] = cell_death['Final list cell death'].tolist()
marker_genes['DNA Repair'] = dna_repair['Final list cell repair'].tolist()

marker_genes['DNA Repair']

#To find which group has the highest expression of marker gene expression I can quickly find the overlap between the ranked expression and the marker genes
cell_annotation = sc.tl.marker_gene_overlap(adata, marker_genes, key='rank_genes_r0.75')
cell_annotation
# This looks good but let's normalize everything and display it as a heatmap
cell_annotation_norm = sc.tl.marker_gene_overlap(adata, marker_genes, key='rank_genes_r0.75', normalize='reference')
sb.heatmap(cell_annotation_norm, cbar=False, annot=True)
#Awesome it looks likt there is one distinct cluster that has higher median expression of some genes in the DNA marker dict we make some additional plots of this with dot plot

#Define a nice colour map for gene expression
colors2 = plt.cm.Reds(np.linspace(0, 1, 128))
colors3 = plt.cm.Greys_r(np.linspace(0.7,0.8,20))
colorsComb = np.vstack([colors3, colors2])
mymap = colors.LinearSegmentedColormap.from_list('my_colormap', colorsComb)

# How to do individual feature maps for #Defa24 #Tff3
sc.pl.umap(adata, color='Defa24', use_raw=False, color_map=mymap)
sc.pl.umap(adata, color='Tff3', use_raw=False, color_map=mymap)



# Check expression of enterocyte markers
#Collate all known markers and get the gene IDs in the data set
ids_cd = np.in1d(adata.var_names, marker_genes['Cell Death'])
ids_dr = np.in1d(adata.var_names, marker_genes['DNA Repair'])
ids_cd
ids_cddr = np.logical_or(ids_cd, ids_dr)

#Calculate the mean expression of known markers for each cell
adata.obs['Cell_Death_and_DNA_Repair_marker_expr'] = adata.X[:,ids_cddr].mean(1)
adata.obs['DNA_Repair_marker_expr'] = adata.X[:,ids_dr].mean(1)
adata.obs['DNA_Repair_marker_expr']
#Plot expression
sc.pl.violin(adata, 'Cell_Death_and_DNA_Repair_marker_expr', groupby='louvain_r0.75', inner='quartile')
sc.pl.violin(adata, 'DNA_Repair_marker_expr', groupby='louvain_r0.75', inner="quart")
sc.pl.umap(adata, color='DNA_Repair_marker_expr', color_map=mymap)

#However the previous section only shows the mean value all expression of all genes in the marker gene set what about expression of indiviudal genes.
#We can look at that using dotplots
#First I will make a var_names object to store the marker genes
var_names = {'DNA Repair': marker_genes['DNA Repair']}
var_names
sc.pl.rank_genes_groups_dotplot(adata, var_names=var_names, standard_scale='var')

# Another was to do this is with the dotplot object functions 
sc.pl.dotplot(adata, marker_genes['DNA Repair'], 'louvain_r0.75', dendrogram=True, standard_scale='var')
dna_repair_missing = ['0610007p08rik', '1110054o05rik', '1810011o10rik', '5730590g19rik', 'Apitd1', 'Bre', 'C77370', 'Ctgf', 'Cyr61', 'Fam175a', 'Fam176a', 'Fgf4', 'Fndc1', 'H2afx', 'H47', 'Lrdd', 'Mkl1', 'Mmd', 'Mum1', 'Phf17', 'Pkm2', 'Pold4', 'Rad51l1', 'Slc9a3r1', 'Supt16h', 'Tcfap4', 'Tdgf1', 'Tmem173', 'Wdr92']
dna_repair = dna_repair[-dna_repair['Final list cell repair'].isin(dna_repair_missing)]


sc.pl.DotPlot(adata, marker_genes['DNA Repair'], groupby='louvain_r0.75', standard_scale='var').show()
dp = sc.pl.dotplot(adata, marker_genes['DNA Repair'], 'louvain_r0.75', dendrogram=True, standard_scale='var', return_fig=True)
dp.add_totals().show()
list(dp)
print(dp)
axes_dict = dp.get_axes()
print(axes_dict)

sc.pl.dotplot(adata, marker_genes['DNA Repair'], 'louvain_r0.75', dendrogram=False, standard_scale='var')



dna_repair_missing = ['1810011o10rik', 'Bre', 'Ctgf', 'Cyr61', 'Fam176a', 'Fgf4', 'Fndc1', 'H47', 'Lrdd', 'Mkl1', 'Mmd', 'Phf17', 'Pkm2', 'Slc9a3r1', 'Tcfap4', 'Tdgf1', 'Tmem173', 'Wdr92']





marker_matches = sc.tl.marker_gene_overlap(adata, marker_genes['DNA Repair'])



