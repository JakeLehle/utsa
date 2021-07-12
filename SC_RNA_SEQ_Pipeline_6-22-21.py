#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 29 17:10:50 2021

@author: jake
"""
#%%
import os
os.chdir("/home/jake/Documents/SC_RNA_SEQ")
#import packages
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
# We need a way to import in R packages using rpy2 so first we import in rpy2 and Ipython
import rpy2.ipython
from rpy2.robjects.packages import importr
# We now use the importr() function to import in R packages and assign them a variable name in python
# The first two packages you have to always bring in are the base are and utility packages that are part of the base R
base = importr('base') 
utils = importr('utils')
scran = importr('scran')
from matplotlib import colors
from matplotlib import rcParams
#So far this is all you need to make the adata object file. Still can't convert it into a SCE file ... I think it needs the key or the barcode added to the file or to have numpy arays converted to something more simple before they can be converteding into an R matrix . There also a possibility it has to do with the colmns and the rows being inverted between the teo data types bt that should be converted automatically for py2rpy.
#%%
import scipy as sp




import combat as c
from gprofiler import gprofiler
import warnings

from rpy2.rinterface import RRuntimeWarning
from gprofiler import gprofiler
import rpy2
import rpy2.rinterface_lib.callbacks
import logging

import anndata2ri

# Ignore R warning messages
#Note: this can be commented out to get more verbose R output
rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)
# Automatically convert rpy2 outputs to pandas dataframes
#pandas2ri.activate()
#anndata2ri.activate()


#What if you need to install a package that isn't already installed in R? Leave this section commented out unless you want to install something.
# Here we use reassign the rpy2.robjects.packages as the function rpackages now
#>>> import rpy2.robjects.packages as rpackages
# We tell python which package we would like to install
#>>> monocle = rpackages.importr('monocle')
# We then select the CRAN mirror to install the package from
#>>> utils.chooseCRANmirror(ind=1)
# Now lets load all of our R packages

RColorBrewer = importr('RColorBrewer')
slingshot = importr('slingshot')
monocle = importr('monocle')
gam = importr('gam')
clusterExperiment = importr('clusterExperiment')
ggplot2 = importr('ggplot2')
plyr = importr('plyr')
MAST = importr('MAST')
basilisk = importr('basilisk')
scRNAseq = importr('scRNAseq')
zellkonverter = importr('zellkonverter')
#%%
#This is just to set up how the figures will display
plt.rcParams['figure.figsize']=(8,8) #rescale figures
sc.settings.verbosity = 3
#sc.set_figure_params(dpi=200, dpi_save=300)
#sc.logging.print_versions()
# Set up data loading
#Data files
# Different sample groups.
sample_strings = ['Duo_M1', 'Duo_M2', 'Jej_M1', 'Jej_M2', 'Il_M1', 'Il_M2']
# Last number in the sample id that is unique all other id numbers are the same. Wonder if this would work with more than 10 groups.
sample_id_strings = ['3', '4', '5', '6', '7', '8']
#Path to the files. Look how the file ends with the id beginning of the file and the only thing missing is the sample id value that is unique for each sample.
file_base = '/home/jake/Documents/SC_RNA_SEQ/data/Haber-et-al_mouse-intestinal-epithelium/GSE92332_RAW/GSM283657'
exp_string = '_Regional_'
#different file types tsv = files with feature and barcode sequences corresponding to row and column indices respectively. mtx
data_file_end = '_matrix.mtx'
barcode_file_end = '_barcodes.tsv'
gene_file_end = '_genes.tsv'
cc_genes_file = '../Macosko_cell_cycle_genes.txt'

# First data set load & annotation
#Parse Filenames
sample = sample_strings.pop(0)
sample_id = sample_id_strings.pop(0)
data_file = file_base+sample_id+exp_string+sample+data_file_end
barcode_file = file_base+sample_id+exp_string+sample+barcode_file_end
gene_file = file_base+sample_id+exp_string+sample+gene_file_end
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
adata.obs['sample'] = [sample]*adata.n_obs
adata.obs['region'] = [sample.split("_")[0]]*adata.n_obs
adata.obs['donor'] = [sample.split("_")[1]]*adata.n_obs
genes.rename(columns={0:'gene_id', 1:'gene_symbol'}, inplace=True)
genes.set_index('gene_symbol', inplace=True)
adata.var = genes
adata.var_names_make_unique()


# Loop to load rest of data sets
for i in range(len(sample_strings)):
    #Parse Filenames
    sample = sample_strings[i]
    sample_id = sample_id_strings[i]
    data_file = file_base+sample_id+exp_string+sample+data_file_end
    barcode_file = file_base+sample_id+exp_string+sample+barcode_file_end
    gene_file = file_base+sample_id+exp_string+sample+gene_file_end
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
    adata_tmp.obs['sample'] = [sample]*adata_tmp.n_obs
    adata_tmp.obs['region'] = [sample.split("_")[0]]*adata_tmp.n_obs
    adata_tmp.obs['donor'] = [sample.split("_")[1]]*adata_tmp.n_obs
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
adata.var_names = [g.split("_")[1] for g in adata.var_names]
adata.var['gene_id'] = [g.split("_")[1] for g in adata.var['gene_id']]

# Check data set annotations
print(adata.obs['region'].value_counts())
print(adata.obs['donor'].value_counts())
print(adata.obs['sample'].value_counts())

# Checking the total size of the data set the first number is the number of single cells in the data set the second is the number of differentially expressed genes.
adata.shape

# Quality control - calculate QC covariates
adata.obs['n_counts'] = adata.X.sum(1) #first sum the counts
adata.obs['log_counts'] = np.log(adata.obs['n_counts']) #now take the log of the n_counts
adata.obs['n_genes'] = (adata.X > 0).sum(1) #find the number of expressed genes
mt_gene_mask = [gene.startswith('mt-') for gene in adata.var_names] #add the gene names
adata.obs['mt_frac'] = adata.X[:, mt_gene_mask].sum(1)/adata.obs['n_counts'] #filter the gene names using the expression > 0

# Quality control - plot QC metrics
#Sample quality plots
t1 = sc.pl.violin(adata, 'n_counts', groupby='sample', size=2, log=True, cut=0) #violin plot of the log caounts of reads for each sample. You can see that worst sample is the jej_M2 due to having the poorest quality of reads. However, there are still enough reads to process these cells. 
t2 = sc.pl.violin(adata, 'mt_frac', groupby='sample') #This shows the fraction of mitochondrial reads (MT frac) this should always be below 20-25%. Again you can see the jej_M2 is again the wost sample in the group but good enough to continue processing.

#Data quality summary plots
p1 = sc.pl.scatter(adata, 'n_counts', 'n_genes', color='mt_frac') #scatterplot of the number of genes vs the number of counts with MT fraction information. You can see in the lower portion of the graph there are cells with high read counts but not many genes the first throught would be to remove these cells as dying outliers but you can see that these cells are purple and have a low MT fraction. We will still probably filter out some of these in the future since they will be difficult to annotate (1000-4000 counts and <~500 genes) just keep this in mind right now.
#Its interesting there are cells in the main clout of points on the left that have high gene counts but also higher fraction of mitochondrial counts. These are actually the cells under stress or dying. When cells are dying there is less mRNA in the nuclues leading to lower overall counts and thus a high fraction of MT mRNA. But cells with high counts high genes and high MT fraction indicate actual biologically relevant cells in the sample with high mitochondrial activity.
p2 = sc.pl.scatter(adata[adata.obs['n_counts']<10000], 'n_counts', 'n_genes', color='mt_frac')
#Thresholding decision: counts. Check to see for the updated use of distplot which will be discontinued in the future.
p3 = sb.distplot(adata.obs['n_counts'], kde=False)
plt.show()
p4 = sb.distplot(adata.obs['n_counts'][adata.obs['n_counts']<4000], kde=False, bins=60)
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
sc.pp.filter_cells(adata, min_counts = 1500)
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
sc.tl.tsne(adata, n_jobs=12) #Note n_jobs works for MulticoreTSNE, but not regular implementation)
sc.tl.umap(adata)
sc.tl.diffmap(adata)
sc.tl.draw_graph(adata) #install fa2 $pip3 install fa2 

sc.pl.pca_scatter(adata, color='n_counts')
sc.pl.tsne(adata, color='n_counts')
sc.pl.umap(adata, color='n_counts')
sc.pl.diffmap(adata, color='n_counts', components=['1,2','1,3'])
sc.pl.draw_graph(adata, color='n_counts')

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
# Perform clustering - using highly variable genes
sc.tl.louvain(adata, key_added='louvain_r1')

#%%
sc.tl.louvain(adata, resolution=0.75, key_added='louvain_r0.75', random_state=10) #play around with resolution to get the 9 total sub-groups

adata.obs['louvain_r0.75'].value_counts()
#Visualize the clustering and how this is reflected by different technical covariates
sc.pl.umap(adata, color=['louvain_r1', 'louvain_r0.75'], palette=sc.pl.palettes.vega_20) #module 'scanpy.plotting.palettes' has no attribute 'default_64'
sc.pl.umap(adata, color=['region', 'n_counts'])
sc.pl.umap(adata, color=['log_counts', 'mt_frac'])
#Calculate marker genes
sc.tl.rank_genes_groups(adata, groupby='louvain_r0.75', key_added='rank_genes_r0.75')
#Plot marker genes
sc.pl.rank_genes_groups(adata, key='rank_genes_r0.75', groups=['0','1','2'], fontsize=12)
sc.pl.rank_genes_groups(adata, key='rank_genes_r0.75', groups=['3','4','5'], fontsize=12)
sc.pl.rank_genes_groups(adata, key='rank_genes_r0.75', groups=['6', '7', '8'], fontsize=12)
#Known marker genes:
marker_genes = dict()
marker_genes['Stem'] = ['Lgr5', 'Ascl2', 'Slc12a2', 'Axin2', 'Olfm4', 'Gkn3']
marker_genes['Enterocyte (Proximal)'] = ['Gsta1','Rbp2','Adh6a','Apoa4','Reg3a','Creb3l3','Cyp3a13','Cyp2d26','Ms4a10','Ace','Aldh1a1','Rdh7','H2-Q2', 'Hsd17b6','Gstm3','Gda','Apoc3','Gpd1','Fabp1','Slc5a1','Mme','Cox7a1','Gsta4','Lct','Khk','Mttp','Xdh','Sult1b1', 'Treh','Lpgat1','Dhrs1','Cyp2c66','Ephx2','Cyp2c65','Cyp3a25','Slc2a2','Ugdh','Gstm6','Retsat','Ppap2a','Acsl5', 'Cyb5r3','Cyb5b','Ckmt1','Aldob','Ckb','Scp2','Prap1']
marker_genes['Enterocyte (Distal)'] = ['Tmigd1','Fabp6','Slc51b','Slc51a','Mep1a','Fam151a','Naaladl1','Slc34a2','Plb1','Nudt4','Dpep1','Pmp22','Xpnpep2','Muc3','Neu1','Clec2h','Phgr1','2200002D01Rik','Prss30','Cubn','Plec','Fgf15','Crip1','Krt20','Dhcr24','Myo15b','Amn','Enpep','Anpep','Slc7a9','Ocm','Anxa2','Aoc1','Ceacam20','Arf6','Abcb1a','Xpnpep1','Vnn1','Cndp2','Nostrin','Slc13a1','Aspa','Maf','Myh14']
marker_genes['Goblet'] = ['Agr2', 'Fcgbp', 'Tff3', 'Clca1', 'Zg16', 'Tpsg1', 'Muc2', 'Galnt12', 'Atoh1', 'Rep15', 'S100a6', 'Pdia5', 'Klk1', 'Pla2g10', 'Spdef', 'Lrrc26', 'Ccl9', 'Bace2', 'Bcas1', 'Slc12a8', 'Smim14', 'Tspan13', 'Txndc5', 'Creb3l4', 'C1galt1c1', 'Creb3l1', 'Qsox1', 'Guca2a', 'Scin', 'Ern2', 'AW112010', 'Fkbp11', 'Capn9', 'Stard3nl', 'Slc50a1', 'Sdf2l1', 'Hgfa', 'Galnt7', 'Hpd', 'Ttc39a', 'Tmed3', 'Pdia6', 'Uap1', 'Gcnt3', 'Tnfaip8', 'Dnajc10', 'Ergic1', 'Tsta3', 'Kdelr3', 'Foxa3', 'Tpd52', 'Tmed9', 'Spink4', 'Nans', 'Cmtm7', 'Creld2', 'Tm9sf3', 'Wars', 'Smim6', 'Manf', 'Oit1', 'Tram1', 'Kdelr2', 'Xbp1', 'Serp1', 'Vimp', 'Guk1', 'Sh3bgrl3', 'Cmpk1', 'Tmsb10', 'Dap', 'Ostc', 'Ssr4', 'Sec61b', 'Pdia3', 'Gale', 'Klf4', 'Krtcap2', 'Arf4', 'Sep15', 'Ssr2', 'Ramp1', 'Calr', 'Ddost']
marker_genes['Paneth'] = ['Gm15284', 'AY761184', 'Defa17', 'Gm14851', 'Defa22', 'Defa-rs1', 'Defa3', 'Defa24', 'Defa26', 'Defa21', 'Lyz1', 'Gm15292', 'Mptx2', 'Ang4']
marker_genes['Enteroendocrine'] = ['Chgb', 'Gfra3', 'Cck', 'Vwa5b2', 'Neurod1', 'Fev', 'Aplp1', 'Scgn', 'Neurog3', 'Resp18', 'Trp53i11', 'Bex2', 'Rph3al', 'Scg5', 'Pcsk1', 'Isl1', 'Maged1', 'Fabp5', 'Celf3', 'Pcsk1n', 'Fam183b', 'Prnp', 'Tac1', 'Gpx3', 'Cplx2', 'Nkx2-2', 'Olfm1', 'Vim', 'Rimbp2', 'Anxa6', 'Scg3', 'Ngfrap1', 'Insm1', 'Gng4', 'Pax6', 'Cnot6l', 'Cacna2d1', 'Tox3', 'Slc39a2', 'Riiad1']
marker_genes['Tuft'] = ['Alox5ap', 'Lrmp', 'Hck', 'Avil', 'Rgs13', 'Ltc4s', 'Trpm5', 'Dclk1', 'Spib', 'Fyb', 'Ptpn6', 'Matk', 'Snrnp25', 'Sh2d7', 'Ly6g6f', 'Kctd12', '1810046K07Rik', 'Hpgds', 'Tuba1a', 'Pik3r5', 'Vav1', 'Tspan6', 'Skap2', 'Pygl', 'Ccdc109b', 'Ccdc28b', 'Plcg2', 'Ly6g6d', 'Alox5', 'Pou2f3', 'Gng13', 'Bmx', 'Ptpn18', 'Nebl', 'Limd2', 'Pea15a', 'Tmem176a', 'Smpx', 'Itpr2', 'Il13ra1', 'Siglecf', 'Ffar3', 'Rac2', 'Hmx2', 'Bpgm', 'Inpp5j', 'Ptgs1', 'Aldh2', 'Pik3cg', 'Cd24a', 'Ethe1', 'Inpp5d', 'Krt23', 'Gprc5c', 'Reep5', 'Csk', 'Bcl2l14', 'Tmem141', 'Coprs', 'Tmem176b', '1110007C09Rik', 'Ildr1', 'Galk1', 'Zfp428', 'Rgs2', 'Inpp5b', 'Gnai2', 'Pla2g4a', 'Acot7', 'Rbm38', 'Gga2', 'Myo1b', 'Adh1', 'Bub3', 'Sec14l1', 'Asah1', 'Ppp3ca', 'Agt', 'Gimap1', 'Krt18', 'Pim3', '2210016L21Rik', 'Tmem9', 'Lima1', 'Fam221a', 'Nt5c3', 'Atp2a3', 'Mlip', 'Vdac3', 'Ccdc23', 'Tmem45b', 'Cd47', 'Lect2', 'Pla2g16', 'Mocs2', 'Arpc5', 'Ndufaf3']
cell_annotation = sc.tl.marker_gene_overlap(adata, marker_genes, key='rank_genes_r0.75')
cell_annotation

cell_annotation_norm = sc.tl.marker_gene_overlap(adata, marker_genes, key='rank_genes_r0.75', normalize='reference')
sb.heatmap(cell_annotation_norm, cbar=False, annot=True)
#Define a nice colour map for gene expression
colors2 = plt.cm.Reds(np.linspace(0, 1, 128))
colors3 = plt.cm.Greys_r(np.linspace(0.7,0.8,20))
colorsComb = np.vstack([colors3, colors2])
mymap = colors.LinearSegmentedColormap.from_list('my_colormap', colorsComb)
#Defa24 #Tff3
sc.pl.umap(adata, color='Defa24', use_raw=False, color_map=mymap)
sc.pl.umap(adata, color='Tff3', use_raw=False, color_map=mymap)
# Check expression of enterocyte markers
#Collate all enterocyte markers and get the gene IDs in the data set
ids_entprox = np.in1d(adata.var_names, marker_genes['Enterocyte (Proximal)'])
ids_entdist = np.in1d(adata.var_names, marker_genes['Enterocyte (Distal)'])
ids_ent = np.logical_or(ids_entprox, ids_entdist)

#Calculate the mean expression of enterocyte markers
adata.obs['Enterocyte_marker_expr'] = adata.X[:,ids_ent].mean(1)

#Plot enterocyte expression
sc.pl.violin(adata, 'Enterocyte_marker_expr', groupby='louvain_r0.75')
sc.pl.umap(adata, color='Enterocyte_marker_expr', color_map=mymap)

#Early enterocyte marker - Arg2
sc.pl.umap(adata, color='Arg2', use_raw=False, color_map=mymap)

sc.pl.violin(adata, groupby='louvain_r0.75', keys='Arg2', use_raw=False)

sc.pl.diffmap(adata, components=['6,9'], color='Arg2', use_raw=False, color_map=mymap)
sc.pl.diffmap(adata, components=['6,9'], color='louvain_r0.75')
sc.pl.violin(adata, 'mt_frac', groupby='louvain_r0.75')
sc.pl.violin(adata, 'log_counts', groupby='louvain_r0.75')
#Check individual stem markers
stem_genes = adata.var_names[np.in1d(adata.var_names, marker_genes['Stem'])]
sc.pl.umap(adata, color=stem_genes[:3], title=stem_genes[:3], color_map=mymap)
sc.pl.umap(adata, color=stem_genes[3:], title=stem_genes[3:], color_map=mymap)

#Check stem marker expression
adata.obs['Stem_marker_expr'] = adata[:,stem_genes].X.mean(1)

sc.pl.violin(adata, 'Stem_marker_expr', groupby='louvain_r0.75')
sc.pl.umap(adata, color='Stem_marker_expr', color_map=mymap)
#Categories to rename
adata.obs['louvain_r0.75'].cat.categories
adata.rename_categories('louvain_r0.75', ['Stem', 'TA', 'Goblet', 'EP (stress)', 'Enterocyte', 'EP (early)', 'Paneth', 'Enteroendocrine', 'Tuft']) #adata.rename_categories('louvain_r0.5', ['TA', 'EP (early)', 'Stem', 'Goblet', 'EP (stress)', 'Enterocyte', 'Paneth', 'Enteroendocrine', 'Tuft'])
adata.obs['louvain_r0.75'].value_counts() #new categories need to have the same number of items as the old categories!
sc.pl.umap(adata, color='louvain_r0.75', size=15, legend_loc='on data')
#%%
#Subcluster enterocytes
sc.tl.louvain(adata, restrict_to=('louvain_r0.75', ['Enterocyte']), resolution=0.2, key_added='louvain_r0.75_entero_sub')
#Show the new clustering
if 'louvain_r0.75_entero_sub_colors' in adata.uns:
    del adata.uns['louvain_r0.75_entero_sub_colors']

sc.pl.umap(adata, color='louvain_r0.75_entero_sub', palette=sc.pl.palettes.godsnot_102)
sc.pl.umap(adata, color='region', palette=sc.pl.palettes.godsnot_102)

#Get the new marker genes
sc.tl.rank_genes_groups(adata, groupby='louvain_r0.75_entero_sub', key_added='rank_genes_r0.75_entero_sub')


#Plot the new marker genes
sc.pl.rank_genes_groups(adata, key='rank_genes_r0.75_entero_sub', groups=['Enterocyte,0','Enterocyte,1','Enterocyte,2'], fontsize=12)

entero_clusts = [clust for clust in adata.obs['louvain_r0.75_entero_sub'].cat.categories if clust.startswith('Enterocyte')]
from scipy import sparse
adata.X = sparse.csc_matrix(adata.X)
for clust in entero_clusts:
    sc.pl.rank_genes_groups_violin(adata, use_raw=False, key='rank_genes_r0.75_entero_sub', groups=[clust], gene_names=adata.uns['rank_genes_r0.75_entero_sub']['names'][clust][90:100])

#%%
#Subset marker gene dictionary to only check for enterocyte markers
marker_genes_entero = {k: marker_genes[k] for k in marker_genes if k.startswith('Enterocyte')}

#Find marker overlap
sc.tl.marker_gene_overlap(adata, marker_genes_entero, key='rank_genes_r0.75_entero_sub', normalize='reference')

#Check enterocyte marker expression
sc.pl.violin(adata[adata.obs['louvain_r0.75']=='Enterocyte'], groupby='louvain_r0.75_entero_sub', keys='Enterocyte_marker_expr')
#Visualize some enterocyte markers
entero_genes = ['Alpi', 'Apoa1', 'Apoa4', 'Fabp1', 'Arg2']
sc.pl.umap(adata, color=entero_genes[:3], title=entero_genes[:3], color_map=mymap)
sc.pl.umap(adata, color=entero_genes[3:], title=entero_genes[3:], color_map=mymap)

sc.pl.diffmap(adata, color='louvain_r0.75_entero_sub', components='3,7')
tmp = adata.obs['louvain_r0.75_entero_sub'].cat.categories

tmp = ['Enterocyte imm. (Distal)' if item == 'Enterocyte,0' else item for item in tmp]
tmp = ['Enterocyte imm. (Proximal)' if item == 'Enterocyte,1' else item for item in tmp]
tmp = ['Enterocyte mature' if item == 'Enterocyte,2' else item for item in tmp]

adata.rename_categories('louvain_r0.75_entero_sub', tmp)

#Subcluster mature enterocytes
sc.tl.louvain(adata, restrict_to=('louvain_r0.75_entero_sub', ['Enterocyte mature']), resolution=0.25, key_added='louvain_r0.75_entero_mat_sub')
#Show the new clustering
if 'louvain_r0.75_entero_mat_sub_colors' in adata.uns:
    del adata.uns['louvain_r0.75_entero_mat_sub_colors']

sc.pl.umap(adata, color='louvain_r0.75_entero_mat_sub', palette=sc.pl.palettes.godsnot_102)
#Get the new marker genes
sc.tl.rank_genes_groups(adata, groupby='louvain_r0.75_entero_mat_sub', key_added='rank_genes_r0.75_entero_mat_sub')
#Plot the new marker genes
sc.pl.rank_genes_groups(adata, key='rank_genes_r0.75_entero_mat_sub', groups=['Enterocyte mature,0','Enterocyte mature,1'], fontsize=12)
entero_mat_clusts = [clust for clust in adata.obs['louvain_r0.75_entero_mat_sub'].cat.categories if clust.startswith('Enterocyte mature')]

for clust in entero_mat_clusts:
    sc.pl.rank_genes_groups_violin(adata, use_raw=False, key='rank_genes_r0.75_entero_mat_sub', groups=[clust], gene_names=adata.uns['rank_genes_r0.75_entero_mat_sub']['names'][clust][90:100])
#%%
#Find marker overlap
sc.tl.marker_gene_overlap(adata, marker_genes_entero, key='rank_genes_r0.75_entero_mat_sub', normalize='reference')
tmp = adata.obs['louvain_r0.75_entero_mat_sub'].cat.categories

tmp = ['Enterocyte mat. (Distal)' if item == 'Enterocyte mature,0' else item for item in tmp]
tmp = ['Enterocyte mat. (Proximal)' if item == 'Enterocyte mature,1' else item for item in tmp]

adata.rename_categories('louvain_r0.75_entero_mat_sub', tmp)
adata.obs['louvain_final'] = adata.obs['louvain_r0.75_entero_mat_sub']
sc.pl.umap(adata, color='louvain_final', palette=sc.pl.palettes.vega_20, legend_loc='on data')
adata.obs['louvain_final'].value_counts()
#Define a variable that stores proximal and distal labels
adata.obs['prox_dist'] = ['Distal' if reg=='Il' else 'Proximal' for reg in adata.obs['region']]
sc.tl.embedding_density(adata, basis='umap', groupby='prox_dist')
adata.obs['prox_dist'].value_counts()
sc.pl.embedding_density(adata, basis='umap', key='umap_density_prox_dist', group='Proximal')
sc.pl.embedding_density(adata, basis='umap', key='umap_density_prox_dist', group='Distal')
adata.obs['louvain_final'].value_counts()
#Subsetting to relevant clusters
clusters_to_include = [g for g in adata.obs['louvain_final'].cat.categories if (g.startswith('Enterocyte') or g.startswith('TA') or g.startswith('Stem') or g.startswith('EP'))]
adata_ent = adata[np.isin(adata.obs['louvain_final'], clusters_to_include),:].copy()

#Subset to highly variable genes
sc.pp.highly_variable_genes(adata_ent, flavor='cell_ranger', n_top_genes=4000, subset=True)

#Recalculating PCA for subset
sc.pp.pca(adata_ent, svd_solver='arpack')
sc.pl.pca(adata_ent)
sc.pl.pca_variance_ratio(adata_ent)
adata_ent.obsm['X_pca'] = adata_ent.obsm['X_pca'][:,0:7]
#%%
#R -i adata_ent

import numpy as np
from scipy import sparse

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, numpy2ri
from rpy2.robjects.conversion import localconverter

ro.r("library(Matrix)")

def dgc_to_csr(r_dgc):
    """Convert (and transpose) a dgCMatrix from R to a csr_matrix in python
    """
    with localconverter(ro.default_converter + pandas2ri.converter):
        X = sparse.csr_matrix(
                (
                    r_dgc.slots["x"], 
                    r_dgc.slots["i"], 
                    r_dgc.slots["p"]
                ),
                shape=tuple(ro.r("dim")(r_dgc))[::-1]
            )
    return X

def csr_to_dgc(csr):
    """Convert (and transpose) a csr matrix from python to a R dgCMatrix (not sure if type is consistent)
    """
    print(csr.shape)
    numeric = ro.r("as.numeric")
    with localconverter(ro.default_converter + ro.numpy2ri.converter):
        X = ro.r("sparseMatrix")(
            i=numeric(csr.indices),
            p=numeric(csr.indptr),
            x=numeric(csr.data),
            index1=False
        )
    return X

for i in range(10):
    X = sparse.rand(1000, 100, density=.1, format="csr")
    assert np.allclose(dgc_to_csr(csr_to_dgc(X)).todense(), X.todense())

#%%
def csr_to_dgc(adata_ent):
    """Convert (and transpose) a csr matrix from python to a R dgCMatrix (not sure if type is consistent)
    """
    print(adata_ent.shape)
    numeric = robjects.r("as.numeric")
    with localconverter(rojects.default_converter + robjects.numpy2ri.converter):
        X = robjects.r("sparseMatrix")(
            i=numeric(adata_ent.indices),
            p=numeric(adata_ent.indptr),
            x=numeric(adata_ent.data),
            index1=False
        )
    return X
from rpy2.robjects.conversion import localconverter
print(adata_ent.shape)
numeric = robjects.r("as.numeric")
with localconverter(robjects.default_converter + robjects.numpy2ri.converter):
    X = robjects.r("sparseMatrix")(
    i=numeric(adata_ent.indices),
    p=numeric(adata_ent.indptr),
    x=numeric(adata_ent.data),
    index1=False)
return X


#%%
adata_ent.obs.dtypes
adata_ent.var.dtypes

robjects.globalenv['adata_ent'] = anndata2ri.py2rpy(adata_ent)




type(adata_ent)
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
#This is what you will need to make non interactive plots
grdevices = importr('grDevices')

#Plot 1
colour_map = robjects.r(f'brewer.pal(20,"Set1")')

grdevices.png(file="/home/jake/Documents/SC_RNA_SEQ/file_1.png", width=512, height=512)
robjects.r('''par(xpd=TRUE)
par(mar=c(4.5,5.5,2,7)
plot(reducedDims(adata_ent)$PCA[,1], reducedDims(adata_ent)$PCA[,2], col=colour_map[colData(adata_ent)$louvain_final], bty="L", xlab="PC1", ylab="PC2")
legend(x=12, y=12, legend=unique(colData(adata_ent)$louvain_final), fill=colour_map[as.integer(unique(colData(adata_ent)$louvain_final))])
''')
grdevices.dev_off()

robjects.r.ls(globalenv)
robjects.globalenv["a"] = 123
print(robjects.r.ls())

print("1:")
adata_ent_start = robjects.r(f'slingshot(adata_ent, clusterLabels = "louvain_final", reducedDim = "PCA", start.clus="Stem")')
print(SlingshotDataSet(adata_ent_start))

print("")
print("2:")
adata_ent_startend <- slingshot(adata_ent, clusterLabels = 'louvain_final', reducedDim = 'PCA', start.clus='Stem', end.clus=c('Enterocyte mat. (Proximal)', 'Enterocyte mat. (Distal)'))
print(SlingshotDataSet(adata_ent_startend))

print("")
print("3:")
adata_ent_simple_startend <- slingshot(adata_ent, clusterLabels = 'louvain_r0.5', reducedDim = 'PCA', start.clus='Stem', end.clus='Enterocyte')
print(SlingshotDataSet(adata_ent_simple_startend))

%%R

#Plot of lineage 1
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plot(reducedDims(adata_ent_startend)$PCA[,c(1,2)], col = colors[cut(adata_ent_startend$slingPseudotime_1,breaks=100)], pch=16, asp = 1, xlab='PC1', ylab='PC2')
lines(slingCurves(adata_ent_startend)$curve1, lwd=2)

#Plot of lineage 2
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plot(reducedDims(adata_ent_startend)$PCA[,c(1,2)], col = colors[cut(adata_ent_startend$slingPseudotime_2,breaks=100)], pch=16, asp = 1, xlab='PC1', ylab='PC2')
lines(slingCurves(adata_ent_startend)$curve2, lwd=2)

#Plot of lineages with clusters visualized
par(xpd=TRUE)
plot(reducedDims(adata_ent_startend)$PCA[,c(1,2)], col = brewer.pal(11,'Set1')[adata_ent$louvain_final], pch=16, asp = 1, bty='L', xlab='PC1', ylab='PC2')
lines(SlingshotDataSet(adata_ent_startend), lwd=2, type='lineages')
legend(x=10, y=20, legend=unique(colData(adata_ent)$louvain_final), fill=brewer.pal(11,'Set1')[as.integer(unique(colData(adata_ent)$louvain_final))])

#Plot of simpler clustering
plot(reducedDims(adata_ent_simple_startend)$PCA[,c(1,2)], col = colors[cut(adata_ent_simple_startend$slingPseudotime_1,breaks=100)], pch=16, asp = 1, xlab='PC1', ylab='PC2')
lines(SlingshotDataSet(adata_ent_simple_startend), lwd=2)
#Subsetting data set - non-batch corrected
cell_mask = np.isin(adata.obs['louvain_final'], clusters_to_include)
adata_ent_nbc = sc.AnnData(X = adata.raw.X[cell_mask,:])
adata_ent_nbc.obs = adata.obs[cell_mask]
adata_ent_nbc.var = adata.var.copy()

#Subset to highly variable genes
sc.pp.highly_variable_genes(adata_ent_nbc, flavor='cell_ranger', n_top_genes=4000, subset=True)
#Recalculating PCA for subset
sc.pp.pca(adata_ent_nbc, svd_solver='arpack')
sc.pl.pca(adata_ent_nbc)
sc.pl.pca_variance_ratio(adata_ent_nbc)
adata_ent_nbc.obsm['X_pca'] = adata_ent_nbc.obsm['X_pca'][:,0:7]
%%R -i adata_ent_nbc

#Plot 1
colour_map = brewer.pal(20,'Set1')
par(xpd=TRUE)
par(mar=c(4.5,5.5,2,11))
plot(reducedDims(adata_ent_nbc)$PCA[,1], reducedDims(adata_ent_nbc)$PCA[,2], col=colour_map[colData(adata_ent_nbc)$louvain_final], bty='L', xlab='PC1', ylab='PC2')
legend(x=24.5, y=0, legend=unique(colData(adata_ent_nbc)$louvain_final), fill=colour_map[as.integer(unique(colData(adata_ent_nbc)$louvain_final))])

#First trajectory: only Stem cells set as root cells
print("1:")
adata_ent_start_nbc <- slingshot(adata_ent_nbc, clusterLabels = 'louvain_final', reducedDim = 'PCA', start.clus='Stem')
print(SlingshotDataSet(adata_ent_start_nbc))

#Second trajectory: Stem cells as root cells and mature enterocytes as end clusters
print("")
print("2:")
adata_ent_startend_nbc <- slingshot(adata_ent_nbc, clusterLabels = 'louvain_final', reducedDim = 'PCA', start.clus='Stem', end.clus=c('Enterocyte mat. (Proximal)', 'Enterocyte mat. (Distal)'))
print(SlingshotDataSet(adata_ent_startend_nbc))

#Third trajectory: Stem cells as root cells and enterocytes as end cluster, non-subclustered data
print("")
print("3:")
adata_ent_simple_startend_nbc <- slingshot(adata_ent_nbc, clusterLabels = 'louvain_r0.5', reducedDim = 'PCA', start.clus='Stem', end.clus='Enterocyte')
print(SlingshotDataSet(adata_ent_simple_startend_nbc))
%%R

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)

#Plot of lineage 1
plot(reducedDims(adata_ent_startend_nbc)$PCA[,c(1,2)], 
     col = colors[cut(adata_ent_startend_nbc$slingPseudotime_1,breaks=100)], 
     pch=16, asp = 1, xlab='PC1', ylab='PC2')
lines(slingCurves(adata_ent_startend_nbc)$curve1, lwd=2)

#Plot of lineage 2
plot(reducedDims(adata_ent_startend_nbc)$PCA[,c(1,2)], 
     col = colors[cut(adata_ent_startend_nbc$slingPseudotime_2,breaks=100)], 
     pch=16, asp = 1, xlab='PC1', ylab='PC2')
lines(slingCurves(adata_ent_startend_nbc)$curve2, lwd=2)

#Combined plot
plot(reducedDims(adata_ent_startend_nbc)$PCA[,c(1,2)], col = 'grey', 
     pch=16, asp = 1, size=0.8, xlab='PC1', ylab='PC2')
points(reducedDims(adata_ent_startend_nbc)$PCA[,c(1,2)], 
       col = colors[cut(adata_ent_startend_nbc$slingPseudotime_1,breaks=100)], 
       pch=16, size=1)
lines(slingCurves(adata_ent_startend_nbc)$curve1, lwd=2)
lines(slingCurves(adata_ent_startend_nbc)$curve2, lwd=2)

#Plot of lineages with clusters visualized
par(xpd=TRUE)
par(mar=c(4.5,5.5,4,1))
plot(reducedDims(adata_ent_startend_nbc)$PCA[,c(1,2)], 
     col = brewer.pal(11,'Set1')[adata_ent_startend_nbc$louvain_final], 
     pch=16, asp = 1, bty='L', xlab='PC1', ylab='PC2')
lines(SlingshotDataSet(adata_ent_startend_nbc), lwd=2, type='lineages')
legend(x=6, y=20, legend=unique(colData(adata_ent_nbc)$louvain_final), fill=brewer.pal(11,'Set1')[as.integer(unique(colData(adata_ent_nbc)$louvain_final))])

#Plot of simpler clustering
plot(reducedDims(adata_ent_simple_startend_nbc)$PCA[,c(1,2)], 
     col = colors[cut(adata_ent_simple_startend_nbc$slingPseudotime_1,breaks=100)], 
     pch=16, asp = 1, xlab='PC1', ylab='PC2')
lines(SlingshotDataSet(adata_ent_simple_startend_nbc), lwd=2)
#Preprocessing for monocle
data_mat_mon = adata.layers['counts'].T
var_mon=adata.var.copy()
obs_mon=adata.obs.copy()
%%R -i data_mat_mon -i obs_mon -i var_mon

#Set up the CellDataSet data structure
pd <- AnnotatedDataFrame(data = obs_mon)
fd <- AnnotatedDataFrame(data = var_mon)
colnames(data_mat_mon) <- rownames(pd)
rownames(data_mat_mon) <- rownames(fd)
ie_regions_cds <- newCellDataSet(cellData=data_mat_mon, phenoData=pd, featureData=fd, expressionFamily=negbinomial.size())

#Normalize the count data
ie_regions_cds <- estimateSizeFactors(ie_regions_cds)

#Calculate dispersions to filter for highly variable genes
ie_regions_cds <- estimateDispersions(ie_regions_cds)


#Filter for Stem, EP, TA, and Enterocytes
cell_types = as.character(pData(ie_regions_cds)$louvain_final)
cell_mask = rep(FALSE, length(cell_types))
cells_to_keep = c("Stem", "EP", "TA", "Enterocyte")
for (item in cells_to_keep) {cell_mask = cell_mask | startsWith(cell_types, item)}
print("Number of cells after filtering:")
print(sum(cell_mask))
print("")

#Filter highly variable genes from our analysis
hvg_mask = fData(ie_regions_cds)$highly_variable
ie_regions_cds <- ie_regions_cds[hvg_mask, cell_mask]

#Do dimensionality reduction
ie_regions_cds <- reduceDimension(ie_regions_cds, norm_method = 'vstExprs', reduction_method='DDRTree', verbose = F, max_components = 7)

#Run for the first time to get the ordering
ie_regions_cds <- orderCells(ie_regions_cds)

#Find the correct root state the corresponds to the 'Stem' cluster
tab1 <- table(pData(ie_regions_cds)$State, pData(ie_regions_cds)$louvain_final)
id = which(colnames(tab1) == 'Stem')
root_name = names(which.max(tab1[,id]))

#Run a second time to get the correct root state that overlaps with Stem cells
ie_regions_cds <- orderCells(ie_regions_cds, root_state=root_name)
%%R -w 1000 -h 800

#Get a nice colour map
custom_colour_map = brewer.pal(length(unique(pData(ie_regions_cds)$louvain_final)),'Paired')

#Find the correct root state that coresponds to the 'Stem' cluster
tab1 <- table(pData(ie_regions_cds)$State, pData(ie_regions_cds)$louvain_final)
id = which(colnames(tab1) == 'Stem')
root_name = names(which.max(tab1[,id]))

# Visualize with our cluster labels
options(repr.plot.width=5, repr.plot.height=4)
plot_complex_cell_trajectory(ie_regions_cds, color_by = 'louvain_final', show_branch_points = T, 
                             cell_size = 2, cell_link_size = 1, root_states = c(root_name)) +
scale_size(range = c(0.2, 0.2)) +
theme(legend.position="left", legend.title=element_blank(), legend.text=element_text(size=rel(1.5))) +
guides(colour = guide_legend(override.aes = list(size=6))) + 
scale_color_manual(values = custom_colour_map)
%%R -w 600 -h 800

#Visualize pseudotime found
options(repr.plot.width=5, repr.plot.height=4)
plot_cell_trajectory(ie_regions_cds, color_by="Pseudotime")
sc.pp.neighbors(adata_ent)
sc.tl.diffmap(adata_ent)
sc.pl.diffmap(adata_ent, components='1,2', color='louvain_final')
sc.pl.diffmap(adata_ent, components='1,3', color='louvain_final')

#Find the stem cell with the highest DC3 value to act as root for the diffusion pseudotime and compute DPT
stem_mask = np.isin(adata_ent.obs['louvain_final'], 'Stem')
max_stem_id = np.argmin(adata_ent.obsm['X_diffmap'][stem_mask,3])
root_id = np.arange(len(stem_mask))[stem_mask][max_stem_id]
adata_ent.uns['iroot'] = root_id

#Compute dpt
sc.tl.dpt(adata_ent)
#Visualize pseudotime over differentiation
sc.pl.diffmap(adata_ent, components='1,3', color='dpt_pseudotime')
%%R

#Set the pseudotime variable
t <- adata_ent_simple_startend$slingPseudotime_1

#Extract the gene expression matrix
Y <- assay(adata_ent_simple_startend)

# fit a GAM with a loess term for pseudotime
#Note: This takes a while
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
})
%%R -w 600 -h 1200

#Select the top 100 most significant genes that change over pseudotime
topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]
heatdata <- assay(adata_ent_simple_startend)[rownames(assay(adata_ent_simple_startend)) %in% topgenes, 
                        order(t, na.last = NA)]

#Scale the data per gene for visualization
heatdata <- t(scale(t(heatdata)))

#Trimm z-score scale
heatdata[heatdata > 3] = 3
heatdata[heatdata < -3] = -3

#Get cluster assignment
heatclus <- adata_ent_simple_startend$louvain_r0.5[order(t, na.last = NA)]

#Set up a clusterExperiment data structure for heatmap visualization
ce <- ClusterExperiment(heatdata, heatclus, transformation = function(x){x})

#Plot the heatmap
plotHeatmap(ce, clusterSamplesData = "orderSamplesValue", visualizeData = 'transformed', fontsize=15)
entero_markers = marker_genes['Enterocyte (Proximal)'] + marker_genes['Enterocyte (Distal)']

%%R -i entero_markers
print(rownames(heatdata)[rownames(heatdata) %in% entero_markers])

%%R

pt1 <- adata_ent_startend_nbc$slingPseudotime_1
clustDat <- adata_ent_startend_nbc$louvain_final

#Subset data to only include cells on lineage 1
clustDat <- clustDat[!is.na(pt1)]
pt1 <- pt1[!is.na(pt1)]
df = data.frame(clusters = clustDat, pt = pt1)

#Bin clusters in same way as pseudotime:
bin_width = 0.5
max_bin = ceiling(max(df$pt)*2)/2
df['bins'] = cut(df$pt, breaks=seq(-bin_width/2, max_bin, bin_width))

#Find dominant cluster in each bin
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

dominant_clusters = sapply(levels(df$bins), function(x){Mode(df$clust[df$bins==x])})
levels(dominant_clusters) <- c(levels(dominant_clusters), 'None')
dominant_clusters[is.na(dominant_clusters)] <- 'None'

#Define colour map
cmap <- brewer.pal(11,'Set1')

#Plot meta-stable states
p <- qplot(adata_ent_startend_nbc$slingPseudotime_1, geom='histogram', main='Cellular density over trajectory', xlab='Pseudotime', binwidth=bin_width, fill=I(cmap[dominant_clusters])) +
theme_bw() +
scale_colour_manual(name="Cluster", values=c('EP (Distal)'=cmap[levels(dominant_clusters)[1]], 'EP (early)'=cmap[levels(dominant_clusters)[2]], EP_early="purple"))
print(p)

#Plot of lineages with clusters visualized
par(xpd=TRUE)
par(mar=c(4.5,5.5,2,1))
plot(reducedDims(adata_ent_startend_nbc)$PCA[,c(1,2)], col = brewer.pal(11,'Set1')[adata_ent$louvain_final], pch=16, asp = 1, bty='L', xlab='PC1', ylab='PC2')
lines(SlingshotDataSet(adata_ent_startend_nbc), lwd=2, type='lineages')
legend(x=10, y=20, legend=unique(colData(adata_ent)$louvain_final), fill=cmap[as.integer(unique(colData(adata_ent)$louvain_final))])
sc.tl.paga(adata, groups='louvain_final')
sc.pl.paga_compare(adata)
sc.pl.paga(adata)
sc.tl.paga(adata, groups='louvain_r0.5')
sc.pl.paga_compare(adata)
sc.pl.paga(adata)
sc.pl.paga_compare(adata, basis='umap')
fig1, ax1 = plt.subplots()
sc.pl.umap(adata, size=40, ax=ax1, show=False)
sc.pl.paga(adata, pos=adata.uns['paga']['pos'], show=False, node_size_scale=10, node_size_power=1, ax=ax1, text_kwds={'alpha':0})
#plt.savefig('./figures/umap_paga_overlay_gut.pdf', dpi=300, format='pdf')
plt.show()
#Create new Anndata object with non-batch corrected data
adata_test = adata.copy()
adata_test.X = adata.raw.X
#Regress out counts and redo pre-processing
sc.pp.regress_out(adata_test, 'n_counts')
sc.pp.pca(adata_test, use_highly_variable=True, svd_solver='arpack')
sc.pp.neighbors(adata_test)
#Recalculate PAGA
sc.tl.paga(adata_test, groups='louvain_r0.5')
sc.pl.paga_compare(adata_test, basis='umap')
sc.pl.paga(adata_test)

#Create new Anndata object for use in MAST with non-batch corrected data as before
adata_test = adata.copy()
adata_test.X = adata.raw.X
adata_test.obs['n_genes'] = (adata_test.X > 0).sum(1) # recompute number of genes expressed per cell
%%R -i adata_test -o ent_de -o paneth_de

#Convert SingleCellExperiment to SingleCellAssay type as required by MAST
sca <- SceToSingleCellAssay(adata_test, class = "SingleCellAssay")

#Scale Gene detection rate
colData(sca)$n_genes = scale(colData(sca)$n_genes)

#Create data subsets for paneth and Enterocyte subpopulations
sca_ent <- subset(sca, with(colData(sca), louvain_r0.5=='Enterocyte'))
sca_paneth <- subset(sca, with(colData(sca), louvain_r0.5=='Paneth'))


#Filter out non-expressed genes in the subsets
print("Dimensions before subsetting:")
print(dim(sca_ent))
print(dim(sca_paneth))
print("")

sca_ent_filt = sca_ent[rowSums(assay(sca_ent)) != 0, ]
sca_paneth_filt = sca_paneth[rowSums(assay(sca_paneth)) != 0, ]

print("Dimensions after subsetting:")
print(dim(sca_ent_filt))
print(dim(sca_paneth_filt))


#Define & run hurdle model - Enterocytes
zlmCond_ent <- zlm(formula = ~prox_dist + donor + n_genes, sca=sca_ent_filt)
summaryCond_ent <- summary(zlmCond_ent, doLRT='prox_distProximal')
summaryDt_ent <- summaryCond_ent$datatable

result_ent <- merge(summaryDt_ent[contrast=='prox_distProximal' & component=='H',.(primerid, `Pr(>Chisq)`)], #P-vals
                 summaryDt_ent[contrast=='prox_distProximal' & component=='logFC', .(primerid, coef)],
                 by='primerid') #logFC coefficients

#Correct for multiple testing (FDR correction) and filtering
result_ent[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')]
ent_de = result_ent[result_ent$FDR<0.01,, drop=F]
ent_de = ent_de[order(ent_de$FDR),]


#Define & run hurdle model - paneth cells
zlmCond_paneth <- zlm(formula = ~prox_dist + donor + n_genes, sca=sca_paneth_filt)
summaryCond_paneth <- summary(zlmCond_paneth, doLRT='prox_distProximal')
summaryDt_paneth <- summaryCond_paneth$datatable

result_paneth <- merge(summaryDt_paneth[contrast=='prox_distProximal' & component=='H',.(primerid, `Pr(>Chisq)`)], #P-vals
                 summaryDt_paneth[contrast=='prox_distProximal' & component=='logFC', .(primerid, coef)],
                 by='primerid') #logFC coefficients

#Correct for multiple testing (FDR correction) and filtering
result_paneth[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')]
paneth_de = result_paneth[result_paneth$FDR<0.01,, drop=F]
paneth_de = paneth_de[order(paneth_de$FDR),]

#Show top 20 differentially expressed genes for enterocytes (up- and down-regulated)
print(ent_de.shape)
ent_de[:20]

#Volcano plot of results
ent_de['-logQ'] = -np.log(ent_de['FDR'])
lowqval_de = ent_de.loc[ent_de['-logQ'] > 200]
other_de = ent_de.loc[ent_de['-logQ'] < 200]

fig, ax = plt.subplots()
sb.regplot(other_de['coef'], other_de['-logQ'], fit_reg=False, scatter_kws={'s':6})
sb.regplot(lowqval_de['coef'], lowqval_de['-logQ'], fit_reg=False, scatter_kws={'s':6})
ax.set_xlabel("log2 FC", fontsize=20)
ax.set_ylabel("-log Q-value", fontsize=20)
ax.tick_params(labelsize=15)

# Label names and positions
x = [i-0.2 for i in lowqval_de['coef']]
y = [i+10 for i in lowqval_de['-logQ']]
labels = lowqval_de['primerid']

# Show only some labels to avoid overcrowding the figure
to_remove = np.where([i < 230 for i in y])[0]
labels = ["" if i in to_remove else lab for i,lab in enumerate(labels) ]

#Move up two labels
y = [y[i]+10 if txt == 'Krt8' else y[i] for i,txt in enumerate(labels)]
y = [y[i]+20 if txt == 'Cd9' else y[i] for i,txt in enumerate(labels)]

#Move down one label
y = [y[i]-20 if txt == 'Phgr1' else y[i] for i,txt in enumerate(labels)]

for i,txt in enumerate(labels):
    ax.annotate(txt, (x[i], y[i]))
plt.show()

#See overlaps with markers expected for proximal and distal enterocytes
prox_de_set = set(ent_de['primerid'][ent_de['coef'] > 0])
dist_de_set = set(ent_de['primerid'][ent_de['coef'] < 0])
print("Fraction of proximal enterocyte markers in up-regulated proximal DE genes: {}".format(len(prox_de_set.intersection(marker_genes['Enterocyte (Proximal)']))/len(marker_genes['Enterocyte (Proximal)'])))
print("Fraction of distal enterocyte markers in up-regulated proximal DE genes: {}".format(len(prox_de_set.intersection(marker_genes['Enterocyte (Distal)']))/len(marker_genes['Enterocyte (Distal)'])))

print()
print("Fraction of proximal enterocyte markers in up-regulated distal DE genes: {}".format(len(dist_de_set.intersection(marker_genes['Enterocyte (Proximal)']))/len(marker_genes['Enterocyte (Proximal)'])))
print("Fraction of distal enterocyte markers in up-regulated distal DE genes: {}".format(len(dist_de_set.intersection(marker_genes['Enterocyte (Distal)']))/len(marker_genes['Enterocyte (Distal)'])))

#Interpretation of differentially expressed genes in paneth cells - g:profiler
gp = GProfiler(return_dataframe=True, user_agent='g:GOSt')

paneth_enrichment = gp.profile(organism='mmusculus', sources=['GO:BP'], user_threshold=0.05,
                               significance_threshold_method='fdr', 
                               background=adata.var_names.tolist(), 
                               query=paneth_de['primerid'].tolist())

#paneth_enrich_results = paneth_enrichment.sort_values('p.value').iloc[:,[2,3,5,6,11]]
paneth_enrich_results = paneth_enrichment.set_index('native').sort_values('p_value').iloc[:,[2,5,7,10,1]]

pd.set_option("display.max_colwidth", 800)
paneth_enrich_results.iloc[:50,:]


from gprofiler_plotting import plot_enrich

plot_enrich(paneth_enrich_results)

#Write to file
adata.write('../data/Haber-et-al_mouse-intestinal-epithelium/Haber_et_al_case_study.h5ad')