import os
import scipy
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import matplotlib.pyplot as plt

scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params('scvelo')
pd.set_option('display.max_columns', 100)
pd.set_option('display.max_rows', 200)
np.set_printoptions(suppress=True)

##=================
## parameters - users need to check and edit
##=================
## Output Directory
## user needs to edit this path
OutDir = "/path/to/Out_Multivelo"

## .loom file - output from Velocyto
## user needs to edit this filename
LoomFile = "/path/to/velocyto_Out/cellsorted_gex_possorted_bam.bam"

## cell annotation file
## with respect to clustering in multimodal data
## output of the script 1_Extract_Seurat_WNN_Out.R
## user needs to edit this filename
AnnotFile = "/path/to/Multimodal_Cluster_Annotation.tsv"

## folder storing the ATAC-seq counts
## user needs to edit
BaseCountDir = "/path/to/Cellranger/count_ARC/outs/"
CountDir = BaseCountDir + "filtered_feature_bc_matrix/"

## peak annotation file from Cellranger output
PeakAnnotFile = BaseCountDir + "atac_peak_annotation.tsv"

## feature linkage file
FeatureLinkageFile = BaseCountDir + "analysis/feature_linkage/feature_linkage.bedpe"

## output of the script 1_Extract_Seurat_WNN_Out.R
## user needs to edit
nn_idx_file = "/path/to/nn_idx.txt"
nn_dist_file = "/path/to/nn_dist.txt"
nn_cells_file = "/path/to/nn_cells.txt"

COUNT_MIN_THR = 1000
COUNT_MAX_THR = 20000
MIN_SHARED_COUNT_THR = 10
N_TOP_GENES_THR = 1000

ATAC_COUNT_MIN_THR = 2000
ATAC_COUNT_MAX_THR = 60000

NUM_PCS_THR = 30
NUM_NEIGHBORS_THR = 50

## target genes whose characteristics is to be plotted
## user needs to edit this file path
TargetGenelistFile = "/path/to/target_genes.txt"

##=================
## Reading in unspliced and spliced counts
## read the RNA part from cellranger output
## filter and select the cells included in the input cell annotation file
##=================
adata_rna = scv.read(LoomFile, cache=True)
adata_rna.obs_names = [x.split(':')[1][:-1] + '-1' for x in adata_rna.obs_names]
## Variable names are not unique. To make them unique, call `.var_names_make_unique`.
adata_rna.var_names_make_unique()
print("\n\n Initial RNA part read - number of rows : " + str(adata_rna.shape[0]))

# Load cell annotations
cell_annot = pd.read_csv(AnnotFile, sep='\t', index_col=0)
print("\n cell_annot - number of rows : " + str(cell_annot.shape[0]))

## Paramita - first do this assignment, and then filter
## the "celltype" field in the annotation file store the cluster information
adata_rna = adata_rna[cell_annot.index,:]
adata_rna.obs['celltype'] = cell_annot['celltype']
print("\n After filtering by cell annotation - adata_rna - number of rows : " + str(adata_rna.shape[0]))

## Paramita - no further filtering is required
if 0:
	sc.pp.filter_cells(adata_rna, min_counts=COUNT_MIN_THR)
	sc.pp.filter_cells(adata_rna, max_counts=COUNT_MAX_THR)

# Top 1000 variable genes are used for downstream analyses.
scv.pp.filter_and_normalize(adata_rna, 
							min_shared_counts=MIN_SHARED_COUNT_THR, 
							n_top_genes=N_TOP_GENES_THR)
print("\n After pp.filter_and_normalize - adata_rna - number of rows : " + str(adata_rna.shape[0]))

# We subset for specific clusters
## Paramita - commented
if 0:
	adata_rna = adata_rna[adata_rna.obs['celltype'].isin(['RG, Astro, OPC', 'IPC', 'V-SVZ', 'Upper Layer', 'Deeper Layer', 'Ependymal cells', 'Subplate'])]
	adata_rna

##=================
## Preprocessing the ATAC counts
##=================
adata_atac = sc.read_10x_mtx(CountDir, var_names='gene_symbols', cache=True, gex_only=False)
adata_atac = adata_atac[:,adata_atac.var['feature_types'] == "Peaks"]

# We aggregate peaks around each gene as well 
# as those that have high correlations with promoter peak or gene expression.
# Peak annotation contains the metadata for all peaks.
# Feature linkage contains pairs of correlated genomic features.
adata_atac = mv.aggregate_peaks_10x(adata_atac, PeakAnnotFile, FeatureLinkageFile, verbose=True)

# Let's examine the total count distribution and remove outliers.
plt.hist(adata_atac.X.sum(1), bins=100, range=(0, 100000));
histfile = OutDir + "/1_histogram_count_distribution.png"
plt.savefig(histfile)

## filter ATAC data
if 0:
	sc.pp.filter_cells(adata_atac, min_counts=ATAC_COUNT_MIN_THR)
	sc.pp.filter_cells(adata_atac, max_counts=ATAC_COUNT_MAX_THR)

# We normalize aggregated peaks with TF-IDF.
mv.tfidf_norm(adata_atac)

##=====================
## Finding shared barcodes and features between RNA and ATAC
##=====================
print("length adata_rna.obs_names : " + str(len(adata_rna.obs_names)))
print("length adata_atac.obs_names : " + str(len(adata_atac.obs_names)))

shared_cells = pd.Index(np.intersect1d(adata_rna.obs_names, adata_atac.obs_names))
shared_genes = pd.Index(np.intersect1d(adata_rna.var_names, adata_atac.var_names))
print("No of shared cells : " + str(len(shared_cells)))
print("No of shared genes : " + str(len(shared_genes)))

# We reload in the raw data and continue with a subset of cells.
adata_rna = scv.read(LoomFile, cache=True)
adata_rna.obs_names = [x.split(':')[1][:-1] + '-1' for x in adata_rna.obs_names]
adata_rna.var_names_make_unique()

adata_rna = adata_rna[shared_cells, shared_genes]
adata_atac = adata_atac[shared_cells, shared_genes]

scv.pp.normalize_per_cell(adata_rna)
scv.pp.log1p(adata_rna)
scv.pp.moments(adata_rna, n_pcs=NUM_PCS_THR, n_neighbors=NUM_NEIGHBORS_THR)

adata_rna.obs['celltype'] = cell_annot.loc[adata_rna.obs_names, 'celltype']
adata_rna.obs['celltype'] = adata_rna.obs['celltype'].astype('category')

# Reorder the categories for color consistency with the manuscript.
if 0:
	all_clusters = ['Upper Layer', 'Deeper Layer', 'V-SVZ', 'RG, Astro, OPC', 'Ependymal cells', 'IPC', 'Subplate']
	adata_rna.obs['celltype'] = adata_rna.obs['celltype'].cat.reorder_categories(all_clusters)

## UMAP computation
scv.tl.umap(adata_rna)

## plot UMAP
scv.pl.umap(adata_rna, color='celltype')
umapfile = OutDir + "/0_UMAP_Celltype.png"
plt.savefig(umapfile)

##==============================
## Smoothing gene aggregagted peaks by neighbors
##==============================
# Write out filtered cells
FilteredCellFile = OutDir + "/filtered_cells.txt"
adata_rna.obs_names.to_frame().to_csv(FilteredCellFile, header=False, index=False)

# Read in Seurat WNN neighbors.
nn_idx = np.loadtxt(nn_idx_file, delimiter=',')
nn_dist = np.loadtxt(nn_dist_file, delimiter=',')
nn_cells = pd.Index(pd.read_csv(nn_cells_file, header=None)[0])

# Make sure cell names match.
np.all(nn_cells == adata_atac.obs_names)

mv.knn_smooth_chrom(adata_atac, nn_idx, nn_dist)

adata_atac


##====================
## Running multi-omic dynamical model

## MultiVelo incorporates chromatin accessibility information into RNA velocity and achieves better lineage predictions.
## The detailed argument list can be shown with "help(mv.recover_dynamics_chrom)".

multiomic_outfile = OutDir + "/multivelo_result.h5ad"
if not os.path.exists(multiomic_outfile):
	# This will take a while. Parallelization is high recommended.
	adata_result = mv.recover_dynamics_chrom(adata_rna, adata_atac, max_iter=5, init_mode="invert", verbose=False, parallel=True, save_plot=False, rna_only=False, fit=True, n_anchors=500, extra_color_key='celltype')
	# Save the result for use later on
	adata_result.write(multiomic_outfile)

# read the existing results
adata_result = sc.read_h5ad(multiomic_outfile)

mv.pie_summary(adata_result)
piefile = OutDir + "/2_piechart_multivelo_summary.png"
plt.savefig(piefile)

mv.switch_time_summary(adata_result)
swtchsummaryfile = OutDir + "/3_boxplot_switch_time_summary.png"
plt.savefig(swtchsummaryfile)

mv.likelihood_plot(adata_result)
likplotfile = OutDir + "/4_likelihood_plot_summary.png"
plt.savefig(likplotfile)

## By default, the velocity genes used for velocity graph is determined as those whose likelihoods are above 0.05. 
## They can be reset with "mv.set_velocity_genes" function upon examining the distributions of variables above if needed.

##=====================
## Computing velocity stream and latent time
##=====================
mv.velocity_graph(adata_result)
mv.latent_time(adata_result)
mv.velocity_embedding_stream(adata_result, basis='umap', color='celltype')
veloembfile = OutDir + "/5_velocity_embedding_stream.png"
plt.savefig(veloembfile)

scv.pl.scatter(adata_result, color='latent_time', color_map='gnuplot', size=20)	#80)
sctrfile = OutDir + "/6_velocity_scatter.png"
plt.savefig(sctrfile)

# ##=====================
# ## examine some genes
# genedf = pd.read_table(TargetGenelistFile, header=None)
# gene_list = genedf[0]	# get the first column - gene names

# ##=================
# ## plot accessbility and expression against gene time.
# ##=================
# # Accessibility/expression by gene time, colored by the four potential states.
# # The solid black curve indicates anchors.
# mv.dynamic_plot(adata_result, gene_list, color_by='state', axis_on=False, frame_on=False)
# outfile = OutDir + "/7_accessibility_or_expression_against_gene_time.png"
# plt.savefig(outfile)

# ##===============
# ## plot velocity against gene time
# ##===============
# # Velocity by gene time, colored by the four potential states.
# # The solid black curve indicates anchors.
# mv.dynamic_plot(adata_result, gene_list, color_by='state', by='velocity', axis_on=False, frame_on=False)
# outfile = OutDir + "/8_velocity_against_gene_time.png"
# plt.savefig(outfile)

# ##===============
# ## accessibility and expression against globally shared latent time.
# ##===============
# # Accessibility/expression by global latent time, colored by cell type assignments.
# # The solid black curve indicates the mean.
# mv.dynamic_plot(adata_result, gene_list, color_by='celltype', gene_time=False, axis_on=False, frame_on=False)
# outfile = OutDir + "/9_accessibility_or_expression_against_latent_time.png"
# plt.savefig(outfile)

# ##===============
# ## Phase portraits on the unspliced-spliced, chromatin-unspliced, 
# ## or 3-dimensional planes can be plotted.
# ##===============
# # Unspliced-spliced phase portraits, colored by celltype.
# mv.scatter_plot(adata_result, gene_list, color_by='celltype', by='us', axis_on=False, frame_on=False)
# outfile = OutDir + "/10_Unspliced_spliced_phase_portraits_by_celltype.png"
# plt.savefig(outfile)
# # Unspliced-spliced phase portraits, colored by log chromatin accessibility.
# # title_more_info shows more information in each subplot title: model, direction, and likelihood.
# mv.scatter_plot(adata_result, gene_list, color_by='c', by='us', cmap='coolwarm', title_more_info=True, axis_on=False, frame_on=False)
# outfile = OutDir + "/11_Unspliced_spliced_phase_portraits_by_log_chromatin_accessibility.png"
# plt.savefig(outfile)
# # Chromatin-unspliced phase portraits, colored by celltype.
# mv.scatter_plot(adata_result, gene_list, color_by='celltype', by='cu', axis_on=False, frame_on=False)
# outfile = OutDir + "/12_Chromatin_Unspliced_phase_portraits_by_celltype.png"
# plt.savefig(outfile)
# # 3D phase portraits, colored by celltype.
# mv.scatter_plot(adata_result, gene_list, color_by='celltype', by='cus', axis_on=False, downsample=2)
# outfile = OutDir + "/13_3D_phase_portraits_by_celltype.png"
# plt.savefig(outfile)

# ##=========================
# ## add velocity arrows to phase portraits to show the predicted directions.
# ##=========================
# mv.scatter_plot(adata_result, gene_list, color_by='celltype', by='us', axis_on=False, frame_on=False, downsample=2, velocity_arrows=True)
# outfile = OutDir + "/14_Unspliced_spliced_velocity_arrows_by_celltype.png"
# plt.savefig(outfile)

# mv.scatter_plot(adata_result, gene_list, color_by='celltype', by='cu', axis_on=False, frame_on=False, downsample=2, velocity_arrows=True)
# outfile = OutDir + "/15_chromatin_Unspliced_velocity_arrows_by_celltype.png"
# plt.savefig(outfile)

# mv.scatter_plot(adata_result, gene_list, color_by='celltype', by='cus', downsample=3, velocity_arrows=True)
# outfile = OutDir + "/16_3D_phase_portraits_velocity_arrows_by_celltype.png"
# plt.savefig(outfile)
