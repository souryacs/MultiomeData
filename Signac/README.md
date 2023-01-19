Multi-omic data (single cell ATAC-seq + single cell RNA-seq) analysis using Signac
===================================================================================

This script "Signac_Code.r" needs to be executed using the command:

	Rscript Signac_Code.r *ATACFile* *scRNASeqFile* *CurrOutDir* *GeneListFile*

	where,

		*ATACFile*: user will provide the scATAC-seq data file (output of CellRanger Arc) of the name "atac_fragments.tsv.gz" as an input. For details of the input scATAC-seq data for the Signac pipeline, refer to the Signac documentation.

		*scRNASeqFile*: user will input the gene expression file output from CellRanger Arc (filtered_feature_bc_matrix.h5) as the second input.

		*CurrOutDir*: output directory to contain all outputs from Signac.

		*GeneListFile*: Optional parameter. File containing a target gene list (one entry per row) whose tracks / features etc. to be analyzed. According to the characteristics of the input data and reference genome, user needs to put the genes of interest whose properties need to be looked at.


	* Note: before execution, first install the required packages (see Installation) and then check the *Code edits / important* section.


Signac Documentation
======================

	User needs to look at the Signac documentation: https://satijalab.org/signac/


Installation
==============

	We recommend users to install R version >= 4.1.0 first, and then install the following R packages:

		Signac, Seurat, ggplot2, JASPAR2020, TFBSTools, rtracklayer, SeuratWrappers, monocle3, Matrix, patchwork, motifmatchr.


	Also, user needs to edit the line 42 of the R script for any reference genome other than mm10.

	

Parameters / edits 
========================

	1) Here we assume that the input scATAC-seq and scRNA-seq data has reference genome mm10.

		The code uses the following R libraries for the mm10 reference genome:
	
			EnsDb.Mmusculus.v79, BSgenome.Mmusculus.UCSC.mm10

		Also, in line 50 of the R script, the current code is:

			annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

		For any other reference genome, 

			1) User needs to install similar R libraries (EnsDb.*, BSgenome.*).

			2) User needs to edit the line 50 of the R script.


	2) User needs to view the parameters section (lines 18 - 89) of the R script and edit them if required.

		** One of the important parameters is the clustering resolution (see line 44) - currently it is set as 0.3 - corresponding cluster and signac object are used for the subsequent analysis. 

		*** We report clustering results for different resolutions (check line 281) so user can pick the resolution accordingly. However, user needs to edit the metadata fields like "SCT_snn_res.0.3"	used throughout the code according to the selected resolution.

	3) The line 47 contains a parameter *Downstream_Analysis* initialized as FALSE. Initially, when the user runs this pipeline, multi-omic data will be processed and the clusters will be generated for various resolutions (see below). 

		** Once user fixes a particular resolution, 
			a) he / she needs to edit this clustering resolution parameter (line 44) - currently it is 0.3
			b) edit the fields of metadata (search the strings like "SCT_snn_res.0.3" throughout the R code)
			c) edit this *Downstream_Analysis* parameter as TRUE, and re-run this pipeline.


Output
========

	1) For the first pass (*Downstream_Analysis* = FALSE), the UMAP plots *Cluster_UMAP_resolution_*.pdf* and corresponding Signac objects *Signac_Res_*.rds* are created. 

		User needs to look at these UMAP plots, in order to select the best resolution, and edit the parameters of lines 44 and 47.

	2) Once the resolution is fixed and *Downstream_Analysis* is TRUE, re-running the pipeline will create other output files. For details of these output files, check the Signac documentation.


Future works
================

	1) Trajectory inference is not incorporated.

	2) Footprint plots are not enabled.


Output files
===============

	All the output files are placed within the mentioned output directory *CurrOutDir*. User needs to go through the detailed documentation of Signac in order to understand individual output files.


Queries
=======

For any query, please e-mail:
Sourya Bhattacharyya

sourya@lji.org

