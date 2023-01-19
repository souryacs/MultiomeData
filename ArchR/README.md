Multi-omic data (single cell ATAC-seq + single cell RNA-seq) analysis using ArchR
===================================================================================

This script "ArchR_Code.r" needs to be executed using the command:

	Rscript ArchR_Code.r *ATACPATH* *scRNASeqObj* *CurrOutDir*

	where,

		*ATACPATH*: user will provide the scATAC-seq data path (output of CellRanger) as an input to this pipeline. For details of the input scATAC-seq data for the ArchR pipeline, refer to the ArchR documentation.

		*scRNASeqObj*: user will input the scRNA-seq processed data (Seurat object in .rds format) as the second input. User can check my uploaded Seurat pipeline regarding processing of scRNA-seq datasets.

		*CurrOutDir*: output directory to contain all outputs from ArchR.

Installation
==============

	We recommend users to install R version >= 4.1.0 first, and then install the following R packages:

	ArchR; ggplot2; rtracklayer

	User also needs to install the "MACS2" peak caller (https://github.com/macs3-project/MACS)

	Also, we suggest users to look at the ArchR documentation:

	brief manual: https://www.archrproject.com/articles/Articles/tutorial.html
	full manual: https://www.archrproject.com/bookdown/index.html


Data Description
=================

	1) Here we assume that the input scATAC-seq and scRNA-seq data has reference genome hg38.

		** For any other reference genome, check the line 6 of the R script and apply corresponding argument to the function "addArchRGenome".

	2) We assume that batch correction is disabled, by setting the parameter batch_correction = FALSE (line 31). User can enable the batch correction, if required.

		** Batch correction is performed with Harmony


Future works
================

	1) The code does not implement the constrained integration between scATAC-seq and scRNA-seq 

		** For details, check the reference manual, and check the code space in lines 231 - 239.

	2) The footprinting option currently plots footprints for a selected set of genes - "GATA1", "CEBPA", "EBF1", "IRF4", "TBX21", "PAX5" (check line 537) - user can edit this gene list.

	3) Trajectory inference is not incorporated.


Output files
===============

	All the output files are placed within the mentioned output directory *CurrOutDir*. User needs to go through the detailed documentation of ArchR in order to understand individual output files.


Queries
=======

For any query, please e-mail:
Sourya Bhattacharyya

sourya@lji.org

