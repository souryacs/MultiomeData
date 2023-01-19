library(Signac)
library(Seurat)
library(ggplot2)
library(JASPAR2020)
library(TFBSTools)
library(rtracklayer)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(patchwork)
# library(cicero)
library(motifmatchr)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)

set.seed(1234)

##=======================
## paramters
##=======================

args <- commandArgs(trailingOnly = TRUE)
fragpath_ATAC <- args[1]
RNA_h5_file <- args[2]
Main_OutDir <- args[3]
if (length(args) > 3) {
  GeneListFile <- args[4]
}

system(paste("mkdir -p", Main_OutDir))

## output log file
sink(paste0(Main_OutDir, "/out_log_cell_Count.txt"))

## Parameters
nCount_RNA_lower <- 1000
nCount_RNA_upper <- 25000
nCount_ATAC_lower <- 1000
nCount_ATAC_upper <- 100000
nucl_signal_upper <- 2
tss_enrich_lower <- 1

## cluster resolution
res_val <- 0.3

## this boolean value indicates if the downstream analysis will take place
Downstream_Analysis <- FALSE

# get gene annotations for mm10
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"

if (length(args) > 3) {
  ## read the target genelist
  GeneListData <- read.table(GeneListFile, header=F, sep="\t", stringsAsFactors=F)
  GeneList <- GeneListData[,1]
  cat(sprintf("\n Input number of genes : %s ", length(GeneList)))
}

QC_OutDir <- paste0(Main_OutDir, '/QC')
system(paste("mkdir -p", QC_OutDir))

Peak_OutDir <- paste0(Main_OutDir, '/Peaks')
system(paste("mkdir -p", Peak_OutDir))

macs2path <- as.character(system("which macs2", intern = TRUE))
cat(sprintf("\n\n ==>> macs2 path : %s ", macs2path))

Cluster_OutDir <- paste0(Main_OutDir, '/Cluster_UMAP')
system(paste("mkdir -p", Cluster_OutDir))

Motif_OutDir <- paste0(Main_OutDir, '/Motif')
system(paste("mkdir -p", Motif_OutDir))

Trajectory_OutDir <- paste0(Main_OutDir, '/Trajectory_Monocle3')
system(paste("mkdir -p", Trajectory_OutDir))

Cicero_OutDir <- paste0(Main_OutDir, '/CoAccessibility_Cicero')
system(paste("mkdir -p", Cicero_OutDir))

Footprint_OutDir <- paste0(Main_OutDir, '/TF_Footprint')
system(paste("mkdir -p", Footprint_OutDir))

tempAnnotationFile <- paste0(Main_OutDir, '/TSS_Positions.bed')

## stores peaks for all clusters
targetpeakfile <- paste0(Peak_OutDir, '/Peaks_Categorized_by_celltype.bed')

SignacObjFile <- paste0(Main_OutDir, '/Signac_Res_', res_val, '.rds')

##=======================
## main code
##=======================

if (file.exists(SignacObjFile) == FALSE) {

  counts <- Read10X_h5(RNA_h5_file)

  # create a Seurat object containing the RNA adata
  sc.obj <- CreateSeuratObject(counts=counts$`Gene Expression`, assay="RNA")
  # create ATAC assay and add it to the object
  sc.obj[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks, sep = c(":", "-"), fragments = fragpath_ATAC, annotation = annotation)
  cat(sprintf("\n\n Read input data sc.obj - number of cells : %s ", ncol(sc.obj)))

  ##================
  ## Quality control
  ##================
  DefaultAssay(sc.obj) <- "ATAC"

  ## Calculate the strength of the nucleosome signal per cell. 
  ## Computes the ratio of fragments between 147 bp and 294 bp (mononucleosome) 
  ## to fragments < 147 bp (nucleosome-free)
  sc.obj <- NucleosomeSignal(sc.obj)

  ## Compute the transcription start site (TSS) enrichment score for each cell, 
  ## as defined by ENCODE: https://www.encodeproject.org/data-standards/terms/.
  sc.obj <- TSSEnrichment(sc.obj)

  ## Visualize QC metrics as a violin plot
  currplot <- VlnPlot(
    object = sc.obj,
    features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
    ncol = 4,
    pt.size = 0
  )
  ggplot2::ggsave(paste0(QC_OutDir, "/QC_Violin_plots.pdf"), plot=currplot, width=10, height=8)

  ## filter out low quality cells
  sc.obj <- subset(
    x = sc.obj,
    subset = nCount_ATAC < nCount_ATAC_upper &
      nCount_RNA < nCount_RNA_upper &
      nCount_ATAC > nCount_ATAC_lower &
      nCount_RNA > nCount_RNA_lower &
      nucleosome_signal < nucl_signal_upper &
      TSS.enrichment > tss_enrich_lower
  )
  # sc.obj
  cat(sprintf("\n\n After filtering - number of cells : %s ", ncol(sc.obj)))

  ##===================
  ## subset ATAC-seq object
  ## use only the standard chromosomes

  ## This is crucial, since the motifs are from all the chromosomes (including random chromosomes)
  ## whereas the ATAC-seq peaks is already defined for the valid chromosomes

  ## subsetting ATAC-seq data is provided here:
  ## https://github.com/timoast/signac/issues/486
  ##===================  
  main.chroms <- standardChromosomes(BSgenome.Mmusculus.UCSC.mm10)
  keep.peaks <- which(as.character(seqnames(granges(sc.obj))) %in% main.chroms)
  sc.obj[["ATAC"]] <- subset(sc.obj[["ATAC"]], features = rownames(sc.obj[["ATAC"]])[keep.peaks])

  ##================
  ## peak calling
  ##================

  ## call peaks using MACS2
  ## Note: these peaks are called using all the cells
  peaks <- CallPeaks(sc.obj, 
    macs2.path = macs2path, 
    outdir = Peak_OutDir, 
    fragment.tempdir = Peak_OutDir)

  # remove peaks on nonstandard chromosomes and in genomic blacklist regions
  peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
  peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)

  ## export the generated peaks (GRanges object)
  export.bed(peaks, paste0(Peak_OutDir, '/Peaks_from_All_cells.bed'))

  # quantify counts in each peak
  macs2_counts <- FeatureMatrix(
    fragments = Fragments(sc.obj),
    features = peaks,
    cells = colnames(sc.obj)
  )

  # ## dump the feature matrix
  # write.table(as.data.frame(macs2_counts), paste0(Peak_OutDir, '/Peaks_from_All_cells_Feature_Matrix.txt'), row.names=F, col.names=T, sep="\t", quote=F, append=F)

  # create a new assay using the MACS2 peak set and add it to the Seurat object
  sc.obj[["peaks"]] <- CreateChromatinAssay(
    counts = macs2_counts,
    # fragments = fragpath,
    fragments = Fragments(sc.obj),
    annotation = annotation
  )

  cat(sprintf("\n\n Number of cells in Seurat object : %s ", ncol(sc.obj)))

  ##===================
  ## Gene expression data processing
  ## also see the UMAPs solely generated from the RNA-seq data
  ##===================

  ## We can normalize the gene expression data using SCTransform, 
  ## and reduce the dimensionality using PCA.
  DefaultAssay(sc.obj) <- "RNA"
  sc.obj <- SCTransform(sc.obj)
  sc.obj <- RunPCA(sc.obj)

  DefaultAssay(sc.obj) <- "SCT"

  ## UMAP visualization
  ## also we have put a "reduction.name" value
  sc.obj <- RunUMAP(object=sc.obj, reduction = 'pca', dims=1:50, reduction.name = "umap.rna", verbose = TRUE)
  sc.obj <- FindNeighbors(sc.obj, reduction = 'pca', dims = 1:50)
  ## test clustering with various resolutions 
  ## for RNA-seq data
  for (res_val in c(0.3, 0.4, 0.5, 0.8, 1, 1.5, 2)) {
    sc.obj <- FindClusters(sc.obj, resolution = res_val, algorithm = 3)
    p1 <- DimPlot(sc.obj, reduction="umap.rna", label = TRUE, label.size = 6) + ggtitle(paste0("RNA UMAP Resolution ", res_val))
    ggplot2::ggsave(paste0(Cluster_OutDir, "/Cluster_UMAP_only_RNA_resolution_", res_val, ".pdf"), plot=p1, width=6, height=6)
  }

  ##===================
  ## DNA accessibility data processing
  ##===================
  DefaultAssay(sc.obj) <- "peaks"

  ## Find most frequently observed features
  sc.obj <- FindTopFeatures(sc.obj, min.cutoff = 5)

  ## Compute the term-frequency inverse-document-frequency
  sc.obj <- RunTFIDF(sc.obj)

  ## dimension reduction using SVD
  sc.obj <- RunSVD(sc.obj)

  ## check https://satijalab.org/signac/articles/pbmc_vignette.html

  ## UMAP visualization for the ATAC-seq data
  sc.obj <- RunUMAP(sc.obj, reduction = 'lsi', dims = 2:40, reduction.name = 'umap.atac', verbose = TRUE)
  sc.obj <- FindNeighbors(object = sc.obj, reduction = 'lsi', dims = 2:40)
  ## test clustering with various resolutions 
  ## for ATAC-seq data
  for (res_val in c(0.3, 0.4, 0.5, 0.8, 1, 1.5, 2)) {
    sc.obj <- FindClusters(object = sc.obj, resolution = res_val, algorithm = 3)
    p1 <- DimPlot(sc.obj, reduction="umap.atac", label = TRUE, label.size = 6) + ggtitle(paste0("ATAC UMAP Resolution ", res_val))
    ggplot2::ggsave(paste0(Cluster_OutDir, "/Cluster_UMAP_only_ATAC_resolution_", res_val, ".pdf"), plot=p1, width=6, height=6)
  }

  ##===================
  ## ***** Clustering using both RNA and ATAC *****
  ## check weighted nearest neighbor analysis
  ## https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html

  ## also check the Signac Vignette
  ## https://satijalab.org/signac/articles/pbmc_multiomic.html
  ##===================
  ## change the default assay
  ## the gene expression data is already normalized by SCTransform
  DefaultAssay(sc.obj) <- "SCT"

  ## build a joint neighbor graph using both assays
  ## this is a Seurat function
  sc.obj <- FindMultiModalNeighbors(
    object = sc.obj,
    reduction.list = list("pca", "lsi"), ## A list of two dimensional reductions, one for each of the modalities to be integrated
    dims.list = list(1:50, 2:40), ## A list containing the dimensions for each reduction to use
    modality.weight.name = "RNA.weight",  ## Variable name to store modality weight in object meta data
    verbose = TRUE
  )

  ## build a joint UMAP visualization
  ## We've omitted assay = "RNA" parameter
  ## use the current default assay
  sc.obj <- RunUMAP(
    object = sc.obj,
    nn.name = "weighted.nn",  ## this is the default value of "weighted.nn.name" parameter (Multimodal neighbor object name) in the FindMultiModalNeighbors function
    reduction.name = "wnn.umap", 
    reduction.key = "wnnUMAP_",
    verbose = TRUE
  )

  ## test clusters with multiple resolution values
  ## we have used weighted nearest neighbor concept
  ## https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html
  for (res_val in c(0.3, 0.4, 0.5, 0.8, 1, 1.5, 2)) {
    ## find clusters from the UMAP
    sc.obj <- FindClusters(object = sc.obj, graph.name = "wsnn", algorithm = 3, resolution = res_val, verbose = TRUE)
    ## print cluster labels
    Idents(object = sc.obj)
    levels(x = sc.obj)
    ## store cell identity classes
    sc.obj[["clustID"]] <- Idents(object = sc.obj)
    ## plot UMAP - check "reduction" parameter
    currplot <- DimPlot(sc.obj, group.by = "clustID", label = TRUE, repel = TRUE, reduction = "wnn.umap", label.size = 6) #+ NoLegend()
    ggplot2::ggsave(paste0(Cluster_OutDir, "/Cluster_UMAP_resolution_", res_val, ".pdf"), plot=currplot, width=6, height=6)
    ## save the Seurat object
    saveRDS(sc.obj, file = paste0(Main_OutDir, '/Signac_Res_', res_val, '.rds'))
  }

} else {
  sc.obj <- readRDS(SignacObjFile)
}

##=====================
## clusterwise cell count
##=====================
if (Downstream_Analysis == TRUE) {

  ## RNA assay
  DefaultAssay(sc.obj) <- "SCT"
  all.clusters <- unique(sc.obj@meta.data[, c("SCT_snn_res.0.3")])
  cat(sprintf("\n\n ==>>> RNA-seq - number of clusters: %s ", length(all.clusters)))
  for (cluster in all.clusters) {
    cat(sprintf("\n Number of cells in (RNA-seq) cluster : %s is %s ", cluster, length(which(sc.obj@meta.data[, c("SCT_snn_res.0.3")] == cluster))))
  }

  ## multimodal assay
  DefaultAssay(sc.obj) <- "SCT"
  all.clusters <- unique(sc.obj@meta.data[, c("wsnn_res.0.3")])
  cat(sprintf("\n\n ==>>> Multimodal - number of clusters: %s ", length(all.clusters)))
  for (cluster in all.clusters) {
    cat(sprintf("\n Number of cells in (multimodal) cluster : %s is %s ", cluster, length(which(sc.obj@meta.data[, c("wsnn_res.0.3")] == cluster))))  
  }

} # end if

##=============================
## marker gene for clusters in RNA-seq data
##=============================
if (Downstream_Analysis == TRUE) {

  DefaultAssay(sc.obj) <- "SCT"
  MarkerGeneOutDir_RNA <- paste0(Main_OutDir, '/MarkerGene_RNA')
  system(paste("mkdir -p", MarkerGeneOutDir_RNA))
  nclust <- length(unique(as.integer(sc.obj@meta.data[, c("SCT_snn_res.0.3")])))
  cat(sprintf("\n RNA-seq - number of clusters: %s ", nclust))
  all.clusters <- (seq(1,nclust) - 1)

  ## marker genes for individual clusters
  for (cluster in all.clusters) {
    if (length(WhichCells(sc.obj, ident=cluster)) > 10) {
      cat(sprintf("\n ---- analyzing one cluster --- %s ", cluster))
      cluster.markers <- FindMarkers(object=sc.obj, group.by="SCT_snn_res.0.3", ident.1=cluster, only.pos=FALSE, min.pct=0.20)
      tmp.file.name <- paste0(MarkerGeneOutDir_RNA, '/MarkerGene_', cluster, '.csv')
      write.csv(x=cluster.markers, file=tmp.file.name)
    }
  } # end cluster loop

  ## marker genes for pairwise clusters
  for (cluster in all.clusters) {
    if (length(WhichCells(sc.obj, ident=cluster)) > 10) {
      CurrOutDir <- paste0(MarkerGeneOutDir_RNA, '/Pairwise_Analysis_cluster_', cluster)
      system(paste("mkdir -p", CurrOutDir))  
      rest.clusters <- all.clusters[all.clusters != cluster] 
      cat(sprintf("\n ---- analyzing one cluster --- %s ", cluster))
      for (next.cluster in all.clusters) {
        if (cluster == next.cluster) {
          next
        }
        if (length(WhichCells(sc.obj, ident=next.cluster)) > 10) {
          cat(sprintf("\n -- processing pair -- cluster : %s next.cluster : %s ", cluster, next.cluster))
          cluster.markers <- FindMarkers(object=sc.obj, group.by="SCT_snn_res.0.3", ident.1=cluster, ident.2=next.cluster, only.pos=FALSE, min.pct=0.20)
          tmp.file.name <- paste0(CurrOutDir, '/', cluster, '_vs_', next.cluster, '.csv')
          write.csv(x=cluster.markers, file=tmp.file.name)
        }
      } # end next.cluster loop
    }
  } # end cluster loop

} # end dummy if

##=============================
## marker gene for clusters in ATAC-seq data
##=============================
if (Downstream_Analysis == TRUE) {

  DefaultAssay(sc.obj) <- "peaks"
  MarkerGeneOutDir_ATAC <- paste0(Main_OutDir, '/MarkerGene_ATAC')
  system(paste("mkdir -p", MarkerGeneOutDir_ATAC))
  nclust <- length(unique(as.integer(sc.obj@meta.data[, c("peaks_snn_res.0.3")])))
  cat(sprintf("\n ATAC-seq - number of clusters: %s ", nclust))
  all.clusters <- (seq(1,nclust) - 1)

  ## marker genes for individual clusters
  for (cluster in all.clusters) {
    if (length(WhichCells(sc.obj, ident=cluster)) > 10) {
      cat(sprintf("\n ---- analyzing one cluster --- %s ", cluster))
      cluster.markers <- FindMarkers(object=sc.obj, group.by="peaks_snn_res.0.3", ident.1=cluster, only.pos=FALSE, min.pct=0.20)
      tmp.file.name <- paste0(MarkerGeneOutDir_ATAC, '/MarkerGene_', cluster, '.csv')
      write.csv(x=cluster.markers, file=tmp.file.name)
    }
  } # end cluster loop

  ## marker genes for pairwise clusters
  for (cluster in all.clusters) {
    if (length(WhichCells(sc.obj, ident=cluster)) > 10) {
      CurrOutDir <- paste0(MarkerGeneOutDir_ATAC, '/Pairwise_Analysis_cluster_', cluster)
      system(paste("mkdir -p", CurrOutDir))  
      rest.clusters <- all.clusters[all.clusters != cluster] 
      cat(sprintf("\n ---- analyzing one cluster --- %s ", cluster))
      for (next.cluster in all.clusters) {
        if (cluster == next.cluster) {
          next
        }
        if (length(WhichCells(sc.obj, ident=next.cluster)) > 10) {
          cat(sprintf("\n -- processing pair -- cluster : %s next.cluster : %s ", cluster, next.cluster))
          cluster.markers <- FindMarkers(object=sc.obj, group.by="peaks_snn_res.0.3", ident.1=cluster, ident.2=next.cluster, only.pos=FALSE, min.pct=0.20)
          tmp.file.name <- paste0(CurrOutDir, '/', cluster, '_vs_', next.cluster, '.csv')
          write.csv(x=cluster.markers, file=tmp.file.name)
        }
      } # end next.cluster loop
    }
  } # end cluster loop

} # end dummy if


##=============================
## marker gene for clusters in Multimodal data
##=============================
if (Downstream_Analysis == TRUE) {

  DefaultAssay(sc.obj) <- "SCT"
  MarkerGeneOutDir_Multimodal <- paste0(Main_OutDir, '/MarkerGene_Multimodal')
  system(paste("mkdir -p", MarkerGeneOutDir_Multimodal))
  nclust <- length(unique(as.integer(sc.obj@meta.data[, c("wsnn_res.0.3")])))
  cat(sprintf("\n Multimodal - number of clusters: %s ", nclust))
  all.clusters <- (seq(1,nclust) - 1)

  ## marker genes for individual clusters
  for (cluster in all.clusters) {
    if (length(WhichCells(sc.obj, ident=cluster)) > 10) {
      cat(sprintf("\n ---- analyzing one cluster --- %s ", cluster))
      cluster.markers <- FindMarkers(object=sc.obj, group.by="wsnn_res.0.3", ident.1=cluster, only.pos=FALSE, min.pct=0.20)
      tmp.file.name <- paste0(MarkerGeneOutDir_Multimodal, '/MarkerGene_', cluster, '.csv')
      write.csv(x=cluster.markers, file=tmp.file.name)
    }
  } # end cluster loop

  ## marker genes for pairwise clusters
  for (cluster in all.clusters) {
    if (length(WhichCells(sc.obj, ident=cluster)) > 10) {
      CurrOutDir <- paste0(MarkerGeneOutDir_Multimodal, '/Pairwise_Analysis_cluster_', cluster)
      system(paste("mkdir -p", CurrOutDir))  
      rest.clusters <- all.clusters[all.clusters != cluster] 
      cat(sprintf("\n ---- analyzing one cluster --- %s ", cluster))
      for (next.cluster in all.clusters) {
        if (cluster == next.cluster) {
          next
        }
        if (length(WhichCells(sc.obj, ident=next.cluster)) > 10) {
          cat(sprintf("\n -- processing pair -- cluster : %s next.cluster : %s ", cluster, next.cluster))
          cluster.markers <- FindMarkers(object=sc.obj, group.by="wsnn_res.0.3", ident.1=cluster, ident.2=next.cluster, only.pos=FALSE, min.pct=0.20)
          tmp.file.name <- paste0(CurrOutDir, '/', cluster, '_vs_', next.cluster, '.csv')
          write.csv(x=cluster.markers, file=tmp.file.name)
        }
      } # end next.cluster loop
    }
  } # end cluster loop

} # end dummy if

##===================
## Call peaks per cluster
##===================
if (Downstream_Analysis == TRUE) {

  DefaultAssay(sc.obj) <- "ATAC"

  if (file.exists(targetpeakfile) == FALSE) {
    peaks_clust <- CallPeaks(object = sc.obj, group.by = "clustID", macs2.path = macs2path, outdir = Peak_OutDir, fragment.tempdir = Peak_OutDir)
    # remove peaks on nonstandard chromosomes and in genomic blacklist regions
    peaks_clust <- keepStandardChromosomes(peaks_clust, pruning.mode = "coarse")
    peaks_clust <- subsetByOverlaps(x = peaks_clust, ranges = blacklist_mm10, invert = TRUE)
    ## export the generated peaks (GRanges object)
    export.bed(peaks_clust, targetpeakfile)  
  }

  ## We can visualize the cell-type-specific MACS2 peak calls alongside the 10x Cellranger peak calls 
  ## (currently being used in the sc.obj object) with the CoveragePlot() function. 
  ## Here the Cellranger peaks are shown in grey and the MACS2 peaks in red:
  if (0) {
    ## we plot near a sample target gene 'CD8A'
    currplot <- CoveragePlot(object = sc.obj, region = "CD8A", ranges = peaks_clust, ranges.title = "MACS2")
    ggplot2::ggsave(paste0(Peak_OutDir, "/Peak_CoveragePlot_CD8A_Gene.pdf"), plot=currplot, width=10, height=8)
  }  

  saveRDS(sc.obj, file = SignacObjFile)

} # end dummy if

##===================
## Get genes linked to peaks
## and peaks linked to genes
##===================
if (Downstream_Analysis == TRUE) {

  DefaultAssay(sc.obj) <- "peaks"

  # first compute the GC content for each peak
  sc.obj <- RegionStats(sc.obj, genome = BSgenome.Mmusculus.UCSC.mm10)

  # link peaks to genes
  ## Note: we use the peaks generated by all the cells
  # sc.obj <- LinkPeaks(object = sc.obj, peak.assay = "peaks", expression.assay = "SCT", genes.use = c("LYZ", "MS4A1"))
  sc.obj <- LinkPeaks(object = sc.obj, peak.assay = "peaks", expression.assay = "SCT")

  ## save this signac object
  saveRDS(sc.obj, file = SignacObjFile)
  cat(sprintf("\n Saved Signac object after link peaks"))

  ## export this links information
  linkinfo_file <- paste0(Peak_OutDir, '/Complete_Link_Info.bed')
  linkinfo <- Links(sc.obj)
  export.bed(linkinfo, linkinfo_file)  

} # end dummy if

##=======
## Use ATAC-seq data to predict gene activity score
## first get the gene activity scores
##=======
if (Downstream_Analysis == TRUE) {

  DefaultAssay(sc.obj) <- "ATAC"

  # compute gene activity score using chromatin accessibility data
  gene.activities <- GeneActivity(sc.obj)

  # add to the Signac object as a new assay
  sc.obj[['GeneAct']] <- CreateAssayObject(counts = gene.activities)

  sc.obj <- NormalizeData(
    object = sc.obj,
    assay = 'GeneAct',
    normalization.method = 'LogNormalize',
    scale.factor = median(sc.obj$nCount_RNA)
  )

  saveRDS(sc.obj, file = SignacObjFile)

} # end dummy if

##===================
## Motif analysis with Signac
## To facilitate motif analysis in Signac, we have create the Motif class to store all the required information, 
## including a list of position weight matrices (PWMs) or position frequency matrices (PFMs) 
## and a motif occurrence matrix.
##===================
if (Downstream_Analysis == TRUE) {

  # DefaultAssay(sc.obj) <- "ATAC"
  DefaultAssay(sc.obj) <- "peaks"

  ## Get a list of motif position frequency matrices from the JASPAR database
  ## check the options
  ## for Human, use species = 'Homo sapiens'
  ## for mouse, use tax_group = 'vertebrates'

  ## mouse
  pfm <- getMatrixSet(
    x = JASPAR2020,
    opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
  )

  cat(sprintf("\n\n After the function getMatrixSet \n\n"))

  # add motif information
  ## To add the DNA sequence motif information required for motif analyses, we can run the AddMotifs() function:
  sc.obj <- AddMotifs(
    object = sc.obj,
    genome = BSgenome.Mmusculus.UCSC.mm10,
    pfm = pfm
  )

  cat(sprintf("\n\n After the function AddMotifs \n\n"))

  ## save the Seurat object
  saveRDS(sc.obj, file = SignacObjFile)

} # end dummy if
  
if (Downstream_Analysis == TRUE) {

  ##=================
  ## dump the TSS positions for all genes
  ##=================
  DefaultAssay(sc.obj) <- "peaks"

  anno1 <- GetTSSPositions(annotation)
  export.bed(anno1, tempAnnotationFile)
  ## read the TSS positions
  ## 1st column: chromosome, 3rd column: position
  TSSDF <- read.table(tempAnnotationFile, header=F, sep="\t", stringsAsFactors=F)

  ##===============
  ## get the complete set of peaks
  ##===============
  DefaultAssay(sc.obj) <- "ATAC"
  peaks_clust <- import(targetpeakfile, format="bed")

  ## get the link information
  DefaultAssay(sc.obj) <- "peaks"
  linkinfo <- Links(sc.obj)

  ## export the linkinfo by constructing a custom data frame
  linkinfo_file_2 <- paste0(Peak_OutDir, '/Complete_Link_Info_2.bed')
  linkinfo_DF <- data.frame(peak=linkinfo$peak, gene=linkinfo$gene, score=linkinfo$score, zscore=linkinfo$zscore, pvalue=linkinfo$pvalue)
  write.table(linkinfo_DF, linkinfo_file_2, row.names=F, col.names=T, sep="\t", quote=F, append=F)

} # end if
 
##===================
## plot genes linked to peaks
## region near genes
##===================
if (Downstream_Analysis == TRUE) {

  cat(sprintf("\n\n Coverage plot for sample genes \n\n"))  

  DefaultAssay(sc.obj) <- "peaks"

  # currently plot a few genes
  GenesLinkedToPeaks <- unique(as.vector(linkinfo$gene))
  plotile <- paste0(Peak_OutDir, "/Plot_Genes_linked_to_peaks.pdf")
  pdf(plotile, width=10, height=8)
  for (i in 1:length(GenesLinkedToPeaks)) {
    currgene <- GenesLinkedToPeaks[i]
    cat(sprintf("\n Coverage plot for the gene : %s ", currgene))  
    ## get the TSS of this gene
    idx <- which(anno1$gene_name == currgene)
    if (length(idx) == 0) {
      next
    }
    currchr <- TSSDF[idx, 1]
    currTSS <- as.integer(TSSDF[idx, 3])
    cat(sprintf("  currchr : %s currTSS : %s ", currchr, currTSS))
    ## scan "linkinfo_DF" to get all the peaks linked to this gene
    idx <- which(linkinfo_DF$gene == currgene)  
    if (length(idx) == 0) {
      cat(sprintf("\n No peaks linked to this gene "))
      if (currgene %in% rownames(sc.obj@assays$SCT@counts)) {
        print(CoveragePlot(object = sc.obj, region = currgene, features = currgene, expression.assay = "SCT", extend.upstream = 5000, extend.downstream = 5000))
      } else {
        print(CoveragePlot(object = sc.obj, region = currgene, extend.upstream = 5000, extend.downstream = 5000))
      }
    } else {
      peakregions <- as.vector(linkinfo_DF[idx, 1])
      peakregions <- gsub("-", "_", peakregions)
      for (j in 1:length(peakregions)) {      
        if (j==1) {
          currchr <- unlist(strsplit(peakregions[j], "_"))[1]
          startpos <- as.integer(unlist(strsplit(peakregions[j], "_"))[2])
          endpos <- as.integer(unlist(strsplit(peakregions[j], "_"))[3])
        } else {
          startpos <- min(startpos, as.integer(unlist(strsplit(peakregions[j], "_"))[2]))
          endpos <- max(endpos, as.integer(unlist(strsplit(peakregions[j], "_"))[3]))
        }
      }
      currregion <- paste0(currchr, "-", (min(startpos, currTSS) - 1000), "-", (max(endpos, currTSS) + 1000))
      cat(sprintf("\n number of peaks linked : %s \n peaks linked to this gene : %s plotting currregion : %s ", length(idx), paste(peakregions, collapse="  "), currregion))
      if (currgene %in% rownames(sc.obj@assays$SCT@counts)) {
        print(CoveragePlot(object = sc.obj, region = currregion, features = currgene, expression.assay = "SCT"))
      } else { 
        print(CoveragePlot(object = sc.obj, region = currregion))
      }
    }
  }
  dev.off()

} # end dummy if

##=================
## plot target genes (if specified)
## region near genes
##=================
if (Downstream_Analysis == TRUE) {

  if (length(args) > 3) {
    cat(sprintf("\n\n Coverage plot for target genes \n\n"))  
    DefaultAssay(sc.obj) <- "peaks"
    plotile <- paste0(Peak_OutDir, "/Plot_TargetGenes_linked_to_peaks.pdf")
    pdf(plotile, width=10, height=8)
    for (i in 1:length(GeneList)) {
      currgene <- GeneList[i]  
      cat(sprintf("\n Coverage plot for the gene : %s ", currgene))  
      ## get the TSS of this gene
      idx <- which(anno1$gene_name == currgene)
      if (length(idx) == 0) {
        next
      }
      currchr <- TSSDF[idx, 1]
      currTSS <- as.integer(TSSDF[idx, 3])
      cat(sprintf("  currchr : %s currTSS : %s ", currchr, currTSS))
      ## scan "linkinfo_DF" to get all the peaks linked to this gene
      idx <- which(linkinfo_DF$gene == currgene)  
      if (length(idx) == 0) {
        cat(sprintf("\n No peaks linked to this gene "))
        if (currgene %in% rownames(sc.obj@assays$SCT@counts)) {
          print(CoveragePlot(object = sc.obj, region = currgene, features = currgene, expression.assay = "SCT", extend.upstream = 5000, extend.downstream = 5000))
        } else {
          print(CoveragePlot(object = sc.obj, region = currgene, extend.upstream = 5000, extend.downstream = 5000))
        }
      } else {
        peakregions <- as.vector(linkinfo_DF[idx, 1])
        peakregions <- gsub("-", "_", peakregions)
        for (j in 1:length(peakregions)) {
          if (j==1) {
            currchr <- unlist(strsplit(peakregions[j], "_"))[1]
            startpos <- as.integer(unlist(strsplit(peakregions[j], "_"))[2])
            endpos <- as.integer(unlist(strsplit(peakregions[j], "_"))[3])
          } else {
            startpos <- min(startpos, as.integer(unlist(strsplit(peakregions[j], "_"))[2]))
            endpos <- max(endpos, as.integer(unlist(strsplit(peakregions[j], "_"))[3]))
          }
        }
        currregion <- paste0(currchr, "-", (min(startpos, currTSS) - 1000), "-", (max(endpos, currTSS) + 1000))
        cat(sprintf("\n number of peaks linked : %s \n peaks linked to this gene : %s plotting currregion : %s ", length(idx), paste(peakregions, collapse="  "), currregion))
        if (currgene %in% rownames(sc.obj@assays$SCT@counts)) {
          print(CoveragePlot(object = sc.obj, region = currregion, features = currgene, expression.assay = "SCT"))
        } else {
          print(CoveragePlot(object = sc.obj, region = currregion))
        }
      }      
    }
    dev.off()
  }
} # end if


##===================
## plot peaks linked to genes
##===================
if (Downstream_Analysis == TRUE) {

  DefaultAssay(sc.obj) <- "peaks"

  ## currently plot a few peaks
  PeaksLinkedToGenes <- unique(as.vector(linkinfo$peak))
  plotile <- paste0(Peak_OutDir, "/Plot_Peaks_linked_to_genes.pdf")
  if (1) {

  pdf(plotile, width=10, height=8)
  for (i in 1:length(PeaksLinkedToGenes)) {
    ## we plot current peaks
    ## and also surrounding 5 Kb regions 
    ranges.show <- StringToGRanges(PeaksLinkedToGenes[i])
    currregion <- PeaksLinkedToGenes[i]
    print(CoveragePlot(object = sc.obj, region = currregion, ranges = ranges.show, extend.upstream = 5000, extend.downstream = 5000))
  }
  dev.off()

} # end dummy if

##=======
## then plot the gene activity scores
## currently plot for a few genes
##=======
if (Downstream_Analysis == TRUE) {

  DefaultAssay(sc.obj) <- 'GeneAct'
  GenesLinkedToPeaks <- unique(as.vector(linkinfo$gene))

  if (1) {
  
    plotfile_geneAct <- paste0(Cluster_OutDir, "/Featureplot_GeneActivityScore_Genes_Linked_to_Peaks.pdf")
    pdf(plotfile_geneAct, width=6, height=6)
    for (i in 1:length(GenesLinkedToPeaks)) {
      currgene <- GenesLinkedToPeaks[i]
      cat(sprintf("\n Featureplot_GeneActivityScore for the linked gene : %s ", currgene))
      ## we plot those genes and surrounding 5 Kb regions
      if (currgene %in% rownames(sc.obj@assays$GeneAct@counts)) {
        print(FeaturePlot(object = sc.obj, features = currgene, pt.size = 0.1, max.cutoff = 'q95', ncol = 1))
      }
    }
    dev.off()

  } # end dummy if

  if (length(args) > 3) {

    plotfile_geneAct <- paste0(Cluster_OutDir, "/Featureplot_GeneActivityScore_SampleGenes.pdf")
    pdf(plotfile_geneAct, width=6, height=6)
    for (i in 1:length(GeneList)) {
      currgene <- GeneList[i]
      ## we plot those genes and surrounding 5 Kb regions
      if (currgene %in% rownames(sc.obj@assays$GeneAct@counts)) {
        cat(sprintf("\n Featureplot_GeneActivityScore for the linked gene : %s ", currgene))
        print(FeaturePlot(object = sc.obj, features = currgene, pt.size = 0.1, max.cutoff = 'q95', ncol = 1))
      }
    }
    dev.off()

  } # end dummy if

} # end if

##===================
## Feature plot of gene expression for a group of genes
## using the RNA assay - for absolute counts
## or using the SCT assay - for normalized counts
##===================
if (Downstream_Analysis == TRUE) {

  DefaultAssay(sc.obj) <- "RNA"

  if (1) {

    plotfile_geneExpr <- paste0(Cluster_OutDir, "/Featureplot_GeneExpr_Genes_Linked_to_Peaks.pdf")
    pdf(plotfile_geneExpr, width=6, height=6)
    for (i in 1:length(GenesLinkedToPeaks)) {
      currgene <- GenesLinkedToPeaks[i]
      cat(sprintf("\n Featureplot_GeneExpr for the linked gene : %s ", currgene))
      ## we plot those genes and surrounding 5 Kb regions
      if (currgene %in% rownames(sc.obj@assays$RNA@counts)) {
        print(FeaturePlot(object = sc.obj, features = currgene, pt.size = 0.1, max.cutoff = 'q95', ncol = 1))
      }
    }
    dev.off()

  } # end dummy if

  if (length(args) > 3) {

    plotfile_geneExpr <- paste0(Cluster_OutDir, "/Featureplot_GeneExpr_SampleGenes.pdf")
    pdf(plotfile_geneExpr, width=6, height=6)
    for (i in 1:length(GeneList)) {
      currgene <- GeneList[i]
      ## we plot those genes and surrounding 5 Kb regions
      if (currgene %in% rownames(sc.obj@assays$RNA@counts)) {
        cat(sprintf("\n Featureplot_GeneExpr for the linked gene : %s ", currgene))
        print(FeaturePlot(object = sc.obj, features = currgene, pt.size = 0.1, max.cutoff = 'q95', ncol = 1))
      }
    }
    dev.off()

  } # end dummy if

} # end if

##===================
## get the peaks enriched for individual clusters
##===================
if (Downstream_Analysis == TRUE) {

  DefaultAssay(sc.obj) <- "peaks"

  ClustList <- sort(unique(sc.obj$clustID))

  ## this directory will store the marker peaks
  system(paste("mkdir -p", paste0(Peak_OutDir, "/Clusterwise_Marker_Peaks")))

  for (i in 1:length(ClustList)) {
     currclust <- ClustList[i]
     cat(sprintf("\n\n Finding peaks enriched in cluster %s ", currclust))
     da_peaks_DF <- FindMarkers(object = sc.obj, ident.1 = currclust, only.pos = TRUE, test.use = 'LR', min.pct = 0.05, latent.vars = 'nCount_peaks')
     markerpeakfile <- paste0(Peak_OutDir, "/Clusterwise_Marker_Peaks/Marker_peaks_cluster_", currclust, ".txt")
     ## write the row names as well
     ## row names contain peak information
     write.table(da_peaks_DF, markerpeakfile, row.names=T, col.names=T, sep="\t", quote=F, append=F) 
  }

} # end dummy if

##===================
## get the peaks enriched for individual pairs of clusters
##===================
if (Downstream_Analysis == TRUE) {

  DefaultAssay(sc.obj) <- "peaks"

  ClustList <- sort(unique(sc.obj$clustID))

  for (i in 1:length(ClustList)) {
    for (j in 1:length(ClustList)) {
      if (i == j) {
        next
      }
      currclust_1 <- ClustList[i]
      currclust_2 <- ClustList[j]
      ## this directory will store the marker peaks
      system(paste("mkdir -p", paste0(Peak_OutDir, "/PairwiseCluster_Marker_Peaks/cluster_", currclust_1)))
      cat(sprintf("\n\n Finding differential peaks between clusters %s and %s ", currclust_1, currclust_2))
      da_peaks_DF <- FindMarkers(object = sc.obj, ident.1 = currclust_1, ident.2 = currclust_2, only.pos = TRUE, test.use = 'LR', min.pct = 0.05, latent.vars = 'nCount_peaks')
      markerpeakfile <- paste0(Peak_OutDir, "/PairwiseCluster_Marker_Peaks/cluster_", currclust_1, "/Marker_peaks_between_cluster_", currclust_1, "_and_", currclust_2, ".txt")
      ## write the row names as well
      ## row names contain peak information
      write.table(da_peaks_DF, markerpeakfile, row.names=T, col.names=T, sep="\t", quote=F, append=F)
    }
  }

} # end dummy if

##==================
## Finding overrepresented motifs
## with respect to the peaks enriched in specific clusters
##==================
if (Downstream_Analysis == TRUE) {

  DefaultAssay(sc.obj) <- "peaks"

  ClustList <- sort(unique(sc.obj$clustID))

  ## this directory will store the marker peaks
  system(paste("mkdir -p", paste0(Motif_OutDir, "/Motifs_Clusterwise_Marker_Peaks")))

  for (i in 1:length(ClustList)) {
    currclust <- ClustList[i]
    cat(sprintf("\n\n Finding motifs corresponding to the peaks enriched in cluster %s ", currclust))
    da_peaks <- read.table(paste0(Peak_OutDir, "/Clusterwise_Marker_Peaks/Marker_peaks_cluster_", currclust, ".txt"), header=T, sep="\t", stringsAsFactors=F, row.names=1)
    if (nrow(da_peaks) == 0) {
      next
    }

    # get top differentially accessible peaks
    top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])
    if (length(top.da.peak) == 0) {
      next
    }

    ## there is an optional routine
    ## choosing a set of background peaks
    ## we do not consider for the moment

    ## test enrichment of motifs with respect to this set of peaks
    motifoutfile <- paste0(Motif_OutDir, '/Motifs_Clusterwise_Marker_Peaks/Motifs_enriched_peaks_cluster_', currclust, '.txt')    
    enriched.motifs <- FindMotifs(object = sc.obj, features = top.da.peak)   
    write.table(enriched.motifs, motifoutfile, row.names=F, col.names=T, sep="\t", quote=F, append=F) 
  
    ## plot the position weight matrices for the motifs, 
    ## so we can visualize the different motif sequences.
    ## we plot the top 10 motifs
    enriched.motifs <- enriched.motifs[order(enriched.motifs$pvalue), ]
    # currplot <- MotifPlot(object = sc.obj, motifs = head(rownames(enriched.motifs)))
    # ggplot2::ggsave(paste0(Motif_OutDir, "/Plot_Motifs_enriched_peaks_cluster_", currclust, ".pdf"), plot=currplot, width=10, height=8)
    pdf(paste0(Motif_OutDir, "/Motifs_Clusterwise_Marker_Peaks/Plot_Motifs_enriched_peaks_cluster_", currclust, ".pdf"), width=10, height=8)
    if (nrow(enriched.motifs) >= 12) {
      print(MotifPlot(object = sc.obj, motifs = rownames(enriched.motifs)[1:6]))
      print(MotifPlot(object = sc.obj, motifs = rownames(enriched.motifs)[7:12]))
    } else if (nrow(enriched.motifs) >= 6) {
      print(MotifPlot(object = sc.obj, motifs = rownames(enriched.motifs)[1:6]))
    } else {
      print(MotifPlot(object = sc.obj, motifs = rownames(enriched.motifs)))
    }
    dev.off()
  }

} # end dummy if

##==================
## Finding overrepresented motifs
## with respect to the differential peaks enriched in specific pair of clusters
##==================
if (Downstream_Analysis == TRUE) {

  DefaultAssay(sc.obj) <- "peaks"

  ClustList <- sort(unique(sc.obj$clustID))
  for (i in 1:length(ClustList)) {
    for (j in 1:length(ClustList)) {
      if (i == j) {
        next
      }    
      currclust_1 <- ClustList[i]
      currclust_2 <- ClustList[j]

      ## this directory will store the marker peaks
      system(paste("mkdir -p", paste0(Motif_OutDir, "/Motifs_PairwiseCluster_Marker_Peaks/cluster_", currclust_1)))

      cat(sprintf("\n\n Finding motifs corresponding to differential peaks between clusters %s and %s ", currclust_1, currclust_2))
      da_peaks <- read.table(paste0(Peak_OutDir, "/PairwiseCluster_Marker_Peaks/cluster_", currclust_1, "/Marker_peaks_between_cluster_", currclust_1, "_and_", currclust_2, ".txt"), header=T, sep="\t", stringsAsFactors=F, row.names=1)
      if (nrow(da_peaks) == 0) {
        next
      }
      
      # get top differentially accessible peaks
      top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])
      if (length(top.da.peak) == 0) {
        next
      }
      
      ## test enrichment of motifs with respect to this set of peaks
      motifoutfile <- paste0(Motif_OutDir, "/Motifs_PairwiseCluster_Marker_Peaks/cluster_", currclust_1, "/Motifs_enriched_peaks_between_cluster_", currclust_1, "_and_", currclust_2, '.txt')
      enriched.motifs <- FindMotifs(object = sc.obj, features = top.da.peak)
      write.table(enriched.motifs, motifoutfile, row.names=F, col.names=T, sep="\t", quote=F, append=F)
      
      ## plot the position weight matrices for the motifs, 
      ## so we can visualize the different motif sequences.
      enriched.motifs <- enriched.motifs[order(enriched.motifs$pvalue), ]
      # currplot <- MotifPlot(object = sc.obj, motifs = head(rownames(enriched.motifs)))
      # ggplot2::ggsave(paste0(Motif_OutDir, "/Plot_Motifs_enriched_peaks_between_cluster_", currclust_1, "_and_", currclust_2, ".pdf"), plot=currplot, width=10, height=8)

      pdf(paste0(Motif_OutDir, "/Motifs_PairwiseCluster_Marker_Peaks/cluster_", currclust_1, "/Plot_Motifs_enriched_peaks_between_cluster_", currclust_1, "_and_", currclust_2, ".pdf"), width=10, height=8)
      if (nrow(enriched.motifs) >= 12) {
        print(MotifPlot(object = sc.obj, motifs = rownames(enriched.motifs)[1:6]))
        print(MotifPlot(object = sc.obj, motifs = rownames(enriched.motifs)[7:12]))
      } else if (nrow(enriched.motifs) >= 6) {
        print(MotifPlot(object = sc.obj, motifs = rownames(enriched.motifs)[1:6]))
      } else {
        print(MotifPlot(object = sc.obj, motifs = rownames(enriched.motifs)))
      }
      dev.off()

    }
  }  

} # end dummy if


