#!/share/apps/R/3.4.3/bin/Rscript

library(Signac)
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
SignacObjFile <- args[1]
OutDir <- args[2]

system(paste("mkdir -p", OutDir))

sc.obj <- readRDS(SignacObjFile)
cat(sprintf("\n\n Number of cells in the original Signac object : %s ", ncol(sc.obj)))

# extract neighborhood graph
nn_idx <- sc.obj@neighbors$weighted.nn@nn.idx
nn_dist <- sc.obj@neighbors$weighted.nn@nn.dist
nn_cells <- sc.obj@neighbors$weighted.nn@cell.names

# save neighborhood graph
write.table(nn_idx, paste0(OutDir, "/nn_idx.txt"), sep = ',', row.names = F, col.names = F, quote = F)
write.table(nn_dist, paste0(OutDir, "/nn_dist.txt"), sep = ',', row.names = F, col.names = F, quote = F)
write.table(nn_cells, paste0(OutDir, "/nn_cells.txt"), sep = ',', row.names = F, col.names = F, quote = F)

## user needs to edit 
## we assume that the clustering is performed with resolution = 0.3
## so we use the fields "SCT_snn_res.0.3", "peaks_snn_res.0.3" and "wsnn_res.0.3" and fetch corresponding information
## for clustering with different resolutions, use the corresponding resolution value and field names

## also write the cell-wise cluster labels 

## RNA-seq cluster
outdf <- data.frame(barcode=colnames(sc.obj), celltype=sc.obj@meta.data$SCT_snn_res.0.3)
write.table(outdf, paste0(OutDir, "/RNA_Cluster_Annotation.tsv"), row.names=F, col.names=T, sep="\t", quote=F, append=F) 

## ATAC-seq cluster
outdf <- data.frame(barcode=colnames(sc.obj), celltype=sc.obj@meta.data$peaks_snn_res.0.3)
write.table(outdf, paste0(OutDir, "/ATAC_Cluster_Annotation.tsv"), row.names=F, col.names=T, sep="\t", quote=F, append=F)

## multimodal cluster
outdf <- data.frame(barcode=colnames(sc.obj), celltype=sc.obj@meta.data$wsnn_res.0.3)
write.table(outdf, paste0(OutDir, "/Multimodal_Cluster_Annotation.tsv"), row.names=F, col.names=T, sep="\t", quote=F, append=F)
