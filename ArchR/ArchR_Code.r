library(ArchR)
library(ggplot2)
library(rtracklayer)
addArchRThreads(threads=16)
## reference genome
addArchRGenome("hg38")

args <- commandArgs(trailingOnly = TRUE)

## path containing the scATAC-seq data from cellranger
ATACPath <- args[1]
scRNASeuratObjFile <- args[2]
CurrOutDir <- args[3]

system(paste("mkdir -p", CurrOutDir))

## output log file
x <- Sys.time()
x <- gsub(" PST", "", x)
x <- gsub("-| |:", "_", x)
textfile <- paste0(CurrOutDir, '/out_', x, '.log')
sink(textfile)

MarkerGeneOutDir <- paste0(CurrOutDir, '/MarkerGenes')
system(paste("mkdir -p", MarkerGeneOutDir))

MarkerPeakOutDir <- paste0(CurrOutDir, '/MarkerPeaks')
system(paste("mkdir -p", MarkerPeakOutDir))

## batch effect correction parameter
batch_correction <- FALSE

##***********************
## read the scRNA-seq data object first
##***********************
seRNA <- readRDS(scRNASeuratObjFile)
cat(sprintf("\n\n *** Number of cells in Seurat RNA seq object : %s ", length(colnames(seRNA))))
## tabulate the number of cells present in each cluster of seRNA
write.table(as.data.frame(table(seRNA@meta.data$seurat_clusters)), paste0(CurrOutDir, '/Cluster_Info_scRNASeq_Seurat.txt'), row.names=F, col.names=T, sep="\t", quote=F, append=F)

##***********************
## process the scATAC-seq data
##***********************
if (file.exists(paste0(CurrOutDir, '/Save-ArchR-Project.rds')) == FALSE) {

	##==============
	## read the input files
	##==============
	inputFiles <- getInputFiles(paths = ATACPath)
	cat(sprintf("\n Read input files : %s ", paste(inputFiles, collapse=" ")))

	##==============
	## create arrow files
	##==============
	## the option addGeneScoreMat = TRUE automatically calculates the Gene scores 
	## (using the chromatin accessibility data) and add them in the arrow files
	ArrowFiles <- createArrowFiles(inputFiles = inputFiles, sampleNames = c('pbmc_granulocyte'), filterTSS = 4, filterFrags = 1000, addTileMat = TRUE, addGeneScoreMat = TRUE)

	##==============
	## inferring doublets
	##==============
	## k Refers to how many cells near a "pseudo-doublet" to count.
	## knnMethod Refers to the embedding to use for nearest neighbor search.
	doubScores <- addDoubletScores(input = ArrowFiles, k = 10, knnMethod = "UMAP", LSIMethod = 1)

	##==============
	## create ArchR project
	##==============
	## copyArrows = TRUE is recommened so that you maintain an unaltered copy for later usage.
	proj <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = CurrOutDir, copyArrows = FALSE)

	cat(sprintf("\n Memory Size = %s MB ", round(object.size(proj) / 10^6, 3)))
	cat(sprintf("\n Available data matrices : %s ", paste(getAvailableMatrices(proj), collapse=" ")))
	cat(sprintf("\n Number of cells of the ArchR project : %s ", nCells(input = proj)))

	## Plotting QC metrics - log10(Unique Fragments) vs TSS enrichment score
	## not required - since this QC plot is done by ArchR itself
	if (0) {
		df <- getCellColData(proj, select = c("log10(nFrags)", "TSSEnrichment"))
		p <- ggPoint(x = df[,1], y = df[,2], colorDensity = TRUE, continuousSet = "sambaNight", xlabel = "Log10 Unique Fragments", ylabel = "TSS Enrichment", xlim = c(log10(500), quantile(df[,1], probs = 0.99)), ylim = c(0, quantile(df[,2], probs = 0.99))) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
		ggsave(filename = paste0(currdir, '/QC_nFrag_vs_TSSEnrichment.pdf'), plot = p)
	}

	##==============
	## filtering doublets
	##==============
	proj <- filterDoublets(proj)
	cat(sprintf("\n\n **** After filtering the doublets - number of cells of the ArchR project : %s ", nCells(input = proj)))

	##==============
	## dimensionality reduction
	##==============
	## creates a "reducedDims" object "IterativeLSI"
	proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI", iterations = 2, clusterParams = list(resolution = c(0.2), sampleCells = 10000, n.start = 10), varFeatures = 25000, dimsToUse = 1:30)

	##==============
	## batch effect correction with harmony
	##==============
	if (batch_correction == TRUE) {
		## creates a new "reducedDims" object called "Harmony"
		proj <- addHarmony(ArchRProj = proj, reducedDims = "IterativeLSI", name = "Harmony", groupBy = "Sample")
	}

	##==============
	## clustering using Seurat's findcluster function
	##==============
	proj <- addClusters(input=proj, reducedDims="IterativeLSI", method="Seurat", name="Clusters", resolution=0.8)

	##==============
	## run UMAP
	##==============
	proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI", name = "UMAP", nNeighbors = 30, minDist = 0.5, metric = "cosine")

	## plot UMAP
	p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
	p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
	plotPDF(p1, p2, name = 'Plot_UMAP_Sample_Clusters_reducedDims_IterativeLSI.pdf', ArchRProj = proj, addDOC = FALSE, width = 8, height = 6)

	if (batch_correction == TRUE) {
		## UMAP embedding for the reduced dimension generated by Harmony
		proj <- addUMAP(ArchRProj = proj, reducedDims = "Harmony", name = "UMAPHarmony", nNeighbors = 30, minDist = 0.5, metric = "cosine")
		p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
		p4 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
		plotPDF(p3, p4, name = 'Plot_UMAP_Sample_Clusters_reducedDims_Harmony.pdf', ArchRProj = proj, addDOC = FALSE, width = 8, height = 6)

		## plot UMAP showing all of UMAP plots
		plotPDF(p1, p2, p3, p4, name = 'Plot_UMAP_Sample_Clusters_BOTH_UMAP_Harmony.pdf', ArchRProj = proj, addDOC = FALSE, width = 10, height = 8)
	}	

	## tabulate the number of cells present in each cluster
	write.table(as.data.frame(table(proj$Clusters)), paste0(CurrOutDir, '/Cluster_Info.txt'), row.names=F, col.names=T, sep="\t", quote=F, append=F)

	##==============
	## Identify marker genes based on gene scores, 
	##==============
	## get the set of clusters
	ClusterSet <- sort(unique(proj$Clusters))
	cat(sprintf("\n==>> set of clusters : %s ", paste(ClusterSet, collapse=" ")))

	## call the getMarkerFeatures() function with useMatrix = "GeneScoreMatrix". 
	## know the cluster-specific features with groupBy = "Clusters" 
	## This function returns a SummarizedExperiment object 
	## containing relevant information on the marker features identified. 
	markersGS <- getMarkerFeatures(ArchRProj = proj, useMatrix = "GeneScoreMatrix", groupBy = "Clusters", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")

	## get a list of DataFrame objects, one for each of our clusters, 
	## containing the relevant marker features using the getMarkers() function:
	markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
	# head(markerList$C1)	## specifying cluster 1
	# head(markerList$C2)	## specifying cluster 2
	for (i in 1:length(ClusterSet)) {
		write.table(as.data.frame(markerList[[ClusterSet[i]]]), paste0(MarkerGeneOutDir, '/marker_genes_', ClusterSet[i], '.txt'), row.names=F, col.names=T, sep="\t", quote=F, append=F)
	}

	## construct the list of marker genes
	## select the top 10 marker genes from each cluster
	markerGenes <- c()
	for (i in 1:length(ClusterSet)) {
		df <- as.data.frame(markerList[[ClusterSet[i]]])
		if (nrow(df) > 10) {
			markerGenes <- c(markerGenes, as.vector(df$name[1:10]))
		} else {
			markerGenes <- c(markerGenes, as.vector(df$name))
		}
	}
	cat(sprintf("\n Number of marker genes : %s ", length(markerGenes)))
	cat(sprintf("\n List of marker genes : %s ", paste(markerGenes, collapse=" ")))

	## To visualize all of the marker features simultaneously, 
	## we can create a heatmap using the markerHeatmap() function, 
	heatmapGS <- markerHeatmap(seMarker = markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25", labelMarkers = markerGenes,transpose = TRUE)

	## save the heatmap
	plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

	## Visualizing Marker Genes on an Embedding
	p <- plotEmbedding(ArchRProj = proj, colorBy = "GeneScoreMatrix", name = markerGenes, embedding = "UMAP",  imputeWeights = NULL, plotAs = "points")	#quantCut = c(0.01, 0.95),

	## Important: To plot a specific gene, we can subset this plot list: like p$CD14
	plotPDF(plotList = p, name = "Plot-UMAP-Marker-Genes-WO-Imputation", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

	##==============
	## Marker Genes Imputation with MAGIC
	##==============
	proj <- addImputeWeights(proj)
	p <- plotEmbedding(ArchRProj = proj, colorBy = "GeneScoreMatrix", name = markerGenes, embedding = "UMAP", imputeWeights = getImputeWeights(proj), plotAs = "points")
	## Important: To plot a specific gene, we can subset this plot list: like p$CD14
	plotPDF(plotList = p, name = "Plot-UMAP-Marker-Genes-W-Imputation", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

	##==============
	## Track Plotting with ArchRBrowser
	##==============
	p <- plotBrowserTrack(ArchRProj = proj, groupBy = "Clusters", geneSymbol = markerGenes, upstream = 50000, downstream = 50000)

	## Important: To plot a track of a specific gene, we can select one from the list.
	## like grid::grid.newpage(); grid::grid.draw(p$CD14)
	plotPDF(plotList = p, name = "Plot-Tracks-Marker-Genes", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

	## save ArchR project
	saveArchRProject(ArchRProj = proj, outputDirectory = CurrOutDir, load = FALSE)	

	##===========
	## launch our local interactive genome browser
	##===========
	## we use the ArchRBrowser() function.
	## currently commented
	# ArchRBrowser(proj)

	##==============
	## unconstrained integration
	##==============
	## parameters "useMatrix", "matrixName" and "reducedDims" have default values
	## parameter "seRNA" contains the input seRNA object
	## parameter "groupRNA" indicates the field in Seurat metadata containing the cluster information
	## parameters "nameCell", "nameGroup" and "nameScore" are field names - "_Un" suffix corresponds to  unconstrained integration
	proj <- addGeneIntegrationMatrix(ArchRProj = proj, useMatrix = "GeneScoreMatrix", matrixName = "GeneIntegrationMatrix", reducedDims = "IterativeLSI", seRNA = seRNA, addToArrow = FALSE, groupRNA = "seurat_clusters", nameCell = "predictedCell_Un", nameGroup = "predictedGroup_Un", nameScore = "predictedScore_Un")

	## match the scATAC-seq clusters with the clusters from the scRNA-seq
	## shows which scRNA-seq cell type is most abundant in each of the scATAC-seq clusters.
	cM <- as.matrix(confusionMatrix(proj$Clusters, proj$predictedGroup_Un))
	write.table(cM, paste0(CurrOutDir, '/scATAC_scRNA_Cluster_Confusion_Matrix.txt'), row.names=T, col.names=T, sep="\t", quote=F, append=F)
	preClust <- colnames(cM)[apply(cM, 1 , which.max)]
	preClust_cM <- cbind(preClust, rownames(cM))
	write.table(preClust_cM, paste0(CurrOutDir, '/scATAC_scRNA_Cluster_Match.txt'), row.names=F, col.names=T, sep="\t", quote=F, append=F)

	##==============
	## constrained integration
	##==============
	## further operations are possible only when we know the cell annotations (or a subset)
	## to be implemented




	##==============
	## Comparing Unconstrained and Constrained Integrations
	##==============
	## create a color palette 
	pal <- paletteDiscrete(values = seRNA@meta.data$seurat_clusters)
	## visualize the integration by overlaying the scRNA-seq cell types on our scATAC-seq data based on the unconstrained integration
	p1 <- plotEmbedding(proj, colorBy = "cellColData", name = "predictedGroup_Un", pal = pal)

	## visualize the integration by overlaying the scRNA-seq cell types on our scATAC-seq data based on the constrained integration
	if (0) {
		## possible when we have results for constrained integration
		## to be implemented
		p2 <- plotEmbedding(proj, colorBy = "cellColData", name = "predictedGroup_Co", pal = pal)
	}

	if (0) {
		## possible when we have results for constrained integration
		## to be implemented
		plotPDF(p1,p2, name = "Plot-UMAP-RNA-Integration.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
	}
	plotPDF(p1, name = "Plot-UMAP-RNA-Integration.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

	##************************
	## adding Pseudo-scRNA-seq profiles for each scATAC-seq cell
	##************************
	## possible when we have results for constrained integration


	##*************************
	## Labeling scATAC-seq clusters with scRNA-seq information
	##*************************
	## possible when we have results for constrained integration


	##*************************
	## Making Pseudo-bulk Replicates
	## need to call peaks
	## it'll be preferable to use scRNA-seq defined clusters
	## here we are just using the clusters defined by the ATAC-seq data
	##*************************
	## parallel processing is returning error
	## set threads = 1
	proj <- addGroupCoverages(ArchRProj = proj, groupBy = "Clusters", threads = 1)

	## save ArchR project
	saveArchRProject(ArchRProj = proj, outputDirectory = CurrOutDir, load = FALSE)

	##*************************
	## Calling Peaks 
	## using the clusters defined by the ATAC-seq data
	##*************************
	## Calling Peaks w/ Macs2
	pathToMacs2 <- findMacs2()
	proj <- addReproduciblePeakSet(ArchRProj = proj, groupBy = "Clusters", pathToMacs2 = pathToMacs2)
	## peak set as a GRanges object
	peakset_MACS2 <- getPeakSet(proj)
	## use export.bed function from "rtracklayer" library
	## to export this set of peaks as a data frame
	export.bed(peakset_MACS2, paste0(CurrOutDir, '/Peaks_MACS2.txt'))

	## Calling Peaks w/ TileMatrix
	proj2 <- addReproduciblePeakSet(ArchRProj = proj, groupBy = "Clusters", peakMethod = "Tiles", method = "p")
	peakset_TILE <- getPeakSet(proj2)
	export.bed(peakset_TILE, paste0(CurrOutDir, '/Peaks_Tile.txt'))
	cat(sprintf("\n\n ** fraction of overlapping peaks between MACS2 and tile method : %s ", length(subsetByOverlaps(getPeakSet(proj2), getPeakSet(proj))) / length(getPeakSet(proj2))))

	## add peak matrix
	proj <- addPeakMatrix(proj)
	getAvailableMatrices(proj)

	##*************************
	## Identifying Marker Peaks with ArchR
	##*************************

	## use the addMarkerFeatures() function in combination with useMatrix = "PeakMatrix".
	## using the clusters defined by the ATAC-seq data
	markersPeaks <- getMarkerFeatures(ArchRProj = proj, useMatrix = "PeakMatrix", groupBy = "Clusters", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
	markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
	
	## dump the list of marker peaks
	## here we are using ATAC-seq clusters	
	ClusterSet <- sort(unique(proj$Clusters))
	for (i in 1:length(ClusterSet)) {
		write.table(as.data.frame(markerList[[ClusterSet[i]]]), paste0(MarkerPeakOutDir, '/marker_peaks_', ClusterSet[i], '.txt'), row.names=F, col.names=T, sep="\t", quote=F, append=F)
	}

	## visualize these marker peaks (or any features output by getMarkerFeatures()) 
	## as a heatmap using the markerHeatmap() function.
	heatmapPeaks <- markerHeatmap(seMarker = markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5", transpose = TRUE)
	plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

	## Marker Peak MA and Volcano Plots
	## for individual cell groups (clusters)
	## here we use ATACseq cluster C1
	pma <- markerPlot(seMarker = markersPeaks, name = "C1", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs="MA")
	pv <- markerPlot(seMarker = markersPeaks, name = "C1", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "Volcano")
	plotPDF(pma, pv, name = "ATAC_Cluster_C1-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

	## Marker Peaks in Browser Tracks
	## Here we specify plotting the GATA1 gene via the geneSymbol parameter.
	## we also specify using marker peaks for the cluster C1
	p <- plotBrowserTrack(ArchRProj = proj, groupBy = "Clusters", geneSymbol = c("GATA1"), features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)["C1"], upstream = 50000, downstream = 50000)
	plotPDF(p, name = "Plot-Tracks-With-Features", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

	##*************************
	## Pairwise Testing Between Groups
	##*************************
	## here we use ATACseq clusters (group by "Clusters") C1 and C2
	markerTest <- getMarkerFeatures(ArchRProj = proj, useMatrix = "PeakMatrix", groupBy = "Clusters", testMethod = "wilcoxon", bias = c("TSSEnrichment", "log10(nFrags)"), useGroups = "C1", bgdGroups = "C2")
	pma <- markerPlot(seMarker = markerTest, name = "C1", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
	pv <- markerPlot(seMarker = markerTest, name = "C1", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
	plotPDF(pma, pv, name = "C1-vs-C2-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

	##*******************
	## Motif Enrichment in Differential Peaks
	## like in the above example - differential peaks between classes C1 and C2
	##*******************

	## add motif annotations
	proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")

	## test differentially accessible peaks between classes C1 and C2 for enrichment of various motifs
	## motifs upregulated in class C1
	motifs_C1_C2_Up_C1 <- peakAnnoEnrichment(seMarker = markerTest, ArchRProj = proj, peakAnnotation = "Motif", cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
	## plot the motifs
	df <- data.frame(TF = rownames(motifs_C1_C2_Up_C1), mlog10Padj = assay(motifs_C1_C2_Up_C1)[,1])
	df <- df[order(df$mlog10Padj, decreasing = TRUE),]
	df$rank <- seq_len(nrow(df))
	ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + geom_point(size = 1) + ggrepel::geom_label_repel(data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), size = 1.5, nudge_x = 2, color = "black") + theme_ArchR() + ylab("-log10(FDR) Motif Enrichment") + xlab("Rank Sorted TFs Enriched") + scale_color_gradientn(colors = paletteContinuous(set = "comet"))

	## motifs downregulated in class C1 (i.e. upregulated in class C2)
	motifs_C1_C2_Up_C2 <- peakAnnoEnrichment(seMarker = markerTest, ArchRProj = proj, peakAnnotation = "Motif", cutOff = "FDR <= 0.1 & Log2FC <= -0.5")
	df <- data.frame(TF = rownames(motifs_C1_C2_Up_C2), mlog10Padj = assay(motifs_C1_C2_Up_C2)[,1])
	df <- df[order(df$mlog10Padj, decreasing = TRUE),]
	df$rank <- seq_len(nrow(df))
	ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + geom_point(size = 1) + ggrepel::geom_label_repel(data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), size = 1.5, nudge_x = 2, color = "black") + theme_ArchR() + ylab("-log10(FDR) Motif Enrichment") + xlab("Rank Sorted TFs Enriched") + scale_color_gradientn(colors = paletteContinuous(set = "comet"))

	plotPDF(ggUp, ggDo, name = "C1-vs-C2-Markers-Motifs-Enriched", width = 5, height = 5, proj = proj, addDOC = FALSE)

	##*******************
	## Motif Enrichment in Marker Peaks
	##*******************

	## perform motif enrichment on our marker peaks identified using getMarkerFeatures()
	enrichMotifs <- peakAnnoEnrichment(seMarker = markersPeaks, ArchRProj = proj, peakAnnotation = "Motif", cutOff = "FDR <= 0.1 & Log2FC >= 0.5")

	## plot these motif enrichments across all cell groups
	heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
	plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

	## enrichment of Encode TF Binding Sites
	proj <- addArchRAnnotations(ArchRProj = proj, collection = "EncodeTFBS")
	enrichEncode <- peakAnnoEnrichment(seMarker = markersPeaks, ArchRProj = proj, peakAnnotation = "EncodeTFBS", cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
	heatmapEncode <- plotEnrichHeatmap(enrichEncode, n = 7, transpose = TRUE)
	plotPDF(heatmapEncode, name = "EncodeTFBS-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

	## enrichment testing with respect to the curated peak calls from bulk ATAC-seq experiments
	## already included in the ArchR data
	proj <- addArchRAnnotations(ArchRProj = proj, collection = "ATAC")
	enrichATAC <- peakAnnoEnrichment(seMarker = markersPeaks, ArchRProj = proj, peakAnnotation = "ATAC", cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
	heatmapATAC <- plotEnrichHeatmap(enrichATAC, n = 7, transpose = TRUE)
	plotPDF(heatmapATAC, name = "ATAC-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

	## enrichment with respect to CODEX TFBSs
	proj <- addArchRAnnotations(ArchRProj = proj, collection = "Codex")
	enrichCodex <- peakAnnoEnrichment(seMarker = markersPeaks, ArchRProj = proj, peakAnnotation = "Codex", cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
	heatmapCodex <- plotEnrichHeatmap(enrichCodex, n = 7, transpose = TRUE)
	plotPDF(heatmapCodex, name = "Codex-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

	##=============
	## There is also provision for custom enrichment
	## with respect to a given set of peaks
	##=============
	EncodePeaks <- c(
	  Encode_K562_GATA1 = "https://www.encodeproject.org/files/ENCFF632NQI/@@download/ENCFF632NQI.bed.gz",
	  Encode_GM12878_CEBPB = "https://www.encodeproject.org/files/ENCFF761MGJ/@@download/ENCFF761MGJ.bed.gz",
	  Encode_K562_Ebf1 = "https://www.encodeproject.org/files/ENCFF868VSY/@@download/ENCFF868VSY.bed.gz",
	  Encode_K562_Pax5 = "https://www.encodeproject.org/files/ENCFF339KUO/@@download/ENCFF339KUO.bed.gz"
	)

	## add a custom annotation to our ArchRProject using the addPeakAnnotation() function.
	## we call our custom annotation “ChIP”.
	proj <- addPeakAnnotations(ArchRProj = proj, regions = EncodePeaks, name = "ChIP")
	## we use this custom annotation to perform the peak annotation enrichment
	enrichRegions <- peakAnnoEnrichment(seMarker = markersPeaks, ArchRProj = proj, peakAnnotation = "ChIP", cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
	heatmapRegions <- plotEnrichHeatmap(enrichRegions, n = 7, transpose = TRUE)
	plotPDF(heatmapRegions, name = "Regions-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

	##******************************
	## ChromVAR Deviatons Enrichment with ArchR
	## chromVAR is designed for predicting enrichment of TF activity on a per-cell basis from sparse chromatin accessibility data
	##******************************

	##===========
	## Motif Deviations
	##===========
	## First, lets make sure we have added motif annotations to our ArchRProject.
	if("Motif" %ni% names(proj@peakAnnotation)){
		proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")
	}
	## add a set of background peaks which are used in computing deviations
	## Background peaks are chosen using the chromVAR::getBackgroundPeaks() function 
	## which samples peaks based on similarity in GC-content 
	## and number of fragments across all samples using the Mahalanobis distance.
	proj <- addBgdPeaks(proj)
	## compute per-cell deviations accross all of our motif annotations
	proj <- addDeviationsMatrix(ArchRProj = proj, peakAnnotation = "Motif", force = TRUE)
	## access these deviations
	plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)
	plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

	## extract a subset of motifs for downstream analysis
	motifs <- c("GATA1", "CEBPA", "EBF1", "IRF4", "TBX21", "PAX5")
	markerMotifs <- getFeatures(proj, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
	markerMotifs

	## %ni% expression which is an ArchR helper function that provides the opposite of %in% from base R.
	markerMotifs <- grep("z:", markerMotifs, value = TRUE)
	markerMotifs <- markerMotifs[markerMotifs %ni% "z:SREBF1_22"]
	markerMotifs

	## plot the distribution of chromVAR deviation scores for each cluster.
	## here the "groupBy" argument uses Clusters; recommended to use "Clusters2", i.e. integrated cluster between scATAC-seq and scRNA-seq
	p <- plotGroups(ArchRProj = proj, groupBy = "Clusters", colorBy = "MotifMatrix", name = markerMotifs, imputeWeights = getImputeWeights(proj))
	plotPDF(p, name = "Plot-Groups-Deviations-w-Imputation", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

	## overlay the z-scores on our UMAP embedding
	p <- plotEmbedding(ArchRProj = proj, colorBy = "MotifMatrix", name = sort(markerMotifs), embedding = "UMAP", imputeWeights = getImputeWeights(proj), plotAs = "points")
	plotPDF(p, name = "Plot-Groups-Deviations-w-Imputation-zscore-UMAP-color-MotifMatrix", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

	## overlay the gene scores for each of these TFs on the UMAP embedding.
	markerRNA <- getFeatures(proj, select = paste(motifs, collapse="|"), useMatrix = "GeneScoreMatrix")
	markerRNA <- markerRNA[markerRNA %ni% c("SREBF1","CEBPA-DT")]
	markerRNA
	p <- plotEmbedding(ArchRProj = proj, colorBy = "GeneScoreMatrix", name = sort(markerRNA), embedding = "UMAP", imputeWeights = getImputeWeights(proj), plotAs = "points")
	plotPDF(p, name = "Plot-Groups-Deviations-w-Imputation-GeneScore-UMAP-color-GeneScoreMatrix", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

	## we previously linked our scATAC-seq data with corresponding scRNA-seq data,
	## plot the linked gene expression for each of these TFs on the UMAP embedding.
	markerRNA <- getFeatures(proj, select = paste(motifs, collapse="|"), useMatrix = "GeneIntegrationMatrix")
	markerRNA <- markerRNA[markerRNA %ni% c("SREBF1","CEBPA-DT")]
	markerRNA
	p <- plotEmbedding(ArchRProj = proj, colorBy = "GeneIntegrationMatrix", name = sort(markerRNA), embedding = "UMAP", continuousSet = "blueYellow", imputeWeights = getImputeWeights(proj), plotAs = "points")
	plotPDF(p, name = "Plot-Groups-Deviations-w-Imputation-GeneExprScore-UMAP-color-GeneIntegrationMatrix", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

	##===============
	## custom deviations - we can add 
	##===============

	## save ArchR project
	saveArchRProject(ArchRProj = proj, outputDirectory = CurrOutDir, load = FALSE)	

} else {
	proj <- loadArchRProject(path = CurrOutDir)

	## top 3 marker genes per cluster
	ClusterSet <- sort(unique(proj$Clusters))
	markerGenes <- c()
	for (i in 1:length(ClusterSet)) {
		df <- read.table(paste0(MarkerGeneOutDir, '/marker_genes_', ClusterSet[i], '.txt'), header=T, sep="\t", stringsAsFactors=F)
		markerGenes <- c(markerGenes, as.vector(df$name[1:3]))
	}
	cat(sprintf("\n List of marker genes (top 3 from each cluster) : %s ", paste(markerGenes, collapse=" ")))

	## get the marker peaks per cluster
	markersPeaks <- getMarkerFeatures(ArchRProj = proj, useMatrix = "PeakMatrix", groupBy = "Clusters", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
	markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
}


##=============
## dump motifs - later merge this code above
##=============
motif_markerpeak_outdir <- paste0(CurrOutDir, '/Annotations/Motifs_Marker_Peaks')
system(paste("mkdir -p", motif_markerpeak_outdir))
enrichMotifs <- peakAnnoEnrichment(seMarker = markersPeaks, ArchRProj = proj, peakAnnotation = "Motif", cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
aa <- assays(enrichMotifs)
AssayList <- names(aa)
MotifList <- rownames(aa[[AssayList[1]]])
for (i in 1:length(ClusterSet)) {
	cat(sprintf("\n\n ==>> Dump motifs for marker peaks - processing cluster : %s ", ClusterSet[i]))
	for (j in 1:length(AssayList)) {
		cat(sprintf("   - processing assay : %s ", AssayList[j]))
		if (j == 1) {
			finalDF <- data.frame(Motif=MotifList, x=aa[[AssayList[j]]][[ClusterSet[i]]])
			colnames(finalDF) <- c('Motif', AssayList[j])
		} else {
			currDF <- data.frame(x=aa[[AssayList[j]]][[ClusterSet[i]]])
			colnames(currDF) <- c(AssayList[j])
			finalDF <- cbind.data.frame(finalDF, currDF)
		}
	}
	write.table(finalDF, paste0(motif_markerpeak_outdir, '/Motif_marker_peak_', ClusterSet[i], '.txt'), row.names=F, col.names=T, sep="\t", quote=F, append=F)
}


##***********************
## Motif Footprinting
##***********************
## obtain the positions of the relevant motifs
motifPositions <- getPositions(proj)
## subset to a few TF motifs that we are interested in
motifs <- c("GATA1", "CEBPA", "EBF1", "IRF4", "TBX21", "PAX5")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]
markerMotifs

## use the group coverages
## it would be better if we use clusters jointly annotated by scRNA-seq and scATAC-seq
## here we are using clusters defined by the ATAC-seq

## already done before
# proj <- addGroupCoverages(ArchRProj = proj, groupBy = "Clusters")

seFoot <- getFootprints(ArchRProj = proj, positions = motifPositions[markerMotifs], groupBy = "Clusters")

##===============
## Normalization of Footprints for Tn5 Bias
##===============
## option 1: subtract Tn5 bias
plotFootprints(seFoot = seFoot, ArchRProj = proj, normMethod = "Subtract", plotName = "Footprints-Subtract-Bias", addDOC = FALSE, smoothWindow = 5)
## option 2: divide by Tn5 bias
plotFootprints(seFoot = seFoot, ArchRProj = proj, normMethod = "Divide", plotName = "Footprints-Divide-Bias", addDOC = FALSE, smoothWindow = 5)

##===============
## Footprinting Without Normalization for Tn5 Bias
##===============
plotFootprints(seFoot = seFoot, ArchRProj = proj, normMethod = "None", plotName = "Footprints-No-Normalization", addDOC = FALSE, smoothWindow = 5)

##===============
## Feature Footprinting
##===============

## already done before
# proj <- addGroupCoverages(ArchRProj = proj, groupBy = "Clusters")

## create TSS insertion profiles without normalization for Tn5 bias
seTSS <- getFootprints(ArchRProj = proj, positions = GRangesList(TSS = getTSS(proj)), groupBy = "Clusters", flank = 2000)

##  plot the TSS insertion profiles for each cell group
plotFootprints(seFoot = seTSS, ArchRProj = proj, normMethod = "None", plotName = "TSS-No-Normalization", addDOC = FALSE, flank = 2000, flankNorm = 100)

##***********************
##  co-accessibility in ArchR
##***********************
proj <- addCoAccessibility(ArchRProj = proj, reducedDims = "IterativeLSI")

## retrieve this co-accessibility information
## returns a DataFrame object if returnLoops = FALSE.
cA <- getCoAccessibility(ArchRProj = proj, corCutOff = 0.5, resolution = 1, returnLoops = FALSE)

## it would be better if we use clusters jointly annotated by scRNA-seq and scATAC-seq
## here we are using clusters defined by the ATAC-seq
p <- plotBrowserTrack(ArchRProj = proj, groupBy = "Clusters", geneSymbol = markerGenes, upstream = 50000, downstream = 50000, loops = getCoAccessibility(proj))
plotPDF(plotList = p, name = "Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

## save ArchR project
saveArchRProject(ArchRProj = proj, outputDirectory = CurrOutDir, load = FALSE)	

##***********************
##  peak-to-gene linkage in ArchR
##***********************
proj <- addPeak2GeneLinks(ArchRProj = proj, reducedDims = "IterativeLSI")
p2g <- getPeak2GeneLinks(ArchRProj = proj, corCutOff = 0.45, resolution = 1, returnLoops = FALSE)
p <- plotBrowserTrack(ArchRProj = proj, groupBy = "Clusters", geneSymbol = markerGenes, upstream = 50000, downstream = 50000, loops = getPeak2GeneLinks(proj))
plotPDF(plotList = p, name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
p <- plotPeak2GeneHeatmap(ArchRProj = proj, groupBy = "Clusters")

# ##***********************
# ##  Trajectory analysis in ArchR
# ##***********************

# ## needs scRNA-seq annotation



# ## save ArchR project
# saveArchRProject(ArchRProj = proj, outputDirectory = CurrOutDir, load = FALSE)	
