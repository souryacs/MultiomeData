#!/bin/bash

##===================
## Running Velocyto on Single Cell data (output from cellranger)
## Output: .loom file (input of multivelo)
##===================

## path to GTF file for the reference genome (assuming mm10)
## user needs to edit
GTFFile='/path/to/CellRanger_ARC/refdata-cellranger-arc-mm10-2020-A-2.0.0/genes/genes.gtf'

## Directory containing cellranger-arc output from the multi-omic data
## user needs to edit them
CellRangerDir='/path/to/cellranger/output'

## output directory for Velocyto
## user needs to edit them
OutDir='velocyto_Out'
mkdir -p $OutDir

## bam file - gene expression - output from cellramger-arc
bamfile=$CellRangerDir'/outs/gex_possorted_bam.bam'

## sorted bam file - gene expression
sortedfile=$CellRangerDir'/outs/cellsorted_gex_possorted_bam.bam'

## we first manually sort the file using samtools
## https://github.com/velocyto-team/velocyto.py/issues/212
sammtools sort -@ 8 -t CB -O BAM -o $sortedfile $bamfile

## unzip the cellranger output barcode file
gunzip $CellRangerDir'/outs/filtered_feature_bc_matrix/barcodes.tsv.gz'

## run velocyto
velocyto run \
    -o $OutDir \
    -b $CellRangerDir'/outs/filtered_feature_bc_matrix/barcodes.tsv' \
    -m $basedir'/Ref/mm10_rmsk.gtf' \
    $bamfile \
    $GTFFile

## re-zip the cellranger output barcode file
gzip $CellRangerDir'/outs/filtered_feature_bc_matrix/barcodes.tsv'



 

