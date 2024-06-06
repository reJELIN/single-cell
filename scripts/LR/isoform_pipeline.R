library(Isosceles)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(scuttle)
library(DEXSeq)
library(RColorBrewer)
library(pheatmap)
library(dittoSeq)
library(Nebulosa)
library(biomaRt)
library(rtracklayer)
library(purrr)
library(ggtranscript)
library(dplyr)
library(magrittr)
library(Seurat)
library(SeuratObject)
library(waffle)
library(gridExtra)
library(grid)
library(cowplot)
library(gridGraphics)
#library(ggarrange)
library(ggpubr)
library(cowplot)
library(RColorBrewer)
library(reshape2)
library(optparse)

option_list <- list(
    make_option("--input_dir_ge", help="Input path of the alignment"),
    make_option("--output_dir_ge", help="Output path")
    make_option("--gtf_file",help="gtf file for the isoform analysis"),
    make_option("--sample_id"),
    make_option("--isoform_transcript_matrix"),
    make_option("--seurat_obj_rda")
)

parser <- OptionParser(usage="Rscript %prog [options]", description = " ", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)

input_dir_ge <- args$options$input_dir_ge
output_dir_ge <- args$option$output_dir_ge
gtf_file <- args$option$gtf_file
sample_id <- args$option$sample_id
isoform_transcript_matrix <- args$option$isoform_transcript_matrix
seurat_obj_rda <- args$option$seurat_obj_rda
transcript_counts <- readRDS()

transcript_counts <- as.data.frame(transcript_counts)

rownames(transcript_counts) <- transcript_counts[,1]
transcript_counts <- transcript_counts[,-1]

load()


isoform_name <- rownames(transcript_counts)
barcodes <- colnames(transcript_counts)

#get barcodes from the seurat object and select them in the transcript matrix
barcodes_sobj <- rownames(sobj@meta.data)
transcript_counts <- transcript_counts[,barcodes_sobj]

#import gtf file in R
my_obj <- import(gtf_file)

gtf <- my_obj %>% dplyr::as_tibble()
gtf$transcript_id <- paste0(gtf$transcript_id,".",gtf$transcript_version)

transcript_id <- my_obj@elementMetadata@listData$transcript_id
transcript_version <- my_obj@elementMetadata@listData$transcript_version
intact_transcript <- paste0(transcript_id,".",transcript_version)
my_obj@elementMetadata@listData$transcript_id <- intact_transcript

df_gtf <- data.frame(transcript_id=transcript_id,
                     transcript_version=transcript_version,
                     transcript_ensembl=intact_transcript) 
#get ride of NA first line
df_gtf <- df_gtf[-1,]
df_gtf_no_duplicated <- df_gtf %>% distinct()

#length(intersect(isoform_name,df_gtf_no_duplicated$transcript_id))

#keep only isoform in nanopore matrice
isoform_df <- df_gtf_no_duplicated[df_gtf_no_duplicated$transcript_id %in% isoform_name,]
rownames(isoform_df) <- isoform_df$transcript_id


isoform_name_df <- data.frame(transcript_id=isoform_name)

merge_isoform_df <- merge(isoform_name_df,isoform_df, by="transcript_id")
