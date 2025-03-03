library(Seurat)
library(SeuratObject)
library(optparse)
library(stringr)

option_list <- list(
    make_option("--input_rda_ge", help="Input path of the rda file"),
    make_option("--output_cluster_ctrl", help="Output path of the text file"),
    make_option("--output_cluster_ofint", help="Output path of the text file"),
    make_option("--cluster_ctrl",help="cluster of control to subset"),
    make_option("--cluster_ofint",help="cluster of interest to subset"))

parser <- OptionParser(usage="Rscript %prog [options]", description = " ", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)

input_rda_ge <- args$options$input_rda_ge
output_cluster_ctrl <- args$options$output_cluster_ctrl
output_cluster_ofint <- args$options$output_cluster_ofint
cluster_ctrl <- args$options$cluster_ctrl
cluster_ofint <- args$options$cluster_ofint

print("loading RDA object")
load(input_rda_ge)

cluster_ctrl <- as.integer(unlist(str_split(cluster_ctrl,",")))
cluster_ofint <- as.integer(unlist(str_split(cluster_ofint,",")))

sub_sobj <- subset(x = sobj, idents = cluster_ctrl)
bc_vec <- unlist(str_split(rownames(sub_sobj@meta.data),"-"))
bc_vec <- bc_vec[bc_vec != "1"]
write.table(bc_vec,file=output_cluster_ctrl,quote = F,col.names = F, row.names = F)

sub_sobj <- subset(x = sobj, idents = cluster_ofint)
bc_vec <- unlist(str_split(rownames(sub_sobj@meta.data),"-"))
bc_vec <- bc_vec[bc_vec != "1"]
write.table(bc_vec,file=output_cluster_ofint,quote = F,col.names = F, row.names = F)