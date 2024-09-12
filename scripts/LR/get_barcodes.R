library(Seurat)
library(SeuratObject)
library(optparse)
library(stringr)

option_list <- list(
    make_option("--input_rda_ge", help="Input path of the rda file"),
    make_option("--output_barcodes_file_txt_ge", help="Output path of the text file")
)

parser <- OptionParser(usage="Rscript %prog [options]", description = " ", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)

input_rda_ge <- args$options$input_rda_ge
output_barcodes_file_txt_ge <- args$options$output_barcodes_file_txt_ge

print("loading RDA object")
load(input_rda_ge)

barcodes_vector=vapply(stringr::str_split(rownames(sobj@meta.data),"-"), `[`, 1, FUN.VALUE=character(1))
print("Writing barcodes table")
write.table(barcodes_vector,file=output_barcodes_file_txt_ge,quote = F,col.names = F, row.names = F)