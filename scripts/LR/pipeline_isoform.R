suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(tibble))
suppressMessages(library(ggplot2))
suppressMessages(library(scuttle))
suppressMessages(library(DEXSeq))
suppressMessages(library(RColorBrewer))
suppressMessages(library(dittoSeq))
suppressMessages(library(Nebulosa))
suppressMessages(library(biomaRt))
suppressMessages(library(rtracklayer))
suppressMessages(library(purrr))
suppressMessages(library(ggtranscript))
suppressMessages(library(dplyr))
suppressMessages(library(magrittr))
suppressMessages(library(Seurat))
suppressMessages(library(SeuratObject))
suppressMessages(library(waffle))
suppressMessages(library(gridExtra))
suppressMessages(library(grid))
suppressMessages(library(cowplot))
suppressMessages(library(gridGraphics))
suppressMessages(library(ggpubr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(reshape2))
suppressMessages(library(optparse))
suppressMessages(library(AnnotationDbi))
suppressMessages(library(presto))
suppressMessages(library(httr))

option_list <- list(
    make_option("--input_rda_ge", help="Input path of the rda file"),
    make_option("--output_rda_ge", help="Output path of the rda file"),
    make_option("--output_dir_ge", help="Output path"),
    make_option("--gtf_file",help="gtf file for the isoform analysis"),
    make_option("--path_isoform_transcript_matrix",help="path to isoform transcript matrix"),
    make_option("--species",help="name of the species")
)

parser <- OptionParser(usage="Rscript %prog [options]", description = "", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)

input_rda_ge <- args$options$input_rda_ge
output_rda_ge <- args$options$output_rda_ge
output_dir_ge <- args$option$output_dir_ge
gtf_file <- args$option$gtf_file
species_ensembl <- args$option$species
isoform_transcript_matrix <- args$option$path_isoform_transcript_matrix

if(species_ensembl == "human"){ensemble2gene <- read.table(file="/mnt/beegfs/pipelines/single-cell/lr_1.3_test/single-cell/common/database/ensembl_hsapiens.txt",header=T,sep="\t")}

if(species_ensembl == "mouse"){ensemble2gene <- read.table(file="/mnt/beegfs/pipelines/single-cell/lr_1.3_test/single-cell/common/database/ensembl_mmusculus.txt",header=T,sep="\t")}

quite_message <- function(message_txt){
    message(message_txt)
    quit("no", -1)
}

load(input_rda_ge)

scmat <- Seurat::Read10X(data.dir = isoform_transcript_matrix)

scmat <- scmat[,rownames(sobj@meta.data)]
isoform_name <- scmat@Dimnames[[1]]
barcodes <- scmat@Dimnames[[2]]

#import gtf file in R
my_obj <- import(gtf_file)

#add concat column transcript_id and transcript_version
#used for plotting
gtf <- my_obj %>% dplyr::as_tibble()
gtf$transcript_id <- paste0(gtf$transcript_id,".",gtf$transcript_version)

#add transcription version and id in the gtf file object
transcript_id <- my_obj@elementMetadata@listData$transcript_id
transcript_version <- my_obj@elementMetadata@listData$transcript_version
transcript_name <- my_obj@elementMetadata@listData$transcript_name
gene_id <- my_obj@elementMetadata@listData$gene_id
gene_name <- my_obj@elementMetadata@listData$gene_name
#concat version and id in order to have ENST isoform
intact_transcript <- paste0(transcript_id,".",transcript_version)
my_obj@elementMetadata@listData$transcript_id <- intact_transcript

#dataframe with 3 columns of transcript info
df_gtf <- data.frame(transcript_id=transcript_id,
                     transcript_version=transcript_version,
                     transcript_ensembl=intact_transcript,
                     transcript_name = transcript_name,
                     gene_id = gene_id,
                     gene_name = gene_name)
                      
#get ride of NA first line
df_gtf <- df_gtf[-1,]
df_gtf_no_duplicated <- df_gtf %>% distinct()
df_gtf_no_duplicated <- na.omit(df_gtf_no_duplicated) 

#length(intersect(isoform_name,df_gtf_no_duplicated$transcript_id))

#keep only isoform in nanopore matrice
isoform_df <- df_gtf_no_duplicated[df_gtf_no_duplicated$transcript_id %in% isoform_name,]
rownames(isoform_df) <- isoform_df$transcript_id

isoform_name_df <- data.frame(transcript_id=isoform_name)
merge_isoform_df <- merge(isoform_name_df,isoform_df, by="transcript_id")
rownames(scmat) <- merge_isoform_df$transcript_name

sobj[['Isoform']] = CreateAssayObject(counts = scmat)

sobj <- Seurat::NormalizeData(object=sobj,normalization.method="LogNormalize",scale.factor = 10000,assay="Isoform")

features_to_scale=sobj@assays$Isoform@counts@Dimnames[[1]]

sobj <- ScaleData(object=sobj,assay="Isoform",vars.to.regress=c(),features = features_to_scale)

###wilcoauc
group_contrasts <- utils::combn(levels(as.factor(sobj@meta.data$seurat_clusters)), 2, simplify = FALSE)
marker_results_list <- lapply(group_contrasts,function(group_contrast){wilcoxauc(X=sobj,group_by="seurat_clusters",groups_use=group_contrast,seurat_assay = 'Isoform',assay="data")}) #TODO test with 10% pct
names(marker_results_list) <- sapply(group_contrasts, paste0,collapse = "__")

#setconfig(config(sslverifypeer = 0L)) or httr::set_config(httr::config(ssl_verifypeer = FALSE))

#set_config(config(ssl_verifypeer = 0L))
#get gene_symbol for each ENST ensembl id version
#mart <- useMart("ensembl","hsapiens_gene_ensembl")
#'www', 'useast', 'asia'.

#ensembl <- useEnsembl(biomart = "ensembl", 
#                   dataset = "hsapiens_gene_ensembl", 
#                   mirror = "useast")
                   
marker_results_list_translate <- lapply(names(marker_results_list),function(contrast_name){
	
	marker_df <- marker_results_list[[contrast_name]] %>% as.data.frame()
	marker_df$transcript_name <- marker_df$feature
	marker_df$transcript_name <- as.character(marker_df$transcript_name)
	
	#ensemble2gene <- getBM(attributes=c("external_gene_name","ensembl_gene_id","ensembl_transcript_id_version"),
    #                   filters = "ensembl_transcript_id_version",
    #                   values = marker_df$feature, 
    #                   mart = ensembl)
                       
    #ensemble2gene$ensembl_transcript_id_version <- as.character(ensemble2gene$ensembl_transcript_id_version)   
     
    #marker_df <- left_join(marker_df,ensemble2gene,by="ensembl_transcript_id_version",relationship = "many-to-many")
    marker_df <- na.omit(marker_df)
    marker_df$contrast <- contrast_name
       
    return(marker_df)})

#bind each list of markers together
marker_df <- do.call(rbind, marker_results_list_translate)
marker_df <- merge(marker_df,df_gtf_no_duplicated, by="transcript_name")
#THIS IS NOT ISOFORM SWITCH distribution, we find DEG between each clusters and keep only DEG genes that have more than one transcripts version
isoswitch_df <- marker_df %>%
    #dplyr::filter(.data$pct_in > 25 | .data$pct_out > 25) %>%
	dplyr::filter(.data$padj <= 0.05) %>%
	dplyr::filter(.data$logFC > 0.1) %>%
	dplyr::group_by(.data$gene_id, .data$contrast) %>%
	dplyr::mutate(n_groups = length(unique(.data$group)),n_transcripts = length(unique(.data$transcript_name))) %>%
	dplyr::ungroup() %>%
	dplyr::filter(.data$n_groups > 1, .data$n_transcripts > 1) %>%
	dplyr::select(-"n_groups", -"n_transcripts") %>%
	dplyr::arrange(.data$contrast, .data$gene_id) %>%
	as.data.frame()
	

path_ftpl <- paste0(output_dir_ge,"isoform/FeaturePlots")
path_report <- paste0(output_dir_ge,"isoform/reports")
dir.create(file.path(path_report), recursive = TRUE)
dir.create(file.path(path_ftpl), recursive = TRUE)

write.table(isoswitch_df,file=paste0(output_dir_ge,"/isoform/",sobj@misc$params$sample.name.GE,"_",sobj@misc$params$clustering$ident,"_isoform_findmarkers_all.txt",sep=","))

Idents(sobj) <- sobj@meta.data$seurat_clusters

cluster_umap <- DimPlot(sobj,reduction=sobj@misc$params$clustering$umap,dims=c(1,2),group.by='seurat_clusters',pt.size=2,label=F) + labs(y='UMAP 2',x='UMAP 1') + ggtitle('') + theme(axis.text = (element_text(size=15,colour='black')),
	     axis.title.x = element_text(size=15),
	     axis.title.y = element_text(size=15),
	     legend.position= "right",
	     legend.text = element_text(size=12),
	     legend.title= element_text(size=12))

for(i_gene in unique(isoswitch_df$gene_name)){

	tmp_df_isoform <- isoswitch_df[isoswitch_df$gene_name == i_gene,]
	
	data <- sobj@assays$Isoform@data[unique(tmp_df_isoform$feature),] %>% t() %>%as.data.frame()

	colnames(data) <- colnames(data)[1:length(colnames(data))]

	data$cluster_labels <- sobj@meta.data$seurat_clusters

	isoform_columns <- grep(i_gene, names(data), value=TRUE)
	
	tmp_df_mean <- setNames(data.frame(matrix(ncol = length(data), nrow = 0)), c(grep(i_gene,names(data), value=TRUE),"cluster_labels"))
	
	for(i in unique(data$cluster_labels)){
	 tmp_row <- colMeans(Filter(is.numeric, data[data$cluster_labels == i,]))
	 tmp_row <- append(tmp_row,i)
	 
	 names(tmp_row) <- c(grep(i_gene,names(data), value=TRUE),"cluster_labels")
	 
	 tmp_row_df <- as.data.frame(t(tmp_row))
	 colnames(tmp_row_df) <- names(tmp_row)
	 
	 tmp_df_mean = rbind(tmp_df_mean,tmp_row_df)
	}
	
	#p_grid_table <- grid.arrange(gridExtra::tableGrob(tmp_df_mean))
	

	longdata = tmp_df_mean %>% 
	  pivot_longer(-cluster_labels)

	p1 <- ggplot(longdata, aes(fill=name, y=as.numeric(value), x=cluster_labels)) + 
	    geom_bar(position="stack", stat="identity") + labs(y="Average logCounts") + theme_classic() + labs(y="Average logcounts", x="clusters", fill="Transcript name") +theme_classic() +
	     theme(axis.text = (element_text(size=15,colour='black')),
	     axis.title.x = element_text(size=15),
	     axis.title.y = element_text(size=15),
	     legend.position= "right",
	     legend.text = element_text(size=12),
	     legend.title= element_text(size=12))


	#create empty dataframe to append row
	tmp_df_prop <- setNames(data.frame(matrix(ncol = length(data), nrow = 0)), c(grep(i_gene,names(data), value=TRUE),"cluster_labels"))

	for(i in unique(data$cluster_labels)){

	 tmp_row <- colSums(Filter(is.numeric, data[data$cluster_labels == i,]))*100/sum(colSums(Filter(is.numeric, data[data$cluster_labels == i,])))
	 tmp_row <- append(tmp_row,i)
	 
	 names(tmp_row) <- c(grep(i_gene,names(data), value=TRUE),"cluster_labels")
	 
	 tmp_row_df <- as.data.frame(t(tmp_row))
	 colnames(tmp_row_df) <- names(tmp_row)
	 
	 tmp_df_prop = rbind(tmp_df_prop,tmp_row_df)
	}

	longdata_prop = tmp_df_prop %>% 
	  pivot_longer(-cluster_labels)

	p2 <- ggplot(longdata_prop, aes(fill = name, 
		                  y = as.numeric(value), 
		                  x = cluster_labels)) + 
	  geom_bar(position="stack",stat="identity")  + labs(y="Percentage (%)", x="clusters", fill="Transcript name") +theme_classic() + theme(axis.text = (element_text(size=15)),
	 			axis.title.x = element_text(size=15,colour="black"),
	 			axis.title.y = element_text(size=15),
	 			legend.position= "right",
	 			legend.text = element_text(size=12),
	 			legend.title= element_text(size=12)) 
		   
	    

	p3 <- DotPlot(object=sobj,assay="Isoform",
	  features=unique(tmp_df_isoform$feature),
	  cols = c("blue", "red"),
	  group.by = "seurat_clusters") + labs(y="clusters", x="transcript name") +theme_classic() + theme(axis.text = (element_text(size=15,angle = 50, vjust = 0.5, hjust=1)),
	 			axis.title.x = element_text(size=15),
	 			axis.title.y = element_text(size=15),
	 			legend.position= "right",
	 			legend.text = element_text(size=12),
	 			legend.title= element_text(size=12)) 

	
	p5  <- FeaturePlot(object=sobj, features = unique(tmp_df_isoform$feature),reduction = sobj@misc$params$clustering$umap,cols = c("grey", "red"))
	
	if(length(unique(tmp_df_isoform$feature)) <= 4){
	pdf(paste0(path_ftpl,"/",i_gene,".pdf"),width=7*length(unique(tmp_df_isoform$feature))-1.5,height=5*length(unique(tmp_df_isoform$feature))-1.5)
	plot(p5)
	dev.off()
	}
	
	if(length(unique(tmp_df_isoform$feature)) > 4){
	pdf(paste0(path_ftpl,"/",i_gene,".pdf"),width=5*length(unique(tmp_df_isoform$feature))-1.5,height=3*length(unique(tmp_df_isoform$feature))-1.5)
	plot(p5)
	dev.off()
	}
	
	if(length(unique(tmp_df_isoform$feature)) > 8){
	pdf(paste0(path_ftpl,"/",i_gene,".pdf"),width=2*length(unique(tmp_df_isoform$feature))-1.5,height=1.5*length(unique(tmp_df_isoform$feature))-1.5)
	plot(p5)
	dev.off()
	}
	
	
	if(length(ls(pattern="violin_to_plot_")) > 0){
	   rm(list=ls(pattern="violin_to_plot_"))
	}
	
	p_violin  <- Seurat::VlnPlot(sobj,assay="Isoform", features = c(unique(tmp_df_isoform$feature)),group.by="seurat_clusters",combine=FALSE)
	
	
	caller <- function(i){p_violin[[i]]}
	
	for(i in seq(1,length(p_violin))){
	if(i == length(p_violin)){
		p_violin[[i]] <- p_violin[[i]] + theme(axis.text = (element_text(size=15)),
		 			axis.title.x = element_text(size=15),
		 			axis.title.y = element_text(size=15),
		 			legend.position= "right",
		 			legend.text = element_text(size=12),
		 			legend.title= element_text(size=12),
		 			plot.title=element_text(size=10))}else{
		p_violin[[i]] <- p_violin[[i]] + theme(axis.text = (element_text(size=15)),
					legend.position= "none",
		 			axis.title.x = element_text(size=15),
		 			axis.title.y = element_text(size=15),
		 			plot.title=element_text(size=10))}
		 			}
	for(i in seq(1,length(p_violin))){
	assign(paste0("violin_to_plot_",i),p_violin[[i]])
	}
	
	
	p_violin_all <- patchwork::wrap_plots(lapply(ls(pattern="violin_to_plot_"),get)) + patchwork::plot_layout(ncol = 3)
 	
	#color palette and transcript biotype name for the fill of ggtranscript plot
	palette_hex <- c("#EE3B3B","#76EEC6","#EE1289","#B23AEE","#1C86EE","#FFF68F","#EEA2AD","#B3EE3A","#2E8B57","#EEEE00","#EE7600")
    transcript_biotype_col <-c("protein_coding","nonsense_mediated_decay","pseudogene","lncRNA","retained_intron","transcribed_unprocessed_pseudogene","non_stop_decay","processed_pseudogene","misc_RNA","TEC","transcribed_unitary_pseudogene")


	vec_gene <- unique(isoswitch_df[isoswitch_df$gene_name == i_gene,]$feature)
	
	
	i_gene_annotation_from_gtf <- gtf %>% dplyr::filter(!is.na(gene_name), gene_name == i_gene)
	 
	 #if biomart gene symbol not found is looking based on ENST transcript id
	 if(dim(i_gene_annotation_from_gtf)[1] == 0){i_gene_annotation_from_gtf <- gtf %>% dplyr::filter(transcript_name %in% vec_gene)}
	  
	  
	  
	i_gene_annotation_from_gtf <- i_gene_annotation_from_gtf %>% 
	  dplyr::select(
	    seqnames,
	    start,
	    end,
	    strand,
	    type,
	    gene_name,
	    transcript_name,
	    transcript_type,
	    transcript_id
	  )
	  
	#i_gene_annotation_from_gtf$transcript_id_gene_name <- paste0(i_gene_annotation_from_gtf$transcript_id,".",i_gene_annotation_from_gtf$gene_name)
	  
	i_exons <- i_gene_annotation_from_gtf %>% dplyr::filter(type == "exon") %>% dplyr::filter(transcript_name %in% vec_gene)
	    
	i_gene_rescaled <- shorten_gaps(
	  i_exons, 
	  to_intron(i_exons, "transcript_name"), 
	  group_var = "transcript_name"
	  )
	
	
	transcript_data_exons <- dplyr::filter(i_gene_rescaled,type == "exon")
	transcript_data_intron <- dplyr::filter(i_gene_rescaled,type == "intron")
	
	p_ggtranscript <-  ggplot(transcript_data_exons,aes(xstart = start,xend = end,y = transcript_name)) +
	 			geom_range(aes(fill = transcript_type)) +
	 			geom_intron(data = transcript_data_intron,
	 			arrow.min.intron.length = 200) +
	 			theme_classic() + theme(axis.text = (element_text(size=15)),
	 			axis.title.x = element_text(size=12),
	 			axis.title.y = element_text(size=12),
	 			axis.text.x = element_text(angle = 50, vjust = 0.5, hjust=1,size=12),
	 			axis.text.y = element_text(angle = 50, vjust = 0.5, hjust=1,size=12),
	 			legend.position= "top",
	 			) + scale_fill_manual(breaks=transcript_biotype_col,values=palette_hex) + ylab("transcript version") + xlab("length") + labs(fill="Transcript type")
	

   	ggarrange_plot <- ggarrange(p1,p2,p3,p_ggtranscript,p_violin_all,cluster_umap,ncol = 3, nrow = 2)
   	
   	nb_transcript_for_plot=unique(tmp_df_isoform$feature)
   	

   	if(length(nb_transcript_for_plot) <= 3){
        pdf(paste0(path_report,"/",i_gene,".pdf"),width=7*3,height=3.5*3)
    	plot(ggarrange_plot)
    	dev.off()
   	}
   	
   	if(length(nb_transcript_for_plot) >= 4){
   	    pdf(paste0(path_report,"/",i_gene,".pdf"),width=7*4,height=3*4)
    	plot(ggarrange_plot)
    	dev.off()
   	}
   	
   	if(length(nb_transcript_for_plot) >= 8){
   	    pdf(paste0(path_report,"/",i_gene,".pdf"),width=1.5*length(nb_transcript_for_plot),height=1.5*length(nb_transcript_for_plot))
        plot(ggarrange_plot)
        dev.off()
   	}

}

sobj@assays[["Isoform"]]@scale.data <- matrix(nrow = 0, ncol = 0)

save(file=output_rda_ge,sobj,compress="bzip2")