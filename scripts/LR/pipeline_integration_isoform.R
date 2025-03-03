suppressMessages(library(optparse))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(magrittr))
suppressMessages(library(Seurat))
suppressMessages(library(SeuratObject))
suppressMessages(library(rtracklayer))
suppressMessages(library(presto))
suppressMessages(library(stringr))
suppressMessages(library(ggtranscript))
suppressMessages(library(ggpubr))

option_list <- list(
    make_option("--input_int_rda", help="Input path of the rda file"),
    make_option("--output_rda_int", help="Output path of the rda file"),
    make_option("--output_dir", help="Output path global"),
    make_option("--gtf_file",help="gtf file for the isoform analysis")
    )
    
#make_option("--species",help="name of the species")

parser <- OptionParser(usage="Rscript %prog [options]", description = " ", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)

input_int_rda <- args$options$input_int_rda
output_rda_int <- args$options$output_rda_int
output_dir_int <- args$options$output_dir
gtf_file <- args$options$gtf_file

load(input_int_rda)

sample_name_list <- unique(sobj@meta.data$orig.ident)
seurat_read_iso <- function(x){
    return(Seurat::Read10X(paste0(output_dir_int,x,'/',x,'/transcript_raw_feature_bc_matrix/')))
}

#load isoform matrix of each sample
iso_X_list <- lapply(sample_name_list, seurat_read_iso)
names(iso_X_list) <- sample_name_list

#get commun isoform between samples
common_iso <- Reduce(intersect, lapply(iso_X_list, function(x){rownames(x)}))
iso_X_list <- lapply(iso_X_list, function(x){x[common_iso,]})

# filter barcordes only present in the seurat object in the isoform matrix
keep_barcodes <- function(x,sobj,iso_mtx_list){
	mtx <- iso_mtx_list[[x]]
	bc <-rownames(sobj@meta.data[sobj@meta.data$orig.ident == x,])
	colnames(mtx) <- paste0(x,"_",colnames(mtx))
	return(mtx[,bc])
}

iso_X_list <- lapply(sample_name_list,keep_barcodes,sobj=sobj,iso_mtx_list=iso_X_list)
names(iso_X_list) <- sample_name_list

# binding each isoform matrix together
scmat_iso <- Reduce(cbind,iso_X_list)

#get isoforme name and barcodes
isoform_name <- scmat_iso@Dimnames[[1]]
barcodes <- scmat_iso@Dimnames[[2]]

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

rownames(scmat_iso) <- merge_isoform_df$transcript_name

sobj[['Isoform']] = CreateAssayObject(counts = scmat_iso)

sobj <- Seurat::NormalizeData(object=sobj,normalization.method="LogNormalize",scale.factor = 10000,assay="Isoform")
 
features_to_scale=sobj@assays$Isoform@counts@Dimnames[[1]]
 
sobj <- ScaleData(object=sobj,assay="Isoform",vars.to.regress=c(),features = features_to_scale)

###wilcoauc
group_contrasts <- utils::combn(levels(as.factor(sobj@meta.data$seurat_clusters)), 2, simplify = FALSE)
marker_results_list <- lapply(group_contrasts,function(group_contrast){wilcoxauc(X=sobj,group_by="seurat_clusters",groups_use=group_contrast,seurat_assay = 'Isoform',assay="data")}) #TODO test with 10% pct
names(marker_results_list) <- sapply(group_contrasts, paste0,collapse = "__")
                   
marker_results_list_translate <- lapply(names(marker_results_list),function(contrast_name){
	
	marker_df <- marker_results_list[[contrast_name]] %>% as.data.frame()
	marker_df$transcript_name <- marker_df$feature
	marker_df$transcript_name <- as.character(marker_df$transcript_name)
    marker_df <- na.omit(marker_df)
    marker_df$contrast <- contrast_name
       
    return(marker_df)})

#bind each list of markers together
marker_df <- do.call(rbind, marker_results_list_translate)
marker_df <- merge(marker_df,df_gtf_no_duplicated, by="transcript_name")
#THIS IS NOT ISOFORM SWITCH distribution, we find DEG between each clusters and keep only DEG genes that have more than one transcripts version
isoswitch_df <- marker_df %>%
	dplyr::filter(.data$padj <= 0.05) %>%
	dplyr::filter(.data$logFC > 0.1) %>%
	dplyr::group_by(.data$gene_id, .data$contrast) %>%
	dplyr::mutate(n_groups = length(unique(.data$group)),n_transcripts = length(unique(.data$transcript_name))) %>%
	dplyr::ungroup() %>%
	dplyr::filter(.data$n_groups > 1, .data$n_transcripts > 1) %>%
	dplyr::select(-"n_groups", -"n_transcripts") %>%
	dplyr::arrange(.data$contrast, .data$gene_id) %>%
	as.data.frame()
	


ouput_dir_dims<-strsplit(output_rda_int,split="isoform")[[1]][1]

path_ftpl <- paste0(ouput_dir_dims,"isoform/FeaturePlots")
path_report <- paste0(ouput_dir_dims,"isoform/reports")
dir.create(file.path(path_report), recursive = TRUE)
dir.create(file.path(path_ftpl), recursive = TRUE)

if("name.int" %in% names(sobj@misc$params)){misc_sample_name=sobj@misc$params$name.int}else{misc_sample_name=sobj@misc$params$name.grp}

write.table(isoswitch_df,file=paste0(ouput_dir_dims,"/isoform/",misc_sample_name,"_",sobj@misc$params$clustering$ident,"_isoform_findmarkers_all.txt"),sep=",")

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

save(file=output_rda_int,sobj,compress="bzip2")