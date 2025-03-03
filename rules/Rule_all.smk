"""
##########################################################################
This function make the input of rule all
##########################################################################
"""
    
def get_targets(STEPS):
    targets = {}
    if "Alignment_countTable_GE" in STEPS:
        targets["Alignment_countTable_GE"]=[
        #multiqc
        expand(os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/QC_reads/{sample_name_ge}_RAW.html"), zip, sample_name_ge=ALIGN_SAMPLE_NAME_GE),
        #alignment
        expand(os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/KALLISTOBUS/run_info.json"), zip, sample_name_ge=ALIGN_SAMPLE_NAME_GE),
        #build_count_matrix
        expand(os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/KALLISTOBUS/{sample_name_ge}.mtx"), zip, sample_name_ge=ALIGN_SAMPLE_NAME_GE),
        expand(os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/KALLISTOBUS/{sample_name_ge}.barcodes.txt"), zip, sample_name_ge=ALIGN_SAMPLE_NAME_GE),
        expand(os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/KALLISTOBUS/{sample_name_ge}.genes.txt"), zip, sample_name_ge=ALIGN_SAMPLE_NAME_GE),
        expand(os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/KALLISTOBUS/Materials_and_Methods.txt"), zip, sample_name_ge=ALIGN_SAMPLE_NAME_GE)
        ]
    if "Alignment_countTable_ADT" in STEPS:
        targets["Alignment_countTable_ADT"]=[
        #multiqc
        expand(os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/QC_reads/{sample_name_adt}_RAW.html"), zip, sample_name_adt=ALIGN_SAMPLE_NAME_ADT),
        #alignment
        expand(os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/run_info.json"), zip, sample_name_adt=ALIGN_SAMPLE_NAME_ADT),
        #build_count_matrix
        expand(os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/{sample_name_adt}.mtx"), zip, sample_name_adt=ALIGN_SAMPLE_NAME_ADT),
        expand(os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/{sample_name_adt}.barcodes.txt"), zip, sample_name_adt=ALIGN_SAMPLE_NAME_ADT),
        expand(os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/{sample_name_adt}.genes.txt"), zip, sample_name_adt=ALIGN_SAMPLE_NAME_ADT),
        expand(os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/Materials_and_Methods.txt"), zip, sample_name_adt=ALIGN_SAMPLE_NAME_ADT)
        ]
    if "Alignment_annotations_TCR_BCR" in STEPS:
        targets["Alignment_annotations_TCR_BCR"]=[
        #multiqc
        expand(os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR + "/{sample_name_tcr_bcr}/QC_reads/{sample_name_tcr_bcr}_RAW.html"), zip, sample_name_tcr_bcr=ALIGN_SAMPLE_NAME_TCR_BCR),
        #alignment_annotations_
        expand(os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR + "/{sample_name_tcr_bcr}/{sample_name_tcr_bcr}_CellRanger/outs/filtered_contig_annotations.csv"), zip, sample_name_tcr_bcr=ALIGN_SAMPLE_NAME_TCR_BCR),
        expand(os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR + "/{sample_name_tcr_bcr}/{sample_name_tcr_bcr}_CellRanger/outs/web_summary.html"), zip, sample_name_tcr_bcr=ALIGN_SAMPLE_NAME_TCR_BCR),
        expand(os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR + "/{sample_name_tcr_bcr}/{sample_name_tcr_bcr}_CellRanger/outs/Materials_and_Methods.txt"), zip, sample_name_tcr_bcr=ALIGN_SAMPLE_NAME_TCR_BCR)
        ]
    if "Droplets_QC_GE" in STEPS:
        targets["Droplets_QC_GE"]=[
        expand(os.path.normpath("{outputqc_droplets_dir_ge}" + "/QC_droplets/" + "{sample_name_ge}_kneeplot.png"), zip, outputqc_droplets_dir_ge=QC_OUTPUT_DIR_GE, sample_name_ge=QC_SAMPLE_NAME_GE) if  str(QC_EMPTYDROPS_RETAIN) == "NULL" else expand(os.path.normpath("{outputqc_droplets_dir_ge}" + "/QC_droplets_retain" + str(QC_EMPTYDROPS_RETAIN) + "/{sample_name_ge}_kneeplot.png"), zip, outputqc_droplets_dir_ge=QC_OUTPUT_DIR_GE, sample_name_ge=QC_SAMPLE_NAME_GE),
        expand(os.path.normpath("{outputqc_droplets_dir_ge}" + "/QC_droplets/" + "{sample_name_ge}_saturation_plot.png"), zip, outputqc_droplets_dir_ge=QC_OUTPUT_DIR_GE, sample_name_ge=QC_SAMPLE_NAME_GE) if  str(QC_EMPTYDROPS_RETAIN) == "NULL" else expand(os.path.normpath("{outputqc_droplets_dir_ge}" + "/QC_droplets_retain" + str(QC_EMPTYDROPS_RETAIN) + "/{sample_name_ge}_saturation_plot.png"), zip, outputqc_droplets_dir_ge=QC_OUTPUT_DIR_GE, sample_name_ge=QC_SAMPLE_NAME_GE),
        expand(os.path.normpath("{outputqc_droplets_dir_ge}" + "/QC_droplets/" + "{sample_name_ge}_QChist.png"), zip, outputqc_droplets_dir_ge=QC_OUTPUT_DIR_GE, sample_name_ge=QC_SAMPLE_NAME_GE) if  str(QC_EMPTYDROPS_RETAIN) == "NULL" else expand(os.path.normpath("{outputqc_droplets_dir_ge}" + "/QC_droplets_retain" + str(QC_EMPTYDROPS_RETAIN) + "/{sample_name_ge}_QChist.png"), zip, outputqc_droplets_dir_ge=QC_OUTPUT_DIR_GE, sample_name_ge=QC_SAMPLE_NAME_GE),
        expand(os.path.normpath("{outputqc_droplets_dir_ge}" + "/QC_droplets/" + "{sample_name_ge}_QC_NON-NORMALIZED.rda"), zip, outputqc_droplets_dir_ge=QC_OUTPUT_DIR_GE, sample_name_ge=QC_SAMPLE_NAME_GE) if  str(QC_EMPTYDROPS_RETAIN) == "NULL" else expand(os.path.normpath("{outputqc_droplets_dir_ge}" + "/QC_droplets_retain" + str(QC_EMPTYDROPS_RETAIN) + "/{sample_name_ge}_QC_NON-NORMALIZED.rda"), zip, outputqc_droplets_dir_ge=QC_OUTPUT_DIR_GE, sample_name_ge=QC_SAMPLE_NAME_GE)
        ]
    if "Filtering_GE" in STEPS:
        targets["Filtering_GE"]=[
            expand(os.path.normpath("{output_filtering_dir_ge}/" + FILTERS_FOLDER + "/DOUBLETSKEPT/" + "{sample_name_ge}_QChist.png"), zip, output_filtering_dir_ge=FILERING_OUTPUT_DIR_GE, sample_name_ge=FILERING_SAMPLE_NAME_GE),
            expand(os.path.normpath("{output_filtering_dir_ge}/" + FILTERS_FOLDER + "/DOUBLETSKEPT/" + "{sample_name_ge}_DOUBLETSKEPT_NON-NORMALIZED.rda"), zip, output_filtering_dir_ge=FILERING_OUTPUT_DIR_GE, sample_name_ge=FILERING_SAMPLE_NAME_GE),
            expand(os.path.normpath("{output_filtering_dir_ge}/" + FILTERS_FOLDER + "/DOUBLETSKEPT/LogNormalize/pca/dims15_res0.8/technical/" + "{sample_name_ge}_technical_MULTI_ALL_uMAPs.png"), zip, output_filtering_dir_ge=FILERING_OUTPUT_DIR_GE, sample_name_ge=FILERING_SAMPLE_NAME_GE)
        ]
        if FILERING_DOUBLET_FILTER_METHOD_NAME == "none":
            targets["Filtering_GE"].append(expand(os.path.normpath("{output_filtering_dir_ge}/" + FILTERS_FOLDER + "/DOUBLETSKEPT/" + "/{sample_name_ge}_stat.txt"), zip, output_filtering_dir_ge=FILERING_OUTPUT_DIR_GE, sample_name_ge=FILERING_SAMPLE_NAME_GE))
        else:
            targets["Filtering_GE"].append(expand(os.path.normpath("{output_filtering_dir_ge}/" + FILTERS_FOLDER + "/DOUBLETSFILTER_" + FILERING_DOUBLET_FILTER_METHOD_NAME + "/{sample_name_ge}_QChist.png"), zip, output_filtering_dir_ge=FILERING_OUTPUT_DIR_GE, sample_name_ge=FILERING_SAMPLE_NAME_GE))
            targets["Filtering_GE"].append(expand(os.path.normpath("{output_filtering_dir_ge}/" + FILTERS_FOLDER + "/DOUBLETSFILTER_" + FILERING_DOUBLET_FILTER_METHOD_NAME + "/{sample_name_ge}_FILTERED_NON-NORMALIZED.rda"), zip, output_filtering_dir_ge=FILERING_OUTPUT_DIR_GE, sample_name_ge=FILERING_SAMPLE_NAME_GE))
            targets["Filtering_GE"].append(expand(os.path.normpath("{output_filtering_dir_ge}/" + FILTERS_FOLDER + "/DOUBLETSFILTER_" + FILERING_DOUBLET_FILTER_METHOD_NAME + "/{sample_name_ge}_stat.txt"), zip, output_filtering_dir_ge=FILERING_OUTPUT_DIR_GE, sample_name_ge=FILERING_SAMPLE_NAME_GE))
    if "Norm_DimRed_Eval_GE" in STEPS:
        targets["Norm_DimRed_Eval_GE"]=[
        expand(os.path.normpath("{output_norm_dimred_dir_ge}" + "/" + NDRE_NORM_VTR + "/" + NDRE_DIMRED_VTR + "/" + "{sample_name_ge}_" + NDRE_NORM_VTR + "_" + NDRE_DIMRED_VTR + ".rda"), zip, output_norm_dimred_dir_ge=NDRE_OUTPUT_DIR_GE, sample_name_ge=NDRE_SAMPLE_NAME_GE)
        ]
    if "Clust_Markers_Annot_GE" in STEPS:
        targets["Clust_Markers_Annot_GE"]=[
        expand(os.path.normpath("{output_clust_markers_annot_dir_ge}" + "/" + CMA_CLUST_FOLDER + "/" + "{sample_name_ge}" + "{complement}" + "_" + str(CMA_KEEP_DIM) + "_" + str(CMA_KEEP_RES) + ".rda"), zip, output_clust_markers_annot_dir_ge = CMA_OUTPUT_DIR_GE, sample_name_ge = CMA_SAMPLE_NAME_GE, complement = CMA_COMPLEMENT)
        ]
    if "Adding_ADT" in STEPS:
        targets["Adding_ADT"]=[
        expand("{add_adt_output}_ADT.rda", add_adt_output = ADD_ADT_OUTPUT)
        ]
    if "Adding_TCR" in STEPS:
        targets["Adding_TCR"]=[
        expand("{add_tcr_output}_TCR.rda", add_tcr_output = ADD_TCR_OUTPUT)
        ]
    if "Adding_BCR" in STEPS:
        targets["Adding_BCR"]=[
        expand("{add_bcr_output}_BCR.rda", add_bcr_output = ADD_BCR_OUTPUT)
        ]
    if "Cerebro" in STEPS:
        targets["Cerebro"]=[
        expand("{cerebro_input_rda_no_extention}{cerebro_complement}", cerebro_input_rda_no_extention = CEREBRO_INPUT_RDA_NO_EXTENTION, cerebro_complement = CEREBRO_COMPLEMENT_CRB)
        ]
    if "Int_Norm_DimRed_Eval_GE" in STEPS:
        targets["Int_Norm_DimRed_Eval_GE"]=[
        expand(os.path.normpath("{output_int_norm_dimred_dir_ge}" + "/GROUPED_ANALYSIS/INTEGRATED/{name_int}/" + INT_NDRE_NORM_VTR + "/" + INT_NDRE_DIMRED_VTR + "/" + "{name_int}_" + INT_NDRE_NORM_VTR + "_" + INT_NDRE_DIMRED_VTR + ".rda"), zip, output_int_norm_dimred_dir_ge=INT_NDRE_OUTPUT_DIR_GE, name_int=INT_NDRE_NAME_INT)
        ]
    if "Int_Clust_Markers_Annot_GE" in STEPS:
        targets["Int_Clust_Markers_Annot_GE"]=[
        expand(os.path.normpath("{output_int_clust_markers_annot_dir_ge}" + "/" + INT_CMA_CLUST_FOLDER + "/" + "{name_int}" + "{int_complement}" + "_" + str(INT_CMA_KEEP_DIM) + "_" + str(INT_CMA_KEEP_RES) + ".rda"), zip, output_int_clust_markers_annot_dir_ge = INT_CMA_OUTPUT_DIR_GE, name_int = INT_CMA_NAME_INT, int_complement = INT_CMA_COMPLEMENT)
        ]
    if "Int_Adding_ADT" in STEPS:
        targets["Int_Adding_ADT"]=[
        expand("{int_add_adt_output}_ADT.rda", int_add_adt_output = INT_ADD_ADT_OUTPUT)
        ]
    if "Int_Adding_TCR" in STEPS:
        targets["Int_Adding_TCR"]=[
        expand("{int_add_tcr_output}_TCR.rda", int_add_tcr_output = INT_ADD_TCR_OUTPUT)
        ]
    if "Int_Adding_BCR" in STEPS:
        targets["Int_Adding_BCR"]=[
        expand("{int_add_bcr_output}_BCR.rda", int_add_bcr_output = INT_ADD_BCR_OUTPUT)
        ]
    if "Grp_Norm_DimRed_Eval_GE" in STEPS:
        targets["Grp_Norm_DimRed_Eval_GE"]=[
        expand(os.path.normpath("{output_grp_norm_dimred_dir_ge}" + "/GROUPED_ANALYSIS/NO_INTEGRATED/{name_grp}/" + GRP_NDRE_NORM_VTR + "/" + GRP_NDRE_DIMRED_VTR + "/" + "{name_grp}_" + GRP_NDRE_NORM_VTR + "_" + GRP_NDRE_DIMRED_VTR + ".rda"), zip, output_grp_norm_dimred_dir_ge=GRP_NDRE_OUTPUT_DIR_GE, name_grp=GRP_NDRE_NAME_GRP)
        ]
    if "Grp_Clust_Markers_Annot_GE" in STEPS:
        targets["Grp_Clust_Markers_Annot_GE"]=[
        expand(os.path.normpath("{output_grp_clust_markers_annot_dir_ge}" + "/" + GRP_CMA_CLUST_FOLDER + "/" + "{name_grp}" + "{grp_complement}" + "_" + str(GRP_CMA_KEEP_DIM) + "_" + str(GRP_CMA_KEEP_RES) + ".rda"), zip, output_grp_clust_markers_annot_dir_ge = GRP_CMA_OUTPUT_DIR_GE, name_grp = GRP_CMA_NAME_GRP, grp_complement = GRP_CMA_COMPLEMENT)
        ]
    if "Grp_Adding_ADT" in STEPS:
        targets["Grp_Adding_ADT"]=[
        expand("{grp_add_adt_output}_ADT.rda", grp_add_adt_output = GRP_ADD_ADT_OUTPUT)
        ]
    if "Grp_Adding_TCR" in STEPS:
        targets["Grp_Adding_TCR"]=[
        expand("{grp_add_tcr_output}_TCR.rda", grp_add_tcr_output = GRP_ADD_TCR_OUTPUT)
        ]
    if "Grp_Adding_BCR" in STEPS:
        targets["Grp_Adding_BCR"]=[
        expand("{grp_add_bcr_output}_BCR.rda", grp_add_bcr_output = GRP_ADD_BCR_OUTPUT)
        ]
    if "Alignment_countTable_LR_GE" in STEPS:
        targets["Alignment_countTable_LR_GE"]=[
        expand(ALIGN_OUTPUT_DIR_GE+"/samplesheet/{sample_name_ge}_samplesheet.csv",sample_name_ge=ALIGN_SAMPLE_NAME),
        expand(ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}_GE/{sample_name_ge}_GE/gene_raw_feature_bc_matrix/matrix.mtx.gz",sample_name_ge=ALIGN_SAMPLE_NAME),
        expand(ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}_GE/{sample_name_ge}_GE/transcript_raw_feature_bc_matrix/matrix.mtx.gz",sample_name_ge=ALIGN_SAMPLE_NAME),
        expand(ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}_GE/{sample_name_ge}_GE/tagged.bam",sample_name_ge=ALIGN_SAMPLE_NAME)
        ]
    if "Isoform_Markers_GE" in STEPS:
        targets["Isoform_Markers_GE"]=[
        #expand("{start}/{sample_name_ge}{end}",zip,start=COMPLEMENT_ISOFORM_OUTPUT_RDA_LIST_START,sample_name_ge=ISOFORM_SAMPLE_NAME_GE,end=COMPLEMENT_ISOFORM_OUTPUT_RDA_LIST_END)
        expand(ISFORM_MTX_INTPUT[0]+"{sample_name_ge}/{sample_name_ge}/{sample_name_ge}_"+CMA_CLUST_FOLDER+"_sobj_isoform.rda",sample_name_ge=ISOFORM_SAMPLE_NAME_GE)
        ]
    if "Variants_Markers_GE" in STEPS:
        targets["Variants_Markers_GE"]=[
        expand(ALIGN_OUTPUT_DIR_GE+ "/{sample_name_ge}/" + FILTERS_FOLDER + "/DOUBLETSFILTER_all/{sample_name_ge}_FILTERED_barcodes.txt",sample_name_ge=ALIGN_SAMPLE_NAME_GE),
        expand(ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/split_chrom/tagged_sort_"+chr_list[0]+".bam",sample_name_ge=ALIGN_SAMPLE_NAME_GE),
        expand(ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/barcodes_batch/{sample_name_ge}_barcodes_list_{batch_number}.txt",sample_name_ge=ALIGN_SAMPLE_NAME_GE,batch_number=batch_number_seq),
        expand(ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/tmp_barcodes/DONE_{sample_name_ge}_barcodes_list_{batch_number}.txt",sample_name_ge=ALIGN_SAMPLE_NAME_GE,batch_number=batch_number_seq),
        expand(ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/tmp_merge/DONE_{sample_name_ge}_barcodes_list_{batch_number}.txt",sample_name_ge=ALIGN_SAMPLE_NAME_GE,batch_number=batch_number_seq),
        expand(ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/tmp_mpileup/DONE_{sample_name_ge}_barcodes_list_{batch_number}.txt",sample_name_ge=ALIGN_SAMPLE_NAME_GE,batch_number=batch_number_seq),
        expand(ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/matrix_SNP_{sample_name_ge}.tsv",sample_name_ge=ALIGN_SAMPLE_NAME_GE,batch_number=batch_number_seq)
        ]
    if "Variants_Calling_GE" in STEPS:
        targets["Variants_Calling_GE"]=[
        expand(ALIGN_OUTPUT_DIR_GE+ "/{sample_name_ge}/" + FILTERS_FOLDER + "/DOUBLETSFILTER_all/{sample_name_ge}_FILTERED_barcodes.txt",sample_name_ge=ALIGN_SAMPLE_NAME_GE),
        expand(ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/filtered_tagged.bam",sample_name_ge=ALIGN_SAMPLE_NAME_GE),
        expand(ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/filtered_sort_tagged.bam",sample_name_ge=ALIGN_SAMPLE_NAME_GE),
        expand(ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/{sample_name_ge}_longshot.vcf",sample_name_ge=ALIGN_SAMPLE_NAME_GE),
        expand(ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/{sample_name_ge}_longshot.vcf.gz",sample_name_ge=ALIGN_SAMPLE_NAME_GE),
        expand(ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/{sample_name_ge}_longshot_filtered.vcf.gz",sample_name_ge=ALIGN_SAMPLE_NAME_GE),
        expand(ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/{sample_name_ge}_longshot_filtered.vcf.gz.tbi",sample_name_ge=ALIGN_SAMPLE_NAME_GE)
        ]
    if "Int_Isoform_Markers" in STEPS:
        targets["Int_Isoform_Markers"]=[expand(INT_ISOFORM_OUPUT_DIR_GLOBAL+"GROUPED_ANALYSIS/INTEGRATED/{sample_int}/tmp/DONE_int_isoform_markers_"+INT_CMA_CLUST_FOLDER+".txt",sample_int=INT_ISOFORM_SAMPLE)]
        
    if "Grp_Isoform_Markers" in STEPS:
        targets["Grp_Isoform_Markers"]=[expand(GRP_ISOFORM_OUPUT_DIR_GLOBAL+"GROUPED_ANALYSIS/NO_INTEGRATED/{sample_grp}/tmp/DONE__grp_isoform_markers_"+GRP_CMA_CLUST_FOLDER+".txt",sample_grp=GRP_ISOFORM_SAMPLE)]
        
    if "Variants_Calling_Comparison_GE" in STEPS:
        targets["Variants_Calling_Comparison_GE"]=[
        expand(ouput_dir_ge+"{sample_name_ge}/{sample_name_ge}/SNV_COMPARE/barcodes_ctrl.txt",sample_name_ge=sample_ge_list),
        expand(ouput_dir_ge+"{sample_name_ge}/{sample_name_ge}/SNV_COMPARE/barcodes_ofint.txt",sample_name_ge=sample_ge_list),
        expand(ouput_dir_ge+"{sample_name_ge}/{sample_name_ge}/SNV_COMPARE/bam/{statut}/filtered_{statut}.bam",sample_name_ge=sample_ge_list,statut=["ctrl","to_test"]),
        expand(ouput_dir_ge+"{sample_name_ge}/{sample_name_ge}/SNV_COMPARE/vcf/{statut}/filtered_{statut}.vcf",sample_name_ge=sample_ge_list,statut=["ctrl","to_test"]),
        expand(vcf_out=ouput_dir_ge+"{sample_name_ge}/{sample_name_ge}/SNV_COMPARE/vcf/{statut}/filtered_{statut}.vcf.gz",sample_name_ge=sample_ge_list,statut=["ctrl","to_test"]),
        expand(ouput_dir_ge+"{sample_name_ge}/{sample_name_ge}/SNV_COMPARE/isec/README.txt",sample_name_ge=sample_ge_list)
        ]
        
    return targets
