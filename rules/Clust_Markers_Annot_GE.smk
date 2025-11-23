"""
##########################################################################
This rule make the clustering, to find markers genes and to apply annotations of cell types in single-cell RNA-seq.
##########################################################################
"""

wildcard_constraints:
    sample_name_ge=".+_GE"

"""
This function allows to determine the input .rda file.
"""
def clust_markers_annot_input_ge(wildcards):
    return dic_CMA_INFO[wildcards.sample_name_ge]['CMA_INPUT_RDA']

"""
This function allows to determine the singularity binding parameters.
"""
def clust_markers_annot_params_sing(wildcards):
    rda_folder = os.path.dirname(dic_CMA_INFO[wildcards.sample_name_ge]['CMA_INPUT_RDA'])
    output_folder = wildcards.output_clust_markers_annot_dir_ge
    concat = " -B " + PIPELINE_FOLDER + ":" + os.path.normpath("/WORKDIR/" + PIPELINE_FOLDER) + " -B " + rda_folder + ":" + os.path.normpath("/WORKDIR/" + rda_folder) + " -B " + output_folder + ":" + os.path.normpath("/WORKDIR/" + output_folder)
    if CMA_MARKFILE != "NULL":
        for markfile in list(dict.fromkeys(CMA_MARKFILE.split(","))):
            markfile = os.path.dirname(markfile)
            concat = concat + " -B " + markfile + ":" + os.path.normpath("/WORKDIR/" + markfile)
    if CMA_METADATA_FILE != "NULL":
        for metadatafile in list(dict.fromkeys(CMA_METADATA_FILE.split(","))):
            metadatafile = os.path.dirname(metadatafile)
            concat = concat + " -B " + metadatafile + ":" + os.path.normpath("/WORKDIR/" + metadatafile)
    if CMA_CUSTOM_SCE_REF != "NULL":
        for custom_cse_ref in list(dict.fromkeys(CMA_CUSTOM_SCE_REF.split(","))):
            custom_cse_ref = os.path.dirname(custom_cse_ref)
            concat = concat + " -B " + custom_cse_ref + ":" + os.path.normpath("/WORKDIR/" + custom_cse_ref)
    if CMA_CUSTOM_MARKERS_REF != "NULL":
        for custom_marker_ref in list(dict.fromkeys(CMA_CUSTOM_MARKERS_REF.split(","))):
            custom_marker_ref = os.path.dirname(custom_marker_ref)
            concat = concat + " -B " + custom_marker_ref + ":" + os.path.normpath("/WORKDIR/" + custom_marker_ref)
    return concat

"""
This rule launches the R script to make the clustering, to find markers genes and to apply annotations of cell types.
"""
rule clust_markers_annot_ge:
    input:
        cma_file = clust_markers_annot_input_ge
    output:
        cma_rda_file = os.path.normpath("{output_clust_markers_annot_dir_ge}" + "/" + CMA_CLUST_FOLDER + "/" + "{sample_name_ge}" + "{complement}" + str(CMA_KEEP_DIM) + "_" + str(CMA_KEEP_RES) + ".rda")
    params:
        sing_bind = clust_markers_annot_params_sing,
        pipeline_folder = os.path.normpath("/WORKDIR/" + PIPELINE_FOLDER),
        input_rda = lambda wildcards, input: os.path.normpath("/WORKDIR/" + input[0]),
        output_folder = os.path.normpath("/WORKDIR/" + "{output_clust_markers_annot_dir_ge}") + "/",
        SING_CMA_MARKFILE = ','.join([os.path.normpath("/WORKDIR/" + x) for x in CMA_MARKFILE.split(',')]) if CMA_MARKFILE != "NULL" else "NULL",
        SING_CMA_METADATA_FILE = ','.join([os.path.normpath("/WORKDIR/" + x) for x in CMA_METADATA_FILE.split(',')]) if CMA_METADATA_FILE != "NULL" else "NULL",
        SING_CMA_CUSTOM_SCE_REF = ','.join([os.path.normpath("/WORKDIR/" + x) for x in CMA_CUSTOM_SCE_REF.split(',')]) if CMA_CUSTOM_SCE_REF != "NULL" else "NULL",
        SING_CMA_CUSTOM_MARKERS_REF = ','.join([os.path.normpath("/WORKDIR/" + x) for x in CMA_CUSTOM_MARKERS_REF.split(',')]) if CMA_CUSTOM_MARKERS_REF != "NULL" else "NULL",
        cma_author_name = CMA_AUTHOR_NAME,
        cma_author_mail = CMA_AUTHOR_MAIL,
        cma_keep_dim = CMA_KEEP_DIM,
        cma_keep_res = CMA_KEEP_RES,
        cma_cfr_minscore = CMA_CFR_MINSCORE,
        cma_sr_minscore = CMA_SR_MINSCORE,
        sing_env = SINGULARITY_ENV
    threads:
        1
    resources:
        mem_mb = lambda wildcards, attempt: 10240 + attempt * 5120,
        time_min = lambda wildcards, attempt: attempt * 120
    envmodules:
        "singularity-ce/4.2.2"
    shell:
        """
        export TMPDIR=$TMPDIR
        TMP_DIR=$(mktemp -d -t sc_pipeline-XXXXXXXXXX) && \
        singularity exec --no-home -B $TMP_DIR:/tmp -B $TMP_DIR:$HOME {params.sing_bind} \
        {params.sing_env} \
        Rscript {params.pipeline_folder}/scripts/pipeline_part4.R \
        --input.rda.ge {params.input_rda} \
        --output.dir.ge {params.output_folder} \
        --author.name {params.cma_author_name} \
        --author.mail {params.cma_author_mail} \
        --nthreads {threads} \
        --pipeline.path {params.pipeline_folder} \
        --markfile  {params.SING_CMA_MARKFILE} \
        --custom.sce.ref {params.SING_CMA_CUSTOM_SCE_REF} \
        --custom.markers.ref {params.SING_CMA_CUSTOM_MARKERS_REF} \
        --keep.dims {params.cma_keep_dim} \
        --keep.res {params.cma_keep_res} \
        --cfr.minscore {params.cma_cfr_minscore} \
        --sr.minscore {params.cma_sr_minscore} \
        --metadata.file {params.SING_CMA_METADATA_FILE} && \
        rm -r $TMP_DIR || rm -r $TMP_DIR
        """
