"""
##########################################################################
This rule make the droplets control-quality of genes expression in single-cell RNA-seq.
##########################################################################
"""

wildcard_constraints:
    sample_name_ge=".+_GE"

"""
This function allows to determine the input alignment folder/files.
"""
def QC_droplets_input_ge(wildcards):
    kallisto_folder = dic_SAMPLE_NAME_GE_INFO[wildcards.sample_name_ge]['QC_INPUT_DIR']
    if "Alignment_countTable_GE" in STEPS:
        mtx_file = os.path.normpath(kallisto_folder + "/" + wildcards.sample_name_ge + ".mtx")
        barcodes_file = os.path.normpath(kallisto_folder + "/" + wildcards.sample_name_ge + ".barcodes.txt")
        genes_file = os.path.normpath(kallisto_folder + "/" + wildcards.sample_name_ge + ".genes.txt")
        files=[mtx_file, barcodes_file, genes_file]
    else:
        files=[kallisto_folder]
    return files

"""
This function allows to determine the input alignment folder for params section.
"""
def QC_params_input_folder(wildcards):
    input_folder = os.path.normpath("/WORKDIR/" + dic_SAMPLE_NAME_GE_INFO[wildcards.sample_name_ge]['QC_INPUT_DIR']) + "/"
    return input_folder

"""
This function allows to determine the singularity binding parameters.
"""
def QC_params_sing(wildcards):
    kallisto_folder = dic_SAMPLE_NAME_GE_INFO[wildcards.sample_name_ge]['QC_INPUT_DIR']
    output_folder = wildcards.outputqc_droplets_dir_ge + "/"
    concat = " -B " + PIPELINE_FOLDER + ":/WORKDIR/" + PIPELINE_FOLDER + " -B " + kallisto_folder + ":" + os.path.normpath("/WORKDIR/" + kallisto_folder) + " -B " + output_folder + ":" + os.path.normpath("/WORKDIR/" + output_folder)
    if QC_MT_FILE != "NULL": concat = concat + " -B " + QC_MT_FILE + ":" + os.path.normpath("/WORKDIR/" + QC_MT_FILE)
    if QC_RB_FILE != "NULL": concat = concat + " -B " + QC_RB_FILE + ":" + os.path.normpath("/WORKDIR/" + QC_RB_FILE)
    if QC_ST_FILE != "NULL": concat = concat + " -B " + QC_ST_FILE + ":" + os.path.normpath("/WORKDIR/" + QC_ST_FILE)
    if QC_TRANSLATION_FILE != "NULL": concat = concat + " -B " + QC_TRANSLATION_FILE + ":" + os.path.normpath("/WORKDIR/" + QC_TRANSLATION_FILE)
    if QC_METADATA_FILE != "NULL":
        for metadatafile in list(dict.fromkeys(QC_METADATA_FILE.split(","))):
            metadatafile = os.path.dirname(metadatafile)
            concat = concat + " -B " + metadatafile + ":" + os.path.normpath("/WORKDIR/" + metadatafile)
    return concat

"""
This rule launches R scipt to read count matrix and perform droplets control-quality.
"""
rule QC_droplets_ge:
    input:
        QC_droplets_input_ge
    output:
        kneeplot_file = os.path.normpath("{outputqc_droplets_dir_ge}" + "/QC_droplets/" + "{sample_name_ge}_kneeplot.png") if  str(QC_EMPTYDROPS_RETAIN) == "NULL" else os.path.normpath("{outputqc_droplets_dir_ge}" + "/QC_droplets_retain" + str(QC_EMPTYDROPS_RETAIN) + "/{sample_name_ge}_kneeplot.png"),
        saturation_file = os.path.normpath("{outputqc_droplets_dir_ge}" + "/QC_droplets/" + "{sample_name_ge}_saturation_plot.png") if  str(QC_EMPTYDROPS_RETAIN) == "NULL" else os.path.normpath("{outputqc_droplets_dir_ge}" + "/QC_droplets_retain" + str(QC_EMPTYDROPS_RETAIN) + "/{sample_name_ge}_saturation_plot.png"),
        QC_hist_unfiltred_file =  os.path.normpath("{outputqc_droplets_dir_ge}" + "/QC_droplets/" + "{sample_name_ge}_QChist.png") if str(QC_EMPTYDROPS_RETAIN) == "NULL" else os.path.normpath("{outputqc_droplets_dir_ge}" +  "/QC_droplets_retain" + str(QC_EMPTYDROPS_RETAIN) + "/{sample_name_ge}_QChist.png"),
        unfiltred_non_norm_rda = os.path.normpath("{outputqc_droplets_dir_ge}" + "/QC_droplets/" + "{sample_name_ge}_QC_NON-NORMALIZED.rda") if  str(QC_EMPTYDROPS_RETAIN) == "NULL" else os.path.normpath("{outputqc_droplets_dir_ge}" + "/QC_droplets_retain" + str(QC_EMPTYDROPS_RETAIN) + "/{sample_name_ge}_QC_NON-NORMALIZED.rda")
    params:
        sing_bind = QC_params_sing,
        pipeline_folder = os.path.normpath("/WORKDIR/" + PIPELINE_FOLDER),
        # input_folder = lambda wildcards, input: os.path.normpath("/WORKDIR/" + input[0]) + "/",
        input_folder = QC_params_input_folder,
        output_folder = os.path.normpath("/WORKDIR/" + "{outputqc_droplets_dir_ge}") + "/",
        SING_QC_MT_FILE = os.path.normpath("/WORKDIR/" + QC_MT_FILE) if QC_MT_FILE != "NULL" else "NULL",
        SING_QC_RB_FILE = os.path.normpath("/WORKDIR/" + QC_RB_FILE) if QC_RB_FILE != "NULL" else "NULL",
        SING_QC_ST_FILE = os.path.normpath("/WORKDIR/" + QC_ST_FILE) if QC_ST_FILE != "NULL" else "NULL",
        SING_QC_TRANSLATION_FILE = os.path.normpath("/WORKDIR", QC_TRANSLATION_FILE) if QC_TRANSLATION_FILE != "NULL" else "NULL",
        SING_QC_METADATA_FILE = ','.join([os.path.normpath("/WORKDIR/" + x) for x in QC_METADATA_FILE.split(',')]) if QC_METADATA_FILE != "NULL" else "NULL",
        qc_species = QC_SPECIES,
        qc_author = QC_AUTHOR_NAME,
        qc_author_mail = QC_AUTHOR_MAIL,
        qc_emptydrop_fdr = QC_EMPTYDROPS_FDR,
        qc_droplets_limit = QC_DROPLETS_LIMIT,
        qc_emptydrop_retain = QC_EMPTYDROPS_RETAIN,
        qc_translation_bool = QC_TRANSLATION_BOOL,
        qc_pcmito_min = QC_PCMITO_MIN,
        qc_pcmito_max = QC_PCMITO_MAX,
        qc_pcribo_min = QC_PCRIBO_MIN,
        qc_pcribo_max = QC_PC_RIBO_MAX,
        qc_min_features = QC_MIN_FEATURES,
        qc_min_counts = QC_MIN_COUNTS,
        qc_min_cells = QC_MIN_CELLS,
        sing_env = SINGULARITY_ENV
    threads:
        2
    resources:
        mem_mb = (lambda wildcards, attempt: min(50000 + attempt * 3072, 40960)), #previous max: 40960 min: 3072
        time_min = (lambda wildcards, attempt: min(attempt * 90, 200)) #previous max: 200
    shell:
        """
        export TMPDIR=$TMPDIR
        TMP_DIR=$(mktemp -d -t sc_pipeline-XXXXXXXXXX) && \
        singularity exec --no-home -B $TMP_DIR:/tmp {params.sing_bind} \
        {params.sing_env} \
        Rscript {params.pipeline_folder}/scripts/pipeline_part1.R \
        --input.dir.ge {params.input_folder} \
        --output.dir.ge {params.output_folder} \
        --sample.name.ge {wildcards.sample_name_ge} \
        --species {params.qc_species} \
        --author.name {params.qc_author} \
        --author.mail {params.qc_author_mail} \
        --nthreads {threads} \
        --pipeline.path {params.pipeline_folder} \
        --emptydrops.fdr {params.qc_emptydrop_fdr} \
        --droplets.limit {params.qc_droplets_limit} \
        --emptydrops.retain {params.qc_emptydrop_retain} \
        --translation {params.qc_translation_bool} \
        --pcmito.min {params.qc_pcmito_min} \
        --pcmito.max {params.qc_pcmito_max} \
        --pcribo.min {params.qc_pcribo_min} \
        --pcribo.max {params.qc_pcribo_max} \
        --min.features {params.qc_min_features} \
        --min.counts {params.qc_min_counts} \
        --min.cells {params.qc_min_cells} \
        --mt.genes.file {params.SING_QC_MT_FILE} \
        --crb.genes.file {params.SING_QC_RB_FILE} \
        --str.genes.file {params.SING_QC_ST_FILE} \
        --translation.file {params.SING_QC_TRANSLATION_FILE} \
        --metadata.file {params.SING_QC_METADATA_FILE} && \
        rm -r $TMP_DIR || rm -r $TMP_DIR
        """
