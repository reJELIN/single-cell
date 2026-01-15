"""
##########################################################################
This rule make the translation of seurat file into cerebro file in single-cell RNA-seq.
##########################################################################
"""

"""
This function allows to determine the input .rda file.
"""
def cerebro_input(wildcards):
    return wildcards.cerebro_input_rda_no_extention + ".rda"

"""
This function allows to determine the singularity binding parameters.
"""
def cerebro_params_sing(wildcards):
    rda_crb_folder = os.path.dirname(wildcards.cerebro_input_rda_no_extention)
    concat = " -B " + PIPELINE_FOLDER + ":" + os.path.normpath("/WORKDIR/" + PIPELINE_FOLDER) + " -B " + rda_crb_folder + ":" + os.path.normpath("/WORKDIR/" + rda_crb_folder)
    if CEREBRO_GMT_FILE != "NULL":
        gmt_folder = os.path.dirname(CEREBRO_GMT_FILE)
        concat = concat + " -B " + gmt_folder + ":" + os.path.normpath("/WORKDIR/" + gmt_folder)
    if CEREBRO_METADATA_FILE != "NULL":
        for metadatafile in list(dict.fromkeys(CEREBRO_METADATA_FILE.split(","))):
            metadatafile = os.path.dirname(metadatafile)
            concat = concat + " -B " + metadatafile + ":" + os.path.normpath("/WORKDIR/" + metadatafile)
    return concat

"""
This rule launches the R script to translate seurat file into cerebro file.
"""
rule cerebro:
    input:
        cerebro_rda_file = cerebro_input
    output:
        cerebro_crb_file = expand("{{cerebro_input_rda_no_extention}}{cerebro_complement}", cerebro_complement = CEREBRO_COMPLEMENT_CRB)
    params:
        sing_bind = cerebro_params_sing,
        pipeline_folder = os.path.normpath("/WORKDIR/" + PIPELINE_FOLDER),
        input_rda = lambda wildcards, input: os.path.normpath("/WORKDIR/" + input[0]),
        SING_CEREBRO_GMT_FILE = os.path.normpath("/WORKDIR/" + CEREBRO_GMT_FILE) if CEREBRO_GMT_FILE != "NULL" else "NULL",
        SING_CEREBRO_METADATA_FILE = ','.join([os.path.normpath("/WORKDIR/" + x) for x in CEREBRO_METADATA_FILE.split(',')]) if CEREBRO_METADATA_FILE != "NULL" else "NULL",
        cerebro_author_name = CEREBRO_AUTHOR_NAME,
        cerebro_author_mail = CEREBRO_AUTHOR_MAIL,
        cerebro_version = CEREBRO_VERSION,
        cerebro_groups = CEREBRO_GROUPS,
        cerebro_remove_other_red = CEREBRO_REMOVE_OTHER_RED,
        cerebro_remove_other_ident = CEREBRO_REMOVE_OTHER_IDENT,
        cerebro_remove_mt = CEREBRO_REMOVE_MT,
        cerebro_remove_crb = CEREBRO_REMOVE_CRB,
        cerebro_remove_str = CEREBRO_REMOVE_STR,
        cerebro_only_pos_de = CEREBRO_ONLY_POS_DE,
        cerebro_remove_custom_de = CEREBRO_REMOVE_CUSTOM_DE,
        singularity_env_cerebro = SINGULARITY_ENV_CEREBRO
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(15000 * attempt , 102400)),
        time_min = (lambda wildcards, attempt: min(attempt * 60, 200))
    shell:
        """
        export TMPDIR=$TMPDIR
        TMP_DIR=$(mktemp -d -t sc_pipeline-XXXXXXXXXX) && \
        singularity exec --contain --home $TMP_DIR:$HOME -B $TMP_DIR:/tmp {params.sing_bind} \
        {params.singularity_env_cerebro} \
        Rscript {params.pipeline_folder}/scripts/pipeline_CEREBRO.R \
        --input.rda.ge {params.input_rda} \
        --author.name {params.cerebro_author_name} \
        --author.mail {params.cerebro_author_mail} \
        --nthreads {threads} \
        --pipeline.path {params.pipeline_folder} \
        --version {params.cerebro_version} \
        --groups {params.cerebro_groups} \
        --remove.other.reductions {params.cerebro_remove_other_red} \
        --remove.other.idents {params.cerebro_remove_other_ident} \
        --remove.mt.genes {params.cerebro_remove_mt} \
        --remove.crb.genes {params.cerebro_remove_crb} \
        --remove.str.genes {params.cerebro_remove_str} \
        --only.pos.DE {params.cerebro_only_pos_de} \
        --remove.custom.DE {params.cerebro_remove_custom_de} \
        --gmt.file {params.SING_CEREBRO_GMT_FILE} \
        --metadata.file {params.SING_CEREBRO_METADATA_FILE}
        """
