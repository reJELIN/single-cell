"""
##########################################################################
This rule make the normalization, dimensions reduction and evaluation in single-cell RNA-seq.
##########################################################################
"""

"""
This function allows to determine the input .rda file.
"""
def norm_dimred_input_ge(wildcards):
    return dic_NDRE_INFO[wildcards.sample_name_ge]['NDRE_INPUT_RDA']

"""
This function allows to determine the singularity binding parameters.
"""
def norm_dimred_params_sing(wildcards):
    rda_folder = os.path.dirname(dic_NDRE_INFO[wildcards.sample_name_ge]['NDRE_INPUT_RDA'])
    output_folder = wildcards.output_norm_dimred_dir_ge
    concat = " -B " + PIPELINE_FOLDER + ":" + os.path.normpath("/WORKDIR/" + PIPELINE_FOLDER) + " -B " + rda_folder + ":" + os.path.normpath("/WORKDIR/" + rda_folder) + " -B " + output_folder + ":" + os.path.normpath("/WORKDIR/" + output_folder)
    if NDRE_METADATA_FILE != "NULL":
        for metadatafile in list(dict.fromkeys(NDRE_METADATA_FILE.split(","))):
            metadatafile = os.path.dirname(metadatafile)
            concat = concat + " -B " + metadatafile + ":" + os.path.normpath("/WORKDIR/" + metadatafile)
    return concat

"""
This rule launches the R script to apply normalization, dimensions reduction and evaluation of parameters.
"""
rule norm_dimred_ge:
    input:
        rda_file = norm_dimred_input_ge
    output:
        ndre_Eval_rda_file = os.path.normpath("{output_norm_dimred_dir_ge}" + "/" + NDRE_NORM_VTR + "/" + NDRE_DIMRED_VTR + "/" + "{sample_name_ge}_" + NDRE_NORM_VTR + "_" + NDRE_DIMRED_VTR + ".rda")
    params:
        sing_bind = norm_dimred_params_sing,
        pipeline_folder = os.path.normpath("/WORKDIR/" + PIPELINE_FOLDER),
        input_rda = lambda wildcards, input: os.path.normpath("/WORKDIR/" + input[0]),
        output_folder = os.path.normpath("/WORKDIR/" + "{output_norm_dimred_dir_ge}") + "/",
        SING_NDRE_METADATA_FILE = ','.join([os.path.normpath("/WORKDIR/" + x) for x in NDRE_METADATA_FILE.split(',')]) if NDRE_METADATA_FILE != "NULL" else "NULL"
    threads:
        4
    resources:
        mem_mb = (lambda wildcards, attempt: min(5120 + attempt * 5120, 61440)),
        time_min = (lambda wildcards, attempt: min(attempt * 120, 200))
    shell:
        """
        export TMPDIR={GLOBAL_TMP}
        TMP_DIR=$(mktemp -d -t sc_pipeline-XXXXXXXXXX) && \
        singularity exec --no-home -B $TMP_DIR:/tmp {params.sing_bind} \
        {SINGULARITY_ENV} \
        Rscript {params.pipeline_folder}/scripts/pipeline_part3.R \
        --input.rda.ge {params.input_rda} \
        --output.dir.ge {params.output_folder} \
        --author.name {NDRE_AUTHOR_NAME} \
        --author.mail {NDRE_AUTHOR_MAIL} \
        --nthreads {threads} \
        --pipeline.path {params.pipeline_folder} \
        --eval.markers {NDRE_EVAL_MARKERS} \
        --features.n {NDRE_FEATURES_N} \
        --norm.method {NDRE_NORM_METHOD} \
        --dimred.method {NDRE_DIMRED_METHOD} \
        --vtr.biases {NDRE_VTR_BIASES} \
        --vtr.scale {NDRE_VTR_SCALE} \
        --dims.max {NDRE_DIM_MAX} \
        --dims.min {NDRE_DIM_MIN} \
        --dims.steps {NDRE_DIM_STEPS} \
        --res.max {NDRE_RES_MAX} \
        --res.min {NDRE_RES_MIN} \
        --res.steps {NDRE_RES_STEPS} \
        --metadata.file {params.SING_NDRE_METADATA_FILE} && \
        rm -r $TMP_DIR || rm -r $TMP_DIR
        """
