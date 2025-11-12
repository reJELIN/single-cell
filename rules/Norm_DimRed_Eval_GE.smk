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
        SING_NDRE_METADATA_FILE = ','.join([os.path.normpath("/WORKDIR/" + x) for x in NDRE_METADATA_FILE.split(',')]) if NDRE_METADATA_FILE != "NULL" else "NULL",
        ndre_author_name = NDRE_AUTHOR_NAME,
        ndre_author_mail = NDRE_AUTHOR_MAIL,
        ndre_eval_markers = NDRE_EVAL_MARKERS,
        ndre_features_n = NDRE_FEATURES_N,
        ndre_norm_method = NDRE_NORM_METHOD,
        ndre_dimred_method = NDRE_DIMRED_METHOD,
        ndre_vtr_biases = NDRE_VTR_BIASES,
        ndre_vtr_scale = NDRE_VTR_SCALE,
        ndre_dim_max = NDRE_DIM_MAX,
        ndre_dim_min = NDRE_DIM_MIN,
        ndre_dim_steps = NDRE_DIM_STEPS,
        ndre_res_max = NDRE_RES_MAX,
        ndre_res_min = NDRE_RES_MIN,
        ndre_res_steps = NDRE_RES_STEPS,
        singularity_env =  SINGULARITY_ENV
    threads:
        4
    resources:
        mem_mb = (lambda wildcards, attempt: min(5120 + attempt * 5120, 61440)),
        time_min = (lambda wildcards, attempt: min(attempt * 120, 200))
    shell:
        """
        export TMPDIR=$TMPDIR
        TMP_DIR=$(mktemp -d -t sc_pipeline-XXXXXXXXXX) && \
        singularity exec --no-home -B $TMP_DIR:/tmp {params.sing_bind} \
        {params.singularity_env} \
        Rscript {params.pipeline_folder}/scripts/pipeline_part3.R \
        --input.rda.ge {params.input_rda} \
        --output.dir.ge {params.output_folder} \
        --author.name {params.ndre_author_name} \
        --author.mail {params.ndre_author_mail} \
        --nthreads {threads} \
        --pipeline.path {params.pipeline_folder} \
        --eval.markers {params.ndre_eval_markers} \
        --features.n {params.ndre_features_n} \
        --norm.method {params.ndre_norm_method} \
        --dimred.method {params.ndre_dimred_method} \
        --vtr.biases {params.ndre_vtr_biases} \
        --vtr.scale {params.ndre_vtr_scale} \
        --dims.max {params.ndre_dim_max} \
        --dims.min {params.ndre_dim_min} \
        --dims.steps {params.ndre_dim_steps} \
        --res.max {params.ndre_res_max} \
        --res.min {params.ndre_res_min} \
        --res.steps {params.ndre_res_steps} \
        --metadata.file {params.SING_NDRE_METADATA_FILE} && \
        rm -r $TMP_DIR || rm -r $TMP_DIR
        """
