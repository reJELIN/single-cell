"""
##########################################################################
This rule make the integration, normalization, dimensions reduction and evaluation in single-cell RNA-seq.
##########################################################################
"""

# wildcard_constraints:
#    name_int = INT_NDRE_NAME_INT

"""
This function allows to determine the input .rda file.
"""
def int_norm_dimred_input_ge_snakemake(wildcards):
    return dic_INT_NDRE_INFO[wildcards.name_int]['INT_NDRE_INPUT_LIST_RDA'].split(',')

def int_norm_dimred_input_ge_R(wildcards):
    return ','.join([os.path.normpath("/WORKDIR/" + x) for x in dic_INT_NDRE_INFO[wildcards.name_int]['INT_NDRE_INPUT_LIST_RDA'].split(',')])

"""
This function allows to determine the singularity binding parameters.
"""
def int_norm_dimred_params_sing(wildcards):
    sing_rda_folder = ''.join([" -B " + os.path.dirname(x) + ":" + os.path.normpath("/WORKDIR/" + os.path.dirname(x)) for x in dic_INT_NDRE_INFO[wildcards.name_int]['INT_NDRE_INPUT_LIST_RDA'].split(',')])
    output_folder = wildcards.output_int_norm_dimred_dir_ge
    concat = " -B " + PIPELINE_FOLDER + ":" + os.path.normpath("/WORKDIR/" + PIPELINE_FOLDER) + sing_rda_folder + " -B " + output_folder + ":" + os.path.normpath("/WORKDIR/" + output_folder)
    if INT_NDRE_METADATA_FILE != "NULL":
        for metadatafile in list(dict.fromkeys(INT_NDRE_METADATA_FILE.split(","))):
            metadatafile = os.path.dirname(metadatafile)
            concat = concat + " -B " + metadatafile + ":" + os.path.normpath("/WORKDIR/" + metadatafile)
    return concat

"""
This rule launches the R script to apply normalization, dimensions reduction and evaluation of parameters.
"""
rule int_norm_dimred_ge:
    input:
        int_rda_file = int_norm_dimred_input_ge_snakemake
    output:
        int_ndre_Eval_rda_file = os.path.normpath("{output_int_norm_dimred_dir_ge}" + "/GROUPED_ANALYSIS/INTEGRATED/{name_int}/" + INT_NDRE_NORM_VTR + "/" + INT_NDRE_DIMRED_VTR + "/" + "{name_int}_" + INT_NDRE_NORM_VTR + "_" + INT_NDRE_DIMRED_VTR + ".rda")
    params:
        sing_int_bind = int_norm_dimred_params_sing,
        pipeline_folder = os.path.normpath("/WORKDIR/" + PIPELINE_FOLDER),
        int_input_rda = int_norm_dimred_input_ge_R,
        int_output_folder = os.path.normpath("/WORKDIR/" + "{output_int_norm_dimred_dir_ge}") + "/",
        SING_INT_NDRE_METADATA_FILE = ','.join([os.path.normpath("/WORKDIR/" + x) for x in INT_NDRE_METADATA_FILE.split(',')]) if INT_NDRE_METADATA_FILE != "NULL" else "NULL",
        INT_SINGULARITY_ENV = INT_SINGULARITY_ENV,
        INT_NDRE_AUTHOR_NAME = INT_NDRE_AUTHOR_NAME,
        INT_NDRE_AUTHOR_MAIL = INT_NDRE_AUTHOR_MAIL,
        INT_NDRE_EVAL_MARKERS = INT_NDRE_EVAL_MARKERS,
        INT_NDRE_MIN_CELLS = INT_NDRE_MIN_CELLS,
        INT_NDRE_INT_METHOD = INT_NDRE_INT_METHOD,
        INT_NDRE_VTR_BATCH = INT_NDRE_VTR_BATCH,
        INT_NDRE_FEATURES_N = INT_NDRE_FEATURES_N,
        INT_NDRE_NORM_METHOD = INT_NDRE_NORM_METHOD,
        INT_NDRE_DIMRED_METHOD = INT_NDRE_DIMRED_METHOD,
        INT_NDRE_VTR_BIASES = INT_NDRE_VTR_BIASES,
        INT_NDRE_VTR_SCALE = INT_NDRE_VTR_SCALE,
        INT_NDRE_DIM_MAX = INT_NDRE_DIM_MAX,
        INT_NDRE_DIM_MIN = INT_NDRE_DIM_MIN,
        INT_NDRE_DIM_STEPS = INT_NDRE_DIM_STEPS,
        INT_NDRE_RES_MAX = INT_NDRE_RES_MAX,
        INT_NDRE_RES_MIN = INT_NDRE_RES_MIN,
        INT_NDRE_RES_STEPS = INT_NDRE_RES_STEPS
    threads:
        4
    resources:
        mem_mb = (lambda wildcards, attempt: min(5120 * len(dic_INT_NDRE_INFO[wildcards.name_int]['INT_NDRE_INPUT_LIST_RDA'].split(',')) + (attempt-1) * 10240, 184320)),
        time_min = (lambda wildcards, attempt: min(60 * len(dic_INT_NDRE_INFO[wildcards.name_int]['INT_NDRE_INPUT_LIST_RDA'].split(',')) + attempt * 120, 4320))
    shell:
        """
        export TMPDIR=$TMPDIR
        TMP_DIR=$(mktemp -d -t sc_pipeline-XXXXXXXXXX) && \
        singularity exec --no-home -B $TMP_DIR:/tmp {params.sing_int_bind} \
        {params.INT_SINGULARITY_ENV} \
        Rscript {params.pipeline_folder}/scripts/Integration_part1.R \
        --input.list.rda {params.int_input_rda} \
        --output.dir.int {params.int_output_folder} \
        --name.int {wildcards.name_int} \
        --author.name {params.INT_NDRE_AUTHOR_NAME} \
        --author.mail {params.INT_NDRE_AUTHOR_MAIL} \
        --nthreads {threads} \
        --pipeline.path {params.pipeline_folder} \
        --eval.markers {params.INT_NDRE_EVAL_MARKERS} \
        --min.cells {params.INT_NDRE_MIN_CELLS} \
        --integration.method {params.INT_NDRE_INT_METHOD}  \
        --vtr.batch {params.INT_NDRE_VTR_BATCH} \
        --features.n {params.INT_NDRE_FEATURES_N} \
        --norm.method {params.INT_NDRE_NORM_METHOD} \
        --dimred.method {params.INT_NDRE_DIMRED_METHOD} \
        --vtr.biases {params.INT_NDRE_VTR_BIASES} \
        --vtr.scale {params.INT_NDRE_VTR_SCALE} \
        --dims.max {params.INT_NDRE_DIM_MAX} \
        --dims.min {params.INT_NDRE_DIM_MIN} \
        --dims.steps {params.INT_NDRE_DIM_STEPS} \
        --res.max {params.INT_NDRE_RES_MAX} \
        --res.min {params.INT_NDRE_RES_MIN} \
        --res.steps {params.INT_NDRE_RES_STEPS} \
        --metadata.file {params.SING_INT_NDRE_METADATA_FILE} && \
        rm -r $TMP_DIR || rm -r $TMP_DIR
        """
