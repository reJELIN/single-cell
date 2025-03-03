##########################################################################
"""
rule adding isoform matrix and results from isoform switching in seurat 
object
"""
##########################################################################

def input_rda_isoform_marker_ge(wildcards):
    return ISOFORM_DICT[wildcards.sample_name_ge]["input_rda"]

def ouput_folder_isoform_marker_ge(wildcards):
    return ISOFORM_DICT[wildcards.sample_name_ge]["output_folder"]

def isoform_mtx_isoform_marker_ge(wildcards):
    return ISOFORM_DICT[wildcards.sample_name_ge]["isoform_mtx"]   


rule isoform_marker_ge:
    input:
	    input_rda_file = input_rda_isoform_marker_ge
    output:
	    output_rda_file = ISFORM_MTX_INTPUT[0]+"{sample_name_ge}/{sample_name_ge}/{sample_name_ge}_"+CMA_CLUST_FOLDER+"_sobj_isoform.rda"
    params:
	    sing_bind = SING_PARAMS_OUTPUT,
	    folder_dir_ge = ouput_folder_isoform_marker_ge,
	    isoform_mtx = isoform_mtx_isoform_marker_ge,
	    gtf_ge = GTF_GE,
	    workflow_dir = PIPELINE_FOLDER,
	    species = ISOFORM_SPECIES
    threads:
	    1
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 10240, 30720)),
        time_min = (lambda wildcards, attempt: min(attempt * 30, 200))
    shell:
        """
        export TMPDIR={GLOBAL_TMP}
        export HOME={HOME}
        TMP_DIR=$(mktemp -d -t sc_pipeline-XXXXXXXXXX) &&
    	singularity exec --no-home -B $TMP_DIR:/tmp -B $TMP_DIR:$HOME -B {params.workflow_dir}/common/database:{params.workflow_dir}/common/database/ {params.sing_bind} {SINGULARITY_ENV_LR} \
    	Rscript {params.workflow_dir}/scripts/LR/pipeline_isoform.R \
    	--input_rda_ge {input.input_rda_file} \
    	--output_rda_ge {output.output_rda_file} \
    	--output_dir_ge {params.folder_dir_ge} \
    	--gtf_file {params.gtf_ge} \
    	--path_isoform_transcript_matrix {params.isoform_mtx} \
    	--species {params.species} && 
    	rm -r $TMP_DIR
        """