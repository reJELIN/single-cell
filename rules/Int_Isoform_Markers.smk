##########################################################################
"""
rule adding isoform matrix and results from isoform switching in seurat 
object
"""
##########################################################################

def input_int_isoform_markers(wildcards):
    return INT_ISOFORM_INPUT_RDA[wildcards.sample_int]

def output_int_isoform_markers(wildcards):
    return INT_ISOFORM_OUTPUT[wildcards.sample_int]
    
rule int_isoform_markers:
    input: 
        input_rda=input_int_isoform_markers
    output:
        done_out=INT_ISOFORM_OUPUT_DIR_GLOBAL+"GROUPED_ANALYSIS/INTEGRATED/{sample_int}/tmp/DONE_int_isoform_markers_"+INT_CMA_CLUST_FOLDER+".txt"
    params:
        sing_bind=SING_PARAMS_OUTPUT_INT,
        ouput_rda=output_int_isoform_markers,
        ouput_dir=INT_ISOFORM_OUPUT_DIR_GLOBAL,
        gtf_int=GTF_INT,
        workflow_dir=PIPELINE_FOLDER,
    resources:
        mem_mb=(lambda wildcards, attempt: min(attempt * 25600, 51200)),
        time_min=(lambda wildcards, attempt: min(attempt * 60, 200))
    shell:
        """
        export TMPDIR={GLOBAL_TMP}
        export HOME={HOME}
        TMP_DIR=$(mktemp -d -t sc_pipeline-XXXXXXXXXX) &&
    	singularity exec --no-home -B $TMP_DIR:/tmp -B $TMP_DIR:$HOME -B {params.workflow_dir}/common/database:{params.workflow_dir}/common/database/ {params.sing_bind} {SINGULARITY_ENV_LR} \
    	Rscript {params.workflow_dir}/scripts/LR/pipeline_integration_isoform.R \
    	--input_int_rda {input.input_rda} \
    	--output_rda_int {params.ouput_rda} \
    	--gtf_file {params.gtf_int} \
    	--output_dir {params.ouput_dir} &&
    	rm -r $TMP_DIR &&
    	touch {output.done_out}
        """

