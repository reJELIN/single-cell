##########################################################################
"""
rule adding isoform matrix and results from isoform switching in seurat 
object
"""
##########################################################################

rule isoform_marker_ge:
    input:
	    input_rda_file = expand("{start}/{sample_name_ge}{end}",zip,start=COMPLEMENT_ISOFORM_INPUT_RDA_LIST_START,sample_name_ge=ISOFORM_SAMPLE_NAME_GE,end=COMPLEMENT_ISOFORM_INPUT_RDA_LIST_END)
    output:
	    output_rda_file = expand("{start}/{sample_name_ge}{end}",zip,start=COMPLEMENT_ISOFORM_OUTPUT_RDA_LIST_START,sample_name_ge=ISOFORM_SAMPLE_NAME_GE,end=COMPLEMENT_ISOFORM_OUTPUT_RDA_LIST_END)
    params:
	    sing_bind = expand("{sing_params}",sing_params=SING_PARAMS_OUTPUT),
	    folder_dir_ge = expand("{folder_output}",folder_output=ISOFORM_OUTPUT_RDA),
	    isoform_mtx = expand("{folder_path}{sample_name_ge}/{sample_name_ge}/transcript_processed_feature_bc_matrix/",zip,folder_path=ISFORM_MTX_INTPUT,sample_name_ge=ISOFORM_SAMPLE_NAME_GE),
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