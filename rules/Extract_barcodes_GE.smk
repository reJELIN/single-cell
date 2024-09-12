rule extract_barcodes_from_rda:
    input:
        rda_file=ALIGN_OUTPUT_DIR_GE+ "/{sample_name_ge}/" + FILTERS_FOLDER + "/DOUBLETSFILTER_all/{sample_name_ge}_FILTERED_NON-NORMALIZED.rda"
    output:
        txt_barcodes=ALIGN_OUTPUT_DIR_GE+ "/{sample_name_ge}/" + FILTERS_FOLDER + "/DOUBLETSFILTER_all/{sample_name_ge}_FILTERED_barcodes.txt"
    params:
        workflow_dir=PIPELINE_FOLDER,
        directory_rda=ALIGN_OUTPUT_DIR_GE+ "/{sample_name_ge}/" + FILTERS_FOLDER + "/DOUBLETSFILTER_all/"
    threads:
        2
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 1024, 2048)),
        time_min = lambda wildcards, attempt: min(attempt * 10, 30)
    shell:
        """
        module load singularity
        singularity exec --no-home -B {params.directory_rda}:{params.directory_rda} -B {params.workflow_dir}:{params.workflow_dir} {params.workflow_dir}/envs/singularity/single_cell.simg Rscript \
        {params.workflow_dir}/scripts/LR/get_barcodes.R \
        --input_rda_ge {input.rda_file} \
        --output_barcodes_file_txt_ge {output.txt_barcodes}
        """