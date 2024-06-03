try:
    min_cells=config["Alignment_countTable_LR_GE"]["min.cells"]
except KeyError:
    min_cells=1
try:
    min_genes=config["Alignment_countTable_LR_GE"]["min.genes"]
except KeyError:
    min_genes=1
try:
    max_mito=config["Alignment_countTable_LR_GE"]["max.mito"]
except KeyError:
    max_mito=1
try:
    species=config["Alignment_countTable_LR_GE"]["species"]
except KeyError:
    species="human"
    
if species == "human":
    ref_genome_dir="/mnt/beegfs/database/bioinfo/cellranger/2020-A/refdata-gex-GRCh38-2020-A"
    
if species == "mouse":
    ref_genome_dir="/mnt/beegfs/database/bioinfo/cellranger/2020-A/refdata-gex-mm10-2020-A"
    
"""
#################################################################################################################
These rule create the sample sheet in csv format that is mandatory in order to run epi2melab wf single-cell
#################################################################################################################
"""
    
rule create_sample_sheet:
    input:
        design_file_path=DESIGN_FILE_GE
    output:
        output_samplesheet=ALIGN_OUTPUT_DIR_GE+"/samplesheet/{sample_name_ge}_samplesheet.csv"
    params:
        output_path=ALIGN_OUTPUT_DIR_GE,
        workflow_dir=PIPELINE_FOLDER,
        sample_id="{sample_name_ge}"
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 256, 2048)),
        time_min = (lambda wildcards, attempt: min(attempt * 1, 5))
    shell:
        """
        python3 \
        {params.workflow_dir}/scripts/LR/create_samplesheet.py \
        --sample_id {params.sample_id} \
        --output_path {params.output_path} \
        --design_file_path {input.design_file_path}
        """

"""
#################################################################################################################
These rule launch the epi2melab wf single-cell workflow
#################################################################################################################
"""

rule alignment_inputs_ge_lr:
    input:
        ref_genome_dir=ref_genome_dir,
        samplesheet=ALIGN_OUTPUT_DIR_GE+"/samplesheet/{sample_name_ge}_samplesheet.csv"
    output:
        genes_counts_matrix=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/{sample_name_ge}.gene_expression.counts.tsv",
        transcript_counts_matrix=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/{sample_name_ge}.transcript_expression.counts.tsv"
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 300000, 356000)),
        time_min = (lambda wildcards, attempt: min(attempt * 2880, 5760))
    threads:
        48
    params:
        output_path=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/",
        workflow_dir=PIPELINE_FOLDER,
        sample_id="{sample_name_ge}",
        min_cells=min_cells,
        min_genes=min_genes,
        max_mito=max_mito,
        fastq_path=lambda wcs: FASTQ_PATH_GE[wcs.sample_name_ge],
    shell:
        """
        export TMPDIR={GLOBAL_TMP}
        bash {params.workflow_dir}/scripts/LR/launcher_epi2me_wf_single_cell_workflow.sh \
        -s {input.samplesheet} \
        -g {params.min_genes} \
        -c {params.min_cells} \
        -m {params.max_mito} \
        -r {input.ref_genome_dir} \
        -o {params.output_path} \
        -n {params.sample_id} \
        -f {params.fastq_path}
        """

"""
#################################################################################################################
These rule compress and rename the gene expression counts matrix in order to make it work with the srsc pipeline
#################################################################################################################
"""

rule compress_and_rename:
    input:
        genes_counts_matrix=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/{sample_name_ge}.gene_expression.counts.tsv"
    output:
        umi_tools=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/{sample_name_ge}_GE.tsv.gz"
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 256, 2048)),
        time_min = (lambda wildcards, attempt: min(attempt * 1, 5))
    threads:
        1
    params:
        output_path=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/"
    shell:
        """
        bzip2 -kz {input.genes_counts_matrix}
        mv {input.genes_counts_matrix}.bz2 {params.output_path}/{sample_name_ge}/{sample_name_ge}_GE.tsv.gz
        """