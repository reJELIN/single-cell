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
    max_mito=100
try:
    species=config["Alignment_countTable_LR_GE"]["species"]
except KeyError:
    species="human"
    
try:
    species=config["species"]
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
These rule launch the epi2melab wf single-cell workflowdata
#################################################################################################################
"""
rule alignment_inputs_ge_lr:
    input:
        ref_genome_dir=ref_genome_dir,
        samplesheet=ALIGN_OUTPUT_DIR_GE+"/samplesheet/{sample_name_ge}_samplesheet.csv"
    output:
        genes_counts_matrix=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}_GE/{sample_name_ge}_GE/gene_raw_feature_bc_matrix/matrix.mtx.gz",
        transcript_counts_matrix=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}_GE/{sample_name_ge}_GE/transcript_raw_feature_bc_matrix/matrix.mtx.gz",
        bam=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}_GE/{sample_name_ge}_GE/tagged.bam"
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 50000, 600000)),
        time_min = (lambda wildcards, attempt: min(attempt * 2880, 10080)),
    envmodules:
        "java/12.0.2",
        "nextflow/21.10.6",
        "singularity/3.6.3",
    threads:
        20
    params:
        output_path=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}_GE/",
        workflow_dir=PIPELINE_FOLDER,
        sample_id="{sample_name_ge}_GE",
        min_cells=min_cells,
        min_genes=min_genes,
        max_mito=max_mito,
        fastq_path=lambda wcs: FASTQ_PATH_GE[wcs.sample_name_ge]
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
        -f {params.fastq_path} \
        -t $TMPDIR
        """
    
#rule epi2me_pipeline:
#    input:
#        ref_genome=ref_genome_dir=ref_genome_dir,
#        samplesheet=ALIGN_OUTPUT_DIR_GE+"/samplesheet/{sample_name_ge}_samplesheet.csv"
#        pipeline_dir="",
#        fastq_dir="",
#    output:
#        file1="",
#        file2="",
#    handover: True
#    envmodules:
#        "java/12.0.2",
#        "nextflow/21.10.6",
#        "singularity/3.6.3",
#    resources:
#        mem_mb="",
#        runtime="",
#        tmpdir="./tmp/",
#    threads: 20
#    params:
#        out_dir=lambda wildcards, output: os.path.commonprefix(output),
#    log:
#        "logs/..."
#    benchmark:
#        "benchmark/..."
#    shell:
#        "nexflot run "
#        "--threads {threads} "
#        "..."
#        "> {log} 2>&1 "