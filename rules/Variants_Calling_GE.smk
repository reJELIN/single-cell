wildcard_constraints:
    sample_name_ge=""

if species_longshot == "human":
    ref_genome_dir="/mnt/beegfs02/database/bioinfo/cellranger/2020-A/refdata-gex-GRCh38-2020-A/fasta/genome.fa"
if species_longshot == "mouse":
    ref_genome_dir="/mnt/beegfs02/database/bioinfo/cellranger/2020-A/refdata-gex-mm10-2020-A/fasta/genome.fa"

rule filter_bam:
    input:
        main_bam=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/tagged.bam",
        txt_barcodes=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/"+FILTERS_FOLDER +"/DOUBLETSFILTER_all/{sample_name_ge}_GE_FILTERED_barcodes.txt"
    output:
        filtered_bam=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/filtered_tagged.bam"
    resources:
        mem_mb = lambda wildcards, attempt: min(attempt * 128, 2048),
        time_min = lambda wildcards, attempt: min(attempt * 120, 360)
    conda:
        "/mnt/beegfs02/pipelines/single-cell/lr_1.3_test/single-cell/envs/conda/5009276213d3fd3f1bcae2865c827914_.yaml"
    threads:
        3
    shell:
        """
        samtools view -@ 3 -D CB:{input.txt_barcodes} {input.main_bam} -h -b -o {output.filtered_bam}
        """
        
rule sort_bam:
    input:
        filtered_bam=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/filtered_tagged.bam"
    output:
        sort_bam=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/filtered_sort_tagged.bam"
    resources:
        mem_mb = lambda wildcards, attempt: min(attempt *5120 , 10240),
        time_min = lambda wildcards, attempt: min(attempt * 120, 360)
    conda:
        "/mnt/beegfs02/pipelines/single-cell/lr_1.3_test/single-cell/envs/conda/5009276213d3fd3f1bcae2865c827914_.yaml"
    threads:
        3
    shell:
        """
        samtools sort -@ 3 {input.filtered_bam} > {output.sort_bam}
        samtools index {output.sort_bam}
        """

rule variant_calling_longshot:
    input:
        filtered_bam=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/filtered_sort_tagged.bam"
    output:
        vcf_longshot_output=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/{sample_name_ge}_longshot.vcf"
    resources:
        mem_mb = lambda wildcards, attempt: min(attempt *10240 , 30720),
        time_min = lambda wildcards, attempt: min(attempt * 360, 2880)
    params:
        input_fa=ref_genome_dir
    conda:
        "/mnt/beegfs02/pipelines/single-cell/lr_1.3_test/single-cell/envs/conda/env_longshot.yml"
    threads:
       1
    shell:
        """
        longshot --bam {input.filtered_bam} --ref {params.input_fa} --out {output.vcf_longshot_output} -y 20 -c 6 -C 10000 --strand_bias_pvalue_cutoff 0.01 --force_overwrite --auto_max_cov
        """

rule compress_index_vcf:
    input:
        vcf_longshot_input=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/{sample_name_ge}_longshot.vcf"
    output:
        vcf_longshot_output_gz=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/{sample_name_ge}_longshot.vcf.gz"
    resources:
        mem_mb = lambda wildcards, attempt: min(attempt *1024 , 3072),
        time_min = lambda wildcards, attempt: min(attempt * 5, 60)
    threads:
        3
    conda:
        "/mnt/beegfs02/pipelines/unofficial-snakemake-wrappers/shared_install/d58a6ebdc411f7854aad53828a4f5722_.yaml"
    shell:
        """
        bcftools view --write-index=tbi --with-header -O z -o {output.vcf_longshot_output_gz} {input.vcf_longshot_input}
        """

rule filering_vcf:
    input:
        vcf_longshot_input_gz=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/{sample_name_ge}_longshot.vcf.gz"
    output:
        vcf_longshot_filtered_output=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/{sample_name_ge}_longshot_filtered.vcf.gz"
    resources:
        mem_mb = lambda wildcards, attempt: min(attempt *1024 , 3072),
        time_min = lambda wildcards, attempt: min(attempt * 5, 60
    threads:
        3
    shell:
        """
        module load vcftools/0.1.16 
        vcftools --gzvcf {input.vcf_longshot_input_gz} --minQ 30 --minDP 10 --recode --stdout | bgzip -c > {output.vcf_longshot_filtered_output}#--remove-indels
        """

rule index_filtered_vcf:
    input:
        vcf_longshot_filtered_input=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/{sample_name_ge}_longshot_filtered.vcf.gz"
    output:
        index_vcf_longshot_filtered_output=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/{sample_name_ge}_longshot_filtered.vcf.gz.tbi"
    resources:
        mem_mb = lambda wildcards, attempt: min(attempt *1024 , 3072),
        time_min = lambda wildcards, attempt: min(attempt * 5, 60
    threads:
        3
    conda:
        "/mnt/beegfs02/pipelines/single-cell/lr_1.3_test/single-cell/envs/conda/d58a6ebdc411f7854aad53828a4f5722_.yaml"
    shell:
        """
        bcftools index -t {input.vcf_longshot_filtered_input} -o {index_vcf_longshot_filtered_output}
        """