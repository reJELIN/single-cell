wildcard_constraints:
    sample_name_ge=".+_GE",
    batch_number_seq=""

#bcftools -q 0 -Q 0 -x -C 0
#ruleorder: split_bams_by_barcodes > extract_barcodes_from_rda > create_batch_list > split_bams_by_barcodes

if VARIANTS_SPECIES == "human":
    ref_genome="/mnt/beegfs02/database/bioinfo/cellranger/2020-A/refdata-gex-GRCh38-2020-A/fasta/genome.fa"
    
if VARIANTS_SPECIES == "mouse":
    ref_genome="/mnt/beegfs02/database/bioinfo/cellranger/2020-A/refdata-gex-mm10-2020-A/fasta/genome.fa"
    
    
#def filtering_rda_path_input_ge(wildcards):
#    return dict_FILTERING_RDA[wildcards.sample_name_ge]

rule split_bam_by_chromosome:
    input:
        main_bam=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/tagged.bam"
    output:
        chr1=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/split_chrom/tagged_sort_"+chr_list[0]+".bam"
    params:
        dir_split_chrom=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/split_chrom",
        chrom_list=chr_list
    threads:
        1

    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 512, 10240)),
        time_min = lambda wildcards, attempt: min(attempt * 360, 1440)
    conda:
        PIPELINE_FOLDER+"/envs/conda/5009276213d3fd3f1bcae2865c827914_.yaml"
    shell:
        """
        mkdir -pv {params.dir_split_chrom}
        for chr in {params.chrom_list}; do
            samtools view -b -o {params.dir_split_chrom}/tagged_$chr.bam {input.main_bam} $chr
            samtools sort {params.dir_split_chrom}/tagged_$chr.bam -o {params.dir_split_chrom}/tagged_sort_$chr.bam
            samtools index -b {params.dir_split_chrom}/tagged_sort_$chr.bam
        done;
        """

rule create_batch_list:
    input:
        barcodes_list=ALIGN_OUTPUT_DIR_GE+ "/{sample_name_ge}/" + FILTERS_FOLDER + "/DOUBLETSFILTER_all/{sample_name_ge}_FILTERED_barcodes.txt",
    output:
        txt_barcodes=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/barcodes_batch/{sample_name_ge}_barcodes_list_{batch_number}.txt"
    params:
        directory_barcodes_batch=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/",
        workflow_dir=PIPELINE_FOLDER,
        sample_id="{sample_name_ge}",
        batch_id="{batch_number}"
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 256, 512)),
        time_min = lambda wildcards, attempt: min(attempt * 10, 60)
    shell:
        """
        python3 {params.workflow_dir}/scripts/LR/split_barcodes.py \
        --barcodes_list {input.barcodes_list} \
        --output_dir {params.directory_barcodes_batch} \
        --sample_id {params.sample_id} \
        --batch_id {params.batch_id}
        """

rule split_bam_per_bc:
    input:
        tagged_bam=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/tagged.bam"
    output:
        end_of_run=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/tmp/end_of_split_bam_per_bc.DONE"
    params:
        datadir=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}"
    threads:
        4
    resources:
        mem_mb = lambda wildcards, attempt: min(attempt * 1024 * 15, 1024*45),
        time_min = lambda wildcards, attempt: min(attempt * 420, 1440)
    conda:
        PIPELINE_FOLDER+"/envs/conda/5009276213d3fd3f1bcae2865c827914_.yaml"
    shell:
        """
        mkdir -p {params.datadir}/tmp
        mkdir -p {params.datadir}/per_cell_bams

        samtools split --max-split 20000 \
        -d CB \
        --threads {threads} \
        -f '{params.datadir}/per_cell_bams/%!.bam' \
        --no-PG \
        {input.tagged_bam}

        samtools index -@ 3 -M {params.datadir}/per_cell_bams/*.bam

        touch {output.end_of_run}
        """

rule mpileup_by_barcodes:
    input:
        end_of_=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/tmp/end_of_split_bam_per_bc.DONE"
         batch=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/barcodes_batch/{sample_name_ge}_barcodes_list_{batch_number}.txt"
    output:
        barcodes_batch_done=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/tmp_mpileup/DONE_{sample_name_ge}_barcodes_list_{batch_number}.txt"
    params:
        input_dir=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/per_cell_bams/",
        output_dir=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/barcodes_pileup/",
        ref_genome=ref_genome,
        bed_file=BED_FILE,
        barcodes_batch_done_dir=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/tmp_mpileup/"
    threads:
        1
    resources:
        mem_mb = lambda wildcards, attempt: min(attempt * 512, 2048),
        runtime = lambda wildcards, attempt: min(attempt * 30, 60)
    conda:
        PIPELINE_FOLDER+"/envs/conda/bam_readcount.yaml"
    shell:
        """
        mkdir -p {params.output_dir}
        for i in `cat {input.batch}`;
        do
            echo $i
            bam-readcount \
            -l {params.bed_file} \
            -f {params.ref_genome} \
            {params.input_dir}$i.bam > {params.output_dir}$i.tsv
        done
        mkdir -p {params.barcodes_batch_done_dir}
        cp {input.batch} {output.barcodes_batch_done}
        """

#            bcftools mpileup {params.input_dir}$i/$i.merge_sort.bam \
#            -R {params.bed_file} \
#            -f {params.ref_genome} \
#            -Ov -o {params.output_dir}$i.vcf \
#            --indels-cns -B -Q1 --max-BQ 60 --delta-BQ 99 -F0.2 \
#            -o15 -e1 -h110 --del-bias 0.4 --indel-bias 0.7 -M40000 \
#            --poly-mqual --seqq-offset 130 --indel-size 80 -a INFO/AD,FORMAT/AD,FORMAT/DP4,FORMAT/DV -A
rule create_snp_matrix:
    input:
        expand(ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/tmp_mpileup/DONE_{sample_name_ge}_barcodes_list_{batch_number}.txt",batch_number=batch_number_seq,sample_name_ge=ALIGN_SAMPLE_NAME_GE)
    output:
        matrix=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/matrix_SNP_{sample_name_ge}.tsv"
    params:
        vcf_input_dir=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/barcodes_pileup/",
        workflow_dir=PIPELINE_FOLDER,
        sample_id="{sample_name_ge}",
        output_dir=ALIGN_OUTPUT_DIR_GE+"/{sample_name_ge}/{sample_name_ge}/"
    threads:
        1
    conda:
       PIPELINE_FOLDER+"/envs/conda/env_pyvcf.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: min(attempt * 1024, 5120),
        runtime = lambda wildcards, attempt: min(attempt * 60, 120)
    shell:
        """
        python3 {params.workflow_dir}/scripts/LR/make_matrix_from_vcf.py \
        --input_directory {params.vcf_input_dir} \
        --output_directory {params.output_dir} \
        --prefix {params.sample_id}
        """
