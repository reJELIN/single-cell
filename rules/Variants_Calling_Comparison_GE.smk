def input_rda(wildcards):
    return design_dict[wildcards.sample_name_ge]["input_rda"]

def cluster_ctrl(wildcards):
    return design_dict[wildcards.sample_name_ge]["clusters_control"]

def cluster_ofint(wildcards):
    return design_dict[wildcards.sample_name_ge]["clusters_of_interest"]
    
if config["species"] == "human":
    ref_genome="/mnt/beegfs02/database/bioinfo/cellranger/2020-A/refdata-gex-GRCh38-2020-A/fasta/genome.fa"
if config["species"] == "mouse":
    ref_genome="/mnt/beegfs02/database/bioinfo/cellranger/2020-A/refdata-gex-mm10-2020-A/fasta/genome.fa"

rule extract_barcodes:
    input:
       rda_file=input_rda
    output:
        barcodes_ctrl=ouput_dir_ge+"{sample_name_ge}/{sample_name_ge}/SNV_COMPARE/barcodes_ctrl.txt",
        barcodes_ofint=ouput_dir_ge+"{sample_name_ge}/{sample_name_ge}/SNV_COMPARE/barcodes_to_test.txt"
    params:
        clusters_of_interest=cluster_ofint,
        clusters_control=cluster_ctrl,
        output_dir_path=ouput_dir_ge,
        workflow_dir=PIPELINE_FOLDER
    resources:
        mem_mb = lambda wildcards, attempt: min(attempt *5120 , 15360),
        time_min = lambda wildcards, attempt: min(attempt * 20, 60)
    threads:
        2
    shell:
        """
        module load singularity
            
        singularity exec --no-home -B {params.output_dir_path}:{params.output_dir_path} -B {params.workflow_dir}:{params.workflow_dir} {params.workflow_dir}/envs/singularity/single_cell.simg Rscript \
        {params.workflow_dir}/scripts/LR/get_barcodes_variants_compare.R \
        --input_rda_ge {input.rda_file} \
        --output_cluster_ctrl {output.barcodes_ctrl} \
        --output_cluster_ofint {output.barcodes_ofint} \
        --cluster_ctrl {params.clusters_control} \
        --cluster_ofint {params.clusters_of_interest}
        """
        

rule subset_bam:
    input:
        barcodes_ctrl=ouput_dir_ge+"{sample_name_ge}/{sample_name_ge}/SNV_COMPARE/barcodes_{statut}.txt"
        bam_in=ouput_dir_ge+"{sample_name_ge}/{sample_name_ge}/tagged.bam"
    output:
        bam_filtered=ouput_dir_ge+"{sample_name_ge}/{sample_name_ge}/SNV_COMPARE/bam/{statut}/filtered_{statut}.bam
    resources:
        mem_mb = lambda wildcards, attempt: min(attempt 1024*15 , 1024*45),
        time_min = lambda wildcards, attempt: min(attempt * 180, 360)
    threads:
        3
    shell:
        """
        module load samtools
        subset-bam \
        --bam {input.bam_in} \
        --cell-barcodes {input.barcodes_ctrl} \
        --bam-tag CB:Z --cores 3 --log-level debug \
        --out-bam {output.bam_filtered}
        
        samtools index -@ 3 -b {output.bam_filtered}
        """

rule run_longshot:
    input:
        ouput_dir_ge+"{sample_name_ge}/{sample_name_ge}/SNV_COMPARE/bam/{statut}/filtered_{statut}.bam
    output:
        ouput_dir_ge+"{sample_name_ge}/{sample_name_ge}/SNV_COMPARE/vcf/{statut}/filtered_{statut}.vcf
    resources:
        mem_mb = lambda wildcards, attempt: min(attempt 1024*15 , 1024*45),
        time_min = lambda wildcards, attempt: min(attempt * 360, 1440)
    params:
        ref_genome=ref_genome
    threads:
        2
    conda:
        PIPELINE_FOLDER+"envs/conda/env_longshot.yml"
    shell:
        """
        longshot --bam {input} \
        --ref {params.ref_genome} \
        --out {output} \
        -y 20 \
        -c 20 \
        -C 10000 \
        -e 50 \
        --force_overwrite
        """

rule compress_vcf:
    input:
        vcf_in=ouput_dir_ge+"{sample_name_ge}/{sample_name_ge}/SNV_COMPARE/vcf/{statut}/filtered_{statut}.vcf
    output:
        vcf_out=ouput_dir_ge+"{sample_name_ge}/{sample_name_ge}/SNV_COMPARE/vcf/{statut}/filtered_{statut}.vcf.gz
    resources:
        mem_mb = lambda wildcards, attempt: min(attempt 1024*2 , 1024*5),
        time_min = lambda wildcards, attempt: min(attempt * 60, 180)
    threads:
        3
    conda:
        PIPELINE_FOLDER+"/envs/conda/d58a6ebdc411f7854aad53828a4f5722_.yaml
    shell:
        """
        bcftools view --write-index=tbi --with-header -O z -o {output.vcf_out} {input.vcf_in}
        """
        
rule bcftools_isec:
    input:
        vcf_gz_to_test=ouput_dir_ge+"{sample_name_ge}/{sample_name_ge}/SNV_COMPARE/vcf/to_test/filtered_to_test.vcf.gz",
        vcf_gz_ctrl=ouput_dir_ge+"{sample_name_ge}/{sample_name_ge}/SNV_COMPARE/vcf/ctrl/filtered_ctrl.vcf.gz"
    output:
        readme=ouput_dir_ge+"{sample_name_ge}/{sample_name_ge}/SNV_COMPARE/isec/README.txt"
    resources:
        mem_mb = lambda wildcards, attempt: min(attempt 1024*2 , 1024*5),
        time_min = lambda wildcards, attempt: min(attempt * 60, 180)
    threads:
        2
    params:
        prefix=ouput_dir_ge+"{sample_name_ge}/{sample_name_ge}/SNV_COMPARE/isec/"
    conda:
        ""
    shell:
        """
        bcftools isec --prefix {params.prefix} {input.vcf_gz_to_test} {input.vcf_gz_ctrl}
        """