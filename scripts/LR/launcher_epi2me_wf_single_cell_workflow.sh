#!/bin/bash

while getopts ":s:o:r:g:c:m:f:n:" option; do
    case $option in
        s) 
            sample_sheet=$OPTARG;;
        o)
            output_dir_path=$OPTARG;;
        r)
            ref_genome_dir=$OPTARG;;
        g)
            min_genes=$OPTARG;;
        c)
            min_cells=$OPTARG;;
        m)
            max_mito=$OPTARG;;
        f)
            fastq_dir_path=$OPTARG;;
        n)
            sample_name=$OPTARG;;
        t)
            tmpdir=$OPTARG;;

   esac
done


module load java/21.0.6-jdk
module load singularity-ce/4.2.2

export TMPDIR=${tmpdir}

echo "TMP DIR: ${TMPDIR}"

#parameters
pipeline_version="v2.0.3"
path_to_pipeline="/mnt/beegfs02/pipelines/wf-scLongReads_Nanopore/"${pipeline_version}"/wf-single-cell/"

mkdir -p ${output_dir_path}

#launch
export NXF_SINGULARITY_CACHEDIR="/mnt/beegfs02/pipelines/wf-scLongReads_Nanopore/"${pipeline_version}"/sing_img/"
nextflow run ${path_to_pipeline} \
    -w ${output_dir_path}/workspace \
    -profile singularity \
    --fastq ${fastq_dir_path} \
    --threads 20 \
    --single_cell_sample_sheet ${sample_sheet} \
    --ref_genome_dir ${ref_genome_dir} \
    --matrix_min_genes ${min_genes} \
    --matrix_min_cells ${min_cells} \
    --matrix_max_mito ${max_mito} \
    --out_dir ${output_dir_path} \
    -resume
