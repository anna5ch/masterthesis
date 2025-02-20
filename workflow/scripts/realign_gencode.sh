#!/bin/bash

# Script to perform genome mapping using minimap2 and samtools for gencode files

THREADS=12 
GENOME="/home/dwissel/data/seq/anna/ref/SIRV4.fasta"
FASTQ_DIR="/home/dwissel/data/seq/kinnex/raw"
OUTPUT_DIR="/home/anna/seq_data/seq/kinnex/sirvs_full_aligned_new"

# create output dir
mkdir -p ${OUTPUT_DIR}

for sample_dir in ${FASTQ_DIR}/day*; do
    #sample=$(basename ${sample_dir} | cut -d'-' -f1,2)
    sample=$(basename ${sample_dir})
    FASTQ="${FASTQ_DIR}/${sample}/sirvs.reads.fastq.gz"
    #FASTQ_R1="${FASTQ_DIR}/${sample}-r1.gencode.fastq.gz"
    #FASTQ_R2="${FASTQ_DIR}/${sample}-r2.gencode.fastq.gz"
    
    # Check if the gencode file exists
    if [ -f "${FASTQ}" ]; then
    #if [ -f "${FASTQ_R1}" ] && [ -f "${FASTQ_R2}" ]; then   
        OUTPUT_BAM="${OUTPUT_DIR}/${sample}/${sample}.sirvs.sorted.bam"

        mkdir -p ${OUTPUT_DIR}/${sample}
        # Run minimap2 and samtools
        # LONG READS:
        if [ -f "${OUTPUT_BAM}" ]; then
            echo "Output BAM file already exists for ${sample}, skipping..."
            continue
        fi

        #echo $OUTPUT_BAM
        minimap2 -ax splice:hq -uf -t ${THREADS} ${GENOME} ${FASTQ} | \
        samtools sort -@ 4 -m4g -o ${OUTPUT_BAM} --write-index -
        
        # SHORT READS:
        #minimap2 -ax sr -t ${THREADS} ${GENOME} ${FASTQ_R1} ${FASTQ_R2} | \
        #samtools sort -@ 4 -m4g -o ${OUTPUT_BAM} --write-index - 


    else
        echo "No gencode.fastq.gz found for ${sample}, skipping..."
    fi
done
