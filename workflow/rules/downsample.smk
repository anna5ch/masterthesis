
rule all:
    input:
        expand(
            config["raw_data_path"] + "{sample}/gencode.fastq.gz",
            sample=[
                "day0-rep1",
                "day0-rep2",
                "day0-rep3",
                "day5-rep1",
                "day5-rep2",
                "day5-rep3",
            ],
        ),
        expand(
            config["raw_data_path_short"] + "{sample}-r1.gencode.fastq.gz",
            sample=[
                "day0-rep1",
                "day0-rep2",
                "day0-rep3",
                "day5-rep1",
                "day5-rep2",
                "day5-rep3",
            ],
        ),
        expand(
            config["raw_data_path_short"] + "{sample}-r2.gencode.fastq.gz",
            sample=[
                "day0-rep1",
                "day0-rep2",
                "day0-rep3",
                "day5-rep1",
                "day5-rep2",
                "day5-rep3",
            ],
        ),


rule subsample_bam:
    input:
        input_name=config["datapath"] + "{sample}.aligned.sorted.bam",
    output:
        output_name=config["datapath"]
        + "downsampled/{sample}."
        + config["dataset"]
        + ".sorted.bam",
    log:
        "logs/align/downsample/downsample_{sample}.log",
    params:
        number_to_sample=5000000,
        seed=42,
        read_name_sam="qname",
        is_transcriptome=False,
    script:
        "../scripts/subsample_bam.R"


rule downsample_fastq:
    input:
        raw_fastq="/home/dwissel/data/seq/kinnex/raw/{sample}/gencode.reads.fastq.gz",
        bam="/home/anna/seq_data/seq/kinnex/gencode_aligned_new/downsampled/{sample}.gencode.sorted.bam",
    output:
        filtered_fastq=config["raw_data_path"] + "{sample}/gencode.fastq.gz",
    params:
        "/home/anna/seq_data/seq/kinnex/gencode_aligned_new/downsampled/{sample}_name.txt",
    threads: 1
    shell:
        """
        samtools view {input.bam} | cut -f 1 | sort | uniq  > {params}
        bash ../scripts/bbmap/filterbyname.sh in={input.raw_fastq} out={output.filtered_fastq} names={params} include=t
        """


rule downsample_fastq_short:
    input:
        raw_fastq1="/home/dwissel/data/seq/kinnex/raw/illumina/FASTQ/{sample}-r1.gencode.fastq.gz",
        raw_fastq2="/home/dwissel/data/seq/kinnex/raw/illumina/FASTQ/{sample}-r2.gencode.fastq.gz",
        bam="/home/anna/seq_data/seq/illumina/gencode_aligned_new/downsampled/{sample}.gencode.sorted.bam",
    output:
        filtered_fastq1=config["raw_data_path_short"] + "{sample}-r1.gencode.fastq.gz",
        filtered_fastq2=config["raw_data_path_short"] + "{sample}-r2.gencode.fastq.gz",
    params:
        "/home/anna/seq_data/seq/illumina/gencode_aligned_new/downsampled/{sample}_name.txt",
    threads: 1
    shell:
        """
        samtools view {input.bam} | cut -f 1 | sort | uniq  > {params}
        bash ../scripts/bbmap/filterbyname.sh in={input.raw_fastq1} in2={input.raw_fastq2} out={output.filtered_fastq1} out2={output.filtered_fastq2} names={params} int=f include=t
        """
