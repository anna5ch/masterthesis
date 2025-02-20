rule starindex:
    input:
        genome=config["genome_fasta"],
        gtf=config["reference_gtf"],
    output:
        output_dir=directory(config["star_index_dir"]),
        STAR_index=config["star_index_dir"] + "SA",
    log:
        "../logs/align/STAR/STAR_index.log",
    params:
        readlength=config["readlength"],
    conda:
        "../envs/star.yaml"
    threads: config["STAR_ncores"]
    shell:
        "echo 'STAR version:\n' > {log}; STAR --version >> {log}; "
        "STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {output.output_dir} --genomeFastaFiles {input.genome} --sjdbGTFfile {input.gtf} --sjdbOverhang {params.readlength}"


rule run_star:
    input:
        index=config["star_index_dir"] + "SA",
        fastq1=config["raw_data_path"]
        + "{sample}-r1."
        + config["dataset"]
        + ".fastq.gz",
        fastq2=config["raw_data_path"]
        + "{sample}-r2."
        + config["dataset"]
        + ".fastq.gz",
    output:
        bam_out=config["datapath_short"]
        + "{sample}."
        + config["dataset"]
        + ".sorted.bam",
    threads: config["STAR_ncores"]
    log:
        "logs/align/STAR/STAR_{sample}.log",
    params:
        tmp_out=config["datapath_short"] + "{sample}/Aligned.sortedByCoord.out.bam",
        out_dir=directory(config["datapath_short"] + "{sample}/"),
        star_dir=directory(config["star_index_dir"]),
        starextraparams=config["additional_star_align"],
    conda:
        "../envs/star.yaml"
    shell:
        """
        echo 'STAR version:\n' > {log}; STAR --version >> {log}; 
        if file {input.fastq1} | grep -Eq 'gzip compressed|gzip'; then
            read_command="--readFilesCommand gunzip -c"
        else
            read_command=""
        fi
        
        STAR --genomeDir {params.star_dir} --readFilesIn {input.fastq1} {input.fastq2} --runThreadN {threads} --outFileNamePrefix {params.out_dir} --outSAMtype BAM SortedByCoordinate {params.starextraparams} $read_command
        mv {params.tmp_out} {output.bam_out}
        """


rule bamindex:
    input:
        bam=config["datapath_short"] + "{sample}." + config["dataset"] + ".sorted.bam",
    output:
        config["datapath_short"] + "{sample}." + config["dataset"] + ".sorted.bam.bai",
    log:
        "../logs/align/samtools_index_{sample}.log",
    conda:
        "../envs/star.yaml"
    shell:
        "echo 'samtools version:\n' > {log}; samtools --version >> {log}; "
        "samtools index {input.bam}"


rule run_minimap2:
    input:
        target=config["genome_fasta"],
        query="/home/dwissel/data/seq/kinnex/raw/{sample}/"
        + config["dataset"]
        + "reads.fastq.gz",
    output:
        "/home/anna/seq_data/seq/kinnex/sirvs_full_aligned_new/{sample}/{sample}."
        + config["dataset"]
        + ".sorted.bam",
    params:
        extra="-ax splice:hq -uf",
        sorting="coordinate",
        sort_extra=f"-m{config['quantify_sort_bam_memory_gb']}g",
    log:
        "logs/align/minimap2/{sample}.stderr",
    threads: 12
    wrapper:
        f"{config['snakemake_wrapper_version']}/bio/minimap2/aligner"


rule make_blastdb:
    input:
        config["blast_fasta"],
    output:
        blast_db="GCA_genomedb",
    shell:
        "makeblastdb -in {input} -dbtype nucl -out GCA_genomed"


rule run_blast:
    input:
        sample_bam="{sample}.bam",
        db="GCA_genomedb",
    output:
        "blast_file.txt",
    shell:
        """
        samtools fasta {input.sample_bam} > {input.sample_bam.replace('.bam', '.fasta')}
        blastn -query day0-rep1.sirvs.fasta -db {input.db} -out {output} -outfmt 6 
        """


rule fastq_to_unaligned_bam:
    input:
        config["raw_data_path"] + "{sample}/" + config["dataset"] + ".reads.fastq.gz",
    output:
        config["path_unmapped"] + "{sample}.bam",
    conda:
        "../envs/base.yaml"
    shell:
        """
        echo "samtools import -0 {input}  -o {output}"
        samtools import -0 {input} -o {output}
        """
