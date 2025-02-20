rule quantify_gffread_extract_tool_transcriptome:
    input:
        annotation="results/collapse/{tool}/{sample}/{sample}.{tool}_output_collapsed.gtf",
    output:
        "results/quantify/{tool}/{sample}_transcriptome.fasta",
    params:
        fasta=config["genome_fasta"],
    threads: 1
    log:
        "logs/quantify/gffread_extract/{tool}_{sample}.stderr",
    conda:
        "../envs/gffread.yaml"
    shell:
        "gffread -w {output} -g {params.fasta} {input.annotation} 2> {log}"


rule quantify_minimap2_map_novel_transcriptome:
    input:
        target="results/quantify/{tool}/{sample}_transcriptome.fasta",
        query=config["raw_data_path"] + "{sample}/" + config["dataset"] + ".fastq.gz",
    output:
        "results/quantify/minimap2/transcriptome/{tool}/{sample}.aligned.sorted.bam",
    params:
        extra="-x map-hifi --for-only --sam-hit-only",
        sorting="coordinate",
        sort_extra=f"-m{config['quantify_sort_bam_memory_gb']}g",
    threads: config["quantify_map_bam_threads"] + 2
    log:
        "logs/quantify/minimap2/{tool}_{sample}.stderr",
    wrapper:
        f"{config['snakemake_wrapper_version']}/bio/minimap2/aligner"


rule quantify_run_oarfish:
    input:
        bams="results/quantify/minimap2/transcriptome/{tool}/{sample}.aligned.sorted.bam",
    output:
        "results/quantify/{tool}/oarfish_{tool}_{sample}/oarfish_{tool}_{sample}.quant",
    params:
        output_path="results/quantify/{tool}/oarfish_{tool}_{sample}/oarfish_{tool}_{sample}",
        n_bins=config["quantify_oarfish_n_bins"],
        filter_group=config["quantify_oarfish_filter_group"],
    threads: config["quantify_threads"]
    log:
        "logs/quantify/oarfish/{tool}_{sample}.stderr",
    conda:
        "../envs/oarfish.yaml"
    shell:
        """
        oarfish --threads {threads} --filter-group {params.filter_group} --model-coverage --bins {params.n_bins} --alignments {input.bams} --output {params.output_path} 2> {log}
        """


rule quantify_summarize_oarfish:
    input:
        input_path="results/quantify/{tool}/oarfish_{tool}_{sample}/oarfish_{tool}_{sample}.quant",
    output:
        "results/quantify/{tool}/oarfish_{tool}_{sample}/oarfish_{sample}_transcript_counts.tsv",
    params:
        sample_name="{sample}",
        gtf_annotation_path=config["reference_gtf"],
    threads: 1
    log:
        "logs/quantify/oarfish/summarize_{tool}_{sample}.stderr",
    conda:
        "../envs/oarfish.yaml"
    shell:
        """python scripts/aggregate_oarfish.py \
                    --sample_paths {input.input_path} \
                    --sample_names {params.sample_name} \
                    --gtf_annotation_path {params.gtf_annotation_path} \
                    --output_path {output} 2> {log}"""


rule quantify_run_salmon:
    input:
        read1=config["raw_data_path_short"]
        + "{sample}-r1."
        + config["dataset"]
        + ".fastq.gz",
        read2=config["raw_data_path_short"]
        + "{sample}-r2."
        + config["dataset"]
        + ".fastq.gz",
        transcriptome_fasta="results/quantify/{tool}/{sample}_transcriptome.fasta",
    output:
        quant="results/quantify/{tool}/salmon_{tool}_{sample}/salmon_{tool}_{sample}.quant",
    params:
        index_dir="results/quantify/salmon/index",
        libtype="A",
        sample_name="{sample}",
        quant="results/quantify/{tool}/salmon_{tool}_{sample}/quant.sf",
    threads: config["quantify_threads"]
    log:
        "logs/quantify/salmon/{tool}_{sample}.stderr",
    conda:
        "../envs/salmon.yaml"
    shell:
        """
        salmon index -t {input.transcriptome_fasta} -i {params.index_dir} -p {threads}
        salmon quant -i {params.index_dir} -l {params.libtype} -1 {input.read1} -2 {input.read2} -p {threads} -o results/quantify/{wildcards.tool}/salmon_{wildcards.tool}_{wildcards.sample} 2> {log}
        python3 scripts/convert_quant.py {params.quant} {output.quant}
        """


rule quantify_summarize_salmon:
    input:
        input_path="results/quantify/{tool}/salmon_{tool}_{sample}/salmon_{tool}_{sample}.quant",
        gtf_annotation_path=config["reference_gtf"],
    output:
        "results/quantify/{tool}/salmon_{tool}_{sample}/salmon_{sample}_transcript_counts.tsv",
    params:
        sample_name="{sample}",
    threads: 1
    log:
        "logs/quantify/salmon/summarize_{tool}_{sample}.stderr",
    conda:
        "../envs/oarfish.yaml"
    shell:
        """python scripts/aggregate_oarfish.py \
                    --sample_paths {input.input_path} \
                    --sample_names {params.sample_name} \
                    --gtf_annotation_path {input.gtf_annotation_path} \
                    --output_path {output} 2> {log}"""
