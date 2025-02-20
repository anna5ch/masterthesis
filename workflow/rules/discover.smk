rule remove_part_annotation_bambu:
    input:
        annotation_gtf=config["annotation_gtf_bambu"],
    output:
        modified_gtf="results/discovery/bambu/modified_annotation.gtf",
        kept_ids="results/discovery/kept_ids.csv",
    params:
        region=config["requested_region_removal"],
        percentage=config["requested_thresh_removal"],
        mid_output="results/discovery/bambu/modified_annotation_bambu.gtf",
    shell:
        """
        python3 scripts/modify_GTF.py {input.annotation_gtf} {params.percentage} {params.region} {params.mid_output} --path_kept_ids {output.kept_ids}
        python3 scripts/modify_GTF2.py {params.mid_output} {output.modified_gtf}
        #grep -Pv '\tgene\t' 
        """


rule remove_part_annotation_isoquant:
    input:
        annotation_gtf=config["annotation_gtf_isoquant"],
        kept_ids="results/discovery/kept_ids.csv",
    output:
        modified_gtf="results/discovery/isoquant/modified_annotation_isoquant.gtf",
    params:
        region=config["requested_region_removal"],
        percentage=config["requested_thresh_removal"],
    shell:
        "python3 scripts/modify_GTF.py {input.annotation_gtf} {params.percentage} {params.region} {output.modified_gtf} --keep_ids {input.kept_ids} "


rule remove_part_annotation_stringtie:
    input:
        annotation_gtf=config["reference_gtf"],
        kept_ids="results/discovery/kept_ids.csv",
    output:
        modified_gtf="results/discovery/stringtie/modified_annotation.gtf",
    params:
        region=config["requested_region_removal"],
        percentage=config["requested_thresh_removal"],
    shell:
        "python3 scripts/modify_GTF.py {input.annotation_gtf} {params.percentage} {params.region} {output.modified_gtf} --keep_ids {input.kept_ids} "


rule remove_part_annotation_flair:
    input:
        annotation_gtf=config["reference_gtf"],
        kept_ids="results/discovery/kept_ids.csv",
    output:
        modified_gtf="results/discovery/flair/modified_annotation.gtf",
    params:
        region=config["requested_region_removal"],
        percentage=config["requested_thresh_removal"],
    shell:
        "python3 scripts/modify_GTF.py {input.annotation_gtf} {params.percentage} {params.region} {output.modified_gtf} --keep_ids {input.kept_ids} "


rule prepare_make_db_files_sirv:
    input:
        input_gtf="results/discovery/isoquant/modified_annotation_isoquant.gtf",
    output:
        output_db="results/discovery/isoquant/modified_annotation_isoquant.db",
    params:
        checklines=config["gffutils_checklines"],
    log:
        "logs/prepare/make_db_files/db.log",
    conda:
        "../envs/gffutils.yaml"
    script:
        "../scripts/create_db_files.py"


rule prepare_make_db_files_gencode:
    input:
        input_gtf="results/prepare/standardize_gtf_files/gencode.v45.primary_assembly.annotation.named.gtf",
    output:
        output_db="results/prepare/make_db_files/gencode.v45.primary_assembly.annotation.named.db",
    params:
        checklines=config["gffutils_checklines"],
    log:
        "logs/prepare/make_db_files/gencode.log",
    conda:
        "../envs/gffutils.yaml"
    script:
        "../scripts/create_db_files.py"


rule discovery_isoseq:
    input:
        unmapped_bams=config["path_unmapped"] + "{sample}.bam",
    output:
        output_bams="results/discovery/isoseq/{sample}/{sample}_output.bam",
        mapped_bams="results/discovery/isoseq/{sample}/{sample}_mapped_output.bam",
        collapsed_gff="results/discovery/isoseq/{sample}/{sample}_collapsed.gff",
        output_gtf="results/discovery/isoseq/{sample}/{sample}.isoseq_output.gtf",
    params:
        #pbmm2 only accepts .fasta/.fa/.fasta.gz/.fa.gz
        #genome_fasta=config["genome_fasta"].replace(".fna", ".fa.gz"),
        genome_fasta="/home/anna/seq_data/seq/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta",
    threads: 24
    benchmark:
        "results/benchmarks/{sample}.isoseq.benchmark.txt"
    log:
        stdout="logs/discovery/isoseq_{sample}.stdout",
        stderr="logs/discovery/isoseq_{sample}.stderr",
    conda:
        "../envs/pigeon.yaml"
    shell:
        """
        isoseq cluster2 {input.unmapped_bams} {output.output_bams} --singletons --num-threads {threads}
        pbmm2 align --preset ISOSEQ --sort {output.output_bams} {params.genome_fasta} {output.mapped_bams} --num-threads {threads}
        isoseq collapse --do-not-collapse-extra-5exons {output.mapped_bams} {output.collapsed_gff}
        gffread {output.collapsed_gff} -T -o {output.output_gtf}
        """


rule discovery_bambu:
    input:
        bam_file=config["datapath"] + "{sample}." + config["dataset"] + ".sorted.bam",
        #modified_gtf="results/discovery/bambu/modified_annotation.gtf", {input.modified_gtf}
    output:
        output_gtf="results/discovery/bambu/{sample}/{sample}.bambu_output.gtf",
    params:
        genome_fasta=config["genome_fasta"],
        output_dir="results/discovery/bambu/{sample}/",
        output_gtf="results/discovery/bambu/{sample}/{sample}."
        + config["dataset"]
        + ".bambu_output.gtf",
    benchmark:
        repeat("results/benchmarks/{sample}.bambu.benchmark.txt", 3)
    log:
        stdout="logs/discovery/bambu_{sample}.stdout",
        stderr="logs/discovery/bambu_{sample}.stderr",
    conda:
        "../envs/bambu.yaml"
    shell:
        """
        Rscript scripts/discovery_bambu.R {input.bam_file} NULL {params.genome_fasta} {params.output_dir} 1> {log.stdout} 2> {log.stderr}
        mv {params.output_gtf} {output.output_gtf}
        """


rule discovery_isoquant:
    input:
        bam_file=config["datapath"] + "{sample}." + config["dataset"] + ".sorted.bam",
        #modified_db="results/discovery/isoquant/modified_annotation_isoquant.db",--genedb {input.modified_db}
    output:
        final_output="results/discovery/isoquant/{sample}/{sample}.isoquant_output.gtf",
    params:
        region=config["requested_region_removal"],
        percentage=config["requested_thresh_removal"],
        genome_fasta=config["genome_fasta"],
        data_type="pacbio_ccs",
        output_file="results/discovery/isoquant/{sample}/{sample}/{sample}.transcript_models.gtf",
    threads: 7
    benchmark:
        "results/benchmarks/{sample}.isoquant.benchmark.txt"
    log:
        stdout="logs/discovery/isoquant_{sample}.stdout",
        stderr="logs/discovery/isoquant_{sample}.stderr",
    conda:
        "../envs/isoquant.yaml"
    # for do de novo discovery, leave out the --genedb option
    shell:
        """
        isoquant.py --reference {params.genome_fasta}  --bam {input.bam_file} --data_type {params.data_type} -o results/discovery/isoquant/{wildcards.sample} --prefix {wildcards.sample} --model_construction_strategy all --report_canonical all --report_novel_unspliced true --polya_requirement never --transcript_quantification unique_only --gene_quantification unique_only --threads {threads} 1> {log.stdout} 2> {log.stderr}
        cp {params.output_file} {output.final_output}
        """


rule discovery_scallop:
    input:
        bam_file=config["datapath_short"]
        + "{sample}."
        + config["dataset"]
        + ".sorted.bam",
    output:
        output_file="results/discovery/scallop/{sample}/{sample}.scallop_output.gtf",
    benchmark:
        repeat("results/benchmarks/{sample}.scallop.benchmark.txt", 3)
    log:
        stdout="logs/discovery/scallop_{sample}.stdout",
        stderr="logs/discovery/scallop_{sample}.stderr",
    conda:
        "../envs/scallop.yaml"
    shell:
        "scallop -i {input.bam_file} -o {output.output_file} --min_transcript_coverage 0 --min_single_exon_coverage 0 --min_transcript_length_base 0 --min_transcript_length_increase 0 --min_mapping_quality 0 --max_num_cigar 10000 --min_bundle_gap 0 --min_num_hits_in_bundle 0 --min_flank_length 0 --min_splice_bundary_hits 0 1> {log.stdout} 2> {log.stderr}"


rule discovery_StringTie2:
    input:
        bam_file=config["datapath_short"]
        + "{sample}."
        + config["dataset"]
        + ".sorted.bam",
        #modified_gtf="results/discovery/stringtie/modified_annotation.gtf", -G {input.modified_gtf}
    output:
        output_file="results/discovery/stringtie/{sample}/{sample}.stringtie_output.gtf",
    params:
        mid_output="results/discovery/stringtie/{sample}/{sample}.stringtie_mid_output.gtf",
    benchmark:
        repeat("results/benchmarks/{sample}.stringtie.benchmark.txt", 3)
    log:
        stdout="logs/discovery/stringtie2_{sample}.stdout",
        stderr="logs/discovery/stringtie2_{sample}.stderr",
    conda:
        "../envs/stringtie2.yaml"
    shell:
        """
        stringtie {input.bam_file} -o {params.mid_output} -m 30 -a 0 -j 0 -t -c 0.001 -s 0 -f 0  1> {log.stdout} 2> {log.stderr}
        python3 scripts/modify_stringtie_GTF.py {params.mid_output} {output.output_file}
        """


rule discovery_StringTie2_long:
    input:
        bam_file=config["datapath"] + "{sample}." + config["dataset"] + ".sorted.bam",
        #modified_gtf="results/discovery/stringtie/modified_annotation.gtf", -G {input.modified_gtf}
    output:
        output_file="results/discovery/stringtie_long/{sample}/{sample}.stringtie_long_output.gtf",
    params:
        mid_output="results/discovery/stringtie_long/{sample}/{sample}.stringtie_mid_output.gtf",
    benchmark:
        repeat("results/benchmarks/{sample}.stringtie_long.benchmark.txt", 3)
    log:
        stdout="logs/discovery/stringtie2_long_{sample}.stdout",
        stderr="logs/discovery/stringtie2_long_{sample}.stderr",
    conda:
        "../envs/stringtie2.yaml"
    shell:
        """
        stringtie {input.bam_file} -o {params.mid_output} -L -m 30 -a 0 -j 0 -t -c 0.001 -s 0 -f 0  1> {log.stdout} 2> {log.stderr}
        python3 scripts/modify_stringtie_GTF.py {params.mid_output} {output.output_file}
        """


rule discovery_StringTie2_mixed:
    input:
        short_bam=config["datapath_short"]
        + "{sample}."
        + config["dataset"]
        + ".sorted.bam",
        long_bam=config["datapath"] + "{sample}." + config["dataset"] + ".sorted.bam",
    output:
        output_file="results/discovery/stringtie_mixed/{sample}/{sample}.stringtie_mixed_output.gtf",
    benchmark:
        repeat("results/benchmarks/{sample}.stringtie2_mixed.benchmark.txt", 3)
    log:
        stdout="logs/discovery/stringtie2_mixed_{sample}.stdout",
        stderr="logs/discovery/stringtie2_mixed_{sample}.stderr",
    conda:
        "../envs/stringtie2.yaml"
    shell:
        "stringtie {input.short_bam} {input.long_bam} --mix -o {output.output_file} 1> {log.stdout} 2> {log.stderr}"


rule discovery_flair:
    input:
        transcriptome_gtf="results/discovery/flair/modified_annotation.gtf",
        bam_file=config["datapath"] + "{sample}." + config["dataset"] + ".sorted.bam",
        raw_reads=config["raw_data_path"]
        + "{sample}/"
        + config["dataset"]
        + ".fastq.gz",
    output:
        final_output="results/discovery/flair/{sample}/{sample}.flair_output.gtf",
    params:
        genome=config["genome_fasta"],
        bed12_output="results/discovery/flair/{sample}/{sample}."
        + config["dataset"]
        + ".sorted.bed",
        corrected_bed_output="results/discovery/flair/{sample}/{sample}_all_corrected.bed",
        mid_output="results/discovery/flair/{sample}/{sample}.isoforms.gtf",
        output_path="results/discovery/flair/{sample}/{sample}",
    threads: config["quantify_map_bam_threads"]
    benchmark:
        repeat("results/benchmarks/{sample}.flair.benchmark.txt", 3)
    log:
        stdout="logs/discovery/flair_{sample}.stdout",
        stderr="logs/discovery/flair_{sample}.stderr",
    conda:
        "../envs/flair.yaml"
    shell:
        """
        bamToBed -bed12 -i {input.bam_file} > {params.bed12_output}
        flair correct --query {params.bed12_output} --gtf {input.transcriptome_gtf} --genome {params.genome} --threads {threads} --output {params.output_path};
        flair collapse -g {params.genome} -q {params.corrected_bed_output} -r {input.raw_reads} --threads {threads} --support 0 --filter ginormous --gtf {input.transcriptome_gtf} --output {params.output_path}
        cp {params.mid_output} {output.final_output}
        """


rule discovery_FLAMES:
    input:
        transcriptome_gtf="results/discovery/bambu/modified_annotation.gtf",
        bam_file=config["datapath"] + "{sample}." + config["dataset"] + ".sorted.bam",
    output:
        final_output="results/discovery/flames/{sample}/{sample}.flames_output.gtf",
    params:
        genome=config["genome_fasta"],
        outdir="results/discovery/flames/{sample}",
        flames_config="config/config_FLAMES.json",
        mid_output="results/discovery/flames/{sample}/isoform_annotated.gff3",
        mid_output_2="results/discovery/flames/{sample}/isoform_annotated.gtf",
        transcriptome_gtf="results/discovery/flames/modified_annotation.gtf",
    threads: 8
    benchmark:
        repeat("results/benchmarks/{sample}.flames.benchmark.txt", 3)
    log:
        stdout="logs/discovery/flames_{sample}.stdout",
        stderr="logs/discovery/flames_{sample}.stderr",
    conda:
        "../envs/flames.yaml"
    shell:
        """
            #gffread {input.transcriptome_gtf} -T -o {params.transcriptome_gtf}
            #gtfsort --input {input.transcriptome_gtf} --output {params.transcriptome_gtf}
            Rscript scripts/discovery_FLAMES.R \
                {params.transcriptome_gtf} \
                {params.genome} \
                {input.bam_file} \
                {params.outdir} \
                {params.flames_config} \
            > {log.stdout} 2> {log.stderr}
            
            # Convert gff3 to gtf 
            gffread {params.mid_output} -T -o {params.mid_output_2}
            python3 scripts/modify_flames_GTF.py {params.mid_output_2} {output.final_output}
            """


rule discovery_mandalorion:
    input:
        fastq=config["raw_data_path"] + "{sample}.fastq.gz",
    output:
        output="results/discovery/mandalorion/{sample}/{sample}.mandalorion_output.gtf",
    params:
        annotation_gtf=config["reference_gtf"],
        genome_fasta=config["genome_fasta"],
        output_path="results/discovery/mandalorion/{sample}/",
        fofn_file="results/discovery/mandalorion/{sample}/{sample}.fofn",
        mid_output="results/discovery/mandalorion/{sample}/Isoforms.filtered.clean.gtf",
        modules="APDF",
    benchmark:
        repeat("results/benchmarks/{sample}.mandalorion.benchmark.txt", 3)
    log:
        stdout="logs/discovery/mandalorion_{sample}.stdout",
        stderr="logs/discovery/mandalorion_{sample}.stderr",
    conda:
        "../envs/mandalorion.yaml"
    shell:
        """
        echo {input.fastq} > {params.fofn_file}
        python3 scripts/Mandalorion/Mando.py -p {params.output_path} -g {params.annotation_gtf} -G {params.genome_fasta} -f {params.fofn_file} -M {params.modules} 1> {log.stdout} 2> {log.stderr}
        mv {params.mid_output} {output.output}
        """
