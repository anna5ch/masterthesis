rule prepare_oarfish_filtering_files:
    input:
        input_gtf="results/collapse/{tool}/{sample}/{sample}.{tool}_output_collapsed.gtf",
        input_quant="results/quantify/{tool}/oarfish_{tool}_{sample}/oarfish_{tool}_{sample}.quant",
        input_tc_file="results/quantify/{tool}/oarfish_{tool}_{sample}/oarfish_{sample}_transcript_counts.tsv",
        transcriptome_fasta="results/quantify/{tool}/{sample}_transcriptome.fasta",
    output:
        nuc_file="results/filtering/{tool}/{sample}/oarfish_{sample}.{tool}_nuc.gtf",
        gffcmp_out="results/filtering/{tool}/{sample}/oarfish_gffcmp_{sample}.{tool}",
        out="results/quantify/{tool}/oarfish_{tool}_{sample}/oarfish_{sample}_transcript_counts_formatted.tsv",
        tracking_out="results/filtering/{tool}/{sample}/oarfish_gffcmp_{sample}.{tool}.tracking",
    params:
        reference=config["reference_gtf"],
        genome_fasta=config["genome_fasta"],
    log:
        stdout="logs/filter/oarfish_prepare/{tool}_{sample}.stdout",
        stderr="logs/filter/oarfish_prepare/{tool}_{sample}.stderr",
    conda:
        "../envs/base.yaml"
    shell:
        """
        bedtools nuc -fi {params.genome_fasta} -bed {input.input_gtf} > {output.nuc_file}
        gffcompare -r {params.reference} {input.input_gtf} -o {output.gffcmp_out}
        python3 scripts/expand_quant.py {input.input_quant} {input.input_gtf} {output.out} 1> {log.stdout} 2> {log.stderr}
        """


rule run_oarfish_feature_matrix:
    input:
        input_gtf="results/collapse/{tool}/{sample}/{sample}.{tool}_output_collapsed.gtf",
        quant_file="results/quantify/{tool}/oarfish_{tool}_{sample}/oarfish_{sample}_transcript_counts_formatted.tsv",
        nuc_file="results/filtering/{tool}/{sample}/oarfish_{sample}.{tool}_nuc.gtf",
        gffcmp_out="results/filtering/{tool}/{sample}/oarfish_gffcmp_{sample}.{tool}",
        kept_ids_file="results/discovery/kept_ids.csv",
    output:
        #filtered_gtf="results/filtering/{tool}/{sample}/oarfish_{sample}.{tool}_output_filtered_{requested_filters}_{requested_threshs}.gtf",
        feature_matrix_file="results/filtering/{tool}/{sample}/oarfish_{sample}.{tool}_output_feature_matrix.tsv",
    log:
        stdout="logs/filter/feature_matrix/{tool}_{sample}.stdout",
        stderr="logs/filter/feature_matrix/{tool}_{sample}.stderr",
    conda:
        "../envs/base.yaml"
    shell:
        """
        python3 scripts/create_feature_matrix.py --transcript_counts {input.quant_file} --gtf {input.input_gtf} --output {output.feature_matrix_file} --ml_label {input.gffcmp_out}.annotated.gtf --tracking {input.gffcmp_out}.tracking --nuc {input.nuc_file} --kept_ids_file {input.kept_ids_file} 1> {log.stdout} 2> {log.stderr}
        """


rule run_oarfish_filtering:
    input:
        feature_matrix_file="results/filtering/{tool}/{sample}/oarfish_{sample}.{tool}_output_feature_matrix.tsv",
        input_gtf="results/collapse/{tool}/{sample}/{sample}.{tool}_output_collapsed.gtf",
        quant_file="results/quantify/{tool}/oarfish_{tool}_{sample}/oarfish_{sample}_transcript_counts_formatted.tsv",
        nuc_file="results/filtering/{tool}/{sample}/oarfish_{sample}.{tool}_nuc.gtf",
        gffcmp_out="results/filtering/{tool}/{sample}/oarfish_gffcmp_{sample}.{tool}",
    output:
        filtered_gtf="results/filtering/{tool}/{sample}/oarfish_{sample}.{tool}_output_filtered_{requested_filters}_{requested_threshs}.gtf",
    log:
        stdout="logs/filter/oarfish/{tool}_{sample}_{requested_filters}_{requested_threshs}.stdout",
        stderr="logs/filter/oarfish/{tool}_{sample}_{requested_filters}_{requested_threshs}.stderr",
    conda:
        "../envs/base.yaml"
    shell:
        """
        extra_flag=""
        # Liste von Tools, für die die Flag gilt
        if [[ "{wildcards.tool}" == "bambu" || "{wildcards.tool}" == "isoquant" || "{wildcards.tool}" == "flames" ]]; then
            extra_flag="--kept_ids results/discovery/kept_ids.csv"
        fi
        python3 scripts/filter_and_create_plots.py --gtf {input.input_gtf} --feature_matrix_file {input.feature_matrix_file} --output {output.filtered_gtf} --ml_label {input.gffcmp_out}.annotated.gtf --tracking {input.gffcmp_out}.tracking $extra_flag 1> {log.stdout} 2> {log.stderr}
        """


rule prepare_salmon_filtering_files:
    input:
        input_gtf="results/collapse/{tool}/{sample}/{sample}.{tool}_output_collapsed.gtf",
        input_quant="results/quantify/{tool}/salmon_{tool}_{sample}/salmon_{tool}_{sample}.quant",
        input_tc_file="results/quantify/{tool}/salmon_{tool}_{sample}/salmon_{sample}_transcript_counts.tsv",
        transcriptome_fasta="results/quantify/{tool}/{sample}_transcriptome.fasta",
    output:
        nuc_file="results/filtering/{tool}/{sample}/salmon_{sample}.{tool}_nuc.gtf",
        gffcmp_out="results/filtering/{tool}/{sample}/salmon_gffcmp_{sample}.{tool}",
        out="results/quantify/{tool}/salmon_{tool}_{sample}/salmon_{sample}_transcript_counts_formatted.tsv",
        tracking_out="results/filtering/{tool}/{sample}/salmon_gffcmp_{sample}.{tool}.tracking",
    params:
        reference=config["reference_gtf"],
        fasta_file=config["genome_fasta"],
    log:
        stdout="logs/filter/salmon_prepare/{tool}_{sample}.stdout",
        stderr="logs/filter/salmon_prepare/{tool}_{sample}.stderr",
    conda:
        "../envs/base.yaml"
    shell:
        """
        bedtools nuc -fi {params.fasta_file} -bed {input.input_gtf} > {output.nuc_file}
        gffcompare -r {params.reference} {input.input_gtf} -o {output.gffcmp_out}
        python3 scripts/expand_quant.py {input.input_quant} {input.input_gtf} {output.out} 1> {log.stdout} 2> {log.stderr}
        """


rule run_salmon_feature_matrix:
    input:
        input_gtf="results/collapse/{tool}/{sample}/{sample}.{tool}_output_collapsed.gtf",
        quant_file="results/quantify/{tool}/salmon_{tool}_{sample}/salmon_{sample}_transcript_counts_formatted.tsv",
        nuc_file="results/filtering/{tool}/{sample}/salmon_{sample}.{tool}_nuc.gtf",
        gffcmp_out="results/filtering/{tool}/{sample}/salmon_gffcmp_{sample}.{tool}",
        kept_ids_file="results/discovery/kept_ids.csv",
    output:
        #filtered_gtf="results/filtering/{tool}/{sample}/salmon_{sample}.{tool}_output_filtered_{requested_filters}_{requested_threshs}.gtf",
        feature_matrix_file="results/filtering/{tool}/{sample}/salmon_{sample}.{tool}_output_feature_matrix.tsv",
    log:
        stdout="logs/filter/feature_matrix/{tool}_{sample}.stdout",
        stderr="logs/filter/feature_matrix/{tool}_{sample}.stderr",
    conda:
        "../envs/base.yaml"
    shell:
        """
        python3 scripts/create_feature_matrix.py --transcript_counts {input.quant_file} --gtf {input.input_gtf} --output {output.feature_matrix_file} --ml_label {input.gffcmp_out}.annotated.gtf --tracking {input.gffcmp_out}.tracking --nuc {input.nuc_file} --kept_ids_file {input.kept_ids_file} 1> {log.stdout} 2> {log.stderr}
        """


rule run_salmon_filtering:
    input:
        feature_matrix_file="results/filtering/{tool}/{sample}/salmon_{sample}.{tool}_output_feature_matrix.tsv",
        input_gtf="results/collapse/{tool}/{sample}/{sample}.{tool}_output_collapsed.gtf",
        quant_file="results/quantify/{tool}/salmon_{tool}_{sample}/salmon_{sample}_transcript_counts_formatted.tsv",
        nuc_file="results/filtering/{tool}/{sample}/salmon_{sample}.{tool}_nuc.gtf",
        gffcmp_out="results/filtering/{tool}/{sample}/salmon_gffcmp_{sample}.{tool}",
    output:
        filtered_gtf="results/filtering/{tool}/{sample}/salmon_{sample}.{tool}_output_filtered_{requested_filters}_{requested_threshs}.gtf",
    log:
        stdout="logs/filter/salmon/{tool}_{sample}_{requested_filters}_{requested_threshs}.stdout",
        stderr="logs/filter/salmon/{tool}_{sample}_{requested_filters}_{requested_threshs}.stderr",
    conda:
        "../envs/base.yaml"
    shell:
        """
        extra_flag=""
        # Liste von Tools, für die die Flag gilt
        if [[ "{wildcards.tool}" == "bambu" || "{wildcards.tool}" == "isoquant" || "{wildcards.tool}" == "flames" ]]; then
            extra_flag="--kept_ids results/discovery/kept_ids.csv"
        fi
        python3 scripts/filter_and_create_plots.py --gtf {input.input_gtf} --feature_matrix_file {input.feature_matrix_file} --output {output.filtered_gtf} --ml_label {input.gffcmp_out}.annotated.gtf --tracking {input.gffcmp_out}.tracking $extra_flag 1> {log.stdout} 2> {log.stderr}
        """
