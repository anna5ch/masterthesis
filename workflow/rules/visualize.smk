rule modify_GTF_new:
    input:
        input_gtf="results/collapse/{tool}/{sample}/{sample}.{tool}_output_collapsed.gtf",
        kept_ids="results/discovery/kept_ids.csv",
    output:
        "results/collapse/{tool}/{sample}/{sample}.{tool}_output_collapsed_only_new.gtf",
    params:
        chromosomes=config["requested_region_removal"],
    shell:
        """
            python3 scripts/filter_gtf_by_chromosome.py {input.input_gtf} {output} --kept_ids {input.kept_ids}
        """


rule gffcmp_only_new:
    input:
        "results/collapse/{tool}/{sample}/{sample}.{tool}_output_collapsed_only_new.gtf",
    output:
        gffcmp="results/filtering/{tool}/{sample}/only_new_{sample}.{tool}",
        tracking="results/filtering/{tool}/{sample}/only_new_{sample}.{tool}.tracking",
    params:
        reference=config["reference_gtf"],
    conda:
        "../envs/base.yaml"
    shell:
        "gffcompare -r {params.reference} {input} -o {output.gffcmp}"


rule create_boxplots:
    input:
        expand(
            "results/filtering/{tool}/{sample}/only_new_{sample}.{tool}",
            tool=config["requested_tools"],
            sample=config["samples"],
        ),
    output:
        "results/visualize/boxplot.png",
    script:
        "../scripts/boxplots.R"


rule modify_feature_matrix_salmon:
    input:
        gtf_file="results/collapse/{tool}/{sample}/{sample}.{tool}_output_collapsed_only_new.gtf",
        feature_matrix="results/filtering/{tool}/{sample}/salmon_{sample}.{tool}_output_feature_matrix.tsv",
    output:
        "results/filtering/{tool}/{sample}/salmon_{sample}.{tool}_output_feature_matrix_only_new.tsv",
    shell:
        "python3 scripts/modify_feature_matrix.py {input.gtf_file} {input.feature_matrix} {output}"


rule modify_feature_matrix_oarfish:
    input:
        gtf_file="results/collapse/{tool}/{sample}/{sample}.{tool}_output_collapsed_only_new.gtf",
        feature_matrix="results/filtering/{tool}/{sample}/oarfish_{sample}.{tool}_output_feature_matrix.tsv",
    output:
        "results/filtering/{tool}/{sample}/oarfish_{sample}.{tool}_output_feature_matrix_only_new.tsv",
    shell:
        "python3 scripts/modify_feature_matrix.py {input.gtf_file} {input.feature_matrix} {output}"


rule plot_tpr_fpr:
    input:
        tracking_files=expand(
            "results/filtering/{tool}/{sample}/oarfish_gffcmp_{sample}.{tool}.tracking",
            sample=config["samples"],
            tool=config["requested_tools_long"],
        ),
        feature_matrix_files=expand(
            "results/filtering/{tool}/{sample}/oarfish_{sample}.{tool}_output_feature_matrix.tsv",
            sample=config["samples"],
            tool=config["requested_tools_long"],
        ),
        tracking_files2=expand(
            "results/filtering/{tool}/{sample}/salmon_gffcmp_{sample}.{tool}.tracking",
            sample=config["samples"],
            tool=config["requested_tools_short"],
        ),
        feature_matrix_files2=expand(
            "results/filtering/{tool}/{sample}/salmon_{sample}.{tool}_output_feature_matrix.tsv",
            sample=config["samples"],
            tool=config["requested_tools_short"],
        ),
    output:
        "results/visualize/tpr_fpr_plot_{sample}.png",
    params:
        num_removed=69,
    script:
        "../scripts/create_TPR_FPR_plots.R"


rule plot_tpr_fpr_only_new:
    input:
        tracking_files=expand(
            "results/filtering/{tool}/{sample}/only_new_{sample}.{tool}.tracking",
            sample=config["samples"],
            tool=config["requested_tools_long"],
        ),
        feature_matrix_files=expand(
            "results/filtering/{tool}/{sample}/oarfish_{sample}.{tool}_output_feature_matrix_only_new.tsv",
            sample=config["samples"],
            tool=config["requested_tools_long"],
        ),
        tracking_files2=expand(
            "results/filtering/{tool}/{sample}/only_new_{sample}.{tool}.tracking",
            sample=config["samples"],
            tool=config["requested_tools_short"],
        ),
        feature_matrix_files2=expand(
            "results/filtering/{tool}/{sample}/salmon_{sample}.{tool}_output_feature_matrix_only_new.tsv",
            sample=config["samples"],
            tool=config["requested_tools_short"],
        ),
    output:
        "results/visualize/only_new/tpr_fpr_plot_{sample}_only_new.png",
    params:
        num_removed=69,
    script:
        "../scripts/create_TPR_FPR_plots_only_new.R"


rule plot_benchmark:
    input:
        expand(
            "results/benchmarks/{sample}.{tool}.benchmark.txt",
            sample=config["samples"],
            tool=config["requested_tools"],
        ),
    output:
        "results/visualize/benchmark/{sample}_benchmark_plots.png",
    params:
        "results/visualize/benchmark/",
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/create_benchmark_plots.py"


rule plot_bam:
    input:
        expand(
            config["datapath"] + "{sample}/{sample}.sirvs.sorted.bam",
            sample=config["samples"],
        ),
    output:
        "results/visualize/read_depth_and_reads.png",
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/create_bam_plots.py"
