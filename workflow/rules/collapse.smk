rule standardize:
    input:
        discovery_gtf="results/discovery/{tool}/{sample}/{sample}.{tool}_output.gtf",
    output:
        standardized_discovery_gtf="results/collapse/{tool}/{sample}/{sample}.{tool}_output_standardized.gtf",
    conda:
        "../envs/gffread.yaml"
    shell:
        "gffread -E {input.discovery_gtf} -T -o- | more > {output.standardized_discovery_gtf}"


rule collapse:
    input:
        discovery_gtf="results/collapse/{tool}/{sample}/{sample}.{tool}_output_standardized.gtf",
    output:
        collapsed_gtf="results/collapse/{tool}/{sample}/{sample}.{tool}_output_collapsed.gtf",
    params:
        prefix="results/collapse/{tool}/{sample}/{sample}.{tool}",
        mid_output="results/collapse/{tool}/{sample}/{sample}.{tool}.combined.gtf",
    conda:
        "../envs/gffread.yaml"
    shell:
        """
        gffcompare -X {input.discovery_gtf} -o {params.prefix}
        python3 scripts/replace_transcript_id_with_oid.py {params.mid_output} {output.collapsed_gtf}
        """
