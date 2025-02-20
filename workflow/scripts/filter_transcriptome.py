import pandas as pd
import re


def filter_transcriptome(expressed_transcripts, gtf_file, output_file):
    transcript_ids = expressed_transcripts["transcript_id"].tolist()
    filtered_gtf = []
    genes = []
    gene_ids_with_transcripts = set()
    with open(gtf_file, "r") as gtf:
        for line in gtf:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if parts[2] == "transcript" or parts[2] == "exon":
                attributes = parts[8]
                transcript_id = None
                for attribute in attributes.split(";"):
                    if "transcript_id" in attribute:
                        transcript_id = attribute.split('"')[1]
                    if "gene_id" in attribute:
                        gene_id = attribute.split('"')[1]
                if transcript_id in transcript_ids:
                    filtered_gtf.append(line)
                    gene_ids_with_transcripts.add(gene_id)
            elif parts[2] == "gene":
                genes.append(line)

    for gene in genes:
        gene_id = re.search(r'gene_id "([^"]+)"', gene).group(1)
        if gene_id in gene_ids_with_transcripts:
            filtered_gtf.append(gene)

    with open(output_file, "w") as output:
        for line in filtered_gtf:
            output.write(line)


expression_data = pd.read_csv(
    "/home/anna/seq_data/seq/expressed_genes_wtc11.tsv", sep="\t"
)

# Filter transcripts with expression in day0/day5 for Illumina and kinnex
expressed_transcripts = expression_data[
    (expression_data["illumina_day0"] == True)
    & (expression_data["kinnex_day0"] == True)
    & (expression_data["illumina_day5"] == True)
    & (expression_data["kinnex_day5"] == True)
]


gtf_file = (
    "/home/dwissel/data/seq/ref/clean/gencode.v45.primary_assembly.annotation.named.gtf"
)
output_file = (
    "/home/anna/seq_data/seq/gencode.v45.primary_assembly.annotation.named.filtered.gtf"
)

filter_transcriptome(expressed_transcripts, gtf_file, output_file)
