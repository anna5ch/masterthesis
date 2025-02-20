import pandas as pd
from gtfparse import read_gtf
import argparse


def parse_gtf(gtf_file):
    # read gtf file into data frame
    gtf_data = read_gtf(gtf_file)
    gtf_df = gtf_data[["transcript_id", "gene_id"]].to_pandas().drop_duplicates()
    return gtf_df


def add_gene_id_to_quant(quant_file, gtf_df, output_file):
    # include gene_id to transcript id of the quant file
    quant_df = pd.read_csv(quant_file, sep="\t")

    if "tname" in quant_df.columns:
        # merge quant_df (transcript_ids and num_reads) with gtf_df(gene_id and transcript_id) on transcript id
        merged_df = quant_df.merge(
            gtf_df, left_on="tname", right_on="transcript_id", how="left"
        )

        # drop redundant cols
        merged_df = merged_df.drop(columns=["tname", "len"])
        merged_df["num_reads"] = merged_df["num_reads"].astype(int)

        # gene_id as the first column, then transcript_id and num_reads
        cols = ["gene_id", "transcript_id", "num_reads"] + [
            col
            for col in merged_df.columns
            if col not in ["gene_id", "transcript_id", "num_reads"]
        ]
        merged_df = merged_df[cols]
    else:
        # merge quant_df (transcript_ids and num_reads) with gtf_df(gene_id and transcript_id) on transcript id
        merged_df = quant_df.merge(
            gtf_df, left_on="Name", right_on="transcript_id", how="left"
        )

        # drop redundant cols
        merged_df = merged_df.drop(columns=["Name", "Length", "EffectiveLength", "TPM"])
        merged_df["num_reads"] = merged_df["NumReads"].astype(int)

        # geneq_id as the first column, then transcript_id and num_reads
        cols = ["gene_id", "transcript_id", "num_reads"] + [
            col
            for col in merged_df.columns
            if col not in ["gene_id", "transcript_id", "num_reads"]
        ]
        merged_df = merged_df[cols]

    # write output
    merged_df.sort_values(["gene_id", "transcript_id"]).to_csv(
        output_file, sep="\t", index=False
    )


if __name__ == "__main__":
    # read arguments
    parser = argparse.ArgumentParser(
        description="Add gene_id to .quant file from oarfish from GTF file."
    )
    parser.add_argument(
        "quant_file", type=str, help="Path to the input .quant file of oarfish"
    )
    parser.add_argument(
        "gtf_file", type=str, help="Path to the standardized/collapsed input GTF file"
    )
    parser.add_argument("output_file", type=str, help="Path to the output .quant file")
    args = parser.parse_args()

    # read gtf file to df
    gtf_df = parse_gtf(args.gtf_file)
    add_gene_id_to_quant(args.quant_file, gtf_df, args.output_file)
