import argparse
from typing import List, Optional
import numpy as np
import pandas as pd
from typeguard import typechecked


@typechecked
def aggregate_oarfish(
    sample_paths: List[str],
    sample_names: List[str],
    gtf_annotation_path: str,
    output_path: str,
    separate_sirvs: Optional[bool] = False,
    sirv_output_path: Optional[str] = None,
) -> None:
    df_dict = {}
    for sample_path, sample_name in zip(sample_paths, sample_names):
        df = pd.read_csv(sample_path, sep="\t", comment="#")
        if sample_name == sample_names[0]:
            df_dict["transcript_id"] = df["tname"].values.tolist()
        df_dict[sample_name] = df["num_reads"].values.tolist()
    df = pd.DataFrame(df_dict)
    gtf_df = pd.read_csv(gtf_annotation_path, sep="\t", header=None, comment="#")
    gtf_df.columns = [
        "seqname",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attribute",
    ]
    gtf_df["transcript_id"] = gtf_df["attribute"].str.extract(
        r'transcript_id "([^"]+)"'
    )
    gtf_df["gene_id"] = gtf_df["attribute"].str.extract(r'gene_id "([^"]+)"')
    gene_transcript_df = gtf_df[["gene_id", "transcript_id"]].copy(deep=True)
    del gtf_df
    df = df.merge(
        gene_transcript_df,
        how="left",
        left_on="transcript_id",
        right_on="transcript_id",
    ).drop_duplicates()

    if separate_sirvs:
        gencode_mask = df.transcript_id.str.contains("ENST")
        df_sirvs = df.loc[np.logical_not(gencode_mask)][
            ["gene_id", "transcript_id"] + sample_names
        ]
        df_gencode = df.loc[gencode_mask][["gene_id", "transcript_id"] + sample_names]
        df_sirvs.to_csv(
            path_or_buf=sirv_output_path, sep="\t", header=True, index=False
        )
    else:
        df_gencode = df[["gene_id", "transcript_id"] + sample_names]

    df_gencode.sort_values(["gene_id", "transcript_id"]).to_csv(
        path_or_buf=output_path, sep="\t", header=True, index=False
    )
    return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--sample_paths",
        nargs="*",
        type=str,
    )
    parser.add_argument(
        "--sample_names",
        nargs="*",
        type=str,
    )
    parser.add_argument("--gtf_annotation_path", type=str)
    parser.add_argument("--output_path", type=str)
    parser.add_argument("--separate_sirvs", type=bool, default=False)
    parser.add_argument("--sirv_output_path", type=str, default="")

    args = parser.parse_args()
    aggregate_oarfish(
        sample_paths=args.sample_paths,
        sample_names=args.sample_names,
        gtf_annotation_path=args.gtf_annotation_path,
        output_path=args.output_path,
        separate_sirvs=args.separate_sirvs,
        sirv_output_path=args.sirv_output_path,
    )
