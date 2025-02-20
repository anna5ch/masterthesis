import pandas as pd
import numpy as np
import argparse
import os
import re
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc, precision_recall_curve
from scipy.special import expit, logit


def read_kept_ids(kept_ids_file):
    kept_ids = set()
    with open(kept_ids_file, "r") as file:
        for line in file:
            line = line.strip()
            if not line.startswith("#") and line:  # Skip comment lines and empty lines
                kept_ids.add(line)
    return kept_ids


def get_label(gffcmp_file, kept_ids_file):
    kept_ids = read_kept_ids(kept_ids_file)
    label_df = pd.read_csv(gffcmp_file, sep="\t", header=None)
    label_df.columns = [
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
    # filter so only "transcript" rows are kept
    label_df = label_df[label_df["feature"] == "transcript"]
    # get relevant info
    label_df["transcript_id"] = label_df["attribute"].str.extract(
        r'transcript_id "([^"]+)"'
    )
    label_df["class_code"] = label_df["attribute"].str.extract(r'class_code "([^"]+)"')
    # if = or c then it its a 1
    label_df["label"] = label_df["class_code"].isin(["=", "c"]).astype(int)

    label_df["label"] = label_df.apply(
        lambda row: 1 if row["transcript_id"] in kept_ids else 0, axis=1
    )

    # Check if no labels were set to 1
    if label_df["label"].sum() == 0:
        # If no matches, set the first label to 1
        label_df.iloc[0, label_df.columns.get_loc("label")] = 1

    label_df = label_df.loc[:, ["transcript_id", "label"]].set_index("transcript_id")
    label_df.name = "label"

    return label_df


def create_feature_matrix(
    gtf_file, transcript_counts_file, gffcomp_file, nuc_file, kept_ids_file
):
    # Store transcript counts as df
    transcript_counts_df = pd.read_csv(transcript_counts_file, sep="\t")

    # Store gtf file as df
    gtf_df = pd.read_csv(gtf_file, sep="\t", comment="#", header=None)

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

    # Extract exons
    gtf_df_exons = gtf_df[gtf_df["feature"].str.contains("exon")].copy()
    gtf_df_exons.columns = [
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
    gtf_df_exons["transcript_id"] = gtf_df["attribute"].str.extract(
        r'transcript_id "([^"]+)"'
    )

    ########################
    #      GC Content      #
    ########################
    gc_content_df = pd.read_csv(nuc_file, sep="\t", comment="#", header=None)
    gc_content_df = gc_content_df[gc_content_df[2].str.contains("transcript")]
    gc_content_df["transcript_id"] = gc_content_df[8].str.extract(
        r'transcript_id "([^"]+)"'
    )
    gc_df = pd.DataFrame(
        {
            "transcript_id": gc_content_df["transcript_id"],
            "gc_content": gc_content_df[10],
        }
    )
    gc_df.set_index("transcript_id", inplace=True)

    feature_matrices = {}

    for day in transcript_counts_df.columns[2::]:
        ########################
        #       Log CPM        #
        ########################
        c = 0.001
        cpm = np.log1p(
            ((transcript_counts_df[day] + c) / (transcript_counts_df[day].sum() + c))
            * 1e6
        )

        ########################
        #   Gene proportion    #
        ########################
        day_counts = transcript_counts_df[["transcript_id", "gene_id", day]]

        gene_counts = day_counts.groupby("gene_id")[day].sum()
        gene_proportion = day_counts[day] / gene_counts[day_counts["gene_id"]].values
        if gene_proportion.isna().values.any():
            # only happens if day_counts is divided by 0
            gene_proportion = gene_proportion.fillna(0)

        gene_proportion_df = pd.DataFrame({"gene_proportion": gene_proportion})
        gene_proportion_df.index = day_counts["transcript_id"]
        gene_proportion_df.name = "gene_proportion"

        ########################
        #  Exonic length   #
        ########################
        transcript_lengths = gtf_df_exons.groupby("transcript_id")[
            ["end", "start"]
        ].apply(lambda group: (group["end"] - group["start"] + 1).sum())
        transcript_lengths.name = "exonic_length"

        # merge dfs with transcript_id as the index and features as columns
        feature_matrix = pd.DataFrame(
            {"transcript_id": transcript_counts_df["transcript_id"], "CPM": cpm}
        )
        feature_matrix.set_index("transcript_id", inplace=True)
        # feature_matrix = pd.merge(
        #    feature_matrix, transcript_lengths, left_index=True, right_index=True
        # )
        feature_matrix = pd.merge(
            feature_matrix, gene_proportion_df, left_index=True, right_index=True
        )
        # feature_matrix = pd.merge(
        #    feature_matrix, gc_df, left_index=True, right_index=True
        # )

        if gffcomp_file is not None:
            label_df = get_label(gffcomp_file, kept_ids_file)
            feature_matrix = pd.merge(
                feature_matrix, label_df, left_index=True, right_index=True
            )

            ########################
            # Machine Learning(LR) #
            ########################

            predicted_df = get_logistic_regression(feature_matrix)
            predicted_probabilities = predicted_df[
                ["transcript_id", "predicted_probability"]
            ]

            # Set transcript_id as index for predicted_probabilities
            predicted_probabilities.set_index("transcript_id", inplace=True)
            feature_matrix = pd.merge(
                feature_matrix,
                predicted_probabilities,
                left_index=True,
                right_index=True,
                how="inner",
            )

        feature_matrices[day] = feature_matrix
        return feature_matrices


def filter_by_value(feature_matrix, filter_option, threshold, output_file):
    # Filter rows where filtering option is above the threshold or NaN
    filtered_transcripts = feature_matrix[
        (feature_matrix[filter_option] > threshold)
        | feature_matrix[filter_option].isna()
    ]
    return filtered_transcripts


def get_logistic_regression(feature_matrix, random_state=42):
    feature_matrix_cleaned = feature_matrix.dropna()

    Y = feature_matrix_cleaned["label"]
    X = feature_matrix_cleaned.drop(columns=["label"])
    if len(Y.unique()) < 2:
        # If there are not at least 2 classes, log a warning and return an empty DataFrame
        print(
            "Warning: Logistic Regression requires at least 2 classes, but only one class is present."
        )
        return pd.DataFrame(), pd.DataFrame()

    model = LogisticRegression(random_state=random_state).fit(X, Y)

    feature_matrix_no_nan = feature_matrix.dropna(subset=X.columns)

    all_predictions = model.predict(feature_matrix_no_nan.drop(columns=["label"]))
    predicted_probabilities = model.predict_proba(
        feature_matrix_no_nan.drop(columns=["label"])
    )[:, 1]

    result_df = pd.DataFrame(
        {
            "transcript_id": feature_matrix_no_nan.index,
            "predicted_probability": predicted_probabilities,
            "predicted_label": all_predictions,
        }
    )

    return result_df


def filter_by_logistic_regression(
    feature_matrix, threshold, output_file, test_size=0.2, random_state=42
):
    feature_matrix_cleaned = feature_matrix.dropna()

    Y = feature_matrix_cleaned["label"]
    X = feature_matrix_cleaned.drop(columns=["label"])
    if len(Y.unique()) < 2:
        # If there are not at least 2 classes, log a warning and return an empty DataFrame
        print(
            "Warning: Logistic Regression requires at least 2 classes, but only one class is present."
        )
        return pd.DataFrame(), pd.DataFrame()

    X_train, X_test, y_train, y_test = train_test_split(
        X, Y, test_size=test_size, random_state=random_state
    )

    model = LogisticRegression(random_state=random_state)
    model.fit(X_train, y_train)

    feature_matrix_no_nan = feature_matrix.dropna(subset=X.columns)

    all_predictions = model.predict(feature_matrix_no_nan.drop(columns=["label"]))
    predicted_probabilities = model.predict_proba(
        feature_matrix_no_nan.drop(columns=["label"])
    )[:, 1]

    result_df = pd.DataFrame(
        {
            "transcript_id": feature_matrix_no_nan.index,
            "predicted_probability": predicted_probabilities,
            "predicted_label": all_predictions,
        }
    )

    filtered_transcripts = result_df[result_df["predicted_probability"] > threshold]

    # filtered_transcripts.to_csv(output_file, header=True, index=False, sep='\t', columns=['transcript_id'])
    return filtered_transcripts


def write_output(gtf_df, filtered_transcripts, output_file):
    merged_df = pd.merge(
        gtf_df,
        filtered_transcripts,
        how="inner",
        left_on="transcript_id",
        right_index=True,
    )
    merged_df.to_csv(output_file, header=True, index=False, sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filtering Transcripts Script")
    parser.add_argument(
        "--transcript_counts",
        type=str,
        help="Path to transcript counts file",
        required=True,
    )
    parser.add_argument("--gtf", type=str, help="Path to GTF file", required=True)
    parser.add_argument("--output", type=str, help="Output file name", required=True)
    parser.add_argument(
        "--ml_label", type=str, help="Path to gffcmp file", required=True
    )
    parser.add_argument(
        "--tracking", type=str, help="Path to gffcmp tracking file", required=True
    )
    parser.add_argument(
        "--nuc", type=str, help="Path to Bedtools NUC file", required=True
    )
    parser.add_argument(
        "--kept_ids_file", type=str, help="Path to kept_ids file", required=True
    )
    args = parser.parse_args()

    transcript_counts_file = args.transcript_counts
    gtf_file = args.gtf
    fm_output = args.output
    ML_label_file = args.ml_label
    tracking_file = args.tracking
    nuc_file = args.nuc
    kept_ids_file = args.kept_ids_file

    # create feature matrix
    print("Create feature matrix")
    feature_matrices = create_feature_matrix(
        gtf_file, transcript_counts_file, ML_label_file, nuc_file, kept_ids_file
    )
    # store feature matrix in new file

    feature_matrices["num_reads"].to_csv(fm_output, sep="\t", index=True)
