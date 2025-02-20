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


def get_label(gffcmp_file):
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
    label_df = label_df.loc[:, ["transcript_id", "label"]].set_index("transcript_id")
    label_df.name = "label"

    return label_df


def parse_feature_matrix_file(feature_matrix_file):
    df = pd.read_csv(feature_matrix_file, sep="\t")
    return df


def parse_gtf_df(gtf_file):
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
    gtf_df["transcript_id"] = gtf_df["attribute"].str.extract(
        r'transcript_id "([^"]+)"'
    )
    return gtf_df


def filter_by_value(feature_matrix, filter_option, threshold):
    # Filter rows where filtering option is above the threshold or NaN
    filtered_transcripts = feature_matrix[
        (feature_matrix[filter_option] > threshold)
        | feature_matrix[filter_option].isna()
    ]
    return filtered_transcripts


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
        right_on="transcript_id",
    )
    merged_df.to_csv(output_file, header=True, index=False, sep="\t")


def parse_tracking_file(tracking_file):
    data = []
    with open(tracking_file) as file:
        for line in file:
            parts = line.strip().split("\t")
            transcript_id = parts[4].strip().split("|")[1]
            match_status = parts[3].strip()
            y_true = 1 if match_status == "=" or match_status == "c" else 0
            data.append((transcript_id, y_true))
    return pd.DataFrame(data, columns=["transcript_id", "y_true"])


def generate_prec_rec_plot(tp_fp_array, value_array, output_file, samplename):
    y_true = np.array(tp_fp_array)
    y_scores = expit(np.array(value_array))

    precision, recall, thresholds = precision_recall_curve(
        y_true, y_scores, pos_label=1
    )
    pr_auc = auc(recall, precision)

    # Plot precision-recall curve
    plt.figure()
    plt.plot(
        recall,
        precision,
        color="darkgreen",
        lw=2,
        label="PR curve (area = %0.2f)" % pr_auc,
    )
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title("Precision-Recall curve for ")
    plt.legend(loc="lower left")
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.savefig(output_file)
    plt.close()


def generate_roc_curve(tp_fp_array, value_array, output_file):
    y_true = np.array(tp_fp_array)
    y_scores = expit(np.array(value_array))

    if np.any(np.isnan(y_true)) or np.any(np.isnan(y_scores)):
        print("One of the arrays contains NaN values.")
        exit()

    # Compute ROC curve and AUC area
    fpr, tpr, thresholds = roc_curve(y_true, y_scores, pos_label=1)
    roc_auc = auc(fpr, tpr)

    # Plot ROC curve
    plt.figure()
    plt.plot(
        fpr, tpr, color="darkorange", lw=2, label="ROC curve (area = %0.2f)" % roc_auc
    )
    plt.plot([0, 1], [0, 1], color="navy", lw=2, linestyle="--")
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("ROC curve")
    plt.legend(loc="lower right")
    plt.savefig(output_file)
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filtering Transcripts Script")
    parser.add_argument("--gtf", type=str, help="Path to GTF file", required=True)
    parser.add_argument(
        "--feature_matrix_file",
        type=str,
        help="Path to feature matrix file",
        required=True,
    )
    parser.add_argument(
        "--threshold", type=float, help="Threshold value", required=False, default=0.5
    )
    parser.add_argument("--output", type=str, help="Output file name", required=True)
    parser.add_argument(
        "--ml_label", type=str, help="Path to gffcmp file", required=True
    )
    parser.add_argument(
        "--tracking", type=str, help="Path to gffcmp tracking file", required=True
    )
    parser.add_argument(
        "--filter_method",
        type=str,
        help="Filtering method, default: gene_proportion",
        required=False,
        default="gene_proportion",
    )
    parser.add_argument(
        "--kept_ids",
        type=str,
        help="Filtered transcripts in kept_ids.csv",
        required=False,
        default=None,
    )
    args = parser.parse_args()
    gtf_file = args.gtf
    threshold = args.threshold
    output = args.output
    ML_label_file = args.ml_label
    tracking_file = args.tracking
    feature_matrix_file = args.feature_matrix_file
    plot_out = args.output.replace("filtered_", "").replace(".gtf", "_plot.png")
    plot_out2 = args.output.replace("filtered_", "").replace(
        ".gtf", "_plot_prec_recall.png"
    )
    filter_method = args.filter_method
    kept_ids = args.kept_ids

    feature_matrix = parse_feature_matrix_file(feature_matrix_file)
    gtf_df = parse_gtf_df(gtf_file)

    output_file = output
    if filter_method in ["gene_proportion", "CPM", "exonic_length", "gc_content"]:
        filtered_transcripts = filter_by_value(feature_matrix, filter_method, threshold)
        write_output(gtf_df, filtered_transcripts, output_file)
        tracking_df = parse_tracking_file(tracking_file)
        merged_df = pd.merge(
            tracking_df,
            feature_matrix,
            on="transcript_id",
            how="inner",
        )
        if kept_ids is not None:
            # read file with kept_ids, to only have new discovered isoforms for evaluation
            ignore_transcripts_df = pd.read_csv(
                kept_ids, header=None, names=["transcript_id"]
            )
            # remove kept_ids:
            filtered_df = merged_df[
                ~merged_df["transcript_id"].isin(ignore_transcripts_df["transcript_id"])
            ]
            generate_roc_curve(
                filtered_df["y_true"], filtered_df[filter_method], plot_out
            )
            generate_prec_rec_plot(
                filtered_df["y_true"],
                filtered_df[filter_method],
                plot_out2,
                filter_method,
            )
        else:
            generate_roc_curve(merged_df["y_true"], merged_df[filter_method], plot_out)
            generate_prec_rec_plot(
                merged_df["y_true"], merged_df[filter_method], plot_out2, filter_method
            )

    elif filter_method == "ml":
        filtered_transcripts = filter_by_logistic_regression(
            feature_matrix, threshold, output_file
        )
        write_output(gtf_df, filtered_transcripts, output_file)
        tracking_df = parse_tracking_file(tracking_file)
        merged_df = pd.merge(
            tracking_df,
            feature_matrix,
            on="transcript_id",
            how="left",
        )
        generate_roc_curve(merged_df["y_true"], merged_df["label"], plot_out)
        generate_prec_rec_plot(
            merged_df["y_true"], merged_df[filter_method], plot_out2, filter_method
        )
    else:
        raise ValueError(f"Invalid filter method: {filter_method}")
