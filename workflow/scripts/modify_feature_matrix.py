import sys
import pandas as pd


def load_gtf_transcript_ids(gtf_file):
    transcript_ids = set()
    with open(gtf_file, "r") as gtf:
        for line in gtf:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if parts[2] == "transcript":
                attributes = parts[8]
                transcript_id = attributes.split('transcript_id "')[1].split('"')[0]
                transcript_ids.add(transcript_id)
    return transcript_ids


def filter_feature_matrix(feature_matrix_file, gtf_transcript_ids, output_file):
    # Load the feature matrix
    df = pd.read_csv(feature_matrix_file, sep="\t")

    # only keep the "novel ones" so the ones which where removed in the beginning
    filtered_df = df[df["transcript_id"].isin(gtf_transcript_ids)]

    filtered_df.to_csv(output_file, sep="\t", index=False)


if __name__ == "__main__":
    gtf_file = sys.argv[1]
    feature_matrix_file = sys.argv[2]
    output_file = sys.argv[3]

    # Load transcript_ids from the GTF file
    gtf_transcript_ids = load_gtf_transcript_ids(gtf_file)

    # Filter the feature matrix and save the output
    filter_feature_matrix(feature_matrix_file, gtf_transcript_ids, output_file)
