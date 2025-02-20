import pandas as pd
import argparse


def convert_file(input_file, output_file):

    df = pd.read_csv(input_file, sep="\t")

    output_df = pd.DataFrame(
        {"tname": df["Name"], "len": df["Length"], "num_reads": df["NumReads"]}
    )

    output_df.to_csv(output_file, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert a file format.")
    parser.add_argument("input", help="Path to the input file")
    parser.add_argument("output", help="Path to the output file")
    args = parser.parse_args()

    convert_file(args.input, args.output)
