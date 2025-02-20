#!/usr/bin/env python3
import re
import argparse


def reformat_gtf(input_file, output_file):
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
                continue

            fields = line.strip().split("\t")
            attributes = fields[8]
            # Extract gene_id and transcript_id
            gene_id_match = re.search(r'gene_id "gene:([^"]+)"', attributes)
            transcript_id_match = re.search(
                r'transcript_id "transcript:([^"]+)"', attributes
            )

            if not gene_id_match or not transcript_id_match:
                continue
            gene_id = gene_id_match.group(1)
            transcript_id = transcript_id_match.group(1)

            # reorder attributes: gene_id first, then transcript_id
            new_attributes = f'gene_id "{gene_id}"; transcript_id "{transcript_id}"'

            fields[8] = new_attributes
            outfile.write("\t".join(fields) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reformat GTF file")
    parser.add_argument("input_file", help="Input GTF file")
    parser.add_argument("output_file", help="Output GTF file")
    args = parser.parse_args()

    reformat_gtf(args.input_file, args.output_file)
