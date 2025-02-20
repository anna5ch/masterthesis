import argparse


def filter_and_write_gtf(file_path, output):
    gene_ids_with_transcripts = set()
    with open(file_path, "r") as file:
        for line in file:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if parts[2] == "transcript":
                info = parts[8]
                for attr in info.split(";"):
                    attr = attr.strip()
                    if "gene_id" in attr:
                        gene_id = attr.split('"')[1]
                        gene_ids_with_transcripts.add(gene_id)
    with open(file_path, "r") as file, open(output, "w") as outfile:
        for line in file:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            info = parts[8]
            for attr in info.split(";"):
                attr = attr.strip()
                if "gene_id" in attr:
                    gene_id = attr.split('"')[1]
                    if gene_id in gene_ids_with_transcripts:
                        outfile.write(line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Remove part of transcripts from a GTF file."
    )
    parser.add_argument("gtf", type=str, help="The GTF file in question")
    parser.add_argument("output_file", type=str, help="Output file path")
    args = parser.parse_args()
    filter_and_write_gtf(args.gtf, args.output_file)
