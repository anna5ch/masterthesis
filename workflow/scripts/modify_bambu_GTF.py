import argparse


def fix_bambu_gtf_file(input_file, output_file):
    try:
        with open(input_file, "r") as infile, open(output_file, "w") as outfile:
            for line in infile:
                if line.startswith("#"):
                    # Skip header lines
                    outfile.write(line)
                    continue

                fields = line.strip().split("\t")
                if len(fields) < 9:
                    outfile.write(line)
                    continue

                attributes = fields[8]

                # Look gene_id with transcript_id in it
                if 'gene_id: "transcript_id' in attributes:
                    parts = attributes.split(";")
                    new_parts = []
                    for part in parts:
                        part = part.strip()
                        if part.startswith("gene_id"):
                            ids = part.split('"')
                            if len(ids) > 2 and " " in ids[1]:
                                # Extract second part of the gene_id which is the gene_id
                                fixed_gene_id = ids[1].split(" ")[1]
                                new_parts.append(f'gene_id "{fixed_gene_id}"')
                            else:
                                new_parts.append(part)
                        else:
                            new_parts.append(part)
                    fields[8] = "; ".join(new_parts).strip() + ";"
                outfile.write("\t".join(fields) + "\n")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description="Fix malformed gene_id entries in a GTF file."
        )
    parser.add_argument("input")
    parser.add_argument("output")
    args = parser.parse_args()

    fix_bambu_gtf_file(args.input, args.output)
