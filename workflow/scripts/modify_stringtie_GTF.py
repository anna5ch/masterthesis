import re
import argparse


def modify_gtf(input_file, output_file):
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9:
                # if new transcript, then no reference_id
                outfile.write(line)
            else:
                attributes = fields[8]
                gene_id = re.search(r'gene_id "([^"]+)"', attributes)
                transcript_id = re.search(r'transcript_id "([^"]+)"', attributes)
                reference_id = re.search(r'reference_id "([^"]+)"', attributes)
                ref_gene_id = re.search(r'ref_gene_id "([^"]+)"', attributes)

                if not reference_id and not ref_gene_id:
                    outfile.write(line)
                    continue

                new_attributes = attributes
                if reference_id:
                    new_transcript_id = reference_id.group(1)
                    new_attributes = re.sub(
                        r'transcript_id "([^"]+)"',
                        f'transcript_id "{new_transcript_id}"',
                        new_attributes,
                    )
                if ref_gene_id:
                    new_gene_id = ref_gene_id.group(1)
                    new_attributes = re.sub(
                        r'gene_id "([^"]+)"', f'gene_id "{new_gene_id}"', new_attributes
                    )

                # Ersetze die Attribute und setze die Zeile zusammen
                fields[8] = new_attributes
                outfile.write("\t".join(fields) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Path to the input file")
    parser.add_argument("output", help="Path to the output file")
    args = parser.parse_args()

    modify_gtf(args.input, args.output)
