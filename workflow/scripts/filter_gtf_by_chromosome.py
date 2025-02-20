import sys


def load_ids_from_file(ids_file):
    with open(ids_file, "r") as file:
        return {line.strip() for line in file}


def filter_gtf_by_chromosome(input_file, output_file, chromosomes):
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
            else:
                chrom = line.split("\t")[0]
                if chrom in chromosomes:
                    outfile.write(line)


def filter_gtf_by_ids(input_file, output_file, kept_ids):
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
            else:
                attributes = line.split("\t")[8]
                transcript_id = next(
                    (
                        attr.split('"')[1]
                        for attr in attributes.split(";")
                        if attr.strip().startswith("transcript_id")
                    ),
                    None,
                )
                if transcript_id not in kept_ids:
                    outfile.write(line)


if __name__ == "__main__":

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Determine if we're filtering by chromosome or IDs
    if sys.argv[3] == "--chromosomes":
        chromosomes = sys.argv[4:]
        filter_gtf_by_chromosome(input_file, output_file, chromosomes)
    elif sys.argv[3] == "--kept_ids":
        ids_file = sys.argv[4]
        kept_ids = load_ids_from_file(ids_file)
        filter_gtf_by_ids(input_file, output_file, kept_ids)
    else:
        print("Usage:")
        print("  To filter by chromosome:")
        print(
            "    python filter_gtf_by_chromosome.py <input_file> <output_file> --chromosomes <chromosome1> <chromosome2> ..."
        )
        print("  To filter by IDs:")
        print(
            "    python filter_gtf_by_chromosome.py <input_file> <output_file> --kept_ids <ids_file>"
        )
