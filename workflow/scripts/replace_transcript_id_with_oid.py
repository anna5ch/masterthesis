import sys


def replace_transcript_id_with_oid(input_file, output_file):
    transcript_to_oid = {}

    # First pass: Create a mapping from transcript_id to oId
    with open(input_file, "r") as infile:
        for line in infile:
            if line.strip() and line.split("\t")[2] == "transcript":
                attributes = line.split("\t")[8]
                transcript_id = next(
                    (
                        attr.split('"')[1]
                        for attr in attributes.split(";")
                        if attr.strip().startswith("transcript_id")
                    ),
                    None,
                )
                o_id = next(
                    (
                        attr.split('"')[1]
                        for attr in attributes.split(";")
                        if attr.strip().startswith("oId")
                    ),
                    None,
                )
                if transcript_id and o_id:
                    transcript_to_oid[transcript_id] = o_id

    # Second pass: Update the GTF file
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            if line.strip():
                fields = line.split("\t")
                attributes = fields[8]

                # Find the transcript_id in the line
                transcript_id = next(
                    (
                        attr.split('"')[1]
                        for attr in attributes.split(";")
                        if attr.strip().startswith("transcript_id")
                    ),
                    None,
                )

                # Replace transcript_id with the corresponding oId, if it exists
                if transcript_id in transcript_to_oid:
                    new_id = transcript_to_oid[transcript_id]
                    attributes = attributes.replace(
                        f'transcript_id "{transcript_id}"', f'transcript_id "{new_id}"'
                    )

                fields[8] = attributes
                line = "\t".join(fields) + "\n"
                outfile.write(line.rstrip() + "\n")


if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    replace_transcript_id_with_oid(input_file, output_file)
