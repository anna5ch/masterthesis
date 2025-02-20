import argparse
import random
import os


def filter_and_write_gtf(
    file_path, output_path, percentage, region, keep_ids_file=None, path_kept_ids=None
):
    transcripts_in_region = {}
    filtered_lines = []
    kept_transcripts = set()
    region_variants = [region, f"chr{region}"] if region != "all" else None
    all_other_transcripts = {}

    # If a keep_ids_file is provided, load the IDs to keep
    if keep_ids_file:
        try:
            with open(keep_ids_file, "r") as file:
                keep_ids = set(
                    line.strip() for line in file if not line.startswith("#")
                )
        except TypeError:
            keep_ids = []
    else:
        keep_ids = None

    with open(file_path, "r") as file:
        for line in file:
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            info = parts[8]
            chromosome = parts[0]
            transcript_id = None
            gene_id = None
            for attr in info.split(";"):
                attr = attr.strip()
                if "gene_id" in attr:
                    gene_id = attr.split('"')[1]
                if "transcript_id" in attr:
                    transcript_id = attr.split('"')[1]
                if transcript_id is not None and gene_id is not None:
                    break
            if (
                transcript_id is None
                and gene_id.startswith("ERCC")
                or chromosome not in [f"chr{i}" for i in range(1, 23)]
            ):
                # or gene_id not in [f"chr{i}" for i in range(1, 23)]
                continue  # Skip this transcript
            if transcript_id:
                if region == "all" or chromosome in region_variants:
                    if chromosome not in transcripts_in_region:
                        transcripts_in_region[chromosome] = []
                    transcripts_in_region[chromosome].append((transcript_id, line))
                else:
                    if transcript_id not in all_other_transcripts:
                        all_other_transcripts[transcript_id] = []
                    all_other_transcripts[transcript_id].append(line)
                    filtered_lines.append(line)
            else:
                filtered_lines.append(line)

    if keep_ids:
        # Keep only the transcripts listed in keep_ids
        for chromosome, transcripts in transcripts_in_region.items():
            for transcript_id, line in transcripts:
                if transcript_id in keep_ids:
                    filtered_lines.append(line)
                    kept_transcripts.add(transcript_id)
    else:
        # Select a subset of transcripts randomly based on the percentage
        unique_transcript_ids = set()

        for transcripts in transcripts_in_region.values():
            for transcript_id, _ in transcripts:
                unique_transcript_ids.add(transcript_id)

        # total_transcripts = len(unique_transcript_ids)
        # num_to_remove = int(total_transcripts * percentage)
        # num_per_chromosome = num_to_remove // len(transcripts_in_region)
        # Uniformly remove transcripts across chromosomes
        num_removed = 0
        for chromosome, transcripts in transcripts_in_region.items():
            unique_ids = list(set([i[0] for i in transcripts]))
            num_to_remove = int(len(unique_ids) * percentage)
            transcripts_to_remove = random.sample(
                unique_ids, min(num_to_remove, len(unique_ids))
            )
            num_removed += len(transcripts_to_remove)
            for transcript_id, line in transcripts:
                if transcript_id not in transcripts_to_remove:
                    filtered_lines.append(line)
                    kept_transcripts.add(transcript_id)

        # for transcript_id in all_other_transcripts.keys():
        # add all other kept ids
        #   kept_transcripts.add(transcript_id)

        kept_transcripts = sorted(kept_transcripts)

        # Save the kept transcript IDs to a file if no keep_ids file was provided
        with open(path_kept_ids, "w") as file:
            file.write(f"# Number of removed transcripts: {num_removed}\n")
            for transcript_id in kept_transcripts:
                file.write(f"{transcript_id}\n")

    # Write the output GTF
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as file:
        # Add line with percentage and region passed as an argument
        file.write(f"# Percentage: {percentage}, Region: {region}\n")
        file.writelines(filtered_lines)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Remove part of transcripts from a GTF file."
    )
    parser.add_argument("gtf", type=str, help="The GTF file in question")
    parser.add_argument(
        "percentage", type=float, help="Percentage of transcripts to remove"
    )
    parser.add_argument(
        "region",
        type=str,
        help="Region to remove transcripts ('all' or chromosome number)",
    )
    parser.add_argument("output_file", type=str, help="Output file path")
    parser.add_argument(
        "--keep_ids",
        type=str,
        help="File containing transcript IDs to keep",
    )
    parser.add_argument(
        "--path_kept_ids",
        type=str,
        help="File path to store kept transcript IDs file",
    )
    args = parser.parse_args()

    filter_and_write_gtf(
        args.gtf,
        args.output_file,
        args.percentage,
        args.region,
        args.keep_ids,
        args.path_kept_ids,
    )
