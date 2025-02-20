import os
from collections import defaultdict


def process_tracking_file(filepath):
    line_count = 0
    unique_ids = set()

    with open(filepath, "r") as f:
        for line in f:
            line_count += 1
            fields = line.strip().split("\t")
            if len(fields) >= 5:
                match_status = fields[3]
                if match_status in ("c", "="):
                    ids = fields[2].split("|")[1]
                    unique_ids.add(ids)

    return line_count, unique_ids


def read_kept_ids(kept_ids_file):
    kept_ids = set()
    with open(kept_ids_file, "r") as file:
        for line in file:
            line = line.strip()
            if not line.startswith("#") and line:  # Skip comment lines and empty lines
                kept_ids.add(line)
    return kept_ids


if __name__ == "__main__":
    base_path = "/home/anna/mas-seq-discovery/workflow/results/filtering"
    kept_ids_file = base_path.replace("filtering", "discovery/kept_ids.csv")
    kept_ids = read_kept_ids(kept_ids_file)
    # "../results/filtering"
    tools = [
        "bambu",
        "isoseq",
        "isoquant",
        "stringtie",
        "scallop",
        "flames",
        "flair",
        "stringtie_long",
    ]
    results = defaultdict(lambda: defaultdict(dict))
    for tool in tools:
        tool_path = os.path.join(base_path, tool)
        if os.path.exists(tool_path):
            for root, _, files in os.walk(tool_path):
                for file in files:
                    if file.endswith(".tracking") and file.startswith("only_new"):
                        # sample name
                        sample = file.split(".")[0].replace("only_new_", "")
                        filepath = os.path.join(root, file)
                        line_count, unique_ids = process_tracking_file(filepath)
                        filtered_unique_ids = unique_ids - kept_ids
                        results[tool][sample]["line_count"] = line_count
                        results[tool][sample]["unique_ids"] = len(filtered_unique_ids)
    for tool, samples in results.items():
        print(f"Tool: {tool}")
        for sample, stats in samples.items():
            print(f"  Sample: {sample}")
            print(f"    Total isoforms: {stats['line_count']}")
            print(f"    Detected IDs (c or =): {stats['unique_ids']}")
