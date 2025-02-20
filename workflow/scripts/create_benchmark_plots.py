import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from collections import defaultdict


def read_benchmark_file(filepath):
    # Extract sample and tool name from the filename
    filename = os.path.basename(filepath)
    sample_name = filename.split(".")[0]  # Get the sample from the filename
    tool_name = filename.split(".")[1]  # Get the tool name

    # Read the file into a Pandas DataFrame
    df = pd.read_csv(
        filepath,
        sep="\t",
        skiprows=1,
        names=[
            "s",
            "h:m:s",
            "max_rss",
            "max_vms",
            "max_uss",
            "max_pss",
            "io_in",
            "io_out",
            "mean_load",
            "cpu_time",
        ],
    )
    df["tool"] = tool_name  # Add tool name as a column
    df["sample"] = sample_name  # Add sample name as a column
    return df


def plot_combined_benchmark(df, sample_name, output_file):
    sns.set(style="whitegrid")

    # Define the metrics to plot
    metrics = {
        "s": "Duration (s)",
        "max_rss": "Max RSS (MB)",
        "cpu_time": "CPU Time (s)",
    }

    # Create a figure with subplots for each metric
    fig, axes = plt.subplots(1, len(metrics), figsize=(18, 6), sharey=False)
    fig.suptitle(f"Benchmark Analysis for Sample: {sample_name}", fontsize=16)

    for ax, (metric, title) in zip(axes, metrics.items()):
        sns.boxplot(
            data=df, x="tool", y=metric, ax=ax, palette=sns.color_palette("dark:#5A9_r")
        )
        # sns.stripplot(data=df, x="tool", y=metric, ax=ax, color=".25", jitter=True, size=5)  # Add individual points
        ax.set_title(title)
        ax.set_xlabel("Tool")
        ax.set_ylabel(title)

    # Adjust layout and save the figure
    plt.tight_layout(rect=[0, 0, 1, 0.96])  # Leave space for the title
    plt.savefig(output_file)
    plt.close()


if __name__ == "__main__":
    # Snakemake input files
    input_files = snakemake.input
    output_path = snakemake.params[0]

    # Group input files by sample name using a hash map
    sample_files = defaultdict(list)
    for filepath in input_files:
        sample_name = os.path.basename(filepath).split(".")[
            0
        ]  # Extract sample name from filename
        sample_files[sample_name].append(filepath)

    # Process each sample and generate plots
    for sample_name, file_list in sample_files.items():
        # Read and combine data for the sample
        dfs = [read_benchmark_file(filepath) for filepath in file_list]
        combined_df = pd.concat(dfs, ignore_index=True)

        # Output file path
        output_file = os.path.join(output_path, f"{sample_name}_benchmark_plots.png")
        # Create combined box plots for the sample
        plot_combined_benchmark(combined_df, sample_name, output_file)
