import pysam
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np   
import pandas as pd


def get_bam_metrics(bam_file):
    samfile = pysam.AlignmentFile(bam_file, "rb")
    # Get the total number of reads in the BAM file
    num_reads = samfile.count()

    # Get read depth (average coverage)
    depth = []
    for pileupcolumn in samfile.pileup():
        depth.append(pileupcolumn.n)
    avg_depth = np.mean(depth) if depth else 0
    samfile.close()
    return num_reads, avg_depth


def plot_read_depth_and_reads(bam_files, output_file):
    # Store the results for each BAM file
    results = []
    for bam_file in bam_files:
        num_reads, avg_depth = get_bam_metrics(bam_file)
        bam_file = bam_file.split("/")[-1].split(".")[0]
        results.append(
            {
                "BAM File": bam_file,
                "Number of Reads": num_reads,
                "Average Depth": avg_depth,
            }
        )

    df = pd.DataFrame(results)

    # Set up the plot
    fig, axes = plt.subplots(1, 2, figsize=(18, 8))

    # Bar plot for Number of Reads
    sns.barplot(
        data=df, x="BAM File", y="Number of Reads", ax=axes[0], color="darkolivegreen"
    )
    axes[0].set_title("Number of Reads per Sample File")
    axes[0].set_xlabel("Sample")
    axes[0].set_ylabel("Number of Reads")
    axes[0].tick_params(axis="x", rotation=45)

    # Bar plot for Average Depth
    sns.barplot(data=df, x="BAM File", y="Average Depth", ax=axes[1], color="bisque")
    axes[1].set_title("Average Read Depth per Sample File")
    axes[1].set_xlabel("Sample")
    axes[1].set_ylabel("Average Depth")
    axes[1].tick_params(axis="x", rotation=45)

    # Adjust layout and save
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()


# Example usage:
if __name__ == "__main__":
    bam_files = snakemake.input
    output_file = snakemake.output[0]

    plot_read_depth_and_reads(bam_files, output_file)
