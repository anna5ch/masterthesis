library(FLAMES)


args <- commandArgs(trailingOnly = TRUE)

transcriptome_gtf <- args[1]
genome <- args[2]
bam_file <- args[3]
outdir <- args[4]
config_file <- args[5]

config <- jsonlite::fromJSON(config_file)

dir.create(outdir)
find_isoform(
  annotation = transcriptome_gtf,
  genome_fa = genome,
  genome_bam = bam_file,
  outdir = outdir,
  config = config
)
