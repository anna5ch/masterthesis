library(bambu)

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 3 || length(args) > 4) {
    stop("Usage: Rscript discovery_bambu.R <reads.bam> <annotations.gtf> <genome.fasta> [output_dir]")
  }

  reads_bam <- args[1]
  annotations <- args[2]
  genome <- args[3]

  if (length(args) == 4) {
    output_dir <- args[4]
    if (!file.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
  } else {
    output_dir <- getwd() # Default to current working directory
  }


  ncore <- 2

  output_file <- sub(".sorted.bam", ".bambu_output.gtf", basename(reads_bam))
  if (length(args) == 4) {
    output_file <- file.path(output_dir, output_file)
  }
  if (annotations != "NULL") {
    bambuAnnotations <- prepareAnnotations(annotations)
  } else {
    bambuAnnotations <- NULL
  }
  se.discoveryOnly <- bambu(reads = reads_bam, annotations = bambuAnnotations, genome = genome, NDR = 1, quant = FALSE, lowMemory = TRUE, opt.discovery = list(min.readCount = 0, min.readFractionByGene = 0, min.exonDistance = 0))
  writeToGTF(se.discoveryOnly, output_file)
}

main()
