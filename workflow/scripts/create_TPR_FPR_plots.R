#!/usr/bin/env Rscript
library(ggplot2)
library(gridExtra)
library(pROC)
library(viridis)
library(PRROC)
library(RColorBrewer)
library(dplyr)


parse_tracking_file <- function(tracking_file) {
  # read tracking file from gffcompare and extract y_true values
  data <- read.table(tracking_file, header = FALSE, stringsAsFactors = FALSE)
  transcript_id <- gsub(".*\\|(.*?)\\|.*", "\\1", data$V5)
  match_status <- data$V4
  y_true <- ifelse(match_status == "=", 1, 0)
  return(data.frame(transcript_id, y_true))
}

read_feature_matrix <- function(feature_matrix_file) {
  # read feature matrix file
  feature_matrix <- read.table(feature_matrix_file, header = TRUE, sep = "\t")
  return(feature_matrix)
}

extract_sample_tool <- function(filename) {
  # get sample name and tool name of file path
  tracking_pattern <- "(?:salmon|oarfish)_gffcmp_(.*)\\.(.*)\\.tracking"
  feature_matrix_pattern <- "(?:salmon|oarfish)_(.*)\\.(.*)_output_feature_matrix\\.tsv"
  base_filename <- basename(filename)

  if (grepl(tracking_pattern, base_filename)) {
    matches <- regmatches(base_filename, regexec(tracking_pattern, base_filename))
    sample <- matches[[1]][2]
    tool <- matches[[1]][3]
    return(list(sample = sample, tool = tool))
  } else if (grepl(feature_matrix_pattern, base_filename)) {
    matches <- regmatches(base_filename, regexec(feature_matrix_pattern, base_filename))
    sample <- matches[[1]][2]
    tool <- matches[[1]][3]
    return(list(sample = sample, tool = tool))
  }
}

generate_plots <- function(tracking_files, feature_matrix_files, num_removed) {
  sample_names <- unique(sapply(tracking_files, function(f) extract_sample_tool(f)$sample))

  for (sample in sample_names) {
    plots <- list()
    pr_plots <- list()
    # Get the tracking and feature matrix files for this sample
    sample_tracking_files <- tracking_files[grepl(sample, tracking_files)]
    sample_feature_matrix_files <- feature_matrix_files[grepl(sample, feature_matrix_files)]

    for (i in seq_along(sample_tracking_files)) {
      tracking_file <- sample_tracking_files[i]
      feature_matrix_file <- sample_feature_matrix_files[i]

      tracking_info <- extract_sample_tool(tracking_file)
      feature_info <- extract_sample_tool(feature_matrix_file)
      if (tracking_info$sample != feature_info$sample ||
        tracking_info$tool != feature_info$tool) {
        stop(paste(
          "Mismatch between tracking file", tracking_file,
          "and feature matrix file", feature_matrix_file
        ))
      }

      tracking_data <- parse_tracking_file(tracking_file)
      feature_matrix <- read_feature_matrix(feature_matrix_file)

      # Merge data
      merged_data <- merge(tracking_data, feature_matrix, by = "transcript_id")
      # Initialize a plot for this tool
      tool_plot <- ggplot() +
        theme_minimal() +
        theme(panel.background = element_rect(fill = "lightgrey"))

      tool_pr_plot <- ggplot() +
        theme_minimal() +
        theme(panel.background = element_rect(fill = "lightgrey"))
      # Create ROC curves for each feature
      features <- names(feature_matrix)[-1] # Exclude 'transcript_id'

      for (feature in features) {
        if (feature != "label" && feature != "y_true" && feature != "exonic_length" && feature != "gc_content") {
          if (length(unique(merged_data$y_true)) == 2) {
            roc_obj <- roc(merged_data$y_true, merged_data[[feature]], quiet = TRUE)
            # df so several lines are in one plot
            roc_data <- data.frame(
              fpr = 1 - roc_obj$specificities,
              tpr = roc_obj$sensitivities,
              feature = feature
            )

            # Add the ROC line for the current feature
            tool_plot <- tool_plot +
              geom_line(data = roc_data, aes(x = fpr, y = tpr, color = feature), linewidth = 0.5) +
              labs(
                title = paste("ROC Curve for Tool:", tracking_info$tool),
                x = "False Positive Rate", y = "True Positive Rate"
              )

            # Precision-Recall curve
            pr_obj <- pr.curve(
              scores.class0 = merged_data[[feature]],
              weights.class0 = as.integer(merged_data$y_true == 1),
              curve = TRUE
            )
            pr_data <- data.frame(
              recall = pr_obj$curve[, 1],
              precision = pr_obj$curve[, 2],
              feature = feature
            )
            # pr_data$recall <- pr_data$recall * sum(feature_matrix$label)/num_removed

            pr_data <- pr_data %>%
              group_by(recall, feature) %>%
              summarise(precision = mean(precision), .groups = "drop")

            tool_pr_plot <- tool_pr_plot +
              geom_line(data = pr_data, aes(x = recall, y = precision, color = feature), linewidth = 0.5) +
              labs(
                title = paste("Precision-Recall Curve for Tool:", tracking_info$tool),
                x = "Recall", y = "Precision",
              ) + ylim(0, NA)
          }
        }
      }

      tool_plot <- tool_plot + scale_color_viridis_d(name = "Feature", option = "cividis")
      tool_pr_plot <- tool_pr_plot + scale_color_viridis_d(name = "Feature", option = "cividis")

      # Store the plot for this tool
      plots[[length(plots) + 1]] <- tool_plot
      pr_plots[[length(pr_plots) + 1]] <- tool_pr_plot
    }

    ncol <- 2
    width <- ncol * 5
    height <- ceiling(length(plots) / ncol) * 3
    # Combine the tool plots in one file for this sample
    output_file <- paste0("results/visualize/tpr_fpr_plot_", sample, ".png")
    ggsave(output_file, do.call(grid.arrange, c(plots, ncol = 2)), width = width, height = height)

    pr_output_file <- paste0("results/visualize/prec_recall_plot_", sample, ".png")
    ggsave(pr_output_file, do.call(grid.arrange, c(pr_plots, ncol = 2)), width = width, height = height)
  }
}


tracking_file <- c(snakemake@input$tracking_files, snakemake@input$tracking_files2)
feature_matrix_file <- c(snakemake@input$feature_matrix_files, snakemake@input$feature_matrix_files2)
num_removed <- snakemake@params$num_removed

# Generate plots
generate_plots(tracking_file, feature_matrix_file, num_removed)
