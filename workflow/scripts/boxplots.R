library(ggplot2)
library(dplyr)
library(reshape2)


calculate_f1 <- function(sens, prec) {
  # calculate F1 score from sensitivity and precision from gffcmp
  if (sens + prec == 0) {
    return(0)
  }
  return(2 * ((sens * prec) / (sens + prec)))
}

extract_metrics <- function(file) {
  # extract sensitivity and precision on transcript level
  lines <- readLines(file)

  # Get tool and sample from the file name
  file_parts <- unlist(strsplit(basename(file), split = "[._]"))
  sample <- file_parts[3]
  tool <- file_parts[4]

  # Look for the Intron chain or Transcript level line in the file
  for (line in lines) {
    if (grepl("Transcript level", line)) {
      # Extract Sensitivity and Precision
      parts <- unlist(strsplit(line, split = "[|]"))
      sens <- as.numeric(trimws(unlist(strsplit(parts[1], ":"))[2]))
      prec <- as.numeric(trimws(parts[2]))

      # Calculate F1 score
      f1 <- calculate_f1(sens / 100, prec / 100)

      # Add the data to the metrics_data data frame
      metrics_data <<- rbind(metrics_data, data.frame(Tool = tool, Sample = sample, F1 = f1, Sensitivity = sens, Precision = prec, stringsAsFactors = FALSE))
    }
  }
}

plot_prec_sens_bar <- function(data) {
  avg_metrics <- data %>%
    group_by(Tool) %>%
    summarise(
      Avg_Sensitivity = mean(Sensitivity),
      Avg_Precision = mean(Precision)
    )

  long_data <- reshape2::melt(avg_metrics,
    id.vars = "Tool",
    measure.vars = c("Avg_Sensitivity", "Avg_Precision"),
    variable.name = "Metric",
    value.name = "Value"
  )

  # Create the bar plot
  bar_plot <- ggplot(long_data, aes(x = Tool, y = Value, fill = Metric)) +
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    theme_minimal(base_size = 15) +
    labs(
      title = "Average Sensitivity and Precision by Tool",
      x = "Tool", y = "Percentage (%)"
    ) +
    scale_fill_manual(
      values = c("Avg_Sensitivity" = "darkolivegreen", "Avg_Precision" = "bisque3"),
      name = "Metric",
      labels = c("Avg Sensitivity", "Avg Precision")
    ) +
    theme(
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      panel.grid.major = element_line(color = "grey80"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  ggsave("results/visualize/precision_sensitivity_barplot.png", plot = bar_plot, width = 10, height = 6, dpi = 300)
}

plot_boxplot_F1 <- function(metrics_data) {
  p <- ggplot(metrics_data, aes(x = Tool, y = F1)) +
    geom_boxplot(fill = "darkolivegreen", color = "black") +
    theme_minimal(base_size = 15) +
    labs(title = "F1 Scores Across Tools (Transcript Level)", x = "Tool", y = "F1 Score") +
    theme(
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      panel.grid.major = element_line(color = "grey80"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  ggsave("results/visualize/boxplot.png", plot = p, width = 10, height = 6, dpi = 300)
}


# input files from snakemake
input_files <- snakemake@input
metrics_data <- data.frame(Tool = character(), Sample = character(), F1 = double(), Sensitivity = double(), Precision = double(), stringsAsFactors = FALSE)

# Loop through all the input files and extract metrics
for (file in input_files) {
  extract_metrics(file)
}

# Plot the boxplot for F1 scores
plot_boxplot_F1(metrics_data)

# plot the sens and prec as barplot
plot_prec_sens_bar(metrics_data)
