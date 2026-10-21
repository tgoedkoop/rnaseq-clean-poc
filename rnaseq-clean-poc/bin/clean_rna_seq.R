#!/usr/bin/env Rscript
# ==========================================================================
#  clean_rna_seq.R
# --------------------------------------------------------------------------
#  Cleans a differential expression results table (e.g. DESeq2 output),
#  classifies each gene's significance, and draws a volcano plot.
#
#  Required columns: log2FoldChange, padj
#
#  Usage:
#    clean_rna_seq.R <input_results.csv> [padj_cutoff] [lfc_cutoff]
#
#  Defaults: padj_cutoff = 0.05, lfc_cutoff = 1
#
#  Outputs (written to the current working directory):
#    cleaned_de_results.csv   - all genes with complete stats + Significance
#    significant_genes.csv    - only the Up/Downregulated genes
#    volcano_plot.pdf         - volcano plot, y-axis = -log10(padj)
# ==========================================================================

# Load required libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(readr))

# 1. Parse command line arguments ------------------------------------------
# FIX: thresholds are now arguments (with defaults) instead of being
# hard-coded in two places, so Nextflow can pass them in from `params`.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: clean_rna_seq.R <input_results.csv> [padj_cutoff] [lfc_cutoff]",
       call. = FALSE)
}

input_file  <- args[1]
padj_cutoff <- if (length(args) >= 2) as.numeric(args[2]) else 0.05
lfc_cutoff  <- if (length(args) >= 3) as.numeric(args[3]) else 1

# Fail early if a threshold was not a valid number
if (is.na(padj_cutoff) || is.na(lfc_cutoff)) {
  stop("padj_cutoff and lfc_cutoff must be numeric.", call. = FALSE)
}

# FIX: clearer output names. The "cleaned" file holds ALL genes;
# a separate file holds only the significant genes.
output_all  <- "cleaned_de_results.csv"
output_sig  <- "significant_genes.csv"
output_plot <- "volcano_plot.pdf"

cat("Thresholds: padj <", padj_cutoff, "| |log2FC| >", lfc_cutoff, "\n")

# 2. Read the DE results ---------------------------------------------------
cat("Reading input data from:", input_file, "\n")
df <- read_csv(input_file, show_col_types = FALSE)

# Check the required columns exist before doing anything else
required_cols <- c("log2FoldChange", "padj")
missing_cols  <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Input is missing required column(s): ",
       paste(missing_cols, collapse = ", "), call. = FALSE)
}

# 3. Clean and classify ----------------------------------------------------
# Remove genes with missing statistics (e.g. DESeq2 sets padj = NA for
# low-count genes removed by independent filtering), then classify them.
cat("Filtering NAs and classifying significance...\n")
cleaned_df <- df %>%
  filter(!is.na(padj) & !is.na(log2FoldChange)) %>%
  mutate(
    Significance = case_when(
      padj < padj_cutoff & log2FoldChange >  lfc_cutoff ~ "Upregulated",
      padj < padj_cutoff & log2FoldChange < -lfc_cutoff ~ "Downregulated",
      TRUE                                              ~ "Not Significant"
    )
  )

n_removed <- nrow(df) - nrow(cleaned_df)
cat("Removed", n_removed, "of", nrow(df), "rows with missing statistics.\n")
print(table(cleaned_df$Significance))

# 4. Save the tables -------------------------------------------------------
cat("Saving cleaned data to:", output_all, "\n")
write_csv(cleaned_df, output_all)

cat("Saving significant genes to:", output_sig, "\n")
write_csv(filter(cleaned_df, Significance != "Not Significant"), output_sig)

# 5. Volcano plot ----------------------------------------------------------
# FIX: the y-axis now uses padj, the same statistic used for classification,
# so the dashed threshold line matches the point colours exactly.
# padj can be exactly 0 (numerical underflow), which would give
# -log10(0) = Inf and silently drop the point. Replace zeros with the
# smallest non-zero padj in the data so those genes plot at the top of
# the existing scale instead of stretching the axis.
min_nonzero_padj <- min(cleaned_df$padj[cleaned_df$padj > 0])
n_zero <- sum(cleaned_df$padj == 0)
if (n_zero > 0) {
  cat("Note:", n_zero, "gene(s) with padj = 0 plotted at the min non-zero padj.\n")
}
plot_df <- cleaned_df %>%
  mutate(neg_log10_padj = -log10(pmax(padj, min_nonzero_padj)))

volcano <- ggplot(plot_df,
                  aes(x = log2FoldChange, y = neg_log10_padj,
                      color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  # Colour-blind-safe palette (Okabe-Ito)
  scale_color_manual(values = c("Upregulated"     = "#D55E00",
                                "Downregulated"   = "#0072B2",
                                "Not Significant" = "grey80")) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
             linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(padj_cutoff),
             linetype = "dashed", color = "black", alpha = 0.5) +
  theme_minimal() +
  labs(
    title    = "RNA-seq Differential Expression",
    subtitle = paste0("Thresholds: |log2FC| > ", lfc_cutoff,
                      ", padj < ", padj_cutoff),
    x = "log2 Fold Change",
    y = "-log10(adjusted p-value)"
  )

ggsave(output_plot, plot = volcano, width = 8, height = 6)  # PDF is vector: dpi not needed

cat("Script completed successfully.\n")
