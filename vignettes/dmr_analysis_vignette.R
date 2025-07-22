#!/usr/bin/env Rscript
# amethyst_vs_minithyst_comparison.R - Comprehensive comparison of DMR analysis
# Tests that minithyst produces identical results to amethyst at each pipeline stage

# Ensure we're in the right directory
if (!dir.exists("./Data")) {
  stop("Data directory not found. Run generate_toy_data.R first.")
}

# Load required packages
suppressPackageStartupMessages({
  library(data.table)
  library(rhdf5)
  library(GenomicRanges)
  library(IRanges)
})

# Source minithyst functions
# source("R/utils.R")
# source("R/indexing.R")
# source("R/windowing.R")

library(minithyst)

data_dir <- "/Users/colemste/Desktop/placeholder_for_downloads/"

h5_file <- paste0(data_dir, "pbmc_vignette.h5")
# h5_paths_file <- paste0(data_dir, "test_h5_paths.tsv")
# metadata_file <- paste0(data_dir, "test_cell_metadata.tsv")
metadata_file <- paste0(data_dir, "pbmc_vignette_cellInfo.txt")
annotation_file <- paste0(data_dir, "pbmc_vignette.annot")
# download.file("https://adeylabopen.s3.us-west-2.amazonaws.com/amethyst/pbmc_vignette.h5", h5_file, method = "curl") # Contains site-level methylation information for each cell
# download.file("https://raw.githubusercontent.com/lrylaarsdam/amethyst/main/vignettes/pbmc_vignette/pbmc_vignette_cellInfo.txt", metadata_file) 
# download.file("https://raw.githubusercontent.com/lrylaarsdam/amethyst/main/vignettes/pbmc_vignette/pbmc_vignette.annot", annotation_file)

# Load toy data
# h5_paths <- fread(h5_paths_file)
metadata <- fread(metadata_file)
annotation <- fread(annotation_file)

metadata <- merge(metadata, annotation, by="V1")
colnames(metadata) <- c("cell_id", "Coverage", "CG_Cov", "pct", "CH_Cov", "something", "cluster_id")
h5_paths <- data.table(barcode = metadata$cell_id, path=h5_file)

cat("Data loaded:\n")
cat(sprintf("- H5 files: %d\n", nrow(h5_paths)))
cat(sprintf("- Cells: %d\n", nrow(metadata)))
cat(sprintf("- Cell types: %s\n", paste(unique(metadata$cluster_id), collapse = ", ")))

# Set deterministic parameters
set.seed(42)
setDTthreads(12)

cat("=== Minithyst vignette ===\n")

#===============================================================================
# SETUP: Load data
#===============================================================================


# Test chromosomes (subset for speed)
test_chrs <- c("chr1", "chr2", "chr3")

# Common analysis parameters
CONTEXT <- "CG"
WINDOW_SIZE <- 500
SMOOTHING_WINDOWS <- 3
MIN_COVERAGE <- 5
GENOME <- "hg38"

cat("\nCommon parameters:\n")
cat(sprintf("- Context: %s\n", CONTEXT))
cat(sprintf("- Window size: %d bp\n", WINDOW_SIZE))
cat(sprintf("- Smoothing windows: %d\n", SMOOTHING_WINDOWS))
cat(sprintf("- Min coverage: %d\n", MIN_COVERAGE))

#===============================================================================
# STAGE 1: Chromosome Indexing
#===============================================================================

cat("\n=== STAGE 1: CHROMOSOME INDEXING COMPARISON ===\n")

# Minithyst indexing
cat("Running minithyst indexing...\n")
mini_chr_index <- create_optimized_index(
  h5_paths = h5_paths,
  type = CONTEXT,
  chr_list = test_chrs,
  threads = 1
)

#===============================================================================
# STAGE 2: Window Smoothing
#===============================================================================

cat("\n=== STAGE 2: WINDOW SMOOTHING ===\n")

# Minithyst window calculation
cat("Running minithyst window calculation...\n")
t_0_mini <- Sys.time()
mini_windows <- calc_smoothed_windows(
  h5_paths = h5_paths,
  chr_index = mini_chr_index,
  metadata = metadata[, .(cell_id, cluster_id)],
  type = CONTEXT,
  step = WINDOW_SIZE,
  smooth = SMOOTHING_WINDOWS,
  group_by = "cluster_id",
  genome = GENOME,
  chrList = test_chrs,
  chrSizes = c(248956422, 242193529, 198295559),
  threads = 1
)
t_1_mini <- Sys.time()
t_delta_mini <- t_1_mini - t_0_mini

#===============================================================================
# STAGE 3: DMR Testing
#===============================================================================

cat("\n=== STAGE 3: DMR TESTING ===\n")

# Minithyst DMR testing
cat("Running minithyst DMR testing...\n")
t_0 <- Sys.time()
mini_dmr_results <- test_dmr(
  sum_matrix = mini_windows$sum_matrix,
  comparisons = NULL,  # All vs all
  min_total = MIN_COVERAGE * 2,
  min_group = MIN_COVERAGE
)
t_1 <- Sys.time()
delta_t_original <- t_1 - t_0
# mini_dmr_results_og <- mini_dmr_results  
mini_dmr_results

#===============================================================================
# STAGE 4: DMR Filtering
#===============================================================================

cat("\n=== STAGE 4: DMR FILTERING ===\n")

# Minithyst filtering
cat("Running minithyst DMR filtering...\n")
mini_dmr_filtered <- filter_dmr(
  dmr_matrix = mini_dmr_results,
  method = "bonferroni",
  p_threshold = 0.05,  # Lenient for toy data
  log_threshold = 1.0,
  filter = TRUE
)

#===============================================================================
# STAGE 5: DMR Collapsing
#===============================================================================

cat("\n=== STAGE 5: DMR COLLAPSING ===\n")

# Minithyst collapsing
cat("Running minithyst DMR collapsing...\n")
mini_dmr_collapsed <- collapse_dmr(
  dmr_filtered = mini_dmr_filtered,
  max_gap = 1000,
  min_length = 1000
)
