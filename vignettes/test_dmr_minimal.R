#!/usr/bin/env Rscript
# test_dmr_minimal.R - Minimal test to verify DMR analysis works

# Load functions
source("./R/dmr_analysis.R")

# First run the debug script to understand H5 structure
cat("=== Running H5 structure debug ===\n")
source("./R/debug_h5_structure.R")

cat("\n\n=== Testing DMR Analysis Steps ===\n")

# Read data
h5_paths <- fread("./Data/h5_paths.tsv")
metadata <- fread("./Data/cell_metadata.tsv")[, .(cell_id, cluster_id)]

# Test with just one cell first
cat("\n1. Testing with single cell...\n")
test_h5 <- h5_paths[1]
test_meta <- metadata[1]

# Try indexing
single_chr_index <- index_chromosomes(
  h5_paths = test_h5,
  type = "CG",
  chr_list = "chr1",
  threads = 1
)

if (length(single_chr_index) > 0) {
  cat("  Success! Found chromosomes:", names(single_chr_index), "\n")
} else {
  cat("  Failed to index chromosomes\n")
}

# Test with all cells
cat("\n2. Testing chromosome indexing with all cells...\n")
chr_index <- index_chromosomes(
  h5_paths = h5_paths,
  type = "CG",
  chr_list = c("chr1", "chr2"),  # Just test two chromosomes
  threads = 1
)

cat("Chromosomes indexed:", names(chr_index), "\n")
cat("Cells per chromosome:\n")
for (chr in names(chr_index)) {
  cat(sprintf("  %s: %d cells\n", chr, nrow(chr_index[[chr]])))
}

# Show index structure
cat("\nIndex structure for chr1:\n")
print(head(chr_index[["chr1"]]))

# Step 2: Test window calculation
cat("\n3. Testing window calculation...\n")
windows <- calc_smoothed_windows(
  h5_paths = h5_paths,
  chr_index = chr_index,
  metadata = metadata,
  type = "CG",
  step = 5000,  # Larger windows for testing
  smooth = 1,   # No smoothing for simplicity
  group_by = "cluster_id",
  genome = "hg38",
  threads = 1
)

cat("Window matrix dimensions:", dim(windows$sum_matrix), "\n")
cat("Column names:", names(windows$sum_matrix), "\n")

# Show a few rows
cat("\nFirst few rows of count matrix:\n")
print(head(windows$sum_matrix))

# Check for data
if (nrow(windows$sum_matrix) == 0) {
  cat("\nERROR: No windows generated. Debugging...\n")
  
  # Check if cells have data
  test_path <- h5_paths$path[1]
  test_barcode <- h5_paths$barcode[1]
  test_data <- h5read(test_path, paste0("CG/", test_barcode, "/1"))
  
  cat("Sample data from first cell:\n")
  print(head(test_data))
  
  cat("\nChromosomes in data:", unique(test_data$chr), "\n")
  
  stop("No windows generated. Check chromosome names and data structure.")
}

# Step 3: Test DMR detection
cat("\n4. Testing DMR detection...\n")
dmr_results <- test_dmr(
  sum_matrix = windows$sum_matrix,
  comparisons = NULL,
  min_total = 5,
  min_group = 2
)

cat("DMR test columns added:", grep("_pval|_logFC", names(dmr_results), value = TRUE), "\n")

# Check for significant results
pval_cols <- grep("_pval$", names(dmr_results), value = TRUE)
for (col in pval_cols) {
  sig_count <- sum(dmr_results[[col]] < 0.05, na.rm = TRUE)
  total_tested <- sum(!is.na(dmr_results[[col]]))
  cat(sprintf("  %s: %d/%d significant (p < 0.05)\n", col, sig_count, total_tested))
}

# Show some results
cat("\nExample results (first 10 windows):\n")
result_cols <- c("window", grep("_pval$|_logFC$", names(dmr_results), value = TRUE))
print(dmr_results[1:min(10, nrow(dmr_results)), ..result_cols])

# Step 4: Test filtering
cat("\n5. Testing result filtering...\n")
dmr_filtered <- filter_dmr(
  dmr_matrix = dmr_results,
  method = "none",  # No correction for testing
  p_threshold = 0.05,
  log_threshold = 0.5,
  filter = TRUE
)

cat("Filtered results:", nrow(dmr_filtered), "DMRs\n")

if (nrow(dmr_filtered) > 0) {
  cat("\nExample filtered results:\n")
  print(head(dmr_filtered))
  
  # Step 5: Test collapsing
  cat("\n6. Testing DMR collapsing...\n")
  dmr_collapsed <- collapse_dmr(
    dmr_filtered = dmr_filtered,
    max_gap = 10000,  # 10kb gap
    min_length = 5000  # 5kb minimum
  )
  
  cat("Collapsed DMRs:", nrow(dmr_collapsed), "\n")
  
  if (nrow(dmr_collapsed) > 0) {
    cat("\nExample collapsed DMRs:\n")
    print(head(dmr_collapsed))
  }
} else {
  cat("\nNo significant DMRs found. This might be due to:\n")
  cat("- Small test dataset\n")
  cat("- Large window size (5000bp)\n")
  cat("- Stringent filtering\n")
}

cat("\n=== Test Complete ===\n")

# Check for known DMRs if available
if (file.exists("known_dmrs.bed")) {
  known_dmrs <- fread("known_dmrs.bed", 
                      col.names = c("chr", "start", "end", "name", "score", "strand"))
  
  cat("\nChecking overlap with known DMRs...\n")
  
  # Check which known DMRs are in our test chromosomes
  known_in_test <- known_dmrs[chr %in% c("chr1", "chr2")]
  cat(sprintf("Known DMRs in test chromosomes: %d\n", nrow(known_in_test)))
  
  if (nrow(known_in_test) > 0 & nrow(dmr_filtered) > 0) {
    # Check for overlaps
    for (i in 1:nrow(known_in_test)) {
      dmr <- known_in_test[i]
      nearby <- dmr_filtered[chr == dmr$chr & 
                               abs((start + end)/2 - (dmr$start + dmr$end)/2) < 50000]
      
      if (nrow(nearby) > 0) {
        cat(sprintf("  Found %d DMRs near %s (%s:%d-%d)\n", 
                    nrow(nearby), dmr$name, dmr$chr, dmr$start, dmr$end))
      }
    }
  }
}

# Save test results
saveRDS(list(
  chr_index = chr_index,
  windows = windows,
  dmr_results = dmr_results,
  dmr_filtered = dmr_filtered
), "./Data/test_results.rds")

cat("\nResults saved to test_results.rds\n")
