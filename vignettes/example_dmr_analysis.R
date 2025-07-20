#!/usr/bin/env Rscript
# example_dmr_analysis.R - Example of running simplified DMR analysis

# Load required functions
source("./R/dmr_analysis.R")
source("./R/dmr_utils.R")

# Set up paths and metadata
# These files are created by generate_toy_data.R
h5_paths <- fread("./Data/h5_paths.tsv")
metadata <- fread("./Data/cell_metadata.tsv")

# For the DMR analysis, we'll use cluster_id as the grouping variable
metadata <- metadata[, .(cell_id, cluster_id, sample)]

# Configure analysis
config <- configure_dmr_analysis(
  genome = "hg38",
  context = "CG",
  window_size = 500,
  smoothing_windows = 3,
  min_coverage = 5
)

print(config)

# Option 1: Run complete pipeline
# This runs all steps and saves intermediate results
results <- run_dmr_pipeline(
  h5_paths = h5_paths,
  metadata = metadata,
  config = config,
  output_dir = "dmr_results",
  threads = 4
)

# Option 2: Run step by step for more control
# This allows inspection and modification at each step

# Step 1: Index chromosomes
chr_index <- index_chromosomes(
  h5_paths = h5_paths,
  type = "CG",
  chr_list = paste0("chr", c(1:22, "X")),  # Exclude Y and mitochondrial
  threads = 4
)

# Check indexing results
sapply(chr_index, nrow)

# Step 2: Calculate smoothed windows
windows <- calc_smoothed_windows(
  h5_paths = h5_paths,
  chr_index = chr_index,
  metadata = metadata,
  type = "CG",
  step = 500,
  smooth = 3,
  group_by = "cluster_id",
  genome = "hg38",
  threads = 4
)

# Inspect window counts
dim(windows$sum_matrix)
dim(windows$pct_matrix)

# Step 3: Test for DMRs
# Option A: All vs all (each group against all others)
dmr_results <- test_dmr(
  sum_matrix = windows$sum_matrix,
  comparisons = NULL,  # NULL means all vs all
  min_total = 10,
  min_group = 5
)

# Option B: Specific comparisons
comparisons <- data.frame(
  name = c("T_vs_B", "sample_effect"),
  A = c("T_cell", "T_cell,B_cell"),
  B = c("B_cell", "T_cell,B_cell"),
  stringsAsFactors = FALSE
)

dmr_results_custom <- test_dmr(
  sum_matrix = windows$sum_matrix,
  comparisons = comparisons,
  min_total = 10,
  min_group = 5
)

# Step 4: Filter results
dmr_filtered <- filter_dmr(
  dmr_matrix = dmr_results,
  method = "bonferroni",  # Conservative for genome-wide testing
  p_threshold = 0.01,
  log_threshold = 1.5,    # Minimum 1.5 log2 fold change
  filter = TRUE
)

# Check filtering results
nrow(dmr_filtered)
table(dmr_filtered$test, dmr_filtered$direction)

# Step 5: Collapse adjacent DMRs
dmr_collapsed <- collapse_dmr(
  dmr_filtered = dmr_filtered,
  max_gap = 2000,      # Merge DMRs within 2kb
  min_length = 2000    # Minimum DMR size of 2kb
)

# Inspect collapsed results
summary_stats <- summarize_dmr(dmr_collapsed)
print(summary_stats)

# Step 6: Annotate with genes (optional)
# Requires GTF file
if (file.exists("gencode.v43.annotation.gtf.gz")) {
  dmr_annotated <- annotate_dmr(
    dmr_collapsed = dmr_collapsed,
    gtf_file = "gencode.v43.annotation.gtf.gz",
    feature_type = "gene"
  )
  
  # Show DMRs with gene annotations
  dmr_annotated[!is.na(gene_names)][1:10]
}

# Step 7: Export results
# Export as BED files for genome browser
export_dmr_bed(
  dmr_results = dmr_collapsed,
  output_file = "dmrs_all.bed",
  track_name = "CG_DMRs",
  separate_by = "direction"  # Creates separate files for hyper/hypo
)

# Create tracks for specific groups
for (group in unique(metadata$cluster_id)) {
  track <- create_dmr_track(windows, group)
  
  # Write UCSC genome browser track
  writeLines(track$header, paste0(group, "_methylation.bedGraph"))
  fwrite(track$data, paste0(group, "_methylation.bedGraph"), 
         sep = "\t", col.names = FALSE, append = TRUE)
}

# Advanced usage examples

# 1. Filter by effect size and significance
strong_dmrs <- dmr_collapsed[abs(dmr_logFC) > 2 & dmr_padj < 0.001]

# 2. Find DMRs specific to one comparison
t_cell_specific <- dmr_collapsed[test == "T_cell" & direction == "hyper"]

# 3. Extract methylation values for specific DMRs
dmr_of_interest <- dmr_collapsed[1]  # First DMR
window_coords <- with(dmr_of_interest, {
  paste0(chr, "_", 
         seq(dmr_start, dmr_end - 500, by = 500), "_",
         seq(dmr_start + 500, dmr_end, by = 500))
})

# Get methylation values for these windows
meth_values <- windows$pct_matrix[rownames(windows$pct_matrix) %in% window_coords, ]

# 4. Compare DMRs between different analyses
# For example, CG vs CH methylation
if (exists("ch_dmr_results")) {
  # Find overlapping DMRs
  cg_ranges <- GRanges(dmr_collapsed$chr, 
                       IRanges(dmr_collapsed$dmr_start, dmr_collapsed$dmr_end))
  ch_ranges <- GRanges(ch_dmr_results$chr,
                       IRanges(ch_dmr_results$dmr_start, ch_dmr_results$dmr_end))
  
  overlaps <- findOverlaps(cg_ranges, ch_ranges)
  
  # DMRs with both CG and CH changes
  dual_dmrs <- unique(queryHits(overlaps))
  dmr_collapsed[dual_dmrs]
}

# 5. Generate summary report
cat("\n=== DMR Analysis Summary ===\n")
cat(sprintf("Total windows tested: %d\n", nrow(windows$sum_matrix)))
cat(sprintf("Significant DMRs: %d\n", nrow(dmr_filtered)))
cat(sprintf("After collapsing: %d\n", nrow(dmr_collapsed)))
cat("\nBy direction:\n")
print(table(dmr_collapsed$direction))
cat("\nBy test:\n")
print(table(dmr_collapsed$test))
cat(sprintf("\nMedian DMR size: %.0f bp\n", median(dmr_collapsed$dmr_length)))
cat(sprintf("Total genomic coverage: %.1f Mb\n", 
            sum(dmr_collapsed$dmr_length) / 1e6))