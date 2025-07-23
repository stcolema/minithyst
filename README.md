# minithyst

An R package for identifying differentially methylated regions (DMRs) from single-cell methylation data. This package implements a subset of the `amethyst` package's functionality, focusing on DMR analysis.

## Overview

minithyst provides functions for DMR analysis using a window-based approach with Fisher's exact tests.

## Method

### Approach
- Genomic regions are divided into fixed-size windows (configurable, default 500bp)
- Optional spatial smoothing by aggregating adjacent windows
- Statistical testing using Fisher's exact test for differential methylation between groups
- Multiple testing correction (user-configurable method)
- Adjacent significant windows can be merged into regions

### Implementation Notes
- Uses HDF5 for data storage and retrieval
- Processes data chromosome-by-chromosome
- Supports parallel processing where applicable

## Installation

```r
# Install dependencies
install.packages(c("data.table", "rhdf5", "GenomicRanges", "IRanges"))

# Install from source
devtools::install_github("stcolema/minithyst")
```

## Input Data Format

### HDF5 Structure
Expected HDF5 file organization:
  ```
/CG/
  /barcode1/
  /1/  -> data with columns: chr, pos, c, t
/barcode2/
  /1/  -> data with columns: chr, pos, c, t
```

Required columns:
  - `chr`: chromosome identifier (character)
- `pos`: genomic position (integer)
- `c`: methylated read count (integer) 
- `t`: unmethylated read count (integer)

### Metadata
data.table with columns:
  - `cell_id`: must match HDF5 barcodes
- `cluster_id`: grouping variable for comparisons

## Step-by-Step Analysis

```r
library(minithyst)

# Index HDF5 files for efficient access
chr_index <- index_chromosomes(
  h5_paths = h5_paths,
  type = "CG",
  chr_list = paste0("chr", 1:22),
  threads = 1
)

# Calculate methylation in genomic windows
windows <- calc_smoothed_windows(
  h5_paths = h5_paths,
  chr_index = chr_index, 
  metadata = metadata,
  type = "CG",
  step = 500,
  smooth = 3,
  group_by = "cluster_id",
  genome = "hg38"
)

# Test for differential methylation
dmr_results <- test_dmr(
  sum_matrix = windows$sum_matrix,
  comparisons = NULL,  # NULL performs all vs all comparisons
  min_total = 10,
  min_group = 5
)

# Filter results
dmr_filtered <- filter_dmr(
  dmr_matrix = dmr_results,
  method = "bonferroni",
  p_threshold = 0.01,
  log_threshold = 1.5,
  filter = TRUE
)

# Merge adjacent windows
dmr_collapsed <- collapse_dmr(
  dmr_filtered = dmr_filtered,
  max_gap = 2000,
  min_length = 2000
)
```

## Custom Comparisons

```r
comparisons <- data.frame(
  name = c("T_vs_B"),
  A = c("T_cell"),
  B = c("B_cell"),
  stringsAsFactors = FALSE
)

dmr_results <- test_dmr(
  sum_matrix = windows$sum_matrix,
  comparisons = comparisons,
  min_total = 10,
  min_group = 5
)
```

## Output Files

The pipeline generates:
  - `count_matrix.tsv`: Methylation counts per window and group
- `percent_matrix.tsv`: Methylation percentages per window and group  
- `dmr_filtered.tsv`: Significant windows
- `dmr_collapsed.tsv`: Merged DMR regions
- `dmrs_hyper.bed`, `dmrs_hypo.bed`: Genome browser tracks
- Various R objects saved as .rds files

## Statistical Considerations

### Method Assumptions
- Fisher's exact test assumes independence of observations within windows
- Multiple testing correction is applied genome-wide
- Spatial smoothing reduces resolution but may improve statistical power

### Parameter Selection
- Window size and smoothing parameters affect sensitivity and resolution
- Minimum count thresholds should be chosen based on data characteristics
- Statistical significance thresholds should account for multiple testing burden

## Implementation Status

This package is under development. Key areas requiring validation:
- Statistical test implementation accuracy
- Memory usage scaling with dataset size
- Computational performance characteristics
- Comparison with established DMR detection methods

## License

GNU General Public License v3.0 - see [LICENSE](LICENSE) file.

---

*This package implements one approach to DMR analysis. Users should validate results using appropriate controls and consider alternative methods for production analyses.*
