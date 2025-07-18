#!/usr/bin/env Rscript
# generate_toy_data.R - Generate toy methylation data for DMR analysis example

library(rhdf5)
library(data.table)

#' Generate methylation data for a single cell
#' 
#' @param n_sites Number of methylation sites to generate
#' @param meth_prob Base probability of methylation
#' @param coverage_mean Mean coverage per site
#' @param chr_sizes Named vector of chromosome sizes
#' @param dmr_regions data.table with DMR regions and effect sizes
generate_cell_methylation <- function(n_sites = 100000,
                                      meth_prob = 0.7,
                                      coverage_mean = 10,
                                      chr_sizes,
                                      dmr_regions = NULL) {
  
  # Distribute sites across chromosomes proportionally
  chr_props <- chr_sizes / sum(chr_sizes)
  sites_per_chr <- round(n_sites * chr_props)
  
  # Generate data for each chromosome
  all_data <- rbindlist(lapply(names(chr_sizes), function(chr) {
    n_chr_sites <- sites_per_chr[chr]
    
    # Generate random positions
    positions <- sort(sample(1:chr_sizes[chr], n_chr_sites, replace = FALSE))
    
    # Generate coverage (Poisson distributed)
    coverage <- rpois(n_chr_sites, lambda = coverage_mean)
    coverage[coverage == 0] <- 1  # Ensure at least 1 read
    
    # Base methylation probability
    site_meth_prob <- rep(meth_prob, n_chr_sites)
    
    # Apply DMR effects if provided
    if (!is.null(dmr_regions)) {
      chr_dmrs <- dmr_regions[chr == chr]
      for (i in seq_len(nrow(chr_dmrs))) {
        dmr <- chr_dmrs[i]
        # Find sites in DMR
        in_dmr <- positions >= dmr$start & positions <= dmr$end
        # Apply effect
        site_meth_prob[in_dmr] <- plogis(qlogis(site_meth_prob[in_dmr]) + dmr$effect)
      }
    }
    
    # Generate methylation calls
    methylated <- rbinom(n_chr_sites, size = coverage, prob = site_meth_prob)
    unmethylated <- coverage - methylated
    
    data.table(
      chr = chr,
      pos = positions,
      c = methylated,
      t = unmethylated
    )
  }))
  
  return(all_data)
}

#' Create H5 file with methylation data
#' 
#' @param filename Output H5 file name
#' @param cell_data List of data.tables with methylation data per cell
#' @param barcodes Character vector of cell barcodes
#' @param contexts Character vector of methylation contexts
create_h5_file <- function(filename, cell_data, barcodes, contexts = c("CG", "CH")) {
  
  # Create H5 file
  h5createFile(filename)
  
  # Add data for each context
  for (context in contexts) {
    # Create context group
    h5createGroup(filename, context)
    
    # Add data for each cell
    for (i in seq_along(barcodes)) {
      barcode <- barcodes[i]
      cell_path <- paste0(context, "/", barcode)
      
      # Create cell group
      h5createGroup(filename, cell_path)
      
      # Write data
      data_path <- paste0(cell_path, "/1")
      
      if (context == "CG") {
        # For CG, use the provided data
        h5write(cell_data[[i]], filename, data_path)
      } else {
        # For CH, generate different data with lower methylation
        ch_data <- copy(cell_data[[i]])
        # CH methylation is typically much lower
        ch_meth_prob <- 0.02
        ch_coverage <- ch_data$c + ch_data$t
        ch_methylated <- rbinom(nrow(ch_data), size = ch_coverage, prob = ch_meth_prob)
        ch_data$c <- ch_methylated
        ch_data$t <- ch_coverage - ch_methylated
        
        h5write(ch_data, filename, data_path)
      }
    }
  }
  
  # Close file
  h5closeAll()
}

# Set random seed for reproducibility
set.seed(42)

# Define parameters
barcodes <- c("ACGCGACGGCTATACCGAAGCGCCTATA", 
              "ACTTCTGCCATCGATGATACGTATGGCA",
              "ACTTCTGCCAGTACCTGAATCTGATGCT",
              "TTCATATCAACTACTGCAAGATCTGAAT")

cell_types <- c("T_cell", "T_cell", "B_cell", "B_cell")

# Chromosome sizes (subset of hg38 for toy data)
chr_sizes <- c(
  chr1 = 248956422,
  chr2 = 242193529,
  chr3 = 198295559,
  chr4 = 190214555,
  chr5 = 181538259
)

# Define some DMR regions with effects
# Positive effect = hypermethylated in T cells
# Negative effect = hypomethylated in T cells
dmr_regions <- data.table(
  chr = c("chr1", "chr1", "chr2", "chr2", "chr3", "chr3"),
  start = c(1000000, 5000000, 2000000, 8000000, 3000000, 7000000),
  end = c(1010000, 5010000, 2010000, 8010000, 3010000, 7010000),
  effect = c(2, -2, 1.5, -1.5, 2.5, -2.5)  # log odds ratio
)

# Generate cell-specific DMR effects
t_cell_dmrs <- dmr_regions
b_cell_dmrs <- copy(dmr_regions)
b_cell_dmrs$effect <- -b_cell_dmrs$effect  # Reverse effects for B cells

# Generate methylation data for each cell
cell_data <- list()

for (i in seq_along(barcodes)) {
  cat(sprintf("Generating data for cell %d (%s)...\n", i, cell_types[i]))
  
  # Use appropriate DMR effects based on cell type
  if (cell_types[i] == "T_cell") {
    cell_dmrs <- t_cell_dmrs
    base_meth <- 0.75  # Slightly higher base methylation
  } else {
    cell_dmrs <- b_cell_dmrs
    base_meth <- 0.70  # Slightly lower base methylation
  }
  
  # Add some cell-to-cell variation
  cell_meth <- base_meth + rnorm(1, 0, 0.02)
  
  # Generate data
  cell_data[[i]] <- generate_cell_methylation(
    n_sites = 50000,  # Fewer sites for toy data
    meth_prob = cell_meth,
    coverage_mean = 15,
    chr_sizes = chr_sizes,
    dmr_regions = cell_dmrs
  )
}

# Create H5 file
output_file <- "./Data/methylation_data.h5"
cat(sprintf("Writing H5 file: %s\n", output_file))

create_h5_file(
  filename = output_file,
  cell_data = cell_data,
  barcodes = barcodes,
  contexts = c("CG", "CH")
)

# Generate metadata file
metadata <- data.table(
  cell_id = barcodes,
  cluster_id = cell_types,
  sample = c("sample1", "sample1", "sample2", "sample2"),
  coverage = sapply(cell_data, function(x) sum(x$c + x$t)),
  n_sites = sapply(cell_data, nrow),
  mean_meth_cg = sapply(cell_data, function(x) sum(x$c) / sum(x$c + x$t))
)

fwrite(metadata, "./Data/cell_metadata.tsv", sep = "\t")

# Generate summary statistics
cat("\n=== Toy Data Generation Summary ===\n")
cat(sprintf("Created H5 file: %s\n", output_file))
cat(sprintf("Number of cells: %d\n", length(barcodes)))
cat(sprintf("Chromosomes included: %s\n", paste(names(chr_sizes), collapse = ", ")))
cat(sprintf("Sites per cell: ~%d\n", mean(metadata$n_sites)))
cat(sprintf("Mean coverage: ~%.1fx\n", mean(metadata$coverage) / mean(metadata$n_sites)))
cat(sprintf("Number of DMRs: %d\n", nrow(dmr_regions)))

# Show methylation levels by cell type
cat("\nMean methylation by cell type:\n")
print(metadata[, .(mean_meth = mean(mean_meth_cg)), by = cluster_id])

# Verify file structure
cat("\n=== H5 File Structure ===\n")
h5ls(output_file)

# Test reading data back
cat("\n=== Testing data retrieval ===\n")
test_barcode <- barcodes[1]
test_data <- h5read(output_file, paste0("CG/", test_barcode, "/1"))
cat(sprintf("Successfully read %d CG sites for barcode %s\n", 
            nrow(test_data), test_barcode))
cat(sprintf("Chromosomes: %s\n", paste(unique(test_data$chr), collapse = ", ")))
cat(sprintf("Mean methylation: %.3f\n", sum(test_data$c) / sum(test_data$c + test_data$t)))

# Create a small example DMR region file for testing
# These correspond to the simulated DMRs
known_dmrs <- data.table(
  chr = dmr_regions$chr,
  start = dmr_regions$start,
  end = dmr_regions$end,
  name = paste0("DMR_", 1:nrow(dmr_regions)),
  score = abs(dmr_regions$effect) * 100,
  strand = ifelse(dmr_regions$effect > 0, "+", "-")
)

fwrite(known_dmrs, "./Data/known_dmrs.bed", sep = "\t", col.names = FALSE)
cat("\nCreated known DMRs file: known_dmrs.bed\n")

# Create example H5 paths file
h5_paths <- data.table(
  barcode = barcodes,
  path = rep(output_file, length(barcodes))
)
fwrite(h5_paths, "./Data/h5_paths.tsv", sep = "\t")

cat("\n=== Files created ===\n")
cat("1. methylation_data.h5 - H5 file with methylation data\n")
cat("2. cell_metadata.tsv - Cell metadata\n")
cat("3. known_dmrs.bed - Known DMR regions\n")
cat("4. h5_paths.tsv - H5 file paths\n")

cat("\nYou can now run the DMR analysis example:\n")
cat("  source('example_dmr_analysis.R')\n")
