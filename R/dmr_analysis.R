# dmr_analysis.R - Simplified DMR Analysis from Single-Cell Methylation Data
# Focused reimplementation of Amethyst DMR functionality

library(data.table)
library(rhdf5)
library(GenomicRanges)
library(IRanges)

# Set data.table threads explicitly to avoid conflicts
setDTthreads(1)  # Will be configured by user

#' Index chromosomes in H5 files
#' 
#' @param h5_paths data.table with columns: barcode, path
#' @param type Character, methylation context (e.g., "CG", "CH")
#' @param chr_list Optional character vector of chromosomes to index
#' @param threads Integer, number of threads for parallel processing
#' @return List of data.tables by chromosome with h5 coordinates
index_chromosomes <- function(h5_paths, type = "CG", chr_list = NULL, threads = 1) {
  
  # Configure threading
  old_threads <- getDTthreads()
  setDTthreads(threads)
  on.exit(setDTthreads(old_threads))
  
  # Process each file
  all_indices <- rbindlist(lapply(seq_len(nrow(h5_paths)), function(i) {
    path <- h5_paths$path[i]
    barcode <- h5_paths$barcode[i]
    
    tryCatch({
      # Read chromosome data
      h5_data <- data.table(h5read(path, name = paste0(type, "/", barcode, "/1")))
      h5_data[, index := .I]
      
      # Filter chromosomes if list provided
      if (!is.null(chr_list)) {
        h5_data <- h5_data[chr %in% chr_list]
      } else {
        # Remove alternative contigs by default
        h5_data <- h5_data[!grepl("_|EBV|M", chr)]
      }
      
      # Calculate start and count for each chromosome
      chr_index <- h5_data[, .(
        cell_id = barcode,
        start = min(index),
        count = .N
      ), by = chr]
      
      return(chr_index)
    }, error = function(e) {
      warning(sprintf("Error indexing barcode %s: %s", barcode, e$message))
      return(NULL)
    })
  }))
  
  # Split by chromosome
  split(all_indices, all_indices$chr)
}

#' Calculate smoothed methylation windows
#' 
#' @param h5_paths data.table with columns: barcode, path
#' @param chr_index List of chromosome indices from index_chromosomes
#' @param metadata data.table with cell metadata, must have grouping column
#' @param type Character, methylation context
#' @param step Integer, window step size in bp
#' @param smooth Integer, number of windows to smooth over
#' @param group_by Character, column name in metadata for grouping
#' @param genome Character, genome build for chromosome sizes
#' @param threads Integer, number of threads
#' @return List with sum_matrix and pct_matrix
calc_smoothed_windows <- function(h5_paths, 
                                  chr_index,
                                  metadata,
                                  type = "CG",
                                  step = 500,
                                  smooth = 3,
                                  group_by = "cluster_id",
                                  genome = "hg38",
                                  threads = 1) {
  
  # Configure threading
  old_threads <- getDTthreads()
  setDTthreads(threads)
  on.exit(setDTthreads(old_threads))
  
  # Get chromosome sizes
  chr_sizes <- get_chromosome_sizes(genome)
  
  # Generate genomic windows
  windows <- generate_genomic_windows(chr_sizes, step)
  
  # Get groups from metadata
  groups <- unique(metadata[[group_by]])
  groups <- groups[!is.na(groups)]
  
  # Process by chromosome
  results_by_chr <- lapply(names(chr_index), function(chr) {
    
    # Get indices for this chromosome
    chr_indices <- chr_index[[chr]]
    
    # Initialize group results
    group_results <- lapply(groups, function(grp) {
      
      # Get member cells
      member_cells <- metadata[get(group_by) == grp, cell_id]
      member_indices <- chr_indices[cell_id %in% member_cells]
      
      if (nrow(member_indices) == 0) return(NULL)
      
      # Read and aggregate data for all members
      member_data <- rbindlist(lapply(seq_len(nrow(member_indices)), function(i) {
        idx <- member_indices[i]
        path <- h5_paths[barcode == idx$cell_id, path]
        
        tryCatch({
          # Read h5 data for this chromosome
          h5_data <- data.table(h5read(path, 
                                       name = paste0(type, "/", idx$cell_id, "/1"),
                                       start = idx$start,
                                       count = idx$count))
          
          # Assign to windows
          h5_data[, window := paste0(chr, "_", 
                                     (pos %/% step) * step, "_",
                                     ((pos %/% step) + 1) * step)]
          
          # Aggregate by window
          h5_data[, .(c = sum(c, na.rm = TRUE), 
                      t = sum(t, na.rm = TRUE)), 
                  by = window]
          
        }, error = function(e) {
          warning(sprintf("Error reading %s: %s", idx$cell_id, e$message))
          return(NULL)
        })
      }))
      
      if (is.null(member_data) || nrow(member_data) == 0) return(NULL)
      
      # Aggregate across all members
      result <- member_data[, .(c = sum(c), t = sum(t)), by = window]
      setnames(result, c("window", paste0(grp, "_c"), paste0(grp, "_t")))
      
      return(result)
    })
    
    # Merge group results
    Reduce(function(x, y) {
      if (is.null(x)) return(y)
      if (is.null(y)) return(x)
      merge(x, y, by = "window", all = TRUE)
    }, group_results)
  })
  
  # Combine chromosomes
  count_matrix <- rbindlist(results_by_chr, fill = TRUE)
  
  # Apply smoothing
  if (smooth > 1) {
    count_matrix <- apply_smoothing(count_matrix, smooth, step)
  }
  
  # Calculate percentage matrix
  pct_matrix <- calculate_pct_matrix(count_matrix, groups)
  
  list(sum_matrix = count_matrix, pct_matrix = pct_matrix)
}

#' Test for differentially methylated regions
#' 
#' @param sum_matrix data.table with counts from calc_smoothed_windows
#' @param comparisons data.frame with test specifications or NULL for all-vs-all
#' @param min_total Minimum total observations across groups
#' @param min_group Minimum observations per group
#' @return data.table with test results
test_dmr <- function(sum_matrix, 
                     comparisons = NULL, 
                     min_total = 3, 
                     min_group = 3) {
  
  # Copy to avoid modifying input
  counts <- copy(sum_matrix)
  
  # Filter by minimum total observations
  count_cols <- grep("_c$|_t$", names(counts), value = TRUE)
  total_counts <- rowSums(counts[, ..count_cols], na.rm = TRUE)
  counts <- counts[total_counts >= min_total]
  
  # Fast Fisher's exact test implementation
  fast_fisher <- function(ctg_table) {
    m <- ctg_table[1, 1] + ctg_table[2, 1]
    n <- ctg_table[1, 2] + ctg_table[2, 2]
    k <- ctg_table[1, 1] + ctg_table[1, 2]
    q <- ctg_table[1, 1]
    
    # Two-tailed test
    p_right <- phyper(q, m, n, k, lower.tail = FALSE) + 0.5 * dhyper(q, m, n, k)
    p_left <- phyper(q - 1, m, n, k, lower.tail = TRUE) + 0.5 * dhyper(q, m, n, k)
    
    2 * min(p_right, p_left, 0.5)
  }
  
  if (is.null(comparisons)) {
    # All vs all comparisons
    groups <- unique(sub("_c$", "", grep("_c$", names(counts), value = TRUE)))
    
    for (grp in groups) {
      # Get member and non-member columns
      member_c <- paste0(grp, "_c")
      member_t <- paste0(grp, "_t")
      
      other_c <- setdiff(grep("_c$", names(counts), value = TRUE), member_c)
      other_t <- setdiff(grep("_t$", names(counts), value = TRUE), member_t)
      
      # Calculate totals
      counts[, `:=`(
        member_c = get(member_c),
        member_t = get(member_t),
        nonmember_c = rowSums(.SD[, ..other_c], na.rm = TRUE),
        nonmember_t = rowSums(.SD[, ..other_t], na.rm = TRUE)
      )]
      
      # Apply minimum group filter
      counts[member_c + member_t < min_group | 
               nonmember_c + nonmember_t < min_group, 
             c("member_c", "member_t", "nonmember_c", "nonmember_t") := NA]
      
      # Apply Fisher's test
      counts[, paste0(grp, "_pval") := apply(.SD, 1, function(x) {
        if (any(is.na(x))) return(NA_real_)
        fast_fisher(matrix(x, nrow = 2, byrow = TRUE))
      }), .SDcols = c("member_c", "member_t", "nonmember_c", "nonmember_t")]
      
      # Calculate log fold change
      counts[, paste0(grp, "_logFC") := log2(
        (member_c / (member_c + member_t)) / 
          (nonmember_c / (nonmember_c + nonmember_t))
      )]
    }
    
    # Clean up temporary columns
    counts[, c("member_c", "member_t", "nonmember_c", "nonmember_t") := NULL]
    
  } else {
    # Custom comparisons
    for (i in seq_len(nrow(comparisons))) {
      name <- comparisons$name[i]
      grp_a <- strsplit(comparisons$A[i], ",")[[1]]
      grp_b <- strsplit(comparisons$B[i], ",")[[1]]
      
      # Calculate totals for each comparison
      a_c_cols <- paste0(grp_a, "_c")
      a_t_cols <- paste0(grp_a, "_t")
      b_c_cols <- paste0(grp_b, "_c")
      b_t_cols <- paste0(grp_b, "_t")
      
      counts[, `:=`(
        member_c = rowSums(.SD[, ..a_c_cols], na.rm = TRUE),
        member_t = rowSums(.SD[, ..a_t_cols], na.rm = TRUE),
        nonmember_c = rowSums(.SD[, ..b_c_cols], na.rm = TRUE),
        nonmember_t = rowSums(.SD[, ..b_t_cols], na.rm = TRUE)
      )]
      
      # Apply tests as above
      counts[member_c + member_t < min_group | 
               nonmember_c + nonmember_t < min_group,
             c("member_c", "member_t", "nonmember_c", "nonmember_t") := NA]
      
      counts[, paste0(name, "_pval") := apply(.SD, 1, function(x) {
        if (any(is.na(x))) return(NA_real_)
        fast_fisher(matrix(x, nrow = 2, byrow = TRUE))
      }), .SDcols = c("member_c", "member_t", "nonmember_c", "nonmember_t")]
      
      counts[, paste0(name, "_logFC") := log2(
        (member_c / (member_c + member_t)) / 
          (nonmember_c / (nonmember_c + nonmember_t))
      )]
    }
    
    counts[, c("member_c", "member_t", "nonmember_c", "nonmember_t") := NULL]
  }
  
  return(counts)
}

#' Filter and reshape DMR results
#' 
#' @param dmr_matrix Results from test_dmr
#' @param method Multiple testing correction method
#' @param p_threshold Maximum adjusted p-value
#' @param log_threshold Minimum absolute log fold change
#' @param filter Logical, whether to filter results
#' @return data.table with filtered results
filter_dmr <- function(dmr_matrix, 
                       method = "bonferroni",
                       p_threshold = 0.01,
                       log_threshold = 1,
                       filter = TRUE) {
  
  results <- copy(dmr_matrix)
  
  # Get test columns
  pval_cols <- grep("_pval$", names(results), value = TRUE)
  logfc_cols <- grep("_logFC$", names(results), value = TRUE)
  
  # Melt to long format
  id_cols <- setdiff(names(results), c(pval_cols, logfc_cols))
  
  results <- melt(results, 
                  id.vars = id_cols,
                  measure.vars = list(pval = pval_cols, logFC = logfc_cols),
                  variable.name = "test")
  
  # Extract test names
  results[, test := gsub("_pval$", "", pval_cols[test])]
  
  # Remove NA values
  results <- results[!is.na(pval)]
  
  # Apply multiple testing correction
  results[, padj := p.adjust(pval, method = method), by = test]
  
  # Add direction
  results[, direction := ifelse(logFC < 0, "hypo", "hyper")]
  
  # Filter if requested
  if (filter) {
    results <- results[padj < p_threshold & abs(logFC) > log_threshold]
  }
  
  return(results)
}

#' Collapse adjacent DMRs
#' 
#' @param dmr_filtered Filtered results from filter_dmr
#' @param max_gap Maximum gap between DMRs to merge
#' @param min_length Minimum length of collapsed DMR
#' @return data.table with collapsed DMRs
collapse_dmr <- function(dmr_filtered, 
                         max_gap = 0,
                         min_length = 1000) {
  
  results <- copy(dmr_filtered)
  
  # Ensure we have the expected columns
  if (!all(c("chr", "start", "end") %in% names(results))) {
    stop("DMR results must have chr, start, end columns")
  }
  
  # Group by test and direction
  setorder(results, test, direction, chr, start)
  
  # Use GenomicRanges for merging
  collapsed <- results[, {
    # Create ranges
    gr <- GRanges(seqnames = chr,
                  ranges = IRanges(start = start - max_gap, 
                                   end = end + max_gap))
    
    # Reduce overlapping ranges
    gr_reduced <- reduce(gr)
    
    # Convert back to data.table
    dt_reduced <- data.table(
      chr = as.character(seqnames(gr_reduced)),
      dmr_start = start(gr_reduced) + max_gap,
      dmr_end = end(gr_reduced) - max_gap
    )
    
    # Calculate length
    dt_reduced[, dmr_length := dmr_end - dmr_start + 1]
    
    # Map back statistics (mean of merged regions)
    dt_reduced[, `:=`(
      dmr_padj = mean(padj),
      dmr_logFC = mean(logFC),
      n_windows = .N
    )]
    
    dt_reduced
  }, by = .(test, direction)]
  
  # Filter by minimum length
  collapsed <- collapsed[dmr_length >= min_length]
  
  return(collapsed)
}

# Helper functions

#' Get chromosome sizes for a genome
get_chromosome_sizes <- function(genome) {
  sizes <- switch(genome,
                  "hg38" = data.table(
                    chr = paste0("chr", c(1:22, "X", "Y")),
                    size = c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 
                             159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
                             114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
                             58617616, 64444167, 46709983, 50818468, 156040895, 57227415)
                  ),
                  "hg19" = data.table(
                    chr = paste0("chr", c(1:22, "X", "Y")),
                    size = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067,
                             159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
                             115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
                             59128983, 63025520, 48129895, 51304566, 155270560, 59373566)
                  ),
                  "mm10" = data.table(
                    chr = paste0("chr", c(1:19, "X", "Y")),
                    size = c(195471971, 182113224, 160039680, 156508116, 151834684, 149736546,
                             145441459, 129401213, 124595110, 130694993, 122082543, 120129022,
                             120421639, 124902244, 104043685, 98207768, 94987271, 90702639,
                             61431566, 171031299, 91744698)
                  ),
                  stop("Unsupported genome. Use 'hg38', 'hg19', or 'mm10'")
  )
  return(sizes)
}

#' Generate genomic windows
generate_genomic_windows <- function(chr_sizes, step) {
  rbindlist(lapply(seq_len(nrow(chr_sizes)), function(i) {
    chr <- chr_sizes$chr[i]
    size <- chr_sizes$size[i]
    starts <- seq(0, size - 1, by = step)
    data.table(
      chr = chr,
      start = starts,
      end = pmin(starts + step, size),
      window = paste0(chr, "_", starts, "_", pmin(starts + step, size))
    )
  }))
}

#' Apply smoothing to count matrix
apply_smoothing <- function(count_matrix, smooth, step) {
  # Extract location info
  count_matrix[, c("chr", "start", "end") := tstrsplit(window, "_", fixed = TRUE, type.convert = TRUE)]
  setorder(count_matrix, chr, start)
  
  # Get count columns
  count_cols <- grep("_c$|_t$", names(count_matrix), value = TRUE)
  
  # Apply rolling sum by chromosome
  count_matrix[, (count_cols) := lapply(.SD, function(x) {
    frollsum(x, n = smooth, align = "center", fill = NA)
  }), .SDcols = count_cols, by = chr]
  
  # Remove rows with NA (edges)
  count_matrix <- count_matrix[complete.cases(count_matrix[, ..count_cols])]
  
  # Update positions to reflect smoothed windows
  count_matrix[, `:=`(
    start = start + (smooth - 1) * step / 2,
    end = end - (smooth - 1) * step / 2,
    window = paste0(chr, "_", start, "_", end)
  )]
  
  count_matrix[, c("chr", "start", "end") := NULL]
  
  return(count_matrix)
}

#' Calculate percentage methylation matrix
calculate_pct_matrix <- function(count_matrix, groups) {
  pct_matrix <- copy(count_matrix)
  
  for (grp in groups) {
    c_col <- paste0(grp, "_c")
    t_col <- paste0(grp, "_t")
    
    pct_matrix[, (grp) := round(100 * get(c_col) / (get(c_col) + get(t_col)), 2)]
    pct_matrix[, c(c_col, t_col) := NULL]
  }
  
  return(pct_matrix)
}

# Example usage:
if (FALSE) {
  # Setup
  h5_paths <- data.table(
    barcode = c("cell1", "cell2", "cell3"),
    path = c("path1.h5", "path2.h5", "path3.h5")
  )
  
  metadata <- data.table(
    cell_id = c("cell1", "cell2", "cell3"),
    cluster_id = c("A", "A", "B")
  )
  
  # Run analysis
  chr_index <- index_chromosomes(h5_paths, type = "CG", threads = 4)
  
  windows <- calc_smoothed_windows(
    h5_paths, chr_index, metadata, 
    type = "CG", step = 500, smooth = 3,
    group_by = "cluster_id", threads = 4
  )
  
  dmr_results <- test_dmr(
    windows$sum_matrix,
    min_total = 10,
    min_group = 5
  )
  
  dmr_filtered <- filter_dmr(
    dmr_results,
    method = "bonferroni",
    p_threshold = 0.01,
    log_threshold = 1.5
  )
  
  dmr_collapsed <- collapse_dmr(
    dmr_filtered,
    max_gap = 2000,
    min_length = 2000
  )
}