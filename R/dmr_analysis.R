# dmr_analysis.R - Simplified DMR Analysis from Single-Cell Methylation Data
# Focused reimplementation of Amethyst DMR functionality

#' @import data.table
#' @import GenomicRanges
#' @importFrom rhdf5 h5read h5ls H5Fopen H5Dopen H5Dget_space H5Sget_simple_extent_dims H5Dget_type H5Tget_class h5closeAll
#' @importFrom IRanges IRanges
#' @importFrom stats phyper dhyper p.adjust
NULL

# Set data.table threads explicitly to avoid conflicts
data.table::setDTthreads(1)  # Will be configured by user

#' Index chromosomes in H5 files
#' 
#' @param h5_paths data.table with columns: barcode, path
#' @param type Character, methylation context (e.g., "CG", "CH")
#' @param chr_list Optional character vector of chromosomes to index
#' @param threads Integer, number of threads for parallel processing
#' @return List of data.tables by chromosome with h5 coordinates
#' @export
index_chromosomes <- function(h5_paths, type = "CG", chr_list = NULL, threads = 1) {
  
  # Configure threading
  old_threads <- data.table::getDTthreads()
  data.table::setDTthreads(threads)
  on.exit(data.table::setDTthreads(old_threads))
  
  # Process each file
  all_indices <- data.table::rbindlist(lapply(seq_len(nrow(h5_paths)), function(i) {
    path <- h5_paths$path[i]
    barcode <- h5_paths$barcode[i]
    
    tryCatch({
      # Build the data path
      data_path <- paste0(type, "/", barcode, "/1")
      
      # Read the full dataset (1D compound array)
      h5_data <- rhdf5::h5read(path, name = data_path)
      
      # Convert to data.table if needed
      if (!data.table::is.data.table(h5_data)) {
        h5_data <- data.table::as.data.table(h5_data)
      }
      
      # Get unique chromosomes and their positions
      chr_summary <- h5_data[, .(
        cell_id = barcode,
        first_pos = min(.I),
        last_pos = max(.I),
        n_sites = .N
      ), by = chr]
      
      # Filter chromosomes if list provided
      if (!is.null(chr_list)) {
        chr_summary <- chr_summary[chr %in% chr_list]
      } else {
        # Remove alternative contigs by default
        chr_summary <- chr_summary[!grepl("_|EBV|M", chr)]
      }
      
      return(chr_summary)
      
    }, error = function(e) {
      warning(sprintf("Error indexing barcode %s: %s", barcode, e$message))
      return(NULL)
    })
  }))
  
  # Split by chromosome
  if (nrow(all_indices) > 0) {
    chr_index <- split(all_indices, all_indices$chr)
  } else {
    chr_index <- list()
  }
  
  return(chr_index)
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
#' @export
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
  old_threads <- data.table::getDTthreads()
  data.table::setDTthreads(threads)
  on.exit(data.table::setDTthreads(old_threads))
  
  # Get chromosome sizes
  chr_sizes <- get_chromosome_sizes(genome)
  
  # Generate genomic windows
  windows <- generate_genomic_windows(chr_sizes, step)
  
  # Get groups from metadata
  groups <- unique(metadata[[group_by]])
  groups <- groups[!is.na(groups)]
  
  # Process by chromosome
  results_by_chr <- lapply(names(chr_index), function(chr_name) {
    
    # Get indices for this chromosome
    chr_indices <- chr_index[[chr_name]]
    
    # Initialize group results
    group_results <- lapply(groups, function(grp) {
      
      # Get member cells
      member_cells <- metadata[get(group_by) == grp, cell_id]
      member_indices <- chr_indices[cell_id %in% member_cells]
      
      if (nrow(member_indices) == 0) return(NULL)
      
      # Read and aggregate data for all members
      member_data <- data.table::rbindlist(lapply(seq_len(nrow(member_indices)), function(i) {
        idx <- member_indices[i]
        path <- h5_paths[barcode == idx$cell_id, path]
        
        tryCatch({
          # Read h5 data
          data_path <- paste0(type, "/", idx$cell_id, "/1")
          h5_data <- rhdf5::h5read(path, name = data_path)
          
          # Convert to data.table
          if (!data.table::is.data.table(h5_data)) {
            h5_data <- data.table::as.data.table(h5_data)
          }
          
          # Filter to this chromosome
          h5_data <- h5_data[chr == chr_name]
          
          if (nrow(h5_data) == 0) return(NULL)
          
          # Assign to windows
          h5_data[, window := paste0(chr, "_", 
                                     (pos %/% step) * step, "_",
                                     ((pos %/% step) + 1) * step)]
          
          # Aggregate by window
          window_summary <- h5_data[, .(
            c = sum(c, na.rm = TRUE), 
            t = sum(t, na.rm = TRUE)
          ), by = window]
          
          return(window_summary)
          
        }, error = function(e) {
          warning(sprintf("Error reading %s: %s", idx$cell_id, e$message))
          return(NULL)
        })
      }))
      
      if (is.null(member_data) || nrow(member_data) == 0) return(NULL)
      
      # Aggregate across all members
      result <- member_data[, .(c = sum(c), t = sum(t)), by = window]
      data.table::setnames(result, c("window", paste0(grp, "_c"), paste0(grp, "_t")))
      
      return(result)
    })
    
    # Merge group results
    merged_result <- Reduce(function(x, y) {
      if (is.null(x)) return(y)
      if (is.null(y)) return(x)
      merge(x, y, by = "window", all = TRUE)
    }, group_results)
    
    return(merged_result)
  })
  
  # Combine chromosomes
  count_matrix <- data.table::rbindlist(results_by_chr, fill = TRUE)
  
  # Check if we got any data
  if (nrow(count_matrix) == 0) {
    stop("No data found. Check that H5 files contain data for the specified chromosomes.")
  }
  
  # Replace NA values with 0 in count columns
  count_cols <- grep("_c$|_t$", names(count_matrix), value = TRUE)
  for (col in count_cols) {
    count_matrix[is.na(get(col)), (col) := 0]
  }
  
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
#' @export
test_dmr <- function(sum_matrix, 
                     comparisons = NULL, 
                     min_total = 3, 
                     min_group = 3) {
  
  # Copy to avoid modifying input
  counts <- data.table::copy(sum_matrix)
  
  # Filter windows with no data
  counts <- counts[!is.na(window)]
  
  # Filter by minimum total observations
  count_cols <- grep("_c$|_t$", names(counts), value = TRUE)
  if (length(count_cols) == 0) {
    stop("No count columns found. Expected columns ending in _c or _t")
  }
  
  # Replace NA with 0 in count columns
  for (col in count_cols) {
    counts[is.na(get(col)), (col) := 0]
  }
  
  total_counts <- rowSums(counts[, ..count_cols], na.rm = TRUE)
  counts <- counts[total_counts >= min_total]
  
  # Fast Fisher's exact test implementation
  fast_fisher <- function(ctg_table) {
    m <- ctg_table[1, 1] + ctg_table[2, 1]
    n <- ctg_table[1, 2] + ctg_table[2, 2]
    k <- ctg_table[1, 1] + ctg_table[1, 2]
    q <- ctg_table[1, 1]
    
    # Two-tailed test
    p_right <- stats::phyper(q, m, n, k, lower.tail = FALSE) + 0.5 * stats::dhyper(q, m, n, k)
    p_left <- stats::phyper(q - 1, m, n, k, lower.tail = TRUE) + 0.5 * stats::dhyper(q, m, n, k)
    
    2 * min(p_right, p_left, 0.5)
  }
  
  if (is.null(comparisons)) {
    # All vs all comparisons
    groups <- unique(sub("_c$", "", grep("_c$", names(counts), value = TRUE)))
    
    for (grp in groups) {
      # Get member and non-member columns
      member_c_col <- paste0(grp, "_c")
      member_t_col <- paste0(grp, "_t")
      
      other_c_cols <- setdiff(grep("_c$", names(counts), value = TRUE), member_c_col)
      other_t_cols <- setdiff(grep("_t$", names(counts), value = TRUE), member_t_col)
      
      # Calculate totals using separate operations
      counts[, member_c := get(member_c_col)]
      counts[, member_t := get(member_t_col)]
      
      if (length(other_c_cols) > 0) {
        counts[, nonmember_c := rowSums(.SD, na.rm = TRUE), .SDcols = other_c_cols]
      } else {
        counts[, nonmember_c := 0]
      }
      
      if (length(other_t_cols) > 0) {
        counts[, nonmember_t := rowSums(.SD, na.rm = TRUE), .SDcols = other_t_cols]
      } else {
        counts[, nonmember_t := 0]
      }
      
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
      
      # Validate columns exist
      if (!all(a_c_cols %in% names(counts))) {
        warning(sprintf("Missing columns for group A in comparison %s", name))
        next
      }
      if (!all(b_c_cols %in% names(counts))) {
        warning(sprintf("Missing columns for group B in comparison %s", name))
        next
      }
      
      # Calculate totals
      if (length(a_c_cols) > 0) {
        counts[, member_c := rowSums(.SD, na.rm = TRUE), .SDcols = a_c_cols]
      } else {
        counts[, member_c := 0]
      }
      
      if (length(a_t_cols) > 0) {
        counts[, member_t := rowSums(.SD, na.rm = TRUE), .SDcols = a_t_cols]
      } else {
        counts[, member_t := 0]
      }
      
      if (length(b_c_cols) > 0) {
        counts[, nonmember_c := rowSums(.SD, na.rm = TRUE), .SDcols = b_c_cols]
      } else {
        counts[, nonmember_c := 0]
      }
      
      if (length(b_t_cols) > 0) {
        counts[, nonmember_t := rowSums(.SD, na.rm = TRUE), .SDcols = b_t_cols]
      } else {
        counts[, nonmember_t := 0]
      }
      
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
#' @export
filter_dmr <- function(dmr_matrix, 
                       method = "bonferroni",
                       p_threshold = 0.01,
                       log_threshold = 1,
                       filter = TRUE) {
  
  results <- data.table::copy(dmr_matrix)
  
  # Parse window coordinates if present
  if ("window" %in% names(results) && !all(c("chr", "start", "end") %in% names(results))) {
    results[, c("chr", "start", "end") := data.table::tstrsplit(window, "_", fixed = TRUE, type.convert = TRUE)]
  }
  
  # Get test columns
  pval_cols <- grep("_pval$", names(results), value = TRUE)
  logfc_cols <- grep("_logFC$", names(results), value = TRUE)
  
  # Melt to long format
  id_cols <- setdiff(names(results), c(pval_cols, logfc_cols))
  
  results <- data.table::melt(results, 
                              id.vars = id_cols,
                              measure.vars = list(pval = pval_cols, logFC = logfc_cols),
                              variable.name = "test")
  
  # Extract test names
  results[, test := gsub("_pval$", "", pval_cols[test])]
  
  # Remove NA values
  results <- results[!is.na(pval)]
  
  # Apply multiple testing correction
  results[, padj := stats::p.adjust(pval, method = method), by = test]
  
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
#' @export
collapse_dmr <- function(dmr_filtered, 
                         max_gap = 0,
                         min_length = 1000) {
  
  results <- data.table::copy(dmr_filtered)
  
  # Ensure we have the expected columns
  if (!all(c("chr", "start", "end") %in% names(results))) {
    stop("DMR results must have chr, start, end columns")
  }
  
  # Group by test and direction
  data.table::setorder(results, test, direction, chr, start)
  
  # Use GenomicRanges for merging
  collapsed <- results[, {
    # Create ranges
    gr <- GenomicRanges::GRanges(seqnames = chr,
                                 ranges = IRanges::IRanges(start = start - max_gap, 
                                                           end = end + max_gap))
    
    # Reduce overlapping ranges
    gr_reduced <- GenomicRanges::reduce(gr)
    
    # Convert back to data.table
    dt_reduced <- data.table::data.table(
      chr = as.character(GenomicRanges::seqnames(gr_reduced)),
      dmr_start = GenomicRanges::start(gr_reduced) + max_gap,
      dmr_end = GenomicRanges::end(gr_reduced) - max_gap
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
#' @keywords internal
get_chromosome_sizes <- function(genome) {
  sizes <- switch(genome,
                  "hg38" = data.table::data.table(
                    chr = paste0("chr", c(1:22, "X", "Y")),
                    size = c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 
                             159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
                             114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
                             58617616, 64444167, 46709983, 50818468, 156040895, 57227415)
                  ),
                  "hg19" = data.table::data.table(
                    chr = paste0("chr", c(1:22, "X", "Y")),
                    size = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067,
                             159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
                             115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
                             59128983, 63025520, 48129895, 51304566, 155270560, 59373566)
                  ),
                  "mm10" = data.table::data.table(
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
#' @keywords internal
generate_genomic_windows <- function(chr_sizes, step) {
  data.table::rbindlist(lapply(seq_len(nrow(chr_sizes)), function(i) {
    chr <- chr_sizes$chr[i]
    size <- chr_sizes$size[i]
    starts <- seq(0, size - 1, by = step)
    data.table::data.table(
      chr = chr,
      start = starts,
      end = pmin(starts + step, size),
      window = paste0(chr, "_", starts, "_", pmin(starts + step, size))
    )
  }))
}

#' Apply smoothing to count matrix
#' @keywords internal
apply_smoothing <- function(count_matrix, smooth, step) {
  # Extract location info
  count_matrix[, c("chr", "start", "end") := data.table::tstrsplit(window, "_", fixed = TRUE, type.convert = TRUE)]
  data.table::setorder(count_matrix, chr, start)
  
  # Get count columns
  count_cols <- grep("_c$|_t$", names(count_matrix), value = TRUE)
  
  # Apply rolling sum by chromosome
  count_matrix[, (count_cols) := lapply(.SD, function(x) {
    as.integer(data.table::frollsum(x, n = smooth, align = "center", fill = NA))
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
#' @keywords internal
calculate_pct_matrix <- function(count_matrix, groups) {
  pct_matrix <- data.table::copy(count_matrix)
  
  for (grp in groups) {
    c_col <- paste0(grp, "_c")
    t_col <- paste0(grp, "_t")
    
    if (c_col %in% names(pct_matrix) && t_col %in% names(pct_matrix)) {
      # Calculate percentage
      total <- pct_matrix[[c_col]] + pct_matrix[[t_col]]
      pct_val <- ifelse(total > 0, 
                        round(100 * pct_matrix[[c_col]] / total, 2),
                        NA_real_)
      
      # Add percentage column
      pct_matrix[, (grp) := pct_val]
      
      # Remove count columns
      pct_matrix[, c(c_col, t_col) := NULL]
    }
  }
  
  return(pct_matrix)
}