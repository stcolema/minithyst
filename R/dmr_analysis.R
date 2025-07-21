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

#' Calculate smoothed methylation windows (Exact Amethyst Replication)
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
#' @param chrList Optional character vector of chromosomes to process
#' @param chrSizes Optional numeric vector of chromosome sizes (must match chrList)
#' @return List with sum_matrix and pct_matrix
calc_smoothed_windows <- function(h5_paths, 
                                  chr_index,
                                  metadata,
                                  type = "CG",
                                  step = 500,
                                  smooth = 3,
                                  group_by = "cluster_id",
                                  genome = "hg38",
                                  threads = 1,
                                  chrList = NULL,
                                  chrSizes = NULL) {
  
  # Configure threading (amethyst approach)
  old_threads <- getDTthreads()
  setDTthreads(1)  # Keep data.table single-threaded
  on.exit(setDTthreads(old_threads))
  
  # Validate inputs (following amethyst)
  if (!is.null(chrList) && is.null(chrSizes)) {
    stop("If a chromosome whitelist is used, you must provide matched chromosome sizes.")
  }
  
  # STEP 1: Generate complete genomic windows (exactly like amethyst)
  generate_windows <- function(chromosome, size) {
    starts <- seq(0, size - 1, by = step)
    ends <- pmin(starts + step, size)
    data.frame(chr = chromosome, start = starts, end = ends)
  }
  
  # Get genome definitions (exactly like amethyst)
  if (is.null(chrList)) {
    if (genome %in% c("hg19", "hg38")) {
      chromosome <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                      "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
                      "chr20", "chr21", "chr22", "chrX", "chrY")
      if (genome == "hg19") {
        size <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431,
                  135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
                  59128983, 63025520, 48129895, 51304566, 155270560, 59373566)
      } else if (genome == "hg38") {
        size <- c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636,
                  138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345,
                  83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 57227415)
      }
    } else if (genome %in% c("mm10", "mm39")) {
      chromosome <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                      "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
                      "chrX", "chrY")
      if (genome == "mm10") {
        size <- c(195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 145441459, 129401213,
                  124595110, 130694993, 122082543, 120129022, 120421639, 124902244, 104043685, 98207768,
                  94987271, 90702639, 61431566, 171031299, 91744698)
      } else if (genome == "mm39") {
        size <- c(195154279, 181755017, 159745316, 156860686, 151758149, 149588044, 144995196, 130127694,
                  124359700, 130530862, 121973369, 120092757, 120883175, 125139656, 104073951, 98008968,
                  95294699, 90720763, 61420004, 169476592, 91455967)
      }
    } else {
      stop("Only hg19, hg38, mm10, or mm39 can be accommodated at this time.")
    }
    chromosome_sizes <- data.frame(chromosome, size)
  } else {
    chromosome_sizes <- data.frame(chromosome = chrList, size = chrSizes)
  }
  
  # Filter to available chromosomes
  available_chrs <- intersect(chromosome_sizes$chromosome, names(chr_index))
  chromosome_sizes <- chromosome_sizes[chromosome_sizes$chromosome %in% available_chrs, ]
  
  # Generate complete genomic window grid (exactly like amethyst)
  genomechunks <- do.call(rbind, lapply(1:nrow(chromosome_sizes), function(i) {
    generate_windows(chromosome_sizes$chromosome[i], chromosome_sizes$size[i])
  }))
  
  setDT(genomechunks)
  genomechunks[, window := paste0(chr, "_", start, "_", end)]
  setkey(genomechunks, chr, start, end)
  
  cat(sprintf("Generated %d genomic windows\n", nrow(genomechunks)))
  
  # STEP 2: Get groups and chromosome lists (exactly like amethyst)
  # Create membership data.frame with proper rownames like amethyst expects
  membership <- data.frame(membership = metadata[[group_by]])
  rownames(membership) <- metadata$cell_id
  groups <- as.list(unique(membership$membership))
  groups <- groups[!is.na(groups)]
  
  chr_groups <- as.list(unique(genomechunks$chr))
  
  # STEP 3: Set up parallel processing (exactly like amethyst)
  if (threads > 1) {
    future::plan(future::multicore, workers = threads)
  }
  
  # STEP 4: Process by chromosome (exactly like amethyst approach)
  by_chr <- list()
  
  for (chr in chr_groups) {
    cat("Processing", chr, "\n")
    
    # Get chromosome index (exactly like amethyst)
    sites <- chr_index[[chr]]
    
    # Process groups in parallel (exactly like amethyst)
    chr_group_results <- furrr::future_map(.x = groups, .f = function(gr) {
      
      # Get member cells (exactly like amethyst)
      member_cells <- rownames(membership[membership$membership == gr, , drop = FALSE])
      
      # Get paths for member cells
      member_paths <- h5_paths[barcode %in% member_cells, path]
      member_barcodes <- h5_paths[barcode %in% member_cells, barcode]
      
      # Process each member cell (exactly like amethyst)
      member_results <- furrr::future_pmap(.l = list(member_paths, member_barcodes), 
                                           .f = function(path, barcode) {
                                             tryCatch({
                                               # Get site info for this cell and chromosome
                                               cell_sites <- sites[cell_id == barcode]
                                               
                                               if (nrow(cell_sites) == 0) {
                                                 return(data.table(window = character(0), c = integer(0), t = integer(0)))
                                               }
                                               
                                               # Read H5 data (exactly like amethyst)
                                               data <- data.table(h5read(path, name = paste0(type, "/", barcode, "/1"),
                                                                         start = cell_sites$start,
                                                                         count = cell_sites$count))
                                               
                                               # Apply amethyst's exact window assignment logic
                                               data[pos %% step == 0, pos := pos + 1L]
                                               data[, window := paste0(chr, "_", 
                                                                       plyr::round_any(pos, step, floor), "_",
                                                                       plyr::round_any(pos, step, ceiling))]
                                               data[, c("chr", "pos", "pct") := NULL]
                                               
                                               # Aggregate by window
                                               data[, .(c = sum(c, na.rm = TRUE), t = sum(t, na.rm = TRUE)), by = window]
                                               
                                             }, error = function(e) {
                                               warning(sprintf("Error processing %s %s: %s", chr, barcode, e$message))
                                               return(data.table(window = character(0), c = integer(0), t = integer(0)))
                                             })
                                           }, .progress = TRUE)
      
      # STEP 5: Aggregate in chunks (exactly like amethyst)
      member_results <- split(member_results, ceiling(seq_along(member_results) / 100))
      member_results <- lapply(member_results, function(chunk) {
        chunk <- rbindlist(chunk, fill = TRUE)
        if (nrow(chunk) > 0) {
          chunk[, .(c = sum(c, na.rm = TRUE), t = sum(t, na.rm = TRUE)), by = window]
        } else {
          data.table(window = character(0), c = integer(0), t = integer(0))
        }
      })
      
      member_results <- rbindlist(member_results, fill = TRUE)
      
      if (nrow(member_results) > 0) {
        member_results <- member_results[, .(c = sum(c, na.rm = TRUE), t = sum(t, na.rm = TRUE)), by = window]
        setnames(member_results, c("window", paste0(gr, "_c"), paste0(gr, "_t")))
      } else {
        member_results <- data.table(window = character(0))
        member_results[, paste0(gr, "_c") := integer(0)]
        member_results[, paste0(gr, "_t") := integer(0)]
      }
      
      return(member_results)
      
    }, .progress = TRUE)
    
    # Merge group results for this chromosome (exactly like amethyst)
    by_chr[[chr]] <- Reduce(function(x, y) merge(x, y, by = "window", all = TRUE, sort = FALSE), 
                            chr_group_results)
    
    cat("Completed", chr, "\n")
  }
  
  # STEP 6: Reset threading and combine results (exactly like amethyst)
  if (threads > 1) {
    future::plan(future::sequential)
    gc()
  }
  
  count_matrix <- rbindlist(by_chr, fill = TRUE)
  
  # STEP 7: Join with genomic chunks (exactly like amethyst)
  count_matrix <- merge(genomechunks, count_matrix, by = "window", all.x = TRUE, sort = FALSE)
  
  # STEP 8: Apply smoothing (exactly like amethyst)
  count_matrix[, window := NULL]
  setorder(count_matrix, chr, start)
  
  # Calculate rolling sums (exactly like amethyst)
  count_cols <- grep("_c$|_t$", names(count_matrix), value = TRUE)
  rolling_sums <- count_matrix[, lapply(.SD, function(x) 
    as.integer(frollsum(x, n = smooth, align = "center", fill = NA, na.rm = TRUE))), 
    .SDcols = count_cols]
  
  # Calculate smoothed boundaries (exactly like amethyst)
  count_matrix[, smooth_start := lapply(.SD, function(x) 
    frollapply(x, n = smooth, FUN = function(y) y[1], align = "center", fill = NA)), 
    .SDcols = "start"]
  count_matrix[, smooth_end := lapply(.SD, function(x) 
    frollapply(x, n = smooth, FUN = function(y) y[smooth], align = "center", fill = NA)), 
    .SDcols = "end"]
  
  # Replace with smoothed values (exactly like amethyst)
  count_matrix[, names(rolling_sums) := rolling_sums]
  
  # Apply amethyst's exact edge filtering
  count_matrix <- count_matrix[shift(chr, 1) == shift(chr, -2)]
  
  # Clean up (exactly like amethyst)
  count_matrix[, c("smooth_start", "smooth_end") := NULL]
  count_matrix <- count_matrix[rowSums(count_matrix[, .SD, .SDcols = -c("chr", "start", "end")], na.rm = TRUE) > 0]
  
  cat(sprintf("Final count matrix: %d windows\n", nrow(count_matrix)))
  
  # STEP 9: Calculate percentage matrix (exactly like amethyst)
  pct_matrix <- copy(count_matrix)
  
  for (gr in groups) {
    c_col <- paste0(gr, "_c")
    t_col <- paste0(gr, "_t")
    
    if (c_col %in% names(pct_matrix) && t_col %in% names(pct_matrix)) {
      # Exactly like amethyst calculation
      pct_matrix[, m := round(pct_matrix[[c_col]] * 100 / 
                                (pct_matrix[[c_col]] + pct_matrix[[t_col]]), 2)]
      setnames(pct_matrix, "m", as.character(gr))
    }
  }
  
  # Remove count columns (exactly like amethyst)
  pct_matrix[, (paste0(groups, "_c")) := NULL]
  pct_matrix[, (paste0(groups, "_t")) := NULL]
  
  cat(sprintf("Final percentage matrix: %d windows\n", nrow(pct_matrix)))
  
  return(list(sum_matrix = count_matrix, pct_matrix = pct_matrix))
}

# Helper function for plyr dependency
if (!requireNamespace("plyr", quietly = TRUE)) {
  round_any <- function(x, accuracy, f = round) {
    f(x / accuracy) * accuracy
  }
} else {
  round_any <- plyr::round_any
}

#' Helper function to get chromosome sizes
get_chromosome_sizes <- function(genome) {
  if (genome == "hg38") {
    data.table(
      chr = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
              "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
              "chr20", "chr21", "chr22", "chrX", "chrY"),
      size = c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636,
               138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345,
               83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 57227415)
    )
  } else if (genome == "hg19") {
    data.table(
      chr = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
              "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
              "chr20", "chr21", "chr22", "chrX", "chrY"),
      size = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431,
               135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
               59128983, 63025520, 48129895, 51304566, 155270560, 59373566)
    )
  } else {
    stop("Only hg19 and hg38 supported. Use 'hg38', 'hg19', 'mm10', or 'mm39'")
  }
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


#' Calculate percentage methylation matrix (FIXED)
#' 
#' @param count_matrix data.table with count data
#' @param groups character vector of group names  
#' @return data.table with percentage methylation
calculate_pct_matrix <- function(count_matrix, groups) {
  pct_matrix <- copy(count_matrix)
  
  # FIXED: Calculate all percentages first, then remove count columns
  for (grp in groups) {
    c_col <- paste0(grp, "_c")
    t_col <- paste0(grp, "_t")
    
    if (c_col %in% names(pct_matrix) && t_col %in% names(pct_matrix)) {
      pct_matrix[, (grp) := round(100 * get(c_col) / (get(c_col) + get(t_col)), 2)]
    }
  }
  
  # Remove count columns AFTER all calculations are done
  count_cols <- grep("_c$|_t$", names(pct_matrix), value = TRUE)
  pct_matrix[, (count_cols) := NULL]
  
  return(pct_matrix)
}