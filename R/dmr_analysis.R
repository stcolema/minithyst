# dmr_analysis.R - Simplified DMR Analysis from Single-Cell Methylation Data
# Focused reimplementation of Amethyst DMR functionality

#' 
#' #' @import data.table
#' #' @import GenomicRanges
#' #' @importFrom rhdf5 h5read h5ls H5Fopen H5Dopen H5Dget_space H5Sget_simple_extent_dims H5Dget_type H5Tget_class h5closeAll
#' #' @importFrom IRanges IRanges
#' #' @importFrom stats phyper dhyper p.adjust
#' NULL
#' 
#' # Set data.table threads explicitly to avoid conflicts
#' data.table::setDTthreads(1) # Will be configured by user
#' 
#' #' Index chromosomes in H5 files
#' #'
#' #' @param h5_paths data.table with columns: barcode, path
#' #' @param type Character, methylation context (e.g., "CG", "CH")
#' #' @param chr_list Optional character vector of chromosomes to index
#' #' @param threads Integer, number of threads for parallel processing
#' #' @return List of data.tables by chromosome with h5 coordinates
#' #' @export
#' index_chromosomes <- function(h5_paths, type = "CG", chr_list = NULL, threads = 1) {
#'   # Configure threading
#'   old_threads <- data.table::getDTthreads()
#'   data.table::setDTthreads(threads)
#'   on.exit(data.table::setDTthreads(old_threads))
#' 
#'   # Process each file
#'   all_indices <- data.table::rbindlist(lapply(seq_len(nrow(h5_paths)), function(i) {
#'     path <- h5_paths$path[i]
#'     barcode <- h5_paths$barcode[i]
#' 
#'     tryCatch(
#'       {
#'         # Build the data path
#'         data_path <- paste0(type, "/", barcode, "/1")
#' 
#'         # Read the full dataset (1D compound array)
#'         h5_data <- rhdf5::h5read(path, name = data_path)
#' 
#'         # Convert to data.table if needed
#'         if (!data.table::is.data.table(h5_data)) {
#'           h5_data <- data.table::as.data.table(h5_data)
#'         }
#' 
#'         # Get unique chromosomes and their positions
#'         chr_summary <- h5_data[, .(
#'           cell_id = barcode,
#'           first_pos = min(.I),
#'           last_pos = max(.I),
#'           n_sites = .N
#'         ), by = chr]
#' 
#'         # Filter chromosomes if list provided
#'         if (!is.null(chr_list)) {
#'           chr_summary <- chr_summary[chr %in% chr_list]
#'         } else {
#'           # Remove alternative contigs by default
#'           chr_summary <- chr_summary[!grepl("_|EBV|M", chr)]
#'         }
#' 
#'         return(chr_summary)
#'       },
#'       error = function(e) {
#'         warning(sprintf("Error indexing barcode %s: %s", barcode, e$message))
#'         return(NULL)
#'       }
#'     )
#'   }))
#' 
#'   # Split by chromosome
#'   if (nrow(all_indices) > 0) {
#'     chr_index <- split(all_indices, all_indices$chr)
#'   } else {
#'     chr_index <- list()
#'   }
#' 
#'   return(chr_index)
#' }
#' 
#' 
#' #' Calculate smoothed methylation windows
#' #' 
#' #' @param h5_paths data.table with columns: barcode, path
#' #' @param chr_index List of chromosome indices from index_chromosomes
#' #' @param metadata data.table with cell metadata, must have grouping column
#' #' @param type Character, methylation context ("CG" or "CH")
#' #' @param step Integer, window step size in bp (default: 500)
#' #' @param smooth Integer, kernel size for smoothing adjacent windows (default: 3)
#' #' @param group_by Character, column name in metadata for grouping (default: "cluster_id")
#' #' @param genome Character, genome build for chromosome sizes (default: "hg38")
#' #' @param threads Integer, number of threads (default: 1)
#' #' @param parallel_strategy Character, either "chromosomes" or "groups" (default: "auto")
#' #' @return List with sum_matrix and/or pct_matrix
#' #' @export
#' calc_smoothed_windows <- function(h5_paths, 
#'                                   chr_index,
#'                                   metadata,
#'                                   type = "CG",
#'                                   step = 500,
#'                                   smooth = 3,
#'                                   group_by = "cluster_id",
#'                                   genome = "hg38",
#'                                   threads = 1,
#'                                   parallel_strategy = "auto",
#'                                   return_sum_matrix = TRUE,
#'                                   return_pct_matrix = TRUE) {
#'   
#'   # Input validation
#'   h5_paths <- data.table::as.data.table(h5_paths)
#'   metadata <- data.table::as.data.table(metadata)
#'   
#'   if (!all(c("barcode", "path") %in% names(h5_paths))) {
#'     stop("h5_paths must contain 'barcode' and 'path' columns")
#'   }
#'   
#'   # Prepare groups
#'   if (!"cell_id" %in% names(metadata)) {
#'     metadata[, cell_id := rownames(metadata)]
#'   }
#'   groups <- unique(metadata[[group_by]])
#'   groups <- groups[!is.na(groups)]
#'   n_groups <- length(groups)
#'   n_chr <- length(chr_index)
#'   
#'   # CRITICAL FIX 1: Choose parallelization strategy intelligently
#'   if (parallel_strategy == "auto") {
#'     # Use groups if we have few chromosomes but many groups
#'     # Use chromosomes if we have many chromosomes but few groups
#'     parallel_strategy <- ifelse(n_groups > n_chr * 2, "groups", "chromosomes")
#'   }
#'   
#'   # CRITICAL FIX 2: Configure threading properly
#'   old_dt_threads <- data.table::getDTthreads()
#'   on.exit(data.table::setDTthreads(old_dt_threads))
#'   
#'   if (threads == 1) {
#'     # Single threaded - keep everything sequential
#'     data.table::setDTthreads(1)
#'   } else if (parallel_strategy == "groups") {
#'     # Parallelize across groups, let data.table work single-threaded
#'     data.table::setDTthreads(1)
#'     future::plan(future::multisession, workers = threads)
#'     on.exit(future::plan(future::sequential), add = TRUE)
#'   } else {
#'     # Parallelize across chromosomes, give remaining threads to data.table
#'     dt_threads <- max(1, ceiling(threads / n_chr))
#'     data.table::setDTthreads(dt_threads)
#'   }
#'   
#'   # Generate genomic windows
#'   genomechunks <- generate_genomic_windows(genome, step)
#'   genomechunks[, window := paste0(chr, "_", start, "_", end)]
#'   data.table::setkey(genomechunks, chr, start, end)
#'   
#'   # CRITICAL FIX 3: Efficient processing function
#'   process_group_for_chromosome <- function(gr, chr, sites, h5_paths, metadata) {
#'     # Get cells in this group
#'     member_cells <- metadata[get(group_by) == gr, cell_id]
#'     if (length(member_cells) == 0) return(NULL)
#'     
#'     # Read all data for this group-chromosome combination at once
#'     group_data_list <- lapply(member_cells, function(cell_id) {
#'       # Match cell to h5 path
#'       path_row <- h5_paths[barcode == cell_id]
#'       if (nrow(path_row) == 0) return(NULL)
#'       
#'       # Get sites for this cell
#'       cell_sites <- sites[cell_id == cell_id]
#'       if (nrow(cell_sites) == 0) return(NULL)
#'       
#'       tryCatch({
#'         # Read H5 data
#'         data <- data.table::data.table(
#'           rhdf5::h5read(path_row$path[1], 
#'                         name = paste0(type, "/", cell_id, "/1"),
#'                         start = cell_sites$start,
#'                         count = cell_sites$count)
#'         )
#'         data[, cell_id := cell_id]
#'         return(data)
#'       }, error = function(e) NULL)
#'     })
#'     
#'     # Combine all cell data for this group
#'     group_data <- data.table::rbindlist(group_data_list, fill = TRUE)
#'     if (is.null(group_data) || nrow(group_data) == 0) return(NULL)
#'     
#'     # Window assignment
#'     group_data[pos %% step == 0, pos := pos + 1L]
#'     group_data[, window := paste0(chr, "_", 
#'                                   round_any(pos, step, floor), "_",
#'                                   round_any(pos, step, ceiling))]
#'     
#'     # Aggregate - this uses data.table's optimized grouping
#'     result <- group_data[, .(c = sum(c, na.rm = TRUE), 
#'                              t = sum(t, na.rm = TRUE)), 
#'                          by = window]
#'     
#'     data.table::setnames(result, c("window", paste0(gr, "_c"), paste0(gr, "_t")))
#'     return(result)
#'   }
#'   
#'   # CRITICAL FIX 4: Process based on chosen strategy
#'   if (parallel_strategy == "chromosomes" || threads == 1) {
#'     # Process chromosome by chromosome
#'     chr_results <- list()
#'     
#'     for (chr in names(chr_index)) {
#'       cat("Processing", chr, "\n")
#'       sites <- chr_index[[chr]]
#'       
#'       # Process all groups for this chromosome
#'       chr_group_results <- lapply(groups, function(gr) {
#'         process_group_for_chromosome(gr, chr, sites, h5_paths, metadata)
#'       })
#'       
#'       # Remove NULLs and merge
#'       chr_group_results <- chr_group_results[!sapply(chr_group_results, is.null)]
#'       if (length(chr_group_results) > 0) {
#'         chr_results[[chr]] <- Reduce(function(x, y) {
#'           base::merge(x, y, by = "window", all = TRUE, sort = FALSE)
#'         }, chr_group_results)
#'       }
#'     }
#'     
#'     count_matrix <- data.table::rbindlist(chr_results, fill = TRUE)
#'     
#'   } else {
#'     # Process groups in parallel across all chromosomes
#'     all_results <- furrr::future_map(groups, function(gr) {
#'       # Process this group for all chromosomes
#'       group_chr_results <- lapply(names(chr_index), function(chr) {
#'         sites <- chr_index[[chr]]
#'         process_group_for_chromosome(gr, chr, sites, h5_paths, metadata)
#'       })
#'       
#'       # Combine chromosomes for this group
#'       data.table::rbindlist(group_chr_results, fill = TRUE)
#'     }, .progress = TRUE)
#'     
#'     # Merge all groups
#'     all_results <- all_results[!sapply(all_results, is.null)]
#'     count_matrix <- Reduce(function(x, y) {
#'       base::merge(x, y, by = "window", all = TRUE, sort = FALSE)
#'     }, all_results)
#'   }
#'   
#'   # CRITICAL FIX 5: Ensure we use dplyr for compatibility
#'   if (requireNamespace("dplyr", quietly = TRUE)) {
#'     genomechunks_df <- as.data.frame(genomechunks)
#'     count_matrix_df <- as.data.frame(count_matrix)
#'     count_matrix <- dplyr::left_join(genomechunks_df, count_matrix_df, by = "window")
#'     count_matrix <- data.table::as.data.table(count_matrix)
#'   } else {
#'     # Fallback to data.table merge
#'     count_matrix <- genomechunks[count_matrix, on = "window"]
#'   }
#'   
#'   # Apply smoothing (rest remains the same)
#'   apply_smoothing(count_matrix, smooth, groups)
#' }
#' 
#' #' Helper to generate genomic windows
#' #' @keywords internal
#' generate_genomic_windows <- function(genome, step) {
#'   genome_info <- get_genome_info(genome)
#'   
#'   window_list <- lapply(seq_along(genome_info$chromosomes), function(i) {
#'     chr <- genome_info$chromosomes[i]
#'     size <- genome_info$sizes[i]
#'     starts <- seq(0, size - 1, by = step)
#'     ends <- pmin(starts + step, size)
#'     data.table::data.table(chr = chr, start = starts, end = ends)
#'   })
#'   
#'   data.table::rbindlist(window_list)
#' }
#' 
#' #' Apply smoothing to count matrix
#' #' @keywords internal  
#' apply_smoothing <- function(count_matrix, smooth, groups) {
#'   count_matrix[, window := NULL]
#'   data.table::setorder(count_matrix, chr, start)
#'   
#'   # Get count columns
#'   count_cols <- grep("_c$|_t$", names(count_matrix), value = TRUE)
#'   
#'   if (length(count_cols) > 0) {
#'     # Apply rolling sum
#'     rolling_sums <- count_matrix[, lapply(.SD, function(x) {
#'       data.table::frollsum(x, n = smooth, align = "center", fill = NA, na.rm = TRUE)
#'     }), .SDcols = count_cols]
#'     
#'     count_matrix[, (count_cols) := rolling_sums]
#'     
#'     # Edge filtering for smooth=3
#'     if (smooth == 3) {
#'       count_matrix <- count_matrix[data.table::shift(chr, 1) == data.table::shift(chr, -2)]
#'     } else {
#'       # General case
#'       edge_offset <- floor(smooth / 2)
#'       for (i in 1:edge_offset) {
#'         count_matrix <- count_matrix[
#'           data.table::shift(chr, i) == chr & 
#'             data.table::shift(chr, -i) == chr
#'         ]
#'       }
#'     }
#'     
#'     # Remove empty windows
#'     non_chr_cols <- setdiff(names(count_matrix), c("chr", "start", "end"))
#'     if (length(non_chr_cols) > 0) {
#'       row_sums <- rowSums(count_matrix[, ..non_chr_cols], na.rm = TRUE)
#'       count_matrix <- count_matrix[row_sums != 0]
#'     }
#'   }
#'   
#'   return(count_matrix)
#' }
#' 
#' 
#' 
#' #' Process a single group for a chromosome
#' #' @keywords internal
#' process_group_chromosome <- function(gr, chr, sites, membership, h5_paths, 
#'                                      type, step, smooth) {
#'   # Get member cells
#'   member_cells <- rownames(membership[membership$membership == gr, , drop = FALSE])
#'   
#'   if (length(member_cells) == 0) {
#'     return(data.table(window = character(0), c = integer(0), t = integer(0)))
#'   }
#'   
#'   # Handle prefix if present
#'   if ("prefix" %in% names(h5_paths)) {
#'     h5_paths$cell_id <- paste0(h5_paths$prefix, h5_paths$barcode)
#'     member_indices <- match(member_cells, h5_paths$cell_id)
#'   } else {
#'     member_indices <- match(member_cells, h5_paths$barcode)
#'   }
#'   
#'   valid_indices <- !is.na(member_indices)
#'   
#'   if (!any(valid_indices)) {
#'     return(data.table(window = character(0), c = integer(0), t = integer(0)))
#'   }
#'   
#'   member_paths <- h5_paths$path[member_indices[valid_indices]]
#'   member_barcodes <- h5_paths$barcode[member_indices[valid_indices]]
#'   member_cells_valid <- member_cells[valid_indices]
#'   
#'   # Process each member cell
#'   member_results <- furrr::future_pmap(
#'     list(member_paths, member_barcodes, member_cells_valid),
#'     function(path, barcode, cell_id) {
#'       tryCatch({
#'         # Get cell sites for this chromosome
#'         cell_sites <- sites[sites$cell_id == cell_id, ]
#'         
#'         if (nrow(cell_sites) == 0) {
#'           return(data.table(window = character(0), c = integer(0), t = integer(0)))
#'         }
#'         
#'         # Read H5 data
#'         data <- data.table(
#'           h5read(path, 
#'                  name = paste0(type, "/", barcode, "/1"),
#'                  start = cell_sites$start,
#'                  count = cell_sites$count)
#'         )
#'         
#'         # Apply window logic (exactly matching amethyst)
#'         data[pos %% step == 0, pos := pos + 1L]
#'         data[, window := paste0(chr, "_", 
#'                                 round_any(pos, step, floor), "_",
#'                                 round_any(pos, step, ceiling))]
#'         data[, c("chr", "pos", "pct") := NULL]
#'         
#'         # Aggregate by window
#'         data[, .(c = sum(c, na.rm = TRUE), t = sum(t, na.rm = TRUE)), by = window]
#'         
#'       }, error = function(e) {
#'         return(data.table(window = character(0), c = integer(0), t = integer(0)))
#'       })
#'     }, .progress = FALSE
#'   )
#'   
#'   # Aggregate results in chunks to avoid memory issues
#'   if (length(member_results) > 100) {
#'     member_results <- split(member_results, ceiling(seq_along(member_results) / 100))
#'     member_results <- lapply(member_results, function(chunk) {
#'       chunk_combined <- rbindlist(chunk, fill = TRUE)
#'       if (nrow(chunk_combined) > 0) {
#'         chunk_combined[, .(c = sum(c, na.rm = TRUE), t = sum(t, na.rm = TRUE)), by = window]
#'       } else {
#'         data.table(window = character(0), c = integer(0), t = integer(0))
#'       }
#'     })
#'   }
#'   
#'   member_results <- rbindlist(member_results, fill = TRUE)
#'   
#'   if (nrow(member_results) > 0) {
#'     member_results <- member_results[, .(c = sum(c, na.rm = TRUE), 
#'                                          t = sum(t, na.rm = TRUE)), by = window]
#'     setnames(member_results, c("window", paste0(gr, "_c"), paste0(gr, "_t")))
#'   } else {
#'     member_results <- data.table(window = character(0))
#'     member_results[, paste0(gr, "_c") := integer(0)]
#'     member_results[, paste0(gr, "_t") := integer(0)]
#'   }
#'   
#'   return(member_results)
#' }
#' 
#' #' Internal helper function for rounding (replaces plyr::round_any)
#' #' @param x numeric vector to round
#' #' @param accuracy number to round to  
#' #' @param f rounding function
#' #' @keywords internal
#' round_any <- function(x, accuracy, f = round) {
#'   f(x / accuracy) * accuracy
#' }
#' 
#' #' Get genome information
#' #' @param genome Character string specifying genome build
#' #' @return List with chromosomes and sizes
#' #' @keywords internal
#' get_genome_info <- function(genome) {
#'   if (genome == "hg19") {
#'     list(
#'       chromosomes = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
#'                       "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", 
#'                       "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
#'                       "chr20", "chr21", "chr22", "chrX", "chrY"),
#'       sizes = c(249250621, 243199373, 198022430, 191154276, 180915260, 
#'                 171115067, 159138663, 146364022, 141213431, 135534747, 
#'                 135006516, 133851895, 115169878, 107349540, 102531392, 
#'                 90354753, 81195210, 78077248, 59128983, 63025520, 
#'                 48129895, 51304566, 155270560, 59373566)
#'     )
#'   } else if (genome == "hg38") {
#'     list(
#'       chromosomes = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
#'                       "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", 
#'                       "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
#'                       "chr20", "chr21", "chr22", "chrX", "chrY"),
#'       sizes = c(248956422, 242193529, 198295559, 190214555, 181538259, 
#'                 170805979, 159345973, 145138636, 138394717, 133797422, 
#'                 135086622, 133275309, 114364328, 107043718, 101991189, 
#'                 90338345, 83257441, 80373285, 58617616, 64444167, 
#'                 46709983, 50818468, 156040895, 57227415)
#'     )
#'   } else if (genome == "mm10") {
#'     list(
#'       chromosomes = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
#'                       "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", 
#'                       "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
#'                       "chrX", "chrY"),
#'       sizes = c(195471971, 182113224, 160039680, 156508116, 151834684, 
#'                 149736546, 145441459, 129401213, 124595110, 130694993, 
#'                 122082543, 120129022, 120421639, 124902244, 104043685, 
#'                 98207768, 94987271, 90702639, 61431566, 171031299, 91744698)
#'     )
#'   } else if (genome == "mm39") {
#'     list(
#'       chromosomes = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
#'                       "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", 
#'                       "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
#'                       "chrX", "chrY"),
#'       sizes = c(195154279, 181755017, 159745316, 156860686, 151758149, 
#'                 149588044, 144995196, 130127694, 124359700, 130530862, 
#'                 121973369, 120092757, 120883175, 125139656, 104073951, 
#'                 98008968, 95294699, 90720763, 61420004, 169476592, 91455967)
#'     )
#'   } else {
#'     stop("Only hg19, hg38, mm10, or mm39 are supported. For other genomes, use chrList and chrSizes parameters.")
#'   }
#' }

#' Helper function to get chromosome sizes
get_chromosome_sizes <- function(genome) {
  if (genome == "hg38") {
    data.table(
      chr = c(
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
        "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
        "chr20", "chr21", "chr22", "chrX", "chrY"
      ),
      size = c(
        248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636,
        138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345,
        83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 57227415
      )
    )
  } else if (genome == "hg19") {
    data.table(
      chr = c(
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
        "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
        "chr20", "chr21", "chr22", "chrX", "chrY"
      ),
      size = c(
        249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431,
        135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
        59128983, 63025520, 48129895, 51304566, 155270560, 59373566
      )
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

  # # Filter windows with no data
  # counts <- counts[!is.na(window)]

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
      # passing_windows <- (counts$member_c + counts$member_t < min_group) |
      #   (counts$nonmember_c + counts$nonmember_t < min_group)
      counts[
        member_c + member_t < min_group |
          nonmember_c + nonmember_t < min_group,
        c("member_c", "member_t", "nonmember_c", "nonmember_t") := NA
      ]

      # Apply Fisher's test
      counts[, paste0(grp, "_pval") := apply(.SD, 1, function(x) {
        if (any(is.na(x))) {
          return(NA_real_)
        }
        fast_fisher(matrix(x, nrow = 2, byrow = TRUE))
      }), .SDcols = c("member_c", "member_t", "nonmember_c", "nonmember_t")]

      # Calculate log fold change
      counts[, paste0(grp, "_logFC") := log2(
        (member_c / (member_c + member_t)) /
          (nonmember_c / (nonmember_c + nonmember_t))
      )]
      
      # Clean up temporary columns
      counts[, c("member_c", "member_t", "nonmember_c", "nonmember_t") := NULL]
      
    }

    
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
      counts[
        member_c + member_t < min_group |
          nonmember_c + nonmember_t < min_group,
        c("member_c", "member_t", "nonmember_c", "nonmember_t") := NA
      ]

      counts[, paste0(name, "_pval") := apply(.SD, 1, function(x) {
        if (any(is.na(x))) {
          return(NA_real_)
        }
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
    variable.name = "test"
  )

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
  collapsed <- results[,
    {
      # Create ranges
      gr <- GenomicRanges::GRanges(
        seqnames = chr,
        ranges = IRanges::IRanges(
          start = start - max_gap,
          end = end + max_gap
        )
      )

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
    },
    by = .(test, direction)
  ]

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
      size = c(
        248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
        159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
        114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
        58617616, 64444167, 46709983, 50818468, 156040895, 57227415
      )
    ),
    "hg19" = data.table::data.table(
      chr = paste0("chr", c(1:22, "X", "Y")),
      size = c(
        249250621, 243199373, 198022430, 191154276, 180915260, 171115067,
        159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
        115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
        59128983, 63025520, 48129895, 51304566, 155270560, 59373566
      )
    ),
    "mm10" = data.table::data.table(
      chr = paste0("chr", c(1:19, "X", "Y")),
      size = c(
        195471971, 182113224, 160039680, 156508116, 151834684, 149736546,
        145441459, 129401213, 124595110, 130694993, 122082543, 120129022,
        120421639, 124902244, 104043685, 98207768, 94987271, 90702639,
        61431566, 171031299, 91744698
      )
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
