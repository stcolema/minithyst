# windowing.R - Smoothed window calculation functions for minithyst package

#' @title Calculate smoothed methylation windows
#' 
#' @description Aggregates methylation data into genomic windows and applies kernel smoothing
#' to reduce sparsity. Cells are pooled within groups before smoothing.
#' 
#' @param h5_paths data.table with columns: barcode, path
#' @param chr_index Optimized index from create_optimized_index()
#' @param metadata data.table with cell metadata including grouping column
#' @param type Character, methylation context ("CG" or "CH")
#' @param step Integer, window step size in bp (default: 500)
#' @param smooth Integer, kernel size for smoothing (default: 3)
#' @param kernel Character or numeric vector specifying kernel type:
#'   "uniform" (default), "gaussian", "triangular", or custom weights
#' @param group_by Character, column name in metadata for grouping (default: "cluster_id")
#' @param genome Character, genome build (default: "hg38")
#' @param threads Integer, number of threads (default: 1)
#' @param chrList Optional character vector of chromosomes to process
#' @param chrSizes Optional numeric vector matching chrList
#' @param return_sum_matrix Logical, return count matrix (default: TRUE)
#' @param return_pct_matrix Logical, return percentage matrix (default: TRUE)
#' @return List with sum_matrix and/or pct_matrix based on parameters
#' @export
#' @examples
#' \dontrun{
#' # Create index first
#' index <- create_optimized_index(h5_paths, type = "CG")
#' 
#' # Calculate smoothed windows
#' result <- calc_smoothed_windows(
#'   h5_paths = h5_paths,
#'   chr_index = index,
#'   metadata = metadata,
#'   group_by = "cluster_id",
#'   threads = 4
#' )
#' }
calc_smoothed_windows <- function(h5_paths, 
                                  chr_index,
                                  metadata,
                                  type = "CG",
                                  step = 500,
                                  smooth = 3,
                                  kernel = "uniform",
                                  group_by = "cluster_id",
                                  genome = "hg38",
                                  threads = 1,
                                  chrList = NULL,
                                  chrSizes = NULL,
                                  return_sum_matrix = TRUE,
                                  return_pct_matrix = TRUE) {
  
  # Input validation
  if (!data.table::is.data.table(h5_paths)) {
    h5_paths <- data.table::as.data.table(h5_paths)
  }
  if (!data.table::is.data.table(metadata)) {
    metadata <- data.table::as.data.table(metadata)
  }
  
  # Validate required columns
  if (!all(c("barcode", "path") %in% names(h5_paths))) {
    stop("h5_paths must contain 'barcode' and 'path' columns")
  }
  if (!group_by %in% names(metadata)) {
    stop(sprintf("Column '%s' not found in metadata", group_by))
  }
  if (!inherits(chr_index, "optimized_methylation_index")) {
    stop("chr_index must be created with create_optimized_index()")
  }
  if (!is.null(chrList) && is.null(chrSizes)) {
    stop("chrSizes must be provided with chrList")
  }
  if (!return_sum_matrix && !return_pct_matrix) {
    stop("At least one output matrix must be requested")
  }
  
  # Ensure cell_id column exists
  if (!"cell_id" %in% names(metadata)) {
    if (!is.null(rownames(metadata))) {
      metadata[, cell_id := rownames(metadata)]
    } else {
      stop("metadata must have 'cell_id' column or row names")
    }
  }
  
  # Get kernel weights
  kernel_weights <- get_kernel_weights(kernel, smooth)
  
  # Configure threading
  old_dt_threads <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(old_dt_threads))
  data.table::setDTthreads(threads)
  
  # Get groups
  groups <- unique(metadata[[group_by]])
  groups <- groups[!is.na(groups)]
  
  message(sprintf("Processing %d groups with %d threads", 
                  length(groups), threads))
  
  # Generate genomic windows
  genomechunks <- generate_genomic_windows(
    genome = genome,
    step = step,
    chrList = chrList,
    chrSizes = chrSizes
  )
  
  # Process data
  count_matrix <- process_with_optimized_index(
    index = chr_index,
    metadata = metadata,
    groups = groups,
    group_by = group_by,
    step = step
  )
  
  count_matrix <- reduce_to_real_columns(count_matrix)
  
  # Merge with genomic windows
  count_matrix <- merge_with_windows(genomechunks, count_matrix)
  
  # Apply kernel smoothing
  count_matrix <- apply_kernel_smoothing(
    count_matrix = count_matrix,
    smooth = smooth,
    kernel = kernel,
    kernel_weights = kernel_weights
  )
  
  message(sprintf("Final count matrix: %d windows", nrow(count_matrix)))
  
  # Calculate percentage matrix if requested
  pct_matrix <- NULL
  if (return_pct_matrix) {
    pct_matrix <- calculate_percentage_matrix(count_matrix, groups)
    message(sprintf("Final percentage matrix: %d windows", nrow(pct_matrix)))
  }
  
  # Return results
  if (return_sum_matrix && return_pct_matrix) {
    return(list(sum_matrix = count_matrix, pct_matrix = pct_matrix))
  } else if (return_sum_matrix) {
    return(count_matrix)
  } else {
    return(pct_matrix)
  }
}

#' Process data using optimized index
#' @keywords internal
process_with_optimized_index <- function(index, metadata, groups, group_by, step) {
  
  # Get all chromosomes from index
  all_chr <- unique(unlist(lapply(index$cell_index, function(x) {
    if (!is.null(x$summary) && nrow(x$summary) > 0) {
      x$summary$chr
    } else {
      character(0)
    }
  })))
  
  if (length(all_chr) == 0) {
    stop("No chromosomes found in index")
  }
  
  all_chr <- sort(all_chr)
  message(sprintf("Found %d chromosomes: %s", length(all_chr), 
                  paste(head(all_chr, 5), collapse = ", ")))
  
  # Process each group-chromosome combination
  process_group_chr <- function(gr, chr) {
    # Get member cells
    member_cells <- metadata[get(group_by) == gr, cell_id]
    if (length(member_cells) == 0) return(NULL)
    
    # Find cells in index
    cells_in_index <- intersect(member_cells, names(index$cell_index))
    if (length(cells_in_index) == 0) return(NULL)
    
    # Check which cells have data for this chromosome
    cells_with_chr <- sapply(cells_in_index, function(cell) {
      cell_info <- index$cell_index[[cell]]
      if (!is.null(cell_info$summary) && nrow(cell_info$summary) > 0) {
        chr %in% cell_info$summary$chr
      } else {
        FALSE
      }
    })
    
    cells_to_process <- cells_in_index[cells_with_chr]
    if (length(cells_to_process) == 0) return(NULL)
    
    # Read data using optimized index
    group_data <- read_with_index(index, cells_to_process, chr)
    if (is.null(group_data) || nrow(group_data) == 0) return(NULL)
    
    # Window assignment
    group_data[pos %% step == 0, pos := pos + 1L]
    group_data[, window := paste0(chr, "_", 
                                  round_any(pos, step, floor), "_",
                                  round_any(pos, step, ceiling))]
    
    # Aggregate
    result <- group_data[, .(c = sum(c, na.rm = TRUE), 
                             t = sum(t, na.rm = TRUE)), 
                         by = window]
    
    data.table::setnames(result, c("window", paste0(gr, "_c"), paste0(gr, "_t")))
    return(result)
  }
  
  # Process all combinations
  all_results <- list()
  pb <- utils::txtProgressBar(max = length(all_chr) * length(groups), style = 3)
  counter <- 0
  
  for (chr in all_chr) {
    for (gr in groups) {
      counter <- counter + 1
      utils::setTxtProgressBar(pb, counter)
      
      result <- process_group_chr(gr, chr)
      if (!is.null(result) && nrow(result) > 0) {
        result_key <- paste0(chr, "_", gr)
        all_results[[result_key]] <- result
      }
    }
  }
  close(pb)
  
  if (length(all_results) == 0) {
    stop("No data was successfully processed")
  }
  
  # Merge all results
  count_matrix <- Reduce(function(x, y) {
    merge(x, y, by = "window", all = TRUE)
  }, all_results)
  
  return(count_matrix)
}

#' Reduce count matrix to the actual columns removing columns created by merging
#' @keywords internal
reduce_to_real_columns <- function(dt) {
  # Extract column names (excluding key columns like 'window')
  cols <- setdiff(names(dt), key(dt))
  
  # Get the prefixes
  prefixes <- sub("\\..*$", "", cols)
  
  # Group columns by prefix
  grouped_cols <- split(cols, prefixes)
  
  # Compute row sums by group
  for (prefix in names(grouped_cols)) {
    cols_rel <- grouped_cols[[prefix]]
    dt[[prefix]] <- rowSums(dt[, ..cols_rel], na.rm = TRUE)
  }
  cols_kept <- c(key(dt), unique(prefixes))
  out <- dt[, ..cols_kept]
  return(out)
}

#' Merge count matrix with genomic windows
#' @keywords internal
merge_with_windows <- function(genomechunks, count_matrix) {
  # Use dplyr if available for exact compatibility
  if (FALSE) { # (requireNamespace("dplyr", quietly = TRUE)) {
    genomechunks_df <- as.data.frame(genomechunks)
    count_matrix_df <- as.data.frame(count_matrix)
    merged <- dplyr::left_join(genomechunks_df, count_matrix_df, by = "window")
    return(data.table::as.data.table(merged))
  } else {
    # Use data.table merge
    
    out <- count_matrix[genomechunks, on = "window"]
    non_chr_cols <- unique(c(c("chr", "start", "end"), names(out)))
    return(out[, ..non_chr_cols])
    # return(genomechunks[count_matrix, on = "window"])
  }
}

#' Apply kernel smoothing to count matrix
#' @keywords internal
apply_kernel_smoothing <- function(count_matrix, smooth, kernel, kernel_weights) {
  
  count_matrix[, window := NULL]
  data.table::setorder(count_matrix, chr, start)
  
  # Get count columns
  count_cols <- grep("_c$|_t$", names(count_matrix), value = TRUE)
  
  if (length(count_cols) > 0) {
    if (kernel == "uniform" || identical(kernel_weights, rep(1, smooth))) {
      # Use optimized frollsum for uniform kernel
      rolling_values <- count_matrix[, lapply(.SD, function(x) {
        data.table::frollsum(x, n = smooth, align = "center", 
                             fill = NA, na.rm = TRUE)
      }), .SDcols = count_cols]
    } else {
      # Use frollapply with custom kernel
      rolling_values <- count_matrix[, lapply(.SD, function(x) {
        data.table::frollapply(x, n = smooth, FUN = function(window_vals) {
          valid_idx <- !is.na(window_vals)
          if (sum(valid_idx) == 0) return(NA_real_)
          
          valid_weights <- kernel_weights[valid_idx]
          valid_vals <- window_vals[valid_idx]
          valid_weights <- valid_weights / sum(valid_weights)
          
          sum(valid_vals * valid_weights)
        }, align = "center", fill = NA)
      }), .SDcols = count_cols]
    }
    
    # Update with smoothed values
    count_matrix[, (count_cols) := rolling_values]
    
    # Remove edge windows
    if (smooth == 3) {
      count_matrix <- count_matrix[data.table::shift(chr, 1) == 
                                     data.table::shift(chr, -2)]
    } else {
      edge_offset <- floor(smooth / 2)
      keep_idx <- rep(TRUE, nrow(count_matrix))
      
      for (i in 1:edge_offset) {
        keep_idx <- keep_idx & 
          (data.table::shift(count_matrix$chr, i) == count_matrix$chr) & 
          (data.table::shift(count_matrix$chr, -i) == count_matrix$chr)
      }
      count_matrix <- count_matrix[keep_idx]
    }
    
    # Remove empty windows
    non_chr_cols <- setdiff(names(count_matrix), c("chr", "start", "end"))
    if (length(non_chr_cols) > 0) {
      row_sums <- rowSums(count_matrix[, ..non_chr_cols], na.rm = TRUE)
      count_matrix <- count_matrix[row_sums != 0]
    }
  }
  
  return(count_matrix)
}

#' Calculate percentage matrix from count matrix
#' @keywords internal
calculate_percentage_matrix <- function(count_matrix, groups) {
  pct_matrix <- data.table::copy(count_matrix)
  
  for (gr in groups) {
    c_col <- paste0(gr, "_c")
    t_col <- paste0(gr, "_t")
    
    if (c_col %in% names(pct_matrix) && t_col %in% names(pct_matrix)) {
      pct_matrix[, (as.character(gr)) := round(
        get(c_col) * 100 / (get(c_col) + get(t_col)), 2
      )]
    }
  }
  
  # Remove count columns
  cols_to_remove <- c(paste0(groups, "_c"), paste0(groups, "_t"))
  cols_to_remove <- cols_to_remove[cols_to_remove %in% names(pct_matrix)]
  if (length(cols_to_remove) > 0) {
    pct_matrix[, (cols_to_remove) := NULL]
  }
  
  return(pct_matrix)
}

#' Generate genomic windows
#' @keywords internal
generate_genomic_windows <- function(genome, step, chrList = NULL, chrSizes = NULL) {
  
  # Get chromosome information
  if (is.null(chrList)) {
    genome_info <- get_genome_info(genome)
    chromosome_sizes <- data.table::data.table(
      chromosome = genome_info$chromosomes,
      size = genome_info$sizes
    )
  } else {
    chromosome_sizes <- data.table::data.table(
      chromosome = chrList,
      size = chrSizes
    )
  }
  
  # Generate windows for each chromosome
  window_list <- lapply(seq_len(nrow(chromosome_sizes)), function(i) {
    chr <- chromosome_sizes$chromosome[i]
    size <- chromosome_sizes$size[i]
    
    starts <- seq(0, size - 1, by = step)
    ends <- pmin(starts + step, size)
    
    data.table::data.table(
      chr = chr,
      start = starts,
      end = ends
    )
  })
  
  genomechunks <- data.table::rbindlist(window_list)
  genomechunks[, window := paste0(chr, "_", start, "_", end)]
  data.table::setkey(genomechunks, chr, start, end)
  
  return(genomechunks)
}

#' Get kernel weights for smoothing
#' @keywords internal
get_kernel_weights <- function(kernel, smooth) {
  if (is.numeric(kernel)) {
    if (length(kernel) != smooth) {
      stop(sprintf("Custom kernel must have length %d", smooth))
    }
    return(kernel * smooth / sum(kernel))
  }
  
  switch(kernel,
         "uniform" = rep(1, smooth),
         "gaussian" = {
           center <- (smooth + 1) / 2
           positions <- 1:smooth
           sigma <- smooth / 4
           weights <- exp(-0.5 * ((positions - center) / sigma)^2)
           weights * smooth / sum(weights)
         },
         "triangular" = {
           center <- (smooth + 1) / 2
           positions <- 1:smooth
           weights <- 1 - abs(positions - center) / center
           weights * smooth / sum(weights)
         },
         stop(sprintf("Unknown kernel: %s. Options: uniform, gaussian, triangular, or numeric vector", kernel))
  )
}

#' Get genome information
#' @keywords internal
get_genome_info <- function(genome) {
  switch(genome,
         "hg19" = list(
           chromosomes = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
                           "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", 
                           "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
                           "chr20", "chr21", "chr22", "chrX", "chrY"),
           sizes = c(249250621, 243199373, 198022430, 191154276, 180915260, 
                     171115067, 159138663, 146364022, 141213431, 135534747, 
                     135006516, 133851895, 115169878, 107349540, 102531392, 
                     90354753, 81195210, 78077248, 59128983, 63025520, 
                     48129895, 51304566, 155270560, 59373566)
         ),
         "hg38" = list(
           chromosomes = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
                           "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", 
                           "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
                           "chr20", "chr21", "chr22", "chrX", "chrY"),
           sizes = c(248956422, 242193529, 198295559, 190214555, 181538259, 
                     170805979, 159345973, 145138636, 138394717, 133797422, 
                     135086622, 133275309, 114364328, 107043718, 101991189, 
                     90338345, 83257441, 80373285, 58617616, 64444167, 
                     46709983, 50818468, 156040895, 57227415)
         ),
         "mm10" = list(
           chromosomes = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
                           "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", 
                           "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
                           "chrX", "chrY"),
           sizes = c(195471971, 182113224, 160039680, 156508116, 151834684, 
                     149736546, 145441459, 129401213, 124595110, 130694993, 
                     122082543, 120129022, 120421639, 124902244, 104043685, 
                     98207768, 94987271, 90702639, 61431566, 171031299, 91744698)
         ),
         "mm39" = list(
           chromosomes = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
                           "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", 
                           "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
                           "chrX", "chrY"),
           sizes = c(195154279, 181755017, 159745316, 156860686, 151758149, 
                     149588044, 144995196, 130127694, 124359700, 130530862, 
                     121973369, 120092757, 120883175, 125139656, 104073951, 
                     98008968, 95294699, 90720763, 61420004, 169476592, 91455967)
         ),
         stop(sprintf("Unsupported genome: %s. Use chrList and chrSizes for custom genomes.", genome))
  )
}

#' Internal rounding function
#' @keywords internal
round_any <- function(x, accuracy, f = round) {
  f(x / accuracy) * accuracy
}