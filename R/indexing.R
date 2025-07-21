# indexing.R - Optimized indexing functions for minithyst package

#' Create an optimized index for efficient H5 file access
#' 
#' Creates a multi-level index that enables batch reading of contiguous genomic
#' regions and efficient window-based queries. This significantly improves I/O
#' performance compared to reading individual sites.
#' 
#' @param h5_paths data.table with columns: barcode, path
#' @param type Character, methylation context (e.g., "CG", "CH")
#' @param chr_list Optional character vector of chromosomes to index
#' @param window_size Integer, size of genomic windows for pre-indexing (default: 10000)
#' @param threads Integer, number of threads for parallel processing
#' @return S3 object of class "optimized_methylation_index"
#' @export
#' @examples
#' \dontrun{
#' h5_paths <- data.table(
#'   barcode = c("cell1", "cell2"), 
#'   path = c("data1.h5", "data2.h5")
#' )
#' index <- create_optimized_index(h5_paths, type = "CG", threads = 4)
#' }
create_optimized_index <- function(h5_paths, 
                                   type = "CG", 
                                   chr_list = NULL,
                                   window_size = 10000,
                                   threads = 1) {
  
  # Ensure data.table format
  if (!data.table::is.data.table(h5_paths)) {
    h5_paths <- data.table::as.data.table(h5_paths)
  }
  
  # Validate inputs
  if (!all(c("barcode", "path") %in% names(h5_paths))) {
    stop("h5_paths must contain 'barcode' and 'path' columns")
  }
  
  # Set up threading
  old_threads <- data.table::getDTthreads()
  data.table::setDTthreads(threads)
  on.exit(data.table::setDTthreads(old_threads))
  
  # Initialize index structure
  index <- list(
    cell_index = list(),
    window_index = list(),
    metadata = list(
      type = type,
      window_size = window_size,
      creation_time = Sys.time(),
      h5_paths = h5_paths
    )
  )
  
  # Process each cell
  message("Building optimized index...")
  pb <- utils::txtProgressBar(min = 0, max = nrow(h5_paths), style = 3)
  
  for (i in seq_len(nrow(h5_paths))) {
    utils::setTxtProgressBar(pb, i)
    
    barcode <- h5_paths$barcode[i]
    path <- h5_paths$path[i]
    
    tryCatch({
      # Read all data at once for indexing
      all_data <- data.table::data.table(
        rhdf5::h5read(path, name = paste0(type, "/", barcode, "/1"))
      )
      all_data[, row_idx := .I]
      
      # Filter chromosomes
      if (!is.null(chr_list)) {
        all_data <- all_data[chr %in% chr_list]
      } else {
        # Remove alternative contigs by default
        all_data <- all_data[!grepl("_|EBV|M", chr)]
      }
      
      if (nrow(all_data) == 0) next
      
      # Create cell-level index
      cell_chr_index <- all_data[, .(
        start_idx = min(row_idx),
        end_idx = max(row_idx),
        count = .N,
        min_pos = min(pos),
        max_pos = max(pos),
        n_blocks = sum(diff(row_idx) > 1) + 1
      ), by = chr]
      
      # Create window-level index
      all_data[, window := (pos %/% window_size) * window_size]
      window_index <- all_data[, .(
        start_idx = min(row_idx),
        end_idx = max(row_idx),
        n_sites = .N
      ), by = .(chr, window)]
      
      # Identify contiguous blocks for efficient reading
      blocks <- all_data[, {
        gaps <- which(diff(row_idx) > 1)
        block_starts <- c(1, gaps + 1)
        block_ends <- c(gaps, .N)
        
        list(
          block_id = seq_along(block_starts),
          block_start_idx = row_idx[block_starts],
          block_end_idx = row_idx[block_ends],
          block_start_pos = pos[block_starts],
          block_end_pos = pos[block_ends],
          block_size = block_ends - block_starts + 1
        )
      }, by = chr]
      
      # Store in index
      index$cell_index[[barcode]] <- list(
        summary = cell_chr_index,
        windows = window_index,
        blocks = blocks,
        path = path
      )
      
    }, error = function(e) {
      warning(sprintf("Error indexing %s: %s", barcode, e$message))
    })
  }
  
  close(pb)
  
  # Build global window registry
  message("Building global window index...")
  
  all_windows <- data.table::rbindlist(
    lapply(names(index$cell_index), function(cell) {
      windows <- index$cell_index[[cell]]$windows
      if (!is.null(windows)) {
        windows[, cell_id := cell]
      }
      windows
    }),
    fill = TRUE
  )
  
  if (nrow(all_windows) > 0) {
    index$window_index <- all_windows[, .(
      cells = list(unique(cell_id)),
      total_sites = sum(n_sites)
    ), by = .(chr, window)]
  }
  
  # Calculate optimal read strategies
  index$read_strategies <- calculate_read_strategies(index)
  
  class(index) <- "optimized_methylation_index"
  return(index)
}

#' Calculate optimal read strategies based on index
#' 
#' Analyzes data fragmentation to determine the most efficient read strategy
#' for each cell-chromosome combination.
#' 
#' @param index Optimized index object
#' @return List of read strategies by cell
#' @keywords internal
calculate_read_strategies <- function(index) {
  strategies <- list()
  
  for (cell in names(index$cell_index)) {
    cell_info <- index$cell_index[[cell]]
    if (is.null(cell_info$blocks)) next
    
    chr_strategies <- data.table::rbindlist(
      lapply(unique(cell_info$blocks$chr), function(chr) {
        blocks <- cell_info$blocks[chr == chr]
        summary <- cell_info$summary[chr == chr]
        
        if (nrow(summary) == 0) return(NULL)
        
        # Calculate read efficiency metrics
        total_sites <- summary$count
        n_blocks <- nrow(blocks)
        span <- summary$end_idx - summary$start_idx + 1
        coverage <- sum(blocks$block_size) / span
        
        # Decide strategy based on coverage and fragmentation
        if (all(coverage > 0.8) || n_blocks <= 3) {
          strategy <- "full_range"
          read_ops <- 1
        } else if (n_blocks <= 10) {
          strategy <- "blocks"
          read_ops <- n_blocks
        } else {
          strategy <- "windows"
          read_ops <- NA
        }
        
        data.table::data.table(
          chr = chr,
          strategy = strategy,
          n_blocks = n_blocks,
          coverage = coverage,
          read_ops = read_ops
        )
      }),
      fill = TRUE
    )
    
    if (nrow(chr_strategies) > 0) {
      strategies[[cell]] <- chr_strategies
    }
  }
  
  return(strategies)
}

#' Read methylation data using optimized index
#' 
#' Efficiently reads H5 data using the pre-computed index, minimizing I/O
#' operations by reading contiguous blocks when possible.
#' 
#' @param index Optimized index from create_optimized_index
#' @param cells Character vector of cell IDs to read
#' @param chr Chromosome to read
#' @param start Optional start position for region queries
#' @param end Optional end position for region queries
#' @return data.table with columns: chr, pos, c, t, pct, cell_id
#' @keywords internal
read_with_index <- function(index, cells, chr_of_interest, start = NULL, end = NULL) {
  
  if (!inherits(index, "optimized_methylation_index")) {
    stop("Index must be created with create_optimized_index()")
  }
  
  # Read data for each cell using optimal strategy
  cell_data_list <- lapply(cells, function(cell) {
    
    cell_info <- index$cell_index[[cell]]
    if (is.null(cell_info)) return(NULL)
    
    # Check if chromosome exists for this cell
    if (is.null(cell_info$summary) || nrow(cell_info$summary) == 0) return(NULL)
    if (!chr_of_interest %in% cell_info$summary$chr) return(NULL)
    
    strategy_info <- index$read_strategies[[cell]]
    if (is.null(strategy_info)) return(NULL)
    
    strategy_row <- strategy_info[chr == chr_of_interest]
    if (nrow(strategy_row) == 0) return(NULL)
    
    path <- cell_info$path
    type <- index$metadata$type
    
    # Choose read method based on strategy
    if (any(strategy_row$strategy == "full_range")) {
      # Read entire chromosome range at once
      summary <- cell_info$summary[chr == chr_of_interest]
      data <- data.table::data.table(
        rhdf5::h5read(path, 
                      name = paste0(type, "/", cell, "/1"),
                      start = summary$start_idx,
                      count = summary$count)
      )
      
    } else if (strategy_row$strategy == "blocks") {
      # Read blocks and combine
      blocks <- cell_info$blocks[chr == chr_of_interest]
      block_data <- lapply(seq_len(nrow(blocks)), function(i) {
        blk <- blocks[i]
        rhdf5::h5read(path,
                      name = paste0(type, "/", cell, "/1"),
                      start = blk$block_start_idx,
                      count = blk$block_size)
      })
      data <- data.table::as.data.table(do.call(rbind, block_data))
      
    } else {
      # Window-based reading for highly fragmented data
      windows_needed <- cell_info$windows[chr == chr_of_interest]
      
      if (!is.null(start) && !is.null(end)) {
        window_start <- (start %/% index$metadata$window_size) * 
          index$metadata$window_size
        window_end <- (end %/% index$metadata$window_size) * 
          index$metadata$window_size
        windows_needed <- windows_needed[window >= window_start & 
                                           window <= window_end]
      }
      
      if (nrow(windows_needed) == 0) return(NULL)
      
      # Read each window's data
      window_data <- lapply(seq_len(nrow(windows_needed)), function(i) {
        win <- windows_needed[i]
        rhdf5::h5read(path,
                      name = paste0(type, "/", cell, "/1"),
                      start = win$start_idx,
                      count = win$end_idx - win$start_idx + 1)
      })
      data <- data.table::as.data.table(do.call(rbind, window_data))
    }
    
    # Filter by position if requested
    if (!is.null(start) || !is.null(end)) {
      if (!is.null(start)) data <- data[pos >= start]
      if (!is.null(end)) data <- data[pos <= end]
    }
    
    if (nrow(data) > 0) {
      data[, cell_id := cell]
    }
    
    return(data)
  })
  
  # Combine all cell data
  cell_data <- data.table::rbindlist(
    cell_data_list[!sapply(cell_data_list, is.null)], 
    fill = TRUE
  )
  
  return(cell_data)
}

#' Save optimized index to disk
#' 
#' @param index Optimized index object
#' @param file Path to save the index (.rds extension recommended)
#' @export
save_index <- function(index, file) {
  if (!inherits(index, "optimized_methylation_index")) {
    stop("Object must be an optimized_methylation_index")
  }
  saveRDS(index, file = file, compress = TRUE)
  message(sprintf("Index saved to %s", file))
  invisible(TRUE)
}

#' Load optimized index from disk
#' 
#' @param file Path to the saved index file
#' @return Optimized index object
#' @export
load_index <- function(file) {
  if (!file.exists(file)) {
    stop(sprintf("File not found: %s", file))
  }
  index <- readRDS(file)
  if (!inherits(index, "optimized_methylation_index")) {
    stop("File does not contain an optimized_methylation_index")
  }
  return(index)
}

#' Print method for optimized index
#' 
#' @param x Optimized index object
#' @param ... Additional arguments (ignored)
#' @export
print.optimized_methylation_index <- function(x, ...) {
  cat("Optimized Methylation Index\n")
  cat("==========================\n")
  cat(sprintf("Type: %s\n", x$metadata$type))
  cat(sprintf("Cells: %d\n", length(x$cell_index)))
  cat(sprintf("Window size: %d bp\n", x$metadata$window_size))
  cat(sprintf("Created: %s\n", format(x$metadata$creation_time)))
  
  # Show chromosomes
  all_chr <- unique(unlist(lapply(x$cell_index, function(cell) {
    if (!is.null(cell$summary)) unique(cell$summary$chr) else NULL
  })))
  cat(sprintf("Chromosomes: %s\n", paste(sort(all_chr), collapse = ", ")))
  
  # Show read strategy summary
  if (length(x$read_strategies) > 0) {
    all_strategies <- data.table::rbindlist(x$read_strategies)
    strategy_summary <- all_strategies[, .N, by = strategy]
    cat("\nRead strategies:\n")
    print(strategy_summary)
  }
  
  invisible(x)
}