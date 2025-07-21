# utils.R - Utility functions for minithyst package

#' Diagnose index and metadata compatibility
#' 
#' Checks for common issues with cell ID matching between index, metadata,
#' and h5_paths. Useful for debugging when calc_smoothed_windows fails.
#' 
#' @param chr_index Optimized index object
#' @param metadata data.table with cell metadata
#' @param h5_paths data.table with h5 file paths
#' @return Invisible list with diagnostic information
#' @export
#' @examples
#' \dontrun{
#' # Run diagnostic if calc_smoothed_windows fails
#' diagnose_index_metadata(index, metadata, h5_paths)
#' }
diagnose_index_metadata <- function(chr_index, metadata, h5_paths) {
  message("=== INDEX DIAGNOSTIC ===")
  
  # Check index type and contents
  if (inherits(chr_index, "optimized_methylation_index")) {
    message(sprintf("Index type: Optimized"))
    message(sprintf("Cells in index: %d", length(chr_index$cell_index)))
    
    if (length(chr_index$cell_index) > 0) {
      message(sprintf("Sample cells: %s", 
                      paste(head(names(chr_index$cell_index), 3), collapse = ", ")))
      
      # Check structure of first cell
      first_cell <- names(chr_index$cell_index)[1]
      cell_info <- chr_index$cell_index[[first_cell]]
      message(sprintf("\nFirst cell (%s) structure:", first_cell))
      message(sprintf("  - summary: %s (%d rows)", 
                      !is.null(cell_info$summary), 
                      if (!is.null(cell_info$summary)) nrow(cell_info$summary) else 0))
      message(sprintf("  - windows: %s (%d rows)", 
                      !is.null(cell_info$windows),
                      if (!is.null(cell_info$windows)) nrow(cell_info$windows) else 0))
      message(sprintf("  - path: %s", 
                      if (!is.null(cell_info$path)) cell_info$path else "NULL"))
    }
  } else {
    message("Index type: Unknown/Invalid")
  }
  
  message("\n=== METADATA DIAGNOSTIC ===")
  message(sprintf("Rows: %d", nrow(metadata)))
  message(sprintf("Columns: %s", paste(names(metadata), collapse = ", ")))
  if ("cell_id" %in% names(metadata)) {
    message(sprintf("Sample cell_ids: %s", 
                    paste(head(metadata$cell_id, 3), collapse = ", ")))
  }
  
  message("\n=== H5_PATHS DIAGNOSTIC ===")
  message(sprintf("Rows: %d", nrow(h5_paths)))
  message(sprintf("Sample barcodes: %s", 
                  paste(head(h5_paths$barcode, 3), collapse = ", ")))
  
  # Check alignment
  message("\n=== ALIGNMENT CHECK ===")
  
  index_cells <- names(chr_index$cell_index)
  metadata_cells <- if ("cell_id" %in% names(metadata)) metadata$cell_id else character(0)
  h5_cells <- h5_paths$barcode
  
  # Check overlaps
  index_metadata_overlap <- length(intersect(index_cells, metadata_cells))
  metadata_h5_overlap <- length(intersect(metadata_cells, h5_cells))
  index_h5_overlap <- length(intersect(index_cells, h5_cells))
  
  message(sprintf("Cells in both index and metadata: %d / %d", 
                  index_metadata_overlap, length(metadata_cells)))
  message(sprintf("Cells in both metadata and h5_paths: %d / %d", 
                  metadata_h5_overlap, length(metadata_cells)))
  message(sprintf("Cells in both index and h5_paths: %d / %d", 
                  index_h5_overlap, length(h5_cells)))
  
  if (index_metadata_overlap == 0) {
    warning("No cells found in both index and metadata! This will cause calc_smoothed_windows to fail.")
  }
  
  invisible(list(
    index_cells = index_cells,
    metadata_cells = metadata_cells,
    h5_cells = h5_cells,
    overlaps = list(
      index_metadata = index_metadata_overlap,
      metadata_h5 = metadata_h5_overlap,
      index_h5 = index_h5_overlap
    )
  ))
}

#' Validate inputs for calc_smoothed_windows
#' 
#' Checks that all required inputs are properly formatted before processing.
#' 
#' @param h5_paths data.table with h5 file paths
#' @param chr_index Optimized index object  
#' @param metadata data.table with cell metadata
#' @param group_by Character, grouping column name
#' @return TRUE if valid, otherwise stops with error
#' @export
validate_inputs <- function(h5_paths, chr_index, metadata, group_by) {
  # Check h5_paths
  if (!data.table::is.data.table(h5_paths)) {
    h5_paths <- data.table::as.data.table(h5_paths)
  }
  if (!all(c("barcode", "path") %in% names(h5_paths))) {
    stop("h5_paths must contain 'barcode' and 'path' columns")
  }
  
  # Check index
  if (!inherits(chr_index, "optimized_methylation_index")) {
    stop("chr_index must be created with create_optimized_index()")
  }
  
  # Check metadata
  if (!data.table::is.data.table(metadata)) {
    metadata <- data.table::as.data.table(metadata)
  }
  if (!group_by %in% names(metadata)) {
    stop(sprintf("Column '%s' not found in metadata", group_by))
  }
  if (!"cell_id" %in% names(metadata)) {
    if (!is.null(rownames(metadata))) {
      warning("Using rownames as cell_id")
      metadata[, cell_id := rownames(metadata)]
    } else {
      stop("metadata must have 'cell_id' column or row names")
    }
  }
  
  # Check alignment
  index_cells <- names(chr_index$cell_index)
  metadata_cells <- metadata$cell_id
  overlap <- length(intersect(index_cells, metadata_cells))
  
  if (overlap == 0) {
    stop("No cells found in both index and metadata. Run diagnose_index_metadata() for details.")
  }
  
  message(sprintf("Validated: %d cells found in both index and metadata", overlap))
  return(TRUE)
}

#' Create example data for testing
#' 
#' Generates a small example dataset for testing minithyst functions.
#' 
#' @param n_cells Integer, number of cells to generate
#' @param n_sites_per_cell Integer, approximate sites per cell
#' @return List with h5_paths, metadata, and temporary h5 file
#' @export
#' @examples
#' \dontrun{
#' # Create test data
#' test_data <- create_test_data(n_cells = 10)
#' 
#' # Use in pipeline
#' index <- create_optimized_index(test_data$h5_paths)
#' result <- calc_smoothed_windows(
#'   h5_paths = test_data$h5_paths,
#'   chr_index = index,
#'   metadata = test_data$metadata
#' )
#' 
#' # Clean up
#' unlink(test_data$h5_file)
#' }
create_test_data <- function(n_cells = 10, n_sites_per_cell = 1000) {
  require(rhdf5)
  
  # Generate cell barcodes
  barcodes <- paste0("cell", seq_len(n_cells))
  
  # Create temporary h5 file
  h5_file <- tempfile(fileext = ".h5")
  rhdf5::h5createFile(h5_file)
  
  # Create groups
  rhdf5::h5createGroup(h5_file, "CG")
  
  # Generate data for each cell
  for (i in seq_len(n_cells)) {
    barcode <- barcodes[i]
    
    # Create cell group
    rhdf5::h5createGroup(h5_file, paste0("CG/", barcode))
    
    # Generate random methylation data
    n_sites <- round(n_sites_per_cell * runif(1, 0.8, 1.2))
    chr_sites <- sample(c("chr1", "chr2", "chr3"), n_sites, replace = TRUE)
    
    # Generate positions for each chromosome
    data_list <- lapply(unique(chr_sites), function(chr) {
      chr_n <- sum(chr_sites == chr)
      positions <- sort(sample(1:1000000, chr_n))
      
      data.frame(
        chr = chr,
        pos = positions,
        c = rbinom(chr_n, 10, 0.7),
        t = rbinom(chr_n, 10, 0.3),
        pct = round(runif(chr_n, 0, 100))
      )
    })
    
    cell_data <- do.call(rbind, data_list)
    cell_data <- cell_data[order(cell_data$chr, cell_data$pos), ]
    
    # Write to h5
    rhdf5::h5write(cell_data, h5_file, paste0("CG/", barcode, "/1"))
  }
  
  rhdf5::h5closeAll()
  
  # Create metadata
  metadata <- data.table::data.table(
    cell_id = barcodes,
    cluster_id = sample(1:3, n_cells, replace = TRUE),
    coverage = round(n_sites_per_cell * runif(n_cells, 0.8, 1.2))
  )
  
  # Create h5_paths
  h5_paths <- data.table::data.table(
    barcode = barcodes,
    path = h5_file
  )
  
  return(list(
    h5_paths = h5_paths,
    metadata = metadata,
    h5_file = h5_file
  ))
}

#' Summarize smoothed window results
#' 
#' Provides summary statistics for the output of calc_smoothed_windows.
#' 
#' @param result Output from calc_smoothed_windows
#' @param type Character, either "sum" or "pct" to specify which matrix
#' @return data.table with summary statistics
#' @export
summarize_windows <- function(result, type = "sum") {
  if (is.list(result)) {
    if (type == "sum" && !is.null(result$sum_matrix)) {
      mat <- result$sum_matrix
    } else if (type == "pct" && !is.null(result$pct_matrix)) {
      mat <- result$pct_matrix
    } else {
      stop("Requested matrix type not found in result")
    }
  } else {
    mat <- result
  }
  
  # Get data columns (exclude chr, start, end)
  data_cols <- setdiff(names(mat), c("chr", "start", "end"))
  
  # Calculate coverage per chromosome
  chr_summary <- mat[, .(
    n_windows = .N,
    total_sites = sum(rowSums(.SD, na.rm = TRUE))
  ), by = chr, .SDcols = data_cols]
  
  # Overall statistics
  overall_summary <- list(
    total_windows = nrow(mat),
    total_chromosomes = length(unique(mat$chr)),
    window_size = mat[2, start] - mat[1, start],
    groups = length(data_cols) / ifelse(type == "sum", 2, 1)
  )
  
  message("=== WINDOW SUMMARY ===")
  message(sprintf("Total windows: %d", overall_summary$total_windows))
  message(sprintf("Chromosomes: %d", overall_summary$total_chromosomes))
  message(sprintf("Window size: %d bp", overall_summary$window_size))
  message(sprintf("Groups: %d", overall_summary$groups))
  message("\nPer-chromosome coverage:")
  print(chr_summary)
  
  invisible(list(
    overall = overall_summary,
    by_chromosome = chr_summary
  ))
}