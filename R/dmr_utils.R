# dmr_utils.R - Utility functions for DMR analysis

#' Annotate DMRs with overlapping genes
#' 
#' @param dmr_collapsed Collapsed DMRs from collapse_dmr
#' @param gtf_file Path to GTF annotation file
#' @param feature_type Character, feature type to overlap (default "gene")
#' @return data.table with gene annotations added
annotate_dmr <- function(dmr_collapsed, gtf_file, feature_type = "gene") {
  
  # Read GTF
  gtf <- rtracklayer::import(gtf_file)
  gtf <- as.data.table(gtf)
  
  # Filter for genes
  genes <- gtf[type == feature_type]
  
  # Create ranges
  dmr_gr <- GRanges(
    seqnames = dmr_collapsed$chr,
    ranges = IRanges(start = dmr_collapsed$dmr_start, 
                     end = dmr_collapsed$dmr_end)
  )
  
  gene_gr <- GRanges(
    seqnames = genes$seqnames,
    ranges = IRanges(start = genes$start, end = genes$end),
    gene_name = genes$gene_name
  )
  
  # Find overlaps
  overlaps <- findOverlaps(dmr_gr, gene_gr)
  
  # Aggregate gene names by DMR
  gene_annotations <- data.table(
    dmr_idx = queryHits(overlaps),
    gene_name = gene_gr$gene_name[subjectHits(overlaps)]
  )
  
  gene_annotations <- gene_annotations[, .(
    gene_names = paste(unique(gene_name), collapse = ", ")
  ), by = dmr_idx]
  
  # Add to results
  dmr_collapsed[gene_annotations$dmr_idx, gene_names := gene_annotations$gene_names]
  
  return(dmr_collapsed)
}

#' Create DMR tracks for visualization
#' 
#' @param smoothed_windows Smoothed windows from calc_smoothed_windows
#' @param group_col Character, which group column to extract
#' @return data.table formatted for genome browser
create_dmr_track <- function(smoothed_windows, group_col) {
  
  track <- smoothed_windows$pct_matrix[, .(window, get(group_col))]
  setnames(track, c("window", "score"))
  
  # Parse window coordinates
  track[, c("chr", "start", "end") := tstrsplit(window, "_", fixed = TRUE, type.convert = TRUE)]
  track[, window := NULL]
  
  # Add track header
  track_name <- paste0("track name='", group_col, "_methylation' ",
                       "description='% methylation for ", group_col, "' ",
                       "visibility=full color=0,0,255 altColor=255,0,0 ",
                       "viewLimits=0:100 autoScale=off")
  
  list(header = track_name, data = track[, .(chr, start, end, score)])
}

#' Export DMRs to BED format
#' 
#' @param dmr_results Filtered or collapsed DMR results
#' @param output_file Path to output BED file
#' @param track_name Name for the track
#' @param separate_by Character, column to separate tracks by (e.g., "direction")
export_dmr_bed <- function(dmr_results, 
                           output_file, 
                           track_name = "DMRs",
                           separate_by = NULL) {
  
  if (!is.null(separate_by)) {
    # Export separate files by category
    categories <- unique(dmr_results[[separate_by]])
    
    for (cat in categories) {
      subset_dmr <- dmr_results[get(separate_by) == cat]
      
      bed_data <- subset_dmr[, .(
        chr = chr,
        start = ifelse("dmr_start" %in% names(subset_dmr), dmr_start, start),
        end = ifelse("dmr_end" %in% names(subset_dmr), dmr_end, end),
        name = paste0(test, "_", cat),
        score = round(abs(dmr_logFC) * 100),  # Scale for visualization
        strand = ifelse(dmr_logFC > 0, "+", "-")
      )]
      
      # Write with track line
      file_name <- sub("\\.bed$", paste0("_", cat, ".bed"), output_file)
      track_color <- ifelse(cat == "hyper", "255,0,0", "0,0,255")
      
      cat(paste0("track name='", track_name, "_", cat, "' ",
                 "description='", track_name, " ", cat, "methylated' ",
                 "color=", track_color, "\n"),
          file = file_name)
      
      fwrite(bed_data, file = file_name, sep = "\t", 
             col.names = FALSE, append = TRUE)
    }
  } else {
    # Export single file
    bed_data <- dmr_results[, .(
      chr = chr,
      start = ifelse("dmr_start" %in% names(dmr_results), dmr_start, start),
      end = ifelse("dmr_end" %in% names(dmr_results), dmr_end, end),
      name = test,
      score = round(abs(dmr_logFC) * 100),
      strand = ifelse(dmr_logFC > 0, "+", "-")
    )]
    
    cat(paste0("track name='", track_name, "' description='", track_name, "'\n"),
        file = output_file)
    
    fwrite(bed_data, file = output_file, sep = "\t", 
           col.names = FALSE, append = TRUE)
  }
}

#' Generate DMR summary statistics
#' 
#' @param dmr_results Any DMR result table
#' @return List with summary statistics
summarize_dmr <- function(dmr_results) {
  
  # Basic counts
  n_total <- nrow(dmr_results)
  
  # By test if available
  if ("test" %in% names(dmr_results)) {
    by_test <- dmr_results[, .N, by = test]
  } else {
    by_test <- NULL
  }
  
  # By direction
  if ("direction" %in% names(dmr_results)) {
    by_direction <- dmr_results[, .N, by = direction]
    
    if ("test" %in% names(dmr_results)) {
      by_test_direction <- dmr_results[, .N, by = .(test, direction)]
    } else {
      by_test_direction <- NULL
    }
  } else {
    by_direction <- NULL
    by_test_direction <- NULL
  }
  
  # Length statistics if collapsed
  if ("dmr_length" %in% names(dmr_results)) {
    length_stats <- dmr_results[, .(
      mean_length = mean(dmr_length),
      median_length = median(dmr_length),
      min_length = min(dmr_length),
      max_length = max(dmr_length),
      total_coverage = sum(dmr_length)
    )]
  } else {
    length_stats <- NULL
  }
  
  # Effect size statistics
  if ("dmr_logFC" %in% names(dmr_results)) {
    effect_stats <- dmr_results[, .(
      mean_abs_logFC = mean(abs(dmr_logFC)),
      median_abs_logFC = median(abs(dmr_logFC)),
      max_abs_logFC = max(abs(dmr_logFC))
    )]
  } else if ("logFC" %in% names(dmr_results)) {
    effect_stats <- dmr_results[, .(
      mean_abs_logFC = mean(abs(logFC)),
      median_abs_logFC = median(abs(logFC)),
      max_abs_logFC = max(abs(logFC))
    )]
  } else {
    effect_stats <- NULL
  }
  
  # Compile results
  list(
    total_dmrs = n_total,
    by_test = by_test,
    by_direction = by_direction,
    by_test_direction = by_test_direction,
    length_statistics = length_stats,
    effect_statistics = effect_stats
  )
}

#' Validate H5 file structure
#' 
#' @param h5_path Path to H5 file
#' @param expected_structure List describing expected structure
#' @return Logical indicating if structure is valid
validate_h5_structure <- function(h5_path, 
                                  expected_structure = list(
                                    contexts = c("CG", "CH"),
                                    required_columns = c("chr", "pos", "c", "t")
                                  )) {
  
  # Check file exists
  if (!file.exists(h5_path)) {
    warning(sprintf("H5 file not found: %s", h5_path))
    return(FALSE)
  }
  
  # List file contents
  h5_contents <- h5ls(h5_path)
  
  # Check for methylation contexts
  contexts_found <- unique(h5_contents$name[h5_contents$group == "/"])
  missing_contexts <- setdiff(expected_structure$contexts, contexts_found)
  
  if (length(missing_contexts) > 0) {
    warning(sprintf("Missing contexts in H5 file: %s", 
                    paste(missing_contexts, collapse = ", ")))
  }
  
  # Check structure for each context
  for (context in intersect(expected_structure$contexts, contexts_found)) {
    # Get barcodes for this context
    barcodes <- unique(h5_contents$name[h5_contents$group == paste0("/", context)])
    
    if (length(barcodes) == 0) {
      warning(sprintf("No barcodes found for context: %s", context))
      next
    }
    
    # Check first barcode structure
    test_path <- paste0(context, "/", barcodes[1], "/1")
    
    tryCatch({
      test_data <- h5read(h5_path, test_path, index = list(1:10))
      
      # Check columns
      missing_cols <- setdiff(expected_structure$required_columns, names(test_data))
      if (length(missing_cols) > 0) {
        warning(sprintf("Missing columns in %s: %s", 
                        test_path, paste(missing_cols, collapse = ", ")))
        return(FALSE)
      }
    }, error = function(e) {
      warning(sprintf("Error reading %s: %s", test_path, e$message))
      return(FALSE)
    })
  }
  
  return(TRUE)
}

#' Configure analysis parameters
#' 
#' @param genome Character, genome build
#' @param context Character, methylation context
#' @param window_size Integer, base window size
#' @param smoothing_windows Integer, number of windows to smooth
#' @param min_coverage Integer, minimum coverage threshold
#' @return List of analysis parameters
configure_dmr_analysis <- function(genome = "hg38",
                                   context = "CG",
                                   window_size = 500,
                                   smoothing_windows = 3,
                                   min_coverage = 5) {
  
  # Validate inputs
  if (!genome %in% c("hg38", "hg19", "mm10", "mm39")) {
    stop("Unsupported genome. Use 'hg38', 'hg19', 'mm10', or 'mm39'")
  }
  
  if (!context %in% c("CG", "CH", "CHG", "CHH")) {
    stop("Unsupported methylation context")
  }
  
  if (window_size < 100 || window_size > 10000) {
    warning("Unusual window size. Typical values are 500-2000bp")
  }
  
  # Create configuration
  config <- list(
    genome = genome,
    context = context,
    window_size = window_size,
    smoothing_windows = smoothing_windows,
    smoothed_window_size = window_size * smoothing_windows,
    min_coverage = min_coverage,
    created = Sys.time()
  )
  
  class(config) <- c("dmr_config", "list")
  return(config)
}

#' Print method for DMR configuration
print.dmr_config <- function(x, ...) {
  cat("DMR Analysis Configuration\n")
  cat("--------------------------\n")
  cat("Genome:", x$genome, "\n")
  cat("Context:", x$context, "\n")
  cat("Window size:", x$window_size, "bp\n")
  cat("Smoothing:", x$smoothing_windows, "windows\n")
  cat("Effective window:", x$smoothed_window_size, "bp\n")
  cat("Min coverage:", x$min_coverage, "\n")
  cat("Created:", format(x$created, "%Y-%m-%d %H:%M:%S"), "\n")
}

# Run complete DMR analysis pipeline
run_dmr_pipeline <- function(h5_paths, 
                             metadata,
                             config,
                             output_dir = ".",
                             threads = 1) {
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Log start
  message(sprintf("Starting DMR analysis: %s", Sys.time()))
  message(sprintf("Cells: %d, Threads: %d", nrow(h5_paths), threads))
  
  # Validate inputs
  if (!all(metadata$cell_id %in% h5_paths$barcode)) {
    stop("Not all cells in metadata are present in h5_paths")
  }
  
  # Step 1: Index chromosomes
  message("Step 1: Indexing chromosomes...")
  chr_index <- index_chromosomes(
    h5_paths = h5_paths,
    type = config$context,
    threads = threads
  )
  
  # Save checkpoint
  saveRDS(chr_index, file.path(output_dir, "chr_index.rds"))
  
  # Step 2: Calculate smoothed windows
  message("Step 2: Calculating smoothed windows...")
  windows <- calc_smoothed_windows(
    h5_paths = h5_paths,
    chr_index = chr_index,
    metadata = metadata,
    type = config$context,
    step = config$window_size,
    smooth = config$smoothing_windows,
    genome = config$genome,
    threads = threads
  )
  
  # Save checkpoint
  saveRDS(windows, file.path(output_dir, "smoothed_windows.rds"))
  fwrite(windows$sum_matrix, file.path(output_dir, "count_matrix.tsv"), sep = "\t")
  fwrite(windows$pct_matrix, file.path(output_dir, "percent_matrix.tsv"), sep = "\t")
  
  # Step 3: Test for DMRs
  message("Step 3: Testing for DMRs...")
  dmr_results <- test_dmr(
    sum_matrix = windows$sum_matrix,
    min_total = config$min_coverage * 2,
    min_group = config$min_coverage
  )
  
  # Save checkpoint
  saveRDS(dmr_results, file.path(output_dir, "dmr_raw_results.rds"))
  
  # Step 4: Filter DMRs
  message("Step 4: Filtering DMRs...")
  dmr_filtered <- filter_dmr(
    dmr_matrix = dmr_results,
    method = "bonferroni",
    p_threshold = 0.01,
    log_threshold = 1.5,
    filter = TRUE
  )
  
  # Save results
  fwrite(dmr_filtered, file.path(output_dir, "dmr_filtered.tsv"), sep = "\t")
  
  # Step 5: Collapse DMRs (only if we have results)
  if (nrow(dmr_filtered) > 0) {
    message("Step 5: Collapsing adjacent DMRs...")
    dmr_collapsed <- collapse_dmr(
      dmr_filtered = dmr_filtered,
      max_gap = config$window_size * 4,
      min_length = config$window_size * 4
    )
    
    # Save results
    fwrite(dmr_collapsed, file.path(output_dir, "dmr_collapsed.tsv"), sep = "\t")
    
    # Export BED files
    export_dmr_bed(
      dmr_results = dmr_collapsed,
      output_file = file.path(output_dir, "dmrs.bed"),
      track_name = paste(config$context, "DMRs"),
      separate_by = "direction"
    )
  } else {
    message("No significant DMRs found after filtering")
    dmr_collapsed <- NULL
  }
  
  # Step 6: Generate summary
  message("Step 6: Generating summary...")
  summary_stats <- summarize_dmr(dmr_filtered)
  
  # Save summary
  saveRDS(summary_stats, file.path(output_dir, "dmr_summary.rds"))
  
  message(sprintf("DMR analysis complete: %s", Sys.time()))
  message(sprintf("Results saved to: %s", output_dir))
  
  # Return all results
  list(
    config = config,
    chr_index = chr_index,
    windows = windows,
    dmr_raw = dmr_results,
    dmr_filtered = dmr_filtered,
    dmr_collapsed = dmr_collapsed,
    summary = summary_stats
  )
}

#' Diagnose DMR analysis issues
#' 
#' @param h5_path Path to H5 file
#' @param barcode Barcode to test
#' @param context Methylation context
#' @return List with diagnostic information
diagnose_dmr_data <- function(h5_path, barcode, context = "CG") {
  
  cat("=== DMR Data Diagnostics ===\n")
  
  # Check file exists
  if (!file.exists(h5_path)) {
    stop("H5 file not found: ", h5_path)
  }
  
  # Check structure
  h5_contents <- h5ls(h5_path)
  cat("\nH5 file structure:\n")
  print(table(h5_contents$group))
  
  # Try to read data
  data_path <- paste0(context, "/", barcode, "/1")
  
  tryCatch({
    # Read first 100 rows
    test_data <- h5read(h5_path, data_path, index = list(1:100, NULL))
    
    cat("\nData preview:\n")
    print(head(test_data))
    
    cat("\nData dimensions:", dim(test_data), "\n")
    cat("Column names:", names(test_data), "\n")
    cat("Data class:", class(test_data), "\n")
    
    # Check data types
    cat("\nColumn types:\n")
    print(sapply(test_data, class))
    
    # Check for required columns
    required <- c("chr", "pos", "c", "t")
    missing <- setdiff(required, names(test_data))
    if (length(missing) > 0) {
      cat("\nWARNING: Missing required columns:", paste(missing, collapse = ", "), "\n")
    }
    
    # Check chromosomes
    if ("chr" %in% names(test_data)) {
      cat("\nChromosomes found:", paste(unique(test_data$chr), collapse = ", "), "\n")
    }
    
    # Check coverage
    if (all(c("c", "t") %in% names(test_data))) {
      coverage <- test_data$c + test_data$t
      cat("\nCoverage statistics:\n")
      cat("  Mean:", mean(coverage), "\n")
      cat("  Median:", median(coverage), "\n")
      cat("  Range:", paste(range(coverage), collapse = "-"), "\n")
    }
    
    return(list(
      success = TRUE,
      data = test_data,
      structure = h5_contents
    ))
    
  }, error = function(e) {
    cat("\nERROR reading data:", e$message, "\n")
    return(list(
      success = FALSE,
      error = e$message,
      structure = h5_contents
    ))
  })
}