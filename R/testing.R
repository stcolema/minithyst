# dmr_analysis.R - Simplified DMR Analysis from Single-Cell Methylation Data
# Focused reimplementation of Amethyst DMR functionality

#' @import data.table
#' @import GenomicRanges
#' @importFrom rhdf5 h5read h5ls H5Fopen H5Dopen H5Dget_space H5Sget_simple_extent_dims H5Dget_type H5Tget_class h5closeAll
#' @importFrom IRanges IRanges
#' @importFrom stats phyper dhyper p.adjust
#' @importFrom future plan multisession sequential
#' @importFrom furrr future_map future_pmap
#' @importFrom dplyr left_join
NULL

#' Fast Fisher's exact test implementation
#' @keywords internal
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

#' Filter counts matrix for DMR testing
#' @keywords internal
filter_count_matrix <- function(sum_matrix, min_total) {
  # Copy to avoid modifying input
  counts <- data.table::copy(sum_matrix)
  
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
  return(counts) 
}

#' Perform tests of each group against all others
#' @keywords internal
compare_all_vs_all <- function(counts, min_group){
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
  return(counts)
}

#' @title Test Differentially Methlated Regions
#' @description Test for differentially methylated regions
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
  counts <- filter_count_matrix(sum_matrix, min_total)
  if (is.null(comparisons)) {
    counts <- compare_all_vs_all(counts, min_group)
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

#' @title Filter DMRs
#' @description Filter DMR results based on multiple testing correction and 
#' log-fold change threshold
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
  results <- results[!is.na(pval) & !is.na(logFC)]

  # Apply multiple testing correction
  results[, padj := stats::p.adjust(pval, method = method)]

  # Add direction
  results[, direction := ifelse(logFC < 0, "hypo", "hyper")]

  # Filter if requested
  if (filter) {
    results <- results[(padj < p_threshold) & (abs(logFC) > log_threshold)]
  }

  return(results)
}

#' @title Collapse DMRs
#' 
#' @description Collapse adjacent DMRs into a single entity
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
  
  if (!all(c("chr", "start", "end") %in% names(results))) {
    stop("DMR results must have chr, start, end columns")
  }
  
  data.table::setorder(results, test, direction, chr, start)
  
  # Add group ID for merging overlapping intervals
  results[, `:=`(
    expanded_start = start - max_gap,
    expanded_end = end + max_gap
  )]
  
  results[, group_id := {
    new_group <- c(TRUE, expanded_start[-1] > cummax(expanded_end[-length(expanded_end)]))
    cumsum(new_group)
  }, by = .(test, direction, chr)]
  
  # Collapse by group
  collapsed <- results[, .(
    dmr_start = max(0, min(expanded_start) + max_gap),
    dmr_end = max(expanded_end) - max_gap,
    dmr_padj = mean(padj),
    dmr_logFC = mean(logFC),
    n_windows = .N,
    dmr_length = max(expanded_end) - min(expanded_start) + 1 - (2 * max_gap)
  ), by = .(test, direction, chr, group_id)]
  
  # Filter and clean
  collapsed <- collapsed[dmr_length >= min_length][, group_id := NULL]
  
  return(collapsed)
}
