#!/usr/bin/env Rscript
# verify_toy_data.R - Verify toy data and visualize DMRs

library(rhdf5)
library(data.table)
library(ggplot2)

# Read metadata
metadata <- fread("cell_metadata.tsv")
known_dmrs <- fread("known_dmrs.bed", 
                    col.names = c("chr", "start", "end", "name", "score", "strand"))

# Read H5 file structure
h5_file <- "methylation_data.h5"
h5_structure <- h5ls(h5_file)

cat("=== H5 File Summary ===\n")
cat(sprintf("Contexts: %s\n", paste(unique(h5_structure$group[h5_structure$group != "/"]), collapse = ", ")))
cat(sprintf("Total datasets: %d\n", sum(h5_structure$otype == "H5I_DATASET")))

# Function to calculate methylation in windows
calculate_window_methylation <- function(h5_file, barcode, context, chr, start, end, window_size = 500) {
  # Read data for the chromosome
  data <- h5read(h5_file, paste0(context, "/", barcode, "/1"))
  data <- as.data.table(data)
  
  # Filter to region of interest
  region_data <- data[chr == chr & pos >= start & pos <= end]
  
  if (nrow(region_data) == 0) {
    return(NULL)
  }
  
  # Calculate windows
  region_data[, window := ((pos - start) %/% window_size) * window_size + start]
  
  # Aggregate by window
  window_meth <- region_data[, .(
    n_sites = .N,
    total_c = sum(c),
    total_t = sum(t),
    meth_pct = 100 * sum(c) / sum(c + t)
  ), by = window]
  
  window_meth[, barcode := barcode]
  return(window_meth)
}

# Examine methylation patterns around DMRs
cat("\n=== Methylation Patterns at Known DMRs ===\n")

# Pick first DMR to examine
dmr <- known_dmrs[1]
region_start <- dmr$start - 5000
region_end <- dmr$end + 5000

cat(sprintf("\nExamining DMR: %s (%s:%d-%d)\n", 
            dmr$name, dmr$chr, dmr$start, dmr$end))

# Calculate methylation for all cells in this region
all_window_meth <- rbindlist(lapply(metadata$cell_id, function(barcode) {
  calculate_window_methylation(h5_file, barcode, "CG", 
                               dmr$chr, region_start, region_end)
}))

# Add cell type information
all_window_meth <- merge(all_window_meth, 
                         metadata[, .(cell_id, cluster_id)], 
                         by.x = "barcode", by.y = "cell_id")

# Average by cell type
avg_window_meth <- all_window_meth[, .(
  mean_meth = mean(meth_pct),
  se_meth = sd(meth_pct) / sqrt(.N),
  n_cells = .N
), by = .(window, cluster_id)]

# Plot
p1 <- ggplot(avg_window_meth, aes(x = window, y = mean_meth, color = cluster_id)) +
  geom_ribbon(aes(ymin = mean_meth - se_meth, ymax = mean_meth + se_meth, 
                  fill = cluster_id), alpha = 0.3) +
  geom_line(size = 1) +
  geom_rect(xmin = dmr$start, xmax = dmr$end, ymin = -Inf, ymax = Inf,
            fill = "gray", alpha = 0.2, color = NA) +
  labs(title = sprintf("Methylation at %s", dmr$name),
       x = sprintf("Position on %s", dmr$chr),
       y = "Mean methylation %",
       color = "Cell type",
       fill = "Cell type") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p1)

# Compare all DMRs
cat("\n=== Methylation Differences at All DMRs ===\n")

# Calculate methylation at all DMRs
dmr_meth <- rbindlist(lapply(1:nrow(known_dmrs), function(i) {
  dmr <- known_dmrs[i]
  
  # Get methylation for each cell
  cell_meth <- rbindlist(lapply(metadata$cell_id, function(barcode) {
    data <- h5read(h5_file, paste0("CG/", barcode, "/1"))
    data <- as.data.table(data)
    
    # Filter to DMR
    dmr_data <- data[chr == dmr$chr & pos >= dmr$start & pos <= dmr$end]
    
    if (nrow(dmr_data) == 0) {
      return(NULL)
    }
    
    data.table(
      dmr_name = dmr$name,
      barcode = barcode,
      n_sites = nrow(dmr_data),
      meth_pct = 100 * sum(dmr_data$c) / sum(dmr_data$c + dmr_data$t)
    )
  }))
  
  return(cell_meth)
}))

# Add cell type
dmr_meth <- merge(dmr_meth, metadata[, .(cell_id, cluster_id)], 
                  by.x = "barcode", by.y = "cell_id")

# Calculate differences
dmr_diff <- dmr_meth[, .(
  t_cell_meth = mean(meth_pct[cluster_id == "T_cell"]),
  b_cell_meth = mean(meth_pct[cluster_id == "B_cell"])
), by = dmr_name]

dmr_diff[, diff := t_cell_meth - b_cell_meth]
dmr_diff[, direction := ifelse(diff > 0, "T_cell_hyper", "T_cell_hypo")]

# Merge with known DMR info
dmr_diff <- merge(dmr_diff, known_dmrs[, .(name, expected_direction = strand)], 
                  by.x = "dmr_name", by.y = "name")
dmr_diff[, expected_direction := ifelse(expected_direction == "+", "T_cell_hyper", "T_cell_hypo")]

print(dmr_diff)

# Plot methylation differences
p2 <- ggplot(dmr_meth, aes(x = dmr_name, y = meth_pct, fill = cluster_id)) +
  geom_boxplot(alpha = 0.7) +
  labs(title = "Methylation at Known DMRs",
       x = "DMR",
       y = "Methylation %",
       fill = "Cell type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p2)

# Verify DMR detection works
cat("\n=== Verification Summary ===\n")
cat(sprintf("Correctly oriented DMRs: %d/%d\n", 
            sum(dmr_diff$direction == dmr_diff$expected_direction),
            nrow(dmr_diff)))

# Show methylation distribution
p3 <- ggplot(dmr_meth, aes(x = meth_pct, fill = cluster_id)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~dmr_name, scales = "free_y") +
  labs(title = "Methylation Distribution at DMRs",
       x = "Methylation %",
       y = "Density",
       fill = "Cell type") +
  theme_minimal()

print(p3)

# Global methylation comparison
cat("\n=== Global Methylation Patterns ===\n")

# Sample random windows across genome
set.seed(123)
random_windows <- rbindlist(lapply(names(table(known_dmrs$chr)), function(chr) {
  data.table(
    chr = chr,
    start = sample(1:50000000, 20) * 1000,  # 20 random 1kb windows
    end = NA
  )
}))
random_windows[, end := start + 1000]

# Calculate methylation at random windows
random_meth <- rbindlist(lapply(1:nrow(random_windows), function(i) {
  win <- random_windows[i]
  
  cell_meth <- rbindlist(lapply(metadata$cell_id, function(barcode) {
    data <- h5read(h5_file, paste0("CG/", barcode, "/1"))
    data <- as.data.table(data)
    
    win_data <- data[chr == win$chr & pos >= win$start & pos <= win$end]
    
    if (nrow(win_data) < 5) {  # Require at least 5 sites
      return(NULL)
    }
    
    data.table(
      window = paste0(win$chr, ":", win$start),
      barcode = barcode,
      meth_pct = 100 * sum(win_data$c) / sum(win_data$c + win_data$t)
    )
  }))
}))

# Add cell type
random_meth <- merge(random_meth, metadata[, .(cell_id, cluster_id)], 
                     by.x = "barcode", by.y = "cell_id")

# Compare to DMR methylation
plot_data <- rbind(
  dmr_meth[, .(meth_pct, cluster_id, type = "DMR")],
  random_meth[, .(meth_pct, cluster_id, type = "Random")]
)

p4 <- ggplot(plot_data, aes(x = cluster_id, y = meth_pct, fill = type)) +
  geom_boxplot(alpha = 0.7) +
  labs(title = "DMR vs Random Region Methylation",
       x = "Cell type",
       y = "Methylation %",
       fill = "Region type") +
  theme_minimal()

print(p4)

# Save plots
ggsave("dmr_methylation_profile.pdf", p1, width = 8, height = 6)
ggsave("dmr_methylation_boxplot.pdf", p2, width = 10, height = 6)
ggsave("dmr_methylation_density.pdf", p3, width = 12, height = 8)
ggsave("dmr_vs_random.pdf", p4, width = 8, height = 6)

cat("\nPlots saved as PDF files\n")