#!/usr/bin/env Rscript
# debug_h5_structure.R - Debug the H5 file structure

library(rhdf5)
library(data.table)

# Check if files exist
if (!file.exists("./Data/methylation_data.h5")) {
  stop("methylation_data.h5 not found. Run generate_toy_data.R first.")
}

h5_file <- "./Data/methylation_data.h5"

cat("=== H5 File Structure Debug ===\n\n")

# List all contents
cat("1. Full H5 contents:\n")
h5_contents <- rhdf5::h5ls(h5_file)
print(h5_contents)

# Get a sample barcode
barcodes <- unique(h5_contents$name[h5_contents$group == "/CG"])
test_barcode <- barcodes[1]

cat("\n2. Testing with barcode:", test_barcode, "\n")

# Try different ways to read the data
data_path <- paste0("CG/", test_barcode, "/1")

cat("\n3. Dataset info:\n")
h5_file_handle <- rhdf5::H5Fopen(h5_file)
h5_dataset <- rhdf5::H5Dopen(h5_file_handle, data_path)
h5_space <- rhdf5::H5Dget_space(h5_dataset)
dims <- rhdf5::H5Sget_simple_extent_dims(h5_space)
cat("  Dimensions:", dims$size, "\n")
cat("  Rank:", dims$rank, "\n")

# Get data type - but don't try to close it separately
h5_type <- rhdf5::H5Dget_type(h5_dataset)
type_class <- rhdf5::H5Tget_class(h5_type)
cat("  Data type class:", type_class, "\n")

# Close all handles at once
rhdf5::h5closeAll()

cat("\n4. Testing different read methods:\n")

# Method 1: Read without index
cat("\n  Method 1 - Read full dataset:\n")
tryCatch({
  data1 <- rhdf5::h5read(h5_file, data_path)
  cat("    Success! Shape:", dim(data1), "\n")
  cat("    Class:", class(data1), "\n")
  cat("    Column names:", names(data1), "\n")
  cat("    First few rows:\n")
  print(head(data1))
}, error = function(e) {
  cat("    Error:", e$message, "\n")
})

# Method 2: Read with NULL index
cat("\n  Method 2 - Read with NULL index:\n")
tryCatch({
  data2 <- rhdf5::h5read(h5_file, data_path, index = NULL)
  cat("    Success! Shape:", dim(data2), "\n")
}, error = function(e) {
  cat("    Error:", e$message, "\n")
})

# Method 3: Read specific rows (1D compound dataset)
cat("\n  Method 3 - Read rows 1-10 (1D index):\n")
tryCatch({
  data3 <- rhdf5::h5read(h5_file, data_path, index = list(1:10))
  cat("    Success! Shape:", dim(data3), "\n")
  print(head(data3))
}, error = function(e) {
  cat("    Error:", e$message, "\n")
})

# Method 4: Read with start and count
cat("\n  Method 4 - Read using start/count:\n")
tryCatch({
  data4 <- rhdf5::h5read(h5_file, data_path, start = 1, count = 10)
  cat("    Success! Shape:", dim(data4), "\n")
  print(head(data4))
}, error = function(e) {
  cat("    Error:", e$message, "\n")
})

cat("\n5. Understanding the data structure:\n")
# Read a small sample to understand structure
sample_data <- rhdf5::h5read(h5_file, data_path, index = list(1:100))
cat("  Sample data structure:\n")
str(sample_data)

# Check column info
if (is.data.frame(sample_data)) {
  cat("\n  Column information:\n")
  for (i in seq_along(names(sample_data))) {
    cat(sprintf("    Column %d: %s (type: %s)\n", 
                i, names(sample_data)[i], class(sample_data[[i]])[1]))
  }
}

cat("\n6. Testing chromosome filtering:\n")
if ("chr" %in% names(sample_data)) {
  chr_counts <- table(sample_data$chr)
  cat("  Chromosome counts in first 100 rows:\n")
  print(chr_counts)
}

# Test reading larger chunk
cat("\n7. Performance test - reading larger chunks:\n")
system.time({
  large_data <- rhdf5::h5read(h5_file, data_path, index = list(1:10000))
})
cat("  Read 10,000 rows successfully\n")

cat("\n=== Summary ===\n")
cat("- Data is stored as 1D COMPOUND datasets\n")
cat("- Each compound element has fields: chr, pos, c, t\n")
cat("- Indexing: use list(rows) for 1D datasets\n")
cat("- Full read is most reliable for this structure\n")

# Clean up
rhdf5::h5closeAll()