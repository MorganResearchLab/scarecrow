#!/usr/bin/env Rscript
library(Matrix)
library(data.table)

args <- commandArgs(trailingOnly=TRUE)
print(args)
counts.file <- args[1]


convert_to_mtx_transposed <- function(input_gz, output_mtx, output_rows, output_cols, chunk_size = 10000) {
  # Open connections
  con_in <- gzfile(input_gz, "r")
  
  # Read header (column names - will become row names in transposed matrix)
  header <- readLines(con_in, n = 1)
  new_row_names <- unlist(strsplit(header, "\t"))[-1]  # These will be rows in output
  
  # First pass: count genes (will become columns in transposed matrix) and non-zero entries
  message("Counting genes and non-zero entries...")
  gene_names <- character()
  total_nonzero <- 0
  while(length(line <- readLines(con_in, n = 1)) > 0) {
    gene_names <- c(gene_names, strsplit(line, "\t")[[1]][1])
    vals <- as.numeric(strsplit(line, "\t")[[1]][-1])
    total_nonzero <- total_nonzero + sum(vals != 0)
  }
  n_genes <- length(gene_names)
  n_cells <- length(new_row_names)
  
  # Write header files
  writeLines(gene_names, output_rows)  # These are now columns in transposed matrix
  writeLines(new_row_names, output_cols)  # These are now rows in transposed matrix
  
  # Reopen input file
  close(con_in)
  con_in <- gzfile(input_gz, "r")
  readLines(con_in, n = 1)  # Skip header
  
  # Open MTX file and write header
  con_mtx <- file(output_mtx, "w")
  writeLines("%%MatrixMarket matrix coordinate real general", con_mtx)
  writeLines(sprintf("%d %d %d", n_cells, n_genes, total_nonzero), con_mtx)
  
  # Process data in chunks
  message("Processing data and writing transposed matrix...")
  pb <- txtProgressBar(min = 0, max = n_genes, style = 3)
  chunk <- list()
  gene_idx <- 0
  
  while(length(line <- readLines(con_in, n = 1)) > 0) {
    gene_idx <- gene_idx + 1
    parts <- strsplit(line, "\t")[[1]]
    vals <- as.numeric(parts[-1])
    
    # Record non-zero entries (cell_id, gene_id, count)
    nonzero <- which(vals != 0)
    if(length(nonzero) > 0) {
      chunk[[length(chunk) + 1]] <- data.frame(
        row = nonzero,
        col = gene_idx,
        val = vals[nonzero]
      )
    }
    
    # Write chunk when it reaches chunk_size
    if(length(chunk) >= chunk_size || gene_idx == n_genes) {
      chunk_df <- do.call(rbind, chunk)
      write.table(chunk_df, con_mtx, 
                  col.names = FALSE, row.names = FALSE,
                  sep = " ", quote = FALSE)
      chunk <- list()
    }
    
    setTxtProgressBar(pb, gene_idx)
  }
  
  close(pb)
  close(con_in)
  close(con_mtx)
  
  message(sprintf("\nConversion complete. Transposed matrix dimensions: %d cells × %d genes", 
                  n_cells, n_genes))
  message(sprintf("Output files:\n- %s (matrix)\n- %s (gene names)\n- %s (cell barcodes)", 
                  output_mtx, output_rows, output_cols))
}

convert_to_mtx_transposed(input_gz = counts.file, 
                          output_mtx = gsub(".tsv.gz", ".mtx", counts.file), 
                          output_rows = gsub(".tsv.gz", ".genes.txt", counts.file),
                          output_cols = gsub(".tsv.gz", ".barcodes.txt", counts.file))

