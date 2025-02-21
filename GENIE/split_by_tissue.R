library(dplyr)

# Function to split data by tissue
split_by_tissue <- function(data_file, translate_file, output_dir) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  message("Reading data table...")
  data_table <- read.csv(data_file, header = TRUE, check.names = FALSE, sep="\t")

  row_names <- paste(data_table[[1]], data_table[[2]], sep = "_")
  rownames(data_table) <- row_names

  data_table <- data_table[, -(1:2)]

  message("Reading translation table...")
  translate_table <- read.csv(translate_file, header = TRUE, check.names = FALSE, sep="\t")
  
  if (!all(c("SAMPID", "SMTS") %in% colnames(translate_table))) {
    stop("The translation table must contain columns named 'SAMPID' and 'SMTS'.")
  }
  
  message("Mapping samples to tissues...")
  column_tissue_map <- translate_table %>% 
    select(SAMPID, SMTS) %>% 
    filter(SAMPID %in% colnames(data_table))
  
  message("Splitting data table by tissue...")
  tissues <- unique(column_tissue_map$SMTS)
  
  for (tissue in tissues) {
    tissue_samples <- column_tissue_map %>% filter(SMTS == tissue) %>% pull(SAMPID)
  
    tissue_data <- data_table[, tissue_samples, drop = FALSE]
  
    row_names_combined <- rownames(tissue_data)
    #split_names <- do.call(rbind, strsplit(row_names_combined, "_"))
  
    tissue_data <- cbind(as.data.frame(row_names_combined, stringsAsFactors = FALSE), tissue_data)
  
    colnames(tissue_data)[1] <- c("Combined")
    rownames(tissue_data) <- NULL
    
    output_file <- file.path(output_dir, paste0("tissue_", gsub(" ", "_", tissue), ".csv"))
    write.csv(tissue_data, output_file, row.names = TRUE)
    message("Saved file: ", output_file)
  }
  
  message("Splitting completed.")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript split_by_tissue.R <data_file> <translate_file> <output_dir>", call. = FALSE)
}

data_file <- args[1]
translate_file <- args[2]
output_dir <- args[3]

split_by_tissue(data_file, translate_file, output_dir)