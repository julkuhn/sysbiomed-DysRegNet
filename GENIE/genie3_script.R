library(GENIE3)
library(stringr)
library(do)

run_GENIE3 <- function(file_path, regulator_file = NULL, normalize = FALSE) {
  # Load the required libraries
  if (!requireNamespace("GENIE3", quietly = TRUE)) {
    stop("The GENIE3 package is required but not installed. Please install it first.")
  }
  
  # Set seed for reproducibility
  set.seed(123)
  
  # Load the expression matrix
  message("Reading expression matrix from file: ", file_path)
  exprMatr <- tryCatch(
    {
      as.matrix(read.table(file_path, header = TRUE, row.names = 1, sep = ","))
    },
    error = function(e) {
      stop("Error reading the file. Please ensure the file is a valid expression matrix: ", e$message)
    }
  )
  
  # Normalize the expression matrix if requested
  if (normalize) {
    message("Normalizing the expression matrix...")
    exprMatr <- log2(exprMatr + 0.001)  
    exprMatr <- t(scale(t(exprMatr), center = TRUE, scale = TRUE))
    message("Normalization completed.")
    name <- file.name(file_path)
    name <- str_extract(name, '.*(?=\\.csv)')
    name <- gsub("[[:space:]]", "", name)
    output_file <- paste("/nfs/data2/dysregnet_gtex/GENIE3_output/normalized_", name, ".csv")
    write.table(exprMatr, file = output_file, sep = ",", row.names = TRUE, quote = FALSE)
    message("Normalization saved at: ", output_file)
  }
  

  # Load candidate regulators if provided
  candidateRegulators <- NULL
  if (!is.null(regulator_file)) {
    message("Reading candidate regulators from file: ", regulator_file)
    candidateRegulators <- tryCatch(
      {
        readLines(regulator_file)
      },
      error = function(e) {
        stop("Error reading the candidate regulators file: ", e$message)
      }
    )
    message("Filtering candidate regulators...")
    candidateRegulators <- intersect(candidateRegulators, rownames(exprMatr))
    message(length(candidateRegulators), " candidate regulators remain after filtering.")
  }
  
  # Run GENIE3 with all available cores
  message("Running GENIE3 using all available cores...")
  nCores <- 50
  if(!is.null(candidateRegulators)){
    weightMat <- GENIE3(exprMatr, regulators = candidateRegulators, nCores = nCores)
  } else {
    weightMat <- GENIE3(exprMatr, nCores = nCores)
  }
  
  # Generate the linked list
  message("Generating the link list...")
  linkedList <- getLinkList(weightMat, reportMax = 100000)
  
  return(linkedList)
}

# Command-line Interface
args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 1) {
  file_path <- args[1]
  regulator_file <- if (length(args) > 1 && args[2] != "") args[2] else NULL
  normalize <- if (length(args) > 2) as.logical(args[3]) else FALSE
  
  if (!is.null(regulator_file)) {
    linkedList <- run_GENIE3(file_path, regulator_file, normalize)
  } else {
    linkedList <- run_GENIE3(file_path, normalize = normalize)
  }

  # Save the linked list to a file
  name <- basename(file_path)
  name <- str_extract(name, '.*(?=\\.csv)')
  name <- gsub("[[:space:]]", "", name)
  output_file <- paste0("/nfs/data2/dysregnet_gtex/results_array/linkedList_output_slurm_", name, ".csv")
  write.table(linkedList, file = output_file, sep = ",", row.names = FALSE, quote = FALSE)
  message("Linked list saved to: ", output_file)

} else {
  message("Usage: Rscript script_name.R <expression_matrix_file> [<regulator_file>] [<normalize>]")
}
