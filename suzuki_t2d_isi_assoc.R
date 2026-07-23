.libPaths(Sys.getenv("R_LIBS_USER"))
library(dplyr)
library(readxl)
library(snakecase)
library(tidyr)

# Specify the path to your Excel file
file <- ("/projects/glostrup-AUDIT/people/rnh585/insulin_sensitivity/gwas/t2d_snps/Suzuki_41586_2024_7019_MOESM3_ESM.xlsx")

# Read the Excel file and perform data preprocessing
suzuki_t2d_snps <- read_excel(file, sheet = "ST4", col_names = TRUE, skip = 2) %>%
  dplyr::slice(-1) %>%
  dplyr::rename_all(to_snake_case) %>% 
  dplyr::rename(
    CHROM = chromosome,
    POS = position_bp_b_37,
    risk_allele = alleles,
    other_allele = `8`,
    p_value = mr_mega_association_p_value
  ) %>%
  dplyr::select(
    locus,
    index_snv,
    CHROM,
    POS,
    risk_allele,
    other_allele,
    p_value
  )%>% 
  tidyr::fill(locus, .direction = "down") %>%
  tidyr::fill(CHROM, .direction = "down") %>%
  dplyr::mutate(locus = trimws(sub(",.*", "", locus))) %>%
  dplyr::filter(!is.na(CHROM), !is.na(POS)) %>%
  dplyr::mutate(CHROM = as.integer(CHROM), POS = as.integer(POS))



# Set directory containing _manhattan.txt files
parent_dir <- "/projects/glostrup-AUDIT/people/rnh585/insulin_sensitivity/gwas/results/inter99/linear_model"

# List of file names in the destination folder
file_names <- list.files(parent_dir, pattern = "*_merged.txt", full.names = TRUE)

# Extract the desired basenames
file_basenames <- gsub(".*/|_merged\\.txt$", "", file_names)

# Function to process each file
process_file <- function(file_basename) {
  # Construct the full file path
  full_file_path <- file.path(parent_dir, paste0(file_basename, "_merged.txt"))
  
  # Define output file name
  output_file_name <- file.path(parent_dir, paste0(file_basename, "_t2d_assoc.txt"))
  
  # Read the file
  data <- read.delim(full_file_path, header = TRUE, sep = "\t", dec = ".") %>%
    dplyr::select(CHROM, POS, ALT, ALT_EFFSIZE, PVALUE)
  
  
  # Merge data based on CHR, POS, columns
  merged_data <- left_join(suzuki_t2d_snps, data, by = c("CHROM", "POS")) %>%
    dplyr::rename_all(to_snake_case)
  
  # Write the result to a new file (overwrite if exists)
  write.table(merged_data, output_file_name, sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat("Processed and overwritten (if existed): file ", file_basename, "\n")  # Log processing
}

# Apply the function to each file
lapply(file_basenames, process_file)

# Run it in the terminal with Rscript

