suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(here)
})

# Source shared constants (phenotype lists, label maps, colour scales)
source(here::here("R/constants.R"), local = FALSE)

# ── Core loading function ─────────────────────────────────────────────────────

.load_one_phenotype <- function(data_dir, phenotype, suffix) {
  path <- file.path(data_dir, paste0(phenotype, "_", suffix, "_t2d_assoc.txt"))
  if (!file.exists(path)) return(NULL)

  index_group <- dplyr::case_when(
    phenotype %in% FASTING_PHENOTYPES      ~ "Fasting",
    phenotype %in% OGTT_0_120_PHENOTYPES   ~ "OGTT,0-120",
    phenotype %in% OGTT_0_30_120_PHENOTYPES ~ "OGTT,0-30-120"
  )

  suppressMessages(
    readr::read_delim(path, delim = "\t", col_names = TRUE, show_col_types = FALSE) %>%
      dplyr::rename(
        nearest_gene      = locus,
        index_variant     = index_snv,
        t_2_d_risk_allele = risk_allele
      ) %>%
      dplyr::filter(!is.na(alt)) %>%
      dplyr::mutate(
        is_index    = unname(IS_INDEX_LABELS[phenotype]),
        index_group = index_group
      )
  )
}

.flip_effect <- function(df) {
  df %>%
    dplyr::mutate(
      beta = dplyr::if_else(t_2_d_risk_allele != alt, -alt_effsize, alt_effsize)
    )
}

.build_gene_snp_ra <- function(df) {
  df %>%
    dplyr::mutate(
      gene_snp_ra = paste(nearest_gene, index_variant, t_2_d_risk_allele, sep = " - ")
    )
}

# ── Public API ────────────────────────────────────────────────────────────────

#' Load and combine all phenotype association files for one analysis type.
#'
#' @param data_dir  Path to the folder containing *_bmi_t2d_assoc.txt files.
#' @param suffix    Either "bmi" or "no_bmi".
#' @return A tidy data frame ready for OJS/plotly.
load_analysis <- function(data_dir, suffix) {
  purrr::map_dfr(ALL_PHENOTYPES, ~ .load_one_phenotype(data_dir, .x, suffix)) %>%
    .flip_effect() %>%
    .build_gene_snp_ra() %>%
    dplyr::select(
      nearest_gene, index_variant, t_2_d_risk_allele, gene_snp_ra,
      is_index, index_group, beta, pvalue
    ) %>%
    dplyr::arrange(nearest_gene, index_variant)
}

#' Convenience wrapper: load both BMI-adjusted and unadjusted.
#'
#' @param data_dir  Path to the folder containing the TSV files.
#' @return Named list with $bmi and $no_bmi data frames.
load_combined_data <- function(data_dir = "data") {
  bmi    <- load_analysis(data_dir, "bmi")
  no_bmi <- load_analysis(data_dir, "no_bmi")

  if (nrow(bmi) == 0 && nrow(no_bmi) == 0) {
    warning(
      "No data files found in '", data_dir, "'. ",
      "Run data/generate_sample.R to create sample data for testing."
    )
  }

  list(bmi = bmi, no_bmi = no_bmi)
}
