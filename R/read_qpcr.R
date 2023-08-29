#' Cleaning operation applied to results and amplification csv files
#'
#' @param x data frame to be cleaned
clean_common <- function(x) {
  x %>%
    janitor::clean_names() %>%
    separate(well_position, into = c('row', 'column'), sep = 1, remove = FALSE) %>%
    mutate(
      # for column, coerce to integer first to remove leading zeros
      across(column, ~as.integer(.x) %>% ordered(levels = 1:12)),
      across(row, ~ordered(.x, levels = LETTERS[1:8])),
    )
}

#' Read qPCR 'Results' or 'Standard Curve Result' table from a CSV file
#'
#' @param file path to file; see `readr::read_csv()`
#'
#' @export
read_qpcr_results_csv <- function(file) {
  # Plain results files (no standard curve calibration) have fewer columns.
  # We'll get the col types spec for just the present columns by seeing which
  # columns are present, and selecting just those files from `cts_results`
  cns <- file %>%
    readr::read_csv(
      comment = '#',
      n_max = 0,
      show_col_types = FALSE,
      progress = FALSE
    ) %>%
    names
  cts <- cts_results
  cts$cols <- cts$cols[intersect(names(cts$cols), cns)]

  file %>%
    readr::read_csv(
      comment = '#',
      col_types = cts,
    ) %>%
    clean_common() %>%
    mutate(
      cq = ifelse(cq == 'Undetermined', NA, cq) %>% as.numeric,
      cq_status = ifelse(is.na(cq), 'Undetermined', 'Determined'),
    )
}

#' Read qPCR 'Amplification Data' table from a CSV file
#'
#' @param file path to file; see `readr::read_csv()`
#'
#' @export
read_qpcr_amplification_csv <- function(file) {
  file %>%
    readr::read_csv(
      comment = '#',
      col_types = cts_amp,
    ) %>%
    clean_common()
}

# Column type specs -----------------------------------------------------------

# Column types for the 'Results' and 'Standard Curve Result' CSV files
cts_results <- readr::cols(
  Well = 'i',
  `Well Position` = 'c',
  Sample = 'c',
  Target = 'c',
  Task = 'c',
  Reporter = 'c',
  Quencher = 'c',
  `Amp Status` = 'c',
  # Cq is either a number or 'Undetermined', so read in as a character
  Cq = 'c',
  `Cq Mean` = 'n',
  `Cq Confidence` = 'n',
  `Cq SD` = 'n',
  `Auto Threshold` = 'l',
  Threshold = 'n',
  `Auto Baseline` = 'l',
  `Baseline Start` = 'i',
  `Baseline End` = 'i',
  Omit = 'l',
  # Additional columns in the 'Standard Curve Result' CSV file
  Quantity = 'n',
  Dye = 'c',
  `Quantity Mean` = 'n',
  `Quantity SD` = 'n',
  # Not sure what the 'Tm' columns are
  # Tm1,
  # Tm2,
  # Tm3,
  # Tm4,
  `Y-Intercept` = 'n',
  R2 = 'n',
  Slope = 'n',
  Efficiency = 'n',
  `Standard Deviation` = 'n',
  `Standard Error` = 'n',
)

# Column types for the amplification curves
cts_amp <- readr::cols(
  Well = 'i',
  `Well Position` = 'c',
  `Cycle Number` = 'i',
  Target = 'c',
  Rn = 'n',
  dRn = 'n',
  Sample = 'c',
  Omit = 'l'
)
