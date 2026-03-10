#' Run Multi–Share-Level Simulation and Compile Error Summaries
#'
#' `sim_combo()` automates a full simulation pipeline over multiple
#' beta-sharing levels (`rho.beta`) and compiles the resulting error‐rate
#' summaries into a single stacked output table.  
#'
#' For each share level, the function:
#' \enumerate{
#'   \item Generates multi-study simulated datasets using `gen_simmba_multistudy()`.
#'   \item Computes summary metrics across repetitions using `sim_sum()`.
#'   \item Extracts and annotates the prediction-error summary table.
#' }
#'
#' All share-level summaries are vertically stacked and written into a single
#' Excel file, allowing convenient comparison of model performance as the
#' strength of cross-study signal sharing varies.
#'
#' @param share_levels A numeric vector specifying the values of `rho.beta`
#'   (share levels) to evaluate. Each value represents the proportion of
#'   cross-study shared signal, typically between 0 (no sharing) and 1
#'   (fully shared).
#'
#' @param N Integer. Sample size per study to pass to
#'   `gen_simmba_multistudy()` via the `nsample` argument.
#'
#' @param S Integer. Number of studies, forwarded to `gen_simmba_multistudy()`
#'   through the `nstudy` argument.
#'
#' @param SNR_SD Numeric. Standard deviation of study-specific SNR used in
#'   `gen_simmba_multistudy()` (controls heterogeneity across studies).
#'
#' @param rep Integer. Number of repetitions (`nrep`) to simulate for each
#'   share level. Also sets the repetition indices passed to `sim_sum()`.
#'
#' @param SNR Numeric. Global signal-to-noise ratio forwarded to
#'   `gen_simmba_multistudy()`. Default is `5`.
#'
#' @param SigAlp Numeric. Variance parameter for study-specific random effects
#'   (`sigma.alpha`) in the simulation generator. Default is `0.1`.
#'
#' @param outfile Character string. File name (including extension `.xlsx`)
#'   for saving the output workbook.  
#'   Default: `"sim_combo_output.xlsx"`.
#'
#' @return A tibble combining prediction-error summaries across all specified
#'   share levels.  
#'
#'   The returned table contains the following columns:
#'   \itemize{
#'     \item `metric` — Error metric name (`overall`, `individual_mean`, `pretrain_mean`).
#'     \item `mean`, `sd`, `mcse`, `ci_lower`, `ci_upper` — Summary statistics
#'       computed by `sim_sum()`.
#'     \item `share_level` — The `rho.beta` value corresponding to each block
#'       of results.
#'   }
#'
#'   The output is also written to an Excel file with a single sheet named
#'   `"error_summary"`.
#'
#' @export

sim_combo <- function(share_levels, N, S, SNR_SD, rep, outcome, 
                      SNR = 5, SigAlp = 0.1,
                      outfile = "sim_combo_output.xlsx") {
  
  outcome <- match.arg(outcome,
                       choices = c("continuous", "binary", "survival"))
  
  # Store stacked results
  combined_error <- list()
  
  # Loop through share levels
  for (i in seq_along(share_levels)) {
    
    share_level <- share_levels[i]
    cat(sprintf("Running (%d/%d): share_level = %s\n", 
                i, length(share_levels), share_level))
    
    # ---- Step 1: Simulation ----
    simulated_studies <- gen_simmba_multistudy(
      nsample = N,
      nstudy = S,
      snr = SNR,
      rho.beta = share_level,
      sigma.alpha = SigAlp,
      tau.snr = SNR_SD,
      outcome.type = outcome,
      nrep = rep,
      seed = 1234
    )
    
    # map outcome -> family
    family <- switch(outcome,
                     continuous = "gaussian",
                     binary     = "binomial",
                     survival   = "cox",
                     stop("Unsupported outcome type: ", outcome)
    )
    
    # ---- Step 2: Summaries ----
    simulation_summary <- sim_sum(
      simulated_studies,
      family = family,
      reps = 1:rep,
      verbose = TRUE
    )
    
    sim_error <- simulation_summary$prediction_summary
    
    # ---- Step 3: Add rownames and share level annotation ----
    sim_error <- tibble::rownames_to_column(sim_error, var = "metric")
    sim_error$share_level <- share_level
    
    # Store
    combined_error[[i]] <- sim_error
  }
  
  # ---- Stack all share-level results vertically ----
  final_error <- bind_rows(combined_error)
  
  # ---- Write to Excel (single sheet) ----
  wb <- createWorkbook()
  addWorksheet(wb, "error_summary")
  writeData(wb, "error_summary", final_error)
  saveWorkbook(wb, outfile, overwrite = TRUE)
  
  cat("Finished! Output saved to:", outfile, "\n")
  
  return(final_error)
}
