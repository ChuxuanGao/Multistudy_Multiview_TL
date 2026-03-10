#' Summarize Simulation Results Across Repetitions
#'
#' `sim_sum()` runs evaluation across one or more repetitions in a simulation
#' object and returns a streamlined summary including prediction error
#' statistics and beta-coefficient comparisons.  
#'
#' The function detects valid repetition indices from `studies_obj`, aligns
#' requested repetitions to those available, and iteratively evaluates each
#' repetition using `get_fit_eval()`.  
#' Summary statistics (mean, SD, Monte Carlo SE, and 95% CIs) are computed for
#' overall prediction errors.
#'
#' @param studies_obj A simulation object containing elements such as
#'   `true_betas`, `trainDat`, `testDat`, and optionally `nstudy`.  
#'   Must be formatted such that each of these components may include
#'   repetition-specific entries named `"Rep_k"` (e.g., `"Rep_1"`, `"Rep_2"`).
#'
#' @param reps Integer vector of repetition indices to run.  
#'   Defaults to `1:rep` (user-defined external variable). The function will
#'   automatically adjust requested reps to those available within
#'   `studies_obj$true_betas` if detectable.
#'
#' @param verbose Logical. If `TRUE` (default), prints progress messages,
#'   including rep-index adjustments and run-progress indication.
#'
#' @return A list with two items:
#' \describe{
#'
#'   \item{`prediction_summary`}{
#'     A data frame summarizing prediction-error metrics across repetitions.
#'     Contains the following columns:
#'     \itemize{
#'       \item `mean` — Mean of the metric across available repetitions.
#'       \item `sd` — Standard deviation across repetitions.
#'       \item `mcse` — Monte Carlo standard error (`sd / sqrt(R)`).
#'       \item `ci_lower` — Lower 95% confidence bound.
#'       \item `ci_upper` — Upper 95% confidence bound.
#'     }
#'
#'     Row names correspond to the prediction-error components:
#'     \itemize{
#'       \item `overall`
#'       \item `individual_mean`
#'       \item `pretrain_mean`
#'     }
#'   }
#'
#'   \item{`beta_comparison_list`}{
#'     A named list of length equal to the number of valid repetitions
#'     evaluated.  
#'     Each element corresponds to `"Rep_k"` and contains the beta-coefficient
#'     comparison table returned by `get_fit_eval()`.
#'   }
#'
#' }
#'
#' @export

sim_sum <- function(studies_obj, family, reps = 1:rep, verbose = TRUE) {
  # ---- detect available reps (unchanged) ----
  available_reps <- NULL
  if (!is.null(studies_obj$true_betas) && length(studies_obj$true_betas) > 0) {
    rep_names <- names(studies_obj$true_betas)
    rep_idx <- suppressWarnings(as.integer(gsub(".*?(\\d+).*", "\\1", rep_names)))
    available_reps <- rep_idx[!is.na(rep_idx)]
  }
  
  if (!is.null(available_reps) && length(available_reps) > 0) {
    requested_reps <- unique(as.integer(reps))
    valid_reps <- intersect(requested_reps, available_reps)
    if (length(valid_reps) == 0)
      stop("Requested reps not available in studies_obj$true_betas.")
    if (verbose && length(valid_reps) < length(reps))
      message("Adjusted reps to available: ", paste(valid_reps, collapse = ", "))
    reps <- sort(valid_reps)
  } else {
    reps <- unique(as.integer(reps))
    if (verbose) message("No rep index detected; using reps as provided.")
  }
  
  # ---- detect nstudy (same code) ----
  nstudy <- NA
  if (!is.null(studies_obj$nstudy)) nstudy <- as.integer(studies_obj$nstudy)
  if (is.na(nstudy) &&
      !is.null(studies_obj$trainDat) &&
      !is.null(studies_obj$trainDat$Rep_1))
    nstudy <- length(grep("^Study_", names(studies_obj$trainDat$Rep_1)))
  
  # ---- containers ----
  pred_errors_list <- list()
  metrics_list <- list()
  beta_comparison_list <- vector("list", length(reps))
  names(beta_comparison_list) <- paste0("Rep_", reps)
  
  # ---- main loop ----
  for (i in seq_along(reps)) {
    r <- reps[i]
    if (verbose) message(sprintf("Running repetition %s (%d of %d)...",
                                 r, i, length(reps)))
    
    studies_obj_rep <- studies_obj
    rep_name_r <- paste0("Rep_", r)
    
    if (!is.null(studies_obj$true_betas) &&
        rep_name_r %in% names(studies_obj$true_betas))
      studies_obj_rep$true_betas$Rep_1 <- studies_obj$true_betas[[rep_name_r]]
    
    if (!is.null(studies_obj$trainDat) &&
        rep_name_r %in% names(studies_obj$trainDat))
      studies_obj_rep$trainDat$Rep_1 <- studies_obj$trainDat[[rep_name_r]]
    
    if (!is.null(studies_obj$testDat) &&
        rep_name_r %in% names(studies_obj$testDat))
      studies_obj_rep$testDat$Rep_1 <- studies_obj$testDat[[rep_name_r]]
    
    eval_out <- tryCatch(
      get_fit_eval(studies_obj_rep, family),
      error = function(e) {
        warning(sprintf(
          "Replication %s failed (treated as non-estimable): %s",
          r, e$message
        ))
        return(NULL)
      }
    )
    
    if (is.null(eval_out)) next
    
    # ---- collect beta_comparison ----
    beta_comparison_list[[i]] <- eval_out$beta_comparison
    
    # ---- collect prediction errors ----
    pe <- eval_out$prediction_errors
    pred_errors_list[[i]] <- c(
      overall         = pe$overall         %||% NA,
      individual_mean = pe$individual_mean %||% NA,
      pretrain_mean   = pe$pretrain_mean   %||% NA
    )
  }
  
  # ---- prediction summary ----
  pred_mat <- do.call(rbind, pred_errors_list)
  pred_df <- as.data.frame(pred_mat)
  
  summarize_vec <- function(x) {
    x <- x[!is.na(x)]
    R <- length(x)
    if (R == 0)
      return(data.frame(mean = NA, sd = NA, mcse = NA, ci_lower = NA, ci_upper = NA))
    m <- mean(x); s <- sd(x); mcse <- s/sqrt(R)
    data.frame(mean = m, sd = s, mcse = mcse,
               ci_lower = m - 1.96 * mcse,
               ci_upper = m + 1.96 * mcse)
  }
  
  prediction_summary <- do.call(
    rbind,
    lapply(colnames(pred_df), function(col) summarize_vec(pred_df[[col]]))
  )
  rownames(prediction_summary) <- colnames(pred_df)
  
  if (verbose) message("sim_sum completed.")
  
  # ---- final streamlined output (ONLY 2 ITEMS) ----
  list(
    prediction_summary = prediction_summary,
    beta_comparison_list = beta_comparison_list
  )
}
