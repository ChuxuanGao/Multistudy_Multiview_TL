# -------------------------
# Helper: safe coefficient extractor
# -------------------------
extract_coef_df <- function(fit_obj, s = "lambda.1se") {
  # Returns a one-column data.frame of coefficients with intercept-like rows removed.
  # fit_obj: object accepted by coef(..., s = s)
  mat <- as.matrix(coef(fit_obj, s = s))
  rn  <- rownames(mat)
  # Remove rows that look like intercepts (case-insensitive)
  keep <- !grepl("intercept", rn, ignore.case = TRUE)
  if (all(!keep)) {
    # If nothing left after removing intercept-like rows, return the full vector
    keep <- rep(TRUE, length(keep))
  }
  df <- as.data.frame(mat[keep, , drop = FALSE])
  rownames(df) <- rn[keep]
  return(df)
}

# -------------------------
# 1) retrieve_true_beta: flexible rep_label; returns data.frame of true betas
#    Columns: shared, indi_1..indi_n, true_1..true_n
# -------------------------
retrieve_true_beta <- function(studies_obj, rep_label = "Rep_1") {
  # rep_label should match the names used in studies_obj$true_betas (e.g. "Rep_1")
  if (is.null(studies_obj$true_betas[[rep_label]])) {
    stop(sprintf("retrieve_true_beta: %s not found in studies_obj$true_betas", rep_label))
  }
  rep_obj <- studies_obj$true_betas[[rep_label]]
  
  # shared
  shared <- rep_obj$beta_shared
  
  # infer nstudy from beta_indiv or beta_s if present
  nstudy <- NA
  if (!is.null(rep_obj$beta_indiv)) {
    nstudy <- length(rep_obj$beta_indiv)
    indiv_list <- rep_obj$beta_indiv
  } else {
    indiv_list <- list()
  }
  if (!is.null(rep_obj$beta_s)) {
    nstudy_pre <- length(rep_obj$beta_s)
    pre_list <- rep_obj$beta_s
    if (is.na(nstudy)) nstudy <- nstudy_pre
  } else {
    pre_list <- list()
  }
  if (is.na(nstudy)) {
    stop("Could not infer number of studies from true betas (no beta_indiv or beta_s found).")
  }
  
  # Build a list of columns: first shared, then indi_1.., then true_1..
  out_cols <- list(shared = shared)
  
  # indi
  for (s in seq_len(nstudy)) {
    colname <- paste0("indi_", s)
    if (!is.null(indiv_list[[s]])) {
      out_cols[[colname]] <- indiv_list[[s]]
    } else {
      # fill with NA vector of appropriate length
      out_cols[[colname]] <- rep(NA, length(shared))
    }
  }
  # true pretrained (beta_s)
  for (s in seq_len(nstudy)) {
    colname <- paste0("true_", s)
    if (!is.null(pre_list[[s]])) {
      out_cols[[colname]] <- pre_list[[s]]
    } else {
      out_cols[[colname]] <- rep(NA, length(shared))
    }
  }
  
  rep_betas <- as.data.frame(do.call(cbind, out_cols), stringsAsFactors = FALSE)
  
  # set rownames from feature metadata if available
  # try trainDat[[rep_label]]$Study_1$feature_metadata$featureID
  try({
    fm <- studies_obj$trainDat[[rep_label]]
    if (!is.null(fm) && "Study_1" %in% names(fm)) {
      feat_ids <- fm$Study_1$feature_metadata$featureID
      if (!is.null(feat_ids) && length(feat_ids) == nrow(rep_betas)) {
        rownames(rep_betas) <- as.vector(feat_ids)
      }
    }
  }, silent = TRUE)
  
  return(rep_betas)
}

# -------------------------
# 2) get_train_test: flexible rep_label; works for arbitrary nstudy
# -------------------------
get_train_test <- function(studies_obj, rep_label = "Rep_1") {
  # expects studies_obj$trainDat[[rep_label]] and studies_obj$testDat[[rep_label]]
  if (is.null(studies_obj$trainDat[[rep_label]])) {
    stop(sprintf("get_train_test: %s not found in studies_obj$trainDat", rep_label))
  }
  if (is.null(studies_obj$testDat[[rep_label]])) {
    stop(sprintf("get_train_test: %s not found in studies_obj$testDat", rep_label))
  }
  
  train_rep <- studies_obj$trainDat[[rep_label]]
  test_rep  <- studies_obj$testDat[[rep_label]]
  
  # detect study names in the rep (expect Study_1, Study_2, ...)
  study_names <- grep("^Study_", names(train_rep), value = TRUE)
  if (length(study_names) == 0) {
    # fallback: use all elements that look like feature_table entries
    study_names <- names(train_rep)
  }
  nstudy <- length(study_names)
  
  # collect train X/Y
  X_train_list <- vector("list", nstudy)
  Y_train_list <- vector("list", nstudy)
  Groups <- integer(0)
  for (i in seq_along(study_names)) {
    sn <- study_names[i]
    # feature_table should exist
    X_i <- t(train_rep[[sn]]$feature_table)
    y_i <- train_rep[[sn]]$sample_metadata$Y
    X_train_list[[i]] <- X_i
    Y_train_list[[i]] <- y_i
    Groups <- c(Groups, rep(i, length(y_i)))
  }
  X_train <- do.call(rbind, X_train_list)
  Y_train <- unlist(Y_train_list, use.names = FALSE)
  
  # collect test X/Y
  X_test_list <- vector("list", nstudy)
  Y_test_list <- vector("list", nstudy)
  Groups_test <- integer(0)
  for (i in seq_along(study_names)) {
    sn <- study_names[i]
    X_i_test <- t(test_rep[[sn]]$feature_table)
    y_i_test <- test_rep[[sn]]$sample_metadata$Y
    X_test_list[[i]] <- X_i_test
    Y_test_list[[i]] <- y_i_test
    Groups_test <- c(Groups_test, rep(i, length(y_i_test)))
  }
  X_test <- do.call(rbind, X_test_list)
  Y_test <- unlist(Y_test_list, use.names = FALSE)
  
  return(list(
    X_train = X_train,
    Y_train = Y_train,
    Groups = Groups,
    X_test = X_test,
    Y_test = Y_test,
    Groups_test = Groups_test,
    nstudy = nstudy,
    study_names = study_names
  ))
}

# -------------------------
# 3) fit_pred_evaluate: fits cv.ptLasso, extracts estimated betas and metrics
#    takes rep_betas (true betas data.frame with standardized names)
# -------------------------
fit_pred_evaluate <- function(X_train, Y_train, Groups,
                              X_test, Y_test, Groups_test,
                              family,
                              rep_betas, nstudy = NULL) {
  
  # Fit
  est  <- cv.ptLasso(X_train, Y_train, Groups, family, alphalist = seq(0, 1, length.out=11))
  pred <- predict(est, X_test, groupstest = Groups_test, ytest = Y_test)
  ptLa <- est$fit

  ## -------------------------------------------------
  ## identify alpha-matching fit
  ## -------------------------------------------------
  alpha_vals <- vapply(ptLa, function(x) x$alpha, numeric(1))
  idx_alpha  <- which.min(abs(alpha_vals - est$alphahat))
  
  if (length(idx_alpha) != 1) {
    stop("Could not uniquely identify ptLa element matching alphahat.")
  }
  
  fit_alpha <- ptLa[[idx_alpha]]
  
  ## -------------------------------------------------
  ## determine nstudy
  ## -------------------------------------------------
  if (is.null(nstudy)) {
    # assume fitind length gives nstudy
    if (!is.null(fit_alpha$fitind)) {
      nstudy <- length(fit_alpha$fitind)
    } else {
      stop("fit_pred_evaluate: cannot infer nstudy from fit object; please supply nstudy.")
    }
  }
  ## -------------------------------------------------
  ## Extract estimated betas
  ## -------------------------------------------------
  est_cols <- list()
  p <- nrow(rep_betas)
  
  # overall/shared
  overall_df <- tryCatch(
    extract_coef_df(fit_alpha$fitoverall, s = "lambda.1se"), 
    error = function(e) data.frame(row = character(0), stringsAsFactors = FALSE))
  # ensure same row order as rep_betas
  # overall_df is one-column; convert to vector
  shared_est <- as.numeric(overall_df[,1, drop = TRUE])
  names(shared_est) <- rownames(overall_df)
  
  shared_est_vec <- rep(NA, p)
  common <- intersect(rownames(rep_betas), names(shared_est))
  if (length(common) > 0) {
    shared_est_vec[match(common, rownames(rep_betas))] <- shared_est[common]
  }
  
  est_cols[["shared_est"]] <- shared_est_vec
  
  # ----- individual fits -----
  for (s in seq_len(nstudy)) {
    fitind_s <- fit_alpha$fitind[[s]]
    ci_df <- extract_coef_df(fitind_s, s = "lambda.1se")
    
    est_vec <- as.numeric(ci_df[, 1, drop = TRUE])
    names(est_vec) <- rownames(ci_df)
    
    out_vec <- rep(NA, p)
    common <- intersect(rownames(rep_betas), names(est_vec))
    if (length(common) > 0) {
      out_vec[match(common, rownames(rep_betas))] <- est_vec[common]
    }
    
    est_cols[[paste0("indi_", s, "_est")]] <- out_vec
  }
  
  # ----- pretrained fits -----
  for (s in seq_len(nstudy)) {
    fitpre_s <- fit_alpha$fitpre[[s]]
    cp_df <- extract_coef_df(fitpre_s, s = "lambda.1se")
    
    est_vec <- as.numeric(cp_df[, 1, drop = TRUE])
    names(est_vec) <- rownames(cp_df)
    
    out_vec <- rep(NA, p)
    common <- intersect(rownames(rep_betas), names(est_vec))
    if (length(common) > 0) {
      out_vec[match(common, rownames(rep_betas))] <- est_vec[common]
    }
    
    est_cols[[paste0("true_", s, "_est")]] <- out_vec
  }
  
  ## -------------------------------------------------
  ## combine estimates
  ## -------------------------------------------------
  rep_betas_est <- as.data.frame(do.call(cbind, est_cols), stringsAsFactors = FALSE)
  rownames(rep_betas_est) <- rownames(rep_betas)
  
  beta_comparison <- cbind(rep_betas, rep_betas_est)
  
  ## -------------------------------------------------
  ## prediction errors
  ## -------------------------------------------------
  
  er_overall  <- if (!is.null(pred$erroverall)) pred$erroverall[[1]] else NA
  er_ind_mean <- if (!is.null(pred$errind)) pred$errind[[1]] else NA
  er_pt_mean  <- if (!is.null(pred$errpre)) pred$errpre[[1]] else NA
  
  ## -------------------------------------------------
  ## metrics
  ## -------------------------------------------------
  
  compute_metrics <- function(beta_true, beta_est) {
    true_nonzero <- as.integer(!is.na(beta_true) & beta_true != 0)
    est_nonzero  <- as.integer(!is.na(beta_est) & beta_est != 0)
    
    TP <- sum(true_nonzero == 1 & est_nonzero == 1, na.rm = TRUE)
    FN <- sum(true_nonzero == 1 & est_nonzero == 0, na.rm = TRUE)
    FP <- sum(true_nonzero == 0 & est_nonzero == 1, na.rm = TRUE)
    TN <- sum(true_nonzero == 0 & est_nonzero == 0, na.rm = TRUE)
    
    Sensitivity <- if ((TP + FN) == 0) NA else TP / (TP + FN)
    Specificity <- if ((TN + FP) == 0) NA else TN / (TN + FP)
    Precision   <- if ((TP + FP) == 0) NA else TP / (TP + FP)
    F1_score    <- if (is.na(Precision) || is.na(Sensitivity) || (Precision + Sensitivity) == 0) NA else 2 * Precision * Sensitivity / (Precision + Sensitivity)
    MSE         <- mean((beta_true - beta_est)^2, na.rm = TRUE)
    
    list(Sensitivity = Sensitivity,
         Specificity = Specificity,
         Precision = Precision,
         F1_score = F1_score,
         MSE = MSE)
  }
  
  # Compute metrics for all standardized true columns
  metric_list <- list()
  true_col_names <- names(rep_betas)
  for (tc in true_col_names) {
    est_col <- paste0(tc, "_est")
    if (! (est_col %in% names(beta_comparison)) ) {
      # try alternative exact match (some user code might have different separators)
      alt_match <- grep(paste0("^", tc, "(_|)est$"), names(beta_comparison), value = TRUE, ignore.case = TRUE)
      if (length(alt_match) >= 1) est_col <- alt_match[1] else est_col <- NA
    }
    if (is.na(est_col)) {
      m <- list(Sensitivity = NA, Specificity = NA, Precision = NA, F1_score = NA, MSE = NA)
    } else {
      m <- compute_metrics(beta_comparison[[tc]], beta_comparison[[est_col]])
    }
    # store under standardized name, e.g., "shared", "indi_1", "true_1"
    metric_list[[tc]] <- m
  }
  
  return(list(
    rep_betas_est   = rep_betas_est,
    beta_comparison = beta_comparison,
    prediction_errors = list(overall = er_overall, individual_mean = er_ind_mean, pretrain_mean = er_pt_mean),
    metrics = metric_list
  ))
}

# -------------------------
# 4) get_fit_eval: wrapper for one repetition (uses the above funcs)
# -------------------------
get_fit_eval <- function(studies_obj, family, rep_label = "Rep_1") {
  # STEP 1: true betas (standardized names)
  rep_betas <- retrieve_true_beta(studies_obj, rep_label = rep_label)
  
  # STEP 2: train/test
  tt <- get_train_test(studies_obj, rep_label = rep_label)
  
  # STEP 3: fit, predict, evaluate
  eval_out <- fit_pred_evaluate(
    X_train = tt$X_train,
    Y_train = tt$Y_train,
    Groups = tt$Groups,
    X_test = tt$X_test,
    Y_test = tt$Y_test,
    Groups_test = tt$Groups_test,
    family = family,
    rep_betas = rep_betas,
    nstudy = tt$nstudy
  )
  
  # combine into tidy return
  return(list(
    true_betas       = rep_betas,
    estimated_betas  = eval_out$rep_betas_est,
    beta_comparison  = eval_out$beta_comparison,
    prediction_errors = eval_out$prediction_errors,
    metrics          = eval_out$metrics,
    response = family
  ))
}
