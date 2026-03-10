cv.multiview.revised <- function(x_list, y,
                                 family = binomial(),
                                 lambda = NULL,
                                 s = c("lambda.1se", "lambda.min"),
                                 rho = seq(0, 1, 0.25),
                                 nfolds = 10,
                                 foldid = NULL,
                                 weights = NULL,
                                 offset = NULL,
                                 penalty.factor = NULL,
                                 type.measure = NULL,
                                 alignment = c("lambda", "fraction"),
                                 trace.it = 0,
                                 keep = FALSE,
                                 ...) {
  
  mv_fit <- function(x_list, y, family, rho, lambda = NULL,
                     weights = NULL, offset = NULL,
                     penalty.factor = NULL, trace.it = 0, ...) {
    
    args <- list(
      x_list = x_list,
      y = y,
      family = family,
      rho = rho,
      weights = weights,
      offset = offset,
      lambda = lambda,
      trace.it = trace.it,
      ...
    )
    if (!is.null(penalty.factor)) args$penalty.factor <- penalty.factor
    
    fit_once <- function(extra = list()) {
      do.call(multiview::multiview, c(args, extra))
    }
    
    # First attempt
    out <- tryCatch(
      fit_once(),
      error = function(e) e
    )
    if (!inherits(out, "error")) return(out)
    
    msg <- conditionMessage(out)
    
    # If we hit the step-size failure, back off to safer settings.
    if (grepl("cannot correct step size", msg, fixed = TRUE)) {
      # Key idea: avoid tiny lambdas + allow more iterations.
      # These are passed through ... to the underlying solver (glmnet-like).
      out2 <- tryCatch(
        fit_once(extra = list(lambda.min.ratio = 0.05, nlambda = 50, maxit = 1e5)),
        error = function(e) e
      )
      if (!inherits(out2, "error")) return(out2)
    }
    
    # Otherwise: rethrow the original error with context
    stop(out)
  }
  
  # ------------------------------------------------------------
  # Binomial: create stratified folds when foldid is NULL
  # and validate that EACH fold has both classes in TRAINING.
  # ------------------------------------------------------------
  make_stratified_foldid <- function(y01, K) {
    y01 <- as.integer(y01)
    if (!all(y01 %in% c(0L, 1L))) stop("For binomial, y must be 0/1.")
    idx0 <- which(y01 == 0L); idx1 <- which(y01 == 1L)
    
    f0 <- sample(rep(seq_len(K), length.out = length(idx0)))
    f1 <- sample(rep(seq_len(K), length.out = length(idx1)))
    
    foldid <- integer(length(y01))
    foldid[idx0] <- f0
    foldid[idx1] <- f1
    foldid
  }
  
  validate_binomial_folds <- function(y01, foldid) {
    K <- max(foldid)
    for (f in seq_len(K)) {
      is_test <- foldid == f
      y_train <- y01[!is_test]
      if (length(unique(y_train)) < 2) {
        stop(sprintf(
          "Invalid CV split for binomial: training set in fold %d has only one class. Use stratified foldid or reduce nfolds.",
          f
        ))
      }
    }
    invisible(TRUE)
  }
  
  ## -----------------------------
  ## local helper: build predmat like cv.multiview()
  ## -----------------------------
  build_predmat_mv <- function(outlist, lambda, x_list, offset, foldid,
                               alignment = c("lambda", "fraction"),
                               family, type = "response", ...) {
    
    alignment <- match.arg(alignment)
    N <- nrow(x_list[[1L]])
    
    if (!is.null(offset)) offset <- drop(offset)
    
    predmat <- matrix(NA_real_, N, length(lambda))
    nfolds <- max(foldid)
    nlambda <- length(lambda)
    
    for (i in seq_len(nfolds)) {
      which <- (foldid == i)
      fitobj <- outlist[[i]]
      x_sub_list <- lapply(x_list, function(x) x[which, , drop = FALSE])
      offset_sub <- if (is.null(offset)) NULL else offset[which]
      
      preds <- switch(
        alignment,
        fraction = predict(fitobj, newx = x_sub_list, newoffset = offset_sub,
                           type = type, ...),
        lambda   = predict(fitobj, newx = x_sub_list, s = lambda, newoffset = offset_sub,
                           type = type, ...)
      )
      
      nlami <- min(ncol(preds), nlambda)
      predmat[which, seq_len(nlami)] <- preds[, seq_len(nlami), drop = FALSE]
      
      # pad to the right with last column if needed (matches glmnet/cv behavior)
      if (nlami < nlambda) {
        predmat[which, (nlami + 1L):nlambda] <- preds[, nlami]
      }
    }
    
    rn <- rownames(x_list[[1L]])
    sn <- paste0("s", seq_len(nlambda))
    dimnames(predmat) <- list(rn, sn)
    attr(predmat, "family") <- family
    predmat
  }
  
  ## -----------------------------
  ## setup
  ## -----------------------------
  alignment <- match.arg(alignment)
  N <- nrow(x_list[[1]])
  if (any(vapply(x_list, function(m) nrow(m), integer(1)) != N))
    stop("All views in x_list must have the same number of rows.")
  
  y <- drop(y)
  if (length(y) != N)
    stop("Length of y must equal nrow(x_list[[1]]).")
  
  if (is.null(weights)) weights <- rep(1, N)
  if (length(weights) != N) stop("Length of weights must equal N.")
  if (!is.null(offset) && length(offset) != N)
    stop("Length of offset must equal N when provided.")
  
  if (is.null(type.measure)) {
    if (family$family == "gaussian") type.measure <- "mse"
    else if (family$family == "binomial") type.measure <- "auc"
    else stop("Default type.measure only defined for gaussian and binomial.")
  }
  
  # tolerate typo "lambda.lse"
  if (length(s) != 1) s <- s[1]
  if (is.character(s) && identical(s, "lambda.lse")) s <- "lambda.1se"
  s <- match.arg(s, choices = c("lambda.1se", "lambda.min"))
  
  rho <- sort(unique(rho))
  
  ## -----------------------------
  ## folds
  ## ----------------------------
  
  if (is.null(foldid)) {
    if (family$family == "binomial") {
      y01 <- as.integer(as.numeric(drop(y)) > 0)
      foldid <- make_stratified_foldid(y01, nfolds)
    } else {
      foldid <- sample(rep(seq_len(nfolds), length.out = N))
    }
  } else {
    if (length(foldid) != N) stop("foldid must have length N (nrow of data).")
    foldid <- as.integer(foldid)
    used <- sort(unique(foldid))
    fold_map <- setNames(seq_along(used), used)
    foldid <- as.integer(fold_map[as.character(foldid)])
    nfolds <- length(unique(foldid))
    if (nfolds < 2) stop("foldid must contain at least 2 unique folds.")
  }
  
  # after renumbering:
  if (family$family == "binomial") {
    y01 <- as.integer(as.numeric(drop(y)) > 0)
    validate_binomial_folds(y01, foldid)
  }
  
  fold_levels <- sort(unique(foldid))
  K <- length(fold_levels)
  
  ## IMPORTANT for build_predmat_mv: it expects foldid labels 1..max(foldid)
  ## We already renumbered foldid to 1..K, so max(foldid) == K.
  
  cv_by_rho <- vector("list", length(rho))
  names(cv_by_rho) <- paste0("rho_", rho)
  
  rho_cvmean <- numeric(length(rho))
  rho_cvse   <- numeric(length(rho))
  rho_lambda_min <- numeric(length(rho))
  rho_lambda_1se <- numeric(length(rho))
  rho_lambda_choice <- numeric(length(rho))
  
  ## -----------------------------
  ## outer loop over rho
  ## -----------------------------
  for (r in seq_along(rho)) {
    
    if (trace.it)
      cat(sprintf("Tuning rho = %.2f (%d/%d)\n", rho[r], r, length(rho)))
    
    ## get lambda path if not supplied
    if (is.null(lambda)) {
      tmp_fit <- mv_fit(
        x_list = x_list, y = y, family = family, rho = rho[r],
        weights = weights, offset = offset,
        penalty.factor = penalty.factor,
        trace.it = 0, ...
      )
      lambda_r <- tmp_fit$lambda
    } else {
      lambda_r <- lambda
    }
    
    # fit per fold (training-only), like cv.multiview()
    outlist <- vector("list", K)
    
    if (trace.it) cat("  Inner folds\n")
    
    for (i in seq_along(fold_levels)) {
      
      f <- fold_levels[i]
      is_test <- (foldid == f)
      
      x_train <- lapply(x_list, function(m) m[!is_test, , drop = FALSE])
      y_train <- if (is.matrix(y)) y[!is_test, , drop = FALSE] else y[!is_test]
      
      w_train <- if (is.null(weights)) NULL else weights[!is_test]
      off_train <- if (is.null(offset)) NULL else offset[!is_test]
      
      outlist[[f]] <- mv_fit(
        x_list = x_train, y = y_train, family = family, rho = rho[r],
        weights = w_train, offset = off_train, lambda = lambda_r,
        penalty.factor = penalty.factor,
        trace.it = 0, ...
      )
    }
    
    # build predmat exactly like cv.multiview() does (using outlist)
    predmat <- build_predmat_mv(
      outlist = outlist,
      lambda  = lambda_r,
      x_list  = x_list,
      offset  = offset,
      foldid  = foldid,
      alignment = alignment,
      family = family,
      type = "response"
    )
    
    # compute fold-wise metrics into cvm (lambda x fold) from predmat
    cvm <- matrix(NA_real_, nrow = length(lambda_r), ncol = K)
    for (i in seq_along(fold_levels)) {
      f <- fold_levels[i]
      test_idx <- which(foldid == f)
      y_test <- if (is.matrix(y)) y[test_idx, , drop = FALSE] else y[test_idx]
      preds_fold <- predmat[test_idx, , drop = FALSE]
      
      if (type.measure == "mse") {
        # y_test is a vector here for gaussian in your use; keep simple
        cvm[, i] <- colMeans((as.numeric(y_test) - preds_fold)^2)
      } else if (type.measure == "auc") {
        cvm[, i] <- apply(preds_fold, 2, function(p) {
          suppressMessages({
            rocobj <- pROC::roc(
              response  = y_test,
              predictor = p,
              levels    = c(0, 1),
              direction = "<",
              quiet     = TRUE
            )
            as.numeric(pROC::auc(rocobj))
          })
        })
      } else {
        stop("Unsupported type.measure in this implementation.")
      }
    }
    
    cvmean <- rowMeans(cvm, na.rm = TRUE)
    cvse   <- apply(cvm, 1, sd, na.rm = TRUE) / sqrt(K)
    
    idx_min <- which.min(cvmean)
    lambda_min <- lambda_r[idx_min]
    
    one_se_rule <- cvmean[idx_min] + cvse[idx_min]
    lambda_1se <- max(lambda_r[cvmean <= one_se_rule])
    
    lambda_choice <- if (s == "lambda.min") lambda_min else lambda_1se
    
    j <- which.min(abs(lambda_r - lambda_choice))
    rho_cvmean[r] <- cvmean[j]
    rho_cvse[r]   <- cvse[j]
    
    rho_lambda_min[r] <- lambda_min
    rho_lambda_1se[r] <- lambda_1se
    rho_lambda_choice[r] <- lambda_choice
    
    cv_by_rho[[r]] <- list(
      rho = rho[r],
      lambda = lambda_r,
      cvm = cvm,
      cvmean = cvmean,
      cvse = cvse,
      lambda.min = lambda_min,
      lambda.1se = lambda_1se,
      lambda.choice = lambda_choice,
      type.measure = type.measure,
      s = s,
      alignment = alignment,
      keep = keep,
      # keep semantics matching cv.multiview():
      fit.preval = if (keep) predmat else NULL,
      foldid = foldid
    )
  }
  
  ## -----------------------------
  ## select rho using rho_cvmean (based on user's s)
  ## -----------------------------
  best_r <- which.min(rho_cvmean)
  rho_choice <- rho[best_r]
  lambda_choice <- rho_lambda_choice[best_r]
  
  ## -----------------------------
  ## refit final model (do NOT pass keep)
  ## -----------------------------
  final_fit <- mv_fit(
    x_list = x_list, y = y, family = family, rho = rho_choice, lambda = lambda_choice,
    weights = weights, offset = offset,
    penalty.factor = penalty.factor,
    trace.it = 0, ...
  )
  
  out <- list(
    multiview.fit = final_fit,
    rho = rho,
    cv_by_rho = cv_by_rho,
    rho.cvmean = rho_cvmean,
    rho.cvse = rho_cvse,
    s = s,
    rho.choice = rho_choice,
    lambda.choice = lambda_choice,
    rho.lambda.min = rho_lambda_min,
    rho.lambda.1se = rho_lambda_1se,
    rho.lambda.choice = rho_lambda_choice,
    type.measure = type.measure,
    foldid = foldid,
    nfolds = length(unique(foldid)),
    keep = keep,
    alignment = alignment
  )
  
  class(out) <- "cv.multiview.revised"
  return(out)
}
