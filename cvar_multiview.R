cvar.multiview <- function(x_list, y,
                           family = gaussian(),
                           alpha = c(0, seq(0.05, 0.95, 0.05), 1),
                           rho = c(0, 0.1, 0.25, 0.5, 1, 5, 10),
                           lambda = NULL,
                           s = c("lambda.1se", "lambda.min"),
                           nfolds = 10,
                           foldid = NULL,
                           weights = NULL,
                           offset = NULL,
                           penalty.factor = NULL,
                           type.measure = NULL,
                           alignment = c("lambda", "fraction"),
                           keep = FALSE,
                           trace.it = 0,
                           seed = NULL,
                           verbose = FALSE,
                           ...) {

  mv_fit <- function(x_list, y, family, alpha, rho, lambda = NULL,
                     weights = NULL, offset = NULL,
                     penalty.factor = NULL, trace.it = 0, ...) {
    args <- list(
      x_list = x_list,
      y = y,
      family = family,
      alpha = alpha,
      rho = rho,
      weights = weights,
      offset = offset,
      lambda = lambda,
      trace.it = trace.it,
      ...
    )
    if (!is.null(penalty.factor)) {
      args$penalty.factor <- penalty.factor
    }

    fit_once <- function(extra = list()) {
      do.call(multiview::multiview, c(args, extra))
    }

    out <- tryCatch(fit_once(), error = function(e) e)
    if (!inherits(out, "error")) {
      return(out)
    }

    msg <- conditionMessage(out)
    if (grepl("cannot correct step size", msg, fixed = TRUE)) {
      out2 <- tryCatch(
        fit_once(extra = list(lambda.min.ratio = 0.05, nlambda = 50, maxit = 1e5)),
        error = function(e) e
      )
      if (!inherits(out2, "error")) {
        return(out2)
      }
    }

    stop(out)
  }

  make_stratified_foldid <- function(y01, K) {
    y01 <- as.integer(y01)
    if (!all(y01 %in% c(0L, 1L))) {
      stop("For binomial, y must be 0/1.")
    }
    idx0 <- which(y01 == 0L)
    idx1 <- which(y01 == 1L)

    f0 <- sample(rep(seq_len(K), length.out = length(idx0)))
    f1 <- sample(rep(seq_len(K), length.out = length(idx1)))

    fid <- integer(length(y01))
    fid[idx0] <- f0
    fid[idx1] <- f1
    fid
  }

  validate_binomial_folds <- function(y01, fid) {
    K <- max(fid)
    for (f in seq_len(K)) {
      is_test <- fid == f
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

  build_predmat_mv <- function(outlist, lambda, x_list, offset, foldid,
                               alignment = c("lambda", "fraction"),
                               family, type = "response", ...) {
    alignment <- match.arg(alignment)
    N <- nrow(x_list[[1L]])

    if (!is.null(offset)) {
      offset <- drop(offset)
    }

    predmat <- matrix(NA_real_, N, length(lambda))
    nlambda <- length(lambda)

    for (f in seq_len(max(foldid))) {
      which <- foldid == f
      fitobj <- outlist[[f]]
      x_sub_list <- lapply(x_list, function(x) x[which, , drop = FALSE])
      offset_sub <- if (is.null(offset)) NULL else offset[which]

      preds <- switch(
        alignment,
        fraction = predict(fitobj, newx = x_sub_list, newoffset = offset_sub,
                           type = type, ...),
        lambda = predict(fitobj, newx = x_sub_list, s = lambda,
                         newoffset = offset_sub, type = type, ...)
      )

      nlami <- min(ncol(preds), nlambda)
      predmat[which, seq_len(nlami)] <- preds[, seq_len(nlami), drop = FALSE]
      if (nlami < nlambda) {
        predmat[which, (nlami + 1L):nlambda] <- preds[, nlami]
      }
    }

    dimnames(predmat) <- list(rownames(x_list[[1L]]), paste0("s", seq_len(nlambda)))
    attr(predmat, "family") <- family
    predmat
  }

  score_fold_metric <- function(y_test, preds_fold, family_name, type.measure) {
    if (type.measure == "mse") {
      return(colMeans((as.numeric(y_test) - preds_fold)^2))
    }

    if (type.measure == "deviance") {
      if (family_name == "gaussian") {
        return(colMeans((as.numeric(y_test) - preds_fold)^2))
      }
      if (family_name == "binomial") {
        eps <- 1e-8
        yv <- as.numeric(y_test)
        return(apply(preds_fold, 2, function(p) {
          pp <- pmin(pmax(as.numeric(p), eps), 1 - eps)
          -2 * mean(yv * log(pp) + (1 - yv) * log(1 - pp))
        }))
      }
      stop("type.measure='deviance' currently implemented for gaussian/binomial only.")
    }

    if (type.measure == "auc") {
      return(apply(preds_fold, 2, function(p) {
        suppressMessages({
          rocobj <- pROC::roc(
            response = y_test,
            predictor = p,
            levels = c(0, 1),
            direction = "<",
            quiet = TRUE
          )
          as.numeric(pROC::auc(rocobj))
        })
      }))
    }

    if (type.measure == "class") {
      yv <- as.numeric(y_test)
      return(apply(preds_fold, 2, function(p) {
        mean((as.numeric(p) >= 0.5) != yv)
      }))
    }

    stop("Unsupported type.measure in cvar.multiview implementation.")
  }

  choose_best_index <- function(vals, metric) {
    if (metric %in% c("auc")) {
      which.max(vals)
    } else {
      which.min(vals)
    }
  }

  choose_1se <- function(lambda, cvmean, cvse, idx_best, metric) {
    if (metric %in% c("auc")) {
      threshold <- cvmean[idx_best] - cvse[idx_best]
      min(lambda[cvmean >= threshold])
    } else {
      threshold <- cvmean[idx_best] + cvse[idx_best]
      max(lambda[cvmean <= threshold])
    }
  }

  alignment <- match.arg(alignment)
  if (!is.null(seed)) {
    set.seed(seed)
  }

  N <- nrow(x_list[[1L]])
  if (any(vapply(x_list, function(m) nrow(m), integer(1)) != N)) {
    stop("All views in x_list must have the same number of rows.")
  }

  y <- drop(y)
  if (length(y) != N) {
    stop("Length of y must equal nrow(x_list[[1]]).")
  }

  if (is.null(weights)) {
    weights <- rep(1, N)
  }
  if (length(weights) != N) {
    stop("Length of weights must equal N.")
  }
  if (!is.null(offset) && length(offset) != N) {
    stop("Length of offset must equal N when provided.")
  }

  if (is.null(type.measure)) {
    if (family$family == "gaussian") {
      type.measure <- "mse"
    } else if (family$family == "binomial") {
      type.measure <- "auc"
    } else {
      type.measure <- "deviance"
    }
  }

  if (length(s) != 1) {
    s <- s[1]
  }
  if (is.character(s) && identical(s, "lambda.lse")) {
    s <- "lambda.1se"
  }
  s <- match.arg(s, choices = c("lambda.1se", "lambda.min"))

  alpha <- sort(unique(alpha))
  rho <- sort(unique(rho))

  if (is.null(foldid)) {
    if (family$family == "binomial") {
      y01 <- as.integer(as.numeric(y) > 0)
      foldid <- make_stratified_foldid(y01, nfolds)
    } else {
      foldid <- sample(rep(seq_len(nfolds), length.out = N))
    }
  } else {
    if (length(foldid) != N) {
      stop("foldid must have length N (nrow of data).")
    }
    foldid <- as.integer(foldid)
    used <- sort(unique(foldid))
    fold_map <- setNames(seq_along(used), used)
    foldid <- as.integer(fold_map[as.character(foldid)])
    nfolds <- length(unique(foldid))
    if (nfolds < 2) {
      stop("foldid must contain at least 2 unique folds.")
    }
  }

  if (family$family == "binomial") {
    validate_binomial_folds(as.integer(as.numeric(y) > 0), foldid)
  }

  fold_levels <- sort(unique(foldid))
  K <- length(fold_levels)

  grid <- expand.grid(alpha = alpha, rho = rho)
  combo_results <- vector("list", nrow(grid))

  score_choice <- numeric(nrow(grid))
  lambda_choice_vec <- numeric(nrow(grid))
  lambda_min_vec <- numeric(nrow(grid))
  lambda_1se_vec <- numeric(nrow(grid))

  for (g in seq_len(nrow(grid))) {
    a <- grid$alpha[g]
    r <- grid$rho[g]

    if (trace.it || verbose) {
      cat(sprintf("Tuning alpha = %.3f, rho = %.3f (%d/%d)\n", a, r, g, nrow(grid)))
    }

    if (is.null(lambda)) {
      tmp_fit <- mv_fit(
        x_list = x_list,
        y = y,
        family = family,
        alpha = a,
        rho = r,
        weights = weights,
        offset = offset,
        penalty.factor = penalty.factor,
        trace.it = 0,
        ...
      )
      lambda_g <- tmp_fit$lambda
    } else {
      lambda_g <- lambda
    }

    outlist <- vector("list", K)
    for (i in seq_along(fold_levels)) {
      f <- fold_levels[i]
      is_test <- foldid == f

      x_train <- lapply(x_list, function(m) m[!is_test, , drop = FALSE])
      y_train <- if (is.matrix(y)) y[!is_test, , drop = FALSE] else y[!is_test]
      w_train <- weights[!is_test]
      off_train <- if (is.null(offset)) NULL else offset[!is_test]

      outlist[[f]] <- mv_fit(
        x_list = x_train,
        y = y_train,
        family = family,
        alpha = a,
        rho = r,
        lambda = lambda_g,
        weights = w_train,
        offset = off_train,
        penalty.factor = penalty.factor,
        trace.it = 0,
        ...
      )
    }

    predmat <- build_predmat_mv(
      outlist = outlist,
      lambda = lambda_g,
      x_list = x_list,
      offset = offset,
      foldid = foldid,
      alignment = alignment,
      family = family,
      type = "response"
    )

    cvm <- matrix(NA_real_, nrow = length(lambda_g), ncol = K)
    for (i in seq_along(fold_levels)) {
      f <- fold_levels[i]
      test_idx <- which(foldid == f)
      y_test <- if (is.matrix(y)) y[test_idx, , drop = FALSE] else y[test_idx]
      preds_fold <- predmat[test_idx, , drop = FALSE]
      cvm[, i] <- score_fold_metric(y_test, preds_fold, family$family, type.measure)
    }

    cvmean <- rowMeans(cvm, na.rm = TRUE)
    cvse <- apply(cvm, 1, sd, na.rm = TRUE) / sqrt(K)

    idx_best <- choose_best_index(cvmean, type.measure)
    lambda_min <- lambda_g[idx_best]
    lambda_1se <- choose_1se(lambda_g, cvmean, cvse, idx_best, type.measure)
    lambda_choice <- if (s == "lambda.min") lambda_min else lambda_1se

    choice_idx <- which.min(abs(lambda_g - lambda_choice))
    score_choice[g] <- cvmean[choice_idx]
    lambda_choice_vec[g] <- lambda_choice
    lambda_min_vec[g] <- lambda_min
    lambda_1se_vec[g] <- lambda_1se

    combo_results[[g]] <- list(
      alpha = a,
      rho = r,
      lambda = lambda_g,
      cvm = cvm,
      cvmean = cvmean,
      cvse = cvse,
      lambda.min = lambda_min,
      lambda.1se = lambda_1se,
      lambda.choice = lambda_choice,
      fit.preval = if (keep) predmat else NULL,
      foldid = foldid,
      type.measure = type.measure,
      s = s,
      alignment = alignment
    )
  }

  best_idx <- choose_best_index(score_choice, type.measure)
  alpha_choice <- grid$alpha[best_idx]
  rho_choice <- grid$rho[best_idx]
  lambda_choice <- lambda_choice_vec[best_idx]

  final_fit <- mv_fit(
    x_list = x_list,
    y = y,
    family = family,
    alpha = alpha_choice,
    rho = rho_choice,
    lambda = lambda_choice,
    weights = weights,
    offset = offset,
    penalty.factor = penalty.factor,
    trace.it = 0,
    ...
  )

  selected <- combo_results[[best_idx]]
  idx_min <- which.min(abs(selected$lambda - selected$lambda.min))
  idx_1se <- which.min(abs(selected$lambda - selected$lambda.1se))

  idx_alpha_choice <- which(grid$alpha == alpha_choice)
  idx_alpha_choice <- idx_alpha_choice[order(grid$rho[idx_alpha_choice])]
  cv_by_rho <- lapply(idx_alpha_choice, function(i) {
    ri <- combo_results[[i]]
    list(
      rho = ri$rho,
      lambda = ri$lambda,
      cvm = ri$cvm,
      cvmean = ri$cvmean,
      cvse = ri$cvse,
      lambda.min = ri$lambda.min,
      lambda.1se = ri$lambda.1se,
      lambda.choice = ri$lambda.choice,
      type.measure = ri$type.measure,
      s = ri$s,
      alignment = ri$alignment,
      keep = keep,
      fit.preval = ri$fit.preval,
      foldid = ri$foldid
    )
  })

  rho_seq <- vapply(cv_by_rho, function(z) z$rho, numeric(1))
  rho_lambda_min <- vapply(cv_by_rho, function(z) z$lambda.min, numeric(1))
  rho_lambda_1se <- vapply(cv_by_rho, function(z) z$lambda.1se, numeric(1))
  rho_lambda_choice <- vapply(cv_by_rho, function(z) z$lambda.choice, numeric(1))
  rho_cvmean <- vapply(cv_by_rho, function(z) {
    j <- which.min(abs(z$lambda - z$lambda.choice))
    z$cvmean[j]
  }, numeric(1))
  rho_cvse <- vapply(cv_by_rho, function(z) {
    j <- which.min(abs(z$lambda - z$lambda.choice))
    z$cvse[j]
  }, numeric(1))

  out <- list(
    call = match.call(),
    multiview.fit = final_fit,
    lambda = selected$lambda,
    cvm = selected$cvmean,
    cvsd = selected$cvse,
    cvup = selected$cvmean + selected$cvse,
    cvlo = selected$cvmean - selected$cvse,
    lambda.min = selected$lambda.min,
    lambda.1se = selected$lambda.1se,
    index = matrix(c(idx_min, idx_1se), ncol = 1,
                   dimnames = list(c("min", "1se"), "Lambda")),
    fit.preval = if (keep) selected$fit.preval else NULL,
    foldid = if (keep) foldid else NULL,
    type.measure = type.measure,
    s = s,
    alpha = alpha,
    rho = rho_seq,
    cv_by_rho = cv_by_rho,
    rho.cvmean = rho_cvmean,
    rho.cvse = rho_cvse,
    rho.lambda.min = rho_lambda_min,
    rho.lambda.1se = rho_lambda_1se,
    rho.lambda.choice = rho_lambda_choice,
    alpha.choice = alpha_choice,
    rho.choice = rho_choice,
    lambda.choice = lambda_choice,
    alpha.rho.grid = grid,
    combo.results = combo_results,
    score.choice = score_choice,
    alpha.rho.choice.score = score_choice[best_idx],
    alpha.rho.lambda.min = lambda_min_vec,
    alpha.rho.lambda.1se = lambda_1se_vec,
    alpha.rho.lambda.choice = lambda_choice_vec,
    keep = keep,
    alignment = alignment,
    nfolds = K,
    family = family
  )

  class(out) <- c("cvar.multiview", "cv.multiview", "cv.glmnet")
  out
}
