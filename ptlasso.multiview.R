ptLasso.multiview1.0 <- function(
    x, y,
    alpha = 0.5,
    family = c("gaussian", "binomial"),
    type.measure = c("default", "mse", "auc", "deviance"),
    overall.lambda = c("lambda.1se", "lambda.min"),
    ind.lambda = c("lambda.1se", "lambda.min"), # NEW
    pre.lambda = c("lambda.1se", "lambda.min"), # NEW
    foldid = NULL, # added
    nfolds = 10,
    verbose = FALSE,
    fitoverall = NULL,
    fitind = NULL,
    penalty.factor = NULL, # added
    group.intercepts = TRUE, # added
    en.alpha = 1, # added
    ...
) {

  # =====================================================
  # STEP 1: Match arguments and defaults
  # =====================================================
  this.call <- match.call()

  family <- match.arg(family)

  type.measure <- match.arg(type.measure)
  if (type.measure == "default") {
    type.measure <- if (family == "gaussian") "mse" else "deviance"
  }

  overall.lambda = match.arg(overall.lambda, c("lambda.1se", "lambda.min"))
  ind.lambda = match.arg(ind.lambda, c("lambda.1se", "lambda.min"))
  pre.lambda = match.arg(pre.lambda, c("lambda.1se", "lambda.min"))

  if (alpha < 0 || alpha > 1)
    stop("alpha must be between 0 and 1.")

  # map character family → glmnet family object
  family_fn <- switch(
    family,
    gaussian = gaussian,
    binomial = binomial
  )

  # =====================================================
  # STEP 2: Input validation (multiview only)
  # =====================================================
  if (!is.list(x) || !is.list(x[[1]]))
    stop("x must be a list of studies, each a list of views (matrices).")

  if (!is.list(y))
    stop("y must be a list with one response vector per study.")

  k <- length(x)

  if (length(y) != k)
    stop("Length of y must match number of studies in x.")

  view_names <- names(x[[1]])
  if (is.null(view_names))
    stop("Each study in x must be a *named* list of views.")

  for (kk in seq_len(k)) {
    if (!identical(names(x[[kk]]), view_names))
      stop("All studies must share identical view names and ordering.")
  }

  # =============================================================
  # STEP 3: Dimension consistency checks + total # of variables p
  # =============================================================
  # - Stops ONLY if inconsistency is found
  # - Computes p (total columns across views) as a single number
  #   (assumes p is the same across studies if checks pass)
  # ------------------------------------------------------------
  .ptlasso_precheck_and_p <- function(x, check_colnames = TRUE) {
    if (!is.list(x) || length(x) == 0) stop("`x` must be a non-empty list of studies (each a named list of view matrices).")

    studies <- names(x); if (is.null(studies) || any(studies == "")) studies <- paste0("Study_", seq_along(x))
    view_ref <- names(x[[1]])
    if (is.null(view_ref) || any(view_ref == "")) stop(sprintf("'%s' must be a named list of views.", studies[1]))

    # reference dims + (optional) feature names from first study
    ref_dims <- vapply(view_ref, function(v) {
      m <- x[[1]][[v]]
      if (!is.matrix(m)) stop(sprintf("'%s' -> '%s' must be a matrix.", studies[1], v))
      c(n = nrow(m), p = ncol(m))
    }, numeric(2))
    ref_cols <- if (check_colnames) lapply(view_ref, function(v) colnames(x[[1]][[v]])) else NULL
    names(ref_cols) <- view_ref

    # check all other studies against reference; stop at first inconsistency
    for (i in seq_along(x)) {
      s <- studies[i]
      xi <- x[[i]]
      if (!is.list(xi) || length(xi) == 0) stop(sprintf("'%s' is not a non-empty list of views.", s))

      v_i <- names(xi)
      if (is.null(v_i) || any(v_i == "")) stop(sprintf("'%s' must be a named list of views.", s))
      if (!setequal(v_i, view_ref)) {
        stop(sprintf(
          "View set mismatch in %s.\nMissing: %s\nExtra: %s",
          s,
          paste(setdiff(view_ref, v_i), collapse = ", "),
          paste(setdiff(v_i, view_ref), collapse = ", ")
        ))
      }

      for (v in view_ref) {
        m <- xi[[v]]
        if (!is.matrix(m)) stop(sprintf("'%s' -> '%s' must be a matrix.", s, v))

        if (nrow(m) != ref_dims["n", v] || ncol(m) != ref_dims["p", v]) {
          stop(sprintf(
            "Dimension mismatch in %s -> %s. Expected (n=%d, p=%d), got (n=%d, p=%d).",
            s, v, ref_dims["n", v], ref_dims["p", v], nrow(m), ncol(m)
          ))
        }

        if (check_colnames) {
          cc <- colnames(m); rc <- ref_cols[[v]]
          if (!identical(cc, rc)) {
            stop(sprintf("Colname mismatch in %s -> %s (feature ordering/content differs from %s).", s, v, studies[1]))
          }
        }
      }
    }

    # minimal p as a single number (from reference study; valid if checks passed)
    p <- sum(ref_dims["p", ])
    return(p)
  }

  p = .ptlasso_precheck_and_p(x, check_colnames = T)

  # =====================================================
  # STEP 4: Validate pretrained model inputs (if supplied)
  # =====================================================
  if (!is.null(fitoverall)) {
    if (!inherits(fitoverall, "cv.multiview") && !inherits(fitoverall, "cv.multiview.revised"))
      stop("fitoverall must be a cv.multiview/cvar.multiview/cv.multiview.revised object.")
  }

  if (!is.null(fitind)) {
    if (length(fitind) != k)
      stop("fitind must contain one model per study.")
    if (!all(sapply(fitind, function(z) inherits(z, "cv.multiview") || inherits(z, "cv.multiview.revised"))))
      stop("All elements of fitind must be cv.multiview/cvar.multiview/cv.multiview.revised objects.")
  }

  # =====================================================
  # STEP 5: dimension calculation
  # =====================================================
  this.call$type.measure <- type.measure

  # k groups = number of studies
  group.levels <- names(y)

  if (is.null(group.levels) || any(group.levels == "")) {
    group.levels <- paste0("Group_", seq_along(x))
    names(x) <- group.levels
    if (!is.null(y)) names(y) <- group.levels
  }

  # group sizes
  n_by_group <- sapply(y, length)
  N_all <- sum(n_by_group)

  # groups_all: length N_all, in the same stacking order as y_all <- unlist(y)
  groups_all <- factor(rep(group.levels, times = n_by_group), levels = group.levels)

  # =====================================================
  # STEP 6: fold construction
  # =====================================================
  # -----------------------------
  # foldid_all: for OVERALL model
  # -----------------------------
  if (is.null(foldid)) {
    # Balanced folds across studies (so each fold gets data from each study)
    foldid_all <- integer(N_all)

    # For each group, assign folds within that group, then paste into foldid_all
    start <- 1L
    for (g in seq_len(k)) {
      n_g <- n_by_group[g]
      # random fold assignment within group
      fold_g <- sample(rep(seq_len(nfolds), length.out = n_g))
      foldid_all[start:(start + n_g - 1L)] <- fold_g
      start <- start + n_g
    }
    # renumber to ensure folds are 1..F actually present
    used <- sort(unique(foldid_all))
    fold_map <- setNames(seq_along(used), used)
    foldid_all <- as.integer(fold_map[as.character(foldid_all)])

    nfolds_all <- length(unique(foldid_all))
  } else {
    # User supplied foldid must match stacked length
    if (length(foldid) != N_all)
      stop("Provided foldid must have length equal to total N across studies (sum(sapply(y, length))).")
    foldid_all <- as.integer(foldid)
    nfolds_all <- length(unique(foldid_all))
  }

  # ---------------------------------------------
  # foldid2: for WITHIN-STUDY individual/pretrained
  # ---------------------------------------------
  foldid2 <- vector("list", k)
  names(foldid2) <- group.levels

  for (g in seq_len(k)) {
    n_g <- n_by_group[g]

    # Don’t force the same number of folds if n_g is small
    nfolds_g <- min(nfolds, n_g)

    fold_g <- sample(rep(seq_len(nfolds_g), length.out = n_g))

    # Ensure labels are 1..nfolds_g
    used <- sort(unique(fold_g))
    fold_map <- setNames(seq_along(used), used)
    foldid2[[g]] <- as.integer(fold_map[as.character(fold_g)])
  }

  # =====================================================
  # STEP 7: Fold-wise, group-wise baseline offset
  # - compute a baseline offset per obs using foldid_all + groups_all
  # =====================================================

  if (is.null(penalty.factor))
    penalty.factor <- rep(1, p)

  overall.pf <- penalty.factor

  # helper: safe logit with smoothing for binomial
  .safe_logit <- function(p) log(p / (1 - p))

  .compute_baseline_offset <- function(y_all, groups_all, foldid_all, family) {
    N_all <- length(y_all)
    off <- numeric(N_all)

    folds <- sort(unique(foldid_all))
    grps  <- levels(groups_all)

    if (family$family == "gaussian") {

      for (f in folds) {
        train <- foldid_all != f
        # compute training mean within each group
        for (g in grps) {
          idx <- which(groups_all == g & !train)          # test indices in this fold+group
          if (length(idx) == 0) next
          y_tr_g <- y_all[train & groups_all == g]
          off[idx] <- mean(y_tr_g)
        }
      }
      return(off)

    } else if (family$family == "binomial") {

      # Jeffreys prior smoothing avoids logit(0) / logit(1)
      a <- 0.5; b <- 0.5

      for (f in folds) {
        train <- foldid_all != f
        for (g in grps) {
          idx <- which(groups_all == g & !train)
          if (length(idx) == 0) next
          y_tr_g <- y_all[train & groups_all == g]

          n1 <- sum(y_tr_g == 1)
          n0 <- sum(y_tr_g == 0)
          p_hat <- (n1 + a) / (n1 + n0 + a + b)

          off[idx] <- .safe_logit(p_hat)
        }
      }
      return(off)

    } else {
      stop("Baseline offset only implemented for gaussian and binomial.")
    }
  }

  # baseline offset is only meaningful if you want group-specific baselines
  # (you can still allow group.intercepts=FALSE to skip it)
  baseline_offset_all <- NULL
  if (group.intercepts) {
    y_all <- unlist(y)
    baseline_offset_all <- .compute_baseline_offset(
      y_all       = y_all,
      groups_all  = groups_all,
      foldid_all  = foldid_all,
      family      = family_fn()
    )
  }

  # =====================================================
  # STEP 8: Fit overall multiview model (cvar.multiview)
  # - No groupIntercept view
  # - Use baseline_offset_all as offset argument
  # =====================================================

  if (is.null(fitoverall)) {

    if (verbose)
      cat("Fitting overall multiview model (with baseline offset)\n")

    # 1) stack each view across studies in SAME order as y_all <- unlist(y)
    x_list_all <- lapply(view_names, function(v) {
      do.call(rbind, lapply(x, function(study) study[[v]]))
    })
    names(x_list_all) <- view_names

    y_all <- unlist(y)

    fitoverall <- cvar.multiview(
      x_list         = x_list_all,
      y              = y_all,
      family         = family_fn(),
      alpha          = 1,
      s = overall.lambda, # NEW user-defined
      type.measure   = type.measure,
      foldid         = foldid_all,
      nfolds         = nfolds_all,
      offset         = baseline_offset_all, # NEW: NULL if group.intercepts=FALSE
      penalty.factor = overall.pf,
      keep = TRUE,
      ...
    )
  }

  # =====================================================
  # STEP 9: Extract offsets, global coefficients, support
  # =====================================================
  # OFFSET
  # get the overall fit model
  fitoverall_fit <- fitoverall$multiview.fit
  lamhat = fitoverall$lambda.choice
  rho_overall = fitoverall$rho.choice #NEW: store cv selected rho

  # OFFSET
  r_idx <- which(fitoverall$rho == fitoverall$rho.choice)

  rho_obj <- fitoverall$cv_by_rho[[r_idx]]
  lam_idx <- which.min(abs(rho_obj$lambda - fitoverall$lambda.choice))

  preval_all <- rho_obj$fit.preval[, lam_idx, drop = FALSE]

  lens <- sapply(y, length)
  preval.offset <- split(preval_all, rep(seq_along(lens), times = lens))
  names(preval.offset) <- names(y)

  # COEFFICIENTS
  coef_vec <- as.numeric(coef(fitoverall_fit, s = lamhat))
  # First element is global intercept
  coef_no_int <- coef_vec[-1]
  # SUPPORT set from real multiview features only
  supall <- which(coef_no_int != 0)

  # =====================================================
  # STEP 10: Fit individual multiview models
  # =====================================================
  if (is.null(fitind)) {
    if (verbose)
      cat("Fitting individual multiview models\n")

    fitind <- vector("list", k)

    for (kk in seq_len(k)) {
      if (verbose)
        cat(sprintf("  Study %d / %d\n", kk, k))

      fitind[[kk]] <- cvar.multiview(
        x_list = x[[kk]],
        y      = y[[kk]],
        alpha  = 1,
        foldid = foldid2[[kk]],
        s = ind.lambda, # NEW user-defined
        type.measure = type.measure,
        family = family_fn(),
        keep = TRUE,
        ...
      )
    }
  }

  # =====================================================
  # STEP 11: Fit pretrained multiview models
  # =====================================================
  if (verbose)
    cat("Fitting pretrained multiview models\n")

  fitpre <- vector("list", k)

  if (alpha == 1) {
    fitpre <- fitind
  } else {
    for (kk in seq_len(k)) {
      if (verbose)
        cat(sprintf("  Pretrained model %d / %d\n", kk, k))

      # extract study-specific offset for pretrain model
      offset = (1 - alpha) * preval.offset[[kk]]
      # compute penalty factor for pretrain model (determined by overall model, thus consistent across models)
      alpha_eff <- max(alpha, 1e-9)
      fac = rep(1 / alpha_eff, p)
      fac[supall] = 1
      pf = penalty.factor * fac

      if ((alpha == 0) & (length(supall) == 0)) {
        almost.zero = 1e-09
        fac = rep(1 / almost.zero, p)
        pf = penalty.factor * fac
      }

      fitpre[[kk]] <- cvar.multiview(
        x_list         = x[[kk]],
        y              = y[[kk]],
        alpha          = 1,
        foldid = foldid2[[kk]],
        s = pre.lambda, # NEW user-defined
        type.measure = type.measure,
        family         =  family_fn(),
        offset         = offset,
        penalty.factor = pf,
        ...
      )
    }
  }

  # =====================================================
  # STEP 12: Output
  # =====================================================
  out <- list(
    call            = this.call,
    k               = k, # number of groups
    N_all           = N_all, # total number of observations in x
    features_all    = p, # total number of features from all views
    alpha           = alpha,
    family          = family,
    support.vars    = supall,
    fitoverall      = fitoverall,
    fitind          = fitind,
    fitpre          = fitpre
  )

  class(out) <- "ptLasso.multiview"
  return(out)
}

