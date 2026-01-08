general_garch_fitting <- function(rt, model_type, orders = c(1, 1),
                        cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald"),
                        drange = c(0, 1),
                        meanspec = mean_spec(), Drange = c(0, 1),
                        nonparspec = locpol_spec(), use_nonpar = FALSE,
                        n_test = 0,
                        start_pars = NULL, LB = NULL, UB = NULL, control = list(),
                        control_nonpar = list(), mean_after_nonpar = FALSE,
                        parallel = TRUE, ncores = max(1, future::availableCores() - 1),
                        trunc = "none", presample = 50, Prange = c(1, 5), skip_vcov = FALSE) {


  p <- orders[[1]]
  q <- orders[[2]]

  cond_d <- match.arg(cond_dist)
  lm <- (model_type %in% c("figarch", "figjrgarch", "fiaparch", "fitgarch"))

  if (lm && !all(orders == c(1, 1))) {
    orders <- c(1, 1)
    message("Only orders p=1 with q=1 are currently supported for FIGARCH, FITGARCH, FIGJR-GARCH and FIAPARCH models.")
  }

  if (use_nonpar) {

    nonpar_out <- run_nonpar(nonparspec = nonparspec,
               control_nonpar = control_nonpar, rt = rt, n_test = n_test,
               mean_after_nonpar = mean_after_nonpar, lm = lm)
    rt <- nonpar_out$rt
    meanspec <- nonpar_out$meanspec
    nonpar_result <- nonpar_out$nonpar_result
    n_test <- nonpar_out$n_test

  }

  initial_split(rt, n_test)

  n <- length(rt_core)
  if (Prange[[1]] < 1) {Prange[[1]] <- 1}
  trunc_i <- trunc

  # Settings for ARMA / FARIMA part in the mean; makes objects available in this
  # environment
  identify_arma_settings(meanspec = meanspec, Drange = Drange)

  if ((lm || lm_arma) && is.character(trunc) && trunc == "none") {
    trunc <- n + presample - 1
  } else if ((lm || lm_arma) && is.numeric(trunc) && trunc >= n + presample) {
    trunc <- n + presample - 1
    message("trunc reduced to n + ", presample, " - 1")
  }

  dfun1 <- fun1_selector(cond_d)
  dfun2 <- fun2_selector(cond_d)

  skew_par <- lookup_table[["skew_par"]][[cond_d]]
  extra_par <- lookup_table[["extra_par"]][[cond_d]]

  sigt_fun_adj <- adj_sigt_fun(
    p_ar = p_ar, q_ma = q_ma, pg0 = pg0, qg0 = qg0,
    lm_arma = lm_arma, incl_mean = incl_mean,
    n_arma_pars = n_arma_pars, grab_fun = arma_grab_fun,
    arma_fun = cmean_fun,
    model_type = model_type
  )

  creator_fun <- switch(
    model_type,
    "garch" = goal_fun_creator_garch,
    "gjrgarch" = goal_fun_creator_gjrgarch,
    "tgarch" = goal_fun_creator_tgarch,
    "aparch" = goal_fun_creator_aparch,
    "figarch" = goal_fun_creator_figarch,
    "figjrgarch" = goal_fun_creator_figjrgarch,
    "fitgarch" = goal_fun_creator_fitgarch,
    "fiaparch" = goal_fun_creator_fiaparch
  )

  base_args <- list(
    p = p,
    q = q,
    p_ar = p_ar,
    q_ma = q_ma,
    incl_mean = incl_mean,
    extra_par = extra_par,
    skew_par = skew_par,
    sig_fun = sigt_fun_adj,
    trunc = trunc,
    presample = presample
  )

  extra_args <- switch(
    model_type,
    "garch" = list(sig2_init = 1),
    "gjrgarch" = {
      delta0 <- 2    # use delta0 = 2
      gamma0 <- 0  # use gamma0 = 0   as suitable initial guesses
      sigd_init <- stats::sd(rt_core)^delta0
      rt_dm <- rt_core - mean(rt_core)
      etransf_init <- mean(abs(rt_dm) - gamma0 * rt_dm)^delta0
      list(sigd_init = sigd_init, etransf_init = etransf_init)
    },
    "tgarch" = {
      delta0 <- 1    # use delta0 = 1
      gamma0 <- 0  # use gamma0 = 0   as suitable initial guesses
      sigd_init <- stats::sd(rt_core)^delta0
      rt_dm <- rt_core - mean(rt_core)
      etransf_init <- mean(abs(rt_dm) - gamma0 * rt_dm)^delta0
      list(sigd_init = sigd_init, etransf_init = etransf_init)
    },
    "aparch" = {
      delta0 <- 2    # use delta0 = 2
      gamma0 <- 0  # use gamma0 = 0   as suitable initial guesses
      sigd_init <- stats::sd(rt_core)^delta0
      rt_dm <- rt_core - mean(rt_core)
      etransf_init <- mean(abs(rt_dm) - gamma0 * rt_dm)^delta0
      list(sigd_init = sigd_init, etransf_init = etransf_init)
    },
    "figarch" = {
      list(presample_val = var(rt_core))
    },
    "figjrgarch" = {
      est_s <- suppressWarnings(gjrgarch(
        rt = rt,
        orders = orders,
        cond_dist = cond_dist,
        meanspec = meanspec,
        use_nonpar = use_nonpar,
        nonparspec = nonparspec,
        parallel = parallel,
        ncores = ncores,
        Prange = Prange,
        n_test = n_test
      ))
      delta0 <- 2
      gamma0 <- est_s@pars[["gamma1"]]
      rt_dm <- rt_core - mean(rt_core)
      list(presample_val = mean(abs(rt_dm) - gamma0 * rt_dm)^delta0)
    },
    "fitgarch" = {
      est_s <- suppressWarnings(tgarch(
        rt = rt,
        orders = orders,
        cond_dist = cond_dist,
        meanspec = meanspec,
        use_nonpar = use_nonpar,
        nonparspec = nonparspec,
        parallel = parallel,
        ncores = ncores,
        Prange = Prange,
        n_test = n_test
      ))
      delta0 <- 1
      gamma0 <- est_s@pars[["gamma1"]]
      rt_dm <- rt_core - mean(rt_core)
      list(presample_val = mean(abs(rt_dm) - gamma0 * rt_dm)^delta0)
    },
    "fiaparch" = {
      est_s <- suppressWarnings(aparch(
        rt = rt,
        orders = orders,
        cond_dist = cond_dist,
        meanspec = meanspec,
        use_nonpar = use_nonpar,
        nonparspec = nonparspec,
        parallel = parallel,
        ncores = ncores,
        Prange = Prange,
        n_test = n_test
      ))
      delta0 <- est_s@pars[["delta"]]
      gamma0 <- est_s@pars[["gamma1"]]
      rt_dm <- rt_core - mean(rt_core)
      list(presample_val = mean(abs(rt_dm) - gamma0 * rt_dm)^delta0)
    }
  )

  goal_fun_s <- function(theta, rt, dfun1, dfun2) {

    do.call(
      what = creator_fun,
      args = c(
        list(x = rt, theta = theta, pdf_fun2 = dfun2),
        base_args,
        extra_args
      )
    )

  }

  # Get cond. distribution restrictions and starting values
  # as objects sval, lval and uval in the environment
  cond_d_vals <- cond_d_par_restr(cond_d)

  # Make further starting parameter values and restrictions available
  # using internally saved table function
  pars_and_restr <- start_pars_and_restr[[model_type]](drange, p, q)
  list2env(pars_and_restr, envir = environment())
  rm(pars_and_restr)


  if (is.null(start_pars) || is.null(LB) || is.null(UB)) {
    mean_rt <- list(NULL, mean(rt_core))[[incl_mean + 1]]
    min_rt <- list(NULL, suppressWarnings(0.2 * min(rt_core - mean_rt) + mean_rt))[[incl_mean + 1]]
    max_rt <- list(NULL, suppressWarnings(0.2 * max(rt_core - mean_rt) + mean_rt))[[incl_mean + 1]]
    skew_start <- list(NULL, 0.98)[[skew_par + 1]]
    skew_low <- list(NULL, 1e-15)[[skew_par + 1]]
    skew_up <- list(NULL, Inf)[[skew_par + 1]]

    set_start_LB_UB(
      mean_val = mean_rt,
      ar_start = ar_start,
      ma_start = ma_start,
      D_start = D_start,
      omega_val = omega_start,
      phi_start = phi_start,
      psi_start = psi_start,
      add_pars_start = add_pars_start,
      d_start = d_spar,
      sval = sval,
      skew_start = skew_start,
      mean_low = min_rt,
      ar_low = ar_lb,
      ma_low = ma_lb,
      D_low = D_low,
      omega_low = omega_low,
      phi_low = phi_low,
      psi_low = psi_low,
      add_low = add_low,
      d_low = d_low,
      lval = lval, skew_low = skew_low,
      mean_up = max_rt,
      ar_up = ar_ub, ma_up = ma_ub,
      D_up = D_up,
      omega_up = omega_up,
      phi_up = phi_up, psi_up = psi_up,
      add_up = add_up,
      d_up = d_up,
      uval = uval, skew_up = skew_up,
      start_pars = start_pars, LB = LB, UB = UB
    )


  }

  constr_mean <- if (pg0) {

    low_ar <- 1 + incl_mean
    up_ar <- low_ar + p_ar - 1

    ineqLB_mean <- rep(1e-6, p_ar)
    ineqUB_mean <- rep(Inf, p_ar)

    # Stationarity constraint
    function(theta) {
      ar <- c(1, -theta[low_ar:up_ar])
      abs(polyroot(ar)) - 1
    }

  } else {

    ineqLB_mean <- NULL
    ineqUB_mean <- NULL

    # Not required if no AR-part
    function(theta) {
      NULL
    }

  }

  constr_vol <- if (lm) {

    low_phi <- 2 + n_arma_pars
    up_phi <- low_phi
    low_beta <- low_phi + 1
    up_beta <- low_beta

    add_s <- switch(
      model_type,
      "figarch" = 1,
      "figjrgarch" = 2,
      "fitgarch" = 2,
      "fiaparch" = 3
    )

    d_idx <- up_beta + add_s


    ineqLB_vol <- c(rep(0, 50), 1e-6)
    ineqUB_vol <- rep(Inf, 51)

    function(theta) {

      phi <- theta[low_phi:up_phi]
      beta <- theta[low_beta:up_beta]
      d <- theta[[d_idx]]

      coef_inf <- ar_infty(ar = phi, ma = -beta, d = d, max_i = 50)[-1]

      c2 <- abs(polyroot(c(1, -beta))) - 1
      c(coef_inf, c2)

    }
  } else {

      if (model_type == "garch") {
        low_phi <- 2 + n_arma_pars
        up_phi <- low_phi - 1 + p
        low_beta <- low_phi + 1
        up_beta <- low_beta - 1 + q

        ineqLB_vol <- c(0)
        ineqUB_vol <- 1 - 1e-6

        function(theta) {
          phi <- theta[low_phi:up_phi]
          beta <- theta[low_beta:up_beta]

          sum(phi) + sum(beta)
        }
      } else if (model_type %in% c("gjrgarch", "tgarch", "aparch")) {
        ineqLB_vol <- NULL
        ineqUB_vol <- NULL

        function(theta) {
          NULL
        }
      }

    }

  # Combined constraint bounds
  ineqLB <- c(ineqLB_mean, ineqLB_vol)
  ineqUB <- c(ineqUB_mean, ineqUB_vol)

  # Combined constraint function
  constr_fun <- if (lm || pg0) {
    function(theta, rt) {

      constr1 <- constr_mean(theta)
      constr2 <- constr_vol(theta)

      c(constr1, constr2)

    }

  } else {

    NULL

  }

  if (is.null(control$trace)) {control$trace <- FALSE}

  if (cond_d %in% c("ald", "sald")) {

    # Makes "P_sel" and "result" available
    ald_fit(
      Prange = Prange,
      parallel = parallel,
      ncores = ncores,
      dfun1 = dfun1,
      dfun2 = dfun2,
      goal_fun_s = goal_fun_s,
      start_pars = start_pars,
      LB = LB, UB = UB,
      ineqfun = constr_fun,
      ineqLB = ineqLB,
      ineqUB = ineqUB,
      rt_core = rt_core,
      control = control
    )

  } else {

    P_sel <- NULL

    goal_fun <- function(theta, rt) {
      goal_fun_s(theta = theta, rt = rt, dfun1 = NULL, dfun2 = dfun2)
    }

    result <- tryCatch(
      expr = {
        suppressWarnings(Rsolnp::solnp(
          pars = start_pars,
          fun = goal_fun,
          LB = LB,
          UB = UB,
          control = control,
          ineqfun = constr_fun,
          ineqLB = ineqLB,
          ineqUB = ineqUB,
          rt = rt_core
        ))
      },
      error = function(e1) {
        stop("Error during optimization. You may want to try different starting parameter and / or solver settings.", call. = FALSE)
      }
    )


  }


  if (result$convergence %in% c(1, 2)) {
    warning("Convergence failed. You may want to try different starting parameter and / or solver settings.", call. = FALSE)
  }

  pars <- result$pars

  # Compute log-likelihood which is shifted from the obtained negative
  # log-likelihood at the optimum; "-tail(result$values, 1)" is the
  # log-likelihood at the optimum for the scaled data, which is too large by
  # "n * log(scale_const)" in comparison to original data
  llhood <- -tail(result$values, 1) - n * log(scale_const)

  beta_names <- paste0("beta", 1:q)

  extra_par_name <- switch(
    cond_d,
    "std" = "df",
    "sstd" = "df",
    "ged" = "shape",
    "sged" = "shape",
    "ald" = "P",
    "sald" = "P",
    character(0)
  )


  delta_name <- if (model_type %in% c("aparch", "fiaparch")) {
    "delta"
  } else {
    character(0)
  }

  gamma_name <- if (model_type %in% c("fiaparch", "figjrgarch", "fitgarch")) {
    "gamma"
  } else if (model_type %in% c("aparch", "gjrgarch", "tgarch")) {
    paste0("gamma", 1:p)
  } else {
    character(0)
  }

  par_names <- c(
    list(NULL, "mu")[[incl_mean + 1]],
    names_ar,
    names_ma,
    D_par_name,
    "omega",
    paste0("phi", 1:p),
    beta_names,
    gamma_name,
    delta_name,
    d_name,
    extra_par_name,
    list(NULL, "skew")[[skew_par + 1]]
  )

  if (cond_d %in% c("ald", "sald")) {

    # Fix dfun2
    dfun2_P <- function(x, mu, sigt, shape, skew) {
      dfun2(x = x, mu = mu, sigt = sigt, shape = P_sel, skew = skew)
    }

    # Fix the final goal function to optimize over
    # using dfun1_P and dfun2_P
    goal_fun <- function(theta, rt) {
      goal_fun_s(theta = theta, rt = rt, dfun1 = NULL, dfun2 = dfun2_P)
    }

  }

  # Compute the conditional standard deviations at the optimum
  final_series <- do.call(
    sigt_fun_adj,
    args = c(list(
      theta = unname(pars),
      x = rt_core,
      dfun = dfun1,
      p = p,
      q = q,
      p_ar = p_ar,
      q_ma = q_ma,
      incl_mean = incl_mean,
      extra_par = extra_par,
      skew_par = skew_par,
      trunc = trunc,
      presample = presample
    ), extra_args)
  )

  delta <- switch(
    model_type,
    "garch" = 2,
    "gjrgarch" = 2,
    "tgarch" = 1,
    "aparch" = pars[[n_arma_pars + 2 * p + q + 2]],
    "figarch" = 2,
    "figjrgarch" = 2,
    "fitgarch" = 1,
    "fiaparch" = pars[[n_arma_pars + 2 * p + q + 2]]
  )

  sc_delta <- scale_const^delta

  type_inp <- c("aparch", "fiaparch")[[lm + 1]]
  # Also updates "pars", if ald or sald
  compute_vcov(
    goal_fun = goal_fun,
    pars = pars,
    rt_core = rt_core,
    scale_const = scale_const,
    incl_mean = incl_mean,
    P_sel = P_sel,
    cond_d = cond_d,
    par_names = par_names,
    sc_delta = sc_delta,
    n_arma_pars = n_arma_pars,
    model_type = type_inp,     # treat as APARCH / FIAPARCH with delta fixed at 2 for rescaling
    skip_vcov = skip_vcov
  )

  rescale_estim(
    pars = pars,
    rt_core_o = rt_core_o,
    sigt = final_series$sigt,
    cmeans = final_series$cmeans,
    scale_const = scale_const,
    n_arma_pars = n_arma_pars,
    incl_mean = incl_mean,
    model_type = type_inp,  # treat as APARCH / FIAPARCH with delta fixed at 2 for rescaling
    sc_delta = sc_delta
  )

  names(pars) <- par_names

  # Makes "aic" and "bic" available
  crit_calc(pars = pars, rt = rt_core_o, llhood = llhood)

  # Apply time series formatting (if relevant) to all
  # estimated series
  series_out <- format_applier_ts(
    rt = train_obs,
    list_of_ts = list(
      "cmeans" = cmeans,
      "sigt" = sigt,
      "etat" = etat
    )
  )

  cmeans <- series_out$cmeans
  sigt <- series_out$sigt
  etat <- series_out$etat

  construct_fun <- switch(
    model_type,
    "garch" = fEGarch_fit_garch,
    "gjrgarch" = fEGarch_fit_gjrgarch,
    "tgarch" = fEGarch_fit_tgarch,
    "aparch" = fEGarch_fit_aparch,
    "figarch" = fEGarch_fit_figarch,
    "figjrgarch" = fEGarch_fit_figjrgarch,
    "fitgarch" = fEGarch_fit_fitgarch,
    "fiaparch" = fEGarch_fit_fiaparch
  )

  list_out <- construct_fun(
    pars = pars,
    se = serrors,
    vcov_mat = vcov_mat,
    rt = train_obs,
    sigt = sigt,
    cmeans = cmeans,
    scale_fun = NULL,
    etat = etat,
    llhood = llhood,
    inf_criteria = c("aic" = aic, "bic" = bic),
    cond_dist = cond_d,
    orders = orders,
    long_memo = lm,
    meanspec = meanspec,
    test_obs = test_obs,
    nonpar_model = NULL,
    trunc = trunc_i
  )


  if (use_nonpar) {
    mu <- rep(nonpar_result$mu, length(rt))
    list_out@rt <- nonpar_result$rt_train
    format_list <- format_applier_ts(
      rt = nonpar_result$rt_train,
      list_of_ts = list(
        scale_fun = nonpar_result$scale_fun,
        cmeans = mu
      )
    )
    list_out@scale_fun <- format_list$scale_fun
    list_out@cmeans <- format_list$cmeans
    list_out@test_obs <- nonpar_result$test_obs
    list_out@sigt <- list_out@sigt * format_list$scale_fun
    list_out@nonpar_model <- nonpar_result$est_nonpar
  }

  list_out

}
