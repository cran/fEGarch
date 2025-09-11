# Split series into training and test set
# and keep and omit potential time series formatting
initial_split <- function(rt, n_test) {

  pf <- parent.frame()

  split_rt <- split_ts(rt, n_test)

  rt_core_o <- zoo::coredata(split_rt$train)
  assign("rt_core_o", value = rt_core_o,
         pos = pf)
  #rt_core_o <- zoo::coredata(split_rt$train)

  scale_const <- stats::sd(rt_core_o)
  assign("scale_const", value = scale_const, pos = pf)
  #scale_const <- sd(rt_core_o)         # for rescaling later
  #rt_core <- rt_core_o / scale_const   # now has sample variance 1
  assign("rt_core", value = rt_core_o / scale_const,
         pos = pf)

  assign("test_obs", value = split_rt$test,
         pos = pf)
  assign("train_obs", value = split_rt$train,
         pos = pf)

  # list(
  #   scale_const = scale_const,
  #   rt_core = rt_core,
  #   rt_core_o = rt_core_o,
  #   test_obs = split_rt$test,
  #   train_obs = split_rt$train
  # )

}

# Compute the negative log-likelihood from
# estimated input series, observations and
# the density function
nllhood_calc <- function(input_list, x, dfun) {

  sigt <- input_list$sigt

  cmeans <- input_list$cmeans
  skew <- input_list$skew
  shape <- input_list$shape

  lhoods <- dfun(x, cmeans, sigt, shape, skew)
  lhoods[lhoods == 0] <- 1e-25   # safety measure to avoid -Inf llhood
  -sum(log(lhoods))

}

# Calculate AIC and BIC from fitted parameter vector length,
# observation number and log-likelihood value
crit_calc <- function(pars, rt, llhood) {

  pf <- parent.frame()
  n <- length(rt)
  k <- length(pars)
  m2llhood <- -2 * llhood

  assign("aic", value = (m2llhood + 2 * k) / n, pos = pf)
  assign("bic", value = (m2llhood + k * log(n)) / n, pos = pf)

}

identify_arma_settings <- function(meanspec, Drange) {

  pf <- parent.frame()

  # ARMA orders
  order_arma <- orders(meanspec)
  p_ar <- order_arma[[1]]
  q_ma <- order_arma[[2]]

  # Include unconditional mean parameter?
  incl_mean <- meanspec@include_mean

  # Long-memory setting (for ARMA-part)
  lm_arma <- long_memo(meanspec)

  # Assign ARMA / FARIMA function
  # and parameter grabbing function for ARMA / FARIMA
  arma_fun <- if (lm_arma) {
    grab_fun <- farima_grab_pars
    D_par_name <- "D"
    D_start <- 0.25
    D_low <- Drange[[1]]
    D_up <- Drange[[2]]
    farima_fit
  } else {
    grab_fun <- arma_grab_pars
    D_par_name <- character(0)
    D_start <- D_low <- D_up <- numeric(0)
    arma_fit
  }

  # Adjust (if ARMA relevant)
  pg0 <- p_ar > 0
  qg0 <- q_ma > 0
  # Check if adjustment of sigt-function required due to
  # conditional mean part
  mean_adj_required <- pg0 || qg0 || lm_arma
  # Number of parameters for the mean specification
  n_arma_pars <- incl_mean + p_ar + q_ma + lm_arma

  # Positions of AR-parameters in parameter vector
  low_ar <- 1 + incl_mean
  up_ar <- incl_mean + p_ar

  assign("p_ar", value = p_ar, pos = pf)
  assign("q_ma", value = q_ma, pos = pf)
  assign("incl_mean", value = incl_mean, pos = pf)
  assign("lm_arma", value = lm_arma, pos = pf)
  assign("arma_grab_fun", value = grab_fun, pos = pf)
  assign("cmean_fun", value = arma_fun, pos = pf)
  assign("pg0", value = pg0, pos = pf)
  assign("qg0", value = qg0, pos = pf)
  assign("mean_adj_required", value = mean_adj_required, pos = pf)
  assign("n_arma_pars", value = n_arma_pars, pos = pf)
  assign("D_par_name", value = D_par_name, pos = pf)
  assign("D_start", value = D_start, pos = pf)
  assign("D_low", value = D_low, pos = pf)
  assign("D_up", value = D_up, pos = pf)
  assign("low_ar", value = low_ar, pos = pf)
  assign("up_ar", value = up_ar, pos = pf)
  assign("ar_start", value = rep(0.05 / p_ar, p_ar), pos = pf)
  assign("ma_start", value = rep(0.05 / q_ma, q_ma), pos = pf)
  assign("ar_lb", value = rep(-2, p_ar), pos = pf)
  assign("ma_lb", value = rep(-2, q_ma), pos = pf)
  assign("ar_ub", value = rep(2, p_ar), pos = pf)
  assign("ma_ub", value = rep(2, q_ma), pos = pf)
  assign("names_ar", value = list(character(0), paste0("ar", 1:p_ar))[[pg0 + 1]], pos = pf)
  assign("names_ma", value = list(character(0), paste0("ma", 1:q_ma))[[qg0 + 1]], pos = pf)

}

# Function to adjust the calculations for the conditional means and conditional
# standard deviations for an ARMA or FARIMA extended volatility model;
# returns the base function for each model, if no additional ARMA / FARIMA
# settings are given
adj_sigt_fun <- function(p_ar, q_ma, pg0, qg0, lm_arma, incl_mean, n_arma_pars,
                         grab_fun, arma_fun,
            model_type = c("garch", "figarch", "aparch", "fiaparch", "gjrgarch",
                      "figjrgarch", "loggarch", "filoggarch", "egarch", "fiegarch",
                      "tgarch", "fitgarch")) {

  check <- pg0 || qg0 || lm_arma

  model_type <- match.arg(model_type)

  # Select the correct base function for the model type
  sigt_fun <- switch(
    model_type,
    "garch" = sigt_garch_R,
    "figarch" = sigt_figarch_R,
    "aparch" = sigt_aparch_R,
    "fiaparch" = sigt_fiaparch_R,
    "gjrgarch" = sigt_gjrgarch_R,
    "figjrgarch" = sigt_figjrgarch_R,
    "tgarch" = sigt_tgarch_R,
    "fitgarch" = sigt_fitgarch_R,
    "loggarch" = sigt_loggarch_short_R,
    "filoggarch" = sigt_loggarch_long_R,
    "egarch" = sigt_egarch_short_R,
    "fiegarch" = sigt_egarch_long_R
  )

  sigt_fun_adj <- if (check) {
    function(theta, x, p, q, p_ar, q_ma,
         incl_mean, extra_par, skew_par, trunc, presample, ...) {

      # Collect input arguments (including also the ellipsis ...)
      # as a named list
      args <- mget(names(formals()), environment())
      args[["..."]] <- NULL
      dot_args <- list(...)
      all_args <- c(args, dot_args)
      rm(dot_args)
      rm(args)

      # Obtain ARMA / FARIMA parameters from theta
      arma_pars <- grab_fun(
        theta = theta,
        incl_mean = incl_mean,
        pg0 = pg0,
        qg0 = qg0,
        p = p_ar,
        q = q_ma
      )
      mu <- arma_pars$mu
      ar <- arma_pars$ar
      ma <- arma_pars$ma

      coef_inf <- if (lm_arma) {
        coefs <- ar_infty(ar = ar, ma = ma, d = arma_pars$D, max_i = trunc)
        coefs[[1]] <- 0
        coefs
      } else {
        0
      }

      # Calculate fitted values following these parameters...
      fitted_vs <- arma_fun(x = x, mu = mu, ma = ma, ar = ar, coef_inf = coef_inf, presample = presample)

      # ... adjust input arguments from intermediate results ...
      all_args[["theta"]] <- theta[-(1:n_arma_pars)]     # Remove the first few elements from theta (those for the parameters relevant for ARMA / FARIMA)
      all_args[["x"]] <- x - fitted_vs                   # Compute conditional mean adjusted values
      all_args[["incl_mean"]] <- FALSE                   # Remove unconditional mean for further steps

      #... and calculate conditional standard deviations from the residuals
      # following the specified GARCH-type
      out <- do.call(sigt_fun, args = all_args)
      out$cmeans <- fitted_vs
      out
    }

  } else {
    sigt_fun             # ...else keep function for volatility only
  }

  sigt_fun_adj

}

# Function to grab conditional distribution parameter restrictions from
# internally saved tables
cond_d_par_restr <- function(cond_d) {

  pf <- parent.frame()

  vals <- lookup_table[["vals"]][[cond_d]]

  assign("sval", value = vals[[1]], pos = pf)
  assign("lval", value = vals[[2]], pos = pf)
  assign("uval", value = vals[[3]], pos = pf)

}

set_start_LB_UB <- function(
    mean_val,
    ar_start, ma_start, D_start,
    omega_val,
    phi_start, psi_start,
    add_pars_start,
    d_start,
    sval, skew_start,
    mean_low,
    ar_low, ma_low, D_low,
    omega_low,
    phi_low,
    psi_low,
    add_low,
    d_low,
    lval, skew_low,
    mean_up,
    ar_up, ma_up, D_up,
    omega_up,
    phi_up, psi_up, add_up, d_up,
    uval, skew_up,
    start_pars, LB, UB
) {

  pf <- parent.frame()

  if (is.null(start_pars)) {
    start_pars <- c(
      mean_val,
      ar_start,
      ma_start,
      D_start,
      omega_val,
      phi_start,
      psi_start,
      add_pars_start,
      d_start,
      sval, skew_start
    )
  }

  if (is.null(LB)) {
    LB <- c(
      mean_low,
      ar_low,
      ma_low,
      D_low,
      omega_low,
      phi_low,
      psi_low,
      add_low,
      d_low,
      lval, skew_low
    )
  }

  if (is.null(UB)) {
    UB <- c(
      mean_up,
      ar_up,
      ma_up,
      D_up,
      omega_up,
      phi_up,
      psi_up,
      add_up,
      d_up,
      uval, skew_up
    )
  }

  assign("start_pars", value = start_pars, pos = pf)
  assign("LB", value = LB, pos = pf)
  assign("UB", value = UB, pos = pf)

}

rescale_estim <- function(pars, rt_core_o, sigt, cmeans, scale_const, n_arma_pars, incl_mean,
                          model_type, ...) {

  pf <- parent.frame()

  list2env(list(...), envir = environment())

  # Adjust parameters according to scaling where necessary
  if (incl_mean) {
    pars[[1]] <- pars[[1]] * scale_const    # mean adjustment
  }
  # omega_sig
  if (model_type %in% c("egarch", "fiegarch", "loggarch", "filoggarch")) {
    pars[[n_arma_pars + 1]] <- pars[[n_arma_pars + 1]] + 2 * log(scale_const)
  } else if (model_type %in% c("aparch", "fiaparch")) {
    pars[[n_arma_pars + 1]] <- pars[[n_arma_pars + 1]] * sc_delta
  }

  sigt <- sigt * scale_const
  cmeans <- cmeans * scale_const

  assign("sigt", value = sigt, pos = pf)
  assign("cmeans", value = cmeans, pos = pf)
  assign("etat", value = (rt_core_o - cmeans) / sigt, pos = pf)
  assign("pars", value = pars, pos = pf)

}

ald_fit <- function(Prange, parallel, ncores, dfun1, dfun2, goal_fun_s, start_pars,
                    LB, UB, ineqfun, ineqLB, ineqUB, rt_core, control) {

  pf <- parent.frame()

  P_range <- Prange   # P-values to consider
  oldplan <- future::plan()    # Save old plan settings
  if (parallel) {
    future::plan(future::multisession, workers = ncores)   # Plan multisession, i.e. parallel, programming
    on.exit(expr = {future::plan(oldplan)}, add = TRUE, after = TRUE)  # Fail-safe: in case of errors, old plan settings are reinstated
  }

  P_seq <- P_range[[1]]:P_range[[2]]

  # Run a loop for each P to check and optimize over the remaining parameters
  fitted_models <- furrr::future_map(
    .x = P_seq,
    .f = function(.x, start_pars, goal_fun_s, LB, UB, ineqfun, ineqLB, ineqUB, control, rt_core, dfun1, dfun2) {

      # For each P, fix dfun1
      dfun1_P <- function(x, shape, skew) {
        dfun1(x = x, shape = .x, skew = skew)
      }

      # For each P, fix dfun2
      dfun2_P <- function(x, mu, sigt, shape, skew) {
        dfun2(x = x, mu = mu, sigt = sigt, shape = .x, skew = skew)
      }

      # Now fix the final goal function to optimize over
      # using dfun1_P and dfun2_P
      goal_fun <- function(theta, rt) {
        goal_fun_s(theta = theta, rt = rt, dfun1 = dfun1_P, dfun2 = dfun2_P)
      }

      est <- tryCatch(
        expr = {suppressWarnings(Rsolnp::solnp(
        pars = start_pars,
        fun = goal_fun,
        LB = LB,
        UB = UB,
        ineqfun = ineqfun,
        ineqLB = ineqLB, #rep(1e-6, p_all),    # abs. values of roots must be outside of unit circle
        ineqUB = ineqUB, #rep(Inf, p_all),
        control = control,
        rt = rt_core
      ))},
        error = function(e1) {
          stop("Error during optimization. You may want to try different starting parameter and / or solver settings.", call. = FALSE)
        }
      )

        est
      }, start_pars = start_pars, goal_fun_s = goal_fun_s, LB = LB, UB = UB,
      ineqfun = ineqfun, ineqLB = ineqLB, ineqUB = ineqUB,
      control = control, rt_core = rt_core,
      dfun1 = dfun1, dfun2 = dfun2,
      .progress = FALSE, .options = furrr::furrr_options(seed = NULL)
    )

    future::plan(oldplan)   # Reinstate old plan settings

    coll_nllh <- vapply(    # From each fitted model, obtain the negative
      X = fitted_models,    # log-likelihood at the corresponding optimum
      FUN = function(.x) {
        tail(.x$values, 1)
      },
      FUN.VALUE = numeric(1)
    )

    assign("P_sel", value = P_seq[[which.min(coll_nllh)]], pos = pf)
    assign("result", value = fitted_models[[which.min(coll_nllh)]], pos = pf)

}

compute_vcov <- function(goal_fun, pars, rt_core, scale_const, incl_mean, P_sel,
                         cond_d, par_names, model_type, ...) {

  pf <- parent.frame()

  list2env(list(...), envir = environment())

  # Compute the (negative) hessian at the optimum using a finite-difference
  # approximation via "hessian" from "numDeriv"; is usually more stable
  # nhess <- tryCatch(
  #   expr = {hessCalc(func = goal_fun, pars = unname(pars), rt = rt_core)},
  #   error = function(e1) {
  #     warning("Unable to obtain standard errors.", call. = FALSE)
  #     matrix(NA, nrow = length(pars), ncol = length(pars))
  #   }
  # )

  nhess <- hessCalc(func = goal_fun, pars = unname(pars), rt = rt_core)

  # Obtain the variance-covariance matrix as the inverse of the
  # negative hessian
  I_mat <- diag(length(pars))    # Transformation matrix to rescale variance-covariance-matrix
  if (incl_mean) {
    I_mat[1, 1] <- scale_const   # For mean; not required for omega_sig, because that transformation is just an added constant
  }
  if (model_type %in% c("aparch", "fiaparch")) {
    n_arma_im1 <- n_arma_pars + 1
    I_mat[n_arma_im1, n_arma_im1] <- sc_delta
  }


  vcov_mat <- I_mat %*% solve(nhess) %*% I_mat   # Retransformation included
  diag_vcov <- diag(vcov_mat)

  if (all(!is.na(diag_vcov)) && any(diag_vcov < 0)) {
    vcov_mat <- matrix(NA, nrow = length(pars), ncol = length(pars))
    diag_vcov <- rep(NA, length(pars))
    warning("Unable to compute Hessian matrix. No standard errors available as consequence.")
  }

  if (cond_d %in% c("ald", "sald")) {
    if (cond_d == "ald") {
      pars <- c(pars, P_sel)
      vcov_mat <- rbind(cbind(vcov_mat, NA), NA)
    } else {
      l_p <- length(pars)
      pars <- c(pars[1:(l_p - 1)], P_sel, pars[[l_p]])
      vcov_c <- cbind(vcov_mat[, 1:(l_p - 1)], NA, vcov_mat[, l_p])
      vcov_mat <- rbind(vcov_c[1:(l_p - 1), ], NA, vcov_c[l_p, ])
    }
    diag_vcov <- diag(vcov_mat)
  }

  serrors <- sqrt(diag_vcov)
  names(serrors) <- par_names
  rownames(vcov_mat) <- par_names
  colnames(vcov_mat) <- par_names

  assign("vcov_mat", value = vcov_mat, pos = pf)
  assign("serrors", value = serrors, pos = pf)
  assign("pars", value = pars, pos = pf)

}
