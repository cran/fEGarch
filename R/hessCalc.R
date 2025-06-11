# Calculation of Hessian matrix using numDeriv::hessian;
# usually produces more stable hessian matrix.
hessCalc <- function(pars, func, control.args = list(
				eps = 1e-4,
				d = 1e-3,
				zero.tol = sqrt(.Machine$double.eps / 7e-7),
				r = 4,
				v = 2,
				show.details = FALSE), rt, adj = TRUE) {


  d_init <- control.args$d
  eps_init <- control.args$eps

  out <- tryCatch(expr = {numDeriv::hessian(func = func, x = pars, method = "Richardson",
                    method.args = control.args, rt = rt)},
                  error = function(e1) {NaN})

  sol <- tryCatch(expr = {solve(out)}, error = function(e1) {-1})

  # Try different settings
  while (control.args$eps > 1e-10 && (any(is.infinite(out)) || any(is.nan(out)) || any(diag(sol) < 0))) {

    control.args$eps <- control.args$eps / 10
    out <- tryCatch(expr = {numDeriv::hessian(func = func, x = pars, method = "Richardson",
              method.args = control.args, rt = rt)},
              error = function(e1) {NaN})

    sol <- tryCatch(expr = {solve(out)}, error = function(e1) {-1})

  }

  control.args$eps <- eps_init
  control.args$d <- d_init

  # Try different settings
  while (control.args$d > 1e-8 && (any(is.infinite(out)) || any(is.nan(out)) || any(diag(sol) < 0))) {

    control.args$d <- control.args$d / 10
    out <- tryCatch(expr= {numDeriv::hessian(func = func, x = pars, method = "Richardson",
              method.args = control.args, rt = rt)},
              error = function(e1) {NaN})

    sol <- tryCatch(expr = {solve(out)}, error = function(e1) {-1})

  }

  if (any(is.infinite(out)) || any(is.nan(out)) || any(is.na(tryCatch(expr = {solve(out)}, error = function(e1) {NA})))) {
    out <- matrix(NA, nrow = length(pars), ncol = length(pars))
    warning("Unable to compute Hessian matrix. No standard errors available as consequence.")
  }

  out
}

