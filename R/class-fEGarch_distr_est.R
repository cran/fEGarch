
setClass("fEGarch_distr_est",
  slots = c(
    pars = "numeric",
    se = "numeric",
    vcov_mat = "ANY",
    x = "ANY",
    llhood = "numeric",
    inf_criteria = "numeric",
    dist = "character",
    fixed = "numeric"
  ),
  prototype = list(
    pars = numeric(1),
    se = numeric(1),
    vcov_mat = matrix(numeric(1)),
    x = numeric(1),
    llhood = numeric(1),
    inf_criteria = numeric(1),
    dist = character(1),
    fixed = numeric(1)
  )
)


fEGarch_distr_est <- function(
    pars,
    se,
    vcov_mat,
    x,
    llhood,
    inf_criteria,
    dist,
    fixed
) {

  new("fEGarch_distr_est",
    pars = pars,
    se = se,
    vcov_mat = vcov_mat,
    x = x,
    llhood = llhood,
    inf_criteria = inf_criteria,
    dist = dist,
    fixed = fixed
  )

}


#'@export
#'@rdname show-methods
#'@aliases show,fEgarch_distr_est-method
setMethod(
  "show",
  "fEGarch_distr_est",
  function(object) {
  x <- object
  par_names <- names(x@pars)
  se <- unname(x@se)
  pars <- unname(x@pars)

  tvals <- pars / se
  atvals <- abs(tvals)

  pvals <- 2 * (1 - pnorm(atvals))

  if (length(pars) > 0) {

    df <- data.frame(
      par = sprintf("%.4f", pars),
      se = sprintf("%.4f", se),
      tval = sprintf("%.4f", tvals),
      pval = sprintf("%.4f", pvals)
    )
    row.names(df) <- par_names

  } else {
    df <- "none"
  }

  fixed <- if (length(object@fixed) > 0) {
    names_fixed <- names(object@fixed)

    fix_mean <- if ("mean" %in% names_fixed) {
      paste0(object@fixed[["mean"]], " (mean)")
    } else {
      ""
    }

    fix_sdev <- if ("sdev" %in% names_fixed) {
      paste0(object@fixed[["sdev"]], " (sdev)")
    } else {
      ""
    }

    int <- if (length(names_fixed) == 2) {
      ", "
    } else {
      ""
    }

    paste0(
      fix_mean,
      int,
      fix_sdev
    )
  } else {
    "none"
  }

  part1 <- paste0(
    "*************************************\n",
    "*         Fitted Distribution       *\n",
    "*************************************\n",
    " \n",
    "Distribution: ", object@dist, "\n",
    "Fixed: ", fixed, "\n",
    " \n",
    "Fitted parameters:\n"
  )

  part2 <- paste0(
    " \n",
    "Information criteria:\n",
    "AIC: ", sprintf("%.4f", unname(x@inf_criteria[[1]])),
    ", BIC: ", sprintf("%.4f", unname(x@inf_criteria[[2]]))
  )

  cat(part1, fill = TRUE)
  if (is.atomic(df) && is.character(df) && length(df) == 1 && df == "none") {
    cat(df, fill = TRUE)
  } else {
    print(df)
  }
  cat(part2, fill = TRUE)
  }
)

#'Methods for Accessing Distribution Estimation Elements
#'
#'Accessors to access the elements of the same name in
#'output objects returned by either \code{\link{distr_est}} or
#'its various wrappers like \code{\link{norm_est}}.
#'
#'@param x an object returned by either \code{\link{distr_est}} or
#'its various wrappers like \code{\link{norm_est}}.
#'
#'@details
#'Convenience methods to access the elements of the same name
#'that can otherwise be accessed via the operator \code{@} within
#'objects that inherit from class \code{"fEGarch_distr_est"}, which covers
#'objects returned by either \code{\link{distr_est}} or
#'its various wrappers like \code{\link{norm_est}}.
#'
#'@return
#'The element within the input object of the same name as the method
#'is returned.
#'
#'@export
#'
#'@rdname accessor_methods_distr_est
#'
#'@aliases inf_criteria,fEGarch_distr_est-method
#'
#'@examples
#'x <- rged_s(4000, shape = 1.5) * 2.1 + 3.3
#'est <- ged_est(x)
#'inf_criteria(est)
#'pars(est)
#'
setMethod("inf_criteria", "fEGarch_distr_est", function(x) {x@inf_criteria})
#'@export
#'@rdname accessor_methods_distr_est
#'@aliases llhood,fEGarch_distr_est-method
setMethod("llhood", "fEGarch_distr_est", function(x) {x@llhood})
#'@export
#'@rdname accessor_methods_distr_est
#'@aliases pars,fEGarch_distr_est-method
setMethod("pars", "fEGarch_distr_est", function(x) {x@pars})
#'@export
#'@rdname accessor_methods_distr_est
#'@aliases se,fEGarch_distr_est-method
setMethod("se", "fEGarch_distr_est", function(x) {x@se})
#'@export
#'@rdname accessor_methods_distr_est
#'@aliases vcov_mat,fEGarch_distr_est-method
setMethod("vcov_mat", "fEGarch_distr_est", function(x) {x@vcov_mat})
