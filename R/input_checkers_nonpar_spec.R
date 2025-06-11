# For class "mean_spec"

poly_checker_locpol <- function(poly_order) {
  !((length(poly_order) == 1) && (poly_order %in% c(1, 3)))
}

kernel_checker_locpol <- function(kernel_order) {
  !((length(kernel_order) == 1) && (kernel_order %in% c(0, 1, 2, 3)))
}

bmethod_checker_locpol <- function(boundary_method) {
  !((length(boundary_method) == 1) && (boundary_method %in% c("shorten", "extend")))
}

bwidth_checker_locpol <- function(bwidth) {
  !((is.null(bwidth)) || (is.numeric(bwidth) && length(bwidth) == 1 && bwidth > 0 && bwidth < 0.5))
}

# Function "lmemo_checker" can be reused from "inpout_checkers_egarch_spec.R".
