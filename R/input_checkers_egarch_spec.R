# For class "egarch_spec"



order_checker <- function(garch_order) {
  !((length(garch_order) == 2) && (garch_order[[1]] > 0 && garch_order[[2]] > 0))
}

power_checker <- function(power_transf) {
  !(length(power_transf) == 2)
}

lmemo_checker <- function(long_memo) {
  !(length(long_memo) == 1)
}

modulus_checker <- function(modulus) {
  !(length(modulus) == 2)
}

model_type_checker <- function(model_type) {
  !((length(model_type) == 1) && (model_type %in% c(1, 2)))
}

distr_checker <- function(distribution) {
  !((length(distribution) == 1) && (distribution %in% c("norm", "snorm", "ged", "ald", "sged", "std", "sstd", "sald")))
}
