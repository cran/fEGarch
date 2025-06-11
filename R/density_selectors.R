fun1_selector <- function(dist) {
  switch(
    dist,
    "norm" = function(x, shape, skew) {pdf_norm_v1(x)},
    "std" = function(x, shape, skew) {pdf_std_v1(x, df = shape)},
    "ged" = function(x, shape, skew) {pdf_ged_v1(x, shape = shape)},
    "ald" = function(x, shape, skew) {pdf_ald_v1(x, P = shape)},
    "snorm" = function(x, shape, skew) {pdf_skew_norm_s(x, skew = skew)},
    "sstd" = function(x, shape, skew) {pdf_skew_sstd_s(x, df = shape, skew = skew)},
    "sged" = function(x, shape, skew) {pdf_skew_sged_s(x, shape = shape, skew = skew)},
    "sald" = function(x, shape, skew) {pdf_skew_sald_s(x, P = shape, skew = skew)}
  )
}

fun2_selector <- function(dist) {

  switch(
    dist,
    "norm" = function(x, mu, sigt, shape, skew) {pdf_norm(x, mu = mu, sigt = sigt)},
    "std" = function(x, mu, sigt, shape, skew) {pdf_std(x, mu = mu, sigt = sigt, df = shape)},
    "ged" = function(x, mu, sigt, shape, skew) {pdf_ged(x, mu = mu, sigt = sigt, shape = shape)},
    "ald" = function(x, mu, sigt, shape, skew) {pdf_ald(x, mu = mu, sigt = sigt, P = shape)},
    "snorm" = function(x, mu, sigt, shape, skew) {pdf_skew_norm_final(x, mu = mu, sigt = sigt, skew = skew)},
    "sstd" = function(x, mu, sigt, shape, skew) {pdf_skew_sstd_final(x, mu = mu, sigt = sigt, df = shape, skew = skew)},
    "sged" = function(x, mu, sigt, shape, skew) {pdf_skew_sged_final(x, mu = mu, sigt = sigt, shape = shape, skew = skew)},
    "sald" = function(x, mu, sigt, shape, skew) {pdf_skew_sald_final(x, mu = mu, sigt = sigt, P = shape, skew = skew)}
  )
}

simfun_selector <- function(dist) {

  switch(
    dist,
    "norm" = function(x, shape, skew) {rnorm_s(x)},
    "std" = function(x, shape, skew) {rstd_s(x, df = shape)},
    "ged" = function(x, shape, skew) {rged_s(x, shape = shape)},
    "ald" = function(x, shape, skew) {rald_s(x, P = shape)},
    "snorm" = function(x, shape, skew) {rsnorm_s(x, skew = skew)},
    "sstd" = function(x, shape, skew) {rsstd_s(x, df = shape, skew = skew)},
    "sged" = function(x, shape, skew) {rsged_s(x, shape = shape, skew = skew)},
    "sald" = function(x, shape, skew) {rsald_s(x, P = shape, skew = skew)}
  )

}
