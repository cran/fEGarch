
setClass("fEGarch_fit",
  slots = c(
    pars = "numeric",
    se = "numeric",
    vcov_mat = "matrix",
    rt = "ANY",
    cmeans = "ANY",
    sigt = "ANY",
    etat = "ANY",
    scale_fun = "ANY",
    llhood = "numeric",
    inf_criteria = "numeric",
    cond_dist = "character",
    orders = "numeric",
    long_memo = "logical",
    meanspec = "mean_spec",
    test_obs = "ANY",
    nonpar_model = "ANY",
    trunc = "ANY"
  ),
  prototype = list(
    pars = numeric(1),
    se = numeric(1),
    vcov_mat = matrix(numeric(1)),
    rt = numeric(1),
    cmeans = numeric(1),
    sigt = numeric(1),
    etat = numeric(1),
    scale_fun = numeric(1),
    llhood = numeric(1),
    inf_criteria = numeric(1),
    cond_dist = character(1),
    orders = numeric(1),
    long_memo = logical(1),
    meanspec = mean_spec(),
    test_obs = numeric(1),
    nonpar_model = numeric(1),
    trunc = numeric(1)
  )
)

setClass("fEGarch_fit_egarch",
  contains = "fEGarch_fit",
  slots = c(
    powers = "numeric",
    modulus = "logical"
  ),
  prototype = list(
    orders = numeric(1),
    modulus = logical(1)
  )
)

setClass("fEGarch_fit_loggarch",
  contains = "fEGarch_fit"
)

setClass("fEGarch_fit_fiaparch",
  contains = "fEGarch_fit"
)

setClass("fEGarch_fit_aparch",
  contains = "fEGarch_fit"
)

setClass("fEGarch_fit_figjrgarch",
  contains = "fEGarch_fit"
)

setClass("fEGarch_fit_gjrgarch",
  contains = "fEGarch_fit"
)

setClass("fEGarch_fit_fitgarch",
  contains = "fEGarch_fit"
)

setClass("fEGarch_fit_tgarch",
  contains = "fEGarch_fit"
)

setClass("fEGarch_fit_garch",
  contains = "fEGarch_fit"
)

setClass("fEGarch_fit_figarch",
  contains = "fEGarch_fit"
)

setClass("fEGarch_fit_rugarch_wrapper",
  contains = "fEGarch_fit",
  slots = c(
    rugarch_model = "ANY",
    model = "character"
  ),
  prototype = list(
    rugarch_model = numeric(1),
    model = character(1)
  )
)

fEGarch_fit_egarch <- function(
    pars,
    se,
    vcov_mat,
    rt,
    cmeans,
    sigt,
    etat,
    scale_fun,
    llhood,
    inf_criteria,
    cond_dist,
    orders,
    powers,
    modulus,
    long_memo,
    meanspec,
    test_obs,
    nonpar_model,
    trunc
) {

  new("fEGarch_fit_egarch",
    pars = pars,
    se = se,
    vcov_mat = vcov_mat,
    rt = rt,
    cmeans = cmeans,
    sigt = sigt,
    etat = etat,
    llhood = llhood,
    inf_criteria = inf_criteria,
    cond_dist = cond_dist,
    orders = orders,
    powers = powers,
    modulus = modulus,
    long_memo = long_memo,
    meanspec = meanspec,
    test_obs = test_obs,
    scale_fun = scale_fun,
    nonpar_model = nonpar_model,
    trunc = trunc
  )

}

fEGarch_fit_loggarch <- function(
    pars,
    se,
    vcov_mat,
    rt,
    sigt,
    cmeans,
    etat,
    llhood,
    inf_criteria,
    cond_dist,
    orders,
    long_memo,
    meanspec,
    test_obs,
    scale_fun,
    nonpar_model,
    trunc

) {

  new("fEGarch_fit_loggarch",
    pars = pars,
    se = se,
    vcov_mat = vcov_mat,
    rt = rt,
    sigt = sigt,
    cmeans = cmeans,
    etat = etat,
    llhood = llhood,
    inf_criteria = inf_criteria,
    cond_dist = cond_dist,
    orders = orders,
    long_memo = long_memo,
    meanspec = meanspec,
    test_obs = test_obs,
    scale_fun = scale_fun,
    nonpar_model = nonpar_model,
    trunc = trunc
  )

}


fEGarch_fit_fiaparch <- function(
    pars,
    se,
    vcov_mat,
    rt,
    sigt,
    cmeans,
    etat,
    llhood,
    inf_criteria,
    cond_dist,
    orders,
    long_memo,
    meanspec,
    test_obs,
    scale_fun,
    nonpar_model,
    trunc
) {

  new("fEGarch_fit_fiaparch",
    pars = pars,
    se = se,
    vcov_mat = vcov_mat,
    rt = rt,
    sigt = sigt,
    cmeans = cmeans,
    etat = etat,
    llhood = llhood,
    inf_criteria = inf_criteria,
    cond_dist = cond_dist,
    orders = orders,
    long_memo = long_memo,
    meanspec = meanspec,
    test_obs = test_obs,
    scale_fun = scale_fun,
    nonpar_model = nonpar_model,
    trunc = trunc
  )

}

fEGarch_fit_aparch <- function(
    pars,
    se,
    vcov_mat,
    rt,
    sigt,
    cmeans,
    etat,
    llhood,
    inf_criteria,
    cond_dist,
    orders,
    long_memo,
    meanspec,
    test_obs,
    scale_fun,
    nonpar_model,
    trunc
) {

  new("fEGarch_fit_aparch",
    pars = pars,
    se = se,
    vcov_mat = vcov_mat,
    rt = rt,
    sigt = sigt,
    cmeans = cmeans,
    etat = etat,
    llhood = llhood,
    inf_criteria = inf_criteria,
    cond_dist = cond_dist,
    orders = orders,
    long_memo = long_memo,
    meanspec = meanspec,
    test_obs = test_obs,
    scale_fun = scale_fun,
    nonpar_model = nonpar_model,
    trunc = trunc
  )

}

fEGarch_fit_figjrgarch <- function(
    pars,
    se,
    vcov_mat,
    rt,
    sigt,
    cmeans,
    etat,
    llhood,
    inf_criteria,
    cond_dist,
    orders,
    long_memo,
    meanspec,
    test_obs,
    scale_fun,
    nonpar_model,
    trunc
) {

  new("fEGarch_fit_figjrgarch",
    pars = pars,
    se = se,
    vcov_mat = vcov_mat,
    rt = rt,
    sigt = sigt,
    cmeans = cmeans,
    etat = etat,
    llhood = llhood,
    inf_criteria = inf_criteria,
    cond_dist = cond_dist,
    orders = orders,
    long_memo = long_memo,
    meanspec = meanspec,
    test_obs = test_obs,
    scale_fun = scale_fun,
    nonpar_model = nonpar_model,
    trunc = trunc
  )

}

fEGarch_fit_gjrgarch <- function(
    pars,
    se,
    vcov_mat,
    rt,
    sigt,
    cmeans,
    etat,
    llhood,
    inf_criteria,
    cond_dist,
    orders,
    long_memo,
    meanspec,
    test_obs,
    scale_fun,
    nonpar_model,
    trunc
) {

  new("fEGarch_fit_gjrgarch",
    pars = pars,
    se = se,
    vcov_mat = vcov_mat,
    rt = rt,
    sigt = sigt,
    cmeans = cmeans,
    etat = etat,
    llhood = llhood,
    inf_criteria = inf_criteria,
    cond_dist = cond_dist,
    orders = orders,
    long_memo = long_memo,
    meanspec = meanspec,
    test_obs = test_obs,
    scale_fun = scale_fun,
    nonpar_model = nonpar_model,
    trunc = trunc
  )

}


fEGarch_fit_fitgarch <- function(
    pars,
    se,
    vcov_mat,
    rt,
    sigt,
    cmeans,
    etat,
    llhood,
    inf_criteria,
    cond_dist,
    orders,
    long_memo,
    meanspec,
    test_obs,
    scale_fun,
    nonpar_model,
    trunc
) {

  new("fEGarch_fit_fitgarch",
    pars = pars,
    se = se,
    vcov_mat = vcov_mat,
    rt = rt,
    sigt = sigt,
    cmeans = cmeans,
    etat = etat,
    llhood = llhood,
    inf_criteria = inf_criteria,
    cond_dist = cond_dist,
    orders = orders,
    long_memo = long_memo,
    meanspec = meanspec,
    test_obs = test_obs,
    scale_fun = scale_fun,
    nonpar_model = nonpar_model,
    trunc = trunc
  )

}

fEGarch_fit_tgarch <- function(
    pars,
    se,
    vcov_mat,
    rt,
    sigt,
    cmeans,
    etat,
    llhood,
    inf_criteria,
    cond_dist,
    orders,
    long_memo,
    meanspec,
    test_obs,
    scale_fun,
    nonpar_model,
    trunc
) {

  new("fEGarch_fit_tgarch",
    pars = pars,
    se = se,
    vcov_mat = vcov_mat,
    rt = rt,
    sigt = sigt,
    cmeans = cmeans,
    etat = etat,
    llhood = llhood,
    inf_criteria = inf_criteria,
    cond_dist = cond_dist,
    orders = orders,
    long_memo = long_memo,
    meanspec = meanspec,
    test_obs = test_obs,
    scale_fun = scale_fun,
    nonpar_model = nonpar_model,
    trunc = trunc
  )

}

fEGarch_fit_garch <- function(
    pars,
    se,
    vcov_mat,
    rt,
    sigt,
    cmeans,
    etat,
    llhood,
    inf_criteria,
    cond_dist,
    orders,
    long_memo,
    meanspec,
    test_obs,
    scale_fun,
    nonpar_model,
    trunc
) {

  new("fEGarch_fit_garch",
    pars = pars,
    se = se,
    vcov_mat = vcov_mat,
    rt = rt,
    sigt = sigt,
    cmeans = cmeans,
    etat = etat,
    llhood = llhood,
    inf_criteria = inf_criteria,
    cond_dist = cond_dist,
    orders = orders,
    long_memo = long_memo,
    meanspec = meanspec,
    test_obs = test_obs,
    scale_fun = scale_fun,
    nonpar_model = nonpar_model,
    trunc = trunc
  )

}


fEGarch_fit_figarch <- function(
    pars,
    se,
    vcov_mat,
    rt,
    sigt,
    cmeans,
    etat,
    llhood,
    inf_criteria,
    cond_dist,
    orders,
    long_memo,
    meanspec,
    test_obs,
    scale_fun,
    nonpar_model,
    trunc
) {

  new("fEGarch_fit_figarch",
    pars = pars,
    se = se,
    vcov_mat = vcov_mat,
    rt = rt,
    sigt = sigt,
    cmeans = cmeans,
    etat = etat,
    llhood = llhood,
    inf_criteria = inf_criteria,
    cond_dist = cond_dist,
    orders = orders,
    long_memo = long_memo,
    meanspec = meanspec,
    test_obs = test_obs,
    scale_fun = scale_fun,
    nonpar_model = nonpar_model,
    trunc = trunc
  )

}

fEGarch_fit_rugarch_wrapper <- function(
    pars,
    se,
    vcov_mat,
    rt,
    sigt,
    cmeans,
    etat,
    llhood,
    inf_criteria,
    cond_dist,
    orders,
    long_memo,
    rugarch_model,
    model,
    meanspec,
    test_obs,
    scale_fun,
    nonpar_model,
    trunc
) {

  new("fEGarch_fit_rugarch_wrapper",
    pars = pars,
    se = se,
    vcov_mat = vcov_mat,
    rt = rt,
    sigt = sigt,
    cmeans = cmeans,
    etat = etat,
    llhood = llhood,
    inf_criteria = inf_criteria,
    cond_dist = cond_dist,
    orders = orders,
    long_memo = long_memo,
    rugarch_model = rugarch_model,
    model = model,
    meanspec = meanspec,
    test_obs = test_obs,
    scale_fun = scale_fun,
    nonpar_model = nonpar_model,
    trunc = trunc
  )

}

#================================================#

