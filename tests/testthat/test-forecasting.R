window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
rt <- window.zoo(SP500, end = "2002-12-31")

test_that("fEGarch forecast in mean works as intended", {

    # Parametric models with mean

    expect_equal({
      est <- loggarch_spec() %>%
        fEGarch(rt, n_test = 250)
      fcast <- predict_roll(est)
      zoo::coredata(fcast@cmeans)
    }, {
      est <- loggarch_spec() %>%
        fEGarch(rt, n_test = 250)
      rep(est@pars[["mu"]], 250)
    })


    expect_equal({
      est <- egarch_spec() %>%
        fEGarch(rt, n_test = 250)
      fcast <- predict_roll(est)
      zoo::coredata(fcast@cmeans)
    }, {
      est <- egarch_spec() %>%
        fEGarch(rt, n_test = 250)
      rep(est@pars[["mu"]], 250)
    })

    expect_equal({
      est <- suppressWarnings(aparch(rt, n_test = 250))
      fcast <- predict_roll(est)
      zoo::coredata(fcast@cmeans)
    }, {
      est <- suppressWarnings(aparch(rt, n_test = 250))
      rep(est@pars[["mu"]], 250)
    })

    expect_equal({
      est <- suppressWarnings(fiaparch(rt, n_test = 250))
      fcast <- predict_roll(est)
      zoo::coredata(fcast@cmeans)
    }, {
      est <- suppressWarnings(fiaparch(rt, n_test = 250))
      rep(est@pars[["mu"]], 250)
    })

    # Parametric models without mean

    expect_equal({
      est <- loggarch_spec() %>%
        fEGarch(rt, n_test = 250, meanspec = mean_spec(include_mean = FALSE))
      fcast <- predict_roll(est)
      zoo::coredata(fcast@cmeans)
    }, {
      rep(0, 250)
    })


    expect_equal({
      est <- egarch_spec() %>%
        fEGarch(rt, n_test = 250, meanspec = mean_spec(include_mean = FALSE))
      fcast <- predict_roll(est)
      zoo::coredata(fcast@cmeans)
    }, {
      rep(0, 250)
    })

    expect_equal({
      est <- suppressWarnings(aparch(rt, n_test = 250, meanspec = mean_spec(include_mean = FALSE)))
      fcast <- predict_roll(est)
      zoo::coredata(fcast@cmeans)
    }, {
      rep(0, 250)
    })

    expect_equal({
      est <- suppressWarnings(fiaparch(rt, n_test = 250, meanspec = mean_spec(include_mean = FALSE)))
      fcast <- predict_roll(est)
      zoo::coredata(fcast@cmeans)
    }, {
      rep(0, 250)
    })

    # Semiparametric models

    expect_equal({
      est <- loggarch_spec() %>%
        fEGarch(rt, n_test = 250, use_nonpar = TRUE)
      fcast <- predict_roll(est)
      zoo::coredata(fcast@cmeans)
    }, {
      rep(mean(utils::head(rt, -250)), 250)
    })


    expect_equal({
      est <- egarch_spec() %>%
        fEGarch(rt, n_test = 250, use_nonpar = TRUE)
      fcast <- predict_roll(est)
      zoo::coredata(fcast@cmeans)
    }, {
      rep(mean(utils::head(rt, -250)), 250)
    })

    expect_equal({
      est <- suppressWarnings(aparch(rt, n_test = 250, use_nonpar = TRUE))
      fcast <- predict_roll(est)
      zoo::coredata(fcast@cmeans)
    }, {
      rep(mean(utils::head(rt, -250)), 250)
    })

    expect_equal({
      est <- suppressWarnings(fiaparch(rt, n_test = 250, use_nonpar = TRUE))
      fcast <- predict_roll(est)
      zoo::coredata(fcast@cmeans)
    }, {
      rep(mean(utils::head(rt, -250)), 250)
    })

})
