window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
rt <- window.zoo(SP500, end = "2002-12-31")

test_that("fEGarch rolling forecast works as intended without refitting", {

    # Parametric models

    expect_equal({
      est <- loggarch_spec() %>%
        fEGarch(rt, n_test = 250)
      fcast <- predict_roll(est, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(0.0002681018, 0.0002681018, 0.0002681018, 0.0002681018,
        0.0002681018, 0.0002681018, 0.0002681018, 0.0002681018,
        0.0002681018, 0.0002681018, 0.0002681018, 0.0002681018,
        0.0002681018, 0.0002681018, 0.0002681018, 0.0118712574,
        0.0119869113, 0.0121606552, 0.0120240618, 0.0120440204,
        0.0151599313, 0.0144787844, 0.013663356, 0.0144471387,
        0.0124818283, 0.0137511904, 0.0136381748, 0.0132417846,
        0.0137930024, 0.0135032147)
    }, tolerance = 1e-02)


    expect_equal({
      est <- egarch_spec() %>%
        fEGarch(rt, n_test = 250)
      fcast <- predict_roll(est, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(-0.0012106234, -0.0012106234, -0.0012106234,
        -0.0012106234, -0.0012106234, -0.0012106234,
        -0.0012106234, -0.0012106234, -0.0012106234,
        -0.0012106234, -0.0012106234, -0.0012106234,
        -0.0012106234, -0.0012106234, -0.0012106234,
        0.0071533667, 0.0066413251, 0.0071996152,
        0.0073997976, 0.0077347649, 0.0136579408,
        0.0135750337, 0.0131028407, 0.0157537715,
        0.0153467836, 0.013350933, 0.0135772983,
        0.013534737, 0.0150148164, 0.0143590828)
    }, tolerance = 1e-02)

    expect_equal({
      est <- suppressWarnings(aparch(rt, n_test = 250))
      fcast <- predict_roll(est, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(-0.0011347475, -0.0011347475, -0.0011347475,
        -0.0011347475, -0.0011347475, -0.0011347475,
        -0.0011347475, -0.0011347475, -0.0011347475,
        -0.0011347475, -0.0011347475, -0.0011347475,
        -0.0011347475, -0.0011347475, -0.0011347475,
        0.0081128659, 0.0076779731, 0.0083380492,
        0.0084323347, 0.0087417991, 0.0146794224,
        0.0140671411, 0.0130785475, 0.0162560362,
        0.0150609143, 0.0136083367, 0.0136214387,
        0.0131758577, 0.0150033972, 0.0139267511)
    }, tolerance = 1e-02)

    expect_equal({
      est <- suppressWarnings(fiaparch(rt, n_test = 250))
      fcast <- predict_roll(est, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(-0.0012603512, -0.0012603512, -0.0012603512,
        -0.0012603512, -0.0012603512, -0.0012603512,
        -0.0012603512, -0.0012603512, -0.0012603512,
        -0.0012603512, -0.0012603512, -0.0012603512,
        -0.0012603512, -0.0012603512, -0.0012603512,
        0.0097863034, 0.0090630439, 0.0094991403,
        0.0099682965, 0.0101150485, 0.0143607021,
        0.0134091387, 0.0120268185, 0.015408134,
        0.0166667716, 0.0135127123, 0.0129465634,
        0.01271435, 0.014672088, 0.0151035367)
    }, tolerance = 1e-02)

    # Parametric models without mean

    expect_equal({
      est <- loggarch_spec() %>%
        fEGarch(rt, n_test = 250, meanspec = mean_spec(include_mean = FALSE))
      fcast <- predict_roll(est, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0.0123187605, 0.0124299726, 0.0125561452,
        0.012323092, 0.0122897684, 0.0148918736,
        0.0141120751, 0.013356193, 0.014213322,
        0.0104304509, 0.0137220165, 0.0135633969,
        0.0130901077, 0.0136885, 0.0134183433)
    }, tolerance = 1e-02)


    expect_equal({
      est <- egarch_spec() %>%
        fEGarch(rt, n_test = 250, meanspec = mean_spec(include_mean = FALSE))
      fcast <- predict_roll(est, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0.0076327411, 0.0071276731, 0.007742262,
        0.0079795505, 0.0083515575, 0.0140184002,
        0.0139050024, 0.0133549105, 0.0160539665,
        0.0155338961, 0.0133250533, 0.0135487596,
        0.0134947298, 0.0150048232, 0.0142858859)
    }, tolerance = 1e-02)

    expect_equal({
      est <- suppressWarnings(aparch(rt, n_test = 250, meanspec = mean_spec(include_mean = FALSE)))
      fcast <- predict_roll(est, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0.0080700191, 0.0076096983, 0.0083408617, 0.0085396191,
        0.0089237842, 0.014709234, 0.0142140729, 0.0131785277,
        0.0162924099, 0.0150603319, 0.0135763438, 0.0136512794,
        0.0133126878, 0.0150977569, 0.013978829)
    }, tolerance = 1e-02)

    expect_equal({
      est <- suppressWarnings(fiaparch(rt, n_test = 250, meanspec = mean_spec(include_mean = FALSE)))
      fcast <- predict_roll(est, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0.0097281939, 0.0089440248, 0.0096375659, 0.010420761,
        0.0108033544, 0.0151202502, 0.0142882662, 0.0127403048,
        0.0159291206, 0.0167495578, 0.0138008048, 0.0134003306,
        0.0134396787, 0.0155426014, 0.015631948)
    }, tolerance = 1e-02)


    expect_equal({
      est <- egarch_spec() %>%
        fEGarch(rt, n_test = 250, use_nonpar = TRUE)
      fcast <- predict_roll(est, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(rep(mean(utils::head(rt, -250)), 15),
        0.0065839069, 0.0060903241, 0.0068959571, 0.0073330375,
        0.0078639467, 0.0105501659, 0.0105929054, 0.0101394326,
        0.0126608854, 0.0122101582, 0.0094958681, 0.0099078573,
        0.0100469181, 0.0115377278, 0.0107846487)
    }, tolerance = 1e-02)

    expect_equal({
      est <- suppressWarnings(aparch(rt, n_test = 250, use_nonpar = TRUE))
      fcast <- predict_roll(est, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(rep(mean(utils::head(rt, -250)), 15),
        0.006929891, 0.0065812373, 0.0073409529, 0.0075926555, 0.0080062076,
        0.0122177972, 0.0117574325, 0.0107753284, 0.0134590614, 0.0122470595,
        0.0109645518, 0.0111114585, 0.01085001, 0.0124681813, 0.011390643)
    }, tolerance = 1e-02)

    expect_equal({
      est <- suppressWarnings(fiaparch(rt, n_test = 250, use_nonpar = TRUE))
      fcast <- predict_roll(est, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(rep(mean(utils::head(rt, -250)), 15),
        0.0073930335, 0.0067632251, 0.0073322214, 0.0083405667, 0.0088772311,
        0.0130886201, 0.0126423295, 0.0112292598, 0.0126671852, 0.013546193,
        0.011818647, 0.0114294269, 0.011811047, 0.0132955867, 0.0135067545)
    }, tolerance = 1e-02)

})

test_that("fEGarch rolling forecast works as intended with refitting", {

    # Parametric models

    expect_equal({
      est <- loggarch_spec() %>%
        fEGarch(rt, n_test = 250)
      fcast <- predict_roll(est, refit_after = 100, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(0.0002681018, 0.0002681018, 0.0002681018, 0.0002681018,
        0.0002681018, -0.0009007016, -0.0009007016, -0.0009007016,
        -0.0009007016, -0.0009007016, -0.0007289457, -0.0007289457,
        -0.0007289457, -0.0007289457, -0.0007289457, 0.0118712574,
        0.0119869113, 0.0121606552, 0.0120240618, 0.0120440204,
        0.0143358164, 0.0134093003, 0.0129413398, 0.0137216026,
        0.0124818474, 0.0145961792, 0.0142512692, 0.013454737,
        0.0141129234, 0.0139055308)
    }, tolerance = 1e-02)


    expect_equal({
      est <- egarch_spec() %>%
        fEGarch(rt, n_test = 250)
      fcast <- predict_roll(est, refit_after = 100, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(-0.0012106234, -0.0012106234, -0.0012106234, -0.0012106234,
        -0.0012106234, -0.0013089776, -0.0013089776, -0.0013089776,
        -0.0013089776, -0.0013089776, -0.0013939419, -0.0013939419,
        -0.0013939419, -0.0013939419, -0.0013939419, 0.0071533667,
        0.0066413251, 0.0071996152, 0.0073997976, 0.0077347649,
        0.0138130977, 0.0137706096, 0.0133237512, 0.015935356,
        0.0155708478, 0.0131230704, 0.0134685118, 0.0135902089,
        0.0149662585, 0.0144344685)
    }, tolerance = 1e-02)

    expect_equal({
      est <- suppressWarnings(aparch(rt, n_test = 250))
      fcast <- predict_roll(est, refit_after = 100, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(-0.0011347475, -0.0011347475, -0.0011347475, -0.0011347475,
        -0.0011347475, -0.0012584898, -0.0012584898, -0.0012584898,
        -0.0012584898, -0.0012584898, -0.0014218601, -0.0014218601,
        -0.0014218601, -0.0014218601, -0.0014218601, 0.0081128659,
        0.0076779731, 0.0083380492, 0.0084323347, 0.0087417991,
        0.0147972664, 0.0142126084, 0.0132290761, 0.0162419305,
        0.0150707592, 0.013952051, 0.0140575973, 0.013665882,
        0.0154222796, 0.0143225252)
    }, tolerance = 1e-02)

    expect_equal({
      est <- suppressWarnings(fiaparch(rt, n_test = 250))
      fcast <- predict_roll(est, refit_after = 100, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(-0.0012603512, -0.0012603512, -0.0012603512, -0.0012603512,
        -0.0012603512, -0.0014161309, -0.0014161309, -0.0014161309,
        -0.0014161309, -0.0014161309, -0.0013408128, -0.0013408128,
        -0.0013408128, -0.0013408128, -0.0013408128, 0.0097863034,
        0.0090630439, 0.0094991403, 0.0099682965, 0.0101150485,
        0.0139026504, 0.0130335987, 0.0117988111, 0.0148467807,
        0.0164970277, 0.0136725943, 0.0137261681, 0.01332404,
        0.0150394628, 0.0140054597)
    }, tolerance = 1e-02)

    # Parametric models without mean

    expect_equal({
      est <- loggarch_spec() %>%
        fEGarch(rt, n_test = 250, meanspec = mean_spec(include_mean = FALSE))
      fcast <- predict_roll(est, refit_after = 100, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0.0123187605, 0.0124299726, 0.0125561452, 0.012323092, 0.0122897684,
        0.0146311107, 0.0139784377, 0.0133284695, 0.0140316241, 0.0107775944,
        0.014645289, 0.014441771, 0.0138832366, 0.0144999844, 0.0141828561)
    }, tolerance = 1e-02)


    expect_equal({
      est <- egarch_spec() %>%
        fEGarch(rt, n_test = 250, meanspec = mean_spec(include_mean = FALSE))
      fcast <- predict_roll(est, refit_after = 100, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0.0076327411, 0.0071276731, 0.007742262, 0.0079795505, 0.0083515575,
        0.0139418436, 0.013920696, 0.0134616053, 0.0159700827, 0.0155891482,
        0.0130560496, 0.0134075025, 0.0135452971, 0.0148772159, 0.0143616049)
    }, tolerance = 1e-02)

    expect_equal({
      est <- suppressWarnings(aparch(rt, n_test = 250, meanspec = mean_spec(include_mean = FALSE)))
      fcast <- predict_roll(est, refit_after = 100, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0.0080700191, 0.0076096983, 0.0083408617, 0.0085396191, 0.0089237842,
        0.014721462, 0.0142726525, 0.0132652041, 0.016174958, 0.0149930434, 0.0139271161,
        0.0140794224, 0.0138225751, 0.015494205, 0.0143773198)
    }, tolerance = 1e-02)

    expect_equal({
      est <- suppressWarnings(fiaparch(rt, n_test = 250, meanspec = mean_spec(include_mean = FALSE)))
      fcast <- predict_roll(est, refit_after = 100, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0097281939, 0.0089440248,
        0.0096375659, 0.010420761, 0.0108033544, 0.0147235622, 0.0139823541,
        0.0126027358, 0.0153525611, 0.0165106421, 0.0138494556, 0.0140249392,
        0.0137830685, 0.0154744842, 0.0143274471)
    }, tolerance = 1e-02)


    expect_equal({
      est <- egarch_spec() %>%
        fEGarch(rt, n_test = 250, use_nonpar = TRUE)
      fcast <- predict_roll(est, parallel = FALSE, refit_after = 100)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(-0.0004435214, -0.0004435214, -0.0004435214, -0.0004435214, -0.0004435214,
        -0.0005152874, -0.0005152874, -0.0005152874, -0.0005152874, -0.0005152874,
        -0.0007104339, -0.0007104339, -0.0007104339, -0.0007104339, -0.0007104339,
        0.0065839069, 0.0060903241, 0.0068959571, 0.0073330375, 0.0078639467,
        0.0262324769, 0.026113558, 0.0256021325, 0.0276263001, 0.0271606675,
        0.0348847219, 0.0352535717, 0.0353867998, 0.0366814724, 0.0360883736)
    }, tolerance = 1e-02)

    expect_equal({
      est <- suppressWarnings(aparch(rt, n_test = 250, use_nonpar = TRUE))
      fcast <- predict_roll(est, parallel = FALSE, refit_after = 100)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(-0.0004435214, -0.0004435214, -0.0004435214, -0.0004435214, -0.0004435214,
        -0.0005152874, -0.0005152874, -0.0005152874, -0.0005152874, -0.0005152874,
        -0.0007104339, -0.0007104339, -0.0007104339, -0.0007104339, -0.0007104339,
        0.006929891, 0.0065812373, 0.0073409529, 0.0075926555, 0.0080062076,
        0.024379877, 0.0235379095, 0.0222659856, 0.0247630431, 0.0233514043,
        0.0268816542, 0.0270505628, 0.0267399069, 0.0285635716, 0.0273526311)
    }, tolerance = 1e-02)

    expect_equal({
      est <- suppressWarnings(fiaparch(rt, n_test = 250, use_nonpar = TRUE))
      fcast <- predict_roll(est, parallel = FALSE, refit_after = 100)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5),
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
      c(-0.0004435214, -0.0004435214, -0.0004435214, -0.0004435214, -0.0004435214,
        -0.0005152874, -0.0005152874, -0.0005152874, -0.0005152874, -0.0005152874,
        -0.0007104339, -0.0007104339, -0.0007104339, -0.0007104339, -0.0007104339,
        0.0073929834, 0.0067631773, 0.0073321801, 0.0083405422, 0.0088772214,
        0.0260863117, 0.0252633269, 0.0236010593, 0.0251350347, 0.0260385369,
        0.0290104597, 0.0285208587, 0.028766602, 0.0302529273, 0.0306257318)
    }, tolerance = 1e-02)

})

test_that("fEGarch rolling forecast output correctly formatted", {

  expect_s4_class(egarch_spec() %>% fEGarch(rt, n_test = 250) %>% predict_roll(), "fEGarch_forecast")
  expect_s4_class(loggarch_spec() %>% fEGarch(rt, n_test = 250) %>% predict_roll(), "fEGarch_forecast")
  expect_s4_class(suppressWarnings(aparch(rt, n_test = 250)) %>% predict_roll(), "fEGarch_forecast")
  expect_s4_class(suppressWarnings(fiaparch(rt, n_test = 250)) %>% predict_roll(), "fEGarch_forecast")

})

test_that("fEGarch rolling forecast works as intended for semiparametric models", {

    expect_equal({
      est <- egarch_spec() %>%
        fEGarch(rt, n_test = 250, use_nonpar = TRUE)
      fcast <- predict_roll(est, refit_after = 100, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
        c(0.00658390688108184, 0.00609032406826754, 0.00689595712887946,
        0.00733303750618234, 0.00786394670990562, 0.0262324769291835,
        0.0261135579790234, 0.0256021324543382,  0.0276263000537599,
        0.0271606675070345, 0.0348847219378428, 0.0352535716784695,
        0.0353867998187214, 0.0366814723704593, 0.0360883735619)
    }, tolerance = 1e-02)

    expect_equal({
      est <- egarch_spec() %>%
        fEGarch(rt, n_test = 250, use_nonpar = TRUE)
      fcast <- predict_roll(est, refit_after = NULL, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5)
      )
    }, {
        rep(-0.0004435214, 15)
    }, tolerance = 1e-06)

    expect_equal({
      est <- suppressWarnings(aparch(rt, n_test = 250, use_nonpar = TRUE))
      fcast <- predict_roll(est, refit_after = 100, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@sigt), 5),
        zoo::coredata(fcast@sigt)[101:105],
        utils::tail(zoo::coredata(fcast@sigt), 5)
      )
    }, {
        c(0.00692989096402155, 0.00658123733303244, 0.00734095285505457,
          0.00759265545800982, 0.0080062076387606,
          0.024379877019696, 0.0235379094738339, 0.0222659856219241,
          0.0247630431377694, 0.0233514042687065, 0.0268816541547651,
          0.0270505627951306, 0.0267399069011834, 0.0285635715857611,
          0.0273526311191019)
    }, tolerance = 1e-02)

    expect_equal({
      est <- suppressWarnings(aparch(rt, n_test = 250, use_nonpar = TRUE))
      fcast <- predict_roll(est, refit_after = NULL, parallel = FALSE)
      c(
        utils::head(zoo::coredata(fcast@cmeans), 5),
        zoo::coredata(fcast@cmeans)[101:105],
        utils::tail(zoo::coredata(fcast@cmeans), 5)
      )
    }, {
        rep(-0.0004435214, 15)
    }, tolerance = 1e-06)

})
