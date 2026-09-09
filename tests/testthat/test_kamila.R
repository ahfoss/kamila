library(kamila)

test_that("KAMILA runs as expected, small data set", {
  suppressWarnings(withr::local_rng_version("3.5.0"))
  set.seed(1)
  res1 <- kamila(
    conVar = data.frame(rnorm(15)),
    catFactor = data.frame(factor(sample(1:4, size = 15, rep = TRUE))),
    numClust = 2,
    numInit = 10,
    maxIter = 25,
    conInitMethod = "runif",
    catBw = 0.025,
    verbose = FALSE,
    calcNumClust = "none",
    numPredStrCvRun = 20,
    predStrThresh = 0.8
  )
  expect_identical(res1$finalMemb, c(2, 1, 2, 1, 1, 2, 1, 1, 1, 2, 1, 1, 2, 2, 1))
})

test_that("KAMILA warns if pred-strength is used with a single numClust", {
  expect_warning(
    kamila(
      conVar = data.frame(rnorm(10), rnorm(10)),
      catFactor = data.frame(
        factor(sample(1:4, size = 10, rep = TRUE)),
        factor(sample(1:4, size = 10, rep = TRUE))
      ),
      numClust = 2,
      numInit = 2,
      maxIter = 2,
      conInitMethod = "runif",
      catBw = 0.025,
      verbose = FALSE,
      calcNumClust = "ps",
      numPredStrCvRun = 2,
      predStrThresh = 0.8
    ),
    "Input parameter numClust is a scalar; the prediction strength"
  )
})

test_that('KAMILA throws error if numClust is length > 1 with calcNumClust=="none"', {
  expect_error(
    kamila(
      conVar = data.frame(rnorm(10), rnorm(10)),
      catFactor = data.frame(
        factor(sample(1:4, size = 10, rep = TRUE)),
        factor(sample(1:4, size = 10, rep = TRUE))
      ),
      numClust = 2:5,
      numInit = 2,
      maxIter = 2,
      conInitMethod = "runif",
      catBw = 0.025,
      verbose = FALSE,
      calcNumClust = "none",
      numPredStrCvRun = 2,
      predStrThresh = 0.8
    ),
    'Input parameter numClust must be length 1 if calcNumClust == "none"'
  )
})

test_that("KAMILA prediction strength works with single-variable catFactor data frame (Issue #14)", {
  set.seed(123)
  res <- kamila(
    conVar = data.frame(rnorm(20), rnorm(20)),
    catFactor = data.frame(factor(sample(1:3, size = 20, replace = TRUE))),
    numClust = 2:3,
    numInit = 3,
    maxIter = 10,
    calcNumClust = "ps",
    numPredStrCvRun = 2
  )
  expect_true(is.list(res))
  expect_true(res$nClust$bestNClust %in% 2:3)
})

test_that("radialKDE and KAMILA handle distance 0 without producing -Inf (Issue #9)", {
  withr::with_seed(42, {
    # 1. Direct radialKDE test with distance 0 across various pdim
    radii <- c(0.1, 0.5, 1.2, 1.8, 2.5)
    for (p in 1:3) {
      rkde <- radialKDE(radii = radii, evalPoints = c(0.0, 0.5, 1.0), pdim = p)
      expect_true(all(rkde$kdes > 0))
      expect_true(all(is.finite(log(rkde$kdes))))
      expect_true(rkde$kdes[1] > 0)
    }

    # 2. kamila with duplicate continuous rows does not yield -Inf log-likelihood
    con_dup <- data.frame(v1 = c(1, 1, 5, 5), v2 = c(1, 1, 5, 5))
    cat_dup <- data.frame(cat = factor(c(1, 1, 2, 2)))
    res_dup <- kamila(
      conVar = con_dup,
      catFactor = cat_dup,
      numClust = 2,
      numInit = 1,
      conInitMethod = "sample"
    )
    expect_true(is.finite(res_dup$finalLogLik))
    expect_equal(length(unique(res_dup$finalMemb)), 2)

    # 3. classifyKamila correctly classifies a point sitting exactly on a cluster centroid
    dat <- genMixedData(
      sampSize = 40,
      nConVar = 2,
      nCatVar = 2,
      nCatLevels = 2,
      nConWithErr = 0,
      nCatWithErr = 0,
      popProportions = c(0.5, 0.5),
      conErrLev = 0.05,
      catErrLev = 0.05
    )
    con_df <- data.frame(dat$conVars)
    cat_df <- data.frame(lapply(as.data.frame(dat$catVars), factor))
    res_model <- kamila(conVar = con_df, catFactor = cat_df, numClust = 2, numInit = 3)

    obs_idx <- which(res_model$finalMemb == 1)[1]
    center1 <- res_model$finalCenters[1, ]
    test_con <- data.frame(X1 = center1[1], X2 = center1[2])
    test_cat <- cat_df[obs_idx, , drop = FALSE]
    pred <- classifyKamila(res_model, newData = list(test_con, test_cat))
    expect_equal(pred, 1)

    # 4. classifyKamila with continuous-only model at centroid
    res_con <- kamila(conVar = con_df, numClust = 2, numInit = 3)
    center1_con <- res_con$finalCenters[1, ]
    pred_con <- classifyKamila(res_con, data.frame(X1 = center1_con[1], X2 = center1_con[2]))
    expect_equal(pred_con, 1)
  })
})

test_that("KAMILA throws clear errors when inputs contain NA values (Issue #3)", {
  conDf <- data.frame(x = rnorm(10), y = rnorm(10))
  catDf <- data.frame(f = factor(rep(c("A", "B"), 5)))

  conDf_na <- conDf
  conDf_na[1, 1] <- NA
  expect_error(
    kamila(conVar = conDf_na, catFactor = catDf, numClust = 2, numInit = 2),
    "conVar contains missing values \\(NA\\)"
  )

  catDf_na <- catDf
  catDf_na[1, 1] <- NA
  expect_error(
    kamila(conVar = conDf, catFactor = catDf_na, numClust = 2, numInit = 2),
    "catFactor contains missing values \\(NA\\)"
  )
})
