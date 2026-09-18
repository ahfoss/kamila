test_that("kamila validates numCores input parameter", {
  conVar <- data.frame(
    x = c(rep(0, 15), rep(10, 15)),
    y = c(rep(0, 15), rep(10, 15)),
    stringsAsFactors = TRUE
  )
  catFactor <- data.frame(
    f = factor(rep(c("A", "B"), each = 15)),
    stringsAsFactors = TRUE
  )

  # Valid integer numCores = 1 does not throw error
  expect_error(
    kamila(
      conVar, catFactor,
      numClust = c(2, 3), numInit = 2,
      calcNumClust = "ps", numPredStrCvRun = 2,
      numCores = 1
    ),
    NA
  )

  # numCores <= 0
  expect_error(
    kamila(
      conVar, catFactor,
      numClust = c(2, 3), numInit = 2,
      calcNumClust = "ps", numPredStrCvRun = 2,
      numCores = 0
    ),
    "Input parameter numCores must be a positive integer or a cluster object."
  )
  expect_error(
    kamila(
      conVar, catFactor,
      numClust = c(2, 3), numInit = 2,
      calcNumClust = "ps", numPredStrCvRun = 2,
      numCores = -2
    ),
    "Input parameter numCores must be a positive integer or a cluster object."
  )

  # numCores non-integer numeric
  expect_error(
    kamila(
      conVar, catFactor,
      numClust = c(2, 3), numInit = 2,
      calcNumClust = "ps", numPredStrCvRun = 2,
      numCores = 1.5
    ),
    "Input parameter numCores must be a positive integer or a cluster object."
  )

  # numCores invalid type
  expect_error(
    kamila(
      conVar, catFactor,
      numClust = c(2, 3), numInit = 2,
      calcNumClust = "ps", numPredStrCvRun = 2,
      numCores = "two"
    ),
    "Input parameter numCores must be a positive integer or a cluster object."
  )

  # numCores vector length > 1
  expect_error(
    kamila(
      conVar, catFactor,
      numClust = c(2, 3), numInit = 2,
      calcNumClust = "ps", numPredStrCvRun = 2,
      numCores = c(1, 2)
    ),
    "Input parameter numCores must be a positive integer or a cluster object."
  )

  # numCores NA / NaN
  expect_error(
    kamila(
      conVar, catFactor,
      numClust = c(2, 3), numInit = 2,
      calcNumClust = "ps", numPredStrCvRun = 2,
      numCores = NA
    ),
    "Input parameter numCores must be a positive integer or a cluster object."
  )
  expect_error(
    kamila(
      conVar, catFactor,
      numClust = c(2, 3), numInit = 2,
      calcNumClust = "ps", numPredStrCvRun = 2,
      numCores = NaN
    ),
    "Input parameter numCores must be a positive integer or a cluster object."
  )

  # numCores is ignored when calcNumClust == 'none'
  expect_error(
    kamila(
      conVar, catFactor,
      numClust = 2, numInit = 2,
      calcNumClust = "none",
      numCores = -5
    ),
    NA
  )
})

test_that("parallel prediction strength runs produce valid output with numCores = 2", {
  dat <- genMixedData(
    80,
    nConVar = 2, nCatVar = 2, nCatLevels = 4,
    nConWithErr = 2, nCatWithErr = 2,
    popProportions = c(0.5, 0.5),
    conErrLev = 0.3, catErrLev = 0.8
  )
  catDf <- data.frame(apply(dat$catVars, 2, factor), stringsAsFactors = TRUE)
  conDf <- data.frame(scale(dat$conVars), stringsAsFactors = TRUE)

  res_par <- kamila(
    conVar = conDf,
    catFactor = catDf,
    numClust = 2:3,
    numInit = 2,
    maxIter = 10,
    calcNumClust = "ps",
    numPredStrCvRun = 4,
    predStrThresh = 0.5,
    numCores = 2
  )

  expect_type(res_par, "list")
  expect_true(res_par$nClust$bestNClust %in% 2:3)
  expect_equal(dim(res_par$nClust$psCvRes), c(2, 4))
  expect_equal(rownames(res_par$nClust$psCvRes), c("2", "3"))
  expect_equal(colnames(res_par$nClust$psCvRes), paste("Run", 1:4))
  expect_equal(names(res_par$nClust$psValues), c("2", "3"))
  expect_equal(length(res_par$nClust$avgPredStr), 2)
  expect_equal(length(res_par$nClust$stdErrPredStr), 2)
  expect_true(all(!is.na(res_par$nClust$psCvRes)))
})

test_that("parallel prediction strength supports user-provided cluster object", {
  dat <- genMixedData(
    60,
    nConVar = 2, nCatVar = 2, nCatLevels = 4,
    nConWithErr = 2, nCatWithErr = 2,
    popProportions = c(0.5, 0.5),
    conErrLev = 0.3, catErrLev = 0.8
  )
  catDf <- data.frame(apply(dat$catVars, 2, factor), stringsAsFactors = TRUE)
  conDf <- data.frame(scale(dat$conVars), stringsAsFactors = TRUE)

  cl <- parallel::makeCluster(2)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  res_cl <- kamila(
    conVar = conDf,
    catFactor = catDf,
    numClust = 2:3,
    numInit = 2,
    maxIter = 10,
    calcNumClust = "ps",
    numPredStrCvRun = 2,
    predStrThresh = 0.5,
    numCores = cl
  )

  expect_true(res_cl$nClust$bestNClust %in% 2:3)
  expect_equal(dim(res_cl$nClust$psCvRes), c(2, 2))

  # Verify cluster is still active (not stopped by kamila)
  test_eval <- parallel::parLapply(cl, 1:2, function(x) x * 2)
  expect_equal(test_eval, list(2, 4))
})

test_that("parallel prediction strength is reproducible with set.seed()", {
  dat <- genMixedData(
    60,
    nConVar = 2, nCatVar = 2, nCatLevels = 4,
    nConWithErr = 2, nCatWithErr = 2,
    popProportions = c(0.5, 0.5),
    conErrLev = 0.3, catErrLev = 0.8
  )
  catDf <- data.frame(apply(dat$catVars, 2, factor), stringsAsFactors = TRUE)
  conDf <- data.frame(scale(dat$conVars), stringsAsFactors = TRUE)

  set.seed(999)
  res1 <- kamila(
    conVar = conDf,
    catFactor = catDf,
    numClust = 2:3,
    numInit = 2,
    maxIter = 10,
    calcNumClust = "ps",
    numPredStrCvRun = 2,
    predStrThresh = 0.5,
    numCores = 2
  )

  set.seed(999)
  res2 <- kamila(
    conVar = conDf,
    catFactor = catDf,
    numClust = 2:3,
    numInit = 2,
    maxIter = 10,
    calcNumClust = "ps",
    numPredStrCvRun = 2,
    predStrThresh = 0.5,
    numCores = 2
  )

  expect_equal(res1$nClust$psCvRes, res2$nClust$psCvRes)
  expect_equal(res1$nClust$psValues, res2$nClust$psValues)
  expect_equal(res1$nClust$bestNClust, res2$nClust$bestNClust)
})

test_that("parallel prediction strength works with numPredStrCvRun = 1 and effectiveCores = 1", {
  dat <- genMixedData(
    60,
    nConVar = 2, nCatVar = 2, nCatLevels = 4,
    nConWithErr = 2, nCatWithErr = 2,
    popProportions = c(0.5, 0.5),
    conErrLev = 0.3, catErrLev = 0.8
  )
  catDf <- data.frame(apply(dat$catVars, 2, factor), stringsAsFactors = TRUE)
  conDf <- data.frame(scale(dat$conVars), stringsAsFactors = TRUE)

  # numPredStrCvRun = 1 with numCores = 2 falls back to sequential without starting a cluster
  res_single_run <- kamila(
    conVar = conDf,
    catFactor = catDf,
    numClust = 2:3,
    numInit = 2,
    maxIter = 10,
    calcNumClust = "ps",
    numPredStrCvRun = 1,
    predStrThresh = 0.5,
    numCores = 2
  )

  expect_equal(dim(res_single_run$nClust$psCvRes), c(2, 1))
  expect_true(res_single_run$nClust$bestNClust %in% 2:3)
})

test_that("parallel prediction strength works with continuous-only and categorical-only data", {
  dat <- genMixedData(
    60,
    nConVar = 2, nCatVar = 2, nCatLevels = 4,
    nConWithErr = 2, nCatWithErr = 2,
    popProportions = c(0.5, 0.5),
    conErrLev = 0.3, catErrLev = 0.8
  )
  catDf <- data.frame(apply(dat$catVars, 2, factor), stringsAsFactors = TRUE)
  conDf <- data.frame(scale(dat$conVars), stringsAsFactors = TRUE)

  # Continuous-only in parallel
  res_con <- kamila(
    conVar = conDf,
    numClust = 2:3,
    numInit = 2,
    maxIter = 10,
    calcNumClust = "ps",
    numPredStrCvRun = 2,
    predStrThresh = 0.5,
    numCores = 2
  )
  expect_equal(dim(res_con$nClust$psCvRes), c(2, 2))

  # Categorical-only in parallel
  res_cat <- kamila(
    catFactor = catDf,
    numClust = 2:3,
    numInit = 2,
    maxIter = 10,
    calcNumClust = "ps",
    numPredStrCvRun = 2,
    predStrThresh = 0.5,
    numCores = 2
  )
  expect_equal(dim(res_cat$nClust$psCvRes), c(2, 2))
})
