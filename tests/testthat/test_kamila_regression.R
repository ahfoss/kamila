local_edition(3)

# Helper function to generate deterministic dataset for regression tests
create_regression_dataset <- function() {
  set.seed(42)
  dat <- genMixedData(
    sampSize = 60,
    nConVar = 3,
    nCatVar = 3,
    nCatLevels = 4,
    nConWithErr = 1,
    nCatWithErr = 1,
    popProportions = c(0.5, 0.5),
    conErrLev = 0.2,
    catErrLev = 0.2
  )
  conDf <- data.frame(dat$conVars)
  catDf <- data.frame(lapply(data.frame(dat$catVars), factor))
  list(conDf = conDf, catDf = catDf)
}

test_that("1. Multivariate gold standard end-to-end regression (P=3, Q=3, K=2)", {
  d <- create_regression_dataset()
  set.seed(100)
  res <- kamila(
    conVar = d$conDf,
    catFactor = d$catDf,
    numClust = 2,
    numInit = 5,
    maxIter = 25,
    conInitMethod = "runif",
    catBw = 0.025,
    verbose = FALSE
  )

  expect_snapshot_value(res$finalMemb, style = "deparse", cran = TRUE)
  expect_snapshot_value(res$numIter, style = "deparse", cran = TRUE)
  expect_snapshot_value(res$finalLogLik, style = "deparse", cran = TRUE)
  expect_snapshot_value(res$finalObj, style = "deparse", cran = TRUE)
  expect_snapshot_value(res$finalCenters, style = "deparse", cran = TRUE)
  expect_snapshot_value(res$finalProbs, style = "deparse", cran = TRUE)

  # Classify holdout data
  set.seed(200)
  testCon <- data.frame(matrix(rnorm(15 * 3), nrow = 15, ncol = 3))
  testCat <- data.frame(
    X1 = factor(sample(1:4, 15, replace = TRUE), levels = levels(d$catDf[[1]])),
    X2 = factor(sample(1:4, 15, replace = TRUE), levels = levels(d$catDf[[2]])),
    X3 = factor(sample(1:4, 15, replace = TRUE), levels = levels(d$catDf[[3]]))
  )
  preds <- classifyKamila(res, list(testCon, testCat))
  expect_snapshot_value(preds, style = "deparse", cran = TRUE)
})

test_that("2. Feature weighting schemes regression", {
  d <- create_regression_dataset()

  # Non-uniform weights
  set.seed(100)
  res_weighted <- kamila(
    conVar = d$conDf,
    catFactor = d$catDf,
    numClust = 2,
    numInit = 5,
    conWeights = c(0.2, 0.8, 1.0),
    catWeights = c(0.4, 0.9, 0.1),
    maxIter = 25,
    conInitMethod = "runif",
    catBw = 0.025
  )

  expect_snapshot_value(res_weighted$finalMemb, style = "deparse", cran = TRUE)
  expect_snapshot_value(res_weighted$finalCenters, style = "deparse", cran = TRUE)
  expect_snapshot_value(res_weighted$finalProbs, style = "deparse", cran = TRUE)
  expect_snapshot_value(res_weighted$finalLogLik, style = "deparse", cran = TRUE)
  expect_snapshot_value(res_weighted$finalObj, style = "deparse", cran = TRUE)

  # Zero weights
  set.seed(100)
  res_zero <- kamila(
    conVar = d$conDf,
    catFactor = d$catDf,
    numClust = 2,
    numInit = 5,
    conWeights = c(1.0, 1.0, 0.0),
    catWeights = c(1.0, 1.0, 0.0),
    maxIter = 25,
    conInitMethod = "runif",
    catBw = 0.025
  )

  expect_snapshot_value(res_zero$finalMemb, style = "deparse", cran = TRUE)
  expect_snapshot_value(res_zero$finalCenters, style = "deparse", cran = TRUE)
  expect_snapshot_value(res_zero$finalLogLik, style = "deparse", cran = TRUE)
  expect_snapshot_value(res_zero$finalObj, style = "deparse", cran = TRUE)
})

test_that("3. Initialization methods and categorical bandwidth variations", {
  d <- create_regression_dataset()

  # conInitMethod = "sample"
  set.seed(300)
  res_sample <- kamila(
    conVar = d$conDf,
    catFactor = d$catDf,
    numClust = 2,
    numInit = 5,
    conInitMethod = "sample",
    maxIter = 25
  )
  expect_snapshot_value(res_sample$finalMemb, style = "deparse", cran = TRUE)
  expect_snapshot_value(res_sample$finalCenters, style = "deparse", cran = TRUE)
  expect_snapshot_value(res_sample$finalLogLik, style = "deparse", cran = TRUE)

  # catBw = 0.005 (light smoothing)
  set.seed(400)
  res_bw0 <- kamila(
    conVar = d$conDf,
    catFactor = d$catDf,
    numClust = 2,
    numInit = 5,
    catBw = 0.005,
    maxIter = 25
  )
  expect_snapshot_value(res_bw0$finalMemb, style = "deparse", cran = TRUE)
  expect_snapshot_value(res_bw0$finalProbs, style = "deparse", cran = TRUE)
  expect_snapshot_value(res_bw0$finalLogLik, style = "deparse", cran = TRUE)

  # catBw = 0.1 (strong smoothing)
  set.seed(400)
  res_bw1 <- kamila(
    conVar = d$conDf,
    catFactor = d$catDf,
    numClust = 2,
    numInit = 5,
    catBw = 0.1,
    maxIter = 25
  )
  expect_snapshot_value(res_bw1$finalMemb, style = "deparse", cran = TRUE)
  expect_snapshot_value(res_bw1$finalProbs, style = "deparse", cran = TRUE)
  expect_snapshot_value(res_bw1$finalLogLik, style = "deparse", cran = TRUE)
})

test_that("4. Prediction strength exact numerical metrics", {
  d <- create_regression_dataset()
  set.seed(500)
  res_ps <- kamila(
    conVar = d$conDf,
    catFactor = d$catDf,
    numClust = 2:3,
    numInit = 3,
    maxIter = 15,
    calcNumClust = "ps",
    numPredStrCvRun = 5,
    predStrThresh = 0.65
  )

  expect_snapshot_value(res_ps$nClust$bestNClust, style = "deparse", cran = TRUE)
  expect_snapshot_value(res_ps$nClust$psValues, style = "deparse", cran = TRUE)
  expect_snapshot_value(res_ps$nClust$avgPredStr, style = "deparse", cran = TRUE)
  expect_snapshot_value(res_ps$nClust$stdErrPredStr, style = "deparse", cran = TRUE)
  expect_snapshot_value(res_ps$nClust$psCvRes, style = "deparse", cran = TRUE)
})

test_that("5. Convergence dynamics, iteration limits, verbose traces, and degenerate handling", {
  d <- create_regression_dataset()

  # maxIter = 2
  set.seed(600)
  res_iter2 <- kamila(
    conVar = d$conDf,
    catFactor = d$catDf,
    numClust = 2,
    numInit = 3,
    maxIter = 2
  )
  expect_snapshot_value(res_iter2$numIter, style = "deparse", cran = TRUE)
  expect_snapshot_value(res_iter2$finalMemb, style = "deparse", cran = TRUE)

  # verbose = TRUE
  set.seed(700)
  res_verb <- kamila(
    conVar = d$conDf,
    catFactor = d$catDf,
    numClust = 2,
    numInit = 2,
    maxIter = 5,
    verbose = TRUE
  )
  expect_snapshot_value(res_verb$verbose$totalLogLikVect, style = "deparse", cran = TRUE)
  expect_snapshot_value(res_verb$verbose$winDistVect, style = "deparse", cran = TRUE)
  expect_snapshot_value(res_verb$verbose$totalDist, style = "deparse", cran = TRUE)
  expect_snapshot_value(res_verb$verbose$objectiveVect, style = "deparse", cran = TRUE)
  expect_snapshot_value(res_verb$verbose$membLongList, style = "deparse", cran = TRUE)

  # Degenerate cluster handling (empty cluster creation triggers fallback)
  set.seed(800)
  conDeg <- data.frame(x = c(rep(0, 10), 100, 200), y = c(rep(0, 10), 100, 200))
  catDeg <- data.frame(f = factor(c(rep("A", 10), "B", "C")))
  res_deg <- kamila(
    conVar = conDeg,
    catFactor = catDeg,
    numClust = 4,
    numInit = 5,
    maxIter = 5
  )
  expect_snapshot_value(res_deg$finalMemb, style = "deparse", cran = TRUE)
  expect_snapshot_value(res_deg$finalCenters, style = "deparse", cran = TRUE)
})

test_that("6. Isolated low-level Rcpp kernel exact numerical outputs", {
  # dptmCpp
  pts <- matrix(c(1.0, 2.0, 3.0, 4.0, 5.0, 6.0), nrow = 3, ncol = 2)
  means <- matrix(c(1.0, 2.0, 2.0, 3.0), nrow = 2, ncol = 2)
  wgts <- c(0.5, 2.0)
  d_val <- dptmCpp(pts, means, wgts)
  expect_snapshot_value(d_val, style = "deparse", cran = TRUE)

  # rowMin, rowMax, rowMaxInds
  test_mat <- matrix(c(-1.5, 3.2, 0.0, 4.1, 3.2, -5.0, 4.1, 0.0, 2.2), nrow = 3, byrow = TRUE)
  expect_snapshot_value(rowMin(test_mat), style = "deparse", cran = TRUE)
  expect_snapshot_value(rowMax(test_mat), style = "deparse", cran = TRUE)
  expect_snapshot_value(rowMaxInds(test_mat), style = "deparse", cran = TRUE)

  # radialKDE
  radii <- c(0.5, 1.2, 1.8, 2.3, 3.0, 3.5, 4.1, 5.0)
  eval_pts <- c(0.0, 1.0, 2.5, 4.0, 6.0)
  kde_res <- radialKDE(radii = radii, evalPoints = eval_pts, pdim = 3, returnFun = FALSE)
  expect_snapshot_value(kde_res$kdes, style = "deparse", cran = TRUE)

  # jointTabSmoothedList & getIndividualLogProbs
  catNum <- matrix(c(1, 2, 1, 2, 1, 1, 2, 2, 3, 1), nrow = 5, ncol = 2)
  memb <- c(1, 1, 2, 2, 2)
  numLev <- c(2, 3)
  j_tabs <- jointTabSmoothedList(catNum, memb, numLev, catBw = 0.05, kk = 2)
  expect_snapshot_value(j_tabs, style = "deparse", cran = TRUE)

  logProbsCond <- lapply(j_tabs, function(xx) log(xx / rowSums(xx)))
  ind_log_probs <- getIndividualLogProbs(catNum, c(0.5, 1.0), logProbsCond)
  expect_snapshot_value(ind_log_probs, style = "deparse", cran = TRUE)
})
