test_that("initMeans and dptmCpp validation and sample method work", {
  conVar <- data.frame(x = 1:10, y = 11:20, stringsAsFactors = TRUE)

  # initMeans method = 'sample'
  m_sample <- initMeans(conVar, method = "sample", numClust = 3)
  expect_equal(dim(m_sample), c(3, 2))

  # initMeans invalid method error
  expect_error(initMeans(conVar, method = "invalid", numClust = 2), "Unrecognized mean initialization method")

  # dptmCpp validation errors
  means <- matrix(1:4, nrow = 2)
  expect_error(dptmCpp(conVar, means[, 1, drop = FALSE], c(1, 1)), "Dimensionality of pts and myMeans must be equal")
  expect_error(dptmCpp(conVar, means, c(1)), "Dimensionality of pts must equal number of weights")
})

test_that("kamila error handling and verbose mode work", {
  conVar <- data.frame(x = rnorm(20), y = rnorm(20), stringsAsFactors = TRUE)
  catFactor <- data.frame(f = factor(rep(c("A", "B"), 10)), stringsAsFactors = TRUE)

  # numClust length != 1 when calcNumClust == 'none'
  expect_error(kamila(conVar, catFactor, numClust = c(2, 3), numInit = 2), "Input parameter numClust must be length 1")

  # Neither dataset specified
  expect_error(kamila(numClust = 2, numInit = 2), "At least one of conVar or catFactor must be specified.")

  # Non-dataframe input
  expect_error(
    kamila(matrix(1:20), catFactor, numClust = 2, numInit = 2),
    "Input dataset conVar must be a dataframe."
  )
  expect_error(
    kamila(conVar, factor(1:20), numClust = 2, numInit = 2),
    "Input dataset catFactor must be a dataframe."
  )

  # 0-column dataframe input
  expect_error(
    kamila(conVar[, 0, drop = FALSE], numClust = 2, numInit = 2),
    "Input dataset conVar must have at least 1 column"
  )
  expect_error(
    kamila(catFactor = catFactor[, 0, drop = FALSE], numClust = 2, numInit = 2),
    "Input dataset catFactor must have at least 1 column"
  )

  # Weight length mismatch
  expect_error(
    kamila(conVar, catFactor, numClust = 2, numInit = 2, conWeights = 1),
    "Length of conWeights must equal number of continuous variables"
  )
  expect_error(
    kamila(conVar, catFactor, numClust = 2, numInit = 2, catWeights = c(1, 1)),
    "Length of catWeights must equal number of categorical variables"
  )

  # Invalid weights
  expect_error(
    kamila(conVar, catFactor, numClust = 2, numInit = 2, conWeights = c(1.5, 0.5)),
    "Weights must be in \\[0,1\\]"
  )
  expect_error(
    kamila(conVar, catFactor, numClust = 2, numInit = 2, conWeights = c(-0.1, 0.5)),
    "Weights must be in \\[0,1\\]"
  )
  expect_error(
    kamila(conVar, catFactor, numClust = 2, numInit = 2, catWeights = 1.5),
    "Weights must be in \\[0,1\\]"
  )
  expect_error(
    kamila(conVar, catFactor, numClust = 2, numInit = 2, catWeights = -0.1),
    "Weights must be in \\[0,1\\]"
  )

  # Mismatched row counts
  expect_error(kamila(conVar, catFactor[1:10, , drop = FALSE], numClust = 2, numInit = 2), "don't match")

  # Invalid calcNumClust option
  expect_error(
    kamila(conVar, catFactor, numClust = 2, numInit = 2, calcNumClust = "invalid"),
    "must be either \"none\" or \"ps\""
  )

  # Verbose mode
  kam_verb <- kamila(conVar, catFactor, numClust = 2, numInit = 2, verbose = TRUE)
  expect_named(
    kam_verb$verbose,
    c("totalLogLikVect", "catLogLikVect", "winDistVect", "totalDist", "objectiveVect", "membLongList")
  )
})

test_that("kamila with prediction strength (calcNumClust == 'ps') works and validates parameters", {
  conVar <- data.frame(x = rnorm(40), y = rnorm(40), stringsAsFactors = TRUE)
  catFactor <- data.frame(f1 = factor(rep(c("A", "B"), 20)), stringsAsFactors = TRUE)

  # Invalid numClust inputs
  expect_error(
    kamila(conVar, catFactor, numClust = c(2, 2), numInit = 2, calcNumClust = "ps"),
    "vector of unique integers"
  )
  expect_error(kamila(conVar, catFactor, numClust = c(2, 30), numInit = 2, calcNumClust = "ps"), "cannot exceed")

  # Scalar numClust warning
  expect_warning(kamila(conVar, catFactor, numClust = 2, numInit = 2, calcNumClust = "ps"), "numClust is a scalar")

  # Invalid predStrThresh
  expect_error(
    kamila(conVar, catFactor, numClust = c(2, 3), numInit = 2, calcNumClust = "ps", predStrThresh = 1.5),
    "must be scalar in \\(0,1\\)"
  )

  # Invalid numPredStrCvRun
  expect_error(
    kamila(conVar, catFactor, numClust = c(2, 3), numInit = 2, calcNumClust = "ps", numPredStrCvRun = 0),
    "must be a positive integer"
  )

  # High threshold warning
  expect_warning(
    kamila(
      conVar, catFactor, numClust = c(2, 3), numInit = 2, calcNumClust = "ps",
      numPredStrCvRun = 2, predStrThresh = 0.999
    ),
    "No cluster size is above prediction strength threshold"
  )

  # Successful PS run
  ps_res <- kamila(
    conVar, catFactor, numClust = c(2, 3), numInit = 2, calcNumClust = "ps",
    numPredStrCvRun = 2, predStrThresh = 0.5
  )
  expect_true(ps_res$nClust$bestNClust %in% c(2, 3))
})

test_that("classifyKamila works and validates inputs", {
  conVar <- data.frame(x = rnorm(30), y = rnorm(30), stringsAsFactors = TRUE)
  catFactor <- data.frame(f1 = factor(rep(c("A", "B"), 15)), stringsAsFactors = TRUE)
  kamObj <- kamila(conVar, catFactor, numClust = 2, numInit = 2)

  # Valid classification
  newData <- list(conVar[1:5, ], catFactor[1:5, , drop = FALSE])
  pred <- classifyKamila(kamObj, newData)
  expect_equal(length(pred), 5)
  expect_true(all(pred %in% c(1, 2)))

  # Invalid list length error
  expect_error(classifyKamila(kamObj, list(conVar[1:5, ])), "must be list of length 2")
  expect_error(classifyKamila(kamObj, list()), "newData list must have length 1 or 2")
  expect_error(
    classifyKamila(kamObj, list(conVar[1:5, ], catFactor[1:5, , drop = FALSE], conVar[1:5, ])),
    "newData list must have length 1 or 2"
  )

  # Dataframe passed to mixed model error
  expect_error(classifyKamila(kamObj, conVar[1:5, ]), "must be list of length 2 for mixed data")

  # Non-dataframe/non-list input error
  expect_error(classifyKamila(kamObj, "invalid"), "must be a data frame or list of data frames")
  expect_error(classifyKamila(kamObj, 1:10), "must be a data frame or list of data frames")

  # Non-dataframe element error
  expect_error(
    classifyKamila(kamObj, list(as.matrix(conVar[1:5, ]), catFactor[1:5, , drop = FALSE])),
    "elements of newData must be data frames"
  )
  expect_error(
    classifyKamila(kamObj, list(conVar[1:5, ], as.matrix(catFactor[1:5, , drop = FALSE]))),
    "elements of newData must be data frames"
  )

  # 0-column dataframe error
  expect_error(
    classifyKamila(kamObj, list(conVar[1:5, 0, drop = FALSE], catFactor[1:5, , drop = FALSE])),
    "must have at least 1 column"
  )
  expect_error(
    classifyKamila(kamObj, list(conVar[1:5, ], catFactor[1:5, 0, drop = FALSE])),
    "must have at least 1 column"
  )

  # Row mismatch error
  expect_error(
    classifyKamila(kamObj, list(conVar[1:5, ], catFactor[1:3, , drop = FALSE])),
    "number of observations in con and cat vars don't match"
  )
})

test_that("myCatKern and sumMatList Rcpp helper function work", {
  # myCatKern
  kdat <- data.frame(v1 = factor(c(1, 2, 1, 2)), v2 = factor(c(1, 1, 2, 2)), stringsAsFactors = TRUE)
  ck_tab <- myCatKern(kdat, bw = 0.1, tabOnly = TRUE)
  expect_equal(dim(ck_tab), c(2, 2))

  ck_pred <- myCatKern(kdat, bw = 0.1, tabOnly = FALSE)
  expect_named(ck_pred, c("preds", "tab"))
  expect_equal(length(ck_pred$preds), 4)

  # sumMatList
  m1 <- matrix(1:4, nrow = 2)
  m2 <- matrix(5:8, nrow = 2)
  s_mat <- sumMatList(list(m1, m2))
  expect_equal(s_mat, m1 + m2)
})

test_that("kamila PS handles small cluster size fallback (clustN < 2)", {
  set.seed(123)
  conVar <- data.frame(x = c(0, 0.1, 0.2, 10, 10.1, 10.2, 20, 20.1))
  catFactor <- data.frame(f = factor(c("A", "A", "A", "B", "B", "B", "C", "C")))

  # numClust = c(2, 4) with 8 observations often generates a test cluster of size 1
  ps_small <- suppressWarnings(
    kamila(
      conVar, catFactor, numClust = c(2, 4), numInit = 2, calcNumClust = "ps",
      numPredStrCvRun = 2, predStrThresh = 0.5
    )
  )
  expect_true(!is.null(ps_small$finalMemb))
})

test_that("KAMILA works natively with continuous-only data", {
  set.seed(42)
  conDf <- data.frame(
    x = c(rnorm(25, mean = 0), rnorm(25, mean = 5)),
    y = c(rnorm(25, mean = 0), rnorm(25, mean = 5))
  )

  # Run continuous-only KAMILA
  res_con <- kamila(conVar = conDf, numClust = 2, numInit = 5, maxIter = 15)
  expect_equal(length(res_con$finalMemb), 50)
  expect_true(all(res_con$finalMemb %in% 1:2))
  expect_equal(dim(res_con$finalCenters), c(2, 2))
  expect_equal(length(res_con$finalProbs), 0)
  expect_true(is.numeric(res_con$finalLogLik))
  expect_equal(res_con$finalObj, res_con$finalLogLik)

  # classifyKamila with data frame input
  pred_df <- classifyKamila(res_con, conDf[1:10, ])
  expect_equal(length(pred_df), 10)
  expect_true(all(pred_df %in% 1:2))

  # classifyKamila with list of length 1
  pred_list <- classifyKamila(res_con, list(conDf[1:10, ]))
  expect_equal(pred_list, pred_df)

  # classifyKamila with list of length 2 (second element NULL)
  pred_list2 <- classifyKamila(res_con, list(conDf[1:10, ], NULL))
  expect_equal(pred_list2, pred_df)

  # classifyKamila errors on column mismatch
  expect_error(
    classifyKamila(res_con, conDf[1:5, 1, drop = FALSE]),
    "number of continuous columns in newData does not match model"
  )

  # Continuous-only prediction strength (ps)
  ps_con <- suppressWarnings(
    kamila(
      conVar = conDf, numClust = 2:3, numInit = 2, maxIter = 10,
      calcNumClust = "ps", numPredStrCvRun = 2, predStrThresh = 0.5
    )
  )
  expect_true(ps_con$nClust$bestNClust %in% 2:3)

  # calcApproxBIC handles continuous-only results
  bic_con <- calcApproxBIC(res_con)
  expect_true(is.numeric(bic_con$criteria))
  expect_equal(bic_con$nCatParm, 0)
  expect_equal(bic_con$nConParm, 4)
})

test_that("KAMILA works natively with categorical-only data", {
  set.seed(42)
  catDf <- data.frame(
    v1 = factor(c(rep("A", 25), rep("B", 25))),
    v2 = factor(c(rep("X", 25), rep("Y", 25)))
  )

  # Run categorical-only KAMILA
  res_cat <- kamila(catFactor = catDf, numClust = 2, numInit = 5, maxIter = 15)
  expect_equal(length(res_cat$finalMemb), 50)
  expect_true(all(res_cat$finalMemb %in% 1:2))
  expect_null(res_cat$finalCenters)
  expect_equal(length(res_cat$finalProbs), 2)
  expect_true(is.numeric(res_cat$finalLogLik))
  expect_equal(res_cat$finalObj, res_cat$finalLogLik)

  # classifyKamila with data frame input
  pred_df <- classifyKamila(res_cat, catDf[1:10, ])
  expect_equal(length(pred_df), 10)
  expect_true(all(pred_df %in% 1:2))

  # classifyKamila with list of length 1
  pred_list <- classifyKamila(res_cat, list(catDf[1:10, ]))
  expect_equal(pred_list, pred_df)

  # classifyKamila with list of length 2 (first element NULL)
  pred_list2 <- classifyKamila(res_cat, list(NULL, catDf[1:10, ]))
  expect_equal(pred_list2, pred_df)

  # classifyKamila errors on column mismatch
  expect_error(
    classifyKamila(res_cat, catDf[1:5, 1, drop = FALSE]),
    "number of categorical columns in newData does not match model"
  )

  # Categorical-only prediction strength (ps)
  ps_cat <- suppressWarnings(
    kamila(
      catFactor = catDf, numClust = 2:3, numInit = 2, maxIter = 10,
      calcNumClust = "ps", numPredStrCvRun = 2, predStrThresh = 0.5
    )
  )
  expect_true(ps_cat$nClust$bestNClust %in% 2:3)

  # calcApproxBIC handles categorical-only results
  bic_cat <- calcApproxBIC(res_cat)
  expect_true(is.numeric(bic_cat$criteria))
  expect_equal(bic_cat$nConParm, 0)
  expect_true(bic_cat$nCatParm > 0)
})
