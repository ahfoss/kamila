test_that("dummyCodeOneVar and dummyCodeFactorDf work as expected", {
  fac <- factor(c("a", "b", "a", "c"))
  d1 <- dummyCodeOneVar(fac)
  expect_equal(dim(d1), c(4, 3))
  expect_equal(colnames(d1), c("a", "b", "c"))
  expect_equal(unname(d1[1, ]), c(1, 0, 0))

  df <- data.frame(v1 = factor(c("x", "y")), v2 = factor(c("1", "2")), stringsAsFactors = TRUE)
  df_dummy <- dummyCodeFactorDf(df)
  expect_equal(dim(df_dummy), c(2, 4))
  expect_equal(colnames(df_dummy), c("v1_x", "v1_y", "v2_1", "v2_2"))

  # Non-factor dataframe error
  df_invalid <- data.frame(v1 = c(1, 2), v2 = factor(c("a", "b")), stringsAsFactors = TRUE)
  expect_error(dummyCodeFactorDf(df_invalid), "must have only factor variables")
})

test_that("squaredEuc, distFromData2Centroid, and withinClusterDist work", {
  expect_equal(squaredEuc(c(1, 2), c(4, 6)), 25)

  dat <- data.frame(x = c(0, 3), y = c(0, 4), stringsAsFactors = TRUE)
  centroid <- data.frame(x = 0, y = 0, stringsAsFactors = TRUE)
  dists <- distFromData2Centroid(dat, centroid, squaredEuc)
  expect_equal(dists, c(0, 25))

  centroids <- data.frame(x = c(0, 3), y = c(0, 4), stringsAsFactors = TRUE)
  memberships <- c(1, 2)
  wdist <- withinClusterDist(dat, centroids, squaredEuc, memberships)
  expect_equal(wdist, 0)

  # Error handling in withinClusterDist
  expect_error(withinClusterDist(dat, NULL, squaredEuc, memberships), "centroids cannot be NULL")
  expect_error(withinClusterDist(dat, data.frame(), squaredEuc, memberships), "Must include at least one centroid")
})

test_that("wkmeans works with factor and numeric data and enforces input validation", {
  set.seed(42)
  con <- data.frame(x = rnorm(20), y = rnorm(20), stringsAsFactors = TRUE)
  cat_fac <- data.frame(f = factor(rep(c("A", "B"), 10)), stringsAsFactors = TRUE)
  cat_num <- dummyCodeFactorDf(cat_fac)

  # With factor catData
  res1 <- wkmeans(conData = con, catData = cat_fac, conWeight = 0.6, nclust = 2)
  expect_s3_class(res1, "kmeans")
  expect_true(!is.null(res1$conCenters))
  expect_true(!is.null(res1$catCenters))

  # With numeric catData
  res2 <- wkmeans(conData = con, catData = cat_num, conWeight = 0.5, nclust = 2)
  expect_s3_class(res2, "kmeans")

  # Error validations
  expect_error(wkmeans(con, data.frame(f = c("A", "B")), 0.5, 2), "must be a data frame with all factor variables or all numeric")
  expect_error(wkmeans(data.frame(x = c("a", "b")), cat_fac, 0.5, 2), "must be data frame with all integer/numeric types")
  expect_error(wkmeans(con, cat_fac, 1.5, 2), "conWeight must be numeric and in \\[0,1\\]")
  expect_error(wkmeans(con, cat_fac, 0.5, 0), "nclust must be a positive integer")
})

test_that("clust2class, myPurity, and macroPrecRec work correctly", {
  clust <- c(1, 1, 2, 2)
  trueClass <- factor(c("A", "A", "B", "B"))
  pred <- clust2class(clust, trueClass)
  expect_equal(length(pred), 4)

  purity <- myPurity(clust, trueClass)
  expect_equal(purity, 1.0)

  mpr <- macroPrecRec(clust, trueClass)
  expect_equal(mpr$macroP, 1.0)
  expect_equal(mpr$macroR, 1.0)

  # Length mismatch errors
  expect_error(clust2class(c(1, 2), trueClass), "must be same length")

  # Level count mismatch in macroPrecRec
  expect_error(
    macroPrecRec(c(1, 1, 1, 1), factor(c("A", "B", "C", "D"))),
    "different number of levels"
  )
})
