test_that("prediction.strength, classify.kmeans, cyclicalCoding, getSubsetCon, and getSubsetMix work", {
  # prediction.strength
  testClust <- c(1, 1, 2, 2)
  classifyTest <- c(1, 1, 2, 2)
  ps <- prediction.strength(testClust, classifyTest)
  expect_equal(ps, 1)

  # classify.kmeans
  kmRes <- list(centers = matrix(c(0, 0, 5, 5), nrow = 2, byrow = TRUE))
  newData <- matrix(c(0.1, 0.1, 4.9, 5.1), nrow = 2, byrow = TRUE)
  classRes <- classify.kmeans(kmRes, newData)
  expect_equal(classRes, c(1, 2))

  # cyclicalCoding
  cyc <- cyclicalCoding(c(1, 2, 3, 4, 5, 6, 7))
  expect_equal(dim(cyc), c(7, 2))

  # getSubsetCon
  mat <- matrix(1:10, nrow = 5)
  expect_equal(dim(getSubsetCon(mat, 1:3)), c(3, 2))

  # getSubsetMix
  mix2 <- list(
    data.frame(x = 1:5, stringsAsFactors = TRUE),
    data.frame(f = factor(1:5), stringsAsFactors = TRUE)
  )
  sub2 <- getSubsetMix(mix2, 1:2)
  expect_equal(length(sub2), 2)
  expect_equal(NROW(sub2[[1]]), 2)

  mix3 <- list(
    data.frame(x = 1:5, stringsAsFactors = TRUE),
    data.frame(f = factor(1:5), stringsAsFactors = TRUE),
    data.frame(c = factor(1:5), stringsAsFactors = TRUE)
  )
  sub3 <- getSubsetMix(mix3, 1:2)
  expect_equal(length(sub3), 3)

  # Invalid list length for getSubsetMix
  expect_error(getSubsetMix(list(1), 1:2), "must be list of length 2 or 3")
})

test_that("kamilaMethod works with list data", {
  set.seed(42)
  dat <- list(
    data.frame(x = rnorm(20), y = rnorm(20), stringsAsFactors = TRUE),
    data.frame(f = factor(rep(c("a", "b"), 10)), stringsAsFactors = TRUE)
  )
  res <- kamilaMethod(dat, k = 2)
  expect_true("cluster" %in% names(res))
  expect_equal(length(res$cluster), 20)
})

test_that("nclust1x and nclustFull run without error", {
  set.seed(123)
  conData <- matrix(rnorm(40), nrow = 20)
  res1x <- nclust1x(conData, kmax = 3, verbose = FALSE)
  expect_named(res1x, c("psVec", "psMatrix"))
  expect_equal(dim(res1x$psMatrix), c(2, 2))

  resFull <- nclustFull(conData, kmax = 3, nrep = 2)
  expect_named(resFull, c("k", "psVec", "scores", "thresh", "kmax", "nrep"))
  expect_true(resFull$k %in% c(2, 3))

  # verbose = TRUE in nclust1x
  expect_output(nclust1x(conData, kmax = 3, verbose = TRUE), "Now starting k =")

  # high psThresh fallback in nclustFull
  resFullHigh <- nclustFull(conData, kmax = 3, psThresh = 0.999, nrep = 2)
  expect_equal(resFullHigh$thresh, 0.999)
})
