library(kamila)

test_that("gmsClust runs as expected, small data set", {
  suppressWarnings(withr::local_rng_version("3.5.0"))
  set.seed(1)
  res1 <- gmsClust(
    conData = data.frame(rnorm(15), rnorm(15)),
    catData = data.frame(
      sample(0:1, size = 15, rep = TRUE),
      sample(0:1, size = 15, rep = TRUE)
    ),
    nclust = 2,
    searchDensity = 3
  )
  expect_equal(res1$results$cluster, c(1, 1, 2, 1, 1, 2, 1, 2, 2, 1, 2, 1, 2, 2, 2))
})

test_that("gmsClust handles zero objective function entries and Qcon/Qcat bounds", {
  conData <- data.frame(x = c(1, 2, 10, 11), y = c(1, 2, 10, 11))
  catData <- data.frame(v1 = c(1, 1, 1, 1), v2 = c(0, 0, 0, 0))

  expect_warning(
    gmsClust(conData, catData, nclust = 2, searchDensity = 3),
    "At least one entry of zero in the objective function"
  )

  conData2 <- data.frame(x = c(-100, 100, -100, 100))
  catData2 <- data.frame(v1 = c(1, 0, 1, 0), v2 = c(0, 1, 0, 1))
  res2 <- suppressWarnings(gmsClust(conData2, catData2, nclust = 2, searchDensity = 10))
  expect_true(length(res2$results$cluster) == 4)
})
