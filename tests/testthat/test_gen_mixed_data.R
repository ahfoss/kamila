test_that("genMixedData works with standard inputs", {
  set.seed(123)
  dat <- genMixedData(
    sampSize = 60,
    nConVar = 2,
    nCatVar = 2,
    nCatLevels = 4,
    nConWithErr = 1,
    nCatWithErr = 1,
    popProportions = c(0.4, 0.6),
    conErrLev = 0.2,
    catErrLev = 0.2
  )

  expect_named(
    dat,
    c("trueID", "trueMus", "conVars", "errVariance", "popProbsNoErr", "popProbsWithErr", "catVars")
  )
  expect_equal(length(dat$trueID), 60)
  expect_equal(dim(dat$conVars), c(60, 2))
  expect_equal(dim(dat$catVars), c(60, 2))
  expect_equal(sort(unique(dat$trueID)), c(1, 2))
})

test_that("genMixedData works with logical nConWithErr vector", {
  set.seed(123)
  dat <- genMixedData(
    sampSize = 40,
    nConVar = 3,
    nCatVar = 2,
    nCatLevels = 4,
    nConWithErr = c(TRUE, FALSE, TRUE),
    nCatWithErr = 1,
    popProportions = c(0.5, 0.5),
    conErrLev = 0.1,
    catErrLev = 0.1
  )
  expect_equal(dim(dat$conVars), c(40, 3))
})

test_that("genMixedData handles edge cases for error variables", {
  set.seed(123)
  # Zero con error variables with nCatWithErr = 1
  dat0 <- genMixedData(
    sampSize = 30,
    nConVar = 2,
    nCatVar = 2,
    nCatLevels = 4,
    nConWithErr = 0,
    nCatWithErr = 1,
    popProportions = c(0.5, 0.5),
    conErrLev = 0.1,
    catErrLev = 0.1
  )
  expect_equal(dim(dat0$catVars), c(30, 2))

  # All cat error variables
  datAll <- genMixedData(
    sampSize = 30,
    nConVar = 2,
    nCatVar = 2,
    nCatLevels = 4,
    nConWithErr = 2,
    nCatWithErr = 2,
    popProportions = c(0.5, 0.5),
    conErrLev = 0.1,
    catErrLev = 0.1
  )
  expect_equal(dim(datAll$catVars), c(30, 2))

  # nConVar = 0
  datNoCon <- genMixedData(
    sampSize = 30,
    nConVar = 0,
    nCatVar = 2,
    nCatLevels = 4,
    nConWithErr = 0,
    nCatWithErr = 1,
    popProportions = c(0.5, 0.5),
    conErrLev = 0.1,
    catErrLev = 0.1
  )
  expect_null(datNoCon$conVars)
  expect_null(datNoCon$trueMus)
  expect_null(datNoCon$errVariance)

  # nCatWithErr = 0
  datNoCatErr <- genMixedData(
    sampSize = 30,
    nConVar = 2,
    nCatVar = 2,
    nCatLevels = 4,
    nConWithErr = 1,
    nCatWithErr = 0,
    popProportions = c(0.5, 0.5),
    conErrLev = 0.1,
    catErrLev = 0.1
  )
  expect_equal(dim(datNoCatErr$catVars), c(30, 2))
})

test_that("genMixedData throws errors on invalid input", {
  # More than 2 populations
  expect_error(
    genMixedData(30, 2, 2, 4, 1, 1, c(0.3, 0.3, 0.4), 0.1, 0.1),
    "More than two populations"
  )

  # Logical nConWithErr length mismatch
  expect_error(
    genMixedData(30, 2, 2, 4, c(TRUE), 1, c(0.5, 0.5), 0.1, 0.1),
    "length of nConWithErr must equal nConVar"
  )

  # nConWithErr > nConVar
  expect_error(
    genMixedData(30, 2, 2, 4, 5, 1, c(0.5, 0.5), 0.1, 0.1),
    "nConWithErr must be less than nConVar"
  )

  # nCatLevels not a multiple of popProportions
  expect_error(
    genMixedData(30, 2, 2, 5, 1, 1, c(0.5, 0.5), 0.1, 0.1),
    "number of categorical levels must"
  )

  # Improper nConWithErr type/length
  expect_error(
    genMixedData(30, 2, 2, 4, c(1, 2), 1, c(0.5, 0.5), 0.1, 0.1),
    "improper input to nConWithErr"
  )
})
