test_that("calcApproxBIC calculates BIC and AIC criteria correctly", {
  resObj <- list(
    finalMemb = c(1, 1, 2, 2, 1, 2),
    finalCenters = matrix(c(0, 0, 1, 1), nrow = 2, ncol = 2),
    finalProbs = list(
      matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, ncol = 2),
      matrix(c(0.9, 0.1, 0.4, 0.6), nrow = 2, ncol = 2)
    ),
    finalLogLik = -15.5
  )

  # BIC = TRUE
  bicRes <- calcApproxBIC(resObj, BIC = TRUE)
  expect_named(bicRes, c("criteria", "k", "logLik", "nConParm", "nCatParm"))
  expect_equal(bicRes$nConParm, 4)
  expect_equal(bicRes$nCatParm, 4)
  expect_equal(bicRes$k, 9)
  expect_equal(bicRes$logLik, -15.5)
  expect_equal(bicRes$criteria, log(6) * 9 - 2 * (-15.5))

  # BIC = FALSE (AIC)
  aicRes <- calcApproxBIC(resObj, BIC = FALSE)
  expect_equal(aicRes$criteria, 2 * 9 - 2 * (-15.5))
})
