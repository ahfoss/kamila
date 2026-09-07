test_that("calcPsCpp computes correct prediction strength proportions and handles edge cases", {
  # Case 1: Standard inputs
  testMemb <- c(1, 1, 1, 2, 2, 2)
  teIntoTr <- c(1, 1, 1, 2, 2, 1)

  # Cluster 1: 3 items (1,1,1) -> all pairs (3 pairs) agree -> 3/3 = 1.0
  # Cluster 2: 3 items (2,2,1) -> pairs: (2,2) agree, (2,1) disagree, (2,1) disagree -> 1/3 = 0.3333333
  res <- calcPsCpp(testMemb, teIntoTr, 2)
  expect_equal(res[1], 1.0)
  expect_equal(res[2], 1 / 3)

  # Case 2: Cluster size < 2 returns NA
  testMemb2 <- c(1, 2, 2)
  teIntoTr2 <- c(1, 1, 2)
  res2 <- calcPsCpp(testMemb2, teIntoTr2, 2)
  expect_true(is.na(res2[1]))
  expect_equal(res2[2], 0.0)

  # Case 3: Empty inputs / zero clusters
  res3 <- calcPsCpp(integer(0), integer(0), 2)
  expect_equal(length(res3), 2)
  expect_true(all(is.na(res3)))

  # Case 4: Invalid numClust
  res4 <- calcPsCpp(testMemb, teIntoTr, 0)
  expect_equal(length(res4), 0)

  # Case 5: Compare against R table calculation on random data
  set.seed(42)
  N <- 300
  K <- 4
  tm <- sample(1:K, N, replace = TRUE)
  tr <- sample(1:K, N, replace = TRUE)

  cpp_res <- calcPsCpp(tm, tr, K)

  # R table reference calculation
  tab <- table(factor(tm, levels = 1:K), tr)
  r_res <- ifelse(
    rowSums(tab) < 2,
    NA_real_,
    rowSums(tab * (tab - 1) / 2) / (rowSums(tab) * (rowSums(tab) - 1) / 2)
  )

  expect_equal(unname(cpp_res), unname(r_res))
})

test_that("calcPsCpp performs identically to the legacy O(N^2) pairwise distance matrix implementation", {
  # Exact replica of the legacy pairwise algorithm previously in kamila.R
  legacy_pairwise_ps <- function(testMemb, teIntoTr, numClust) {
    numInTest <- length(testMemb)
    if (numInTest == 0 || numClust <= 0) {
      return(rep(NA_real_, max(0, numClust)))
    }
    testIndList <- mapply(
      x = 1:numClust,
      function(x) which(testMemb == x),
      SIMPLIFY = FALSE
    )
    dMat <- matrix(NaN, nrow = numInTest, ncol = numInTest)
    if (numInTest > 1) {
      for (i in 1:(numInTest - 1)) {
        for (j in (i + 1):numInTest) {
          dMat[i, j] <- teIntoTr[i] == teIntoTr[j]
        }
      }
    }
    psProps <- rep(0, numClust)
    for (cl in 1:numClust) {
      clustN <- length(testIndList[[cl]])
      if (clustN > 1) {
        for (i in 1:(clustN - 1)) {
          for (j in (i + 1):clustN) {
            psProps[cl] <- psProps[cl] + dMat[testIndList[[cl]][i], testIndList[[cl]][j]]
          }
        }
      }
      if (clustN < 2) {
        psProps[cl] <- NA_real_
      } else {
        psProps[cl] <- psProps[cl] / (clustN * (clustN - 1)) * 2
      }
    }
    return(psProps)
  }

  # Test across diverse parameter configurations
  test_configs <- list(
    list(N = 10, K = 2),
    list(N = 30, K = 3),
    list(N = 100, K = 5),
    list(N = 250, K = 4),
    list(N = 400, K = 6)
  )

  set.seed(2026)
  for (cfg in test_configs) {
    N <- cfg$N
    K <- cfg$K
    for (rep in 1:3) {
      testMemb <- sample(1:K, N, replace = TRUE)
      teIntoTr <- sample(1:K, N, replace = TRUE)

      res_cpp <- calcPsCpp(testMemb, teIntoTr, K)
      res_legacy <- legacy_pairwise_ps(testMemb, teIntoTr, K)

      expect_identical(length(res_cpp), length(res_legacy))
      expect_equal(res_cpp, res_legacy, tolerance = 1e-12)
    }
  }

  # Edge case: empty clusters and singleton clusters
  testMemb_edge <- c(1, 2, 2, 4, 4, 4, 4) # cluster 3 is empty, cluster 1 is singleton
  teIntoTr_edge <- c(1, 2, 1, 3, 3, 4, 3)
  res_cpp_edge <- calcPsCpp(testMemb_edge, teIntoTr_edge, 4)
  res_legacy_edge <- legacy_pairwise_ps(testMemb_edge, teIntoTr_edge, 4)
  expect_equal(res_cpp_edge, res_legacy_edge, tolerance = 1e-12)

  # Edge case: all observations in a cluster agree completely (PS = 1.0)
  testMemb_agree <- c(1, 1, 1, 1, 2, 2)
  teIntoTr_agree <- c(2, 2, 2, 2, 1, 1)
  res_cpp_agree <- calcPsCpp(testMemb_agree, teIntoTr_agree, 2)
  res_legacy_agree <- legacy_pairwise_ps(testMemb_agree, teIntoTr_agree, 2)
  expect_equal(res_cpp_agree, c(1.0, 1.0))
  expect_equal(res_cpp_agree, res_legacy_agree, tolerance = 1e-12)

  # Edge case: all observations in a cluster disagree completely (PS = 0.0)
  testMemb_disagree <- c(1, 1, 1)
  teIntoTr_disagree <- c(1, 2, 3)
  res_cpp_disagree <- calcPsCpp(testMemb_disagree, teIntoTr_disagree, 1)
  res_legacy_disagree <- legacy_pairwise_ps(testMemb_disagree, teIntoTr_disagree, 1)
  expect_equal(res_cpp_disagree, 0.0)
  expect_equal(res_cpp_disagree, res_legacy_disagree, tolerance = 1e-12)
})
