library(kamila)
context('Modha-Spangler clustering')

test_that('gmsClust runs as expected, small data set', {
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(1)
  res1 <- gmsClust(
    conData = data.frame(rnorm(15),rnorm(15)),
    catData = data.frame(sample(0:1,size=15,rep=T),
			 sample(0:1,size=15,rep=T)),
    nclust = 2,
    searchDensity = 3
  )
  expect_equal(res1$results$cluster, c(1,1,2,1,1,2,1,2,2,1,2,1,2,2,2))
})
