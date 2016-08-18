library(kamila)
context('KAMILA clustering')

test_that('KAMILA runs as expected, small data set', {
  set.seed(1)
  res1 <- kamila(
    conVar = data.frame(rnorm(15)),
    catFactor = data.frame(factor(sample(1:4,size=15,rep=T))),
    numClust = 2,
    numInit = 10,
    maxIter = 25,
    conInitMethod = 'runif',
    catBw = 0.025,
    verbose = FALSE,
    calcNumClust = 'none',
    numPredStrCvRun = 20,
    predStrThresh = 0.8
  )
  expect_identical(res1$finalMemb, c(2,1,2,1,1,2,1,1,1,2,1,1,2,2,1))
})

test_that('KAMILA warns if pred-strength is used with a single numClust',{
  expect_warning(
    kamila(
      conVar = data.frame(rnorm(10),rnorm(10))
      ,catFactor = data.frame(factor(sample(1:4,size=10,rep=T)),
                              factor(sample(1:4,size=10,rep=T)))
      ,numClust = 2
      ,numInit = 2
      ,maxIter = 2
      ,conInitMethod = 'runif'
      ,catBw = 0.025
      ,verbose = FALSE
      ,calcNumClust = 'ps'
      ,numPredStrCvRun = 2
      ,predStrThresh = 0.8
    ),
    'Input parameter numClust is a scalar; the prediction strength'
  )
})

test_that('KAMILA throws error if numClust is length > 1 with calcNumClust=="none"',{
  expect_error(
    kamila(
      conVar = data.frame(rnorm(10),rnorm(10))
      ,catFactor = data.frame(factor(sample(1:4,size=10,rep=T)),
                              factor(sample(1:4,size=10,rep=T)))
      ,numClust = 2:5
      ,numInit = 2
      ,maxIter = 2
      ,conInitMethod = 'runif'
      ,catBw = 0.025
      ,verbose = FALSE
      ,calcNumClust = 'none'
      ,numPredStrCvRun = 2
      ,predStrThresh = 0.8
    ),
    'Input parameter numClust must be length 1 if calcNumClust == "none"'
  )
})
