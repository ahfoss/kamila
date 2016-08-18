library(kamila)
context('MEDEA helper functions')

test_that('setArgs selects from list appropriately', {
  expect_equal(setArgs(list(a=1,b=2,c=3),name='b',default=4), 2)
  expect_equal(setArgs(list(a=1,b=2,c=3),name='c',default=4), 3)
})

test_that('setArgs defaults correctly when item is not in list',{
  expect_equal(setArgs(list(a=1,b=2,c=3),name='d',default=4),4)
})

test_that('clustOneVar partitions variables as expected', {
  set.seed(1)
  expect_equivalent(clustOneVar(1:100,k=2, nstart=10, iter.max=20), rep(1:2,each=50))
  set.seed(1)
  expect_equivalent(
    clustOneVar(sort(rnorm(1000)),k=10, nstart=10, iter.max=20),
    rep(c(10,5,1,3,7,8,2,6,9,4),times=c(38,67,97,103,158,149,163,117,76,32)))
})
