
library(kamila)
context('Test medeaWgts')

test_that('medeaWgts partitions variables as expected',{
  
  set.seed(1)
  data1 <-  data.frame((1:100)/50,factor(sample(1:4,size=100,replace=T)))
  
  set.seed(1)
  med_output <- suppressWarnings(medeaWgts(
    data1,
    discretizationMethod='kmeans',
    associationMethod='ari',
    verbosity=T
  ))$verbose$cluster
  
  set.seed(1)
  ref_output <- suppressWarnings(kmeans(data1[,1],centers=4,nstart=10,iter.max=200))$cluster
  
  expect_equivalent(med_output, ref_output)
})

# add test from actual data set