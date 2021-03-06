% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kamila.R
\name{classifyKamila}
\alias{classifyKamila}
\title{Classify new data into existing KAMILA clusters}
\usage{
classifyKamila(obj, newData)
}
\arguments{
\item{obj}{An output object from the kamila function.}

\item{newData}{A list of length 2, with first element a data frame of continuous variables, and second element a data frame of categorical factors.}
}
\value{
An integer vector denoting cluster assignments of the new data points.
}
\description{
A function that classifies a new data set into existing KAMILA clusters
using the output object from the kamila function.
}
\details{
A function that takes obj, the output from the kamila function, and newData,
a list of length 2, where the first element is a data frame of continuous
variables, and the second element is a data frame of categorical factors.
Both data frames must have the same format as the original data used
to construct the kamila clustering.
}
\examples{
# Generate toy data set
set.seed(1234)
dat1 <- genMixedData(400, nConVar = 2, nCatVar = 2, nCatLevels = 4,
  nConWithErr = 2, nCatWithErr = 2, popProportions = c(.5,.5),
  conErrLev = 0.2, catErrLev = 0.2)
# Partition the data into training/test set
trainingIds <- sample(nrow(dat1$conVars), size = 300, replace = FALSE)
catTrain <- data.frame(apply(dat1$catVars[trainingIds,], 2, factor), stringsAsFactors = TRUE)
conTrain <- data.frame(scale(dat1$conVars)[trainingIds,], stringsAsFactors = TRUE)
catTest <- data.frame(apply(dat1$catVars[-trainingIds,], 2, factor), stringsAsFactors = TRUE)
conTest <- data.frame(scale(dat1$conVars)[-trainingIds,], stringsAsFactors = TRUE)
# Run the kamila clustering procedure on the training set
kamilaObj <- kamila(conTrain, catTrain, numClust = 2, numInit = 10)
table(dat1$trueID[trainingIds], kamilaObj$finalMemb)
# Predict membership in the test data set
kamilaPred <- classifyKamila(kamilaObj, list(conTest, catTest))
table(dat1$trueID[-trainingIds], kamilaPred)
}
\references{
Foss A, Markatou M; kamila: Clustering Mixed-Type Data in R and Hadoop. Journal of Statistical Software, 83(13). 2018. doi: 10.18637/jss.v083.i13
}
