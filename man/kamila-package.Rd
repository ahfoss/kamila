\name{kamila-package}
\alias{kamila-package}
\docType{package}
\title{
Clustering for mixed continuous and categorical data sets
}
\description{
A collection of methods for clustering mixed type data, including KAMILA (KAy-means for MIxed LArge data) and a flexible implementation of Modha-Spangler clustering
}
\details{
\tabular{ll}{
Package: \tab kamila\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2015-10-06\cr
License: \tab GPL-3\cr
}

}
\author{
Alex Foss and Marianthi Markatou

Maintainer: Alex Foss <alexanderhfoss@gmail.com>
}
\references{
AH Foss, M Markatou, B Ray, and A Heching (in press). A semiparametric method for clustering mixed data. Machine Learning, DOI: 10.1007/s10994-016-5575-7.

DS Modha and S Spangler (2003). Feature weighting in k-means clustering. Machine Learning 52(3), 217-237.
}

\examples{
\dontrun{
# import and format a mixed-type data set
data(Byar, package='clustMD')
Byar$logSpap <- log(Byar$Serum.prostatic.acid.phosphatase)
conInd <- c(5,6,8:10,16)
conVars <- Byar[,conInd]
conVars <- data.frame(scale(conVars))

catVarsFac <- Byar[,-c(1:2,conInd,11,14,15)]
catVarsFac[] <- lapply(catVarsFac, factor)
catVarsDum <- dummyCodeFactorDf(catVarsFac)

# Modha-Spangler clustering with kmeans default Hartigan-Wong algorithm
gmsResHw <- gmsClust(conVars, catVarsDum, nclust = 3)

# Modha-Spangler clustering with kmeans Forgy-Lloyd algorithm
# NOTE searchDensity should be >= 10 for optimal performance:
# this is just a syntax demo
gmsResLloyd <- gmsClust(conVars, catVarsDum, nclust = 3,
  algorithm = "Lloyd", searchDensity = 3)

# KAMILA clustering
kamRes <- kamila(conVars, catVarsFac, numClust=3, numInit=10)

# Plot results
ternarySurvival <- factor(Byar$SurvStat)
levels(ternarySurvival) <- c('Alive','DeadProst','DeadOther')[c(1,2,rep(3,8))]
plottingData <- cbind(
  conVars,
  catVarsFac,
  KamilaCluster = factor(kamRes$finalMemb),
  MSCluster = factor(gmsResHw$results$cluster))
plottingData$Bone.metastases <- ifelse(
  plottingData$Bone.metastases == '1', yes='Yes',no='No')

# Plot Modha-Spangler/Hartigan-Wong results
msPlot <- ggplot(
  plottingData,
  aes(
    x=logSpap,
    y=Index.of.tumour.stage.and.histolic.grade,
    color=ternarySurvival,
    shape=MSCluster))
plotOpts <- function(pl) (pl + geom_point() + 
  scale_shape_manual(values=c(2,3,7)) + geom_jitter())
plotOpts(msPlot)

# Plot KAMILA results
kamPlot <- ggplot(
  plottingData,
  aes(
    x=logSpap,
    y=Index.of.tumour.stage.and.histolic.grade,
    color=ternarySurvival,
    shape=KamilaCluster))
plotOpts(kamPlot)
}
}
