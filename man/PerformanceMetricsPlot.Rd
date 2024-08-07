% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PerformanceMetrics.R
\name{PerformanceMetricsPlot}
\alias{PerformanceMetricsPlot}
\title{Performance Metrics Plot}
\usage{
PerformanceMetricsPlot(Y.pred, Y.test, family = "binomial")
}
\arguments{
\item{Y.pred}{A vector of predicted values. For SELPCCA, this is a vector of
predicted probabilities. For SIDA and SIDANet, this is a vector of predicted class.}

\item{Y.test}{A vector of observed test values.}

\item{family}{A string to denote the family for which  metrics should be provided.
Options are "gaussian", "binomial".}
}
\value{
A ggplot object
}
\description{
Creates a ggplot showing either an AUC curve (family == "binomial") or a scatter of predicted vs. observed values (family == "gaussian")
}
\details{
For a binary outcome, plots an AUC curve. For a continuous (gaussian) outcome, plots a scatter of observed vs. predicted values.
}
\examples{
\dontrun{
data(sidaData)
Xdata=sidaData[[1]]
Y=sidaData[[2]]
Xtestdata=sidaData[[3]]
Ytest=sidaData[[4]]
##---- call cross validation
 mycv=cvSIDA(Xdata,Y,withCov=FALSE,plotIt=FALSE, Xtestdata=Xtestdata,Ytest=Ytest,
             isParallel=FALSE,ncores=NULL,gridMethod='RandomSearch',
             AssignClassMethod='Joint',nfolds=5,ngrid=8,standardize=TRUE,
            maxiteration=20, weight=0.5,thresh=1e-03)
#check output
 test.error=mycv$sidaerror
 test.correlation=mycv$sidacorrelation
 optTau=mycv$optTau
 hatalpha=mycv$hatalpha
 #train metrics
 Y.pred=mycv$PredictedClass.train-1 #to get this in 0 and 1
 Y.train=Y-1 #to get this in 0 and 1
 train.metrics=PerformanceMetrics(Y.pred,Y.train,family='binomial')

 print(train.metrics)
 #obtain predicted class
 Y.pred=mycv$PredictedClass-1 #to get this in 0 and 1
 Ytest.in=Ytest-1 #to get this in 0 and 1
 PerformanceMetricsPlot(Y.pred,Ytest.in,family='binomial')
 }
}
\seealso{
\code{\link{cvSIDA}} \code{\link{selpscca.pred}} \code{\link{predict.SELPCCA}}\code{\link{PerformanceMetrics}}
}
