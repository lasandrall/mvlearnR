% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GenericPlots.R
\name{VarImportancePlot}
\alias{VarImportancePlot}
\title{Variable Importance Plot}
\usage{
VarImportancePlot(object)
}
\arguments{
\item{object}{the output from SIDA, SIDANet, and SELPCCA methods}
}
\value{
A graph of the absolute loadings for variables selected. The variables
are normalized to the variable with the largest weight.
}
\description{
Wrapper function to visualize loadings for variables selected
by SIDA, SIDANet, and SELPCCA methods.
}
\examples{
##---- load SIDA data
data("sidaData")
Xdata <- sidaData[[1]]
Y <- sidaData[[2]]
Xtestdata <- sidaData[[3]]
Ytest <- sidaData[[4]]
##---- call cross validation
mycv=cvSIDA(Xdata,Y,withCov=FALSE,plotIt=FALSE, Xtestdata=Xtestdata,Ytest=Ytest)
##----  Obtain variable importance plot
VarImportancePlot(mycv)
}
