% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sidaPlots.R
\name{CorrelationPlots}
\alias{CorrelationPlots}
\title{Correlation Plots}
\usage{
CorrelationPlots(
  Xtestdata = Xtestdata,
  Ytest = Ytest,
  hatalpha = hatalpha,
  color.palette = NULL
)
}
\arguments{
\item{Xtestdata}{A list with each entry containing views of size
\eqn{ntest \times p_d}, where \eqn{d = 1, \dots, D}.Rows are samples
 and columns are variables. Can use testing or training data}

\item{Ytest}{\eqn{ntest \times 1} vector of class membership.}

\item{hatalpha}{A list of estimated sparse discriminant vectors for each view.}

\item{color.palette}{character vector of length K (number of classes), specifying the colors to use for the classes, respectively.
Defaults to shades of blue and orange (color.BlueOrange). Other option includes red and green combinations (color.GreenRed)}
}
\value{
\item{NULL}{}
}
\description{
Plots for visualizing correlation between estimated
discriminant vectors for pairwise data.
}
\details{
The function will return correlation plot(s).
}
\examples{
 #call sida
data(sidaData)
##---- call sida algorithm to estimate discriminant vectors, and predict on testing data

Xdata=sidaData[[1]]
Y=sidaData[[2]]
Xtestdata=sidaData[[3]]
Ytest=sidaData[[4]]

#call sidatunerange to get range of tuning parameter
ngrid=10
mytunerange=sidatunerange(Xdata,Y,ngrid,standardize=TRUE,weight=0.5,withCov=FALSE)

# an example with Tau set as the lower bound
Tau=c(mytunerange$Tauvec[[1]][1], mytunerange$Tauvec[[2]][1])
mysida=sida(Xdata,Y,Tau,withCov=FALSE,Xtestdata=Xtestdata,Ytest=Ytest,AssignClassMethod='Joint',
            plotIt=FALSE, standardize=TRUE,maxiteration=20,weight=0.5,thresh= 1e-03)

test.error=mysida$sidaerror
test.correlation=mysida$sidacorrelation

#estimated discriminant vectors and predicted class
hatalpha=mysida$hatalpha

predictedClass=mysida$PredictedClass

##----plot discriminant and correlation plots
#---------Correlation plot
mycorrplot=CorrelationPlots(Xtestdata,Ytest,mysida$hatalpha)

}
\references{
Sandra E. Safo, Eun Jeong Min, and Lillian Haine (2022), Sparse Linear Discriminant
Analysis for Multi-view Structured Data, Biometrics.
}
\seealso{
\code{\link{cvSIDA}} \code{\link{sidatunerange}}
\code{\link{DiscriminantPlots}}
}