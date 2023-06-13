

#' @title Data example for SIDA
#' @description Simulated data to demonstrate the use of SIDA.
#' @format A list with 4 elements:
#' \describe{
#'   \item{\code{Xdata}}{A list with each entry containing two views
#'    of training data with dimension 160 × 2000 each. Rows are
#'    samples and columns are variables.}
#'   \item{\code{Y}}{160 × 1 vector of training class membership.
#'   There are two classes each with size 80.}
#'   \item{\code{Xtestdata}}{A list with each entry containing two views
#'    of testing data with dimension 320 × 2000 each. Rows are samples
#'     and columns are variables.}
#'   \item{\code{Ytest}}{ 320 × 1 vector of testing class membership.
#'   There are two classes each with size 160.}
#'}
#'
#' @references
#' Sandra E. Safo, Eun Jeong Min, and Lillian Haine (2019), Sparse Linear Discriminant
#' Analysis for Multi-view Structured Data, submitted.
"sidaData"




#' @title Data example for SIDANet
#' @description Simulated data to demonstrate the use of SIDANet.
#' @format A list with 6 elements:
#' \describe{
#'   \item{\code{XdataNet}}{A list with each entry containing two views
#'   of training data with dimension 240 × 1000 each. Rows are samples
#'   and columns are variables.}
#'   \item{\code{YNet}}{240 × 1 vector of training class membership.
#'   There are three classes each with size 80.}
#'   \item{\code{XtestdataNet}}{A list with each entry containing
#'   two views of testing data with dimension 480 × 1000 each. Rows
#'   are samples and columns are variables.}
#'   \item{\code{YtestNet}}{ 480 × 1 vector of testing class membership.
#'   There are three classes each with size 160.}
#'   \item{\code{myedges}}{A list with each entry containing a 36 × 2
#'    matrix of edge information for each view. Assumes variable
#'    1 is connected to variables 2 to 10, variable 11 is connected
#'    to variables 12 to 20, variable 21 is connected to variables 22
#'     to 30 and variable 31 is connected to variables 32 to 40.
#'     All remaining variables are singletons.}
#'   \item{\code{myedgeweight}}{A list with each entry containing
#'   edgeweight. In this example, views 1 and 2 have edge weights so the
#'   Laplacian of a weighted graph will be used.}
#'}
#'
#' @references
#' Sandra E. Safo, Eun Jeong Min, and Lillian Haine (2019), Sparse Linear Discriminant
#' Analysis for Multi-view Structured Data, submitted.
"sidanetData"


#' @title Data example for SELPscca
#' @description Simulated data with one true canonical correlation
#' vectors for first and second datasets. The first 20 and 15
#' variables are nonzero (i.e., signal variables) in the first
#' canonical correlation vectors for the first and second datasets
#' respectively.
#' @format A list with 7 elements:
#' \describe{
#'   \item{\code{Xdata1}}{A matrix of size \eqn{80 \times 200} for first dataset.
#'   Rows are samples and columns are variables.}
#'   \item{\code{Xdata2}}{A matrix of size \eqn{80 \times 150} for second dataset.
#'   Rows are samples and columns are variables.}
#'   \item{\code{Xtestdata1}}{A matrix of size \eqn{400 \times 200} for first dataset.
#'   Rows are samples and columns are variables.}
#'   \item{\code{Xtestdata2}}{A matrix of size \eqn{400 \times 150} for second dataset.
#'   Rows are samples and columns are variables.}
#'   \item{\code{TrueAlpha}}{The first canonical correlation vector for Xdata1.}
#'   \item{\code{TrueBeta}}{The first canonical correlation vector for Xdata2.}
#'   \item{\code{TrueCorr}}{The first canonical correlation coefficient.}
#'}
#'
#' @references
#' Sandra E. Safo, Jeongyoun Ahn, Yongho Jeon, and Sungkyu Jung (2018),
#' Sparse Generalized Eigenvalue Problem with Application to Canonical
#' Correlation Analysis for Integrative Analysis of Methylation and Gene
#' Expression Data. Biometrics
"selpData"


#' @title Multiomics data pertaining to COVID-19
#' @description RNA Sequencing (RNASeq) and Proteomics data pertaining to COVID-19. Clinical data
#' are also available. Please refer to Overmyer et.al (2021) for a description of the data and Lipman et.al (2022) for how data
#' were pre-processed.
#' @format A list with 3 elements:
#' \describe{
#'   \item{\code{COVIDData[[1]]}}{Proteomics data. A data frame of size \eqn{120 \times 264}.
#'   Rows are samples and columns are variables.}
#'   \item{\code{COVIDData[[2]]}}{RNASeq data. A data frame of size \eqn{120 \times 5800}.
#'   Rows are samples and columns are variables.}
#'   \item{\code{COVIDData[[2]]}}{Clinical and demographic data. A data frame of size \eqn{120 \times 18}.
#'   Rows are samples and columns are variables.}
#'}
#'
#' @references
#' Multi-omic analysis reveals enriched pathways
#' associated with COVID-19 and COVID-19 severity. PLOS ONE, 17(4)
#' Overmyer, K.A., Shishkova, E., Miller, I.J., Balnis, J., Bernstein, M.N., Peters-Clarke, T.M., Meyer, J.G., Quan,
#' Q., Muehlbauer, L.K., Trujillo, E.A., et al.: Large-scale multi-omic analysis of covid-19 severity. Cell systems
#' 12(1), 23–40 (2021)
"COVIDData"
