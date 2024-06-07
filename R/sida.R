

#' @title Sparse Integrative Discriminant Analysis for Multi-View Data
#'
#' @description Performs sparse integrative discriminant analysis of
#' multi-view data to 1) obtain discriminant vectors that are associated
#' and optimally separate subjects into different classes 2) estimate
#' misclassification rate, and total correlation coefficient. Allows for
#' the inclusion of other covariates which are not penalized in the
#' algorithm. It is recommended to use cvSIDA to choose best tuning parameter.
#'
#' @param Xdata A list with each entry containing training views of size \eqn{n \times p_d},
#' where \eqn{d =1, \dots , D} views. Rows are samples and columns are variables. If
#' covariates are available, they should be included as a separate view, and
#' set as the last dataset. For binary or categorical covariates (assumes no
#' ordering), we suggest the use of indicator variables.
#' @param Y \eqn{n \times 1} vector of class membership. Numeric, coded as 1, 2, ....
#' @param Tau \eqn{d \times 1} vector of tuning parameter. It is recommended to
#' use sidatunerange to obtain lower and upper bounds for the tuning parameters
#' since too large a tuning parameter will result in a trivial solution vector
#' (all zeros) and too small may result in non-sparse vectors.
#' @param withCov TRUE or FALSE if covariates are available. If TRUE, please set
#' all covariates as one dataset and should be the last dataset. For binary and
#' categorical variables, use indicator matrices/vectors. Default is FALSE.
#' @param Xtestdata A list with each entry containing testing views of size
#' \eqn{ntest \times p_d}, where \eqn{d = 1, \dots, D}. Rows are samples and columns are variables.
#' The order of the list should be the same as the order for the training data,
#' Xdata. Use if you want to predict on a testing dataset. If no Xtestdata, set to NULL.
#' @param Ytest \eqn{ntest \times 1} vector of test class membership. If no testing data
#' provided, set to NULL.
#' @param AssignClassMethod Classification method. Either Joint or Separate. Joint uses
#' all discriminant vectors from D datasets to predict class membership. Separate predicts
#' class membership separately for each dataset. Default is Joint.
#' @param plotIt TRUE or FALSE. If TRUE, produces discriminants and correlation plots.
#' Default is FALSE.
#' @param standardize TRUE or FALSE. If TRUE, data will be normalized to have mean zero
#' and variance one for each variable. Default is TRUE.
#' @param maxiteration Maximum iteration for the algorithm if not converged.Default is 20.
#' @param weight Balances separation and association. Default is 0.5.
#' @param thresh Threshold for convergence. Default is 0.001.
#'
#' @details The function will return several R objects, which can be assigned to a variable.
#' To see the results, use the “$" operator.
#'
#' @return The output is a list containing the following components.
#'  \item{sidaerror}{Estimated classication error. If testing data provided, this will
#'  be test classification error, otherwise, training error}
#'  \item{sidacorrelation}{Sum of pairwise RV coefficients. Normalized to be within 0 and 1, inclusive.}
#'  \item{hatalpha}{A list of estimated sparse discriminant vectors for each view.}
#'  \item{PredictedClass}{ Predicted class. If AssignClassMethod=’Separate’, this will be a
#'  \eqn{ntest \times D} matrix, with each column the predicted class for each data.}
#'
#'
#' @seealso \code{\link{cvSIDA}}  \code{\link{sidatunerange}}
#'
#' @references
#' Sandra E. Safo, Eun Jeong Min, and Lillian Haine (2022), Sparse Linear Discriminant
#' Analysis for Multi-view Structured Data, Biometrics.
#'
#' @import RSpectra
#' @import foreach
#' @import doParallel
#' @import parallel
#' @importFrom graphics legend par plot points
#' @importFrom stats aggregate density quantile
#' @importFrom utils combn
#' @importFrom graphics lines
#' @importFrom Matrix Matrix
#' @importFrom CVXR norm diag
#'
#' @export
#' @examples
#' \dontrun{
#'  #call sida
#' data(sidaData)
#' ##---- call sida algorithm to estimate discriminant vectors, and predict on testing data
#'
#' Xdata=sidaData[[1]]
#' Y=sidaData[[2]]
#' Xtestdata=sidaData[[3]]
#' Ytest=sidaData[[4]]
#'
#' #call sidatunerange to get range of tuning parameter
#' ngrid=10
#' mytunerange=sidatunerange(Xdata,Y,ngrid,standardize=TRUE,weight=0.5,withCov=FALSE)
#'
#' # an example with Tau set as the lower bound
#' Tau=c(mytunerange$Tauvec[[1]][1], mytunerange$Tauvec[[2]][1])
#' mysida=sida(Xdata,Y,Tau,withCov=FALSE,Xtestdata=Xtestdata,Ytest=Ytest,AssignClassMethod='Joint',
#'             plotIt=FALSE, standardize=TRUE,maxiteration=20,weight=0.5,thresh= 1e-03)
#'
#' test.error=mysida$sidaerror
#' test.correlation=mysida$sidacorrelation
#'
#' #estimated discriminant vectors and predicted class
#' hatalpha=mysida$hatalpha
#'
#' predictedClass=mysida$PredictedClass
#' }
sida=function(Xdata=Xdata,Y=Y,Tau=Tau,withCov=FALSE,
              Xtestdata=NULL,Ytest=NULL,
              AssignClassMethod='Joint',plotIt=FALSE,
              standardize=TRUE,maxiteration=20,weight=0.5,
              thresh= 1e-03){

  XdataOrig=Xdata
  XtestdataOrig=Xtestdata
  YOrig=Y
  YtestOrig=Ytest
  #check inputs
  dsizes=lapply(Xdata, function(x) dim(x))
  n=dsizes[[1]][1]
  nsizes=lapply(Xdata, function(x) dim(x)[1])

  if(all(nsizes!=nsizes[[1]])){
    stop('The datasets  have different number of observations')
  }

  #check data
  if (is.list(Xdata)) {
    D = length(Xdata)
  } else {
    stop("Input data should be a list")
  }

  # #check if minimum of Y is 0, if so shift everything by 1
  # if(min(YOrig)==1){
  #   Y=YOrig+ 1
  #   if(!is.null(Ytest)){
  #     Ytest=YtestOrig+1
  #   }
  # }

  #check if testing data are provided. If not, will set training data as testing data.
  if(is.null(Xtestdata)){
    Xtestdata=Xdata
    Ytest=Y
  }

  if(is.null(AssignClassMethod)){
    AssignClassMethod='Joint'
  }

  if(is.null(withCov)){
    withCov=FALSE
  }

  if(is.null(plotIt)){
    plotIt=FALSE
  }

  if(is.null(standardize)){
    standardize=TRUE
  }

  if(is.null(maxiteration)){
    maxiteration=20
  }

  if(is.null(weight)){
    weight=0.5
  }

  if(is.null(thresh)){
    thresh=1e-03
  }


  #standardize if true
  if(standardize==TRUE){
    Xdata=lapply(Xdata,function(x)scale(x,center=TRUE,scale=TRUE))
    Xtestdata=lapply(Xtestdata,function(x)scale(x,center=TRUE,scale=TRUE))
  }


  nK=length(unique(as.vector(Y))) -1


  #norm function for convergence
  normdiff=function(xnew,xold){
    ndiff=CVXR::norm(xnew-xold,'f')^2 / CVXR::norm(xold,'f')^2
  }
  #initialize
  iter=0
  diffalpha=1
  reldiff=1

  mynsparse=myfastIDAnonsparse(Xdata,Y,weight)
  myalpha=mynsparse$myalphaoldmat
  #while convergence is not met
  while(iter < maxiteration && min(reldiff,max(diffalpha))> thresh){
    iter=iter+1
    #cat("current iteration is", iter, "\n")

    myalphaold=myalpha
    mysidainner=sidainner(Xdata,Y,mynsparse$sqrtminvmat,myalphaold,mynsparse$tildealphamat, mynsparse$tildelambda,Tau,weight,withCov)

    myalpha=mysidainner$hatalpha
    nz=sapply(1:D, function(i) list(colSums(myalpha[[i]]!=0)))
    nz=cbind(c(do.call(rbind,nz)))
    if(any(nz==0)){
      myalpha=myalphaold
      break
    }
    diffalpha=mapply(normdiff, myalpha, myalphaold)
    sumnormdiff=sum(sapply(1:D, function(i) CVXR::norm(myalpha[[i]]-myalphaold[[i]],'f')^2 ))
    sumnormold=sum(sapply(1:D, function(i) CVXR::norm(myalphaold[[i]],'f')^2    ))
    reldiff=sumnormdiff/sumnormold
  }

  #classification
  if(AssignClassMethod=='Joint'){
    myclassify=sidaclassify(myalpha,Xtestdata,Xdata,Y,AssignClassMethod='Joint')
    sidaerror=sum(myclassify$PredictedClass!=Ytest)/length(Ytest)
  }else if(AssignClassMethod=='Separate'){
    myclassify=sidaclassify(myalpha,Xtestdata,Xdata,Y,AssignClassMethod='Separate')

    sidaerror=sapply(1:(nK+1), function(x)  sum(myclassify$PredictedClass[,x]!=Ytest)/length(Ytest) )
  }


  #sum pairwise correlations of training data
  ss=list()
  #sum pairwise correlations
  for(d in 1:D){
    dd=setdiff(seq(1, D, by= 1),d)
    #correlations
    sumCorr2=0;
    for (jj in 1:length(dd)){
      j=dd[jj];
      X1=Xtestdata[[d]]%*%myalpha[[d]]
      X2=Xtestdata[[j]]%*%myalpha[[j]]
      X1=scale(X1, center=TRUE,scale=FALSE)
      X2=scale(X2, center=TRUE,scale=FALSE)
      X1X2=t(X1)%*%X2/dim(X1)[1]
      X1X1=t(X1)%*%X1/dim(X1)[1]
      X2X2=t(X2)%*%X2/dim(X2)[1]
      sumcorr3=sum(CVXR::diag(X1X2%*%t(X1X2)))/(sqrt(sum(CVXR::diag(X1X1%*%X1X1)))*sqrt(sum(CVXR::diag(X2X2%*%X2X2))))
      sumCorr2=sumCorr2+sumcorr3
    }
    ss[[d]]=sumCorr2/length(dd)
  }

  sidacorrelation=sum(do.call(rbind,ss))/D

  #Produce discriminant and correlation plot if plotIt=T
  if(plotIt==TRUE){
    DiscriminantPlots(Xtestdata,Ytest,myalpha)
    CorrelationPlots(Xtestdata,Ytest,myalpha)
  }else{
    myDiscPlot=NULL
    myCorrPlot=NULL
  }
  result=list(sidaerror=sidaerror,sidacorrelation=sidacorrelation,hatalpha=myalpha,
              PredictedClass=myclassify$PredictedClass)
  class(result)="SIDA"
  return(result)

}



#' @title Cross validation for Sparse Integrative Discriminant Analysis for Multi-View Data
#'
#' @description Performs nfolds cross validation to select optimal tuning
#' parameters for sida based on training data, which are then used with
#' the training or testing data to predict class membership. Allows for
#' inclusion of covariates which are not penalized. If you want to apply
#' optimal tuning parameters to testing data, you may also use sida.
#'
#'
#' @param Xdata A list with each entry containing training views of size \eqn{n \times p_d},
#' where \eqn{d =1, \dots , D} views. Rows are samples and columns are variables. If
#' covariates are available, they should be included as a separate view, and
#' set as the last dataset. For binary or categorical covariates (assumes no
#' ordering), we suggest the use of indicator variables.
#' @param Y \eqn{n \times 1} vector of class membership. Numeric, coded as 1, 2, ....
#' @param withCov TRUE or FALSE if covariates are available. If TRUE, please set
#' all covariates as one dataset and should be the last dataset. For binary and
#' categorical variables, use indicator matrices/vectors. Default is FALSE.
#' @param plotIt TRUE or FALSE. If TRUE, produces discriminants and correlation
#'  plots. Default is FALSE.
#' @param Xtestdata A list with each entry containing testing views of size
#' \eqn{ntest \times p_d}, where \eqn{d = 1, \dots, D}. Rows are samples and columns are variables.
#' The order of the list should be the same as the order for the training data,
#' Xdata. Use if you want to predict on a testing dataset. If no Xtestdata, set to NULL.
#' @param Ytest \eqn{ntest \times 1} vector of test class membership. If no testing data
#' provided, set to NULL.
#' @param isParallel TRUE or FALSE for parallel computing. Default is TRUE
#' @param ncores Number of cores to be used for parallel computing. Only used
#' if isParallel=TRUE. If isParallel=TRUE and ncores=NULL, defaults to half the
#' size of the number of system cores.
#' @param gridMethod GridSearch or RandomSearch. Optimize tuning parameters over
#' full grid or random grid. Default is RandomSearch.
#' @param AssignClassMethod Classification method. Either Joint or Separate. Joint uses
#' all discriminant vectors from \eqn{D} datasets to predict class membership. Separate predicts
#' class membership separately for each dataset. Default is Joint.
#' @param nfolds Number of cross validation folds. Default is 5.
#' @param ngrid Number of grid points for tuning parameters. Default is 8 for each
#' view if \eqn{D =2}. If \eqn{D > 2}, default is 5.
#' @param standardize TRUE or FALSE. If TRUE, data will be normalized to have mean zero
#' and variance one for each variable. Default is TRUE.
#' @param maxiteration Maximum iteration for the algorithm if not converged.Default is 20.
#' @param weight Balances separation and association. Default is 0.5.
#' @param thresh Threshold for convergence. Default is 0.001.
#'
#' @details The function will return several R objects, which can be assigned to a variable.
#' To see the results, use the “$" operator.
#'
#' @return A list with the following components:
#'  \item{sidaerror}{Estimated classification error. If testing data provided, this will
#'  be test classification error, otherwise, training error}
#'  \item{sidacorrelation}{Sum of pairwise RV coefficients. Normalized to be within 0 and 1, inclusive.}
#'  \item{hatalpha}{A list of estimated sparse discriminant vectors for each view.}
#'  \item{PredictedClass}{ Predicted class. If AssignClassMethod=’Separate’, this will be a
#'  \eqn{ntest \times D} matrix, with each column the predicted class for each data.}
#'  \item{PredictedClass.train}{Predicted class for train data. If AssignClassMethod=’Separate’, this will
#'  \eqn{ntrain \times D} matrix, with each column the predicted class for each data.}
#'  \item{optTau}{Optimal tuning parameters for each view, not including covariates, if available.}
#'  \item{gridValues}{Grid values used for searching optimal tuning parameters.}
#'  \item{AssignClassMethod}{Classification method used. Joint or Separate.}
#'  \item{gridMethod}{Grid method used. Either GridSearch or RandomSearch}
#'
#'
#' @seealso \code{\link{sida}}
#'
#' @references
#' Sandra E. Safo, Eun Jeong Min, and Lillian Haine (2022) , Sparse Linear
#' Discriminant Analysis for Multi-view Structured Data, Biometrics
#'
#' @importFrom foreach %dopar% foreach
#' @import RSpectra
#' @importFrom igraph spectrum decompose
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster
#' @importFrom CVXR diag
#'
#' @export
#' @examples
#' \dontrun{
#'  #call sida
#' data(sidaData)
#' ##---- call sida algorithm to estimate discriminant vectors, and predict on testing data
#'
#' Xdata=sidaData[[1]]
#' Y=sidaData[[2]]
#' Xtestdata=sidaData[[3]]
#' Ytest=sidaData[[4]]
#'
#' ##---- call cross validation
#' mycv=cvSIDA(Xdata,Y,withCov=FALSE,plotIt=FALSE, Xtestdata=Xtestdata,Ytest=Ytest,
#'             isParallel=FALSE,gridMethod='RandomSearch',
#'             AssignClassMethod='Joint',nfolds=5,ngrid=8,standardize=TRUE,
#'             maxiteration=20, weight=0.5,thresh=1e-03)
#'
#' #check output
#' test.error=mycv$sidaerror
#' test.correlation=mycv$sidacorrelation
#' optTau=mycv$optTau
#' hatalpha=mycv$hatalpha
#'
#' #Obtain more performance metrics (applicable to two classes only)
#'  #train metrics
#'  Y.pred=mycv$PredictedClass.train-1 #to get this in 0 and 1
#'  Y.train=Y-1 #to get this in 0 and 1
#'  train.metrics=PerformanceMetrics(Y.pred,Y.train,family='binomial')
#'
#'  print(train.metrics)
#'  #obtain predicted class
#'  Y.pred=mycv$PredictedClass-1 #to get this in 0 and 1
#'  Ytest.in=Ytest-1 #to get this in 0 and 1
#'  test.metrics=PerformanceMetrics(Y.pred,Ytest.in,family='binomial')
#'  print(test.metrics)
#'  }

cvSIDA=function(Xdata=Xdata,Y=Y,withCov=FALSE,plotIt=FALSE,
                Xtestdata=NULL,Ytest=NULL,isParallel=TRUE,ncores=NULL,
                gridMethod='RandomSearch',AssignClassMethod='Joint',
                nfolds=5,ngrid=8,standardize=TRUE,maxiteration=20,
                weight=0.5,thresh=1e-03){

  starttimeall=Sys.time()

  XdataOrig=Xdata
  XtestdataOrig=Xtestdata
  YOrig=Y
  YtestOrig=Ytest

  #check inputs for training data
  dsizes=lapply(Xdata, function(x) dim(x))
  n=dsizes[[1]][1]
  nsizes=lapply(Xdata, function(x) dim(x)[1])

  if(all(nsizes!=nsizes[[1]])){
    stop('The datasets  have different number of observations')
  }


  #check data
  if (is.list(Xdata)) {
    D = length(Xdata)
    if(D==1){
      stop("There should be at least two datasets")
    }
  } else {
    stop("Input data should be a list")
  }



  #If testing data are not provided, the default is to use training data
  if(is.null(Xtestdata)){
    Xtestdata=Xdata
    Ytest=Y
  }

  #check inputs for testing data
  ntestsizes=lapply(Xtestdata, function(x) dim(x)[1])
  if(all(ntestsizes!=ntestsizes[[1]])){
    stop('The testing datasets  have different number of observations')
  }


  if(is.null(withCov)){
    withCov=FALSE
  }

  if(is.null(plotIt)){
    plotIt=FALSE
  }

  if(is.null(standardize)){
    standardize=TRUE
  }


  #standardize if true
  if(standardize==TRUE){
    Xdata=lapply(Xdata,function(x)scale(x,center=TRUE,scale=TRUE))
    Xtestdata=lapply(Xtestdata,function(x)scale(x,center=TRUE,scale=TRUE))
  }

  if(is.null(gridMethod)){
    gridMethod='RandomSearch'
  }

  if(is.null(AssignClassMethod)){
    AssignClassMethod='Joint'
  }

  if(is.null(isParallel)){
    isParallel=TRUE
  }

  if(is.null(nfolds)){
    nfolds=5
  }

  if(is.null(ngrid)){
    ngrid=8
  }


  if(is.null(maxiteration)){
    maxiteration=20
  }

  if(is.null(weight)){
    weight=0.5
  }

  if(is.null(thresh)){
    thresh=1e-03
  }

  set.seed(1234)
  nK=length(unique(as.vector(Y))) -1

  nc=length(unique(as.vector(Y)))
  Nn=mat.or.vec(nc,1)
  foldid=list()
  for(i in 1:nc)
  {
    Nn[i]=sum(Y==i)
    mod1=Nn[i]%%nfolds
    if(mod1==0){
      foldid[[i]]=sample(c(rep(1:nfolds,times=floor(Nn[i])/nfolds)),Nn[i])
    }else if(mod1> 0){
      foldid[[i]]=sample(c(rep(1:nfolds,times=floor(Nn[i])/nfolds), 1:(Nn[i]%%nfolds)),Nn[i])
    }
  }

  foldid=unlist(foldid)

  #obtain tuning range common to all K


  if(withCov==TRUE){
    #Dnew=D-1
    Dnew=D
  }else if(withCov==FALSE){
    Dnew=D
  }

  if(Dnew>2){
    ngrid=5
  }
  starttimetune=Sys.time()
  print('Getting tuning grid values')
  myTauvec=sidatunerange(Xdata,Y,ngrid,standardize,weight,withCov)
  endtimetune=Sys.time()
  print('Completed at time')
  print(endtimetune-starttimetune)

  #define the grid
  mygrid=expand.grid(do.call(cbind,myTauvec))
  gridcomb=dim(mygrid)[1]
  if(gridMethod=='RandomSearch'){
    if(Dnew==2){
      ntrials=floor(0.2*gridcomb)}
    else if(Dnew>2){
      ntrials=floor(0.15*gridcomb)
    }
    mytune=sample(1:gridcomb, ntrials, replace = FALSE)
    gridValues=mygrid[mytune,]
  }else if(gridMethod=='GridSearch'){
    gridValues=mygrid
  }


  starttimeCV=Sys.time()
  CVOut=matrix(0, nfolds, nrow(gridValues))
  #cross validation
  if(isParallel==TRUE){
    cat("Begin", nfolds,"-folds cross-validation", "\n")
    doParallel::registerDoParallel()
    if(is.null(ncores)){

      # chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
      #
      # if (nzchar(chk) && chk == "TRUE") {
      #   # use 2 cores in CRAN/Travis/AppVeyor
      #   ncores <- 2L
      # } else {
      #   # use all cores in devtools::test()
      #   #ncores=parallel::detectCores()
      #   ncores=ceiling(ncores/2)
      # }

      ncores=parallel::detectCores()
      ncores=ceiling(ncores/2)
      }
    cl=parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    CVOut=matrix(0, nrow(gridValues), nfolds)
    mycv=foreach::foreach(i = 1:nrow(gridValues), .combine='rbind',.export=c('sida',
              'sidainner','myfastinner','myfastIDAnonsparse','mysqrtminv',
              'sidaclassify', 'sidatunerange','DiscriminantPlots',
              'CorrelationPlots'),.packages=c('CVXR','RSpectra')) %dopar% {
      mTau=sapply(1:D, function(itau) list(t(gridValues[,itau][i])))
      #cat("Begin CV-fold", i, "\n")
      CVOut[i,]= sapply(1:nfolds, function(j){
        testInd=which(foldid==j)
        testX=lapply(Xdata, function(x) x[testInd,])
        testY=Y[testInd]
        trainX=lapply(Xdata, function(x) x[-testInd,])
        trainY=Y[-testInd]
        mysida=sida(trainX,trainY,mTau,withCov,
                    Xtestdata=testX,testY,AssignClassMethod='Joint',
                    standardize,maxiteration,weight,thresh, plotIt=F)
        return(min(mysida$sidaerror))
      } )
    }
    CVOut=t(mycv)
    parallel::stopCluster(cl)
  }else if(isParallel==FALSE){
    CVOut=matrix(0, nfolds, nrow(gridValues))
    for (i in 1:nfolds){
      testInd=which(foldid==i)
      testX=lapply(Xdata, function(x) x[testInd,])
      testY=Y[testInd]
      trainX=lapply(Xdata, function(x) x[-testInd,])
      trainY=Y[-testInd]

      cat("Begin CV-fold", i, "\n")
      CVOut[i,]= sapply(1:nrow(gridValues), function(itau){
        mTau=sapply(1:D, function(d) list(t(gridValues[itau,][d])))
        mysida=sida(trainX,trainY,mTau,withCov,Xtestdata=testX,testY,
                    AssignClassMethod='Joint',standardize,maxiteration,
                    weight,thresh, plotIt=F)

        return(min(mysida$sidaerror))
      } )
    }
  }

  endtimeCV=Sys.time()
  print('Cross-validation completed at time')
  print(endtimeCV-starttimeCV)

  print('Getting Results......')
  #compute average classification error
  minEorrInd=max(which(colMeans(CVOut)==min(colMeans(CVOut))))
  optTau=gridValues[ minEorrInd,]

  #Apply on testing data
  moptTau=sapply(1:Dnew, function(i) list(t(gridValues[minEorrInd,][i])))
  mysida=sida(Xdata=Xdata,Y=Y,Tau=moptTau,withCov,Xtestdata=Xtestdata,
              Ytest=Ytest,AssignClassMethod,standardize,maxiteration,
              weight,thresh,plotIt=F)

  #Apply on training data
  mysidaTrain=sida(Xdata=Xdata,Y=Y,Tau=moptTau,withCov,Xtestdata=Xdata,
              Ytest=Y,AssignClassMethod,standardize,maxiteration,weight,thresh,
              plotIt=F)

  ss=list()
  #sum pairwise RV coefficients for testing data
  for(d in 1:D){
    dd=setdiff(seq(1, D, by= 1),d)
    #correlations
    sumCorr2=0;
    for (jj in 1:length(dd)){
      j=dd[jj];
      X1=Xtestdata[[d]]%*%mysida$hatalpha[[d]]
      X2=Xtestdata[[j]]%*%mysida$hatalpha[[j]]
      X1=scale(X1, center=TRUE,scale=FALSE)
      X2=scale(X2, center=TRUE,scale=FALSE)
      X1X2=t(X1)%*%X2/dim(X1)[1]
      X1X1=t(X1)%*%X1/dim(X1)[1]
      X2X2=t(X2)%*%X2/dim(X2)[1]
      sumcorr3=sum(CVXR::diag(X1X2%*%t(X1X2)))/(sqrt(sum(CVXR::diag(X1X1%*%X1X1)))*sqrt(sum(CVXR::diag(X2X2%*%X2X2))))
      sumCorr2=sumCorr2+sumcorr3
    }
    ss[[d]]=sumCorr2/length(dd)
  }

  sidacorrelation=sum(do.call(rbind,ss))/D

  ss=list()
  #sum pairwise RV coefficients for training data
  for(d in 1:D){
    dd=setdiff(seq(1, D, by= 1),d)
    #correlations
    sumCorr2=0;
    for (jj in 1:length(dd)){
      j=dd[jj];
      X1=Xdata[[d]]%*%mysida$hatalpha[[d]]
      X2=Xdata[[j]]%*%mysida$hatalpha[[j]]
      X1=scale(X1, center=TRUE,scale=FALSE)
      X2=scale(X2, center=TRUE,scale=FALSE)
      X1X2=t(X1)%*%X2/dim(X1)[1]
      X1X1=t(X1)%*%X1/dim(X1)[1]
      X2X2=t(X2)%*%X2/dim(X2)[1]
      sumcorr3=sum(CVXR::diag(X1X2%*%t(X1X2)))/(sqrt(sum(CVXR::diag(X1X1%*%X1X1)))*sqrt(sum(CVXR::diag(X2X2%*%X2X2))))
      sumCorr2=sumCorr2+sumcorr3
    }
    ss[[d]]=sumCorr2/length(dd)
  }

  sidacorrelation.train=sum(do.call(rbind,ss))/D

  #Produce discriminant and correlation plot if plotIt=T
  if(plotIt==TRUE){
    DiscriminantPlots(Xtestdata,Ytest,mysida$hatalpha)
    CorrelationPlots(Xtestdata,Ytest,mysida$hatalpha)
  }

  #print out some results
  cat("Estimated Test Classification Error is", mysida$sidaerror, "\n")
  cat("Estimated Train Classification Error is", mysidaTrain$sidaerror, "\n")

  cat("Estimated Test Correlation is", sidacorrelation, "\n")
  cat("Estimated Train Correlation is", sidacorrelation.train, "\n")

  #number of selected variables
  mysum=matrix(NA,nrow=D,ncol=1)
  hatalpha.temp=list()
  hatalpha=mysida$hatalpha
  for(d in 1:D){
    if(dim(hatalpha[[d]])[2] == 1){
      mysum[d,1]=sum(hatalpha[[d]]!=0)
    }else if (dim(hatalpha[[d]])[2] >1){
      hatalpha.temp[[d]]=rowSums(abs(hatalpha[[d]]))
      mysum[d,1]=sum(hatalpha.temp[[d]]!=0)
    }
    #cat("Number of nonzero coefficients in view", d, "is", sum(mysida$hatalpha[[d]]!=0), "\n")
    cat("Number of nonzero coefficients in view", d, "is", as.vector(mysum[d,1]), "\n")
  }

  endtimeall=Sys.time()
  print("Total time used is")
  print(endtimeall-starttimeall)

  result=list(CVOut=CVOut,sidaerror=mysida$sidaerror,sidacorrelation=sidacorrelation,sidaerror.train=mysidaTrain$sidaerror,
              sidacorrelation.train=sidacorrelation.train,
              hatalpha=mysida$hatalpha,PredictedClass=mysida$PredictedClass,
              PredictedClass.train=mysidaTrain$PredictedClass,
              optTau=moptTau,gridValues=gridValues, AssignClassMethod=AssignClassMethod,
              gridMethod=gridMethod,
              InputData=XdataOrig)
  class(result)="SIDA"

  return(result)
}



#' @title Tuning parameter grid values for sida
#'
#' @description Sida function to provide tuning parameter grid values
#'  for each view, not including covariates, if available. It is
#'  recommended to use this to get lower and upper bounds of tuning
#'  parameters for each view that can be used in sida. This function is
#'  called by cvSIDA to select optimal tuning parameters.
#'
#' @param Xdata A list with each entry containing training views of size \eqn{n \times p_d},
#' where \eqn{d =1, \dots , D} views. Rows are samples and columns are variables. If
#' covariates are available, they should be included as a separate view, and
#' set as the last dataset. For binary or categorical covariates (assumes no
#' ordering), we suggest the use of indicator variables.
#' @param Y \eqn{n \times 1} vector of class membership. Numeric, coded as 1, 2, ....
#' @param ngrid Number of grid points for tuning parameters. Default is 10 for each
#' view if \eqn{D =2}. If \eqn{D > 2}, default is 5.
#' @param standardize TRUE or FALSE. If TRUE, data will be normalized to have mean zero
#' and variance one for each variable. Default is TRUE.
#' @param weight Balances separation and association. Default is 0.5.
#' @param withCov TRUE or FALSE if covariates are available. If TRUE, please set
#' all covariates as one dataset and should be the last dataset. For binary and
#' categorical variables, use indicator matrices/vectors. Default is FALSE.
#'
#' @details The function will return an R object with grid values for each data,
#'  not including covariates, if available. To see the results, use the “$" operator.
#'
#' @return An R object containing the following information:
#' \item{Tauvec}{grid values for each data, not including covariates, if available.}
#'
#'
#' @seealso \code{\link{sida}}
#'
#' @references
#' Sandra E. Safo, Eun Jeong Min, and Lillian Haine (2022) , Sparse Linear
#' Discriminant Analysis for Multi-view Structured Data, Biometrics.
#'
#' @importFrom igraph spectrum decompose
#'
#' @export
#' @examples
#' \dontrun{
#'  #call sida
#' data(sidaData)
#' ##---- call sida algorithm to estimate discriminant vectors, and predict on testing data
#'
#' Xdata=sidaData[[1]]
#' Y=sidaData[[2]]
#' Xtestdata=sidaData[[3]]
#' Ytest=sidaData[[4]]
#'
#' #call sidatunerange to get range of tuning parameter
#' ngrid=10
#' mytunerange=sidatunerange(Xdata,Y,ngrid,standardize=TRUE,weight=0.5,withCov=FALSE)
#'
#' # an example with Tau set as the lower bound
#' Tau=c(mytunerange$Tauvec[[1]][1], mytunerange$Tauvec[[2]][1])
#' mysida=sida(Xdata,Y,Tau,withCov=FALSE,Xtestdata=Xtestdata,Ytest=Ytest,AssignClassMethod='Joint',
#'             plotIt=FALSE, standardize=TRUE,maxiteration=20,weight=0.5,thresh= 1e-03)
#'
#' test.error=mysida$sidaerror
#' test.correlation=mysida$sidacorrelation
#'
#' #estimated discriminant vectors and predicted class
#' hatalpha=mysida$hatalpha
#'
#' predictedClass=mysida$PredictedClass
#'#obtain more performance metrics (applicable to two classes)
#'
#'  #train metrics
#'  Y.pred=mysida$PredictedClass.train-1 #to get this in 0 and 1
#'  Y.train=Y-1 #to get this in 0 and 1
#'  train.metrics=PerformanceMetrics(Y.pred,Y.train,family='binomial')
#'  print(train.metrics)
#'
#'  #obtain test predicted class
#'  Y.pred=mysida$PredictedClass-1 #to get this in 0 and 1
#'  Ytest.in=Ytest-1 #to get this in 0 and 1
#'  test.metrics=PerformanceMetrics(Y.pred,Ytest.in,family='binomial')
#'  print(test.metrics)
#'  }

sidatunerange=function(Xdata,Y,ngrid=10,standardize=TRUE,
                       weight=0.5,withCov=TRUE){

  #check size of each data

  dsizes=lapply(Xdata, function(x) dim(x))
  n=dsizes[[1]][1]
  p=lapply(Xdata, function(x) dim(x)[2])
  D=length(dsizes)

  if(is.null(withCov)){
    withCov=FALSE
  }

  if(is.null(ngrid)){
    ngrid=8
  }

  if(is.null(standardize)){
    standardize=TRUE
  }

  if(is.null(weight)){
    weight=0.5
  }

  # if(withCov==TRUE){
  #   D=D-1
  # }

  nK=length(unique(as.vector(Y))) -1


  #standardize if true
  if(standardize==TRUE){
    Xdata=lapply(Xdata,function(x)scale(x,center=TRUE,scale=TRUE))
  }


  #obtain nonsparse solutions
  mynsparse=myfastIDAnonsparse(Xdata,Y,weight)
  myfinner=myfastinner(Xdata,Y,mynsparse$sqrtminvmat,mynsparse$myalphaoldmat,mynsparse$tildealphamat, weight)


  #obtain upper and lower bounds
  ubx=lapply(myfinner$SepAndAssocd, function(x) CVXR::norm(x,'i')/1.2)
  lbx=lapply(1:D, function(x) 1.2*sqrt(log(p[[x]])/n)*ubx[[x]])
  ubx=lapply(1:D, function(x) ubx[[x]])

  #tuning range for each data
  Taugrid=list()
  cc=lapply(1, function(x1,x2)  cbind(lbx,ubx))

  cc=as.matrix(do.call(rbind,cc))
  for(d in 1:D){
    Taugrid[[d]]=seq(as.numeric(cc[d,1]),as.numeric(cc[d,2]),length.out=(ngrid+1))
  }

  # up to 25% sparsity

  myperx=lapply(Taugrid, function(x) quantile(x[1:ngrid], c(.1, .15, .2, .25, .35, .45), type=5))#similar to matlab
  myperx2=do.call(rbind,myperx)
  for(loc in 1:6){
    mTaux=sapply(1:D, function(i) list(t(myperx2[i,loc])))
    myres=sidainner(Xdata,Y,mynsparse$sqrtminvmat,mynsparse$myalphaoldmat,mynsparse$tildealphamat, mynsparse$tildelambda,mTaux,weight,withCov)
    nnz=sapply(1:D, function(i) list(colSums(myres$hatalpha[[i]]!=0)/dsizes[[i]][2]))
    nnz=cbind(c(do.call(rbind,nnz)))
    #print(nnz)
    if(all(nnz<=0.25)){
      break
    }
  }
  #final grid
  Tauvec=sapply(1:D, function(i) list(seq(as.numeric(t(myperx2[i,loc])),as.numeric(ubx[[i]]),len=(ngrid+1))))
  Tauvec=sapply(1:D, function(x) list(Tauvec[[x]][1:ngrid]))

  result=list(Tauvec=Tauvec)
  return(result)
}


#' @title Classification approach for SIDA and SIDANet
#'
#' @description Performs classification using nearest centroid on separate or
#' combined estimated discriminant vectors, and predicts class membership.
#'
#' @param hatalpha A list of estimated sparse discriminant vectors for each
#' view. This may be obtained from sida or cvSIDA.
#' @param Xtestdata A list with each entry containing testing views of size
#' \eqn{ntest \times p_d}, where \eqn{d =1, \dots, D} views. Rows are samples and
#' columns are variables. The order of the list should be the same as the
#' order for the training data, Xdata. If covariates are available, they
#' should be included as a separate view, and set as the last dataset.
#' For binary or categorical covariates (assumes no ordering), we suggest
#' the use of indicator variables. If you want to obtain training error,
#' set as Xdata.
#' @param Xdata A list with each entry containing training views of size
#'  \eqn{n \times p_d}, where \eqn{d = 1, \dots, D} views. Rows are samples and
#'  columns are variables. If covariates are available, they should be included
#'  as a separate view, and set as the last dataset. For binary or categorical
#'  covariates (assumes no ordering), we suggest the use of indicator variables.
#' @param Y \eqn{n \times 1} vector of class membership. Same size as the number of training samples.
#' Numeric, coded as 1, 2, ....
#' @param AssignClassMethod Classification method. Either Joint or Separate.
#' Joint uses all discriminant vectors from D datasets to predict class membership.
#' Separate predicts class membership separately for each dataset. Default is Joint.
#' @param standardize TRUE or FALSE. If TRUE, data will be normalized to have mean
#' zero and variance one for each variable. Default is TRUE.
#'
#'
#' @return An R object containing the following information:
#' \item{PredictedClass}{Predicted class. If AssignClassMethod=’Separate’, this will
#' be a \eqn{ntest × D} matrix, with each column the predicted class for each data.}
#' \item{AssignClassMethod}{Classification method used. Either Joint or Separate.}
#'
#' @seealso \code{\link{cvSIDA}} \code{\link{sida}}
#'
#' @references
#' Sandra E. Safo, Eun Jeong Min, and Lillian Haine (2022) , Sparse Linear
#' Discriminant Analysis for Multi-view Structured Data, Biometrics.
#'
#' @export
#' @examples
#' \dontrun{
#'  #call sida
#' data(sidaData)
#' ##---- call sida algorithm to estimate discriminant vectors, and predict on testing data
#'
#' Xdata=sidaData[[1]]
#' Y=sidaData[[2]]
#' Xtestdata=sidaData[[3]]
#' Ytest=sidaData[[4]]
#'
#' #call sidatunerange to get range of tuning paramater
#' ngrid=10
#'
#' mytunerange=sidatunerange(Xdata,Y,ngrid,standardize=TRUE,weight=0.5,withCov=FALSE)
#' # an example with Tau set as the lower bound
#' Tau=c(mytunerange$Tauvec[[1]][1], mytunerange$Tauvec[[2]][1])
#' mysida=sida(Xdata,Y,Tau,withCov=FALSE,Xtestdata=Xtestdata,Ytest=Ytest)
#' #classification with combined estimated vectors
#' mysida.classify.Joint=sidaclassify(mysida$hatalpha,Xtestdata,Xdata,Y,
#'                                    AssignClassMethod='Joint')
#' mysida.PredClass.Joint=mysida.classify.Joint$PredictedClass
#' #classification with separate estimated vectors
#' mysida.classify.Separate=sidaclassify(mysida$hatalpha,Xtestdata,Xdata,Y,
#'                                      AssignClassMethod='Separate')
#' mysida.PredClass.Separate=mysida.classify.Separate$PredictedClass
#' }
sidaclassify=function(hatalpha=hatalpha,Xtestdata=Xtestdata,Xdata=Xdata,Y=Y,AssignClassMethod='Joint', standardize=TRUE){

  #hatalpha is a list of d estimated SIDA vectors, for each view
  #Y is a vector of training observations
  if(is.null(AssignClassMethod)){
    AssignClassMethod='Joint'
  }

  if(is.null(standardize)){
    standardize=TRUE
  }

  #standardize if true
  if(standardize==TRUE){
    Xdata=lapply(Xdata,function(x)scale(x,center=TRUE,scale=TRUE))
    Xtestdata=lapply(Xtestdata,function(x)scale(x,center=TRUE,scale=TRUE))
  }

  D=length(Xtestdata)
  #classification
  Projtest=sapply(1:D, function(i) list(Xtestdata[[i]]%*%hatalpha[[i]]))
  Projtrain=sapply(1:D, function(i) list(Xdata[[i]]%*%hatalpha[[i]]))


  nc= length(unique(as.vector(Y)))
  ntest=dim(Projtest[[1]])[1]

  PredclassSeparate=list()
  if(AssignClassMethod=='Separate'){
    for(d in 1:D){
      ProjXtestdatad=Projtest[[d]]
      ProjXdatad=Projtrain[[d]]
      ProjXdata1d=cbind(Y,ProjXdatad)
      Projmv=stats::aggregate(ProjXdata1d[,-1],list(ProjXdata1d[,1]),mean)
      distv=list()
      jrep=list()
      for(j in 1: nc){
        rProjm=matrix( rep(Projmv[j,-1],times= ntest), ncol=ncol(ProjXdatad), byrow=TRUE)
        #euclidean distance
        sqdiff=(ProjXtestdatad-as.numeric(rProjm))^2
        dist1=rowSums(sqdiff)^0.5
        jrep[[j]]=j*rep(1,times=ntest)
        distv[[j]]=dist1
      }
      distv=do.call(cbind,distv)
      dim(distv)=c(nc*nrow(distv),1)

      jrep=do.call(cbind,jrep)
      dim(jrep)=c(nc*nrow(jrep),1)

      distv=cbind(jrep, distv)

      #The following code outputs the assigned class
      rdistvX1=matrix(distv[,-1], nrow=ntest,ncol=nc)
      minX1=apply(rdistvX1,1,min)
      minind=which(rdistvX1==minX1, arr.ind=TRUE) #minimum indices
      predclassX1=minind[order(minind[,1]),2]
      PredclassSeparate[[d]]=predclassX1
    }
    Predclass=do.call(cbind,PredclassSeparate)
  }
  else if(AssignClassMethod=='Joint') {
    #classification for joint
    ProjtestJoint=do.call(cbind,Projtest)
    ProjtrainJoint=do.call(cbind, Projtrain)

    ntest=dim(ProjtestJoint)[1]
    ProjtrainJointX1=cbind(Y,ProjtrainJoint)
    Projmv=stats::aggregate(ProjtrainJointX1[,-1],list(ProjtrainJointX1[,1]),mean)

    distv=list()
    jrep=list()
    for(j in 1: nc){
      rProjm=matrix( rep(Projmv[j,-1],times= ntest), ncol=ncol(ProjtrainJoint), byrow=TRUE)
      #euclidean distance
      sqdiff=(ProjtestJoint-as.numeric(rProjm))^2
      dist1=rowSums(sqdiff)^0.5
      jrep[[j]]=j*rep(1,times=ntest)
      distv[[j]]=dist1
    }
    distv=do.call(cbind,distv)
    dim(distv)=c(nc*nrow(distv),1)

    jrep=do.call(cbind,jrep)
    dim(jrep)=c(nc*nrow(jrep),1)

    distv=cbind(jrep, distv)


    #The following code outputs the assigned class
    rdistvX1=as.data.frame(matrix(distv[,-1], nrow=ntest,ncol=nc))
    minX1=apply(as.matrix(rdistvX1),1,min)
    minind=which(rdistvX1==minX1, arr.ind=TRUE) #minimum indices
    Predclass=minind[order(minind[,1]),2]
    #  Predclass=predclassJ
  }
  result=list(PredictedClass=Predclass, AssignClassMethod=AssignClassMethod)
  return(result)
}
