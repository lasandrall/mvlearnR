


#' @title 2-step supervised SELPCCA
#'
#' @description Performs n-fold cross validation to select optimal
#' tuning parameters for SELPCCA based on training data. Then uses the
#' results to build a GLM or survival model for a pre-specified outcome.
#'
#' @param Xdata1 A matrix of size n × p for first dataset. Rows are
#'        samples and columns are variables.
#' @param Xdata2 A matrix of size n × q for second dataset. Rows are
#'        samples and columns are variables.
#' @param Y A vector of size n for the outcome. Continuous outcomes
#'        do not have to be centered or scaled. If family="survival",
#'        Y is a vector of size n indicating the time at which
#'        the event occurred or the observation was censored. See
#'        'event' for more information on how to use function for
#'        a survival outcome.
#' @param fitselpCCA The output of cvselpscca() function or multiplescca(). If
#'        NULL, the algorithm will fit a cvselpscca model.
#' @param family A string to denote the type of prediction model to build.
#'        Options are "gaussian", "binomial", "poisson", or "survival".
#'        When family="survival", a proportional Cox model will be fitted.
#'        Otherwise a generalized linear model will be used.
#' @param event A vector of size n needed when family="survival" to denote whether
#'        or not the event of interest occurred at timepoint Y. Let
#'        event=NULL when family does not equal "survival".
#' @param model.separately A boolean to denote whether or not to use separate
#'        prediction models for Xdata1 and Xdata2. When
#'        model.separately=FALSE, a single model will be fit using
#'        the output for both datasets.
#' @param ncancorr  Number of canonical correlation vectors. Default is 1.
#' @param CovStructure Covariance structure to use in estimating sparse
#'        canonical correlation vectors. Either "Iden" or "Ridge".
#'        Iden assumes the covariance matrix for each dataset
#'        is identity. Ridge uses the sample covariance for each
#'        dataset. See reference article for more details.
#' @param isParallel TRUE or FALSE for parallel computing. Default is TRUE.
#' @param ncores Number of cores to be used for parallel computing. Only
#'        used if isParallel=TRUE. If isParallel=TRUE and ncores=NULL,
#'        defaults to half the size of the number of system cores.
#' @param nfolds Number of cross validation folds. Default is 5.
#' @param ngrid Number of grid points for tuning parameters. Default
#'        is 10 for each dataset.
#' @param standardize TRUE or FALSE. If TRUE, data will be normalized to
#'        have mean zero and variance one for each variable. Note
#'        that this only standardizes Xdata1 and Xdata2. Y will not
#'        be standardized. Default is TRUE.
#' @param thresh Threshold for convergence. Default is 0.0001.
#' @param maxiteration Maximum iteration for the algorithm if not converged.
#' Default is 20.
#' @param showProgress A boolean for whether or not the function
#'        should display text output at various stages in
#'        the function to indicate progress. Default is TRUE.
#'
#'
#' @details The function will return several R objects, which can be assigned
#'  to a variable. To see the results, use the “$" operator.
#'
#' @return The output is a list containing the following components.
#' \item{selp.fit}{The output of the cvselpscca() function.}
#' \item{mod.fit}{The output of the glm() or coxph() regression model.}
#' \item{data.matrix}{The data matrix that was used to build the regression model.}
#' \item{family}{The type of outcome specified.}
#'
#' @seealso \code{\link{cvselpscca}}
#'
#' @references
#' Sandra E. Safo, Jeongyoun Ahn, Yongho Jeon, and Sungkyu Jung (2018),
#'  Sparse Generalized Eigenvalue Problem with Application to Canonical
#'  Correlation Analysis for Integrative Analysis of Methylation and Gene
#'  Expression Data. Biometrics
#' @export
#' @examples
#' ##---- read in data
#' data(sidaData)
#'
#' Xdata1=sidaData[[1]][[1]]
#' Xdata2=sidaData[[1]][[2]]
#' Y=sidaData[[2]]-1
#'
#' myresult=selpscca.pred(Xdata1, Xdata2, Y,fitselpCCA=NULL, family="binomial",
#'              event=NULL,model.separately=FALSE, ncancorr=1,
#'              CovStructure="Iden", isParallel=TRUE, ncores=NULL,
#'              nfolds=5, ngrid=10, standardize=TRUE,thresh=0.0001,
#'              maxiteration=20, showProgress=T)
#'
#' #check output
#' train.correlation=myresult$selp.fit$maxcorr
#' optTau=myresult$selp.fit$optTau
#' hatalpha=myresult$selp.fit$hatalpha
#' hatbeta=myresult$selp.fit$hatbeta
#' predictionModel=summary(myresult$mod.fit)
#'
#'##Performance metrics
#'##Train Performance Metrics
#'newPredictions=predict(myresult, newdata=Xdata1, newdata2=Xdata2, type="response")
#'Y.pred=newPredictions$pred.mod #predicted probabilities
#'Y.train=Y
#'train.metrics=PerformanceMetrics(Y.pred,Y.train,family='binomial',isPlot=TRUE)
#'print(train.metrics)
#'
#'##Test Performance Metrics
#'Y.test=sidaData[[4]]-1
#'newPredictions=predict(myresult, newdata=Xtestdata1, newdata2=Xtestdata2, type="response")
#'Y.pred=newPredictions$pred.mod #predicted probabilities
#'test.metrics=PerformanceMetrics(Y.pred,Y.test,family='binomial',isPlot=TRUE)
#'print(test.metrics)


selpscca.pred <- function(Xdata1, Xdata2, Y, fitselpCCA=NULL, family="gaussian",
                          event=NULL,
                          model.separately=FALSE,
                          ncancorr=1, CovStructure="Iden",
                          isParallel=TRUE, ncores=NULL, nfolds=5,
                          ngrid=10, standardize=TRUE,thresh=0.0001,
                          maxiteration=20,
                          showProgress=T){

#Fit SELPscca
  if(showProgress){cat("Fitting SELPCCA Model \n")}
  if(is.null(fitselpCCA)){
    fit <- cvselpscca(Xdata1=Xdata1,Xdata2=Xdata2,
                               ncancorr=ncancorr,CovStructure=CovStructure,
                               isParallel=isParallel,ncores=ncores,nfolds=nfolds,
                               ngrid=ngrid,standardize=standardize,
                               thresh=thresh,maxiteration=maxiteration)
  }else{
   fit<-fitselpCCA
   ncancorr=dim(fit$hatalpha)[2]
  }


  if(showProgress){cat("Fitting Prediction Model \n")}
  if(model.separately==F){
    #Part 1: Create Data matrix
    selp.dat <- data.frame(Y=Y)
    Xdata1=apply(as.matrix(Xdata1), 2, as.numeric)
    Xdata2=apply(as.matrix(Xdata2), 2, as.numeric)
    scoresX1=Xdata1%*% fit$hatalpha
    scoresX2=Xdata2%*% fit$hatbeta
    selp.dat <- cbind(selp.dat, scoresX1, scoresX2)
    names(selp.dat) <- c("Y", paste0("X1_", 1:ncancorr), paste0("X2_", 1:ncancorr))
    if(family=="survival"){
      selp.dat$event <- event
    }

    #Part 2: Fit Model
    if(family=="gaussian"){
      mod.fit <- stats::lm(Y ~ ., data=selp.dat)
    }else if(family=="binomial"){
      mod.fit <- stats::glm(factor(Y) ~ ., data=selp.dat, family = stats::binomial)
    }else if(family == "poisson"){
      mod.fit <- stats::glm(Y ~ ., data=selp.dat, family = stats::poisson)
    }else if(family == "survival"){
      mod.fit <- survival::coxph(survival::Surv(Y, event) ~ ., data=selp.dat)
    }else{
      stop(paste0("Model ", family, " doesn't exist"))
    }
    all.mod.fits <- mod.fit
    all.selp.dat <- selp.dat

  }else if(model.separately){
    #Part 1: Create Data matrices
    selp.dat1 <- selp.dat2 <- data.frame(Y=Y)
    #matrix 1
    Xdata1=apply(as.matrix(Xdata1), 2, as.numeric)
    scoresX1=Xdata1%*% fit$hatalpha
    selp.dat1 <- cbind(selp.dat1, scoresX1)
    names(selp.dat1) <- c("Y", paste0("X1_", 1:ncancorr))
    #matrix 2
    Xdata2=apply(as.matrix(Xdata2), 2, as.numeric)
    scoresX2=Xdata2%*% fit$hatbeta
    selp.dat2 <- cbind(selp.dat2, scoresX2)
    names(selp.dat2) <- c("Y", paste0("X2_", 1:ncancorr))
    #add event if survival
    if(family=="survival"){
      selp.dat1$event <- event
      selp.dat2$event <- event
    }

    #Part 2: Fit Models
    if(family=="gaussian"){
      mod.fit1 <- stats::lm(Y ~ ., data=selp.dat1)
      mod.fit2 <- stats::lm(Y ~ ., data=selp.dat2)
    }else if(family=="binomial"){
      mod.fit1 <- stats::glm(factor(Y) ~ ., data=selp.dat1, family = stats::binomial)
      mod.fit2 <- stats::glm(factor(Y) ~ ., data=selp.dat2, family = stats::binomial)
    }else if(family == "poisson"){
      mod.fit1 <- stats::glm(Y ~ ., data=selp.dat1, family = stats::poisson)
      mod.fit2 <- stats::glm(Y ~ ., data=selp.dat2, family = stats::poisson)
    }else if(family == "survival"){
      mod.fit1 <- survival::coxph(survival::Surv(Y, event) ~ ., data=selp.dat1)
      mod.fit2 <- survival::coxph(survival::Surv(Y, event) ~ ., data=selp.dat2)
    }else{
      stop(paste0("Model ", family, " doesn't exist"))
    }
    all.mod.fits <- list(mod.fit1, mod.fit2)
    all.selp.dat <- list(selp.dat1, selp.dat2)
  }



  result <- list(selp.fit=fit,
                 mod.fit=all.mod.fits,
                 data.matrix=all.selp.dat,
                 family=family,
                 InputData=list(Xdata1,Xdata2),
                 method="selpscca.pred")

  class(result) <- "SELPCCA"

  return(result)
}




#' @title Prediction for out-of-sample data for SELPCCA predict
#'
#' @description A wrapper function to obtain the canonical variates for
#' an out-of-sample dataset based on a fitted SELPCCA model and then use
#' that information to predict Y based on the fitted GLM or Cox model.
#'
#' @param object  A fitted model of class SELPCCA
#' @param newdata A matrix of size \eqn{n \times p} for the first dataset. Rows are
#'        samples and columns are variables.
#' @param newdata2 A matrix of size \eqn{n \times q} for the second dataset. Rows are
#'        samples and columns are variables.
#' @param type See predict.glm() and predict.coxph() for type options and defaults.
#'
#' @return An object containing the output from predict.glm() or predict.coxph()
#'
#' @seealso \code{\link{cvSIDA}}  \code{\link{sidatunerange}}
#'
#' @export
#' @examples
#' ##---- read in data
#' data(sidaData)
#'
#' Xdata1=sidaData[[1]][[1]]
#' Xdata2=sidaData[[1]][[2]]
#' Xtestdata1=sidaData[[3]][[1]]
#' Xtestdata2=sidaData[[3]][[2]]
#' Y=sidaData[[2]]-1
#'
#'myresult=selpscca.pred(Xdata1, Xdata2, Y,fitselpCCA=NULL, family="binomial",
#'                       event=NULL,model.separately=FALSE, ncancorr=1,
#'                       CovStructure="Iden", isParallel=TRUE, ncores=NULL,
#'                       nfolds=5, ngrid=10, standardize=TRUE,thresh=0.0001,
#'                       maxiteration=20, showProgress=T)
#'
#' #check output
#' train.correlation=myresult$selp.fit$maxcorr
#' optTau=myresult$selp.fit$optTau
#' hatalpha=myresult$selp.fit$hatalpha
#' hatbeta=myresult$selp.fit$hatbeta
#' predictionModel=summary(myresult$mod.fit)
#'
#'##Performance metrics
#'##Train Performance Metrics
#'newPredictions=predict(myresult, newdata=Xdata1, newdata2=Xdata2, type="response")
#'Y.pred=newPredictions$pred.mod #predicted probabilities
#'Y.train=Y
#'train.metrics=PerformanceMetrics(Y.pred,Y.train,family='binomial',isPlot=TRUE)
#'print(train.metrics)
#'
#'##Test Performance Metrics
#'Y.test=sidaData[[4]]-1
#'newPredictions=predict(myresult, newdata=Xtestdata1, newdata2=Xtestdata2, type="response")
#'Y.pred=newPredictions$pred.mod #predicted probabilities
#'test.metrics=PerformanceMetrics(Y.pred,Y.test,family='binomial',isPlot=TRUE)
#'print(test.metrics)

predict.SELPCCA <- function(object, newdata, newdata2, type="response"){
  newdata=apply(as.matrix(newdata), 2, as.numeric)
  newdata2=apply(as.matrix(newdata2), 2, as.numeric)
  scoresX1=newdata %*% object$selp.fit$hatalpha
  scoresX2=newdata2 %*% object$selp.fit$hatbeta

  if(length(object$mod.fit)==2){
    dat <- list()
    dat[[1]] <- as.data.frame(scoresX1)
    dat[[2]] <- as.data.frame(scoresX2)
    if(object$family == "survival"){
      dat[[1]]$event <- 0
      dat[[2]]$event <- 0
    }
    names(dat[[1]]) <- names(object$data.matrix[[1]])[-1]
    names(dat[[2]]) <- names(object$data.matrix[[2]])[-1]

    pred.mod <- list()
    if(object$family != "survival"){
      pred.mod[[1]] <- stats::predict(object$mod.fit[[1]], dat[[1]],
                               type=type)
      pred.mod[[2]] <- stats::predict(object$mod.fit[[2]], dat[[2]],
                               type=type)
    }else{
      if(type=="response"){type="lp"}
      pred.mod[[1]] <- stats::predict(object$mod.fit[[1]], dat, type=type)
      pred.mod[[2]] <- stats::predict(object$mod.fit[[2]], dat, type=type)
    }
  }else{
    dat <- as.data.frame(cbind(scoresX1, scoresX2))
    if(object$family == "survival"){
      dat$event <- 0
    }
    names(dat) <- names(object$data.matrix)[-1]

    if(object$family != "survival"){
      pred.mod <- stats::predict(object$mod.fit, dat,
                          type=type)
    }else{
      if(type=="response"){type="lp"}
      pred.mod <- stats::predict(object$mod.fit, dat, type=type)
    }
  }

  result=list(pred.mod=pred.mod,
              Xdata1=newdata, Xdata2=newdata2)
  class(result)="SELP-Predict"
  return(result)
}


