

#' @title Cross validation for Sparse Canonical Correlation Analysis
#'
#' @description Performs n-fold cross validation to select optimal
#' tuning parameters for SELPCCA based on training data. If you want
#' to apply optimal tuning parameters to testing data, you may also
#' use multiplescca.
#'
#' @param Xdata1 A matrix of size \eqn{n \times p} for first dataset.
#' Rows are samples and columns are variables.
#' @param Xdata2 A matrix of size \eqn{n \times q} for second dataset.
#' Rows are samples and columns are variables.
#' @param ncancorr Number of canonical correlation vectors. Default is 1.
#' @param CovStructure Covariance structure to use in estimating sparse
#' canonical correlation vectors. Either "Iden" or "Ridge". Iden assumes
#' the covariance matrix for each dataset is identity. Ridge uses the
#' sample covariance for each dataset. See reference article for more details.
#' @param isParallel TRUE or FALSE for parallel computing. Default is TRUE.
#' @param ncores Number of cores to be used for parallel computing. Only used
#' if isParallel=TRUE. If isParallel=TRUE and ncores=NULL, defaults to half
#' the size of the number of system cores.
#' @param nfolds Number of cross validation folds. Default is 5.
#' @param ngrid Number of grid points for tuning parameters. Default is 10
#' for each dataset.
#' @param standardize TRUE or FALSE. If TRUE, data will be normalized to have
#'  mean zero and variance one for each variable. Default is TRUE.
#' @param maxiteration Maximum iteration for the algorithm if not converged.
#' Default is 20.
#' @param thresh Threshold for convergence. Default is 0.0001.
#'
#' @details The function will return several R objects, which can be assigned
#'  to a variable. To see the results, use the “$" operator.
#'
#' @return The output is a list containing the following components.
#' \item{hatalpha}{Estimated sparse canonical correlation vectors for first dataset.}
#' \item{hatbeta}{Estimated sparse canonical correlation vectors for second dataset.}
#' \item{CovStructure}{Covariance structure used in estimating sparse canonical
#' correlation vectors. Ei- ther "Iden" or "Ridge".}
#' \item{optTau}{Optimal tuning parameters for each dataset.}
#' \item{maxcorr}{Estimated canonical correlation coefficient.}
#' \item{tunerange}{Grid values for each dataset used for searching optimal
#' tuning paramters.}
#'
#' @seealso \code{\link{multiplescca}}
#'
#' @references
#' Sandra E. Safo, Jeongyoun Ahn, Yongho Jeon, and Sungkyu Jung (2018),
#'  Sparse Generalized Eigenvalue Problem with Application to Canonical
#'  Correlation Analysis for Integrative Analysis of Methylation and Gene
#'  Expression Data. Biometrics
#' @export
#' @examples
#' ##---- read in data
#' data(selpData)
#'
#' Xdata1=selpData[[1]]
#' Xdata2=selpData[[2]]
#'
#' ##---- call cross validation to estimate first canonical correlation vectors
#' ncancorr=1
#' mycv=cvselpscca(Xdata1=Xdata1,Xdata2=Xdata2,ncancorr=ncancorr,CovStructure="Iden",
#'                 isParallel=TRUE,ncores=NULL,nfolds=5,ngrid=10,
#'                 standardize=TRUE,thresh=0.0001,maxiteration=20)
#'
#' #check output
#' train.correlation=mycv$maxcorr
#' optTau=mycv$optTau
#' hatalpha=mycv$hatalpha
#' hatbeta=mycv$hatbeta
#'
#' #obtain correlation plot using training data
#' scoresX1=Xdata1%*% hatalpha
#' scoresX2=Xdata2%*% hatbeta
#' plot(scoresX1, scoresX2,lwd=3,
#'      xlab=paste(
#'        "First Canonical correlation variate for dataset", 1),
#'      ylab=paste("First Canonical correlation variate for dataset", 2),
#'      main=paste("Correlation plot for datasets",1, "and" ,2, ",", "\u03C1 =", mycv$maxcorr))
#'
#' #obtain correlation plot using testing data
#' Xtestdata1=selpData[[3]]
#' Xtestdata2=selpData[[4]]
#' scoresX1=Xtestdata1%*%hatalpha
#' scoresX2=Xtestdata2%*%hatbeta
#' mytestcorr=round(abs(cor(Xtestdata1%*%hatalpha,Xtestdata2%*%hatbeta)),3)
#'
#' plot(scoresX1, scoresX2,lwd=3,xlab=paste(
#'   "First Canonical correlation variate for dataset", 1),
#'   ylab=paste("First Canonical correlation variate for dataset", 2),
#'   main=paste("Correlation plot for datasets",1, "and" ,2, ",", "\u03C1 =", mytestcorr))
#'
cvselpscca=function(Xdata1=Xdata1,Xdata2=Xdata2,ncancorr=ncancorr,CovStructure="Iden",isParallel=TRUE,ncores=NULL,nfolds=5,ngrid=10,standardize=TRUE,thresh=1e-04,maxiteration=20){

  Xdata1Orig=Xdata1
  Xdata2Orig=Xdata2
  if(is.null(Xdata2)){
    stop('There must be two datasets')
  }

  nX=dim(Xdata1)[1]
  nY=dim(Xdata2)[1]
  n=nX
  if(nX!=nY){
    stop('Xdata1 and Xdata2 have different number of observations')
  }
  #set defaults

  if(is.null(ncancorr)){
    ncancorr=1
  }

  if(is.null(CovStructure)){
    CovStructure="Iden"
  }
  if(is.null(isParallel)){
    isParallel=TRUE
  }

  if(is.null(nfolds)){
    nfolds=5
  }

  if(is.null(ngrid)){
    ngrid=10
  }


  if(is.null(maxiteration)){
    maxiteration=20
  }

  if(is.null(thresh)){
    thresh=1e-04
  }

  if(is.null(standardize)){
    standardize=TRUE
  }

  if(standardize==TRUE){
    Xdata1=scale(Xdata1)
    Xdata2=scale(Xdata2)
  }

  Tauxvec=list()
  Tauyvec=list()
  tunerange=cvtunerange(Xdata1,Xdata2,ncancorr,CovStructure)
  Tauxrange=tunerange$TauX1range
  Tauyrange=tunerange$TauX2range

  for(j in 1:ncancorr){
    Tauxvec[[j]]=10^(seq(log(Tauxrange[j,1])/log(10),log(Tauxrange[j,2])/log(10), (log(Tauxrange[j,2])/log(10)-log(Tauxrange[j,1])/log(10))/ngrid))
    Tauyvec[[j]]=10^(seq(log(Tauyrange[j,1])/log(10),log(Tauyrange[j,2])/log(10), (log(Tauyrange[j,2])/log(10)-log(Tauyrange[j,1])/log(10))/ngrid))
  }
  Tauxvec=matrix(unlist(Tauxvec),nrow=length(Tauxvec),byrow=TRUE)
  Tauyvec=matrix(unlist(Tauyvec),nrow=length(Tauyvec),byrow=TRUE)

  myoptTau=list()

  set.seed(1)
  if((n%%nfolds)>=1){foldid=sample(c(rep(1:nfolds,floor(n/nfolds)),1:(n%%nfolds)),n,replace=FALSE)}else
  {foldid=sample(rep(1:nfolds,floor(n/nfolds)),n,replace=FALSE)}


  diff_alpha=1
  diff_beta=1
  iter=0
  reldiff=1

  #ncores=c()
  while((iter<maxiteration) && min(reldiff, max(diff_alpha,diff_beta))> thresh){
    iter=iter+1
    print(paste("Current iteration is", iter))
    mycorr2=list()
    if(dim(Tauyvec)[1]==1){
      Tauy=t(apply(t(Tauyvec[,1:ngrid]),1,stats::median))
    }else{Tauy=as.matrix(apply(Tauyvec[,1:ngrid],1,stats::median))}

    if(isParallel==TRUE){ #if using parallel

      doParallel::registerDoParallel()
      if(is.null(ncores)){
        ncores=parallel::detectCores()
        ncores=ceiling(ncores/2)}
      cl=parallel::makeCluster(ncores)
      doParallel::registerDoParallel(cl)

      mycorr2=foreach::foreach(jj=1:(dim(Tauxvec)[2]-1),.export=c('myfastnonsparsecca','mysqrtminv','selpscca','minv','fastsvd'),.packages=c('CVXR')) %dopar%{
        mycorr=list()
        if(dim(Tauxvec)[1]==1){
          Taux=t(Tauxvec[,jj])
        }else{Taux=as.matrix(Tauxvec[,jj])}

        for(ii in 1:nfolds){
          which.row=foldid==ii
          myfastcca=myfastnonsparsecca(Xdata1[-which(which.row),],Xdata2[-which(which.row),],CovStructure)
          tildeA=myfastcca$tildeA
          tildeB=myfastcca$tildeB
          tilderho=myfastcca$tilderho
          Ux=myfastcca$Ux
          Sigma12r=myfastcca$Sigma12r
          Uy=myfastcca$Uy

          if(iter==1){
            myalphaold=as.matrix(tildeA[,1:ncancorr])
            mybetaold=as.matrix(tildeB[,1:ncancorr])
            tilderhoold=tilderho[1:ncancorr]
          }

          myselps=selpscca(Xdata1[-which(which.row),],Xdata2[-which(which.row),],mybetaold,myalphaold, tilderhoold,ncancorr,Ux,Sigma12r,Uy,Taux,Tauy,CovStructure)
          myhatalpha=myselps$myhatalpha
          myhatbeta=myselps$myhatbeta

          myUtest=Xdata1[which(which.row),]%*%myhatalpha
          myVtest=Xdata2[which(which.row),]%*%myhatbeta

          myUtrain=Xdata1[-which(which.row),]%*%myhatalpha
          myVtrain=Xdata2[-which(which.row),]%*%myhatbeta

          mycorrTrain=prod(abs(diag(stats::cor(myUtrain,myVtrain))))
          mycorrTest=prod(abs(diag(stats::cor(myUtest,myVtest))))

          mycorr[[ii]]=c(t(Taux), ii ,mycorrTrain, mycorrTest)
        }

        mycorr=matrix(unlist(mycorr),nrow=length(mycorr),byrow=TRUE)
        isnanTrain=is.nan(mycorr[,ncancorr+2])
        isnanTest=is.nan(mycorr[,ncancorr+3])
        return(c(t(Taux),mean(mycorr[!isnanTrain,ncancorr+2]), mean(mycorr[!isnanTest,ncancorr+3])))
      }
      mycorr2=matrix(unlist(mycorr2),nrow=length(mycorr2),byrow=TRUE)
      mydiffmean=abs(mycorr2[,ncancorr+1]-mycorr2[,ncancorr+2])
      if(length(stats::na.omit(mydiffmean[mydiffmean!=0]))==0){#If all are zeros
        optTaux=Tauxvec[,1]
      }else{
        row=max(which((mydiffmean==min(stats::na.omit(mydiffmean[mydiffmean!=0]))),arr.ind=TRUE))
        optTaux=Tauxvec[,row]
      }

      mycorr2=list()
      Taux=as.matrix(optTaux)
      mycorr2=foreach::foreach(jj=1:(dim(Tauyvec)[2]-1),.export=c('myfastnonsparsecca','mysqrtminv','selpscca','minv','fastsvd'),.packages=c('CVXR')) %dopar%{
        mycorr=list()
        if(dim(Tauyvec)[1]==1){
          Tauy=t(Tauyvec[,jj])
        }else{Tauy=as.matrix(Tauyvec[,jj])}

        for(ii in 1:nfolds){
          which.row=foldid==ii
          myfastcca=myfastnonsparsecca(Xdata1[-which(which.row),],Xdata2[-which(which.row),],CovStructure)
          tildeA=myfastcca$tildeA
          tildeB=myfastcca$tildeB
          tilderho=myfastcca$tilderho
          Ux=myfastcca$Ux
          Sigma12r=myfastcca$Sigma12r
          Uy=myfastcca$Uy

          if(iter==1){
            myalphaold=as.matrix(tildeA[,1:ncancorr])
            mybetaold=as.matrix(tildeB[,1:ncancorr])
            tilderhoold=tilderho[1:ncancorr]
          }

          myselps=selpscca(Xdata1[-which(which.row),],Xdata2[-which(which.row),],mybetaold,myalphaold, tilderhoold,ncancorr,Ux,Sigma12r,Uy,Taux,Tauy,CovStructure)
          myhatalpha=myselps$myhatalpha
          myhatbeta=myselps$myhatbeta

          myUtest=Xdata1[which(which.row),]%*%myhatalpha
          myVtest=Xdata2[which(which.row),]%*%myhatbeta

          myUtrain=Xdata1[-which(which.row),]%*%myhatalpha
          myVtrain=Xdata2[-which(which.row),]%*%myhatbeta

          mycorrTrain=prod(abs(diag(stats::cor(myUtrain,myVtrain))))
          mycorrTest=prod(abs(diag(stats::cor(myUtest,myVtest))))

          mycorr[[ii]]=c(t(Tauy), ii ,mycorrTrain, mycorrTest)
        }
        mycorr=matrix(unlist(mycorr),nrow=length(mycorr),byrow=TRUE)
        isnanTrain=is.nan(mycorr[,ncancorr+2])
        isnanTest=is.nan(mycorr[,ncancorr+3])
        return(c(t(Tauy),mean(mycorr[!isnanTrain,ncancorr+2]), mean(mycorr[!isnanTest,ncancorr+3])))
      }
      mycorr2=matrix(unlist(mycorr2),nrow=length(mycorr2),byrow=TRUE)
      mydiffmean=abs(mycorr2[,ncancorr+1]-mycorr2[,ncancorr+2])
      if(length(stats::na.omit(mydiffmean[mydiffmean!=0]))==0){#If all are zeros
        optTauy=Tauyvec[,1]
      }else{
        row=max(which((mydiffmean==min(stats::na.omit(mydiffmean[mydiffmean!=0]))),arr.ind=TRUE))
        optTauy=Tauyvec[,row]}

      #stop cluster
      parallel::stopCluster(cl)
    }else{ #if not using parallel
      for(jj in 1:(dim(Tauxvec)[2]-1)){
        mycorr=list()
        if(dim(Tauxvec)[1]==1){
          Taux=t(Tauxvec[,jj])
        }else{Taux=as.matrix(Tauxvec[,jj])}

        for(ii in 1:nfolds){
          which.row=foldid==ii
          myfastcca=myfastnonsparsecca(Xdata1[-which(which.row),],Xdata2[-which(which.row),],CovStructure)
          tildeA=myfastcca$tildeA
          tildeB=myfastcca$tildeB
          tilderho=myfastcca$tilderho
          Ux=myfastcca$Ux
          Sigma12r=myfastcca$Sigma12r
          Uy=myfastcca$Uy

          if(iter==1){
            myalphaold=as.matrix(tildeA[,1:ncancorr])
            mybetaold=as.matrix(tildeB[,1:ncancorr])
            tilderhoold=tilderho[1:ncancorr]
          }

          myselps=selpscca(Xdata1[-which(which.row),],Xdata2[-which(which.row),],mybetaold,myalphaold, tilderhoold,ncancorr,Ux,Sigma12r,Uy,Taux,Tauy,CovStructure)
          myhatalpha=myselps$myhatalpha
          myhatbeta=myselps$myhatbeta

          myUtest=Xdata1[which(which.row),]%*%myhatalpha
          myVtest=Xdata2[which(which.row),]%*%myhatbeta

          myUtrain=Xdata1[-which(which.row),]%*%myhatalpha
          myVtrain=Xdata2[-which(which.row),]%*%myhatbeta

          mycorrTrain=prod(abs(diag(stats::cor(myUtrain,myVtrain))))
          mycorrTest=prod(abs(diag(stats::cor(myUtest,myVtest))))

          mycorr[[ii]]=c(t(Taux), ii ,mycorrTrain, mycorrTest)
        }

        mycorr=matrix(unlist(mycorr),nrow=length(mycorr),byrow=TRUE)
        isnanTrain=is.nan(mycorr[,ncancorr+2])
        isnanTest=is.nan(mycorr[,ncancorr+3])
        mycorr2[[jj]]=c(t(Taux),mean(mycorr[!isnanTrain,ncancorr+2]), mean(mycorr[!isnanTest,ncancorr+3]))
      }
      mycorr2=matrix(unlist(mycorr2),nrow=length(mycorr2),byrow=TRUE)
      mydiffmean=abs(mycorr2[,ncancorr+1]-mycorr2[,ncancorr+2])
      if(length(stats::na.omit(mydiffmean[mydiffmean!=0]))==0){#If all are zeros
        optTaux=Tauxvec[,1]
      }else{
        row=max(which((mydiffmean==min(stats::na.omit(mydiffmean[mydiffmean!=0]))),arr.ind=TRUE))
        optTaux=Tauxvec[,row]
      }

      mycorr2=list()
      Taux=as.matrix(optTaux)
      for(jj in 1:(dim(Tauyvec)[2]-1)){
        mycorr=list()
        if(dim(Tauyvec)[1]==1){
          Tauy=t(Tauyvec[,jj])
        }else{as.matrix(Tauyvec[,jj])}

        for(ii in 1:nfolds){
          which.row=foldid==ii
          myfastcca=myfastnonsparsecca(Xdata1[-which(which.row),],Xdata2[-which(which.row),],CovStructure)
          tildeA=myfastcca$tildeA
          tildeB=myfastcca$tildeB
          tilderho=myfastcca$tilderho
          Ux=myfastcca$Ux
          Sigma12r=myfastcca$Sigma12r
          Uy=myfastcca$Uy

          if(iter==1){
            myalphaold=as.matrix(tildeA[,1:ncancorr])
            mybetaold=as.matrix(tildeB[,1:ncancorr])
            tilderhoold=tilderho[1:ncancorr]
          }

          myselps=selpscca(Xdata1[-which(which.row),],Xdata2[-which(which.row),],mybetaold,myalphaold, tilderhoold,ncancorr,Ux,Sigma12r,Uy,Taux,Tauy,CovStructure)
          myhatalpha=myselps$myhatalpha
          myhatbeta=myselps$myhatbeta

          myUtest=Xdata1[which(which.row),]%*%myhatalpha
          myVtest=Xdata2[which(which.row),]%*%myhatbeta

          myUtrain=Xdata1[-which(which.row),]%*%myhatalpha
          myVtrain=Xdata2[-which(which.row),]%*%myhatbeta

          mycorrTrain=prod(abs(diag(stats::cor(myUtrain,myVtrain))))
          mycorrTest=prod(abs(diag(stats::cor(myUtest,myVtest))))

          mycorr[[ii]]=c(t(Tauy), ii ,mycorrTrain, mycorrTest)
        }
        mycorr=matrix(unlist(mycorr),nrow=length(mycorr),byrow=TRUE)
        isnanTrain=is.nan(mycorr[,ncancorr+2])
        isnanTest=is.nan(mycorr[,ncancorr+3])
        mycorr2[[jj]]=c(t(Tauy),mean(mycorr[!isnanTrain,ncancorr+2]), mean(mycorr[!isnanTest,ncancorr+3]))
      }
      mycorr2=matrix(unlist(mycorr2),nrow=length(mycorr2),byrow=TRUE)
      mydiffmean=abs(mycorr2[,ncancorr+1]-mycorr2[,ncancorr+2])
      if(length(stats::na.omit(mydiffmean[mydiffmean!=0]))==0){#If all are zeros
        optTauy=Tauyvec[,1]
      }else{
        row=max(which((mydiffmean==min(stats::na.omit(mydiffmean[mydiffmean!=0]))),arr.ind=TRUE))
        optTauy=Tauyvec[,row]}
    }

    Tauy=as.matrix(optTauy)
    myfastcca=myfastnonsparsecca(Xdata1 ,Xdata2,CovStructure)
    tildeA=myfastcca$tildeA
    tildeB=myfastcca$tildeB
    tilderho=myfastcca$tilderho

    myalphaold=as.matrix(tildeA[,1:ncancorr])
    mybetaold=as.matrix(tildeB[,1:ncancorr])
    tilderhoold=tilderho[1:ncancorr]

    if(iter>1){
      myalphaold=myalpha
      mybetaold=mybeta
      tilderhoold=mytilderho

      myalpha2old=myalpha2
      mybeta2old=mybeta2
    }

    myselpscca= selpscca(Xdata1,Xdata2,mybetaold, myalphaold, tilderhoold,ncancorr,myfastcca$Ux,myfastcca$Sigma12r,myfastcca$Uy,Taux,Tauy,CovStructure)
    myoptTau[[iter]]=cbind(iter*matrix(1,nrow=ncancorr,ncol=1),as.matrix(optTaux),as.matrix(optTauy))

    myalpha=myselpscca$myhatalpha
    mybeta=myselpscca$myhatbeta
    mytilderho=myselpscca$mycorrmat

    myalpha2=myselpscca$myalphaicon
    mybeta2=myselpscca$myalphaicon

    if( min(colSums(abs(myalpha)))==0 || min(colSums(abs(mybeta)))==0){
      myalpha=myalphaold
      mybeta=mybetaold
      break
    }

    myalpha2=myselpscca$myalphaicon
    mybeta2=myselpscca$myalphaicon

    diff_alpha=norm(myalpha-myalphaold,'f')^2 / norm(myalphaold,'f')^2 #
    diff_beta=norm(mybeta-mybetaold,'f')^2 / norm(mybetaold,'f')^2 #
    sumnew=norm(myalpha-myalphaold,'f')^2 + norm(mybeta-mybetaold,'f')^2 #
    sumold=norm(myalphaold,'f')^2+norm(mybetaold,'f')^2 #
    reldiff=sumnew/sumold #
  }



  print("Applying optimal tuning parameter on whole data")

  for(i in 1:length(myoptTau)){
    if(i==1){temp=matrix(unlist(myoptTau[[i]]),nrow=ncancorr,byrow=FALSE)}else{
      temp=rbind(temp,matrix(unlist(myoptTau[[i]]),nrow=ncancorr,byrow=FALSE))}}

  myoptTau=temp

  mymulticca=multiplescca(Xdata1,Xdata2,ncancorr,myoptTau,CovStructure,standardize,maxiteration,thresh)

  result=list(optTau=myoptTau,hatalpha=mymulticca$myalpha,
              hatbeta=mymulticca$mybeta,maxcorr=mymulticca$estCorr,CovStructure=CovStructure,
              tunerange=matrix(c(Tauxvec,Tauyvec),nrow = 2, byrow=TRUE), method="cvselpscca", InputData=list(Xdata1Orig,Xdata2Orig))

  class(result)="SELPCCA"
  return(result)
}


#' @title Tuning parameter range
#'
#' @description Obtain upper and lower bounds of tuning parameters
#' for each canonical correlation vector. It is recommended to use
#' cvselpscca to choose optimal tuning parameters for each dataset.
#'
#' @param Xdata1 A matrix of size \eqn{n \times p} for first dataset.
#' Rows are samples and columns are variables.
#' @param Xdata2 A matrix of size \eqn{n \times q} for second dataset.
#'  Rows are samples and columns are variables.
#' @param ncancorr Number of canonical correlation vectors. Default is one.
#' @param CovStructure Covariance structure to use in estimating sparse
#' canonical correlation vectors. Either "Iden" or "Ridge". Iden assumes
#' the covariance matrix for each dataset is identity. Ridge uses the
#' sample covariance for each dataset. See reference article for more details.
#' @param standardize TRUE or FALSE. If TRUE, data will be normalized to have
#' mean zero and variance one for each variable. Default is TRUE.
#'
#' @details The function will return tuning ranges for sparse estimation of
#' canonical correlation vectors. To see the results, use the “$" operator.
#'
#' @return The output is a list containing the following components.
#' \item{TauX1range}{A \eqn{ncancorr \times 2} matrix of upper and lower bounds of tuning
#' parameters for each canonical correlation vector for first dataset.}
#' \item{TauX2range}{A \eqn{ncancorr \times 2} matrix of upper and lower bounds of tuning
#'  parameters for each canonical correlation vector for second dataset.}
#'
#' @seealso \code{\link{multiplescca}}  \code{\link{cvselpscca}}
#'
#' @references
#' Sandra E. Safo, Jeongyoun Ahn, Yongho Jeon, and Sungkyu Jung (2018),
#'  Sparse Generalized Eigenvalue Problem with Application to Canonical
#'  Correlation Analysis for Integrative Analysis of Methylation and Gene
#'  Expression Data. Biometrics
#' @export
#' @examples
#' ##---- read in data
#' data(selpData)
#'
#' Xdata1=selpData[[1]]
#' Xdata2=selpData[[2]]
#'
#'   ##---- estimate first canonical correlation vectors
#' ncancorr=1
#'
#' #use cvtunerange for range of tuning parameters
#' mytunerange=cvtunerange(Xdata1=Xdata1,Xdata2=Xdata2,ncancorr=ncancorr,
#'                         CovStructure="Iden",standardize=TRUE)
#' print(mytunerange)
#'
#' #Fix Tau for first and second datasets as 1.1 and 1.0 respectively
#' Tau=matrix(c(1,1.2,1),nrow=1)
#' mysparsevectors=multiplescca(Xdata1=Xdata1,Xdata2=Xdata2,ncancorr=ncancorr,
#'                              Tau=Tau, CovStructure="Iden",standardize=TRUE,
#'                              maxiteration=20, thresh=0.0001)
#'
#' #example with two canonical correlation vectors
#' #use cvselpscca to obtain optimal tuning parameters
#' mycv=cvselpscca(Xdata1=Xdata1,Xdata2=Xdata2,ncancorr=ncancorr,
#'                 CovStructure="Iden",isParallel=TRUE,ncores=NULL,nfolds=5,
#'                 ngrid=10, standardize=TRUE,thresh=0.0001,maxiteration=20)
#'
#' Tau=mycv$optTau
#' mysparsevectors=multiplescca(Xdata1=Xdata1,Xdata2=Xdata2,ncancorr=ncancorr,
#'                              Tau=Tau, CovStructure="Iden",standardize=TRUE,maxiteration=20,
#'                              thresh=0.0001)
cvtunerange=function(Xdata1=Xdata1,Xdata2=Xdata2,ncancorr=ncancorr,CovStructure="Iden",standardize=TRUE){

  #set defaults
  #check data sizes

  p=dim(Xdata1)[2]
  n=dim(Xdata2)[1]
  q=dim(Xdata2)[2]

  if(dim(Xdata1)[1]!=dim(Xdata2)[1]){
    stop('The datasets  have different number of observations')
  }


  if(is.null(ncancorr)){
    ncancorr=1
  }

  if(is.null(CovStructure)){
    CovStructure="Iden"
  }

  if(is.null(standardize)){
    standardize=TRUE
  }

  if(standardize==TRUE){
    Xdata1=scale(Xdata1)
    Xdata2=scale(Xdata2)
  }


  myfastcca=myfastnonsparsecca(Xdata1,Xdata2,CovStructure)
  tildeA=myfastcca$tildeA
  tildeB=myfastcca$tildeB
  Ux=myfastcca$Ux
  Sigma12r=myfastcca$Sigma12r
  Uy=myfastcca$Uy

  mymaxlambdaalpha=matrix(rep(NA,ncancorr),nrow=ncancorr,ncol=1)
  mymaxlambdabeta=matrix(rep(NA,ncancorr),nrow=ncancorr,ncol=1)

  for(j in 1:ncancorr){
    mymaxlambdaalpha[j,1]= norm(Ux%*%Sigma12r%*%t(Uy)%*%tildeB[,j],"I")
    mymaxlambdabeta[j,1] = norm(Uy%*%t(Sigma12r)%*%t(Ux)%*%tildeA[,j],"I")
  }

  if(CovStructure=="Iden"){
    lb=mymaxlambdaalpha
    TauX1range=cbind(sqrt(log(p)/n)*lb, mymaxlambdaalpha/1.2)
    TauX2range=cbind(sqrt(log(q)/n)*lb, mymaxlambdabeta/1.2)
  }else if(CovStructure=="Ridge"){
    lb= mymaxlambdaalpha/3
    TauX1range=cbind(sqrt(log(p)/n)*lb, mymaxlambdaalpha/3)
    TauX2range=cbind(sqrt(log(q)/n)*lb, mymaxlambdabeta/3)
  }
  return(list(TauX1range=TauX1range,TauX2range=TauX2range))
}



#' @title Sparse canonical correlation vectors for fixed tuning parameters
#'
#' @description Performs n-fold cross validation to select optimal
#' tuning parameters for SELPCCA based on training data. If you want
#' to apply optimal tuning parameters to testing data, you may also
#' use multiplescca.
#'
#' @param Xdata1 A matrix of size \eqn{n \times p} for first dataset.
#' Rows are samples and columns are variables.
#' @param Xdata2 A matrix of size \eqn{n \times q} for second dataset.
#' Rows are samples and columns are variables.
#' @param ncancorr Number of canonical correlation vectors. Default is 1.
#' @param Tau A vector of matrix of fixed tuning parameters for each dataset.
#' @param CovStructure Covariance structure to use in estimating sparse
#' canonical correlation vectors. Either "Iden" or "Ridge". Iden assumes
#' the covariance matrix for each dataset is identity. Ridge uses the
#' sample covariance for each dataset. See reference article for more details.
#' @param standardize TRUE or FALSE. If TRUE, data will be normalized to have
#'  mean zero and variance one for each variable. Default is TRUE.
#' @param maxiteration Maximum iteration for the algorithm if not converged.
#' Default is 20.
#' @param thresh Threshold for convergence. Default is 0.0001.
#'
#' @details The function will return several R objects, which can be assigned
#'  to a variable. To see the results, use the “$" operator.
#'
#' @return The output is a list containing the following components.
#' \item{hatalpha}{Estimated sparse canonical correlation vectors for first dataset.}
#' \item{hatbeta}{Estimated sparse canonical correlation vectors for second dataset.}
#' \item{maxcorr}{Estimated canonical correlation coefficient.}
#'
#' @seealso \code{\link{cvselpscca}}  \code{\link{cvtunerange}}
#'
#' @references
#' Sandra E. Safo, Jeongyoun Ahn, Yongho Jeon, and Sungkyu Jung (2018),
#'  Sparse Generalized Eigenvalue Problem with Application to Canonical
#'  Correlation Analysis for Integrative Analysis of Methylation and Gene
#'  Expression Data. Biometrics
#' @export
#' @examples
#' ##---- read in data
#' data(selpData)
#'
#' Xdata1=selpData[[1]]
#' Xdata2=selpData[[2]]
#'
#'   ##---- estimate first canonical correlation vectors
#' ncancorr=1
#'
#' #use cvtunerange for range of tuning parameters
#' mytunerange=cvtunerange(Xdata1=Xdata1,Xdata2=Xdata2,ncancorr=ncancorr,
#'                         CovStructure="Iden",standardize=TRUE)
#' print(mytunerange)
#'
#' #Fix Tau for first and second datasets as 1.1 and 1.0 respectively
#' Tau=matrix(c(1,1.2,1),nrow=1)
#' mysparsevectors=multiplescca(Xdata1=Xdata1,Xdata2=Xdata2,ncancorr=ncancorr,
#'                              Tau=Tau, CovStructure="Iden",standardize=TRUE,
#'                              maxiteration=20, thresh=0.0001)
#'
#' #example with two canonical correlation vectors
#' #use cvselpscca to obtain optimal tuning parameters
#' mycv=cvselpscca(Xdata1=Xdata1,Xdata2=Xdata2,ncancorr=ncancorr,
#'                 CovStructure="Iden",isParallel=TRUE,ncores=NULL,nfolds=5,
#'                 ngrid=10, standardize=TRUE,thresh=0.0001,maxiteration=20)
#'
#' Tau=mycv$optTau
#' mysparsevectors=multiplescca(Xdata1=Xdata1,Xdata2=Xdata2,ncancorr=ncancorr,
#'                              Tau=Tau, CovStructure="Iden",standardize=TRUE,maxiteration=20,
#'                              thresh=0.0001)
multiplescca=function(Xdata1=Xdata1,Xdata2=Xdata2,ncancorr=ncancorr,Tau=Tau,CovStructure="Iden",standardize=TRUE,maxiteration=20,thresh=1e-04){

  Xdata1Orig=Xdata1
  Xdata2Orig=Xdata2

  if(dim(Xdata1)[1]!=dim(Xdata2)[1]) stop("Xdata1 and Xdata2 have different number of observations")

  if(is.null(ncancorr)){
    ncancorr=1
  }

  if(is.null(CovStructure)){
    CovStructure="Iden"
  }

  if(is.null(standardize)){
    standardize=TRUE
  }

  if(is.null(Tau)){
    stop('Tuning parameter (s) missing. You must provide one' )
  }

  if(standardize==TRUE){
    Xdata1=scale(Xdata1)
    Xdata2=scale(Xdata2)
  }

  if(is.null(maxiteration)){
    maxiteration=20
  }

  if(is.null(thresh)){
    thresh=1e-04
  }

  myfastcca=myfastnonsparsecca(Xdata1,Xdata2,CovStructure)
  tildealpha=myfastcca$tildeA
  tildebeta=myfastcca$tildeB
  tilderho=myfastcca$tilderho
  Ux=myfastcca$Ux
  Sigma12r=myfastcca$Sigma12r
  Uy=myfastcca$Uy

  if(max(Tau[,1])==1){
    maxiteration=maxiteration
    Tau=do.call(rbind, replicate(maxiteration, Tau, simplify=FALSE))
    Tau[,1]=matrix(rep(1:maxiteration,each=ncancorr),nrow= 1)
  }else{
    maxiteration=max(Tau[,1])
  }


  mybeta=as.matrix(tildebeta[,1:ncancorr])
  myalpha=as.matrix(tildealpha[,1:ncancorr])
  myrho=tilderho[1:ncancorr]
  iter=0
  diffalpha=1
  diffbeta=1
  reldiff=1
  while((iter<maxiteration) && min(reldiff, max(diffalpha,diffbeta))>thresh){
    iter=iter+1
    print(paste("Current Iteration Is:", iter))
    mybetaold=mybeta
    myalphaold=myalpha
    tilderhoold=myrho

    myselps=selpscca(Xdata1,Xdata2,mybetaold,myalphaold,tilderhoold,ncancorr,Ux,Sigma12r,Uy,as.matrix(Tau[Tau[,1]==iter,2]),as.matrix(Tau[Tau[,1]==iter,3]),CovStructure)
    myalpha=myselps$myhatalpha
    mybeta=myselps$myhatbeta
    myrho=myselps$mycorrmat

    if( min(colSums(abs(myalpha)))==0 || min(colSums(abs(mybeta)))==0){
      myalpha=myalphaold
      mybeta=mybetaold
      break
    }

    diffalpha=norm(myalpha-myalphaold,'f')^2 / norm(myalphaold,'f')^2
    diffbeta=norm(mybeta-mybetaold,'f')^2 / norm(mybetaold,'f')^2
    sumnew=norm(myalpha-myalphaold,'f')^2 + norm(mybeta-mybetaold,'f')^2
    sumold=norm(myalphaold,'f')^2+norm(mybetaold,'f')^2
    reldiff=sumnew/sumold
  }

  maxcorr=round(diag(abs(stats::cor(Xdata1%*%myalpha,Xdata2%*%mybeta))),3)
  hatalpha=myalpha
  hatbeta=mybeta

  NonZeroAlpha=colSums(hatalpha!=0)
  NonZeroBeta=colSums(hatbeta!=0)
  print(paste(c("Number of non-zero Xdata1 or View 1:", NonZeroAlpha),collapse=" "))
  print(paste(c("Number of non-zero Xdata2 or View 2:", NonZeroBeta),collapse=" "))
  print(paste(c("Corr(Xdata1*alpha,Xdata2*beta):", t(maxcorr)),collapse=" "))
  print(paste("Sparse CCA CovStructure used is:", CovStructure))

  result=list(estCorr=maxcorr,myalpha=hatalpha,mybeta=hatbeta,method="multiplescca", InputData=list(Xdata1Orig,Xdata2Orig))
  class(result)="SELPCCA"

  return(result)
}
