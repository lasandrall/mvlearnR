

IDANonsparse=function(Xdata, Y,weight){

  #   %--------------------------------------------------------------------------
  #   %IDANonsparse.R: function to obtain solution to integrative lda problem
  # without sparsity. Use this function if you want to perform joint integration
  # and separation without sparse constraints
  # %--------------------------------------------------------------------------
  #
  #

  #check data
  if (is.list(Xdata)) {
    D = length(Xdata)
  } else {
    stop("Input data should be a list")
  }


  #define weights
  w1=weight;
  w2=2*(1-weight)/(D*(D-1))

  Y=as.vector(Y)

  #nc=max(unique(as.vector(Y)))
  nc=length(unique(as.vector(Y)))

  Crxd=list()
  Sbx=list()
  myalphaold1=list()
  myalphaold2=list()
  myalphaoldmat=list()
  rmyalphaoldmat=list()
  sqrtminvmat=list()
  tildealphamat=list()
  tildelambda=list()

  for (d in 1:D){
    myX=as.matrix(Xdata[[d]])
    n=dim(myX)[1]
    p=dim(myX)[2]

    myX=t(myX)

    #if n is less than p, project into n space
    if(n<p){
      mysvd=svd(myX);
      Ux1=mysvd$u;
      V=mysvd$v;
      W=diag(mysvd$d)
      R=W%*%t(V)

      rdata2=cbind(Y,t(R))
      rdata=t(rdata2)
      mrd=aggregate(rdata2[,-1],list(rdata2[,1]),mean)
      mr=rowMeans(rdata[-1,])

      #nc=max(unique(Y))
      #nc=length(unique(as.vector(Y)))
      C=list()
      for(i in 1:nc)
      {
        C[[i]]=rdata2[rdata2[,1]==i,-1] - matrix(rep(t(mrd[mrd[,1]==i,-1]),times=sum(Y==i)) ,ncol=ncol(rdata2[,-1]),byrow=TRUE)
      }
      C=as.matrix(do.call(rbind,C))
      Swrx=t(C)%*%C /(n-1)

      #Crx=R-t(matrix(rep(mr,times=n),ncol=ncol(R),byrow=TRUE))
      Crx=R-rowMeans(R)
      Srx=Crx%*%t(Crx)/(n-1)

      Sbrx=Srx-Swrx
      Sbrx=Sbrx + t(Sbrx)
      Sbx[[d]]=Sbrx
      lambda=sqrt(log(p)/n)
      Strx=Swrx + lambda*diag(n)
      Strx=Strx + t(Strx)
    }else{
      #work in p space
      R=t(myX)
      rdata2=cbind((as.data.frame(Y)[,1]),R)
      rdata=t(rdata2)
      mrd=aggregate(rdata2[,-1],list(rdata2[,1]),mean)
      mr=rowMeans(rdata[-1,])

      #nc=max(unique(Y))
      #nc=length(unique(as.vector(Y)))
      C=list()
      for(i in 1:nc)
      {
        C[[i]]=rdata2[rdata2[,1]==i,-1] - matrix(rep(t(mrd[mrd[,1]==i,-1]),times=sum(Y==i)) ,ncol=ncol(rdata2[,-1]),byrow=TRUE)
      }
      C=as.matrix(do.call(rbind,C))
      Swrx=t(C)%*%(C) /(n-1)

      #Crx=R-t(matrix(rep(mr,times=n),ncol=ncol(R),byrow=TRUE))
      Crx=(t(R)-matrix(rep(colMeans(R),n),nrow=dim(R)[2],ncol=dim(R)[1]))
      Srx=(Crx)%*%t(Crx)/(n-1)

      Sbrx=Srx-Swrx
      Sbrx=Sbrx + t(Sbrx)
      Sbx[[d]]=Sbrx
      lambda=sqrt(log(p)/n)
      Strx=Swrx + lambda*diag(p)
      Strx=Strx + t(Strx)

    }



    Crxd[[d]]=Crx;

    #Set mybetaold and myalphaold as LDA solutions

    sqrtminv= mysqrtminv(Strx)$sqrtminv;
    sqrtminvmat[[d]]=sqrtminv;


    #myeigen=Re(eigs(sqrtminv%*%Sbrx%*%sqrtminv,nc-1))
    #myeigen=eigen(sqrtminv%*%Sbrx%*%sqrtminv,symmetric=TRUE)
    myeigen=eigs_sym(sqrtminv%*%Sbrx%*%sqrtminv,nc-1,which="LM")
    if(n < p){
      myalphaold1[[1]]=Ux1%*%myeigen$vectors
    }else{
      myalphaold1[[1]]=myeigen$vectors
    }
    myalphaoldmat[[d]]=do.call(rbind,lapply(myalphaold1, function(x) x/norm(x,'2')))
    rmyalphaoldmat[[d]]=myeigen$vectors
  }





  #nonsparse solution to integrative LDA



  for(d in 1:D){
    dd=setdiff(seq(1, D, by= 1),d)
    #cross-covariance
    rSumassociation=0;
    for (jj in 1:length(dd)){
      j=dd[jj];
      myalphaold=rmyalphaoldmat[[j]];
      Sdj=Crxd[[d]]%*%t(Crxd[[j]])/(n-1);
      rassociation=Sdj%*%myalphaold%*%t(myalphaold)%*%t(Sdj)
      rSumassociation=rSumassociation + rassociation + t(rassociation);
    }
    #solution to integrative LDA
    myinteig=eigs_sym(sqrtminvmat[[d]]%*%( w1*Sbx[[d]] +  w2*rSumassociation)%*%sqrtminvmat[[d]],nc-1,which="LM")
    myalphaold2[[1]]=myinteig$vectors
    tildealphamat[[d]]=do.call(rbind,lapply(myalphaold2, function(x) x/norm(x,'2')))
    tildelambda[[d]]=myinteig$values
  }

  result=list(tildealphamat=tildealphamat,tildelambda=tildelambda, myalphaoldmat=myalphaoldmat,sqrtminvmat=sqrtminvmat);
  return(result)
}
