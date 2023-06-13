
#' @title Correlation Plots
#'
#' @description Plots for visualizing correlation between estimated
#' discriminant vectors for pairwise data.
#'
#' @param Xtestdata A list with each entry containing views of size
#' \eqn{ntest \times p_d}, where \eqn{d = 1, \dots, D}.Rows are samples
#'  and columns are variables. Can use testing or training data
#' @param Ytest \eqn{ntest \times 1} vector of class membership.
#' @param hatalpha A list of estimated sparse discriminant vectors for each view.
#' @param color.palette  character vector of length K (number of classes), specifying the colors to use for the classes, respectively.
#' Defaults to shades of blue and orange (color.BlueOrange). Other option includes red and green combinations (color.GreenRed)
#'
#' @details The function will return correlation plot(s).
#'
#' @return
#' \item{NULL}{}
#'
#' @seealso \code{\link{cvSIDA}} \code{\link{sidatunerange}}
#' \code{\link{DiscriminantPlots}}
#'
#' @references
#' Sandra E. Safo, Eun Jeong Min, and Lillian Haine (2022), Sparse Linear Discriminant
#' Analysis for Multi-view Structured Data, Biometrics.
#' @import ggplot2
#' @export
#' @examples
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
#'
#' ##----plot discriminant and correlation plots
#' #---------Correlation plot
#' mycorrplot=CorrelationPlots(Xtestdata,Ytest,mysida$hatalpha)
#'
CorrelationPlots=function(Xtestdata=Xtestdata,Ytest=Ytest,hatalpha=hatalpha,color.palette=NULL){

  if(is.null(color.palette)){
    color.palette=color.BlueOrange(length(unique(Ytest)))
  }
  dsizes=lapply(Xtestdata, function(x) dim(x))
  D=length(dsizes)
  mycomb=t(utils::combn(D, 2, FUN = NULL, simplify = TRUE))

  # nc=unique(Ytest)
  # nnc=length(unique(Ytest))
  # mycol=mat.or.vec(length(Ytest),1)
  # mypch=mat.or.vec(length(Ytest),1)
  # for(j in 1:length(mycol)){
  #   mycol[j]=nc[Ytest[j]]
  #   mypch[j]=nc[Ytest[j]]
  # }

  for(d in 1:dim(mycomb)[1]){
    dd=mycomb[d,]
    Scoresd=Xtestdata[[dd[1]]]%*%hatalpha[[dd[1]]]
    Scoresj=Xtestdata[[dd[2]]]%*%hatalpha[[dd[2]]]
    Scores=cbind.data.frame(Scoresd[,1], Scoresj[,1], Ytest)
    colnames(Scores)=c("Disc1", "Disc2", "Class")

    #calculate RV coefficient
    X1=Xtestdata[[dd[1]]]%*%hatalpha[[dd[1]]]
    X2=Xtestdata[[dd[2]]]%*%hatalpha[[dd[2]]]
    X1=scale(X1, center=TRUE,scale=FALSE)
    X2=scale(X2, center=TRUE,scale=FALSE)
    X1X2=t(X1)%*%X2/dim(X1)[1]
    X1X1=t(X1)%*%X1/dim(X1)[1]
    X2X2=t(X2)%*%X2/dim(X2)[1]
    RVCoeff=sum(diag(X1X2%*%t(X1X2)))/(sqrt(sum(diag(X1X1%*%X1X1)))*sqrt(sum(diag(X2X2%*%X2X2))))
    RVCoeff=round(RVCoeff,digits=2)


    # plot(Scores[,1], Scores[,2],col=mycol,lwd=3,pch=mypch
    #      ,xlab=paste(
    #        "First Discriminant Score for View", dd[1]),ylab=paste("First Discriminant Score for View", dd[2]),xaxt="n",
    #      yaxt="n", main=paste("Correlation plot for views",dd[1], "and" ,dd[2],  "\u03C1 =", RVCoeff))
    # graphics::legend("topleft",bty = "n",cex=0.8,legend=c(nc) ,col = 1:max(nc), pch=1:max(nc),title="Class")

    Classes=factor(Ytest)
    print(
    ggplot2::ggplot(Scores, aes(Scores[,1], Scores[,2])) +
      geom_point(aes(shape=Classes,colour = Classes),size=4) +
      xlab(paste(
        "First Discriminant Score for View", dd[1])) +
      ylab(paste("First Discriminant Score for View", dd[2])) +
      ggtitle(paste("Correlation plot for views",dd[1], "and" ,dd[2],",", "\u03C1 =", RVCoeff)) +
      scale_colour_manual(values=color.palette) +
      scale_fill_manual(values = color.palette) +
      theme_bw() +
      theme(axis.title = element_text(face="bold"))+
      theme(axis.text = element_blank())+
      theme(axis.ticks.x = element_blank())+
      theme(axis.ticks.y=element_blank())+
      theme(axis.ticks = element_blank())
    )



  }
  return(invisible(NULL))


}



#' @title Discriminant Plots
#'
#' @description Plots discriminant scores for visualizing class separation
#'
#' @param Xtestdata A list with each entry containing views of size
#' \eqn{ntest \times p_d}, where \eqn{d = 1, \dots, D}. Rows are samples and
#' columns are variables. Can use testing or training data.
#' @param Ytest \eqn{ntest \times 1} vector of class membership.
#' @param hatalpha A list of estimated sparse discriminant vectors for each view.
#' @param color.palette  character vector of length K (number of classes), specifying the colors to use for the classes, respectively.
#' Defaults to shades of blue and orange (color.BlueOrange). Other option includes red and green combinations (color.GreenRed)
#' @details The function will return discriminant plots.
#'
#' @return
#' \item{NULL}{}
#'
#' @seealso \code{\link{cvSIDA}} \code{\link{sidatunerange}}
#' \code{\link{CorrelationPlots}}
#'
#' @references
#' Sandra E. Safo, Eun Jeong Min, and Lillian Haine (2023), Sparse Linear Discriminant
#' Analysis for Multi-view Structured Data, Biometrics.
#'
#' @import ggplot2
#' @export
#' @examples
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
#'
#' ##----plot discriminant plots
#' #---------Discriminant plot
#' mydisplot=DiscriminantPlots(Xtestdata,Ytest,mysida$hatalpha)
#'
DiscriminantPlots=function(Xtestdata=Xtestdata,Ytest=Ytest,hatalpha=hatalpha,color.palette=NULL){



  dsizes=lapply(Xtestdata, function(x) dim(x))
  D=length(dsizes)

  nc=unique(Ytest)
  if(length(unique(Ytest))==2){
    if(is.null(color.palette)){
    color.palette=color.BlueOrange(2)
    }
    #graphics::par(mfrow=c(1,D))
    for(d in 1:D){
      myScores=Xtestdata[[d]]%*%hatalpha[[d]]
      ss21=myScores[Ytest==nc[1],]
      ss22=myScores[Ytest==nc[2],]
      dd1=stats::density(ss21);
      dd2=stats::density(ss22);
      plot(dd1,xlim=c(min(myScores[,1]-0.7),max(myScores[,1]+0.7)),col=color.palette[1],cex.axis=1.5,
           cex.lab=1.5,lwd=2.5,xlab="First Discriminant Score ",main=paste("Discriminant Plot for View",d),ylim=c(0,max(dd1$y,dd2$y)))
      graphics::lines(dd2,xlim=c(min(myScores[,1]-0.7),max(myScores[,1]+0.7)),col=color.palette[2],cex.axis=1,cex.lab=1,lwd=2.5)
      graphics::points(cbind(ss21,rep(0,length(ss21))),col=color.palette[1],pch=3,cex=4)
      graphics::points(cbind(ss22,rep(0,length(ss22))),col=color.palette[2],pch=1,cex=4)
      #legend("topleft", legend = c(paste("Class", nc[1]),paste("Class", nc[2])), col = c("red","black"),cex=1.2,lty=1,lwd=2.5)
      graphics::legend("topleft", inset=c(0,0),legend=c(nc), col = color.palette,pch=c(3,1),title="Class")

    }
    #graphics::par(mfrow=c(1,1))
  }else if(length(unique(Ytest))>2){
    # graphics::par(mfrow=c(1,D))
    # nnc=length(unique(Ytest))
    # mycol=mat.or.vec(length(Ytest),1)
    # mypch=mat.or.vec(length(Ytest),1)
    # for(j in 1:length(mycol)){
    #   mycol[j]=nc[Ytest[j]]
    #   mypch[j]=nc[Ytest[j]]
    #
    # }
    if(is.null(color.palette)){
      color.palette=color.BlueOrange(length(unique(Ytest)))
    }
    Classes=factor(Ytest)
    for(d in 1:D){
      #Scores=cbind.data.frame(Ytest,Xtestdata[[d]]%*%cbind(hatalpha[[d]],hatalpha[[d]]))
      Scores=cbind.data.frame(Ytest,Xtestdata[[d]]%*%hatalpha[[d]])
      # plot(Scores[,2], Scores[,3],col=mycol,lwd=2.5,pch=mypch,
      #      xlab=paste(
      #        "First Discriminant Score for View", d),
      #      ylab=paste("Second Discriminant Score for View", d),
      #      xaxt="n",
      #      yaxt="n",
      #      main=paste("Discriminant Plot for View",d),
      #      cex.lab=1.5,cex.axis=1.5,cex.main=1.5,cex.sub=1.5)
      # graphics::par(xpd=TRUE)
      # graphics::legend("topright",bty = "n",legend=c(nc),col = 1:max(nc), pch=1:max(nc),title="Class")




      #dev.new()
      print(
        ggplot2::ggplot(Scores[,-1], aes(Scores[,2], Scores[,3])) +
          geom_point(aes(shape=Classes,colour = Classes),size=4) +
          xlab(paste(
            "First Discriminant Score for View", d)) +
          ylab(paste("Second Discriminant Score for View", d)) +
          ggtitle(paste("Discriminant Plot for View",d)) +
          scale_colour_manual(values=color.palette) +
          scale_fill_manual(values = color.palette) +
          theme_bw() +
          stat_ellipse(aes(Scores[,2], Scores[,3], color=Classes, fill = Classes),type = "norm", geom = "polygon", alpha = 0.1)+
          guides(colour = guide_legend(override.aes = list(size=5))) +
          scale_size_continuous(range=c(10,15))+
          theme_bw() +
          theme(axis.title = element_text(face="bold"))+
          theme(axis.text = element_blank())+
          theme(axis.ticks.x = element_blank())+
          theme(axis.ticks.y=element_blank())+
          theme(axis.ticks = element_blank())
        )
      Sys.sleep(2)

    }
    #graphics::par(mfrow=c(1,1))
  }


  return(invisible(NULL))

}


#' @title Loadings Plots
#'
#' @description Plots discriminant and canonical vectors  to visualize how
#' selected variables contribute to the first and second discriminant (for SIDA and SIDANet)
#' or canonical correlation (for SELPCCA) vectors.  Variables farther from the origin and close to first or second axis
#' have higher impact on first or second discriminant/canonical vectors, respectively.
#' Variables farther from the origin and between both first and second axes have similar higher contributions to the
#' first and second discriminant/canonical correlation vectors. In both situations, for SIDA and SIDANet, this suggests that
#' these variables contribute more to the separation of classes and association of views. For SELPCCA, this suggests that
#' these variables contribute more to the association between the two views. This plot can
#' only be generated for classification and association problems with 3 or more classes (SIDA and SIDANet),
#' or for CCA problems with two or more canonical correlation vectors requested (i.e. ncancorr > 1 for SELPCCA).
#'
#' @param object the output from SIDA, SIDANet, and SELPCCA methods
#' @param color.line color to use for plotting direction vectors. Default is "darkgray".
#' @param keep.loadings numeric, specifying how many variables to represent on loadings plot. This is useful
#' in situations where the number of variables selected is large, and could clutter the plot. If this number is more
#' than the variables selected, it will be set to the maximum number of variables selected for each view.
#' Default is plotting all selected variables.
#' @details The function will return loading plots, one for each view.
#'
#' @return
#' \item{NULL}{}
#'
#' @seealso \code{\link{cvSIDA}} \code{\link{DiscriminantPlots}}
#' \code{\link{CorrelationPlots}}
#'
#' @references
#' Sandra E. Safo, Eun Jeong Min, and Lillian Haine (2023), Sparse Linear Discriminant
#' Analysis for Multi-view Structured Data, Biometrics.
#' Sandra E. Safo, Jeongyoun Ahn, Yongho Jeon, and Sungkyu Jung (2018),
#'  Sparse Generalized Eigenvalue Problem with Application to Canonical
#'  Correlation Analysis for Integrative Analysis of Methylation and Gene
#'  Expression Data. Biometrics
#'
#' @import ggplot2
#' @import ggthemes
#' @export
#' @examples
#'data("sidanetData")
#'Xdata <- sidanetData[[1]]
#'Y <- sidanetData[[2]] #class membership already coded as 1,2,...
#'Xtestdata <- sidanetData[[3]]
#'Ytest <- sidanetData[[4]] #class membership already coded as 1,2,...
#'#edge information
#'myedges=sidanetData[[5]]
#'myedgeweight=sidanetData[[6]]
#'##---- call cross validation
#'mycvsidanet=cvSIDANet(Xdata,Y,myedges,myedgeweight,withCov=FALSE,plotIt=FALSE,Xtestdata=Xtestdata,
#'                      Ytest=Ytest)
#'LoadingsPlots(mycvsidanet,keep.loadings=c(3,3))

LoadingsPlots=function(object,color.line="darkgray",keep.loadings=NULL){
  if(class(object)=="SIDA" | class(object)=="SIDANet"){
    hatalpha=object$hatalpha
  }else if( class(object)=="SELPCCA"){

    if(object$method=="selpscca.pred"){
      L=dim(hatalpha[[1]])[2]
      hatalpha=list(object$selp.fit$hatalpha,object$selp.fit$hatbeta)
      for(j in 1:L){
        hatalpha[[j]]=qr.Q(qr(hatalpha[[j]]))
      }
    }else{ hatalpha=list(object$hatalpha,object$hatbeta)
    L=dim(hatalpha[[1]])[2]
    for(j in 1:L){
      hatalpha[[j]]=qr.Q(qr(hatalpha[[j]]))
    }}

  }

  if(is.null(color.line)){
    color.line="darkgray"
  }

  D <- length(hatalpha)

  L=dim(hatalpha[[1]])[2]


  ncomp=dim(hatalpha[[1]])[2]
  for(jj in 1:D){

    X1=as.data.frame(object$InputData[[jj]])
    if(ncomp == 1){
      if(class(object)=="SIDA" | class(object)=="SIDANet"){
        stop("Loadings plot not applicable with one discriminant vector" , call. = FALSE)
      }else if(class(object)=="SELPCCA"){
        stop("Loadings plot not applicable with one CCA vector " , call. = FALSE)
      }
    }else if (ncomp >1){
      hatalpha1=rowSums(abs(hatalpha[[jj]]))
    }
      hatalpha2=hatalpha1[order(hatalpha1,decreasing=TRUE)]
      col1=order(hatalpha1,decreasing=TRUE)

      X1var.Ind=which(as.matrix(hatalpha1)!=0, arr.ind = TRUE)

      if(!is.null(keep.loadings)){

        if(keep.loadings[[jj]]>sum(hatalpha1!=0)){
          warning("keep.loadings is greater than maximum number of variables selected, setting to this maximum")
          keep.loadings[[jj]]=sum(hatalpha1!=0)
        }
        # X1var.Ind1=which(as.data.frame(hatalpha1)!=0, arr.ind = TRUE)[,1]
        #
        #
        # X1var.Ind=colnames(X1[,X1var.Ind1])
        # which(as.data.frame(hatalpha1[col1[1:keep.loadings[[jj]]]])!=0, arr.ind = TRUE)
        myloadings=as.data.frame(scale(hatalpha[[jj]][col1[1:keep.loadings[jj]],],center=FALSE, scale=FALSE))

        #X1var.ind[,1][col1[1:keep.loadings[[jj]]]]
        X1var.ind=colnames(X1)[col1[1:keep.loadings[[jj]]]]
        #X1var.Ind1=(col1[X1var.Ind[,1]])[1:keep.loadings[[jj]]]
        var.names=sub("\\;.*", "", X1var.ind)

      }else if(is.null(keep.loadings)){
        X1var.Ind=which(as.matrix(hatalpha1)!=0, arr.ind = TRUE)
        #X1var.Ind=which(as.matrix(hatalpha1)!=0, arr.ind = TRUE)
        X1var.Ind1=X1var.Ind[,1]
        var.names=sub("\\;.*", "", colnames(X1[,X1var.Ind1]))
        myloadings=as.data.frame(scale(hatalpha[[jj]][X1var.Ind1,],center=FALSE, scale=FALSE))
      }

    # X1var.Ind=X1var.Ind[,1]
    # var.names=sub("\\;.*", "", colnames(X1[,X1var.Ind]))

     #myloadings=as.data.frame(scale(hatalpha[[jj]][col1[1:keep.loadings[jj]],],center=FALSE, scale=FALSE))
     if(keep.loadings[[jj]]==1){
       myloadings=as.data.frame(t(scale(myloadings,center=FALSE,scale=FALSE)))
     }

  print(
    ggplot2::ggplot(myloadings, aes(myloadings[,1], myloadings[,2])) +
      geom_text(aes(label = var.names), size = 4) +
      geom_segment(data = myloadings[,1:2], aes(x=0, y=0, xend=myloadings[,1], yend=myloadings[,2]),
                   linewidth=1,
                   arrow=arrow(length=unit(0.3, "cm")), color = color.line) +
      xlab(
        if(class(object)=="SIDA" | class(object)=="SIDANet"){
          paste("First Discriminant Vector for View", jj)
        }else if(class(object)=="SELPCCA"){
          paste("First Canonical Vector for View", jj)
        }
      ) +
      ylab(if(class(object)=="SIDA" | class(object)=="SIDANet"){
        paste("Second Discriminant Vector for View", jj)
      }else if(class(object)=="SELPCCA"){
        paste("Second Canonical Vector for View", jj)
      }
      ) +
      ggtitle("Loading Plot for View", jj) + theme_bw() +
      theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      geom_vline(xintercept = 0, color = "darkgrey") +
      geom_hline(yintercept = 0, color = "darkgrey")+
      ggthemes::theme_stata(scheme="s2manual")+
      xlim(-1,1)
  )
  Sys.sleep(2)

  }

#  return(invisible(list(SimilirarityMatrix=Mjk, ViewCombinations=mycomb)))
}




#' @title Biplots for Discriminant Scores or Canonical Correlation Variates for each View
#'
#' @description Biplots  to visualize discriminant scores/ canonical variates and how
#' selected variables contribute to the first and second discriminant (for SIDA and SIDANet)
#' or canonical correlation (for SELPCCA) vectors.  Variables farther from the origin and close to first or second axis
#' have higher impact on first or second discriminant/canonical vectors, respectively.
#' Variables farther from the origin and between both first and second axes have similar higher contributions to the
#' first and second discriminant/canonical correlation vectors. In both situations, for SIDA and SIDANet, this suggests that
#' these variables contribute more to the separation of classes and association of views. For SELPCCA, this suggests that
#' these variables contribute more to the association between the two views. This plot can
#' only be generated for classification and association problems with 3 or more classes (SIDA and SIDANet),
#' or for CCA problems with two or more canonical correlation vectors requested (i.e. ncancorr > 1 for SELPCCA).
#'
#' @param object the output from SIDA, SIDANet, and SELPCCA methods
#' @param Y a vector of class membership for grouping canonical correlatoin variates and discriminant scores.
#' @param color.palette  character vector of length K (number of classes), specifying the colors to use for the classes, respectively.
#' Defaults to shades of blue and orange (color.BlueOrange). Other option includes red and green combinations (color.GreenRed)
#' @param keep.loadings numeric vector of length D (number of views), specifying how many variables
#' to represent on loadings plot for each view. This is useful
#' in situations where the number of variables selected is large, and could clutter the plot. If this number is more
#' than the variables selected, it will be set to the maximum number of variables selected for each view.
#' Default is plotting all selected variables.
#' @details The function will return loading plots, one for each view.
#'
#' @return
#' \item{NULL}{}
#'
#' @seealso \code{\link{cvSIDA}} \code{\link{DiscriminantPlots}}
#' \code{\link{CorrelationPlots}}
#'
#' @references
#' Sandra E. Safo, Eun Jeong Min, and Lillian Haine (2023), Sparse Linear Discriminant
#' Analysis for Multi-view Structured Data, Biometrics.
#' Sandra E. Safo, Jeongyoun Ahn, Yongho Jeon, and Sungkyu Jung (2018),
#'  Sparse Generalized Eigenvalue Problem with Application to Canonical
#'  Correlation Analysis for Integrative Analysis of Methylation and Gene
#'  Expression Data. Biometrics
#'
#' @import ggplot2
#' @import ggthemes
#' @export
#' @examples
#'data("sidanetData")
#'Xdata <- sidanetData[[1]]
#'Y <- sidanetData[[2]] #class membership already coded as 1,2,...
#'Xtestdata <- sidanetData[[3]]
#'Ytest <- sidanetData[[4]] #class membership already coded as 1,2,...
#'#edge information
#'myedges=sidanetData[[5]]
#'myedgeweight=sidanetData[[6]]
#'##---- call cross validation
#'mycvsidanet=cvSIDANet(Xdata,Y,myedges,myedgeweight,withCov=FALSE,plotIt=FALSE,Xtestdata=Xtestdata,
#'                      Ytest=Ytest)
#'WithinViewBiplot(mycvsidanet,Y, color.palette=NULL,keep.loadings=c(3,3))

WithinViewBiplot=function(object,Y, color.palette=NULL,keep.loadings=NULL){
  if(class(object)=="SIDA" | class(object)=="SIDANet"){
    hatalpha=object$hatalpha
  }else if( class(object)=="SELPCCA"){
    if(object$method=="selpscca.pred"){
      hatalpha=list(object$selp.fit$hatalpha,object$selp.fit$hatbeta)
      L=dim(hatalpha[[1]])[2]
      for(j in 1:L){
        hatalpha[[j]]=qr.Q(qr(hatalpha[[j]]))
      }
    }else{ hatalpha=list(object$hatalpha,object$hatbeta)
    L=dim(hatalpha[[1]])[2]
    for(j in 1:L){
      hatalpha[[j]]=qr.Q(qr(hatalpha[[j]]))
    }}

  }

  if(is.null(color.palette)){
    color.palette=color.BlueOrange(length(unique(Y)))
  }

  D <- length(hatalpha)

  L=dim(hatalpha[[1]])[2]



  # for(j in 1:L){
  #   hatalpha[[j]]=qr.Q(qr(hatalpha[[j]]))
  # }


  Classes=factor(Y)

  ncomp=dim(hatalpha[[1]])[2]
  for(jj in 1:D){

    X1=as.data.frame(object$InputData[[jj]])
    if(ncomp == 1){
      if(class(object)=="SIDA" | class(object)=="SIDANet"){
        stop("Loadings plot not applicable with one discriminant vector" , call. = FALSE)
      }else if(class(object)=="SELPCCA"){
        stop("Loadings plot not applicable with one CCA vecto " , call. = FALSE)
      }
    }else if (ncomp >1){
      hatalpha1=rowSums(abs(hatalpha[[jj]]))
    }

    hatalpha2=hatalpha1[order(hatalpha1,decreasing=TRUE)]
    col1=order(hatalpha1,decreasing=TRUE)

    X1var.Ind=which(as.matrix(hatalpha1)!=0, arr.ind = TRUE)

   if(!is.null(keep.loadings)){
    if(keep.loadings[[jj]]>sum(hatalpha1!=0)){
      warning("keep.loadings is greater than maximum number of variables selected, setting to this maximum")
      keep.loadings[[jj]]=sum(hatalpha1!=0)
    }

    # X1var.Ind1=(col1[X1var.Ind[,1]])[1:keep.loadings[[jj]]]
    # var.names=sub("\\;.*", "", colnames(X1)[X1var.Ind1])

    myloadings=as.data.frame(scale(hatalpha[[jj]][col1[1:keep.loadings[jj]],],center=FALSE, scale=FALSE))

    #X1var.ind[,1][col1[1:keep.loadings[[jj]]]]
    X1var.ind=colnames(X1)[col1[1:keep.loadings[[jj]]]]
    #X1var.Ind1=(col1[X1var.Ind[,1]])[1:keep.loadings[[jj]]]
    var.names=sub("\\;.*", "", X1var.ind)



   }else if(is.null(keep.loadings)){
     X1var.Ind=which(as.matrix(hatalpha1)!=0, arr.ind = TRUE)
     X1var.Ind1=X1var.Ind[,1]
     var.names=sub("\\;.*", "", colnames(X1)[X1var.Ind1])
     myloadings=as.data.frame(scale(hatalpha[[jj]][X1var.Ind1,],center=FALSE, scale=FALSE))

   }


    #X1var.Ind=X1var.Ind[,1]
    #var.names=sub("\\;.*", "", colnames(X1[,X1var.Ind]))

    #Scores=cbind.data.frame(Y,as.matrix(X1[,X1var.Ind])%*%hatalpha[[jj]][X1var.Ind,]) #selected variables
    Scores=cbind.data.frame(Y,as.matrix(X1)%*%hatalpha[[jj]])
    if(keep.loadings[[jj]]==1){
      myloadings=as.data.frame(t(myloadings))
    }else{
      myloadings=as.data.frame(myloadings)
    }

    mxmax=max(abs(Scores[,2]))+20
    mymax=max(abs(Scores[,3]))+20

    print(
      ggplot2::ggplot() +
        geom_point(data = Scores[,-1], aes(shape=Classes,x = Scores[,2], y = Scores[,3], colour = Classes),size=4) +
        stat_ellipse(aes(Scores[,2], Scores[,3], color=Classes, fill = Classes),type = "norm", geom = "polygon", alpha = 0.1)+
        geom_segment(data = myloadings[,1:2], aes(x=0, y=0, xend=mxmax*myloadings[,1], yend=mymax*myloadings[,2]),
                     linewidth=1,
                     arrow=arrow(length=unit(0.3, "cm")), color = "black") +
        scale_y_continuous(sec.axis = sec_axis(~./26)) +
        geom_text(data = myloadings, aes(x=mxmax*myloadings[,1], y=mymax*myloadings[,2], label = var.names),
                  size = 3.5, hjust = "outward", vjust = "outward") +
        scale_colour_manual(values=color.palette) +
        scale_fill_manual(values = color.palette) +
        xlab(paste(
          "First Discriminant Score for View ", jj)) +
        ylab(paste("Second Discriminant Score for View ", jj)) +
        ggtitle(
          if(class(object)=="SIDA" | class(object)=="SIDANet"){
            paste("SIDA Biplot for View ",jj)
          }else if(class(object)=="SELPCCA"){
            paste("SELPCCA Biplot for View ",jj)
          }
        )+
        theme_bw() +
        theme(axis.line = element_line(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        geom_vline(xintercept = 0, color = "darkgrey") +
        geom_hline(yintercept = 0, color = "darkgrey") +
        ggthemes::theme_stata(scheme="s2manual")+
        xlim(min(Scores[,2])-20,max(Scores[,2])+20)
    )
    Sys.sleep(2)

  }

  #  return(invisible(list(SimilirarityMatrix=Mjk, ViewCombinations=mycomb)))
}


#' @title Biplots for Discriminant Scores or Canonical Correlation Variates between pairs of views
#'
#' @description Biplots  to visualize discriminant scores/ canonical variates  and how
#' selected variables contribute to the first and second discriminant (for SIDA and SIDANet)
#' or canonical correlation (for SELPCCA) vectors.  Variables farther from the origin and close to first or second axis
#' have higher impact on first or second discriminant/canonical vectors, respectively.
#' Variables farther from the origin and between both first and second axes have similar higher contributions to the
#' first and second discriminant/canonical correlation vectors. In both situations, for SIDA and SIDANet, this suggests that
#' these variables contribute more to the separation of classes and association of views. For SELPCCA, this suggests that
#' these variables contribute more to the association between the two views. This plot can
#' only be generated for classification and association problems with 3 or more classes (SIDA and SIDANet),
#' or for CCA problems with two or more canonical correlation vectors requested (i.e. ncancorr > 1 for SELPCCA).
#'
#' @param object the output from SIDA, SIDANet, and SELPCCA methods
#' @param Y a vector of class membership for grouping canonical correlatoin variates and discriminant scores.
#' @param color.palette  character vector of length K (number of classes), specifying the colors to use for the classes, respectively.
#' Defaults to shades of blue and orange (color.BlueOrange). Other option includes red and green combinations (color.GreenRed)
#' @param keep.loadings numeric vector of length D (number of views), specifying how many variables
#' to represent on loadings plot for each view. This is useful
#' in situations where the number of variables selected is large, and could clutter the plot. If this number is more
#' than the variables selected, it will be set to the maximum number of variables selected for each view.
#' Default is plotting all selected variables.
#' @details The function will return loading plots, one for each view.
#'
#' @return
#' \item{NULL}{}
#'
#' @seealso \code{\link{cvSIDA}} \code{\link{DiscriminantPlots}}
#' \code{\link{CorrelationPlots}}
#'
#' @references
#' Sandra E. Safo, Eun Jeong Min, and Lillian Haine (2023), Sparse Linear Discriminant
#' Analysis for Multi-view Structured Data, Biometrics.
#' Sandra E. Safo, Jeongyoun Ahn, Yongho Jeon, and Sungkyu Jung (2018),
#'  Sparse Generalized Eigenvalue Problem with Application to Canonical
#'  Correlation Analysis for Integrative Analysis of Methylation and Gene
#'  Expression Data. Biometrics
#'
#' @import ggplot2
#' @import ggthemes
#' @export
#' @examples
#'data("sidanetData")
#'Xdata <- sidanetData[[1]]
#'Y <- sidanetData[[2]] #class membership already coded as 1,2,...
#'Xtestdata <- sidanetData[[3]]
#'Ytest <- sidanetData[[4]] #class membership already coded as 1,2,...
#'#edge information
#'myedges=sidanetData[[5]]
#'myedgeweight=sidanetData[[6]]
#'##---- call cross validation
#'mycvsidanet=cvSIDANet(Xdata,Y,myedges,myedgeweight,withCov=FALSE,plotIt=FALSE,Xtestdata=Xtestdata,
#'                      Ytest=Ytest)
#'BetweenViewBiplot(mycvsidanet, Y,keep.loadings=c(3,3) )

BetweenViewBiplot=function(object,Y, color.palette=NULL,keep.loadings=c(20,30)){
  if(class(object)=="SIDA" | class(object)=="SIDANet"){
    hatalpha=object$hatalpha
  }else if( class(object)=="SELPCCA"){
    if(object$method=="selpscca.pred"){
      hatalpha=list(object$selp.fit$hatalpha,object$selp.fit$hatbeta)
      L=dim(hatalpha[[1]])[2]
      for(j in 1:L){
        hatalpha[[j]]=qr.Q(qr(hatalpha[[j]]))
      }
    }else{
      hatalpha=list(object$hatalpha,object$hatbeta)
      L=dim(hatalpha[[1]])[2]
    for(j in 1:L){
      hatalpha[[j]]=qr.Q(qr(hatalpha[[j]]))
    }}

  }
  if(is.null(color.palette)){
    color.palette=color.BlueOrange(length(unique(Y)))
  }

  D <- length(hatalpha)

  L=dim(hatalpha[[1]])[2]



  # for(j in 1:L){
  #   hatalpha[[j]]=qr.Q(qr(hatalpha[[j]]))
  # }

  mycomb=combn(D,2) #pairwise combination
  ncomp=dim(hatalpha[[1]])[2]
  myloadings=list()

  Classes=factor(Y)

  ncomp=dim(hatalpha[[1]])[2]
  for(jj in 1:dim(mycomb)[2]){

    X1=as.data.frame(object$InputData[[mycomb[1,jj]]])
    X2=as.data.frame(object$InputData[[mycomb[2,jj]]])

    if(ncomp == 1){
      if(class(object)=="SIDA" | class(object)=="SIDANet"){
        stop("Loadings plot not applicable with one discriminant vector" , call. = FALSE)
      }else if(class(object)=="SELPCCA"){
        stop("Loadings plot not applicable with one CCA vecto " , call. = FALSE)
      }
    }else if (ncomp >1){
      hatalpha1=rowSums(abs(hatalpha[[mycomb[1,jj]]]))
      col1=order(hatalpha1,decreasing=TRUE)
      hatalpha2=rowSums(abs(hatalpha[[mycomb[2,jj]]]))
      col2=order(hatalpha2,decreasing=TRUE)

    }


    if(!is.null(keep.loadings)){
      if(keep.loadings[[mycomb[1,jj]]]>sum(hatalpha1!=0)){
        warning("keep.loadings is greater than maximum number of variables selected, setting to this maximum")
        keep.loadings[[mycomb[1,jj]]]=sum(hatalpha1!=0)
      }

      if(keep.loadings[[mycomb[2,jj]]]>sum(hatalpha2!=0)){
        warning("keep.loadings is greater than maximum number of variables selected, setting to this maximum")
        keep.loadings[[mycomb[2,jj]]]=sum(hatalpha2!=0)
      }
      #for one view
      mycolnames=sub("\\;.*", "", colnames(object$InputData[[mycomb[1,jj]]]))
      var1.names <- mycolnames[col1[1:keep.loadings[[mycomb[1,jj]]]]]
      #for another view
      mycolnames=sub("\\;.*", "", colnames(object$InputData[[mycomb[2,jj]]]))
      var2.names <- mycolnames[col2[1:keep.loadings[[mycomb[2,jj]]]]]



      myloadings[[1]]=as.data.frame(scale(hatalpha[[mycomb[1,jj]]][col1[1:keep.loadings[[mycomb[1,jj]]]],],center=FALSE,scale=FALSE))
      myloadings[[2]]=as.data.frame(scale(hatalpha[[mycomb[2,jj]]][col2[1:keep.loadings[[mycomb[2,jj]]]],],center=FALSE,scale=FALSE))


    }else if(is.null(keep.loadings)){
      X1var.Ind=which(as.matrix(hatalpha1)!=0, arr.ind = TRUE)
      X1var.Ind=X1var.Ind[,1]
      X2var.Ind=which(as.matrix(hatalpha2)!=0, arr.ind = TRUE)
      X2var.Ind=X2var.Ind[,1]
      var1.names=sub("\\;.*", "", colnames(X1[,X1var.Ind]))
      var2.names=sub("\\;.*", "", colnames(X2[,X2var.Ind]))

      myloadings[[1]]=as.data.frame(scale(hatalpha[[mycomb[1,jj]]][X1var.Ind,],center=FALSE, scale=FALSE))
      myloadings[[2]]=as.data.frame(scale(hatalpha[[mycomb[2,jj]]][X2var.Ind,],center=FALSE,scale=FALSE))

    }

    if(length(intersect(var1.names,var2.names))!=0){
      var1.names <- paste(var1.names,mycomb[1,jj],sep="_")
      var2.names <- paste(var2.names,mycomb[2,jj],sep="_")
    }

    U=as.matrix(X1)%*%hatalpha[[mycomb[1,jj]]]
    V=as.matrix(X2)%*%hatalpha[[mycomb[2,jj]]]
    Z=U+V
    Scores=cbind.data.frame(Y,Z)

    if(keep.loadings[[mycomb[1,jj]]]==1){
      myloadings[[1]]=as.data.frame(t(myloadings[[1]]))
    }

    if(keep.loadings[[mycomb[2,jj]]]==1){
      myloadings[[2]]=as.data.frame(t(myloadings[[2]]))
    }

    mxmax=max(abs(Scores[,2]))+35
    mymax=max(abs(Scores[,3]))+35

    print(
      ggplot2::ggplot() +
        geom_point(data = Scores[,-1], aes(shape=Classes,x = Scores[,2], y = Scores[,3], colour = Classes),size=4) +
        stat_ellipse(aes(Scores[,2], Scores[,3], color=Classes, fill = Classes),type = "norm", geom = "polygon", alpha = 0.1)+
        geom_segment(data = myloadings[[1]][,1:2], aes(x=0, y=0, xend=mxmax*myloadings[[1]][,1], yend=mymax*myloadings[[1]][,2],                                                       ),
                     linewidth=1,
                     arrow=arrow(length=unit(0.3, "cm")), color = "black") +
        geom_text(data = myloadings[[1]], aes(x=mxmax*myloadings[[1]][,1], y=mymax*myloadings[[1]][,2], label = var1.names),
                  size = 4, hjust = "outward", vjust = "outward") +
        geom_segment(data = myloadings[[2]][,1:2], aes(x=0, y=0, xend=mxmax*myloadings[[2]][,1], yend=mymax*myloadings[[2]][,2],
                                                       ),
                     linewidth=2,linetype=2,
                     arrow=arrow(length=unit(0.3, "cm")), color = "red") +
        geom_text(data = myloadings[[2]], aes(x=mxmax*myloadings[[2]][,1], y=mymax*myloadings[[2]][,2],
                                              label = var2.names),
                  size = 4, hjust = "outward", vjust = "outward") +
        scale_y_continuous(sec.axis = sec_axis(~./26)) +
        scale_colour_manual(values=color.palette) +
        scale_fill_manual(values = color.palette) +
        # scale_linetype_manual("view1", values=c("view1"=2))+
        # scale_linetype_manual("view2", values=c("view1"=1))+
        xlab(paste(
          "First Discriminant Score")) +
        ylab(paste("Second Discriminant Score")) +
        ggtitle(
          if(class(object)=="SIDA" | class(object)=="SIDANet"){
            paste("SIDA Biplot for Views ",mycomb[1,jj], "and", mycomb[2,jj])
          }else if(class(object)=="SELPCCA"){
            paste("SELPCCA Biplot for Views ",mycomb[1,jj], "and", mycomb[2,jj])
          }
        )+
        theme_bw() +
        theme(axis.line = element_line(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        geom_vline(xintercept = 0, color = "darkgrey") +
        geom_hline(yintercept = 0, color = "darkgrey") +
        ggthemes::theme_stata(scheme="s2manual")+
        xlim(min(Scores[,2])-20,max(Scores[,2])+20)
    )
    Sys.sleep(2)

  }

  return(invisible(list(var1names=var1.names, var2names=var2.names)))
}


