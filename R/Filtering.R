

#' Supervised Filtering
#'
#' @description Performs univariate supervised filtering on multi-source
#'  data. A separate model will be fit for each feature within each view of data
#'  and all features with p-values less than the specified threshold will be retained.
#'
#' @param X A list containing all data sources. Each row must represent a
#' subject and each column represents a feature.
#' @param Y An outcome vector of length equal to the number of rows in each view of X.
#' @param method Options are "linear" for linear regression, "logistic" for logistic
#' regression, "t.test" for a 2-sample unpaired T-test, or "kw" for a Kruskal-Wallis
#' test. Default is "linear".
#' @param padjust Boolean on whether or not to adjust pvalue for multiple testing.
#' Default is "F".
#' @param adjmethod Options are "holm", "hochberg",  "hommel", "bonferroni", "BH"
#' "BY","fdr","none". Default is "BH" if padjust is True.
#' @param thresh P-value threshold to determine which features to keep after filtering.
#' Default will keep all features with a p-value < 0.05.
#' @param center Boolean on whether or not to center the features prior to filtering.
#' @param scale Boolean on whether or not to scale the features to have variance 1 prior to filtering.
#' @param standardize Boolean on whether or not to center and scale the features to have mean 0 and variance 1 prior to filtering.
#' @param log2TransForm Boolean on whether or not to log2 transform the features prior to filtering. Will return an error if TRUE but data
#' have negative values.
#' @param Xtest Optional list containing test data. If included, filtering will be performed
#' only on the training data, X, but Xtest will be subsetted to the same group of features.
#'
#' @return A list containing the following (and others):
#'   \item{X}{List of the filtered X data}
#'   \item{Y}{Vector of the outcome}
#'   \item{Xtest}{List of the subsetted Xtest data}
#'   \item{method}{Method used for filtering}
#'   \item{pval.mat}{Dataset containing the calculated p-values for each feature, coefficients, and whether significant.}
#'
#' @importFrom umap umap
#' @importFrom stats binomial glm kruskal.test lm p.adjust prcomp t.test var
#'
#' @export
#' @examples
#' ##---- read in data
#' data(sidaData)
#'
#' Xdata=sidaData[[1]]
#' Y=sidaData[[2]]
#'
#' data.red=filter_supervised(Xdata, Y,  method="t.test", padjust=FALSE,adjmethod=NULL, thresh=0.05,
#' center=FALSE, scale=FALSE, standardize=FALSE, log2TransForm=FALSE, Xtest=NULL)
#'
#' ##-----Plot Result via UMAP
#' umapPlot(data.red)

#############################################################################################################
# Authors:
#   Elise Palzer, Unviersity of Minnesota
# Sandra E. Safo, Unversity of Minnesota
# created: 2023
#
# Copyright (C) 2023
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#############################################################################################################

filter_supervised <- function(X, Y, method="linear", padjust=FALSE,adjmethod="BH", thresh=0.05,
                              center=FALSE, scale=FALSE, standardize=FALSE, log2TransForm=FALSE, Xtest=NULL){
  if(is.null(adjmethod)){
    adjmethod="BH"
  }
  XOrig=X
  XtestOrig=Xtest

  for(i in 1:length(X)){

    if(log2TransForm){
      #check data
      if(min(X[[1]])<0){
        stop("'negative values in data, cannot log2-transform'",
             call. = FALSE)
      }
      X[[i]] <- apply(X[[i]], 2, log2)
      if(length(XtestOrig)>0){
        Xtest[[i]] <- apply(Xtest[[i]], 2, log2)
      }
    }

    if(standardize){
      mean.vec <- apply(X[[i]], 2, mean, na.rm=TRUE)
      var.vec <- apply(X[[i]], 2, function(x) sqrt(var(x,na.rm=TRUE)))
      X[[i]] <- t( (t(X[[i]]) - mean.vec) /var.vec)
      if(length(XtestOrig)>0){
        Xtest[[i]] <- t( (t(Xtest[[i]])-mean.vec)/var.vec)
      }
      #X[[i]]=apply(X[[i]], 2, function(x) scale(x))
    }


    if(center){
      mean.vec <- apply(X[[i]], 2, mean,na.rm=TRUE)
      X[[i]] <- t(t(X[[i]]) - mean.vec )
      if(length(XtestOrig)>0){
        Xtest[[i]] <- t(t(Xtest[[i]]) - mean.vec)
      }
    }
    if(scale){
      var.vec <- apply(X[[i]], 2, function(x) sqrt(var(x,na.rm=TRUE)))
      X[[i]] <- t(t(X[[i]])/var.vec)
      if(length(XtestOrig)>0){
        Xtest[[i]] <- t(t(Xtest[[i]])/var.vec)
      }
    }

    # if(standardize){
    #   mean.vec <- apply(XOrig[[i]], 2, mean)
    #   var.vec <- apply(XOrig[[i]], 2, function(x) sqrt(var(x)))
    #   X[[i]] <- t( (t(XOrig[[i]]) - mean.vec) /var.vec)
    #   # if(length(XtestOrig)>0){
    #   #   Xtest[[i]] <- t(t(Xtest[[i]])/var.vec)
    #   # }
    #   #X[[i]]=apply(X[[i]], 2, function(x) scale(x))
    # }



  }

  coef.mat <- pval.mat <- pval.adj.mat <- red.mat <- mean.mat <- X.red <- Xtest.red <- list()
  for(i in 1:length(X)){
    if(method == "linear"){
      coef.mat[[i]] <- apply(X[[i]], 2, function(x)
        lm(as.numeric(Y) ~ as.numeric(x))$coefficients[2])
      pval.mat[[i]] <- apply(X[[i]], 2, function(x)
        summary(lm(as.numeric(Y) ~ as.numeric(x)))$coefficients[2,4])
    }else if(method == "logistic"){
      coef.mat[[i]] <- apply(X[[i]], 2, function(x)
        glm(factor(Y) ~ as.numeric(x), family = binomial())$coefficients[2])
      pval.mat[[i]] <- apply(X[[i]], 2, function(x)
        summary(glm(factor(Y) ~ as.numeric(x), family = binomial()))$coefficients[2,4])
    }else if(method == "t.test"){
      coef.mat[[i]] <- apply(X[[i]], 2, function(x)
        t.test(as.numeric(x) ~ factor(Y))$estimate[1]-t.test(as.numeric(x) ~ factor(Y))$estimate[2])
      pval.mat[[i]] <- apply(X[[i]], 2, function(x)
        t.test(as.numeric(x) ~ factor(Y))$p.value)
      #coef.mat[[i]] <- pval.mat[[i]]*NA
    }else if(method == "kw"){
      coef.mat[[i]] <- apply(X[[i]], 2, function(x)
        kruskal.test(as.numeric(x) ~ factor(Y))$statistic)
      pval.mat[[i]] <- apply(X[[i]], 2, function(x)
        kruskal.test(as.numeric(x) ~ factor(Y))$p.value)
    }else{
      cat("Warning: Method does not exist")
      quit()
    }



    if(padjust==TRUE){
      pval.adj.mat[[i]] <-p.adjust(pval.mat[[i]], method=adjmethod)
        red.mat[[i]] <- 1:ncol(X[[i]]) %in% which(pval.adj.mat[[i]] < thresh)
        X.red[[i]] <- X[[i]][,red.mat[[i]]]
        pval.mat[[i]]=pval.adj.mat[[i]]
    }else{
      red.mat[[i]] <- 1:ncol(X[[i]]) %in% which(pval.mat[[i]] < thresh)
      X.red[[i]] <- X[[i]][,red.mat[[i]]]
    }

    if(length(Xtest)>0){
      Xtest.red[[i]] <- Xtest[[i]][,red.mat[[i]]]
    }
  }

  # #label for ttest
  # if(method=="t.test"){
  #   che=as.data.frame(t.test(as.numeric(X[[1]][,1]) ~ factor(Y))$estimate)
  #   mean.diff.label=rownames(che)[1]-rownames(che)[2]
  # }

  temp <- data.frame(Coef = NULL,
                     Pval = NULL,
                     Keep = NULL,
                     View = NULL)
  t1=list()
  for(i in 1:length(X)){

    t1[[i]] <- data.frame(Coef = coef.mat[[i]],
                          Pval = pval.mat[[i]],
                          Keep = red.mat[[i]],
                          View = i)
    temp <- rbind(temp,t1[[i]])


    # if(method=="t.test"){
    #   colnames(temp)[1]="Mean Difference"
    # }else if(method=="linear"){
    #   colnames(temp)[1]="Coef"
    # }else if(method=="logistic"){
    #   colnames(temp)[1]="Log ORs"
    # }else if(method=="kw"){
    #   colnames(temp)[1]="Coef"
    # }


  }

  #print results for top 10 significant variables by views
  for(i in 1:length(X)){
    mydata2=temp[temp$View==i,]
    mydata2.Sorted= mydata2[order(mydata2[,2]),]
    mydatathresh=mydata2.Sorted[mydata2.Sorted[,2]< thresh, ]
    if(length(mydatathresh)==0){
      print(paste0("No variable is significnat for View ",i))
      print(mydatathresh)
    }else{
      print(paste0("Printing top 10 results for  significant variables for View ",i))
      print(mydatathresh[1:10,])
    }
  }

  result <- list(X=X.red, Y=Y,
                 Xtest=Xtest.red,
                 method=method,
                 pval.mat=temp,
                 pval.matlist=t1,
                 significant.thresh=thresh,
                 padjust=padjust,
                 adjmethod=adjmethod,
                 X_Original=XOrig,
                 Xtest_Original=XtestOrig,
                 Center=center,
                 Scale=scale,
                 Log2Transform=log2TransForm,
                 Standardize=standardize)
  return(result)
}


#' Unsupervised Filtering
#'
#' @description Performs univariate unsupervised filtering on multi-source
#'  data. A separate model will be fit for each feature within each view of data
#'  and all features with p-values less than the specified threshold will be retained.
#'
#' @param X A list containing all data sources. Each row must represent a
#' subject and each column represents a feature.
#' @param method Options are "variance" which will keep the \code{pct.keep} percent of
#'  features with the highest variance, and "IQR", which will keep the features with the
#'  median amount of variance (+/- pct.keep/2). Default is "variance".
#' @param pct.keep Percent of variables to keep in each view of data. Default is 10\%.
#' @param center Boolean on whether or not to center the features after filtering.
#' @param scale Boolean on whether or not to scale the features after filtering.
#' @param standardize Boolean on whether or not to center and scale the features to have mean 0 and variance 1 after filtering.
#' @param log2TransForm Boolean on whether or not to log2 transform the features prior to filtering. Will return an error if TRUE but data
#' have negative values.
#' @param Xtest Optional list containing test data. If included, filtering will be performed
#' only on the training data, X, but Xtest will be subsetted to the same group of features.
#'
#' @return A list containing the following
#'   \item{X}{List of the filtered X data}
#'   \item{Xtest}{List of the subsetted Xtest data}
#'   \item{method}{Method used for filtering}
#'   \item{var.mat}{Dataset containing the calculated mean and variances for each feature.}
#'
#' @importFrom stats binomial glm kruskal.test lm p.adjust prcomp t.test var
#' @export
#' @examples
#' ##---- read in data
#' data(sidaData)
#'
#' Xdata=sidaData[[1]]
#'
#' data.red=filter_unsupervised(Xdata,  method="variance", pct.keep=10,
#'                     center=FALSE, scale=FALSE, standardize=FALSE,  log2TransForm=FALSE,
#'                      Xtest=NULL)

#############################################################################################################
# Authors:
#   Elise Palzer, University of Minnesota
#   Sandra E. Safo, University of Minnesota
# created: 2023
#
# Copyright (C) 2023
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#############################################################################################################


filter_unsupervised <- function(X, method="variance", pct.keep=10,
                              center=FALSE, scale=FALSE,standardize=FALSE,  log2TransForm=FALSE,
                              Xtest=NULL){
  var.mat <- mean.mat <- red.mat <- X.red <- Xtest.red <- list()



  for(i in 1:length(X)){

    #A boolean on whether or not to log2() transform all X datasets prior to filtering.
    if(log2TransForm){
      if(min(X[[1]])<0){
        stop("'negative values in data, cannot log2-transform'",
             call. = FALSE)
      }
      X[[i]] <- apply(X[[i]], 2, log2)
      if(length(Xtest)>0){
        Xtest[[i]] <- apply(Xtest[[i]], 2, log2)
      }

    }



    var.mat[[i]] <- apply(X[[i]], 2, function(x) var(as.numeric(x),na.rm=TRUE))
    mean.mat[[i]] <- apply(X[[i]], 2, function(x) mean(as.numeric(x),na.rm=TRUE))
  }
  if(method == "variance"){
    for(i in 1:length(X)){
      keep.var <- which(var.mat[[i]] > quantile(var.mat[[i]], 1-(pct.keep/100)))
      red.mat[[i]] <- 1:ncol(X[[i]]) %in% keep.var
      X.red[[i]] <- X[[i]][,keep.var]
      if(length(Xtest) > 0){
        Xtest.red[[i]] <- Xtest[[i]][,keep.var]
      }
    }
  }else if(method == "IQR"){
    for(i in 1:length(X)){
      bad.var <- which(var.mat[[i]] > quantile(var.mat[[i]], .5+(pct.keep/200)) |
                         var.mat[[i]] < quantile(var.mat[[i]], .5-(pct.keep/200)) )
      red.mat[[i]] <- 1:ncol(X[[i]]) %in% bad.var
      red.mat[[i]] <- (red.mat[[i]] == FALSE)
      X.red[[i]] <- X[[i]][,-bad.var]
      if(length(Xtest) > 0){
        Xtest.red[[i]] <- Xtest[[i]][,-bad.var]
      }
    }
  }
  #Plot
  temp <- data.frame(Mean = NULL,
                     Variance = NULL,
                     Keep = NULL,
                     View = NULL)
  for(i in 1:length(X)){
    t1 <- data.frame(Mean = mean.mat[[i]],
                     Variance = var.mat[[i]],
                     Keep = red.mat[[i]],
                     View = i)
    temp <- rbind(temp,t1)
  }

  # for(i in 1:length(X.red)){
  #   X.red.temp=X.red
  #   Xtest.red.temp=Xtest.red
  #
  #   # if(log2TransForm){
  #   #   X.red[[i]] <- apply(X.red.temp[[i]], 2, log2)
  #   #   if(length(Xtest.red.temp)>0){
  #   #     Xtest.red[[i]] <- apply(Xtest.red.temp[[i]], 2, log2)
  #   #   }
  #
  #   }
    if(standardize){
      mean.vec <- apply(X.red[[i]], 2, mean, na.rm=TRUE)
      var.vec <- apply(X.red[[i]], 2, function(x) sqrt(var(x,na.rm=TRUE)))
      X.red[[i]] <- t( (t(X.red[[i]]) - mean.vec) /var.vec)
      if(length(Xtest)>0){
        Xtest.red[[i]] <- t( (t(Xtest.red[[i]])-mean.vec)/var.vec)
      }

    }

    if(center){
      mean.vec <- apply(X.red[[i]], 2, mean,na.rm=TRUE)
      X.red[[i]] <- t(t(X.red[[i]]) - mean.vec)
      if(length(Xtest)>0){
        Xtest.red[[i]] <- t(t(Xtest.red[[i]]) - mean.vec)
      }
    }
    if(scale){
      var.vec <- apply(X.red[[i]], 2, function(x) sqrt(var(x,na.rm=TRUE)))
      X.red[[i]] <- t(t(X.red[[i]])/var.vec)
      if(length(Xtest)>0){
        Xtest.red[[i]] <- t(t(Xtest.red[[i]])/var.vec)
      }
    }





  result <- list(X=X.red,
                 Xtest=Xtest.red,
                 method=method,
                 var.mat=temp,
                 Scale=scale,
                 Center=center,
                 Log2Transform=log2TransForm,
                 Standardize=standardize)
  return(result)
}


#' @param X A list containing all data sources. Each row must represent a
#' subject and each column represents a feature.
#' @param Y An outcome vector of length equal to the number of rows in each view of X.
#' @param method Options are "linear" for linear regression, "logistic" for logistic
#' regression, "t.test" for a 2-sample unpaired T-test, or "kw" for a Kruskal-Wallis
#' test. Default is "linear".
#' @param padjust Boolean on whether or not to adjust pvalue for multiple testing.
#' Default is "F".
#' @param adjmethod Options are "holm", "hochberg",  "hommel", "bonferroni", "BH"
#' "BY","fdr","none". Default is "BH" if padjust is True.
#' @param thresh P-value threshold to determine which features to keep after filtering.
#' Default will keep all features with a p-value < 0.05.
#' @param center Boolean on whether or not to center the features prior to filtering.
#' @param scale Boolean on whether or not to scale the features to have variance 1 prior to filtering.
#' @param standardize Boolean on whether or not to center and scale the features to have mean 0 and variance 1 prior to filtering.
#' @param log2TransForm Boolean on whether or not to log2 transform the features prior to filtering. Will return an error if TRUE but data
#' have negative values.
#' @param Xtest Optional list containing test data. If included, filtering will be performed
#' only on the training data, X, but Xtest will be subsetted to the same group of features.
#' @param CovAdjust matrix of covariates with which to adjust. default = NULL
#' @importFrom stats lm glm p.adjust  t.test kruskal.test
filter_supervised_covadjust <- function(X, Y, method="linear", padjust=F,adjmethod="BH", thresh=0.05,
                                        center=F, scale=F, standardize=F, log2TransForm=F, Xtest=NULL, CovAdjust=NULL){

  if(length(thresh)==1){
    thresh=matrix(thresh,nrow=length(X),ncol=1)
  }

  if(is.null(adjmethod)){
    adjmethod="BH"
  }
  XOrig=X
  XtestOrig=Xtest

  for(i in 1:length(X)){

    if(log2TransForm){
      #check data
      if(min(X[[1]])<0){
        stop("'negative values in data, cannot log2-transform'",
             call. = FALSE)
      }
      X[[i]] <- apply(X[[i]], 2, log2)
      if(length(XtestOrig)>0){
        Xtest[[i]] <- apply(Xtest[[i]], 2, log2)
      }
    }

    if(standardize){
      mean.vec <- apply(X[[i]], 2, mean, na.rm=TRUE)
      var.vec <- apply(X[[i]], 2, function(x) sqrt(var(x,na.rm=TRUE)))
      X[[i]] <- t( (t(X[[i]]) - mean.vec) /var.vec)
      if(length(XtestOrig)>0){
        Xtest[[i]] <- t( (t(Xtest[[i]])-mean.vec)/var.vec)
      }
      #X[[i]]=apply(X[[i]], 2, function(x) scale(x))
    }


    if(center){
      mean.vec <- apply(X[[i]], 2, mean,na.rm=TRUE)
      X[[i]] <- t(t(X[[i]]) - mean.vec )
      if(length(XtestOrig)>0){
        Xtest[[i]] <- t(t(Xtest[[i]]) - mean.vec)
      }
    }
    if(scale){
      var.vec <- apply(X[[i]], 2, function(x) sqrt(var(x,na.rm=TRUE)))
      X[[i]] <- t(t(X[[i]])/var.vec)
      if(length(XtestOrig)>0){
        Xtest[[i]] <- t(t(Xtest[[i]])/var.vec)
      }
    }

    # if(standardize){
    #   mean.vec <- apply(XOrig[[i]], 2, mean)
    #   var.vec <- apply(XOrig[[i]], 2, function(x) sqrt(var(x)))
    #   X[[i]] <- t( (t(XOrig[[i]]) - mean.vec) /var.vec)
    #   # if(length(XtestOrig)>0){
    #   #   Xtest[[i]] <- t(t(Xtest[[i]])/var.vec)
    #   # }
    #   #X[[i]]=apply(X[[i]], 2, function(x) scale(x))
    # }



  }

  coef.mat <- pval.mat <- pval.adj.mat <- red.mat <- mean.mat <- X.red <- Xtest.red <- mydata <- list()
  for(i in 1:length(X)){
    if(method == "linear"){
      if (!is.null(CovAdjust)){
        mydata[[i]]=cbind.data.frame(Y,X[[i]], CovAdjust)
        merged_cov <- paste(colnames(CovAdjust),collapse=' + ')
        colY=colnames(mydata[[1]])[1]
        fmt <- c(paste(c(paste(c(colY,"~ %s"), collapse = ' '),merged_cov), collapse = ' + '))
        fo.strings <- as.data.frame(outer(fmt, colnames(X[[i]]), sprintf))
        coef.mat[[i]]=apply(fo.strings, 2, function(formulas)
          lm(formulas,data = mydata[[i]])$coefficients[2])
        pval.mat[[i]]= apply(fo.strings, 2, function(formulas)
          summary(lm(formulas, data = mydata[[i]]))$coefficients[2,4])
      }else if(is.null(CovAdjust)){
        coef.mat[[i]] <- apply(X[[i]], 2, function(x)
          lm(as.numeric(Y) ~ as.numeric(x))$coefficients[2])
        pval.mat[[i]] <- apply(X[[i]], 2, function(x)
          summary(lm(as.numeric(Y) ~ as.numeric(x)))$coefficients[2,4])
      }
    }else if(method == "logistic"){
      if (!is.null(CovAdjust)){
        mydata[[i]]=cbind.data.frame(Y,X[[i]], CovAdjust)
        merged_cov <- paste(colnames(CovAdjust),collapse=' + ')
        colY=colnames(mydata[[1]])[1]
        fmt <- c(paste(c(paste(c(colY,"~ %s"), collapse = ' '),merged_cov), collapse = ' + '))
        fo.strings <- as.data.frame(outer(fmt, colnames(X[[i]]), sprintf))
        coef.mat[[i]]=apply(fo.strings, 2, function(formulas)
          glm(formulas, family=binomial(link='logit'), data = mydata[[i]])$coefficients[2])
        che=glm(Y ~ A0A075B6H9 + Age + Sex,family=binomial(link='logit'), data = mydata[[i]] )$coefficients[2]
        pval.mat[[i]]= apply(fo.strings, 2, function(formulas)
          summary(glm(formulas, family=binomial(link='logit'),data = mydata[[i]]))$coefficients[2,4])

        # for(j in 1:264){
        #   print(j)
        #   che=glm(fo.strings[j],family=binomial(link='logit'), data = mydata[[i]] )$coefficients[2]
        # }


      }else if(is.null(CovAdjust)){
        coef.mat[[i]] <- apply(X[[i]], 2, function(x)
          glm(factor(Y) ~ as.numeric(x), family = binomial())$coefficients[2])
        pval.mat[[i]] <- apply(X[[i]], 2, function(x)
          summary(glm(factor(Y) ~ as.numeric(x), family = binomial()))$coefficients[2,4])
      }
    }else if(method == "t.test"){
      coef.mat[[i]] <- apply(X[[i]], 2, function(x)
        t.test(as.numeric(x) ~ factor(Y))$estimate[1]-t.test(as.numeric(x) ~ factor(Y))$estimate[2])
      pval.mat[[i]] <- apply(X[[i]], 2, function(x)
        t.test(as.numeric(x) ~ factor(Y))$p.value)
      #coef.mat[[i]] <- pval.mat[[i]]*NA
    }else if(method == "kw"){
      coef.mat[[i]] <- apply(X[[i]], 2, function(x)
        kruskal.test(as.numeric(x) ~ factor(Y))$statistic)
      pval.mat[[i]] <- apply(X[[i]], 2, function(x)
        kruskal.test(as.numeric(x) ~ factor(Y))$p.value)
    }else{
      cat("Warning: Method does not exist")
      quit()
    }
    if(padjust==TRUE){
      pval.adj.mat[[i]] <-p.adjust(pval.mat[[i]], method=adjmethod)
      red.mat[[i]] <- 1:ncol(X[[i]]) %in% which(pval.adj.mat[[i]] < thresh[i])
      X.red[[i]] <- X[[i]][,red.mat[[i]]]
      pval.mat[[i]]=pval.adj.mat[[i]]
    }else{
      red.mat[[i]] <- 1:ncol(X[[i]]) %in% which(pval.mat[[i]] < thresh[i])
      X.red[[i]] <- X[[i]][,red.mat[[i]]]
    }

    if(length(Xtest)>0){
      Xtest.red[[i]] <- Xtest[[i]][,red.mat[[i]]]
    }
  }

  # #label for ttest
  # if(method=="t.test"){
  #   che=as.data.frame(t.test(as.numeric(X[[1]][,1]) ~ factor(Y))$estimate)
  #   mean.diff.label=rownames(che)[1]-rownames(che)[2]
  # }

  temp <- data.frame(Coef = NULL,
                     Pval = NULL,
                     Keep = NULL,
                     View = NULL)
  t1=list()
  for(i in 1:length(X)){

    t1[[i]] <- data.frame(Coef = coef.mat[[i]],
                          Pval = pval.mat[[i]],
                          Keep = red.mat[[i]],
                          View = i)
    temp <- rbind(temp,t1[[i]])



    # if(method=="t.test"){
    #   colnames(temp)[1]="Mean Difference"
    # }else if(method=="linear"){
    #   colnames(temp)[1]="Coef"
    # }else if(method=="logistic"){
    #   colnames(temp)[1]="Log ORs"
    # }else if(method=="kw"){
    #   colnames(temp)[1]="Coef"
    # }


  }

  #print results for top 10 significant variables by views
  for(i in 1:length(X)){
    mydata2=temp[temp$View==i,]
    mydata2.Sorted= mydata2[order(mydata2[,2]),]
    mydatathresh=mydata2.Sorted[mydata2.Sorted[,2]< thresh[i], ]
    if(length(mydatathresh)==0){
      print(paste0("No variable is significnat for View ",i))
      print(mydatathresh)
    }else{
      print(paste0("Printing top 10 results for  significant variables for View ",i))
      print(mydatathresh[1:10,])
    }
  }

  result <- list(X=X.red, Y=Y,
                 Xtest=Xtest.red,
                 method=method,
                 #pval.mat=temp,
                 pval.mat=t1,
                 significant.thresh=thresh,
                 padjust=padjust,
                 adjmethod=adjmethod,
                 X_Original=XOrig,
                 Xtest_Original=XtestOrig,
                 Center=center,
                 Scale=scale,
                 Log2Transform=log2TransForm,
                 Standardize=standardize,
                 CovAdjust=CovAdjust)
  return(result)
}




