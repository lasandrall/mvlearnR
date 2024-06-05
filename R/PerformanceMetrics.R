
#############################################################################################################
# Authors:
#  Sandra Safo, Unviersity of Minnesota
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


  #' Performance Metrics
  #'
  #' @description Estimates performance metrics for a predicted model. Currently works
  #'  for binary and continuous outcomes
  #'
  #'
  #' @param Y.pred A vector of predicted values. For SELPCCA, this is a vector of
  #'              predicted probabilities. For SIDA and SIDANet, this is a vector of predicted class.
  #' @param Y.test A vector of observed test values.
  #' @param family A string to denote the family for which  metrics should be provided.
  #'        Options are "gaussian", "binomial".
  #' @details For a binary outcome, we provide the following metrics:
  #'          "Accuracy","Error rate","Sensitivity", "Specificity", "Matthews Correlation Coefficient (MCC)",
  #'          "Balanced Accuracy","Balanced Error Rate", "F1 Score", "False.Discovery.Rate", and
  #'          "Positive Predictive Value".
  #'
  #'          For a  continuous outcome, we provide the following metrics:
  #'          "Mean Squared Error","Root Mean.Squared Error","Relative Squared Error",
  #'.         "Root Relative Squared.Error", "Root Absolute Error", "Mean Absolute Error".
  #' @return An output of performance metrics:
  #'   \item{Metrics}{A table of estimated metrics}
  #' @seealso \code{\link{cvSIDA}} \code{\link{selpscca.pred}} \code{\link{predict.SELPCCA}}
  #'
  #' @import ROCit
  #' @importFrom ROCit rocit
  #' @export
  #'
  #' @examples
  #' data(sidaData)
  ##---- call sida algorithm to estimate discriminant vectors, and predict on testing data
  #' Xdata=sidaData[[1]]
#' Y=sidaData[[2]]
#' Xtestdata=sidaData[[3]]
#'Ytest=sidaData[[4]]
#' ##---- call cross validation
#'  mycv=cvSIDA(Xdata,Y,withCov=FALSE,plotIt=FALSE, Xtestdata=Xtestdata,Ytest=Ytest,
#'              isParallel=FALSE,ncores=NULL,gridMethod='RandomSearch',
#'              AssignClassMethod='Joint',nfolds=5,ngrid=8,standardize=TRUE,
#'             maxiteration=20, weight=0.5,thresh=1e-03)
#' #check output
#'  test.error=mycv$sidaerror
#'  test.correlation=mycv$sidacorrelation
#'  optTau=mycv$optTau
#'  hatalpha=mycv$hatalpha
#'  #train metrics
#'  Y.pred=mycv$PredictedClass.train-1 #to get this in 0 and 1
#'  Y.train=Y-1 #to get this in 0 and 1
#'  train.metrics=PerformanceMetrics(Y.pred,Y.train,family='binomial',isPlot=FALSE)
#'
#'  print(train.metrics)
#'  #obtain predicted class
#'  Y.pred=mycv$PredictedClass-1 #to get this in 0 and 1
#'  Ytest.in=Ytest-1 #to get this in 0 and 1
#'  test.metrics=PerformanceMetrics(Y.pred,Ytest.in,family='binomial')
#'  print(test.metrics)
PerformanceMetrics = function(Y.pred, Y.test, family='binomial'){
  if (!(family %in% c("binomial","gaussian"))){
    stop('metrics available only for binary and gaussian outcomes')
  }
  if (family=='binomial'){
    
    Y.pred.class=round(Y.pred)
    TP=sum(Y.pred.class[Y.test==1]==1)
    TN=sum(Y.pred.class[Y.test==0]==0)
    FP=sum(Y.pred.class[Y.test==0]==1)
    FN=sum(Y.pred.class[Y.test==1]==0)
    
    ACC=(TP + TN)/(TP+TN +FP+FN)
    Error.rate=1-ACC
    
    BAC=0.5*(TP/(TP + FN) + TN/(TN + FP))
    
    BER= 0.5*(FN/(TP + FN) + FP/(TN + FP))
    
    F1= 2*TP/(2*TP + FP +FN)
    
    FDR=FP/(TP + FP)
    
    A=TP + FP
    B=TP+ FN
    C=TN + FP
    D=TN + FN
    if(min(c(A, B,C, D))==0){
      MCC=( (TP * TN) - (FP * FN))/1
    }else{
      MCC= (TP * TN - FP * FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    }
    
    PPV=TP/(TP + FP)
    SENS= TP/(TP + FN)
    SPEC= TN/(TN + FP)
    metrics=rbind.data.frame(ACC,Error.rate,SENS, SPEC, MCC, BAC,BER, F1, FDR, PPV)
    rownames(metrics)=c("Accuracy","Error.rate","Sensitivity", "Specificity", "Matthews.Correlation.Coefficient",
                        "Balanced.Accuracy","Balanced.Error.Rate", "F1.Score", "False.Discovery.Rate", "Positive.Predictive.Value")
    colnames(metrics)="Metrics"
    
  }else if(family=='gaussian'){
    MSE= mean((Y.pred-Y.test)^2)
    RMSE = sqrt(mean((Y.pred-Y.test)^2))
    
    RRSE = sqrt(sum((Y.pred-Y.test)^2)/sum((Y.test - mean(Y.test))^2)) #root relative squared error
    
    MAE= mean(abs(Y.pred-Y.test)) #mean absolute error
    RAE= sum(abs(Y.pred-Y.test)) / sum(abs(Y.test-mean(Y.test))) #relative absolute error
    
    RSE=sum((Y.pred-Y.test)^2)/sum( (Y.test-mean(Y.test))^2) #relative squared error
    
    metrics=rbind.data.frame(MSE, RMSE, RSE, RRSE, RAE, MAE )
    rownames(metrics)=c("Mean.Squared.Error","Root.Mean.Squared.Error","Relative.Squared.Error",
                        "Root.Relative.Squared.Error", "Root.Absolute.Error", "Mean.Absolute.Error")
    colnames(metrics)="Metrics"
    
  }
  return(metrics)
}


#' Performance Metrics Plot
#'
#' @description Creates a ggplot showing either an AUC curve (family == "binomial") or a scatter of predicted vs. observed values (family == "gaussian")
#'
#'
#' @param Y.pred A vector of predicted values. For SELPCCA, this is a vector of
#'              predicted probabilities. For SIDA and SIDANet, this is a vector of predicted class.
#' @param Y.test A vector of observed test values.
#' @param family A string to denote the family for which  metrics should be provided.
#'        Options are "gaussian", "binomial".
#' @details For a binary outcome, plots an AUC curve. For a continuous (gaussian) outcome, plots a scatter of observed vs. predicted values. 
#' @return A ggplot object
#' @seealso \code{\link{cvSIDA}} \code{\link{selpscca.pred}} \code{\link{predict.SELPCCA}}\code{\link{PerformanceMetrics}}
#'
#' @import ROCit
#' @importFrom ROCit rocit
#' @export
#'
#' @examples
#' data(sidaData)
##---- call sida algorithm to estimate discriminant vectors, and predict on testing data
#' Xdata=sidaData[[1]]
#' Y=sidaData[[2]]
#' Xtestdata=sidaData[[3]]
#'Ytest=sidaData[[4]]
#' ##---- call cross validation
#'  mycv=cvSIDA(Xdata,Y,withCov=FALSE,plotIt=FALSE, Xtestdata=Xtestdata,Ytest=Ytest,
#'              isParallel=FALSE,ncores=NULL,gridMethod='RandomSearch',
#'              AssignClassMethod='Joint',nfolds=5,ngrid=8,standardize=TRUE,
#'             maxiteration=20, weight=0.5,thresh=1e-03)
#' #check output
#'  test.error=mycv$sidaerror
#'  test.correlation=mycv$sidacorrelation
#'  optTau=mycv$optTau
#'  hatalpha=mycv$hatalpha
#'  #train metrics
#'  Y.pred=mycv$PredictedClass.train-1 #to get this in 0 and 1
#'  Y.train=Y-1 #to get this in 0 and 1
#'  train.metrics=PerformanceMetrics(Y.pred,Y.train,family='binomial',isPlot=FALSE)
#'
#'  print(train.metrics)
#'  #obtain predicted class
#'  Y.pred=mycv$PredictedClass-1 #to get this in 0 and 1
#'  Ytest.in=Ytest-1 #to get this in 0 and 1
#'  PerformanceMetricsPlot(Y.pred,Ytest.in,family='binomial')
PerformanceMetricsPlot = function(Y.pred, Y.test, family = "binomial"){
  if (!(family %in% c("binomial","gaussian"))){
    stop('metrics available only for binary and gaussian outcomes')
  }
  
  if (family == "binomial"){
    ROCit_obj <- ROCit::rocit(score=Y.pred,class=Y.test)
    maxIndex <- which.max(ROCit_obj$TPR - ROCit_obj$FPR)
    FPR <- 1 - ROCit_obj$FPR[maxIndex]
    TPR <- ROCit_obj$TPR[maxIndex]
    cutoff <- ROCit_obj$Cutoff[maxIndex]
    
    data.frame(FPR = ROCit_obj$FPR, TPR = ROCit_obj$TPR) %>%
      ggplot(aes(x = FPR, y = TPR))+
      geom_line(color = "red")+
      geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
      annotate("text", x = 0.75, y = 0.28, label = paste0("AUC: ",
                                                          round(ROCit_obj$AUC,4)))+
      annotate("text", x = 0.75, y = 0.2, 
               label = paste0("Optimal Cutoff: ", round(cutoff,4)))+
      annotate("text", x = 0.75, y = 0.12, 
               label = paste0("Sensitivity: ",round(TPR,4)))+
      annotate("text", x = 0.75, y = 0.04, 
               label = paste0("Specificity: ",round(FPR,4)))+
      theme_bw()+
      labs(x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)")
  }else if (family == "gaussian"){
    df=cbind.data.frame(Y.pred, Y.test)
    colnames(df)=c("Predicted.Y", "Observed.Y")
    ggplot(df, aes(x=Observed.Y, y=Predicted.Y)) +
        geom_point() +
        geom_smooth(method=lm) +#add linear trend line
        ylab(paste("Predicted Values")) +
        xlab(paste("Observed Values"))+
      theme_bw()
  }
}


