

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
  #' @param isPlot Boolean on whether or not to generate a plot. If family is binomial,
  #'        an AUC plot is generated. Currently applicable to SELPCCA. Set isPlot to FALSE
  #'        for SIDA and SIDANet. If isPlot is TRUE and predictions
  #'        are from SIDA and SIDANet, AUC estimate coincides with balanced accuracy.
  #'        If family is gaussian, a scatter plot of observed vs
  #'        predicted is generate.
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
  #' @import ggplot2
  #' @import ROCit
  #' @import ggthemes
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
#'              isParallel=TRUE,ncores=NULL,gridMethod='RandomSearch',
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
#'  test.metrics=PerformanceMetrics(Y.pred,Ytest.in,family='binomial',isPlot=FALSE)
#'  print(test.metrics)


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


PerformanceMetrics = function(Y.pred, Y.test, family='binomial',isPlot=TRUE){
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


    if(isPlot==TRUE){
      # rocobj=pROC::roc(Y.test,Y.pred)
      # AUC.value=round(pROC::auc(rocobj),3)
      # print(
      #   pROC::ggroc(rocobj, legacy.axes = TRUE,alpha = 0.5, colour = "red", linetype = 2, linewidth=2)
      #   + xlab("FPR") +
      #     ylab("TPR") +
      #   ggplot2::geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="black", linetype=1)+
      #   ggplot2::ggtitle("ROC Plot") +
      #   ggplot2::annotate("text", x = 0.2, y = 0.5, label =paste0("AUC =", AUC.value), parse = TRUE)
      # # )
      # )
      #print(
        ROCit_obj <- ROCit::rocit(score=Y.pred,class=Y.test)
        ROCit_plot <- plot(ROCit_obj, legend=F, YIndex=F, col="red")
        text(0.75, 0.28, paste0("AUC: ",
                                round(ROCit_obj$AUC,4)))
        text(0.75, 0.2, paste0("Optimal Cutoff: ",
                               round(ROCit_plot$`optimal Youden Index point`[4],4)))
        text(0.75, 0.12, paste0("Sensitivity: ",
                                round(ROCit_plot$`optimal Youden Index point`[3],4)))
        text(0.75, 0.04, paste0("Specificity: ",
                                round(1-ROCit_plot$`optimal Youden Index point`[2],4)))
      #)

    }


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

         if(isPlot==TRUE){
           df=cbind.data.frame(Y.pred, Y.test)
           colnames(df)=c("Predicted.Y", "Observed.Y")
           print(
           ggplot2::ggplot(df, aes(x=Observed.Y, y=Predicted.Y)) +
             geom_point() +
             geom_smooth(method=lm) +#add linear trend line
             ylab(paste("Predicted Values")) +
             xlab(paste("Observed Values")) )  }
      }else{
        print('metrics available only for binary and gaussian outcomes')
      }
  return(metrics)

    }


