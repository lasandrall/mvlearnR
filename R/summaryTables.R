#' Deviance Table for predicted, supervised data.
#'
#' @description Produced a deviance table for supervised predictions. 
#'
#' @param fit the output from the filter.supervised() function
#' @param Xtestdata1 a matrix or data.frame with the first view's test X data
#' @param Xtestdata2 a matrix or data.frame with the second view's test X data
#' @param Ytest a vector of test response values
#'
#' @return A data.frame
#'
#' @export
#'
#' @examples
#'
devianceTable = function(fit, Xtestdata1, Xtestdata2, Ytest){
  method = class(fit)
  
  if (!(method %in% c("SIDA", "SIDANet", "SELPCCA"))){
    stop("fit must be the result of cvsida, sidanet, or cvselpcca")
  }
  # results from SIDA or SIDANet are come directly from the fit
  # so do those and return early
  if(method=="SIDA"){
    return(
      data.frame(`SIDA Test Classification Error` = round(fit$sidaerror,4),
                 `SIDA Correlation` = round(fit$sidacorrelation,4))
    )
  }else if(method=="SIDANet"){
    return(
      data.frame(`SIDANet Test Classification Error` = round(fit$sidaneterror,4),
                      `SIDANet Correlation` = round(fit$sidanetcorrelation,4))
    )
  }
  ### Otherwise we're dealing with SELPCCA which is more complicated
  family = fit$family
  modeled_separately = fit$model.separately
  temp.test <- list(Xtestdata1, Xtestdata2)
  for(i in 1:length(temp.test)){
    temp.test[[i]]  <- matrix(as.numeric(as.character(unlist(temp.test[[i]]))), ncol=ncol(temp.test[[i]]))
  }
  if(family=="gaussian"){
    fit.pred <- predict(fit, newdata=temp.test[[1]], 
                        newdata2=temp.test[[2]])
    if(modeled_separately){
      df2a <- data.frame(loglik = round(as.numeric(logLik(fit$mod.fit[[1]])), 3),
                         Deviance = round(deviance(fit$mod.fit[[1]]),3),
                         AIC = round(AIC(fit$mod.fit[[1]]),3),
                         BIC = round(BIC(fit$mod.fit[[1]]),3),
                         McFadden_R2 = round(as.numeric(pscl::pR2(fit$mod.fit[[1]])['McFadden']),3),
                         `Train MSE`=round(sum((fit$mod.fit[[1]]$fitted.values - dat$Y)^2)/length(dat$Y),3),
                         `Test MSE`= round(sum((fit.pred[[1]] - dat$Ytest)^2)/length(fit.pred[[2]] ),3))
      df2b <- data.frame(loglik = round(as.numeric(logLik(fit$mod.fit[[2]])), 3),
                         Deviance = round(deviance(fit$mod.fit[[2]]),3),
                         AIC = round(AIC(fit$mod.fit[[2]]),3),
                         BIC = round(BIC(fit$mod.fit[[2]]),3),
                         McFadden_R2 = round(as.numeric(pscl::pR2(fit$mod.fit[[2]])['McFadden']),3),
                         `Train MSE`=round(sum((fit$mod.fit[[2]]$fitted.values - dat$Y)^2)/length(dat$Y),3),
                         `Test MSE`= round(sum((fit.pred[[2]] - dat$Ytest)^2)/length(fit.pred[[2]]),3))
      df2 <- rbind(df2a, df2b)
    }else{
      df2 <- data.frame(loglik = round(as.numeric(logLik(fit$mod.fit)), 3),
                        Deviance = round(deviance(fit$mod.fit),3),
                        AIC = round(AIC(fit$mod.fit),3),
                        BIC = round(BIC(fit$mod.fit),3),
                        McFadden_R2 = round(as.numeric(pscl::pR2(fit$mod.fit)['McFadden']),3),
                        `Train MSE`=round(sum((fit$mod.fit$fitted.values - dat$Y)^2)/length(dat$Y),3),
                        `Test MSE`= round(sum((fit.pred$pred.mod - dat$Ytest)^2)/length(fit.pred$pred.mod ),3))
    }
    df2
  }else if(family=="binomial"){
    fit.pred <- predict(fit, newdata=temp.test[[1]], 
                        newdata2=temp.test[[2]], type="response")
    if(modeled_separately){
      df2a <- data.frame(loglik = round(as.numeric(logLik(fit$mod.fit[[1]])),3),
                         Deviance =round( deviance(fit$mod.fit[[1]]),3),
                         AIC = round(AIC(fit$mod.fit[[1]]),3),
                         BIC = round(BIC(fit$mod.fit[[1]]),3),
                         McFadden_R2 = round(as.numeric(pscl::pR2(fit$mod.fit[[1]])['McFadden']),3),
                         `Train Misclassification Rate`=round(sum(abs(round(fit$mod.fit[[1]]$fitted.values) - fit$data.matrix$Y))/length(fit$data.matrix$Y),3),
                         `Test Misclassification Rate`= round(sum(abs(round(fit.pred[[1]]) - Ytest))/length(fit.pred[[1]] ),3))
      df2b <- data.frame(loglik = round(as.numeric(logLik(fit$mod.fit[[2]])),3),
                         Deviance =round( deviance(fit$mod.fit[[2]]),3),
                         AIC = round(AIC(fit$mod.fit[[2]]),3),
                         BIC = round(BIC(fit$mod.fit[[2]]),3),
                         McFadden_R2 = round(as.numeric(pscl::pR2(fit$mod.fit[[2]])['McFadden']),3),
                         `Train Misclassification Rate`=round(sum(abs(round(fit$mod.fit[[2]]$fitted.values) - fit$data.matrix$Y))/length(fit$data.matrix$Y),3),
                         `Test Misclassification Rate`= round(sum(abs(round(fit.pred[[2]]) - Ytest))/length(fit.pred[[2]] ),3))
      df2 <- rbind(df2a, df2b)
    }else{
      df2 <- data.frame(loglik = round(as.numeric(logLik(fit$mod.fit)),3),
                        Deviance =round( deviance(fit$mod.fit),3),
                        AIC = round(AIC(fit$mod.fit),3),
                        BIC = round(BIC(fit$mod.fit),3),
                        McFadden_R2 = round(as.numeric(pscl::pR2(fit$mod.fit)['McFadden']),3),
                        `Train Misclassification Rate`=round(sum(abs(round(fit$mod.fit$fitted.values) - fit$data.matrix$Y))/length(fit$data.matrix$Y),3),
                        `Test Misclassification Rate`= round(sum(abs(round(fit.pred$pred.mod) - Ytest))/length(fit.pred$pred.mod ),3))
    }
    df2
  }else if(fit$family=="poisson"){
    fit.pred <- predict(fit, newdata=temp.test[[1]], 
                        newdata2=temp.test[[2]], type="response")
    if(modeled_separately){
      df2a <- data.frame(loglik = round(as.numeric(logLik(fit$mod.fit[[1]])),3),
                         Deviance = round(deviance(fit$mod.fit[[1]]),3),
                         AIC = round(AIC(fit$mod.fit[[1]]),3),
                         BIC = round(BIC(fit$mod.fit[[1]]),3),
                         McFadden_R2 = round(as.numeric(pscl::pR2(fit$mod.fit[[1]])['McFadden']),3),
                         `Train MSE`=round(sum((round(fit$mod.fit[[1]]$fitted.values) - fit$data.matrix$Y)^2)/length(fit$data.matrix$Y),3),
                         `Test MSE`= round(sum((round(fit.pred[[1]]) - test)^2)/length(fit.pred[[1]] ),3))
      df2b <- data.frame(loglik = round(as.numeric(logLik(fit$mod.fit[[2]])),3),
                         Deviance = round(deviance(fit$mod.fit[[2]]),3),
                         AIC = round(AIC(fit$mod.fit[[2]]),3),
                         BIC = round(BIC(fit$mod.fit[[2]]),3),
                         McFadden_R2 = round(as.numeric(pscl::pR2(fit$mod.fit[[2]])['McFadden']),3),
                         `Train MSE`=round(sum((round(fit$mod.fit[[2]]$fitted.values) - fit$data.matrix$Y)^2)/length(fit$data.matrix$Y),3),
                         `Test MSE`= round(sum((round(fit.pred[[2]]) - Ytest)^2)/length(fit.pred[[2]] ),3))
      df2 <- rbind(df2a, df2b)
    }else{
      df2 <- data.frame(loglik = round(as.numeric(logLik(fit$mod.fit)),3),
                        Deviance = round(deviance(fit$mod.fit),3),
                        AIC = round(AIC(fit$mod.fit),3),
                        BIC = round(BIC(fit$mod.fit),3),
                        McFadden_R2 = round(as.numeric(pscl::pR2(fit$mod.fit)['McFadden']),3),
                        `Train MSE`=round(sum((round(fit$mod.fit$fitted.values) - fit$data.matrix$Y)^2)/length(fit$data.matrix$Y),3),
                        `Test MSE`= round(sum((round(fit.pred$pred.mod) - Ytest)^2)/length(fit.pred$pred.mod ),3))
    }
    df2
  }else if(fit$family=="survival"){
    if(modeled_separately){
      df2a <- data.frame(loglik =round( as.numeric(logLik(fit$mod.fit[[1]])),3),
                         Concordance = round(summary(fit$mod.fit[[1]])$concordance[1],3),
                         Concordance_SE = round(summary(fit$mod.fit[[1]])$concordance[2],3))
      df2b <- data.frame(loglik =round( as.numeric(logLik(fit$mod.fit[[2]])),3),
                         Concordance = round(summary(fit$mod.fit[[2]])$concordance[1],3),
                         Concordance_SE = round(summary(fit$mod.fit[[2]])$concordance[2],3))
      df2 <- rbind(df2a, df2b)
    }else{
      df2 <- data.frame(loglik =round( as.numeric(logLik(fit$mod.fit)),3),
                        Concordance = round(summary(fit$mod.fit)$concordance[1],3),
                        Concordance_SE = round(summary(fit$mod.fit)$concordance[2],3)
      )
    }
    df2
  }
}
