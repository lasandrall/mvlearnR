#modified to include covariates for adjustment in models
filter.supervised.covadjust <- function(X, Y, method="linear", padjust=F,adjmethod="BH", thresh=0.05,
                              center=F, scale=F, standardize=F, log2TransForm=F, Xtest=NULL, CovAdjust=NULL){
  # if(is.null(adjmethod)){
  #   adjmethod="BH"
  # }
  # XOrig=X
  # XtestOrig=Xtest
  #
  # for(i in 1:length(X)){
  #
  #   #convert data to numeric
  #   X[[i]]=apply(as.matrix(X[[i]]), 2, as.numeric)
  #
  #   if(length(XtestOrig)>0){
  #     Xtest[[i]] <- apply(as.matrix(Xtest[[i]]), 2, as.numeric)
  #   }
  #
  #   if(log2TransForm){
  #     X[[i]] <- apply(X[[i]], 2, log2)
  #     if(length(XtestOrig)>0){
  #       Xtest[[i]] <- apply(Xtest[[i]], 2, log2)
  #     }
  #   }
  #
  #   if(standardize){
  #     mean.vec <- apply(X[[i]], 2, mean, na.rm=TRUE)
  #     var.vec <- apply(X[[i]], 2, function(x) sqrt(var(x,na.rm=TRUE)))
  #     X[[i]] <- t( (t(X[[i]]) - mean.vec) /var.vec)
  #     if(length(XtestOrig)>0){
  #       Xtest[[i]] <- t( (t(Xtest[[i]])-mean.vec)/var.vec)
  #     }
  #     #X[[i]]=apply(X[[i]], 2, function(x) scale(x))
  #   }
  #
  #
  #   if(center){
  #     mean.vec <- apply(X[[i]], 2, mean,na.rm=TRUE)
  #     X[[i]] <- t(t(X[[i]]) - mean.vec )
  #     if(length(XtestOrig)>0){
  #       Xtest[[i]] <- t(t(Xtest[[i]]) - mean.vec)
  #     }
  #   }
  #   if(scale){
  #     var.vec <- apply(X[[i]], 2, function(x) sqrt(var(x,na.rm=TRUE)))
  #     X[[i]] <- t(t(X[[i]])/var.vec)
  #     if(length(XtestOrig)>0){
  #       Xtest[[i]] <- t(t(Xtest[[i]])/var.vec)
  #     }
  #   }
  #
  #   # if(standardize){
  #   #   mean.vec <- apply(XOrig[[i]], 2, mean)
  #   #   var.vec <- apply(XOrig[[i]], 2, function(x) sqrt(var(x)))
  #   #   X[[i]] <- t( (t(XOrig[[i]]) - mean.vec) /var.vec)
  #   #   # if(length(XtestOrig)>0){
  #   #   #   Xtest[[i]] <- t(t(Xtest[[i]])/var.vec)
  #   #   # }
  #   #   #X[[i]]=apply(X[[i]], 2, function(x) scale(x))
  #   # }
  #
  #
  #
  # }
  #
  # coef.mat <- pval.mat <- pval.adj.mat <- red.mat <- mean.mat <- X.red <- Xtest.red <- list()
  # for(i in 1:length(X)){
  #   if(!is.null(CovAdjust)){
  #     mydata=cbind.data.frame(X[[i]],CovAdjust)
  #   }else{
  #     CovAdjust=1
  #   }
  #   if(method == "linear"){
  #     coef.mat[[i]] <- apply(X[[i]], 2, function(x)
  #       lm(as.numeric(Y) ~ as.numeric(x) + CovAdjust,data=mydata)$coefficients[2])
  #     pval.mat[[i]] <- apply(X[[i]], 2, function(x)
  #       summary(lm(as.numeric(Y) ~ as.numeric(x)) + CovAdjust,data=mydata)$coefficients[2,4])
  #   }else if(method == "logistic"){
  #     coef.mat[[i]] <- apply(X[[i]], 2, function(x)
  #       glm(factor(Y) ~ as.numeric(x) , family = binomial())$coefficients[2])
  #     pval.mat[[i]] <- apply(X[[i]], 2, function(x)
  #       summary(glm(factor(Y) ~ as.numeric(x) ,family = binomial()))$coefficients[2,4])
  #   }else if(method == "t.test"){
  #     coef.mat[[i]] <- apply(X[[i]], 2, function(x)
  #       t.test(as.numeric(x) ~ factor(Y))$estimate[1]-t.test(as.numeric(x) ~ factor(Y))$estimate[2])
  #     pval.mat[[i]] <- apply(X[[i]], 2, function(x)
  #       t.test(as.numeric(x) ~ factor(Y))$p.value)
  #     #coef.mat[[i]] <- pval.mat[[i]]*NA
  #   }else if(method == "kw"){
  #     coef.mat[[i]] <- apply(X[[i]], 2, function(x)
  #       kruskal.test(as.numeric(x) ~ factor(Y))$statistic)
  #     pval.mat[[i]] <- apply(X[[i]], 2, function(x)
  #       kruskal.test(as.numeric(x) ~ factor(Y))$p.value)
  #   }else{
  #     cat("Warning: Method does not exist")
  #     quit()
  #   }
  #
  #
  #
  #   if(padjust==TRUE){
  #     pval.adj.mat[[i]] <-p.adjust(pval.mat[[i]], method=adjmethod)
  #     red.mat[[i]] <- 1:ncol(X[[i]]) %in% which(pval.adj.mat[[i]] < thresh)
  #     X.red[[i]] <- X[[i]][,red.mat[[i]]]
  #     pval.mat[[i]]=pval.adj.mat[[i]]
  #   }else{
  #     red.mat[[i]] <- 1:ncol(X[[i]]) %in% which(pval.mat[[i]] < thresh)
  #     X.red[[i]] <- X[[i]][,red.mat[[i]]]
  #   }
  #
  #   if(length(Xtest)>0){
  #     Xtest.red[[i]] <- Xtest[[i]][,red.mat[[i]]]
  #   }
  # }
  #
  # # #label for ttest
  # # if(method=="t.test"){
  # #   che=as.data.frame(t.test(as.numeric(X[[1]][,1]) ~ factor(Y))$estimate)
  # #   mean.diff.label=rownames(che)[1]-rownames(che)[2]
  # # }
  #
  # temp <- data.frame(Coef = NULL,
  #                    Pval = NULL,
  #                    Keep = NULL,
  #                    View = NULL)
  # for(i in 1:length(X)){
  #   t1 <- data.frame(Coef = coef.mat[[i]],
  #                    Pval = pval.mat[[i]],
  #                    Keep = red.mat[[i]],
  #                    View = i)
  #   temp <- rbind(temp,t1)
  #
  #
  #   # if(method=="t.test"){
  #   #   colnames(temp)[1]="Mean Difference"
  #   # }else if(method=="linear"){
  #   #   colnames(temp)[1]="Coef"
  #   # }else if(method=="logistic"){
  #   #   colnames(temp)[1]="Log ORs"
  #   # }else if(method=="kw"){
  #   #   colnames(temp)[1]="Coef"
  #   # }
  #
  #
  # }
  #
  # #print results for top 10 significant variables by views
  # for(i in 1:length(X)){
  #   mydata2=temp[temp$View==i,]
  #   mydata2.Sorted= mydata2[order(mydata2[,2]),]
  #   mydatathresh=mydata2.Sorted[mydata2.Sorted[,2]< thresh, ]
  #   if(length(mydatathresh)==0){
  #     print(paste0("No variable is significnat for View ",i))
  #     print(mydatathresh)
  #   }else{
  #     print(paste0("Printing top 10 results for  significnat variables for View ",i))
  #     print(mydatathresh[1:10,])
  #   }
  # }
  #
  # #i should add padjust to output
  # result <- list(X=X.red, Y=Y,
  #                Xtest=Xtest.red,
  #                method=method,
  #                pval.mat=temp,
  #                significant.thresh=thresh,
  #                adjmethod=adjmethod,
  #                X_Original=XOrig,
  #                Xtest_Original=XtestOrig,
  #                Center=center,
  #                Scale=scale,
  #                Log2Transform=log2TransForm,
  #                Standardize=standardize)
  # return(result)

  # filter.supervised <- function(X, Y, method="linear", padjust=FALSE,adjmethod="BH", thresh=0.05,
  #                               center=FALSE, scale=FALSE, standardize=FALSE, log2TransForm=FALSE, Xtest=NULL){

  #allows to make thresh a vector of length X data
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




