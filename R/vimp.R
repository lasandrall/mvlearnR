vimp=function(fit,method){

  graph.dat <- list()

  if(method=="SIDA"){
    k <- length(fit$hatalpha)
    hatalpha=fit$hatalpha
    mysum=matrix(NA,nrow=k,ncol=1)
    for(i in 1:k){
      mysum[i,1]=sum(fit$hatalpha[[i]]!=0)
    }
    topk=min(20,min(mysum))
    hatalpha.temp=list()
    for(i in 1:k){
      # hatalpha.temp[[i]]=rowMeans(abs(hatalpha[[i]]))
      # col1<-order(abs(hatalpha.temp[[i]]), decreasing = T)[1:topk]
      # mycolnames=sub("\\;.*", "", colnames(fit$InputData[[i]]))
      # graph.dat[[i]] <- data.frame(mycolnames[col1])
      # #graph.dat[[i]] <- data.frame(name=order(abs(fit$hatalpha[[i]]), decreasing = T)[1:20])
      # graph.dat[[i]]$val <- hatalpha.temp[[i]][col1]
      if(dim(hatalpha[[i]])[2] > 1){
        hatalpha.temp[[i]]=rowMeans(abs(hatalpha[[i]]))
        col1<-order(abs(hatalpha.temp[[i]]), decreasing = T)[1:topk]
        mycolnames=sub("\\;.*", "", colnames(as.data.frame(fit$InputData[[i]])))
        graph.dat[[i]] <- data.frame(mycolnames[col1])
        #graph.dat[[i]] <- data.frame(name=order(abs(fit$hatalpha[[i]]), decreasing = T)[1:20])
        graph.dat[[i]]$val <- hatalpha.temp[[i]][col1]
        #graph.dat[[i]]$abs_val <- rowMeans(abs(graph.dat[[i]]$val))
        graph.dat[[i]]$abs_val <- abs(graph.dat[[i]]$val)
      }else if(dim(hatalpha[[i]])[2] ==1){
        col1<-order(abs(hatalpha[[i]]), decreasing = T)[1:topk]
        mycolnames=sub("\\;.*", "", colnames(as.data.frame(fit$InputData[[i]])))
        graph.dat[[i]] <- data.frame(mycolnames[col1])
        graph.dat[[i]]$val <- fit$hatalpha[[i]][col1]
        graph.dat[[i]]$abs_val <- abs(graph.dat[[i]]$val)
      }
      #graph.dat[[i]]$abs_val <- rowMeans(abs(graph.dat[[i]]$val))
      graph.dat[[i]]$nri <- graph.dat[[i]]$abs_val/graph.dat[[i]]$abs_val[1]
      colnames(graph.dat[[i]])[1]="name"
    }
    p <- list()
    for(i in 1:k){
      graph.dat[[i]]$name <- factor(graph.dat[[i]]$name,
                                    levels=graph.dat[[i]]$name[topk:1])
      print( ggplot(graph.dat[[i]], aes(name, graph.dat[[i]]$nri)) +
        geom_bar(stat="identity", fill="#7A0019") + theme_bw() + theme(legend.position="none") +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
        coord_flip() + ylab("Normalized relative importance") + xlab('') +
        theme(axis.title = element_text(face="bold"))+
        theme(axis.text = element_text(face="bold"))+
        ggtitle(paste0("Variable importance for View ", i))
      )
      Sys.sleep(2)
    }
    }else if(method=="SELP"){

    hatalpha=fit$hatalpha
    hatbeta=fit$hatbeta
    mysum=matrix(NA,nrow=2,ncol=1)
    mysum[1,1]=sum(fit$hatalpha[,1]!=0)
    mysum[2,1]=sum(fit$hatbeta[,1]!=0)
    topk=min(20,min(mysum))
      for(i in 1){
        #graph.dat[[1]]=data.frame(name=order(abs(hatalpha[,i]), decreasing = T)[1:20])
        col1<-order(abs(hatalpha[,i]), decreasing = T)[1:topk]
        #col1<-(name=order(abs(hatalpha[,i]), decreasing = T)[1:20])
        mycolnames=sub("\\;.*", "", colnames(as.data.frame(fit$InputData[[i]])))
        graph.dat[[1]] <- data.frame(mycolnames[col1])
        #print(graph.dat[[1]][1:5,])
        #print(hatalpha[,i][graph.dat[[i]]$name])
        graph.dat[[1]]$val <- hatalpha[,i][col1]
        graph.dat[[1]]$abs_val <- abs(graph.dat[[i]]$val)
        graph.dat[[1]]$nri <- graph.dat[[i]]$abs_val/graph.dat[[i]]$abs_val[1]
        colnames(graph.dat[[1]])[1]="name"
      }
      for(i in 1){
        #col1=data.frame(name=order(abs(hatbeta[,i]), decreasing = T)[1:20])
        col1=order(abs(hatbeta[,i]), decreasing = T)[1:topk]
        mycolnames=sub("\\;.*", "", colnames(as.data.frame(fit$InputData[[2]])))
        graph.dat[[2]] <- data.frame(mycolnames[col1])
        graph.dat[[2]]$val <- hatbeta[,i][col1]
        graph.dat[[2]]$abs_val <- abs(graph.dat[[2]]$val)
        graph.dat[[2]]$nri <- graph.dat[[2]]$abs_val/graph.dat[[2]]$abs_val[1]
        colnames(graph.dat[[2]])[1]="name"
      }
      p <- list()

      for(i in 1:2){
        graph.dat[[i]]$name <- factor(graph.dat[[i]]$name,
                                      levels=graph.dat[[i]]$name[topk:1])
        # graph.dat[[i]]$name <- factor(colnames(graph.dat[[i]]),
        #                               levels=colnames(graph.dat[[i]])[20:1])
        print(
        ggplot(graph.dat[[i]], aes(name, graph.dat[[i]]$nri)) +
          geom_bar(stat="identity", fill="#7A0019") + theme_bw() + theme(legend.position="none") +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
          coord_flip() + ylab("Normalized relative importance") + xlab('') +
          theme(axis.title = element_text(face="bold"))+
          theme(axis.text = element_text(face="bold"))+
          ggtitle(paste0("Variable importance for View ", i, " CCA Vector 1"))
         )

        Sys.sleep(5)
      }
  #
    }
  #return(invisible(NULL))
}


