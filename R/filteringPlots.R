#' UMAP Plot
#'
#' @description Wrapper function to plot a UMAP of the results after
#' supervised filtering. See "umap" R package for more details on the
#' method.
#'
#' @param filtering_fit the output from the filter.supervised() function
#' @param filteredData Boolean on whether to plot UMAP on filtered or original data.
#' Default is filtered data.
#'
#' @return A graph of the UMAP
#'
#' @import umap
#' @import ggplot2
#' @importFrom stats prcomp
#'
#' @export
#'
#' @examples
#' ##---- read in data
#' data(sidaData)
#'
#' Xdata=sidaData[[1]]
#' Y=sidaData[[2]]
#'
#' data.red=filter.supervised(Xdata, Y,  method="t.test", padjust=TRUE, thresh=0.05,
#'                  center=FALSE, scale=FALSE, Xtest=NULL)
#'
#' ##-----Plot Result via UMAP
#' umapPlot(data.red)
umapPlot <- function(object,filteredData=TRUE){
  view.pca<-view.umap<-list()
  ##If supervised with binary/categorical outcome, plot 2 UMAPs
  #plot 1 = UMAP on filtered data
  #Seems Y has to be numeric
  
  for(i in 1:length(object$X)){
    
    if(filteredData==FALSE){
      mydata=object$X_Original[[i]]
    }else{
      mydata=object$X[[i]]
    }
    
    view.umap[[i]] <- umap(cbind(mydata, object$Y))
    
    #PCA
    t1 <- prcomp(mydata, rank.=10)
    view.pca[[i]] <- umap(cbind(t1$x, object$Y))
    
    my.cols <- c(rgb(0,0,0,alpha=0.8), #grey
                 rgb(1,0,0,alpha=0.8), #red
                 rgb(0,0,1,alpha=0.8), #blue
                 rgb(0,1,0,alpha=0.8)) #green
    print(
      ggplot2::ggplot(as.data.frame(view.umap[[i]]$layout), aes(view.umap[[i]]$layout[,1], view.umap[[i]]$layout[,2], col=factor(object$Y))) +
        geom_point() + theme_bw() + xlab("UMAP Component 1") +
        ylab("UMAP Component 2") +
        ggtitle( paste("View", i, "- UMAP on filtered data")) +
        scale_colour_manual(values=my.cols) +
        theme(axis.title = element_text(face="bold"))+
        theme(axis.text = element_text(face="bold"))+
        guides(color = guide_legend(title = "Outcome"))
    )
    print(
      ggplot2::ggplot(as.data.frame(view.pca[[i]]$layout), aes(view.pca[[i]]$layout[,1], view.pca[[i]]$layout[,2], col=factor(object$Y))) +
        geom_point() + theme_bw() + xlab("UMAP Component 1") +
        ylab("UMAP Component 2") +
        ggtitle(paste("View", i, "- UMAP on PCA filtered data")) +
        theme(axis.title = element_text(face="bold"))+
        theme(axis.text = element_text(face="bold"))+
        scale_colour_manual(values=my.cols) +
        guides(color = guide_legend(title = "Outcome"))
    )
    Sys.sleep(5)
    
  }
  
}

#' Volcano Plot
#'
#' @description Wrapper function for volcano plots of the results after
#' supervised filtering.
#'
#' @param object the output from the filter.supervised() function
#'
#' @return A graph of the volcano plot
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' ##---- read in data
#' data(COVID)
#'
#' #make omics data numeric
#' Proteomics= apply(as.matrix(COVIDData[[1]]), 2, as.numeric)
#' RNASeq= apply(as.matrix(COVIDData[[2]]), 2, as.numeric)
#' Clinical= COVIDData[[3]]
#' X=list(Proteomics, RNASeq)
#' Y=Clinical$DiseaseStatus.Indicator
#'
#' data.red=filter.supervised(X, Y,  method="t.test", padjust=TRUE,adjmethod="BH",
#' thresh=0.05,center=TRUE, scale=TRUE, Xtest=NULL)
#'
#' ##-----Volcano Plot of Result
#' volcanoPlot(data.red)
volcanoPlot <- function(fit, plotIt = TRUE){
  results = fit[["pval.mat"]]
  results$VariableName = sub("\\;.*", "", row.names(results))
  xlab = case_when(
    fit$method == "logistic" ~ ("Log Odds Ratio"),
    fit$method == "linear" ~ ("Means"),
    fit$method == "t.test" ~ ("Mean Difference"),
    fit$method == "kw" ~ ("KW Test Statistic")
  )
  ylab = case_when(
    fit$padjust ~ "-log[10]*padj",
    !fit$padjust ~ "-log[10]*p"
  )
  
  plots = lapply(unique(results$View),
                 FUN = function(view){
                   sub_results = results %>%
                     dplyr::filter(View == view)
                   labeled_names = sub_results %>% 
                     dplyr::filter(Keep) %>%
                     group_by(Coef < 0) %>%
                     arrange(Pval, Coef) %>%
                     slice_head(n = 15) %>%
                     pull(VariableName)
                   
                   sub_results %>%
                     mutate(Name = ifelse(VariableName %in% labeled_names,
                                          VariableName, NA)) %>%
                     ggplot(aes(x = Coef, y = -log10(Pval),
                                color = Keep,
                                label = Name))+
                     geom_hline(yintercept = -log10(0.05),
                                linetype = "dashed")+
                     geom_point()+
                     geom_text(color = "black",
                               check_overlap = TRUE, size=4,
                               hjust=0, nudge_x = -0.25, fontface="bold")+
                     scale_color_manual(values = c("#999999", "#0072B2", "#56B4E9")) +
                     guides(colour = guide_legend(override.aes = list(size=5))) +
                     theme_bw() +
                     theme(axis.title = element_text(face="bold"))+
                     theme(axis.text = element_text(face="bold"))+
                     guides(colour = guide_legend(override.aes = list(size=5))) +
                     xlab(xlab)+
                     ylab(
                       if(fit$padjust){
                         expression("-log"[10]*"padj")
                       }else{
                         expression("-log"[10]*"p")
                       }
                     )+
                     ggtitle(paste0("Volcano plot for view ", view))
                 }
  )
  
  if (plotIt){
    gridExtra::grid.arrange(grobs = plots, nrow = 1)
  }else{
    return(plots)
  }
  
}
