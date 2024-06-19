#' UMAP Plot
#'
#' @description Wrapper function to plot a UMAP of the results after
#' supervised filtering. See "umap" R package for more details on the
#' method.
#'
#' @param fit the output from the filter_supervised() function
#' @param useFilteredData Boolean on whether to plot UMAP on filtered or original data.
#' Default is filtered data.
#' @param usePrincipleComponents boolean on whether to apply PCA first 
#' @param plotIt boolean, whether to print the result or just return it (default = TRUE)
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
#' data.red=filter_supervised(Xdata, Y,  method="t.test", padjust=TRUE, thresh=0.05,
#'                  center=FALSE, scale=FALSE, Xtest=NULL)
#'
#' ##-----Plot Result via UMAP
#' umapPlot(data.red)
umapPlot <- function(fit,
                     useFilteredData=TRUE,
                     usePrincipleComponents = TRUE,
                     plotIt = TRUE){
  view.pca<-view.umap<-list()
  ##If supervised with binary/categorical outcome, plot 2 UMAPs
  #plot 1 = UMAP on filtered data
  #Seems Y has to be numeric
  
  plots = lapply(1:length(fit$X),
                 FUN = function(i){
                   
                   if(!useFilteredData){
                     mydata=fit$X_Original[[i]]
                   }else{
                     mydata=fit$X[[i]]
                   }
                   
                   view.umap[[i]] <- umap(cbind(mydata, fit$Y))
                   this_title =  paste("View", i, "- UMAP on filtered data")
                   if (usePrincipleComponents){
                     t1 <- prcomp(mydata, rank.=10)
                     view.pca[[i]] <- umap(cbind(t1$x, fit$Y))
                     view.umap[[i]] = view.pca[[i]]
                     this_title =  paste("View", i, "- UMAP on PCA on filtered data")
                   }
                   
                  
                   my.cols <- c(rgb(0,0,0,alpha=0.8), #grey
                                rgb(1,0,0,alpha=0.8), #red
                                rgb(0,0,1,alpha=0.8), #blue
                                rgb(0,1,0,alpha=0.8)) #green
                   ggplot2::ggplot(as.data.frame(view.umap[[i]]$layout), 
                                        aes(view.umap[[i]]$layout[,1], view.umap[[i]]$layout[,2], 
                                            col=factor(fit$Y))) +
                     geom_point() + theme_bw() + xlab("UMAP Component 1") +
                     ylab("UMAP Component 2") +
                     ggtitle(this_title) +
                     scale_colour_manual(values=my.cols) +
                     theme(axis.title = element_text(face="bold"))+
                     theme(axis.text = element_text(face="bold"))+
                     guides(color = guide_legend(title = "Outcome"))
                 })
  if (plotIt){
    gridExtra::grid.arrange(grobs = plots, nrow = 1)
    return()
  }else{
    return(plots)
  }
}

#' Volcano Plot
#'
#' @description Wrapper function for volcano plots of the results after
#' supervised filtering.
#'
#' @param fit the output from the filter_supervised() function
#' @param plotIt boolean, whether to print the result (TRUE) or just return it
#'
#' @return A graph of the volcano plot
#' @import ggplot2
#' @importFrom dplyr "%>%"
#' @importFrom grDevices rgb
#'
#' @export
#'
#' @examples
#' ##---- read in data
#' data(COVIDData)
#'
#' #make omics data numeric
#' Proteomics= apply(as.matrix(COVIDData[[1]]), 2, as.numeric)
#' RNASeq= apply(as.matrix(COVIDData[[2]]), 2, as.numeric)
#' Clinical= COVIDData[[3]]
#' X=list(Proteomics, RNASeq)
#' Y=Clinical$DiseaseStatus.Indicator
#'
#' data.red=filter_supervised(X, Y,  method="t.test", padjust=TRUE,adjmethod="BH",
#' thresh=0.05,center=TRUE, scale=TRUE, Xtest=NULL)
#'
#' ##-----Volcano Plot of Result
#' volcanoPlot(data.red)
volcanoPlot <- function(fit, plotIt = TRUE){
  results = fit[["pval.mat"]]
  results$VariableName = sub("\\;.*", "", row.names(results))
  xlab = dplyr::case_when(
    fit$method == "logistic" ~ ("Log Odds Ratio"),
    fit$method == "linear" ~ ("Means"),
    fit$method == "t.test" ~ ("Mean Difference"),
    fit$method == "kw" ~ ("KW Test Statistic")
  )
  ylab = dplyr::case_when(
    fit$padjust ~ "-log[10]*padj",
    !fit$padjust ~ "-log[10]*p"
  )
  
  plots = lapply(unique(results$View),
                 FUN = function(view){
                   sub_results = results %>%
                     dplyr::filter(View == view)
                   labeled_names = sub_results %>% 
                     dplyr::filter(Keep) %>%
                     dplyr::group_by(Coef < 0) %>%
                     dplyr::arrange(Pval, Coef) %>%
                     dplyr::slice_head(n = 15) %>%
                     dplyr::pull(VariableName)
                   
                   sub_results %>%
                     dplyr::mutate(Name = ifelse(VariableName %in% labeled_names,
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
