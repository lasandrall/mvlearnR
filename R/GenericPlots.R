
#' @title Variable Importance Plot
#'
#' @description Wrapper function to visualize loadings for variables selected
#' by SIDA, SIDANet, and SELPCCA methods.
#'
#' @param object the output from SIDA, SIDANet, and SELPCCA methods
#'
#' @return A graph of the absolute loadings for variables selected. The variables
#' are normalized to the variable with the largest weight.
#'
#' @export
#'
#' @examples
#' ##---- load SIDA data
#' data("sidaData")
#' Xdata <- sidaData[[1]]
#' Y <- sidaData[[2]]
#' Xtestdata <- sidaData[[3]]
#' Ytest <- sidaData[[4]]
#' ##---- call cross validation
#' mycv=cvSIDA(Xdata,Y,withCov=FALSE,plotIt=FALSE, Xtestdata=Xtestdata,Ytest=Ytest)
#' ##----  Obtain variable importance plot
#' VarImportancePlot(mycv)

#############################################################################################################
# Authors:
#   Sandra E. Safo, Unviersity of Minnesota
#   Elise Palzer, Unviersity of Minnesota
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


VarImportancePlot <- function(object){
  if(class(object)=="SIDA" | class(object)=="SIDANet"){
    method="SIDA"
    myobject=object

    hatalpha=object$hatalpha
    k <- length(object$hatalpha)
    mysum=matrix(NA,nrow=k,ncol=1)

    hatalpha.temp=list()
    for(i in 1:k){
      if(dim(hatalpha[[i]])[2] == 1){
         mysum[i,1]=sum(hatalpha[[i]]!=0)
      }else if (dim(hatalpha[[i]])[2] >1){
          hatalpha.temp[[i]]=rowSums(abs(hatalpha[[i]]))
          mysum[i,1]=sum(hatalpha.temp[[i]]!=0)
      }
    }
    topk=min(20,min(mysum))

    dfView =list()
    hatalpha.temp=list()
    df=list()
    for(i in 1:k){
      if(dim(hatalpha[[i]])[2] > 1){
        hatalpha.temp[[i]]=rowMeans(abs(hatalpha[[i]]))
        col1<-order(abs(hatalpha.temp[[i]]), decreasing = T)
        mycolnames=sub("\\;.*", "", colnames(as.data.frame(object$InputData[[i]])))
        col1name <- data.frame(mycolnames[col1])
        col3 <- hatalpha.temp[[i]][col1]
        col4 <- col3/col3[1]
        view.no=data.frame(rep(i,dim(col1name)[1]))
        df[[i]] <- cbind.data.frame(view.no, col1name,  data.frame(col3), data.frame(col4))
        df2=df[[i]]
        # df <- rbind.data.frame(df, cbind.data.frame(i, col1name,   col3, col4))
        # df2 <- as.data.frame(df)
        colnames(df2) <- c("View", "Variable Name",
                            "Mean Loading",
                           "Normalized Relative Importance")
        dfView[[i]] <- cbind.data.frame(df2[,1:2],round(df2[,c(3:4)],3))
        print(dfView[[i]][1:topk,])
      }else if(dim(hatalpha[[i]])[2] ==1){
        col1<-order(abs(hatalpha[[i]]), decreasing = T)
        mycolnames=sub("\\;.*", "", colnames(as.data.frame(object$InputData[[i]])))
        col1name <- data.frame(mycolnames[col1])
        col2=hatalpha[[i]][col1]
        col3 <- abs(hatalpha[[i]][col1])
        col4 <- col3/col3[1]
        view.no=data.frame(rep(i,dim(col1name)[1]))
        df[[i]] <- cbind.data.frame(view.no, col1name, data.frame(col2), data.frame(col3), data.frame(col4))
        df2=df[[i]]
        colnames(df2) <- c("View", "Variable Name",
                           "Loading", "Absolute Loading",
                           "Normalized Relative Importance")
        dfView[[i]] <- cbind.data.frame(df2[,1:2],round(df2[,c(3:5)],3))
        print(dfView[[i]][1:topk,])
      }

    }
      # #graph.dat[[i]]$abs_val <- rowMeans(abs(graph.dat[[i]]$val))
      # graph.dat[[i]]$nri <- graph.dat[[i]]$abs_val/graph.dat[[i]]$abs_val[1]
      # colnames(graph.dat[[i]])[1]="name"
      # df <- NULL
      # col1 <- order(abs(object$hatalpha[[i]]), decreasing = T)
      # col2 <- object$hatalpha[[i]][col1]
      # col1name=colnames(object$InputData[[i]])[col1]
      # col3 <- abs(col2)
      # col4 <- col3/col3[1]
      # #df <- rbind(df, cbind(i, col1, col2, col3, col4))

      # df <- rbind.data.frame(df, cbind.data.frame(i, col1name, col2, col3, col4))
      # df2 <- as.data.frame(df)
      # colnames(df2) <- c("View", "Variable Name",
      #                    "Loading", "Absolute Loading",
      #                    "Normalized Relative Importance")
      # dfView[[i]] <- cbind.data.frame(df2[,1:2],round(df2[,c(3:5)],3))
      # print(dfView[[i]][topk:1,])
      # print(length(df))
      # df2 <- do.call(rbind.data.frame,df)
      # print(df2[1:200,])
      # colnames(df2) <- c("View", "Variable Name",
      #                    "Loading", "Absolute Loading",
      #                    "Normalized Relative Importance")
      # dfView <- cbind.data.frame(df2[,1:2],round(df2[,c(3:5)],3))
      # print(dfView[1:topk,])



    # df <- as.data.frame(df)
    # names(df) <- c("View", "Variable Name",
    #                "Loading", "Absolute Loading",
    #                "Normalized Relative Importance")
    # df <- round(df,3)
    # df

    # df2 <- as.data.frame(df)
    # colnames(df2) <- c("View", "Variable Name",
    #                    "Loading", "Absolute Loading",
    #                    "Normalized Relative Importance")
    # dfView <- cbind.data.frame(df2[,1:2],round(df2[,c(3:5)],3))
    # print(dfView[1:20,])

    #print(df)



    #get graphs
    print(vimp(myobject,method))

    result=list(dfView)

  }else if(class(object)=="SELP-Predict"|class(object)=="SELPCCA"){
    method="SELP"
    if(object$method=="cvselpscca"|object$method=="multiplescca"){
      myobject=object
      k=ncol(object$hatalpha)
      #col1 <- order(abs(object$hatalpha[,i]), decreasing = T)
      #col2 <- object$hatalpha[,i][col1]
    }else if(object$method=="selpscca.pred"){
      myobject=object$selp.fit
      k <- ncol(object$selp.fit$hatalpha)
      #col1 <- order(abs(object$selp.fit$hatalpha[,i]), decreasing = T)
      #col2 <- object$selp.fit$hatalpha[,i][col1]
  }
    df <- NULL
    for(i in 1:k){

      if(object$method=="cvselpscca"|object$method=="multiplescca"){
        col1 <- order(abs(object$hatalpha[,i]), decreasing = T)
        col2 <- object$hatalpha[,i][col1]
      }else if(object$method=="selpscca.pred"){
        col1 <- order(abs(object$selp.fit$hatalpha[,i]), decreasing = T)
        col2 <- object$selp.fit$hatalpha[,i][col1]
      }

      #col1 <- order(abs(object$selp.fit$hatalpha[,i]), decreasing = T)
      col1name=colnames(as.data.frame(object$InputData[[1]]))[col1]
      #print(col1name[1:5])
      #col2 <- object$selp.fit$hatalpha[,i][col1]
      col3 <- abs(col2)
      col4 <- col3/col3[1]
      df <- rbind.data.frame(df, cbind.data.frame(i, col1name, col2, col3, col4))


    }
    df2 <- as.data.frame(df)
        colnames(df2) <- c("View 1 Canonical Correlation Vector", "Variable Name",
                   "Loading", "Absolute Loading",
                   "Normalized Relative Importance")
    dfView1 <- cbind.data.frame(df2[,1:2],round(df2[,c(3:5)],3))
    for(j in 1:k){
      ithcomp=dfView1[dfView1[,1]==j,]
      print(ithcomp[1:10,])
    }

    #print(dfView1[1:20,])

    if(object$method=="cvselpscca"|object$method=="multiplescca"){
      k=ncol(object$hatbeta)
      #col1 <- order(abs(object$hatbeta[,i]), decreasing = T)
      #col2 <- object$hatbeta[,i][col1]
    }else if(object$method=="selpscca.pred"){
      k <- ncol(object$selp.fit$hatbeta)
      #col1 <- order(abs(object$selp.fit$hatbeta[,i]), decreasing = T)
      #col2 <- object$selp.fit$hatbeta[,i][col1]
    }
    df <- NULL
    for(i in 1:k){

      if(object$method=="cvselpscca"|object$method=="multiplescca"){
        #k=object$hatbeta
        col1 <- order(abs(object$hatbeta[,i]), decreasing = T)
        col2 <- object$hatbeta[,i][col1]
      }else if(object$method=="selpscca.pred"){
        #k <- ncol(object$selp.fit$hatbeta)
        col1 <- order(abs(object$selp.fit$hatbeta[,i]), decreasing = T)
        col2 <- object$selp.fit$hatbeta[,i][col1]
      }

      #col1 <- order(abs(object$selp.fit$hatbeta[,i]), decreasing = T)
      col1name=colnames(as.data.frame(object$InputData[[2]]))[col1]
      #print(col1name[1:5])
      #col2 <- object$selp.fit$hatbeta[,i][col1]
      col3 <- abs(col2)
      col4 <- col3/col3[1]
      df <- rbind.data.frame(df, cbind.data.frame(i, col1name, col2, col3, col4))
    }
    df2 <- as.data.frame(df)
    colnames(df2) <- c("View 2 Canonical Vector", "Variable Name",
                       "Loading", "Absolute Loading",
                       "Normalized Relative Importance")
    dfView2 <- cbind.data.frame(df2[,1:2],round(df2[,c(3:5)],3))

    for(j in 1:k){
      ithcomp=dfView2[dfView2[,1]==j,]
      print(ithcomp[1:10,])
    }
    #print(dfView2[1:20,])


    #get graphs
    print(vimp(myobject,method))

    result=list(dfView1,dfView2)
    return(invisible(result))
  }
}


#' @title Network visualization of selected variables from integrative analysis methods
#'
#' @description Wrapper function to visualize  graph of similarity matrix for selected variables.
#' We estimate pairwise similarity matrix using low-dimensional representations of our
#' sparse integrative analysis methods (selpcca, sida, sidanet). We follow ideas in González et al. 2012
#' to create bipartite graph (bigraph) where variables or nodes from one view are connected to variables
#' or nodes from another view. We construct the bigraph from a pairwise similarity matrix obtained
#' from the outputs of our integrative analysis methods.  We estimate the similarity
#' score between a pair of selected variables from two views  by calculating the inner product of each selected variable
#' and the sum of canonical variates (for SELPCCA) or discriminant vectors (for SIDA, SIDANet) for the pairs of views.
#' As noted in González et al. 2012, the entries in the similarity matrix is a robust approximation of the Pearson correlation between pairs of variables
#' and the two views under consideration. This network graph
#' has potential to shed light on the complex associations between pairs of views.
#' @param object the output from SIDA, SIDANet, and SELPCCA methods
#' @param cutoff a vector containing  numeric values between 0 and 1 of similarity cutoff to use when
#' generating graphs. Length of vector is number of pairwise data combinations. Variable pairs with high
#' similarity measure may be of interest. The relevance of the associations can be explored by changing the cutoff.
#' This can also be used to reduce the size of the graph, for dense network. Default is 0.5 meaning that graph will only be generated
#' for variable pairs with similarity value greater than 0.5 for each data pair.
#' @param color.node vector of length two, specifying the colors of nodes for pairs of views. Defaults
#' to white and yellow.
#' @param lty.edge character vector of length 2, specifying the line type for edges with positive and negative weights, respectively.
#' Can be one of "solid", "dashed", "dotted", "dotdash", "longdash" and "twodash".
#' See igraph package for more details. Defaults to c("solid", "dashed"), where positive weights are solid lines, and negative weights are dashed lines.
#' @param show.edge.labels boolen indicating whether or not to show weights as edge labels.
#' @param vertex.frame.color a character string of color to use as frame for nodes. Defaults to "red".
#' @param layout.fun a function, specifying how the vertices will be placed on the graph. Refer to igraph package using help(layout) for more details.
#' Default is layout.fruchterman.reingold.
#' @param show.color.key boolen indicating whether or not to show color key on plot. Defaults to TRUE.
#' Positive weights or similarity values (correlations) are indicated as red and negative values are indicated as green.
#' @param save should the plot be saved? If so, choose one of these options: "jpeg", "tiff", "png" or "pdf"
#' @param name.save character string for the name of the file to be saved.
#' @details The function will return D R objects, where D is the number of views. To see the results, use the “$" operator.
#' @return A network graph for variables selected. Each list will contain similarity matrix, cutoff used, and indices of pairings.
#' \item{networkGraph}{a graph object for each pair of views (if more than two views) that can be interrogated in cytoscape}
#' \item{SimilarityMatrix}{ the similarity matrix used for generating the network for each pair of views}
#' \item{cutoff}{the cutoff used when generating network}
#' \item{pairs}{The pairs of views for which network(s) were generated}
#' @references
#' Elise Palzer and Sandra E. Safo 2023. Submitted
#'González I., Lê Cao K-A., Davis, M.J. and Déjean, S. (2012).
#'Visualising associations between paired omics data sets. J. Data Mining 5:19.
#'
#' @import igraph
#' @import grDevices
#' @import methods
#' @import graphics
#' @importFrom grDevices dev.cur dev.off jpeg tiff pdf
#' @importFrom stats cor
#' @importFrom utils combn
#' @importFrom grDevices dev.new colorRamp colorRampPalette
#' @importFrom methods is
#' @importFrom graphics strwidth strheight layout image box axis title
#' @export
#'
#' @examples
#' ##---- load SIDA data
#' data("sidaData")
#' Xdata <- sidaData[[1]]
#' Y <- sidaData[[2]]
#' Xtestdata <- sidaData[[3]]
#' Ytest <- sidaData[[4]]
#' ##---- call cross validation
#' mycv=cvSIDA(Xdata,Y,withCov=FALSE,plotIt=FALSE, Xtestdata=Xtestdata,Ytest=Ytest)
#' ##----  Obtain relevance network
#' networkPlot(mycv,cutoff=0.7)


################################################################################
#Part of this function was borrowed from the mixOmics package and modified for mvlearnR
################################################################################
networkPlot=function(object, cutoff = NULL,
                           color.node = NULL,
                           lty.edge = c("solid","dashed"),
                           show.edge.labels = FALSE,
                           show.color.key = TRUE,
                           vertex.frame.color="red",
                           layout.fun = NULL,
                           save = NULL,
                           name.save = NULL){
  #call inner function to calculate similarity matrix

  mymat= networkplotinner(object)
  myComb=mymat$ViewCombinations
  nComb=dim(mymat$ViewCombinations)[2]

  if(is.null(cutoff)){
    cutoff2=0.5*as.vector.data.frame(matrix(1,nrow=nComb))
  }else if(length(cutoff)!=nComb){
    warning("'number of cutoff is less than the data combinations- setting others to 0.5'",
            call. = FALSE)
    cutoff2=c(cutoff,0.5*as.vector.data.frame(matrix(1,nrow=nComb-length(cutoff))))
  }else{
    cutoff2=cutoff
  }



  #code for generating network graph follow network plot in mixOmics.
  #nComb=dim(mymat$ViewCombinations)[2]
  res=list()
  for(jj in 1:nComb){

    mat=mymat$SimilirarityMatrix[[jj]]
    cutoff=cutoff2[jj]
    if(cutoff>max(abs(mat))){
      warning("'cutoff is greater than largest value in similarity matrix- setting to largest entry'",
              call. = FALSE)
      cutoff=max(abs(mat)-0.01)
    }

    if (!is.null(save)) {
      if (!save %in% c("jpeg", "tiff", "png", "pdf"))
        stop("'save' must be one of 'jpeg', 'png', 'tiff' or 'pdf'.",
             call. = FALSE)
    }
    if (!is.null(name.save)) {
      if (!is.character(name.save) || length(name.save) > 1)
        stop("'name.save' must be a character.", call. = FALSE)
    }else {
      if (!is.null(save)){
        #name.save = paste0("network_", gsub(".", "_", deparse(substitute(mat)),
        #                                    fixed = TRUE))

        name.save = paste0("network plot for views ", mycomb[1,jj], " and ", mycomb[2,jj])
      }

    }
    if (!is.null(save)) {
      while (dev.cur() > 1) dev.off()
      if (save == "jpeg")
        jpeg(filename = paste0(name.save, ".jpeg"), res = 600,
             width = 4000, height = 4000)
      if (save == "png")
        jpeg(filename = paste0(name.save, ".png"), res = 600,
             width = 4000, height = 4000)
      if (save == "tiff")
        tiff(filename = paste0(name.save, ".tiff"), res = 600,
             width = 4000, height = 4000)
      if (save == "pdf")
        pdf(file = paste0(name.save, ".pdf"))
    }

    p = nrow(mat)
    q = ncol(mat)
    row.names.plot = TRUE
    row.names = rownames(mat)
    col.names = colnames(mat)

    if (row.names.plot == TRUE) {
      row.names.plot = row.names
    }else {
      row.names.plot = rep("", p)
    }
    col.names.plot = TRUE

    if (col.names.plot == TRUE) {
      col.names.plot = col.names
    }else {
      col.names.plot = rep("", q)
    }

    if(is.null(vertex.frame.color)){
      vertex.frame.color="red"
    }
    if (is.null(color.node))
      color.node = c("white", "yellow")

    if (!is.list(color.node)) {
      if (!is.vector(color.node) || length(color.node) !=
          2)
        stop("'color.node' must be a vector of length 2.",
             call. = FALSE)
    }else {
      stop("'color.node' must be a vector of length 2.",
           call. = FALSE)
    }

    shape.node = c("circle", "rectangle")
    lwd.edge = 1
    cex.edge.label = 1

    if (length(lty.edge) == 1)
      lty.edge = c(lty.edge, lty.edge)
    choices = c("solid", "dashed", "dotted", "dotdash", "longdash",
                "twodash", "blank")
    lty.edge = choices[pmatch(lty.edge, choices, duplicates.ok = TRUE)]

    if (length(lwd.edge) == 1)
      lwd.edge = c(lwd.edge, lwd.edge)
    if (length(lwd.edge) != 2 || any(!is.finite(lwd.edge)) ||
        any(lwd.edge <= 0))
      stop("'lwd.edge' must be positive.")
    if (!is.logical(show.edge.labels))
      stop("'show.edge.labels' must be a logical constant (TRUE or FALSE).",
           call. = FALSE)
    if (!is.finite(cex.edge.label) || cex.edge.label < 0 || length(cex.edge.label) >
        1)
      stop("'cex.edge.label' must be a non-negative numerical value.",
           call. = FALSE)
    if (!is.logical(show.color.key))
      stop("'show.color.key' must be a logical constant (TRUE or FALSE).",
           call. = FALSE)

    # if (length(keysize) != 2 || any(!is.finite(keysize)))
    #   stop("'keysize' must be a numeric vector of length 2.",
    #        call. = FALSE)
    # if (length(keysize.label) != 1 || any(!is.finite(keysize)))
    #   stop("'keysize' must be a numeric vector of length 1.",
    #        call. = FALSE)

    keysize = c(1, 1)
    keysize.label = 1

    if (!is.null(layout.fun) && !is(layout.fun, "function"))
      stop("'layout.fun' must be a valid layout function.",
           call. = FALSE)

    w = as.vector(t(mat))

    node.X = row.names
    node.Y = col.names
    nodes = data.frame(name = c(node.X, node.Y), group = c(rep("x",
                                                               p), rep("y", q)))
    node.X = rep(node.X, each = q)
    node.Y = rep(node.Y, p)


    color.edge=color.GreenRed(100)
    cex.node.name=1

    relations = data.frame(from = node.X, to = node.Y, weight = w)
    id = bin.color(w, cutoff = cutoff, breaks = NULL, col = color.edge,
                   symkey = TRUE)
    col.id = id$bin
    color.edge = id$col[col.id]
    idx = (abs(w) >= cutoff)
    relations = relations[idx, ]
    color.edge = color.edge[idx]
    myplot = graph.data.frame(relations, directed = FALSE, vertices = nodes)
    V(myplot)$label.color = "black"
    V(myplot)$label.family = "sans"

    V(myplot)$label = c(row.names.plot, col.names.plot)
    V(myplot)$color = color.node[1]
    V(myplot)$color[V(myplot)$group == "y"] = color.node[2]
    V(myplot)$shape = shape.node[1]
    V(myplot)$shape[V(myplot)$group == "y"] = shape.node[2]

    if (show.edge.labels)
      E(myplot)$label = round(E(myplot)$weight, 2)

    E(myplot)$label.color = "black"
    E(myplot)$color = color.edge
    E(myplot)$lty = lty.edge[1]
    E(myplot)$lty[E(myplot)$weight < 0] = lty.edge[2]
    E(myplot)$width = lwd.edge[1]
    E(myplot)$width[E(myplot)$weight < 0] = lwd.edge[2]


    myplot = delete.vertices(myplot, which(degree(myplot) == 0)) #delete vertices with no edges
    lwid = c(keysize[1], 4)
    lhei = c(keysize[2], 4)
    lmat = matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE)
    nc = length(id$col)
    x = seq(0, 1, length = nc + 2)
    z.mat = seq(0, 1, length = nc + 1)
    z.mat = matrix(z.mat, ncol = 1)

    if ((id$lim[1] < -cutoff) & (id$lim[2] < cutoff)) {
      xv = c(0, x[nc + 1])
      lv = round(c(id$lim[1], -cutoff), 2)
      col = c(id$col, "white")
    }
    if ((id$lim[1] > -cutoff) & (id$lim[2] > cutoff)) {
      xv = c(x[2], 1)
      lv = round(c(cutoff, id$lim[2]), 2)
      col = c("white", id$col)
    }
    if ((id$lim[1] < -cutoff) & (id$lim[2] > cutoff)) {
      idn = max(which(id$breaks < 0))
      idp = min(which(id$breaks > 0))
      xv = c(0, x[idn + 1], x[idp], 1)
      lv = round(c(id$lim[1], -cutoff, cutoff, id$lim[2]),
                 2)
      col = c(id$col[1:idn], "white", id$col[(idn + 1):nc])
    }
    nn = vcount(myplot)
    V(myplot)$label.cex = min(2.5 * cex.node.name/log(nn), 1)
    E(myplot)$label.cex = min(2.25 * cex.edge.label/log(nn), 1)


    cex0 = 2 * V(myplot)$label.cex
    def.par = par(no.readonly = TRUE)
    dev.new()
    par(pty = "s", mar = c(0, 0, 0, 0), mfrow = c(1, 1))
    plot(1:100, 1:100, type = "n", axes = FALSE, xlab = "", ylab = "")
    cha = V(myplot)$label
    cha = paste("", cha, "")
    xh = strwidth(cha, cex = cex0) * 1.5
    yh = strheight(cha, cex = cex0) * 3
    V(myplot)$size = xh
    V(myplot)$size2 = yh
    dev.off()
    if (is.null(layout.fun)) {
      l = layout.fruchterman.reingold(myplot, weights = (1 - abs(E(myplot)$weight)))
      #l = layout.fruchterman.reingold(myplot, weights = E(myplot)$weight)
    }else {
      l = layout.fun(myplot)
    }

    if (isTRUE(show.color.key)) {

      layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
      par(mar = c(3, 2, 2, 2), cex = 0.75)
      image(z.mat, col = col, xaxt = "n", yaxt = "n")
      box()
      par(usr = c(0, 1, 0, 1))
      axis(1, at = xv, labels = lv, cex.axis = keysize.label)
      title("Color key", font.main = 1, cex.main = keysize.label)
      par(def.par)
      par(new = TRUE)
    }
    par(pty = "s", mar = c(0, 0, 0, 0), mfrow = c(1, 1))
    #igraph::plot.igraph(myplot, layout = l)
    igraph::plot.igraph(myplot, layout = l, vertex.size=20, vertex.size2=15, vertex.frame.width=2,
                        vertex.label.font=2,edge.size=4,edge.label.font=2,
                        edge.width=2,vertex.frame.color=vertex.frame.color)




    # igraph::plot.igraph(myplot[[i]],vertex.size=20,vertex.color="yellow", vertex.frame.width=3,
    #                     vertex.shape="square",vertex.label.font=2,vertex.label.font=2,vertex.frame.color="red",
    #                     edge.size=4, cex.edge.label=20, edge.label.font=2, edge.width=2, edge.color="black")


    par(def.par)

    res[[jj]] = list(NetworkGraph = myplot)

    mat[abs(mat) < cutoff] = 0
    res[[jj]]$SimilirarityMatrix = mat
    res[[jj]]$cutoff = cutoff
    res[[jj]]$pairs =t(myComb[,jj])

    # if (!is.null(save))
    #    dev.off()


  }
  return(invisible(res))

}

################################################################################
#This function was borrowed from the mixOmics package and modified for mvlearnR
################################################################################
bin.color=function (mat, cutoff, breaks, col, symkey)
{
  if(cutoff>max(abs(mat))){
    warning("'cutoff is greater than largest value in similarity matrix- setting to largest entry'",
            call. = FALSE)
    cutoff=max(abs(mat))
  }

  if (isTRUE(symkey)) {
    max.mat = max(abs(mat))
    min.mat = -max.mat
  }
  else {
    max.mat = max(mat)
    min.mat = min(mat)
  }
  if (missing(breaks) || is.null(breaks)) {
    if (is(col, "function"))
      breaks = 32
    else breaks = length(col)
  }
  if (length(breaks) == 1) {
    if (isTRUE(symkey)) {
      if ((breaks/2) - trunc(breaks/2) != 0)
        stop("'breaks' must be a even number if 'symkey = TRUE'",
             call. = FALSE)
      if (cutoff == 0) {
        breaks = c(seq(min.mat, max.mat, length = breaks +
                         1))
      }
      else {
        nb = breaks/2
        breaks = c(seq(min.mat, -cutoff, length = nb +
                         1), 0, seq(cutoff, max.mat, length = nb + 1))
        id = which(breaks == 0)
        breaks = breaks[-c(id - 1, id + 1)]
      }
    }
    else {
      breaks = breaks + 1
      if ((min.mat < -cutoff) & (max.mat < cutoff))
        breaks = seq(min.mat, -cutoff, length = breaks)
      if ((min.mat > -cutoff) & (max.mat > cutoff))
        breaks = seq(cutoff, max.mat, length = breaks)
      if ((min.mat < -cutoff) & (max.mat > cutoff)) {
        if (cutoff == 0) {
          breaks = c(seq(min.mat, max.mat, length = breaks))
        }
        else {
          long = max.mat - min.mat - 2 * cutoff
          bin = long/breaks
          breaks = seq(cutoff, -min.mat, by = bin)
          o = order(breaks, decreasing = TRUE)
          breaks = c(-breaks[o], 0, seq(cutoff, max.mat,
                                        by = bin))
          id = which(breaks == 0)
          breaks = breaks[-c(id - 1, id + 1)]
        }
      }
    }
  }
  ncol = length(breaks) - 1
  if (is(col, "function"))
    col = col(ncol)
  if (length(breaks) != length(col) + 1)
    stop("must have one more break than colour", call. = FALSE)
  min.breaks = min(breaks)
  max.breaks = max(breaks)
  mat[mat < min.breaks] = min.breaks
  mat[mat > max.breaks] = max.breaks
  bin = .bincode(as.double(mat), as.double(breaks), TRUE, TRUE)
  return(invisible(list(bin = bin, col = col, breaks = breaks,
                        lim = c(min.mat, max.mat))))
}

#############################################################################################################
# Authors:
#   Sandra E. Safo, Unviersity of Minnesota
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



networkplotinner=function(object){

  if(class(object)=="SIDA" | class(object)=="SIDANet"){
    hatalpha=object$hatalpha
    #L=dim(hatalpha[[1]])[2]
    L=length(hatalpha)
    for(j in 1:L){
      hatalpha[[j]]=qr.Q(qr(hatalpha[[j]]))
    }
  }else if( class(object)=="SELPCCA"){
    if(object$method=="selpscca.pred"){
      hatalpha=list(object$selp.fit$hatalpha,object$selp.fit$hatbeta
                    )
      L=dim(hatalpha[[1]])[2]
      for(j in 1:L){
        hatalpha[[j]]=qr.Q(qr(hatalpha[[j]]))
      }

    }else{ hatalpha=list(object$hatalpha,object$hatbeta)
    L=dim(hatalpha[[1]])[2]
    for(j in 1:L){
      hatalpha[[j]]=qr.Q(qr(hatalpha[[j]]))
    }
    }

  }



  Mjk=list()
  D <- length(hatalpha)
  #obtain variables that are selected
  hatalpha.temp=list()
  InputData.temp=list()

  #L=dim(hatalpha[[1]])[2]



  mycomb=utils::combn(D,2) #pairwise combination
  ncomp=dim(hatalpha[[1]])[2]

  for(jj in 1:dim(mycomb)[2]){

    X1=as.data.frame(object$InputData[[mycomb[1,jj]]])
    X2=as.data.frame(object$InputData[[mycomb[2,jj]]])


    if(ncomp == 1){
      #hatalpha.temp[[i]]=abs(hatalpha[[i]]!=0)
      hatalpha1=hatalpha[[mycomb[1,jj]]]
      hatalpha2=hatalpha[[mycomb[2,jj]]]
    }else if (ncomp >1){
      #hatalpha.temp[[i]]=rowSums(abs(hatalpha[[i]]))
      hatalpha1=rowSums(abs(hatalpha[[mycomb[1,jj]]]))
      hatalpha2=rowSums(abs(hatalpha[[mycomb[2,jj]]]))
    }

    X1var.Ind=which(as.matrix(hatalpha1)!=0, arr.ind = TRUE)
    X1var.Ind=X1var.Ind[,1]
    colname1=sub("\\;.*", "", colnames(X1[,X1var.Ind]))

    UL=as.matrix(X1[,X1var.Ind])%*%hatalpha[[mycomb[1,jj]]][X1var.Ind,]
    #UL=as.matrix(X1)%*%hatalpha[[mycomb[1,jj]]]

    #variable selected for second view
    X2var.Ind=which(as.matrix(hatalpha2)!=0, arr.ind = TRUE)
    X2var.Ind=X2var.Ind[,1]
    colname2=sub("\\;.*", "", colnames(X2[,X2var.Ind]))
    VL=as.matrix(X2[,X2var.Ind])%*%hatalpha[[mycomb[2,jj]]][X2var.Ind,]
    #VL=as.matrix(X2)%*%hatalpha[[mycomb[2,jj]]]
    #bisect
    ZL=UL + VL


    X11corr=   stats::cor(as.matrix(X1[,X1var.Ind]), ZL,use = "pairwise" )
    X21corr=   stats::cor(as.matrix(X2[,X2var.Ind]), ZL, use = "pairwise" )


    rownames(X11corr)=colname1
    rownames(X21corr)=colname2
    if(length(intersect(colname1,colname2))!=0){
      rownames(X11corr) <- paste(rownames(X11corr),mycomb[1,jj],sep="_")
      rownames(X21corr) <- paste(rownames(X21corr),mycomb[2,jj],sep="_")
    }
    simmat=as.matrix(X11corr)%*%as.matrix(t(X21corr))
    if(max(abs(simmat))>1){
      Mjk[[jj]]=simmat/(max(abs(simmat))+0.001)
    }else{
      Mjk[[jj]]=as.matrix(X11corr)%*%as.matrix(t(X21corr))
    }
  }



  return(invisible(list(SimilirarityMatrix=Mjk, ViewCombinations=mycomb)))
}



#############################################################################################################
# # This code was borrowed from the mixOmics Package
#Authors:
#   Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
#
# created: 2015
# last modified:
#
# Copyright (C) 2015
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


#-- green-black-red gradient colors --#
#-------------------------------------#
color.GreenRed =
  function (n, alpha = 1)
  {
    #-- checking general input parameters --------------------------------------#
    #---------------------------------------------------------------------------#

    #-- n
    if (length(n) > 1 || !is.finite(n))
      stop("'n' must be an integer positive value.", call. = FALSE)

    if (n < 1)
      stop("'n' must be an integer positive value.", call. = FALSE)

    #-- alpha
    if (length(alpha) > 1 || !is.finite(alpha))
      stop("'alpha' must be an numeric value in the range [0, 1].", call. = FALSE)

    if (alpha < 0 || alpha > 1)
      stop("'alpha' must be an numeric value in the range [0, 1].", call. = FALSE)

    alpha = round(255 * alpha)

    #-- end checking --#
    #------------------#

    ramp = colorRampPalette(c("green", "darkgreen", "black", "darkred", "red"))
    ramp = ramp(101)
    green = ramp[1:43]
    red = ramp[59:101]
    ramp = colorRamp(c(green, "black", red), space = "Lab")

    rgb(ramp(seq(0, 1, length = n)), alpha = alpha, maxColorValue = 255)
  }


#motivated by color.GreenRed Palette
#modified by Sandra Safo
color.BlueOrange =
  function (n, alpha = 1)
  {
    #-- checking general input parameters --------------------------------------#
    #---------------------------------------------------------------------------#

    #-- n
    if (length(n) > 1 || !is.finite(n))
      stop("'n' must be an integer positive value.", call. = FALSE)

    if (n < 1)
      stop("'n' must be an integer positive value.", call. = FALSE)

    #-- alpha
    if (length(alpha) > 1 || !is.finite(alpha))
      stop("'alpha' must be an numeric value in the range [0, 1].", call. = FALSE)

    if (alpha < 0 || alpha > 1)
      stop("'alpha' must be an numeric value in the range [0, 1].", call. = FALSE)

    alpha = round(255 * alpha)

    #-- end checking --#
    #------------------#

    ramp = colorRampPalette(c("#56B4E9", "blue", "#999999", "darkorange", "orange"))
    ramp = ramp(101)
    blue = ramp[1:43]
    orange = ramp[59:101]
    ramp = colorRamp(c(blue, "#999999", orange), space = "Lab")

    rgb(ramp(seq(0, 1, length = n)), alpha = alpha, maxColorValue = 255)
  }

