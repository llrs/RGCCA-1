#' Create a biplot from a rgcca object. 
#' The first two principal components of the superbloc are used as a common space to plot the individuals and the variables.
#' Variables coordinates are computed by a multiple linear regression of the variables on the principal components
#' Require packages ggplot and purrr
#' 
#' @param object A RGCCA object
#' @param superbloc matrix. The superbloc used for the rgcca
#' @param superbloc_pos Numeric. Position of the superbloc in the rgcca object (default is last bloc)
#' @param var_names Vector of characters. Variables to be plotted in the biplot (default is all the variables of the superbloc (column))
#' @param ind_names Vector of characters. Names of the individuals (default is the rownames of the superbloc)
#' @param ind_col Vector. Colors of the individuals in the biplot if specified
#' @param var_col Vector. Colors of the individuals in the biplot if specified
#' @param scale Logical. True if the superbloc needs to be scaled for the regression
#' @param expend Numeric. To scale the ratio individual/variable
#' @param show_r_squared Logical. Should the r.squared coefficient be shown
#' @param arrows Logical. Should arrows be plotted to the variable coordinates
#' @param main String of character. The main title of the biplot
#' 
#' @return a ggplot biplot 
#' 
#' @examples
#'
#' #############
#' # Example 1 #
#' #############
#' rm(list=ls())
#' library(ggplot2)
#' library(RGCCA)
#' library(purrr)
#' ############ Chargement et mise en place des données ###############################
#' data(Russett)
#' X_agric =as.matrix(Russett[,c("gini","farm","rent")])
#' X_ind = as.matrix(Russett[,c("gnpr","labo")])
#' X_polit = as.matrix(Russett[ , c("demostab", "demoinst","dictator")])
#' A = list(X_agric, X_ind, X_polit)
#' 
#' 
#' ########### rgcca avec superbloc ###################################################
#' superbloc = cbind(A[[1]],A[[2]])
#' B = list(X_agric, X_ind, superbloc, X_polit)
#' C = matrix(c(0, 0, 1, 0,
#'              0, 0, 1, 0,
#'              1, 1, 0, 1,
#'              0, 0, 1, 0), 4, 4)
#' B = lapply(B,scale2)
#' fit.rgcca = rgcca(B, C, tau = c(1,1,1,1), scheme = "factorial", scale = FALSE, ncomp = c(1,1,2,1))
#' ind_col = factor(apply(Russett[, 9:11], 1, which.max),labels =c("demostab", "demoinst", "dictator"))
#' rgcca.biplot(fit.rgcca, superbloc, ind_col = ind_col)
#' 
#' 
#' 
#'@export rgcca.biplot

library(ggplot2)
library(ggrepel)

rgcca.biplot = function(object, 
                        superbloc, 
                        superbloc_pos = length(object$Y), 
                        var = colnames(superbloc),
                        var_labels = var, 
                        var_col = rep("red",length(var_labels)),
                        ind = rownames(superbloc),
                        ind_labels = ind,
                        ind_col = rep(1,length(ind_names)),
                        ind_size = 1,
                        scale = TRUE,
                        expand = 1,
                        show_r_square = FALSE,
                        main = NULL,
                        arrows = TRUE,
                        shorten_label = T,
                        target = NULL){
  if (dim(object$Y[[superbloc_pos]])[2] < 2){
    return("Two components on the superbloc are required for the biplot")
  } else {
    Y = object$Y[[superbloc_pos]]
  }
  if (scale == TRUE) {superbloc = scale2(superbloc)}
  r = list()
  r.square = list()
  for (i in 1:length(var)){
    r[[i]] = lm(superbloc[,var[i]] ~ Y[,1] + Y[,2])
    # print(paste(var_labels[i],"Adjusted R-squared :",summary(r[[i]])$adj.r.squared))
    r.square[[i]] = summary(r[[i]])$adj.r.squared
    r[[i]] = r[[i]]$coefficient[2:3]
  }
  r = reduce(r,rbind)
  r.square = reduce(r.square,rbind)
  rownames(r) = rownames(r.square) = var_labels
  if (show_r_square) {print(as.matrix(r.square[order(r.square,decreasing = TRUE),1]))}
  
  unsigned.range <- function(x)
    c(-abs(min(x, na.rm=TRUE)), abs(max(x, na.rm=TRUE)))
  rangx1 <- unsigned.range(Y[, 1L])
  rangx2 <- unsigned.range(Y[, 2L])
  rangy1 <- unsigned.range(r[, 1L])
  rangy2 <- unsigned.range(r[, 2L])
  
  if (shorten_label == T){
    library(stringr)
    shorten_label = function(x){
      if (str_length(x) >10){
        return(str_c(str_sub(x,1,10),"..."))
      } else {
        return(x)
      }
    }
    var_labels = vapply(var_labels,shorten_label,"a")
  }
  
  
  
  ratio <- max(rangy1/rangx1, rangy2/rangx2)/expand
  r=r*ratio
  df1 = data.frame(ind,Y,ind_col)
  df2 = data.frame(y1 = r[,1], y2 = r[,2])
  g = ggplot(data = df1, aes(comp1,comp2))+
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = 0) + 
      # geom_text_repel(aes(colour = ind_col, label = ind_labels)) +
      geom_point(aes(colour = ind_col), size = ind_size) +
      # geom_label_repel(label = round(100*ind,1)) +
      scale_colour_gradient(low = "yellow", high = "brown") +
      # geom_label_repel(label = round(100*atrophy,1),label.size = 0.25,size = 2) +
      # geom_text(data = df2, aes(y1,y2, label = var_labels, colour = var_col)) +
      geom_text_repel(data = df2, aes(y1,y2, label = var_labels), colour = var_col) +
      ggtitle(main)    
  
  if (arrows == TRUE){
    g = g + geom_segment(data = df2, aes(x = 0, y = 0, xend = y1, yend = y2), 
                         arrow=arrow(length=unit(0.2,"cm")), colour = var_col) 
  } 
  if (!is.null(target)){
    g = g + geom_segment(data = target, aes(x = 0, y = 0, xend = x, yend = y),
                             arrow=arrow(length=unit(0.2,"cm")), colour = 1)
    g = g + geom_text_repel(data = target, aes(x,y, label = rownames(target)), colour = 1)
  }
  g
  return(g)
}
