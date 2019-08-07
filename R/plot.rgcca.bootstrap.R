plot.bootstrap.rgcca = function(object, nb_bloc = length(object$A)){
  if (class(object) != "rgcca.bootstrap"){
    return("object is not rgcca.bootstrap !")
  }
  for (i in 1:nb_bloc){
    color = rep("green3", nrow(object$CI_sel[[i]]))
    Sel     = barplot(object$CI_sel[[i]][, 1], col = color, 
                      ylim = c(min(0, min(object$CI_sel[[i]][, 2])), max(0, object$CI_sel[[i]][, 3])),
                      las = 1,  omi = c(50, 4, 4, 4))
    segments(Sel,object$CI_sel[[i]][, 2], Sel, object$CI_sel[[i]][, 3])
    segments(Sel-0.1, object$CI_sel[[i]][, 2], Sel+0.1, object$CI_sel[[i]][, 2])
    segments(Sel-0.1, object$CI_sel[[i]][, 3], Sel+0.1, object$CI_sel[[i]][, 3])
    title(main = paste("Selected variables in bloc", i),cex.main = 3)
  }
}


