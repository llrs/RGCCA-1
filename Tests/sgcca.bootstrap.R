#' Return variables that are most often selected by a sgcca with permutated individuals. The function uses parallelisation on windows (parallel package). This first version is built for two explicatives blocks

#' @param A A list that contains the J blocks of variables \mathbf{X_1}, \mathbf{X_2}, ..., \mathbf{X_J}
#' @param B A number: how many permutation the function must perform
#' @param c1 A vector sparcity parameters used by the sgcca
#' @param C A design matrix that describes the relationships between blocks (default: complete design)
#' @param tau tau is either a 1 \times J vector or a \mathrm{max}(ncomp) \times J matrix, and contains the values of the shrinkage parameters (default: tau = 1, for each block and each dimension). If tau = "optimal" the shrinkage paramaters are estimated for each block and each dimension using the Schafer and Strimmer (2005) analytical formula . If tau is a 1\times J numeric vector, tau[j] is identical across the dimensions of block \mathbf{X}_j. If tau is a matrix, tau[k, j] is associated with \mathbf{X}_{jk} (kth residual matrix for block j)
#' @param scheme The value is "horst", "factorial" or "centroid" (default: "centroid")
#' @param top Number of variables to return
#' @param plot A logical, should a plot of top variables be displayed
#' @return A list containing for each block a matrix with the top selected variables and their occurences, ordered by occurence
#' @value \item {top_var} list with variables selected according to the top parameter for each bloc
#' @value \item {count} list with the percentage of times each variable was selected for each bloc


#  out= list(top_var = top_var, count = count_sort)

bootstrap_sgcca = function(A,nb_boot,c1,
                           C = 1 - diag(length(A)),
                           tau = rep(1,length(A)), 
                           scheme = "centroid",
                           top = 15,
                           tol = .Machine$double.eps,
                           plot = FALSE){
  J = length(A) - 1
  trials = seq(1,nb_boot)
  
  n_core = parallel::detectCores()
  cl = makeCluster(n_core-1)
  e = environment()
  parallel::clusterExport(cl,c("A","c1","C","trials","scheme","J","tol"),envir = e)  #Variables utilisées dans le cluster
  parallel::clusterEvalQ(cl,library(RGCCA))
  
  boot_sgcca = function(trial){
    ind = sample(NROW(A[[1]]),replace = TRUE)
    Bscr = lapply(A, function(x) x[ind,])
    res = sgcca(Bscr, C,c1 = c1, ncomp = rep(1,length(A)),
                scheme = scheme, scale = TRUE,tol = tol)
    STAB = res$a
    return(STAB)
  }
  
  results = parLapply(cl,trials,boot_sgcca) 
  res = list()    #Pour concatener les résultats dans une liste de trois blocs. au final 1 élement par block, 1 variable par colonne et la valeur prise pour chaque bootstrap en ligne
  for (i in 1:J){
    res[[i]] = results[[1]][[i]]
    for (j in 2:nb_boot){
      res[[i]] = cbind(res[[i]],results[[j]][[i]])
    }
  }
  
  parallel::stopCluster(cl)
  
  count = lapply(res, function(x) apply(x,1,function(y) sum(y!=0)))
  count_sort = lapply(count, function(x) sort(x,decreasing = TRUE))
  
  if (plot){
    g=list()
    Stab = list()
    for (j in 1:J){
      Stab[[j]] = data.frame(Variables = names(count_sort[[j]])[1:top],Count = count_sort[[j]][1:top]/nb_boot)
      g[[j]] = ggplot(Stab[[j]], aes(x = reorder(Variables, Count), y = Count))+
        coord_flip() +
        geom_bar(stat = "identity")+
        ggtitle(paste0("Stability selection for block ", j, ", nb_boot = ",nb_boot))+
        xlab("Metabolites")+
        ylab("Occurence")
    }
    lapply(g,function(x) print(x))
  }
  top_var = lapply(count_sort, function(x) names(x)[1:top])
  out= list(top_var = top_var, count = count_sort)
  return(out)
}
