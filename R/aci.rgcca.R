#' Run a fixed number of RGCCA with bootstrap
#' 
#' 
#' @param object A RGCCA object
#' @param A A list that contains the J blocks of variables \mathbf{X_1}, \mathbf{X_2}, ..., \mathbf{X_J}.
#' @param B A numeric value. Number of bootstrap to be done
#' @param alpha A numeric value between 0 and 1. Extreme values to be removed to select variables
#' @param plot A character string. "all" for plotting every variables, "selected" for plotting only impacting variables
#' @param ndim A numeric value. The component to consider for the weight of the variables
#' @param verbose (mac-linux only)
#' @param nb_cores (mac-linux only) An integer. Number of cores to use for the parallelisation. 
#' (In windows it takes all but one)
#' 
#' @return \item {a_boot} weigth vectors for each iteration of the bootstrap
#' @return \item {CI} CI coefficients: for each variable the mean, the min and max without extreme values
#' @return \item {CI_sel} CI coefficients only for selected variables (with consistent direction across simulations)
#' @return \item {inner_relation} 
#' 
#' return(list(a_boot = Astar, CI = mat, CI_sel = mat_sel, inner_relation = inner_relation))
#' 
aci.rgcca = function(object, 
                     A, 
                     B, 
                     alpha = 0.05, 
                     ndim = 1, 
                     verbose = FALSE, 
                     plot = "selected", 
                     nb_cores = 1) {
  n = NROW(A[[1]])
  p = sapply(A, ncol)
  J = length(A)
  
  ##############################################
  # Initialization of the outer weight vectors #
  ##############################################
  if(any(object$ncomp!=1)){
    Yinit <- sapply(object$Y, function(x) x[, ndim])
  }else{
    Yinit <- sapply(object$Y, cbind)
  }  
  W = lapply(object$a, function(x) x[, ndim])
  
  ########################################
  # Construction of the Boostrap samples #
  if (Sys.info()["sysname"] == "Windows"){
    n_cores = parallel::detectCores()
    cl = parallel::makeCluster(n_cores-1) #7 coeurs
    e = environment()

    parallel::clusterExport(cl,c("object", "A","n","J","ndim","W","B", "bootstrap"),envir = e)  #Variables utilis?es dans le cluster
    parallel::clusterEvalQ(cl, library(RGCCA))
    # tryCatch({parallel::clusterEvalQ(cl, load_all("RGCCA/R/."))},
    #          error = function(e) {warning("error : please check if you are in the directory containing RGCCA functions")})
    # parallel::clusterEvalQ(cl, load_all("RGCCA/R/."))
    boot_b  = parallel::parLapply(cl,1:B, function(z) bootstrap(n = n, J = J, A = A, object = object, W = W, ndim = ndim))
    # tryCatch({boot_b  = parallel::parLapply(cl,
    #                                          1:B, 
    #                                          function(z) bootstrap(n = n, 
    #                                                                J = J, 
    #                                                                A = A, 
    #                                                                object = object, 
    #                                                                W = W, 
    #                                                                ndim = ndim))},
    #          error = function(e) {warning("an error occured with sgcca.crit"); return(NULL) })
    parallel::stopCluster(cl)
    # boot_b  = lapply(1:B, function(z) bootstrap(n = n, J = J, A = A, object = object, W = W, ndim = ndim))
  } else {
    if (verbose){
      boot_b  = pbmcapply::pbmclapply(1:B, function(z) bootstrap(n = n, J = J, A = A, object = object, W = W, ndim = ndim), mc.cores = nb_cores)
    }else{
      boot_b  = parallel::mclapply(1:B, function(z) bootstrap(n = n, J = J, A = A, object = object, W = W, ndim = ndim), mc.cores = nb_cores)
    }
  }
  ########################################
  # Construction of the Coefficeints'CI  #
  ########################################
  Astar = NULL
  J = ifelse(dim(A[[J]])[2] == 1,J-1,J)  # Avoid an error if the bloc to explain has only one dimension (a vector)
  for (i in 1:J){
    Astar[[i]] = sapply(1:B, FUN = function(x) boot_b[[x]][[2]][[i]])
  }
  
  M1 = lapply(Astar, function(w) apply(w, 1,  function(x) c(mean(x), sd(x))))
  
  mat  = list()
  tail = qnorm(1-alpha/(2))
  for (j in 1:J){
    mat[[j]]           <- cbind(W[[j]], M1[[j]][1, ]-tail*M1[[j]][2, ], 
                                M1[[j]][1, ]+tail*M1[[j]][2, ],
                                pnorm(0, mean = abs(M1[[j]][1,]), sd = M1[[j]][2,]),
                                pnorm(0, mean = abs(M1[[j]][1,]), sd = M1[[j]][2,])*length(M1[[j]][1,]))
    rownames(mat[[j]]) <- colnames(A[[j]])
    colnames(mat[[j]]) <- c("Initial weights","Lower Bound","Upper Bound","p(X)>0","p(X)>0 corrected bonf")
  }
  
  ########################################
  # Construction of the Correlations'CI  #
  ########################################
  J = length(A)
  
  MAT_COR = t(sapply(1:B, FUN = function(x) boot_b[[x]][[1]]))
  M2      = apply(MAT_COR, 2,  function(x) c(mean(x), sd(x)))
  
  inner_relation           = matrix(0, 3, (J*(J-1))/2)
  connection               = cbind(rep(paste0("X", 1:J), J), rep(paste0("X", 1:J), each = J))
  colnames(inner_relation) = matrix(paste0(connection[, 1], "-", connection[, 2]), J, J)[upper.tri(diag(J))]
  rownames(inner_relation) = c("Initial Coorelation","Lower Bound","Upper Bound")
  inner_relation[1, ]      = (cor(Yinit)*object$C)[upper.tri(diag(J))]
  inner_relation[2, ]      = (matrix(M2[1, ]-tail*M2[2, ], J, J)*object$C)[upper.tri(diag(J))]
  inner_relation[3, ]      = (matrix(M2[1, ]+tail*M2[2, ], J, J)*object$C)[upper.tri(diag(J))]
  inner_relation           = t(inner_relation)
  
  
  ########################################
  # Construction of a vizualisation      #
  ########################################
  J = ifelse(dim(A[[J]])[2] == 1,J-1,J)
  mat = lapply(mat, function(x) x[order(x[, 1], decreasing = TRUE), ])
  mat_sel = lapply(mat, function(x) x[which(x[,2]/x[,3]>0),])
  if (plot == "all"){
    par(cex = .8)
    for (j in 1:J){
      color = rep("red", nrow(mat[[j]])) ; color[which(mat[[j]][, 2]/mat[[j]][, 3]>0)] = "green3"
      r     = barplot(mat[[j]][, 1], col = color, ylim = c(min(0, min(mat[[j]][, 2])), max(0, mat[[j]][, 3])), las = 2,  omi = c(50, 4, 4, 4))
      segments(r, mat[[j]][, 2], r, mat[[j]][, 3])
      segments(r-0.1, mat[[j]][, 2], r+0.1, mat[[j]][, 2])
      segments(r-0.1, mat[[j]][, 3], r+0.1, mat[[j]][, 3])
    }
  }else if (plot == "selected"){
    par(cex = .8)
    for (j in 1:J){
      color = rep("green3", nrow(mat_sel[[j]]))
      Sel     = barplot(mat_sel[[j]][, 1], col = color, 
                        ylim = c(min(0, min(mat_sel[[j]][, 2])), max(0, mat_sel[[j]][, 3])),
                        las = 1,  omi = c(50, 4, 4, 4))
      segments(Sel, mat_sel[[j]][, 2], Sel, mat_sel[[j]][, 3])
      segments(Sel-0.1, mat_sel[[j]][, 2], Sel+0.1, mat_sel[[j]][, 2])
      segments(Sel-0.1, mat_sel[[j]][, 3], Sel+0.1, mat_sel[[j]][, 3])
    }
  }
  out = list(a_boot = Astar, CI = mat, CI_sel = mat_sel, inner_relation = inner_relation, A = A)
  class(out) = "rgcca.bootstrap"
  return(out)
}
