
#' This function does a fuzzy clusterwise RGCCA. The user choose a number of cluster to consider. The function will then build one model by cluster, using for the cluster c, the blocks \eqn{\mathbf{U_c^mX_1}, \mathbf{U_c^mX_2}, ..., \mathbf{U_c^mX_J}}, where \eqn{\mathbf{U_c}} is a diagonal matrix, and \eqn{\sum_cu_{ci}=1}, m denote the user prescribed fuzzy weight scalar, often called the “fuzzifier” (Bezdek, 1974a).
#' \eqn{\mathbf{U_c}} are randomly initialized. The functions has two steps that are repeated until convergence: 
#' First step : With \eqn{\mathbf{U_c}} fixed, builds the RGCCA models for each cluster.
#' Second step : With RGCCA weights fixed, optimize \eqn{\mathbf{U_c}}


#' @param A  A list that contains the \eqn{J} blocks of variables \eqn{\mathbf{X_1}, \mathbf{X_2}, ..., \mathbf{X_J}}.
#' @param n_cl Number of cluster of the model
#' @param init_cw Choice of initialisation of \eqn{\mathbf{U_c}}. "Random" or an integer. If integer, the function will choose from a number of random \eqn{\mathbf{U_c}} the ones that maximises RGCCA at the first iteration.
#' @param m An integer, it denotes the user prescribed fuzzy weight scalar, often called the “fuzzifier” (Bezdek, 1974a)
#' @param max_iter An integer, maximum number of iteration in the algorithm if convergence is not reached.
#' @param C  A design matrix that describes the relationships between blocks (default: complete design).
#' @param tau tau is either a \eqn{1 \times J} vector or a \eqn{\mathrm{max}(ncomp) \times J} matrix, and contains the values 
#' of the shrinkage parameters (default: tau = 1, for each block and each dimension).
#' If tau = "optimal" the shrinkage paramaters are estimated for each block and each dimension using the Schafer and Strimmer (2005)
#' analytical formula . If tau is a \eqn{1\times J} numeric vector, tau[j] is identical across the dimensions of block \eqn{\mathbf{X}_j}. 
#' If tau is a matrix, tau[k, j] is associated with \eqn{\mathbf{X}_{jk}} (\eqn{k}th residual matrix for block \eqn{j})
#' @param ncomp  A \eqn{1 \times J} vector that contains the numbers of components for each block (default: rep(1, length(A)), which gives one component per block.)
#' @param scheme The value is "horst", "factorial", "centroid" or the g function (default: "centroid").
#' @param scale  If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
#' @param verbose  If verbose = TRUE, the progress will be report while computing (default: TRUE).
#' @param init The mode of initialization to use in RGCCA algorithm. The alternatives are either by Singular Value Decompostion ("svd") or random ("random") (Default: "svd").
#' @param bias A logical value for biaised or unbiaised estimator of the var/cov (default: bias = TRUE).
#' @param tol The stopping value for convergence.
#' @param sameBlockWeight A logical value indicating if the different blocks should have the same weight in the analysis (default, sameBlockWeight=TRUE)
#' @param returnA A logical value indicating if the A list should be return as a result (default, sameBlockWeight=FALSE)

#' @return \item{fit.rgcca.cw}{}
#' @return \item{max_step}{}
#' @return \item{U}{}
#' @references Tenenhaus A. and Tenenhaus M., (2011), Regularized Generalized Canonical Correlation Analysis, Psychometrika, Vol. 76, Nr 2, pp 257-284.
#' @references Tenenhaus A. et al., (2013), Kernel Generalized Canonical Correlation Analysis, submitted.
#' @references Schafer J. and Strimmer K., (2005), A shrinkage approach to large-scale covariance matrix estimation and implications for functional genomics. Statist. Appl. Genet. Mol. Biol. 4:32.
#' @title Regularized Generalized Canonical Correlation Analysis (RGCCA) 
#' @examples
#' 
#' 
#' #############
#' # Example 1 #
#' #############
#' library(RGCCA)
#' C = matrix(c(0,1,1,
#'              1,0,1,
#'              1,1,0),3,3)
#' X_agric =as.matrix(Russett[,c("gini","farm","rent")])
#' X_ind = as.matrix(Russett[,c("gnpr","labo")])
#' X_polit = as.matrix(Russett[ , c("ecks", "death")])
#' 
#' 
#' lab = as.vector(apply(Russett[, 9:11], 1, which.max))
#' ind_demoins = which(lab==2)
#' A = list(X_agric, X_ind, X_polit)
#' A = lapply(A, function(x) x[-ind_demoins,])
#' A = lapply(A, function(x) scale2(x,scale=TRUE, bias = TRUE))
#' A = lapply(A, function(x) x/sqrt(NCOL(x)))
#' n = dim(A[[1]])[1]    # Number of individuals
#' n_cl = 2
#' 
#' res = cw.rgcca.hwang(A = A,
#'                      n_cl = n_cl, 
#'                      init_cw = 1000,
#'                      # U = list(diag(1,nrow(A[[1]]))),
#'                      max_iter = 300, 
#'                      m = 2, 
#'                      C = 1 - diag(length(A)), 
#'                      tau = rep(0, length(A)), 
#'                      ncomp = rep(1,length(A)), 
#'                      scheme = "horst", 
#'                      scale = FALSE, 
#'                      init = "svd",
#'                      bias = TRUE, 
#'                      tol = 1e-08, 
#'                      verbose = FALSE)
#' 
#' cl = apply(sapply(res$U,diag),1,which.max)
#' 
#' result.rgcca = rgcca(A = A,
#'                      C = C,
#'                      tau = rep(0,3),
#'                      ncomp = rep(1,3),
#'                      scheme = "horst",
#'                      scale = F,
#'                      verbose = F)
#'
#' lab = lab[-ind_demoins]
#' lab[which(lab == 3)] = 2
#' y = apply(as.matrix(rownames(Russett)),2, function(x) substr(x,1,2))
#' y = y[-ind_demoins]
#' plot(result.rgcca$Y[[1]], result.rgcca$Y[[2]], col = "white", 
#'      xlab = "Y1 (Agric. inequality)", ylab = "Y2 (Industrial Development)")
#' text(result.rgcca$Y[[1]], result.rgcca$Y[[2]], rownames(Russett)[-ind_demoins], col = lab, cex = .7)
#' text(result.rgcca$Y[[1]], result.rgcca$Y[[2]], adj = c(0.1,-1), y, col = cl, cex = .7)
#'
#' table(lab,cl)


cw.rgcca.hwang = function(A,
                          n_cl, 
                          init_cw = "random",
                          max_iter = 200, 
                          m = 2, 
                          C = 1 - diag(length(A)), 
                          tau = rep(1, length(A)), 
                          ncomp = rep(1,length(A)), 
                          scheme = "centroid", 
                          scale = TRUE, 
                          init = "svd",
                          bias = TRUE, 
                          tol = 1e-08, 
                          verbose = TRUE){
  U = init.cw(n_cl = n_cl,
              n_ind = nrow(A[[1]]), 
              init_cw = init_cw,
              C = C,
              tau = tau, 
              ncomp = ncomp, 
              scheme = scheme, 
              scale = scale, 
              init = init,
              bias = bias, 
              tol = tol, 
              verbose = verbose)
  x=list()
  crit = list()
  UU = U
  n_ind = nrow(A[[1]])
  for (u in 1:max_iter){
    fit.rgcca = cw.rgcca(A =A, 
                         U = U, 
                         C = C,
                         tau = tau, 
                         ncomp = ncomp,
                         scheme = scheme,
                         scale = scale,
                         init = init,
                         bias = bias, 
                         tol = tol, 
                         verbose = verbose)
    crit[[u]] = fit.rgcca$crit
    e = lapply(1:n_cl,function(x) rep(0,n_ind))
    for (i in 1:n_ind){
      for (c in 1:n_cl){
        for (k in 1:length(A)){
          for (j in 1:length(A)){
            x[[j+length(A)*(k-1)]] = C[j,k]*(A[[k]]%*%fit.rgcca$fit.rgcca[[c]]$a[[k]]-A[[j]]%*%fit.rgcca$fit.rgcca[[c]]$a[[j]])
            # if (j!=k){x[[j+length(A)*(k-1)]] = A[[k]]%*%fit.rgcca[[c]]$a[[k]]-A[[j]]%*%fit.rgcca[[c]]$a[[j]]}
            # else{x[[j+length(A)*(k-1)]] = rep(0,n)}
          }
        }
        e[[c]][i] = Reduce("+",lapply(x, function(z) z[i]^2))
      }
      for (c in 1:n_cl){
        UU[[c]][i,i] = 1/(Reduce("+",lapply(e, function(x) (e[[c]][i]/x[i])^(1/(m-1)))))
      }
    }
    iter = abs(UU[[1]]-U[[1]])
    if (max(iter) < tol){
      print(paste("Convergence at step",u))
      umax = u
      return(list(fit.rgcca.cw =fit.rgcca, max_step = u, U = UU))
    } else {
      u=u+1
      U = UU
    }
    # Debug
    # if (any(is.na(iter))){
    #   print(paste("Stop at step",u))
    #   umax = u
    # }
  }
  return(list(fit.rgcca.cw =fit.rgcca, max_step = u, U = UU,crit = ))
}
