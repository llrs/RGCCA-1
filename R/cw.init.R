
# Internal function used by cw.rgcca.hwang. It initiates the matrix U. 
# If init_cw = "random" is selected, only one set of matrix are generated.
# If init_cw = n with n an integer, n random set of matrix are generated, the one that maximise the criteria of clusterwise RGCCA.

init.cw = function(n_cl,
                   n_ind, 
                   init_cw = "random",
                   C = 1 - diag(length(A)),
                   tau = rep(1, length(A)), 
                   ncomp = rep(1,length(A)), 
                   scheme = "centroid", 
                   scale = TRUE, 
                   init = "svd",
                   bias = TRUE, 
                   tol = 1e-08, 
                   verbose = TRUE){
  if (init_cw == "random"){
    U = lapply(1:n_cl, function(x) diag(runif(n_ind)))
    SumDiag = Reduce("+",U)
    U = lapply(U, function(x) x/diag(SumDiag))
    return(U)
  }
  else if(init_cw == "cw" ){}
  else {
    UU = list()
    for (u in 1:init_cw){
      U = lapply(1:n_cl, function(x) diag(runif(n)))
      SumDiag = Reduce("+",U)
      UU[[u]] = lapply(U, function(x) x/diag(SumDiag))
    }
    crit = lapply(UU, function(x) cw.rgcca(A =A, U = x, C = C,tau = tau, ncomp = ncomp,
                                           scheme = scheme ,scale = scale,init = init,
                                           bias = bias, tol = tol, verbose = verbose)$crit)
    return(U = UU[[which.max(crit)]])
  }
}

