par.cw.rgcca.hwang = function(A,
                       n_cl, 
                       init_cw = "random",
                       U_ref = NULL, 
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
  res = cw.rgcca.hwang(A = A,
                       n_cl = n_cl, 
                       init_cw = init_cw,
                       max_iter = max_iter, 
                       m = m, 
                       C = C, 
                       tau = tau, 
                       ncomp = ncomp, 
                       scheme = scheme, 
                       scale = scale, 
                       init = init,
                       bias = bias, 
                       tol = tol, 
                       verbose = verbose)
  if (!is.null(U_ref)){
    U1 = sapply(res$U,diag)
    M = cor(U_ref,U1)
    UU1 = U1
    for (i in 1:n_cl){
      UU1[,i] = U1[,which.max(M[,i])]
    }
    return(U = UU1)
  } else{
    U = res$U
    return(U = U)
  }
}