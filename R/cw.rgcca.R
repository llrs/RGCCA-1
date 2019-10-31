# Internal function. Computes clusterwise RGCCA by building one model per cluster, cluster defined by diagonal matrix \eqn{\mathbf{U_c}} for c from 1 to n_cl (number of cluster), and \eqn{\sum_cu_{ci}=1}.

cw.rgcca = function(A,
                    U =list(diag(1,nrow(A[[1]]))), 
                    m=2, 
                    C = 1 - diag(length(A)), 
                    tau = rep(1, length(A)), 
                    ncomp = rep(1,length(A)), 
                    scheme = "centroid", 
                    scale = TRUE, 
                    init = "svd",
                    bias = TRUE, 
                    tol = 1e-08, 
                    verbose = TRUE){
  fit.rgcca = lapply(U, function(x) rgcca(A = lapply(A, function(z) x**(m/2)%*%z),C = C,tau = tau,ncomp = ncomp,
                                          scheme = scheme,scale = scale,init = init,bias = bias,tol = tol,verbose = verbose))
  return(list(fit.rgcca = fit.rgcca, crit_cl = sum(sapply(fit.rgcca, function(x) x$crit[length(x$crit)]))))
}
