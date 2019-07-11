#' Run through a set of constraint parameters c1s to select the best with permutation
#' 
#' @param A A list that contains the J blocks of variables \mathbf{X_1}, \mathbf{X_2}, ..., \mathbf{X_J}
#' @param c1s A matrix containing sets of constraint variables, one row by set. If null, sgcca.permute takes 10 sets between min values ($1/sqrt(ncol)$) and 1
#' @param nperm Number of permutation tested for each set of constraint
#' @param C A design matrix that describes the relationships between blocks (default: complete design)
#' @param scheme The value is "horst", "factorial" or "centroid" (default: "centroid")
#' @param plot A logical, should a plot of coeffi
#' @value A list containing :
#' @value \item {pval}
#' @value \item {zstat}
#' @value \item {bestpenalties}
#' @value \item {permcrit} 
#' @value \item {crit}

sgcca.permute.cors = function(A,
                         c1s = NULL,
                         nperm = 20,
                         C = 1 - diag(length(A)),
                         ncomp = rep(1,length(A)),
                         scheme = "centroid",
                         plot = FALSE,
                         tol = .Machine$double.eps){
  if (is.null(c1s)) {
    c1s <- matrix(NA, nrow = 10, ncol = length(A))
    for (k in 1:length(A)) {
      c1s[, k] <- pmin(1/(seq(0.8, 0.1, len = 10) * 
                                  sqrt(ncol(A[[k]]))), 1)
    }
  }
  crits = rep(0,NROW(c1s))
  for (i in 1:NROW(c1s)){
    out = sgcca(A, C, c1 = c1s[i,], ncomp = ncomp, scheme = scheme)
    crits[i] = sum(C*cor(Reduce("cbind",out$Y)))
  }
  if (Sys.info()["sysname"] == "Windows"){
    n_cores = parallel::detectCores()
    e=environment()
    cl = parallel::makeCluster(n_cores)
    parallel::clusterExport(cl,c("A","c1s","nperm","C","ncomp","scheme","out","crit","crits","tol"),envir = e)
    parallel::clusterEvalQ(cl, library(devtools))
    tryCatch({parallel::clusterEvalQ(cl, load_all("RGCCA/R/."))},
             error = function(e) {warning("error : probably an issue with the localisation of RGCCA functions")})
    permcrit = tryCatch({parallel::parSapply(cl,1:nperm, 
                                   function(x) sgcca.cors(A = A, 
                                                          C = C, 
                                                          c1s = c1s, 
                                                          ncomp = ncomp, 
                                                          scheme = scheme, 
                                                          tol = tol))},
                        error = function(e) {warning("an error occured with sgcca.crit"); return(NULL) })
    parallel::stopCluster(cl)
    if (is.null(permcrit)){return(NULL)} 
  } else {
    permcrit = parallel::mlcapply(1:nperm, 
                                   function(x) sgcca.cors(A = A, 
                                                          C = C, 
                                                          c1s = c1s, 
                                                          ncomp = ncomp, 
                                                          scheme = scheme, 
                                                          tol = tol))
  }
  pvals = zs = matrix(NA,nrow = NROW(c1s), ncol = NCOL(c1s)+1)
  # chaque ligne : [critere 1, critere 2, valeur pvals ou zs]
  for (i in 1:NROW(c1s)){
      pvals[i,] = c(c1s[i,],mean(permcrit[i,]>=crits[i]))
      zs[i,] = c(c1s[i,],(crits[i] - mean(permcrit[i,])) / (sd(permcrit[i,])))# + 0.05))
  }
  
  # zs <- c(zs, (cors[i] - mean(permcors[, i]))/(sd(permcors[,i]) + 0.05))
  
  bestpenalties = c1s[which.max(zs[,length(A)+1]),1:length(A)]
  sgcca.best = sgcca(A, C = C, c1 = bestpenalties, ncomp = ncomp, scheme = scheme)
  out = list(pvals = pvals, zstat = zs , 
             bestpenalties =  bestpenalties,
             sgcca.best = sgcca.best,
             permcrit = permcrit, crit = crits )
  if (plot){
    par(mfrow=c(2,1), mar = c(2,1,2,1))
    plot(1:NROW(c1s),out$crit, ylim = c(0, max(crits,permcrit)))
    for (i in 1:nperm){
      points(1:NROW(c1s),out$permcrit[,i], col = "green")
    }
    for (j in 1:(NROW(c1s)-1)){
      segments(j,out$crit[j],j+1,out$crit[j+1])
    }
    
    plot(1:NROW(c1s),out$zstat[,length(A)+1])
    for (j in 1:(NROW(c1s)-1)){
      segments(j,out$zstat[j,length(A)+1],j+1,out$zstat[j+1,length(A)+1])
    }
    
  }
  return(out)
}
