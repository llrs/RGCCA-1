# Fonction interne utilisee par sgcca.permute

sgcca.crit = function(A, C, c1s, ncomp, scheme, tol, crit = crit){
  for (k in 1:length(A)){
    A[[k]] = A[[k]][sample(1:nrow(A[[k]])),]
  }
  for (i in 1:NROW(c1s)){
    out = sgcca(A = A, C = C, c1 = c1s[i,], ncomp = ncomp, scheme = scheme, tol = tol)
    crit[i]=out$crit
  }
  return(crit)
}