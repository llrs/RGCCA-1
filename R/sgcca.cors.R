# Fonction interne utilisee par sgcca.permute

sgcca.cors = function(A, C, c1s, ncomp, scheme, tol){
  cors = rep(0,NROW(c1s))
  for (k in 1:length(A)){
    A[[k]] = A[[k]][sample(1:nrow(A[[k]])),]
  }
  for (i in 1:NROW(c1s)){
    out = sgcca(A = A, C = C, c1 = c1s[i,], ncomp = ncomp, scheme = scheme, tol = tol)
    Y=Reduce("cbind",out$Y)
    cors[i] = sum(C*cor(Y))
  }
  return(cors)
}