#### Pour prendre l'intersection des individus des blocs consideres, en retirant les outliers et les batch effects si specifies



rgcca.entry = function(X,rownames = NULL, labels = NULL, batch = NULL, outliers = NULL) {
  J = length(X)
  A = list()
  if (is.null(rownames)){
    rownames = rownames(X[[1]])
    for (j in 2:J){
      rownames = intersect(rownames,rownames(X[[j]]))
    }
  } else {rownames = rownames}
  if (!is.null(outliers)){
    rownames = setdiff(rownames,outliers)
  } 
  A = lapply(X, function(x) x[rownames,])
  if (!is.null(labels)){
    if (!is.null(batch)){
      B = A
      modcombat <- model.matrix(~1, data =  A[[J]])
      for (k in batch){
        b = labels[rownames,k]
        for (j in 1:J){
          B[[j]] = t(B[[j]])
          B[[j]] <- ComBat(dat = B[[j]], batch = b, mod = modcombat, par.prior = TRUE, prior.plots = FALSE)
          B[[j]] = t(B[[j]])
        }
      }
    } else {B = NULL}
  } else {B = NULL}
  return(list(Blocks = A, Blocks.corr = B, rownames = rownames)) 
}