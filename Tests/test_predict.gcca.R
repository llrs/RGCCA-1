rm(list = ls())


library(devtools)#to load package
library(CMA)#to generate learning sets
load_all("~/Documents/THESE/Implementation/RGCCA/.")

#######################################
# example 3: RGCCA  in MCCV with hand #  
#######################################

#Starts by loading data
data("Russett")
N       = nrow(Russett)
lab     = apply(Russett[, 9:11], 1, which.max)
Ytest   = matrix(0, N, 3)
X_agric = as.matrix(Russett[ , c("gini","farm","rent")])
X_ind   = as.matrix(Russett[ , c("gnpr","labo")])
X_polit = as.matrix(Russett[ , c("demostab", "dictatur")])
A       = list(X_agric, X_ind, X_polit)
#Define the design matrix (output = C) 
C            = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
result.rgcca = rgcca(A, C, tau = rep(1, 3), ncomp = rep(1, 3), scheme = "factorial", verbose = TRUE)

rate       = 0.8
niter      = 300 
FOLD       = GenerateLearningsets(n = N, method = "MCCV", niter = niter, ntrain = rate*N)@learnmatrix
nb_to_pred = N - dim(FOLD)[2]
Ytest      = matrix(0, (nb_to_pred)*dim(FOLD)[1], 3)
for (i in 1:nrow(FOLD)){
  ind  = seq(nb_to_pred*(i-1)+1, i*nb_to_pred, by = 1)
  B    = lapply(A, function(x) x[FOLD[i, ], ])
  B    = lapply(B, scale2)
  resB = rgcca(B, C, tau = rep(1, 3), scheme = "factorial", scale = FALSE, verbose = FALSE, scale_size_bloc = T)
  #  look for potential conflicting sign among components within the loo loop.
  for (k in 1:length(B)){
    if (cor(result.rgcca$a[[k]], resB$a[[k]]) >= 0) 
      resB$a[[k]] = resB$a[[k]] else resB$a[[k]] = -resB$a[[k]]
  }
  Btest         = lapply(A, function(x) x[-FOLD[i, ], ])
  Btest[[1]]    = (Btest[[1]] - t(replicate(nb_to_pred, attr(B[[1]],"scaled:center")))) / (t(replicate(nb_to_pred, attr(B[[1]],"scaled:scale")))) / sqrt(NCOL(B[[1]]))
  Btest[[2]]    = (Btest[[2]] - t(replicate(nb_to_pred, attr(B[[2]],"scaled:center")))) / (t(replicate(nb_to_pred, attr(B[[2]],"scaled:scale")))) / sqrt(NCOL(B[[2]]))
  Btest[[3]]    = (Btest[[3]] - t(replicate(nb_to_pred, attr(B[[3]],"scaled:center")))) / (t(replicate(nb_to_pred, attr(B[[3]],"scaled:scale")))) / sqrt(NCOL(B[[3]]))
  Ytest[ind, 1] = Btest[[1]]%*%resB$a[[1]]
  Ytest[ind, 2] = Btest[[2]]%*%resB$a[[2]]
  Ytest[ind, 3] = Btest[[3]]%*%resB$a[[3]]
}

rm(list=setdiff(ls(), c("Ytest", "FOLD", "N", "A", "C", "result.rgcca", "nb_to_pred", "lab")))
source('~/Documents/THESE/Implementation/RGCCA/R/predict.gcca.R')

library(devtools)
load_all("~/Documents/THESE/Implementation/RGCCA/.")

###############################################
# example 3: RGCCA  in MCCV with predict.gcca #  
###############################################

Ytest_newFunction = matrix(0, dim(Ytest)[1], dim(Ytest)[2])
for (i in 1:nrow(FOLD)){
  ind  = seq(nb_to_pred*(i-1)+1, i*nb_to_pred, by = 1)
  B    = lapply(A, function(x) x[FOLD[i, ], ])
  resB = rgcca(B, C, tau = rep(1, 3), scheme = "factorial", scale = T, verbose = FALSE, scale_size_bloc = T)
  #  look for potential conflicting sign among components within the loo loop.
  for (k in 1:length(B)){
    if (cor(result.rgcca$a[[k]], resB$a[[k]]) >= 0) 
    {resB$a[[k]] = resB$a[[k]]; resB$astar[[k]] = resB$astar[[k]]} else {resB$a[[k]] = -resB$a[[k]]; resB$astar[[k]] = -resB$astar[[k]]}
  }
  Btest                    = lapply(A, function(x) x[-FOLD[i, ], ])
  names(resB$A)            = names(Btest) = c("Agr", "Ind", "Pol")
  Ytest_newFunction[ind, ] = Reduce("cbind", predict.gcca(object = resB, newA = Btest, type = "classification", 
                                                          new_scaled = F, bloc_to_pred = "Pol", y.train = lab[FOLD[i, ]], y.test = lab[-FOLD[i, ]], 
                                                          scale_size_bloc = T, fit = "lda")$pred)
}

identical(Ytest_newFunction, Ytest)

#Agregate predictions
individuals = c(apply(FOLD, 1, function(x) (1:N)[-x]))
SUMUP       = apply(X = Ytest, 2, function(x) tapply(X = x, INDEX = individuals, FUN = mean))

#Plot results
data("Russett")
plot(result.rgcca$Y[[1]], result.rgcca$Y[[2]], col = "white", xlab = "Y1 (Agric. inequality)", ylab = "Y2 (Ind. Development)")
text(result.rgcca$Y[[1]], result.rgcca$Y[[2]], rownames(Russett), col = lab, cex = .7)
text(SUMUP[, 1], SUMUP[, 2], substr(rownames(Russett), 1, 1), col = lab, cex = .7)



