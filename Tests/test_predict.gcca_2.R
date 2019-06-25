rm(list = ls())


library(devtools)#to load package
library(CMA)#to generate learning sets
load_all("~/Documents/THESE/Implementation/RGCCA/.")

#Starts by loading data
data("Russett")
N       = nrow(Russett)
lab     = apply(Russett[, 9:11], 1, which.max)
X_agric = as.matrix(Russett[ , c("gini","farm","rent")])
X_ind   = as.matrix(Russett[ , c("gnpr","labo")])
X_polit = as.matrix(Russett[ , c("demostab", "dictatur")])
A       = list(X_agric, X_ind, X_polit)
#Define the design matrix (output = C) 
C            = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), 3, 3)
result.rgcca = rgcca(A, C, tau = rep(1, 3), ncomp = rep(1, 3), scheme = "factorial", verbose = TRUE)

rate       = 0.6
niter      = 15 
# FOLD       = GenerateLearningsets(n = N, method = "MCCV", niter = niter, ntrain = rate*N)@learnmatrix
FOLD       = GenerateLearningsets(y = as.factor(lab), method = "CV", fold = 6, strat = T, niter = niter)@learnmatrix
nb_to_pred = N - dim(FOLD)[2]
Ytest      = matrix(0, (nb_to_pred)*dim(FOLD)[1], 3)
SCORE      = matrix(0, nrow(FOLD))


for (i in 1:nrow(FOLD)){
  ind  = seq(nb_to_pred*(i-1)+1, i*nb_to_pred, by = 1)
  B    = lapply(A, function(x) x[FOLD[i, ], ])
  resB = rgcca(B, C, tau = rep(1, 3), scheme = "factorial", scale = T, verbose = FALSE, scale_size_bloc = T)
  #  look for potential conflicting sign among components within the loo loop.
  for (k in 1:length(B)){
    if (cor(result.rgcca$a[[k]], resB$a[[k]]) >= 0) 
    {resB$a[[k]] = resB$a[[k]]; resB$astar[[k]] = resB$astar[[k]]} else {resB$a[[k]] = -resB$a[[k]]; resB$astar[[k]] = -resB$astar[[k]]}
  }
  Btest         = lapply(A, function(x) x[-FOLD[i, ], ])
  names(resB$A) = names(Btest) = c("Agr", "Ind", "Pol")
  PRED          = predict.gcca(object = resB, newA = Btest, type = "classification", new_scaled = F, bloc_to_pred = "Pol", 
                               y.train = lab[FOLD[i, ]], y.test = lab[-FOLD[i, ]], scale_size_bloc = T, fit = "lda")
  SCORE[i, ]   = PRED$score
}

mean(SCORE)
hist(SCORE)



