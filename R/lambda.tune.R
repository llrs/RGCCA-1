lambda.tune <- function(A.train = A.train, A.test = A.test, C, i, method="CV", metric = "auc", type = "classification",
                        scaled = T, scaled_size_bloc = T, tau_Y = 0, lambda.grid){
    # x is a list of B (p_k by n) matrices of predictors
    # y dataframe of survival times and censoring, 
    # D is a B by B binary design matrix describing connections between blocks
    # trainmat matrix of training sets indices
    # outer_it outerloop iteration
    # lambda.grid is a matrix containing all combinations of lambdas over blocks to test
    # TODO : propose several different metrics : inner_cv_metric=c("auc", "F", "wrec")
    
    cat("### CV iteration : ", i, "\n")
    
    library(MASS)
    library(CMA)
    library(survival)
    #  library(multiblox)
    source("ALS_multiblox_NewtRaph.R")
    source("multiblox.predict.R")
    source("multiblox.score.R")
    source("functions.R")
    
    N <- nrow(A.train[[1]]) # nb of individuals
    
    beta.train <- eta.train <- eta.test <- NULL
    pred.score <- res <- model <- cox.model <- cv <- NULL
    
    ### Beta initialization without any intercept
    beta0 <- NULL
    B     <- length(A.train)
    for (b in 1:B){
      beta0[[b]] <- matrix(2, nrow=(ncol(A.train[[b]])), ncol=1)
    }
    
    ### pathwise warm restart              
    for (l in 1:nrow(lambda.grid)){
      cat("lambdas : ", l, " -> ")
      for (bid in 1:length(lambda.grid[l,])) {
        cat(" ", lambda.grid[l,bid], sep="")
      }
      cat("\n")
      
      ##2) Estimating the betas on the training set x, I, R, D, lambda.opt = 0, eps = 0.001, max.iter = 10000, beta.init = NULL
      res[[l]]<- rgcca(A = A.train, C = C, tau = c(lambda.grid[l,], tau_Y), ncomp = rep(1, length(A.train)), scheme = scheme, 
                       scale = !scaled, scale_size_bloc = !scaled_size_bloc , init = "svd", bias = T, tol = 1e-8, verbose = F)

      beta.train <- res[[l]]$a
      
      if(metric=="deviance"){
        ### deviance of partial likelihood
        cv[[l]] <- predict.gcca(object = res[[l]], newA = A.test, type = type, new_scaled = !scaled, bloc_to_pred = bloc_to_pred, 
                                scale_size_bloc = !scaled_size_bloc, fit = metric)
      }else if(metric=="exactConcordanceIndex"){
        ### predict
        pred <- multiblox.predict(model.train=res[[l]], x.train=X.train, y.train=y.train, x.test=X.test) 
        
        ##5) get and accumulate the score
        pred.score[[l]] <- multiblox.score2(pred$ychapo, as.matrix(y.test), metric=metric)$perf
      }
      
      beta0 <- res[[l]]$beta # update beta0 for warm restart
    }
    return(list(res=res, pred.score=pred.score, cv=cv, metric=metric))
  }
