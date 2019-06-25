#
# ancillary.sgcca.surv.cv for potential parallel run
#
require(snow)
require(MASS)

ancillary.sgcca.lm.cv <- function(Xlist, C, ncomp, scale, method, y,
                                  trainmat_ind, penalty, scheme="horst")
{
  # cat("Processing: ", head(trainmat_ind), "... for ", as.numeric(penalty), "\n")
  #sanity
  num_block <- length(Xlist)
  if (method=="LOOCV") {
    xtest <- lapply(Xlist, function(mm) t(as.matrix(mm[-trainmat_ind,])))
  } else {
    xtest <- lapply(Xlist, function(mm) { as.matrix(mm[-trainmat_ind,])})
    
  }
  if (scale) {
    
    xxlist.train <- lapply(Xlist, function(mm) scale2(mm[trainmat_ind, ]))
    scl_fun <- function(data, scaled) {
      scale(data, center = attr(scaled, "scaled:center"),
            scale = attr(scaled, "scaled:scale")) }
    xxlist.test <- mapply(scl_fun, xtest, xxlist.train, SIMPLIFY=FALSE)
  } else {
    xxlist.train <- lapply(Xlist, function(mm) mm[trainmat_ind, ])
    xxlist.test <- lapply(Xlist, function(mm) mm[-trainmat_ind, ])
  }
  
  sgcca.res <- sgcca(A = xxlist.train, C = C, c1 = penalty, ncomp = ncomp, scheme=scheme)
  DataTrain <- as.data.frame(Reduce("cbind", sgcca.res$Y))
  colnames(DataTrain) <- c(paste("Y", 1:(ncol(DataTrain) - 1), sep = ""), "y")
  # DataTrain <- data.frame(DataTrain, y=sgcca.res$Y[[num_block]][, 1])
  reslm <- lm(y~.,data=DataTrain)
  
  xxtest <- mapply(function(data, a){ data%*%a}, xxlist.test, sgcca.res$a, SIMPLIFY=F)

  DataTest <- as.data.frame(Reduce("cbind", xxtest))
  colnames(DataTest) <- c(paste("Y", 1:(ncol(DataTest) - 1), sep = ""), "y")
  # DataTest <- data.frame(DataTest, y=xxtest[[num_block]][, 1])
  
  ychapo <- predict(reslm, DataTest)
  if (any(is.na(ychapo))) warning("NA in predictions.")
  mse <- mean((y[-trainmat_ind] - ychapo)**2)
  return(list(mse=mse))
}

#
# sgcca.lm.cv
#
#' cross validation for the sGCCA-Survival Analysis.
#'
#' Functions used to compute the optimal parameters used in a context of sGCCA
#' multiblock classification. The main function call an ancillary function
#' that performs the most internal CV loop
#'
#' @param Xlist List of the blocks. As for now, the last item of the list is the
#'              block corresponding to the class (Factor y describing the classes)
#' @param C     Matrix describing the design, supposed to be symetrical
#' @param y     Dependent variable
#' @param nfold Integer, number of folds for the CV frame. Default is the number of n of
#'              measurements (which is the LOO framework)
#' @param ncomp Integer vector, vector of number of components, should be in
#'              1:(blocks column numbers) for each block
#' @param scale Boolean, should the data be scaled internally? Default to FALSE.
#' @param params Data.frame containing the tabulated values of penalizing
#'              parameters to assess. A Default value is arbitrary.
#' @param cl    cluster descriptor. If the snow package is installed and
#'              if cl is not NULL, the CV computations will be performed on
#'              the cluster  described by cl. Default is NULL
#' @return optim vector of numeric with optimal penalty values.
#' @examples \dontrun{
#' library(gliomaData)
#' library(sgcca)
#' #read data
#' data(ge_cgh_locIGR)
#' Xlist <- ge_cgh_locIGR$multiblocks
#' #set design
#' C <- matrix(c(0,0,1,0,0,1,1,1,0),3,3)
#' # gridded parameters
#' P <-  expand.grid(c11 = c(3, 4, 5), c12 = c(3, 4, 5), c13 = 10)
#' # no cluster
#' cl <- NULL
#' opt <- sgcca.lm.cv(Xlist=Xlist,C,ncomp=c(1, 1, 1), scale=T, params=P, nfold=5, cl)
#' #
#' #
#' library(gliomaData)
#' library(sgcca)
#' #read data
#' data(ge_cgh_locIGR)
#' Xlist <- ge_cgh_locIGR$multiblocks
#' #set design
#' C <- matrix(c(0,0,1,0,0,1,1,1,0),3,3)
#' # gridded parameters
#' P <-  expand.grid(c11 = c(3, 4, 5), c12 = c(3, 4, 5), c13 = 10)/10
#' # no cluster
#' cl <- NULL
#' opt <- sgcca.lm.cv(Xlist=Xlist, C, ncomp=c(1, 1, 1), scale=T, params=P, nfold=5, cl)
#' #
#' # in case you get a bunch of available nodes
#' Host209425 <-list(host = "is209425",
#'          snowlib = "/tmp/usr/x86_64-pc-linux-gnu-library/2.15")
#' Host209426 <-list(host = "localhost")
#' cl <- makeCluster(c(rep(list(Options209425), 6),
#'                     rep(list(Options209426), 4)), type = "SOCK")
#' clusterEvalQ(cl, {library(sgcca)})
#' opt <- sgcca.lm.cv(Xlist=Xlist,C,ncomp=c(1, 1, 1), scale=T, params=P, nfold=5, cl)
#' }
#' @export sgcca.lm.cv
sgcca.lm.cv <- function (Xlist, C,  y, ncomp=rep(1, length(Xlist)),
                           nfold = nrow(Xlist[[1]]), scale = FALSE,
                           params=expand.grid(c11 = c(0.2, 0.3, 0.4),
                                              c12 = c(0.2, 0.3, 0.4),
                                              c13 = 1),
                           scheme="horst",
                           cl=NULL)
{
  num_samples <- nrow(Xlist[[1]])
  out <- y
  
  if (nfold == num_samples) {
    trainmat <- GenerateLearningsets(num_samples, method = "LOOCV")@learnmatrix
    attr(trainmat, "method") <- "LOOCV"
  } else {
    trainmat <- GenerateLearningsets(num_samples, method = "CV", fold = nfold)@learnmatrix
    attr(trainmat, "method") <- "CV"
  }
  
  #flatten the loop on folds and the internloop on grided parameters
  # explicitely create the list of all configurations of CV and GParam
  tmp = expand.grid(fold=1:nfold, grid_param=1:nrow(params))
  flattened_folds <- trainmat[tmp[,'fold'],]
  flattened_folds <- split(flattened_folds,1:nrow(flattened_folds) )
  flattened_grid_params <- params[tmp[,'grid_param'], ]
  flattened_grid_params <- split(flattened_grid_params, 1:nrow(flattened_grid_params))
  
  if (is.null(cl)) {
    pb <- txtProgressBar(style = 3)
    obj = mapply(ancillary.sgcca.lm.cv,  # apply this function
                 trainmat_ind=flattened_folds,#for this fold and par
                 penalty=flattened_grid_params,
                 MoreArgs=list(Xlist, C, ncomp, scale, #on this
                               y=out, method=attr(trainmat,"method"), scheme=scheme)
    )
    rescv <- unlist(obj)
    close(pb)
  } else {
    tm = snow.time(
      obj <- clusterMap(cl,ancillary.sgcca.lm.cv,  # apply this function
                        trainmat_ind=flattened_folds,#for this fold and par
                        penalty=flattened_grid_params,
                        MoreArgs=list(Xlist, C, ncomp, scale, #on this
                                      y=out, method=attr(trainmat,"method"), scheme=scheme)
      )
    )
    print(tm)
    rescv <- unlist(obj) # clusterMap does not SIMPLIFY unlike mapply
  }
  # reshape : the lines are the folds and the col are the diff gridded params
  rescv = matrix(rescv, nrow=nfold)
  return(opt = unlist(params[which.min(colSums(rescv, na.rm=T)), ]))
}

