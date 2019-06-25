#
# ancillary.sgcca.cor.cv for potential parallel run
#

require(snow)
require(MASS)
ancillary.sgcca.cor.cv <- function(Xlist, C, ncomp, scale, method, yLabel,
                                   trainmat_ind, penalty, scheme="horst")
{
  get.upper.tri <- function(m) m[upper.tri(m)]
  
  require(MASS)
  require(sgcca)
  cat("Processing: ", head(trainmat_ind), "... for ", as.numeric(penalty), "\n")
  #sanity
  num_block <- length(Xlist)
  if (method=="LOOCV") {
    xtest <- lapply(Xlist, function(mm) t(as.matrix(mm[-trainmat_ind,])))
  } else {
    xtest <- lapply(Xlist, function(mm) { as.matrix(mm[-trainmat_ind,])})
  }
  if (scale) {
    xxlist.train <- lapply(Xlist, function(mm) scale2(mm[trainmat_ind, ],scale=FALSE))
    scl_fun <- function(data, scaled) {
      scale(data, center = attr(scaled, "scaled:center"))}#,
##### Caution : code modified            scale = attr(scaled, "scaled:scale")) }
    #       xxlist.train <- lapply(Xlist, function(mm) myscale(mm[trainmat_ind, ]))
    #       scl_fun <- function(data, scaled) {
    #                 scale(data, center = attr(scaled, "scaled:center"),
    #                             scale = attr(scaled, "std")) }
    
    xxlist.test <- mapply(scl_fun, xtest, xxlist.train, SIMPLIFY=FALSE)
  } else {
    xxlist.train <- lapply(Xlist, function(mm) mm[trainmat_ind, ])
    xxlist.test <- lapply(Xlist, function(mm) mm[-trainmat_ind, ])
  }
  sgcca.res <- sgcca(A = xxlist.train, C = C, c1 = penalty, ncomp = ncomp,scheme=scheme,
                     scale=FALSE)
  xxtest <- mapply(function(data, a){ data%*%a}, xxlist.test, sgcca.res$astar)
  if (method=="LOOCV") {
    return( list(xxtest=xxtest) )
  }else{
    sumcor<- sum(get.upper.tri(cor(xxtest)*C)**2)
    return( list(sumcor=sumcor) )
  }
}

#
# sgcca.cor.cv
#
#' cross validation for the sGCCA using correlations to evaluate models
#'
#' Functions used to compute the optimal parameters used in a context of sGCCA
#' multiblock classification. The main function call an ancillary function
#' that performs the most internal CV loop
#'
#' @param Xlist List of the blocks. As for now, the last item of the list is the
#'              block corresponding to the class (Factor y describing the classes)
#' @param C     Matrix describing the design, supposed to be symetrical
#' @param nfold Integer, number of folds for the CV frame. Default is the number of n of
#'              measurements (which is the LOO framework)
#' @param ncomp Integer vector, vector of number of components, should be in
#'              1:(blocks column numbers) for each block
#' @param scale Boolean, should the data be scaled internally? Default to FALSE.
#' @param params Data.frame containing the tabulated values of penalizing
#'              parameters to assess. A Default value is arbitrary.
#' @param scheme Character ("horst", "factorial", "centroid")
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
#' opt <- sgcca.cor.cv(Xlist=Xlist,C,ncomp=c(1, 1, 1), scale=T, params=P, nfold=5, cl)
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
#' P <-  expand.grid(c11 = c(3, 4, 5), c12 = c(3, 4, 5), c13 = 10)
#' # no cluster
#' cl <- NULL
#' opt <- sgcca.cor.cv(Xlist=Xlist,C,ncomp=c(1, 1, 1), scale=T, params=P, nfold=5, cl)
#' #
#' # in case you get a bunch of available nodes
#' Host209425 <-list(host = "is209425",
#'          snowlib = "/tmp/usr/x86_64-pc-linux-gnu-library/2.15")
#' Host209426 <-list(host = "localhost")
#' cl <- makeCluster(c(rep(list(Options209425), 6),
#'                     rep(list(Options209426), 4)), type = "SOCK")
#' clusterEvalQ(cl, {library(sgcca)})
#' opt <- sgcca.cor.cv(Xlist=Xlist,C,ncomp=c(1, 1, 1), scale=T, params=P, nfold=5, cl)
#' }
#' @export sgcca.cor.cv
sgcca.cor.cv <- function (Xlist, C,  ncomp=rep(1, length(Xlist)),
                          nfold = nrow(Xlist[[1]]), scale = FALSE,
                          params=expand.grid(c11 = c(3, 4, 5),
                                             c12 = c(3, 4, 5),
                                             c13 = 10),
                          scheme="horst",
                          cl=NULL)
{
  num_samples <- nrow(Xlist[[1]])
  outcome <- tail(Xlist, 1)[[1]]
  for(k in 1:ncol(outcome)) outcome[, k] <- k*outcome[, k]
  out <- apply(outcome, 1, sum)
  
  if (nfold == num_samples) {
    trainmat <- GenerateLearningsets(num_samples, out, method = "LOOCV")@learnmatrix
    attr(trainmat, "method") <- "LOOCV"
    get.upper.tri <- function(m) m[upper.tri(m)]
  } else {
    trainmat <- GenerateLearningsets(num_samples, out, method = "CV", fold = nfold)@learnmatrix
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
    obj = mapply(ancillary.sgcca.cor.cv,  # apply this function
                 trainmat_ind=flattened_folds,#for this fold and par
                 penalty=flattened_grid_params,
                 MoreArgs=list(Xlist, C, ncomp, scale, #on this
                               yLabel=out, method=attr(trainmat,"method"), scheme=scheme)
    )
    if (attr(trainmat,"method") == "LOOCV"){
      obj_all    = data.frame(Reduce("cbind", obj))
      indx_split = factor(rep(1:nrow(params), each = nfold))
      var_params = split(names(obj_all), indx_split)
      obj_params = lapply(var_params, function(x) obj_all[, x])
      rescv      = lapply(obj_params, function(x) sum(get.upper.tri(cor(t(x))*C)**2))
      rescv      = unlist(rescv)
    }else{
      rescv <- unlist(obj)
    }
    close(pb)
  } else {
    tm = snow.time(
      obj <- clusterMap(cl,ancillary.sgcca.cor.cv,  # apply this function
                        trainmat_ind=flattened_folds,#for this fold and par
                        penalty=flattened_grid_params,
                        MoreArgs=list(Xlist, C, ncomp, scale, #on this
                                      yLabel=out, method=attr(trainmat,"method"), scheme=scheme)
      )
    )
    print(tm)
    
    if (attr(trainmat,"method") == "LOOCV"){
      obj_all    = data.frame(Reduce("cbind", obj))
      indx_split = factor(rep(1:nrow(params), each = nfold))
      var_params = split(names(obj_all), indx_split)
      obj_params = lapply(var_params, function(x) obj_all[, x])
      rescv      = lapply(obj_params, function(x) sum(get.upper.tri(cor(t(x))*C)**2))
      rescv      = unlist(rescv)
    }else{
      rescv <- unlist(obj) # clusterMap does not SIMPLIFY unlike mapply
    }
  }
  # reshape : the lines are the folds and the col are the diff gridded params
  if (attr(trainmat,"method") == "LOOCV"){
    print(rescv)
    return(opt = unlist(params[which.max(rescv), ]))
  }else{
    rescv <- matrix(rescv, nrow=nfold)
    print(rescv)
    return(opt = unlist(params[which.max(colSums(rescv, na.rm=T)), ])) 
  }
}
