require(methods)
#' learningsets  class for support in sgcca package.
#'
#' Some details about this class and my plans for it in the body.
#'
# @slot coucou A matrix of dimension \code{niter x ntrain}
#'
#' \describe{
#'     \item{\code{learnmatrix}:}{ A matrix of dimension \code{niter x
#'    ntrain}. Each row contains the indices of those observations
#'       representing the learningset for one iteration. If \code{method =
#'    CV}, zeros appear due to rounding issues.}
#'     \item{\code{method}:}{The method used to generate the \code{learnmatrix},
#'                           s.\code{\link{GenerateLearningsets}}}
#'     \item{\code{ntrain}:}{Number of observations in one learning set.If
#'       \code{method = CV}, this number is not attained for all
#'       iterations, due to rounding issues.}
#'     \item{\code{iter}:}{Number of iterations (different learningsets)
#'       that are stored in \code{learnmatrix}.}
#'   }
#' @name learningsets-class
#' @rdname learningsets-class
#' @exportClass learningsets
setClass(Class="learningsets", 
         representation(learnmatrix="matrix", method="character", ntrain="numeric",
                        iter="numeric")
         )

roundvector <- function(x, maxint){
  fx <- floor(x)
  aftercomma <- x-fx
  roundorder <- order(aftercomma, decreasing=TRUE)
  i <- 1
  while(sum(fx) < maxint){ 
    fx[roundorder[i]] <- ceiling(x[roundorder[i]])
    i <- i+1
  }
  return(fx)
}

rowswaps <- function(blocklist){
  
  cols <- length(blocklist)
  fold <- nrow(blocklist[[1]])
  learnmatrix <- blocklist[[1]]
  for(i in 2:cols) learnmatrix <- cbind(learnmatrix, blocklist[[i]])
  rs <- rowSums( learnmatrix == 0)
  allowedzeros <- ceiling(sum(rs)/fold)
  indmatrix <-  matrix(rep(1:fold, each=cols), nrow=fold, byrow=TRUE) 
  while(any(rs > allowedzeros)){
    indmatrix <- replicate(cols, sample(1:fold))
    temp2list <- blocklist
    for(i in 1:cols) temp2list[[i]] <- blocklist[[i]][indmatrix[,i], ]
    learnmatrix <- temp2list[[1]]
    for(i in 2:cols) learnmatrix <- cbind(learnmatrix, temp2list[[i]])
    rs <- rowSums( learnmatrix == 0)
  }
  return(indmatrix)
}


#' Generate a learningsets manager a helper for the CV, LOOCV, etc processings
#' 
#' Create an object that interact with various CV schemes. It deals with
#' multiblocks objects.
#'
#' @param n numeric, number of folds chosen
#' @param y factoc vector to create balanced os stratified CV subsets
#' @param method character in c("LOOCV", "CV", "MCCV", "bootstrap")
#' @param fold Integer
#' @param niter Integer
#' @param ntrain Integer
#' @param strat boolean for stratification
#  @return a learningsets-class object
#' @references Martin Slawski, Anne-Laure Boulesteix and Christoph Bernau. (2009). CMA: Synthesis of microarray-based classification. R package version 1.14.0.
#' @export GenerateLearningsets
GenerateLearningsets <- function (n, y, method = c("LOOCV", "CV", "MCCV", "bootstrap"), 
          fold = NULL, niter = NULL, ntrain = NULL, strat = FALSE) 
{
  if (!missing(n)) {
    if (length(n) != 1 || n < 0) 
      stop("'n' must be a positive integer ! \n")
    n <- as.integer(n)
    if (!is.null(fold) && n <= fold) 
      stop("'n' is too small \n")
    if (!is.null(ntrain) && n <= ntrain) 
      stop("'n' is too small \n")
  }
  if (missing(n) & missing(y)) 
    stop("At least one of 'n' or 'y' mus be given \n")
  if (!missing(y)) 
    n <- length(y)
  method <- match.arg(method, c("LOOCV", "CV", "MCCV", "bootstrap"))
  if (!is.element(method, eval(formals(GenerateLearningsets)$method))) 
    stop("method must be one of 'LOOCV', 'CV', 'MCCV', 'bootstrap' \n")
  if (strat & missing(y)) 
    stop("If 'strat=TRUE', 'y' (class memberships) must be given \n")
  if (method == "MCCV") {
    if (is.null(niter) | is.null(ntrain)) 
      stop("With the MCCV method, arguments niter and ntrain should be given.")
    if (strat) {
      taby <- table(y)
      prop <- taby/sum(taby)
      classize <- roundvector(prop * ntrain, ntrain)
      if (any(classize < 1)) 
        stop("Generation of learningsets failed, one or several classes are too small. \n")
      indlist <- sapply(names(taby), function(z) which(y == 
        z), simplify = FALSE)
      learnmatrix <- matrix(nrow = niter, ncol = ntrain)
      lower <- cumsum(c(1, classize[-length(classize)]))
      upper <- cumsum(classize)
      for (i in 1:length(indlist)) learnmatrix[, lower[i]:upper[i]] <- t(replicate(niter, 
                                                                                   sample(indlist[[i]], classize[i], replace = FALSE)))
    }
    else learnmatrix <- t(replicate(niter, sample(n, ntrain, 
                                                  replace = FALSE)))
  }
  if (method == "CV") {
    if (is.null(niter)) 
      niter <- 1
    if (is.null(fold)) 
      stop("With the CV method, argument 'fold' must be given.")
    if (!strat) {
      if (fold == n) 
        method <- "LOOCV"
      else {
        size <- n/fold
        learnmatrix <- matrix(0, niter * fold, n - floor(size))
        size.int <- floor(size)
        size.vector <- rep(size.int, fold)
        if (size.int != size) 
          size.vector[1:((size - size.int) * fold)] <- size.vector[1:((size - 
            size.int) * fold)] + 1
        group.index <- c()
        for (j in 1:fold) group.index <- c(group.index, 
                                           rep(j, size.vector[j]))
        for (i in 1:niter) {
          group.index <- group.index[sample(n, n, replace = FALSE)]
          for (j in 1:fold) {
            whichj <- which(group.index == j)
            learnmatrix[j + (i - 1) * fold, 1:length(whichj)] <- whichj
          }
        }
        learnmatrix <- learnmatrix[, 1:max(size.vector), 
                                   drop = FALSE]
        if (size.int != size) 
          learnmatrix <- t(apply(learnmatrix, 1, function(z) setdiff(0:n, 
                                                                     z)))
        if (size.int == size) 
          learnmatrix <- t(apply(learnmatrix, 1, function(z) setdiff(1:n, 
                                                                     z)))
      }
    }
    else {
      taby <- table(y)
      prop <- taby/sum(taby)
      siz <- n - floor(n/fold)
      classize <- roundvector(prop * siz, siz)
      if (any(taby < fold)) 
        stop("Generation of learningsets failed, one or several classes are smaller than the number of folds. \n")
      indlist <- sapply(names(taby), function(z) which(y == 
        z), simplify = FALSE)
      templist <- vector(mode = "list", length = length(indlist))
      for (i in 1:length(indlist)) {
        outp <- do.call(GenerateLearningsets, args = list(n = taby[i], 
                                                          method = "CV", niter = niter, fold = fold))@learnmatrix
        templist[[i]] <- t(apply(outp, 1, function(z) ifelse(z == 
          0, 0, indlist[[i]][z])))
      }
      topass <- lapply(templist, function(z) z[1:fold, 
                                               , drop = FALSE])
      swaporder <- rowswaps(topass)
      nrep <- 1
      while (nrep < niter) {
        swaporder <- rbind(swaporder, swaporder[1:fold, 
                                                , drop = FALSE] + fold * nrep)
        nrep <- nrep + 1
      }
      for (i in 1:length(templist)) templist[[i]] <- templist[[i]][swaporder[, 
                                                                             i], ]
      learnmatrix <- templist[[1]]
      for (i in 2:length(indlist)) learnmatrix <- cbind(learnmatrix, 
                                                        templist[[i]])
    }
  }
  if (method == "LOOCV") 
    learnmatrix <- matrix(rep(1:n, each = n - 1), nrow = n)
  if (method == "bootstrap") {
    if (is.null(niter)) 
      stop("If 'method=bootstrap', the argument 'niter' must be given. \n")
    if (!strat) 
      learnmatrix <- t(replicate(niter, sample(n, replace = TRUE)))
    else {
      taby <- table(y)
      if (any(taby) < 1) 
        stop("Generation of learningsets failed, one or several classes are too small. \n")
      indlist <- sapply(names(taby), function(z) which(y == 
        z), simplify = FALSE)
      learnmatrix <- matrix(nrow = niter, ncol = n)
      lower <- cumsum(c(1, taby[-length(taby)]))
      upper <- cumsum(taby)
      for (i in 1:length(indlist)) {
        learnmatrix[, lower[i]:upper[i]] <- t(replicate(niter, 
                                                        sample(indlist[[i]], taby[i], replace = TRUE)))
      }
    }
  }
  if (strat & is.element(method, c("CV", "MCCV", "bootstrap"))) 
    method <- paste("stratified", method)
  new("learningsets", learnmatrix=learnmatrix, method=method,
      ntrain=ncol(learnmatrix), iter=nrow(learnmatrix))
  #return(list(learnmatrix = learnmatrix, method = method,
  #            ntrain = ncol(learnmatrix), iter = nrow(learnmatrix)))
}

  