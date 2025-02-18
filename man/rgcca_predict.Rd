% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rgcca_predict.R
\name{rgcca_predict}
\alias{rgcca_predict}
\title{Predict RGCCA}
\usage{
rgcca_predict(
  rgcca_res,
  newA,
  model = "regression",
  fit = "lm",
  bloc_to_pred = NULL,
  new_scaled = TRUE,
  regress_on = "block"
)
}
\arguments{
\item{rgcca_res}{A fitted RGCCA object (see  \code{\link[RGCCA]{rgcca}})}

\item{newA}{A list of either a dataframe/matrix or a vector giving the blocks to be predicted}

\item{model}{A character corresponding to the model of prediction among : regression or classification}

\item{fit}{A character giving the function used to compare the trained and the tested models}

\item{bloc_to_pred}{A character giving the block to predicted (must be the same name among train and test set)}

\item{new_scaled}{A boolean scaling the blocks to predict}

\item{regress_on}{A boolean indicating if the regression is performed on the blocks (by default) or on the components}
}
\description{
Predict a new block from a RGCCA
}
\examples{
data("Russett")
blocks = list(
agriculture = Russett[, 1:3],
industry = Russett[, 4:5],
politic = Russett[, 6:11]
)
C = connection = matrix(c(0, 0, 1,
0, 0, 1,
1, 1, 0),
3, 3)
object1 = rgcca(blocks, connection = C, tau = c(0.7,0.8,0.7),
    ncomp = c(3,2,4), superblock = FALSE, response = 3)
A = lapply(object1$call$blocks, function(x) x[1:32,])
object = rgcca(A, connection = C, tau = c(0.7,0.8,0.7),
    ncomp = c(3,2,4), scale = FALSE, scale_block = FALSE, superblock = FALSE, response = 3)
newA = lapply(object1$call$blocks, function(x) x[-c(1:32),])
newA = lapply( newA, function(x) x[, sample(1:NCOL(x))] )
newA = sample(newA, length(newA))
bloc_to_pred = "industry"
to_pred_train = kmeans(A[[bloc_to_pred]], 3)$cluster
to_pred_test = kmeans(newA[[bloc_to_pred]], 3)$cluster
political_regime = factor(apply(Russett[, 9:11], 1, which.max),
labels =c("demostab", "demoinst", "dictator"))
res  = rgcca_predict(object, newA, bloc_to_pred = "industry")
# ( res  = rgcca_predict(object, newA, "regression", "cor", "industry") )
# rgcca_predict(object, newA, bloc_to_pred = "industry", fit = "r2")$score # 0.5586036
res  = rgcca_predict(object, newA)
library(MASS)
}
