# Generate 3 data sets so that first 25 features are correlated across
# the data sets...

library(PMA)
u <- matrix(rnorm(50),ncol=1)
v1 <- matrix(c(rep(.5,25),rep(0,75)),ncol=1)
v2 <- matrix(c(rep(1,25),rep(0,25)),ncol=1)
v3 <- matrix(c(rep(.5,25),rep(0,175)),ncol=1)

x1 <- u%*%t(v1) + matrix(rnorm(50*100),ncol=100)
x2 <- u%*%t(v2) + matrix(rnorm(50*50),ncol=50)
x3 <- u%*%t(v3) + matrix(rnorm(50*200),ncol=200)

xlist <- lapply(list(x1, x2, x3),scale2)

# Run MultiCCA.permute w/o specifying values of tuning parameters to
# try.
# The function will choose the lambda for the ordered data set.
# Then permutations will be used to select optimal sum(abs(w)) for
# standard data sets.
# We assume that x1 is standard, x2 is ordered, x3 is standard:

set.seed(987654321)
perm.out <- MultiCCA.permute(xlist, type=c("standard"),nperm = 30)#, "ordered", "standard")
                                          
print(perm.out)
plot(perm.out)
out.cca <- MultiCCA(xlist, type=c("standard"),
                penalty=perm.out$bestpenalties, ncomponents=1, ws=perm.out$ws.init)
print(out.cca)
# Or if you want to specify tuning parameters by hand:
# this time, assume all data sets are standard:
# perm.out <- MultiCCA.permute(xlist, type="standard",
                             # penalties=cbind(c(1.1,1.1,1.1),c(2,3,4),c(5,7,10)), ws=perm.out$ws.init)
# print(perm.out)
# plot(perm.out)

# Making use of the fact that the features are ordered:
# out <- MultiCCA(xlist, type="ordered", penalty= perm.out$bestpenalties)
par(mfrow=c(3,1))
PlotCGH(out.cca$ws[[1]], chrom=rep(1,ncol(x1)))
PlotCGH(out.cca$ws[[2]], chrom=rep(2,ncol(x2)))
PlotCGH(out.cca$ws[[3]], chrom=rep(3,ncol(x3)))


### Sélection 2D
A = xlist
c1s = perm.out$penalties
for (i in 1:nrow(c1s)){ 
  c1s[i,] = c1s[i,] / sqrt(ncol(xlist[[i]]))
}
nperm = 30
C = 1-diag(3)
# C = matrix(c(0,0,1,0,0,1,1,1,0),3,3)

perm.sgcca = sgcca.permute(xlist,c1s = t(c1s),C = C,nperm = nperm,scheme = "horst",plot = TRUE)

out.sgcca = sgcca(A, C = C, c1 = perm.sgcca$bestpenalties, scheme = "horst")
par(mfrow=c(3,1))
PlotCGH(out.sgcca$a[[1]], chrom=rep(1,ncol(x1)))
PlotCGH(out.sgcca$a[[2]], chrom=rep(2,ncol(x2)))
PlotCGH(out.sgcca$a[[3]], chrom=rep(3,ncol(x3)))

out.sgcca$a
################## 
# Comparaison pas complétement valable car sggca utilise cov() pour le critére et cca utilise cor()





####################  Construction d'un dataset de données corrélées:
library(mvtnorm)
library(MASS)
sigma = matrix(c(1,0.2,0.8,
                 0.2,1,0.7,
                 0.8,0.7,1),3,3, byrow = TRUE)

X = mvrnorm(50, rep(0,3),Sigma = sigma)

v1 <- matrix(c(rep(.5,25),rep(0,75)),ncol=1)
v2 <- matrix(c(rep(1,25),rep(0,25)),ncol=1)
v3 <- matrix(c(rep(.5,25),rep(0,175)),ncol=1)


x1 <- X[,1]%*%t(v1) + 0.5*matrix(rnorm(50*100),ncol=100)
x2 <- X[,2]%*%t(v2) + 0.5*matrix(rnorm(50*50),ncol=50)
x3 <- X[,3]%*%t(v3) + 0.5*matrix(rnorm(50*200),ncol=200)

xlist <- lapply(list(x1, x2, x3),scale2)
set.seed(987654321)
nperm = 50

## Selection via PMA package, crit = correlation
perm.cca <- MultiCCA.permute(xlist, type=c("standard"),nperm = nperm)#, "ordered", "standard")

print(perm.cca)
plot(perm.cca)

out.cca <- MultiCCA(xlist, type=c("standard"),
                    penalty=perm.cca$bestpenalties, ncomponents=1, ws=perm.out$ws.init)
par(mfrow=c(3,1), mar = rep(2,4))
PlotCGH(out.cca$ws[[1]], chrom=rep(1,ncol(x1)))
PlotCGH(out.cca$ws[[2]], chrom=rep(2,ncol(x2)))
PlotCGH(out.cca$ws[[3]], chrom=rep(3,ncol(x3)))



# c1s = perm.out$penalties
# for (i in 1:nrow(c1s)){ 
#   c1s[i,] = c1s[i,] / sqrt(ncol(A[[i]]))
# }
C = matrix(c(0,1,1,
             1,0,1,
             1,1,0),3,3)

t1 = Sys.time()
perm.sgcca = sgcca.permute(xlist,C = C,nperm = nperm,scheme = "centroid",plot = TRUE, tol = 1e-5)
T2 = Sys.time()-t1

out.sgcca = sgcca(xlist, C = C, c1 = perm.sgcca$bestpenalties, scheme = "centroid")
par(mfrow=c(3,1), mar = rep(2,4))
PlotCGH(out.sgcca$a[[1]], chrom=rep(1,ncol(x1)))
PlotCGH(out.sgcca$a[[2]], chrom=rep(2,ncol(x2)))
PlotCGH(out.sgcca$a[[3]], chrom=rep(3,ncol(x3)))



tm = microbenchmark(test_cca = MultiCCA(xlist, type=c("standard"),penalty = perm.cca$bestpenalties),
                    test_sgcca = rgcca(xlist,C = C,scheme = "centroid", tol = 1e-5,verbose = FALSE,scale_size_bloc = F),
                    times = 20)


