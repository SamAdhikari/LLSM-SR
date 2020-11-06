##load required packages
library(mvtnorm)
library(igraph)
library(Rcpp)
library(RcppArmadillo)

##source r scripts to fit the model
source('BetaupdateRWCOVSR.R')
source('InterceptupdateRWCOVSR.R')
source('llsmRWCOVSR.R')
source('logitInverseCovSR.R')
source('MCMCsampleRWCOVSR.R')
source('SRupdate.R')
source('ZupdateRWCOVSR.R')
source('Zprior.R')
source('betaprior.R')
source('getIndicesYY.R')
source('SigmaUpdateRW.R')
source('adjust.my.tune.R')
source('which.suck.R')
sourceCpp('Likelihoods.cpp')
source('updateMu.R')


##Function to generate networks
logitInverse = function(d,Xi,Beta,Si,Rj)
{
  return(1/(1+exp(d-sum(Beta*Xi)-Si-Rj)))
}

inverselogit_LSM = function(Beta,distMat,XX,S,R)
{
  D = distMat
  prob = array(NA,dim=c(nrow(D),nrow(D)))
  for(i in 1:nrow(D)){
    for(j in i:nrow(D)){
      prob[i,j] = logitInverse(D[i,j],XX[i,j,],Beta,S[i],R[j])
      prob[j,i] = logitInverse(D[j,i],XX[j,i,],Beta,S[j],R[i])
    }
  }
  return(prob)
}

genDD = function(ZZ)
{
  dd <- as.matrix(dist(ZZ,diag=TRUE,upper=TRUE))
  diag(dd) = 0
  return(dd)
}

getProb_LSM = function(distMat,Beta,XX,S,R)
{
  nn = dim(distMat)[1]
  if(is.null(XX)){
    XX = array(0,dim=c(nn,nn,1)) }
  Prob = inverselogit_LSM(Beta=Beta,distMat,XX,S,R)    
  return(Prob)
}

### Set parameters to generate longitudinal networks
TT = 3
n = 25
Sigma = 0.5
Beta = array(c(0.1,0,0.5,1,-1.5,0.1),dim=c(2,TT))
XX = list(length=TT)
for(tt in 1:TT){
  XX[[tt]] = array(NA,dim=c(n,n,2))
  XX[[tt]][,,1] =  t(sapply(1:n,function(x)rnorm(n,0,1)))
  XX[[tt]][,,2] = t(sapply(1:n,function(x)rnorm(n,0,1)))
}

Z = list(length=TT)
YY = list(length = TT)

Z[[1]] = t(sapply(1:n,function(x)rnorm(2,0,Sigma)))
for(kk in 2:(TT-1)){
  Z[[kk]] = Z[[kk-1]] + rnorm(2,0,Sigma)
}
Z[[TT]] = Z[[TT-1]] + rnorm(2,0,2*Sigma)

D = lapply(1:TT,function(x)genDD(ZZ=Z[[x]]))

S = rnorm(n,0.1,0.5)
R = rnorm(n,0,.1)

for(tt in 1:TT){
   Prob = getProb_LSM(Beta = Beta[,tt],distMat = D[[tt]],XX=XX[[tt]],S=S,R=R)
   YY[[tt]] = array(NA,dim=c(n,n))
   for(jj in 1:n){
  YY[[tt]][,jj] = rbinom(n,1,Prob[,jj])
   }
  diag(YY[[tt]]) = 0
  dimnames(YY[[tt]])[[1]] = dimnames(YY[[tt]])[[2]] = 1:n
} 


##plot the networks
par(mfrow=c(2,3)); lapply(1:TT,function(x)plot.igraph(graph.adjacency(YY[[x]])))


####################
#Fit the model

burnin = 500
thin = 10
niter = 10000

ptm = proc.time()
LLSMfitRWCOVSR =  llsmRWCOVSR(Y=YY,X=XX,initialVals = NULL, 
                              priors = NULL, tune = NULL, 
                              tuneIn = TRUE, dd=2, niter=5000)

ptm2 = proc.time() - ptm





