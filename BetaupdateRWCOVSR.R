BetaupdateRWCOVSR <-
  function(Intercept,llikAll,MuBeta,VarBeta,tune,acc,Y,Z,TT,X,Beta,nn,dd,pp,SS,RR)  
{
    for(tt in 1:TT){
        BetaNew = Beta[,tt]
        RRtt = RR[which(RR[,1]%in%rownames(Z[[tt]])),2]
        SStt = SS[which(SS[,1]%in%rownames(Z[[tt]])),2]
        #print(RRtt)
        for(kk in 1:length(Beta[,tt])){
            BetaNew[kk] = Beta[kk,tt] + tune[kk,tt]*rnorm(1,0,1)
          #  print(BetaNew)
            #compute loglikelihood at proposed value
            #   loglik =likelihoodCovSR(Y=Y[[x]], Z=Z[[x]], intercept=Intercept, XX=X[[x]], Beta=Beta[,x],SS=SStt,RR=RRtt) 
            
        #    llikNew = likelihoodCovSR(Y=Y[[tt]], Z=Z[[tt]], 
            #                          intercept=Intercept, XX=X[[tt]], Beta=BetaNew,SS=SStt,RR=RRtt) 
            
          llikNew =  FullLogLikCOVSR(YY = Y[[tt]], ZZ = Z[[tt]],XX = X[[tt]],
                                     Beta = BetaNew, SS = SStt,RR = RRtt,
                                     intercept = Intercept,
                                     nn = nn[tt],dd = dd, pp = pp)
         #   print(llikNew)
         #  print(llikAll[tt])
            #log prior at current value
            priorOld = betaprior(Beta[kk,tt], MuBeta, VarBeta)
            #log prior at new value
            priorNew = betaprior(BetaNew[kk],MuBeta,VarBeta)
            #logratio
            logratio = llikNew - llikAll[tt] + priorNew - priorOld
            if(is.finite(logratio)){
            if(logratio > log(runif(1,0,1))){
                Beta[kk,tt] = BetaNew[kk]
                acc[kk,tt] = acc[kk,tt] + 1
                llikAll[tt] = llikNew
            }
            }else{BetaNew[kk]=Beta[kk,tt]}
        }
    }
        return(list(Beta=Beta,acc = acc, llik = llikAll))
}
