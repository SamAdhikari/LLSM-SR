MCMCsampleRWCOVSR <-
function(niter,Y,Z,X,Intercept,Beta,SS,RR,TT,dd,nn,pp,MuBeta, 
         MuInt,VarInt,VarZ,VarSS, VarRR, priorV,
                      accZ,accInt,accBeta,accSS,accRR,
         tuneBeta,tuneZ,tuneInt,tuneSS,tuneRR,A,B,gList)
{
    
    ZFinal = list()
    BetaFinal = array(NA,dim=c(pp,TT,niter))
    InterceptFinal = rep(NA,niter)
    Likelihood = rep(NA,niter)
    ZVarFinal = list()    
    SSFinal = array(NA,dim=c(dim(SS)[1],niter))
    RRFinal =array(NA,dim=c(dim(RR)[1],niter))
    llikOld = list()
    length(llikOld) = TT 
    RRpos = lapply(1:TT,function(x) which(RR[,1]%in%rownames(Z[[x]])))
    SSpos = lapply(1:TT,function(x) which(SS[,1]%in%rownames(Z[[x]])))
    
    for(iter in 1:niter){
        #update Z
        llikOld = lapply(1:TT,function(x){
            RRtt = RR[RRpos[[x]],2]
            SStt = SS[SSpos[[x]],2]
            return(sapply(1:nn[x],function(y){
                likelihoodiCOVSR(ii=y,dd=dd,nn=nn[x],pp=pp,Yt=Y[[x]],Xt=X[[x]],
                Zt=Z[[x]],SS=SStt,RR=RRtt,intercept=Intercept,Beta=Beta[,x])
            }))
            # likelihoodiCovSR(ii=y,Y=Y[[x]],XX=X[[x]],
            #                  Z=Z[[x]],intercept=Intercept,Beta=Beta[,x],SS=SStt,RR=RRtt)}))
        })
        Zupdt = ZupdateRWCOVSR(Y=Y,Z=Z,TT=TT,X=X,
                        Intercept=Intercept,Beta=Beta,dd=dd,
                        var=VarZ,llikOld=llikOld,acc=accZ,
                        tune=tuneZ,nn=nn,pp=pp,SS=SS,RR=RR)        
        Z = Zupdt$Z
        accZ = Zupdt$acc
     #  print(Z)
     SSupdt = SSupdate(SS=SS,RR=RR,YY=Y,ZZ=Z, XX=X,TT=TT,Beta=Beta,intercept=Intercept,tune=tuneSS,acc=accSS,
     nn=nn,pp=pp,dd=dd, priorVar = VarSS)
     SS = SSupdt$SS
     accSS = SSupdt$acc
     oldlikelihood = SSupdt$llik
     
     RRupdt = RRupdate(SS=SS,RR=RR,YY=Y,ZZ=Z,XX=X,TT=TT,Beta=Beta,intercept=Intercept,tune=tuneRR,acc=accRR,
     nn=nn,pp=pp,dd=dd,oldlikelihood=oldlikelihood, priorVar = VarRR)
     RR = RRupdt$RR
     accRR = RRupdt$acc
        #update beta 
	 llikBeta = sapply(1:TT,function(x){
	    RRtt = RR[RRpos[[x]],2]
	    SStt = SS[SSpos[[x]],2]
        loglik = FullLogLikCOVSR(YY=Y[[x]],ZZ=Z[[x]],XX=X[[x]],Beta=Beta[,x],SS=SStt,RR=RRtt,
	        intercept=Intercept,nn=nn[x],dd=dd,pp=pp)
            return(loglik)})
	
        Betaupdt = BetaupdateRWCOVSR(Intercept=Intercept,llikAll=llikBeta,
                                   MuBeta = MuBeta, VarBeta = VarInt,
                              tune = tuneBeta, acc = accBeta,Y=Y,Z=Z,TT=TT,
                              X=X,Beta=Beta,nn=nn,dd=dd,pp=pp,SS=SS,RR=RR)
        Beta = Betaupdt$Beta
        accBeta = Betaupdt$acc
        #update intercept
        llikAll = sum(Betaupdt$llik)
        #        llikAll = sum(sapply(1:TT,function(x){
        #    RRtt = RR[RRpos[[x]],2]
        #    SStt = SS[SSpos[[x]],2]
        #    loglik = FullLogLikCOVSR(YY=Y[[x]],ZZ=Z[[x]],XX=X[[x]],Beta=Beta[,x],SS=SStt,RR=RRtt,
        #                             intercept=Intercept,nn=nn[x],dd=dd,pp=pp)
        #    return(loglik)}))
        
         Intupdt = InterceptupdateRWCOVSR(Intercept=Intercept,llikAll=llikAll,
                            MuBeta=MuInt,VarBeta=VarInt,tune=tuneInt,
                           acc=accInt,Y=Y,Z=Z,TT=TT,X=X,Beta=Beta,nn=nn,dd=dd,pp=pp,SS=SS,RR=RR)

         Intercept = Intupdt$Intercept
         accInt = Intupdt$acc
         llikAll = Intupdt$llikAll
        #Update Sender and Receiver effects
         MuBeta = sapply(1:pp, function(x) updateMu(betaK = Beta[x, ], TT = TT, 
                                                   priorM = 0, priorV = priorV, SigmaSq = VarInt))
          
        #        VarZ = SigmaUpdateRW(Aprior=A,Bprior=B,Z=Z,nn=nn,dd=dd,TT=TT,gList=gList)
        #STORE UPDATES
        SSFinal[,iter] = SS[,2]
        RRFinal[,iter] = RR[,2]
        InterceptFinal[iter] = Intercept
        BetaFinal[,,iter] = Beta 
        ZFinal[[iter]] = Z
        Likelihood[iter] = llikAll
        ZVarFinal[[iter]] = VarZ	
        print(iter)
        
        ##save the chain after every 10,000 draws
        if(iter == 5000| iter == 10000 | iter == 20000 | iter == 30000 |iter == 60000){
            draws = list(Z=ZFinal,Intercept=InterceptFinal,Beta=BetaFinal,RR=RRFinal,SS=SSFinal,
                        Likelihood=Likelihood,VarZ=ZVarFinal )
            save(draws, file = paste('PartialChain',iter,'.RData',sep = ''))
        }
        
    }
    draws = list(Z=ZFinal,Intercept=InterceptFinal,Beta=BetaFinal,RR=RRFinal,SS=SSFinal,
                 Likelihood=Likelihood,VarZ=ZVarFinal )
    accZ = lapply(1:TT,function(x)accZ[[x]]/niter)
    accInt = accInt/niter
    accBeta = accBeta/niter
    accSS = accSS/niter
    accRR = accRR/niter
    acc = list(accZ=accZ,accInt=accInt,accBeta=accBeta, accSS = accSS, accRR=accRR)
    return(list(draws=draws,acc=acc))
}
