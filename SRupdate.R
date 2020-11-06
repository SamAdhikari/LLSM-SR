#Update sender and receiver effect for each individual 
#by Metropolis Hastings algorithm
SSupdate = function(SS,RR,YY,ZZ, XX,TT,Beta,intercept,tune,acc,nn,pp,dd, priorVar){
    RRtt = list()
    nameList = list()
    for(tt in 1:TT){
        nameList[[tt]] = rownames(ZZ[[tt]])
        RRtt[[tt]] = RR[which(RR[,1]%in%nameList[[tt]]), 2]
    }
    NewSS = SS
    
    oldlikelihood = rep(NA,dim(SS)[1])
    for(ii in 1:dim(SS)[1]){
        oldlikelihood[ii] = sum(sapply(1:TT,function(tt){
            kk = which(rownames(ZZ[[tt]]) == SS[ii,1])
            if(length(kk)>0){
                SStt = SS[which(SS[,1]%in%rownames(ZZ[[tt]])),2]
                return(likelihoodiCOVSR(ii=kk,dd=dd,nn=nn[tt],
                pp=pp,Yt=YY[[tt]],Xt=XX[[tt]],
                Zt=ZZ[[tt]],SS=SStt,RR=RRtt[[tt]],
                intercept=intercept,Beta=Beta[,tt]))
            }else(return(0)) }) )
    }
    
 #   print(NewSS)
 #   print(SS)
    for(ii in 1:dim(SS)[1]){
        NewSS[ii,2] = SS[ii,2] + tune[ii]*rnorm(1,0,1)
        # oldlikelihood = sum(sapply(1:TT,function(tt){
        #    SStt = SS[which(SS[,1]%in%rownames(ZZ[[tt]])),2]
        #    kk = which( rownames(ZZ[[tt]]) == SS[ii,1])
        #    if(length(kk)>0){return(likelihoodiCOVSR(ii=kk,dd=dd,nn=nn[tt],pp=pp,Yt=YY[[tt]],Xt=XX[[tt]],
        #        Zt=ZZ[[tt]],SS=SStt,RR=RRtt[[tt]],intercept=intercept,Beta=Beta[,tt]))
        #        }else(return(0)) }) )
        
        oldprior = dnorm(SS[ii,2], 0, sqrt(priorVar), log=TRUE)
    
        newlikelihood =  sum(sapply(1:TT,function(tt){
            
            kk = which(rownames(ZZ[[tt]])== SS[ii,1])
            if(length(kk)>0){
                NewSStt = NewSS[which(NewSS[,1]%in%rownames(ZZ[[tt]])),2]
              return(likelihoodiCOVSR(ii=kk,dd=dd,nn=nn[tt],pp=pp,Yt=YY[[tt]],Xt=XX[[tt]],
                                      Zt=ZZ[[tt]],SS=NewSStt,RR=RRtt[[tt]],intercept=intercept,Beta=Beta[,tt]))
              
               }else(return(0))}))
        newprior = dnorm(NewSS[ii,2], 0, sqrt(priorVar), log=TRUE)
    
        logRatio = newlikelihood - oldlikelihood[ii] + newprior - oldprior
        if(is.finite(logRatio)){
        if(logRatio > log(runif(1))){
            SS[ii,2] = NewSS[ii,2]
            oldlikelihood[ii] = newlikelihood
            acc[ii] = acc[ii]+1
        }
        }else{
            NewSS[ii,2] = SS[ii,2]        
        }
    }
 #   print(SS)
    return(list(SS=SS,acc=acc,llik = oldlikelihood))
 }

RRupdate = function(SS,RR,YY,ZZ,XX,TT,Beta,intercept,tune,acc,nn,pp,dd,oldlikelihood, priorVar)
  {
    SStt = list()
    for(tt in 1:TT){
        SStt[[tt]] = SS[which(SS[,1] %in% rownames(ZZ[[tt]])),2]
    }
    NewRR = RR
    for(ii in 1:length(RR[,2]))
      {
        NewRR[ii,2] = RR[ii,2] +tune[ii]*rnorm(1,0,1)
        
        # oldlikelihood = sum(sapply(1:TT,function(tt){
        #    RRtt = RR[which(RR[,1]%in%rownames(ZZ[[tt]])),2]
        #     kk = which(rownames(ZZ[[tt]])== RR[ii,1])
        #     if(length(kk)>0){
        #              return(likelihoodiCOVSR(ii=kk,dd=dd,nn=nn[tt],pp=pp,Yt=YY[[tt]],Xt=XX[[tt]],
        #                              Zt=ZZ[[tt]],SS=SStt[[tt]],RR=RRtt,intercept=intercept,Beta=Beta[,tt]))
        #      }else(return(0))}))
        
        oldprior = dnorm(RR[ii,2], 0, sqrt(priorVar), log=TRUE)
        
        newlikelihood =  sum(sapply(1:TT,function(tt){
            kk = which(rownames(ZZ[[tt]])==RR[ii,1])
            if(length(kk)>0){
            NewRRtt = NewRR[which(NewRR[,1]%in%rownames(ZZ[[tt]])),2]
            return(likelihoodiCOVSR(ii=kk,dd=dd,nn=nn[tt],pp=pp,Yt=YY[[tt]],Xt=XX[[tt]],
                                      Zt=ZZ[[tt]],SS=SStt[[tt]],RR=NewRRtt,intercept=intercept,Beta=Beta[,tt]))
        }else(return(0))}))
        
        newprior = dnorm(NewRR[ii,2], 0, sqrt(priorVar), log=TRUE)
      #  print(newlikelihood)
      #  print(oldlikelihood)
        
            logRatio = newlikelihood - oldlikelihood[ii] + newprior - oldprior
            if(is.finite(logRatio)){
            if(logRatio > log(runif(1))){
                RR[ii,2] = NewRR[ii,2]
                oldlikelihood[ii] = newlikelihood
                acc[ii] = acc[ii] +1
            }
            }else{
                NewRR[ii,2] = RR[ii,2] 
            }
            }
        return(list(RR=RR,acc=acc))
    }
        
