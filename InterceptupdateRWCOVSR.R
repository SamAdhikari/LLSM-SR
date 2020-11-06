InterceptupdateRWCOVSR <-
function(Intercept,llikAll,MuBeta,
            VarBeta,tune,acc,Y,Z,TT,X,Beta,nn,dd,pp,SS,RR)
{
    #propose new value for intercept
    IntNew = Intercept + tune*rnorm(1,0,1)
    #compute loglikelihood at proposed value
    llikNew =sum(sapply(1:TT,function(x){
        RRtt = RR[which(RR[,1]%in%rownames(Z[[x]])),2]
        SStt = SS[which(SS[,1]%in%rownames(Z[[x]])),2]
        #loglik =likelihoodCovSR(Y=Y[[x]],Z=Z[[x]],intercept=IntNew, XX=X[[x]],
        #                        Beta=Beta[,x],SS=SStt,RR=RRtt) 
        loglik = FullLogLikCOVSR(YY=Y[[x]],ZZ=Z[[x]],XX=X[[x]],Beta=Beta[,x],SS=SStt,RR=RRtt,
                        intercept=IntNew,nn=nn[x],dd=dd,pp=pp)
        return(loglik)}))
    #log prior at current value
    priorOld = betaprior(Intercept, MuBeta, VarBeta)
    #log prior at new value
    priorNew = betaprior(IntNew,MuBeta,VarBeta)
    #logratio
    logratio = llikNew - llikAll + priorNew - priorOld
    if(is.finite(logratio)){
    if(logratio > log(runif(1,0,1))){
        Intercept = IntNew
        acc = acc + 1
        llikAll = llikNew
    }
    }
    return(list(Intercept=Intercept,acc = acc,llikAll = llikAll))
}
