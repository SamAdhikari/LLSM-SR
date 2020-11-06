#' @title Function to run MCMC sampler for the LLSM-RW model
#'
#' @description
#' \code{llsmRWCOV} runs MCMC sampler for the LLSM-RW model.
#'
#' @details
#' \code{llsmRW} runs MCMC sampler for the LLSM-RW model and returns samples from the posteriors chains of the parameters,
#' the posterior likelihood at the accpeted parameters, a list of acceptance rates from the metropolis hastings sampling, 
#' and a list of the tuning values if \code{tuneIn} is set to TRUE
#'
#' @param Y A list of sociomatrix for observed networks
#' @param initialVals A list of values for initializing the chain for \code{intercept} and \code{ZZ}. Default is set to NULL, 
#' when random initialization is used. 
#' @param priors A list of parameters for prior distribution specified as \code{MuBeta}, \code{VarBeta}, \code{VarZ}, 
#' \code{A} and \code{B}
#' If set to NULL, default priors is used
#' @param tune A list of tuning parameters. If set to NULL, default values are used.
#' @param tuneIn Logical option to specify whether to auto tune the chain or not. Default is \code{TRUE}
#' @param dd Dimension of the latent space
#' @param niter Number of MCMC iterations to run
#'
#' @aliases llsmRW
#' @export


llsmRWCOVSR <-
function(Y,X=NULL,initialVals = NULL, priors = NULL, tune = NULL, 
                      tuneIn = TRUE, dd, niter)
{    
    nn = sapply(1:length(Y),function(x) nrow(Y[[x]]))
    TT = length(Y) #number of time steps
    if(is.null(X)){
	pp = 1}else(pp = dim(X[[1]])[3])
    gList = getIndicesYY(Y,TT,nn)$gg
    
    C = lapply(1:TT,function(tt){
        diag(nn[tt]) - (1/nn[tt]) * array(1, dim = c(nn[tt],nn[tt]))})
    Z0 = lapply(1:TT,function(tt){
        g = graph.adjacency(Y[[tt]]);
        ss = shortest.paths(g);
        ss[ss > 4] = 4;
        Z0 = cmdscale(ss,k = dd);
        dimnames(Z0)[[1]] = dimnames(Y[[tt]])[[1]];
        return(Z0)})
    Z00 = lapply(1:TT,function(tt)C[[tt]]%*%Z0[[tt]])
 
    print(Z00)
    
    if(is.null(X)){
        XX= lapply(1:TT,function(x)array(0,dim=c(nn[x],nn[x],1))) 
     }else{
	    XX = list()
	    for(tt in 1:TT){
		XX[[tt]] = array(0,dim=c(nn[tt]*pp,nn[tt]))
		a = 1
		b = nn[tt]
		for(ll in 1:pp){
		    XX[[tt]][a:b,] = X[[tt]][,,ll]
		    a = b + 1
		    b = b + nn[tt]
 		} }    }

    #Priors
    if(is.null(priors)){
        MuInt= 0 
        VarInt = 1000
        VarZ = diag(150,dd)
        A = 100
        B = 150
        VarSS = 1
        VarRR = 1
        priorV = 1000
    }else{
        if(class(priors) != 'list')(stop("priors must be of class list, if not NULL"))
        MuInt = priors$MuBeta
        VarInt = priors$VarBeta
        VarZ = priors$VarZ
        A = priors$A
        B = priors$B
        VarSS = priors$VarSS
        VarRR = priors$VarRR
        priorV = priors$priorV
    }
    
    ##starting values
    if(is.null(initialVals)){
#         Z0 = list()
#         for(i in 1:TT){  
#             #     ZZ = t(replicate(nn[i],rnorm(dd,0,1)))
#             ZZ = array(NA,dim=c(nn[i],dd))    
#             Z0[[i]] = ZZ		 
#         }
        Intercept0  = rnorm(1,0,1)
        # Intercept0 = 0
        Beta0 = sapply(1:TT,function(x)rnorm(pp,0,1))
        if(pp == 1){
            Beta0 = t(matrix(Beta0))
        }
        MuBeta = rnorm(pp,0,1)
        uniqueName = unique(unlist(lapply(1:TT,function(tt) dimnames(Y[[tt]])[[1]])))
        SS0 =data.frame('Names'=(uniqueName),
                        'SS'= rnorm(length(uniqueName),0,1) )
        RR0 =data.frame('Names'=uniqueName,
                        'RR'= rnorm(length(uniqueName),0,1) )
        #  print(SS0)
        print("Starting Values Set")
    }else{
        if(class(initialVals)!= 'list')(stop("initialVals must be of class list, if not NULL"))
        Z0 = initialVals$Z
        Intercept0 = initialVals$intercept
        Beta0 = initialVals$Beta
        SS0 = initialVals$SS0
        RR0 = initialVals$RR0
        MuBeta = initialVals$MuBeta
        }
   
    ###tuning parameters#####
    if(is.null(tune)){
        a.number = 5
        tuneInt = 1
        tuneBeta = array(1,dim=c(pp,TT))
        tuneSS = tuneRR = rep(1,dim(RR0)[1])
	if(pp == 1){tuneBeta = t(matrix(tuneBeta))}
        tuneZ =  lapply(1:TT, function(x) rep(1.2,nn[x]))          
    } else{
        if(class(tune) != 'list')(stop("tune must be of class list, if not NULL"))
        a.number = 1
        tuneInt = tune$tuneInt
        tuneBeta = tune$tuneBeta
        tuneZ = tune$tuneZ
        tuneSS = tune$tuneSS
        tuneRR = tune$tuneRR
    }
    accZ = lapply(1:TT,function(x)rep(0,nn[x]))
    accInt = 0
    accBeta = array(0,dim=c(pp,TT))
    accSS = accRR = rep(0,dim(RR0)[1])
        ###tuning the Sampler####
    do.again = 1
    tuneX = 1
    if(tuneIn == TRUE){
        #  while(do.again ==1){
            print('Tuning the Sampler')
            for(counter in 1:a.number ){                
                rslt = MCMCsampleRWCOVSR(niter = 700,Y=Y, Z = Z0, X = XX,
                                         Intercept = Intercept0,
                                  Beta = Beta0, SS = SS0, RR = RR0, 
                                  TT = TT, dd = dd, nn = nn, pp = pp,
                                  MuBeta = MuBeta,
                                  MuInt = MuInt, VarInt = VarInt,
                                  VarZ = VarZ, VarSS = VarSS, VarRR = VarRR,
                                  priorV = priorV,
                                  accZ = accZ, accInt = accInt,
                                  accBeta = accBeta, accSS = accSS, accRR = accRR,
                                  tuneBeta = tuneBeta, tuneRR = tuneRR, tuneSS = tuneSS,
                                  tuneZ = tuneZ, tuneInt = tuneInt, 
                                  A = A, B = B, gList=gList)
                tuneZ = lapply(1:TT,function(x)adjust.my.tune(tuneZ[[x]],
                                                              rslt$acc$accZ[[x]],2))
               tuneInt = adjust.my.tune(tuneInt,rslt$acc$accInt, 1)
                tuneSS = adjust.my.tune(tuneSS,rslt$acc$accSS,1)
                tuneRR = adjust.my.tune(tuneRR,rslt$acc$accRR,1)
                tuneBeta = sapply(1:TT,function(x){sapply(1:pp,function(y){
			adjust.my.tune(tuneBeta[y,x],rslt$acc$accBeta[y,x],1)})})
		if(pp ==1){tuneBeta = t(matrix(tuneBeta))}
                print(paste('TuneDone = ',tuneX))
                tuneX = tuneX+1
            }
            # extreme = lapply(1:TT,function(x)which.suck(rslt$acc$accZ[[x]],2))
            # do.again = max(sapply(extreme, length)) > 50
            # }
        print("Tuning is finished")
        save(rslt,file = 'MCMCfit_tuning.RData')
    }
    
    rslt = MCMCsampleRWCOVSR(niter = niter,Y=Y,Z=Z0,X=XX,Intercept=Intercept0,
                                       Beta=Beta0,SS=SS0, RR=RR0,
                                      TT=TT,dd=dd,nn=nn,pp=pp,
                                       MuBeta = MuBeta,
                                       MuInt=MuInt,VarInt=VarInt,
                                       VarZ=VarZ,VarSS = VarSS, VarRR = VarRR,
                                       priorV = priorV,
                                       accZ=accZ,accInt=accInt,
                                       accBeta=accBeta,accSS=accSS,accRR=accRR,
                                       tuneBeta=tuneBeta,tuneRR=tuneRR,tuneSS=tuneSS,
                                       tuneZ=tuneZ,tuneInt=tuneInt,
                             A=A,B=B,gList=gList) 
    ##Procrustean transformation of latent positions
    #     #######################################

    Ztransformed = lapply(1:niter, function(ii) {lapply(1:TT,
                        function(tt){z= rslt$draws$Z[[ii]][[tt]];
                                 z = C[[tt]]%*%z;
                                 pr = t(Z00[[tt]])%*% z;
                                 ssZ = svd(pr);
                                 tx = ssZ$v%*%t(ssZ$u);
                                 zfinal = z%*%tx;
                                 dimnames(zfinal)[[1]] = dimnames(rslt$draws$Z[[ii]][[tt]])[[1]]
                                 return(zfinal)})})    
    rslt$draws$ZZ = Ztransformed
    rslt$call = match.call()
    rslt$tune = list(tuneZ = tuneZ, tuneInt = tuneInt,tuneBeta=tuneBeta,
                     tuneSS = tuneSS,tuneRR=tuneRR)
    class(rslt) = 'LLSM'
    rslt       
}
