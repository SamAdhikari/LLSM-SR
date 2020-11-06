Zprior <-
function(Z,MuZ,VarZ)
{
    return(dmvnorm(Z,mean = MuZ, sigma = VarZ,log=TRUE))
}
