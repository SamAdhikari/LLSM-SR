logitInverseCOVSR <-
function(intercept,d,Xi,Beta,Si,Rj)
{
    return(1/(1+exp(d-intercept-sum(Beta*Xi)-Si-Rj)))
}
