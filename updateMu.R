#sample mean
updateMu = 
function(betaK, TT, priorM, priorV, SigmaSq ){
  denom =  (SigmaSq + TT*priorV) / priorV*SigmaSq
  postM = (priorM/priorV + sum(betaK)/SigmaSq)/ denom
  postV = 1/denom
  muK = rnorm(1, mean = postM, sd = sqrt(postV))
}
