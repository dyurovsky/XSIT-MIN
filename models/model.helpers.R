# Computes BIC Given a LogLikelihood, #Paramaters, and Data Set Size
bic <- function(lp, params, n) {-2*lp +params*log(n)}

# Computes Strength in Memory for Statistical Accumulation Model
#
# Input
# gamma: Initial Encoding Strength Parameter
# lambda: Memory Decay Parameter
# nrefs: Number of Referents
# int: Number of Intervening Trials
#
# Output
# a list of (Same trial, Switch trial, New Label Trial)
accumulation.pred <- function(gamma,lambda,nrefs,int) {
  sn <- (gamma/nrefs)*int^-lambda
  on <- sn
  wn <- on/4
  
  return(c(sn,on,wn))
}

# Computes Strength in Memory for Single Referent Model
#
# Input
# gamma: Initial Encoding Strength Parameter
# lambda: Memory Decay Parameter
# nrefs: Number of Referents
# int: Number of Intervening Trials
#
# Output
# a list of (Same trial, Switch trial, New Label Trial)
singleref.pred <- function(gamma,lambda,nrefs,int) {
  sn <- gamma*int^-lambda
  on <- 0
  wn <- 0
  
  return(c(sn,on,wn))
}

# Computes Strength in Memory for Integrated Model
#
# Input
# sigma: Attention Allocation Parameter
# gamma: Initial Encoding Strength Parameter
# lambda: Memory Decay Parameter
# nrefs: Number of Referents
# int: Number of Intervening Trials
#
# Output
# a list of (Same trial, Switch trial, New Label Trial)
integrated.pred <- function(sigma,gamma,lambda,nrefs,int) {
  sn <- gamma*sigma*int^-lambda
  on <- gamma*((1-sigma)/(nrefs-1))*int^-lambda
  wn <- on/4
  
  return(c(sn,on,wn))
}

# Computes Prop. Correct at Test from Memory Strength
#
# Input
# preds: a list of(Same Trial encoding, Switch Trial encoding, 
#                     New Label Trial encoding)
# nrefs: Number of Referents
#
# Output
# a list of (Same trial, Switch trial, New Label Trial)
prop.from.pred <- function(preds,nrefs) {
  
  ps <- mean(rbinom(100000,1,preds[1] + (1-preds[1])/nrefs))
  po <- mean(rbinom(100000,1,preds[2] + (1-preds[2])/nrefs))
  wn <- mean(rbinom(100000,1,(1-preds[3])/nrefs) - preds[3])
  
  return(c(ps,po,wn))
}