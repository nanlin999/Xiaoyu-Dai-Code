library(Rcpp)
# mcf: calculate MCF value for each test
sourceCpp("/Users/xiaoyudai/Documents/Paper2/Rcode/mcf.cpp")


##############################################################################


##############################################################################
BH <- function(p.org, alpha = 0.1){
  # alpha: nominal FDR level
  
  # number of test
  n = length(p.org)
  
  # search for s
  s = alpha
  while(n*s/max(1,sum(p.org<=s)) > alpha){
    s = s - 0.0001
  }
  # print(n*s/max(1,sum(p.org<=s)))
  # index of rejected tests
  return(which(p.org <= s))
}
##############################################################################

# Storey-BH method
##############################################################################
SBH <- function(p.org, alpha, lambda){
  # alpha:      nominal FDR level
  # lambda: threshold to estimate pi0
  
  # search for s
  s = alpha
  while(s/(1-lambda)*(1+sum(p.org>lambda))/max(1,sum(p.org<=s)) > alpha){
    s = s - 0.0001
  }
  # print(s/(1-lambda)*(1+sum(p.org>lambda))/max(1,sum(p.org<=s)))
  # index of rejected tests
  return(which(p.org <= s))
}
##############################################################################

# Selective SeqStep (Barber & Candes, 2015)
##############################################################################
SS <- function(p.org, alpha, s){
  # alpha: nominal FDR level
  # s: threshold for each test
  
  k = length(p.org)
  while(s/(1-s)*(1+sum(p.org[1:k]>s))/max(1,sum(p.org[1:k]<s)) > alpha){
    k = k - 1
    if(k<0){
      stop('No suitable k found since s is too large!')
    }
  }
  # print(s/(1-s)*(1+sum(p.org[1:k]>s))/max(1,sum(p.org[1:k]<s)))
  return(which(p.org[1:k]<s))
}
##############################################################################

# Adaptive SeqStep (Barber & Candes, 2015)
##############################################################################
AS <- function(p.org, alpha, s, lambda){
  # alpha:      nominal FDR level
  # s:      threshold for each test
  # lambda: threshold to estimate pi0
  
  k = length(p.org)
  while(s/(1-lambda)*(1+sum(p.org[1:k]>lambda))/max(1,sum(p.org[1:k]<s)) > alpha){
    k = k - 1
    if(k<0){
      stop('No suitable k found since s is too large!')
    }
  }
  # print(s/(1-lambda)*(1+sum(p.org[1:k]>lambda))/max(1,sum(p.org[1:k]<s))) # FDR.hat
  return(which(p.org[1:k]<s))
}
##############################################################################


##############################################################################
MCF <- function(p.org, p.next, alpha, s = 0.1*q, lambda = 0.5, method = 'BH'){
  # general methods combined with MCF modification to adopt for discrete tests
  # alpha: nominal FDR level
  # s: threhold for each test (needed for 'SS', 'AS')
  # lambda: threshold to estimate pi0 (needed for 'SBH', 'SS')
  
  # possible method types
  if(!method %in% c('BH', 'SBH', 'SS', 'AS')){
    stop("unrecognized method")
  }
  
  if(method == 'BH'){
    rej = BH.MCF(p.org, p.next, alpha)
  } else if (method == 'SBH'){
    rej = SBH.MCF(p.org, p.next, alpha, lambda)
  } else if (method == 'SS'){
    rej = SS.MCF(p.org, p.next, alpha, s)
  } else {
    rej = AS.MCF(p.org, p.next, alpha, s, lambda)
  }
  
  return(rej)
}
##############################################################################


##############################################################################
BH.MCF <- function(p.org, p.next, alpha){
  # samples of randomized p-value
  n.test = length(p.org)
  n.rep <- 1000
  randp.sample <- rep(0,n.rep*n.test)
  for(i in 0:(n.rep-1)){
    randp = runif(n.test,p.next,p.org)
    randp.sample[(i*n.test+1):((i+1)*n.test)] <- randp
  }
  # determin s
  s = alpha
  while(s/(sum(randp.sample<=s)/(n.test*n.rep)) > alpha){
    s = s - 0.0001
    
  }
  # print(s/(sum(randp.sample<=s)/(n.test*n.rep)))
  # reject based on MCF
  R = round(sum(randp.sample<=s)/n.rep) # number of rejections
  # calculating MCF r's
  r <- mcf(p.org, p.next, s)
  rej = sort(r, decreasing = T, index.return = T)$ix[1:R]
  return(rej)
}
##############################################################################


##############################################################################
SBH.MCF <- function(p.org, p.next, alpha, lambda){
  # samples of randomized p-value
  n.test = length(p.org)
  n.rep <- 100
  randp.sample <- rep(0,n.rep*n.test)
  for(i in 0:(n.rep-1)){
    randp = runif(n.test,p.next,p.org)
    randp.sample[(i*n.test+1):((i+1)*n.test)] <- randp
  }
  # determin s
  s = alpha
  while(s/(1-lambda) * (1+sum(randp.sample>lambda)/n.rep) / (sum(randp.sample<=s)/n.rep) > alpha){
    s = s - 0.0001
  }
  # print(s/(1-lambda) * (1+sum(randp.sample>lambda)/n.rep) / (sum(randp.sample<=s)/n.rep))
  # reject based on MCF
  R = round(sum(randp.sample<=s)/n.rep) # number of rejections
  # calculating MCF r's
  r <- mcf(p.org, p.next, s)
  rej = sort(r, decreasing = T, index.return = T)$ix[1:R]
  return(rej)
}
##############################################################################


##############################################################################
SS.MCF <- function(p.org, p.next, alpha, s){
  # fdr estimator
  n = length(p.org)
  n.rep = 10
  randp = matrix(rep(0,n.rep*n), nrow = n)
  for(i in 1:n.rep){
    randp[,i] = runif(n, p.next, p.org)
  }
  
  fdr.est <- function(p.org, p.next, s, k){
    randp.sample = as.vector(randp[1:k,])
    A = sum(randp.sample > s)/n.rep
    R = max(sum(randp.sample <= s)/n.rep, 1)
    fdr.est = s/(1-s)*(1+A)/R
    return(fdr.est)
  }
  
  # determine k
  k = length(p.org)
  fdr.tmp = fdr.est(p.org, p.next, s, k)
  while(fdr.tmp > alpha){
    k = k-1 
    fdr.tmp = fdr.est(p.org, p.next, s, k)
    # print(c(k,fdr.tmp))
    if(k<0){
      stop('No suitable k found since s is too large!')
    }
  }
  # determine rejetions
  randp.sample = as.vector(randp[1:k,])
  R = round(sum(randp.sample<=s)/n.rep)
  r = mcf(p.org[1:k], p.next[1:k], s)
  rej = sort(r, decreasing = T, index.return = T)$ix[1:R]
  return(rej) 
}
##############################################################################


##############################################################################
AS.MCF <- function(p.org, p.next, alpha, s, lambda){
  # fdr estimator
  n = length(p.org)
  n.rep = 10
  randp = matrix(rep(0,n.rep*n), nrow = n)
  for(i in 1:n.rep){
    randp[,i] = runif(n, p.next, p.org)
  }
  
  fdr.est <- function(p.org, p.next, s, lambda, k){
    randp.sample = as.vector(randp[1:k,])
    A = sum(randp.sample > lambda)/n.rep
    R = max(sum(randp.sample <= s)/n.rep, 1)
    fdr.est = s/(1-lambda)*(1+A)/R
    return(fdr.est)
  }
  # determine k
  k = n
  fdr.tmp = fdr.est(p.org, p.next, s, lambda, k)
  while(fdr.tmp > alpha){
    k = k-1 
    fdr.tmp = fdr.est(p.org, p.next, s, lambda, k)
    # print(c(k,fdr.tmp))
    if(k<0){
      stop('No suitable k found since s is too large!')
    }
  }
  # determine rejetions
  randp.sample = as.vector(randp[1:k,])
  R = round(sum(randp.sample<=s)/n.rep)
  r = mcf(p.org[1:k], p.next[1:k], s)
  rej = sort(r, decreasing = T, index.return = T)$ix[1:R]
  return(rej) 
}
##############################################################################

















