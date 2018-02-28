library(Rcpp)
sourceCpp("/Users/xiaoyudai/Documents/Paper/Rcode/mcf.cpp")

MCF = function(p.org, p.next, alpha , n.rep=1000){
  # MCF-based multiple testing procedure.
  # References: 
  #  Dai, X., Lin, N., Li, D., Wang, T.(2016) A non-randomized procedure for discrete multiple testing based on randomized test. (submitted)
  # Input:
  #   p.org:  raw p-value for each test
  #   p.next: the next possible smaller p-value, i.e. p^-
  #   alpha:  nominal FDR level
  #   n.rep:  number of samples of randomized p-values, to estimate null proportion pi_0 and alternative distribution F_1 
  # Returns:
  #   pi0: estimated null proportion  
  #   Rej: rejection result, a vector of length m (number of tests): 
  #         1: rejected
  #         0: accepted
  #
  
  # number of tests
  n.test = length(p.org)

  # generate and pool randomized p-values
  n.ecdf <- n.rep*n.test
  randp.ecdf <- rep(0,n.ecdf)
  sum.pi0 = 0
  for(i in 0:(n.rep-1)){
    randp = runif(n.test,p.next,p.org)
    randp.ecdf[(i*n.test+1):((i+1)*n.test)] <- randp
    # Storey's estimator with lambda = 0.5
    storey.lambda = 0.5
    sum.pi0 = sum.pi0 + (sum(randp > storey.lambda)+1) / ((1-storey.lambda)*5000)
  }
  
  ### average pi0 ##################
  pi0 = sum.pi0/n.rep
  
  
  # find lambda.star
  a = 0.0000001
  b = 0.9999999
  lambda.star = 0
  num.rej = 0
  pFDR = 0
  while( abs(pFDR-alpha)>0.0000001 ){
    lambda.star = (a+b)/2
    cdf <- sum(randp.ecdf<lambda.star)/length(randp.ecdf)
    pFDR = pi0*lambda.star/cdf  
    num.rej = round(n.test*cdf) 
    print(c(a,b,lambda.star, pFDR, num.rej))
    if(pFDR > alpha){
      b=lambda.star
    } else {
      a=lambda.star
    }  
  }
  
  # find MCF values corresponding to lambda.star
  mcf.value <- mcf(p.org, p.next, lambda.star)
  REJ <- rep(0,n.test)
  rej.idx = sort(mcf.value, decreasing = T, index.return = T)$ix[1:num.rej]
  REJ[rej.idx]=1
  
  return(list(REJ=REJ,pi0=pi0))
}



# example
# two-sided Fisher's Exact Test
p.org = rep(0,5000)
p.next = rep(0,5000)
for(j in 1:5000){
  n1 = rpois(1, lambda = 20) + 1
  n2 = rpois(1, lambda = 25) + 1
  a <- rbinom(1, n1, 0.5)
  # 1-4500 true null, 4501-5000 true non-null
  if(j < 4501){
    b <- rbinom(1, n2, 0.5)
  } else {
    b <- rbinom(1, n2, 0.8)
  }
  l <- max(0,a+b-n2)
  u <- min(a+b,n1)
  prob.all <- dhyper(l:u,n1,n2,(a+b))
  prob.obs <- dhyper(a,n1,n2,(a+b))
  p.org[j] <- sum(prob.all[which(prob.all<=prob.obs)])
  p.next[j] <- sum(prob.all[which(prob.all<prob.obs)])
}
p.org[p.org>1] <- 1 #remove rounding mistake
p.next[p.next>1] <- 1

result = MCF(p.org,p.next,alpha=0.1,n.rep = 1000)
names(result)
result$pi0

