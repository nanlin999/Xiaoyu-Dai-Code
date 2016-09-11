library(Rcpp)
sourceCpp("/Users/xiaoyudai/Documents/multiple-testing/Rcode/Rcpp/MCF.cpp")

library(cp4p)


MT = function(p.org, p.next, alpha , n.rep=1000, pi0.method="pounds"){
  n.test = length(p.org)

  # generate and pool randomized p-values
  n.ecdf <- n.rep*n.test
  randp.ecdf <- rep(0,n.ecdf)
  for(i in 0:(n.rep-1)){
    randp.ecdf[(i*n.test+1):((i+1)*n.test)] <- runif(n.test,p.next,p.org)
  }
  
  # estimate pi0
  methods = c("abh", "st.spline", "st.boot", "langaas", "histo", "pounds", "jiang", "slim")
  if(pi0.method %in% methods){
    pi0 = estim.pi0(randp.ecdf,pi0.method = pi0.method)$pi0
  }else{
    warning("Warning: pi0 method doesn't find")
    pi0 = 1
  }
  
  
  # find lambda.star
  a = 0.0000001
  b = 0.1
  lambda.star = 0
  num.rej = 0
  pFDR = 0
  while( abs(pFDR-alpha)>0.0000001 ){
    lambda.star = (a+b)/2
    cdf <- sum(randp.ecdf<lambda.star)/length(randp.ecdf)
    pFDR = pi0*lambda.star/cdf  
    num.rej = round(n.test*cdf) 
    if(pFDR > alpha){
      b=lambda.star
    }
    if(pFDR < alpha){
      a=lambda.star
    }  
  }
  
  # find MCF values corresponding to lambda.star
  mcf <- MCF(p.org, p.next, lambda.star)
  REJ <- rep(0,n.test)
  REJ[which(mcf>sort(mcf)[n.test-num.rej])]=1
  
  return(list(REJ=REJ,pi0=pi0))
}

result = MT(p.org,p.next,alpha=0.1,pi0.method = "histo")
result$pi0
