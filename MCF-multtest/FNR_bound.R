FNR_bound_bt <- function(lambda1) {
  N <- rpois(n.test, lambda1) + 1 # total of the binomial
  P.all <- c(rep(0.5,4500), 0.5 + sample(c(-1,1),500,replace = TRUE)*runif(500,min=0.2,max=0.5))

  p.org <- rep(0,n.test) 
  p.next <- rep(0,n.test)
  p.min <- rep(0,n.test)
  X <- rep(0,n.test)
  
  ### get pvalues ####################
  for(j in 1:n.test){
    n <- N[j]
    a <- rbinom(1,n,P.all[j])
    X[j] = a
    p.org[j] <- binom.test(a,n,p=0.5,alternative="two.sided",conf.level = 0.95)$p.value
    if(a>0&&a<n){
      p.next.temp1 = binom.test(a-1,n,p=0.5,alternative="two.sided",conf.level = 0.95)$p.value
      p.next.temp2 = binom.test(a+1,n,p=0.5,alternative="two.sided",conf.level = 0.95)$p.value
      if(p.next.temp1<p.org[j] & p.next.temp2<p.org[j]){
        p.next[j] = max(p.next.temp1, p.next.temp2)
      } else {
        p.next[j] = min(p.next.temp1, p.next.temp2)
      }
    }else{
      p.next[j] <- 0
    }
    p.min[j] <- binom.test(0,n,p=0.5,alternative="two.sided",conf.level = 0.95)$p.value
  }
  p.org[p.org>1] <- 1 #remove rounding mistake
  p.next[p.next>1] <- 1
  p.min[p.min>1] <- 1
  
  n.rep <- 1000
  n.ecdf <- n.rep*n.test
  randp.ecdf <- rep(0,n.ecdf)
  sum.pi0 = 0
  for(i in 0:(n.rep-1)){
    randp = runif(n.test,p.next,p.org)
    randp.ecdf[(i*n.test+1):((i+1)*n.test)] <- randp
    # print(i)
  }
  pi0 = 0.9
  pi1 = 0.1
  
  lambda.all <- c(seq(0.00001,0.01,by=0.00003), seq(0.01,0.03,by=0.0001))
  pFDR.all <- rep(0,length(lambda.all))
  pFNR_bound.all <- rep(0,length(lambda.all))
  
  for(i in 1:length(lambda.all)){
    l = lambda.all[i]
    cdf = sum(randp.ecdf<l)/length(randp.ecdf)
    pFDR.all[i] = pi0*l/cdf
    pFNR_bound.all[i] = pi1*(1-(cdf-pi0*l)/pi1)/(1-cdf)
    # print(i)
  }
  
  alpha.all = seq(0.001,0.2,by=0.002)
  n.alpha = length(alpha.all)
  pFNR_bound = rep(0,n.alpha)
  
  for(i.alpha in 1:n.alpha){
    ### get lambda.star #############
    alpha = alpha.all[i.alpha]
    
    if(pFDR.all[1]<alpha){
      idx.lambda.star = max(which(pFDR.all<alpha))
    }else{
      idx.lambda.star = 1
    }
    pFNR_bound[i.alpha] = pFNR_bound.all[idx.lambda.star]
    # print(c(alpha,idx.lambda.star,pFDR.all[idx.lambda.star],pFNR_bound[i.alpha]))
  }
  return(list(alpha=alpha.all, pFNR_bound=pFNR_bound))
}



FNR_bound_fet <- function(lambda1, lambda2) {
  n.test = 5000
  N1 <- rpois(n.test, lambda1) # total of each binomial
  N2 <- rpois(n.test, lambda2)
  p1.Null <- runif(n.test*0.9, min=0.1,max=0.9)
  p2.Null <- p1.Null
  p1.Nonnull <- 0.5 * runif(n.test*0.1,min=0.1,max=0.9) 
  p2.Nonnull <- p1.Nonnull + runif(n.test*0.1,min=0.2,max=0.5)
  p1 <- c(p1.Null, p1.Nonnull)
  p2 <- c(p2.Null, p2.Nonnull)

  
  ### generate each binomial ##############################
  p.org <- rep(0,n.test) 
  p.next <- rep(0,n.test)
  p.min <- rep(0,n.test)
  Z1 <- rep(0,n.test)
  Z2 <- rep(0,n.test)
  Z3 <- rep(0,n.test)
  Z4 <- rep(0,n.test)
  
  ### get pvalues ####################
  for(j in 1:n.test){
    n1 = N1[j]
    n2 = N2[j]
    a <- rbinom(1,n1,p1[j])
    b <- rbinom(1,n2,p2[j])
    Z1[j] = a
    Z2[j] = n1-a
    Z3[j] = b
    Z4[j] = n2-b
    l <- max(0,a+b-n2)
    u <- min(a+b,n1)
    prob.all <- dhyper(l:u,n1,n2,(a+b))
    prob.obs <- dhyper(a,n1,n2,(a+b))
    p.org[j] <- sum(prob.all[which(prob.all<=prob.obs)])
    p.next[j] <- sum(prob.all[which(prob.all<prob.obs)])
    p.min[j] <- min(prob.all)
  }
  p.org[p.org>1] <- 1 #remove rounding mistake
  p.next[p.next>1] <- 1
  p.min[p.min>1] <- 1
  
  n.rep <- 1000
  n.ecdf <- n.rep*n.test
  randp.ecdf <- rep(0,n.ecdf)
  sum.pi0 = 0
  for(i in 0:(n.rep-1)){
    randp = runif(n.test,p.next,p.org)
    randp.ecdf[(i*n.test+1):((i+1)*n.test)] <- randp
    # print(i)
  }
  pi0 = 0.9
  pi1 = 0.1
  
  lambda.all <- c(seq(0.00001,0.01,by=0.00003), seq(0.01,0.03,by=0.0001))
  pFDR.all <- rep(0,length(lambda.all))
  pFNR_bound.all <- rep(0,length(lambda.all))
  
  for(i in 1:length(lambda.all)){
    l = lambda.all[i]
    cdf = sum(randp.ecdf<l)/length(randp.ecdf)
    pFDR.all[i] = pi0*l/cdf
    pFNR_bound.all[i] = pi1*(1-(cdf-pi0*l)/pi1)/(1-cdf)
    # print(i)
  }
  
  alpha.all = seq(0.001,0.2,by=0.002)
  n.alpha = length(alpha.all)
  pFNR_bound = rep(0,n.alpha)
  
  for(i.alpha in 1:n.alpha){
    ### get lambda.star #############
    alpha = alpha.all[i.alpha]
    
    if(pFDR.all[1]<alpha){
      idx.lambda.star = max(which(pFDR.all<alpha))
    }else{
      idx.lambda.star = 1
    }
    pFNR_bound[i.alpha] = pFNR_bound.all[idx.lambda.star]
    # print(c(alpha,idx.lambda.star,pFDR.all[idx.lambda.star],pFNR_bound[i.alpha]))
  }
  return(list(alpha=alpha.all, pFNR_bound=pFNR_bound))
}

