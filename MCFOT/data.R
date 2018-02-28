simu.FET <- function(n.test, Pi1, mu){
  # simulate ordered multiple testing data from FET
  # null:     p1 = p2 ~ U[0,1]
  # non-null: p1 ~ U[0,0.5], p2 = p1 + U[0.2,0.5]
  # Input:
  #   Pi1: function of pi1
  #   mu: mean of Poisson dist to get n_i
  # Output:
  #   p.org: raw p-value
  #   p.next: next smaller possible p-value
  #   H: indicator of whether null (0) or non-null (1)
  
  p.org = rep(0,n.test)
  p.next = rep(0,n.test)
  H = rep(0,n.test)
  for(i in 1:n.test){
    pi1 = Pi1(i/n.test)
    pi0 = max(0,1-pi1)
    H[i] = sample(x = 0:1, size = 1, prob = c(pi0, pi1))
    if(H[i] == 0){
      p1 = runif(1,min=0,max=1)
      p2 = p1
    } else {
      p1 = runif(1,min=0,max=0.5)
      p2 = p1 + runif(1,min=0.2,max=0.5)
    }
    
    n1 = rpois(1, mu)
    n2 = rpois(1, mu)
    a <- rbinom(1,n1,p1)
    b <- rbinom(1,n2,p2)
    l <- max(0,a+b-n2)
    u <- min(a+b,n1)
    prob.all <- dhyper(l:u,n1,n2,(a+b))
    prob.obs <- dhyper(a,n1,n2,(a+b))
    p.org[i] <- sum(prob.all[which(prob.all<=prob.obs)])
    p.next[i] <- sum(prob.all[which(prob.all<prob.obs)])
  }
  return(list(p.org=p.org, p.next=p.next, H=H))
}




