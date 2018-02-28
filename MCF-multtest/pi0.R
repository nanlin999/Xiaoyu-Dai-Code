
### define function ##############
library(cp4p)
PI = function(p.org, p.next, n.rep=1000, pi0.method="pounds"){
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
  return(pi0)
}








### data set-up ###############
lambda1 <- 20 # or 20
lambda2 <- 20 # or 20

n.test = 5000
p.org <- rep(0,n.test) 
p.next <- rep(0,n.test)
p.min <- rep(0,n.test)
Z1 <- rep(0,n.test)
Z2 <- rep(0,n.test)
Z3 <- rep(0,n.test)
Z4 <- rep(0,n.test)

N <- rep(0,n.test)
X <- rep(0,n.test)


for(i in 1:4500){
  ### Fisher's Exact Test ########################
  # n1 <- rpois(1,lambda1)
  # n2 <- rpois(1,lambda2)
  # p.same <- runif(1,min=0,max=1)
  # a <- rbinom(1,n1,p.same)
  # b <- rbinom(1,n2,p.same)
  # Z1[i] = a
  # Z2[i] = n1-a
  # Z3[i] = b
  # Z4[i] = n2-b
  # l <- max(0,a+b-n2)
  # u <- min(a+b,n1)
  # prob.all <- dhyper(l:u,n1,n2,(a+b))
  # prob.obs <- dhyper(a,n1,n2,(a+b))
  # p.org[i] <- sum(prob.all[which(prob.all<=prob.obs)])
  # p.next[i] <- sum(prob.all[which(prob.all<prob.obs)])
  # p.min[i] <- min(prob.all)
  
  ### binomial test ####################
  n1 <- rpois(1,lambda1)
  a <- rbinom(1,n1,0.5)
  X[i] = a
  N[i] = n1
  p.org[i] <- binom.test(a,n1,p=0.5,alternative="two.sided",conf.level = 0.95)$p.value
  if(a>0&&a<n1){
    p.next[i] <- min(binom.test(a-1,n1,p=0.5,alternative="two.sided",conf.level = 0.95)$p.value,
                     binom.test(a+1,n1,p=0.5,alternative="two.sided",conf.level = 0.95)$p.value)
  }else{
    p.next[i] <- 0
  }
  p.min[i] <- binom.test(0,n1,p=0.5,alternative="two.sided",conf.level = 0.95)$p.value
}
for(i in 4501:5000){
  ### Fisher's Exact Test ########################
  # n1 <- rpois(1,lambda1)
  # n2 <- rpois(1,lambda2)
  # potential.p1 <- runif(1,min=0.5,max=1)
  # potential.p2 <- runif(1,min=0,max=0.5)
  # p1 <- max(potential.p1,potential.p2)
  # p2 <- min(potential.p1,potential.p2)
  # p.same <- runif(1,min=0,max=1)
  # p1 <- p.same*0.5+0.5
  # p2 <- p.same*0.5
  # a <- rbinom(1,n1,p1)
  # b <- rbinom(1,n2,p2)
  # Z1[i] = a
  # Z2[i] = n1-a
  # Z3[i] = b
  # Z4[i] = n2-b
  # l <- max(0,a+b-n2)
  # u <- min(a+b,n1)
  # prob.all <- dhyper(l:u,n1,n2,(a+b))
  # prob.obs <- dhyper(a,n1,n2,(a+b))
  # p.org[i] <- sum(prob.all[which(prob.all<=prob.obs)])
  # p.next[i] <- sum(prob.all[which(prob.all<prob.obs)])
  # p.min[i] <- min(prob.all)
  
  ### binomial test ####################
  n1 <- rpois(1,lambda1)
  p.alter <- 0.5 + sample(c(-1,1),1)*runif(1,min=0.3,max=0.5)
  a <- rbinom(1,n1,p.alter)
  X[i] = a
  N[i] = n1
  p.org[i] <- binom.test(a,n1,p=0.5,alternative="two.sided",conf.level = 0.95)$p.value
  if(a>0&&a<n1){
    p.next[i] <- min(binom.test(a-1,n1,p=0.5,alternative="two.sided",conf.level = 0.95)$p.value,
                     binom.test(a+1,n1,p=0.5,alternative="two.sided",conf.level = 0.95)$p.value)
  }else{
    p.next[i] <- 0
  }
  p.min[i] <- binom.test(0,n1,p=0.5,alternative="two.sided",conf.level = 0.95)$p.value
}
p.org[p.org>1] <- 1 #remove rounding mistake
p.next[p.next>1] <- 1
p.min[p.min>1] <- 1

tmp.Chen.FET = data.frame(cbind(Z1,Z3,Z1+Z2,Z3+Z4))
tmp.Chen.BT = data.frame(cbind(X,N-X))


ptm <- proc.time()
PI(p.org,p.next,n.rep=1000,pi0.method="pounds")
proc.time() - ptm

ptm <- proc.time()
PI(p.org,p.next,n.rep=1000,pi0.method="histo")
proc.time() - ptm

ptm <- proc.time()
PI(p.org,p.next,n.rep=1000,pi0.method="st.spline")
proc.time() - ptm

alpha=0.05

ptm <- proc.time()
k <- 1
m <- 1
while(k<=n.test){
  m <- sum(p.min < alpha/k)
  if(m<=k) break
  k = k+1
}
R <- which(p.min < alpha/k)
temp.p <- p.org[R]
estim.pi0(temp.p,pi0.method = "pounds")$pi0
proc.time() - ptm

ptm <- proc.time()
k <- 1
m <- 1
while(k<=n.test){
  m <- sum(p.min < alpha/k)
  if(m<=k) break
  k = k+1
}
R <- which(p.min < alpha/k)
temp.p <- p.org[R]
estim.pi0(temp.p,pi0.method = "st.spline")$pi0
proc.time() - ptm

ptm <- proc.time()
k <- 1
m <- 1
while(k<=n.test){
  m <- sum(p.min < alpha/k)
  if(m<=k) break
  k = k+1
}
R <- which(p.min < alpha/k)
temp.p <- p.org[R]
estim.pi0(temp.p,pi0.method = "histo")$pi0
proc.time() - ptm

ptm <- proc.time()
mt.Chen = GeneralizedFDREstimators(data=tmp.Chen.BT,
                                   Test= "Binomial Test",FET_via = "IndividualMarginals",
                                   FDRlevel=alpha, lambda=0.5, epsilon=1)
mt.Chen$Generalized_Estimator$pi0Est
proc.time() - ptm