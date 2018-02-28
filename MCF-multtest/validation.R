
### data set-up ###############
lambda1 <- 15 # or 15

### fixed across repeated simulations ########################
n.test = 5000
# N <- rpois(n.test, lambda1) + 1 # total of the binomial
# P.all <- c(rep(0.5,4500), 0.5 + runif(500,min=0.2,max=0.5))

alpha = 0.1

p.org <- rep(0,n.test) 
p.next <- rep(0,n.test)
p.min <- rep(0,n.test)
X <- rep(0,n.test)

### get pvalues ####################
for(j in 1:n.test){
  n <- N[j]
  a <- rbinom(1,n,P.all[j])
  X[j] = a
  p.org[j] <- binom.test(a,n,p=0.5,alternative="greater",conf.level = 0.95)$p.value
  if(a>0&&a<n){
    p.next.temp1 = binom.test(a-1,n,p=0.5,alternative="greater",conf.level = 0.95)$p.value
    p.next.temp2 = binom.test(a+1,n,p=0.5,alternative="greater",conf.level = 0.95)$p.value
    if(p.next.temp1<p.org[j] & p.next.temp2<p.org[j]){
      p.next[j] = max(p.next.temp1, p.next.temp2)
    } else {
      p.next[j] = min(p.next.temp1, p.next.temp2)
    }
  }else{
    p.next[j] <- 0
  }
  p.min[j] <- binom.test(0,n,p=0.5,alternative="greater",conf.level = 0.95)$p.value
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
  storey.lambda = 0.5
  sum.pi0 = sum.pi0 + (sum(randp > storey.lambda)+1) / ((1-storey.lambda)*5000)
  randp.ecdf[(i*n.test+1):((i+1)*n.test)] <- randp
  # print(i)
}

### calculate pi0 ##################
pi0 = sum.pi0/n.rep
print(paste('pi0 = ', pi0))

lambda.all <- c(seq(0.00001,0.01,by=0.00003), seq(0.01,0.03,by=0.0001))
pFDR.all <- rep(0,length(lambda.all))
for(i in 1:length(lambda.all)){
  l = lambda.all[i]
  cdf = sum(randp.ecdf<l)/length(randp.ecdf)
  pFDR.all[i] = pi0*l/cdf
  # print(i)
}

if(pFDR.all[1]<alpha){
  idx.lambda.star = max(which(pFDR.all<alpha))
}else{
  idx.lambda.star = 1
}
lambda.star = lambda.all[idx.lambda.star]


l = lambda.star
cdf = sum(randp.ecdf<l)/length(randp.ecdf)
R <- round(n.test*cdf)
r <- rep(0,n.test)
for(s in 1:n.test){
  if(l>p.org[s]){
    r[s] <- 1}else if(l<p.next[s]){
      r[s] <- 0}else{
        r[s] <- (l-p.next[s])/(p.org[s]-p.next[s])
      }
}
rej.mcf <- sort(r, decreasing = F, index.return = T)$ix[(n.test-R+1):n.test]
FDP = sum(rej.mcf<4501)/length(rej.mcf)
FDP


#######################################################

lambda.star
lambda = lambda.star
pi0
pi1 = 1-pi0
pi1
u = sum(randp.ecdf<lambda)/length(randp.ecdf)
F1_lamba = (u-pi0*lambda)/pi1


W_i = rep(0,5000)
Pr_0 = rep(0,5000)
Pr_1 = rep(0,5000)
Pr_w = rep(0,5000)
A_lambda = rep(0,5000)
B_lambda = rep(0,5000)


for(i in 1:5000){
  n = N[i]
  p = P.all[i]
   
  # pr(x), x=0,...,n
  pr_x_all = dbinom(n:0, n, 0.5)
  # possible p-values
  p_all = cumsum(pr_x_all)
  # prob of possible p-values
  pr_p = dbinom(n:0, n, p)
  # how many possible p-values under lambda
  id_a = sum(p_all<lambda)
  # 
  if(id_a == 0){
    a_lambda = 0
    Pr_1[i] = 0
  } else {
    a_lambda = p_all[id_a]
    Pr_1[i] = sum(pr_p[1:id_a])
  }
  b_lambda = p_all[id_a+1]
  #
  W_i[i] = (lambda - a_lambda)/(b_lambda - a_lambda)
  Pr_w[i] = pr_p[id_a+1]
  #
  if(id_a == (n+1)){
    Pr_0[i] = 0
  } else {
    Pr_0[i] = sum(pr_p[(id_a+2):(n+1)])
  }
  A_lambda[i] = a_lambda
  B_lambda[i] = b_lambda
}

# lambda
# mean((W_i*Pr_w + Pr_1)[1:4500])
# 
# F1_lambda
# mean((W_i*Pr_w + Pr_1)[4501:5000])


####################################

q = 0.2
id_greater_q = which(W_i>=q)
id_smaller_q = which(W_i<q)
temp.u = (sum(Pr_1[id_smaller_q]) + sum((Pr_1+Pr_w)[id_greater_q]))/5000


q = 0.0001
id_greater_q = which(W_i>=q)
id_smaller_q = which(W_i<q)
temp.u = (sum(Pr_1[id_smaller_q]) + sum((Pr_1+Pr_w)[id_greater_q]))/5000
temp.u

while(temp.u > u){
  q = q+0.0001
  id_greater_q = which(W_i>=q)
  id_smaller_q = which(W_i<q)
  temp.u = (sum(Pr_1[id_smaller_q]) + sum((Pr_1+Pr_w)[id_greater_q]))/5000
  print(c(q,temp.u))
}

id_greater_q_h0 = id_greater_q[which(id_greater_q<4501)]
id_greater_q_h1 = id_greater_q[which(id_greater_q>4500)]



rho = (sum(Pr_1[1:4500]) + sum(Pr_w[id_greater_q_h0]) )/4500 # 0.00675, 0.00675
lambda                                                       # 0.00802, 0.00787

pi0*lambda/u # nominal FDR level = 0.1: 0.0999553, 0.09974428 
pi0*rho/u    # estimate FDR:            0.0841412, 0.08560531
FDP          # actural  FDP:            0.0888889, 0.08732394

(sum((Pr_w*W_i)[4501:500])-sum(Pr_w[id_greater_q_h1]))/4500 # A-B = -0.001738494, -0.00166716
rho - lambda                                                # C-D = -0.001268855, -0.00111559













