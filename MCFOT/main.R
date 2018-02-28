source("~/Documents/Paper2/Rcode/methods.R")
source("~/Documents/Paper2/Rcode/data.R")
source("~/Documents/Paper2/Rcode/performance.R")
source("~/Documents/Paper2/Rcode/AT_randp.R")

# non-null proportion
Pi1 <- function(t){
  gamma = 0.2 #
  b = 10
  return(gamma*exp(-b*t)*b/(1-exp(-b)))
}

# generate data
n.test = 5000
mu = 25
data.temp = simu.FET(n.test, Pi1, mu)
p.org = data.temp$p.org
p.next = data.temp$p.next
H = data.temp$H
# plot(H[1:2000],cex=.1)

perf = matrix(nrow = 10, ncol = 4)
perf[,1] = c('BH','SBH','SS','AS','BH.MCF','SBH.MCF','SS.MCF','AS.MCF','AT','AT.MCF')
colnames(perf) = c('method','R','fdr','power')

q = 0.1
s = 0.1*q
lambda = 0.5
for(i in 1:nrow(perf)){
  mtd = perf[i,1]
  if (mtd == 'BH'){
    rej = BH(p.org, q)
  } else if (mtd == 'SBH') {
    rej = SBH(p.org, q, lambda)
  } else if (mtd == 'SS') {
    rej = SS(p.org, q, s) 
  } else if (mtd == 'AS') {
    rej = AS(p.org, q, s, lambda)
  } else if (mtd == 'BH.MCF') {
    rej = BH.MCF(p.org, p.next, q)
  } else if (mtd == 'SBH.MCF') {
    rej = SBH.MCF(p.org, p.next, q, lambda)
  } else if (mtd == 'SS.MCF') {
    rej = SS.MCF(p.org, p.next, q, s)
  } else if (mtd == 'AS.MCF') {
    rej = AS.MCF(p.org, p.next, q, s, lambda)
  } else if (mtd == 'AT') {
    k.hat = HingeExp(p.org*0.99,alpha=q,C=2)
    rej = 1:(k.hat)
  } else if (mtd == 'AT.MCF') {
    k.hat = HingeExp_randp(p.org*0.99,p.next*0.99,alpha=q,C=2)
    rej = 1:(k.hat)
  }
  # number of rjections, fdr, statistical power
  tmp = performance(rej, H)
  perf[i,2] = tmp$R
  perf[i,3] = tmp$fdr
  perf[i,4] = tmp$power
}

write.csv(perf, '~/Documents/Paper2/Simu/perf1.csv', quote = F, row.names = F)






for(j in 1:3){
  n.test = 5000
  mu = 25
  data.temp = simu.FET(n.test, Pi1, mu)
  p.org = data.temp$p.org
  p.next = data.temp$p.next
  H = data.temp$H
  
  perf = matrix(nrow = 8, ncol = 4)
  perf[,1] = c('BH','SBH','SS','AS','BH.MCF','SBH.MCF','SS.MCF','AS.MCF')
  colnames(perf) = c('method','R','fdr','power')
  
  q = 0.1
  s = 0.1*q
  lambda = 0.5
  for(i in 1:8){
    mtd = perf[i,1]
    if (mtd == 'BH'){
      rej = BH(p.org, q)
    } else if (mtd == 'SBH') {
      rej = SBH(p.org, q, lambda)
    } else if (mtd == 'SS') {
      rej = SS(p.org, q, s) 
    } else if (mtd == 'AS') {
      rej = AS(p.org, q, s, lambda)
    } else if (mtd == 'BH.MCF') {
      rej = BH.MCF(p.org, p.next, q)
    } else if (mtd == 'SBH.MCF') {
      rej = SBH.MCF(p.org, p.next, q, lambda)
    } else if (mtd == 'SS.MCF') {
      rej = SS.MCF(p.org, p.next, q, s)
    } else if (mtd == 'AS.MCF') {
      rej = AS.MCF(p.org, p.next, q, s, lambda)
    }
    # number of rjections, fdr, statistical power
    tmp = performance(rej, H)
    perf[i,2] = tmp$R
    perf[i,3] = tmp$fdr
    perf[i,4] = tmp$power
  }
  
  write.csv(perf, paste0('~/Documents/Paper2/Simu/perf',j,'.csv'), quote = F, row.names = F)
}


perf = read.csv('~/Documents/Paper2/Simu/perf.csv')
perf1 = read.csv('~/Documents/Paper2/Simu/perf1.csv')
perf2 = read.csv('~/Documents/Paper2/Simu/perf2.csv')
perf3 = read.csv('~/Documents/Paper2/Simu/perf3.csv')

tmp1 = (perf[,2:4] + perf1[,2:4] + perf2[,2:4] + perf3[,2:4])/4
result = cbind(perf[,1],tmp1)



