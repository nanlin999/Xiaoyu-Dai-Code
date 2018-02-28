source("~/Documents/Paper2/Rcode/methods.R")
source("~/Documents/Paper2/Rcode/data.R")
source("~/Documents/Paper2/Rcode/performance.R")
source("~/Documents/Paper2/Rcode/AT_randp.R")
source("~/Documents/Paper2/Rcode/multiplot.R")
library(ggplot2)

methods = c('SS','SS-MCF','AS','AS-MCF','AT','AT-MCF')
n.method = length(methods)
n.simu = 100
FDP = matrix(rep(0,n.method*n.simu), nrow = n.method)
POWER = matrix(rep(0,n.method*n.simu), nrow = n.method)
row.names(FDP) = methods
row.names(POWER) = methods

# non-null proportion
Pi1 <- function(t){
  gamma = 0.05 #
  b = 50
  return(gamma*exp(-b*t)*b/(1-exp(-b)))
}

# parameter configuration

# nominal FDR level
alpha = 0.1
# s for SS, AS, SS.MCF and AS.MCF
s = 0.5 * alpha
# lambda for AS and AS.MCF
lambda = 0.5


for(i.simu in 1:n.simu){
  # generate data
  n.test = 5000
  mu = 25
  data.temp = simu.FET(n.test, Pi1, mu)
  p.org = data.temp$p.org
  p.next = data.temp$p.next
  H = data.temp$H
  
  for(mtd in methods){
    if (mtd == 'SS'){
      rej = SS(p.org, alpha, s)
    } else if (mtd == 'SS-MCF') {
      rej = SS.MCF(p.org, p.next, alpha, s)
    } else if (mtd == 'AS') {
      rej = AS(p.org, alpha, s, lambda) 
    } else if (mtd == 'AS-MCF') {
      rej = AS.MCF(p.org, p.next, alpha, s, lambda)
    } else if (mtd == 'AT') {
      k_hat = HingeExp(p.org*0.99,alpha,C=2)
      rej = 1:k_hat
    } else if (mtd == 'AT-MCF') {
      k_hat = HingeExp_randp(p.org*0.99,p.next*0.99,alpha,C=2)
      rej = 1:k_hat
    } else {
      stop("unrecognized method!")
    }
    perf = performance(rej, H)
    # perf[i,2] = tmp$R
    FDP[mtd,i.simu] = perf$fdr
    POWER[mtd,i.simu] = perf$power
  }
}

write.table(FDP, file = '~/Documents/Paper2/Simu/FDP/fdp_gamma005b50.txt', col.names = F, quote = F)
write.table(POWER, file = '~/Documents/Paper2/Simu/POWER/power_gamma005b50.txt', col.names = F, quote = F)


# box plot

FDP = read.table('~/Documents/Paper2/Simu/FDP/fdp_gamma005b50.txt', row.names = 1)
POWER =  read.table('~/Documents/Paper2/Simu/POWER/power_gamma005b50.txt', row.names = 1)

FDP.dat <- data.frame(methods = factor(rep(methods, each=n.simu)), FDP = as.vector(t(FDP)))
FDP.dat$methods=factor(FDP.dat$methods , levels=levels(FDP.dat$methods)[c(5,1,3,6,2,4)])

FDP.plt <- ggplot(FDP.dat, aes(x=methods, y=FDP)) + 
  geom_boxplot() +
  coord_cartesian(ylim = c(0,0.15)) +
  geom_hline(yintercept=0.1) + 
  xlab(expression(paste(gamma,' = 0.05, b = 50'))) 
FDP.plt

POWER.dat <- data.frame(methods = factor(rep(methods, each=n.simu)), POWER = as.vector(t(POWER)))
POWER.dat$methods=factor(POWER.dat$methods , levels=levels(POWER.dat$methods)[c(5,1,3,6,2,4)])

POWER.plt <- ggplot(POWER.dat, aes(x=methods, y=POWER)) + 
  geom_boxplot() + 
  coord_cartesian(ylim = c(0,1)) + 
  xlab(expression(paste(gamma,' = 0.05, b = 50')))
  # geom_hline(yintercept=0.1)
POWER.plt

a4 = FDP.plt
b4 = POWER.plt


# combine plots
pdf("~/Documents/Paper2/Simu/FDP/fdp.pdf",width=12,height=8)
multiplot(a1,a3,a2,a4, cols = 2 )
dev.off()

pdf("~/Documents/Paper2/Simu/POWER/power.pdf",width=12,height=8)
multiplot(b1,b3,b2,b4, cols = 2 )
dev.off()

