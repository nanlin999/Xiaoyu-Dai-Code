library(qvalue)  
library(Rcpp)
sourceCpp("/Users/xiaoyudai/Documents/multiple-testing/Rcode/Rcpp/MCF.cpp")

### CD4 VS CD8 ###

ORG <- read.table("/Users/xiaoyudai/Documents/multiple-testing/Comparison/CD8_ES/CD8_ES.txt", header=F)


Z1 = ORG[,4] +ORG[,6]
Z2 = ORG[,5] +ORG[,7] - (ORG[,4] +ORG[,6])
Z3 = ORG[,8] +ORG[,10]
Z4 = ORG[,9] +ORG[,11] - (ORG[,8] +ORG[,10])


# remove small counts
counts_cutoff = 15
idx = which((Z1+Z2)<counts_cutoff|(Z3+Z4)<counts_cutoff)
idx = which((Z1+Z2)<counts_cutoff|(Z3+Z4)<counts_cutoff)
V = ORG[-idx,]
Z1 = Z1[-idx]
Z2 = Z2[-idx]
Z3 = Z3[-idx]
Z4 = Z4[-idx]

n.test = nrow(V)
p.org <- rep(0,n.test)
p.next <- rep(0,n.test)

for(i in 1:n.test){
  n1 <- Z1[i] + Z2[i]
  n2 <- Z3[i] + Z4[i]
  a <- Z1[i]
  b <- Z3[i]  
  l <- max(0,a+b-n2)
  u <- min(a+b,n1)
  prob.all <- dhyper(l:u,n1,n2,(a+b))
  prob.obs <- dhyper(a,n1,n2,(a+b))
  p.org[i] <- sum(prob.all[which(prob.all<=prob.obs)])
  p.next[i] <- sum(prob.all[which(prob.all<prob.obs)])
  print(i)
}
# rounding problems
p.org[which(p.org>1)] = 1
p.next[which(p.next>1)] = 1


MT <- qvalue(p.org) ## 5 mins
pi0 = MT$pi0 # 1
QVALUE <- MT$qvalues
# sum(QVALUE<0.1)


pi1 <- 1-pi0

n.rep <- 10
n.ecdf <- n.rep*n.test
randp.ecdf <- rep(0,n.ecdf)
for(i in 0:(n.rep-1)){
  randp.ecdf[(i*n.test+1):((i+1)*n.test)] <- runif(n.test,p.next,p.org)
  print(i)  
}

sourceCpp("/Users/xiaoyudai/Documents/multiple-testing/Rcode/Rcpp/MCF.cpp")
### under FDR nominal level 0.1 ########################
a = 0.00001
b = 0.1
l = 0
R = 0
pFDR = 0
leng = length(randp.ecdf)
while( abs(pFDR-0.1)>0.0000001 ){
  l = (a+b)/2
  cdf <- sum(randp.ecdf<l)/leng
  pFDR = (1-pi1)*l/cdf  
  pFDR
  if(pFDR > 0.1){
    b=l
  }
  if(pFDR < 0.1){
    a=l
  }  
  print(c(a,b,l,pFDR))
}
print(c(a,b,l,pFDR))
mcf <- MCF(p.org, p.next, l)
rej.0.1 <- which(mcf >= quantile(mcf, probs = 1-cdf, type = 1))
# rej.0.1 <- which(mcf>sort(mcf)[n.test-R])
rejC.0.1 <- which(QVALUE<0.1)

### under FDR nominal level 0.05 ###################
a = 0.00001
b = 0.1
l = 0
R = 0
pFDR = 0
leng = length(randp.ecdf)
while( abs(pFDR-0.05)>0.0000001 ){
  l = (a+b)/2
  cdf <- sum(randp.ecdf<l)/leng
  pFDR = (1-pi1)*l/cdf  
  R = round(n.test*cdf) 
  pFDR
  if(pFDR > 0.05){
    b=l
  }
  if(pFDR < 0.05){
    a=l
  }  
  print(c(a,b,l,pFDR))
}
print(c(a,b,l,pFDR))
mcf <- MCF(p.org, p.next, l)
rej.0.05 <- which(mcf >= quantile(mcf, probs = 1-cdf, type = 1))
rejC.0.05 <- which(QVALUE<0.05)


REJ010 <- rep(0,n.test)
REJC010 <- rep(0,n.test)
REJ010[rej.0.1] <- 1
REJC010[rejC.0.1] <- 1
REJ005 <- rep(0,n.test)
REJC005 <- rep(0,n.test)
REJ005[rej.0.05] <- 1
REJC005[rejC.0.05] <- 1

c(sum(REJ010), sum(REJC010), sum(REJ010*REJC010))
c(sum(REJ005), sum(REJC005), sum(REJ005*REJC005))





ALL = data.frame(cbind(V,p.org,p.next,REJ010,REJC010,REJ005,REJC005))
colnames(ALL) = c('chrome','start','end','A1','A2','A3','A4','B1','B2','B3','B4', 'porg','pnext', 'MCF010','Qvalue010', 'MCF005','Qvalue005')
write.table(ALL,"/Users/xiaoyudai/Documents/multiple-testing/Comparison/CD4_CD8/ALL.txt",row.name=F,quote=FALSE)


### MCF to find DMR
library(DSS)
ALL <- read.table("/Users/xiaoyudai/Documents/multiple-testing/Comparison/CD8_ES/ALL.txt", header=T)
### under FDR nominal level 0.1 ##############
MCF = ALL[,14]
Qvalue = ALL[,15]

data = data.frame(matrix(rep(0,11*nrow(ALL)),ncol=11))
colnames(data) = c('chr','pos','mu1','mu2','diff','diff.se','stat','phi1','phi2','pval','fdr')

data1 = data
data1[,1:2] = ALL[,1:2]
data1[,10] = 1-MCF
data1[,3] = (ALL[,4]+ALL[,6])/(ALL[,5]+ALL[,7])
data1[,4] = (ALL[,8]+ALL[,10])/(ALL[,9]+ALL[,11])
data1[,5] = data1[,3]-data1[,4]
# dmrs1 <- callDMR(data1, p.threshold=0.01,minlen=1000,dis.merge=100,pct.sig=0.5,minCG = 20)
dmrs1 <- callDMR(data1, p.threshold=0.01,minlen=500,dis.merge=100,pct.sig=0.5)
dim(dmrs1)

data2 = data
data2[,1:2] = ALL[,1:2]
data2[,10] = 1-Qvalue
data2[,3] = (ALL[,4]+ALL[,6])/(ALL[,5]+ALL[,7])
data2[,4] = (ALL[,8]+ALL[,10])/(ALL[,9]+ALL[,11])
data2[,5] = data2[,3]-data2[,4]
# dmrs2 <- callDMR(data2, p.threshold=0.01,minlen=1000,dis.merge=100,pct.sig=0.5,minCG = 20)
dmrs2 <- callDMR(data2, p.threshold=0.01,minlen=500,dis.merge=100,pct.sig=0.5)
dim(dmrs2)

write.table(dmrs1[,1:3],"/Users/xiaoyudai/Documents/multiple-testing/Comparison/bedtools/CD4_ES/MCF010.txt",
            row.names=F, quote=FALSE, col.names=F, sep = "\t")
write.table(dmrs2[,1:3],"/Users/xiaoyudai/Documents/multiple-testing/Comparison/bedtools/CD4_ES/Q010.txt",
            row.name=F, quote=FALSE, col.names=F, sep = "\t")

write.table(dmrs1[,1:8],"/Users/xiaoyudai/Documents/multiple-testing/Comparison/CD4_CD8/MCF010.txt",row.name=F,quote=FALSE)
write.table(dmrs2[,1:8],"/Users/xiaoyudai/Documents/multiple-testing/Comparison/CD4_CD8/Q010.txt",row.name=F,quote=FALSE)

### under FDR nominal level 0.05 ##############
 
write.table(dmrs1[,1:3],"/Users/xiaoyudai/Documents/multiple-testing/Comparison/bedtools/CD4_CD8/MCF005.txt",
            row.names=F, quote=FALSE, col.names=F, sep = "\t")
write.table(dmrs2[,1:3],"/Users/xiaoyudai/Documents/multiple-testing/Comparison/bedtools/CD4_CD8/Q005.txt",
            row.name=F, quote=FALSE, col.names=F, sep = "\t")


write.table(dmrs1[,1:8],"/Users/xiaoyudai/Documents/multiple-testing/Comparison/CD4_ES/MCF005.txt",row.name=F,quote=FALSE)
write.table(dmrs2[,1:8],"/Users/xiaoyudai/Documents/multiple-testing/Comparison/CD4_ES/Q005.txt",row.name=F,quote=FALSE)




