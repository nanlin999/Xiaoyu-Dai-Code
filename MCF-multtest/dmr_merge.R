ALL <- read.table("/Users/xiaoyudai/Documents/multiple-testing/Comparison/CD4_CD8/ALL.txt", header=T)

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
dmrs1 <- callDMR(data1, p.threshold=0.01,minlen=1000,dis.merge=100,pct.sig=0.7, minCG=20 )
dim(dmrs1)

data2 = data
data2[,1:2] = ALL[,1:2]
data2[,10] = 1-Qvalue
data2[,3] = (ALL[,4]+ALL[,6])/(ALL[,5]+ALL[,7])
data2[,4] = (ALL[,8]+ALL[,10])/(ALL[,9]+ALL[,11])
data2[,5] = data2[,3]-data2[,4]
dmrs2 <- callDMR(data2, p.threshold=0.01,minlen=1500,dis.merge=100,pct.sig=0.7, minCG=20)
dim(dmrs2)

write.table(dmrs1[,1:8],"/Users/xiaoyudai/Documents/multiple-testing/Comparison/CD8_ES/MCF010.txt",row.name=F,quote=FALSE)
write.table(dmrs2[,1:8],"/Users/xiaoyudai/Documents/multiple-testing/Comparison/CD8_ES/Q010.txt",row.name=F,quote=FALSE)

dmrs1 <- read.table("/Users/xiaoyudai/Documents/multiple-testing/Comparison/CD4_ES/MCF005.txt",header  = T)
dmrs2 <- read.table("/Users/xiaoyudai/Documents/multiple-testing/Comparison/CD4_ES/Q005.txt",header = T)

write.table(dmrs1[,1:3],"/Users/xiaoyudai/Documents/multiple-testing/Comparison/bedtools/CD8_ES/MCF010.txt",
            col.names = F, row.name=F,quote=FALSE,sep='\t')
write.table(dmrs2[,1:3],"/Users/xiaoyudai/Documents/multiple-testing/Comparison/bedtools/CD8_ES/Q010.txt",
            col.names = F, row.name=F,quote=FALSE,sep='\t')


inter <- read.table("/Users/xiaoyudai/Documents/multiple-testing/Comparison/bedtools/CD8_ES/inter.txt", header=T)
additional <- read.table("/Users/xiaoyudai/Documents/multiple-testing/Comparison/bedtools/CD8_ES/additional.txt", header=T)

dim(inter)
dim(additional)
dim(dmrs1)
