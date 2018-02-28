source("~/Documents/Paper3/Rcode/FDR-HMM/HMM-FDR/mt.hmm.R")
library(SDMTools)
library(DSS)
library(discretization)
library(infotheo)
# library(smbinning)


medip_1 = read.table('/Users/xiaoyudai/Documents/Paper3/Data/chr18/MeDIP_1.bed')
medip_2 = read.table('/Users/xiaoyudai/Documents/Paper3/Data/chr18/MeDIP_2.bed')
mre_1 = read.table('/Users/xiaoyudai/Documents/Paper3/Data/chr18/MRE_1.bed')
mre_2 = read.table('/Users/xiaoyudai/Documents/Paper3/Data/chr18/MRE_2.bed')


colnames(Brain_density) = c('chr', 'start', 'end', 'brain_density')
colnames(Brain_ratio) = c('chr', 'start', 'end', 'brain_ratio')
colnames(Skin_density) = c('chr', 'start', 'end', 'skin_density')
colnames(Skin_ratio) = c('chr', 'start', 'end', 'skin_ratio')

merge1 = merge(Brain_density, Brain_ratio, by = c('chr','start','end'))
merge2 = merge(Skin_density, Skin_ratio, by = c('chr','start','end'))

bis_all = merge(merge1, merge2, by = c('chr','start','end'))
write.table(bis_all, file = '/Users/xiaoyudai/Documents/Paper3/Data/train/bis_all_new.txt', quote = F, sep = " ", col.names = T, row.names = F)
bis_all = read.table('/Users/xiaoyudai/Documents/Paper3/Data/train/bis_all.txt', header = T)

dirwrite <- "/Users/xiaoyudai/Documents/Paper3/Data/chr18/crf/"

allcpgfile <- "/Users/xiaoyudai/Documents/Paper3/Data/chr18/density_1.bed"

density_1 = read.table(allcpgfile, header = F)
write.table(density_1, allcpgfile, quote = F, sep = " ", col.names = F, row.names = F)

writefile <- paste(dirwrite, "MeDIP_bin_1.bed", sep = "")
countMeDIPbin('/Users/xiaoyudai/Documents/Paper3/Data/chr18/MeDIP_1.bed', writefile = writefile, binlength = 2)

writefile <- paste(dirwrite, "MRE_bin_1.bed", sep = "")
countMREbin('/Users/xiaoyudai/Documents/Paper3/Data/chr18/MRE_1.bed', writefile = writefile, binlength = 2)

writefile <- paste(dirwrite, "MeDIP_bin_2.bed", sep = "")
countMeDIPbin('/Users/xiaoyudai/Documents/Paper3/Data/chr18/MeDIP_2.bed', writefile = writefile, binlength = 2)

writefile <- paste(dirwrite, "MRE_bin_2.bed", sep = "")
countMREbin('/Users/xiaoyudai/Documents/Paper3/Data/chr18/MRE_2.bed', writefile = writefile, binlength = 2)


MRE_2 = read.table('/Users/xiaoyudai/Documents/Paper3/Data/chr18/crf/MRE_bin_2.bed', header = T)
colnames(MRE_2) = c('chr', 'start', 'end', 'skin_mre')


MRE_2$chr = 'chr18'
bis_all = merge(bis_all, MRE_2, by = c('chr','start','end'))

bis_all = read.table('/Users/xiaoyudai/Documents/Paper3/Data/train/bis_all.txt', header = T)

dat1 <- all[,c('chr','start','brain_density','brain_ratio')]
colnames(dat1) <- c('chr','pos','N','X')
dat1$X = round(dat1$N * dat1$X)

dat2 <- all[,c('chr','start','skin_density','skin_ratio')]
colnames(dat2) <- c('chr','pos','N','X')
dat2$X = round(dat2$N * dat2$X)

BSobj <- makeBSseqData( list(dat1, dat2), c("C1", "N1") )
dmlTest.sm <- DMLtest(BSobj, group1=c("C1"), group2=c("N1"), smoothing=TRUE)
colnames(dmlTest.sm)[2] = 'start'
# dmls <- callDML(dmlTest.sm, p.threshold=0.001)


all = read.table('/Users/xiaoyudai/Documents/Paper3/Data/chr18/crf/all.bed', header = F)
colnames(all) = c('chr','start','end','brain_density','brain_ratio','skin_density','skin_ratio','brain_mre','brain_medip','skin_mre','skin_medip')
all = merge(all, dmlTest.sm, by = c('chr','start'))
all = all[order(all$start),]

idx = which(all$brain_density > 10 & all$skin_density > 10)
length(idx)
all = all[idx,]

# impute <- function(x) {
#   res = x  
#   id0 = which(x==0)
#   for (i in id0) {
#     print(i)
#     id1 = max(which(x[1:i] != 0))
#     id2 = min(which(x[i+1:length(x)] != 0))
#     res[i] = (x[id1] + x[id2]) / 2
#   }
#   return(res)
# }
# 
# all$brain_medip = impute(all$brain_medip)
# all$skin_medip = impute(all$skin_medip)
# write.table(all, file = '/Users/xiaoyudai/Documents/Paper3/Data/chr18/crf/all_impute.txt', quote = F, sep = " ", col.names = T, row.names = F)
# 
# all$brain_mre = impute(all$brain_mre)
# all$skin_mre = impute(all$skin_mre)
# 
# write.table(all, file = '/Users/xiaoyudai/Documents/Paper3/Data/chr18/crf/all_impute.txt', quote = F, sep = " ", col.names = T, row.names = F)
# all = read.table('/Users/xiaoyudai/Documents/Paper3/Data/chr18/crf/all_impute.txt', header = T)

n_train_seq = 100
l_train_seq = 290
label.threshold=0.0001
n_clusters=30

# source('~/Documents/Paper3/Data/code/split_with_augmented_features.r')
res = split_with_augmented_features(all, n_train_seq, l_train_seq, label.threshold, n_clusters)
train_crf = res$train_crf
test_crf = res$test_crf
all = res$all
train_rnn = res$train_rnn
test_rnn = res$test_rnn

write.table(train_crf, file = '/Users/xiaoyudai/Documents/Paper3/Data/chr18/crf/train_crf', quote = F, sep = " ", col.names = F, row.names = F)
write.table(test_crf, file = '/Users/xiaoyudai/Documents/Paper3/Data/chr18/crf/test_crf', quote = F, sep = " ", col.names = F, row.names = F)
write.table(train_rnn, file = '/Users/xiaoyudai/Documents/Paper3/Data/chr18/crf/train_rnn', quote = F, sep = "\t", col.names = F, row.names = F)
write.table(test_rnn, file = '/Users/xiaoyudai/Documents/Paper3/Data/chr18/crf/test_rnn', quote = F, sep = "\t", col.names = F, row.names = F)


result = read.table('/Users/xiaoyudai/Documents/Paper3/Data/chr18/crf/result')
pred = as.character(result[,ncol(result)-2])
pred = as.numeric(substr(pred,1,1))
truth = test_crf$label
accuracy(truth, pred)
confusion.matrix(truth, pred)


prob0 = result[,ncol(result)-1]
prob0 = as.character(prob0)
prob0 = as.numeric(substr(prob0,3,10))

prob1 = result[,ncol(result)]
prob1 = as.character(prob1)
prob1 = as.numeric(substr(prob1,3,10)) 
# accuracy(truth, prob1)
# confusion.matrix(truth, prob1)

lsi.crf = prob0/prob1

q = 0.03
res.crf <- mt.hmm(lsi.crf, q)$de
sum(res.crf)

test_all <- all[,c('chr','start','mu1','mu2','diff','diff.se','stat','phi1','phi2','pval','fdr')]
colnames(test_all)[2] = 'pos'
test_all[,'pval'] = 1 - res.crf
dmr_crf = callDMR(test_all, delta=0, p.threshold=0.01, minlen = 100, pct.sig = 0.5, minCG = 2)
dim(dmr_crf)


qval = read.table('~/Documents/Paper3/Data/chr18/q_Brain_Skin_chr18.bed', header = T)
q_mnm = 0.01
# dmr_mnm = qval[which(qval$qvalue < q_mnm & qval$chrSt>test_all[1,'pos']),]
dmr_mnm = qval[which(qval$qvalue < q_mnm),]
dim(dmr_mnm)




a = dmlTest.sm
a$pval = 1
dim(callDML(a, delta=0, p.threshold=0.0001))



accuracy(truth, res.crf)
confusion.matrix(truth, res.crf)

accuracy(res.MnM, truth)
confusion.matrix(res.MnM, truth)

mt_eval(res.crf, truth)
mt_eval(res.MnM, truth)
confusion.matrix(res.crf, res.MnM)


mt_eval <- function(pred, true){
  false_positive = sum(pred==1 & true==0)
  false_negative = sum(pred==0 & true==1)
  fdr = false_positive / sum(pred==1)
  fnr = false_negative / sum(pred==0)
  return(c(fdr, fnr))
}


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






















qval = read.table('~/Documents/Paper3/Data/chr18/q_Brain_Skin_chr18.bed', header = T)
res.MnM <- rep(0, 295966)
res.MnM[qval[,'qvalue'] < q] = 1

res.MnM = res.MnM[100001:295966]
res.MnM = res.MnM[!res.MnM == ' ']

sum(res.MnM)

length(res.MnM)















label[1000:3000]

m = 0
l = 0
for (i in 1:295966) {
  if(label[i] == 1){
    l = l+1
  } else {
    l = 0
  }
  m = max(m,l)
}

which(label==1)


























