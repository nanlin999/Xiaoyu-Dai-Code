source("~/Documents/Paper3/Rcode/FDR-HMM/HMM-FDR/mt.hmm.R")
library(SDMTools)
library(DSS)
library(discretization)
library(infotheo)
library(ggplot2)
library(pROC)
library("GenomicRanges")
library(methylMnM)

brain_cpg = read.table('~/Documents/Paper3/Data/Brain_ES_chr18/brain_cpg.bed')
brain_density = read.table('~/Documents/Paper3/Data/Brain_ES/brain_density.bed')
es_bs = read.table('~/Documents/Paper3/Data/Brain_ES/es_bs.bed')

colnames(brain_cpg) = c('chr', 'start', 'end', 'brain_cpg')
colnames(brain_density) = c('chr', 'start', 'end', 'brain_density')
colnames(es_bs) = c('chr', 'start', 'end', 'es_methy', 'es_total')

all_bs = merge(brain_cpg, brain_density, by=c('chr','start','end'))
all_bs = merge(all_bs, es_bs, by=c('chr','start','end'))
all_bs = all_bs[order(all_bs$start),]
all_bs$brain_methy = round(all_bs$brain_cpg*all_bs$brain_density)

write.table(all_bs[,c('chr','start','end','brain_methy','brain_density','es_methy','es_total')], 
            file = '/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES/all_bs.bed', quote = F, sep = " ", col.names = F, row.names = F)




binlength = 2000

dirwrite <- "/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_chr18/bin2000/"

writefile <- paste(dirwrite, "brain_medip_bin.bed", sep = "")
countMeDIPbin('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_chr18/brain_medip.bed', writefile = writefile, binlength = binlength)

writefile <- paste(dirwrite, "brain_mre_bin.bed", sep = "")
countMREbin('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_chr18/brain_mre.bed', writefile = writefile, binlength = binlength)

writefile <- paste(dirwrite, "es_medip_bin.bed", sep = "")
countMeDIPbin('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_chr18/es_medip.bed', writefile = writefile, binlength = binlength)

writefile <- paste(dirwrite, "es_mre_bin.bed", sep = "")
countMREbin('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_chr18/es_mre.bed', writefile = writefile, binlength = binlength)

writefile = paste(dirwrite, "cpgbin.bed", sep = "")
countcpgbin('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_chr18/all_bs.bed', writefile = writefile, binlength = binlength)

file <- '~/Documents/Paper3/Data/TriMRE_frags.bed'
file1 <- '~/Documents/Paper3/Data/Brain_ES_chr18/all_bs.bed'
allcpgfile <- paste(dirwrite, "cpgbin.bed", sep = "")
five_Mre_CpGsite <- read.table(file, header = FALSE, as.is = TRUE)
four_Mre_CpGsite <- five_Mre_CpGsite[five_Mre_CpGsite[, 4] != "ACGT", ]
mrecpg.site <- four_Mre_CpGsite[four_Mre_CpGsite[, 4] != "CGCG", ]
writefile <- paste(dirwrite, "three_mre_cpg_bin.bed", sep = "")
countMREcpgbin(mrecpg.site, file.allcpgsite = file1, file.bin = allcpgfile,
               writefile = writefile, binlength = binlength)



datafile1 <- paste(dirwrite, "brain_medip_bin1.bed", sep = "")
datafile2 <- paste(dirwrite, "es_medip_bin1.bed", sep = "")
datafile3 <- paste(dirwrite, "brain_mre_bin1.bed", sep = "")
datafile4 <- paste(dirwrite, "es_mre_bin1.bed", sep = "")
datafile <- c(datafile1, datafile2, datafile3, datafile4)
cpgfile <- paste(dirwrite, "cpgbin1.bed", sep = "")
mrecpgfile <- paste(dirwrite, "three_mre_cpg_bin1.bed", sep = "")
writefile <- paste(dirwrite, "pval_Brain_es_chr18.bed", sep = "")
reportfile <- paste(dirwrite, "report_Brain_es_chr18.txt", sep = "")
MnM.test(file.dataset = datafile, file.cpgbin = cpgfile,
         file.mrecpgbin = mrecpgfile, writefile = writefile, reportfile = reportfile,
         mreratio = 3/7, method = "XXYY", psd = 2, mkadded = 1, a = 1e-16,
         cut = 100, top = 500)


datafile <- paste(dirwrite, "pval_Brain_es_chr18.bed", sep = "")
writefile <- paste(dirwrite, "q_Brain_es_chr18.bed", sep = "")
reportfile <- paste(dirwrite, "report_q_Brain_es_chr18.bed", sep = "")
MnM.qvalue(datafile, writefile, reportfile)

qval = read.table(writefile, header = T)

idx = floor(all_mm1[,'start']/binlength)+1
all_mm1[,'mm_q_2000'] = qval$qvalue[idx]

pred_mm = ifelse(all_mm1$mm_q_2000<0.01, 1, 0)
truth = all_mm1$label
confusion.matrix(truth, pred_mm) # 3386 7570



qval = read.table('~/Documents/Paper3/Data/Brain_ES_chr18/bin2000/q_Brain_es_chr18.bed', header = T)

idx = floor(all_mm1[,'start']/2000)+1
all_mm1$Medip1_2000 = qval$Medip1[idx]
all_mm1$Medip2_2000 = qval$Medip2[idx]
all_mm1$MRE1_2000 = qval$MRE1[idx]
all_mm1$MRE2_2000 = qval$MRE2[idx]



all_mm1 = all_mm1[,c(1:27,59:64)]
colnames(all_mm1)



file <- '~/Documents/Paper3/Data/TriMRE_frags.bed'
file1 <- '~/Documents/Paper3/Data/Brain_ES/all_bs.bed'
allcpgfile <- '~/Documents/Paper3/Data/Brain_ES/all_bs.bed'
five_Mre_CpGsite <- read.table(file, header = FALSE, as.is = TRUE)
four_Mre_CpGsite <- five_Mre_CpGsite[five_Mre_CpGsite[, 4] != "ACGT", ]
mrecpg.site <- four_Mre_CpGsite[four_Mre_CpGsite[, 4] != "CGCG", ]
writefile <- paste(dirwrite, "three_mre_cpg.bed", sep = "")
countMREcpgbin(mrecpg.site, file.allcpgsite = file1, # file.bin = allcpgfile,
               writefile = writefile, binlength = 2)



all_bs = read.table('~/Documents/Paper3/Data/Brain_ES_chr18/all_bs.bed')
all_bs$V1 = 18
write.table(all_bs, '~/Documents/Paper3/Data/Brain_ES_chr18/all_bs1.bed', quote = F, sep = "\t", col.names = F, row.names = F)

all = read.table('~/Documents/Paper3/Data/Brain_ES_chr18/all_bs5.bed')
colnames(all) = c('chr','start','end','brain_methy','brain_total','es_methy','es_total','brain_medip','brain_mre','es_medip','es_mre')

idx = which(all$brain_total>10 & all$es_total>10)
all = all[idx,]

dat1 <- all[,c('chr','start','brain_total','brain_methy')]
colnames(dat1) <- c('chr','pos','N','X')

dat2 <- all[,c('chr','start','es_total','es_methy')]
colnames(dat2) <- c('chr','pos','N','X')

BSobj <- makeBSseqData( list(dat1, dat2), c("C1", "N1") )
dmlTest.sm <- DMLtest(BSobj, group1=c("C1"), group2=c("N1"), smoothing=TRUE)
colnames(dmlTest.sm)[2] = 'start'

all = merge(all, dmlTest.sm, by = c('chr','start'))
all = all[order(all$start),]

mre_cpg = read.table('~/Documents/Paper3/Data/Brain_ES_chr18/three_mre_cpg.bed', header = T)
colnames(mre_cpg) = c('chr','start','end','is_mre_cpg')
mre_cpg$chr = 18

all1 = merge(all, mre_cpg, by = c('chr','start','end'), all.x = TRUE)
idx = is.na(all1$is_mre_cpg)
all1[idx,'is_mre_cpg'] = 0

write.table(all1, '~/Documents/Paper3/Data/Brain_ES_chr18/crf/all1.bed', quote = F, sep = " ", col.names = T, row.names = F)
all = read.table('~/Documents/Paper3/Data/Brain_ES/crf/all_bs10.bed', header = T)


# idx = which(all$brain_total>10 & all$es_total>10 & all$is_mre_cpg==1 & all$brain_medip>2 & all$es_medip>2)
# length(idx)
# all = all[idx,]

#################################################
#################################################
#################################################
#################################################
#################################################

all_raw = read.table('~/Documents/Paper3/Data/Brain_ES_chr18/crf/all_raw', header = T)
add = all_raw[which(all_raw$brain_total>10 & all_raw$es_total>10),]
add = add[order(add$start),c(8:12)]

all_mm = cbind(all_gf, add)
all_mm = all_mm[1:323583,]

write.table(all_mm1, '~/Documents/Paper3/Data/Brain_ES_chr18/crf/all_mm', quote = F, sep = " ", col.names = T, row.names = F)
write.table(all_mm2, '~/Documents/Paper3/Data/Skin_ES_chr18/crf/all_mm', quote = F, sep = " ", col.names = T, row.names = F)

# all_gf = read.table('~/Documents/Paper3/Data/Brain_ES_chr18/crf/all_gf.bed', header = T)

all_mm1 = read.table('~/Documents/Paper3/Data/Brain_ES_chr18/crf/all_mm', header = T)
all_mm2 = read.table('~/Documents/Paper3/Data/Skin_ES_chr18/crf/all_mm', header = T)

raw_1 = read.table('~/Documents/Paper3/Data/Brain_ES_chr18/all_bs5.bed')
raw_2 = read.table('~/Documents/Paper3/Data/Skin_ES_chr18/all_bs5')

colnames(raw_1) = c("chr", "start", "end", "brain_methy", "brain_total", "es_methy", "es_total", "brain_medip", "brain_mre", "es_medip", "es_mre")
colnames(raw_2) = c("chr", "start", "end", "skin_methy", "skin_total", "es_methy", "es_total", "skin_medip", "skin_mre", "es_medip", "es_mre")

merge1 = merge(all_mm1, raw_1, by = c('chr', 'start', 'end', 'brain_methy', "brain_total", "es_methy", "es_total"))
merge2 = merge(all_mm2, raw_2, by = c('chr', 'start', 'end', 'skin_methy', "skin_total", "es_methy", "es_total"))

merge1 = merge1[order(merge1$start),]
merge2 = merge2[order(merge2$start),]

colnames(merge1)[8:11] = c('brain_medip', 'brain_mre', 'es_medip', 'es_mre')
merge1$brain_medip = merge1$brain_medip.y
merge1$brain_mre = merge1$brain_mre.y
merge1$es_medip = merge1$es_medip.y
merge1$es_mre = merge1$es_mre.y

colnames(merge2)[8:11] = c('skin_medip', 'skin_mre', 'es_medip', 'es_mre')
merge2$skin_medip = merge2$skin_medip.y
merge2$skin_mre = merge2$skin_mre.y
merge2$es_medip = merge2$es_medip.y
merge2$es_mre = merge2$es_mre.y

write.table(merge1[,c(1:63)], '~/Documents/Paper3/Data/Brain_ES_chr18/crf/all_mm', quote = F, sep = " ", col.names = T, row.names = F)
write.table(merge2[,c(1:63)], '~/Documents/Paper3/Data/Skin_ES_chr18/crf/all_mm', quote = F, sep = " ", col.names = T, row.names = F)

###################################################

all_mm1 = read.table('~/Documents/Paper3/Data/Brain_ES_chr18/crf/all_mm_1', header = T)
all_mm2 = read.table('~/Documents/Paper3/Data/Skin_ES_chr18/crf/all_mm', header = T)



feature_names <- c("chr",             "start",           "end",             "pval",            "fdr",             "is_mre_cpg",     
                   "label",           "dist",            "rmsk_LINE",       "rmsk_LowSimple",  "rmsk_SINE",       "rmsk_DNA",        "rmsk_other",     
                   "utr3",            "utr5",            "cpgi",            "intron",          "exon",            "refGene",         "ratio_medip",    
                   "ratio_mre",       "dif",             "ratio_medip_10",  "ratio_mre_10",    "dif_10",          "ratio_medip_50",  "ratio_mre_50",   
                   "dif_50",          "ratio_medip_200", "ratio_mre_200",   "dif_200",         "ratio_medip_1k",  "ratio_mre_1k",    "dif_1k",         
                   "mm_q_500",  "mm_q_10",     "mm_q_50",    "mm_q_200",       "mm_q_1000",   "mm_q_2000" )

feature_names1 = c("chr",             "start",           "end",             "pval",            "fdr",             "is_mre_cpg",     
                   "label",           "dist",            "rmsk_LINE",       "rmsk_LowSimple",  "rmsk_SINE",       "rmsk_DNA",        "rmsk_other",     
                   "utr3",            "utr5",            "cpgi",            "intron",          "exon",            "refGene",          
                   "brain_medip",     "brain_medip_10",  "brain_medip_50",  "brain_mre",       "brain_mre_10",    "brain_mre_50",  
                   "es_medip",        "es_medip_10",     "es_medip_50",     "es_mre",          "es_mre_10",       "es_mre_50",
                   "mm_q_500",  "mm_q_10",     "mm_q_50",    "mm_q_200",       "mm_q_1000",   "mm_q_2000" )

feature_names2 = c("chr",             "start",           "end",             "pval",            "fdr",             "is_mre_cpg",     
                   "label",           "dist",            "rmsk_LINE",       "rmsk_LowSimple",  "rmsk_SINE",       "rmsk_DNA",        "rmsk_other",     
                   "utr3",            "utr5",            "cpgi",            "intron",          "exon",            "refGene",          
                   "skin_medip",      "skin_medip_10",   "skin_medip_50",   "skin_mre",        "skin_mre_10",     "skin_mre_50",  
                   "es_medip",        "es_medip_10",     "es_medip_50",     "es_mre",          "es_mre_10",       "es_mre_50",
                   "mm_q_500",  "mm_q_10",     "mm_q_50",    "mm_q_200",       "mm_q_1000",   "mm_q_2000" )


# comb_a = all_mm1[,feature_names]
# comb_b = all_mm2[,feature_names]


comb_a = all1[which(all1$chr=='chr18'),]
comb_b = all1[which(all1$chr=='chr17'),]

comb_a = all_mm1[,c(1:3,12:57)]
comb_b = all_mm2[,c(1:3,12:57)]


cl = colnames(comb_a)
colnames(comb_a) = NA
colnames(comb_b) = NA

all_mm = rbind(comb_a, comb_b)
colnames(all_mm) = cl

# all_split = read.table('~/Documents/Paper3/Data/Brain_ES_chr18/crf/all_split_1000.bed', header = T, blank.lines.skip = F)


label.threshold = 0.001
n_clusters = 200
rnn_sub_length = 10


# source('~/Documents/Paper3/Data/code/split_with_augmented_features.r')
# res = split_with_augmented_features(all_gf, label.threshold, n_clusters, genome_feature)
# res = create_cluster_feature(all_split, n_clusters, every_n)
res = create_cluster_feature(all_split, n_clusters)

train_crf = res$train_crf
test_crf = res$test_crf
# all = res$all
# write.table(all, file = '/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES/crf/all_sm.bed', quote = F, sep = " ", col.names = T, row.names = F)

# train_rnn = res$train_rnn
# test_rnn = res$test_rnn

# idx = which(train_crf[1:nrow(all_gf),]$label != 1)
# train_tmp = train_crf[1:nrow(all_gf),]
# train_tmp[idx,] = c(' ')
# for(i in 1:10) {
#   print(i)
#   train_crf = rbind(train_crf, train_tmp)
# }

train_crf[which(is.na(train_crf$ratio_medip_class)), ] = ' '
test_crf[which(is.na(test_crf$ratio_medip_class)), ] = ' '


write.table(train_crf, file = '/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_chr18/crf/train_crf', quote = F, sep = " ", col.names = F, row.names = F)
write.table(test_crf, file = '/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_chr18/crf/test_crf', quote = F, sep = " ", col.names = F, row.names = F)

train_idx = which(train_crf$label[1:nrow(all_mm)] != ' ')
test_mlp = all_mm[, c(14,16:27,44:63,15)]
for(i in c(2,14:33)) {
  test_mlp[,i] = scaling(test_mlp[,i])
}
train_mlp = test_mlp[train_idx,]

write.table(train_mlp, file = '/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_chr18/crf/train_mlp', quote = F, sep = "\t", col.names = F, row.names = F)
write.table(test_mlp, file = '/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_chr18/crf/test_mlp', quote = F, sep = "\t", col.names = F, row.names = F)

write.table(train_rnn, file = '/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_chr18/crf/train_rnn', quote = F, sep = "\t", col.names = F, row.names = F)
write.table(test_rnn, file = '/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_chr18/crf/test_rnn', quote = F, sep = "\t", col.names = F, row.names = F)


result = read.table('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_chr18/crf/result')
pred = as.character(result[,35])
pred = as.numeric(substr(pred,1,1))
truth = as.numeric(test_crf$label)
truth = truth[which(truth!=' ')]
accuracy(truth, pred)
confusion.matrix(truth, pred)


truth = ifelse(all1$fdr < 0.001, 1, 0)
mean(truth)

pred_mm = ifelse(all_17$mm_q<0.1, 1, 0)
# truth = all_mm1$label
# test_idx = floor(323583/5*4):323583
accuracy(truth, pred_mm)
confusion.matrix(truth, pred_mm) # 3386 7570  # 8331 9364
                                 # 632  1943  # 1047 2612


train_idx = which(train_crf$label[1:nrow(all_mm)] != ' ')
test_idx = which(train_crf$label[1:nrow(all_mm)] == ' ')
accuracy(truth[train_idx], pred[train_idx])
confusion.matrix(truth[train_idx], pred[train_idx])
accuracy(truth[test_idx], pred[test_idx])
confusion.matrix(truth[test_idx], pred[test_idx])


p_train = train_crf$label[1:nrow(all)]
p_train[which(p_train == ' ')] = 0
auc(truth, as.numeric(p_train))


test_seq = c(1:(train_seq[1,1]-1))
for(i in 1:(nrow(train_seq)-1) ){
  test_seq = c(test_seq, train_seq[i,2]:train_seq[i+1,1])
}
test_seq = c(test_seq, (train_seq[n_train_seq, 2]+1):nrow(test_crf))
truth = truth[test_seq]
pred = pred[test_seq]


# pred = mnm_label + pred
# pred = pred%%2
# truth = all_sm$label

mean(abs(pred-truth))
hist(abs(pred-truth))
mean( (pred-truth) < 2 & (pred-truth)> -2 )
pred_label = ifelse(pred > 8, 1, 0)
truth_label = ifelse(truth > 8, 1, 0)
accuracy(truth_label, pred_label)
confusion.matrix(truth_label, pred_label)

prob0 = result[,ncol(result)-1]
prob0 = as.character(prob0)
prob0 = as.numeric(substr(prob0,3,10))

prob1 = result[,ncol(result)]
prob1 = as.character(prob1)
prob1 = as.numeric(substr(prob1,3,10)) 

auc(truth, prob1)

auc(truth, 1-all_17$mm_q) # pval: 0.5901 fdr: 0.6083
                               # pval: 0.5386

hist(prob0, breaks=100)
hist(prob0[which(prob0<0.9)], breaks=20)
summary(prob0[which(prob0 < 0.9)])


test_idx = which(train_crf$label[1:nrow(all_gf)] == ' ')
prob0_test = prob0[test_idx]
hist(prob0_test, breaks=100)
hist(prob0_test[which(prob0_test<0.9)], breaks=20)


lsi.crf = prob0/prob1
hist(lsi.crf[which(lsi.crf<500)])

lsi.crf = prob0
q = 0.1
res.crf <- mt.hmm(lsi.crf, q)$de
sum(res.crf)
1 - sum(res.crf * truth) / sum(res.crf)
c(sum(res.crf) - sum(res.crf * truth), sum(res.crf * truth))


# test_all <- all_gf[,c('chr','start','mu1','mu2','diff','diff.se','stat','phi1','phi2','pval','fdr')]

test_all <- all_mm2[,c('chr','start')]
test_all[,'mu1'] = 0
test_all[,'mu2'] = 0
test_all[,'diff'] = 0
test_all[,'diff.se'] = 0
test_all[,'stat'] = 0
test_all[,'phi1'] = 0
test_all[,'phi2'] = 0
test_all[,'pval'] = all_mm2$pval
test_all[,'fdr'] = all_mm2$fdr

colnames(test_all)[2] = 'pos'

test_all[,'pval'] = 1 - res.crf
dmr_crf = callDMR(test_all, delta=0, p.threshold=0.01, minlen = 300, pct.sig = 0.7, minCG = 10)
dim(dmr_crf)


qval = read.table('~/Documents/Paper3/Data/Skin_ES_chr18/q_Skin_ES_chr18.bed', header = T)
hist(qval$qvalue)
hist(qval$qvalue[which(qval$qvalue < 0.9)])
q_mnm = 0.1
dmr_mnm = qval[which(qval$qvalue < q_mnm & qval$chrSt>test_all[1,'pos']),]
# dmr_mnm = qval[which(qval$qvalue < q_mnm),]
dim(dmr_mnm)


idx = floor(all_sm[,'start']/500)+1
p_mnm = 1 - qval$qvalue[idx]
p_mnm[which(is.na(p_mnm))] = 0
auc(truth, p_mnm)




dmr_mnm = qval[which(qval$qvalue < 0.6), ]
mnm_label = rep(0, nrow(all_sm))
for (i in 1:nrow(dmr_mnm)) {
  start = dmr_mnm[i,'chrSt']
  end = dmr_mnm[i,'chrEnd']
  idx = which(all_sm$start >= start & all_sm$end <= end)
  mnm_label[idx] = 1
}
# truth = ifelse(all_sm$fdr < 0.001, 1, 0)
truth = test_crf$label
accuracy(truth, mnm_label)
confusion.matrix(truth, mnm_label)





result_mlp = read.table('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_chr18/crf/result_mlp')
pred = result_mlp$V2 > 0.5

accuracy(truth, pred)
confusion.matrix(truth, pred)
accuracy(truth[train_idx], pred[train_idx])
confusion.matrix(truth[train_idx], pred[train_idx])
accuracy(truth[test_idx], pred[test_idx])
confusion.matrix(truth[test_idx], pred[test_idx])



prob0.mlp = result_mlp$V1
prob1.mlp = result_mlp$V2

lsi.mlp = prob0.mlp

q = 0.1
res.mlp <- mt.hmm(lsi.mlp, q)$de
sum(res.mlp)

test_all[,'pval'] = 1 - res.mlp
dmr_mlp = callDMR(test_all, delta=0, p.threshold=0.01, minlen = 100, pct.sig = 0.7, minCG = 5)
dim(dmr_mlp)




prob1.rnn = read.table('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_chr18/crf/result_rnn')
prob1.rnn = prob1.rnn$V1
prob0.rnn = 1 - prob1.rnn

lsi.rnn = prob0.rnn/prob1.rnn

q = 0.0001
res.rnn = rep(0, nrow(test_crf))
res.rnn[1:nrow(result_rnn)] <- mt.hmm(lsi.rnn, q)$de
sum(res.rnn)

test_all <- all_sm[round(0.2*nrow(all)):nrow(all), 
                   c('chr','start','mu1','mu2','diff','diff.se','stat','phi1','phi2','pval','fdr')]
colnames(test_all)[2] = 'pos'
test_all[,'pval'] = 1 - res.rnn
dmr_rnn = callDMR(test_all, delta=0, p.threshold=0.01, minlen = 100, pct.sig = 0.5, minCG = 1)
dim(dmr_rnn)







