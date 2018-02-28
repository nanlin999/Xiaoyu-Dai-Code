all2 = read.table('~/Documents/Paper3/Data/Brain_ES_all/all2', sep = '\t', header = F)

all3 = all2[,c(1,2,3,4,8,12,13)]
all3$V4 = round(all3$V4 * all3$V8)
colnames(all3) = c('chr','start','end','brain_methy','brain_total','es_methy','es_total')
write.table(all3, '~/Documents/Paper3/Data/Brain_ES_all/all_raw', quote = F, sep = " ", col.names = T, row.names = F)
# write.table(all3, '~/Documents/Paper3/Data/Brain_ES_all/all3', quote = F, sep = "\t", col.names = F, row.names = F)


all_raw = read.table('~/Documents/Paper3/Data/Brain_ES_all/all_raw')


brain_mre = read.table('~/Documents/Paper3/Data/Brain/HuFGM02_BrainGerminalMatrix_MeDIP_A02758.bed', sep = '\t', header = F)
es_mre = read.table('~/Documents/Paper3/Data/H1ES/B1_H1Es_MRE_HS1052.bed', sep = '\t', header = F)

new = es_mre[,c(1,2,3)]
new$V3 = new$V2 + 2
new$V1 = as.character(new$V1)

for(i in 1:22){
  new[which(new$V1 == i),1] = paste('chr', i, sep = '')
}
new[which(new$V1 == 'X'),1] = 'chrX'
new[which(new$V1 == 'Y'),1] = 'chrY'
unique(new$V1)

new$V2 = format(new$V2, scientific = FALSE, trim = T)
new$V3 = format(new$V3, scientific = FALSE, trim = T)

write.table(new, '~/Documents/Paper3/Data/Brain_ES_all/brain_mre', quote = F, sep = "\t", col.names = F, row.names = F)
write.table(new, '~/Documents/Paper3/Data/Brain_ES_all/es_mre', quote = F, sep = "\t", col.names = F, row.names = F)

all_raw = read.table('~/Documents/Paper3/Data/Brain_ES_all/all7', sep = '\t', header = F)
colnames(all_raw) = c('chr','start','end','brain_methy','brain_total','es_methy','es_total','brain_medip','brain_mre','es_medip','es_mre')
write.table(all_raw, '~/Documents/Paper3/Data/Brain_ES_all/all_raw', quote = F, sep = " ", col.names = T, row.names = F)


id = which(all_raw$brain_total > 10 & all_raw$es_total > 10)

all = all_raw[id,]


dat1 <- all[,c('chr','start','brain_total','brain_methy')]
colnames(dat1) <- c('chr','pos','N','X')

dat2 <- all[,c('chr','start','es_total','es_methy')]
colnames(dat2) <- c('chr','pos','N','X')

BSobj <- makeBSseqData( list(dat1, dat2), c("C1", "N1") )
dmlTest.sm <- DMLtest(BSobj, group1=c("C1"), group2=c("N1"), smoothing=TRUE)
colnames(dmlTest.sm)[2] = 'start'

all1 = merge(all, dmlTest.sm, by = c('chr','start'))

write.table(all1, file = '/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_all/all_gf.bed', quote = F, sep = " ", col.names = T, row.names = F)


all[,'label'] = ifelse(all$fdr < 0.001, 1, 0)


all_gf = read.table('~/Documents/Paper3/Data/Brain_ES_all/crf/all_gf1.bed', header = T)


all_gf$brain_medip = all_gf$brain_medip/quantile(all_gf[which(all_gf$brain_medip>0),'brain_medip'], 0.75) * 10
all_gf$brain_mre = all_gf$brain_mre/quantile(all_gf[which(all_gf$brain_mre>0),'brain_mre'], 0.75) * 10
all_gf$es_medip = all_gf$es_medip/quantile(all_gf[which(all_gf$es_medip>0),'es_medip'], 0.75) * 10
all_gf$es_mre = all_gf$es_mre/quantile(all_gf[which(all_gf$es_mre>0),'es_mre'], 0.75) * 10

all_gf[,'dist'] = all_gf$start - c(all_gf$start[1], all_gf$start[1:nrow(all_gf)-1])

all_gf$brain_medip_10 = ma(all_gf$brain_medip, 10)
all_gf$brain_medip_50 = ma(all_gf$brain_medip, 50)
all_gf$brain_medip_200 = ma(all_gf$brain_medip, 200)
all_gf$brain_medip_1k = ma(all_gf$brain_medip, 1000)

all_gf$es_medip_10 = ma(all_gf$es_medip, 10)
all_gf$es_medip_50 = ma(all_gf$es_medip, 50)
all_gf$es_medip_200 = ma(all_gf$es_medip, 200)
all_gf$es_medip_1k = ma(all_gf$es_medip, 1000)

all_gf$brain_mre_10 = ma(all_gf$brain_mre, 10)
all_gf$brain_mre_50 = ma(all_gf$brain_mre, 50)
all_gf$brain_mre_200 = ma(all_gf$brain_mre, 200)
all_gf$brain_mre_1k = ma(all_gf$brain_mre, 1000)

all_gf$es_mre_10 = ma(all_gf$es_mre, 10)
all_gf$es_mre_50 = ma(all_gf$es_mre, 50)
all_gf$es_mre_200 = ma(all_gf$es_mre, 200)
all_gf$es_mre_1k = ma(all_gf$es_mre, 1000)

all_gf$ratio_medip = abs(all_gf$brain_medip - all_gf$es_medip) / min(all_gf$brain_medip+1, all_gf$es_medip+1)
all_gf$ratio_medip_10 = abs(all_gf$brain_medip_10 - all_gf$es_medip_10) / min(all_gf$brain_medip_10+1, all_gf$es_medip_10+1)
all_gf$ratio_medip_50 = abs(all_gf$brain_medip_50 - all_gf$es_medip_50) / min(all_gf$brain_medip_50+1, all_gf$es_medip_50+1)
all_gf$ratio_medip_200 = abs(all_gf$brain_medip_200 - all_gf$es_medip_200) / min(all_gf$brain_medip_200+1, all_gf$es_medip_200+1)
all_gf$ratio_medip_1k = abs(all_gf$brain_medip_1k - all_gf$es_medip_1k) / min(all_gf$brain_medip_1k+1, all_gf$es_medip_1k+1)

all_gf$ratio_mre = abs(all_gf$brain_mre - all_gf$es_mre) / min(all_gf$brain_mre+1, all_gf$es_mre+1)
all_gf$ratio_mre_10 = abs(all_gf$brain_mre_10 - all_gf$es_mre_10) / min(all_gf$brain_mre_10+1, all_gf$es_mre_10+1)
all_gf$ratio_mre_50 = abs(all_gf$brain_mre_50 - all_gf$es_mre_50) / min(all_gf$brain_mre_50+1, all_gf$es_mre_50+1)
all_gf$ratio_mre_200 = abs(all_gf$brain_mre_200 - all_gf$es_mre_200) / min(all_gf$brain_mre_200+1, all_gf$es_mre_200+1)
all_gf$ratio_mre_1k = abs(all_gf$brain_mre_1k - all_gf$es_mre_1k) / min(all_gf$brain_mre_1k+1, all_gf$es_mre_1k+1)

all_gf$dif = abs((all_gf$brain_mre+1)*(all_gf$es_medip+1) - (all_gf$es_mre+1)*(all_gf$brain_medip+1))
all_gf$dif_10 = abs((all_gf$brain_mre_10+1)*(all_gf$es_medip_10+1) - (all_gf$es_mre_10+1)*(all_gf$brain_medip_10+1))
all_gf$dif_50 = abs((all_gf$brain_mre_50+1)*(all_gf$es_medip_50+1) - (all_gf$es_mre_50+1)*(all_gf$brain_medip_50+1))
all_gf$dif_200 = abs((all_gf$brain_mre_200+1)*(all_gf$es_medip_200+1) - (all_gf$es_mre_200+1)*(all_gf$brain_medip_200+1))
all_gf$dif_1k = abs((all_gf$brain_mre_1k+1)*(all_gf$es_medip_1k+1) - (all_gf$es_mre_1k+1)*(all_gf$brain_medip_1k+1))


write.table(all1, file = '/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_all/all_gf.bed', quote = F, sep = " ", col.names = T, row.names = F)


###################################
all1 = read.table('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_all/all_gf.bed', header = T)
all2 = read.table('/Users/xiaoyudai/Documents/Paper3/Data/Skin_ES_all/all_gf.bed', header = T)

all1 = all1[,c(1:11,19:66)]


all1_chr18 = all1[which(all1$chr == 'chr18'), ]
all2_chr18 = all2[which(all2$chr == 'chr18'), ]

for(i in 28:59) {
  ratio = quantile(all1_chr18[which(all1_chr18[,i]>0),i], 0.75) / quantile(all2_chr18[which(all2_chr18[,i]>0),i], 0.75)
  all2_chr18[,i] = all2_chr18[,i] * ratio
}

colnames(all1_chr18) = NA
colnames(all2_chr18) = NA

all_gf = rbind(all1_chr18, all2_chr18)
colnames(all_gf) = colnames(all1)

label.threshold = 0.001
n_clusters = 1000
rnn_sub_length = 5


res = create_cluster_feature(all_gf, n_clusters)

train_crf = res$train_crf
test_crf = res$test_crf


# train_crf = train_crf[,c(1:16,20:21,25:26,30)]
# test_crf = test_crf[,c(1:16,20:21,25:26,30)]
write.table(train_crf, file = '/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_all/crf/train_crf', quote = F, sep = " ", col.names = F, row.names = F)
write.table(test_crf, file = '/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_all/crf/test_crf', quote = F, sep = " ", col.names = F, row.names = F)

result = read.table('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_all/crf/result')
pred = as.character(result[,31])
pred = as.numeric(substr(pred,1,1))
truth = as.numeric(test_crf$label)
truth = truth[which(truth!=' ')]
accuracy(truth, pred)
confusion.matrix(truth, pred)


prob0 = result[,ncol(result)-1]
prob0 = as.character(prob0)
prob0 = as.numeric(substr(prob0,3,10))

prob1 = result[,ncol(result)]
prob1 = as.character(prob1)
prob1 = as.numeric(substr(prob1,3,10)) 

auc(truth, prob1)




##################################################

binlength = 500

dirwrite <- "/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_all/mm/"

writefile <- paste(dirwrite, "brain_medip_bin.bed", sep = "")
countMeDIPbin('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_all/mm/brain_medip1', writefile = writefile, binlength = binlength)

writefile <- paste(dirwrite, "brain_mre_bin.bed", sep = "")
countMREbin('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_all/mm/brain_mre1', writefile = writefile, binlength = binlength)

writefile <- paste(dirwrite, "es_medip_bin.bed", sep = "")
countMeDIPbin('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_all/mm/es_medip1', writefile = writefile, binlength = binlength)

writefile <- paste(dirwrite, "es_mre_bin.bed", sep = "")
countMREbin('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_all/mm/es_mre1', writefile = writefile, binlength = binlength)

writefile = paste(dirwrite, "cpgbin.bed", sep = "")
countcpgbin('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_all/mm/all_bs', writefile = writefile, binlength = binlength)

file <- '~/Documents/Paper3/Data/TriMRE_frags.bed'
file1 <- '~/Documents/Paper3/Data/Brain_ES_all/mm/all_bs'
allcpgfile <- paste(dirwrite, "cpgbin.bed", sep = "")
five_Mre_CpGsite <- read.table(file, header = FALSE, as.is = TRUE)
four_Mre_CpGsite <- five_Mre_CpGsite[five_Mre_CpGsite[, 4] != "ACGT", ]
mrecpg.site <- four_Mre_CpGsite[four_Mre_CpGsite[, 4] != "CGCG", ]
writefile <- paste(dirwrite, "three_mre_cpg_bin.bed", sep = "")
countMREcpgbin(mrecpg.site, file.allcpgsite = file1, file.bin = allcpgfile,
               writefile = writefile, binlength = binlength)



brain_medip_bin = read.table('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_all/mm/brain_medip_bin.bed', header = T)
brain_mre_bin = read.table('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_all/mm/brain_mre_bin.bed', header = T)
es_medip_bin = read.table('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_all/mm/es_medip_bin.bed', header = T)
es_mre_bin = read.table('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_all/mm/es_mre_bin.bed', header = T)
cpg_bin = read.table('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_all/mm/cpgbin.bed', header = T)
three_mre_cpg_bin = read.table('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_all/mm/three_mre_cpg_bin.bed', header = T)

colnames(brain_medip_bin) = c('V1','V2','V3','brain_medip_bin')
colnames(brain_mre_bin) = c('V1','V2','V3','brain_mre_bin')
colnames(es_medip_bin) = c('V1','V2','V3','es_medip_bin')
colnames(es_mre_bin) = c('V1','V2','V3','es_mre_bin')
colnames(cpgbin) = c('V1','V2','V3','cpgbin')
colnames(three_mre_cpg_bin) = c('V1','V2','V3','three_mre_cpg_bin')

all = merge(brain_medip_bin, brain_mre_bin, by = c('V1','V2','V3'), sort = F)
all = merge(all, es_medip_bin, by = c('V1','V2','V3'), sort = F)
all = merge(all, es_mre_bin, by = c('V1','V2','V3'), sort = F)
all = merge(all, cpg_bin, by = c('V1','V2','V3'), sort = F)
all = merge(all, three_mre_cpg_bin, by = c('V1','V2','V3'), sort = F)

for(name in c('brain_medip_bin', 'brain_mre_bin', 'es_medip_bin', 'es_mre_bin', 'cpgbin', 'three_mre_cpg_bin') ){
  x = all[,c('V1','V2','V3',name)]
  path = paste('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_all/mm/', name, '1.bed', sep = '')
  colnames(x)[4] = 'V4'
  write.table(x, path, quote = F, sep = "\t", col.names = T, row.names = F)
}


datafile1 <- paste(dirwrite, "brain_medip_bin1.bed", sep = "")
datafile2 <- paste(dirwrite, "es_medip_bin1.bed", sep = "")
datafile3 <- paste(dirwrite, "brain_mre_bin1.bed", sep = "")
datafile4 <- paste(dirwrite, "es_mre_bin1.bed", sep = "")
datafile <- c(datafile1, datafile2, datafile3, datafile4)
cpgfile <- paste(dirwrite, "cpgbin1.bed", sep = "")
mrecpgfile <- paste(dirwrite, "three_mre_cpg_bin1.bed", sep = "")
writefile <- paste(dirwrite, "pval_Brain_es_all.bed", sep = "")
reportfile <- paste(dirwrite, "report_Brain_es_all.txt", sep = "")
MnM.test(file.dataset = datafile, file.cpgbin = cpgfile,
         file.mrecpgbin = mrecpgfile, writefile = writefile, reportfile = reportfile,
         mreratio = 3/7, method = "XXYY", psd = 2, mkadded = 1, a = 1e-16,
         cut = 100, top = 500)


datafile <- paste(dirwrite, "pval_Brain_es_all.bed", sep = "")
writefile <- paste(dirwrite, "q_Brain_es_all.bed", sep = "")
reportfile <- paste(dirwrite, "report_q_Brain_es_all.bed", sep = "")
MnM.qvalue(datafile, writefile, reportfile)

qval = read.table('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_all/mm/q_Brain_es_all.bed', header = T)


for( chrome in unique(all1$chr)) {
  print(chrome)
  id_all1 = which(all1$chr == chrome)
  id_qval = which(qval$chr == chrome)
  sub_qval = qval[id_qval,]
  
  idx = floor(all1[id_all1,'start']/500)+1
  
  all1[id_all1, 'mm_Medip1'] = sub_qval[idx, 'Medip1']
  all1[id_all1, 'mm_Medip2'] = sub_qval[idx, 'Medip2']
  all1[id_all1, 'mm_MRE1'] = sub_qval[idx, 'MRE1']
  all1[id_all1, 'mm_MRE2'] = sub_qval[idx, 'MRE2']
  all1[id_all1, 'mm_q'] = sub_qval[idx, 'qvalue']
}

write.table(all1, '~/Documents/Paper3/Data/Brain_ES_all/mm/all', quote = F, sep = " ", col.names = T, row.names = F)

all1 = read.table('~/Documents/Paper3/Data/Brain_ES_all/mm/all', header = T)








