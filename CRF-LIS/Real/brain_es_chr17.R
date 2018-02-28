brain_cpg = read.table('~/Documents/Paper3/Data/Brain_ES_chr17/brain_cpg')
brain_density = read.table('~/Documents/Paper3/Data/Brain_ES_chr17/brain_density')
es_bs = read.table('~/Documents/Paper3/Data/Brain_ES_chr17/es_methy')

colnames(brain_cpg) = c('chr', 'start', 'end', 'brain_cpg')
colnames(brain_density) = c('chr', 'start', 'end', 'brain_density')
colnames(es_bs) = c('chr', 'start', 'end', 'es_methy', 'es_total')

all_bs = merge(brain_cpg, brain_density, by=c('chr','start','end'))
all_bs = merge(all_bs, es_bs, by=c('chr','start','end'))
all_bs = all_bs[order(all_bs$start),]
all_bs$brain_methy = round(all_bs$brain_cpg*all_bs$brain_density)

write.table(all_bs[,c('chr','start','end','brain_methy','brain_density','es_methy','es_total')], 
            file = '/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_chr17/all_bs', quote = F, sep = "\t", col.names = F, row.names = F)


binlength = 500

dirwrite <- "/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_chr17/"

writefile <- paste(dirwrite, "brain_medip_bin.bed", sep = "")
countMeDIPbin('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_chr17/brain_medip1', writefile = writefile, binlength = binlength)

writefile <- paste(dirwrite, "brain_mre_bin.bed", sep = "")
countMREbin('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_chr17/brain_mre1', writefile = writefile, binlength = binlength)

writefile <- paste(dirwrite, "es_medip_bin.bed", sep = "")
countMeDIPbin('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_chr17/es_medip1', writefile = writefile, binlength = binlength)

writefile <- paste(dirwrite, "es_mre_bin.bed", sep = "")
countMREbin('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_chr17/es_mre1', writefile = writefile, binlength = binlength)

writefile = paste(dirwrite, "cpgbin.bed", sep = "")
countcpgbin('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_chr17/all_bs', writefile = writefile, binlength = binlength)

file <- '~/Documents/Paper3/Data/TriMRE_frags.bed'
file1 <- '~/Documents/Paper3/Data/Brain_ES_chr17/all_bs'
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
writefile <- paste(dirwrite, "pval_Brain_es_chr17.bed", sep = "")
reportfile <- paste(dirwrite, "report_Brain_es_chr17.txt", sep = "")
MnM.test(file.dataset = datafile, file.cpgbin = cpgfile,
         file.mrecpgbin = mrecpgfile, writefile = writefile, reportfile = reportfile,
         mreratio = 3/7, method = "XXYY", psd = 2, mkadded = 1, a = 1e-16,
         cut = 100, top = 500)


datafile <- paste(dirwrite, "pval_Brain_es_chr17.bed", sep = "")
writefile <- paste(dirwrite, "q_Brain_es_chr17.bed", sep = "")
reportfile <- paste(dirwrite, "report_q_Brain_es_chr17.bed", sep = "")
MnM.qvalue(datafile, writefile, reportfile)

qval = read.table('~/Documents/Paper3/Data/Brain_ES_chr17/q_Brain_es_chr17.bed', header = T)


idx = floor(all[,'start']/500)+1
all$mm_brain_medip = qval$Medip1[idx]
all$mm_es_medip = qval$Medip2[idx]
all$mm_brain_mre = qval$MRE1[idx]
all$mm_es_mre = qval$MRE2[idx]
all$mm_q = qval$qvalue[idx]

write.table(all, '~/Documents/Paper3/Data/Brain_ES_chr17/crf/all', quote = F, sep = " ", col.names = T, row.names = F)


all <- read.table('~/Documents/Paper3/Data/Brain_ES_chr17/all_bs4', header = F)
colnames(all) = colnames(all_18)[1:11]

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

all = all[,c(1:11,19,20)]
all[,'is_mre_cpg'] = 0

all[,'label'] = ifelse(all$fdr < 0.001, 1, 0)

all[,'dist'] = all$start - c(all$start[1], all$start[1:nrow(all)-1])










file <- '~/Documents/Paper3/Data/TriMRE_frags.bed'
five_Mre_CpGsite <- read.table(file, header = FALSE, as.is = TRUE)
five_Mre_CpGsite <- five_Mre_CpGsite[which(five_Mre_CpGsite[,1]=='chr17'),]
four_Mre_CpGsite <- five_Mre_CpGsite[five_Mre_CpGsite[, 4] != "ACGT", ]
mrecpg.site <- four_Mre_CpGsite[four_Mre_CpGsite[, 4] != "CGCG", ]



mre_cpg = read.table('~/Documents/Paper3/Data/Brain_ES_chr18/three_mre_cpg.bed', header = T)
colnames(mre_cpg) = c('chr','start','end','is_mre_cpg')
mre_cpg$chr = 18

all1 = merge(all, mre_cpg, by = c('chr','start','end'), all.x = TRUE)
idx = is.na(all1$is_mre_cpg)
all1[idx,'is_mre_cpg'] = 0

write.table(all1, '~/Documents/Paper3/Data/Brain_ES_chr18/crf/all1.bed', quote = F, sep = " ", col.names = T, row.names = F)
all = read.table('~/Documents/Paper3/Data/Brain_ES/crf/all_bs10.bed', header = T)



all_17 = read.table('~/Documents/Paper3/Data/Brain_ES_chr17/crf/all', header = T)
all_18 = read.table('~/Documents/Paper3/Data/Brain_ES_chr18/crf/all_mm_1', header = T)

all_ref = rbind(all_18, all_17)


summary(all_17$es_medip)







