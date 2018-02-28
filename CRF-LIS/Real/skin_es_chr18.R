skin_cpg = read.table('~/Documents/Paper3/Data/Skin_ES_chr18/skin_cpg')
skin_density = read.table('~/Documents/Paper3/Data/Skin_ES_chr18/skin_density')
es_bs = read.table('~/Documents/Paper3/Data/Skin_ES_chr18/es_bs.bed')

colnames(skin_cpg) = c('chr', 'start', 'end', 'skin_cpg')
colnames(skin_density) = c('chr', 'start', 'end', 'skin_density')
colnames(es_bs) = c('chr', 'start', 'end', 'es_methy', 'es_total')

all_bs = merge(skin_cpg, skin_density, by=c('chr','start','end'))
all_bs = merge(all_bs, es_bs, by=c('chr','start','end'))
all_bs = all_bs[order(all_bs$start),]
all_bs$skin_methy = round(all_bs$skin_cpg*all_bs$skin_density)

write.table(all_bs[,c('chr','start','end','skin_methy','skin_density','es_methy','es_total')], 
            file = '/Users/xiaoyudai/Documents/Paper3/Data/Skin_ES_chr18/all_bs', quote = F, sep = " ", col.names = F, row.names = F)



binlength = 2000

dirwrite <- "/Users/xiaoyudai/Documents/Paper3/Data/Skin_ES_chr18/bin2000/"

writefile <- paste(dirwrite, "skin_medip_bin.bed", sep = "")
countMeDIPbin('/Users/xiaoyudai/Documents/Paper3/Data/Skin_ES_chr18/skin_medip', writefile = writefile, binlength = binlength)

writefile <- paste(dirwrite, "skin_mre_bin.bed", sep = "")
countMREbin('/Users/xiaoyudai/Documents/Paper3/Data/Skin_ES_chr18/skin_mre', writefile = writefile, binlength = binlength)

writefile <- paste(dirwrite, "es_medip_bin.bed", sep = "")
countMeDIPbin('/Users/xiaoyudai/Documents/Paper3/Data/Skin_ES_chr18/es_medip.bed', writefile = writefile, binlength = binlength)

writefile <- paste(dirwrite, "es_mre_bin.bed", sep = "")
countMREbin('/Users/xiaoyudai/Documents/Paper3/Data/Skin_ES_chr18/es_mre.bed', writefile = writefile, binlength = binlength)

writefile = paste(dirwrite, "cpgbin.bed", sep = "")
countcpgbin('/Users/xiaoyudai/Documents/Paper3/Data/Skin_ES_chr18/all_bs', writefile = writefile, binlength = binlength)


file <- '~/Documents/Paper3/Data/TriMRE_frags.bed'
file1 <- '~/Documents/Paper3/Data/Skin_ES_chr18/all_bs'
allcpgfile <- paste(dirwrite, "cpgbin.bed", sep = "")
five_Mre_CpGsite <- read.table(file, header = FALSE, as.is = TRUE)
four_Mre_CpGsite <- five_Mre_CpGsite[five_Mre_CpGsite[, 4] != "ACGT", ]
mrecpg.site <- four_Mre_CpGsite[four_Mre_CpGsite[, 4] != "CGCG", ]
writefile <- paste(dirwrite, "three_mre_cpg_bin.bed", sep = "")
countMREcpgbin(mrecpg.site, file.allcpgsite = file1, file.bin = allcpgfile,
               writefile = writefile, binlength = binlength)


datafile1 <- paste(dirwrite, "skin_medip_bin1.bed", sep = "")
datafile2 <- paste(dirwrite, "es_medip_bin1.bed", sep = "")
datafile3 <- paste(dirwrite, "skin_mre_bin1.bed", sep = "")
datafile4 <- paste(dirwrite, "es_mre_bin1.bed", sep = "")
datafile <- c(datafile1, datafile2, datafile3, datafile4)
cpgfile <- paste(dirwrite, "cpgbin1.bed", sep = "")
mrecpgfile <- paste(dirwrite, "three_mre_cpg_bin1.bed", sep = "")
writefile <- paste(dirwrite, "pval_Skin_ES_chr18.bed", sep = "")
reportfile <- paste(dirwrite, "report_Skin_ES_chr18.txt", sep = "")
MnM.test(file.dataset = datafile, file.cpgbin = cpgfile,
         file.mrecpgbin = mrecpgfile, writefile = writefile, reportfile = reportfile,
         mreratio = 3/7, method = "XXYY", psd = 2, mkadded = 1, a = 1e-16,
         cut = 100, top = 500)

datafile <- paste(dirwrite, "pval_Skin_ES_chr18.bed", sep = "")
writefile <- paste(dirwrite, "q_Skin_ES_chr18.bed", sep = "")
reportfile <- paste(dirwrite, "report_q_Skin_ES_chr18.bed", sep = "")
MnM.qvalue(datafile, writefile, reportfile)

qval = read.table(writefile, header = T)

idx = floor(all_mm2[,'start']/binlength)+1
all_mm2[,'mm_q_2000'] = qval$qvalue[idx]

pred_mm = ifelse(all_mm2$mm_q_2000<0.001, 1, 0)
truth = all_mm2$label
confusion.matrix(truth, pred_mm) # 3386 7570


qval = read.table('~/Documents/Paper3/Data/Skin_ES_chr18/bin2000/q_Skin_ES_chr18.bed', header = T)

idx = floor(all_mm2[,'start']/2000)+1
all_mm2$Medip1_2000 = qval$Medip1[idx]
all_mm2$Medip2_2000 = qval$Medip2[idx]
all_mm2$MRE1_2000 = qval$MRE1[idx]
all_mm2$MRE2_2000 = qval$MRE2[idx]



all_mm2 = all_mm2[,c(1:27,59:64)]
colnames(all_mm2)






all_bs = read.table('~/Documents/Paper3/Data/Skin_ES_chr18/all_bs')
all_bs$V1 = 18
write.table(all_bs, '~/Documents/Paper3/Data/Skin_ES_chr18/all_bs1', quote = F, sep = "\t", col.names = F, row.names = F)

skin_mre = read.table('~/Documents/Paper3/Data/Skin_ES_chr18/skin_mre')
skin_mre$V3 = skin_mre$V2 + 2
write.table(skin_mre[,c(1,2,3)], '~/Documents/Paper3/Data/Skin_ES_chr18/skin_mre1', quote = F, sep = "\t", col.names = F, row.names = F)

es_mre = read.table('~/Documents/Paper3/Data/Skin_ES_chr18/es_mre.bed')
es_mre$V3 = es_mre$V2 + 2
write.table(es_mre[,c(1,2,3)], '~/Documents/Paper3/Data/Skin_ES_chr18/es_mre1', quote = F, sep = "\t", col.names = F, row.names = F)


all = read.table('~/Documents/Paper3/Data/Skin_ES_chr18/all_bs5')
colnames(all) = c('chr','start','end','skin_methy','skin_total','es_methy','es_total','skin_medip','skin_mre','es_medip','es_mre')
write.table(all, '~/Documents/Paper3/Data/Skin_ES_chr18/crf/all_raw', quote = F, sep = " ", col.names = T, row.names = F)

all1 = read.table('~/Documents/Paper3/Data/Skin_ES_chr18/crf/all_raw', header = T)
idx = which(all$skin_total>10 & all$es_total>10)
all = all[idx,]

dat1 <- all[,c('chr','start','skin_total','skin_methy')]
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

write.table(all1[,c(1:11,19,20,21)], '~/Documents/Paper3/Data/Skin_ES_chr18/crf/all_gf.bed', quote = F, sep = " ", col.names = T, row.names = F)
all = read.table('~/Documents/Paper3/Data/BraiSkin_ES_chr18n_ES/crf/all_gf.bed', header = T)





