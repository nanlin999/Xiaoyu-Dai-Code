skin_mre = read.table('~/Documents/Paper3/Data/Skin/Skin03_Keratinocyte_MRE_A13914.bed', sep = '\t', header = F)

new = skin_mre[,c(1,2,3)]
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

write.table(new, '~/Documents/Paper3/Data/Skin_ES_all/skin_mre', quote = F, sep = "\t", col.names = F, row.names = F)

all2 = read.table('~/Documents/Paper3/Data/Skin_ES_all/all2', sep = '\t', header = F)
all3 = all2[,c(1,2,3,4,8,12,13)]
all3$V4 = round(all3$V4 * all3$V8)
write.table(all3, '~/Documents/Paper3/Data/Skin_ES_all/all3', quote = F, sep = "\t", col.names = F, row.names = F)


all_raw = read.table('~/Documents/Paper3/Data/Skin_ES_all/all7', sep = '\t', header = F)
colnames(all_raw) = c('chr','start','end','skin_methy','skin_total','es_methy','es_total','skin_medip','skin_mre','es_medip','es_mre')
write.table(all_raw, '~/Documents/Paper3/Data/Skin_ES_all/all_raw', quote = F, sep = " ", col.names = T, row.names = F)

id = which(all_raw$skin_total > 10 & all_raw$es_total > 10)
all = all_raw[id,]

dat1 <- all[,c('chr','start','skin_total','skin_methy')]
colnames(dat1) <- c('chr','pos','N','X')

dat2 <- all[,c('chr','start','es_total','es_methy')]
colnames(dat2) <- c('chr','pos','N','X')

BSobj <- makeBSseqData( list(dat1, dat2), c("C1", "N1") )
dmlTest.sm <- DMLtest(BSobj, group1=c("C1"), group2=c("N1"), smoothing=TRUE)
colnames(dmlTest.sm)[2] = 'start'

all1 = merge(all, dmlTest.sm, by = c('chr','start'))

write.table(all1, file = '/Users/xiaoyudai/Documents/Paper3/Data/Skin_ES_all/all_gf.bed', quote = F, sep = " ", col.names = T, row.names = F)

all_gf[,'label'] = ifelse(all_gf$fdr < 0.001, 1, 0)

all_gf = all_gf[,c(1:11,19:21)]

##  add genomic feature..... ###


all_gf = read.table('~/Documents/Paper3/Data/Skin_ES_all/crf/all_gf1.bed', header = T)


all_gf$skin_medip = all_gf$skin_medip/quantile(all_gf[which(all_gf$skin_medip>0),'skin_medip'], 0.75) * 10
all_gf$skin_mre = all_gf$skin_mre/quantile(all_gf[which(all_gf$skin_mre>0),'skin_mre'], 0.75) * 10
all_gf$es_medip = all_gf$es_medip/quantile(all_gf[which(all_gf$es_medip>0),'es_medip'], 0.75) * 10
all_gf$es_mre = all_gf$es_mre/quantile(all_gf[which(all_gf$es_mre>0),'es_mre'], 0.75) * 10

all_gf[,'dist'] = all_gf$start - c(all_gf$start[1], all_gf$start[1:nrow(all_gf)-1])

all_gf$skin_medip_10 = ma(all_gf$skin_medip, 10)
all_gf$skin_medip_50 = ma(all_gf$skin_medip, 50)
all_gf$skin_medip_200 = ma(all_gf$skin_medip, 200)
all_gf$skin_medip_1k = ma(all_gf$skin_medip, 1000)

all_gf$es_medip_10 = ma(all_gf$es_medip, 10)
all_gf$es_medip_50 = ma(all_gf$es_medip, 50)
all_gf$es_medip_200 = ma(all_gf$es_medip, 200)
all_gf$es_medip_1k = ma(all_gf$es_medip, 1000)

all_gf$skin_mre_10 = ma(all_gf$skin_mre, 10)
all_gf$skin_mre_50 = ma(all_gf$skin_mre, 50)
all_gf$skin_mre_200 = ma(all_gf$skin_mre, 200)
all_gf$skin_mre_1k = ma(all_gf$skin_mre, 1000)

all_gf$es_mre_10 = ma(all_gf$es_mre, 10)
all_gf$es_mre_50 = ma(all_gf$es_mre, 50)
all_gf$es_mre_200 = ma(all_gf$es_mre, 200)
all_gf$es_mre_1k = ma(all_gf$es_mre, 1000)

all_gf$ratio_medip = abs(all_gf$skin_medip - all_gf$es_medip) / min(all_gf$skin_medip+1, all_gf$es_medip+1)
all_gf$ratio_medip_10 = abs(all_gf$skin_medip_10 - all_gf$es_medip_10) / min(all_gf$skin_medip_10+1, all_gf$es_medip_10+1)
all_gf$ratio_medip_50 = abs(all_gf$skin_medip_50 - all_gf$es_medip_50) / min(all_gf$skin_medip_50+1, all_gf$es_medip_50+1)
all_gf$ratio_medip_200 = abs(all_gf$skin_medip_200 - all_gf$es_medip_200) / min(all_gf$skin_medip_200+1, all_gf$es_medip_200+1)
all_gf$ratio_medip_1k = abs(all_gf$skin_medip_1k - all_gf$es_medip_1k) / min(all_gf$skin_medip_1k+1, all_gf$es_medip_1k+1)

all_gf$ratio_mre = abs(all_gf$skin_mre - all_gf$es_mre) / min(all_gf$skin_mre+1, all_gf$es_mre+1)
all_gf$ratio_mre_10 = abs(all_gf$skin_mre_10 - all_gf$es_mre_10) / min(all_gf$skin_mre_10+1, all_gf$es_mre_10+1)
all_gf$ratio_mre_50 = abs(all_gf$skin_mre_50 - all_gf$es_mre_50) / min(all_gf$skin_mre_50+1, all_gf$es_mre_50+1)
all_gf$ratio_mre_200 = abs(all_gf$skin_mre_200 - all_gf$es_mre_200) / min(all_gf$skin_mre_200+1, all_gf$es_mre_200+1)
all_gf$ratio_mre_1k = abs(all_gf$skin_mre_1k - all_gf$es_mre_1k) / min(all_gf$skin_mre_1k+1, all_gf$es_mre_1k+1)

all_gf$dif = abs((all_gf$skin_mre+1)*(all_gf$es_medip+1) - (all_gf$es_mre+1)*(all_gf$skin_medip+1))
all_gf$dif_10 = abs((all_gf$skin_mre_10+1)*(all_gf$es_medip_10+1) - (all_gf$es_mre_10+1)*(all_gf$skin_medip_10+1))
all_gf$dif_50 = abs((all_gf$skin_mre_50+1)*(all_gf$es_medip_50+1) - (all_gf$es_mre_50+1)*(all_gf$skin_medip_50+1))
all_gf$dif_200 = abs((all_gf$skin_mre_200+1)*(all_gf$es_medip_200+1) - (all_gf$es_mre_200+1)*(all_gf$skin_medip_200+1))
all_gf$dif_1k = abs((all_gf$skin_mre_1k+1)*(all_gf$es_medip_1k+1) - (all_gf$es_mre_1k+1)*(all_gf$skin_medip_1k+1))


write.table(all_gf, file = '/Users/xiaoyudai/Documents/Paper3/Data/Skin_ES_all/all_gf.bed', quote = F, sep = " ", col.names = T, row.names = F)



























