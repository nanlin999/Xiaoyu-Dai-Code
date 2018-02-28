refGene <- read.table('~/Documents/Paper3/Data/genFeature/refGene/refGene.txt', header = F)
idx = which(refGene$V3 == 'chr17')
refGene = refGene[idx,c(3,5,6)]
refGene = refGene[order(refGene$V5),]
refGene_gr = with(refGene, GRanges(V3, IRanges(start = V5, end = V6)))

rmsk <- read.table('~/Documents/Paper3/Data/genFeature/rmsk/rmsk.txt', header = F)
idx = which(rmsk$V6 == 'chr17')
rmsk = rmsk[idx,]
rmsk_gr = with(rmsk, GRanges(V6, IRanges(start = V7, end = V8)))

cpgi <- read.table('~/Documents/Paper3/Data/genFeature/cpgi/cpgIslandExt.txt', header = F)
idx = which(cpgi$V2 == 'chr17')
cpgi = cpgi[idx,c(2,3,4)]
cpgi = cpgi[order(cpgi$V3),]
cpgi_gr = with(cpgi, GRanges(V2, IRanges(start = V3, end = V4)))


intron <- read.table('~/Documents/Paper3/Data/genFeature/intron/intron', header = F)
idx = which(intron$V1 == 'chr17')
intron = intron[idx,c(1,2,3)]
intron = intron[order(intron$V2),]
intron_gr = with(intron, GRanges(V1, IRanges(start = V2, end = V3)))

exon <- read.table('~/Documents/Paper3/Data/genFeature/exon/exon', header = F)
idx = which(exon$V1 == 'chr17')
exon = exon[idx,c(1,2,3)]
exon = exon[order(exon$V2),]
exon_gr = with(exon, GRanges(V1, IRanges(start = V2, end = V3)))

utr3 <- read.table('~/Documents/Paper3/Data/genFeature/utr3/utr3', header = F)
idx = which(utr3$V1 == 'chr17')
utr3 = utr3[idx,c(1,2,3)]
utr3 = utr3[order(utr3$V2),]
utr3_gr = with(utr3, GRanges(V1, IRanges(start = V2, end = V3)))

utr5 <- read.table('~/Documents/Paper3/Data/genFeature/utr5/utr5', header = F)
idx = which(utr5$V1 == 'chr17')
utr5 = utr5[idx,c(1,2,3)]
utr5 = utr5[order(utr5$V2),]
utr5_gr = with(utr5, GRanges(V1, IRanges(start = V2, end = V3)))



rmsk_LINE = rmsk[which(rmsk$V12 == 'LINE'), ]
rmsk_LowSimple = rmsk[which(rmsk$V12 == 'Low_complexity' |  rmsk$V12 == 'Simple_repeat'), ]
# rmsk_LTR = rmsk[which(rmsk$V12 == 'LTR'), ]
rmsk_SINE = rmsk[which(rmsk$V12 == 'SINE'), ]
rmsk_DNA = rmsk[which(rmsk$V12 == 'DNA'), ]
rmsk_other = rmsk[which(rmsk$V12 == 'Other'), ]


rmsk_LINE_gr = with(rmsk_LINE, GRanges(V6, IRanges(start = V7, end = V8)))
rmsk_LowSimple_gr = with(rmsk_LowSimple, GRanges(V6, IRanges(start = V7, end = V8)))
# rmsk_LTR_gr = with(rmsk_LTR, GRanges(V6, IRanges(start = V7, end = V8)))
rmsk_SINE_gr = with(rmsk_SINE, GRanges(V6, IRanges(start = V7, end = V8)))
rmsk_DNA_gr = with(rmsk_DNA, GRanges(V6, IRanges(start = V7, end = V8)))
rmsk_other_gr = with(rmsk_other, GRanges(V6, IRanges(start = V7, end = V8)))



name = 'refGene'
between = as.data.frame(intersect(all_gr, refGene_gr)) # change this as well!
between = between[,c(1,2,3)]
between[,name] = 1
colnames(between)[1] = 'chr'

all_gf = merge(all_gf, between, by = c('chr', 'start', 'end'), all.x=T)
idx = which(is.na(all_gf[,name]))
all_gf[idx,name] = 0
head(all_gf)
table(all_gf[,name])



name = 'notAnno'
idx = which(rowSums(all_gf[,c(15:26)]) == 0)
vec_notAnno = rep(0,nrow(all_gf))
vec_notAnno[idx] = 1
all_gf[,name] = vec_notAnno
head(all_gf)
table(all_gf[,name])

write.table(all_gf, '~/Documents/Paper3/Data/Skin_ES_chr18/crf/all_gf.bed', quote = F, sep = " ", col.names = T, row.names = F)

all_gf = read.table('~/Documents/Paper3/Data/Brain_ES_chr18/crf/all_gf.bed', header = T)
all_gr = with(all_gf, GRanges(chr, IRanges(start = start, end = end)))

all_sm = all_new[which(all_new$is.rsmk == 1), c(1:107)]

cnames = colnames(all_gf)

cpg.bed = read.table('~/Documents/Paper3/software/methylCRF/mscse-methylcrf-fd50dc97b171/chr18/cpg.bed')
cpg.bed = cpg.bed[,c(1:4)]
colnames(cpg.bed) = c('chr','start','end','id')

core_promoter = read.table('~/Documents/Paper3/software/methylCRF/mscse-methylcrf-fd50dc97b171/chr18/core_promoter_cpg.bin')
colnames(core_promoter) = c('id','core_prometor')
  
distal_promoter = read.table('~/Documents/Paper3/software/methylCRF/mscse-methylcrf-fd50dc97b171/chr18/distal_promoter_cpg.bin')
colnames(distal_promoter) = c('id','distal_promoter')

proximal_promoter = read.table('~/Documents/Paper3/software/methylCRF/mscse-methylcrf-fd50dc97b171/chr18/proximal_promoter_cpg.bin')
colnames(proximal_promoter) = c('id','proximal_promoter')

cpgi = read.table('~/Documents/Paper3/software/methylCRF/mscse-methylcrf-fd50dc97b171/chr18/cpgi_cpg.bin')
colnames(cpgi) = c('id','cpgi')

cpgi_1k = read.table('~/Documents/Paper3/software/methylCRF/mscse-methylcrf-fd50dc97b171/chr18/cpgi_1kshore_cpg.bin')
colnames(cpgi_1k) = c('id','cpgi_1k')

cpgi_2k = read.table('~/Documents/Paper3/software/methylCRF/mscse-methylcrf-fd50dc97b171/chr18/cpgi_2kshore_cpg.bin')
colnames(cpgi_2k) = c('id','cpgi_2k')

gene_body = read.table('~/Documents/Paper3/software/methylCRF/mscse-methylcrf-fd50dc97b171/chr18/gene_body_cpg.bin')
colnames(gene_body) = c('id','gene_body')

rmsk_DNA = read.table('~/Documents/Paper3/software/methylCRF/mscse-methylcrf-fd50dc97b171/chr18/rmsk_DNA')
colnames(rmsk_DNA) = c('id','rmsk_DNA')

rmsk_RNA = read.table('~/Documents/Paper3/software/methylCRF/mscse-methylcrf-fd50dc97b171/chr18/rmsk_RNA')
colnames(rmsk_RNA) = c('id','rmsk_RNA')

rmsk_LINE = read.table('~/Documents/Paper3/software/methylCRF/mscse-methylcrf-fd50dc97b171/chr18/rmsk_LINE')
colnames(rmsk_LINE) = c('id','rmsk_LINE')

rmsk_LTR = read.table('~/Documents/Paper3/software/methylCRF/mscse-methylcrf-fd50dc97b171/chr18/rmsk_LTR')
colnames(rmsk_LTR) = c('id','rmsk_LTR')

rmsk_SINE = read.table('~/Documents/Paper3/software/methylCRF/mscse-methylcrf-fd50dc97b171/chr18/rmsk_SINE')
colnames(rmsk_SINE) = c('id','rmsk_SINE')

rmsk_lowsimple = read.table('~/Documents/Paper3/software/methylCRF/mscse-methylcrf-fd50dc97b171/chr18/rmsk_lowsimple')
colnames(rmsk_lowsimple) = c('id','rmsk_lowsimple')

rmsk_other = read.table('~/Documents/Paper3/software/methylCRF/mscse-methylcrf-fd50dc97b171/chr18/rmsk_other')
colnames(rmsk_other) = c('id','rmsk_other')

intron = read.table('~/Documents/Paper3/software/methylCRF/mscse-methylcrf-fd50dc97b171/chr18/intron')
colnames(intron) = c('id','intron')

exon = read.table('~/Documents/Paper3/software/methylCRF/mscse-methylcrf-fd50dc97b171/chr18/exon')
colnames(exon) = c('id','exon')

utr3 = read.table('~/Documents/Paper3/software/methylCRF/mscse-methylcrf-fd50dc97b171/chr18/utr3')
colnames(utr3) = c('id','utr3')

utr5 = read.table('~/Documents/Paper3/software/methylCRF/mscse-methylcrf-fd50dc97b171/chr18/utr5')
colnames(utr5) = c('id','utr5')

notAnno = read.table('~/Documents/Paper3/software/methylCRF/mscse-methylcrf-fd50dc97b171/chr18/notAnno')
colnames(notAnno) = c('id','notAnno')




gc = read.table('~/Documents/Paper3/software/methylCRF/mscse-methylcrf-fd50dc97b171/chr18/gdata.tbl')
colnames(gc) = c('id','gc10','gc75','gc250','cpgprev','cpgratio','phastCons46way')


all_ref = merge(cpg.bed, gc, by=c('id'))
all_ref = all_ref[,c(2:7,9:10)]
all_gf = merge(all_gf, all_ref, by=c('chr','start','end'))
all_gf = all_gf[order(all_gf$start),]

head(all_gf)

# all_ref$chr = 18

write.table(all_gf, '~/Documents/Paper3/Data/Brain_ES_chr18/crf/all_gf.bed', quote = F, sep = " ", col.names = T, row.names = F)

all_ref = all_ref[,c(2:4,12:17)]
all_new = merge(all_gf, all_ref, by=c('chr','start','end'))
all_gf1 = all_new[order(all_new$start),]

all_gf1 = all_new[,c(1:23,25,26,60:66,28:58)]
colnames(all_gf1)[29] = 'cpgi'
write.table(all_gf, '~/Documents/Paper3/Data/Brain_ES_chr18/crf/all_gf.bed', quote = F, sep = " ", col.names = T, row.names = F)


head(all_gf)
x = c(0, all_gf$start[2:nrow(all_gf)] - all_gf$end[1:(nrow(all_gf)-1)])
all_gf$dist = x
all_gf1 = all_gf[,c(1:66,68,69)]


all_gf$skin_medip = all_gf$skin_medip/quantile(all_gf[which(all_gf$skin_medip>0),'skin_medip'], 0.75) * 10
all_gf$skin_mre = all_gf$skin_mre/quantile(all_gf[which(all_gf$skin_mre>0),'skin_mre'], 0.75) * 10
all_gf$es_medip = all_gf$es_medip/quantile(all_gf[which(all_gf$es_medip>0),'es_medip'], 0.75) * 10
all_gf$es_mre = all_gf$es_mre/quantile(all_gf[which(all_gf$es_mre>0),'es_mre'], 0.75) * 10

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


idx = floor(all_gf[,'start']/500)+1
all_gf$mm_brain_medip = qval$Medip1[idx]
all_gf$mm_es_medip = qval$Medip2[idx]
all_gf$mm_brain_mre = qval$MRE1[idx]
all_gf$mm_es_mre = qval$MRE2[idx]
all_gf$mm_q = qval$qvalue[idx]

all_mm = all_gf[,c(1:11,19:70)]
write.table(all_mm, '~/Documents/Paper3/Data/Skin_ES_chr18/crf/all_mm', quote = F, sep = " ", col.names = T, row.names = F)










