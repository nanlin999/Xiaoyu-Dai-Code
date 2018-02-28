library("GenomicRanges")
library(ggplot2)

# write.table(dmlTest.sm, '~/Documents/Paper3/Data/chr18/crf/dmlTestsm', quote = F, sep = " ", col.names = T, row.names = F)
# wgbs_all <- all_sm[round(0.2*nrow(all)):nrow(all), c('chr','start','end','label')]
# wgbs_all$label <- (wgbs_all$fdr < label.threshold) * 1
# wgbs_all$label <- test_crf$label
gr_all <- with(all_mm2, GRanges(chr, IRanges(start = start, end = end), label=label))


range_crf = dmr_crf[,1:3]
range_mnm = dmr_mnm[,1:3]
range_mlp = dmr_mlp[,1:3]
range_rnn = dmr_rnn[,1:3]

gr_crf = with(range_crf, GRanges(chr, IRanges(start = start, end = end)))
gr_mnm = with(range_mnm, GRanges(chr, IRanges(start = chrSt, end = chrEnd)))
gr_mlp = with(range_mlp, GRanges(chr, IRanges(start = start, end = end)))
gr_rnn = with(range_rnn, GRanges(chr, IRanges(start = start, end = end)))

# dif_crf = as.data.frame(setdiff(gr_crf, gr_mnm))
# dif_mnm = as.data.frame(setdiff(gr_mnm, gr_crf))
gr_crf = as.data.frame(gr_crf)
gr_mnm = as.data.frame(gr_mnm)
gr_mlp = as.data.frame(gr_mlp)
gr_rnn = as.data.frame(gr_rnn)


crf_methy_dif = c()
for (i in 1:nrow(gr_crf)) {
  start = gr_crf[i,'start']
  end = gr_crf[i,'end']
  idx = which(all_raw$start >= start & all_raw$end <= end)
  if (length(idx) > 0) {
    methy_dif = abs( sum(all_raw[idx,'skin_methy'])/sum(all_raw[idx,'skin_total']) -
                       sum(all_raw[idx,'es_methy'])/sum(all_raw[idx,'es_total']) )
  } else {
    methy_dif = 0
  }
  crf_methy_dif = c(crf_methy_dif, methy_dif)
}

mm_methy_dif = c()
for (i in 1:nrow(gr_mnm)) {
  start = gr_mnm[i,'start']
  end = gr_mnm[i,'end']
  idx = which(all_raw$start >= start & all_raw$end <= end)
  if (length(idx) > 0) {
    methy_dif = abs( sum(all_raw[idx,'skin_methy'])/sum(all_raw[idx,'skin_total']) -
                       sum(all_raw[idx,'es_methy'])/sum(all_raw[idx,'es_total']) )
  } else {
    methy_dif = 0
  }
  mm_methy_dif = c(mm_methy_dif, methy_dif)
}

mlp_methy_dif = c()
for (i in 1:nrow(gr_mlp)) {
  start = gr_mlp[i,'start']
  end = gr_mlp[i,'end']
  idx = which(all_raw$start >= start & all_raw$end <= end)
  if (length(idx) > 0) {
    methy_dif = abs( sum(all_raw[idx,'brain_methy'])/sum(all_raw[idx,'brain_total']) -
                       sum(all_raw[idx,'es_methy'])/sum(all_raw[idx,'es_total']) )
  } else {
    methy_dif = 0
  }
  mlp_methy_dif = c(mlp_methy_dif, methy_dif)
}


crf_labelratio = c()
crf_cpgnum = c()
for (i in 1:nrow(gr_crf)) {
  start = gr_crf[i,'start']
  end = gr_crf[i,'end']
  idx = which(all_mm2$start >= start & all_mm2$end <= end)
  if (length(idx) > 0) {
    labelratio = mean(all_mm2[idx,'label'])
  } else {
    labelratio = 0
  }
  crf_labelratio = c(crf_labelratio, labelratio)
  crf_cpgnum = c(crf_cpgnum, length(idx))
}

mnm_labelratio = c()
mnm_cpgnum = c()
for (i in 1:nrow(gr_mnm)) {
  start = gr_mnm[i,'start']
  end = gr_mnm[i,'end']
  idx = which(all_mm2$start >= start & all_mm2$end <= end)
  if (length(idx) > 0) {
    labelratio = mean(all_mm2[idx,'label'])
  } else {
    labelratio = 0
  }
  mnm_labelratio = c(mnm_labelratio, labelratio)
  mnm_cpgnum = c(mnm_cpgnum, length(idx))
}

mlp_labelratio = c()
mlp_cpgnum = c()
for (i in 1:nrow(gr_mlp)) {
  start = gr_mlp[i,'start']
  end = gr_mlp[i,'end']
  idx = which(all_mm$start >= start & all_mm$end <= end)
  if (length(idx) > 0) {
    labelratio = mean(all_mm[idx,'label'])
  } else {
    labelratio = 0
  }
  mlp_labelratio = c(mlp_labelratio, labelratio)
  mlp_cpgnum = c(mlp_cpgnum, length(idx))
}

rnn_labelratio = c()
rnn_cpgnum = c()
for (i in 1:nrow(gr_rnn)) {
  start = gr_rnn[i,'start']
  end = gr_rnn[i,'end']
  idx = which(wgbs_all$start >= start & wgbs_all$end <= end)
  if (length(idx) > 0) {
    labelratio = mean(wgbs_all[idx,'label'])
  } else {
    labelratio = 0
  }
  rnn_labelratio = c(rnn_labelratio, labelratio)
  rnn_cpgnum = c(rnn_cpgnum, length(idx))
}

DF <- rbind(data.frame(method='CRF', obs=crf_methy_dif),
            data.frame(method='MnM', obs=mm_methy_dif)),
            data.frame(method='MLP', obs=mlp_methy_dif))
DF$method <- as.factor(DF$method)
ggplot(DF, aes(x=obs, fill=method)) +
  geom_histogram(binwidth=0.1, colour="black", position="dodge")


DF <- rbind(data.frame(method='CRF', obs=crf_labelratio),
            data.frame(method='MnM', obs=mnm_labelratio)),
            data.frame(method='MLP', obs=mlp_labelratio))
DF$method <- as.factor(DF$method)
ggplot(DF, aes(x=obs, fill=method)) +
  geom_histogram(binwidth=0.1, colour="black", position="dodge") + 
  xlab('true DMC ratio')


df <- data.frame(
  value = c(crf_methy_dif,mm_methy_dif,mlp_methy_dif),
  group = c(rep('CRF',length(crf_methy_dif)), rep('MnM',length(mm_methy_dif)),
            rep('DNN',length(mlp_methy_dif)))
)
ggplot(df, aes(factor(group), value, fill=group)) + 
  geom_boxplot() + 
  coord_cartesian(ylim = c(0,1)) +
  xlab('methylation ratio difference (DMR)')


df <- data.frame(
  value = c(crf_labelratio,mnm_labelratio,mlp_labelratio),
  group = c(rep('CRF',length(crf_labelratio)), rep('MnM',length(mnm_labelratio)),
            rep('DNN',length(mlp_labelratio)))
)
ggplot(df, aes(factor(group), value, fill=group)) + 
  geom_boxplot() + 
  coord_cartesian(ylim = c(0,1)) +
  xlab('label ratio (DMR)')

df <- data.frame(
  value = c(dmr_crf$length,dmr_mlp$length),
  group = c(rep('CRF',nrow(dmr_crf)), rep('DNN',nrow(dmr_mlp)))
)
ggplot(df, aes(factor(group), value, fill=group)) + 
  geom_boxplot() + 
  coord_cartesian(ylim = c(0,1000)) +
  geom_hline(yintercept = 500, linetype=2) +
  xlab('DMR length')


hist(dmr_crf$length, breaks=100)
hist(dmr_mlp$length, breaks=100)



hist(crf_cpgnum, breaks=100)
hist(mnm_cpgnum, breaks=100)
hist(mlp_cpgnum, breaks=100)

c(sum(crf_cpgnum), sum(crf_cpgnum*crf_labelratio))
sum(crf_cpgnum*crf_labelratio)/sum(crf_cpgnum)

c(sum(mnm_cpgnum), sum(mnm_cpgnum*mnm_labelratio))
sum(mnm_cpgnum*mnm_labelratio)/sum(mnm_cpgnum)

c(sum(mlp_cpgnum), sum(mlp_cpgnum*mlp_labelratio))
sum(mlp_cpgnum*mlp_labelratio)/sum(mlp_cpgnum)

c(sum(rnn_cpgnum), sum(rnn_cpgnum*rnn_labelratio))
sum(rnn_cpgnum*rnn_labelratio)/sum(rnn_cpgnum)

sum(gr_crf$width)
sum(gr_mnm$width)
sum(gr_mlp$width)























for (name in all_feature_names) {
  count = table(all$label, all[,name])
  barplot(count, xlab=name)
  print(name)
  print(count)
}











dif_mnm_labelratio = c()
dif_mnm_cpgnum = c()
for (i in 1:nrow(x)) {
  start = x[i,'start']
  end = x[i,'end']
  idx = which(wgbs_all$start >= start & wgbs_all$end <= end)
  labelratio = mean(wgbs_all[idx,'label'])
  dif_mnm_labelratio = c(dif_mnm_labelratio, labelratio)
  dif_mnm_cpgnum = c(dif_mnm_cpgnum, length(idx))
}

sum(dif_mnm_cpgnum)



dif_crf_pval = c()
for (i in 1:nrow(dif_crf)) {
  start = dif_crf[i,'start']
  end = dif_crf[i,'end']
  idx = which(wgbs_all$start >= start & wgbs_all$end <= end)
  pval = mean(wgbs_all[idx,'pval'])
  dif_crf_pval = c(dif_crf_pval, pval)
}

dif_mnm_pval = c()
for (i in 1:nrow(dif_mnm)) {
  start = dif_mnm[i,'start']
  end = dif_mnm[i,'end']
  idx = which(wgbs_all$start >= start & wgbs_all$end <= end)
  pval = mean(wgbs_all[idx,'pval'])
  dif_mnm_pval = c(dif_mnm_pval, pval)
}

boxplot(dif_mnm_pval)
boxplot(dif_crf_pval)




x = c(1,2,3,4,5,6,7)
y = c(7)

kk = kmeans(x, centers = 3)

fitted(kk)





mnm_labelratio = c()
n_cpg = c()
for (i in 1:nrow(dmr_mnm)) {
  start = dmr_mnm[i,'chrSt']
  end = dmr_mnm[i,'chrEnd']
  idx = which(wgbs_all$start >= start & wgbs_all$end <= end)
  labelratio = mean(wgbs_all[idx,'label'])
  n_cpg = c(n_cpg, length(idx))
  if (!is.na(labelratio)) {
    mnm_labelratio = c(mnm_labelratio, labelratio)
    if(labelratio == 1){
      print(length(idx))
    }
  }
}






mre = read.table('~/Documents/Paper3/Data/chr18/MRE_1.bed')
colnames(mre) = c('chr', 'start', 'end')
mre$chr = 'chr18'
gr_mre = with(mre, GRanges(chr, IRanges(start = start, end = end)))

cpg = read.table('~/Documents/Paper3/Data/chr18/density_1.bed')
colnames(cpg) = c('chr', 'start', 'end')
cpg$chr = 'chr18'
gr_cpg = with(cpg, GRanges(chr, IRanges(start = start, end = end)))


dif = intersect(gr_mre, gr_cpg)
dim(dif)
