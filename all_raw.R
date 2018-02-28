
all_bs = read.table('~/Documents/Paper3/Data/Brain_ES_chr18/all_bs.bed')
colnames(all_bs) = c('chr','start','end','brain_methy','brain_total','es_methy','es_total')
idx = floor(all_bs[,'start']/500)+1
all_bs$mm_brain_medip = qval$Medip1[idx]
all_bs$mm_es_medip = qval$Medip2[idx]
all_bs$mm_brain_mre = qval$MRE1[idx]
all_bs$mm_es_mre = qval$MRE2[idx]
all_bs$mm_q = qval$qvalue[idx]

write.table(all_bs, '~/Documents/Paper3/Data/Brain_ES_chr18/crf/all_raw', quote = F, sep = " ", col.names = T, row.names = F)

all_raw = read.table('~/Documents/Paper3/Data/Brain_ES_chr18/crf/all_raw', header = T)
head(all_raw)

pred_train = train_crf$label
pred_train[which(pred_train == ' ')] = 0
pred_train = as.numeric(pred_train)

auc(truth, 1 - all_mm$mm_q) # 0.5632
auc(truth, pred_train) # 0.5989
auc(truth, prob1) # 0.8468
auc(truth, prob1.mlp) # 0.7454
auc(truth[1:length(prob1.rnn)], prob1.rnn) # 0.6877

auc(truth[idx_test], prob1[idx_test]) # 0.8039
auc(truth[idx_test], prob1.mlp[idx_test]) # 0.



lsi.crf = prob0/prob1

q = 0.1

res.crf <- mt.hmm(lsi.crf, q)$de
res.mlp <- mt.hmm(lsi.mlp, q)$de
res.mnm = ifelse(all_mm$mm_q < q, 1, 0)

res.crf <- mt.hmm(lsi.crf, 0.01)$de
res.mlp <- mt.hmm(lsi.mlp, 5e-8)$de
res.mnm = ifelse(all_mm$mm_q < 0.7, 1, 0)
c(sum(res.crf), sum(res.mlp), sum(res.mnm))

res.crf <- mt.hmm(lsi.crf, 0.001)$de
res.mlp <- mt.hmm(lsi.mlp, 1e-12)$de
res.mnm = ifelse(all_mm$mm_q < 0.012, 1, 0)
c(sum(res.crf), sum(res.mlp), sum(res.mnm))

confusion.matrix(truth, res.mnm)
confusion.matrix(truth, res.crf)
confusion.matrix(truth, res.mlp)

dif <- abs(all_mm$brain_methy/all_mm$brain_total - all_mm$es_methy/all_mm$es_total)

DF <- rbind(data.frame(method='CRF', obs=dif[which(res.crf==1)]),
            data.frame(method='MnM', obs=dif[which(res.mnm==1)]),
            data.frame(method='DNN', obs=dif[which(res.mlp==1)]))
DF$method <- as.factor(DF$method)
ggplot(DF, aes(x=obs, fill=method)) +
  geom_histogram(binwidth=0.1, colour="black", position="dodge") +
  xlab('methylation ratio difference')

df <- data.frame(
  value = c(dif[which(res.crf==1)],dif[which(res.mnm==1)],dif[which(res.mlp==1)]),
  group = c(rep('CRF',sum(res.crf)), rep('MnM',sum(res.mnm)),
            rep('DNN',sum(res.mlp)))
)
ggplot(df, aes(factor(group), value, fill=group)) + 
  geom_boxplot() + 
  coord_cartesian(ylim = c(0,1)) +
  xlab('methylation ratio difference')

 


idx_ovlp = which(res.mnm*res.crf == 1)

length(idx_ovlp)

sum(truth[idx_ovlp] == 1)
mean(truth[idx_ovlp] == 1)










confusion.matrix(truth[idx_test], res.mnm[idx_test])
confusion.matrix(truth[idx_test], res.crf[idx_test])






