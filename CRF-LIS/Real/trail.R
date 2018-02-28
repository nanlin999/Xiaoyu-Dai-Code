range_crf = dmr_crf[,1:3]
range_mnm = dmr_mnm[,1:3]
range_crf$chr = 'chr18'
gr_crf = with(range_crf, GRanges(chr, IRanges(start = start, end = end)))
gr_mnm = with(range_mnm, GRanges(chr, IRanges(start = chrSt, end = chrEnd)))
gr_crf = as.data.frame(gr_crf)
gr_mnm = as.data.frame(gr_mnm)
crf_labelratio = c()
crf_cpgnum = c()
for (i in 1:nrow(gr_crf)) {
  start = gr_crf[i,'start']
  end = gr_crf[i,'end']
  idx = which(wgbs_all$start >= start & wgbs_all$end <= end)
  if (length(idx) > 0) {
    labelratio = mean(wgbs_all[idx,'label'])
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
  idx = which(wgbs_all$start >= start & wgbs_all$end <= end)
  if (length(idx) > 0) {
    labelratio = mean(wgbs_all[idx,'label'])
  } else {
    labelratio = 0
  }
  mnm_labelratio = c(mnm_labelratio, labelratio)
  mnm_cpgnum = c(mnm_cpgnum, length(idx))
}

ovlp_labelratio = c()
ovlp_cpgnum = c()
for (i in 1:nrow(ovlp)) {
  start = ovlp[i,'start']
  end = ovlp[i,'end']
  idx = which(wgbs_all$start >= start & wgbs_all$end <= end)
  if (length(idx) > 0) {
    labelratio = mean(wgbs_all[idx,'label'])
  } else {
    labelratio = 0
  }
  ovlp_labelratio = c(ovlp_labelratio, labelratio)
  ovlp_cpgnum = c(ovlp_cpgnum, length(idx))
}


DF <- rbind(data.frame(method='CRF', obs=crf_labelratio),
            data.frame(method='MnM', obs=mnm_labelratio),
            data.frame(method='ovlp_mnm', obs=ovlp_labelratio))
DF$method <- as.factor(DF$method)
ggplot(DF, aes(x=obs, fill=method)) +
  geom_histogram(binwidth=0.1, colour="black", position="dodge")

c(sum(crf_cpgnum), sum(crf_cpgnum*crf_labelratio))
sum(crf_cpgnum*crf_labelratio)/sum(crf_cpgnum)

c(sum(mnm_cpgnum), sum(mnm_cpgnum*mnm_labelratio))
sum(mnm_cpgnum*mnm_labelratio)/sum(mnm_cpgnum)



hist(crf_cpgnum, breaks=0:155, xlim=c(0,40))
hist(mnm_cpgnum, breaks=0:40, add=T, col='grey', alpha=0.1)



id_crf = which(res.crf == 1)
ref = abs(all_gf$brain_methy/all_gf$brain_total - all_gf$es_methy/all_gf$es_total)


boxplot(ref[id_crf])
boxplot(mg$ref, add=T)
boxplot(mg1$ref)



df <- data.frame(
  value = c(ref[id_crf],mg$ref,mg1$ref),
  group = c(rep('crf_only',length(id_crf)), rep('both',length(mg$ref)),
            rep('mnm_only',length(mg1$ref)))
)

df

library(ggplot2)
ggplot(df, aes(factor(group), value)) + geom_boxplot() + coord_cartesian(ylim = c(0,1))
                                                                         


id = which((all_gf$brain_medip-all_gf$es_medip)>50)
id = which(all_gf$ratio_medip>14)

hist(all_gf$brain_medip[which(all_gf$brain_medip<50)],breaks = 100)

all_gf[id,c(4:12,65)]

x = c()
y = c()
for(value in sort(unique(all_gf$brain_medip))) {
  id = which(all_gf$brain_medip == value)
  x = c(x, value)
  y = c(y, mean(all_gf[id,'brain_methy']/all_gf[id,'brain_total']))
  # print(c(x, mean(all_gf[id,'brain_methy']/all_gf[id,'brain_total'])))
}
plot(y)

x = c()
y = c()
for(value in sort(unique(all_gf$es_medip))) {
  id = which(all_gf$es_medip == value)
  x = c(x, value)
  y = c(y, mean(all_gf[id,'es_methy']/all_gf[id,'es_total']))
  # print(c(x, mean(all_gf[id,'brain_methy']/all_gf[id,'brain_total'])))
}
plot(y)

train_feature_vec = all_gf$brain_medip
cluster_label = EqualFreq2(all_gf$pval, n = 10)
feature_vec = train_feature_vec


rank_convert <- function(x) {
  res = rep(1, length(x))
  for(value in unique(x)){
    res = res + ifelse(x>value, 1, 0)
  }
  return(res) 
}


x_cluster = clustering(train_feature_vec, cluster_label, feature_vec)

x_cluster1 = rank_convert(discretize(all_gf$brain_medip, n = 100)[,1])

x_cluster1 = all_gf$brain_mre + 1
id = which(all_gf$brain_mre > 0)
x_cluster1[id] = rank_convert(discretize(all_gf$brain_mre[id], n = 100)[,1]) + 1

# x_cluster1[id] = EqualFreq2(all_gf$brain_mre[id], n = 14) + 1
x_cluster1[id] = discretize(all_gf$brain_mre[id], n = 30)[,1] + 1


x = c()
y = c()
for(value in sort(unique(x_cluster1))) {
  id = which(x_cluster1 == value)
  x = c(x, value)
  y = c(y, mean(all_gf[id,'brain_methy']/all_gf[id,'brain_total']))
}
plot(x,y)

x_cluster2 = EqualFreq2(all_gf$es_medip, n = 20)
x_cluster2 = discretize(all_gf$es_medip, n = 50)[,1]

x = c()
y = c()
for(value in sort(unique(x_cluster2))) {
  id = which(x_cluster2 == value)
  x = c(x, value)
  y = c(y, mean(all_gf[id,'es_methy']/all_gf[id,'es_total']))
  # print(c(x, mean(all_gf[id,'brain_methy']/all_gf[id,'brain_total'])))
}
plot(y)


x_cluster1 = rank_convert(discretize(all_gf$brain_medip, n = 100)[,1])
x_cluster2 = rank_convert(discretize(all_gf$es_medip, n = 100)[,1])

feature = abs(x_cluster1 - x_cluster2)/min(x_cluster1, x_cluster2)
x = c()
y = c()
for(value in sort(unique(feature))) {
  id = which(feature == value)
  x = c(x, value)
  # y = c(y, mean(all_gf[id,'pval']))
  y = c(y, mean(abs(all_gf[id,'es_methy']/all_gf[id,'es_total'] - 
                      all_gf[id,'brain_methy']/all_gf[id,'brain_total'])))
  # y = c(y, mean(all_gf[id,'label']))
}
plot(y)










gr_crf = with(range_crf, GRanges(chr, IRanges(start = start, end = end)))
gr_mnm = with(range_mnm, GRanges(chr, IRanges(start = chrSt, end = chrEnd)))
ovlp = as.data.frame(findOverlaps(gr_crf, gr_mnm))

ovlp = subsetByOverlaps(gr_crf, gr_mnm)
ovlp = as.data.frame(subsetByOverlaps(gr_mnm, gr_crf))

cpg_mnm = subsetByOverlaps(gr_all, gr_mnm)
cpg_mnm = as.data.frame(cpg_mnm)
colnames(cpg_mnm)[1] = 'chr'


cpg_crf = cbind(all_gf[,c('chr','start','end')],res.crf)
cpg_crf = cpg_crf[which(cpg_crf$res.crf==1),]
cpg_crf$chr = 'chr18'
gr_cpg_crf = with(cpg_crf, GRanges(chr, IRanges(start = start, end = end)))

ovlp = as.data.frame(subsetByOverlaps(gr_cpg_crf, gr_mnm))

truth = all_gf[,c('chr','start','end','label')]
truth$chr = 'chr18'

truth[,'ref'] = ref

mg = merge(truth, ovlp, by=c('chr','start','end'))

mg1 = merge(truth, cpg_mnm, by=c('chr','start','end'))

mean(mg$label)
nrow(mg)


colnames(all_sm)

feature_names = colnames(all_sm)[23:64]

for( name in colnames(all_sm)[23:64]) {
  print(name)
  all_sm[,name] = all_sm[,name]/quantile(all_sm[which(all_sm[,name] != 0),name], 0.5) * 10
  print(summary(all_sm[,name]))
}

summary(all_sm)

dt = all_sm[1:1000,c(feature_names, 'fdr')]

kk = kmeans(dt, centers = 10)
table(kk$cluster)

cluster = kk$cluster


w = 1
while(test_order(dt, kk$cluster) == FALSE) {
  w = w+1
  dt$fdr = w * dt$fdr 
  kk = kmeans(dt, centers = 10)
  cluster = kk$cluster
  print(w)
}

test_order <- function(dt, cluster) {
  res = 0
  for(i in 1:9) {
    idx_prev = which(cluster == i)
    idx_after = which(cluster == (i+1))
    max_prev = max(dt$fdr[idx_prev])
    min_after = min(dt$fdr[idx_after])
    
    if(max_prev > min_after) {
      res = 1
    }
  }
  
  if(res == 0){
    return(TRUE)
  } else {
    return(FALSE)
  }
}






for(i in 1:ncol(train_crf)){
  print(length(unique(train_crf[,i])))
}





y = discretize(x, disc = 'equalfreq', nbins = 2000)
nrow(unique(y))
length(unique(x))

EqualFreq2 <- function(x,n){
  nx <- length(x)
  nrepl <- floor(nx/n)
  nplus <- sample(1:n,nx - nrepl*n)
  nrep <- rep(nrepl,n)
  nrep[nplus] <- nrepl+1
  x[order(x)] <- rep(seq.int(n),nrep)
  x
}

y = EqualFreq2(x, 1000)
length(unique(y))
table(y)



all_gf_ref = all_gf





all_gf$brain_medip = rank_convert(discretize(all_gf$brain_medip, n = 100)[,1])
all_gf$es_medip = rank_convert(discretize(all_gf$es_medip, n = 100)[,1])

id = which(all_gf$brain_mre > 0)
all_gf$brain_mre[which(all_gf$brain_mre == 0)] = 1
all_gf$brain_mre[id] = rank_convert(discretize(all_gf$brain_mre[id], n = 100)[,1]) + 1

id = which(all_gf$es_mre > 0)
all_gf$es_mre[which(all_gf$es_mre == 0)] = 1
all_gf$es_mre[id] = rank_convert(discretize(all_gf$es_mre[id], n = 100)[,1]) + 1

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

all_gf$ratio_medip = abs(all_gf$brain_medip - all_gf$es_medip) / min(all_gf$brain_medip, all_gf$es_medip)
all_gf$ratio_medip_10 = abs(all_gf$brain_medip_10 - all_gf$es_medip_10) / min(all_gf$brain_medip_10, all_gf$es_medip_10)
all_gf$ratio_medip_50 = abs(all_gf$brain_medip_50 - all_gf$es_medip_50) / min(all_gf$brain_medip_50, all_gf$es_medip_50)
all_gf$ratio_medip_200 = abs(all_gf$brain_medip_200 - all_gf$es_medip_200) / min(all_gf$brain_medip_200, all_gf$es_medip_200)
all_gf$ratio_medip_1k = abs(all_gf$brain_medip_1k - all_gf$es_medip_1k) / min(all_gf$brain_medip_1k, all_gf$es_medip_1k)

all_gf$ratio_mre = abs(all_gf$brain_mre - all_gf$es_mre) / min(all_gf$brain_mre, all_gf$es_mre)
all_gf$ratio_mre_10 = abs(all_gf$brain_mre_10 - all_gf$es_mre_10) / min(all_gf$brain_mre_10, all_gf$es_mre_10)
all_gf$ratio_mre_50 = abs(all_gf$brain_mre_50 - all_gf$es_mre_50) / min(all_gf$brain_mre_50, all_gf$es_mre_50)
all_gf$ratio_mre_200 = abs(all_gf$brain_mre_200 - all_gf$es_mre_200) / min(all_gf$brain_mre_200, all_gf$es_mre_200)
all_gf$ratio_mre_1k = abs(all_gf$brain_mre_1k - all_gf$es_mre_1k) / min(all_gf$brain_mre_1k, all_gf$es_mre_1k)

all_gf$dif = abs(all_gf$brain_mre*all_gf$es_medip - all_gf$es_mre*all_gf$brain_medip)
all_gf$dif_10 = abs(all_gf$brain_mre_10*all_gf$es_medip_10 - all_gf$es_mre_10*all_gf$brain_medip_10)
all_gf$dif_50 = abs(all_gf$brain_mre_50*all_gf$es_medip_50 - all_gf$es_mre_50*all_gf$brain_medip_50)
all_gf$dif_200 = abs(all_gf$brain_mre_200*all_gf$es_medip_200 - all_gf$es_mre_200*all_gf$brain_medip_200)
all_gf$dif_1k = abs(all_gf$brain_mre_1k*all_gf$es_medip_1k - all_gf$es_mre_1k*all_gf$brain_medip_1k)


all_gf$label = ifelse(all_gf$fdr<0.001, 1, 0)

# feature = discretize(all_gf$brain_mre, n = 20)[,]
feature = discretize(all_gf$ratio_mre_50, n = 100)[,1]
x = c()
y = c()
for(value in sort(unique(feature))) {
  id = which(feature == value)
  x = c(x, value)
  # y = c(y, mean(all_gf[id,'es_methy']/all_gf[id,'es_total']))
  # y = c(y, mean(abs(all_gf[id,'es_methy']/all_gf[id,'es_total'] - all_gf[id,'brain_methy']/all_gf[id,'brain_total'])))
  y = c(y, mean(all_gf[id,'fdr']))
}
plot(y)


write.table(all_gf, '~/Documents/Paper3/Data/Brain_ES_chr18/crf/all_gf.bed', quote = F, sep = " ", col.names = T, row.names = F)

colnames(all_split) = colnames(all_gf)
id = which(!is.na(all$label))
all_split[id,] = all_gf



idx = which(all_gf$dist > 1000)

all_split = all_gf
all_split[idx,] = ' '

idx = c(1, idx)
new = data.frame(matrix(data=c(' '), nrow = nrow(all_gf)+length(idx), ncol = ncol(all_gf)))
for(i in 1:(length(idx)-1)){
  start = idx[i]
  end = idx[i+1] - 1
  new[(start+i):(end+i),] = all_gf[start:end,]
  print(i)
}

all_split = all_gf
i = 1
j = 1
while(i < nrow(all_split)) {
  if(as.numeric(all_split[i,'dist']) > 1000) {
    all_split = rbind(all_split[1:(i-1),], ' ', all_split[i:nrow(all_split),])
    print(j)
    j = j+1
    i = i+1
  }
  i = i+1
}





idx = (1:64716)*5 + 2

summary(all_mm$start[idx+1] - all_mm$end[idx-1])

n = floor((max(all_mm$end)-min(all_mm$start)) / 100)
rand_idx <- rand_parts(min(all_mm$start):max(all_mm$end), n, l=100)

rand_idx$
rand_gr <- with(wgbs_all, GRanges(chr, IRanges(start = start, end = end)))


test_crf = all[, c(class_feature_names,'label')]

set.seed(4321)
train_idx = sample(1:nrow(test_crf), floor(nrow(test_crf)/5), replace = F)
train_crf = test_crf
train_crf[setdiff(1:nrow(test_crf),train_idx),] = ' '

y = train_crf$label
n_seq = 0
for(i in 1:(length(y)-1) ) {
  if(y[i] == ' ' & y[i+1] != ' '){
    n_seq = n_seq + 1
  }
  print(i)
}
print(n_seq)

for(seed in 4321:4421) {
  train_idx = sample(1:nrow(test_crf), floor(nrow(test_crf)/5), replace = F)
  train_crf = test_crf
  train_crf[setdiff(1:nrow(test_crf),train_idx),] = ' '
  
  y = train_crf$label
  n_seq = 0
  for(i in 1:(length(y)-1) ) {
    if(y[i] == ' ' & y[i+1] != ' '){
      n_seq = n_seq + 1
    }
  }
  print(c(seed, n_seq))
}







label = all1$label

dif = label[2:13406403] - label[1:13406402]

sum(dif==1)/sum(label==0) # 0 -> 1: 0.02285902
sum(dif==-1)/sum(label==1) # 1 -> 0: 0.2550832








