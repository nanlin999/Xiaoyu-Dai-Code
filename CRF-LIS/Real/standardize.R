all_mm1$brain_medip = all_mm1$brain_medip/quantile(all_mm1[which(all_mm1$brain_medip>0),'brain_medip'], 0.75) * 10
all_mm1$brain_mre = all_mm1$brain_mre/quantile(all_mm1[which(all_mm1$brain_mre>0),'brain_mre'], 0.75) * 10
all_mm1$es_medip = all_mm1$es_medip/quantile(all_mm1[which(all_mm1$es_medip>0),'es_medip'], 0.75) * 10
all_mm1$es_mre = all_mm1$es_mre/quantile(all_mm1[which(all_mm1$es_mre>0),'es_mre'], 0.75) * 10

set.seed(1234)
train_idx = sample(1:nrow(all_mm1), floor(nrow(all_mm1)/5), replace = F)

# methods: 1 = CAIM, 2 = CACC, 3 = Ameva
cut_point_medip1 <- disc.Topdown(cbind(all_mm1$brain_medip[train_idx], floor(all_mm1$pval*30)[train_idx]), method = 1)$cutp[[1]]
cut_point_mre1 <- disc.Topdown(cbind(all_mm1$brain_mre[train_idx], floor(all_mm1$pval*30)[train_idx]), method = 1)$cutp[[1]]
cut_point_medip2 <- disc.Topdown(cbind(all_mm1$es_medip[train_idx], floor(all_mm1$pval*30)[train_idx]), method = 1)$cutp[[1]]
cut_point_mre2 <- disc.Topdown(cbind(all_mm1$es_medip[train_idx], floor(all_mm1$pval*30)[train_idx]), method = 1)$cutp[[1]]


###############################################

all_mm1$brain_medip = as.vector(uniq_rank(all_mm1$brain_medip))
all_mm1$brain_mre = as.vector(uniq_rank(all_mm1$brain_mre))
all_mm1$es_medip = as.vector(uniq_rank(all_mm1$es_medip))
all_mm1$es_mre = as.vector(uniq_rank(all_mm1$es_mre))

all_mm1$brain_medip = all_mm1$brain_medip/quantile(all_mm1[which(all_mm1$brain_medip>0),'brain_medip'], 0.75) * 10
all_mm1$brain_mre = all_mm1$brain_mre/quantile(all_mm1[which(all_mm1$brain_mre>0),'brain_mre'], 0.75) * 10
all_mm1$es_medip = all_mm1$es_medip/quantile(all_mm1[which(all_mm1$es_medip>0),'es_medip'], 0.75) * 10
all_mm1$es_mre = all_mm1$es_mre/quantile(all_mm1[which(all_mm1$es_mre>0),'es_mre'], 0.75) * 10

all_mm1$brain_medip = clustering(all_mm1$brain_medip, cut_point_medip1)
all_mm1$brain_mre = clustering(all_mm1$brain_mre, cut_point_mre1)
all_mm1$es_medip = clustering(all_mm1$es_medip, cut_point_medip2)
all_mm1$es_mre = clustering(all_mm1$es_mre, cut_point_mre2)

all_mm1$brain_medip = all_mm1$brain_medip/quantile(all_mm1[which(all_mm1$brain_medip>0),'brain_medip'], 0.75) * 10
all_mm1$brain_mre = all_mm1$brain_mre/quantile(all_mm1[which(all_mm1$brain_mre>0),'brain_mre'], 0.75) * 10
all_mm1$es_medip = all_mm1$es_medip/quantile(all_mm1[which(all_mm1$es_medip>0),'es_medip'], 0.75) * 10
all_mm1$es_mre = all_mm1$es_mre/quantile(all_mm1[which(all_mm1$es_mre>0),'es_mre'], 0.75) * 10

#####################################

all_mm2$skin_medip = as.vector(uniq_rank(all_mm2$skin_medip))
all_mm2$skin_mre = as.vector(uniq_rank(all_mm2$skin_mre))
all_mm2$es_medip = as.vector(uniq_rank(all_mm2$es_medip))
all_mm2$es_mre = as.vector(uniq_rank(all_mm2$es_mre))

all_mm2$skin_medip = all_mm2$skin_medip/quantile(all_mm2[which(all_mm2$skin_medip>0),'skin_medip'], 0.75) * 10
all_mm2$skin_mre = all_mm2$skin_mre/quantile(all_mm2[which(all_mm2$skin_mre>0),'skin_mre'], 0.75) * 10
all_mm2$es_medip = all_mm2$es_medip/quantile(all_mm2[which(all_mm2$es_medip>0),'es_medip'], 0.75) * 10
all_mm2$es_mre = all_mm2$es_mre/quantile(all_mm2[which(all_mm2$es_mre>0),'es_mre'], 0.75) * 10

all_mm2$skin_medip = clustering(all_mm2$skin_medip, cut_point_medip1)
all_mm2$skin_mre = clustering(all_mm2$skin_mre, cut_point_mre1)
all_mm2$es_medip = clustering(all_mm2$es_medip, cut_point_medip2)
all_mm2$es_mre = clustering(all_mm2$es_mre, cut_point_mre2)

all_mm2$skin_medip = all_mm2$skin_medip/quantile(all_mm2[which(all_mm2$skin_medip>0),'skin_medip'], 0.75) * 10
all_mm2$skin_mre = all_mm2$skin_mre/quantile(all_mm2[which(all_mm2$skin_mre>0),'skin_mre'], 0.75) * 10
all_mm2$es_medip = all_mm2$es_medip/quantile(all_mm2[which(all_mm2$es_medip>0),'es_medip'], 0.75) * 10
all_mm2$es_mre = all_mm2$es_mre/quantile(all_mm2[which(all_mm2$es_mre>0),'es_mre'], 0.75) * 10


#############################################

all_mm1$brain_medip_10 = as.numeric(ma(all_mm1$brain_medip, 10))
all_mm1$brain_medip_50 = as.numeric(ma(all_mm1$brain_medip, 50))
all_mm1$brain_medip_200 = as.numeric(ma(all_mm1$brain_medip, 200))
all_mm1$brain_medip_1k = as.numeric(ma(all_mm1$brain_medip, 1000))

all_mm1$es_medip_10 = as.numeric(ma(all_mm1$es_medip, 10))
all_mm1$es_medip_50 = as.numeric(ma(all_mm1$es_medip, 50))
all_mm1$es_medip_200 = as.numeric(ma(all_mm1$es_medip, 200))
all_mm1$es_medip_1k = as.numeric(ma(all_mm1$es_medip, 1000))

all_mm1$brain_mre_10 = as.numeric(ma(all_mm1$brain_mre, 10))
all_mm1$brain_mre_50 = as.numeric(ma(all_mm1$brain_mre, 50))
all_mm1$brain_mre_200 = as.numeric(ma(all_mm1$brain_mre, 200))
all_mm1$brain_mre_1k = as.numeric(ma(all_mm1$brain_mre, 1000))

all_mm1$es_mre_10 = as.numeric(ma(all_mm1$es_mre, 10))
all_mm1$es_mre_50 = as.numeric(ma(all_mm1$es_mre, 50))
all_mm1$es_mre_200 = as.numeric(ma(all_mm1$es_mre, 200))
all_mm1$es_mre_1k = as.numeric(ma(all_mm1$es_mre, 1000))


list_name <- c('brain_medip_10', 'brain_medip_50', 'brain_medip_200', 'brain_medip_1k', 
               'es_medip_10', 'es_medip_50', 'es_medip_200', 'es_medip_1k',
               'brain_mre_10', 'brain_mre_50', 'brain_mre_200', 'brain_mre_1k',
               'es_mre_10', 'es_mre_50', 'es_mre_200', 'es_mre_1k')

for(name in list_name) {
  all_mm1[,name] = all_mm1[,name]/quantile(all_mm1[which(all_mm1[,name]>0),name], 0.75) * 10
}

add_on = 1

all_mm1$ratio_medip = abs(all_mm1$brain_medip - all_mm1$es_medip) / min(all_mm1$brain_medip+add_on, all_mm1$es_medip+add_on)
all_mm1$ratio_medip_10 = abs(all_mm1$brain_medip_10 - all_mm1$es_medip_10) / min(all_mm1$brain_medip_10+add_on, all_mm1$es_medip_10+add_on)
all_mm1$ratio_medip_50 = abs(all_mm1$brain_medip_50 - all_mm1$es_medip_50) / min(all_mm1$brain_medip_50+add_on, all_mm1$es_medip_50+add_on)
all_mm1$ratio_medip_200 = abs(all_mm1$brain_medip_200 - all_mm1$es_medip_200) / min(all_mm1$brain_medip_200+add_on, all_mm1$es_medip_200+add_on)
all_mm1$ratio_medip_1k = abs(all_mm1$brain_medip_1k - all_mm1$es_medip_1k) / min(all_mm1$brain_medip_1k+add_on, all_mm1$es_medip_1k+add_on)

all_mm1$ratio_mre = abs(all_mm1$brain_mre - all_mm1$es_mre) / min(all_mm1$brain_mre+add_on, all_mm1$es_mre+add_on)
all_mm1$ratio_mre_10 = abs(all_mm1$brain_mre_10 - all_mm1$es_mre_10) / min(all_mm1$brain_mre_10+add_on, all_mm1$es_mre_10+add_on)
all_mm1$ratio_mre_50 = abs(all_mm1$brain_mre_50 - all_mm1$es_mre_50) / min(all_mm1$brain_mre_50+add_on, all_mm1$es_mre_50+add_on)
all_mm1$ratio_mre_200 = abs(all_mm1$brain_mre_200 - all_mm1$es_mre_200) / min(all_mm1$brain_mre_200+add_on, all_mm1$es_mre_200+add_on)
all_mm1$ratio_mre_1k = abs(all_mm1$brain_mre_1k - all_mm1$es_mre_1k) / min(all_mm1$brain_mre_1k+add_on, all_mm1$es_mre_1k+add_on)

all_mm1$dif = abs((all_mm1$brain_mre+add_on)*(all_mm1$es_medip+add_on) - (all_mm1$es_mre+add_on)*(all_mm1$brain_medip+add_on))
all_mm1$dif_10 = abs((all_mm1$brain_mre_10+add_on)*(all_mm1$es_medip_10+add_on) - (all_mm1$es_mre_10+add_on)*(all_mm1$brain_medip_10+add_on))
all_mm1$dif_50 = abs((all_mm1$brain_mre_50+add_on)*(all_mm1$es_medip_50+add_on) - (all_mm1$es_mre_50+add_on)*(all_mm1$brain_medip_50+add_on))
all_mm1$dif_200 = abs((all_mm1$brain_mre_200+add_on)*(all_mm1$es_medip_200+add_on) - (all_mm1$es_mre_200+add_on)*(all_mm1$brain_medip_200+add_on))
all_mm1$dif_1k = abs((all_mm1$brain_mre_1k+add_on)*(all_mm1$es_medip_1k+add_on) - (all_mm1$es_mre_1k+add_on)*(all_mm1$brain_medip_1k+add_on))


all_mm2$skin_medip_10 = as.numeric(ma(all_mm2$skin_medip, 10))
all_mm2$skin_medip_50 = as.numeric(ma(all_mm2$skin_medip, 50))
all_mm2$skin_medip_200 = as.numeric(ma(all_mm2$skin_medip, 200))
all_mm2$skin_medip_1k = as.numeric(ma(all_mm2$skin_medip, 1000))

all_mm2$es_medip_10 = as.numeric(ma(all_mm2$es_medip, 10))
all_mm2$es_medip_50 = as.numeric(ma(all_mm2$es_medip, 50))
all_mm2$es_medip_200 = as.numeric(ma(all_mm2$es_medip, 200))
all_mm2$es_medip_1k = as.numeric(ma(all_mm2$es_medip, 1000))

all_mm2$skin_mre_10 = as.numeric(ma(all_mm2$skin_mre, 10))
all_mm2$skin_mre_50 = as.numeric(ma(all_mm2$skin_mre, 50))
all_mm2$skin_mre_200 = as.numeric(ma(all_mm2$skin_mre, 200))
all_mm2$skin_mre_1k = as.numeric(ma(all_mm2$skin_mre, 1000))

all_mm2$es_mre_10 = as.numeric(ma(all_mm2$es_mre, 10))
all_mm2$es_mre_50 = as.numeric(ma(all_mm2$es_mre, 50))
all_mm2$es_mre_200 = as.numeric(ma(all_mm2$es_mre, 200))
all_mm2$es_mre_1k = as.numeric(ma(all_mm2$es_mre, 1000))


list_name <- c('skin_medip_10', 'skin_medip_50', 'skin_medip_200', 'skin_medip_1k', 
               'es_medip_10', 'es_medip_50', 'es_medip_200', 'es_medip_1k',
               'skin_mre_10', 'skin_mre_50', 'skin_mre_200', 'skin_mre_1k',
               'es_mre_10', 'es_mre_50', 'es_mre_200', 'es_mre_1k')

for(name in list_name) {
  all_mm2[,name] = all_mm2[,name]/quantile(all_mm2[which(all_mm2[,name]>0),name], 0.75) * 10
}


all_mm2$ratio_medip = abs(all_mm2$skin_medip - all_mm2$es_medip) / min(all_mm2$skin_medip+add_on, all_mm2$es_medip+add_on)
all_mm2$ratio_medip_10 = abs(all_mm2$skin_medip_10 - all_mm2$es_medip_10) / min(all_mm2$skin_medip_10+add_on, all_mm2$es_medip_10+add_on)
all_mm2$ratio_medip_50 = abs(all_mm2$skin_medip_50 - all_mm2$es_medip_50) / min(all_mm2$skin_medip_50+add_on, all_mm2$es_medip_50+add_on)
all_mm2$ratio_medip_200 = abs(all_mm2$skin_medip_200 - all_mm2$es_medip_200) / min(all_mm2$skin_medip_200+add_on, all_mm2$es_medip_200+add_on)
all_mm2$ratio_medip_1k = abs(all_mm2$skin_medip_1k - all_mm2$es_medip_1k) / min(all_mm2$skin_medip_1k+add_on, all_mm2$es_medip_1k+add_on)

all_mm2$ratio_mre = abs(all_mm2$skin_mre - all_mm2$es_mre) / min(all_mm2$skin_mre+add_on, all_mm2$es_mre+add_on)
all_mm2$ratio_mre_10 = abs(all_mm2$skin_mre_10 - all_mm2$es_mre_10) / min(all_mm2$skin_mre_10+add_on, all_mm2$es_mre_10+add_on)
all_mm2$ratio_mre_50 = abs(all_mm2$skin_mre_50 - all_mm2$es_mre_50) / min(all_mm2$skin_mre_50+add_on, all_mm2$es_mre_50+add_on)
all_mm2$ratio_mre_200 = abs(all_mm2$skin_mre_200 - all_mm2$es_mre_200) / min(all_mm2$skin_mre_200+add_on, all_mm2$es_mre_200+add_on)
all_mm2$ratio_mre_1k = abs(all_mm2$skin_mre_1k - all_mm2$es_mre_1k) / min(all_mm2$skin_mre_1k+add_on, all_mm2$es_mre_1k+add_on)

all_mm2$dif = abs((all_mm2$skin_mre+add_on)*(all_mm2$es_medip+add_on) - (all_mm2$es_mre+add_on)*(all_mm2$skin_medip+add_on))
all_mm2$dif_10 = abs((all_mm2$skin_mre_10+add_on)*(all_mm2$es_medip_10+add_on) - (all_mm2$es_mre_10+add_on)*(all_mm2$skin_medip_10+add_on))
all_mm2$dif_50 = abs((all_mm2$skin_mre_50+add_on)*(all_mm2$es_medip_50+add_on) - (all_mm2$es_mre_50+add_on)*(all_mm2$skin_medip_50+add_on))
all_mm2$dif_200 = abs((all_mm2$skin_mre_200+add_on)*(all_mm2$es_medip_200+add_on) - (all_mm2$es_mre_200+add_on)*(all_mm2$skin_medip_200+add_on))
all_mm2$dif_1k = abs((all_mm2$skin_mre_1k+add_on)*(all_mm2$es_medip_1k+add_on) - (all_mm2$es_mre_1k+add_on)*(all_mm2$skin_medip_1k+add_on))




##########################################

idx1 = 1:nrow(all_mm1)
idx2 = nrow(all_mm1) + 1:nrow(all_mm2)

for(i in 26:49) {
  # all_mm[,i] = all_mm[,i]/quantile(all_mm[which(all_mm[,i]>0),i], 0.75) * 10
  all_mm[idx1,i] = all_mm[idx1,i]/quantile(all_mm[which(all_mm[idx1,i]>0),i], 0.75) * 10
  all_mm[idx2,i] = all_mm[idx2,i]/quantile(all_mm[which(all_mm[idx2,i]>0),i], 0.75) * 10
}


for(binlength in c('10', '50', '200', '500', '1000', '200')) {
  name1 = paste('ratio_medip_', binlength, sep='')
  name2 = paste('ratio_mre_', binlength, sep='')
  name3 = paste('dif_', binlength, sep='')
  
  medip1 = paste('Medip1_', binlength, sep='')
  medip2 = paste('Medip2_', binlength, sep='')
  mre1 = paste('MRE1_', binlength, sep='')
  mre2 = paste('MRE2_', binlength, sep='')
  
  all_mm[,name1] = abs(all_mm[,medip1] - all_mm[,medip2]) / pmin(all_mm[,medip1], all_mm[,medip2])
  all_mm[,name2] = abs(all_mm[,mre1] - all_mm[,mre2]) / pmin(all_mm[,mre1], all_mm[,mre2])
  all_mm[,name3] = all_mm[,medip1] * all_mm[,mre2] -  all_mm[,medip2] * all_mm[,mre1]
  
  all_mm[which(is.na(all_mm[,name2])),name2] = 0
}








