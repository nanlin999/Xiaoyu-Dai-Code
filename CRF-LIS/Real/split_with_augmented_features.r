



split_with_augmented_features <- function(all, label.threshold=0.01, n_clusters=10, genome_feature) {
  
  # label = rep(0, nrow(all))
  # idx = which( (all$brain_methy/all$brain_total - all$es_methy/all$es_total) > 0.3 )
  # label[idx] = 1
  all[,'label'] = ifelse(all$fdr < 0.001, 1, 0)
  
  # label = ifelse(all$fdr < label.threshold, 1, 0)
  # label = discretize(all$fdr, disc = 'equalfreq', nbins = 9)
  # all[,'label'] = label

  
  # feature_names = colnames(all)[c(16,49:68)]
  feature_names = colnames(all)[c(28,45:64)]  
  print(feature_names)
  # class_feature_names = colnames(all)[c(14,16:27)]                      
  class_feature_names = colnames(all)[c(15:27)]  
  
  # all[,'fdr_class'] = discretize(all$fdr, disc = 'equalfreq', nbins = 10)
  
  # train_seq = rand_parts(seq(1,nrow(all)), n_train_seq, l_train_seq)
  
  # train = c()
  # for (i in 1:nrow(train_seq)){
  #   print(i)
  #   start = train_seq[i,1]
  #   end = train_seq[i,2]
  #   train = rbind(train, all[start:end, c(feature_names,'label')])
  # }
  
  for (name in feature_names) {
    # all[,name] = kmeans(cbind(scaling(all[,name]), kmeans_label), n_clusters)$cluster
    # print(name)
    # all[,name] = clustering(scaling(all[,name]), cluster_label, train.size)
    print(name)
    new_name = paste(name, '_class', sep = "")
    # if (grepl('mre', name)) {
    #     all[,new_name] = discretize(all[,name], disc = 'equalwidth', nbins = n_clusters)
    #   } else {
    #     all[,new_name] = discretize(all[,name], disc = 'equalfreq', nbins = n_clusters)
    #   }
    # all[,new_name] = clustering(train[,name], train[,'fdr_class'], all[,name])
    # all[,new_name] = clustering(train[,name], train[,'label'], all[,name])
    # all[,new_name] = discretize(all[,name], disc = 'equalfreq', nbins = n_clusters)
    all[,new_name] = EqualFreq2(all[,name], n = n_clusters)
    print(length(unique(all[,new_name])))
    class_feature_names = c(class_feature_names, new_name)
  }
  all_feature_names = c(feature_names, class_feature_names)
  # print(all_feature_names)
  
  
  # train = c()
  # for (i in 1:nrow(train_seq)){
  #   print(i)
  #   start = train_seq[i,1]
  #   end = train_seq[i,2]
  #   train = rbind(train, all[start:end, c(class_feature_names,'label')], c(' '))
  # }
  # 
  
  # every_n = 5
  # idx_split = c(1 : (round(nrow(all)/every_n)-1) ) * every_n
  # train = all[, c(class_feature_names,'label')]
  # 
  # for(i in (every_n/5):(every_n-1)){
  #   train[idx_split+i,] = c(' ')
  # }
  # train = rbind(train, all[1:1000, c(class_feature_names,'label')])
  # 
  # train = all[1:64000, c(class_feature_names,'label')]
  
  # train = all[1:round(nrow(all)/5), c('dist', class_feature_names,'label')]
  
  # test = all
  # idx = 1
  # while(idx < (nrow(test)+1)) {
  #   if(test[idx,'cpgprev'] > 750) {
  #     test = rbind(test[1:idx,], ' ', test[(idx+1):nrow(test),] )
  #     idx = idx+1
  #   }
  #   idx = idx+1
  #   print(idx)
  # }
  
  # all_split = all_gf[,c(14:32,49:68)]
  # i = 1
  # j = 1
  # while(i < nrow(all_split)) {
  #   if(as.numeric(all_split[i,'dist']) > 1000) {
  #     all_split = rbind(all_split[1:(i-1),], ' ', all_split[i:nrow(all_split),])
  #     print(j)
  #     j = j+1
  #     i = i+1
  #   }
  #   i = i+1
  # }
  # 
  # write.table(all_split, file = '/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_chr18/crf/all_split_1000.bed', quote = F, sep = " ", col.names = T, row.names = F)
  
  
  
  # test[which(test$cpgprev>750),] = ' '
  test = all_split[, c(class_feature_names,'label')]
  
  every_n = 10
  idx_split = c(1 : (round(nrow(all)/every_n)-1) ) * every_n
  train = test

  for(i in (every_n/5):(every_n-1)){
    train[idx_split+i,] = c(' ')
  }
  train = rbind(train, test[1:1000, ])

  
  idx = which(train$label != ' ' & (!is.na(train$label)) )
  train_rnn = train[idx,]
  
  idx = which(test$label != ' ' & (!is.na(test$label)) )
  test_rnn = test[idx,]
  
  
  # train[idx_0,] = c(' ')
  # test[idx_0,] = c(' ')
  
  return(list(train_crf = train,
              test_crf = test,
              all = all,
              train_rnn=train_rnn,
              test_rnn = test_rnn))
  
  # split_ratio = 0.2
  # 
  # return(list(train_crf = all[1:round(split_ratio*nrow(all)), c(class_feature_names,'label')], 
  #             test_crf = all[round(split_ratio*nrow(all)):nrow(all), c(class_feature_names,'label')], 
  #             all = all, 
  #             train_rnn= all[1:round(split_ratio*nrow(all)), c(8:11,21,23:64,22,20)], 
  #             test_rnn = all[round(split_ratio*nrow(all)):nrow(all), c(8:11,21,23:64,22,20)]))
  
}

ma <- function(x,n=5){
  res = filter(x,rep(1/n,n), sides=2)
  idx = which(is.na(res))
  res[idx] = x[idx]
  return(res)
}


scaling <- function(x) {
  x = as.numeric(x)
  # return(x / sqrt(sum(x^2))) 
  x = (x - min(x))/(max(x)-min(x))
  return(x)
}


rand_parts <- function(seq, n, l) {
  set.seed(1234)
  indices = seq(1, length(seq) - (l - 1) * n)
  result = c()
  offset = 0
  for (i in sort(sample(indices, n))) {
    i = i + offset
    result = rbind(result, c(i, i+l-1))
    offset = offset + l -1
  }
  return(result)
}

sm <- function(all, window_size){
  brain_medip_sm = c()
  brain_mre_sm = c()
  es_medip_sm = c()
  es_mre_sm = c()
  for (i in 1:nrow(all)) {
    print(i)
    low = all[i,'start'] - window_size/2
    up = all[i,'end'] + window_size/2
    idx = which(all$start > low & all$end < up)
    brain_medip_sm = c(brain_medip_sm, mean(all[idx,'brain_medip']))
    brain_mre_sm = c(brain_mre_sm, mean(all[idx,'brain_mre']))
    es_medip_sm = c(es_medip_sm, mean(all[idx,'es_medip']))
    es_mre_sm = c(es_mre_sm, mean(all[idx,'es_mre']))
  }
  return(list(brain_medip_sm=brain_medip_sm,
              brain_mre_sm=brain_mre_sm,
              es_medip_sm=es_medip_sm,
              es_mre_sm=es_mre_sm))
}

EqualFreq2 <- function(x,n){
  nx <- length(x)
  nrepl <- floor(nx/n)
  nplus <- sample(1:n,nx - nrepl*n)
  nrep <- rep(nrepl,n)
  nrep[nplus] <- nrepl+1
  x[order(x)] <- rep(seq.int(n),nrep)
  x
}

uniq_rank <- function(x){
  x_unique <- unique(x)
  x_ranks <- rank(x_unique)
  res = x_ranks[match(x,x_unique)]
  return(res)
}



create_cluster_feature <- function(all, n_clusters) {
  
 
  feature_names = colnames(all)[c(16,44:63)] # all_mm
  # feature_names = colnames(all)[c(16,8:11,28:29,32:33,36:37,40:41,59:63)] 
  # feature_names = colnames(all)[c(8,20:25,50:64)] 
  # feature_names = colnames(all)[c(16,20:49)] 
  print(feature_names)
  # class_feature_names = colnames(all)[c(6,9:19)]  # all_mm                   
  class_feature_names = colnames(all)[c(14,17:27)]  # all_mm                   
  # class_feature_names = colnames(all)[c(15:27)]  # all_mm                   
  
  # all$label = ifelse(all$pval < 0.001, 1, 0)
  
  all$label = ifelse(all$fdr < 0.001, 1, 0)
  
  # for(name in feature_names){
  #   hist(all[,name], xlab = name)
  # }
  
  # all[which(is.na(all$label)),] = ' '
  # idx = which(all$label != ' ')
  
  # for(name in c('brain_methy', 'brain_total', 'es_methy', 'es_total')){
  #   all[,name] = as.numeric(all[, name])
  # }
  
  # all$label[idx] = ifelse( abs(all$brain_methy[idx]/all$brain_total[idx] - all$es_methy[idx]/all$es_total[idx])>0.3, 1, 0)
  
  # prob_gap = c(rep(0.001, 20), 0.001)
  # set.seed(1234)
  train_idx = sample(1:323583, floor(323583/5), replace = F)
  # train_idx = 323581 + sample(1:526818, floor(526818/5), replace = F)
  
  # pval_label = floor(all[train_idx,'pval']*10)

  for (i in 1:length(feature_names)) {
    name = feature_names[i]
    print(name)
    all[, name] = as.numeric(all[, name])
    new_name = paste(name, '_class', sep = "")
    # all[idx,new_name] = EqualFreq2(all[idx,name], n = n_clusters)
    # print(length(unique(all[,name])))
    all[,new_name] = discretize(all[,name], disc = 'equalfreq', nbins = 8000)
    
    # cut_point <- c(-1e-5, as.numeric(quantile(all[train_idx, name], probs = seq(0, 1, prob_gap[i]))), 1e20)
    # cut_point = sort(unique(cut_point))
    
    # if(grepl('mm_q', name)) {
    #   cut_point = c(-1e-20, 0.000000000001, 0.0000001, 0.0001, 0.01, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1e20)
    #   # cut_point = c(-1e-20, 0.000000000001, 0.0000001, 0.0001, 0.01, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, seq(0.9,1,by=0.0001), 1e20)
    #   # cut_point <- c(-1e-20, as.numeric(quantile(all[train_idx, name], probs = seq(0, 1, 0.00001))), 1e20)
    #   # cut_point = sort(unique(cut_point))
    # } else {
    #   cut_point <- c(-1e-20, as.numeric(quantile(all[train_idx, name], probs = seq(0, 1, 0.00001))), 1e20)
    #   cut_point = sort(unique(cut_point))
    # }
    # 
    # all[,new_name] = cut(all[,name], breaks = cut_point, labels = FALSE)
    # all[which(is.na(all[,new_name])), new_name] = max(all[,new_name]) + 1
    
    print(length(unique(all[,new_name])))
    # print(table(all[,new_name]))
    # plot(table(all[,new_name]))
    class_feature_names = c(class_feature_names, new_name)
  }
  
  # train = all[1:70000,c(class_feature_names,'label')]
  # test = all[70001:nrow(all), c(class_feature_names,'label')]
  
  
  # train_crf = all[, c(class_feature_names,'label')]
  # test_crf = all[, c(class_feature_names,'label')]
  
  # gap_train = which( (all$chr != 'chr1' & all$chr != 'chr2') | all$dist > 1000)
  # gap_test = which(all$chr != 'chr18' | all$dist > 1000)
 
  # train_crf[gap_train,] = ' '
  # test_crf[gap_test,] = ' '
  
  # train_crf = all[1:323583, c(class_feature_names,'label')]
  # test_crf = all[323584:nrow(all), c(class_feature_names,'label')]
  
  train_crf = all[, c(class_feature_names,'label')]
  test_crf = all[, c(class_feature_names,'label')]
  
  set.seed(1234)
  train_idx = sample(1:nrow(train_crf), floor(nrow(train_crf)/5), replace = F)
  train_crf[setdiff(1:nrow(train_crf),train_idx),] = ' '
  
  # test_crf = all[1:323583, c(class_feature_names,'label')]
  
  write.table(train_crf, file = '/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_chr18/crf/train_crf', quote = F, sep = " ", col.names = F, row.names = F)
  write.table(test_crf, file = '/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_chr18/crf/test_crf', quote = F, sep = " ", col.names = F, row.names = F)
  
  
  
  train_crf = all[1:nrow(comb_a), c(class_feature_names,'label')]
  train_crf[setdiff(1:nrow(comb_a),train_idx),] = ' '
  test_crf = all[nrow(comb_a):nrow(all), c(class_feature_names,'label')]
  
  
  train_crf = all[1:floor(323583/5*2), c(class_feature_names,'label')]
  test_crf = all[floor(323583/5*2):nrow(all), c(class_feature_names,'label')]
  
  
  test_crf = all[323582:nrow(all), c(class_feature_names,'label')]
  write.table(test_crf, file = '/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_chr18/crf/test_crf', quote = F, sep = " ", col.names = F, row.names = F)
  
  
  
  
  test_crf[sample(1:nrow(test_crf), floor(nrow(test_crf)/5)*4, replace = F),] = ' '
  
  
  train_crf = all[323584:500000, c(class_feature_names,'label')]
  test_crf = all[500000:nrow(all), c(class_feature_names,'label')]
  
  all[which(all$dist>750),] = ' '
  
  train_crf = all[1:64000, c(class_feature_names,'label')]
  test_crf = all[1:nrow(all), c(class_feature_names,'label')]
  
  set.seed(1234)
  test_crf = all[, c(class_feature_names,'label')]
  train_idx = sample(1:nrow(test_crf), floor(nrow(test_crf)/5), replace = F)
  train_crf = test_crf
  train_crf[setdiff(1:nrow(test_crf),train_idx),] = ' '
  test_crf[train_idx,] = ' '
  
  
  train_idx = which(train_crf$label[1:nrow(all)] != ' ')
  test_mlp = all[, c(14,16:27,44:63,15)]
  for(i in c(2,14:33)) {
    test_mlp[,i] = scaling(test_mlp[,i])
  }
  train_mlp = test_mlp[train_idx,]
  
  
  rnn_sub_length = 10
  
  n_sub = floor(nrow(all)/rnn_sub_length)
  test_rnn = all[1:(rnn_sub_length*n_sub), c(14,16:27,44:63,15)]
  for(i in c(2,14:33)) {
    test_rnn[,i] = scaling(test_rnn[,i])
  }
  # set.seed(1234)
  # train_start = sample(0:(n_sub-1), floor(n_sub/5), replace = F) * rnn_sub_length
  # train_idx = c()
  # for(i in 1:rnn_sub_length) {
  #   train_idx = c(train_idx, train_start + i)
  # }
  
  train_start = c(0:(n_sub-1)) * rnn_sub_length
  train_idx = c()
  for(i in 1:(rnn_sub_length/5)) {
    train_idx = c(train_idx, train_start + i)
  }
  
  train_idx = sort(train_idx)
  train_rnn = test_rnn[train_idx,]
  
  write.table(train_rnn, file = '/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_chr18/crf/train_rnn', quote = F, sep = "\t", col.names = F, row.names = F)
  write.table(test_rnn, file = '/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES_chr18/crf/test_rnn', quote = F, sep = "\t", col.names = F, row.names = F)
  
  
  
  
  
  return(list(train_crf = train_crf,
              test_crf = test_crf,
              train_mlp=train_mlp,
              test_mlp = test_mlp,
              train_rnn=train_rnn,
              test_rnn = test_rnn))
              # all = all,
              # train_rnn=train_rnn,
              # test_rnn = test_rnn))
}









