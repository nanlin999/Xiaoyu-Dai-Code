ensemble_result = data.frame(matrix(NA, nrow = nrow(all_gf)))



result_cpgi = read.table('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES/crf/result_cpgi')
prob1_cpgi = result_cpgi[,ncol(result_cpgi)]
prob1_cpgi = as.character(prob1_cpgi)
prob1_cpgi = as.numeric(substr(prob1_cpgi,3,10)) 
idx_1 = which(all[, 'is.cpgi'] == 1)
res = rep(NA,nrow(all))
res[idx_1] = prob1_cpgi
ensemble_result$prob1_cpgi = res

result_exoniphy = read.table('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES/crf/result_exoniphy')
prob1_exoniphy = result_exoniphy[,ncol(result_exoniphy)]
prob1_exoniphy = as.character(prob1_exoniphy)
prob1_exoniphy = as.numeric(substr(prob1_exoniphy,3,10)) 
idx_1 = which(all[, 'is.exoniphy'] == 1)
res = rep(NA,nrow(all))
res[idx_1] = prob1_exoniphy
ensemble_result$prob1_exoniphy = res

result_exon = read.table('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES/crf/result_exon')
prob1_exon = result_exon[,ncol(result_exon)]
prob1_exon = as.character(prob1_exon)
prob1_exon = as.numeric(substr(prob1_exon,3,10)) 
idx_1 = which(all[, 'is.exon'] == 1)
res = rep(NA,nrow(all))
res[idx_1] = prob1_exon
ensemble_result$prob1_exon = res

result_intron = read.table('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES/crf/result_intron')
prob1_intron = result_intron[,ncol(result_intron)]
prob1_intron = as.character(prob1_intron)
prob1_intron = as.numeric(substr(prob1_intron,3,10)) 
idx_1 = which(all[, 'is.intron'] == 1)
res = rep(NA,nrow(all))
res[idx_1] = prob1_intron
ensemble_result$prob1_intron = res

result_nestrepeat = read.table('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES/crf/result_nestrepeat')
prob1_nestrepeat = result_nestrepeat[,ncol(result_nestrepeat)]
prob1_nestrepeat = as.character(prob1_nestrepeat)
prob1_nestrepeat = as.numeric(substr(prob1_nestrepeat,3,10)) 
idx_1 = which(all[, 'is.nestrepeat'] == 1)
res = rep(NA,nrow(all))
res[idx_1] = prob1_nestrepeat
ensemble_result$prob1_nestrepeat = res

result_simplerepeat = read.table('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES/crf/result_simplerepeat')
prob1_simplerepeat = result_simplerepeat[,ncol(result_simplerepeat)]
prob1_simplerepeat = as.character(prob1_simplerepeat)
prob1_simplerepeat = as.numeric(substr(prob1_simplerepeat,3,10)) 
idx_1 = which(all[, 'is.simplerepeat'] == 1)
res = rep(NA,nrow(all))
res[idx_1] = prob1_simplerepeat
ensemble_result$prob1_simplerepeat = res

result_refgene = read.table('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES/crf/result_refgene')
prob1_refgene = result_refgene[,ncol(result_refgene)]
prob1_refgene = as.character(prob1_refgene)
prob1_refgene = as.numeric(substr(prob1_refgene,3,10)) 
idx_1 = which(all[, 'is.refGene'] == 1)
res = rep(NA,nrow(all))
res[idx_1] = prob1_refgene
ensemble_result$prob1_refgene = res

result_rsmk = read.table('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES/crf/result_rsmk')
prob1_rsmk = result_rsmk[,ncol(result_rsmk)]
prob1_rsmk = as.character(prob1_rsmk)
prob1_rsmk = as.numeric(substr(prob1_rsmk,3,10)) 
idx_1 = which(all[, 'is.rsmk'] == 1)
res = rep(NA,nrow(all))
res[idx_1] = prob1_rsmk
ensemble_result$prob1_rsmk = res

result_utr3 = read.table('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES/crf/result_utr3')
prob1_utr3 = result_utr3[,ncol(result_utr3)]
prob1_utr3 = as.character(prob1_utr3)
prob1_utr3 = as.numeric(substr(prob1_utr3,3,10)) 
idx_1 = which(all[, 'is.utr3'] == 1)
res = rep(NA,nrow(all))
res[idx_1] = prob1_utr3
ensemble_result$prob1_utr3 = res

result_utr5 = read.table('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES/crf/result_utr5')
prob1_utr5 = result_utr3[,ncol(result_utr5)]
prob1_utr5 = as.character(prob1_utr5)
prob1_utr5 = as.numeric(substr(prob1_utr5,3,10)) 
idx_1 = which(all[, 'is.utr5'] == 1)
res = rep(NA,nrow(all))
res[idx_1] = prob1_utr5
ensemble_result$prob1_utr5 = res

result_non = read.table('/Users/xiaoyudai/Documents/Paper3/Data/Brain_ES/crf/result_non')
prob1_non = result_non[,ncol(result_non)]
prob1_non = as.character(prob1_non)
prob1_non = as.numeric(substr(prob1_non,3,10)) 
idx_1 = which( (all$is.cpgi + all$is.exon + all$is.exoniphy + all$is.intron + all$is.refGene + all$is.rsmk 
                + all$is.simplerepeat + all$is.nestrepeat) == 0)
res = rep(NA,nrow(all))
res[idx_1] = prob1_non
ensemble_result$prob1_non = res


head(ensemble_result)

ensemble_result = ensemble_result[,c(2:12)]

write.table(ensemble_result[,c(2:12)], '~/Documents/Paper3/Data/Brain_ES/crf/ensenmble_result', quote = F, sep = " ", col.names = T, row.names = F)

prob1 = rowMeans(ensemble_result, na.rm = T)
prob0 = 1-prob1

truth = all_gf$label
accuracy(truth[idx], prob1[idx])
confusion.matrix(truth[idx], prob1[idx])




















