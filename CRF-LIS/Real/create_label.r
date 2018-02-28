Brain_ratio = read.table('~/Documents/Paper3/Data/Brain/HuFGM02_BrainGerminalMatrix_Bisulfite-Seq_A04698_CpG.bedGraph')
Brain_density = read.table('~/Documents/Paper3/Data/Brain/HuFGM02_BrainGerminalMatrix_Bisulfite-Seq_A04698_density.bedGraph')

Brain_ratio = Brain_ratio[Brain_ratio$V1 == 'chr18',]
Brain_density = Brain_density[Brain_density$V1 == 'chr18',]


Skin_ratio = read.table('~/Documents/Paper3/Data/Skin/E058.fm.bedGraph')
Skin_density = read.table('~/Documents/Paper3/Data/Skin/E058.rc.bedGraph')

Skin_ratio = Skin_ratio[Skin_ratio$V1 == 'chr18',]
Skin_density = Skin_density[Skin_density$V1 == 'chr18',]

write.table(Brain_ratio, file = '/Users/xiaoyudai/Documents/Paper3/Data/train/tmp.txt', quote = F, sep = " ", col.names = F, row.names = F)
write.table(Brain_density, file = '/Users/xiaoyudai/Documents/Paper3/Data/train/tmp1.txt', quote = F, sep = " ", col.names = F, row.names = F)
write.table(Skin_ratio, file = '/Users/xiaoyudai/Documents/Paper3/Data/train/tmp2.txt', quote = F, sep = " ", col.names = F, row.names = F)
write.table(Skin_density, file = '/Users/xiaoyudai/Documents/Paper3/Data/train/tmp3.txt', quote = F, sep = " ", col.names = F, row.names = F)


summary(as.numeric(Brain_ratio$V2))

V2 = seq(0,78019500,500)
V3 = V2 + 500
V1 = rep('chr1', length(V2))
V4 = rep(0, length(V2))
V5 = rep(0, length(V2))
V6 = rep(0, length(V2))
V7 = rep(0, length(V2))

dt = data.frame(V1, V2, V3, V4, V5, V6, V7)

for (i in 1:length(V2)) {
  print(i)
  min = dt$V2[i]
  max = dt$V3[i]
  
  idx = which(Brain_ratio$V2 >= min & Brain_ratio$V3 < max)
  r1 = Brain_ratio[idx,'V4']
  d1 = Brain_density[idx,'V4']
  dt[i,'V4'] = round(sum(r1*d1))
  dt[i,'V5'] = sum(d1)
  
  idx = which(Skin_ratio$V2 >= min & Skin_ratio$V3 < max)
  r1 = Skin_ratio[idx,'V4']
  d1 = Skin_density[idx,'V4']
  dt[i,'V6'] = round(sum(r1*d1))
  dt[i,'V7'] = sum(d1)
}
  
write.table(dt, '~/Documents/Paper3/Data/train/dt.txt')















