library(edgeR)
library(statmod)
library(methylMnM)


MRE_1 = read.table('~/Documents/Paper3/Data/Brain/HuFGM02_BrainGerminalMatrix_MRE_A02761.bed')
MRE_1 = MRE_1[MRE_1$V1 == '18', ]

MeDIP_1 = read.table('~/Documents/Paper3/Data/Brain/HuFGM02_BrainGerminalMatrix_MeDIP_A02758.bed') 
MeDIP_1 = MeDIP_1[MeDIP_1$V1 == '18', ]

MRE_2 = read.table('~/Documents/Paper3/Data/Skin/Skin03_Keratinocyte_MRE_A13914.bed')
MRE_2 = MRE_2[MRE_2$V1 == '18', ]

MeDIP_2 = read.table('~/Documents/Paper3/Data/Skin/Skin03_Keratinocyte_MeDIP_A14101.bed') 
MeDIP_2 = MeDIP_2[MeDIP_2$V1 == '18', ]


dirwrite <- "/Users/xiaoyudai/Documents/Paper3/Data/trial/"

write.table(MeDIP_1, file = '/Users/xiaoyudai/Documents/Paper3/Data/trial/MeDIP_1.bed', quote = F, sep = " ", col.names = F, row.names = F)
writefile <- paste(dirwrite, "MeDIP_bin_1.bed", sep = "")
countMeDIPbin('/Users/xiaoyudai/Documents/Paper3/Data/trial/MeDIP_1.bed', writefile = writefile, binlength = 500)

write.table(MRE_1, file = '/Users/xiaoyudai/Documents/Paper3/Data/trial/MRE_1.bed', quote = F, sep = " ", col.names = F, row.names = F)
writefile <- paste(dirwrite, "MRE_bin_1.bed", sep = "")
countMREbin('/Users/xiaoyudai/Documents/Paper3/Data/trial/MRE_1.bed', writefile = writefile, binlength = 500)

# Special case: Skin_MeDIP too large
# write.table(MeDIP_2, file = '/Users/xiaoyudai/Documents/Paper3/Data/trial/MeDIP_2.bed', quote = F, sep = " ", col.names = F, row.names = F)
writefile <- paste(dirwrite, "MeDIP_bin_2.bed", sep = "")
countMeDIPbin('/Users/xiaoyudai/Documents/Paper3/Data/trial/Skin_MeDIP_chr18.bed', writefile = writefile, binlength = 500)

write.table(MRE_2, file = '/Users/xiaoyudai/Documents/Paper3/Data/trial/MRE_2.bed', quote = F, sep = " ", col.names = F, row.names = F)
writefile <- paste(dirwrite, "MRE_bin_2.bed", sep = "")
countMREbin('/Users/xiaoyudai/Documents/Paper3/Data/trial/MRE_2.bed', writefile = writefile, binlength = 500)



density_1 = read.table('~/Documents/Paper3/Data/Brain/HuFGM02_BrainGerminalMatrix_Bisulfite-Seq_A04698_density.bedGraph')
density_1 = density_1[density_1$V1 == 'chr18', ]
write.table(density_1, file = '/Users/xiaoyudai/Documents/Paper3/Data/trial/density_1.bed', quote = F, sep = " ", col.names = F, row.names = F)
writefile = paste(dirwrite, "cpgbin.bed", sep = "")
countcpgbin('/Users/xiaoyudai/Documents/Paper3/Data/trial/density_1.bed', writefile = writefile, binlength = 500)

remove(density_1, MeDIP_1, MeDIP_2, MRE_1, MRE_2)

file <- '~/Documents/Paper3/Data/TriMRE_frags.bed'
file1 <- '~/Documents/Paper3/Data/Brain/HuFGM02_BrainGerminalMatrix_Bisulfite-Seq_A04698_density.bedGraph'
allcpgfile <- paste(dirwrite, "cpgbin.bed", sep = "")
five_Mre_CpGsite <- read.table(file, header = FALSE, as.is = TRUE)
four_Mre_CpGsite <- five_Mre_CpGsite[five_Mre_CpGsite[, 4] != "ACGT", ]
mrecpg.site <- four_Mre_CpGsite[four_Mre_CpGsite[, 4] != "CGCG", ]
writefile <- paste(dirwrite, "three_mre_cpg.bed", sep = "")
countMREcpgbin(mrecpg.site, file.allcpgsite = file1, file.bin = allcpgfile,
               writefile = writefile, binlength = 500)



datafile1 <- paste(dirwrite, "MeDIP_bin_1_tmp.bed", sep = "")
datafile2 <- paste(dirwrite, "MeDIP_bin_2_tmp.bed", sep = "")
datafile3 <- paste(dirwrite, "MRE_bin_1.bed", sep = "")
datafile4 <- paste(dirwrite, "MRE_bin_2.bed", sep = "")
datafile <- c(datafile1, datafile2, datafile3, datafile4)
cpgfile <- paste(dirwrite, "cpgbin_tmp.bed", sep = "")
mrecpgfile <- paste(dirwrite, "three_mre_cpg_tmp.bed", sep = "")
writefile <- paste(dirwrite, "pval_Brain_Skin_chr18.bed", sep = "")
reportfile <- paste(dirwrite, "report_Brain_Skin_chr18.txt", sep = "")
MnM.test(file.dataset = datafile, file.cpgbin = cpgfile,
         file.mrecpgbin = mrecpgfile, writefile = writefile, reportfile = reportfile,
         mreratio = 3/7, method = "XXYY", psd = 2, mkadded = 1, a = 1e-16,
         cut = 100, top = 500)


datafile <- paste(dirwrite, "pval_Brain_Skin_chr18.bed", sep = "")
writefile <- paste(dirwrite, "q_Brain_Skin_chr18.bed", sep = "")
reportfile <- paste(dirwrite, "report_q_Brain_Skin_chr18.bed", sep = "")
MnM.qvalue(datafile, writefile, reportfile)

qval = read.table(writefile, header = T)

qval[which(as.numeric(qval$qvalue) < 0.01), ]



