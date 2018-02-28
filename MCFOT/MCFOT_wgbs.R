
D1 <- read.table('~/Documents/Paper2/Real/WGBS/CD4_ES/ALL.txt', header = T) # CD4 vs ES cells, main comparison
D2 <- read.table('~/Documents/Paper2/Real/WGBS/CD8_ES/ALL.txt', header = T) # CD8 vs ES cells, used for ordering

# for trial
# C1 <- D1[which(D1$chrome == 'chr1'),c(1,2,3,12,13)]
# C2 <- D2[which(D2$chrome == 'chr1'),c(1,2,3,12,13)]
C1 <- D1[1:100000,c(1,2,3,12,13)]
C2 <- D2[1:100000,c(1,2,3,12,13)]


all <- merge(C1,C2, by='start')

idx = sort(all$porg.y,index.return = T, decreasing = F)$ix
new <- all[idx,]
p.org <- new[,"porg.x"]
p.next <- new[,"pnext.x"]

p.org[which(p.org==1)] = 0.999999999 # avoid conflict when pval=1 in HingeExp function



head(all)

plot(p.org,cex=.1)
plot(p.next,cex=.1,col='red')

pdf('/Users/xiaoyudai/Documents/Paper2/Real/WGBS/p1.pdf',width=9,height=5.4)
plot(p.org,cex=.1,xlab='',ylab='p-value')
dev.off()


source('~/Documents/Paper2/Rcode/AT_randp.R')
source('~/Documents/Paper2/Rcode/methods.R')


### Run all methods
# Here we apply multiple comparison methods to one of the unordered sequences of p-values, "pvals_lowdose_ttest" or "pvals_lowdose_permutation_test". We also apply accumulation test methods to the ordered sequence of pvalues, "signed_pvals_reordered".
methods=c(
  'HingeExp (C=2)', 'HingeExp (C=2) with MCF',
  'SeqStep (C=2)', 'SeqStep (C=2) with MCF',
  'ForwardStop', 'ForwardStop with MCF',
  'SepStep', 'SepStep with MCF',
  'Adaptive SepStep', 'Adaptive SepStep with MCF')

# We run methods at a range of desired FDR levels alpha=0.01,0.02,...,0.90. For each method, we evaluate its performance by counting the number of discoveries (i.e. number of rejected hypotheses) that it is able to make, at each alpha.
max_alpha=.2; alphalist = seq(0.01,max_alpha,by=0.01); num_alpha=length(alphalist);

NumRej = matrix(0,length(methods),num_alpha)
for(i in 1:num_alpha){
  alpha=alphalist[i]
  NumRej[which(methods=='HingeExp (C=2)'),i] = HingeExp(p.org, alpha=alpha, C=2)
  NumRej[which(methods=='HingeExp (C=2) with MCF'),i] = HingeExp_randp(p.org, p.next, alpha=alpha, C=2)
  NumRej[which(methods=='SeqStep (C=2)'),i] = SeqStep(p.org, alpha=alpha, C=2)
  NumRej[which(methods=='SeqStep (C=2) with MCF'),i] = SeqStep_randp(p.org, p.next, alpha=alpha, C=2)
  NumRej[which(methods=='ForwardStop'),i] = ForwardStop(p.org, alpha=alpha)
  NumRej[which(methods=='ForwardStop with MCF'),i] = ForwardStop_randp(p.org, p.next, alpha=alpha)
  NumRej[which(methods=='SepStep'),i] = length(SS(p.org, alpha, alpha))
  NumRej[which(methods=='SepStep with MCF'),i] = length(SS.MCF(p.org, p.next, alpha, alpha))
  NumRej[which(methods=='Adaptive SepStep'),i] = length(AS(p.org, alpha, alpha, 0.5))
  NumRej[which(methods=='Adaptive SepStep with MCF'),i] = length(AS.MCF(p.org, p.next, alpha, alpha, 0.5))
  print(i)
}

write.table(NumRej, file = '~/Documents/Paper2/Real/WGBS/res1.txt', col.names = F, quote = F)


# Note that for ForwardStop and HingeExp, we must shift the p-values slightly away from 1 to ensure that we do not get values of infinity.

### Plot results
# Here we plot the results for each method.
cols=c(rep('black',2),rep('green',2),rep('black',2),rep('red',2),rep('blue',2))
ltys=rep(c(1,3),5)
lwds=c(rep(1,10))
pchs=c(20,20,6,6,18,18,4,4,17,17)

pdf('/Users/xiaoyudai/Documents/Paper2/Real/WGBS/res1.pdf',width=10,height=7.4)
n = 60000
plot(0:1,0:1,type='n',xlim=range(alphalist),ylim=c(0,n),xlab=expression(paste('Target FDR level ',alpha)),ylab='# of discoveries',axes=FALSE)
axis(side=1,at=0:10/10)
axis(side=2)
alpha_pt=(1:(10*max_alpha))*10
for(i in length(methods):1){
  points(alphalist,NumRej[i,],type='l',col=cols[i],lty=ltys[i],lwd=lwds[i])
  points(alphalist[alpha_pt],NumRej[i,alpha_pt],col=cols[i],pch=pchs[i],lwd=lwds[i])
}
legend(0,n,methods,col=cols,lty=ltys,lwd=lwds,pch=pchs,seg.len=3,cex=0.9)
box(); 
dev.off()
