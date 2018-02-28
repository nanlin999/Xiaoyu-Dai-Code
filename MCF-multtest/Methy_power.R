mcf <- read.table("/Users/xiaoyudai/Documents/multiple-testing/Comparison/CD4_CD8/MCF005.txt",head=T)
q <- read.table("/Users/xiaoyudai/Documents/multiple-testing/Comparison/CD4_CD8/Q005.txt",head=T)

sum(ALL[,17])
sum(ALL[,16][which(ALL[,17]==1)]==0)

library(cp4p)
library(qvalue)
library(ggplot2)

################################## CD4 vs CD8 ##################################

ALL <- read.table("/Users/xiaoyudai/Documents/multiple-testing/Comparison/CD4_CD8/ALL.txt",head=T)
head(ALL)
p.org <- ALL[,12]
p.next <- ALL[,13]
n.test = nrow(ALL)

alpha.all <- seq(0.001,0.2,by=0.001)
power.MCF <- rep(0,length(alpha.all))
power.q <- rep(0,length(alpha.all))

### generate randomized p-values ###
n.rep <- 10
n.ecdf <- n.rep*n.test
randp.ecdf <- rep(0,n.ecdf)
for(i in 0:(n.rep-1)){
  randp.ecdf[(i*n.test+1):((i+1)*n.test)] <- runif(n.test,p.next,p.org)
  print(i)  
}
leng = length(randp.ecdf)

###### q-value method #############
MT <- qvalue(p.org)
QVALUE <- MT$qvalues

######
pi0 = MT$pi0
pi1 <- 1-pi0
#############

for(i.simu in 1:length(alpha.all)){
  alpha = alpha.all[i.simu]
  ####### MCF method ############
  a = 0.00000001
  b = 1
  l = 0
  R = 0
  pFDR = 0
  while( abs(pFDR-alpha)>0.0000001 ){
    l = (a+b)/2
    cdf <- sum(randp.ecdf<l)/leng
    pFDR = (1-pi1)*l/cdf  
    R = round(n.test*cdf) 
    if(pFDR > alpha){
      b=l
    }
    if(pFDR < alpha){
      a=l
    }  
    print(c(a,b,l,pFDR))
  }
  power.MCF[i.simu] = R
  ###### q-value method #############
  power.q[i.simu] = sum(QVALUE<alpha)
  
  print(i.simu)
}


delta <- (ALL[,4]+ALL[,6])/(ALL[,5]+ALL[,7]) - (ALL[,8]+ALL[,10])/(ALL[,9]+ALL[,11])
############# alpha = 0.1 ####################
idx.mcf <- which(ALL[,14]==1)
idx.q <- which(ALL[,15]==1)
idx.dif <- setdiff(idx.mcf,idx.q)
hist(delta[idx.dif])
dat <- data.frame(xx = c(delta[idx.q],delta[idx.dif]),yy = c(rep(0,length(idx.q)),rep(1,length(idx.dif))))
p.delta.11 <- ggplot(dat,aes(x=xx)) + 
          geom_histogram(data=subset(dat,yy == '0'),fill = "red", alpha = 0.2, binwidth=0.05) +
          geom_histogram(data=subset(dat,yy == '1'),fill = "blue", alpha = 0.2, binwidth=0.05) +
          xlab("Methylation level difference") +
          ylab("CD4_CD8") +
          ggtitle(expression(paste('nominal FDR level ',alpha,' = 0.1')))
p.delta.11
############# alpha = 0.05 ####################
idx.mcf <- which(ALL[,16]==1)
idx.q <- which(ALL[,17]==1)
idx.dif <- setdiff(idx.mcf,idx.q)
dat <- data.frame(xx = c(delta[idx.q],delta[idx.dif]),yy = c(rep(0,length(idx.q)),rep(1,length(idx.dif))))
p.delta.12 <- ggplot(dat,aes(x=xx)) + 
  geom_histogram(data=subset(dat,yy == '0'),fill = "red", alpha = 0.2, binwidth=0.05) +
  geom_histogram(data=subset(dat,yy == '1'),fill = "blue", alpha = 0.2, binwidth=0.05) +
  xlab("Methylation level difference") +
  ylab(" ") +
  ggtitle(expression(paste('nominal FDR level ',alpha,' = 0.05')))
p.delta.12

n = length(alpha.all)
power = data.frame(rep(alpha.all,2), c(power.MCF,power.q), 
                 c(rep('MCF',n),rep('q-value',n)))
colnames(power) = c('alpha', 'power', 'method')

p1 <- ggplot(power, aes(x = alpha, y = power, linetype=method)) + 
  geom_line() + 
  xlab("nominal FDR level") + 
  ylab("number of significant CpGs")
p1


################################# CD4 vs ES ###############################

ALL <- read.table("/Users/xiaoyudai/Documents/multiple-testing/Comparison/CD4_ES/ALL.txt",head=T)
head(ALL)
p.org <- ALL[,12]
p.next <- ALL[,13]
n.test = nrow(ALL)

alpha.all <- seq(0.001,0.2,by=0.001)
power.MCF <- rep(0,length(alpha.all))
power.q <- rep(0,length(alpha.all))

### generate randomized p-values ###
n.rep <- 10
n.ecdf <- n.rep*n.test
randp.ecdf <- rep(0,n.ecdf)
for(i in 0:(n.rep-1)){
  randp.ecdf[(i*n.test+1):((i+1)*n.test)] <- runif(n.test,p.next,p.org)
  print(i)  
}
leng = length(randp.ecdf)

###### q-value method #############
MT <- qvalue(p.org)
QVALUE <- MT$qvalues

######
pi0 = MT$pi0
pi1 <- 1-pi0
#############

for(i.simu in 1:length(alpha.all)){
  alpha = alpha.all[i.simu]
  ####### MCF method ############
  a = 0.00000001
  b = 1
  l = 0
  R = 0
  pFDR = 0
  while( abs(pFDR-alpha)>0.0000001 ){
    l = (a+b)/2
    cdf <- sum(randp.ecdf<l)/leng
    pFDR = (1-pi1)*l/cdf  
    R = round(n.test*cdf) 
    if(pFDR > alpha){
      b=l
    }
    if(pFDR < alpha){
      a=l
    }  
    print(c(a,b,l,pFDR))
  }
  power.MCF[i.simu] = R
  ###### q-value method #############
  power.q[i.simu] = sum(QVALUE<alpha)
  
  print(i.simu)
}


delta <- (ALL[,4]+ALL[,6])/(ALL[,5]+ALL[,7]) - (ALL[,8]+ALL[,10])/(ALL[,9]+ALL[,11])
############# alpha = 0.1 ####################
idx.mcf <- which(ALL[,14]==1)
idx.q <- which(ALL[,15]==1)
idx.dif <- setdiff(idx.mcf,idx.q)
hist(delta[idx.dif])
dat <- data.frame(xx = c(delta[idx.q],delta[idx.dif]),yy = c(rep(0,length(idx.q)),rep(1,length(idx.dif))))
p.delta.21 <- ggplot(dat,aes(x=xx)) + 
  geom_histogram(data=subset(dat,yy == '0'),fill = "red", alpha = 0.2, binwidth=0.05) +
  geom_histogram(data=subset(dat,yy == '1'),fill = "blue", alpha = 0.2, binwidth=0.05) +
  xlab("Methylation level difference") +
  ylab("CD4_ES")
p.delta.21
############# alpha = 0.05 ####################
idx.mcf <- which(ALL[,16]==1)
idx.q <- which(ALL[,17]==1)
idx.dif <- setdiff(idx.mcf,idx.q)
dat <- data.frame(xx = c(delta[idx.q],delta[idx.dif]),yy = c(rep(0,length(idx.q)),rep(1,length(idx.dif))))
p.delta.22 <- ggplot(dat,aes(x=xx)) + 
  geom_histogram(data=subset(dat,yy == '0'),fill = "red", alpha = 0.2, binwidth=0.05) +
  geom_histogram(data=subset(dat,yy == '1'),fill = "blue", alpha = 0.2, binwidth=0.05) +
  xlab("Methylation level difference") +
  ylab(" ")
p.delta.22

n = length(alpha.all)
power = data.frame(rep(alpha.all,2), c(power.MCF,power.q), 
                   c(rep('MCF',n),rep('q-value',n)))
colnames(power) = c('alpha', 'power', 'method')

p2 <- ggplot(power, aes(x = alpha, y = power, linetype=method)) + 
  geom_line() + 
  xlab("nominal FDR level") + 
  ylab("number of significant CpGs")
p2




################################# CD8 vs ES ###############################

ALL <- read.table("/Users/xiaoyudai/Documents/multiple-testing/Comparison/CD8_ES/ALL.txt",head=T)
head(ALL)
p.org <- ALL[,12]
p.next <- ALL[,13]
n.test = nrow(ALL)

alpha.all <- seq(0.001,0.2,by=0.001)
power.MCF <- rep(0,length(alpha.all))
power.q <- rep(0,length(alpha.all))

### generate randomized p-values ###
n.rep <- 10
n.ecdf <- n.rep*n.test
randp.ecdf <- rep(0,n.ecdf)
for(i in 0:(n.rep-1)){
  randp.ecdf[(i*n.test+1):((i+1)*n.test)] <- runif(n.test,p.next,p.org)
  print(i)  
}
leng = length(randp.ecdf)

###### q-value method #############
MT <- qvalue(p.org)
QVALUE <- MT$qvalues

######
pi0 = MT$pi0
pi1 <- 1-pi0
#############

for(i.simu in 1:length(alpha.all)){
  alpha = alpha.all[i.simu]
  ####### MCF method ############
  a = 0.00000001
  b = 1
  l = 0
  R = 0
  pFDR = 0
  while( abs(pFDR-alpha)>0.0000001 ){
    l = (a+b)/2
    cdf <- sum(randp.ecdf<l)/leng
    pFDR = (1-pi1)*l/cdf  
    R = round(n.test*cdf) 
    if(pFDR > alpha){
      b=l
    }
    if(pFDR < alpha){
      a=l
    }  
    print(c(a,b,l,pFDR))
  }
  power.MCF[i.simu] = R
  ###### q-value method #############
  power.q[i.simu] = sum(QVALUE<alpha)
  
  print(i.simu)
}


delta <- (ALL[,4]+ALL[,6])/(ALL[,5]+ALL[,7]) - (ALL[,8]+ALL[,10])/(ALL[,9]+ALL[,11])
############# alpha = 0.1 ####################
idx.mcf <- which(ALL[,14]==1)
idx.q <- which(ALL[,15]==1)
idx.dif <- setdiff(idx.mcf,idx.q)
hist(delta[idx.dif])
dat <- data.frame(xx = c(delta[idx.q],delta[idx.dif]),yy = c(rep(0,length(idx.q)),rep(1,length(idx.dif))))
p.delta.31 <- ggplot(dat,aes(x=xx)) + 
  geom_histogram(data=subset(dat,yy == '0'),fill = "red", alpha = 0.2, binwidth=0.05) +
  geom_histogram(data=subset(dat,yy == '1'),fill = "blue", alpha = 0.2, binwidth=0.05) +
  xlab("Methylation level difference") +
  ylab("CD8_ES")
p.delta.31
############# alpha = 0.05 ####################
idx.mcf <- which(ALL[,16]==1)
idx.q <- which(ALL[,17]==1)
idx.dif <- setdiff(idx.mcf,idx.q)
dat <- data.frame(xx = c(delta[idx.q],delta[idx.dif]),yy = c(rep(0,length(idx.q)),rep(1,length(idx.dif))))
p.delta.32 <- ggplot(dat,aes(x=xx)) + 
  geom_histogram(data=subset(dat,yy == '0'),fill = "red", alpha = 0.2, binwidth=0.05) +
  geom_histogram(data=subset(dat,yy == '1'),fill = "blue", alpha = 0.2, binwidth=0.05) +
  xlab("Methylation level difference") +
  ylab(" ")
p.delta.32

n = length(alpha.all)
power = data.frame(rep(alpha.all,2), c(power.MCF,power.q), 
                   c(rep('MCF',n),rep('q-value',n)))
colnames(power) = c('alpha', 'power', 'method')

p3 <- ggplot(power, aes(x = alpha, y = power, linetype=method)) + 
  geom_line() + 
  xlab("nominal FDR level") + 
  ylab("number of significant CpGs")
p3



############### multiple plot ###################################
library(gridExtra)
require(grid)

grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}


pdf("/Users/xiaoyudai/Documents/Paper/Real/power.pdf",width=8,height=12)

grid_arrange_shared_legend(p1,p2,p3)

dev.off()


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



pdf("/Users/xiaoyudai/Documents/Paper/Real/delta.pdf",width=8,height=8)

multiplot(p.delta.11,p.delta.21,p.delta.31, p.delta.12,p.delta.22,p.delta.32, cols=2)

dev.off()










