ALL <- read.table("/Users/xiaoyudai/Documents/multiple-testing/Comparison/CD4_CD8/ALL.txt", header=T)
ALL1 <- read.table("/Users/xiaoyudai/Documents/multiple-testing/Comparison/CD4_ES/ALL.txt", header=T)
ALL2 <- read.table("/Users/xiaoyudai/Documents/multiple-testing/Comparison/CD8_ES/ALL.txt", header=T)

MCF = ALL[,14]
Q = ALL[,15]
rej.mcf = which(MCF==1)
rej.q = which(Q==1)
extra = setdiff(rej.mcf,rej.q)
both.p.org = ALL[rej.q,12]
both.p.next = ALL[rej.q,13]
both.dif = both.p.org - both.p.next
extra.p.org = ALL[extra,12]
extra.p.next = ALL[extra,13]
extra.dif = extra.p.org - extra.p.next
dat1 <- data.frame(xx = c(both.dif,extra.dif),yy = c(rep(0,length(both.dif)), rep(1,length(extra.dif))))
p1 <- ggplot(dat1,aes(x=xx)) + 
        geom_histogram(data=subset(dat,yy == 0),fill = "red", alpha = 0.2) +
        geom_histogram(data=subset(dat,yy == 1),fill = "blue", alpha = 0.4) + 
        xlab(expression(p- p^'-')) +
        ylab("CD4_CD8") +
        ggtitle(expression(paste('nominal FDR level ',alpha,' = 0.1')))


MCF = ALL1[,14]
Q = ALL1[,15]
rej.mcf = which(MCF==1)
rej.q = which(Q==1)
extra = setdiff(rej.mcf,rej.q)
both.p.org = ALL1[rej.q,12]
both.p.next = ALL1[rej.q,13]
both.dif = both.p.org - both.p.next
extra.p.org = ALL1[extra,12]
extra.p.next = ALL1[extra,13]
extra.dif = extra.p.org - extra.p.next
dat2 <- data.frame(xx = c(both.dif,extra.dif),yy = c(rep(0,length(both.dif)), rep(1,length(extra.dif))))
p2 <- ggplot(dat2,aes(x=xx)) + 
  geom_histogram(data=subset(dat,yy == 0),fill = "red", alpha = 0.2) +
  geom_histogram(data=subset(dat,yy == 1),fill = "blue", alpha = 0.4) + 
  coord_cartesian(ylim = c(0, 1000000)) + 
  xlab(expression(p- p^'-')) +
  ylab("CD4_ES")


MCF = ALL2[,14]
Q = ALL2[,15]
rej.mcf = which(MCF==1)
rej.q = which(Q==1)
extra = setdiff(rej.mcf,rej.q)
both.p.org = ALL2[rej.q,12]
both.p.next = ALL2[rej.q,13]
both.dif = both.p.org - both.p.next
extra.p.org = ALL2[extra,12]
extra.p.next = ALL2[extra,13]
extra.dif = extra.p.org - extra.p.next
dat3 <- data.frame(xx = c(both.dif,extra.dif),yy = c(rep(0,length(both.dif)), rep(1,length(extra.dif))))
p3 <- ggplot(dat3,aes(x=xx)) + 
  geom_histogram(data=subset(dat,yy == 0),fill = "red", alpha = 0.2) +
  geom_histogram(data=subset(dat,yy == 1),fill = "blue", alpha = 0.4) + 
  coord_cartesian(ylim = c(0, 1000000)) + 
  xlab(expression(p- p^'-')) +
  ylab("CD8_ES")



MCF = ALL[,16]
Q = ALL[,17]
rej.mcf = which(MCF==1)
rej.q = which(Q==1)
extra = setdiff(rej.mcf,rej.q)
both.p.org = ALL[rej.q,12]
both.p.next = ALL[rej.q,13]
both.dif = both.p.org - both.p.next
extra.p.org = ALL[extra,12]
extra.p.next = ALL[extra,13]
extra.dif = extra.p.org - extra.p.next
dat4 <- data.frame(xx = c(both.dif,extra.dif),yy = c(rep(0,length(both.dif)), rep(1,length(extra.dif))))
p4 <- ggplot(dat4,aes(x=xx)) + 
  geom_histogram(data=subset(dat,yy == 0),fill = "red", alpha = 0.2) +
  geom_histogram(data=subset(dat,yy == 1),fill = "blue", alpha = 0.4) + 
  xlab(expression(p- p^'-')) +
  ylab(" ") +
  ggtitle(expression(paste('nominal FDR level ',alpha,' = 0.05')))


MCF = ALL1[,16]
Q = ALL1[,17]
rej.mcf = which(MCF==1)
rej.q = which(Q==1)
extra = setdiff(rej.mcf,rej.q)
both.p.org = ALL1[rej.q,12]
both.p.next = ALL1[rej.q,13]
both.dif = both.p.org - both.p.next
extra.p.org = ALL1[extra,12]
extra.p.next = ALL1[extra,13]
extra.dif = extra.p.org - extra.p.next
dat5 <- data.frame(xx = c(both.dif,extra.dif),yy = c(rep(0,length(both.dif)), rep(1,length(extra.dif))))
p5 <- ggplot(dat5,aes(x=xx)) + 
  geom_histogram(data=subset(dat,yy == 0),fill = "red", alpha = 0.2) +
  geom_histogram(data=subset(dat,yy == 1),fill = "blue", alpha = 0.4) + 
  coord_cartesian(ylim = c(0, 1000000)) + 
  xlab(expression(p- p^'-')) +
  ylab(" ")


MCF = ALL2[,16]
Q = ALL2[,17]
rej.mcf = which(MCF==1)
rej.q = which(Q==1)
extra = setdiff(rej.mcf,rej.q)
both.p.org = ALL2[rej.q,12]
both.p.next = ALL2[rej.q,13]
both.dif = both.p.org - both.p.next
extra.p.org = ALL2[extra,12]
extra.p.next = ALL2[extra,13]
extra.dif = extra.p.org - extra.p.next
dat6 <- data.frame(xx = c(both.dif,extra.dif),yy = c(rep(0,length(both.dif)), rep(1,length(extra.dif))))
p6 <- ggplot(dat6,aes(x=xx)) + 
  geom_histogram(data=subset(dat,yy == 0),fill = "red", alpha = 0.2) +
  geom_histogram(data=subset(dat,yy == 1),fill = "blue", alpha = 0.4) + 
  coord_cartesian(ylim = c(0, 1000000)) + 
  xlab(expression(p- p^'-')) +
  ylab(" ")




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

pdf("/Users/xiaoyudai/Documents/Paper/Real/pvalueDif.pdf",width=8,height=8)

# multiplot(p1,p3,p5,p2,p4,p6, cols=2)
multiplot(p1,p2,p3,p4,p5,p6, cols = 2)

dev.off()










