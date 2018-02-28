library(DSS)
library(ggplot2)

ALL <- read.table("/Users/xiaoyudai/Documents/multiple-testing/Comparison/CD4_CD8/ALL.txt", header=T)
ALL1 <- read.table("/Users/xiaoyudai/Documents/multiple-testing/Comparison/CD4_ES/ALL.txt", header=T)
ALL2 <- read.table("/Users/xiaoyudai/Documents/multiple-testing/Comparison/CD8_ES/ALL.txt", header=T)



z.prop = function(x1,x2,n1,n2){
  numerator = (x1/n1) - (x2/n2)
  p.common = (x1+x2) / (n1+n2)
  denominator = sqrt(p.common * (1-p.common) * (1/n1 + 1/n2))
  z.prop.ris = numerator / denominator
  return(z.prop.ris)
}

f <- function(ALL, alpha, text_ylab = ' ', text_title = ' '){
  if(alpha == 0.05){
    rej.mcf = which(ALL$MCF005==1)
    rej.q = which(ALL$Qvalue005==1)
  } else {
    rej.mcf = which(ALL$MCF010==1)
    rej.q = which(ALL$Qvalue010==1)
  }
  extra = setdiff(rej.mcf, rej.q)
  both = rej.q
  
  X1.extra = ALL[extra,4] + ALL[extra,6]
  N1.extra = ALL[extra,5] + ALL[extra,7] 
  X2.extra = ALL[extra,8] + ALL[extra,10]
  N2.extra = ALL[extra,9] + ALL[extra,11]
  
  X1.both = ALL[both,4] + ALL[both,6]
  N1.both = ALL[both,5] + ALL[both,7] 
  X2.both = ALL[both,8] + ALL[both,10]
  N2.both = ALL[both,9] + ALL[both,11]
  
  Z.extra = z.prop(X1.extra, X2.extra, N1.extra, N2.extra)
  Z.both = z.prop(X1.both, X2.both, N1.both, N2.both)
  
  
  
 
  dat <- data.frame(xx = c(Z.both,Z.extra),yy = c(rep(0,length(Z.both)),rep(1,length(Z.extra))))
  p <- ggplot(dat,aes(x=xx)) + 
    geom_histogram(data=subset(dat,yy == '0'),fill = "red", alpha = 0.2, binwidth=0.05) +
    geom_histogram(data=subset(dat,yy == '1'),fill = "blue", alpha = 0.4, binwidth=0.05) +
    coord_cartesian(xlim = c(-10, 10)) +
    xlab("Methylation level difference") +
    ylab(text_ylab) + 
    ggtitle(text_title)
  
  return(p)
}


p1 = f(ALL, alpha = 0.05, text_ylab = 'CD4_CD8', text_title = expression(paste('nominal FDR level ',alpha,' = 0.05')))
p4 = f(ALL, alpha = 0.1, text_title = expression(paste('nominal FDR level ',alpha,' = 0.1')))
p2 = f(ALL1, alpha = 0.05, text_ylab = 'CD4_ES')
p5 = f(ALL1, alpha = 0.1)
p3 = f(ALL2, alpha = 0.05, text_ylab = 'CD8_ES')
p6 = f(ALL2, alpha = 0.1)


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

pdf("~/Desktop/trial.pdf",width=8,height=8)

# multiplot(p1,p3,p5,p2,p4,p6, cols=2)
multiplot(p1,p2,p3,p4,p5,p6, cols = 2)

dev.off()






