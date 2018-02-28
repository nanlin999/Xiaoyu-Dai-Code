### data set-up ###############
lambda1 <- 15 # or 15

### fixed across repeated simulations ########################
n.test = 50000
N <- rpois(n.test, lambda1) + 1 # total of the binomial
P.all <- c(rep(0.5,n.test*0.9), 0.5 + sample(c(-1,1),n.test*0.1,replace = TRUE)*runif(n.test*0.1,min=0.2,max=0.5))

p.org <- rep(0,n.test) 
p.next <- rep(0,n.test)
p.min <- rep(0,n.test)
X <- rep(0,n.test)

### get pvalues ####################
for(j in 1:n.test){
  n <- N[j]
  a <- rbinom(1,n,P.all[j])
  X[j] = a
  p.org[j] <- binom.test(a,n,p=0.5,alternative="two.sided",conf.level = 0.95)$p.value
  if(a>0&&a<n){
    p.next.temp1 = binom.test(a-1,n,p=0.5,alternative="two.sided",conf.level = 0.95)$p.value
    p.next.temp2 = binom.test(a+1,n,p=0.5,alternative="two.sided",conf.level = 0.95)$p.value
    if(p.next.temp1<p.org[j] & p.next.temp2<p.org[j]){
      p.next[j] = max(p.next.temp1, p.next.temp2)
    } else {
      p.next[j] = min(p.next.temp1, p.next.temp2)
    }
  }else{
    p.next[j] <- 0
  }
  p.min[j] <- binom.test(0,n,p=0.5,alternative="two.sided",conf.level = 0.95)$p.value
}
p.org[p.org>1] <- 1 #remove rounding mistake
p.next[p.next>1] <- 1
p.min[p.min>1] <- 1

### our MCF-based method #######################
l = 0.1
# R <- round(n.test*cdf)
r <- rep(0,n.test)
for(s in 1:n.test){
  if(l>p.org[s]){
    r[s] <- 1}else if(l<p.next[s]){
      r[s] <- 0}else{
        r[s] <- (l-p.next[s])/(p.org[s]-p.next[s])
      }
}


dt = data.frame(mcf = r,
                cdf = c(rep('G0^m', n.test*0.9),
                         rep('G1^m', n.test*0.1)))

plt1 <- ggplot(dt, aes(mcf, linetype = cdf)) + 
    stat_ecdf() + 
    ylab("Binomial Test")
plt1


r.smooth = r
for(i in 1:length(r.smooth)) {
  if(r.smooth[i]!=0 & r.smooth[i]!=1) {
    r.smooth[i] = runif(1, min=0, max=1)
  }
}

dt = data.frame(mcf = r.smooth,
                cdf = c(rep('G0^* ', length(r.smooth)*0.9),
                         rep('G1^* ', length(r.smooth)*0.1)))

plt2 <- ggplot(dt, aes(mcf, linetype = cdf)) + 
    stat_ecdf() + 
    ylab("")
plt2
  


### data set-up ###############
lambda1 <- 20 # or 20
lambda2 <- 20 # or 20

n.test = 50000
N1 <- rpois(n.test, lambda1) # total of each binomial
N2 <- rpois(n.test, lambda2)
p1.Null <- runif(n.test*0.9, min=0.1,max=0.9)
p2.Null <- p1.Null
p1.Nonnull <- 0.5 * runif(n.test*0.1,min=0.1,max=0.9) 
p2.Nonnull <- p1.Nonnull + runif(n.test*0.1,min=0.2,max=0.5)
p1 <- c(p1.Null, p1.Nonnull)
p2 <- c(p2.Null, p2.Nonnull)

p.org <- rep(0,n.test) 
p.next <- rep(0,n.test)
p.min <- rep(0,n.test)
Z1 <- rep(0,n.test)
Z2 <- rep(0,n.test)
Z3 <- rep(0,n.test)
Z4 <- rep(0,n.test)

### get pvalues ####################
for(j in 1:n.test){
  n1 = N1[j]
  n2 = N2[j]
  a <- rbinom(1,n1,p1[j])
  b <- rbinom(1,n2,p2[j])
  Z1[j] = a
  Z2[j] = n1-a
  Z3[j] = b
  Z4[j] = n2-b
  l <- max(0,a+b-n2)
  u <- min(a+b,n1)
  prob.all <- dhyper(l:u,n1,n2,(a+b))
  prob.obs <- dhyper(a,n1,n2,(a+b))
  p.org[j] <- sum(prob.all[which(prob.all<=prob.obs)])
  p.next[j] <- sum(prob.all[which(prob.all<prob.obs)])
  p.min[j] <- min(prob.all)
}
p.org[p.org>1] <- 1 #remove rounding mistake
p.next[p.next>1] <- 1
p.min[p.min>1] <- 1

### our MCF-based method #######################
l = 0.1
# R <- round(n.test*cdf)
r <- rep(0,n.test)
for(s in 1:n.test){
  if(l>p.org[s]){
    r[s] <- 1}else if(l<p.next[s]){
      r[s] <- 0}else{
        r[s] <- (l-p.next[s])/(p.org[s]-p.next[s])
      }
}

dt = data.frame(mcf = r,
                cdf = c(rep('G0^m', n.test*0.9),
                        rep('G1^m', n.test*0.1)))

plt3 <- ggplot(dt, aes(mcf, linetype = cdf)) + 
  stat_ecdf() + 
  ylab("Fisher's Exact Test") 
plt3


r.smooth = r
for(i in 1:length(r.smooth)) {
  if(r.smooth[i]!=0 & r.smooth[i]!=1) {
    r.smooth[i] = runif(1, min=0, max=1)
  }
}

dt = data.frame(mcf = r.smooth,
                cdf = c(rep('G0^* ', length(r.smooth)*0.9),
                        rep('G1^* ', length(r.smooth)*0.1)))

plt4 <- ggplot(dt, aes(mcf, linetype = cdf)) + 
  stat_ecdf() + 
  ylab("")
plt4

pdf("/Users/xiaoyudai/Documents/Paper/Simu_rep/Condition_check.pdf",width=12,height=8)
multiplot(plt1, plt3, plt2, plt4, cols=2)
dev.off()









# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
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

