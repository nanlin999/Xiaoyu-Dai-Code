mse.setup1.20 = c()

for(storey.lambda in seq(0,1,0.05)){
  
  ### data set-up ###############
  lambda1 <- 20 # or 20
  lambda2 <- 20 # or 20
  
  ### fixed across repeated simulations ########################
  n.test = 5000
  N1 <- rpois(n.test, lambda1) # total of each binomial
  N2 <- rpois(n.test, lambda2)
  p1.Null <- runif(4500, min=0,max=1)
  p2.Null <- p1.Null
  p1.Nonnull <- 0.5 * runif(500,min=0,max=1) 
  p2.Nonnull <- p1.Nonnull + runif(500,min=0.2,max=0.5)
  p1 <- c(p1.Null, p1.Nonnull)
  p2 <- c(p2.Null, p2.Nonnull)
  
  ############
  N1.temp = N1
  N2.temp = N2
  p1.temp = p1
  p2.temp = p2
  #####################
  
  ### repeated simulations ##################
  alpha.all = seq(0.001,0.2,by=0.002)
  n.alpha = length(alpha.all)
  n.rep.simu = 100
  
  pi0.all = c()
  
  for(i.simu in 1:n.rep.simu){
    ### generate each binomial ##############################
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
    
    ### get samples of randomized pvalues #########################
    n.rep <- 1000
    n.ecdf <- n.rep*n.test
    randp.ecdf <- rep(0,n.ecdf)
    sum.pi0 = 0
    for(i in 0:(n.rep-1)){
      randp = runif(n.test,p.next,p.org)
      sum.pi0 = sum.pi0 + (sum(randp > storey.lambda)) / ((1-storey.lambda)*5000)
      randp.ecdf[(i*n.test+1):((i+1)*n.test)] <- randp
      # print(i)
    }
    
    ### calculate pi0 ##################
    pi0 = sum.pi0/n.rep
    pi0.all = c(pi0.all, pi0)
  }
  
  mse = mean((pi0.all - 0.9)^2)
  print(mse)
  mse.setup1.20 = c(mse.setup1.20, mse)
}

mse.setup2.15 = c()

for(storey.lambda in seq(0,1,0.05)){
  
  lambda1 <- 15 # or 15
  
  ### fixed across repeated simulations ########################
  n.test = 5000
  N <- rpois(n.test, lambda1) + 1 # total of the binomial
  P.all <- c(rep(0.5,4500), 0.5 + sample(c(-1,1),500,replace = TRUE)*runif(500,min=0.2,max=0.5))
  
  #####
  N.temp = N
  P.temp = P.all
  ### repeated simulations ##################
  alpha.all = seq(0.001,0.2,by=0.002)
  n.alpha = length(alpha.all)
  n.rep.simu = 100
  
  pi0.all = c()
  
  for(i.simu in 1:n.rep.simu){
    ### generate each binomial ##############################
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
    
    ### get samples of randomized pvalues #########################
    ### get samples of randomized pvalues #########################
    n.rep <- 1000
    n.ecdf <- n.rep*n.test
    sum.pi0 = 0
    for(i in 0:(n.rep-1)){
      randp = runif(n.test,p.next,p.org)
      # sum.pi0 = sum.pi0 + (sum(randp > storey.lambda)+1) / ((1-storey.lambda)*5000)
      sum.pi0 = sum.pi0 + (sum(randp > storey.lambda)) / ((1-storey.lambda)*5000)
      # print(i)
    }
    
    ### calculate pi0 ##################
    pi0 = sum.pi0/n.rep
    print(pi0)
    pi0.all = c(pi0.all, pi0)
  }
  mse = mean((pi0.all - 0.9)^2)
  print(mse)
  mse.setup2.15 = c(mse.setup2.15, mse)
}

x = rep(seq(0,1,0.05),50)
y = mse.setup1.25

plot(x,y, cex=.1)
fit2<-smooth.spline(x,y)
lines(fit2)


m = matrix(mse.setup1.25, nrow=50, byrow = T)

mm = apply(m,2,mean)

plot(seq(0,1,0.05), mse.setup2.15)
fit2<-smooth.spline(seq(0,1,0.05),mse.setup1.25,cv = TRUE)
lines(fit2)

dt1 = data.frame(x = seq(0,1,0.05), y=mse.setup1.25)
p1 <- ggplot(data=dt1, aes(x=x, y=y)) +
    geom_line() + 
    xlab("t") + 
    ylab("mean squred error") +
    ggtitle(expression(paste('            setup1: ',mu[1],' = ',mu[2],' = 25')))
    

dt2 = data.frame(x = seq(0,1,0.05), y=mse.setup1.20)
p2 <- ggplot(data=dt2, aes(x=x, y=y)) +
  geom_line() + 
  xlab("t") + 
  ylab("mean squred error") +
  ggtitle(expression(paste('            setup1: ',mu[1],' = ',mu[2],' = 20')))

dt3 = data.frame(x = seq(0,1,0.05), y=mse.setup2.15)
p3 <- ggplot(data=dt3, aes(x=x, y=y)) +
  geom_line() + 
  xlab("t") + 
  ylab("mean squred error") +
  ggtitle(expression(paste('            setup2: ',mu,' = 15')))

dt4 = data.frame(x = seq(0,1,0.05), y=mse.setup2.10)
p4 <- ggplot(data=dt3, aes(x=x, y=y)) +
  geom_line() + 
  xlab("t") + 
  ylab("mean squred error") +
  ggtitle(expression(paste('            setup2: ',mu,' = 10')))


pdf("/Users/xiaoyudai/Documents/Paper/Simu_rep/pi/pi09.pdf",width=8,height=6)

multiplot(p1, p2, p3, p4, cols=2)

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
