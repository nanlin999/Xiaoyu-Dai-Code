library(fdrDiscreteNull)
library(cp4p)
library(ggplot2)

### data set-up ###############
lambda1 <- 10 # or 15

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

result = matrix(nrow = n.alpha*n.rep.simu, ncol = 15)
colnames(result) = c('FDR.q', 'FNR.q', 'Power.q',
                     'FDR.mcf', 'FNR.mcf', 'Power.mcf',
                     'FDR.Gilbert', 'FNR.Gilbert', 'Power.Gilbert',
                     'FDR.Chen', 'FNR.Chen', 'Power.Chen',
                     'FDR.mid', 'FNR.mid', 'Power.mid')


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
  randp.ecdf <- rep(0,n.ecdf)
  sum.pi0 = 0
  for(i in 0:(n.rep-1)){
    randp = runif(n.test,p.next,p.org)
    storey.lambda = 0.5
    sum.pi0 = sum.pi0 + (sum(randp > storey.lambda)+1) / ((1-storey.lambda)*5000)
    randp.ecdf[(i*n.test+1):((i+1)*n.test)] <- randp
    # print(i)
  }
  
  ### calculate pi0 ##################
  pi0 = sum.pi0/n.rep
  print(paste('pi0 = ', pi0))
  
  ### calculate pFDR.all in order to get lambda.star ###################
  lambda.all <- c(seq(0.00001,0.01,by=0.00003), seq(0.01,0.03,by=0.0001))
  pFDR.all <- rep(0,length(lambda.all))
  for(i in 1:length(lambda.all)){
    l = lambda.all[i]
    cdf = sum(randp.ecdf<l)/length(randp.ecdf)
    pFDR.all[i] = pi0*l/cdf
    # print(i)
  }
  
  ### dataframe for Chen and Doerge
  tmp.Chen.BT = as.matrix(cbind(X,N-X))
  ### one sample of mid-pvalues ###################
  p.mid = (p.org+p.next)/2
  ### one sample of randomized pvalues #################
  p.rand = runif(n.test, min=p.next, max=p.org)
  
  # lambda.all.q <- c(seq(0.00001,0.01,by=0.00003), seq(0.01,0.03,by=0.0001))
  # pFDR.all.q <- rep(0,length(lambda.all.q))
  # pi0.q = (sum(p.rand>0.5)+1)/(5000*0.5)
  # for(i in 1:length(lambda.all)){
  #   l = lambda.all.q[i]
  #   cdf = sum(p.rand<l)/length(p.rand)
  #   pFDR.all[i] = pi0.q*l/cdf
  #   # print(i)
  # }
  
  
  ### Comparing differnt multiple testing methods ########################## 
  for(i.alpha in 1:n.alpha){
    # id in the result
    id = (i.simu-1) * n.alpha + i.alpha
    
    ### get lambda.star #############
    alpha = alpha.all[i.alpha]
    
    if(pFDR.all[1]<alpha){
      idx.lambda.star = max(which(pFDR.all<alpha))
    }else{
      idx.lambda.star = 1
    }
    lambda.star = lambda.all[idx.lambda.star]
    
    # if(pFDR.all.q[1]<alpha){
    #   idx.lambda.star.q = max(which(pFDR.all.q<alpha))
    # }else{
    #   idx.lambda.star.q = 1
    # }
    # lambda.star.q = lambda.all.q[idx.lambda.star.q]
    
    ### q-value method ##########################
    rej.q <- which(p.rand<lambda.star)
    # FDR.q[i.simu,i.alpha] <- sum(rej.q<4501)/length(rej.q)
    # FNR.q[i.simu,i.alpha] = sum((1:5000)[-rej.q]>4500)/(5000-length(rej.q))
    result[id,'FDR.q'] = sum(rej.q<4501)/length(rej.q)
    result[id,'FNR.q'] = sum((1:5000)[-rej.q]>4500)/(5000-length(rej.q))
    result[id,'Power.q'] = sum(rej.q>4500)/500
    
    
    
    ### our MCF-based method #######################
    l = lambda.star
    cdf = sum(randp.ecdf<l)/length(randp.ecdf)
    # R <- round(n.test*cdf)
    r <- rep(0,n.test)
    for(s in 1:n.test){
      if(l>p.org[s]){
        r[s] <- 1}else if(l<p.next[s]){
          r[s] <- 0}else{
            r[s] <- (l-p.next[s])/(p.org[s]-p.next[s])
          }
    }
    
    rej.mcf <- which(r >= quantile(r, probs = 1-cdf, type = 1))
    # rej.mcf <- sort(r, decreasing = F, index.return = T)$ix[(n.test-R+1):n.test]
    # rej.mcf <- which(r>=sort(r)[n.test-R])
    # FDR.mcf[i.simu,i.alpha] <- sum(rej.mcf<4501)/length(rej.mcf)
    # FNR.mcf[i.simu,i.alpha] = sum((1:5000)[-rej.mcf]>4500)/(5000-length(rej.mcf))
    result[id,'FDR.mcf'] = sum(rej.mcf<4501)/length(rej.mcf)
    result[id,'FNR.mcf'] = sum((1:5000)[-rej.mcf]>4500)/(5000-length(rej.mcf))
    result[id,'Power.mcf'] = sum(rej.mcf>4500)/500
    
    #### Gilbert's method #######################
    k <- 1
    m <- 1
    while(k<=n.test){
      m <- sum(p.min < alpha/k)
      if(m<=k) break
      k = k+1
    }
    R <- which(p.min < alpha/k)
    temp.p <- p.org[R]
    indx <- (1:n.test)[R]
    temp.indx <- max(which(sort(temp.p)<(alpha*seq(1:m)/m)))
    temp.rej <- which(temp.p<sort(temp.p)[temp.indx+1])
    rej.Gilbert <- indx[temp.rej]
    # FDR.Gilbert[i.simu,i.alpha] = sum(rej.Gilbert<4501)/length(rej.Gilbert)
    # FNR.Gilbert[i.simu,i.alpha] = sum((1:5000)[-rej.Gilbert]>4500)/(5000-length(rej.Gilbert))
    result[id,'FDR.Gilbert'] = sum(rej.Gilbert<4501)/length(rej.Gilbert)
    result[id,'FNR.Gilbert'] = sum((1:5000)[-rej.Gilbert]>4500)/(5000-length(rej.Gilbert))
    result[id,'Power.Gilbert'] = sum(rej.Gilbert>4500)/500
    
    
    ### Chen and Doerge's method ##########################
    mt.Chen = GeneralizedEstimatorsGrouped(data=tmp.Chen.BT,
                                           test_in = "Binomial Test",FET_via_in = "IndividualMarginals",
                                           grpby = 'quantileOfRowTotal', ngrp_in = 3,
                                           FDRlevel_in = alpha) # , lambda_in = 0.5, epsilon_in = 1)
    rej.Chen = mt.Chen$aBH$IndicesOfDiscoveries

    # FDR.Chen[i.simu,i.alpha] = sum(rej.Chen<4501)/length(rej.Chen)
    # FNR.Chen[i.simu,i.alpha] = sum((1:5000)[-rej.Chen]>4500)/(5000-length(rej.Chen))
    result[id,'FDR.Chen'] = sum(rej.Chen<4501)/length(rej.Chen)
    result[id,'FNR.Chen'] = sum((1:5000)[-rej.Chen]>4500)/(5000-length(rej.Chen))
    result[id,'Power.Chen'] = sum(rej.Chen>4500)/500
    
    
    ### BH algorithm on mid-pvalues
    idx.mid = max(which(sort(p.mid)<(alpha*seq(1:n.test)/n.test)))
    rej.mid = which(p.mid<sort(p.mid)[idx.mid+1])
    # FDR.mid[i.simu,i.alpha] = sum(rej.mid<4501)/length(rej.mid)
    # FNR.mid[i.simu,i.alpha] = sum((1:5000)[-rej.mid]>4500)/(5000-length(rej.mid))
    result[id,'FDR.mid'] = sum(rej.mid<4501)/length(rej.mid)
    result[id,'FNR.mid'] = sum((1:5000)[-rej.mid]>4500)/(5000-length(rej.mid))
    result[id,'Power.mid'] = sum(rej.mid>4500)/500
    
    # 
    print(c(i.simu, i.alpha))
  }
  
  write.csv(result, file = '~/Documents/Paper/Simu_rep/BT/result.csv', quote = F, row.names = F)
}



FDR.mcf = matrix(rep(0,n.alpha*n.rep.simu), ncol=n.alpha)
FDR.q = matrix(rep(0,n.alpha*n.rep.simu), ncol=n.alpha)
FDR.Gilbert = matrix(rep(0,n.alpha*n.rep.simu), ncol=n.alpha)
FDR.Chen = matrix(rep(0,n.alpha*n.rep.simu), ncol=n.alpha)
FDR.mid = matrix(rep(0,n.alpha*n.rep.simu), ncol=n.alpha)
Power.mcf = matrix(rep(0,n.alpha*n.rep.simu), ncol=n.alpha)
Power.q = matrix(rep(0,n.alpha*n.rep.simu), ncol=n.alpha)
Power.Gilbert = matrix(rep(0,n.alpha*n.rep.simu), ncol=n.alpha)
Power.Chen = matrix(rep(0,n.alpha*n.rep.simu), ncol=n.alpha)
Power.mid = matrix(rep(0,n.alpha*n.rep.simu), ncol=n.alpha)
FNR.mcf = matrix(rep(0,n.alpha*n.rep.simu), ncol=n.alpha)
FNR.q = matrix(rep(0,n.alpha*n.rep.simu), ncol=n.alpha)
FNR.Gilbert = matrix(rep(0,n.alpha*n.rep.simu), ncol=n.alpha)
FNR.Chen = matrix(rep(0,n.alpha*n.rep.simu), ncol=n.alpha)
FNR.mid = matrix(rep(0,n.alpha*n.rep.simu), ncol=n.alpha)

for(i.simu in 1:n.rep.simu){
  for(i.alpha in 1:n.alpha){
    FDR.mcf[i.simu,i.alpha] = result[(i.simu-1) * n.alpha + i.alpha,'FDR.mcf']
    FDR.q[i.simu,i.alpha] = result[(i.simu-1) * n.alpha + i.alpha,'FDR.q']
    FDR.Gilbert[i.simu,i.alpha] = result[(i.simu-1) * n.alpha + i.alpha,'FDR.Gilbert']
    FDR.Chen[i.simu,i.alpha] = result[(i.simu-1) * n.alpha + i.alpha,'FDR.Chen']
    FDR.mid[i.simu,i.alpha] = result[(i.simu-1) * n.alpha + i.alpha,'FDR.mid']
    Power.mcf[i.simu,i.alpha] = result[(i.simu-1) * n.alpha + i.alpha,'Power.mcf']
    Power.q[i.simu,i.alpha] = result[(i.simu-1) * n.alpha + i.alpha,'Power.q']
    Power.Gilbert[i.simu,i.alpha] = result[(i.simu-1) * n.alpha + i.alpha,'Power.Gilbert']
    Power.Chen[i.simu,i.alpha] = result[(i.simu-1) * n.alpha + i.alpha,'Power.Chen']
    Power.mid[i.simu,i.alpha] = result[(i.simu-1) * n.alpha + i.alpha,'Power.mid']
    FNR.mcf[i.simu,i.alpha] = result[(i.simu-1) * n.alpha + i.alpha,'FNR.mcf']
    FNR.q[i.simu,i.alpha] = result[(i.simu-1) * n.alpha + i.alpha,'FNR.q']
    FNR.Gilbert[i.simu,i.alpha] = result[(i.simu-1) * n.alpha + i.alpha,'FNR.Gilbert']
    FNR.Chen[i.simu,i.alpha] = result[(i.simu-1) * n.alpha + i.alpha,'FNR.Chen']
    FNR.mid[i.simu,i.alpha] = result[(i.simu-1) * n.alpha + i.alpha,'FNR.mid']
  }
}


FDR.mcf.mean.1 = apply(FDR.mcf,2,mean,na.rm=TRUE)
FDR.q.mean.1 = apply(FDR.q,2,mean,na.rm=TRUE)
FDR.mid.mean.1 = apply(FDR.mid,2,mean,na.rm=TRUE)
FDR.Chen.mean.1 = apply(FDR.Chen,2,mean,na.rm=TRUE)
FDR.Gilbert.mean.1 = apply(FDR.Gilbert,2,mean,na.rm=TRUE)
Power.mcf.mean.1 = apply(Power.mcf,2,mean,na.rm=TRUE)
Power.q.mean.1 = apply(Power.q,2,mean,na.rm=TRUE)
Power.mid.mean.1 = apply(Power.mid,2,mean,na.rm=TRUE)
Power.Chen.mean.1 = apply(Power.Chen,2,mean,na.rm=TRUE)
Power.Gilbert.mean.1 = apply(Power.Gilbert,2,mean,na.rm=TRUE)
FNR.mcf.mean.1 = apply(FNR.mcf,2,mean,na.rm=TRUE)
FNR.q.mean.1 = apply(FNR.q,2,mean,na.rm=TRUE)
FNR.mid.mean.1 = apply(FNR.mid,2,mean,na.rm=TRUE)
FNR.Chen.mean.1 = apply(FNR.Chen,2,mean,na.rm=TRUE)
FNR.Gilbert.mean.1 = apply(FNR.Gilbert,2,mean,na.rm=TRUE)

FDR.mcf.sd.1 = apply(FDR.mcf,2,sd,na.rm=TRUE)
FDR.q.sd.1 = apply(FDR.q,2,sd,na.rm=TRUE)
FDR.mid.sd.1 = apply(FDR.mid,2,sd,na.rm=TRUE)
FDR.Chen.sd.1 = apply(FDR.Chen,2,sd,na.rm=TRUE)
FDR.Gilbert.sd.1 = apply(FDR.Gilbert,2,sd,na.rm=TRUE)


####################################################################
####################################################################
### data set-up ###############
lambda1 <- 15 # or 10

### fixed across repeated simulations ########################
n.test = 5000
N <- rpois(n.test, lambda1) # total of the binomial
P.all <- c(rep(0.5,4500), 0.5 + sample(c(-1,1),500,replace = TRUE)*runif(500,min=0.3,max=0.5))

### repeated simulations ##################
alpha.all = seq(0.001,0.2,by=0.002)
n.alpha = length(alpha.all)
n.rep.simu = 100

result2 = matrix(nrow = n.alpha*n.rep.simu, ncol = 15)
colnames(result2) = c('FDR.q', 'FNR.q', 'Power.q',
                      'FDR.mcf', 'FNR.mcf', 'Power.mcf',
                      'FDR.Gilbert', 'FNR.Gilbert', 'Power.Gilbert',
                      'FDR.Chen', 'FNR.Chen', 'Power.Chen',
                      'FDR.mid', 'FNR.mid', 'Power.mid')


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
  n.rep <- 1000
  n.ecdf <- n.rep*n.test
  randp.ecdf <- rep(0,n.ecdf)
  sum.pi0 = 0
  for(i in 0:(n.rep-1)){
    randp = runif(n.test,p.next,p.org)
    storey.lambda = 0.5
    sum.pi0 = sum.pi0 + (sum(randp > storey.lambda)) / ((1-storey.lambda)*5000)
    randp.ecdf[(i*n.test+1):((i+1)*n.test)] <- randp
    # print(i)
  }
  
  ### calculate pi0 ##################
  pi0 = sum.pi0/n.rep
  print(paste('pi0 = ', pi0))
  
  ### calculate pFDR.all in order to get lambda.star ###################
  lambda.all <- c(seq(0.00001,0.01,by=0.00003), seq(0.01,0.03,by=0.0001))
  pFDR.all <- rep(0,length(lambda.all))
  for(i in 1:length(lambda.all)){
    l = lambda.all[i]
    cdf = sum(randp.ecdf<l)/length(randp.ecdf)
    pFDR.all[i] = pi0*l/cdf
    # print(i)
  }
  
  ### dataframe for Chen and Doerge
  tmp.Chen.BT = as.matrix(cbind(X,N-X))
  ### one sample of mid-pvalues ###################
  p.mid = (p.org+p.next)/2
  ### one sample of randomized pvalues #################
  p.rand = runif(n.test, min=p.next, max=p.org)
  
  ### Comparing differnt multiple testing methods ########################## 
  for(i.alpha in 1:n.alpha){
    # id in the result
    id = (i.simu-1) * n.alpha + i.alpha
    
    ### get lambda.star #############
    alpha = alpha.all[i.alpha]
    if(pFDR.all[1]<alpha){
      idx.lambda.star = max(which(pFDR.all<alpha))
    }else{
      idx.lambda.star = 1
    }
    lambda.star = lambda.all[idx.lambda.star]
    
    ### q-value method ##########################
    rej.q <- which(p.rand<lambda.star)
    # FDR.q[i.simu,i.alpha] <- sum(rej.q<4501)/length(rej.q)
    # FNR.q[i.simu,i.alpha] = sum((1:5000)[-rej.q]>4500)/(5000-length(rej.q))
    result2[id,'FDR.q'] = sum(rej.q<4501)/length(rej.q)
    result2[id,'FNR.q'] = sum((1:5000)[-rej.q]>4500)/(5000-length(rej.q))
    result2[id,'Power.q'] = sum(rej.q>4500)/500
    
    
    
    ### our MCF-based method #######################
    l = lambda.star
    cdf = sum(randp.ecdf<l)/length(randp.ecdf)
    # R <- round(n.test*cdf)
    r <- rep(0,n.test)
    for(s in 1:n.test){
      if(l>p.org[s]){
        r[s] <- 1}else if(l<p.next[s]){
          r[s] <- 0}else{
            r[s] <- (l-p.next[s])/(p.org[s]-p.next[s])
          }
    }
    
    rej.mcf <- which(r >= quantile(r, probs = 1-cdf, type = 1))
    # rej.mcf <- sort(r, decreasing = F, index.return = T)$ix[(n.test-R+1):n.test]
    # FDR.mcf[i.simu,i.alpha] <- sum(rej.mcf<4501)/length(rej.mcf)
    # FNR.mcf[i.simu,i.alpha] = sum((1:5000)[-rej.mcf]>4500)/(5000-length(rej.mcf))
    result2[id,'FDR.mcf'] = sum(rej.mcf<4501)/length(rej.mcf)
    result2[id,'FNR.mcf'] = sum((1:5000)[-rej.mcf]>4500)/(5000-length(rej.mcf))
    result2[id,'Power.mcf'] = sum(rej.mcf>4500)/500
    
    #### Gilbert's method #######################
    k <- 1
    m <- 1
    while(k<=n.test){
      m <- sum(p.min < alpha/k)
      if(m<=k) break
      k = k+1
    }
    R <- which(p.min < alpha/k)
    temp.p <- p.org[R]
    indx <- (1:n.test)[R]
    temp.indx <- max(which(sort(temp.p)<(alpha*seq(1:m)/m)))
    temp.rej <- which(temp.p<sort(temp.p)[temp.indx+1])
    rej.Gilbert <- indx[temp.rej]
    # FDR.Gilbert[i.simu,i.alpha] = sum(rej.Gilbert<4501)/length(rej.Gilbert)
    # FNR.Gilbert[i.simu,i.alpha] = sum((1:5000)[-rej.Gilbert]>4500)/(5000-length(rej.Gilbert))
    result2[id,'FDR.Gilbert'] = sum(rej.Gilbert<4501)/length(rej.Gilbert)
    result2[id,'FNR.Gilbert'] = sum((1:5000)[-rej.Gilbert]>4500)/(5000-length(rej.Gilbert))
    result2[id,'Power.Gilbert'] = sum(rej.Gilbert>4500)/500
    
    
    ### Chen and Doerge's method ##########################
    mt.Chen = GeneralizedEstimatorsGrouped(data=tmp.Chen.BT,
                                       test_in = "Binomial Test", ngrp_in = 3,
                                       FDRlevel_in =alpha, lambda_in = 0.5, epsilon_in =1)
    rej.Chen = mt.Chen$Gen
    # FDR.Chen[i.simu,i.alpha] = sum(rej.Chen<4501)/length(rej.Chen)
    # FNR.Chen[i.simu,i.alpha] = sum((1:5000)[-rej.Chen]>4500)/(5000-length(rej.Chen))
    result2[id,'FDR.Chen'] = sum(rej.Chen<4501)/length(rej.Chen)
    result2[id,'FNR.Chen'] = sum((1:5000)[-rej.Chen]>4500)/(5000-length(rej.Chen))
    result2[id,'Power.Chen'] = sum(rej.Chen>4500)/500
    
    
    ### BH algorithm on mid-pvalues
    idx.mid = max(which(sort(p.mid)<(alpha*seq(1:n.test)/n.test)))
    rej.mid = which(p.mid<sort(p.mid)[idx.mid+1])
    # FDR.mid[i.simu,i.alpha] = sum(rej.mid<4501)/length(rej.mid)
    # FNR.mid[i.simu,i.alpha] = sum((1:5000)[-rej.mid]>4500)/(5000-length(rej.mid))
    result2[id,'FDR.mid'] = sum(rej.mid<4501)/length(rej.mid)
    result2[id,'FNR.mid'] = sum((1:5000)[-rej.mid]>4500)/(5000-length(rej.mid))
    result2[id,'Power.mid'] = sum(rej.mid>4500)/500
    
    # 
    print(c(i.simu, i.alpha))
  }
  
  write.csv(result2, file = '~/Documents/Paper/Simu_rep/BT/result2.csv', quote = F, row.names = F)
}


FDR.mcf = matrix(rep(0,n.alpha*n.rep.simu), ncol=n.alpha)
FDR.q = matrix(rep(0,n.alpha*n.rep.simu), ncol=n.alpha)
FDR.Gilbert = matrix(rep(0,n.alpha*n.rep.simu), ncol=n.alpha)
FDR.Chen = matrix(rep(0,n.alpha*n.rep.simu), ncol=n.alpha)
FDR.mid = matrix(rep(0,n.alpha*n.rep.simu), ncol=n.alpha)
Power.mcf = matrix(rep(0,n.alpha*n.rep.simu), ncol=n.alpha)
Power.q = matrix(rep(0,n.alpha*n.rep.simu), ncol=n.alpha)
Power.Gilbert = matrix(rep(0,n.alpha*n.rep.simu), ncol=n.alpha)
Power.Chen = matrix(rep(0,n.alpha*n.rep.simu), ncol=n.alpha)
Power.mid = matrix(rep(0,n.alpha*n.rep.simu), ncol=n.alpha)
FNR.mcf = matrix(rep(0,n.alpha*n.rep.simu), ncol=n.alpha)
FNR.q = matrix(rep(0,n.alpha*n.rep.simu), ncol=n.alpha)
FNR.Gilbert = matrix(rep(0,n.alpha*n.rep.simu), ncol=n.alpha)
FNR.Chen = matrix(rep(0,n.alpha*n.rep.simu), ncol=n.alpha)
FNR.mid = matrix(rep(0,n.alpha*n.rep.simu), ncol=n.alpha)

for(i.simu in 1:n.rep.simu){
  for(i.alpha in 1:n.alpha){
    FDR.mcf[i.simu,i.alpha] = result2[(i.simu-1) * n.alpha + i.alpha,'FDR.mcf']
    FDR.q[i.simu,i.alpha] = result2[(i.simu-1) * n.alpha + i.alpha,'FDR.q']
    FDR.Gilbert[i.simu,i.alpha] = result2[(i.simu-1) * n.alpha + i.alpha,'FDR.Gilbert']
    FDR.Chen[i.simu,i.alpha] = result2[(i.simu-1) * n.alpha + i.alpha,'FDR.Chen']
    FDR.mid[i.simu,i.alpha] = result2[(i.simu-1) * n.alpha + i.alpha,'FDR.mid']
    Power.mcf[i.simu,i.alpha] = result2[(i.simu-1) * n.alpha + i.alpha,'Power.mcf']
    Power.q[i.simu,i.alpha] = result2[(i.simu-1) * n.alpha + i.alpha,'Power.q']
    Power.Gilbert[i.simu,i.alpha] = result2[(i.simu-1) * n.alpha + i.alpha,'Power.Gilbert']
    Power.Chen[i.simu,i.alpha] = result2[(i.simu-1) * n.alpha + i.alpha,'Power.Chen']
    Power.mid[i.simu,i.alpha] = result2[(i.simu-1) * n.alpha + i.alpha,'Power.mid']
    FNR.mcf[i.simu,i.alpha] = result2[(i.simu-1) * n.alpha + i.alpha,'FNR.mcf']
    FNR.q[i.simu,i.alpha] = result2[(i.simu-1) * n.alpha + i.alpha,'FNR.q']
    FNR.Gilbert[i.simu,i.alpha] = result2[(i.simu-1) * n.alpha + i.alpha,'FNR.Gilbert']
    FNR.Chen[i.simu,i.alpha] = result2[(i.simu-1) * n.alpha + i.alpha,'FNR.Chen']
    FNR.mid[i.simu,i.alpha] = result2[(i.simu-1) * n.alpha + i.alpha,'FNR.mid']
  }
}



FDR.mcf.mean.2 = apply(FDR.mcf,2,mean,na.rm=TRUE)
FDR.q.mean.2 = apply(FDR.q,2,mean,na.rm=TRUE)
FDR.mid.mean.2 = apply(FDR.mid,2,mean,na.rm=TRUE)
FDR.Chen.mean.2 = apply(FDR.Chen,2,mean,na.rm=TRUE)
FDR.Gilbert.mean.2 = apply(FDR.Gilbert,2,mean,na.rm=TRUE)
Power.mcf.mean.2 = apply(Power.mcf,2,mean,na.rm=TRUE)
Power.q.mean.2 = apply(Power.q,2,mean,na.rm=TRUE)
Power.mid.mean.2 = apply(Power.mid,2,mean,na.rm=TRUE)
Power.Chen.mean.2 = apply(Power.Chen,2,mean,na.rm=TRUE)
Power.Gilbert.mean.2 = apply(Power.Gilbert,2,mean,na.rm=TRUE)
FNR.mcf.mean.2 = apply(FNR.mcf,2,mean,na.rm=TRUE)
FNR.q.mean.2 = apply(FNR.q,2,mean,na.rm=TRUE)
FNR.mid.mean.2 = apply(FNR.mid,2,mean,na.rm=TRUE)
FNR.Chen.mean.2 = apply(FNR.Chen,2,mean,na.rm=TRUE)
FNR.Gilbert.mean.2 = apply(FNR.Gilbert,2,mean,na.rm=TRUE)

FDR.mcf.sd.2 = apply(FDR.mcf,2,sd,na.rm=TRUE)
FDR.q.sd.2 = apply(FDR.q,2,sd,na.rm=TRUE)
FDR.mid.sd.2 = apply(FDR.mid,2,sd,na.rm=TRUE)
FDR.Chen.sd.2 = apply(FDR.Chen,2,sd,na.rm=TRUE)
FDR.Gilbert.sd.2 = apply(FDR.Gilbert,2,sd,na.rm=TRUE)




### ploting #####################
### set up ###################
library(ggplot2)
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



### FDR upper bound ####################

################# start #############################

n = length(alpha.all)
FDR = data.frame(rep(alpha.all,5), c(FDR.mcf.mean.1,FDR.q.mean.1,FDR.mid.mean.1,FDR.Chen.mean.1,FDR.Gilbert.mean.1), 
                 c(rep('MCF',n), rep('Habiger',n), rep('mid',n), rep('Chen',n), 
                   rep('Gilbert',n)))
colnames(FDR) = c('alpha', 'fdr', 'method')
Power = data.frame(rep(alpha.all,5), c(Power.mcf.mean.1,Power.q.mean.1,Power.mid.mean.1,Power.Chen.mean.1,Power.Gilbert.mean.1), 
                   c(rep('MCF',n), rep('Habiger',n), rep('mid',n), rep('Chen',n), 
                     rep('Gilbert',n)))
colnames(Power) = c('alpha', 'power', 'method')
pFNR_bound = FNR_bound_bt(lambda1=10)$pFNR_bound
FNR = data.frame(rep(alpha.all,3), c(FNR.mcf.mean.1,FNR.q.mean.1,pFNR_bound), 
                 c(rep('MCF',n), rep('Habiger',n), rep('pFDR_bound',n)))
colnames(FNR) = c('alpha', 'fnr', 'method')
FDR.sd = data.frame(rep(alpha.all,2), c(FDR.mcf.sd.1,FDR.q.sd.1), 
                    c(rep('MCF',n), rep('Habiger',n)))
colnames(FDR.sd) = c('alpha', 'sd', 'method')


p1 <- ggplot(FDR, aes(x = alpha, y = fdr, linetype=method, size=method)) + 
  geom_line() + 
  xlab("nominal FDR level") + 
  ylab("true FDR") +
  coord_cartesian(xlim = c(0, 0.15), ylim = c(0, 0.2)) +
  scale_linetype_manual(values=c(6,5,2,1,3)) +
  scale_size_manual(values=c(0.3,0.3,0.6,0.6,0.3)) +
  ggtitle(expression(paste(mu,' = 10')))
p1
p2 <- ggplot(Power, aes(x = alpha, y = power, linetype=method, size=method)) + 
  geom_line() + 
  xlab("nominal FDR level") + 
  ylab("true Power") +
  coord_cartesian(xlim = c(0, 0.2), ylim = c(0, 0.5)) +
  scale_linetype_manual(values=c(6,5,2,1,3)) + 
  scale_size_manual(values=c(0.3,0.3,0.6,0.6,0.3))
p2
p3 <- ggplot(FNR, aes(x = alpha, y = fnr, linetype=method, size=method)) + 
  geom_line() + 
  xlab("nominal FDR level") + 
  ylab("true FNR") +
  coord_cartesian(xlim = c(0, 0.15), ylim = c(0.05, 0.1)) +
  scale_linetype_manual(values=c(2,1,3)) +
  scale_size_manual(values=c(0.6,0.6,0.6)) +
  ggtitle(expression(paste(mu,' = 10')))
p3
p4 <- ggplot(FDR.sd, aes(x = alpha, y = sd, linetype=method, size=method)) + 
  geom_line() + 
  xlab("nominal FDR level") + 
  ylab("standard deviation of FDR") +
  coord_cartesian(xlim = c(0, 0.15), ylim = c(0, 0.04)) +
  scale_linetype_manual(values=c(2,1)) +
  scale_size_manual(values=c(0.3,0.3))
p4


n = length(alpha.all)
FDR = data.frame(rep(alpha.all,5), c(FDR.mcf.mean.2,FDR.q.mean.2,FDR.mid.mean.2,FDR.Chen.mean.2,FDR.Gilbert.mean.2), 
                 c(rep('MCF',n), rep('Habiger',n), rep('mid',n), rep('Chen',n), 
                   rep('Gilbert',n)))
colnames(FDR) = c('alpha', 'fdr', 'method')
Power = data.frame(rep(alpha.all,5), c(Power.mcf.mean.2,Power.q.mean.2,Power.mid.mean.2,Power.Chen.mean.2,Power.Gilbert.mean.2), 
                   c(rep('MCF',n), rep('Habiger',n), rep('mid',n), rep('Chen',n), 
                     rep('Gilbert',n)))
colnames(Power) = c('alpha', 'power', 'method')
pFNR_bound = FNR_bound_bt(lambda1=15)$pFNR_bound
FNR = data.frame(rep(alpha.all,3), c(FNR.mcf.mean.2,FNR.q.mean.2,pFNR_bound), 
                 c(rep('MCF',n), rep('Habiger',n), rep('pFDR_bound',n)))
colnames(FNR) = c('alpha', 'fnr', 'method')
FDR.sd = data.frame(rep(alpha.all,2), c(FDR.mcf.sd.2,FDR.q.sd.2), 
                    c(rep('MCF',n), rep('Habiger',n)))
colnames(FDR.sd) = c('alpha', 'sd', 'method')


p5 <- ggplot(FDR, aes(x = alpha, y = fdr, linetype=method, size=method)) + 
  geom_line() + 
  xlab("nominal FDR level") + 
  ylab("true FDR") +
  coord_cartesian(xlim = c(0, 0.15), ylim = c(0, 0.2)) +
  scale_linetype_manual(values=c(6,5,2,1,3)) +
  scale_size_manual(values=c(0.3,0.3,0.6,0.6,0.3)) +
  ggtitle(expression(paste(mu,' = 15')))
p5
p6 <- ggplot(Power, aes(x = alpha, y = power, linetype=method, size=method)) + 
  geom_line() + 
  xlab("nominal FDR level") + 
  ylab("true Power") +
  coord_cartesian(xlim = c(0, 0.2), ylim = c(0, 0.9)) +
  scale_linetype_manual(values=c(6,5,2,1,3)) + 
  scale_size_manual(values=c(0.3,0.3,0.6,0.6,0.3))
p6
p7 <- ggplot(FNR, aes(x = alpha, y = fnr, linetype=method, size=method)) + 
  geom_line() + 
  xlab("nominal FDR level") + 
  ylab("true FNR") +
  coord_cartesian(xlim = c(0, 0.15), ylim = c(0, 0.1)) +
  scale_linetype_manual(values=c(2,1,3)) +
  scale_size_manual(values=c(0.6,0.6,0.6)) +
  ggtitle(expression(paste(mu,' = 15')))
p7
p8 <- ggplot(FDR.sd, aes(x = alpha, y = sd, linetype=method, size=method)) + 
  geom_line() + 
  xlab("nominal FDR level") + 
  ylab("standard deviation of FDR") +
  coord_cartesian(xlim = c(0, 0.15), ylim = c(0, 0.025)) +
  scale_linetype_manual(values=c(2,1)) +
  scale_size_manual(values=c(0.3,0.3))
p8
### save the plot ####################
pdf("/Users/xiaoyudai/Documents/Paper/Simu_rep/BT/bt.pdf",width=8,height=12)
grid_arrange_shared_legend(p1,p5,p2, p6,p4,p8)
dev.off()

pdf("/Users/xiaoyudai/Documents/Paper/Simu_rep/BT/bt_FNR.pdf",width=8,height=12)
grid_arrange_shared_legend(p3,p7)
dev.off()




