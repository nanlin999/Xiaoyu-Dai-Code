library(CRF)

NUM = 100 # number of nodes in each CRF
n.train = 10

mc.train = matrix(0, nrow = n.train, ncol = NUM)
hmm.train = matrix(0, nrow = n.train, ncol = NUM)
for(i in 1:n.train){
  dt = rdata.hmm(NUM, pii, A, f0, 1, f1)
  mc.train[i,] = dt$s + 1   # status: 1 or 2
  hmm.train[i,] = dt$o
}


n.nodes <- NUM
n.states <- 2

adj <- matrix(0, n.nodes, n.nodes)
for (i in 1:(n.nodes-1))
{
  adj[i, i+1] <- 1
}


mrf.new <- make.crf(adj, n.states)
mrf.new <- make.features(mrf.new)
mrf.new <- make.par(mrf.new, 4)
mrf.new$node.par[1,1,1] <- 1
for (i in 1:mrf.new$n.edges)
{
  mrf.new$edge.par[[i]][1,1,1] <- 2
  mrf.new$edge.par[[i]][1,2,1] <- 3
  mrf.new$edge.par[[i]][2,1,1] <- 4
}
mrf.new <- train.mrf(mrf.new, mc.train)
mrf.new$par
mrf.new$node.pot <- mrf.new$node.pot / rowSums(mrf.new$node.pot)
mrf.new$edge.pot[[1]] <- mrf.new$edge.pot[[1]] / rowSums(mrf.new$edge.pot[[1]])
mrf.new$node.pot[1,]
mrf.new$edge.pot[[1]]


crf.new <- make.crf(adj, n.states)
crf.new <- make.features(crf.new, 2, 1)
crf.new <- make.par(crf.new, 5)
crf.new$node.par[1,1,1] <- 1
for (i in 1:crf.new$n.edges)
{
  crf.new$edge.par[[i]][1,1,] <- 2
  crf.new$edge.par[[i]][1,2,] <- 3
  crf.new$edge.par[[i]][2,1,] <- 4
}
crf.new$node.par[,1,2] <- 5
# crf.new$node.par[,1,3] <- 6


hmm.nf <- lapply(1:dim(hmm.train)[1], function(i) matrix(1, crf.new$n.nf, crf.new$n.nodes))
for (i in 1:dim(hmm.train)[1])
{
  hmm.nf[[i]][2, ] <- log(dnorm(hmm.train[i,],f0[1],f0[2])/dnorm(hmm.train[i,],f1[1],f1[2])) # log(f0/f1)
  # hmm.nf[[i]][2, ] <- log(dnorm(hmm.train[i,],0,1))
  # hmm.nf[[i]][2, ] <- log(dnorm(hmm.train[i,],1,1))
}
hmm.ef <- lapply(1:dim(hmm.train)[1], function(i) matrix(1, crf.new$n.ef, crf.new$n.edges))

crf.new <- train.crf(crf.new, mc.train, hmm.nf, hmm.ef)
crf.new$par

crf.new$node.pot <- crf.new$node.pot / rowSums(crf.new$node.pot)
crf.new$edge.pot[[1]] <- crf.new$edge.pot[[1]] / rowSums(crf.new$edge.pot[[1]])
crf.new$node.pot[1,] # initial prob
crf.new$edge.pot[[1]] # estimated traisition matrix

hmm.infer <- matrix(0, nrow=dim(hmm.train)[1], ncol=dim(hmm.train)[2])
for (i in 1:dim(hmm.train)[1])
{
  crf.new <- crf.update(crf.new, hmm.nf[[i]], hmm.ef[[i]])
  hmm.infer[i,] <- decode.chain(crf.new)
}
sum(hmm.infer != mc.train)/length(mc.train)




# test 

# n.test = 100
# 
# mc.test = matrix(0, nrow = n.test, ncol = NUM)
# hmm.test = matrix(0, nrow = n.test, ncol = NUM)
# for(i in 1:n.test){
#   dt = rdata.hmm(NUM, pii, A, f0, 1, f1)
#   mc.test[i,] = dt$s + 1   # status: 1 or 2
#   hmm.test[i,] = dt$o
# }
# 
# test.nf <- lapply(1:dim(hmm.train)[1], function(i) matrix(1, crf.new$n.nf, crf.new$n.nodes))
# for (i in 1:dim(hmm.test)[1])
# {
#   test.nf[[i]][2, ] <- log(dnorm(hmm.test[i,],1,1)/dnorm(hmm.test[i,],0,1)) # log(f1/f0)
# }
# test.ef <- lapply(1:dim(hmm.test)[1], function(i) matrix(1, crf.new$n.ef, crf.new$n.edges))
# 
# hmm.infer <- matrix(0, nrow=dim(hmm.test)[1], ncol=dim(hmm.test)[2])
# for (i in 1:dim(hmm.test)[1])
# {
#   crf.new <- crf.update(crf.new, hmm.nf[[i]], hmm.ef[[i]])
#   hmm.infer[i,] <- decode.chain(crf.new)
# }


# sum(hmm.infer != mc.test)

# res = infer.chain(crf.new)




# the HMM data
n.test = 1000 # num of tests
set.seed(4321)
rdata1<-rdata.hmm(n.test, pii, A, f0, 1, f1)
# the observed values
x1<-rdata1$o
# the unobserved states
theta1<-rdata1$s
# the EM algorithm
em.res1<-em1.hmm(x1, maxiter=500) # f0 is always N(0,1)
# HMM procedure
lsi.pi<-em.res1$lf
pi.res<-mt.hmm(lsi.pi, q)
# the decision rule
pi.de<-pi.res$de



hmm.test = matrix(x1, nrow = n.test/NUM, ncol = NUM, byrow = T)

test.nf <- lapply(1:dim(hmm.test)[1], function(i) matrix(1, crf.new$n.nf, crf.new$n.nodes))
for (i in 1:dim(hmm.test)[1])
{
  test.nf[[i]][2, ] <- log(dnorm(hmm.test[i,],f0[1],f0[2])/dnorm(hmm.test[i,],f1[1],f1[2])) # log(f0/f1)
}
test.ef <- lapply(1:dim(hmm.test)[1], function(i) matrix(1, crf.new$n.ef, crf.new$n.edges))

hmm.infer <- matrix(0, nrow=dim(hmm.test)[1], ncol=dim(hmm.test)[2])
for (i in 1:dim(hmm.test)[1])
{
  crf.new <- crf.update(crf.new, test.nf[[i]], test.ef[[i]])
  hmm.infer[i,] <- infer.chain(crf.new)$node.bel[,1]
  # hmm.infer[i,] <- decode.chain(crf.new)
}
# sum(hmm.infer!=(theta1+1))

lsi.crf <- as.vector(t(hmm.infer))
crf.res<-mt.hmm(lsi.crf, q)
# the decision rule
crf.de<-crf.res$de
#

1 - sum(crf.de*theta1)/sum(crf.de)

pi.fdp = 1 - sum(pi.de*theta1)/sum(pi.de)
crf.fdp = 1 - sum(crf.de*theta1)/sum(crf.de)




idx0 <- which(theta1==0)
idx1 <- which(theta1==1)


plot(lsi.pi[idx1],cex=.5)

hist(lsi.pi[idx1])
hist(lsi.pi[idx0],add=T,col='red')

hist(lsi.crf[idx1])
hist(lsi.crf[idx0],add=T)


pi.res<-mt.hmm(lsi.pi, 0.)
# the decision rule
pi.de<-pi.res$de
pi.infer <- rep(0,1000)
pi.infer[pi.de] = 1
sum(pi.infer!=theta1)







