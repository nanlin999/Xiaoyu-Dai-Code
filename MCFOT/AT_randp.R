##############################################################################
## These functions implement 
##   the Accumulation Test methods from the paper:
## Ang Li & Rina Foygel Barber,
##   "Accumulation tests for FDR control
##      in ordered hypothesis testing"
## Available from http://arxiv.org/abs/1505.07352
## (Several methods from other papers also implemented,
##        as noted below - see citations in paper)
##############################################################################

##############################################################################
## HingeExp method,
##    i.e. an accumulation test with the HingeExp function:
##       h(p) = C * log(1/(C*(1-p))) * (p>1-1/C)
##############################################################################
create_HingeExp_function = function(C=2){
  function(x){C*log(1/(C*(1-x)))*(x>1-1/C)}
}
HingeExp = function(pvals,alpha=0.2,C=2){
  AT(pvals,create_HingeExp_function(C),alpha=alpha)
}
HingeExp_randp = function(p.org,p.next,alpha=0.2,C=2){
  AT_randp(p.org,p.next,n_rand=100,create_HingeExp_function(C),alpha=alpha)
}

##############################################################################


##############################################################################
## ForwardStop method (G'Sell et al 2013),
##    i.e. an accumulation test with the ForwardStop function:
##       h(p) = log(1/(1-p))
##############################################################################
create_ForwardStop_function = function(){
  function(x){log(1/(1-x))}
}
ForwardStop = function(pvals,alpha=0.2){
  AT(pvals,create_ForwardStop_function(),alpha=alpha)
}
ForwardStop_randp = function(p.org,p.next,alpha=0.2){
  AT_randp(p.org,p.next,n_rand=100,create_ForwardStop_function(),alpha=alpha)
}
##############################################################################


##############################################################################
## SeqStep method (Barber&Candes 2015),
##    i.e. an accumulation test with the step function:
##       h(p) = C * (p>1-1/C)
##############################################################################
create_SeqStep_function = function(C=2){
  function(x){C*(x>1-1/C)}
}
SeqStep = function(pvals,alpha=0.2,C=2){
  AT(pvals,create_SeqStep_function(C),alpha=alpha)
}
SeqStep_randp = function(p.org,p.next,alpha=0.2,C=2){
  AT_randp(p.org,p.next,n_rand=100,create_SeqStep_function(C),alpha=alpha)
}
####################################################################################


##############################################################################
## SeqStep+ method (Barber&Candes 2015),
##    i.e. an accumulation test with the step function:
##       h(p) = C * (p>1-1/C)
##         & with the conservative correction
##              for estimating FDR
##############################################################################
SeqStepPlus = function(pvals,alpha=0.2,C=2){
  AT(pvals,create_SeqStep_function(C),alpha=alpha,numerator_plus=C,denominator_plus=1)
}
SeqStepPlus_randp = function(p.org,p.next,alpha=0.2,C=2){
  AT_randp(p.org,p.next,n_rand=100,create_SeqStep_function(C),alpha=alpha,numerator_plus=C,denominator_plus=1)
}
##############################################################################


##############################################################################
## Accumulation test for a generic function "hfun"
# Accumulation Test (Li&Barber, 2016)
# ref: https://www.stat.uchicago.edu/~rina/accumulationtests.html
##############################################################################
AT = function(p.org, h_function, alpha=0.1, numerator_plus=0, denominator_plus=0){
  n.test = length(p.org)
  
  # h_funtion at each raw p-value
  h_eval = unlist(lapply(p.org,h_function))
  
  k_hat = max(0, which((numerator_plus + cumsum(h_eval))/(denominator_plus + 1:n.test) <= alpha))
  return(k_hat)
}
##############################################################################

# Accumulative Test based on randomized p-values
##############################################################################
AT_randp <- function(p.org, p.next, n_rand = 100, h_function, alpha=0.1, numerator_plus=0, denominator_plus=0){
  # n_rand: number of sampels of randomized p-value
  
  n.test = length(p.org)
  
  # h_funtion at each raw p-value
  h_eval = rep(0,n.test)
  for(i in 1:n_rand){
    p.rand = runif(n.test,p.next,p.org)
    h_eval = h_eval + unlist(lapply(p.rand,h_function))
  }
  h_eval = h_eval/n_rand
  
  k_hat = max(0, which((numerator_plus + cumsum(h_eval))/(denominator_plus + 1:n.test) <= alpha))
  return(k_hat)
}


##############################################################################

