# MCF-based discrete multiple testing function for ordered multiple testing problem
# now adaptive to 4 procedures: 
#   Storey-BH (SBH): Storey, J.D. (2002) Adirect approach to false discovery rate.
#   Adaptive Selective SeqStep (AS): Lei, L., Fithian, W. (2016) Power of Ordered Hypothesis Testing.
#   Accumulation Test (AT): Li A., Barber, R.F. (2016) Accumulation Tests for FDR Control in Ordered Hypothesis Testing.
#   Weighted BH (WBH): Li A., Barber, R.F. (2016) Multiple Testing with the Structure Adaptive Benjamini-Hochberg Algorithm.

MCF_main <- function(p_org, p_next, alpha = 0.2, method = 'SBH', sample_num = 100, lambda = 0.5, threshould_AT = 0.05, hfun = NULL, weight = NULL){
  # Inputs:
  #   p_org: raw p-values
  #   p_next: largest possible p-value smaller than raw p-value. 0 if not exist
  #   alpha: nominal FDR level. vector is not acceptable
  #   method: one of 4 procedures that MCF is adaptive to: 'SBH', 'AS', 'AT', 'WBH'
  #   sample_num: number of samples for randomized p-values in the MCF-based method
  #   lambda: cutoff to estimate the null proportion pi_0, only used when method is 'SBH' or 'AS'
  #   threshould_AS: pre-specified threshould for each test, only used when method is 'AS'
  #   hfun: accumulative function, only used when mehtod is 'AT'
  #   weight: weight for each test, only used when method is 'WBH'
  #
  # Output:
  #   matrix of decisions: 0 if accepted, 1 if rejected
  
  if(method == 'SBH'){
    output = MCF_SBH(p_org, p_next, alpha, sample_num, lambda)
  } else if(method == 'AS'){
    output = MCF_AS(p_org, p_next, alpha, sample_num, lambda, threshould_AS)
  } else if(method == 'AT'){
    output = MCF_AT(p_org, p_next, alpha, sample_num, hfun)
  } else if(method == 'WBH'){
    output = MCF_WBH(p_org, p_next, alpha, sample_num, weight)
  } else {
    stop('unrecognized method')
  }
  
  return(output)
}



MCF_BH <- function(p_org, p_next, alpha, sample_num, lambda){
  n_test = length(p_org)
  sample_randomized_p = runif(n_test*sample_num, min = p_next, max = p_org)
  s = alpha
  if(s/(1-lambda)*(1+sum(sample_randomized_p>lambda))/(sum(sample_randomized_p<=lambda)) >= alpha){
    s = s - alpha/10000
  }
  
  mcf_value = 
}















