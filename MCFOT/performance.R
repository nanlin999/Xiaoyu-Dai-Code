performance <- function(rej, H){
  R = length(rej) # total number of rejections
  V = sum(H[rej]==0) # false rejections
  S = R-V # true rejections
  fdr = V/R
  power = S/sum(H==1)
  
  return(list(R=R, fdr=fdr, power=power))
}