#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// input: x: original p-values
//        y: next possible p-values
//        lambda: threshould on randomized p-values
// output: MCF for each test


// [[Rcpp::export]]
NumericVector mcf(NumericVector x, NumericVector y, double lambda) {
  int n = x.size();
  NumericVector r(n);
  for(int i=0; i<n; i++){
    if(lambda > x[i]){
      r[i] = 1;
    }else if(lambda < y[i]){
      r[i] = 0;
    }else{
      r[i] = (lambda-y[i])/(x[i]-y[i]);
    }  
  }
  
    
  return r;
}
