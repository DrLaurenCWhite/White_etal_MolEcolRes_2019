gdfb=function(mres, my.d=3){
  xx=cbind(coefficients(mres), coefficients(mres) + t(apply(X=dfbeta(mres), MARGIN=2, FUN=range))) #t transposes the data
  colnames(xx)=c("orig", "min", "max")
  round(xx, digits=my.d) #looks fine, not much variation from the original.
  return(xx)
}