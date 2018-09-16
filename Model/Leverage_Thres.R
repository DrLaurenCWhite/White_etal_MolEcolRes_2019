
LevThres=function(mod.res){
  cookie=max(as.vector(influence(mod.res)$hat))
  eins=(2*length(coefficients(mod.res)))/length(residuals(mod.res))
  zwei=(3*length(coefficients(mod.res)))/length(residuals(mod.res))
  thres=c(eins, zwei)
  worry=cookie>thres
  df=data.frame(thres,worry)
  return(df)
}