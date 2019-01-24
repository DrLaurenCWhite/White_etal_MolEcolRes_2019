# library(lme4)
# library(car)
# library(gtools)
# library(pscl)
# library(rela)
# library(MuMIn)
# 
# setwd("~/Dropbox/LinearModsCours/DataFiles")
# 
# rm(list=ls()) #clears the workspace. Removes any defined objects
# 
# ##IMPORT DATA
# xdata=read.table(file="02_water_cons.txt", header=T, sep="\t")
# str(xdata)
# 
# ##Fit Model
# res=lm(data=xdata, liters.per.hour~number.participants) #fits the model response~predictor(s)
# 
# #What is the max cook's distance?
# max(cooks.distance(res))
# 
# attributes(res)
# length(res$residuals)

#Does the max cook's distance exceed any of the three thresholds?
CooksThres=function(model.res, alpha=0.05){
  cookie=max(abs(cooks.distance(model.res)))
  eins=1
  zwei=4/(nobs(model.res))
  drei=qf(p=1-alpha, df1=length(model.res$coefficients), df2=(length(model.res$residuals)-length(model.res$coefficients)), lower.tail=T)
  thres=c(eins, zwei, drei)
  worry=cookie>thres
  df=data.frame(thres,worry)
  return(df)
}


