rm(list = ls())

library(data.table)
library(car)
source("./diagnostic_fcns.r")
source("./CooksDistance_Thres.R")
source("./GetDFBeta.R")

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}



All=read.csv("Dropout_Overseq_Data_all_withoutXY.csv")
All$Rounds.of.Capture=as.factor(All$Rounds.of.Capture)


###########################USING UNIQUE MAPPED READS###############################

#DROPOUT
#Inspect Predictors =CapRounds, PercEnd and SeqEffort
hist(All$Unique.Exome.Reads, breaks=30) #Fine

table(All$Rounds.of.Capture) 
plot(All$Rounds.of.Capture, All$Unique.Exome.Reads)


#Insepct Response
hist(All$DropOut, breaks=50) #Maybe some influential cases?


res=lm(data=All, DropOut~Rounds.of.Capture+Unique.Exome.Reads)
null=lm(data=All, DropOut~Unique.Exome.Reads)
anova(res, null)
summary(res)
confint(res)


diagnostics.plot(res)
max(abs(dffits(res))) #threshold is 2

#Cooks distance
max(cooks.distance((res)))
CooksThres(res)

#leverage
max(as.vector(influence(res)$hat))
lev.thresh(res) #too high above the threshold? Some influential cases?

#Collinearity
vif(res)

#Check whether influential cases are a problem
lev=as.vector(influence(res)$hat)
l.thresh=lev.thresh(res)
sel.data=subset(All, lev<l.thresh)
sel.res=lm(data=All, DropOut~Rounds.of.Capture+Unique.Exome.Reads)
diagnostics.plot(sel.res)
summary(sel.res) #nah, basically the same





#OVERSEQ
#Inspect Predictors =CapRounds, PercEnd and SeqEffort
hist(All$Unique.Exome.Reads, breaks=30) #Fine


table(All$Rounds.of.Capture) 
plot(All$Rounds.of.Capture, All$Unique.Exome.Reads)


#Insepct Response
hist(All$OverSeq, breaks=50)


res=lm(data=All, OverSeq~Rounds.of.Capture+Unique.Exome.Reads)
null=lm(data=All, OverSeq~Unique.Exome.Reads)
anova(res, null)
summary(res)
confint(res)


diagnostics.plot(res)
max(abs(dffits(res))) #threshold is 2

#Cooks distance
max(cooks.distance((res)))
CooksThres(res)

#leverage
max(as.vector(influence(res)$hat))
lev.thresh(res) #too high above the threshold?

#Collinearity (without interaction)
vif(res)

#Check whether influential cases are a problem
lev=as.vector(influence(res)$hat)
l.thresh=lev.thresh(res)
sel.data=subset(All, lev<l.thresh)
sel.res=lm(data=All, OverSeq~Rounds.of.Capture+Unique.Exome.Reads)
diagnostics.plot(sel.res)
summary(sel.res) #nah, basically the same

