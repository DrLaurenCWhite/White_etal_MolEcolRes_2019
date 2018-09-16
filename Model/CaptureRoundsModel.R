rm(list = ls())

library(car)
source("./CooksDistance_Thres.R")
source("./GetDFBeta.R")
source("./Leverage_Thres.R")


#Import Data. Set Rounds of Capture as a factor
fox = read.csv("ModelData.csv", header=TRUE, na.strings=c(NA))
str(fox)
fox$CapRounds=as.factor(fox$CapRounds)


#Inspect Predictors = Rounds of Capture (CapRounds), Percent Subject DNA (PercEnd) and Sequencing effort (SeqEffort)
hist(fox$Shot.PercEnd, breaks=30) #A little skewed, try log transform
min(fox$Shot.PercEnd) # Minimum is above zero
hist(log(fox$Shot.PercEnd), breaks=30) #Better
#fox$SampPercEnd.log=log(fox$Shot.PercEnd) <- Already done

hist(fox$Prod.Reads)

table(fox$CapRounds)
plot(fox$CapRounds, fox$Shot.PercEnd)
plot(fox$CapRounds, fox$Prod.Reads)

#z.tranform covariates
#fox$SampPercEnd.log.z=scale(fox$SampPercEnd.log) <- Already done
#fox$Prod.Reads.z=scale(fox$Prod.Reads) <- Already done

#Implement model
res=lm(data=fox, Unique.Exome.Reads~CapRounds*SampPercEnd.log.z*Prod.Reads.z)


#Check assumptions
#Normality of residuals
qqnorm(residuals(res))
qqline(residuals(res))  # Maybe some influential cases.

max(abs(dffits(res))) # Typical threshold is 2. Looks fine

#Cooks distance
max(cooks.distance((res)))
CooksThres(res) #Above one of three typical thresholds

#leverage
max(as.vector(influence(res)$hat))
LevThres(res) #Above both our thresholds. Possibly some overly influential cases

#Collinearity (without interaction)
red=lm(data=fox, Unique.Exome.Reads~CapRounds+SampPercEnd.log.z+Prod.Reads.z)
vif(red) #Looks fine. Max VIF=1.13


#Check whether influential cases are a problem
lev.thresh<-function(model.res){
  k=length(coefficients(model.res))
  n=length(residuals(model.res))
  return(2*(k+1)/n)
}
l.thresh=lev.thresh(res) #Get leverage threshold
lev=as.vector(influence(res)$hat) 
sel.data=subset(fox, lev<l.thresh) # Extract data excluding potentially influential cases
sel.res=lm(data=sel.data, Unique.Exome.Reads~CapRounds*SampPercEnd.log.z*Prod.Reads.z) #Run reduced model without potentially influential cases
summary(sel.res) #Compare reduced model to full-data model
summary(res) #Very similar results, especially for three-way interaction of interest. Keep full-data model

#Full-null model comparison
null=lm(data=fox, Unique.Exome.Reads~1)
anova(res, null) #Full-null comparison. Is significant.

#Finally, look at model output. 
summary(res)

coefficients(res)
confint(res)
