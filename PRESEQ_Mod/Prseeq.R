#for a clean workspace
rm(list = ls())

#load the usual culprits
require(magrittr); require(data.table); require(plyr); require(parallel); require(dplyr); library(ggplot2)


#Import Preseq Output into list. Edit dataframe names
import.multiple.csv.files=function(mypath,mypattern)
{
  tmp.list.1<-list.files(mypath, pattern=mypattern)
  tmp.list.2<-list(length=length(tmp.list.1))
  for (i in 1:length(tmp.list.1)){tmp.list.2[[i]]<-read.csv(tmp.list.1[i], header=TRUE)}
  names(tmp.list.2)<-tmp.list.1
  tmp.list.2
}
setwd("./SampleExomeProjection")
preseqfiles =import.multiple.csv.files("./", mypattern="_Exome_projection.csv$")
setwd("../")
y=as.vector(NULL)
for (i in names(preseqfiles)) {
  x=sub("_Exome_projection.csv","", i)
  y=c(y, x)
  names(preseqfiles)=c(y)
}

#Get Mapping Rate infomation
trot=read.csv("MappingRate.csv")
trot <- as.data.table(trot)

#flatten the list items into a single data.table
data <- lapply(preseqfiles %>% names , function(d) as.data.table(preseqfiles[[d]])[ , Sample := d]  ) %>% bind_rows
#attach the mapping rate info to the preseq data and calculate adjusted total
data <-  trot[data,on="Sample"][,adjusted_total_reads := TOTAL_READS / Enriched.TargetMapRate][]


#write.csv(data, "adjustedtotalReads.csv") #<-already done


	#Plot Preseq Curves for ALL libraries

#Function to get scientific notation for plots
fancy_scientific <- function(l) {
  l <- format(l, scientific = TRUE)
  l <- gsub("^(.*)e", "'\\1'e", l)
  l <- gsub("e", "%*%10^", l)
  parse(text=l)
}

#convert rounds to factor
data$CapRounds=as.factor(data$CapRounds)

#Colour by Rounds of Capture
ggplot(data, aes(group = Sample , x = adjusted_total_reads ,  y = EXPECTED_DISTINCT, col=CapRounds)) +
  geom_line(size=0.5) +
  scale_y_continuous(labels = fancy_scientific, limits=c(0, 1e7)) +
  scale_x_continuous(labels=c(0, 5, 10, 15, 20),limits = c(0, 20000000)) +
  xlab("Total Reads Sequenced") + ylab("Predicted Number of\nUnique Exome Reads") +
  theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
        panel.background = element_blank(),
        axis.line = element_line(colour="black"),
        panel.grid.major = element_line(colour="#f0f0f0"),
        axis.text = element_text(size=15),
        axis.title = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

#Colour by Percent Subject DNA
ggplot(data, aes( group = Sample , x = adjusted_total_reads ,  y = EXPECTED_DISTINCT, col=Shot.PercEnd)) +
  geom_line(size=0.5) + 
  scale_y_continuous(labels = fancy_scientific, limits=c(0, 1e7)) +
  scale_x_continuous(labels=c(0, 5, 10, 15, 20),limits = c(0, 20000000)) +
  scale_color_gradient(name="Percent\nSubject\nDNA", trans="log", breaks=c(3,30), low="red", high="blue") +
  xlab("Total Reads Sequenced (Millions of Reads)") + ylab("Predicted Number of\nUnique Exome Reads") +
  theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
        panel.background = element_blank(),
        axis.line = element_line(colour="black"),
        panel.grid.major = element_line(colour="#f0f0f0"),
        axis.text = element_text(size=15),
        axis.title = element_text(size=15),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))


#a function that returns a linearly interpolated y at point x, given the x and y coordinates of two other points
interpolate_y_at_x_from_two_neighbours <- function(neighbour_xs,neighbour_ys,x){
  slope <- (neighbour_ys[2]-neighbour_ys[1])/(neighbour_xs[2]-neighbour_xs[1])
  dist_from <- x - neighbour_xs[1]
  y_from <- neighbour_ys[1]
  y_from + (slope * dist_from)
}

#must be adjusted for each output. If there are future versions of this, this will be automated.
x_query <- 15e6 #or 20e6 or 17e6 or 14e6 or whatever 

#isolate points either side of the interpolation point (x_query)
caught <- data[ , .(dist_from_xquery = x_query - adjusted_total_reads , adjusted_total_reads , EXPECTED_DISTINCT , nearest = x_query - adjusted_total_reads ) , by="Sample"]
catcher <- data[,.N,Sample][,N := NULL][ , dist_from_xquery := 0][]

#select datapoints just below (a) and above (b) the query (for each sample)
b <- caught[catcher,,on=c("Sample","dist_from_xquery"),roll=T]
a <- caught[catcher,,on=c("Sample","dist_from_xquery"),roll=-Inf]
near <- bind_rows(a,b); rm(a,b)

near[order(Sample)]

#perform the interpolation
#contains a harmless self reference that might generate an ignorable warning
interp <- near[, exp_distinct_at_xquery := interpolate_y_at_x_from_two_neighbours(neighbour_xs = adjusted_total_reads , neighbour_ys = EXPECTED_DISTINCT , x = x_query ) ,Sample]

#Plot points that we've just interpolated
ggplot(data, aes( group = Sample , x = adjusted_total_reads ,  y = EXPECTED_DISTINCT )) + geom_line(size=.2) + geom_point(data=interp,aes(x=x_query , y = exp_distinct_at_xquery)) + xlim(c(0,1e8))
 
#save the data 
write.csv(interp[,.(Sample,exp_distinct_at_xquery)], "expected_at_15M.csv") #or name as appropriate for x_query


