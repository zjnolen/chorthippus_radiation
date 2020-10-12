## Read data. Each column is topology where rows are gene trees used
weightings_noerrortrees <- "weightings_md10first_saurabh"
weights <- read.table(weightings_noerrortrees, header = T)
#head(weights)
# Retrieve the names of the topologies
topoNames<- names(weights)

## Normalize rows so weights add to sum 1
weights_ratio <- weights / apply(weights, 1, sum)
#weights_ratio[,1]
## Exclude any rows where data is missing
good_rows = which(is.na(apply(weights,1,sum)) == F)
weights <- weights[good_rows,]

## Find the most common topology
# average_weightings<- apply(weights, 2, mean)
# average_weightings_dataf<- as.data.frame(average_weightings)
# which.max(average_weightings_dataf$average_weightings)
# average_weightings_dataf[8376,]

average_weightings_ratio<- apply(weights_ratio, 2, mean)
average_weightings_ratio_dataf<- as.data.frame(average_weightings_ratio)
which.max(average_weightings_ratio_dataf$average_weightings_ratio)


x <- as.character("rest")
x <- as.character(rep(x, 10395))
type<- as.data.frame(x)
type$x <- as.character(type$x)

topoNames<- as.data.frame(topoNames)
average_weightings_ratio_dataf<- cbind(topoNames, average_weightings_ratio_dataf)
average_weightings_ratio_dataf<- cbind(average_weightings_ratio_dataf, type)

bp<- ggplot(average_weightings_ratio_dataf, aes(x=reorder(topoNames,-average_weightings_ratio), y=average_weightings_ratio, fill=x))+ geom_bar(stat = "identity")+scale_fill_manual(values = c("rest"="gray"))+theme_void()
bp+theme(legend.position = "none")

average_weightings_ratio_dataf[[8376,3]] <- "top"
average_weightings_ratio_dataf[[501,3]] <- "top"
average_weightings_ratio_dataf[[4806,3]] <- "top"

geography<- c(7263,7269,7270,7308,7314,7315,7353,7359,7360,7366,7367,7368,7369,7370,7371,1194,1197,1198,1209,1212,1213,1224,1227,1228,1231,1232,1233,1234,1235,1236,2979,2982,2983,2994,2997,2998,3009,3012,3013,3016,3017,3018,3019,3020,3021)
#length(geography)

for(i in 1:length(geography)){
  average_weightings_ratio_dataf[[geography[i],3]] <- "geography"
}

## Find the totals
(average_weightings_ratio_dataf[8376,2]+average_weightings_ratio_dataf[501,2]+average_weightings_ratio_dataf[4806,2])*100  ## 0.5161998

for(i in 1:length(geography)){
  geo_sum<- sum(average_weightings_ratio_dataf[geography[i],2])
}
geo_sum*100   ## 0.007926024
#barplot(average_weightings_ratio)
average_weightings_ratio_dataf[8376,]

### There is three topology which are seen more common in the barplot, which one of them ??

### First :Topology 8376

## Find the second and third most common topology

# average_weightings_ordered<- average_weightings_dataf[with(average_weightings_dataf, order(-average_weightings)),]
# average_weightings_ordered_dataf <- as.data.frame(average_weightings_ordered)

average_weightings_ratio_ordered<- average_weightings_ratio_dataf[with(average_weightings_ratio_dataf, order(-average_weightings_ratio)),]
average_weightings_ratio_ordered_dataf <- as.data.frame(average_weightings_ratio_ordered)

#barplot(average_weightings_ratio_ordered)
library(ggplot2)
barplot<- ggplot(average_weightings_ratio_dataf, aes(x=reorder(topoNames,-average_weightings_ratio), y=average_weightings_ratio, fill=x))+ geom_bar(stat = "identity")+scale_fill_manual(values = c("top"="#FF0000", "geography"="blue", "rest"="white"))+theme_void()
barplot+theme(legend.position = "none")

average_weightings_ratio_dataf$average_weightings_ratio
which(average_weightings_ratio_dataf$average_weightings_ratio>=0.0017863482)
which(average_weightings_ratio_dataf$average_weightings_ratio>=0.0016893039)
which(average_weightings_ratio_dataf$average_weightings_ratio>=0.0016863455)
## second :Topology 501 and third : Topology 4806
## What is the average weightings
average_weightings_ratio_dataf[8376,2]*100   ## 0.1786348
average_weightings_ratio_dataf[4806,]*100   ## 0.1689304
average_weightings_ratio_dataf[501,]*100    ## 0.1686346

### Show with distribution plot (bar plot) only the first 3 and geographic topologies

#average_weightings_ratio_dataf<- cbind(topoNames, average_weightings_ratio_dataf)
ss_av_weightings_ratio_dataf<- average_weightings_ratio_dataf[c(8376,4806,501,7263,7269,7270,7308,7314,7315,7353,7359,7360,7366,7367,7368,7369,7370,7371,1194,1197,1198,1209,1212,1213,1224,1227,1228,1231,1232,1233,1234,1235,1236,2979,2982,2983,2994,2997,2998,3009,3012,3013,3016,3017,3018,3019,3020,3021),]

#ss_av_weightings_ratio_ordered<- ss_av_weightings_ratio_dataf[with(ss_av_weightings_ratio_dataf, order(ss_av_weightings_ratio)),]
#ss_av_weightings_ratio_ordered_dataf <- as.data.frame(ss_av_weightings_ratio_ordered)

barplot<- ggplot(ss_av_weightings_ratio_dataf, aes(x=reorder(topoNames,-average_weightings_ratio), y=average_weightings_ratio, fill=x))+ geom_bar(stat = "identity")+scale_fill_manual(values = c("top"="#FF0000", "geography"="blue"))+theme_void()
barplot+theme(legend.position = "none")

which(average_weightings_ratio_dataf$average_weightings_ratio>=9.577279e-05)


### The results from the second runs of datasets ###
### Mindepth 10 first run data second twisst

## Read data. Each column is topology where rows are gene trees used
weightings_noerrortrees <- "weightings_md10first_saurabh_2"
weights_2 <- read.table(weightings_noerrortrees, header = T)
#head(weights)
# Retrieve the names of the topologies
topoNames<- names(weights_2)
## Normalize rows so weights add to sum 1
weights_ratio_2 <- weights_2 / apply(weights_2, 1, sum)
#weights_ratio_2[,1]
## Exclude any rows where data is missing
good_rows = which(is.na(apply(weights_2,1,sum)) == F)
weights_2 <- weights_2[good_rows,]
## No bad row

## Find the most common topology

average_weightings_ratio_2<- apply(weights_ratio_2, 2, mean)
average_weightings_ratio_dataf_2<- as.data.frame(average_weightings_ratio_2)
which.max(average_weightings_ratio_dataf_2$average_weightings_ratio_2)
barplot(average_weightings_ratio_2)


x <- as.character("rest")
x <- as.character(rep(x, 10395))
type<- as.data.frame(x)
type$x <- as.character(type$x)

topoNames<- as.data.frame(topoNames)
average_weightings_ratio_dataf_2<- cbind(topoNames, average_weightings_ratio_dataf_2)
average_weightings_ratio_dataf_2<- cbind(average_weightings_ratio_dataf_2, type)

bp_2<- ggplot(average_weightings_ratio_dataf_2, aes(x=reorder(topoNames,-average_weightings_ratio_2), y=average_weightings_ratio_2, fill=x))+ geom_bar(stat = "identity")+scale_fill_manual(values = c("rest"="gray"))+theme_void()
bp_2+theme(legend.position = "none")

average_weightings_ratio_dataf_2[[8376,3]] <- "top"
average_weightings_ratio_dataf_2[[501,3]] <- "top"
average_weightings_ratio_dataf_2[[4806,3]] <- "top"

geography<- c(7263,7269,7270,7308,7314,7315,7353,7359,7360,7366,7367,7368,7369,7370,7371,1194,1197,1198,1209,1212,1213,1224,1227,1228,1231,1232,1233,1234,1235,1236,2979,2982,2983,2994,2997,2998,3009,3012,3013,3016,3017,3018,3019,3020,3021)
#length(geography)

for(i in 1:length(geography)){
  average_weightings_ratio_dataf_2[[geography[i],3]] <- "geography"
}


### There is again three topology which are seen more common in the barplot, which one of them ??

### First :Topology 501

## Find the second and third most common topology


average_weightings_ratio_ordered_2<- average_weightings_ratio_dataf_2[with(average_weightings_ratio_dataf_2, order(-average_weightings_ratio_2)),]
average_weightings_ratio_ordered_dataf_2 <- as.data.frame(average_weightings_ratio_ordered_2)
barplot(average_weightings_ratio_ordered_2)

which(average_weightings_ratio_dataf_2$average_weightings_ratio_2>=0.0017821183) ## 501
which(average_weightings_ratio_dataf_2$average_weightings_ratio_2>=0.0017565757) ## 8376
which(average_weightings_ratio_dataf_2$average_weightings_ratio_2>=1.63638e-03) ## 4806

## second :Topology 8376 and third : Topology 4806
## What is the average weightings
average_weightings_ratio_dataf_2[8376,]*100   ## 0.1756576
average_weightings_ratio_dataf_2[4806,]*100   ## 0.1636389
average_weightings_ratio_dataf_2[501,]*100    ## 0.1782118


ss_av_weightings_ratio_2<- average_weightings_ratio_dataf_2[c(8376,4806,501,7263,7269,7270,7308,7314,7315,7353,7359,7360,7366,7367,7368,7369,7370,7371,1194,1197,1198,1209,1212,1213,1224,1227,1228,1231,1232,1233,1234,1235,1236,2979,2982,2983,2994,2997,2998,3009,3012,3013,3016,3017,3018,3019,3020,3021),]


barplot_2<- ggplot(ss_av_weightings_ratio_2, aes(x=reorder(topoNames,-average_weightings_ratio_2), y=average_weightings_ratio_2, fill=x))+ geom_bar(stat = "identity")+scale_fill_manual(values = c("top"="#FF0000", "geography"="blue"))+theme_void()
barplot_2+theme(legend.position = "none")


### The second techical replicate of mindepth 10
## First run of Twisst
setwd("~/Documents/master thesis project/twisst/mindepth_10")
## Read data. Each column is topology where rows are gene trees used
weightings_noerrortrees <- "weightings_md10second_saurabh"
weights_second <- read.table(weightings_noerrortrees, header = T)
#head(weights_second)
# Retrieve the names of the topologies
topoNames<- names(weights_second)
## Normalize rows so weights add to sum 1
weights_ratio_second <- weights_second / apply(weights_second, 1, sum)
#weights_ratio_second[,1]
## Exclude any rows where data is missing
good_rows = which(is.na(apply(weights_ratio_second,1,sum)) == F)
weights_ratio_second <- weights_ratio_second[good_rows,]

## Find the most common topology

average_weightings_ratio_second<- apply(weights_ratio_second, 2, mean)
average_weightings_ratio_dataf_second<- as.data.frame(average_weightings_ratio_second)
which.max(average_weightings_ratio_dataf_second$average_weightings_ratio_second)
barplot(average_weightings_ratio_second)

### There is three topology which are seen more common in the barplot, which one of them ??

### First :Topology 813


x <- as.character("rest")
x <- as.character(rep(x, 10395))
type<- as.data.frame(x)
type$x <- as.character(type$x)

topoNames<- as.data.frame(topoNames)
average_weightings_ratio_dataf_second<- cbind(topoNames, average_weightings_ratio_dataf_second)
average_weightings_ratio_dataf_second<- cbind(average_weightings_ratio_dataf_second, type)

bp_second<- ggplot(average_weightings_ratio_dataf_second, aes(x=reorder(topoNames,-average_weightings_ratio_second), y=average_weightings_ratio_second, fill=x))+ geom_bar(stat = "identity")+scale_fill_manual(values = c("rest"="gray"))+theme_void()
bp_second+theme(legend.position = "none")

average_weightings_ratio_dataf_second[[8376,3]] <- "top"
average_weightings_ratio_dataf_second[[501,3]] <- "top"
average_weightings_ratio_dataf_second[[813,3]] <- "top"

geography<- c(7263,7269,7270,7308,7314,7315,7353,7359,7360,7366,7367,7368,7369,7370,7371,1194,1197,1198,1209,1212,1213,1224,1227,1228,1231,1232,1233,1234,1235,1236,2979,2982,2983,2994,2997,2998,3009,3012,3013,3016,3017,3018,3019,3020,3021)
#length(geography)

for(i in 1:length(geography)){
  average_weightings_ratio_dataf_second[[geography[i],3]] <- "geography"
}

## Find the second and third most common topology

average_weightings_ratio_ordered_second<- average_weightings_ratio_dataf_second[with(average_weightings_ratio_dataf_second, order(-average_weightings_ratio_second)),]
average_weightings_ratio_ordered_dataf_second <- as.data.frame(average_weightings_ratio_ordered_second)
barplot(average_weightings_ratio_ordered_second)


which(average_weightings_ratio_dataf_second$average_weightings_ratio_second>=1.452586e-03)
which(average_weightings_ratio_dataf_second$average_weightings_ratio_second>=1.398628e-03)
which(average_weightings_ratio_dataf_second$average_weightings_ratio_second>=1.35161e-03)

## second :Topology 501 and third : Topology 8376

ss_av_weightings_ratio_second<- average_weightings_ratio_dataf_second[c(8376,813,501,7263,7269,7270,7308,7314,7315,7353,7359,7360,7366,7367,7368,7369,7370,7371,1194,1197,1198,1209,1212,1213,1224,1227,1228,1231,1232,1233,1234,1235,1236,2979,2982,2983,2994,2997,2998,3009,3012,3013,3016,3017,3018,3019,3020,3021),]

barplot_second<- ggplot(ss_av_weightings_ratio_second, aes(x=reorder(topoNames,-average_weightings_ratio_second), y=average_weightings_ratio_second, fill=x))+ geom_bar(stat = "identity")+scale_fill_manual(values = c("top"="#FF0000", "geography"="blue"))+theme_void()
barplot_second+theme(legend.position = "none")

## Second run of Twisst

## Read data. Each column is topology where rows are gene trees used
weightings_noerrortrees <- "weightings_md10second_saurabh_2"
weights_second_2 <- read.table(weightings_noerrortrees, header = T)
#head(weights_second)
# Retrieve the names of the topologies
topoNames<- names(weights_second_2)
## Normalize rows so weights add to sum 1
weights_ratio_second_2 <- weights_second_2 / apply(weights_second_2, 1, sum)
#weights_ratio_second[,1]
## Exclude any rows where data is missing
good_rows = which(is.na(apply(weights_ratio_second_2,1,sum)) == F)
weights_ratio_second_2 <- weights_ratio_second_2[good_rows,]

## Find the most common topology

average_weightings_ratio_second_2<- apply(weights_ratio_second_2, 2, mean)
average_weightings_ratio_dataf_second_2<- as.data.frame(average_weightings_ratio_second_2)
which.max(average_weightings_ratio_dataf_second_2$average_weightings_ratio_second_2)
barplot(average_weightings_ratio_second_2)

### There is three topology which are seen more common in the barplot, which one of them ??

### First :Topology 8376
x <- as.character("rest")
x <- as.character(rep(x, 10395))
type<- as.data.frame(x)
type$x <- as.character(type$x)

topoNames<- as.data.frame(topoNames)
average_weightings_ratio_dataf_second_2<- cbind(topoNames, average_weightings_ratio_dataf_second_2)
average_weightings_ratio_dataf_second_2<- cbind(average_weightings_ratio_dataf_second_2, type)

bp_second_2<- ggplot(average_weightings_ratio_dataf_second_2, aes(x=reorder(topoNames,-average_weightings_ratio_second_2), y=average_weightings_ratio_second_2, fill=x))+ geom_bar(stat = "identity")+scale_fill_manual(values = c("rest"="gray"))+theme_void()
bp_second_2+theme(legend.position = "none")

average_weightings_ratio_dataf_second_2[[8376,3]] <- "top"
average_weightings_ratio_dataf_second_2[[501,3]] <- "top"
average_weightings_ratio_dataf_second_2[[669,3]] <- "top"

geography<- c(7263,7269,7270,7308,7314,7315,7353,7359,7360,7366,7367,7368,7369,7370,7371,1194,1197,1198,1209,1212,1213,1224,1227,1228,1231,1232,1233,1234,1235,1236,2979,2982,2983,2994,2997,2998,3009,3012,3013,3016,3017,3018,3019,3020,3021)
#length(geography)

for(i in 1:length(geography)){
  average_weightings_ratio_dataf_second_2[[geography[i],3]] <- "geography"
}

## Find the second and third most common topology

average_weightings_ratio_ordered_second_2<- average_weightings_ratio_dataf_second_2[with(average_weightings_ratio_dataf_second_2, order(-average_weightings_ratio_second_2)),]
average_weightings_ratio_ordered_dataf_second_2 <- as.data.frame(average_weightings_ratio_ordered_second_2)
barplot(average_weightings_ratio_ordered_second_2)

which(average_weightings_ratio_dataf_second_2$average_weightings_ratio_second_2>=1.399765e-03)
which(average_weightings_ratio_dataf_second_2$average_weightings_ratio_second_2>=1.398153e-03)
which(average_weightings_ratio_dataf_second_2$average_weightings_ratio_second_2>=1.315452e-03)

## second :Topology 501 and third : Topology 669


ss_av_weightings_ratio_second_2<- average_weightings_ratio_dataf_second_2[c(8376,669,501,7263,7269,7270,7308,7314,7315,7353,7359,7360,7366,7367,7368,7369,7370,7371,1194,1197,1198,1209,1212,1213,1224,1227,1228,1231,1232,1233,1234,1235,1236,2979,2982,2983,2994,2997,2998,3009,3012,3013,3016,3017,3018,3019,3020,3021),]

barplot_second_2<- ggplot(ss_av_weightings_ratio_second_2, aes(x=reorder(topoNames,-average_weightings_ratio_second_2), y=average_weightings_ratio_second_2, fill=x))+ geom_bar(stat = "identity")+scale_fill_manual(values = c("top"="#FF0000", "geography"="blue"))+theme_void()
barplot_second_2+theme(legend.position = "none")
