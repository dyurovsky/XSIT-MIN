# Clear All Previous Variables
rm(list=ls())

# Get Lab Version of Useful R Functions
source('~/Projects/Other/Ranalysis/useful.R')

library(rjson)
library(chron)

library(dplyr)

# Make a List of All Files to Read
e1.results <- list.files(path = "../data/raw/exp1/", pattern = '*.results',
                         all.files = FALSE)
e1.results <- sapply(e1.results,function(x){paste("exp1/",x,sep="")},
                     USE.NAMES=FALSE)
e1.results <- data.frame(File=e1.results,exp="Exp 1",stringsAsFactors=FALSE)

e1.results$expt <- "Exp1"

e2.results <- list.files(path = "../data/raw/exp2/", pattern = '*.results',
                         all.files = FALSE)
e2.results <- sapply(e2.results,function(x){paste("exp2/",x,sep="")},
                     USE.NAMES=FALSE)
e2.results <- data.frame(File=e2.results,exp="Exp 2",stringsAsFactors=FALSE)
all.results <- merge(e1.results,e2.results,all=TRUE)


# Prepare Matrix to Hold all of the Data
all.data <- as.data.frame(matrix(ncol = 0, nrow = 0))

# Main Loop Through All Conditions -- SLOW
for(f in 1:nrow(all.results)) {

  data <- read.table(paste("../data/raw/", all.results$File[f],sep=""),
                     sep="\t",header=TRUE, stringsAsFactors=FALSE)

  
  long.data <- as.data.frame(matrix(ncol = 0, nrow = 20*nrow(data)))
  c <- 1 #participant counter
  
  for (i in 1:nrow(data)) {
    d <- fromJSON(as.character(data$Answer.data[i])) #participant data
    
    for (j in 1:length(d)) {
      
      #Grab Relevant Fields from the JSON Mess
      long.data$subid[c] <- data$workerid[i]
      long.data$submit.date[c] <-  paste(word(data$assignmentsubmittime[i],
                                              start=2,end=3),
                                         word(data$assignmentsubmittime[i]
                                              ,start=-1L))
      long.data$submit.time[c] <-  word(data$assignmentsubmittime[i]
                                              ,start=4)
      long.data$trial.num[c] <- j
      long.data$itemNum[c] <- d[[j]]$itemNum
      long.data$trialType[c] <- d[[j]]$trialType
      long.data$samePos[c] <- d[[j]]$samePos
      long.data$chosen[c] <- d[[j]]$chosen
      long.data$chosenIdx[c] <- d[[j]]$chosen_idx
      long.data$kept[c] <- d[[j]]$kept
      long.data$keptIdx[c] <- d[[j]]$kept_idx
      long.data$rt[c] <- d[[j]]$rt

      c <- c + 1 #increment participant counter
    }
  }
  
  #get Num Refs/Interval from File Names
  subs <- strsplit(all.results$File[f],'\\.');
  subs <- strsplit(unlist(subs)[1],'_');
  long.data$numPic <- unlist(subs)[4];
  interval <- unlist(subs)[5];
  long.data$interval <- interval;
  
  long.data$test <- 0
  long.data$exp = all.results[f,"exp"]
  long.data$subid = paste(long.data$subid,long.data$exp)
  
  #grab only the continuation trials
  # this is super ugly
  if (interval == 0) {
    long.data$test[long.data$trial.num==6 |
                     long.data$trial.num==8 |
                     long.data$trial.num==10 |
                     long.data$trial.num==12 |
                     long.data$trial.num==14 |
                     long.data$trial.num==16 |
                     long.data$trial.num==18 |
                     long.data$trial.num==20] <- 1;
  } else if(interval == 1){
    long.data$test[long.data$trial.num==7 |
                     long.data$trial.num==8 |
                     long.data$trial.num==11 |
                     long.data$trial.num==12 |
                     long.data$trial.num==15 |
                     long.data$trial.num==16 |
                     long.data$trial.num==19 |
                     long.data$trial.num==20] <- 1;
  } else if(interval == 3){
    long.data$test[long.data$trial.num==9 |
                     long.data$trial.num==10 |
                     long.data$trial.num==11 |
                     long.data$trial.num==12 |
                     long.data$trial.num==17 |
                     long.data$trial.num==18 |
                     long.data$trial.num==19 |
                     long.data$trial.num==20] <- 1;  
  } else if(interval == 7){
    long.data$test[long.data$trial.num==13 |
                     long.data$trial.num==14 |
                     long.data$trial.num==15 |
                     long.data$trial.num==16 |
                     long.data$trial.num==17 |
                     long.data$trial.num==18 |
                     long.data$trial.num==19 |
                     long.data$trial.num==20] <- 1;  
  }
  
  all.data <- rbind(all.data,long.data);
}

# Compute Day/Time of Each Hit for Excluding Multiples
all.data$day.and.time <- chron(dates = all.data$submit.date,
                               times = all.data$submit.time,
                               format=c("mon d y","h:m:s"))

# Sort Data Chronologically 
all.data <- all.data[with(all.data,order(subid,day.and.time)),]

# Find Subjects who Participated More Than Once
drop.subs <- all.data %.%
  select(subid) %.%
  group_by(subid) %.%
  summarise(drop = n() > 20) %.%
  filter(drop)

#Grab Earliest HIT for Each Participant
all.drops <- matrix(0,nrow(all.data))
for(subid in drop.subs$subid) {
  rows <- as.integer(all.data$subid == subid)
  all.drops[rows & (cumsum(rows) > 20)] <- 1
}
all.data <- subset(all.data,!all.drops)

#Convert Numerical Codes to English-Readable
all.data$subid <- as.factor(all.data$subid)
all.data[with(all.data,exp=="Exp 2" & trialType==2),"trialType"] <- 3
all.data$trialType <- factor(all.data$trialType, 
                             labels = c('Same','Switch','New Label'))

#Grab The Test Trials
test.data <- all.data[all.data$test ==  1,]
test.data$correct <- test.data$chosen == test.data$kept

# Grab Example Trials for Exclusion
example.data <- all.data[nchar(all.data$kept) == 0,]

# Exclude for Getting Examples Wrong
include.subs <- example.data %.%
  select(subid,chosen) %.%
  group_by(subid) %.%
  summarise(include = sum(chosen == "squirrel" | chosen=="tomato") == 4) %.%
  filter(include)

test.data <- merge(test.data,include.subs)
test.data <- test.data[with(test.data, order(exp,subid,trial.num)),]
keep.data <- filter(test.data,include)

# Mark the First Trial for Potential Analysis
keep.data$first.trial <- FALSE
keep.data$first.trial[with(keep.data,(interval==0 & trial.num==6) |
                             (interval==1 & trial.num==7) |
                             (interval==3 & keep.data$trial.num==9) |
                             (interval==7 & trial.num==13))] <- TRUE


keep.data$numPicN <- as.numeric(keep.data$numPic)
keep.data$intervalN <- as.numeric(keep.data$interval)

#Renumber subids for anonymity
keep.data$subid <- with(keep.data,
                        factor(subid,levels = subid[seq(1,length(subid),8)],
                               labels = 1:n.unique(subid)))

write.csv(filter(keep.data,exp=="Exp 1"),'../data/exp1_long.csv')
write.csv(filter(keep.data,exp=="Exp 2"),'../data/exp2_long.csv')
