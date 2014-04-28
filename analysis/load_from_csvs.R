####################### Read in Data ######################
e2 <- read.csv("../data/exp1_long.csv")
e3 <- read.csv("../data/exp2_long.csv")

e2$exp <- "Exp 1"
e3$exp <- "Exp 2"

d <- merge(e2,e3,all=TRUE)
d$expt <- factor(d$exp)

#### DATA MUNGING ######
all.data <- d %.%
  select(exp,subid,trialType,correct,numPicN,intervalN) %.%
  mutate(
    trialType = factor(trialType,levels=c("Same","Switch","New Label")),
    intervalN = intervalN + 1,
    log.numPic = log2(numPicN),
    log.interval = log2(intervalN)) %.%
  group_by(exp,intervalN,numPicN,trialType,subid)

#### DESCRIPTIVES ######
mss <- all.data %.%
  summarise(
    sums = sum(correct), 
    correct = mean(correct)) %.%
  arrange(exp,trialType,intervalN,numPicN)

ms <- mss %.% #By-condition, across subjects
  summarise(
    prop = mean(correct),
    cih = ci.high(correct),
    cil = ci.low(correct),
    n = n()) %.%
  arrange(exp,trialType,intervalN,numPicN)

#Add English-Readable columns for nicer graphing
ms$graph.trialType <- factor(ms$trialType,
                             labels=c("Same","Switch","New\nLabel"))
ms$graph.numPic <- sapply(ms$numPic,function(x){paste(x," Referents")})