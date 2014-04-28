# clear all previous variables
rm(list=ls())

# Get Lab Version of Useful R Functions
source('~/Projects/Other/Ranalysis/useful.R')

# Read in Experimental data from .csvs
# all.data: trial-level
# mss: subject-level
# ms: condition-level
source('load_from_csvs.R')

# Load Libraries for Data Manipulation and Graphing
library(directlabels)
library(stringr)
library(dplyr)
library(xtable)

e1.data <- filter(all.data,exp=="Exp 1")
e2.data <- filter(all.data,exp=="Exp 2")

#### LMERS ######
e1.3way.lmer <- glmer(correct ~ log.numPic * log.interval * 
                                      trialType + (trialType | subid), 
                                    family="binomial", 
                                    data = e1.data)

e2.2way.lmer <- glmer(correct ~ log.numPic * trialType + 
                        log.interval * trialType + (trialType | subid), 
                      family="binomial", 
                      data = e2.data)
#### TABLES ####

#TABLE S1
e1.tab <- as.data.frame(summary(e1.3way.lmer)$coef)
e1.tab$Predictor <- c("Intercept","Log(Referents)","Log(Interval)",
                      "Switch Trial","Log(Referents)*Log(Interval)",
                      "Log(Referents)*Switch Trial",
                      "Log(Interval)*Switch Trial",
                      "Log(Referents)*Log(Interval)*Switch Trial")
rownames(e1.tab) <- NULL
e1.tab <- e1.tab[,c(5,1:4)]
e1.tab$stars <- sapply(e1.tab[,5],getstars)
names(e1.tab)[6] <- ""

names(e1.tab)[4:5] <- c("$z$ value","$p$ value")

print(xtable(e1.tab,
             align = c("l","l","r","r","r","r","l"),
             label = "tab:exp1_reg"),
      include.rownames=FALSE,hline.after=c(0,nrow(e1.tab)),
      sanitize.text.function=function(x){x})
         

#TABLE S2
e2.tab <- as.data.frame(summary(e2.2way.lmer)$coef)
e2.tab <- e2.tab[c(1,2,4,3,5,6),] #reorder for consistency
e2.tab$Predictor <- c("Intercept","Log(Referents)","Log(Interval)",
                      "New Label Trial","Log(Referents)*New Label Trial",
                      "Log(Interval)*New Label Trial")
rownames(e2.tab) <- NULL
e2.tab <- e2.tab[,c(5,1:4)]
e2.tab$stars <- sapply(e2.tab[,5],getstars)
names(e2.tab)[6] <- ""

names(e2.tab)[4:5] <- c("$z$ value","$p$ value")

print(xtable(e2.tab,
             align = c("l","l","r","r","r","r","l"),
             label = "tab:exp2_reg"),
      include.rownames=FALSE,hline.after=c(0,nrow(e2.tab)),
      sanitize.text.function=function(x){x})

#### STATS #########

## CHI-SQUARE TESTS##
#get a distribution on number of correct responses in each condition
ss <- aggregate(sums ~ exp + trialType + numPicN + intervalN , data=mss,
                FUN = function(x) {hist(x,breaks=c(-1,0,1,2,3,4),
                                        plot=FALSE)$counts})
ss$n <- summarise(mss,n=n())$n


#get the null distribution expected under a binomial guessing model
nulls <- sapply(c(1/2,1/3,1/4,1/8),
                FUN=function(x){dbinom(c(0,1,2,3,4),4,x)})
nulls <- matrix(nulls,nrow=16*ncol(nulls),ncol=nrow(nulls),byrow=TRUE)

for(i in 1:nrow(ss)){
  #simulate.p=TRUE is probably better but doesn't produce dfs.
  #sumulate.p=FALSE generates warnings, but produces all of the same results
  chisq <- chisq.test(ss[i,5],p=nulls[i,],simulate.p=FALSE)
  ss$chisq[i] <- chisq$statistic
  ss$pval[i] <- chisq$p.value
  ss$df[i] <- chisq$parameter

}

ms$graph.trialType <- factor(ms$trialType,
                             labels=c("Same","Switch","New\nLabel"))
ms$graph.numPic <- sapply(ms$numPic,function(x){paste(x," Referents")})


quartz(width=10,height=5.5,title = "Experiment 1 Data")
ggplot(ms, aes(x=intervalN, y=prop, colour=trialType,label=graph.trialType))+
  facet_grid(expt ~ graph.numPic) +
  geom_pointrange(aes(ymin = prop-cih,
                      ymax = prop+cih),
                  size = .8) +
  geom_hline(aes(yintercept=1/numPicN),lty=2)  +
  scale_y_continuous(limits = c(0,1),breaks=c(0,.25,.5,.75,1),
                     name = "Prop. Choosing Repeated Referent") +
  scale_x_continuous(limits=c(.9,10.3), breaks=c(1, 2, 4, 8),
                     name = "Intervening Trials") + 
  theme_bw(base_size=18) +
  theme(legend.position="none", axis.title.y=element_text(vjust=0.27)) +
  #scale_color_manual(values=c("#F8766D","#00BA38","#619CFF")) 
  scale_color_grey(start=.2,end=.5) +
  geom_dl(method=list("last.qp",cex=1,hjust=-.15)) 

ggsave("exp1_2_data.pdf",dpi=600,title="Experiments 1 and 2")


## WILCOXON TESTS##
ws <- as.data.frame(aggregate(correct ~ numPicN + intervalN + trialType , data=mss,
                              FUN=function(x) {wilcox.test(x,mu=.125, exact=FALSE, 
                                                           correct=TRUE)})$correct)

ms$W <- ws$statistic
ms$p <- ws$p.value

ts <- ddply(mss, .(numPicN,intervalN,trialType), function(x){t.test(x$correct,1/x$numPic)$statistic})
ts$ps <- ddply(mss, .(numPicN,intervalN,trialType), function(x){t.test(x$correct,1/x$numPic)$p.value})$V1



