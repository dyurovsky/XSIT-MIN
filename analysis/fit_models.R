rm(list=ls())

#get lab version of useful R functions
source('~/Projects/Other/Ranalysis/useful.R')

# Read in Experimental data from .csvs
# all.data: trial-level
# mss: subject-level
# ms: condition-level
source('load_from_csvs.R')

library(directlabels)
library(dplyr)
library(xtable)

library(rstan)
set_cppo("fast")  # for best running speed

# Data for Just Experiment 1
e1.mss <- filter(mss,exp=="Exp 1")
e1.ms <- filter(ms,exp=="Exp 1")

#### DESCRIPTIVES ######
mss <- all.data %.%
  group_by(intervalN,numPicN,trialType,subid,add=FALSE) %.%
  summarise(
    sums = sum(correct), 
    correct = mean(correct)) %.%
  arrange(trialType,intervalN,numPicN)

ms <- mss %.% #By-condition, across subjects
  summarise(
    prop = mean(correct),
    cih = ci.high(correct),
    cil = ci.low(correct),
    n = n()) %.%
  arrange(trialType,intervalN,numPicN)

#Add English-Readable columns for nicer graphing
ms$graph.trialType <- factor(ms$trialType,
                             labels=c("Same","Switch","New\nLabel"))
ms$graph.numPic <- sapply(ms$numPic,function(x){paste(x," Referents")})

###################### Set Up Data Structures for Model #######################
REF.TYPES <- 4
NUM.REFS <- c(2,3,4,8)
INT.TYPES <- 4
NUM.INTS <- c(1,2,4,8)
NUM.TRIALS <- 4
E1.COND.TYPES <- 2
COND.TYPES <- 3

# By-subject number correct in each condition
e1.ms_s <- filter(e1.mss,trialType=="Same")$sums

ms_s <- filter(mss,trialType=="Same")$sums
ms_o <- filter(mss,trialType=="Switch")$sums
ms_w <- filter(mss,trialType=="New Label")$sums

# Number of subjects in each condition
e1.n_s <- filter(e1.ms,trialType=="Same")$n

n_s <- filter(ms,trialType=="Same")$n
n_o <- filter(ms,trialType=="Switch")$n
n_w <- filter(ms,trialType=="New Label")$n

# Make a NumRefs x NumIntervals matrix for each condition
s.index.mat <- function (n,begin.flag) {
  if(begin.flag=="begin")
    mat <- head(c(0,cumsum(n))+1,-1)
  else
    mat <- cumsum(n)  
  mat <- t(matrix(mat,nrow=REF.TYPES,ncol=INT.TYPES))
 
}

e1.s_b <- s.index.mat(e1.n_s,'begin')
e1.s_e <- s.index.mat(e1.n_s,'end')

s_b <- s.index.mat(n_s,'begin')
s_e <- s.index.mat(n_s,'end')
o_b <- s.index.mat(n_o,'begin')
o_e <- s.index.mat(n_o,'end')
w_b <- s.index.mat(n_w,'begin')
w_e <-s.index.mat(n_w,'end')

# Package Experiment 1 up for RStan
e1.dat <- list(Trials = NUM.TRIALS,
               PicConds = REF.TYPES,
               IntConds = INT.TYPES,
               NumPic = NUM.REFS,
               Int = NUM.INTS,
               P_s = e1.n_s,
               P_o = n_o,
               B_s = e1.s_b,
               B_o = o_b,
               E_s = e1.s_e,
               E_o = o_e,
               S = e1.ms_s,
               O = ms_o)

# Package All data up for RStan
# Not Used in This Analysis, but in principle...
both.dat <- list(Trials = NUM.TRIALS,
                 PicConds = REF.TYPES,
                 IntConds = INT.TYPES,
                 NumPic = NUM.REFS,
                 Int = NUM.INTS,
                 P_s = n_s,
                 P_o = n_o,
                 P_w = n_w,
                 B_s = s_b,
                 B_o = o_b,
                 B_w = w_b,
                 E_s = s_e,
                 E_o = o_e,
                 E_w = w_e,
                 S = ms_s,
                 O = ms_o,
                 W = ms_w)

######################### Fit Parameters for Exp 1 First ######################

# Statistical Accumulation Model
fit.accumulation <- stan(file="../models/accumulation.e1.stan.R", 
                         data = e1.dat, iter = 1000, chains = 1)

# Single Referent Tracking Model
fit.singleref <- stan(file="../models/singleref.e1.stan.R", 
                      data = e1.dat, iter = 1000, chains = 1)

# Integrated Model
fit.integrated <- stan(file="../models/integrated.e1.stan.R", data = e1.dat, 
                       iter = 1000, chains = 1)

# Extract Inferred Parameter Values
accumulation.lambda <- mean(extract(fit.accumulation,"lambda")$lambda)
accumulation.gamma <- mean(extract(fit.accumulation,"gamma")$gamma)

singleref.lambda <- mean(extract(fit.singleref,"lambda")$lambda)
singleref.gamma <- mean(extract(fit.singleref,"gamma")$gamma)

integrated.sigma <- mean(extract(fit.integrated,"sigma")$sigma)
integrated.lambda <- mean(extract(fit.integrated,"lambda")$lambda)
integrated.gamma <- mean(extract(fit.integrated,"gamma")$gamma)

############################ Compute Fit Statistics ###########################
# Read in Helper Functions for Running Models
source('../models/model.helpers.R')

# Log Likelihoods and BICs
accumulation.lp <- mean(extract(fit.accumulation,"lp__")$lp__)
accumulation.bic <- bic(accumulation.lp,2,(sum(e1.n_s,n_o)))

singleref.lp <- mean(extract(fit.singleref,"lp__")$lp__)
singleref.bic <- bic(singleref.lp,2,(sum(e1.n_s,n_o)))

integrated.lp <- mean(extract(fit.integrated,"lp__")$lp__)
integrated.bic <- bic(integrated.lp,3,(sum(e1.n_s,n_o)))

# Compute Predictions for Each Model in Each Experimental Condition
# Used to Compute r^2
accumulation.ps <- NULL
singleref.ps <- NULL
integrated.ps <- NULL

for (i in 1:INT.TYPES) {
  for (n in 1:REF.TYPES) {
    
    # Compute Strength in Memory given Parameters
    accumulation.preds <- accumulation.pred(accumulation.gamma,
                                            accumulation.lambda,
                                            NUM.REFS[n],NUM.INTS[i])
    
    singleref.preds <- singleref.pred(singleref.gamma, singleref.lambda,
                                            NUM.REFS[n],NUM.INTS[i])
    
    integrated.preds <- integrated.pred(integrated.sigma, integrated.gamma,
                                        integrated.lambda,
                                        NUM.REFS[n],NUM.INTS[i])
    
    # Compute Prop. Correct given Strength in Memory
    accumulation.ps <- c(accumulation.ps,
                         prop.from.pred(accumulation.preds,NUM.REFS[n]))
    singleref.ps <- c(singleref.ps,
                      prop.from.pred(singleref.preds,NUM.REFS[n]))
    integrated.ps <- c(integrated.ps,
                       prop.from.pred(integrated.preds,NUM.REFS[n]))
  }
}

# Reorder Predictions to be Condition x NumRefs x Interval
nrows <- COND.TYPES*REF.TYPES*INT.TYPES
e1.nrows <-  E1.COND.TYPES*REF.TYPES*INT.TYPES
accumulation.ps <- c(accumulation.ps[seq(1,nrows,3)],
                     accumulation.ps[seq(2,nrows,3)],
                     accumulation.ps[seq(3,nrows,3)])

singleref.ps <- c(singleref.ps[seq(1,nrows,3)],
                  singleref.ps[seq(2,nrows,3)],
                  singleref.ps[seq(3,nrows,3)])

integrated.ps <- c(integrated.ps[seq(1,nrows,3)],
                   integrated.ps[seq(2,nrows,3)],
                   integrated.ps[seq(3,nrows,3)])

# Attach to experimental data
e1.ms$accumulation <- head(accumulation.ps,e1.nrows)
e1.ms$singleref <- head(singleref.ps,e1.nrows)
e1.ms$integrated <- head(integrated.ps,e1.nrows)

ms$accumulation <- accumulation.ps
ms$singleref <- singleref.ps
ms$integrated <- integrated.ps

# Compute Experiment 1 R2s
e1.accumulation.r2 <- with(e1.ms,cor(prop,accumulation))^2
e1.singleref.r2 <- with(e1.ms,cor(prop,singleref))^2
e1.integrated.r2 <- with(e1.ms,cor(prop,integrated))^2

accumulation.r2 <- with(ms,cor(prop,accumulation))^2
singleref.r2 <- with(ms,cor(prop,singleref))^2
integrated.r2 <- with(ms,cor(prop,integrated))^2

############################### Visualizations ################################
#TABLE S2

fit.table <- data.frame(Model=c("Statistical Accumulation", 
                                "Single Referent", "Integrated"))
fit.table$"Log Likelihood" <- c(accumulation.lp,singleref.lp,integrated.lp)
fit.table$BIC <- c(accumulation.bic,singleref.bic,integrated.bic)
fit.table$"E1 $r^{2}$" <- c(e1.accumulation.r2,e1.singleref.r2,
                            e1.integrated.r2)
fit.table$"E1+2 $r^{2}$" <- c(accumulation.r2,singleref.r2,
                              integrated.r2)

# Table of Model Fit Statistics (Paper Table 1)
print(xtable(fit.table,
             align = c("l","l","r","r","r","r"),
             label = "tab:model",
             digits = c(0,0,0,0,2,2)),
      include.rownames=FALSE,hline.after=c(0,nrow(fit.table)),
      sanitize.text.function=function(x){x})

# Graph Integrated Model Fit (Paper Figure 5)
quartz(width=10,height=4,title = "Integrated Model Model Fit")
ggplot(ms, aes(x=intervalN, y=prop, colour=trialType,
               label=graph.trialType)) +
  facet_grid( ~ graph.numPic) +
  geom_pointrange(aes(ymin = prop-cih,
                      ymax = prop+cih),
                  size = .8) +
  geom_hline(aes(yintercept=1/numPicN),lty=2)  +
  scale_y_continuous(limits = c(0,1),breaks=c(0,.25,.5,.75,1),
                     name = "Prop. Choosing Repeated Referent") +
  scale_x_continuous(limits=c(.9,10.3), breaks=c(1, 2, 4, 8),
                     name = "Intervening Trials") + 
  theme_bw(base_size=18) + 
  theme(legend.position="none",
        axis.title.y=element_text(hjust=.6,vjust=0.3)) +
  scale_color_manual(values=c("dark gray","#808080","#333333","#626262" )) +
  geom_line(aes(group=trialType,y=integrated,colour="dark gray"),
            lty=1,size=1.5) +
  geom_dl(method = list("last.qp",cex=1,hjust=-.15))

ggsave("e1and2fit.pdf",dpi=600,title="Two Process Model Fit")
