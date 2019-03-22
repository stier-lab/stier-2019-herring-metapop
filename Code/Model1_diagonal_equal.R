## ==================
##  Load Libraries
## ==================

library(ggplot2)
library(reshape2)
library(R2jags)
library(coda)
library(gridExtra)
library(MASS)
library(Hmisc)
library(MCMCvis)

setwd("~/github/herring_synchrony_hg/Code")
source('multiplot.R')
source('theme_acs.R')


#test 

## ==================
##  Time and Sites
## ==================
years = seq(1950,2015)
nYears = length(years)
nSites = 11

## ==================
##  Load PDO time series
## ==================

pdo<-read.csv("pdo.csv") 
pdosummer<-subset(pdo,month %in% c(3,4,5,6)) #average march<-june
pdoxb<-c(tapply(pdosummer$Value,list(pdosummer$year),mean))
pdo2<-pdoxb[97:162] #1940-2015

## ==================
##  Load Spawn Index time series
## ==================

x <- read.csv("HG_Spawn_Survey_1940_2015.csv")
x <- x[c(1,3,13,14,15,16)]
x$presabs <-ifelse(x$SHI>0,1,0)

x2 <- x[,c(1,2,4)]
w <- reshape(x2, 
             timevar   = "section_name",
             idvar     = c("year"),
             direction = "wide")[,-1]

w[w==0] <- NA#replace zeros with NAs
Y= as.matrix(w)
Y = Y[-c(1:10),-c(4,7)] #drop site 4 and 7
Y = log(Y)
logSHI<-Y

ym<-melt(Y)
ym<-data.frame(ym,rep(seq(1:nrow(Y)),ncol(Y)))
names(ym)<-c("crap","site","logSHI","time")

## ==================
##  Load Catch time series
## ==================

c <- read.csv("herring_catch_local2015.csv")
yr_c<-unique(c$Year)
c$presabs <-ifelse(c$TotalCatch>0,1,0)

unique(c[,c('Section','Name')]) #look at sections and names

c<-drop.levels(subset(c,Section %in% c(1,2,3,5,6,12,21,22,23,24,25)))#subset out Cartwright Sound (4), Masset (11)
ctab<-data.frame(tapply(c$CatchJan_Apr,list(c$Year,c$Name),sum)) #just spring catch

data.frame(colnames(Y),colnames(ctab)) #mismatch column names
ctab2<-ctab[,c(11,7,8,2,5,6,3,9,1,4,10)] #re order so catch table matches spawn table 
data.frame(colnames(logSHI),colnames(ctab2)) #double check column orders match

colnames(logSHI)<-colnames(ctab2)

ctab2<-as.matrix(ctab2)
logcatch<-log(ctab2+1)

## ==================
##  Catch Priors and Design matrix 
## ==================

#make a prior matrix that sets the upper bound of the uniform distribution for prior Pc to be 0.95 when it's a non-zero and 0.01 when it's zero
nonz        <- which(logcatch>0)
c_prior_sim <- matrix(0.95,nrow=nrow(logcatch),ncol=ncol(logcatch))
design.c    <- ifelse(logcatch>0,1,0)

## ==================
## Estimating Catch Rates Only in Sites where Catch Reported
## ==================

#Identify Which Row-Column Combinations Where Catch>0
INDEX      <- NULL  
INDEX.zero <- NULL  
for(i in 1:nrow(logcatch)){
  
  if(length(logcatch[i,][logcatch[i,]>0])>0){ 
    temp   <- data.frame(row = rep(i,length(which(logcatch[i,]>0))),col = which(logcatch[i,]>0))
    INDEX <- rbind(INDEX,temp)
  }
  
  if(length(logcatch[i,][logcatch[i,]==0])>0){ 
    temp.2 <- data.frame(row = rep(i,length(which(logcatch[i,]==0))),col = which(logcatch[i,]==0))
    INDEX.zero <- rbind(INDEX.zero,temp.2)
  }
}

nIndex      <- nrow(INDEX) #nubmer of row-column combinations with catch>0
nIndex.zero <- nrow(INDEX.zero) #nubmer of row-column combinations with catch>0

Catch.dummy <- logcatch
Catch.dummy[Catch.dummy > 0] <- 1 #put 1 in places where >0 catch

##read in INDEX, nIndex, and Catch.dummy as data in jags.data


###read in a vector of q to index two q's
#surface 1950-1987 
length(1950:1987)
#scuba 1988-present 
length(1988:2015)

q_idx<-c(rep(1,38),rep(2,28))


# q_idx<-matrix(0,ncol=2,nrow=nYears)
# q_idx[c(rep(1,38),rep(0,28)),1]
# q_idx[,2]<-c(rep(0,38),rep(1,28))

## ==================
## JAGS CODE to fit model
## ==================

#Begin JAGS code
jagsscript = cat("
                 model {  
                
##########################
#OBSERVATION MODEL PRIORS and LIKELIHOODs
##########################


tauR2  ~ dgamma(2,2); #this is the estimated variance prior for Spawn Index to Spawn Biomass Conversion

#Catchability Conversion from Spawn Index to Spawn Biomass 

for(q in 1:2) {
  log.q[q] ~ dnorm(0,10) 
  }


# Prop. Catch (pc) Prior  (USES INDEX to find row-column combinations) for prior 

for(k in 1:nIndex){
         Pc.logit[k] ~ dnorm(-1.38625,2)
         Pc[k] = exp(Pc.logit[k]) / (exp(Pc.logit[k]) + 1)

}
                 
for(i in 1:nYears) {
  for(j in 1:nSites) {
    Y[i,j] ~  dnorm((X[i,j]+log.q[q_idx[i]]),1/tauR2); #obs eq for spawn index. a normal dist mean Xq LIKELIHOOD STATEMENT*
  }
}


for(k in 1:nIndex){
    Pc.mat[INDEX[k,1],INDEX[k,2]] <-  Pc[k] 
    tl[INDEX[k,1],INDEX[k,2]]     <- Z[INDEX[k,1],INDEX[k,2]] + log(Pc[k])   #total catch
    ctab[INDEX[k,1],INDEX[k,2]]   ~ dnorm((tl[INDEX[k,1],INDEX[k,2]]),1/0.0001) #distribution of catch data  LIKELIHOOD STATEMENT ****  
}

for(k in 1:nIndex.zero){
    Pc.mat[INDEX.zero[k,1],INDEX.zero[k,2]] <- 0
}


##########################
#PROCESS MODEL PRIORS
##########################

# POPULATION GROWTH parameters U: 
Umu   ~  dnorm(0,1); #average population growth

#COVARIATE
pdoz2   <- 1
pdocoef ~ dnorm(0,1/pdoz2); #estimated impact of pdo of herring

#VARINCE MATRIX - Diagonal and Equal - same variance all sites 

# For indepenent and equal variances
sigma2 ~ dgamma(0.001,0.001); #set initial


for(j in 1:nSites) {
  sigma2.all[j] <- sigma2; #replicate initial value for all sigma2
}
    

# Estimate the initial state vector of population abundances

for(j in 1:nSites) {
      Z[1,j] ~ dnorm(5,1/100); # vague normal prior for first time step #changed from inverse 3/31
      X[1,j] <- Z[1,j] + log(1-Pc.mat[1,j]) # changed to reflect handwritten model 3/31
  }


#  Delta Prior

for(i in 1:(nYears)){
  for(j in 1:nSites) {
    delta[i,j] ~ dnorm(0,1/sigma2.all[j]) #changed from inverse 3/31
  }
}

      
for(i in 2:nYears){
  for(j in 1:nSites) {

  #fishing by site and global pdo estimates and U estimates 
  Z[i,j] <- X[(i-1),j]+Umu+pdocoef*pdo[i-1]+delta[(i-1),j]

  #Estimate the state X with catch by site 
  X[i,j] <-Z[i,j]+log(1-Pc.mat[i,j]);

 }
}

} 
 ",file="diag_equal_design_c_noUsig_2q_pc_norm.RData.txt")


#data going into the model
jags.data = list("Y"=logSHI, "nYears"=nYears,"nSites"=nSites,"pdo"=pdo2,"ctab"=logcatch,
                 "INDEX"=INDEX,"INDEX.zero"=INDEX.zero, "nIndex"=nIndex,"nIndex.zero"=nIndex.zero,"q_idx"=q_idx) # named list

jags.params=c("X","sigma2","Umu","delta","Z","pdocoef","log.q","Pc","tauR2")


model.loc="diag_equal_design_c_noUsig_2q_pc_norm.RData.txt" # name of the txt file

n.chains = 4
n.burnin = 500000
n.thin = 500
n.iter = 1000000

#number of recorded mcmc
runL <- n.chains*(n.iter-n.burnin)/n.thin

inits = NULL

for(i in 1:n.chains){
 inits[[i]]    <- list(
   "Umu" = rnorm(1,0,0.1),
   #"Usig2" = runif(1,0,1),
   "pdocoef" = rnorm(1,0,0.1),
   #"log.q" = rnorm(1,0,0.1),
   "tauR2" = runif(1,0,1)
   ) 
}

model = jags(jags.data, inits=inits,parameters.to.save=jags.params,
            model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin, n.thin=n.thin, n.iter=n.iter, DIC=TRUE)

#attach.jags(model)

mDIC<-model$BUGSoutput$DIC 
pD<-model$BUGSoutput$pD
devi<-model$BUGSoutput$deviance

## ==================
## Save/Load model output
## ==================

setwd("~/Dropbox/Projects/In Progress/Herring_Haida_Gwaii/Code")
save("model","mDIC","inits",
     file="diag_equal_design_c_noUsig_2q_pcnorm2.RData")

setwd("~/Dropbox/Projects/In Progress/Herring_Haida_Gwaii/Herring-GitHub/Herring-Metapopulation-HG/Code")
load("diag_equal_design_c_noUsig_2q_pcnorm2.RData")

write.csv(MCMCsummary(model),"AllParams_pc_norm.csv")

#subset out which didn't converge
# tdf<- data.frame(MCMCsummary(model))
# write.csv(tdf,"AllParams_pc_norm.csv")
# tdf2<- subset(tdf, Rhat > 1.1)
# write.csv(tdf2,"AllParams_pc_norm_problems.csv")

#plot chains of poorly converged 
# tvec<- rownames(tdf2)
# MCMCtrace(model,params=c('Pc\\[6\\]'),ISB=FALSE)
# MCMCtrace(model,params=c('Pc\\[43\\]'),ISB=FALSE)
# MCMCtrace(model,params=c('X\\[2\\]'),ISB=FALSE)


MCMCsummary(model,params='sigma2',n.eff=TRUE)
MCMCpstr(model,params='Umu',func=mean)
MCMCpstr(model,params='Umu',func=function(x) quantile(x,probs=c(0.025,0.975)))
MCMCplot(model,params='Umu')

jags.params=c("X","sigma2","Umu","delta","Z","pdocoef","log.q","Pc","tauR2")


MCMCtrace(model,params='Pc',filename="Pc_chains_posteriors.pdf")
MCMCtrace(model,params='delta',filename="delta_chains_posteriors.pdf")
MCMCtrace(model,params='X',filename="X_chains_posteriors.pdf")
MCMCtrace(model,params='Umu',filename="Umu_chains_posteriors.pdf")
MCMCtrace(model,params='sigma2',filename="Sigma2_chains_posteriors.pdf")
MCMCtrace(model,params='pdocoef',filename="pdocoef_chains_posteriors.pdf")
MCMCtrace(model,params='tauRr2',filename="tauR2_chains_posteriors.pdf")
MCMCtrace(model,params='log.q',filename="logq_chains_posteriors.pdf")




