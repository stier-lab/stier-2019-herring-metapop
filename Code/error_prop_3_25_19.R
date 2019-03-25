#load a bunch of packages 
library(devtools)
library(ggplot2)
library(reshape2)
library(tidyr)
library(gdata)
library(reshape2)
library(Hmisc)
library(coda)
library(cowplot)
library(here)
library(scales)
library(ggpubr)


## ========================Load Model Output================================##
load("diag_equal_design_c_noUsig_2q_pcnorm2.RData")


createMcmcList = function(model) {
  McmcArray = as.array(model$BUGSoutput$sims.array)
  McmcList = vector("list",length=dim(McmcArray)[2])
  for(i in 1:length(McmcList)) McmcList[[i]] = as.mcmc(McmcArray[,i,])
  McmcList = mcmc.list(McmcList)
  return(McmcList)
}


myList<-createMcmcList(model) #mcmc output 



## ========================years, sites, and mcmcm chains================================##

years = seq(1950,2015)
nYears = length(years)
nSites = 11

n.chains = 4
n.burnin = 500000
n.thin = 500
n.iter = 1000000

#number of recorded mcmc
runL <- n.chains*(n.iter-n.burnin)/n.thin

## ========================PDO TIME SERIES ================================##


#load PDO raw data
pdo<-read.csv(here("data","pdo.csv"))
pdosummer<-subset(pdo,month %in% c(3,4,5,6)) #average march<-june
pdoxb<-c(tapply(pdosummer$Value,list(pdosummer$year),mean))
pdo2<-pdoxb[97:162] #1940-2015

pdosummer<-subset(pdo,month %in% c(3,4,5,6)) #average april<-june
pdoxb<-c(tapply(pdosummer$Value,list(pdosummer$year),mean))
pdo2<-pdoxb[87:162] #1940-20135
pdo3<-pdo2[11:76] #1950-2015
plot.ts(pdo3)

pdocoef<-melt(model$BUGSoutput$sims.array[,,"pdocoef"])$value

hist(pdocoef)

distpdoci<-smedian.hilow(pdocoef,conf.int=0.95)
distpdoiqr<-smedian.hilow(pdocoef,conf.int=0.5)

pdodf<-data.frame("pdo"=pdo3,
                  "medianci"=distpdoci[1]*pdo3,
                  "maxci"=distpdoci[3]*pdo3,
                  "minci"=distpdoci[2]*pdo3,
                  "maxiqr"=distpdoiqr[3]*pdo3,
                  "miniqr"=distpdoiqr[2]*pdo3
)

pdodf$year<-as.numeric(rownames(pdodf))
pdodf2<-data.frame(site=rep("pdocoef",nrow(pdodf)),pdodf)
pdodf2$hc<-ifelse(pdodf2$medianci>0,"hot","cold")


pdotsgg<-ggplot(pdodf2,aes(x=year,y=medianci))+
  geom_hline(yintercept=0)+
  geom_ribbon(aes(ymin=minci,ymax=maxci),fill="grey70",alpha=0.2)+
  geom_ribbon(aes(ymin=miniqr,ymax=maxiqr),fill="grey")+
  geom_line(lty=2)+
  geom_point(aes(colour=hc))+
  theme_pubr()+
  ylab("PDO Effect")+
  xlab("Year")+
  scale_x_continuous(breaks=bk)


print(pdotsgg)






## ========================Spawn and Fishing TIME SERIES ================================##

x <- read.csv(here("data","HG_Spawn_Survey_1940_2015.csv"))
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


c <- read.csv(here("data","herring_catch_local2015.csv"))
yr_c<-unique(c$Year)
c$presabs <-ifelse(c$TotalCatch>0,1,0)

unique(c[,c('Section','Name')]) #look at sections and names

c<-drop.levels(subset(c,Section %in% c(1,2,3,5,6,12,21,22,23,24,25)))#subset out Cartwright Sound (4), Masset (11)
ctab<-data.frame(tapply(c$CatchJan_Apr,list(c$Year,c$Name),sum)) #just spring catch

data.frame(colnames(Y),colnames(ctab)) #mismatch column names
ctab2<-ctab[,c(11,7,8,2,5,6,3,9,1,4,10)] #re order so catch table matches spawn table 
data.frame(colnames(logSHI),colnames(ctab2)) #double check column orders match

colnames(Y)<-colnames(ctab2)

ctab2<-as.matrix(ctab2)
logcatch<-log(ctab2+1)



## ========================Biomass by subpopulation witrh error propagation ================================##


xpop<-melt(model$BUGSoutput$sims.array[,,"X"])


#plot out individual stocklets with halo confidence intervals 


edf1<-melt(model$BUGSoutput$sims.array[,,colnames(myList[[1]])[which(colnames(myList[[1]])=="X[1,1]"):which(colnames(myList[[1]])=="X[66,11]")]])
names(edf1)<-c("num","chain","response","value")
edf1$chain<-factor(edf1$chain)

edf1$group<-factor(sort(rep(seq(1:nSites),runL)))
edf1$year2<-sort(rep(1950:2015,runL*nSites))
edf1$year<-sort(rep(seq(1:nYears),runL*nSites))
edf1$section<-sort(rep(colnames(Y),nYears*runL))


tsmed<-tapply(edf1$value,list(edf1$response),median)
tsupper<-melt(tapply(edf1$value,list(edf1$response),quantile,probs=c(0.95)))
tslower<-melt(tapply(edf1$value,list(edf1$response),quantile,probs=c(0.05)))

xdat<-data.frame(exp(tsmed),exp(tsupper[,2]),exp(tslower[,2]))
colnames(xdat)<-c("median","upper","lower")
xdat$year<-rep(1950:2015,nSites)

secvec<-colnames(Y)
sec<-c()

for(i in 1:11){
  temp<-rep(secvec[i],66)
  sec<-c(sec,temp)
}

xdat$section<-sec

ggplot(xdat,aes(x=year,y=median))+
  geom_hline(yintercept=0)+
  geom_ribbon(aes(ymin=lower,ymax=upper),fill="dodgerblue",alpha=0.2)+
  # geom_ribbon(aes(ymin=miniqr,ymax=maxiqr),fill="grey")+
  geom_line(colour="grey")+
  # geom_point(aes(colour=hc))+
  # theme_acs()+
  ylab("Estimated Biomass")+
  xlab("Year")+
  facet_wrap(~section,scales="free")+
  scale_x_continuous(breaks=bk)+
  scale_y_continuous(labels=comma)




## ========================Fishing Figures ================================##


#average estimated proportion caught by fishing 
tpc<-model$BUGSoutput$mean$Pc


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


emat <- matrix(0,nrow=nYears,ncol=nSites)

for(i in 1:156){
  
  emat[INDEX[i,1],INDEX[i,2]]<-tpc[i]
  
}


colnames(emat)<-colnames(Y)
pc_tab<-melt(emat)
colnames(pc_tab) <- c("year2","section","pc")
pc_tab$year<-seq(1950,2015,1)

arch<-data.frame(tapply(pc_tab$pc,list(pc_tab$year),mean))
arch$year<-as.numeric(rownames(arch))
names(arch)<-c("pc","year")

pc_tab2<-subset(pc_tab,pc>0)

sec<-data.frame(tapply(pc_tab2$pc,list(pc_tab2$year),mean))
sec$year<-as.numeric(rownames(sec))
names(sec)<-c("pc","year")


fish<-merge(arch,sec,by=("year"),all=TRUE)
colnames(fish)<-c("year","archipelago","stocklet")
fish[is.na(fish)]<-0

fish2<-melt(fish,id.vars="year")
colnames(fish2)<-c("year","var","pc")

bk<-seq(1950,2015,by=10)

#Arhipelago Versus Stocklet 

ggplot(fish2,aes(x=year,y=pc))+
  geom_line(aes(lty=var,colour=var))+
  theme_pubr()+
  scale_y_continuous(limits=c(0,1),breaks=c(0,0.2,0.4,0.6,0.8,1.0))+
  scale_x_continuous(breaks=bk)+
  ylab("Proportion Biomass Caught (F)")


#pull out upper and lower quantiles for fishing x year combinations for 156 places where fishing happened

#first Identify Which Row-Column Combinations Where Catch>0
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


model$BUGSoutput$sims.array[,,colnames(myList[[1]])[which(colnames(myList[[1]])=="Pc[1]"):which(colnames(myList[[1]])=="Pc[156]")]]

edf1<-melt(model$BUGSoutput$sims.array[,,colnames(myList[[1]])[which(colnames(myList[[1]])=="Pc[1]"):which(colnames(myList[[1]])=="Pc[156]")]])
names(edf1)<-c("num","chain","response","value")
edf1$chain<-factor(edf1$chain)

#estimate means for each of the diffrent Pcs
pcmed <-tapply(edf1$value,list(edf1$response),median)
pcupper<-tapply(edf1$value,list(edf1$response),quantile,probs=c(0.95))
pclower<-tapply(edf1$value,list(edf1$response),quantile,probs=c(0.05))


#propogate the medians and quantiles into three separate matrices in their right row/column combos
medINDEX <- Catch.dummy
upperINDEX <- Catch.dummy
lowerINDEX <- Catch.dummy

for(i in 1:156){
  medINDEX[INDEX[i,1],INDEX[i,2]] <- pcmed[i]
  upperINDEX[INDEX[i,1],INDEX[i,2]] <- pcupper[i]
  lowerINDEX[INDEX[i,1],INDEX[i,2]] <- pclower[i]
  
}


fdf<-  data.frame(
  melt(medINDEX),
  melt(upperINDEX)$value,
  melt(lowerINDEX)$value
)

names(fdf) <- c("year","site","median","upper","lower")


ggplot(fdf,aes(x=year,y=median))+
  geom_line()+
  geom_ribbon(aes(ymin=lower,ymax=upper),fill="blue",alpha=0.2)+
  facet_wrap(~site)+
  theme_pubr()+
  ylab("Proportion Caught")

fdf_fe<-subset(fdf,fdf$median>0)

sort(tapply(fdf_fe$median,list(fdf_fe$site),mean))*100

#propogate error from iterations of MCMC 
fa <- model$BUGSoutput$sims.array[,,colnames(myList[[1]])[which(colnames(myList[[1]])=="Pc[1]"):which(colnames(myList[[1]])=="Pc[156]")]]
fdf2 <- melt(fa)
names(fdf2)<- c("iter","chain","site","pc")
fdf2 <- spread(fdf2, key = site, value = pc)

fnames<-unlist(dimnames(fa))
fnames2<-c("iter","chain",fnames)
fdf3 <- fdf2[, fnames2]


#insert mean of each of the values into a matrix at right location

#melt to long form
#estiamte the mean and sd across all sites including zeros and excluding zeros
#save to a new vector of means and sds and then summarize with figure 

archmean <- Catch.dummy
emat <-matrix(nrow=66,ncol=4000) #matrix for catch at arhcipelago per year
emat2 <-matrix(nrow=66,ncol=4000) #matrix for catch average per year excluding sites where no catch occurred

for(i in 1:4000){
  
  evec <- fdf3[i,][-c(1:2)] 
  
  for(j in 1:156){
    archmean[INDEX[j,1],INDEX[j,2]] <- unlist(c(evec[j]))
  }
  
  temp <- melt(archmean) #including sites that didn't get fishes 
  names(temp) <- c("year","site","pc")
  emat[,i] <- tapply(temp$pc,list(temp$year),mean)
  
  temp2 <- subset(temp,pc>0) #excluding sites that didn't get fished 
  yr <- data.frame("year" = 1950:2015)
  t<-tapply(temp2$pc,list(temp2$year),mean) 
  y<-as.numeric(names(t))
  tdf <- data.frame("pc"=t,"year"=y)
  tdf2 <-merge(tdf,yr, by = "year",all=TRUE)
  emat2[,i] <- tdf2$pc
  
}


rownames(emat)<-1950:2015
memat <- melt(emat)

tapply(memat$value,list(memat$Var1),mean)

pcmed <-tapply(memat$value,list(memat$Var1),median)
pcupper<-tapply(memat$value,list(memat$Var1),quantile,probs=c(0.95))
pclower<-tapply(memat$value,list(memat$Var1),quantile,probs=c(0.05))

fdf4 <- data.frame("median"=pcmed,"upper"=pcupper,"lower"=pclower,"year"=1950:2015)

ggplot(fdf4,aes(x=year,y=median))+
  geom_line()+
  geom_ribbon(aes(ymin=lower,ymax=upper),fill="blue",alpha=0.2)+
  theme_pubr()


rownames(emat2)<-1950:2015
memat2 <- melt(emat2)
memat2[is.na(memat2)]<-0

pcmed <-tapply(memat2$value,list(memat2$Var1),median)
pcupper<-tapply(memat2$value,list(memat2$Var1),quantile,probs=c(0.95),na.rm=TRUE)
pclower<-tapply(memat2$value,list(memat2$Var1),quantile,probs=c(0.05),na.mr=TRUE)

fdf5 <- data.frame("median"=pcmed,"upper"=pcupper,"lower"=pclower,"year"=1950:2015)

ggplot(fdf4,aes(x=year,y=median))+
  geom_line()+
  geom_ribbon(aes(ymin=lower,ymax=upper),fill="blue",alpha=0.2)+
  theme_pubr()


pcmerge<-merge(fdf4,fdf5, by = "year",all=TRUE)


colnames(pcmerge) <- c("year","archomed","archupper",
                       "archlower","submed","subupper","sublower")


#*** Fig 4b 
ggplot(pcmerge)+
  geom_line(aes(x=year,y=archomed))+
  geom_ribbon(aes(x=year,ymin=archlower,ymax=archupper),fill="blue",alpha=0.2)+
  geom_line(aes(x=year,y=submed),lty=2)+
  geom_ribbon(aes(x=year,ymin=sublower,ymax=subupper),fill="red",alpha=0.2)+
  theme_pubr()+
  xlab("Year")+
  ylab("Proportion Biomass Caught")+
  scale_y_continuous(limits=c(0,0.6),breaks=c(0,0.2,0.4,0.6))+
  scale_x_continuous(breaks=seq(1950,2020,by=10))


