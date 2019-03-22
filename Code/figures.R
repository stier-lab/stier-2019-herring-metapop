## ==================
##  Load relevant packages
## ==================
library(devtools)
library(ggplot2)
library(reshape2)
library(tidyr)
library(gdata)
library(reshape2)
library(Hmisc)
library(coda)
library(cowplot)
library(ecofolio)
library (synchrony)
library(here)
library(scales)

# install.packages(c("plyr", "reshape", "MuMIn", "robustbase", "devtools"))
# devtools::install_github("ecofolio", username="seananderson")


## ========================LOAD Relevant Code================================##

# setwd("~/Dropbox/Projects/In Progress/Herring_Haida_Gwaii/Herring-GitHub/Herring-Metapopulation-HG/Code")
source('code/theme_publication.R')
source('code/multiplot.R')
  load("code/diag_equal_design_c_noUsig_2q_pcnorm2.RData")
  

## ========================LOAD Relevant Functions================================## 

createMcmcList = function(model) {
  McmcArray = as.array(model$BUGSoutput$sims.array)
  McmcList = vector("list",length=dim(McmcArray)[2])
  for(i in 1:length(McmcList)) McmcList[[i]] = as.mcmc(McmcArray[,i,])
  McmcList = mcmc.list(McmcList)
  return(McmcList)
}

stat_sum_single <- function(fun, geom="point", ...) {
  stat_summary(fun.y=fun, geom=geom, size = 1, ...)
}

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

# Function to order correlations by strength of correlation
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

bk<-seq(1950,2015,by=10)


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

## ========================LOAD Relevant Data================================##

#reconstruct MCMC ouput
myList<-createMcmcList(model) #mcmc output 

## ==================
##  Time and Sites
## ==================
years = seq(1950,2015)
nYears = length(years)
nSites = 11


## ==================
##  Load PDO time series and subset to 1940-2015
## ==================
pdo<-read.csv("data/pdo.csv") 
pdosummer<-subset(pdo,month %in% c(3,4,5,6)) #average march<-june
pdoxb<-c(tapply(pdosummer$Value,list(pdosummer$year),mean))
pdo2<-pdoxb[97:162] #1940-2015


## ==================
##  Load Spawn Index time series
## ==================
x <- read.csv("data/HG_Spawn_Survey_1940_2015.csv")
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
c <- read.csv("data/herring_catch_local2015.csv")
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

## ========================FIGURES================================##


## ==================
##  Fig 1 Spawn Index At Different Sites***
## ==================

s1<-subset(x, year > 1949)
is.na(s1$SHI)<-!s1$SHI #change zero to NA, which is accurate 
s2<-subset(s1,section_name %in% c("Louscoone Inlet","Juan Perez Sound","Rennell Sound","Skidegate Inlet","Skincuttle Inlet","Englefield Bay"))


ggplot(s2,aes(x=year,y=SHI,group=section_name))+
  geom_line(lty=1)+
  geom_point(shape=21,fill="white")+
  #scale_y_log10()+
  # theme_acs()+
  #ggtitle("Herring By Spawn Area")+
  facet_wrap(~section_name,scales="free_y")+
  theme(legend.position="none")+
  scale_x_continuous(limits=c(1950,2014),breaks=c(1950,1970,1990,2010))+
  # scale_y_log10(limits=c(100,1000000),
  #               breaks=c(10^2,10^4,10^6),
  #               labels = trans_format("log10", math_format(10^.x)))+  
  scale_y_continuous(labels=comma)+
  theme_Publication()+
  ylab("Spawn Index")+
  xlab("Year")

## ==================
##  Fig S1 Time series of data availability 
## ==================

s3<-subset(s1,section_name !="Cartwright Sound")
s3<-drop.levels(subset(s3,section_name !="Masset Inlet"))

unique(s3$section_name)

plot.spawn<-ggplot(s3,aes(x=year,y=section_name,fill=factor(presabs)))+
              geom_tile(colour="black")+
              scale_fill_manual(values=c("#0C5790","#BC3C05"))+
              theme_Publication()+
              theme(legend.position="none")+
              scale_x_continuous(limits=c(1950,2015),breaks=c(1950,1960,1970,1980,1990,2000,2010,2020))+
              ylab("Subpopulation")+
              xlab("")

#problem is you've got tasu in there fro fishing but subsetted out something else uyou didn't want to ! 
ctab3<-ctab2
# unique(s3$section_name)
# s3$section_name[c(1,2,3)]
# 
# data.frame(colnames(ctab3),unique(s3$section_name))

colnames(ctab3)<-unique(s3$section_name)
ctab3<-melt(ctab3)
ctab3$presabs<-ifelse(ctab3$value>0,1,0)
names(ctab3)<- c("year","section_name","catch","presabs")

ctab3$section_name <- factor(ctab3$section_name, levels = c(
  "Cumshewa Inlet", "Englefield Bay",          
"Juan Perez Sound","Laskeek Bay",             
"Louscoone Inlet","Naden Harbour",           
"Port Louis","Rennell Sound",           
"Skidegate Inlet","Skincuttle Inlet",        
"Tasu Sound & Gowgaia Bay"
))

plot.fish<-ggplot(ctab3,aes(x=year,y=section_name,fill=factor(presabs)))+
        geom_tile(colour="black")+
        scale_fill_manual(values=c("#0C5790","#BC3C05"))+
        theme_Publication()+
        theme(legend.position="none")+
        scale_x_continuous(limits=c(1950,2015),breaks=c(1950,1960,1970,1980,1990,2000,2010,2020))+
        ylab("")+
        xlab("")

plot_grid(plot.spawn,plot.fish,labels="AUTO")


## ==================
##  Fig 2 Scaled estimated spawning biomass (X) from model and cross correlation***
## ==================

tempmatx<-scale(model$BUGSoutput$mean$X)
colnames(tempmatx) <- colnames(Y)
temp2x<-melt(tempmatx)
colnames(temp2x) <-c("time","site","x")
#temp2$site<-as.numeric(temp2$site)
temp2x$year<-years
temp2x$year<-as.numeric(temp2x$year)

#drop uncertain sites
temp2x<-drop.levels(subset(temp2x,site!="Tasu.Sound...Gowgaia.Bay"))
temp2x<-drop.levels(subset(temp2x,site!="Naden.Harbour" ))

t3x<-data.frame("x"=tapply(temp2x$x,list(temp2x$year),mean))
t3x$year<-years
t3x$year<-as.numeric(t3x$year)
t3x$site<-12
t3x$max<-tapply(temp2x$x,list(temp2x$time),max)
t3x$min<-tapply(temp2x$x,list(temp2x$time),min)
names(t3x)<-c("x","year","site","max","min")

bk<-seq(1950,2015,by=10)

# bold colors 
xscalegg<-ggplot()+
  #geom_ribbon(data=t3,aes(x=year,ymin=min,ymax=max),colour="grey60")+
  geom_hline(yintercept=0,colour="grey")+
  geom_line(data=temp2x,aes(x=year,y=x,colour=factor(site),lty=factor(site),group=site),size=1)+
  geom_line(data=t3x,aes(x=year,y=x),colour="grey",size=3,alpha=0.75)+
  # geom_vline(xintercept=1967,colour="grey",lty=2)+ #this is the closure of the reduction
  #  geom_vline(xintercept=1994,colour="grey",lty=2)+ #this is the big recent closure
  theme_Publication()+
  scale_x_continuous(breaks=bk)+
  # scale_colour_brewer(palette="Spectral")+
  # scale_colour_colorblind()+
  # scale_colour_brewer(colours=rainbow(11))+
  scale_colour_viridis_d()+
  xlab("Year")+
  ylab("Scaled Estimated Biomass (B)")


print(xscalegg)



#*** remove the two sites tasu and naden from this plot 

## ==================
##  Fig 2bc ccf of pre-fishery biomass (Z) histogram and tile plot
## ==================

#COVARIANCE OF PREDICTED STOCK BIOMASS
xmat<-model$BUGSoutput$mean$X #can look at Z here perhaps to look at pre fishery biomass correlation
colnames(xmat)<-colnames(Y)

xmat<- xmat[,-c(1,6)]

covmat_x <-matrix(NA,9,9)

for(i in 1:9){
  for(j in 1:9){
    tmp<-data.frame(xmat[,i],xmat[,j])
    tmp<-tmp[complete.cases(tmp),]
    tmp<-subset(tmp,tmp[,1]>0.011 | tmp[,2]>0.011) #subset out zeros in one or both columns 
    s1<-ccf(rank(tmp[,1]),rank(tmp[,2]),lag.max=0)
    covmat_x[i,j]<-s1$acf[,,1]
  }
}


covmat_x[lower.tri(covmat_x)] <-NA
diag(covmat_x) <- NA

colnames(covmat_x)<-colnames(Y)[-c(1,6)]
rownames(covmat_x)<-colnames(Y)[-c(1,6)]

cdf<-melt(covmat_x)
cdf2<- cdf[complete.cases(cdf),]
# cdf$X1<-substr(cdf$Var1,5,28)
# cdf$X1<-substr(cdf$Var2,5,28)

cdf$X2<-factor(cdf$Var2, levels = colnames(Y)[-c(1,6)])
cdf$X1<-factor(cdf$Var1, levels = colnames(Y)[-c(1,6)])


xccfgg<-ggplot(cdf2,aes(x=value,fill=..x..))+
  geom_histogram(colour="black")+
  geom_vline(xintercept=0,lty=2)+
  #scale_x_continuous(limits=c(-1,1))+
  scale_fill_gradient2(low="#0C5790",high="#BC3C05",limit=c(-1,1))+
  ylab("Frequency")+
  xlab("Pairwise Cross Correlaiton of Estimated Biomass (X)")+
  theme_Publication()



#All different groups
ccfmat<-ggplot(cdf,aes(x=Var1,y=Var2,fill=value,colour=value))+
  geom_tile()+
  scale_fill_gradient2(low="#0C5790",high="#BC3C05",limit=c(-1,1))+
  scale_colour_gradient2(low="#0C5790",high="#BC3C05",limit=c(-1,1))+
  theme_Publication()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="top")+
  xlab("")+
  ylab("")

multiplot(xccfgg,ccfmat,cols=2)

bottom_row <- plot_grid(xccfgg,ccfmat,labels=c('B','C'))
plot_grid(xscalegg,bottom_row,labels=c('A',''),align="h",ncol=1)


## ==================
##  Fig S3 Estimated biomass by subpopulation
## ==================

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
  

## ==================
##  Fig 3 Population Growth (U) and PDO effect  
## ==================

####Umu and Ui estimates 
Umudf<-melt(model$BUGSoutput$sims.array[,,"Umu"])
Umudf$Var3<-"Umu"
Umudf<-Umudf[,c(1,2,4,3)]
umuci<-smedian.hilow(Umudf$value,conf.int=0.95)
umuci[4]<-"Umu"
names(umuci)<-c("Median","Lower","Upper","ttt")
umudf2<-data.frame(t(data.frame(umuci)))
umudf2$Median<-as.numeric(as.character(umudf2$Median))
umudf2$Lower<-as.numeric(as.character(umudf2$Lower))
umudf2$Upper<-as.numeric(as.character(umudf2$Upper))

umuci2<-quantile(Umudf$value,c(0.25,0.75))
umudf2$Upper<-umuci2[2]
umudf2$Lower<-umuci2[1]

umuplot<-ggplot(umudf2,aes(x=ttt,y=Median,ymin=Lower,ymax=Upper))+
  geom_pointrange()+
  theme_Publication()+
  scale_y_continuous(limits=c(-0.1,0.1))+
  ylab("Intrinsic Growth Rate")

print(umuplot)


###PDO Coef Effect
pdo<-read.csv("data/pdo.csv")
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


# ggplot(pdodf2,aes(x=site,y=medianci))+
#   geom_boxplot()+
#   geom_point(aes(colour=medianci))+
#   scale_colour_gradient2(low="#BC3C05",high="#0C5790")+
#   coord_flip()+
#   theme_Publication()+
#   scale_y_continuous(limits=c(-0.10,0.2),breaks=c(-0.1,0,0.1,0.2))


pdobw<-ggplot(pdodf2,aes(x=site,y=medianci))+
  geom_boxplot()+
  geom_point(aes(colour=medianci))+
  scale_colour_gradient2(low="#BC3C05",mid="grey",high="#0C5790")+
  # coord_flip()+
  theme_Publication()+
  theme(legend.position="right")+
  scale_y_continuous(limits=c(-0.1,0.1),breaks=c(-0.1,0,0.1))+
  ylab("PDO Effect")

plot_grid(umuplot,pdobw,labels="AUTO") #*****

# multiplot(umuplot,pdobw,cols=2)


pdotsgg<-ggplot(pdodf2,aes(x=year,y=medianci))+
  geom_hline(yintercept=0)+
  geom_ribbon(aes(ymin=minci,ymax=maxci),fill="grey70",alpha=0.2)+
  geom_ribbon(aes(ymin=miniqr,ymax=maxiqr),fill="grey")+
  geom_line(lty=2)+
  geom_point(aes(colour=hc))+
  theme_Publication()+
  ylab("PDO Effect")+
  xlab("Year")+
  scale_x_continuous(breaks=bk)



print(pdotsgg)



## ==================
##  Fig 4 Proportion of fished sites and strength of fishing at multiple scales
## ==================

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
  theme_Publication()+
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
  theme_Publication()+
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
  theme_Publication()


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
  theme_Publication()


pcmerge<-merge(fdf4,fdf5, by = "year",all=TRUE)


colnames(pcmerge) <- c("year","archomed","archupper",
                       "archlower","submed","subupper","sublower")


#*** Fig 4b 
ggplot(pcmerge)+
  geom_line(aes(x=year,y=archomed))+
  geom_ribbon(aes(x=year,ymin=archlower,ymax=archupper),fill="blue",alpha=0.2)+
  geom_line(aes(x=year,y=submed),lty=2)+
  geom_ribbon(aes(x=year,ymin=sublower,ymax=subupper),fill="red",alpha=0.2)+
  theme_Publication()+
  xlab("Year")+
  ylab("Proportion Biomass Caught")+
  scale_y_continuous(limits=c(0,0.6),breaks=c(0,0.2,0.4,0.6))+
  scale_x_continuous(breaks=seq(1950,2020,by=10))
  


secfish<-ggplot(pc_tab2,aes(x=section,y=pc))+
  geom_hline(yintercept=mean(tapply(pc_tab2$pc,list(pc_tab2$section),mean)))+
  geom_boxplot()+
  geom_point()+
  coord_flip()+
  theme_Publication()+
  ylab("Fishing Effect (Pc)")

print(secfish)

tapply(pc_tab2$pc,list(pc_tab2$section),mean)


## ==================
##  Fig 5 Realized Population Growth Rate
## ==================

dmat<-model$BUGSoutput$mean$delta
colnames(dmat)<-colnames(Y)
dmat2<- exp(data.frame(dmat+median(Umudf$value)+pdodf$medianci))
dmat2$year<-c(1950:2015)
dmat2$bin<-ifelse(dmat2$year<1995,"1950-1994","1995-2015")

dmat3<-melt(dmat2,id.vars=c("bin","year"))


dmat3<-drop.levels(subset(dmat3,variable!="SHI.Tasu.Sound...Gowgaia.Bay"))
dmat3<-drop.levels(subset(dmat3,variable!="SHI.Naden.Harbour"))

#plot realized growth rate by Stocklet

ggplot(dmat3,aes(x=year,y=value))+
  geom_line(aes(colour=variable))+
  # geom_point(aes(colour=variable))+
  geom_hline(yintercept=1,colour="grey")+
  # scale_colour_brewer(palette="Spectral")+
  scale_colour_colorblind()+
  facet_wrap(~variable,scales="free")+
  ylab("Realized Population Growth exp(U+pdoeff+delta)")+
  theme(legend.position="none")+
  theme_Publication()


#plot realized growth rate for two times #*** Fig 5
ggplot(dmat3,aes(x=bin,y=value,group=variable,colour=variable,lty=variable))+
  #stat_summary(fun.data=mean_cl_normal,geom="linerange")+
  #stat_summary(fun.y=mean,geom="point")+
  stat_sum_single(median, geom="line")+
  theme_Publication()+
  # scale_colour_brewer(palette="Spectral")+
  # scale_colour_colorblind()+
  scale_colour_viridis_d()+
  # theme(legend.position="none")+
  #facet_wrap(~variable,scales="free")+
  geom_hline(yintercept=1,colour="grey")+
  xlab("Time Period")+
  ylab("Realized Growth Rate")



## ========================Portfolio================================##




## ==================
##  Is there a portfolio effect?
## ==================

#estimate the portfolio effect through a)ratio of CVs, b)mean-variance, and c) synchrony index

#plot out the time series of the 9 core sites estimated biomass 
xmat3<-exp(data.frame(model$BUGSoutput$mean$X))
colnames(xmat3)<-colnames(Y)
xmat3$year<-seq(1950,2015,by=1)
xmat3<-xmat3[,-c(1,6)]
xmatlong<-melt(xmat3,id.vars="year")

ggplot(xmatlong, aes(year, value, colour = variable))+
  geom_line()+
  # scale_colour_brewer(palette="Spectral")+
  scale_colour_colorblind()+
  theme_Publication()

xmat4<-xmat3[,-10] #remove year column

#estimate cv of average subpopulation and of 

cvarch<-sd(rowSums(xmat4))/mean(rowSums(xmat4))

cvsubpop<-tapply(xmatlong$value,list(xmatlong$variable),sd)/tapply(xmatlong$value,list(xmatlong$variable),mean)

mean(cvsubpop)


#estimate whether Taylor's power law deviates from 2, 
#if so need to use mean-variance, if not can use CV ratio

fit_taylor(xmat4) #estimates at v close t 2


plot_mv(xmat4, show = "linear", ci = TRUE)

#mean-variance
pe_mv(xmat4, ci = TRUE)
pe_avg_cv(xmat4, ci = TRUE, boot_reps = 1000)

#synchrony Loreau and deMazencourt estiamte for full eries
synchrony(xmat4)


## ==================
##  Fig 6a change in portfolio through time - process variance (delta)
## ==================

#delta through time 
tempmat<-model$BUGSoutput$mean$delta
colnames(tempmat) <- colnames(Y) #<-seq(1,11,1)
temp2<-melt(tempmat)
colnames(temp2) <-c("time","site","delta")
#temp2$site<-as.numeric(temp2$site)
temp2$year<-years
temp2$year<-as.numeric(temp2$year)
temp2<-drop.levels(subset(temp2,site!="Tasu.Sound...Gowgaia.Bay"))
temp2<-drop.levels(subset(temp2,site!="Naden.Harbour"))

t3<-data.frame("delta"=tapply(temp2$delta,list(temp2$time),mean))
t3$year<-years
t3$year<-as.numeric(t3$year)
t3$site<-12
t3$max<-tapply(temp2$delta,list(temp2$time),max)
t3$min<-tapply(temp2$delta,list(temp2$time),min)
names(t3)<-c("delta","year","site","max","min")

#***Fig 6a
dp<-ggplot()+
  # geom_ribbon(data=t3,aes(x=year,ymin=min,ymax=max),colour="grey60")+
  geom_line(data=temp2,aes(x=year,y=delta,colour=factor(site),lty=factor(site),group=site))+
  geom_line(data=t3,aes(x=year,y=delta),colour="gray30",alpha=0.75,size=2)+
  geom_hline(yintercept=0)+
  #geom_line(data=t3,aes(x=year,y=max))+
  #geom_line(data=t3,aes(x=year,y=min))+
  #scale_colour_continuous(low="green",high="red")+
  geom_vline(xintercept=1967,colour="grey",lty=2)+ #this is the closure of the reduction 
  geom_vline(xintercept=1994,colour="grey",lty=2)+ #this is the big recent closure 
  theme_Publication()+
  # scale_colour_brewer(palette="Spectral")+
  # scale_colour_colorblind()+
  scale_colour_viridis_d()+
  scale_x_continuous(breaks=bk)+
  xlab("Year")+
  ylab("Detrended Population Performance (Delta)")
print(dp)



## ==================
##  Fig 6b change in portfolio through time - different portfolio metrics 
## ==================

#way 2 & way 3 - estimate how ratio of  archipelago CV to stocklet CV changes through time
#also estimate mean-variance andersen approach 

##### loop through 10 year moving window analaysis 

xmat<-exp(data.frame(model$BUGSoutput$mean$X)) #raw biomass
colnames(xmat)<-colnames(Y)
xmat$year<-seq(1950,2015,by=1)

#make a matrix with 10 year windows 
yearstring <- seq(1950,2015,by=1)

em<- matrix(0,ncol=2,nrow=56) 
em[1,1] <- 1950
em[1,2] <- 1960

for(i in 1:55){
  
  em[i+1,1]<-yearstring[i]+1 
  em[i+1,2]<-yearstring[i]+11
  
}


em2<-em
em2<-cbind(em2,rep(0,56),rep(0,56),rep(0,56),rep(0,56))
cvvec<-rep(0,9)

for(y in 1:56){
  
  print(y)
  
  temp<- xmat[xmat$year>=em[y,1] & xmat$year<=em[y,2],]
  temp2<- temp[,-c(1,6,12)] #subtract tasu, naden, and year column
  cvvec<-rep(0,9)
  
  for(i in 1:9){
    cvvec[i] <- sd(temp2[,i])/mean(temp2[,i]) #cv for each subpop for this subset of time
  }
  
  em2[y,3] <- mean(cvvec) #average cv of subpopulations for that subset of time 
  em2[y,4] <- sd(rowSums(temp2))/mean(rowSums(temp2)) #cv of total archipelago population for that subset of time
  em2[y,5] <- pe_mv(temp2) #throwing weird error here, i think an indexing thing 
  em2[y,6] <- synchrony(temp2)
}


ratio<- em2
ratio<-data.frame(ratio)
colnames(ratio)<-c("start","finish","subCV","archCV","anders_mv","LM_Index")
ratio$ratio<-ratio$subCV/ratio$archCV
ratio$LM_Index=1-ratio$LM_Index

#CV
plot(ratio$start,ratio$ratio,type="l",xlab="Start of 10 yr moving window date",ylab="ratio of subpopCV to arch CV")

#Anderson metric
plot(ratio$start,ratio$anders_mv,type="l",xlab="Start of 10 yr moving window date",ylab="Anderson mean-variance metric")

#loreau and deMazancourt syncrhony index
plot(ratio$start,ratio$LM_Index,type="l",xlab="Start of 10 yr moving window date",ylab="Loreau and de Mazancourt Index")

ratio2<-melt(ratio[,-c(3,4)],id.vars=c("start","finish"))

ggplot(ratio2,aes(x=start,y=value))+
  geom_line()+
  facet_wrap(~variable,ncol=1,scales="free")+
  theme_Publication()


ratio3<-subset(ratio2,variable!="anders_mv")
plot_sync <-ggplot(ratio3,aes(x=start,y=value))+
  geom_line()+
  facet_wrap(~variable,ncol=1,scales="free")+
  theme_Publication()

#could add some variance to this by bootstrapping populations 

# ####### Generally archipelago CV is lower than sub population 
# ratio2<- ratio[,-c(5:6)]
# ratio3 <-melt(ratio2,id.vars=c("start","finish"))
# 
# 
# #plot out the different subpopulation and archipelago
# ggplot(ratio3,aes(x=start,y=value))+
#   geom_line(aes(colour=variable))+
#   xlab("year")+
#   ylab("CV")+
#   scale_x_continuous(breaks=bk)+
#   theme_Publication()
# 

plot_grid(dp,plot_sync,labels="AUTO",ncol=1) #***** Fig6abc



## ==================
##  Fig X change in portfolio through time - Supplement for X and delta 
## ==================

tempmat<-model$BUGSoutput$mean$X
colnames(tempmat)<-colnames(Y)
tempmat<-data.frame(exp(tempmat))
tempmat<-tempmat[,-c(1,6)]
tempmat$year<- 1950:2015
tempmat<-melt(tempmat,id.vars=c("year"))
names(tempmat)<- c("year","site","biomass")

archbiomass <- tapply(tempmat$biomass,list(tempmat$year),sum)

#way 1 - estimate the cross correlation of the subpopulation biomass through time (Xi)


xmat<-data.frame(model$BUGSoutput$mean$X)
colnames(xmat)<-colnames(Y)
xmat$year<-seq(1950,2015,by=1)

#make a matrix with 10 year windows 
yearstring <- seq(1950,2015,by=1)

em_x<- matrix(0,ncol=2,nrow=56) 
em_x[1,1] <- 1950
em_x[1,2] <- 1960

for(i in 1:55){
  
  em_x[i+1,1]<-yearstring[i]+1 
  em_x[i+1,2]<-yearstring[i]+11
  
}

em2_x <-em_x
em2_x<-cbind(em2_x,rep(0,56),rep(0,56))

#loop through each 10 year window estimate the average cross correlation 
for(y in 1:56){
  
  print(y)
  
  temp<- xmat[xmat$year>=em_x[y,1] & xmat$year<=em_x[y,2],] #pull out years of interest
  temp2<- temp[,-c(1,6,12)] #subtract naden and tasu and year column
  
  covmat<-matrix(NA,9,9)
  
  for(i in 1:9){
    for(j in 1:9){
      tmp<-data.frame(temp2[,i],temp2[,j])
      tmp<-tmp[complete.cases(tmp),]
      s1<-ccf(rank(tmp[,1]),rank(tmp[,2]),lag.max=0) 
      covmat[i,j]<-s1$acf[,,1]
    }
  }
  
  covmat[lower.tri(covmat)] <-NA
  diag(covmat) <- NA
  
  em2_x[y,3] <-   mean(covmat,na.rm=T)
  em2_x[y,4] <-   sd(covmat,na.rm=T)
  
}

colnames(em2_x)<-c("start","finish","mean","sd")

em3_x<-data.frame(em2_x)
em4_x<-melt(em3_x,id.vars=c("start","finish"))

write.csv(em4_x,"cross_correlation_moving_window_X.csv")
em4_x <- read.csv("code/cross_correlation_moving_window_X.csv")

ggplot(em4_x,aes(x=start,y=value))+
  # geom_smooth(aes(colour=variable),se=F)+
  # geom_smooth(aes(colour=variable),method="lm")+
  geom_line(aes(colour=variable))+
  # theme_acs()+
  facet_wrap(~variable,ncol=1)+
  xlab("year")+
  ylab("Cross Correlation of Estimated Biomass 10yr mving window")+
  scale_x_continuous(breaks=bk)




#######################################################
##way 5 - look at how CCF of delta changed through time 
#######################################################



dmat<-data.frame(model$BUGSoutput$mean$delta)
colnames(dmat)<-colnames(Y)
dmat$year<-seq(1950,2015,by=1)


d <- seq(1950,2015,by=1)

em_d<- matrix(0,ncol=2,nrow=56) 
em_d[1,1] <- 1950
em_d[1,2] <- 1960

for(i in 1:55){
  
  em_d[i+1,1]<-d[i]+1 
  em_d[i+1,2]<-d[i]+11
  
}


em2_d<-em_d
em2_d<-cbind(em2_d,rep(0,56),rep(0,56))

for(y in 1:56){
  
  print(y)
  
  temp<- dmat[dmat$year>=em[y,1] & dmat$year<=em_d[y,2],]
  temp2<- temp[,-c(1,6,12)] #remove tasu and naden 
  
  covmat<-matrix(NA,9,9)
  
  for(i in 1:9){
    for(j in 1:9){
      tmp<-data.frame(temp2[,i],temp2[,j])
      tmp<-tmp[complete.cases(tmp),]
      s1<-ccf(rank(tmp[,1]),rank(tmp[,2]),lag.max=0) 
      covmat[i,j]<-s1$acf[,,1]
    }
  }
  
  covmat[lower.tri(covmat)] <-NA
  diag(covmat) <- NA
  
  em2_d[y,3] <-   mean(covmat,na.rm=T)
  em2_d[y,4] <-   sd(covmat,na.rm=T)
  
}

colnames(em2_d)<-c("start","finish","mean","sd")

em3_d<-data.frame(em2_d)
em4_d<-melt(em3_d,id.vars=c("start","finish"))

write.csv(em4_d,"code/cross_correlation_moving_window_d.csv")
em4_d <-read.csv("code/cross_correlation_moving_window_d.csv")
bk<-seq(1950,2015,by=10)


ggplot(em4_d,aes(x=start,y=value))+
  # geom_smooth(aes(colour=variable),se=F)+
  # geom_smooth(aes(colour=variable),method="lm")+
  geom_line(aes(colour=variable))+
  facet_wrap(~variable,ncol=1,scales="free")+
  xlab("year")+
  ylab("mean CCF of delta ")+
  scale_x_continuous(breaks=bk)+
  theme_Publication()



#mean ccf of realized population growth rate 




dmat<-data.frame(model$BUGSoutput$mean$delta)
dmat<-data.frame(dmat+median(Umudf$value)+pdodf$medianci)

colnames(dmat)<-colnames(Y)
dmat$year<-seq(1950,2015,by=1)


d <- seq(1950,2015,by=1)

em<- matrix(0,ncol=2,nrow=56) 
em[1,1] <- 1950
em[1,2] <- 1960

for(i in 1:55){
  
  em[i+1,1]<-d[i]+1 
  em[i+1,2]<-d[i]+11
  
}


em2<-em
em2<-cbind(em2,rep(0,56),rep(0,56))

for(y in 1:56){
  
  print(y)
  
  temp<- dmat[dmat$year>=em[y,1] & dmat$year<=em[y,2],]
  temp2<- temp[,-c(1,6,12)] #remove tasu and naden 
  
  covmat<-matrix(NA,9,9)
  
  for(i in 1:9){
    for(j in 1:9){
      tmp<-data.frame(temp2[,i],temp2[,j])
      tmp<-tmp[complete.cases(tmp),]
      s1<-ccf(rank(tmp[,1]),rank(tmp[,2]),lag.max=0) 
      covmat[i,j]<-s1$acf[,,1]
    }
  }
  
  covmat[lower.tri(covmat)] <-NA
  diag(covmat) <- NA
  
  em2[y,3] <-   mean(covmat,na.rm=T)
  em2[y,4] <-   sd(covmat,na.rm=T)
  
}

colnames(em2)<-c("start","finish","mean","sd")

em3<-data.frame(em2)
em4<-melt(em3,id.vars=c("start","finish"))
em4_d2 <- em4

write.csv(em4_d2,"cross_correlation_moving_window_realized_growth.csv")
em4_d2 <- read.csv("code/cross_correlation_moving_window_realized_growth.csv")


ggplot(em4_d2,aes(x=start,y=value))+
  # geom_smooth(aes(colour=variable),se=F)+
  # geom_smooth(aes(colour=variable),method="lm")+
  geom_line(aes(colour=variable))+
  facet_wrap(~variable,ncol=1)+
  xlab("year")+
  ylab("mean CCF of realized population growth (delta+umu+pdoeff) ")+
  scale_x_continuous(breaks=bk)+
  theme_Publication()


#plotted in the paper is realized population growth ccf but i find that confusing 
#i think the SD of the deltas works betters 
