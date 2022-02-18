library(ChainLadder)
library(robustgam)
library(robustgam)
library(reshape2)
library(ggplot2)

setwd("")
source('functions.R')

##############
fam=poisson()
set.seed(1)

reserveALL=NULL

for(s in 1:1000){
cum<-read.csv(paste("./cum_cf/cum_CF_",s,".csv",sep=""))[,-1]
if(any(cum2incr(cum)<1)) next
ncol=ncol(cum)
colnames(cum)=1:ncol
data=cutoff(cum2incr(cum),ncol)
train=data$train
test=data$test

###fit models on the original data

#fit glm
glmFit <- glmReserve(incr2cum(train), mse.method = "bootstrap", nsim=2)

#fit gam
gamFit <- gamReserve(incr2cum(train), fam, K=2, mse.method = "bootstrap", nsim=2)

#fit robustgam
rgamFit <- rgamReserve(incr2cum(train), fam, M=10, K=2, mse.method = "bootstrap", nsim=2)


####randomly add one outlier in the data

long=as.data.frame(train)
long=long[complete.cases(long), ]
ind=c(ncol,nrow(long))
indout=sample(c(1:nrow(long))[-ind],1)

a=as.numeric(long[indout,][1])
b=as.numeric(long[indout,][2])
train.cont=train
train.cont[a,b]=10*train.cont[a,b]
train.out=incr2cum(train.cont)

###fit models on the contaminated data

#fit glm
glmFit.out <- glmReserve(train.out, mse.method = "bootstrap", nsim=2)

#fit gam
gamFit.out <- gamReserve(train.out, fam, K=2, mse.method = "bootstrap", nsim=2)

#fit robustgam
rgamFit.out <- rgamReserve(train.out, fam, M=10, K=2, mse.method = "bootstrap", nsim=2)

###store the results

#prediction on each incremental data
glm.pred=lowertri(cum2incr(glmFit$FullTriangle))
gam.pred=lowertri(cum2incr(gamFit$FullTriangle))
rgam.pred=lowertri(cum2incr(rgamFit$FullTriangle))
glm.out.pred=lowertri(cum2incr(glmFit.out$FullTriangle))
gam.out.pred=lowertri(cum2incr(gamFit.out$FullTriangle))
rgam.out.pred=lowertri(cum2incr(rgamFit.out$FullTriangle))


###summary output table
re.glm=tail(glmFit$summary[,4],1)
re.gam=tail(gamFit$summary[,4],1)
re.rgam=tail(rgamFit$summary[,4],1)
re.glm.out=tail(glmFit.out$summary[,4],1)
re.gam.out=tail(gamFit.out$summary[,4],1)
re.rgam.out=tail(rgamFit.out$summary[,4],1)
re.true=sum(lowertri(test),na.rm=TRUE)

reserve=c(re.glm,re.gam,re.rgam,re.glm.out,re.gam.out,re.rgam.out,re.true)
reserveALL=rbind(reserveALL, reserve)
}

est=reserveALL[,1:6]
true=reserveALL[,7]
result=apply(est, 2, function(x){sqrt(sum((x-true)^2)/length(true))})
write.csv(reserveALL, file="reserveALL.csv")


####produce table

reserveALL=read.csv("reserveALL.csv")[,2:8]
est=reserveALL[,1:7]
true=reserveALL[,7]
mean=apply(est, 2, mean)
mspe=apply(est, 2, function(x){sqrt(sum((x-true)^2)/length(true))})
mape=apply(est, 2, function(x){sum(abs(x-true))/length(true)})
percent=apply(est, 2, quantile, 0.995)

table=cbind(mean, percent, mape, mspe)
rownames(table)=c("GLM","GAM","rGAM", "GLM","GAM","rGAM", "True")
colnames(table)=c("Mean","99.5 Quantile","MAPE","RMSPE")

library(xtable)
xtable(table, digits=0, rownames=T)