###functions for claim reserving using robust gam###

#first derivative of the Huber's loss function
huber_psi=function(x, M){
  y = x
  y = y*(abs(x)<M)+M*sign(x)*(abs(x)>=M)
  return(y)
}

#Huber's loss function
huber_loss=function(x, M){
  y = x^2/2
  y = y*(abs(x)<M)+(M*(abs(x)-0.5*M))*(abs(x)>=M)
  return(y)
}

#calculate effective degree freedom of a robustgam

edf.rgam=function(model, y, M){
  r=y-model$fitted.values
  weight=diag(huber_psi(r, M)/r)
  A=t(model$B)%*%weight%*%(model$B)
  edf=sum(diag(solve(A+model$sD)%*%A))
  return(edf)
}


#function to transform triangle into dataframe long format

tran.triagle=function(tr.incr){
  lda <- as.data.frame(tr.incr, origin = names(dimnames(tr.incr))[1], dev = names(dimnames(tr.incr))[2])
  names(lda)[1:3] <- c("origin", "dev", "value")
  #lda <- transform(lda, origin = factor(origin, levels = dimnames(tr.incr)[[1]]))
  lda <- transform(lda, origin = as.numeric(origin))
  ldaFit <- subset(lda, !is.na(lda$value))
  ldaOut <- subset(lda, is.na(lda$value))
  return(list(ldaFit=ldaFit,ldaOut=ldaOut))
}


#function to perform claim reserving using gam, input data should be cum triangle

gamReserve=function(triangle, fam, K=5, mse.method = "bootstrap", nsim = 1000){
  
  data=triangle
  lda=tran.triagle(cum2incr(data))
  ldaFit=lda$ldaFit
  ldaOut=lda$ldaOut
  modelfit <- gam(value ~ s(origin, bs="cr", k=K, m=3)+s(dev, bs="cr", k=K, m=3), data=ldaFit, family=fam) 
  yp <- predict(modelfit, newdata = ldaOut, type = "response")
  resMeanAy <- tapply(yp, factor(ldaOut$origin), sum)
  resMeanTot <- sum(resMeanAy)
  S.E=0
  
  
  if (mse.method == "bootstrap") {
    nO <- nrow(ldaFit)
    nP <- nrow(ldaOut)
    bias <- sqrt(nO/df.residual(modelfit))
    resMeanAyB <-  matrix(0, nsim, length(resMeanAy))
    sims.fitted <- matrix(0, nsim, nO)
    sims.predict <- matrix(0, nsim, nP)
    mu <- modelfit$fitted.values
    mup <- sqrt(mu)
    phi <- sum(((ldaFit$value-mu)/mup)^2)/df.residual(modelfit)
    
    ind=c(ncol(data),nO)
    
    for (i in 1:nsim) {
      ybad <- 1
      while (ybad) {
        res.P=(ldaFit$value-mu)/mup
        rn=res.P
        rn[-ind]<- sample(res.P[-ind], nO-length(ind), replace = TRUE)
        yB <- rn * mup + mu
        if (all(yB >= 0)) 
          ybad <- 0
      }
      yB=round(yB)
      modelfitB <- gam(yB ~ s(origin, bs="cr", k=K, m=3)+s(dev, bs="cr", k=K, m=3), data=ldaFit, family=fam)
      ymB <- predict(modelfitB, newdata = ldaOut, type = "response")
      resMeanAyB[i, ] <- as.numeric(tapply(ymB, factor(ldaOut$origin), sum))
      sims.fitted[i, ] <- modelfitB$fitted.values
      sims.predict[i, ] <- ymB
    }
    
    dimnames(resMeanAyB)[[2]] <- as.character(sort(unique(ldaOut$origin)))
    
    S.E.Ay=sqrt((bias*apply(resMeanAyB, 2, sd))^2+phi*resMeanAy)
    S.E.Tot=sqrt((bias*sd(rowSums(resMeanAyB)))^2+phi*resMeanTot)
    
    S.E <- c(S.E.Ay, S.E.Tot)
  }
  
  ##produce summary
  IBNR <- round(c(resMeanAy, resMeanTot))
  CV <- S.E/IBNR
  Latest <- getLatestCumulative(data)
  Latest <- Latest[-(1:(length(Latest) - length(unique(ldaOut$origin))))]
  Latest <- c(Latest, total = sum(Latest))
  Ultimate <- Latest + IBNR
  resDf <- data.frame(Latest = Latest, Dev.To.Date = Latest/Ultimate, 
                      Ultimate = Ultimate, IBNR = IBNR, S.E = S.E, CV = CV)
  row.names(resDf) <- names(Latest)
  ldaOut$value <- round(yp)
  FullTriangle.inc <- as.triangle(rbind(ldaFit, ldaOut))
  FullTriangle.cum <- incr2cum(FullTriangle.inc)
  out <- list(summary = resDf, res.P=res.P, FullTriangle.inc = FullTriangle.inc, FullTriangle = FullTriangle.cum, model = modelfit, sims.fitted = matrix(0), sims.predict = matrix(0),
              sims.reserve.mean = matrix(0))
  if (mse.method == "bootstrap") {
    out$sims.fitted<- sims.fitted
    out$sims.predict<- sims.predict
    out$sims.reserve.mean <- resMeanAyB
  }
  class(out) <- "gamReserve"
  return(out)
}


#function to perform claim reserving using robust gam, input data should be cum triangle

rgamReserve=function(triangle, fam, M=1.6, K=5, trim=0.99, mse.method = "bootstrap", s=10, nsim = 1000, sp.min=1e-6, sp.max=1e-2){
  
  data=triangle
  lda=tran.triagle(cum2incr(data))
  ldaFit=lda$ldaFit
  ldaOut=lda$ldaOut
  robustfit.gic <- robustgam.GIC(cbind(ldaFit$origin,ldaFit$dev), ldaFit$value, family=fam, p=3, K=K, c=M, show.msg=FALSE, len=10, smooth.basis='cr',sp.min=sp.min, sp.max=sp.max)
  modelfit <- robustfit.gic$optim.fit
  Xnew=data.frame(X1=ldaOut$origin,X2=ldaOut$dev)
  yp <- pred.robustgam(modelfit, Xnew, type = "response")$predict.values
  resMeanAy <- tapply(yp, factor(ldaOut$origin), sum)
  resMeanTot <- sum(resMeanAy)
  df.modelfit=length(ldaFit$value)-edf.rgam(modelfit, ldaFit$value, M)
  S.E=0
  
  
  if (mse.method == "bootstrap") {
    nO <- nrow(ldaFit)
    nP <- nrow(ldaOut)
    bias <- sqrt(nO/df.modelfit)
    resMeanAyB <-  matrix(0, nsim, length(resMeanAy))
    sims.fitted <- matrix(0, nsim, nO)
    sims.predict <- matrix(0, nsim, nP)
    mu <- modelfit$fitted.values
    mup <- sqrt(mu)
    phi <- TMSE((ldaFit$value-mu)/mup, trim)/(df.modelfit-(length(ldaFit$value)-as.integer(trim*length(ldaFit$value))))
    ind=c(ncol(data),nO)
      
    for (i in 1:nsim) {
      ybad <- 1
      while (ybad) {
        res.P=(ldaFit$value-mu)/mup
        rn=res.P
        rn[-ind]<- sample(res.P[-ind], nO-length(ind), replace = TRUE)
        yB <- rn * mup + mu
        if (all(yB >= 0)) 
          ybad <- 0
      }
      yB=round(yB)
      modelfitB <- robustgam(cbind(ldaFit$origin,ldaFit$dev), yB, family=fam, p=3, K=K, c=M, show.msg=FALSE, sp=robustfit.gic$optim.sp, smooth.basis='cr')
      ymB <- pred.robustgam(modelfitB, Xnew, type = "response")$predict.values
      resMeanAyB[i, ] <- as.numeric(tapply(ymB, factor(ldaOut$origin), sum))
      sims.fitted[i, ] <- modelfitB$fitted.values
      sims.predict[i, ] <- ymB
    }
    
    dimnames(resMeanAyB)[[2]] <- as.character(sort(unique(ldaOut$origin)))
    
    S.E.Ay=sqrt((bias*apply(resMeanAyB, 2, sd))^2+phi*resMeanAy)
    S.E.Tot=sqrt((bias*sd(rowSums(resMeanAyB)))^2+phi*resMeanTot)
    
    S.E <- c(S.E.Ay, S.E.Tot)
  }
 
  if (mse.method == "sbootstrap") {
    nO <- nrow(ldaFit)
    nP <- nrow(ldaOut)
    bias <- sqrt(nO/df.modelfit)
    resMeanAyB <-  matrix(0, nsim, length(resMeanAy))
    sims.fitted <- matrix(0, nsim, nO)
    sims.predict <- matrix(0, nsim, nP)
    mu <- modelfit$fitted.values
    mup <- sqrt(mu)
    phi <- TMSE((ldaFit$value-mu)/mup, trim)/(df.modelfit-(length(ldaFit$value)-as.integer(trim*length(ldaFit$value))))
    ind=c(ncol(data),nO)
    
    for (i in 1:nsim) {
      ybad <- 1
      while (ybad) {
        res.P=(ldaFit$value-mu)/mup
        rn=res.P
        rn[-ind]<- stratboot(res.P[-ind],s=s)
        yB <- rn * mup + mu
        if (all(yB >= 0)) 
          ybad <- 0
      }
      yB=round(yB)
      modelfitB <- robustgam(cbind(ldaFit$origin,ldaFit$dev), yB, family=fam, p=3, K=K, c=M, show.msg=FALSE, sp=robustfit.gic$optim.sp, smooth.basis='cr')
      ymB <- pred.robustgam(modelfitB, Xnew, type = "response")$predict.values
      resMeanAyB[i, ] <- as.numeric(tapply(ymB, factor(ldaOut$origin), sum))
      sims.fitted[i, ] <- modelfitB$fitted.values
      sims.predict[i, ] <- ymB
    }
    
    dimnames(resMeanAyB)[[2]] <- as.character(sort(unique(ldaOut$origin)))
    
    S.E.Ay=sqrt((bias*apply(resMeanAyB, 2, sd))^2+phi*resMeanAy)
    S.E.Tot=sqrt((bias*sd(rowSums(resMeanAyB)))^2+phi*resMeanTot)
    
    S.E <- c(S.E.Ay, S.E.Tot)
  }
  
   
  ##produce summary
  IBNR <- round(c(resMeanAy, resMeanTot))
  CV <- S.E/IBNR
  Latest <- getLatestCumulative(data)
  Latest <- Latest[-(1:(length(Latest) - length(unique(ldaOut$origin))))]
  Latest <- c(Latest, total = sum(Latest))
  Ultimate <- Latest + IBNR
  resDf <- data.frame(Latest = Latest, Dev.To.Date = Latest/Ultimate, 
                      Ultimate = Ultimate, IBNR = IBNR, S.E = S.E, CV = CV)
  row.names(resDf) <- names(Latest)
  ldaOut$value <- round(yp)
  FullTriangle.inc <- as.triangle(rbind(ldaFit, ldaOut))
  FullTriangle.cum <- incr2cum(FullTriangle.inc)
  out <- list(summary = resDf, res.P=res.P, FullTriangle.inc = FullTriangle.inc, FullTriangle = FullTriangle.cum, model = modelfit, sims.fitted = matrix(0), sims.predict = matrix(0),
              sims.reserve.mean = matrix(0))
  
    out$sims.fitted<- sims.fitted
    out$sims.predict<- sims.predict
    out$sims.reserve.mean <- resMeanAyB
  
  class(out) <- "rgamReserve"
  return(out)
}


###function to perform stratified bootstrap
stratboot=function(res, s=11){

order.res<-order(abs(res))
n=length(res)
S=list()
ns=ceiling(n/s)
for(i in 1:s){
  S[[i]]=order.res[(ns*(i-1)+1):(ns*i)]
}
S=lapply(S, function(x) x[!is.na(x)])

strat.sample=NULL
for(j in 1:s){
  strat.sample=c(strat.sample,sample(S[[j]],length(S[[j]]),replace=TRUE))
}
strat.sample=sample(strat.sample)
return(res[strat.sample])
}


###cutoff triangle, the input is incremental claim

cutoff=function(tri,cut=7){
  test=tri[1:cut,1:cut]
  train=test
  for(i in 1:nrow(train)){
    for (j in 1:ncol(train)){
    if((i+j)<(cut+2)){
      train[i,j]=test[i,j]
    } else {
      train[i,j]=NA
    }
    }
  }
  return(list(train=as.triangle(train), test=as.triangle(test)))
}



###lower triangle

lowertri=function(tri){
  trinew=matrix(NA, nrow(tri), ncol(tri))
  for(i in 1:nrow(tri)){
    for (j in 1:ncol(tri)){
      if((i+j)>nrow(tri)+1){
        trinew[i,j]=tri[i,j]
      } 
    }
  }
  return(trinew)
}


uptri=function(data, nrow=10, ncol=10){
  origin=NULL
  dev=NULL
  for(i in 1:nrow){
  origin=c(origin, 1:(ncol-i+1))
  }
  for(j in 1:ncol){
    dev=c(dev, rep(j,ncol-j+1))
  } 
  datanew=data.frame(value=data, origin, dev)
  trinew=as.triangle(datanew, "origin", "dev", "value")
  return(trinew)
}

TMSE=function(error, trim=0.95){
errorsq=sort(error^2)[1:as.integer(trim*(length(error^2)))]
sum(errorsq)
}





exp.rgam=function(triangle,nex=5, M=50){
lda=tran.triagle(cum2incr(triangle))
ldaFit=lda$ldaFit
ldaOut=lda$ldaOut

###fit robust gam
robustfit.gic <- robustgam.GIC(cbind(ldaFit$origin,ldaFit$dev), ldaFit$value, family=fam, p=3, K=5, c=M, show.msg=FALSE, len=10, smooth.basis='cr')
robustfit <- robustfit.gic$optim.fit

#plot(nonrobustfit$fitted.values[ldaFit$origin==1])

###expolation length
nex=5
extra=data.frame(origin=rep(unique(lda$ldaFit$origin),nex),dev=rep((max(lda$ldaFit$origin)+1):(max(lda$ldaFit$origin)+nex),each=length(unique(lda$ldaFit$origin))),value=NA)
ldaOut=rbind(ldaOut,extra)

Xnew=data.frame(X1=ldaOut$origin,X2=ldaOut$dev)
yp <- pred.robustgam(robustfit, Xnew, type = "response")$predict.values

ldaOut$value <- round(yp)
ldaFit$value=robustfit$fitted.values
ldawhole=rbind(ldaFit,ldaOut)
expTriagle=as.triangle(ldawhole, origin="origin", dev="dev", "value")
return(expTriagle)
}




exp.gam=function(triangle,nex=5){
  K=5
  lda=tran.triagle(cum2incr(triangle))
  ldaFit=lda$ldaFit
  ldaOut=lda$ldaOut
  
  ###fit robust gam
  modelfit <- gam(value ~ s(origin, bs="cr", k=K, m=3)+s(dev, bs="cr", k=K, m=3), data=ldaFit, family=fam) 
  
  
  #plot(nonrobustfit$fitted.values[ldaFit$origin==1])
  
  ###expolation length
  nex=5
  extra=data.frame(origin=rep(unique(lda$ldaFit$origin),nex),dev=rep((max(lda$ldaFit$origin)+1):(max(lda$ldaFit$origin)+nex),each=length(unique(lda$ldaFit$origin))),value=NA)
  ldaOut=rbind(ldaOut,extra)
  
  Xnew=data.frame(origin=ldaOut$origin,dev=ldaOut$dev)
  yp <- predict(modelfit, newdata = Xnew, type = "response")
  
  ldaOut$value <- round(yp)
  ldaFit$value=modelfit$fitted.values
  ldawhole=rbind(ldaFit,ldaOut)
  expTriagle=as.triangle(ldawhole, origin="origin", dev="dev", "value")
  return(expTriagle)
}


