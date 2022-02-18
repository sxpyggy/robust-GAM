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

###fit models on the original data

#fit glm
glmFit <- glmReserve(GenIns, mse.method = "bootstrap", nsim=1000)

#fit gam
gamFit <- gamReserve(GenIns, fam, K=5, mse.method = "bootstrap", nsim=1000)

#fit robustgam
rgamFit <- rgamReserve(GenIns, fam, M=50, K=5, mse.method = "bootstrap", nsim=1000)

##robust bootstrap
rgamFit.s <- rgamReserve(GenIns, fam, M=50, K=5, mse.method = "sbootstrap", s=15, nsim=1000)



####add outlier

###good choice 82, 28, 44
a=2
b=8
ori=cum2incr(GenIns)
cont=ori
cont[a,b]=10*cont[a,b]
#cont[b,a]=10*cont[b,a]
#cont[5,5]=10*cont[5,5]
GenIns.out=incr2cum(cont)

###fit models on the contaminated data

#fit glm
glmFit.out <- glmReserve(GenIns.out, mse.method = "bootstrap", nsim=1000)

#fit gam
gamFit.out <- gamReserve(GenIns.out, fam, K=5, mse.method = "bootstrap", nsim=1000)

#fit robustgam
rgamFit.out <- rgamReserve(GenIns.out, fam, M=50, K=5, mse.method = "bootstrap", nsim=1000)

##robust bootstrap
rgamFit.s.out <- rgamReserve(GenIns.out, fam, M=50, K=5, mse.method = "sbootstrap", s=15, nsim=1000)


###store the results

##tables

###summary output table

library(xtable)

tablefinal <- data.frame(cbind(glmFit$summary[,c(4,5,6)],gamFit$summary[,c(4,5,6)],rgamFit$summary[,c(4,5,6)]))
colnames(tablefinal)=c(rep(c("Estimated reserves","Prediction error","Prediction error ratio"),3))

xtable(tablefinal, digits=c(0,rep(c(0,0,2),3)), rownames=T)


table=rbind(tail(glmFit.out$summary[,c(4,5,6)], 1),tail(gamFit.out$summary[,c(4,5,6)], 1),tail(rgamFit.out$summary[,c(4,5,6)], 1),tail(rgamFit.s.out$summary[,c(4,5,6)], 1))
rownames(table)=c("GLM","GAM","rGAM","rGAM.S")
colnames(table)=c("Estimated reserves","Prediction error","Prediction error ratio")

xtable(table, digits=c(0,rep(c(0,0,2),1)), rownames=T)


##graphs
res.glm=scale(residuals(glmFit$model))
res.gam=scale(gamFit$res.P)
res.rgam=scale(rgamFit$res.P)

resdata=data.frame(x=rep(1:length(res.glm), 3), Method=as.factor(rep(c(1,2,3),each=length(res.glm))), value=c(res.glm, res.gam, res.rgam))
levels(resdata$Method)=c("GLM","GAM","rGAM")


p0 <-ggplot(resdata, aes(x=x, y=value)) + 
  geom_point(aes(colour =Method), size=2)+
  scale_color_manual(values=c("grey","#F8766D", "#619CFF"))+
  facet_wrap(~ Method, scales = "free_x") +
  labs(y= "Standardised Residuals", x = "")+
  theme_bw()+
  theme(legend.position = "None", text=element_text(size=16), axis.text = element_text(size=16), strip.text.x = element_text(size = 16),panel.spacing.x=unit(0.5, "lines"))

ggsave("residual.pdf", plot = p0, width = 30, height = 10, units = "cm")


#prediction on each incremental data
glm.pred=lowertri(cum2incr(glmFit$FullTriangle))
gam.pred=lowertri(cum2incr(gamFit$FullTriangle))
rgam.pred=lowertri(cum2incr(rgamFit$FullTriangle))
glm.out.pred=lowertri(cum2incr(glmFit.out$FullTriangle))
gam.out.pred=lowertri(cum2incr(gamFit.out$FullTriangle))
rgam.out.pred=lowertri(cum2incr(rgamFit.out$FullTriangle))

plotdata=rbind(melt(glm.out.pred/glm.pred),melt(gam.out.pred/gam.pred),melt(rgam.out.pred/rgam.pred))
colnames(plotdata)=c("Acc", "Dev", "Ratio")
plotdata$Method=as.factor(rep(c(1,2,3),each=100))
levels(plotdata$Method)=c("GLM","GAM","rGAM")
plotdata$Acc=as.factor(plotdata$Acc)
plotdata$Dev=as.factor(plotdata$Dev)

p1 <-ggplot(plotdata, aes(Dev, Acc)) + 
  geom_tile(aes(fill = Ratio), colour = "white") +
  geom_point(aes(x=b, y=11-a), colour="black", size=2)+
  facet_wrap(~ Method, scales = "free_x") +
  scale_fill_gradient2(low = "white", mid="steelblue", high = "red", midpoint=1, na.value="white", breaks=seq(1, as.integer(max(plotdata$Ratio, na.rm=TRUE)), by=2))+
  scale_y_discrete(limits = rev(levels(plotdata$Acc)))+
  scale_x_discrete(position = "bottom")+
  labs(y= "Accident Period", x = "Development Period")+
  theme_bw()+
  theme(text=element_text(size=16), axis.text = element_text(size=16), strip.text.x = element_text(size = 16),panel.spacing.x=unit(0.5, "lines"))

ggsave("ratio.pdf", plot = p1, width = 30, height = 10, units = "cm")


#column factor
fac.glm= attr(ata(uptri(glmFit$model$fitted.values)), "vwtd")
fac.gam= attr(ata(uptri(gamFit$model$fitted.values)), "vwtd")
fac.rgam= attr(ata(uptri(rgamFit$model$fitted.values)), "vwtd")
fac.glm.out= attr(ata(uptri(glmFit.out$model$fitted.values)), "vwtd")
fac.gam.out= attr(ata(uptri(gamFit.out$model$fitted.values)), "vwtd")
fac.rgam.out= attr(ata(uptri(rgamFit.out$model$fitted.values)), "vwtd")

noout=data.frame(Method=as.factor(rep(1:3, each=9)),value=c(fac.glm, fac.gam, fac.rgam), Dev=as.factor(rep(2:10, 3)))
out=data.frame(Method=as.factor(rep(1:3, each=9)),value=c(fac.glm.out, fac.gam.out, fac.rgam.out), Dev=as.factor(rep(2:10, 3)))
levels(noout$Method)=c("GLM","GAM","rGAM")
levels(out$Method)=c("GLM","GAM","rGAM")

p2=ggplot(noout, aes(x=Dev, y=value, group=(Method))) +
  geom_line(aes(linetype=Method, color=Method),size=1)+
  geom_point(aes(shape=Method, color=Method),size=4)+
  scale_color_manual(values=c("grey","#F8766D", "#619CFF"))+
  ggtitle("") + 
  labs(y = "Column Effect",x= "Development Period") +
  ylim(0,6)+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=16), axis.title.x=element_text(size=16), axis.title.y=element_text(size=16),legend.title=element_text(size=16), 
        legend.text=element_text(size=16))

ggsave("columnnoout.pdf", plot = p2 , width = 20, height = 15, units = "cm")


p3=ggplot(out, aes(x=Dev, y=value, group=(Method))) +
  geom_line(aes(linetype=Method, color=Method),size=1)+
  geom_point(aes(shape=Method, color=Method),size=4)+
  scale_color_manual(values=c("grey","#F8766D", "#619CFF"))+
  ggtitle("") + 
  labs(y = "Column Effect",x= "Development Period") +
  ylim(0,6)+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=16), axis.title.x=element_text(size=16), axis.title.y=element_text(size=16),legend.title=element_text(size=16), 
        legend.text=element_text(size=16))

ggsave("columnout.pdf", plot = p3 , width = 20, height = 15, units = "cm")



#boot distribution
bootglm.out=apply(glmFit.out$sims.reserve.mean,1,sum)
bootgam.out=apply(gamFit.out$sims.reserve.mean,1,sum)
bootrgam.out=apply(rgamFit.out$sims.reserve.mean,1,sum)
bootrgam.s.out=apply(rgamFit.s.out$sims.reserve.mean,1,sum)



tmp1=data.frame(x=c(density(bootglm.out)$x,density(bootgam.out)$x,density(bootrgam.out)$x,density(bootrgam.s.out)$x),y=c(density(bootglm.out)$y,density(bootgam.out)$y,density(bootrgam.out)$y,density(bootrgam.s.out)$y),model=rep(rep(c("GLM","GAM","rGAM","rGAM.S"),each=512),1))
plotdist=tmp1
plotdist$model=factor(plotdist$model,levels=c("GLM","GAM","rGAM","rGAM.S"))
p4=ggplot(plotdist,aes(x = x, y = y, group = model, col=model)) +
geom_line(size=1) +
scale_color_manual(values=c("grey","#F8766D", "#619CFF","purple"))+
xlab("Total Claim Reserve") + ylab("Density")+ theme_bw()+   
guides(col=guide_legend(title="Method"))+
theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=16), axis.title.x=element_text(size=16), axis.title.y=element_text(size=16),legend.title=element_text(size=16), 
      legend.text=element_text(size=16))
ggsave("bootdist1.pdf", plot = p4 , width = 20, height = 15, units = "cm")


tmp2=data.frame(x=c(density(bootrgam.out)$x,density(bootrgam.s.out)$x),y=c(density(bootrgam.out)$y,density(bootrgam.s.out)$y),model=rep(rep(c("rGAM","rGAM.S"),each=512),1))
plotdist=tmp2
plotdist$model=factor(plotdist$model,levels=c("rGAM","rGAM.S"))
p5=ggplot(plotdist,aes(x = x, y = y, group = model, col=model)) +
  geom_line(size=1) +
  geom_vline(xintercept=c(quantile(bootrgam.out,0.995), quantile(bootrgam.s.out,0.995)), colour = c("#619CFF","purple"), linetype="dotted", size=1.2)+
  scale_color_manual(values=c("#619CFF","purple"))+
  xlab("Total Claim Reserve") + ylab("Density")+ theme_bw()+   
  guides(col=guide_legend(title="Method"))+
  theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=16), axis.title.x=element_text(size=16), axis.title.y=element_text(size=16),legend.title=element_text(size=16), 
        legend.text=element_text(size=16))
ggsave("bootdist2.pdf", plot = p5 , width = 20, height = 15, units = "cm")


#outstanding reserve by accident period

noout=data.frame(Method=as.factor(rep(1:3, each=9)),value=c(glmFit$summary["IBNR"][-10,], gamFit$summary["IBNR"][-10,], rgamFit$summary["IBNR"][-10,]), Dev=as.factor(rep(2:10, 3)))
out=data.frame(Method=as.factor(rep(1:3, each=9)),value=c(glmFit.out$summary["IBNR"][-10,], gamFit.out$summary["IBNR"][-10,], rgamFit.out$summary["IBNR"][-10,]), Dev=as.factor(rep(2:10, 3)))
levels(noout$Method)=c("GLM","GAM","rGAM")
levels(out$Method)=c("GLM","GAM","rGAM")

p6=ggplot(noout, aes(x=Dev, y=value, group=(Method))) +
  geom_line(aes(linetype=Method, color=Method),size=1)+
  scale_color_manual(values=c("grey","#F8766D", "#619CFF"))+
  ggtitle("") + 
  labs(y = "Outstanding Claim Reserve",x= "Accident Period") +
  ylim(0,1.7e+07)+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=16), axis.title.x=element_text(size=16), axis.title.y=element_text(size=16),legend.title=element_text(size=16), 
        legend.text=element_text(size=16))

ggsave("accnoout.pdf", plot = p6 , width = 20, height = 15, units = "cm")


p7=ggplot(out, aes(x=Dev, y=value, group=(Method))) +
  geom_line(aes(linetype=Method, color=Method),size=1)+
  scale_color_manual(values=c("grey","#F8766D", "#619CFF"))+
  ggtitle("") + 
  labs(y = "Outstanding Claim Reserve",x= "Accident Period") +
  ylim(0,1.7e+07)+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=16), axis.title.x=element_text(size=16), axis.title.y=element_text(size=16),legend.title=element_text(size=16), 
        legend.text=element_text(size=16))

ggsave("accout.pdf", plot = p7 , width = 20, height = 15, units = "cm")



#column factor with extrapolation
nex=5
expTriagle.gam=exp.gam(GenIns)
expTriagle.gam.out=exp.gam(GenIns.out)
expTriagle.rgam=exp.rgam(GenIns)
expTriagle.rgam.out=exp.rgam(GenIns.out)


fac.glm= attr(ata(uptri(glmFit$model$fitted.values)), "vwtd")
fac.gam= attr(ata(expTriagle.gam), "vwtd")
fac.rgam= attr(ata(expTriagle.rgam), "vwtd")
fac.glm.out= attr(ata(uptri(glmFit.out$model$fitted.values)), "vwtd")
fac.gam.out= attr(ata(expTriagle.gam.out), "vwtd")
fac.rgam.out= attr(ata(expTriagle.rgam.out), "vwtd")

noout=data.frame(Method=as.factor(rep(c(rep(1,9), rep(2,9+nex), rep(3,9+nex)),2)),value=c(fac.glm, fac.gam, fac.rgam), Dev=as.factor(c(rep(c(2:10, 2:15, 2:15), 2))))
out=data.frame(Method=as.factor(rep(c(rep(1,9), rep(2,9+nex), rep(3,9+nex)),2)),value=c(fac.glm.out, fac.gam.out, fac.rgam.out), Dev=as.factor(c(rep(c(2:10, 2:15, 2:15), 2))))
levels(noout$Method)=c("GLM","GAM","rGAM")
levels(out$Method)=c("GLM","GAM","rGAM")

p2=ggplot(noout, aes(x=Dev, y=value, group=(Method))) +
  geom_line(aes(linetype=Method, color=Method),size=1)+
  geom_point(aes(shape=Method, color=Method),size=4)+
  scale_color_manual(values=c("grey","#F8766D", "#619CFF"))+
  ggtitle("") + 
  labs(y = "Column Effect",x= "Development Period") +
  ylim(0,6)+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=16), axis.title.x=element_text(size=16), axis.title.y=element_text(size=16),legend.title=element_text(size=16), 
        legend.text=element_text(size=16))

ggsave("columnnooutext.pdf", plot = p2 , width = 20, height = 15, units = "cm")


p3=ggplot(out, aes(x=Dev, y=value, group=(Method))) +
  geom_line(aes(linetype=Method, color=Method),size=1)+
  geom_point(aes(shape=Method, color=Method),size=4)+
  scale_color_manual(values=c("grey","#F8766D", "#619CFF"))+
  ggtitle("") + 
  labs(y = "Column Effect",x= "Development Period") +
  ylim(0,6)+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=16), axis.title.x=element_text(size=16), axis.title.y=element_text(size=16),legend.title=element_text(size=16), 
        legend.text=element_text(size=16))

ggsave("columnoutext.pdf", plot = p3 , width = 20, height = 15, units = "cm")


