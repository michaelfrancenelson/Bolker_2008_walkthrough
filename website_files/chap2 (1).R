###################################################
### chunk number 1: 
###################################################
source("chapskel.R")


###################################################
### chunk number 2: 
###################################################
require(emdbook)
data(SeedPred)
data(SeedPred_wide)
set.seed(1001)
library(lattice)
## set trellis to white fg, gray fill
trellis.par.set(color=FALSE)
## palette(gray(seq(0,0.9,len=10)))


###################################################
### chunk number 3: 
###################################################
long = subset(SeedPred,as.numeric(station)<8 & dist==25 & date<as.Date("1999-04-04"))[,
                                                c("station","date","dist","species","seeds")]
rownames(long) = 1:nrow(long)
head(long)


###################################################
### chunk number 4: 
###################################################
wide = head(subset(SeedPred_wide,as.numeric(station)<8 & dist==25)[,c(1,3,2,4,6)])
rownames(wide) = 1:nrow(wide)
wide


###################################################
### chunk number 5: 
###################################################
reshape(wide,direction="long",timevar="date",varying=4:5)
reshape(long,direction="wide",timevar="date",idvar=c("station","dist","species"))


###################################################
### chunk number 6: 
###################################################
library(reshape)
recast(wide,formula=...~.,id.var=c("station","dist","species"))
recast(long,formula=station+dist+species~...,id.var=c("station","dist","species","date"))


###################################################
### chunk number 7:  eval=FALSE
###################################################
## data = read.table("mydata.dat",header=TRUE)


###################################################
### chunk number 8: 
###################################################
1:5
c("yes","no","maybe")


###################################################
### chunk number 9: 
###################################################
SeedPred[[2]]
SeedPred[["species"]]
SeedPred$species


###################################################
### chunk number 10: 
###################################################
SeedPred[,2]
SeedPred[,"species"]


###################################################
### chunk number 11: 
###################################################
SeedPred[1:10,]


###################################################
### chunk number 12: 
###################################################
summary(SeedPred[,1:4])


###################################################
### chunk number 13: 
###################################################
str(SeedPred)


###################################################
### chunk number 14: 
###################################################
class(SeedPred)
sapply(SeedPred,class)


###################################################
### chunk number 15: 
###################################################
oldwid = options(width=75)


###################################################
### chunk number 16: 
###################################################
head(SeedPred)


###################################################
### chunk number 17: 
###################################################
options(width=65)


###################################################
### chunk number 18:  eval=FALSE
###################################################
## table(SeedPred$station,SeedPred$species)


###################################################
### chunk number 19: 
###################################################
t1 = table(SeedPred$station,SeedPred$species)
head(t1)


###################################################
### chunk number 20:  eval=FALSE
###################################################
## SeedPred[24,"species"] = "mmu"


###################################################
### chunk number 21: 
###################################################
attach(SeedPred,warn=FALSE)
SeedPred_10 = subset(SeedPred,dist==10)
SeedPred_25 = subset(SeedPred,dist==25)
s10_means = tapply(SeedPred_10$seeds,list(SeedPred_10$date,SeedPred_10$species),mean,na.rm=TRUE)
s10_dates = unique(SeedPred_10$date)
s25_means = tapply(SeedPred_25$seeds,list(SeedPred_25$date,SeedPred_25$species),mean,na.rm=TRUE)
s25_dates = unique(SeedPred_25$date)


###################################################
### chunk number 22: 
###################################################
op = par(cex=1.5,las=1,bty="l",lwd=2)
matplot(s10_dates,s10_means,type="b",axes=FALSE,col="black",pch=1:8,lty=1,log="y",
        xlab="",ylab="Mean seeds remaining",ylim=c(0.05,15),cex=0.5,
        xlim=c(10665,10950))
mtext("Date",side=1,line=2.25,cex=1.5)
matlines(s25_dates,s25_means,col="gray",pch=1:8,lty=1,type="b",cex=0.5)
axis(side=2,at=c(0.05,0.5,5))
axis.Date(side=1,s10_dates)
box()
par(xpd=NA)
legend(10700,15,c("10 m","25 m"),col=c("black","gray"),pch=16,ncol=2)
#t(outer(levels(SeedPred$species),c("10 m","25 m"),paste,sep=": ")),pch=rep(1:8,each=2),
#       col=rep(c("black","gray"),8),ncol=4,cex=0.5)
y.adj = c(abz=0,cd=0,cor=-0.2,dio=-0.25,mmu=0,pol=0,psd=-0.03,uva=0.2)
text(rep(10930,length(levels(species))),
     s10_means[nrow(s10_means),]+y.adj,levels(species),cex=0.5,adj=0)
points(rep(10954,length(levels(species))),
       s10_means[nrow(s10_means),]+y.adj,pch=1:8,cex=0.5,lwd=1)
par(op)


###################################################
### chunk number 23:  eval=FALSE
###################################################
## bwplot(s10_means ~ species|cut(as.numeric(date),5),data=SeedPred)


###################################################
### chunk number 24: 
###################################################
op = par(mfrow=c(1,2))
par(cex=1.5,las=1,bty="l",lwd=2)
par(mar=c(5,4,2,0)+0.1,mgp=c(2.5,1,0))
plot(jitter(SeedPred$available),jitter(SeedPred$taken),xlab="Seeds available",ylab="Seeds taken",lwd=1)
##
corner.label2("a","topleft",inset=0.025)
library(plotrix)
par(mar=c(5,2,2,2)+0.1)
with(subset(SeedPred,available>0),
     sizeplot(available,taken,scale=0.6,
              xlim=c(0.3,6),ylim=c(-0.5,5.5),axes=FALSE,
              xlab="Seeds available",ylab="",col="gray",pow=0.4))
t1 = with(subset(SeedPred,available>0),table(available,taken))
axis(side=1,at=1:5)
axis(side=2,at=0:5,labels=rep("",6))
box()
par(xpd=NA); corner.label2("b","topleft",inset=0.025); par(xpd=FALSE)
text(row(t1)[t1>0],col(t1)[t1>0]-1,t1[t1>0],cex=0.5)


###################################################
### chunk number 25: 
###################################################
op = par(cex=1.5,las=1,bty="l",lwd=2,mgp=c(2.5,1,0))
cols = gray.colors(6,start=0.6)
barplot(t(log10(t1+1)),beside=TRUE,legend=FALSE,
        xlab="",
        ylab="log10(1+# observations)",axes=FALSE,ylim=c(0,4),
        col=cols)
mtext(side=1,"Number of seeds available",line=2.25,cex=1.5)
axis(side=2,at=0:3)
par(xpd=NA)
legend(0,4.25,0:5,fill=cols,
       ncol=6,x.intersp=0.1,bty="n",cex=0.75)
text(15,4.3,"Number taken",xpd=TRUE,adj=1,cex=0.75)
par(op)


###################################################
### chunk number 26:  eval=FALSE
###################################################
## library(scatterplot3d)
## avail = row(t1)[t1>0]
## y = col(t1)[t1>0]-1
## z = log10(t1[t1>0])
## scatterplot3d(-avail,-y,z,type="h",
##               angle=50,pch=16)
## ## or with RGL
## library(rgl)
## open3d(FOV=1)
## load("seedgridpos.RData")
## par3d(seedgridpos)
## plot3d(avail,y,z,lit=TRUE,
##        col.pt="gray",xlab="",ylab="",zlab="",type="s",
##        axes=FALSE,
##        size=0.5,
##        zlim=c(0,4))
## plot3d(avail,y,z,add=TRUE,type="h",
##        axes=FALSE,box=FALSE,size=4,col=gray(0.2))
## ## don't draw bbox ...
## axes3d(edges=c("x+-","y--","z++"))
## ## axes3d()
## grid3d(c("x+","y-","z"))
## text3d(8.75,5,2,"Frequency",adj=0.5,size=2)
## par3d(ignoreExtent=TRUE)
## text3d(3,6.5,-0.6,"Available",adj=0.5,size=2)
## text3d(-1,3.5,-0.6,"Taken",adj=0.5,size=2)
## ## r1 = rgl.pop("bboxdeco")
## rgl.bbox(alpha=1,col=NA,front="cull",back="cull",
##          xat=1,xlab="",yat=-1,ylab="",zat=1,zlab="") ## HACK
## rgl.material(color="black",alpha=1)
## rgl.postscript(file="seed3d.eps")
## rgl.postscript(file="seed3d.pdf",fmt="pdf")
## rgl.close()


###################################################
### chunk number 27: 
###################################################
frac.taken = SeedPred$taken/SeedPred$available


###################################################
### chunk number 28: 
###################################################
mean.frac.by.avail = tapply(frac.taken,available,mean,na.rm=TRUE)


###################################################
### chunk number 29: 
###################################################
n.by.avail = table(available)
sd.by.avail = tapply(frac.taken,available,sd,na.rm=TRUE)
se.by.avail = sd.by.avail/sqrt(n.by.avail)


###################################################
### chunk number 30: 
###################################################
library(gplots,warn.conflicts=FALSE)
op = par(cex=1.5,las=1,bty="l",lwd=2,mgp=c(3,1,0))
b = barplot2(mean.frac.by.avail[-1],plot.ci=TRUE,
  ci.l=mean.frac.by.avail[-1]-se.by.avail[-1],
  ci.u=mean.frac.by.avail[-1]+se.by.avail[-1],xlab="",
  ylab="Fraction taken")
mtext("Number of seeds available",side=1,line=2,cex=1.5)
par(op)


###################################################
### chunk number 31: 
###################################################
mean.frac.by.avail.sp = tapply(frac.taken,list(available,species),mean,na.rm=TRUE)
mean.frac.by.avail.sp = na.omit(mean.frac.by.avail.sp)
barplot(mean.frac.by.avail.sp,beside=TRUE)


###################################################
### chunk number 32: 
###################################################
nz = subset(SeedPred,taken>0)


###################################################
### chunk number 33: 
###################################################
barchart(table(nz$available,nz$species,nz$taken),stack=FALSE)


###################################################
### chunk number 34: 
###################################################
## trellis.par.set(bar.fill=new.bar.fill)
trellis.par.set(canonical.theme(color=FALSE))
nz = subset(SeedPred,taken>0)
print(barchart(table(nz$available,nz$species,nz$taken),stack=FALSE,
               xlab="Frequency"))
##op = par(mfrow=c(3,3),mar=c(2,2,1,1))
##frac.taken.by.species = split(frac.taken,Species)
##invisible(sapply(frac.taken.by.species,hist,xlab="",ylab="",main="",
##       col="gray")) ## ,ylim=c(0,15)))
##par(op)


###################################################
### chunk number 35: 
###################################################
detach(SeedPred)


###################################################
### chunk number 36: 
###################################################
data(ReedfrogPred)
gcols = gray(c(0.9,0.7,0.4))
op = par(cex=1.5,las=1,bty="l",lwd=2)
b = boxplot(propsurv~size*density*pred,data=ReedfrogPred,axes=FALSE,
  col=rep(rep(gcols,each=2),2),ylab="Proportion surviving")
axis(side=2)
axis(side=1,labels=FALSE,at=1:12)
staxlab(side=1,at=1:12,labels=gsub("\\.[^.]*$","",b$names),nlines=2)
text(3.5,0.5,"no pred",cex=1.5)
text(9.5,1,"pred",cex=1.5,xpd=NA)
legend(1,0.2,fill=gcols,c(10,25,35),ncol=3,cex=0.7,bty="n")
text(1,0.23,"density",adj=0)
box()
par(op)                                                                  


###################################################
### chunk number 37:  eval=FALSE
###################################################
## bwplot(propsurv~density|pred*size,data=ReedfrogPred,horizontal=FALSE)


###################################################
### chunk number 38: 
###################################################
library(splines)
data(ReedfrogFuncresp)
data(ReedfrogSizepred)
ltys = c(1,3,2,4)
op = par(mfrow=c(1,2))
par(cex=1.5,las=1,bty="l",lwd=2,mgp=c(2.5,1,0),mar=c(5,3.4,4,2)+0.1)
attach(ReedfrogFuncresp,warn=FALSE)
plot(Initial,Killed,xlim=c(0,120),  ylab="Number killed",
     xlab="Initial density",axes=FALSE)
axis(side=2)
axis(side=1,at=seq(0,120,by=40))
box()
lines(lowess(Initial,Killed),lty=ltys[1])
lines(predict(smooth.spline(Initial,Killed,df=5),x=0:100),lty=ltys[2])
##lm1 = lm(Killed ~ ns(Initial, df = 5), data = ReedfrogSizepred)
##lm2 = lm(Killed ~ ns(Initial), data = ReedfrogSizepred)
##lmq = lm(Killed ~ Initial+I(Initial^2), data = ReedfrogSizepred)
##lines(predict(lm1,newdata=data.frame(Initial=1:100)),lty=2)
## lines(predict(lmq,newdata=data.frame(Initial=1:100)),lty=4)
meanvals = tapply(Killed,Initial,mean)
lines(unique(Initial),meanvals,lty=ltys[3])
## abline(lm(Killed ~ Initial, data = ReedfrogSizepred),lty=3)
corner.label2("a","topleft",inset=0.025)
par(xpd=NA)
legend(55,15,
       c("lowess","spline","means"),
       lty=ltys,
       bty="n",xjust=0,yjust=1,
       y.intersp=0.8)
par(xpd=FALSE)
detach(ReedfrogFuncresp)
attach(ReedfrogSizepred,warn=FALSE)
## TEST ME: FONT
sizeplot(TBL,Kill,xlab="Tadpole Size (TBL in mm)",ylab="Number killed",
         xlim=c(0,40),scale=1)
corner.label2("b","topleft",inset=0.025)
##lines(lowess(TBL,Kill))
## lines(lowess(TBL,Kill,f=0.6))
##lm1 = lm(Kill ~ ns(TBL, df = 3), data = ReedfrogSizepred)
## lines(predict(lm1,newdata=data.frame(TBL=0:40)))
meanvals = tapply(Kill,TBL,mean)
lines(unique(TBL),meanvals,lty=ltys[3])
detach(ReedfrogSizepred)


###################################################
### chunk number 39: 
###################################################
data(DamselRecruitment)
data(DamselRecruitment_sum)
attach(DamselRecruitment,warn=FALSE)
init.dens = init/area*1000
surv.dens = surv/area*1000
detach(DamselRecruitment)
op = par(lwd=2,cex=1.5,las=1,bty="l",mgp=c(2.5,1,0)) 
par(lwd=1)
plot(init.dens,surv.dens,cex=0.5,
  xlab=expression(paste("Initial density",({}/0.1*m^2))),
         ylab="Recruit density (6 months)",log="x",axes=FALSE)
par(lwd=2)
axis(side=2,at=c(0,5,10))
axis(side=1,at=c(0.5,5,50,500))
box()
plotCI(DamselRecruitment_sum$settler.den,DamselRecruitment_sum$surv.den,
 DamselRecruitment_sum$SE,
   add=TRUE,pch=16,col="darkgray",gap=0)
lines(lowess(init.dens,surv.dens))
par(xpd=NA)
legend("topleft",c("actual","target","lowess"),
       cex=1,
       pt.cex=c(0.5,1,1),
       pch=c(1,16,NA),lty=c(NA,1,1),
       lwd=c(1,2,2),bty="n",
       merge=TRUE,col=c("black","darkgray","black"))
par(op)


###################################################
### chunk number 40: 
###################################################
data(DamselSettlement)
attach(DamselSettlement,warn=FALSE)
op = par(lwd=2,cex=1.5,las=1,bty="l",mar=c(4,4,2,1)+0.1) 
hist(density[density<200],main="",breaks=c(0,seq(1,201,by=4)),col="gray",
        xlab="",
        ylab="Probability density")
box()
lines(density(density[density<200],from=0))
mtext(side=1,expression(paste("Settler density",({}/0.1*m^2))),
        at=100,line=2.5,cex=1.5) 
detach(DamselSettlement)
par(op)


###################################################
### chunk number 41:  eval=FALSE
###################################################
## bwplot(log10(1+density)~pulse|site,data=DamselSettlement,
##               horizontal=FALSE)


###################################################
### chunk number 42:  eval=FALSE
###################################################
## library(reshape)
## x2 = melt(DamselSettlement,measure.var="density")
## x3 = cast(x2,pulse+obs~...)


###################################################
### chunk number 43: 
###################################################
library(reshape)
x2 = melt(DamselSettlement,measure.var="density")
x3 = cast(x2,pulse+obs~...)
print(splom(log10(1+x3[,3:5]),groups=x3$pulse,pch=as.character(1:6),col=1,cex=2,
            xlab=""))
##  pairs(log10(1+x3[,3:5]),cex.labels=1)


###################################################
### chunk number 44: 
###################################################
gobydat = read.csv("GobySurvival.csv",
                colClasses=
                c(rep("factor",4),
                  rep("numeric",4)))


###################################################
### chunk number 45: 
###################################################
attach(gobydat)


###################################################
### chunk number 46: 
###################################################
meansurv = (d1+d2)/2


###################################################
### chunk number 47: 
###################################################
dens.cat = ifelse(density>median(density),"high","low")
dens.cat = factor(dens.cat,levels=c("low","high"))
qual.cat = ifelse(qual>median(qual),"high","low")
qual.cat = factor(qual.cat,levels=c("low","high"))


###################################################
### chunk number 48: 
###################################################
print(xyplot(jitter(meansurv,factor=2)~jitter(density,2)|qual.cat,
             xlab=list(label="Density",cex=1.5),
             ylab=list(label="Mean survival time",cex=1.5),
             scales=list(cex=1.5),
             panel=function(x,y,...) {
               panel.xyplot(x,y,...)
               panel.lmline(x,y)
             }))


###################################################
### chunk number 49: 
###################################################
survtab = table(meansurv); survtab


###################################################
### chunk number 50: 
###################################################
csurvtab = cumsum(rev(survtab));csurvtab


###################################################
### chunk number 51: 
###################################################
csurvtab = rev(csurvtab)


###################################################
### chunk number 52:  eval=FALSE
###################################################
## survtab/csurvtab


###################################################
### chunk number 53: 
###################################################
round(survtab/csurvtab,2)


###################################################
### chunk number 54: 
###################################################
intcat = interaction(qual.cat,dens.cat)
cattab = table(intcat)
survtab = table(meansurv,intcat)
## reverse order
survtab = survtab[nrow(survtab):1,]
## get cumulative sum: this is number surviving until day x
csurvtab = apply(survtab,2,cumsum)
## divide by total number in category
cnsurvtab = sweep(csurvtab,2,cattab,"/")
days = as.numeric(rownames(csurvtab))
fracmort = survtab/csurvtab


###################################################
### chunk number 55: 
###################################################
op = par(mfrow=c(1,2))
par(mar=c(5,4,2,2.5)+0.1)
par(lwd=2,cex=1.5,las=1,bty="l",mgp=c(2.5,1,0))
matplot(days[days<20],fracmort[days<20,],xlab="",
        ylab="Proportional mortality",type="l",col=1)
mtext(side=1,"Time (days)",cex=1.5,line=2)
corner.label2("a","topleft",inset=0.025)
matplot(days,cnsurvtab,type="s",xlab="",
        ylab="Fraction surviving",xlim=c(0,40),
        log="y",col=1)
mtext(side=1,"Time (days)",cex=1.5,line=2)
corner.label2("b","topleft",inset=0.025,bg="white")
par(xpd=NA)
legend(c(5,par("usr")[2]),c(0.3,1.05),
       c("high qual/high density",
         "low qual/low density",
         "high qual/low density",
         "low qual/high density"),
       lty=c(4,1,2,3),col=1,cex=0.75,bty="n")
par(xpd=FALSE)
par(op)


###################################################
### chunk number 56: 
###################################################
detach(gobydat)


###################################################
### chunk number 57: 
###################################################
dat_10 = read.csv("duncan_10m.csv",na.strings="?",comment="")
dat_25 = read.csv("duncan_25m.csv",na.strings="?",comment="")


###################################################
### chunk number 58:  eval=FALSE
###################################################
## dat_10 = dat_10[1:159,]
## dat_25 = dat_25[,1:39]


###################################################
### chunk number 59: 
###################################################
library(reshape)
dat_10_melt = melt(dat_10,id.var=1:2)


###################################################
### chunk number 60: 
###################################################
date_10 = paste(dat_10_melt[,3],"1999",sep=".")


###################################################
### chunk number 61: 
###################################################
dat_10_melt[,3] = as.Date(date_10,format="X%d.%b.%Y")


###################################################
### chunk number 62: 
###################################################
names(dat_10_melt) = c("station","species","date","seeds")


###################################################
### chunk number 63: 
###################################################
dat_25_melt = melt(dat_25,id.var=1:2)
date_25 = paste(dat_25_melt[,3],"1999",sep=".")
dat_25_melt[,3] = as.Date(date_25,format="X%d.%b.%Y")
names(dat_25_melt) = c("station","species","date","seeds")


###################################################
### chunk number 64: 
###################################################
split_10 = split(dat_10_melt,dat_10_melt$station)


###################################################
### chunk number 65: 
###################################################
for (i in 1:length(split_10)) {
  x = split_10[[i]]
  tcum = as.numeric(x$date-x$date[1])
  tint = as.numeric(c(NA,diff(x$date)))
  taken = c(NA,-diff(x$seeds))
  available = c(NA,x$seeds[-nrow(x)])
  split_10[[i]] = data.frame(x,tcum,tint,taken,available)
}


###################################################
### chunk number 66: 
###################################################
dat_10 = do.call("rbind",split_10)


###################################################
### chunk number 67: 
###################################################
split_25 = split(dat_25_melt,dat_25_melt$station)
for (i in 1:length(split_25)) {
  x = split_25[[i]]
  tcum = as.numeric(x$date-x$date[1])
  tint = as.numeric(c(NA,diff(x$date)))
  taken = c(NA,-diff(x$seeds))
  available = c(NA,x$seeds[1:(nrow(x)-1)])
  split_25[[i]] = data.frame(x,tcum,tint,taken,available)
}
dat_25 = do.call("rbind",split_25)


###################################################
### chunk number 68: 
###################################################
dat_10 = data.frame(dat_10,dist=rep(10,nrow(dat_10)))
dat_25 = data.frame(dat_25,dist=rep(25,nrow(dat_25)))
SeedPred = rbind(dat_10,dat_25)


###################################################
### chunk number 69: 
###################################################
SeedPred$station = factor(SeedPred$station)
SeedPred$dist = factor(SeedPred$dist)


###################################################
### chunk number 70: 
###################################################
SeedPred = SeedPred[,c("station","dist","species","date","seeds",
                         "tcum","tint","taken","available")]
SeedPred_wide = reshape(SeedPred[order(SeedPred$date),],
                         direction="wide",timevar="date",idvar=c("station","dist","species"),
                         drop=c("tcum","tint","taken","available"))


###################################################
### chunk number 71: 
###################################################
rm("available","taken","tcum","tint","t1")  ## clean up


###################################################
### chunk number 72:  eval=FALSE
###################################################
## ## don't need to do this every time
## save("SeedPred",file="SeedPred.rda")
## save("SeedPred_wide",file="SeedPred_wide.rda")


###################################################
### chunk number 73: 
###################################################
attach(SeedPred)


###################################################
### chunk number 74: 
###################################################
SeedPred_10 = subset(SeedPred,dist==10)
SeedPred_25 = subset(SeedPred,dist==25)


###################################################
### chunk number 75: 
###################################################
s10_means = tapply(SeedPred_10$seeds,list(SeedPred_10$date,SeedPred_10$species),mean,na.rm=TRUE)
s25_means = tapply(SeedPred_25$seeds,list(SeedPred_25$date,SeedPred_25$species),mean,na.rm=TRUE)


###################################################
### chunk number 76: 
###################################################
matplot(s10_means,log="y",type="b",col=1,pch=1:8,lty=1)
matlines(s25_means,type="b",col="gray",pch=1:8,lty=1)


###################################################
### chunk number 77: 
###################################################
plot(jitter(SeedPred$available),jitter(SeedPred$taken))


###################################################
### chunk number 78: 
###################################################
library(plotrix)
sizeplot(SeedPred$available,SeedPred$taken,scale=0.5,pow=0.5,xlim=c(-2,6),ylim=c(-2,5))
t1 = table(SeedPred$available,SeedPred$taken)
text(row(t1)-1,col(t1)-1,t1)


###################################################
### chunk number 79:  eval=FALSE
###################################################
## library(gplots)
## balloonplot(t1)


###################################################
### chunk number 80: 
###################################################
plot(t1)


###################################################
### chunk number 81: 
###################################################
mosaicplot(~available+taken,data=SeedPred)


###################################################
### chunk number 82: 
###################################################
barplot(t(log10(t1+1)),beside=TRUE,xlab="Available",ylab="log10(1+# observations)")


###################################################
### chunk number 83: 
###################################################
barplot(t(t1+1),log="y",beside=TRUE,xlab="Available",ylab="1+# observations")


###################################################
### chunk number 84: 
###################################################
mean.frac.by.avail = tapply(frac.taken,available,mean,na.rm=TRUE)
n.by.avail = table(available)
se.by.avail = tapply(frac.taken,available,sd,na.rm=TRUE)/sqrt(n.by.avail)
barplot2(mean.frac.by.avail,plot.ci=TRUE,
  ci.l=mean.frac.by.avail-se.by.avail,
  ci.u=mean.frac.by.avail+se.by.avail,xlab="Number available",ylab="Fraction taken")


###################################################
### chunk number 85: 
###################################################
library(plotrix)
frac.taken = SeedPred$taken/SeedPred$available
mean.frac.by.avail.by.species = tapply(frac.taken,list(available,species),mean,na.rm=TRUE)
n.by.avail.by.species = table(available,species)
se.by.avail.by.species = tapply(frac.taken,list(available,species),sd,na.rm=TRUE)/sqrt(n.by.avail.by.species)
b = barplot(mean.frac.by.avail.by.species,beside=TRUE)
plotCI(b,mean.frac.by.avail.by.species,se.by.avail.by.species,add=TRUE,pch=".",gap=FALSE)


###################################################
### chunk number 86: 
###################################################
avail = row(t1)[t1>0]
taken = col(t1)[t1>0]-1
freq = log10(t1[t1>0])


###################################################
### chunk number 87: 
###################################################
library(scatterplot3d)
scatterplot3d(-avail,-taken,freq,type="h",
              angle=50,pch=16)


###################################################
### chunk number 88: 
###################################################
library(rgl)
plot3d(avail,taken,freq,lit=TRUE,
       col.pt="gray",type="s",
       size=0.5,
       zlim=c(0,4))


###################################################
### chunk number 89: 
###################################################
plot3d(avail,taken,freq,add=TRUE,type="h",size=4,col=gray(0.2))
grid3d(c("x+","y-","z"))


###################################################
### chunk number 90: 
###################################################
rgl.close()


###################################################
### chunk number 91: 
###################################################
histogram(~frac.taken|species,xlab="Fraction taken")


###################################################
### chunk number 92: 
###################################################
op = par(mfrow=c(3,3))
for (i in 1:length(levels(species))) {
  hist(frac.taken[species==levels(species)[i]],xlab="Fraction taken",main="",
       col="gray")
}
par(op)


###################################################
### chunk number 93: 
###################################################
detach(SeedPred)


###################################################
### chunk number 94: 
###################################################
data(ReedfrogPred)
data(ReedfrogFuncresp)
data(ReedfrogSizepred)


###################################################
### chunk number 95: 
###################################################
graycols = rep(rep(gray(c(0.4,0.7,0.9)),each=2),2)
boxplot(propsurv~size*density*pred,data=ReedfrogPred,col=graycols)


###################################################
### chunk number 96: 
###################################################
attach(ReedfrogFuncresp,warn=FALSE)


###################################################
### chunk number 97: 
###################################################
plot(Initial,Killed,xlim=c(0,100),ylab="Number killed",xlab="Initial density")


###################################################
### chunk number 98: 
###################################################
lines(lowess(Initial,Killed))


###################################################
### chunk number 99: 
###################################################
meanvals = tapply(Killed,Initial,mean)
densvals = unique(Initial)
lines(densvals,meanvals,lty=3)


###################################################
### chunk number 100: 
###################################################
lms = smooth.spline(Initial,Killed,df = 5)


###################################################
### chunk number 101: 
###################################################
ps = predict(lms,x=0:100)
lines(ps,lty=2)


###################################################
### chunk number 102: 
###################################################
library(splines)
lm1 = lm(Killed ~ ns(Initial, df = 5), data = ReedfrogSizepred)
p1 = predict(lm1,newdata=data.frame(Initial=1:100))
lines(p1,lty=2)


###################################################
### chunk number 103: 
###################################################
lm2 = lm(Killed ~ Initial, data = ReedfrogSizepred)
lmq = lm(Killed ~ Initial+I(Initial^2), data = ReedfrogSizepred)


###################################################
### chunk number 104: 
###################################################
detach(ReedfrogFuncresp)


###################################################
### chunk number 105: 
###################################################
data(DamselRecruitment)
data(DamselRecruitment_sum)
attach(DamselRecruitment)
attach(DamselRecruitment_sum)


###################################################
### chunk number 106: 
###################################################
plot(init.dens,surv.dens,log="x")
plotCI(settler.den,surv.den,SE,
   add=TRUE,pch=16,col="darkgray",gap=0)
lines(lowess(init.dens,surv.dens))


###################################################
### chunk number 107: 
###################################################
detach(DamselRecruitment)
detach(DamselRecruitment_sum)


###################################################
### chunk number 108: 
###################################################
attach(DamselSettlement)
hist(density[density<200],breaks=c(0,seq(1,201,by=4)),col="gray",
        xlab="",
        ylab="Prob. density")
lines(density(density[density<200],from=0))


###################################################
### chunk number 109: 
###################################################
hist(log(1+density))
hist(density[density>0],breaks=50)


###################################################
### chunk number 110: 
###################################################
h1 = hist(density,breaks=c(0,seq(1,201,by=4),500),plot=FALSE)
b= barplot(h1$counts,space=0)
axis(side=1,at=b,labels=h1$mids)


###################################################
### chunk number 111: 
###################################################
bwplot(log10(1+density)~pulse|site,data=DamselSettlement,
              horizontal=FALSE)


###################################################
### chunk number 112: 
###################################################
densityplot(~density,groups=site,data=DamselSettlement,xlim=c(0,100))
bwplot(density~site,horizontal=FALSE,data=DamselSettlement)
bwplot(density~site|pulse,horizontal=FALSE,data=DamselSettlement)
bwplot(log10(1+density)~site|pulse,data=DamselSettlement,
             panel=panel.violin,
             horizontal=FALSE)
boxplot(density~site*pulse)


###################################################
### chunk number 113: 
###################################################
library(reshape)
x2 = melt(DamselSettlement,measure.var="density")
x3 = cast(x2,pulse+obs~...)


###################################################
### chunk number 114: 
###################################################
pairs(log10(1+x3[,3:5]))


###################################################
### chunk number 115: 
###################################################
splom(log10(1+x3[,3:5]),groups=x3$pulse,pch=as.character(1:6),col=1)


###################################################
### chunk number 116: 
###################################################
detach(DamselSettlement)


###################################################
### chunk number 117: 
###################################################
attach(gobydat)
xyplot(jitter(meansurv,factor=2)~jitter(density,2)|qual.cat,
             xlab="Density",ylab="Mean survival time")


###################################################
### chunk number 118: 
###################################################
panel1 = function(x,y) {
  panel.xyplot(x,y)
  panel.lmline(x,y)
}


###################################################
### chunk number 119: 
###################################################
xyplot(jitter(meansurv,factor=2)~jitter(density,2)|qual.cat,
             xlab="Density",ylab="Mean survival time",
             panel=panel1)
detach(gobydat)


###################################################
### chunk number 120: 
###################################################
intcat = interaction(qual.cat,dens.cat)
cattab = table(intcat)


###################################################
### chunk number 121: 
###################################################
survtab = table(meansurv,intcat)


###################################################
### chunk number 122: 
###################################################
survtab = survtab[nrow(survtab):1,]
csurvtab = apply(survtab,2,cumsum)


###################################################
### chunk number 123: 
###################################################
cnsurvtab = sweep(csurvtab,2,cattab,"/")


###################################################
### chunk number 124: 
###################################################
fracmort = survtab/csurvtab


###################################################
### chunk number 125: 
###################################################
days = as.numeric(rownames(csurvtab))


###################################################
### chunk number 126: 
###################################################
matplot(days,cnsurvtab,type="s",xlab="Time (days)",
        ylab="Proportion of cohort surviving",
        log="y")


