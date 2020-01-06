###################################################
### chunk number 1: 
###################################################
library(bbmle)
library(emdbook)
library(plotrix)
library(lattice)
library(Hmisc)
library(MASS)
source("chapskel.R")


###################################################
### chunk number 2: 
###################################################
sizefun = function(S,phi=20,eps=1,beta=1) {
   exp(eps*(phi-S))/(1+exp(beta*eps*(phi-S)))
}
maxxval = function(phi=20,eps=1,beta=1) {
  phi+log(beta-1)/(beta*eps)
}
maxyval = function(beta) {
  (beta-1)^(-1/beta)/(1+1/(beta-1))
}


###################################################
### chunk number 3: 
###################################################
op = par(cex=1.2,mar=c(4,4,0,1)+0.1,
  mgp=c(2.5,0.75,0),las=1,bty="l")
curve(sizefun(x),from=5,to=40,xlab="Prey size",ylab="Predation risk")
curve(sizefun(x,beta=1.1),add=TRUE,lty=2)
curve(sizefun(x,beta=3),add=TRUE,lty=3)
curve(sizefun(x,beta=1.1,eps=5),add=TRUE,lty=4)
text(x=c(19.5,5,6.5,28),y=c(0.9,0.5,0.08,0.5),
     c(expression(list(beta==1,epsilon==1)),expression(list(beta==1.1,epsilon==1)),
       expression(list(beta==1.1,epsilon==5)),
       expression(list(beta==3,epsilon==1))),
     adj=0)
arrows(27.5,0.49,21.5,0.4)
par(op)


###################################################
### chunk number 4: 
###################################################
data(ReedfrogSizepred)
attach(ReedfrogSizepred)


###################################################
### chunk number 5: 
###################################################
modlogist = function(x,eps,beta,phi) {
 exp(eps*(phi-x))/(1+exp(beta*eps*(phi-x)))
}
##ricker = function(x,a,b) {
##  a*x*exp(-b*x)
##}
powricker = function(x,a,b,alpha) {
  b*(x/a*exp(1-x/a))^alpha
}
tricker = function(x,a,b,t,min=0.0001) {
  ifelse(x<t,min,b*((x-t)/a*exp(1-(x-t)/a)))
}


###################################################
### chunk number 6: 
###################################################
NLL.modlogist = function(eps,beta,phi) {
  p.pred = modlogist(TBL,eps,beta,phi)
  -sum(dbinom(Kill,size=10,prob=p.pred,log=TRUE))
}
NLL.modlogist.bb = function(eps,beta,phi,theta) {
  p.pred = modlogist(TBL,eps,beta,phi)
  -sum(dbetabinom(Kill,size=10,prob=p.pred,theta=theta,log=TRUE))
}
##NLL.ricker = function(a,b) {
##  p.pred = ricker(TBL,a,b)
##  -sum(dbinom(Kill,size=10,prob=p.pred,log=TRUE))
##}
NLL.powricker = function(a,b,alpha) {
  p.pred = powricker(TBL,a,b,alpha)
  -sum(dbinom(Kill,size=10,prob=p.pred,log=TRUE))
}
NLL.tricker = function(a,b,t) {
  p.pred = tricker(TBL,a,b,t)
  -sum(dbinom(Kill,size=10,prob=p.pred,log=TRUE))
}


###################################################
### chunk number 7: 
###################################################
FSP.modlogist = mle2(NLL.modlogist,
                    start=list(eps=5,beta=1.1,phi=15))
FSP.modlogist.bb = mle2(NLL.modlogist.bb,
                       start=as.list(c(coef(FSP.modlogist),list(theta=1000))),
                       control=list(parscale=c(1,1,1,1000)))


###################################################
### chunk number 8:  eval=FALSE
###################################################
## confint(FSP.modlogist.bb,which="theta")


###################################################
### chunk number 9: 
###################################################
FSP.modlogist2 = mle2(NLL.modlogist,
                     start=list(eps=0.357,beta=8.99,phi=9.75),
                     control=list(maxit=1000))


###################################################
### chunk number 10: 
###################################################
rbind(coef(FSP.modlogist),coef(FSP.modlogist2))


###################################################
### chunk number 11: 
###################################################
##FSP.ricker = mle2(NLL.ricker,start=list(a=0.4,b=0.3))
FSP.powricker = mle2(NLL.powricker,start=list(a=0.4,b=0.3,alpha=1))
FSP.tricker = mle2(NLL.tricker,start=list(a=0.4,b=0.3,t=8))
##


###################################################
### chunk number 12: 
###################################################
alpha.ci = confint(FSP.powricker,parm="alpha",quietly=TRUE)


###################################################
### chunk number 13: 
###################################################
a1 = AICtab(FSP.powricker,FSP.modlogist,FSP.tricker,FSP.modlogist2,
       weights=TRUE,sort=TRUE)
a1$AIC = round(a1$AIC,1)
a1$weight = round(a1$weight,3)
attr(a1,"row.names") = c("modified logistic (fit 2)",
                          "truncated Ricker",
                          "modified logistic (fit 1)",
                          "generalized Ricker")
latex(a1,file="",title="",table.env=FALSE)


###################################################
### chunk number 14: 
###################################################
op = par(lwd=2,bty="l",las=1,mgp=c(3,1,0),cex=1.5)
sizeplot(TBL,Kill,xlab="Tadpole size\n(total body length in mm)",
         ylab="",
         xlim=c(0,40),scale=1,ylim=c(0,9))
mtext(side=2,"Number killed",line=2,cex=1.5,las=0)
with(as.list(coef(FSP.modlogist)),curve(modlogist(x,eps,beta,phi)*10,
                                    add=TRUE))
## with(as.list(coef(FSP.ricker)),curve(ricker(x,a,b)*10,add=TRUE,lty=2))
with(as.list(coef(FSP.powricker)),curve(powricker(x,a,b,alpha)*10,add=TRUE,lty=2))
with(as.list(coef(FSP.tricker)),curve(tricker(x,a,b,t)*10,add=TRUE,lty=3))
with(as.list(coef(FSP.modlogist2)),curve(modlogist(x,eps,beta,phi)*10,
                                    add=TRUE,lty=4))
legend("topright",c("modified logistic",
                    "generalized Ricker",
                    "truncated Ricker",
                    "modified logistic #2"),cex=0.7,
       lty=1:4,bty="n")


###################################################
### chunk number 15: 
###################################################
rogers.pred = function(N0,a,h,P,T) {
   N0 - lambertW(a*h*N0*exp(-a*(P*T-h*N0)))/(a*h)
}
holling2.pred = function(N0,a,h,P,T) {
  a*N0*P*T/(1+a*h*N0)
}


###################################################
### chunk number 16:  eval=FALSE
###################################################
## curve(rogers.pred(x,a=1,h=0.2,P=1,T=1),from=0,to=60,
##   ylab="Number eaten/unit time",xlab="Initial number",ylim=c(0,5),
##   main="Predation: a=1, h=0.2")
## curve(rogers.pred(x,a=1,h=0.2,P=1,T=5)/5,add=TRUE,lty=2,from=0)
## curve(rogers.pred(x,a=1,h=0.2,P=1,T=0.2)*5,add=TRUE,lty=3,from=0)
## curve(rogers.pred(x,a=1,h=0.2,P=1,T=10)/10,add=TRUE,lty=4,from=0)
## curve(holling2.pred(x,a=1,h=0.2),add=TRUE,lty=1,lwd=2,from=0)
## abline(h=5)
## legend(30,2,
##    c(paste("Rogers, T=",c(0.2,1,5,10),sep=""),
##     "Holling type II"),lwd=c(rep(1,4),2),lty=c(3,1,2,4,1))


###################################################
### chunk number 17: 
###################################################
data(ReedfrogFuncresp)
attach(ReedfrogFuncresp)


###################################################
### chunk number 18: 
###################################################
NLL.rogers = function(a,h,T,P) {
  if (a<0 || h<0) return(NA)
  prop.exp = rogers.pred(Initial,a,h,P,T)/Initial
  -sum(dbinom(Killed,prob=prop.exp,size=Initial,log=TRUE))
}
NLL.holling2 = function(a,h,P=1,T=1) {
  -sum(dbinom(Killed,prob=a*T*P/(1+a*h*Initial),size=Initial,log=TRUE))
}


###################################################
### chunk number 19: 
###################################################
FFR.rogers = mle2(NLL.rogers,start=list(a=0.012,h=0.84),
  data=list(T=14,P=3))
FFR.holling2 = mle2(NLL.holling2,start=list(a=0.012,h=0.84),
  data=list(T=14,P=3))


###################################################
### chunk number 20: 
###################################################
a2 = AICtab(FFR.rogers,FFR.holling2,sort=TRUE,weights=TRUE)
a2$AIC = round(a2$AIC,1)
a2$weight = round(a2$weight,3)
attr(a2,"row.names") = c("Holling type II",
                          "Rogers")
latex(a2,file="",title="",table.env=FALSE)


###################################################
### chunk number 21: 
###################################################
esttab = signif(rbind(coef(FFR.rogers),coef(FFR.holling2)),3)
dimnames(esttab) = list(c("Rogers","Holling type II"),
          c("$a$","$h$"))
latex(esttab,file="",title="",table.env=FALSE)


###################################################
### chunk number 22: 
###################################################
op = par(lwd=2,bty="l",las=1,mgp=c(3,1,0),cex=1.5)
plot(Initial,Killed,xlim=c(0,100))
with(as.list(coef(FFR.rogers)),
     curve(rogers.pred(x,a,h,T=14,P=3),add=TRUE))
with(as.list(coef(FFR.rogers)),
     curve(a*14*3*x/(1+a*h*x),add=TRUE,lty=3))
with(as.list(coef(FFR.holling2)),
     curve(a*14*3*x/(1+a*h*x),add=TRUE,lty=2))
legend("topleft",
       c("Rogers","Rogers (no depletion)",
         "Holling"),
       lty=c(1,3,2),bty="n")
par(op)


###################################################
### chunk number 23: 
###################################################
xpars.Funcresp = list(T=14,P=3,vol=220,area=1.2*0.8,size=12.8)
xpars.Sizepred = list(T=3,P=2,vol=25,area=pi*0.16^2,initprey=10)


###################################################
### chunk number 24: 
###################################################
n.Funcresp = nrow(ReedfrogFuncresp)
n.Sizepred = nrow(ReedfrogSizepred)
combInit = c(ReedfrogFuncresp$Initial,rep(xpars.Sizepred$initprey,n.Sizepred))
combSize = c(rep(xpars.Funcresp$size,n.Funcresp),ReedfrogSizepred$TBL)
combKilled = c(ReedfrogFuncresp$Killed,ReedfrogSizepred$Kill)
combP = rep(c(xpars.Funcresp$P/xpars.Funcresp$area,
               xpars.Sizepred$P/xpars.Sizepred$area),
             c(n.Funcresp,n.Sizepred))
combT = rep(c(xpars.Funcresp$T,xpars.Sizepred$T),
             c(n.Funcresp,n.Sizepred))


###################################################
### chunk number 25: 
###################################################
prop.eaten = function(N0,S,h,P,T,eps,beta,phi,minprop=.Machine$double.eps) {
  a = modlogist(S,eps=eps,beta=beta,phi=phi)
  N.eaten = rogers.pred(N0,a=a,h=h,P=P,T=T)
  prop = N.eaten/N0
  prop[prop<=0] = minprop
  prop[prop>=1] = 1-minprop
  prop
}


###################################################
### chunk number 26: 
###################################################
NLL.rogerscomb = function(a,h,eps,beta,phi,T=combT,P=combP) {
  if (h<0) return(NA)
  prob = prop.eaten(combInit,combSize,h,P,T,eps,beta,phi)
  dprob = dbinom(combKilled,prob=prob,size=combInit,log=TRUE)
  -sum(dprob)
}


###################################################
### chunk number 27: 
###################################################
startvals = c(list(h=coef(FFR.rogers)["h"]),as.list(coef(FSP.modlogist)))


###################################################
### chunk number 28: 
###################################################
FPcomb = mle2(NLL.rogerscomb,start=startvals,
              method="Nelder-Mead")
confint(FPcomb,method="quad")


###################################################
### chunk number 29:  eval=FALSE
###################################################
## FPci2 = try(confint(FPcomb),silent=TRUE)
## startvals2 = list(h=0.88254539,eps=-0.01068388,
##   beta=0.93631205,phi=383.67804375)
## FPcomb2 = mle2(NLL.rogerscomb,start=startvals2,
##   method="Nelder-Mead",control=list(parscale=c(1,1,1,500)))


###################################################
### chunk number 30:  eval=FALSE
###################################################
## cat(h,eps,beta,phi,"\n")


###################################################
### chunk number 31:  eval=FALSE
###################################################
## if (any(!is.finite(prob))) cat("NAs:",h,eps,beta,phi,"\n")


###################################################
### chunk number 32:  eval=FALSE
###################################################
## if (any(!is.finite(dprob))) {
##   browser()
## }


###################################################
### chunk number 33: 
###################################################
FPcomb = mle2(NLL.rogerscomb,start=startvals,
              method="L-BFGS-B",
              lower=c(0.7,0.5,1,14),
              upper=c(1.8,2.25,2,20),
              control=list(parscale=c(1,1,1,10)))
FPcomb.ci = confint(FPcomb)


###################################################
### chunk number 34: 
###################################################
NLL.rogerscomb.bb = function(a,h,eps,beta,phi,
  theta,T=combT,P=combP) {
  if (h<0) return(NA)
  prob = prop.eaten(combInit,combSize,h,P,T,eps,beta,phi)
  dprob = dbetabinom(combKilled,prob=prob,size=combInit,theta=theta,log=TRUE)
  -sum(dprob)
}
FPcomb.bb = mle2(NLL.rogerscomb.bb,start=c(startvals,list(theta=100)),
              method="L-BFGS-B",
              lower=c(0.7,0.5,1,14,0.1),
              upper=c(1.8,2.25,2,20,Inf),
              control=list(parscale=c(1,1,1,10,100)))
theta.ci = confint(FPcomb.bb,parm="theta")


###################################################
### chunk number 35:  eval=FALSE
###################################################
## c1 = coef(FSP.modlogist)
## FSP.expprop.mean = modlogist(12.8,c1["eps"],c1["beta"],c1["phi"])
## FSP.exppropvar = deltavar(exp(eps*(phi-12.8))/(1+exp(beta*eps*(phi-12.8))),
##          meanval=coef(FSP.modlogist),Sigma=vcov(FSP.modlogist))
## FSP.expprop.deltaci = FSP.expprop.mean+c(-1.96,1.96)*sqrt(FSP.exppropvar)


###################################################
### chunk number 36: 
###################################################
c1 = coef(FSP.modlogist)
FSP.expprop.mean = modlogist(12.8,c1["eps"],c1["beta"],c1["phi"])


###################################################
### chunk number 37: 
###################################################
c2 = coef(FPcomb)
FP.expprop.mean = prop.eaten(N0=10,S=12.8,c2["h"],
  P=2/0.08,T=3,eps=c2["eps"],beta=c2["beta"],phi=c2["phi"])


###################################################
### chunk number 38: 
###################################################
set.seed(1001)
FSP.expprop.pars = mvrnorm(5000,mu=c1,Sigma=vcov(FSP.modlogist))
FSP.expprop.val = numeric(5000)
for (i in 1:5000) {
  FSP.expprop.val[i] = modlogist(12.8,FSP.expprop.pars[i,1],
                                  FSP.expprop.pars[i,2],FSP.expprop.pars[i,3])
}
FSP.expprop.ppi = quantile(FSP.expprop.val,c(0.025,0.975))


###################################################
### chunk number 39:  eval=FALSE
###################################################
## FP.exppropvar = deltavar(prop.eaten(N0=10,S=12.8,h,P=2/0.08,T=3,eps,beta,phi),
##   meanval=coef(FPcomb),Sigma=vcov(FPcomb))
## FP.expprop.deltaci = FP.expprop.mean+c(-1.96,1.96)*sqrt(FP.exppropvar)


###################################################
### chunk number 40: 
###################################################
FP.expprop.pars = mvrnorm(5000,mu=c2,Sigma=vcov(FPcomb))
FP.expprop.val = numeric(5000)
for (i in 1:5000) {
  FP.expprop.val[i] = prop.eaten(N0=10,S=12.8,P=2/0.08,T=3,
                                  h=FP.expprop.pars[i,"h"],
                                  eps=FP.expprop.pars[i,"eps"],
                                  beta=FP.expprop.pars[i,"beta"],
                                  phi=FP.expprop.pars[i,"phi"])
}
FP.expprop.ppi = quantile(FP.expprop.val,c(0.025,0.975))


###################################################
### chunk number 41: 
###################################################
##esttab = rbind(c(FSP.expprop.mean,NA,NA),
##               c(NA,FSP.expprop.deltaci),
##                c(NA,FSP.expprop.ppi),
##                c(FP.expprop.mean,NA,NA),
##                c(NA,FP.expprop.deltaci),
##                c(NA,FP.expprop.ppi))
## dimnames(esttab) = list(c("size-pred",
##                            "delta method",
##                            "PPI",
##                            "combined",
##                            "delta method",
##                            "PPI"),
##                          c("mean","low","high"))
esttab = rbind(c(FSP.expprop.mean,FSP.expprop.ppi),
                c(FP.expprop.mean,FP.expprop.ppi))
dimnames(esttab) = list(c("size-pred","combined"),
                         c("mean","low","high"))
latex(round(esttab,3),file="",title="",table.env=FALSE)


###################################################
### chunk number 42: 
###################################################
spredarea = with(xpars.Sizepred,P*T/area)
frpredarea = with(xpars.Funcresp,P*T/area)
spredvol = with(xpars.Sizepred,P*T/vol)
frpredvol = with(xpars.Funcresp,P*T/vol)


###################################################
### chunk number 43: 
###################################################
svec = seq(0,40,by=0.25)
FP.expprop = prop.eaten(N0=10,S=svec,c2["h"],
  P=2/0.08,T=3,eps=c2["eps"],beta=c2["beta"],phi=c2["phi"])
FSP.expprop = with(as.list(coef(FSP.modlogist)),
  modlogist(svec,eps,beta,phi))
FP.expprop.pars = mvrnorm(5000,mu=c2,Sigma=vcov(FPcomb))
FP.expprop.vals = 
  t(apply(FP.expprop.pars,1,
        function(x) {
          with(as.list(x),prop.eaten(N0=10,S=svec,P=2/0.08,T=3,
                                     h=h,
                                  eps=eps,
                                  beta=beta,
                                  phi=phi))}))
FP.expprop.pars = mvrnorm(5000,mu=c2,Sigma=vcov(FPcomb))
FSP.expprop.vals = 
  t(apply(FSP.expprop.pars,1,
        function(x) {
          with(as.list(x),modlogist(svec,eps,beta,phi))}))
FSP.env = t(apply(FSP.expprop.vals,2,quantile,c(0.025,0.975)))
FP.env = t(apply(FP.expprop.vals,2,quantile,c(0.025,0.975)))


###################################################
### chunk number 44: 
###################################################
op = par(lwd=2,bty="l",las=1,mgp=c(3,1,0),cex=1.5)
svec = seq(0,40,by=0.25)
with(ReedfrogSizepred,sizeplot(TBL,Kill/10,ylim=c(0,1),
                               xlab="Total body length",
                               ylab="Proportion eaten"))
lines(svec,FP.expprop)
lines(svec,FSP.expprop,lty=2)
abline(v=12.8,col="gray")
## matlines(svec,FSP.env,lty=2,col="gray")
## matlines(svec,FP.env,lty=1,col="gray")

## c3 = coef(FPcomb2)
## FP.expprop2 = prop.eaten(N0=10,S=svec,c3["h"],
##   P=2/0.08,T=3,eps=c3["eps"],beta=c3["beta"],phi=c3["phi"])
## lines(svec,FP.expprop2,lty=3)
## bit = 0.3
## plotCI(12.8-bit,FP.expprop.mean,li=FP.expprop.ppi[1],ui=FP.expprop.ppi[2],
##        add=TRUE,pch=NA)
## plotCI(12.8+bit,FSP.expprop.mean,li=FSP.expprop.ppi[1],ui=FSP.expprop.ppi[2],
##        add=TRUE,pch=NA,slty=2)
legend("topright",c("size-pred.","combined"),lty=2:1,bty="n")
par(op)


###################################################
### chunk number 45: 
###################################################
detach(ReedfrogSizepred)
detach(ReedfrogFuncresp)


###################################################
### chunk number 46: 
###################################################
library(emdbookx)
data(GobySurvival)
attach(GobySurvival)


###################################################
### chunk number 47: 
###################################################
day1 = d1-1
day2 = ifelse(d2==70,Inf,d2-1)


###################################################
### chunk number 48:  eval=FALSE
###################################################
## day1 = pmax(day1,.Machine$double.eps)


###################################################
### chunk number 49: 
###################################################
xmax = 1
op = par(cex=1.5,lwd=2,
  mgp=c(2.5,0.75,0),las=1,bty="l")
par(mfrow=c(1,2))
shapes = c(1,0.5,0.2)
curve(pweibull(x,shape=shapes[1],scale=0.5,lower.tail=FALSE),from=0.001,to=xmax,log="y",
      ylab="Fraction surviving",xlab="Time",ylim=c(0.1,1),axes=FALSE)
axis(side=1,at=c(0,0.5,1))
axis(side=2)
box()
curve(pweibull(x,shape=shapes[2],scale=0.5,lower.tail=FALSE),add=TRUE,lty=2)
curve(pweibull(x,shape=shapes[3],scale=0.5,lower.tail=FALSE),add=TRUE,lty=3)
## abline(h=exp(-1),col=2)
## abline(v=0.5,col=2)
text(rep(xmax+0.1,3),pweibull(xmax,shape=shapes,scale=0.5,lower.tail=FALSE),
     paste("shape=",shapes,sep=""),
     adj=0,xpd=NA)
text(corner.loc(x=1,y=1),adj=1,"scale=0.5",cex=1.5)
par(mar=c(5,3,4,4)+0.1)
scales = c(0.5,1,2)
curve(pweibull(x,shape=0.5,scale=scales[1],lower.tail=FALSE),from=0.001,to=xmax,
      log="y",
      axes=FALSE,ylim=c(0.1,1),xlab="Time",ylab="")
box()
axis(side=1,at=c(0,0.5,1))
curve(pweibull(x,shape=0.5,scale=scales[2],lower.tail=FALSE),add=TRUE,lty=2)
curve(pweibull(x,shape=0.5,scale=scales[3],lower.tail=FALSE),add=TRUE,lty=3)
text(rep(xmax+0.1,3),pweibull(xmax,shape=0.5,scale=scales,lower.tail=FALSE),
     paste("scale=",scales,sep=""),
     adj=0,xpd=NA)
## abline(h=exp(-1))
text(corner.loc(x=1,y=1),adj=1,"shape=0.5",cex=1.5)
par(op)


###################################################
### chunk number 50: 
###################################################
op = par(lwd=2,bty="l",las=1,mgp=c(3,1,0),cex=1.5)
dens.cat = factor(ifelse(density>median(density),"high","low"),
                   levels=c("low","high"))
qual.cat = factor(ifelse(qual>median(qual),"high","low"),
                   levels=c("low","high"))
intcat = interaction(qual.cat,dens.cat)
cattab = table(intcat)
## when (average) did they die?
meansurv = (d1+d2)/2
## tabulate
morttab = table(meansurv,intcat)
## reverse order
morttab = morttab[nrow(morttab):1,]
## get cumulative sum: this is number surviving until day x
csurvtab = apply(morttab,2,cumsum)
## divide by total number in category
cnsurvtab = sweep(csurvtab,2,cattab,"/")
days = as.numeric(rownames(csurvtab))
matplot(days,cnsurvtab,type="s",xlab="Time (days)",
        ylab="Proportion of cohort surviving",
        log="y",col=1,xlim=c(0,40))
legend("topright",
       c("high quality/high density",
         "low quality/low density",
         "high quality/low density",
         "low quality/high density"),bty="n",
       lty=c(4,1,2,3),cex=0.9)
par(op)


###################################################
### chunk number 51:  eval=FALSE
###################################################
## numsurv = diag(table(d1,d2))
## days = c(1,4,8,11,70)
## csurv = nrow(X)-c(0,cumsum(numsurv))
## 1-exp(log(csurv[2:5]/csurv[1:4])/diff(days))


###################################################
### chunk number 52: 
###################################################
NLL.GS.xqdi = function(lscale0,lscale.q,lscale.d,lscale.i,
                          lscale.x2,lscale.x3,lscale.x4,lscale.x5,
                          lshape) {
  lscalediff = c(0,lscale.x2,lscale.x3,lscale.x4,lscale.x5)
  scale=exp(lscale0+lscalediff[exper]+
    lscale.q*qual+(lscale.d+lscale.i*qual)*density)
  shape = exp(lshape)
  -sum(log(pweibull(day2,shape,scale)-
          pweibull(day1,shape,scale)))
}


###################################################
### chunk number 53: 
###################################################
totmeansurv = mean((d1+d2)/2)
startvals.GS = list(lscale0=log(totmeansurv),
                  lscale.x2=0,lscale.x3=0,lscale.x4=0,lscale.x5=0,
                  lscale.q=0,lscale.d=0,lscale.i=0,
                  lshape=0)
GS.xqdi = mle2(NLL.GS.xqdi,startvals.GS)


###################################################
### chunk number 54:  eval=FALSE
###################################################
## ## this also works!
## dicweib = function(x,shape,scale,log=FALSE) {
##   if (is.matrix(x)) {
##     day1 = x[,1]
##     day2 = x[,2]
##   } else {
##     day1 = x[1]
##     day2 = x[2]
##   }
##   v = log(pweibull(day2,shape,scale)-pweibull(day1,shape,scale))
##   if (log) v else exp(v)
## }
## 
## fexper = factor(exper)
## mle2(cbind(day1,day2)~dicweib(exp(shape),exp(scale)),
##      parameters=list(scale~fexper+qual*density),
##      start=list(scale=log(totmeansurv),shape=0))


###################################################
### chunk number 55: 
###################################################
summary(GS.xqdi)


###################################################
### chunk number 56: 
###################################################
GS.xqd = mle2(NLL.GS.xqdi,startvals.GS,fixed=list(lscale.i=0))


###################################################
### chunk number 57: 
###################################################
GS.xd = mle2(NLL.GS.xqdi,startvals.GS,
             fixed=list(lscale.i=0,lscale.q=0))
GS.x = mle2(NLL.GS.xqdi,startvals.GS,
            fixed=list(lscale.i=0,lscale.q=0,lscale.d=0))


###################################################
### chunk number 58: 
###################################################
anova(GS.xqdi,GS.xqd,GS.xd,GS.x)


###################################################
### chunk number 59: 
###################################################
GS.qdi = mle2(NLL.GS.xqdi,startvals.GS,fixed=list(lscale.x2=0,lscale.x3=0,lscale.x4=0,lscale.x5=0))
GS.qd = mle2(NLL.GS.xqdi,startvals.GS,
             fixed=list(lscale.i=0,
               lscale.x2=0,lscale.x3=0,lscale.x4=0,lscale.x5=0))

GS.xq = mle2(NLL.GS.xqdi,startvals.GS,
            fixed=list(lscale.i=0,lscale.d=0))
GS.d = mle2(NLL.GS.xqdi,startvals.GS,
            fixed=list(lscale.i=0,lscale.x2=0,lscale.x3=0,lscale.x4=0,
              lscale.x5=0,lscale.q=0))
GS.q = mle2(NLL.GS.xqdi,startvals.GS,
            fixed=list(lscale.i=0,
              lscale.x2=0,lscale.x3=0,lscale.x4=0,lscale.x5=0,lscale.d=0))
GS.0 = mle2(NLL.GS.xqdi,startvals.GS,
            fixed=list(lscale.i=0,lscale.x2=0,lscale.x3=0,lscale.x4=0,
              lscale.x5=0,lscale.d=0,lscale.q=0))


###################################################
### chunk number 60: 
###################################################
a1 = AICtab(GS.xqdi,GS.xqd,GS.qdi,GS.xq,GS.xd,GS.qd,GS.x,GS.q,GS.d,GS.0,delta=TRUE,weights=TRUE)
a2 = AICctab(GS.xqdi,GS.xqd,GS.qdi,GS.xq,GS.xd,GS.qd,GS.x,GS.q,GS.d,GS.0,
              nobs=nrow(GobySurvival),delta=TRUE)
b1 = BICtab(GS.xqdi,GS.xqd,GS.qdi,GS.xq,GS.xd,GS.qd,GS.x,GS.q,GS.d,GS.0,
              nobs=nrow(GobySurvival),delta=TRUE)
atab = cbind(a1$df,a1$dAIC,a1$weight,a2$dAICc,b1$dBIC)
rownames(atab) = gsub("GS\\.","",attr(a1,"row.names"))
## need to double backslashes
colnames(atab) = c("params","$\\Delta \\mbox{AIC}$","AIC weights","$\\Delta \\mbox{AIC}_c$","$\\Delta \\mbox{BIC}$")
atab = round(atab[order(atab[,2]),],3)


###################################################
### chunk number 61: 
###################################################
latex(round(atab,2),file="",table.env=FALSE,title="model")


###################################################
### chunk number 62: 
###################################################
data(SeedPred)


###################################################
### chunk number 63: 
###################################################
SeedPred = na.omit(subset(SeedPred,available>0))
attach(SeedPred)


###################################################
### chunk number 64: 
###################################################
nz = subset(SeedPred,taken>0)


###################################################
### chunk number 65:  eval=FALSE
###################################################
## barchart(table(nz$taken,nz$available,nz$dist,nz$species),stack=FALSE)
## barchart(table(nz$taken,nz$species,nz$dist,nz$available),stack=FALSE)
## barchart(table(nz$species,nz$available,nz$dist,nz$taken),stack=FALSE)
## barchart(table(nz$available,nz$dist,nz$taken),stack=FALSE)
## barchart(table(nz$available,nz$species,nz$taken),stack=FALSE)


###################################################
### chunk number 66: 
###################################################
barchart(table(available,dist,taken),stack=FALSE)


###################################################
### chunk number 67:  eval=FALSE
###################################################
## tcumfac = cut(nz$tcum,breaks=c(0,20,40,60,180))
## barchart(table(nz$available,tcumfac,nz$taken),stack=FALSE)
## barchart(table(available,tcumfac,taken),stack=FALSE)


###################################################
### chunk number 68: 
###################################################
pcomb = table(nz$available,nz$taken)
pcomb = sweep(pcomb,1,rowSums(pcomb),"/")
trellis.par.set(canonical.theme(color=FALSE))
print(barchart(pcomb[-1,],stack=FALSE,auto.key=list(x=0.8,y=0.8,corner=c(1,1),
                                        title="Taken",reverse.rows=TRUE),
               xlab="Frequency",ylab="Seeds available"))


###################################################
### chunk number 69: 
###################################################
dzibinom = function(x,prob,size,zprob,log=FALSE) {
  logv = log(1 - zprob) + 
    dbinom(x, prob = prob, size = size, log = TRUE)
  logv = ifelse(x == 0, log(zprob + exp(logv)), logv)
  if (log) logv  else exp(logv)
}
dzibb = function(x,size,prob,theta,zprob,log=FALSE) {
  logv = ifelse(x>size,
    NA,log(1 - zprob) + 
    dbetabinom(x, prob = prob, size = size, theta=theta, log = TRUE))
  logv = ifelse(x == 0, log(zprob + exp(logv)), logv)
  if (log) logv  else exp(logv)
}


###################################################
### chunk number 70: 
###################################################
SP.zibb = mle2(taken~dzibb(size=available,prob,theta,plogis(logitzprob)),
                start=list(prob=0.5,theta=1,logitzprob=0))
print(SP.zibb)


###################################################
### chunk number 71: 
###################################################
cov2cor(vcov(SP.zibb))


###################################################
### chunk number 72: 
###################################################
SP.bb = mle2(taken~dbetabinom(prob,theta,size=available),
              start=list(prob=0.5,theta=1))


###################################################
### chunk number 73: 
###################################################
logLik(SP.bb)-logLik(SP.zibb)


###################################################
### chunk number 74: 
###################################################
SP.zib = mle2(taken~dzibinom(size=available,prob=p,
                              zprob=plogis(logitzprob)),
               start=list(p=0.2,logitzprob=0))


###################################################
### chunk number 75: 
###################################################
AICtab(SP.zib,SP.zibb,SP.bb,sort=TRUE,weights=TRUE)


###################################################
### chunk number 76: 
###################################################
comb = table(taken,available)
pcomb = sweep(comb,2,colSums(comb),"/")


###################################################
### chunk number 77: 
###################################################
mtab = matrix(0,nrow=6,ncol=5)
for (N in 1:5) {
  cvals = coef(SP.zibb)
  mtab[1:(N+1),N] = dzibb(0:N,size=N,prob=cvals["prob"],
                             theta=cvals["theta"],
                             zprob=plogis(cvals["logitzprob"]))
}


###################################################
### chunk number 78: 
###################################################
mtab2 = matrix(0,nrow=6,ncol=5)
for (av in 1:5) {
  mtab2[1:(av+1),av] = dbetabinom(0:av,prob=coef(SP.bb)["prob"],
                             theta=coef(SP.bb)["theta"],size=av)

}
mtab3 = matrix(0,nrow=6,ncol=5)
for (av in 1:5) {
  mtab3[1:(av+1),av] = 
    dzibinom(0:av,prob=coef(SP.zib)["p"],
             zprob=plogis(coef(SP.zib)["logitzprob"]),
             size=av)

}


###################################################
### chunk number 79: 
###################################################
  tmpf = function(m) {
    m2 = m[-1,]
    sweep(m2,2,colSums(m2),"/")
  }
op = par(mfrow=c(2,2),mar=c(2,2,2,2))
##mosaicplot(tmpf(pcomb),main="data")
barplot(tmpf(pcomb),main="data",axes=FALSE)
barplot(tmpf(mtab),main="Z-I beta-binomial",axes=FALSE)
barplot(tmpf(mtab2),main="beta-binomial",axes=FALSE)
barplot(tmpf(mtab3),main="Z-I binomial",axes=FALSE)
par(op)


###################################################
### chunk number 80: 
###################################################
pval = numeric(5)
for (N in 1:5) {
  obs = comb[1:(N+1),N]
  prob = mtab[1:(N+1),N]
  pval[N] = chisq.test(obs,p=prob)$p.value
}


###################################################
### chunk number 81: 
###################################################
latex(t(matrix(format.pval(pval,digits=2,eps=0.001),
               dimnames=list(1:5,NULL))),file="",title="",table.env=FALSE)


###################################################
### chunk number 82: 
###################################################
SP.bb.dist = mle2(taken~dbetabinom(prob,size=available,theta),
  parameters=list(prob~dist-1,theta~dist-1),start=as.list(coef(SP.bb)))


###################################################
### chunk number 83: 
###################################################
anova(SP.bb,SP.bb.dist)


###################################################
### chunk number 84: 
###################################################
startvals =list(lprob=qlogis(coef(SP.bb.dist)["prob.dist10"]),
  ltheta=log(coef(SP.bb.dist)["theta.dist10"]))


###################################################
### chunk number 85: 
###################################################
SP.bb.dist2 = mle2(taken~dbetabinom(plogis(lprob),exp(ltheta),size=available),
  parameters=list(lprob~dist,ltheta~dist),start=startvals)


###################################################
### chunk number 86: 
###################################################
summary(SP.bb.dist2)


###################################################
### chunk number 87: 
###################################################
SP.bb.probdist = mle2(taken~dbetabinom(plogis(lprob),exp(ltheta),size=available),
  parameters=list(lprob~dist),start=startvals)


###################################################
### chunk number 88: 
###################################################
anova(SP.bb,SP.bb.probdist,SP.bb.dist)
AICtab(SP.bb,SP.bb.probdist,SP.bb.dist,sort=TRUE,weights=TRUE)


###################################################
### chunk number 89: 
###################################################
c1 = coef(SP.bb.probdist)
plogis(c(c1[1],c1[1]+c1[2]))


###################################################
### chunk number 90: 
###################################################
SP.bb.sp = mle2(taken~dbetabinom(plogis(lprob),exp(ltheta),size=available),
  parameters=list(lprob~species,ltheta~species),start=startvals)


###################################################
### chunk number 91: 
###################################################
s1 = summary(SP.bb.sp)
s1@coef = s1@coef[,-3] ## drop z value
printCoefmat(s1@coef,has.Pvalue=TRUE)


###################################################
### chunk number 92: 
###################################################
SP.bb.probsp = mle2(taken~dbetabinom(plogis(lprob),exp(ltheta),size=available),
  parameters=list(lprob~species),start=startvals)


###################################################
### chunk number 93: 
###################################################
anova(SP.bb.sp,SP.bb.probsp,SP.bb)
AICtab(SP.bb.sp,SP.bb.probsp,SP.bb,sort=TRUE,weights=TRUE)


###################################################
### chunk number 94:  eval=FALSE
###################################################
## tmpfun = calc_mle2_function(taken~dbetabinom(plogis(lprob),exp(ltheta),size=available),
##   parameters=list(lprob~species-1),start=startvals)$fn
## pars = coef(SP.bb.probsp)
## pars[2:8] = pars[2:8]+pars[1]
## names(pars) = NULL
## do.call("tmpfun",as.list(pars))


###################################################
### chunk number 95: 
###################################################
SP.bb.probsp0 = mle2(taken~dbetabinom(plogis(lprob),exp(ltheta),size=available),
  parameters=list(lprob~species-1),start=startvals,
  method="L-BFGS-B",lower=rep(-10,9),upper=rep(10,9))


###################################################
### chunk number 96:  eval=FALSE
###################################################
##  lprob0 = coef(SP.bb.probsp)["lprob0"]
## lprobdiff = lprob0+coef(SP.bb.probsp)[2:nsp]
## probvals = c(lprob0,lprob0+lprobdiff)
## predprob = plogis(probvals)


###################################################
### chunk number 97:  eval=FALSE
###################################################
## sdvals = coef(summary(SP.bb.probsp))[,"Std. Error"]
## sdvals0 = sdvals["lprob0"]
## sdvalsdiff = sqrt(sdvals[1]^2+sdvals[2:nsp]^2)
## sdvals = c(sdvals0,sdvalsdiff)


###################################################
### chunk number 98:  eval=FALSE
###################################################
## ci.lo = plogis(probvals-1.96*sdvals)
## ci.hi = plogis(probvals+1.96*sdvals)


###################################################
### chunk number 99: 
###################################################
predprob = plogis(coef(SP.bb.probsp0))[1:8]
SP.bb.ci = plogis(confint(SP.bb.probsp0,method="quad"))[1:8,]


###################################################
### chunk number 100: 
###################################################
op = par(cex=1.5,mar=c(4,4,0,1)+0.1,bty="l",mgp=c(2.5,0.75,0),las=1)
data(SeedPred_mass)
plotCI(SeedPred_mass,predprob,li=SP.bb.ci[,1],ui=SP.bb.ci[,2],pch=NA,gap=0.012,
       xlab="Seed mass",ylab="Removal probability (p)")
textrec(SeedPred_mass,predprob,levels(SeedPred$species),cex=0.75,
        bg="white",expand=1.2)
text(SeedPred_mass,predprob,levels(SeedPred$species),cex=0.75)
box()
par(op)


###################################################
### chunk number 101:  eval=FALSE
###################################################
## NLL.bb.probsptrans = function(p) {
##   nsp = 8
##   lprob0 = p[1]
##   lprobdiff = p[2:nsp]
##   lprobtr = p[nsp+1]
##   lprobtrdiff = p[(nsp+2):(2*nsp)]
##   lprob = c(lprob0,lprob0+lprobdiff,
##              lprob0+lprobtr,lprob0+lprobtr+lprobdiff+lprobtrdiff)
##   probvec = plogis(lprob)
##   ltheta0 = p[2*nsp+1]
##   prob = probvec[interaction(SeedPred$species,SeedPred$dist)]
##   theta = exp(ltheta0)
##   -sum(dbetabinom(taken,size=available,prob=prob,theta=theta,log=TRUE))
## }
## parnames(NLL.bb.probsptrans) = c(parnames(NLL.bb.sp)[1:nsp],
##                            paste("tr",parnames(NLL.bb.sp)[1:nsp],sep=""),
##                            "ltheta")
## svec.probtr = c(svec[1:nsp],rep(0,nsp),svec[nsp+1])
## names(svec.probtr) = parnames(NLL.bb.probsptrans)
## SP.bb.probsptrans = mle2(NLL.bb.probsptrans,start=svec.probtr,vecpar=TRUE)


###################################################
### chunk number 102: 
###################################################
SP.bb.probspdist = mle2(taken~dbetabinom(plogis(lprob),exp(ltheta),size=available),
  parameters=list(lprob~species*dist),start=startvals,
  method="L-BFGS-B",lower=rep(-10,9),upper=rep(5,9))


###################################################
### chunk number 103: 
###################################################
AICtab(SP.bb,SP.bb.probsp,SP.bb.probspdist,SP.bb.sp,SP.bb.probdist,SP.bb.dist,
       weights=TRUE,sort=TRUE)


###################################################
### chunk number 104:  eval=FALSE
###################################################
## mean.prop.taken = tapply(taken/available,tint,mean,na.rm=TRUE)
## sd.prop.taken = tapply(taken/available,tint,sd,na.rm=TRUE)
## n.tint = table(tint)
## se.prop.taken = sd.prop.taken/sqrt(n.tint)


###################################################
### chunk number 105: 
###################################################
par(mfrow=c(1,2))
op = par(cex=1.5,
  mar=c(4,4,0,1)+0.1,bty="l",mgp=c(2.5,0.75,0),
  bty="l",las=1)
nint = table(tint)
mtaken = tapply(taken/available,tint,mean,na.rm=TRUE)
setaken = tapply(taken/available,tint,sd,na.rm=TRUE)/sqrt(nint)
par(mfrow=c(1,2))
plotCI(3:10,mtaken,setaken,xlab="Interval (days)",ylab="Proportion taken")
corner.label2("a","topleft",inset=0.025)
plot(date,tint,xlab="Date",ylab="Interval (days)")
corner.label2("b","topleft",inset=0.025)
par(op)


###################################################
### chunk number 106:  eval=FALSE
###################################################
## NLL.bb.probspdate = function(p) {
##   nsp = 8
##   lprob0 = p[1]
##   lprobdiff = p[2:nsp]
##   lprob = c(lprob0,lprob0+lprobdiff)
##   probvec = plogis(lprob)
##   ltheta0 = p[nsp+1]
##   date = p[nsp+2]
##   prob = probvec[SeedPred$species]*exp(-SeedPred$tcum*date)
##   theta = exp(ltheta0)
##   -sum(dbetabinom(taken,size=available,prob=prob,theta=theta,log=TRUE))
## }
## parnames(NLL.bb.probspdate) = c(parnames(NLL.bb.sp)[1:nsp],"ltheta","date")
## svec.probspdate = c(svec[1:(nsp+1)],0.001)
## names(svec.probspdate) = parnames(NLL.bb.probspdate)
## SP.bb.probspdate = mle2(NLL.bb.probspdate,start=svec.probspdate,vecpar=TRUE)


###################################################
### chunk number 107: 
###################################################
SP.bb.probspdate = mle2(taken~dbetabinom(plogis(lprob)*exp(-tcum*date),
  exp(ltheta),size=available),
  parameters=list(lprob~species),
  start=c(startvals,date=0),
  method="L-BFGS-B",lower=c(rep(-10,9),0),upper=c(rep(5,9),2))


###################################################
### chunk number 108: 
###################################################
SP.bb.probdate = mle2(taken~dbetabinom(plogis(lprob)*exp(-tcum*date),
  size=available,exp(ltheta)),
  start=c(startvals,date=0),
  method="L-BFGS-B",lower=c(rep(-10,2),0),upper=c(rep(5,2),2))


###################################################
### chunk number 109: 
###################################################
op = par(cex=1.5,mar=c(4,4,0,1)+0.1,bty="l",mgp=c(2,0.75,0),las=1)
dates = unique(date)
ndate = table(date)
mtaken = tapply(taken/available,date,mean,na.rm=TRUE)
sdtaken = tapply(taken/available,date,sd,na.rm=TRUE)/sqrt(ndate)
plotCI(dates,mtaken,sdtaken,xlab="Date",ylab="",ylim=c(0.001,0.4))
mtext(side=2,"Proportion taken",line=2.5,cex=1.5,las=0)
t0 = as.numeric(SeedPred$date[1])
p0 = plogis(coef(SP.bb.probdate)["lprob"])
dpar = coef(SP.bb.probdate)["date"]
curve(p0*exp(-dpar*(x-t0)),add=TRUE)
par(op)


###################################################
### chunk number 110: 
###################################################
detach(SeedPred)


