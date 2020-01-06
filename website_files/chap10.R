###################################################
### chunk number 1: 
###################################################
library(emdbook)
library(Hmisc)
library(bbmle)
library(MASS)
source("chapskel.R")


###################################################
### chunk number 2: 
###################################################
nlikfun = function(a,c,d) {
  k = c*exp(d*flowers)
  -sum(dnbinom(seedlings,mu=a,size=k,log=TRUE))
}


###################################################
### chunk number 3:  eval=FALSE
###################################################
## mle2(seedling~dbinom(mu=a,size=c*exp(d*flowers)),...)


###################################################
### chunk number 4: 
###################################################
## interesting but perhaps a waste of time here:
## actually not very much difference between the
## x-dependent k and constant k models.  Also,
## the tail of the seedling distribution seems too
## thin (too few points outside upper tail?? not
## really: 2 vs 4.1 expected above 97.5,
##         7 vs 8.2 expected about 90)
## some kind of Q-Q plot would be possible, but some work ...
data(Lily_sum)
Lily_sum <- subset(Lily_sum,seedlings < 50)
attach(Lily_sum)
likfun <- function(a,b,c,d) {
  mu <- a*exp(b*flowers)
  k <- c*exp(d*flowers)
  -sum(dnbinom(seedlings,mu=mu,size=k,log=TRUE))
}
library(bbmle)
lily.mk <- mle2(minuslogl=likfun,start=list(a=1,b=-0.1,c=1,d=-0.1))
lily.m <- mle2(minuslogl=likfun,start=list(a=1,b=-0.1,c=1,d=-0.1),
           fixed=list(d=0))
lily.k <- mle2(minuslogl=likfun,start=list(a=1,b=-0.1,c=1,d=-0.1),
           fixed=list(b=0))
lily.0 <- mle2(minuslogl=likfun,start=list(a=1,b=-0.1,c=1,d=-0.1),
           fixed=list(b=0,d=0))
## AIC table:
AICtab(lily.0,lily.k,lily.m,lily.mk,sort=TRUE)
anova(lily.mk,lily.k,lily.0) 
pchisq(2*(-logLik(lily.m)+logLik(lily.mk)),1,lower.tail=FALSE)
if (FALSE) {
  library(lattice)
  densityplot(~seedlings|equal.count(flowers),from=0,data=X)
  plot(flowers,seedlings)
  with(as.list(coef(lily.mk)),curve(a*exp(b*x),add=TRUE))
  with(as.list(coef(lily.mk)),curve(qnbinom(0.975,mu=a*exp(b*x),size=c*exp(d*x)),add=TRUE,lty=2,type="s"))
  with(as.list(coef(lily.mk)),curve(qnbinom(0.025,mu=a*exp(b*x),size=c*exp(d*x)),add=TRUE,lty=2,type="s"))
  with(as.list(coef(lily.mk)),curve(qnbinom(0.95,mu=a*exp(b*x),size=c*exp(d*x)),add=TRUE,lty=3,type="s"))
  with(as.list(coef(lily.mk)),curve(qnbinom(0.05,mu=a*exp(b*x),size=c*exp(d*x)),add=TRUE,lty=3,type="s"))
  with(as.list(coef(lily.m)),curve(a*exp(b*x),add=TRUE,col=2))
  with(as.list(coef(lily.m)),curve(qnbinom(0.975,mu=a*exp(b*x),size=c*exp(d*x)),add=TRUE,lty=2,type="s",col=2))
  with(as.list(coef(lily.m)),curve(qnbinom(0.025,mu=a*exp(b*x),size=c*exp(d*x)),add=TRUE,lty=2,type="s",col=2))
} ## eventually have to weed out this junk, leave it in for now
## "logistic slicing"



###################################################
### chunk number 5: 
###################################################
op=par(lwd=2,bty="l",las=1,cex=1.5,mgp=c(2.5,1,0),mar=c(5,4,2,3.5)+0.1)
plot(flowers,jitter(seedlings),xlab="Flowers",ylab="Seedlings",axes=FALSE,lwd=1)
axis(side=1,at=seq(0,120,by=40))
axis(side=2)
box()
bw.exp <- 2  ## how much extra smoothing?
k1 = kde2d(flowers,seedlings,n=100,h=bw.exp*c(bandwidth.nrd(flowers),bandwidth.nrd(seedlings)))
## renormalize k1 by y
k2z <- sweep(k1$z,1,rowSums(k1$z),"/")
means <- rowSums(sweep(k2z,2,k1$y,"*"))
## convert from 
k1$z <- t(apply(k1$z,1,cumsum)) ## need to t()
k1$z <- sweep(k1$z,1,(k1$z[,ncol(k1$z)]),"/")
contour(k1,add=TRUE,col="darkgray",levels=c(0.9,0.95,0.975),labcex=1)
lines(k1$x,means,col="darkgray")
text(107,6.5,"mean",col="darkgray",cex=0.9)
with(as.list(coef(lily.k)),curve(a*exp(b*x),add=TRUE))
with(as.list(coef(lily.k)),curve(qnbinom(0.975,mu=a*exp(b*x),size=c*exp(d*x)),add=TRUE,type="s"))
with(as.list(coef(lily.k)),curve(qnbinom(0.95,mu=a*exp(b*x),size=c*exp(d*x)),add=TRUE,type="s"))
with(as.list(coef(lily.k)),curve(qnbinom(0.9,mu=a*exp(b*x),size=c*exp(d*x)),add=TRUE,type="s"))
## xloc <- 105
par(xpd=NA)
xloc <- par("usr")[2]
qlocs <- c(0.9,0.95,0.975)
tweak=c(-1,0,0.5)
ylocs <- with(as.list(coef(lily.k)),qnbinom(qlocs,mu=a,size=c*exp(d*xloc)))+tweak## +1
text(xloc,coef(lily.k)["a"],"mean",pos=4,cex=0.9) ##+0.5,"mean")
text(xloc,ylocs,paste("q(",as.character(qlocs),")",sep=""),cex=0.9,pos=4)
par(xpd=FALSE)
###
## nonparametric density estimation check -- far more elaborate than necessary
if (FALSE) {  ## add marginal density to plot: how sparse is density at low end?
  par(new=TRUE)
  plot(density(flowers,from=0,adjust=0.3),col=2,ylim=c(0,0.2),xlab="",ylab="",ann=FALSE,
       main="",axes=FALSE)
}
par(op)


###################################################
### chunk number 6: 
###################################################
## HACK/FIXME: temporarily set width narrower
oldwid = options(width=55)


###################################################
### chunk number 7: 
###################################################
data(FirDBHFec)
firdata=na.omit(FirDBHFec)
gnls(TOTCONES~a*DBH^b,data=firdata,weights=varPower(form=~DBH),
     start=list(a=1,b=1))


###################################################
### chunk number 8: 
###################################################
## restore
options(oldwid)


###################################################
### chunk number 9: 
###################################################
M = matrix(c(1,0.9,-0.9,0.9,1,0.9,-0.9,0.9,1),nrow=3)
e1 = eigen(M)$values


###################################################
### chunk number 10:  eval=FALSE
###################################################
## -dmvnorm(z,mu,Sigma=V,log=TRUE)


###################################################
### chunk number 11:  eval=FALSE
###################################################
## Lambda = exp(mvrnorm(1,mu=mu,Sigma=V))
## Y = rpois(length(mu),Lambda)


###################################################
### chunk number 12:  eval=FALSE
###################################################
## genf = function() {rnorm(6,sd=1)+rep(rnorm(2,sd=2),each=3)}
## m = t(replicate(1000,genf()))
## cov2cor(var(m))


###################################################
### chunk number 13: 
###################################################
set.seed(1001)
x.true = rgamma(1000,shape=3,scale=10)
x.obs = rnorm(1000,mean=x.true,sd=10)


###################################################
### chunk number 14:  eval=FALSE
###################################################
## ## all superseded by batch file?
## if (!file.exists("gammanorm1.RData")) {
##   time0 = system.time(m1 <- mle2(minuslogl=getdist,start=list(shape=3,scale=10,sd=1),
##                       data=list(dat=x.obs,debug=TRUE),
##                       method="Nelder-Mead"))
##   time0b = system.time(m1.ci <- confint(m2))
##   save("m1","m2","x.obs","x.true","time0","time0b",file="gammanorm1.RData")
## }
## load("gammanorm1.RData")


###################################################
### chunk number 15: 
###################################################
load("gammanorm-batch.RData")
S4trans <- function(obj1) {
  obj2 <- new(class(obj1))
  for (i in slotNames(obj1)) {
    t1 <- try(slot(obj2,i) <- slot(obj1,i))
  }  
  obj2
}
m1 <- S4trans(m1)


###################################################
### chunk number 16: 
###################################################
citab <- round(data.frame(true=c(3,10,10),coef(m1),confint(m1,method="quad"),
                          m1.ci),2)
colnames(citab) <- c("true","MLE","Quadratic\n2.5\\%","Quadratic\n97.5\\%",
                     "Profile\n2.5\\%","Profile\n97.5\\%")
rownames(citab) <- c("shape ($a$)","scale ($s$)","$\\sigma$") 
latex(citab,file="",table.env=FALSE,title="")


###################################################
### chunk number 17: 
###################################################
op=par(lwd=2,bty="l",las=1,cex=1.5,mgp=c(2.5,1,0),mar=c(5,5.5,2,2)+0.1)
hist(x.obs,prob=TRUE,col="lightgray",breaks=50,main="",ylab="",
     xlab="Observed growth rate",ylim=c(0,0.03))
mtext(side=2,"Probability density",at=0.015,line=3.5,cex=1.5,las=0)
box()
##lines(density(x.true),lwd=2)
curve(dgamma(x,shape=3,scale=10),lty=1,lwd=2,add=TRUE)
with(as.list(coef(m1)),curve(dgamma(x,shape=shape,scale=scale),lty=2,lwd=2,add=TRUE))
curve(dgamma(x,shape=3,scale=10),lty=1,lwd=2,add=TRUE)
curve(dnorm(x,mean=0,sd=10)/2,lty=1,lwd=2,col="darkgray",add=TRUE)
curve(dnorm(x,mean=0,sd=coef(m1)["sd"])/2,lty=2,lwd=2,col="darkgray",add=TRUE)
text(c(35,-9),c(0.029,0.023),c("growth","error"),
     col=c("black","gray"))
legend("topright",c("true","estimated"),lty=1:2,lwd=2,bty="n",cex=0.8)
par(op)


###################################################
### chunk number 18: 
###################################################
library(R2WinBUGS)
library(coda)


###################################################
### chunk number 19:  eval=FALSE
###################################################
## ## hidden code for making bugs work under Wine in my Linux system
## Sys.setenv(WINE="/opt/cxoffice/bin/wine")
## bugs.old <- bugs
## bugs <- function(...) { bugs.old(...,useWINE=TRUE) }


###################################################
### chunk number 20: 
###################################################
pos.obs = x.obs[x.obs>0]
f1 = fitdistr(pos.obs,"gamma")


###################################################
### chunk number 21: 
###################################################
neg.obs = x.obs[x.obs<0]
bineg.obs = c(neg.obs,-neg.obs)
tau0 = 1/var(bineg.obs)


###################################################
### chunk number 22: 
###################################################
tstart = pmax(x.obs,1)


###################################################
### chunk number 23: 
###################################################
clist = as.list(coef(f1))
params = list(tau=tau0,sh=clist$shape,rate=clist$rate,x.true=tstart)
## FIXME: fix perturb.params!
#inits = perturb.params(params,list(tau=c(0.5,1.5),sh=c(0.5,1.5)),mult=TRUE)
inits = rep(list(params),5)
inits[[2]]$tau=0.5*tau0
inits[[3]]$tau=1.5*tau0
inits[[4]]$sh=0.5*clist$shape
inits[[5]]$sh=1.5*clist$shape


###################################################
### chunk number 24: 
###################################################
N = length(x.obs)
data = list("x.obs","N")


###################################################
### chunk number 25: 
###################################################
minval = which.min(x.obs)
minx = x.obs[minval]


###################################################
### chunk number 26: 
###################################################
(minval = which.min(x.obs))
(medval = which.min(abs(x.obs-median(x.obs))))
parameters = c("sh","rate","tau","x.true[670]","x.true[675]")


###################################################
### chunk number 27:  eval=FALSE
###################################################
## time2 <- system.time(gn1.bugs <- bugs(data,inits,parameters.to.save=parameters,
##   model.file="gammanorm.bug",n.chains=length(inits)))


###################################################
### chunk number 28:  eval=FALSE
###################################################
## gn1.bugs = bugs(data,inits,parameters.to.save=parameters,
##                  model.file="gammanorm.bug",n.chains=length(inits))


###################################################
### chunk number 29:  eval=FALSE
###################################################
## save("gn1.bugs","time2",file="gn1-bugs.RData")


###################################################
### chunk number 30: 
###################################################
load("gn1-bugs.RData")


###################################################
### chunk number 31: 
###################################################
s = gn1.bugs$summary[1:5,c("mean","2.5%","97.5%")]
s = data.frame(true=c(3,1/10,1/100,round(x.true[minval],1),
                 round(x.true[medval],1)),s)
colnames(s) <- c("true","mean","2.5\\%","97.5\\%")
rownames(s) <- c("shape","rate","$\\tau$",
                 "min val",
                 "median val")
latex(round(s,3),file="",title="",table.env=FALSE)


###################################################
### chunk number 32:  eval=FALSE
###################################################
## ## FIXME: come up with a way to set "which" in coda plots (less clumsily)
## m1bugs <- as.mcmc(gn1.bugs)
## m2bugs <- lapply(m1bugs,function(x)x[,c(1:5,8)])
## m2bugs <- lapply(m2bugs,function(x) { y <- x; colnames(y)[4:5] <- c("median.point","min.point"); y})
## class(m2bugs) <- "mcmc.list"
## trellis.par.set(canonical.theme(color=FALSE))
## print(densityplot(m2bugs,trace=FALSE))


###################################################
### chunk number 33: 
###################################################
rm(bugs,bugs.old)  ## clean up user defs


###################################################
### chunk number 34:  eval=FALSE
###################################################
## ## binom-Poisson mixture
## rbind(table(factor(rbinom(10000,size=rpois(10000,lambda=1),prob=0.8),levels=0:8))/10000,
##       table(factor(rbinom(10000,size=rpois(10000,lambda=2),prob=0.4),levels=0:8))/10000,
## round(dpois(0:8,lambda=0.8),4))


###################################################
### chunk number 35: 
###################################################
tmpf = function(eps,shape,scale,sd,x) {
    exp(dnorm(eps,mean=0,sd=sd,log=TRUE)+dgamma(x-eps,shape=shape,scale=scale,log=TRUE))
  }


###################################################
### chunk number 36: 
###################################################
tmpf(1,shape=3,scale=10,sd=1,x=x.obs[1])


###################################################
### chunk number 37: 
###################################################
i1 = integrate(f=tmpf,lower=-Inf,upper=Inf,shape=3,scale=10,sd=1,x=x.obs[1]); i1$value


###################################################
### chunk number 38: 
###################################################
tmpf2 = function(x,shape,scale,sd) {
    integrate(f=tmpf,lower=-Inf,upper=Inf,shape=shape,scale=scale,sd=sd,x=x)$value
  }
getdist = function(shape,scale,sd,dat,debug=FALSE) {
  v = -sum(log(sapply(dat,tmpf2,shape=shape,scale=scale,sd=sd)))
  if (debug) cat(shape,scale,sd,v,"\n")
  v
}


###################################################
### chunk number 39: 
###################################################
getdist(shape=3,scale=10,sd=1,dat=x.obs)


###################################################
### chunk number 40:  eval=FALSE
###################################################
## m1 = mle2(minuslogl=getdist,start=list(shape=3,scale=10,sd=1),
##                       data=list(dat=x.obs),method="Nelder-Mead")
## m1.ci = confint(m2)


