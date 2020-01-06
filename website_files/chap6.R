###################################################
### chunk number 1: 
###################################################
library(emdbook)
library(plotrix) ## for corner.loc
## try(detach("package:nlme"))
library(bbmle)
library(lattice)
library(MASS)
library(Hmisc,warn=FALSE)
source("chapskel.R")


###################################################
### chunk number 2: 
###################################################
binomNLL1 = function(p,k,N) {
  -sum(dbinom(k,prob=p,size=N,log=TRUE))
}


###################################################
### chunk number 3: 
###################################################
data(ReedfrogPred)
x = subset(ReedfrogPred,pred=="pred" & density==10 & size=="small")
k = x$surv


###################################################
### chunk number 4: 
###################################################
op=par(mfrow=c(1,2),lwd=2,bty="l",las=1,cex=1.5)
p.mle <- sum(k)/40
mlik <- exp(sum(dbinom(k,size=10,prob=p.mle,log=TRUE)))
pvec <- seq(0,1,length=100)
liks <- sapply(pvec,function(x)sum(dbinom(k,size=10,prob=x,log=TRUE)))
plot(pvec,exp(liks),type="l",xlab="Predation probability\nper capita",
     ylab="Likelihood",log="y",ylim=c(1e-20,1),axes=FALSE)
axis(side=1,at=seq(0,1,by=0.25))
## if axTicks worked ...
loc2 = seq(-20,0,by=5)
axis(side=2,at=10^loc2,
     do.call("c",lapply(loc2,
            function(z)substitute(expression(10^x),list(x=z)))))
abline(v=p.mle,lty=2)
abline(h=mlik,lty=2)
par(xpd=NA)
text(p.mle,10^(par("usr")[4]+0.5),
     expression(hat(p)==0.75)) 
text(0,mlik*10,expression(L[max]==5.1 %*% 10^-4),adj=0)
box()
par(xpd=TRUE)
plot(table(factor(k,levels=0:10))/4,xlab="# of successes",
     ylab="Probability")
points(0:10,dbinom(0:10,size=10,prob=0.75),pch=16,col="darkgray")
par(op)


###################################################
### chunk number 5: 
###################################################

ptry <- function() {
  x <- rbinom(4,prob=0.75,size=10)
  sum(dbinom(x,prob=0.75,size=10,log=TRUE))
}
rprobs <- replicate(1000,ptry())
binom.pval <- sum(rprobs<log(mlik))/1000


###################################################
### chunk number 6: 
###################################################
set.seed(1001)
npts <- 50
r.true <- 0.008 ## 5% for every 10 m
d.true <- 8
k.true <- 2
poisdatamean= d.true*exp(-r.true*20)
poisdata <- rpois(npts,poisdatamean)
nbdata <- rnbinom(npts,mu=d.true*exp(-r.true*20),size=k.true)
pdata.x <- runif(npts,0,100)
det.y <- d.true*exp(-r.true*pdata.x)
pdata.y <-  rpois(npts,det.y)


###################################################
### chunk number 7: 
###################################################
O1=optim(fn=binomNLL1,par=c(p=0.5),N=10,k=k,
  method="BFGS")


###################################################
### chunk number 8: 
###################################################
O1$par
exp(-O1$value)


###################################################
### chunk number 9: 
###################################################
library(bbmle)
m1=mle2(minuslogl=binomNLL1,start=list(p=0.5),data=list(N=10,k=k))
m1


###################################################
### chunk number 10: 
###################################################
mle2(k~dbinom(prob=p,size=10),start=list(p=0.5))


###################################################
### chunk number 11: 
###################################################
rm(k,x)


###################################################
### chunk number 12: 
###################################################
data(MyxoTiter_sum)
myxdat = subset(MyxoTiter_sum,grade==1)


###################################################
### chunk number 13: 
###################################################
gammaNLL1 = function(shape,scale) {
  -sum(dgamma(myxdat$titer,shape=shape,scale=scale,log=TRUE))
}


###################################################
### chunk number 14: 
###################################################
gm = mean(myxdat$titer)
cv = var(myxdat$titer)/mean(myxdat$titer)


###################################################
### chunk number 15:  eval=FALSE
###################################################
## m3 = mle2(gammaNLL1,
##   start=list(shape=gm/cv,scale=cv))


###################################################
### chunk number 16: 
###################################################
m3 = mle2(gammaNLL1,
  start=list(shape=45.8,scale=0.151))


###################################################
### chunk number 17: 
###################################################
m3


###################################################
### chunk number 18:  eval=FALSE
###################################################
## m3 = mle2(myxdat$titer~dgamma(shape,scale=scale),
##   start=list(shape=gm/cv,scale=cv))


###################################################
### chunk number 19: 
###################################################
f1 =fitdistr(myxdat$titer,"gamma")


###################################################
### chunk number 20: 
###################################################
op=par(mfrow=c(1,2),lwd=2,bty="l",las=1,cex=1.5)
v = curve3d(gammaNLL1(x,y),from=c(25,0.05),to=c(85,0.3),
        sys3d="contour",
  xlab="Shape",ylab="Scale",drawlabels=FALSE,axes=FALSE)
axis(side=2)
axis(side=1,at=c(30,50,70))
box()
contour(v$x,v$y,v$z,add=TRUE,levels=seq(0,190,by=20),col="gray",
        drawlabels=FALSE)
## persp3d(v$x,v$y,-v$z,col="blue")
cmle <- coef(m3) 
points(cmle[1],cmle[2])
text(cmle[1],cmle[2]+0.02,"MLE")
plot(density(myxdat$titer),main="",xlab="Virus titer",
     ylab="Probability density",zero.line=FALSE,col="darkgray")
n <- nrow(myxdat)
points(myxdat$titer,runif(n,0,0.02))
curve(dgamma(x,shape=cmle["shape"],scale=cmle["scale"]),add=TRUE)
curve(dnorm(x,mean(myxdat$titer),sd(myxdat$titer)),lty=2,add=TRUE)
text(7.8,0.45,"density",col="darkgray",adj=0)
text(4.2,0.38,"Gamma")
arrows(4.43,0.36,5.9,0.32,angle=15)
par(xpd=NA)
text(8.5,0.33,"normal",adj=0)
arrows(8.4,0.33,7.7,0.33,angle=15)
par(op)


###################################################
### chunk number 21: 
###################################################
x = subset(ReedfrogPred,pred=="pred" & density==10 & size=="small")
k = x$surv


###################################################
### chunk number 22: 
###################################################
a=121
b=81


###################################################
### chunk number 23: 
###################################################
op=par(lwd=2,bty="l",las=1,cex=1.5)
curve(dbeta(x,shape1=sum(k)+1,shape2=40-sum(k)+1),
      xlab="Predation probability\nper capita",
      ylab="Probability density",
      from=0,to=1,ylim=c(0,13))
curve(dbeta(x,shape1=sum(k),shape2=40-sum(k)),
      add=TRUE,lwd=3)
curve(dbeta(x,shape1=1,shape2=1),col="darkgray",add=TRUE,
      type="s")
curve(dbeta(x,shape1=121,shape2=81),
      add=TRUE,lty=2,col="darkgray")
curve(dbeta(x,shape1=151,shape2=91),
      add=TRUE,lty=2,n=200)
tlocs <- list(x=c(0.44,0.13,0.82,0.75,0.95),
              y=c(10,3.1,10.8,8,5.9))
text(tlocs$x,tlocs$y,c("prior\n(121,81)",
                       "prior\n(1,1)",
                       "posterior\n(151,111)",
                       "posterior\n(31,11)",
                       "scaled\nlikelihood"),
     cex=0.6)
alocs <- list(x=c(0.151,0.464,0.734,0.720,0.924,
                  0.18, 0.547,0.656,0.725,0.843),
              y=c(2.3,9.047,10.806,7.241,4.833,
                  1.02,7.195,9.973,5.898,3.212))
arrows(alocs$x[1:5],alocs$y[1:5],alocs$x[6:10],alocs$y[6:10],
       angle=15,col=rep(c("darkgray","black"),c(2,3)))
par(op)


###################################################
### chunk number 24: 
###################################################
rm(k)


###################################################
### chunk number 25: 
###################################################
prior.as = function(a,s) {
  dgamma(a,shape=0.01,scale=100)*
  dgamma(s,shape=0.1,scale=10)
}
unscaled.posterior = function(a,s) {
  prior.as(a,s)*exp(-gammaNLL1(shape=a,scale=s))
}


###################################################
### chunk number 26: 
###################################################
## (log) prior for a
aprior = function(a) dgamma(a,shape=0.01,scale=100,log=TRUE)
## (log) prior for s
sprior = function(s) dgamma(s,shape=0.1,scale=10,log=TRUE)
## library(adapt)
## calc mean value (and find mode)
postfun = function(p,log=FALSE) {
  a <- p[1]
  s <- p[2]
  v <- aprior(a)+sprior(s)+sum(dgamma(myxdat$titer,shape=a,scale=s,log=TRUE))
  if (log) v else exp(v)
}
prifun = function(p) {
  a <- p[1]
  s <- p[2]
  exp(aprior(a)+sprior(s))
}
unsc.post <- curve3d(postfun(c(x,y)),from=c(10,0.07),
                     n=c(91,101),
                     to=c(100,0.5),sys3d="none")
avec <- unsc.post$x
svec <- unsc.post$y
unsc.post$z <- unsc.post$z/sum(unsc.post$z)
amat <- avec[row(unsc.post$z)]
smat <- svec[col(unsc.post$z)]
amean <- sum(amat*unsc.post$z)
smean <- sum(smat*unsc.post$z)
amarg <- rowSums(unsc.post$z)
smarg <- colSums(unsc.post$z)
amean2 <- sum(amarg*avec)
smean2 <- sum(smarg*svec)
amode <- amat[unsc.post$z==max(unsc.post$z)]       
smode <- smat[unsc.post$z==max(unsc.post$z)]       
amode2 <- avec[which.max(amarg)]
smode2 <- svec[which.max(smarg)]


###################################################
### chunk number 27: 
###################################################
op=par(lwd=2,bty="l",las=1,cex=1.5,mar=c(4,4,2,2)+0.1)
nf <- layout(matrix(c(2,1,0,3), 2, 2, byrow=TRUE),
             widths=c(1,3),heights=c(3,1))
## don't know why this needs to be repeated ... but it does
par(cex=1.5)
contour(unsc.post$x,unsc.post$y,
        unsc.post$z,
        levels=10^seq(-10,-2),drawlabels=FALSE,
        xlab="Shape",ylab="",lwd=1,col="darkgray")
mtext("Scale",side=2,at=0.3,line=2.5,cex=1.5,las=0)
points(amean,smean,pch=1)
points(amode,smode,pch=2)
arrows(c(41,56),c(0.27,0.18),
       c(45,48.5),c(0.175,0.15),angle=15)
text(c(41,56),c(0.28,0.19),
     c("mean","mode"))
par(mar=c(4,1,2,2)+0.1)
plot(-smarg,svec,type="l",axes=FALSE,xlab="",ylab="")
axis(side=4)
axis(side=1,at=c(-0.04,0),labels=c(0.04,0))
## have to hand-draw box
u <- par("usr")
segments(c(u[1],u[2]),
         c(u[3],u[3]),
         c(u[2],u[2]),
         c(u[3],u[4]),lwd=2)
par(mar=c(2,4,0.5,2)+0.1)
plot(avec,amarg,type="l",xlab="",ylab="",axes=FALSE,
     ylim=c(0,0.04))
axis(side=1)
axis(side=2,at=c(0,0.04))
box()
par(op)


###################################################
### chunk number 28:  eval=FALSE
###################################################
## ## junk
## zz <- log(unsc.post$z)
## zz[!is.finite(zz)] <- -750
## ## library(rgl); persp3d(unsc.post$x,unsc.post$y,zz,col="red")
## unsc.post.log <- curve3d(postfun(c(x,y),log=TRUE),
##                          from=c(1,0.01),to=c(100,1),n=c(71,71),
##                          sys3d="rgl",col="blue")
## persp3d(unsc.post.log$x,unsc.post.log$y,exp(unsc.post.log$z),col="red")
## x = unsc.post.log$x
## y = unsc.post.log$y
## image(x,y,log(-unsc.post.log$z))
## curve(0.003056+6.808135/x,add=TRUE)
## contour(x,y,-unsc.post.log$z,levels=49,add=TRUE)
## peaks = which(-unsc.post.log$z<48,arr.ind=TRUE)
## points(x[peaks[,1]],y[peaks[,2]],col="blue",cex=0.5,pch=16)
## ##
## avec = 1:100
## svec = 0.003056+6.808135/avec
## g = mapply(gammaNLL1,avec,svec)
## plot(g,ylim=c(35,50))
## g2 = -mapply(function(x,y)postfun(c(x,y),log=TRUE),avec,svec)
## plot(g2,ylim=c(46,50))
## ##
## plot2 <- curve3d(log(-postfun(c(x,y),log=TRUE)),
##                  from=c(20,0.1),to=c(60,0.2),n=c(71,71),
##                          sys3d="image")
## contour(plot2$x,plot2$y,plot2$z,levels=seq(3.8,3.9,by=0.01),
##         add=TRUE)
## curve(0.003056+6.808135/x,add=TRUE,col="blue",lwd=2)
## 
## numfuna0 = function(a,s,x) {
##   exp(aprior(a)+sprior(s)+dgamma(x,shape=a,scale=s,log=TRUE)+log(a))
## }
## numfuns0 = function(p) {
##   a <- p[1]
##   s <- p[2]
##   exp(aprior(a)+sprior(s)+dgamma(x,shape=a,scale=s,log=TRUE)+log(s))
## }
## numfunas0 = function(a,s,x) {
##   exp(aprior(a)+sprior(s)+dgamma(x,shape=a,scale=s,log=TRUE)+log(s)+log(a))
## }
## numfuna <- function(a) {
##   sapply(a,numfuna0)
## }
## numfuns <- function(s) {
##   sapply(a,numfuns0)
## }
## numa = adapt(2,lower=c(0,0),upper=c(200,200),minpts=10000,functn=numfuns0)
## subdiv=10000
## limits=c(0,15)
## bval2 =integrate(numfun,lower=limits[1],upper=limits[2],subdivisions=subdiv)$value/
##   integrate(denomfun,lower=limits[1],upper=limits[2],subdivisions=subdiv)$value


###################################################
### chunk number 29: 
###################################################
data(ReedfrogFuncresp)
attach(ReedfrogFuncresp)
binomNLL2 = function(p,N,k) {
  a = p[1]
  h = p[2]
  predprob = a/(1+a*h*N)
  -sum(dbinom(k,prob=predprob,size=N,log=TRUE))
}
O2 = optim(fn=binomNLL2,par=c(a=0.5,h=0.0125),N=Initial,k=Killed)
parnames(binomNLL2) = c("a","h")
m2 = mle2(binomNLL2,start=c(a=0.5,h=0.0125),data=list(N=Initial,k=Killed))
gammaNLL2 = function(a,b,shape) {
  meantiter = a*myxdat$day*exp(-b*myxdat$day)
  -sum(dgamma(myxdat$titer,shape=shape,scale=meantiter/shape,log=TRUE))
}
m4 = mle2(gammaNLL2,start=list(a=1,b=0.2,shape=50),
  method="Nelder-Mead")
detach(ReedfrogFuncresp)


###################################################
### chunk number 30: 
###################################################
op=par(mfrow=c(1,2),lwd=2,bty="l",las=1,cex=1.5)
## frog data
plot(ReedfrogFuncresp$Initial,ReedfrogFuncresp$Killed,
     xlab="Initial density",ylab="Number killed")
with(as.list(coef(m2)),
     {
       curve(a*x/(1+a*h*x),add=TRUE)
       curve(qbinom(0.025,size=floor(x),prob=a/(1+a*h*x)),lty=2,add=TRUE,
             type="s")
       curve(qbinom(0.975,size=floor(x),prob=a/(1+a*h*x)),lty=2,add=TRUE,
             type="s")
     })
plot(myxdat$day,myxdat$titer,
     xlab="Day since infection",ylab="Virus titer",xlim=c(0,10),
     ylim=c(0,9))
with(as.list(coef(m4)),
     {
       curve(a*x*exp(-b*x),add=TRUE)
       curve(qgamma(0.025,shape=shape,scale=a*x*exp(-b*x)/shape),
                    lty=2,add=TRUE)
       curve(qgamma(0.975,shape=shape,scale=a*x*exp(-b*x)/shape),
                    lty=2,add=TRUE)
     })
par(op)


###################################################
### chunk number 31: 
###################################################
binomNLL2 = function(p,N,k) {
  a = p[1]
  h = p[2]
  predprob = a/(1+a*h*N)
  -sum(dbinom(k,prob=predprob,size=N,log=TRUE))
}


###################################################
### chunk number 32: 
###################################################
data(ReedfrogFuncresp)
attach(ReedfrogFuncresp)
O2 = optim(fn=binomNLL2,par=c(a=0.5,h=0.0125),N=Initial,k=Killed)


###################################################
### chunk number 33: 
###################################################
parnames(binomNLL2) = c("a","h")
m2 = mle2(binomNLL2,start=c(a=0.5,h=0.0125),data=list(N=Initial,k=Killed))
m2


###################################################
### chunk number 34: 
###################################################
z = as.list(coef(m2))
prob = with(z,a/(1+a*h*50))
pval = pbinom(5,size=50,prob=prob)


###################################################
### chunk number 35: 
###################################################
gammaNLL2 = function(a,b,shape) {
  meantiter = a*myxdat$day*exp(-b*myxdat$day)
  -sum(dgamma(myxdat$titer,shape=shape,scale=meantiter/shape,log=TRUE))
}


###################################################
### chunk number 36: 
###################################################
m4 = mle2(gammaNLL2,start=list(a=1,b=0.2,shape=50),
  method="Nelder-Mead")
m4


###################################################
### chunk number 37: 
###################################################
mle2(titer~dgamma(shape,scale=a*day*exp(-b*day)/shape),
     start=list(a=1,b=0.2,shape=50),data=myxdat,
     method="Nelder-Mead")


###################################################
### chunk number 38: 
###################################################
library(R2WinBUGS)


###################################################
### chunk number 39: 
###################################################
initial <- ReedfrogFuncresp$Initial
killed <- ReedfrogFuncresp$Killed
n <- length(initial)
inits <- list(list(a=0.5,h=0.02),list(a=1,h=0.015),
              list(a=0.1,h=0.5))
frogpred1.bugs <- bugs(data=list("initial","killed","n"),
                       inits,parameters.to.save=c("a","h"),
                       model.file="frogpred.bug",
                       n.chains=length(inits),n.iter=6000,n.thin=1)


###################################################
### chunk number 40: 
###################################################
library(R2WinBUGS)


###################################################
### chunk number 41: 
###################################################
titer = myxdat$titer
day = myxdat$day
n = length(titer)


###################################################
### chunk number 42: 
###################################################
inits <- list(list(a=4,b=0.2,shape=90),
              list(a=1,b=0.1,shape=50),
              list(a=8,b=0.4,shape=150))


###################################################
### chunk number 43: 
###################################################
myxo1.bugs <- bugs(data=list("titer","day","n"),
                       inits,parameters.to.save=c("a","b","shape"),
                       model.file="myxo1.bug",
                   n.chains=length(inits),n.iter=3000)
  


###################################################
### chunk number 44: 
###################################################
detach(ReedfrogFuncresp)  ## why??
mvals <- signif(myxo1.bugs$summary[,"mean"],3)


###################################################
### chunk number 45: 
###################################################
x = subset(ReedfrogPred,pred=="pred" & density==10 & size=="small")
k = x$surv
op=par(mfrow=c(1,2),lwd=2,bty="l",las=1,cex=1.5)
pchs = c(1,2,5)
plot(pvec,-liks,type="l",xlab="Predation probability\nper capita (p)",
     ylab="Negative log-likelihood",ylim=c(5,30))
p.conf <- confint(m1,quietly=TRUE)
v <- c(p.mle,p.conf)
h <- sapply(c(p.mle,p.conf),binomNLL1,k=k,N=10)
points(v,h,pch=pchs)
abline(v=v[2:3],lty=3)
abline(h=h[2],lty=3)
corner.label("a")
####
plot(table(factor(k,levels=0:10))/4,xlab="Tadpoles eaten",
     ylab="Probability",lwd=4,col="gray")
points(0:10,dbinom(0:10,size=10,prob=0.75),pch=pchs[1],
       type="b")
points(0:10,dbinom(0:10,size=10,prob=p.conf[1]),pch=pchs[2],
       type="b")
points(0:10,dbinom(0:10,size=10,prob=p.conf[2]),pch=pchs[3],
       type="b")
text(c(3.4,7.2,10.6),c(0.18,0.31,0.34),
     paste("p=",round(c(p.conf[1],0.75,p.conf[2]),2),sep=""),
     xpd=NA,cex=0.8)
corner.label("b")
par(op)


###################################################
### chunk number 46:  eval=FALSE
###################################################
## attach(ReedfrogFuncresp)
## Lvec = seq(5,9,by=0.1) 
## likfun = function(lambda) { -sum(dpois(poisdata,lambda,log=TRUE)) }
## likcurve = sapply(Lvec,likfun)
## m0 = mle2(minuslogl=likfun,start=list(lambda=4))
## prof0 = profile(m2)
## c95 = confint(prof0,level=0.95)
## c99 = confint(prof0,level=0.99)
## op=par(mfrow=c(1,2))
## plot(table(factor(poisdata,levels=0:14))/length(poisdata),
##      ylab="Probability",xlab="Clutch size")
## points(0:14,dpois(0:14,coef(m0)),type="b",pch=16,lwd=2,add=TRUE)
## points(0:14,dpois(0:14,poisdatamean),type="b",pch=1,lty=2,lwd=2,add=TRUE)
## points(0:14,dpois(0:14,c99[2]),type="b",pch=2,lty=3,add=TRUE)
## points(0:14,dpois(0:14,c99[1]),type="b",pch=3,lty=4,add=TRUE)
## legend(8,0.2,
##        c("best estimate",
##          "true",
##          "1%","99%"),
##        pch=c(16,2,2,3),
##        lwd=c(2,2,1,1),
##        lty=1:4)
## plot(Lvec,likcurve,type="l",xlab=expression(paste("Estimated "),lambda),
##      ylab="Negative log likelihood")
## abline(h=-logLik(m0))
## abline(v=coef(m0),lwd=2)
## abline(v=c95,lty=2)
## abline(v=confint(prof0,level=0.99),lty=3)
## abline(h=qchisq(c(0.95,0.99),df=1)/2-logLik(m0),lty=c(2,3))
## abline(v=poisdatamean,lty=2,lwd=2)
## points(c(coef(m0),poisdatamean,c99[1],c99[2]),
##        sapply(c(coef(m0),poisdatamean,c99[1],c99[2]),likfun),
##        pch=c(16,1,2,3))
## par(op)
## detach(ReedfrogFuncresp)


###################################################
### chunk number 47: 
###################################################
rvec <- seq(0,0.025,length=75)
dvec <- seq(4,15,length=75)
prof <- profile(m1)
## will generate errors -- OK
## profwide <- profile(m1,zmax=8)
rm(x)
rm(k)


###################################################
### chunk number 48: 
###################################################
attach(ReedfrogFuncresp)
bsurf = curve3d(binomNLL2(c(x,y),N=ReedfrogFuncresp$Initial,k=ReedfrogFuncresp$Killed),
  sys="none",from=c(0.3,0.002),to=c(0.75,0.03),n=c(91,91))
prof <- profile(m2)
conflim95 <- confint(prof,level=0.95)
conflim99 <- confint(prof,level=0.99)
m2.p <- coef(m2)
p2 <- profile(m2,which="a")
p2h <- profile(m2,which="h")
detach(ReedfrogFuncresp)


###################################################
### chunk number 49: 
###################################################
op=par(lwd=2,bty="l",las=1,cex=1.5,mar=c(5,5,4,2)+0.1)
image(bsurf$x,bsurf$y,log(bsurf$z),xlab="Attack rate (a)",
      ylab="",col=gray((20:0)/20))
mtext("Handling time (h)",side=2,at=0.017,line=3.5,cex=1.5,las=0)
points(m2.p[1],m2.p[2],pch=16)
#contour(bsurf$x,bsurf$y,bsurf$z,add=TRUE,drawlabels=FALSE,col="gray")
contour(bsurf$x,bsurf$y,bsurf$z,levels=-logLik(m2)+qchisq(0.95,2)/2,add=TRUE,
        drawlabels=FALSE)
contour(bsurf$x,bsurf$y,bsurf$z,levels=-logLik(m2)+qchisq(0.95,1)/2,add=TRUE,lty=3,
        drawlabels=FALSE)
lines(prof@profile$a$par.vals[,"a"],prof@profile$a$par.vals[,"h"],lty=2)
lines(prof@profile$h$par.vals[,"a"],prof@profile$h$par.vals[,"h"],lty=4)
abline(v=conflim95["a",],lty=3,col="darkgray")
abline(h=conflim95["h",],lty=3,col="darkgray")
abline(h=coef(m2)["h"],lty=5,col="darkgray")
text(c(0.65,0.74),c(0.028,0.025),c("h","a"))
text(c(0.5,0.59),c(0.008,0.013),c("univariate","bivariate"),adj=0,cex=0.75)
arrows(0.495,0.008,0.475,0.009,angle=25,length=0.1)
pt3 = prof@profile$a$par.vals[2,]
pt3 = c(pt3,binomNLL2(pt3,N=ReedfrogFuncresp$Initial,k=ReedfrogFuncresp$Killed))
pt4 = c(pt3[1],coef(m2)["h"])
pt4 = c(pt4,binomNLL2(pt4,N=ReedfrogFuncresp$Initial,k=ReedfrogFuncresp$Killed))
points(pt3["a"],pt3["h"],pch=2)
points(pt4["a"],pt4["h"],pch=8)
text(0.3,coef(m2)["h"]+0.0015,"slice",adj=0)
par(op)


###################################################
### chunk number 50: 
###################################################
profpic <- function(prof,param,trueval,best,which,alpha=c(0.95,0.99),
                    legend=TRUE,legpos=NULL,ylab="Negative log likelihood",...) {
  prof1 = prof@profile[[param]]
  nll = (prof1$z)^2/2+prof@summary@m2logL/2
  plot(prof1$par.vals[,param],nll,type="l",xlab=param,ylab=ylab,...)
  abline(h=-logLik(best))
  abline(v=coef(best)[param],lwd=2)
  if (!missing(trueval)) abline(v=trueval,lwd=2,lty=2)
  nalpha = length(alpha)
  for (i in seq(along=alpha)) {
    crit <- qchisq(alpha[i],1)/2
    abline(h=-logLik(m1)+crit,lty=i+1)
    conflim <- confint(prof,parm=param,level=alpha[i])
    abline(v=conflim,lty=i+1)
  }
  if (legend) {
    if (is.null(legpos)) legpos <- corner.loc(1,1,xoff=0.2)
    legend(legpos$x,legpos$y,
         c(paste("alpha=",alpha,sep=""),
           "true",
           "observed"),
           lty=c(1:nalpha,2,1),
           lwd=c(rep(1,nalpha),2,2),bg="white")
  }
}


###################################################
### chunk number 51: 
###################################################
sqrprofplot <- function (x, levels, conf = c(99, 95, 90, 80, 50)/100, nseg = 50,
          plot.confstr = FALSE, confstr = NULL, absVal = TRUE, sqrVal=FALSE, add = FALSE,
          col.minval="green", lty.minval=2,
          col.conf="magenta", lty.conf=2,
          col.prof="blue", lty.prof=1,
          xlabs=nm, ylab="score", xlim, ylim, ...)
{
    ## Plot profiled likelihood
    ## Based on profile.nls (package stats)
    obj <- x@profile
    confstr <- NULL
    if (missing(levels)) {
        levels <- sqrt(qchisq(pmax(0, pmin(1, conf)), 1))
        confstr <- paste(format(100 * conf), "%", sep = "")
    }
    if (any(levels <= 0)) {
        levels <- levels[levels > 0]
        warning("levels truncated to positive values only")
    }
    if (is.null(confstr)) {
        confstr <- paste(format(100 * pchisq(levels^2, 1)), "%", sep = "")
    }
    mlev <- max(levels) * 1.05
    nm <- names(obj)
    ##    opar <- par(mar = c(5, 4, 1, 1) + 0.1)
    for (i in seq(along = nm)) {
      ## <FIXME> This does not need to be monotonic
      ## cat("**",i,obj[[i]]$par.vals[,i],obj[[i]]$z,"\n")
      ## browser()
      yvals <- obj[[i]]$par.vals[,nm[i],drop=FALSE]
      sp <- splines::interpSpline(yvals, obj[[i]]$z,
                                  na.action=na.omit)
      bsp <-splines:: backSpline(sp)
      ## </FIXME>
      x0 <- predict(bsp,0)$y
      if (missing(xlim)) xlim <- predict(bsp, c(-mlev, mlev))$y
      if (is.na(xlim[1]))
        xlim[1] <- min(yvals)
      if (is.na(xlim[2]))
        xlim[2] <- max(yvals)
      if (sqrVal) {
        if (!add) {
          if (missing(ylim)) ylim <- c(0,mlev^2)
          plot(I(z^2) ~ par.vals[, i], data = obj[[i]], xlab = xlabs[i],
               ylim = ylim, xlim = xlim, ylab = ylab,
               type = "n", ...)
        }
        avals <- rbind(as.data.frame(predict(sp)),
                       data.frame(x = drop(yvals),
                                  y = obj[[i]]$z))
        avals$y <- avals$y^2
        lines(avals[order(avals$x), ], col = col.prof, lty=lty.prof)
      } else {
        if (absVal) {
          if (!add) {
            if (missing(ylim)) ylim <- c(0,mlev)
            plot(abs(z) ~ par.vals[, i], data = obj[[i]], xlab = nm[i],
                 ylim = ylim, xlim = xlim, ylab = expression(abs(z)),
                 type = "n", ...)
          }
          avals <- rbind(as.data.frame(predict(sp)),
                         data.frame(x = yvals,
                                    y = obj[[i]]$z))
          avals$y <- abs(avals$y)
          lines(avals[order(avals$x), ], col = col.prof, lty=lty.prof)
        } else {
          if (!add) {
            if (missing(ylim)) ylim <- c(-mlev,mlev)
            plot(z ~ par.vals[, i], data = obj[[i]], xlab = nm[i],
                 ylim = ylim, xlim = xlim, ylab = expression(z),
                 type = "n", ...)
          }
          lines(predict(sp), col = col.prof, lty=lty.prof)
        }
      }
      abline(v = x0, h=0, col = col.minval, lty = lty.minval)
      for (i in 1:length(levels)) {
        lev <- levels[i]
        confstr.lev <- confstr[i]
        ## Note: predict may return NA if we didn't profile
        ## far enough in either direction. That's OK for the
        ## "h" part of the plot, but the horizontal line made
        ## with "l" disappears.
        pred <- predict(bsp, c(-lev, lev))$y
        ## horizontal
        if (absVal || sqrVal) levs=rep(lev,2) else levs=c(-lev,lev)
        if (sqrVal) lines(pred, levs^2, type = "h", col = col.conf, lty = lty.conf)
        else lines(pred, levs, type = "h", col = col.conf, lty = 2)
        ## vertical
        pred <- ifelse(is.na(pred), xlim, pred)
        if (sqrVal) lines(pred, rep(lev^2,2), type = "l", col = col.conf, lty = lty.conf)
        else if (absVal) lines(pred, rep(lev, 2), type = "l", col = col.conf, lty = lty.conf)
        else {
          lines(c(x0,pred[2]), rep(lev, 2), type = "l", col = col.conf, lty = lty.conf)
          lines(c(pred[1],x0), rep(-lev, 2), type = "l", col = col.conf, lty = lty.conf)
        }
        if (plot.confstr) {
          if (sqrVal) text(labels=confstr.lev,x=x0,y=lev^2,col=col.conf)
          else text(labels=confstr.lev,x=x0,y=lev^2,col=col.conf)
        }
      }
    }
    ## par(opar)
  }

op=par(lwd=2,bty="l",las=1,cex=1.5)
## attach(ReedfrogFuncresp)
sqrprofplot(p2,sqrVal=TRUE,axes=FALSE,conf=0.95,col.conf="gray",
     col.prof="black",col.minval=NA,ylim=c(0,40),xlim=c(0.3,0.75),
     xlab="Attack rate (a)",
     ylab="Negative log-likelihood")
axis(side=1)
## convert scale (deviance) 
##  scale = 2*(L-Lmin)
##  L = scale/2+Lmin
ylocs <- c(47,seq(50,75,by=5))
axis(side=2,at=2*(ylocs+logLik(m2)),labels=ylocs)
box()
aslice = sapply(bsurf$x,
  function(a)binomNLL2(c(a,coef(m2)["h"]),
                       N=ReedfrogFuncresp$Initial,
                       k=ReedfrogFuncresp$Killed))
lines(bsurf$x,2*(aslice+logLik(m2)),lty=3,xpd=NA)
points(coef(m2)["a"],0,pch=16)
points(pt3["a"],2*(pt3[3]+as.numeric(logLik(m2))),pch=2)
points(pt4["a"],2*(pt4[3]+as.numeric(logLik(m2))),pch=8)
par(xpd=NA)
text(c(0.36,0.72),c(28,11),c("slice","profile"),cex=1,
     adj=0)
par(xpd=FALSE)
##detach(ReedfrogFuncresp)
par(op)


###################################################
### chunk number 52: 
###################################################
op=par(lwd=2,bty="l",las=1,cex=1.5)
## frog data
par(mar=c(5,4,2,4)+0.1)
plot(ReedfrogFuncresp$Initial,ReedfrogFuncresp$Killed,
     xlab="Initial density",ylab="Number killed")
tmpf <- function(pars,...) {
  with(as.list(pars),curve(a*x/(1+a*h*x),add=TRUE,...))
}
tmpf(coef(m2))
tmpf(pt3[1:2],lty=2)
tmpf(pt4[1:2],lty=3)
avals <- c(coef(m2)["a"],pt3["a"],pt4["a"])
hvals <- c(coef(m2)["h"],pt3["h"],pt4["h"])
hts <- avals*104/(1+avals*hvals*100)
points(rep(110,3),hts,pch=c(16,2,8),xpd=NA)
text(rep(114,3),hts,c("MLE","profile","slice"),xpd=NA,adj=0)


###################################################
### chunk number 53:  eval=FALSE
###################################################
## library(rgl)
## bsurf2 = bsurf
## bsurf2$z = pmin(bsurf2$z,55)
## lev <- -logLik(m2)+qchisq(0.95,c(1,2))/2
## cL = contourLines(bsurf$x,bsurf$y,bsurf$z,
##   levels = lev)
## ll <- -logLik(m2)
## persp3d(bsurf$x,bsurf$y,bsurf2$z,col="gray",alpha=0.8)
## grid3d(side="x-")
## grid3d(side="z-")
## grid3d(side="y-")
## ## add contours?
## lines3d(cL[[1]]$x,cL[[1]]$y,rep(lev[1],length(cL[[1]]$x)),col="red",size=2)
## lines3d(cL[[2]]$x,cL[[2]]$y,rep(lev[2],length(cL[[2]]$x)),col="blue",size=2)
## lines3d(prof@profile$a$par.vals[,"a"],prof@profile$a$par.vals[,"h"],
##         ll+prof@profile$a$z^2/2,size=2,col="green")
## lines3d(prof@profile$h$par.vals[,"a"],prof@profile$h$par.vals[,"h"],
##         ll+prof@profile$h$z^2/2,size=2,col="yellow")
## ## etc.: various bugs to work out


###################################################
### chunk number 54: 
###################################################
source("sqrprofplot.R")
op=par(mfrow=c(1,2),lwd=2,bty="l",las=1,cex=1.5)
par(mar=c(5,4,2,4)+0.1,mgp=c(2,1,0))
sqrprofplot(p2,sqrVal=TRUE,col.minval="gray",col.conf="gray",
     col.prof="black",conf=c(0.95,0.99),
     xlab="Attack rate (a)",
     ylab=expression(paste(Delta,"Negative log-likelihood")),axes=FALSE,
     ylim=c(0,8))
par(bty="u")
box()
axis(side=1)
axis(side=2,at=c(0,4,8),labels=c(0,2,4))
axis(side=4,at=qchisq(c(0.95,0.99),df=1),
     labels=c(expression(frac(chi[1]^2*(0.95),2)),
       expression(frac(chi[1]^{2}*(0.99),2))),cex.axis=0.8)
hts = qchisq(c(0.95,0.99),1)
arrows(conflim95[1,1],hts[1],conflim95[1,2],hts[1],code=3)
text(coef(m2)["a"],hts[1]+0.3,"95%")
arrows(conflim99[1,1],hts[2],conflim99[1,2],hts[2],code=3)
text(coef(m2)["a"],hts[2]+0.3,"99%")
sqrprofplot(p2h,sqrVal=TRUE,col.minval="gray",col.conf="gray",
     col.prof="black",conf=c(0.95,0.99),
     xlab="Handling time (h)",
     ylab=expression(paste(Delta,"Negative log-likelihood")),axes=FALSE,
     ylim=c(0,8))
box()
axis(side=1)
axis(side=2,at=c(0,4,8),labels=c(0,2,4))
axis(side=4,at=qchisq(c(0.95,0.99),df=1),
     labels=c(expression(frac(chi[1]^2*(0.95),2)),
       expression(frac(chi[1]^{2}*(0.99),2))),cex.axis=0.8)
arrows(conflim95[2,1],hts[1],conflim95[2,2],hts[1],code=3)
text(coef(m2)["h"],hts[1]+0.3,"95%")
arrows(conflim99[2,1],hts[2],conflim99[2,2],hts[2],code=3)
text(coef(m2)["h"],hts[2]+0.3,"99%")
par(op)


###################################################
### chunk number 55: 
###################################################
alpha=c(0.95,0.95,0.99,0.99)
cval = qchisq(alpha,1)/2
cval2 = -logLik(m2)+cval
var = rep(c("$a$","$h$"),2)
vals <- rbind(conflim95,conflim99)
rownames(vals) <- NULL
ctab = data.frame(alpha,cval,cval2,var,vals)
ctab[,-4] = signif(ctab[,-4],3)
ctab[2,1:3]=ctab[4,1:3]=NA
colnames(ctab) <- c("$\\alpha$","$\\frac{\\chi_1^2(\\alpha)}{2}$",
                              "$-\\llik+\\frac{\\chi_1^2(\\alpha)}{2}$",
                              "variable","lower","upper")
latex(ctab,file="",table.env=FALSE,title="",rowname=NULL)



###################################################
### chunk number 56: 
###################################################
attach(ReedfrogFuncresp)


###################################################
### chunk number 57: 
###################################################
p2 = profile(m2)
confint(p2)


###################################################
### chunk number 58: 
###################################################
detach(ReedfrogFuncresp)


###################################################
### chunk number 59: 
###################################################
binomNLL2.a = function(p,N,k,a) {
  h = p[1]
  p = a/(1+a*h*N)
  -sum(dbinom(k,prob=p,size=N,log=TRUE))
}


###################################################
### chunk number 60: 
###################################################
avec = seq(0.3,0.8,length=100)
aprof = numeric(100)
for (i in 1:100) {
  aprof[i] = optim(binomNLL2.a,
         par=0.02,k=ReedfrogFuncresp$Killed,N=ReedfrogFuncresp$Initial,
         a = avec[i],
         method="BFGS")$value
}


###################################################
### chunk number 61: 
###################################################
prof.lower = aprof[1:which.min(aprof)]
prof.avec = avec[1:which.min(aprof)]


###################################################
### chunk number 62: 
###################################################
approx(prof.lower,prof.avec,xout=-logLik(m2)+qchisq(0.95,1)/2)


###################################################
### chunk number 63: 
###################################################
x = subset(ReedfrogPred,pred=="pred" & density==10 & size=="small")
k = x$surv
op=par(lwd=2,bty="l",las=1,cex=1.5)
curve(dbeta(x,shape1=sum(k)+1,shape2=40-sum(k)+1),
      xlab="Predation probability\nper capita",
      ylab="Probability density",
      from=0.4,to=1,yaxs="i")
c1 = tcredint("beta",list(shape1=sum(k)+1,shape2=40-sum(k)+1),
  verbose=TRUE)
v = with(as.list(c1),seq(lower,upper,length=100))
w = dbeta(v,shape1=sum(k)+1,shape2=40-sum(k)+1)
polygon(c(v,rev(v)),c(w,rep(0,length(w))),col="gray")
curve(dbeta(x,shape1=sum(k)+1,shape2=40-sum(k)+1),add=TRUE)
abline(h=c1["p"],lty=2)
qs = qbeta(c(0.025,0.975),shape1=sum(k)+1,shape2=40-sum(k)+1)
v2 = seq(0.4,qs[1],length=100)
w2 = dbeta(v2,shape1=sum(k)+1,shape2=40-sum(k)+1)
polygon(c(v2,rev(v2)),c(w2,rep(0,length(w2))),density=10)
v3 = seq(qs[2],1,length=100)
w3 = dbeta(v3,shape1=sum(k)+1,shape2=40-sum(k)+1)
polygon(c(v3,rev(v3)),c(w3,rep(0,length(w3))),density=10)
text(0.75,2.1,"95%\ncredible\ninterval")
text(0.5,1.4,"2.5% tails")
arrows(c(0.5,0.5),c(1.2,1.2),c(0.58,0.88),
       c(0.27,0.30),angle=15)
par(op)
rm(x)
rm(k)


###################################################
### chunk number 64: 
###################################################
post1 = with(frogpred1.bugs$sims.list,kde2d(a,h,n=40))
dx = diff(post1$x[1:2])
dy = diff(post1$y[1:2])
sz = sort(post1$z)
c1 = cumsum(sz)*dx*dy
c95 = approx(c1,sz,xout=0.05)
dens.h = density(frogpred1.bugs$sims.list$h,
  from=post1$y[1],to=post1$y[length(post1$y)],n=length(post1$y))
dens.a = density(frogpred1.bugs$sims.list$a,
  from=post1$x[1],to=post1$x[length(post1$x)],n=length(post1$x))
frog.coda <- lump.mcmc.list(as.mcmc.bugs(frogpred1.bugs))
cred.h = HPDinterval(frog.coda[,"h"])
cred.a = HPDinterval(frog.coda[,"a"])
## plot(sz,c1,type="l")
## abline(v=qval$y)
## abline(h=qval$x)
## sum(sz[sz>qval$y])*dx*dy
##  get contour lines, discard edgy ones
cl1 = contourLines(post1$x,post1$y,post1$z,level=c95)[[5]]
mean.a <- mean(frogpred1.bugs$sims.list$a)
mean.h <- mean(frogpred1.bugs$sims.list$h)
wmode <- which(post1$z==max(post1$z),TRUE)
mode.a <- post1$x[wmode[1]]
mode.h <- post1$y[wmode[2]]


###################################################
### chunk number 65: 
###################################################
op=par(lwd=2,bty="l",las=1,cex=1.5)
nf <- layout(matrix(c(2,1,0,3), 2, 2, byrow=TRUE),
             widths=c(1,3),heights=c(3,1))
par(cex=1.5)
plot(cl1$x,cl1$y,type="l",
     xlab="Attack rate",ylab="",lwd=1,
     xlim=range(post1$x),ylim=range(post1$y))
mtext("Handling time",side=2,at=0.018,line=3,cex=1.5,las=0)
points(mean.a,mean.h,lwd=1)
points(mode.a,mode.h,pch=2,lwd=1)
points(coef(m2)[1],coef(m2)[2],pch=3,lwd=1)
text(c(0.55,0.53,0.49),c(0.0176,0.0149,0.0174),
     c("mean","mode","MLE"),adj=0,cex=0.7)
contour(bsurf$x,bsurf$y,bsurf$z,levels=-logLik(m2)+qchisq(0.95,2)/2,
        add=TRUE,
        drawlabels=FALSE,lty=2)
legend("topright",
       c("bivariate credible region",
         "bivariate confidence region"),
       lty=1:2,
       lwd=2,bty="n")
oldmar = par(mar=c(5.1,0.5,4.1,0.5))$mar
plot(-dens.h$y,dens.h$x,type="l",axes=FALSE,xlab="",ylab="",xaxs="i",
     xlim=c(-81,0))
axis(side=4,labels=FALSE)
axis(side=1,at=c(-80,0),labels=c(80,0))
dh <- subset(data.frame(x=dens.h$x,y=dens.h$y),
             dens.h$x>cred.h[,"lower"] & dens.h$x<cred.h[,"upper"])
polygon(c(-dh$y,rep(0,length(dh$x))),c(dh$x,rev(dh$x)),col="gray")
## have to hand-draw box
u <- par("usr")
segments(c(u[1],u[2]),
         c(u[3],u[3]),
         c(u[2],u[2]),
         c(u[3],u[4]),lwd=2)
par(mar=c(2,4,0,1)+0.1)
plot(dens.a$x,dens.a$y,type="l",xlab="",ylab="",axes=FALSE,
     ylim=c(0,6),yaxs="i",lwd=2)
axis(side=1)
axis(side=2,at=c(0,6))
box()
da <- subset(data.frame(x=dens.a$x,y=dens.a$y),
             dens.a$x>cred.a[,"lower"] & dens.a$x<cred.a[,"upper"])
polygon(c(da$x,rev(da$x)),c(da$y,rep(0,length(da$x))),col="gray")
par(op)


###################################################
### chunk number 66: 
###################################################
op=par(lwd=2,bty="l",las=1,cex=1.5)
bsurf3 =  curve3d(binomNLL2(c(x,y),N=ReedfrogFuncresp$Initial,k=ReedfrogFuncresp$Killed),
  sys="contour",from=c(0.3,0.001),to=c(0.82,0.04),n=c(51,51),
  levels=-logLik(m2)+qchisq(c(0.8,0.995),2)/2,
  drawlabels=FALSE,
  xlab="Attack rate (a)",ylab="Handling time (h)")
library(ellipse)
lines(ellipse(vcov(m2),centre=coef(m2),level=0.8),col="gray")
lines(ellipse(vcov(m2),centre=coef(m2),level=0.995),col="gray")
legend("topleft",c("profile","information"),lty=1,col=c("black","gray"),
       bty="n")
text(c(0.65,0.72),c(0.027,0.036),c("80%","99.5%"))
par(op)


###################################################
### chunk number 67:  eval=FALSE
###################################################
## m = matrix(c(0.5,1,1,2),nrow=2)
## banana <- function(x,y,rot=FALSE,theta=pi/4){
##   if (rot) {
##     newx = cos(theta)*x+sin(theta)*y
##     newy = -sin(theta)*x+cos(theta)*y
##     x=newx
##     y=newy
##   }
##   100*(y-x^2)^2+(1-x)^2
## } 
## m1 = mle2(banana,start=list(x=1,y=1),data=list(rot=TRUE))
## cov2cor(vcov(m1))
## curve3d(banana(x,y,rot=TRUE),from=c(-1,-1),to=c(2,2),sys3d="contour",
##         levels=0:10)
## points(0,sqrt(2))


###################################################
### chunk number 68: 
###################################################
gourd <- function(x,y,s1=0,s2=3,a=1,s1exp=1,rot=FALSE,theta=pi/4){
  if (rot) {
    newx = cos(theta)*x+sin(theta)*y
    newy = -sin(theta)*x+cos(theta)*y
    x=newx
    y=newy
  }
  a*((exp(s1+s1exp*y)*x)^2+(s2*y)^2)
} 
fun1 = function(x,y) gourd(x,y,s1exp=0,a=0.5)
fun2 = function(x,y) gourd(x,y,s1exp=0,a=0.5,rot=TRUE)
fun3 = function(x,y) gourd(x,y,s1exp=2,a=1)
fun4 = function(x,y) gourd(x,y,s1exp=2,a=1.2,rot=TRUE)
c1 = curve3d(fun1,from=c(-3,-2),to=c(3,2),
          sys3d="none")
c2 = curve3d(fun2,from=c(-3,-2),to=c(3,2),
          sys3d="none")
c3 = curve3d(fun3,from=c(-3,-2),to=c(3,2),
          sys3d="none",n=c(71,71))
c4 = curve3d(fun4,from=c(-3,-2),to=c(3,2),
          sys3d="none",n=c(91,91))
tmpf <- function(fun0,c0,heights=c(-2.3,-1.7,-2),
                 xpos = 2.35,
                 quadlty=1,quadcol="darkgray",
                 legend=FALSE) {
  ## 2D (univariate) confidence region
  contour(c0$x,c0$y,c0$z,
          levels=qchisq(0.95,1)/2,
          drawlabels=FALSE,axes=FALSE,
          xlab="",ylab="",
          xlim=c(-3,3.4),
          ylim=c(-2.4,1.8))
  box(col="gray")
  m = mle2(fun0,start=list(x=0,y=0))
  p = profile(m)
  s = slice(m)
  sliceconf.x = approx(s@profile$x$z,
    s@profile$x$par.vals[,"x"],
    xout = qnorm(c(0.025,0.975)))$y
  profconf.x = confint(p)["x",]
  profconf.y = approx(p@profile$x$z,
    y=p@profile$x$par.vals[,"y"],
    xout = qnorm(c(0.025,0.975)))$y
  ellconf.x = confint(m,method="quad")["x",]
  v = vcov(m)
  slope = v[1,2]/v[1,1]
  ellconf.y = ellconf.x*slope
  lines(ellipse(vcov(m),centre=coef(m),
                t=sqrt(qchisq(0.95,1))),lty=quadlty,
        col=quadcol)
  ## redraw
  contour(c0$x,c0$y,c0$z,
          levels=qchisq(0.95,1)/2,
          drawlabels=FALSE,add=TRUE)
  lines(p@profile$x$par.vals[,"x"],p@profile$x$par.vals[,"y"],lty=2)
  abline(a=0,b=slope)
  points(ellconf.x,ellconf.y,pch=2)
  points(profconf.x,profconf.y,pch=5)
  points(sliceconf.x,rep(0,2),pch=8)
  points(ellconf.x,rep(heights[1],2),pch=2)
  segments(ellconf.x[1],heights[1],ellconf.x[2],heights[1],pch=2)
  points(profconf.x,rep(heights[2],2),pch=5)
  segments(profconf.x[1],heights[2],profconf.x[2],heights[2])
  points(sliceconf.x,rep(heights[3],2),pch=8)
  segments(sliceconf.x[1],heights[3],sliceconf.x[2],heights[3])
  text(rep(xpos,3),heights,c("quad","profile","slice"),cex=0.75,adj=0)
  if (legend) {
    legend("topleft",
           c("conf. region",
             "quadratic",
             "profile"),
           lwd=2,
           lty=c(1,quadlty,2),
           col=c(1,quadcol,1),
           bty="n",cex=0.75)
  }
}


###################################################
### chunk number 69: 
###################################################
op=par(mfrow=c(2,2),lwd=2,bty="l",las=1,cex=1.5,
  mar=c(0,0,0,0))
tmpf(fun1,c1)
tmpf(fun2,c2)
tmpf(fun3,c3,legend=TRUE)
tmpf(fun4,c4)
par(op)


###################################################
### chunk number 70:  eval=FALSE
###################################################
## ## a = initial slope
## ## h = 1/asymptote
## ## half-maximum=
## binomNLL3 = function(p,N,k) {
##   a = p[1]
##   h = p[2]*p[1]
##   p = a/(1+a*h*N)
##   -sum(dbinom(k,prob=p,size=N,log=TRUE))
## }
## parnames(binomNLL3) <- c("a","ha")
## m5 = mle2(binomNLL3,start=c(a=0.5,ha=0.01),data=list(N=ReedfrogFuncresp$Initial,k=ReedfrogFuncresp$Killed))
## curve3d(binomNLL3(c(x,y),N=ReedfrogFuncresp$Initial,k=ReedfrogFuncresp$Killed),
##         from=c(0.3,0.001),to=c(0.7,0.05),sys3d="contour",levels=40:60)


###################################################
### chunk number 71: 
###################################################
data(FirDBHFec)
X = na.omit(FirDBHFec[,c("TOTCONES","DBH","WAVE_NON")])
X$TOTCONES = round(X$TOTCONES)
attach(X)


###################################################
### chunk number 72: 
###################################################
## all kinds of computation so we can draw the figure up front
nbNLL.ab = function(a.w,b.w,a.n,b.n,k) {
  wcode = as.numeric(WAVE_NON)
  a = c(a.n,a.w)[wcode]
  b  = c(b.n,b.w)[wcode]
  predcones = a*DBH^b
  -sum(dnbinom(TOTCONES,mu=predcones,size=k,log=TRUE))
}
nbNLL.0 = function(a,b,k) {
  predcones = a*DBH^b
  -sum(dnbinom(TOTCONES,mu=predcones,size=k,log=TRUE))
}
nbNLL.a = function(a.n,a.w,b,k) {
  wcode = as.numeric(WAVE_NON)
  a = c(a.n,a.w)[wcode]
  b  = b
  predcones = a*DBH^b
  -sum(dnbinom(TOTCONES,mu=predcones,size=k,log=TRUE))
}
nbNLL.b = function(a,b.n,b.w,k) {
  wcode = as.numeric(WAVE_NON)
  b  = c(b.n,b.w)[wcode]
  predcones = a*DBH^b
  -sum(dnbinom(TOTCONES,mu=predcones,size=k,log=TRUE))
}
nbNLL.abk = function(a.n,a.w,b.n,b.w,k.n,k.w) {
  wcode = as.numeric(WAVE_NON)
  a = c(a.n,a.w)[wcode]
  b  = c(b.n,b.w)[wcode]
  k = c(k.n,k.w)[wcode]
  predcones = a*DBH^b
  -sum(dnbinom(TOTCONES,mu=predcones,size=k,log=TRUE))
}
nbfit.0 = mle2(nbNLL.0,start=list(a=1,b=1,k=1))
a = coef(nbfit.0)["a"]
b = coef(nbfit.0)["b"]
k = coef(nbfit.0)["k"]
nbfit.ab = mle2(nbNLL.ab,
  start=list(a.n=a,a.w=a,b.n=b,b.w=b,k=k))
nbfit.a = mle2(nbNLL.a,
  start=list(a.n=a,a.w=a,b=b,k=k))
nbfit.abk = mle2(nbNLL.abk,
  start=list(a.n=a,a.w=a,b.n=b,b.w=b,
    k.n=k,k.w=k))
nbNLL.ak = function(a.n,a.w,b,k.n,k.w) {
  wcode = as.numeric(WAVE_NON)
  a = c(a.n,a.w)[wcode]
  k = c(k.n,k.w)[wcode]
  predcones = a*DBH^b
  -sum(dnbinom(TOTCONES,mu=predcones,size=k,log=TRUE))
}
nbfit.ak = mle2(nbNLL.ak,
  start=list(a.n=a,a.w=a,b=b,k.n=k,k.w=k))
nbNLL.bk = function(a,b.n,b.w,k.n,k.w) {
  wcode = as.numeric(WAVE_NON)
  b  = c(b.n,b.w)[wcode]
  k = c(k.n,k.w)[wcode]
  predcones = a*DBH^b
  -sum(dnbinom(TOTCONES,mu=predcones,size=k,log=TRUE))
}
nbfit.bk = mle2(nbNLL.bk,
  start=list(a=a,b.n=b,b.w=b,k.n=k,k.w=k))
nbNLL.k = function(a,b,k.n,k.w) {
  wcode = as.numeric(WAVE_NON)
  k = c(k.n,k.w)[wcode]
  predcones = a*DBH^b
  -sum(dnbinom(TOTCONES,mu=predcones,size=k,log=TRUE))
}
nbfit.k = mle2(nbNLL.k,
  start=list(a=a,b=b,k.n=k,k.w=k))
nbNLL.b = function(a,b.n,b.w,k) {
  wcode = as.numeric(WAVE_NON)
  b  = c(b.n,b.w)[wcode]
  predcones = a*DBH^b
  -sum(dnbinom(TOTCONES,mu=predcones,size=k,log=TRUE))
}
nbfit.b = mle2(nbNLL.b,
  start=list(a=a,b.n=b,b.w=b,k=k))


###################################################
### chunk number 73: 
###################################################
op=par(lwd=2,bty="l",las=1,cex=1.5)
plot(TOTCONES~DBH,pch=as.numeric(WAVE_NON),cex=0.8,
     xlab="Size (DBH)",ylab="Fecundity (total cones)")
par(xpd=NA)
legend(c(4,10),
       c(250,350),pch=c(1:2,NA),
       lty=c(2,3,1),
       c("nonwave","wave","combined"),bty="n",merge=FALSE)
par(xpd=FALSE)
with(as.list(coef(nbfit.ab)), {
  curve(a.n*x^b.n,add=TRUE,lty=2)
  curve(a.w*x^b.w,add=TRUE,lty=3)
})
with(as.list(coef(nbfit.0)), {
  curve(a*x^b,add=TRUE,lty=1)
})
par(op)


###################################################
### chunk number 74: 
###################################################
detach(X)


###################################################
### chunk number 75: 
###################################################
tmpf = function(x,a=0.4,b=0.1,c=2,d=1) {
  (a+b*x+c*x^2)*exp(-d*x)
}
dpars = formals(tmpf)[-1]
set.seed(1005)
npts = 10
x = runif(npts,min=1,max=7)
y_det = tmpf(x)
y = y_det+rnorm(npts,sd=0.35)
y2 = y_det+rnorm(npts,sd=0.35)
ymax=2
n1 = nls(y~tmpf(x,a,b,c,d),start=list(a=0.4,b=0.1,c=2,d=1))
n2 = nls(y~a*x*exp(-b*x),start=list(a=1,b=0.5))
p0 = rep(mean(y),length(y))
p1 = predict(n1)
xvec = seq(0,7,length=150)
p1vec = predict(n1,newdata=list(x=xvec))
p2vec = predict(n2,newdata=list(x=xvec))
p2 = predict(n2)
calc_r2 = function(y) {
  s0 = sum((y-p0)^2)
  s1 = sum((y-p1)^2)
  s2 = sum((y-p2)^2)
  c(1-s1/s0,1-s2/s0,s0,s1,s2,which.min(c(s0,s1,s2)))
}
r2.0 = calc_r2(y)
r2.1 = calc_r2(y2)
r2vec = t(replicate(500,calc_r2(y_det+rnorm(npts,sd=0.35))))
r2vec_mean = colMeans(r2vec)[1:5]
tv = table(r2vec[,6])


###################################################
### chunk number 76: 
###################################################
op = par(mfrow=c(1,2))
par(lwd=2,bty="l",las=1,cex=1.5, mgp=c(2.5,1,0))
par(mar=c(5,4,2,0.5)+0.1)
plot(x,y,ylim=c(0,max(c(y,ymax))),xlim=c(0,7),axes=FALSE,
      xlab="",ylab="")
points(x[7],y[7],pch=16)
abline(h=p0)
axis(side=1,at=c(0,3,6))
axis(side=2)
box()
tcol = "gray"
curve(tmpf,add=TRUE,from=0,col=tcol)
lines(xvec,p1vec,lty=2)
lines(xvec,p2vec,lty=3)
par(xpd=NA)
legend(c(3.1,8.6),c(1.2,2),
       c("constant",
         "Ricker",
         "gen Ricker",
         "true"),
       col=rep(c("black",tcol),c(3,1)),
       lty=c(1,3,2,1),
       bty="n",
       cex=0.75)     
par(xpd=FALSE)
par(mar=c(5,1,2,3.5)+0.1)
plot(x,y2,ylim=c(0,max(c(y,ymax))),xlim=c(0,7),axes=FALSE,
     xlab="",ylab="")
points(x[7],y2[7],pch=16)
curve(tmpf,add=TRUE,from=0,col=tcol)
abline(h=p0)
xvec = seq(0,7,length=150)
lines(xvec,predict(n1,newdata=list(x=xvec)),lty=2)
lines(xvec,predict(n2,newdata=list(x=xvec)),lty=3)
axis(side=1,at=c(0,3,6))
axis(side=2,labels=FALSE)
box()


###################################################
### chunk number 77: 
###################################################
data(FirDBHFec)
X = na.omit(FirDBHFec[,c("TOTCONES","DBH","WAVE_NON")])
X$TOTCONES = round(X$TOTCONES)


###################################################
### chunk number 78: 
###################################################
nbfit.0 = mle2(TOTCONES~dnbinom(mu=a*DBH^b,size=k),
  start=list(a=1,b=1,k=1),data=X)


###################################################
### chunk number 79: 
###################################################
start.ab = as.list(coef(nbfit.0))
nbfit.ab = mle2(TOTCONES~dnbinom(mu=a*DBH^b,size=k),
  start=start.ab, data=X,
  parameters=list(a~WAVE_NON,b~WAVE_NON))


###################################################
### chunk number 80:  eval=FALSE
###################################################
## start.ab3=list(a.WAVE_NONn=0.3,
##     a.WAVE_NONw=0.3,
##     b.WAVE_NONn=2.3,
##     b.WAVE_NONw=2.3,k=1.5)
## nbfit.ab2 = mle2(TOTCONES~dnbinom(mu=a*DBH^b,size=k),
##   start=start.ab, data=X,
##   parameters=list(a~WAVE_NON-1,b~WAVE_NON-1))
## f2 = calc_mle2_function(TOTCONES~dnbinom(mu=a*DBH^b,size=k),
##                    start=start.ab,data=X,
##                    parameters=list(a~WAVE_NON-1,b~WAVE_NON-1))


###################################################
### chunk number 81: 
###################################################
attach(X)
nbNLL.ab = function(a.w,b.w,a.n,b.n,k) {
  wcode = as.numeric(WAVE_NON)
  a = c(a.n,a.w)[wcode]
  b  = c(b.n,b.w)[wcode]
  predcones = a*DBH^b
  -sum(dnbinom(TOTCONES,mu=predcones,size=k,log=TRUE))
}


###################################################
### chunk number 82:  eval=FALSE
###################################################
## k=c(k.n,k.w)[wcode]


###################################################
### chunk number 83: 
###################################################
nbNLL.a = function(a.n,a.w,b,k) {
  wcode = as.numeric(WAVE_NON)
  a = c(a.n,a.w)[wcode]
  b  = b
  predcones = a*DBH^b
  -sum(dnbinom(TOTCONES,mu=predcones,size=k,log=TRUE))
}
nbNLL.b = function(a,b.n,b.w,k) {
  wcode = as.numeric(WAVE_NON)
  b  = c(b.n,b.w)[wcode]
  predcones = a*DBH^b
  -sum(dnbinom(TOTCONES,mu=predcones,size=k,log=TRUE))
}
nbNLL.abk = function(a.n,a.w,b.n,b.w,k.n,k.w) {
  wcode = as.numeric(WAVE_NON)
  a = c(a.n,a.w)[wcode]
  b  = c(b.n,b.w)[wcode]
  k = c(k.n,k.w)[wcode]
  predcones = a*DBH^b
  -sum(dnbinom(TOTCONES,mu=predcones,size=k,log=TRUE))
}


###################################################
### chunk number 84: 
###################################################
poisNLL.0 <- function(a,b) {
  predcones = a*DBH^b
  -sum(dpois(TOTCONES,lambda=predcones,log=TRUE))
}
poisfit.0 = mle2(poisNLL.0,start=list(a=1,b=1))


###################################################
### chunk number 85: 
###################################################
nbfit.a = mle2(nbNLL.a,
  start=list(a.n=a,a.w=a,b=b,k=k))
nbfit.abk = mle2(nbNLL.abk,
  start=list(a.n=a,a.w=a,b.n=b,b.w=b,
    k.n=k,k.w=k))


###################################################
### chunk number 86: 
###################################################
## redo all with formula interface
nbfit.0 = mle2(TOTCONES~dnbinom(mu=a*DBH^b,size=k),
  start=list(a=1,b=1,k=1),data=X)
start.ab = as.list(coef(nbfit.0))
nbfit.a = mle2(TOTCONES~dnbinom(mu=a*DBH^b,size=k),
  start=start.ab, data=X,
  parameters=list(a~WAVE_NON))
nbfit.b = mle2(TOTCONES~dnbinom(mu=a*DBH^b,size=k),
  start=start.ab, data=X,
  parameters=list(b~WAVE_NON))
nbfit.abk = mle2(TOTCONES~dnbinom(mu=a*DBH^b,size=k),
  start=start.ab, data=X,
  parameters=list(a~WAVE_NON-1,b~WAVE_NON-1,k~WAVE_NON-1))
nbfit.ak = mle2(TOTCONES~dnbinom(mu=a*DBH^b,size=k),
  start=start.ab, data=X,
  parameters=list(a~WAVE_NON,k~WAVE_NON))
nbfit.bk = mle2(TOTCONES~dnbinom(mu=a*DBH^b,size=k),
  start=start.ab, data=X,
  parameters=list(a~WAVE_NON,k~WAVE_NON))


###################################################
### chunk number 87: 
###################################################
anova(nbfit.0,nbfit.a,nbfit.ab)


###################################################
### chunk number 88:  eval=FALSE
###################################################
## nbNLL.ak = function(a.n,a.w,b,k.n,k.w) {
##   wcode = as.numeric(WAVE_NON)
##   a = c(a.n,a.w)[wcode]
##   k = c(k.n,k.w)[wcode]
##   predcones = a*DBH^b
##   -sum(dnbinom(TOTCONES,mu=predcones,size=k,log=TRUE))
## }
## nbfit.ak = mle2(nbNLL.ak,
##   start=list(a.n=a,a.w=a,b=b,k.n=k,k.w=k))
## nbNLL.bk = function(a,b.n,b.w,k.n,k.w) {
##   wcode = as.numeric(WAVE_NON)
##   b  = c(b.n,b.w)[wcode]
##   k = c(k.n,k.w)[wcode]
##   predcones = a*DBH^b
##   -sum(dnbinom(TOTCONES,mu=predcones,size=k,log=TRUE))
## }
## nbfit.bk = mle2(nbNLL.bk,
##   start=list(a=a,b.n=b,b.w=b,k.n=k,k.w=k))
## nbNLL.k = function(a,b,k.n,k.w) {
##   wcode = as.numeric(WAVE_NON)
##   k = c(k.n,k.w)[wcode]
##   predcones = a*DBH^b
##   -sum(dnbinom(TOTCONES,mu=predcones,size=k,log=TRUE))
## }
## nbfit.k = mle2(nbNLL.k,
##   start=list(a=a,b=b,k.n=k,k.w=k))
## nbNLL.b = function(a,b.n,b.w,k) {
##   wcode = as.numeric(WAVE_NON)
##   b  = c(b.n,b.w)[wcode]
##   predcones = a*DBH^b
##   -sum(dnbinom(TOTCONES,mu=predcones,size=k,log=TRUE))
## }
## nbfit.b = mle2(nbNLL.b,
##   start=list(a=a,b.n=b,b.w=b,k=k))


###################################################
### chunk number 89: 
###################################################
detach(X)


###################################################
### chunk number 90: 
###################################################
## model nesting plot
## 6 parameters: a.{w,n},b.{w,n},k.{w,n}
## 5 parameters -- 3 possibilities
## 4 parameters -- 3 possibilities
## 3 parameters -- 1 possibility
##
## points; segments; text indicating parameters; LRT p-values
nlevel = 4
ylocs = rev((1:nlevel)/(nlevel+1))
xlocs = lapply(c(1,3,3,1),function(x) (1:x)/(x+1)-0.1)
txt = list(
  list(
       list("all parameters","equal")),
  list(
       list(expression(a[w]!=a[n])),
       list(expression(b[w]!=b[n])),
       list(expression(k[w]!=k[n]))),
  list(
       list(expression(a[w]!=a[n]),
            expression(b[w]!=b[n])),
       list(expression(a[w]!=a[n]),
            expression(k[w]!=k[n])),
       list(expression(b[w]!=b[n]),
            expression(k[w]!=k[n]))),
       list(
            list("all parameters","different")))
models = list(
  list(nbfit.0),
  list(nbfit.a,nbfit.b,nbfit.k),
  list(nbfit.ab,nbfit.ak,nbfit.bk),
  list(nbfit.abk))


###################################################
### chunk number 91: 
###################################################
#,width=8,height=8>>=
op=par(lwd=2,bty="l",las=1,cex=1.5,mar=c(0,0,0,0))
spc = 0.03 ## y-spacing: was 0.025
ht=0.014  ## rectangle height? was 0.012
wid=0.08
cex1 = 0.7
plot(0:1,0:1,axes=FALSE,type="n",ann=FALSE)
for (i in 1:nlevel) {
  ## draw segments
  p = 0.05
  if (i<nlevel) {
    for (j in 1:length(txt[[i]])) {
      for (k in 1:length(txt[[i+1]])) {
        x1=xlocs[[i]][[j]]
        x2=xlocs[[i+1]][[k]]
        y1=ylocs[i]
        y2=ylocs[i+1]
        segments(x1*(1-p)+x2*p,
                 y1*(1-p)+y2*p,
                 x1*p+x2*(1-p),
                 y1*p+y2*(1-p),lwd=3)
      }
    }
  }   
}
for (i in 1:nlevel) {
  for (j in 1:length(txt[[i]])) {
    yoff = seq(0,by=-spc,length=length(txt[[i]][[j]])+1)
    for (k in 1:length(txt[[i]][[j]])) {
      rect(xlocs[[i]][j]-wid,
           ylocs[i]+yoff[k]-ht,
           xlocs[[i]][j]+wid,
           ylocs[i]+yoff[k]+ht,
           col="white",border="white")
      text(xlocs[[i]][j],ylocs[i]+yoff[k],txt[[i]][[j]][[k]],
           cex=cex1)
    }
    rect(xlocs[[i]][j]-wid,
         ylocs[i]+yoff[length(yoff)]-ht,
         xlocs[[i]][j]+wid,
         ylocs[i]+yoff[length(yoff)]+ht,
         col="white",border="white")
    text(xlocs[[i]][j],ylocs[i]+yoff[length(yoff)],
         paste("D=",
               ## hack!!
               ## can't find another way to get this consistent with the ANOVA table in text
               if (i==3 && j==1) 2271.3 else round(-2*logLik(models[[i]][[j]]),1),
##               signif(-2*logLik(models[[i]][[j]]),5),
##               round(-2*logLik(models[[i]][[j]]),1),
               sep=""),
         cex=cex1)
  }
}
text(1,ylocs,adj=1,
     labels=paste(3:6,"parameters"),cex=0.7)


###################################################
### chunk number 92: 
###################################################
poisfit.ab = mle2(TOTCONES~dpois(a*DBH^b),
  start=list(a=1,b=1), data=X,
  parameters=list(a~WAVE_NON,b~WAVE_NON))
anova(poisfit.ab,nbfit.ab)


###################################################
### chunk number 93: 
###################################################
a1 <- AICtab(nbfit.0,nbfit.a,nbfit.b,
             nbfit.k,
             nbfit.ab,
             nbfit.ak,nbfit.bk,nbfit.abk,delta=TRUE)
a2 <- AICctab(nbfit.0,nbfit.a,nbfit.b,
              nbfit.k,
              nbfit.ab,nbfit.ak,nbfit.bk,nbfit.abk,delta=TRUE,
          nobs=nrow(X))
b1 <- BICtab(nbfit.0,nbfit.a,nbfit.b,
             nbfit.k,
             nbfit.ab,
    nbfit.ak,nbfit.bk,nbfit.abk,delta=TRUE,nobs=nrow(X))
atab <- cbind(a1$df,a1$dAIC,a2$dAICc,b1$dBIC)
rownames(atab) <- attr(a1,"row.names")
## need to double backslashes
colnames(atab) <- c("params","$\\Delta \\mbox{AIC}$","$\\Delta \\mbox{AIC}_c$","$\\Delta \\mbox{BIC}$")


###################################################
### chunk number 94: 
###################################################
latex(round(atab,2),file="",table.env=FALSE,title="model")


###################################################
### chunk number 95: 
###################################################
lnormfit.ab = mle2(TOTCONES+0.001~dlnorm(meanlog=b*log(DBH)+log(a),
  sdlog=sdlog),
  start=list(a=1,b=1,sdlog=0.1),
  data=X,
  parameters=list(a~WAVE_NON,b~WAVE_NON),method="Nelder-Mead")
gammafit.ab = mle2(TOTCONES+0.001~dgamma(scale=a*DBH^b/shape,
  shape=shape),
  start=list(a=1,b=1,shape=2),
  data=X,
  parameters=list(a~WAVE_NON,b~WAVE_NON))


###################################################
### chunk number 96: 
###################################################
s0tab = AICtab(poisfit.ab,gammafit.ab,lnormfit.ab,nbfit.ab,delta=TRUE,sort=TRUE)
stab <- cbind(s0tab$AIC,s0tab$df,s0tab$dAIC)
stab=round(stab,1)
rownames(stab) <- attr(s0tab,"row.names")
## need to double backslashes
colnames(stab) = c("AIC","df","$\\Delta\\mbox{AIC}$")
latex(stab,file="",table.env=FALSE,title="",rowname=c("Neg. binom.","Gamma","Log-normal","Poisson"))


###################################################
### chunk number 97:  eval=FALSE
###################################################
## ## checking Poisson-gamma parameterization of neg binomial
## x = rnbinom(100,mu=1,size=0.1)
## n = length(x)
## inits = list(list(mu=2,k=1))
## nbinom.bugs <- bugs(data=list("x","n"),
##                     inits,parameters.to.save=c("mu","k"),
##                     model.file="negbinom.bug",
##                     n.chains=1,n.iter=2000)


###################################################
### chunk number 98: 
###################################################
DBH <- X$DBH
cones <- X$TOTCONES
n <- length(DBH)



###################################################
### chunk number 99:  eval=FALSE
###################################################
## inits <- list(list(a=0,b=0,k=1),list(a=2,b=-2,k=1),
##               list(a=-2,b=2,k=1))
## firfec0.bugs <- bugs(data=list("DBH","cones","n"),
##                        inits,parameters.to.save=c("a","b","k"),
##                        model.file="firfec0.bug",
##                        n.chains=length(inits),n.iter=9000)
## inits <- list(list(a=c(0,0),b=c(0,0),k=1),
##               list(a=c(2,0),b=c(-2,0),k=1),
##               list(a=c(0,2),b=c(0,-2),k=1))
## grp <- as.numeric(X$WAVE_NON)
## firfec1.bugs <- bugs(data=list("DBH","cones","n","grp"),
##                        inits,parameters.to.save=c("a","b","k"),
##                        model.file="firfec1.bug",
##                        n.chains=length(inits),n.iter=6000)
## inits <- list(list(a=c(0,0),b=c(0,0),k=c(1,1)),
##               list(a=c(2,0),b=c(-2,0),k=c(1,1)),
##               list(a=c(0,2),b=c(0,-2),k=c(1,1)))
## firfec2.bugs <- bugs(data=list("DBH","cones","n","grp"),
##                        inits,parameters.to.save=c("a","b","k"),
##                        model.file="firfec2.bug",
##                        n.chains=length(inits),n.iter=6000)


###################################################
### chunk number 100: 
###################################################
load("firfec-batch.RData")
## load("firfec-batch2.RData")
## load("firfec-batch3.RData")
rm(DBH)


###################################################
### chunk number 101: 
###################################################
flist = list(firfec0.bugs,firfec1.bugs,firfec2.bugs)
x0 = firfec0.bugs$sims.list$deviance
x1 = firfec1.bugs$sims.list$deviance
x2 = firfec2.bugs$sims.list$deviance
mx = min(x2)
## calc (offset) likelihood
x0A = exp(-(x0-mx)/2)
x1A = exp(-(x1-mx)/2)
x2A = exp(-(x2-mx)/2)
## harmonic mean estimator
x0B = 1/mean(1/x0A)
x1B = 1/mean(1/x1A)
x2B = 1/mean(1/x2A)
## plot(density(log(x0A),to=0),xlim=c(-90,0),ylim=c(0,0.04))
## lines(density(log(x1A),to=0),col=2)
## lines(density(log(x2A),to=0),col=4)
postprob = c(x0B,x1B,x2B)/(x0B+x1B+x2B)
## bfac = c(f01=x0B/x1B,f02=x0B/x2B,f12=x1B/x2B)
bdiff <- -2*(log(postprob)-max(log(postprob)))


###################################################
### chunk number 102: 
###################################################
## postprob
## AIC(nbfit.0,nbfit.ab,nbfit.abk,delta=TRUE,
##    weights=TRUE)
b1 = BICtab(nbfit.0,nbfit.ab,
  nbfit.abk,delta=TRUE,nobs=nrow(X))
## BIC(nbfit.0,nbfit.ab,
##    nbfit.abk,nobs=nrow(X))


###################################################
### chunk number 103: 
###################################################
lprior.0 <- function(a,b,k) {
  sum(dnorm(c(a,b),sd=10,log=TRUE))+dunif(k,0.1,5,log=TRUE)
}
lprior.ab <- function(a.w,b.w,a.n,b.n,k) {
  sum(dnorm(c(a.w,a.n,b.w,b.n),sd=10,log=TRUE))+
    dunif(k,0.1,5,log=TRUE)
}
lprior.abk <- function(a.w,b.w,a.n,b.n,k.w,k.n) {
  sum(dnorm(c(a.w,a.n,b.w,b.n),sd=10,log=TRUE))+
    sum(dunif(c(k.w,k.n),0.1,5,log=TRUE))
}
if (FALSE) {
  ## need to rebuild mle/check/tweak
  b0 <- bayesfactor(nbfit.0,log=TRUE,logprior=lprior.0)
  b1 <- bayesfactor(nbfit.ab,log=TRUE,logprior=lprior.ab)
  b2 <- bayesfactor(nbfit.abk,log=TRUE,logprior=lprior.abk)
  bvec <- 2*c(b0,b1,b2)
  bvec <- bvec-min(bvec)
}


###################################################
### chunk number 104: 
###################################################
attach(X)
nblpost.0 <- function(a,b,k) {
  nbNLL.0(a,b,k)-lprior.0(a,b,k)
}
nblpost.ab <- function(a.n,a.w,b.n,b.w,k) {
  nbNLL.ab(a.w,a.n,b.w,b.n,k)-lprior.ab(a.w,a.n,b.w,b.n,k)
}
nblpost.abk<- function(a.n,a.w,b.n,b.w,k.n,k.w) {
  nbNLL.abk(a.w,a.n,b.w,b.n,k.w,k.n)-lprior.abk(a.w,a.n,b.w,b.n,k.w,k.n)
}
## redo fit with -1
nbfit.ab = mle2(TOTCONES~dnbinom(mu=a*DBH^b,size=k),
 start=start.ab, data=X,
 parameters=list(a~WAVE_NON-1,b~WAVE_NON-1))
start2 = as.list(coef(nbfit.ab))
names(start2) <- names(formals(nblpost.ab)) 
start3 = as.list(coef(nbfit.abk))
names(start3) <- names(formals(nblpost.abk)) 
p0 <- mle2(minuslogl=nblpost.0,start=list(a=1,b=1,k=1))
p.ab <- mle2(minuslogl=nblpost.ab,start=start2)
p.abk <- mle2(minuslogl=nblpost.abk,start=start3)
detach(X)
marglik <- function(obj,method="laplace",log=FALSE) {
  v <- vcov(obj)
  d <- nrow(v)
  if (method=="laplace") {
    logdet = c(determinant(v,logarithm=TRUE)$modulus)
    r <- d/2*log(2*pi) + 1/2*logdet+c(logLik(obj))
  } else if (method=="adapt") {
    require(adapt)
    c1 <- coef(obj)
    sd <- sqrt(diag(v))
    lower <- c1-3*sd
    upper <- c1+3*sd
    minval <- 0
    fn <- function(p,log=FALSE) {
      par <- as.list(p)
      names(par) <- names(c1)
      m <- obj@minuslogl
      x <- do.call("m",par)-minval
      if (log) x else exp(x)
    }
    adapt(d,lower,upper,functn=fn)
  }
  if (log) r else exp(r)
}
mlik2 <- sapply(list(p0,p.ab,p.abk),marglik,log=TRUE)
mdiff <- -mlik2-min(-mlik2)


###################################################
### chunk number 105: 
###################################################
apptab = round(cbind(bdiff,mdiff,b1$dBIC),1)
dimnames(apptab) = list(c("null","$a$, $b$ differ","$a$, $b$, $k$ differ"),
          c("harmonic mean","Laplace","BIC"))
latex(apptab,file="",title="",table.env=FALSE)


###################################################
### chunk number 106: 
###################################################
DICvec <- sapply(flist,"[[","DIC")
pDvec <- sapply(flist,"[[","pD")


###################################################
### chunk number 107:  eval=FALSE
###################################################
## op=par(lwd=2,bty="l",las=1,cex=1.5,mgp=c(2.3,1,0))
## curve(qchisq(0.95,x)/2,from=1,to=10,
##       xlab="Extra parameters",
##       ylab="Critical log-likelihood\ndifference")
## curve(1*x,add=TRUE,lty=2)
## curve(x*log(10)/2,add=TRUE,lty=3)
## curve(x*log(100)/2,add=TRUE,lty=4)
## legend("bottomright",
##        lty=1:4,
##        c("LRT","AIC","BIC, nobs=10",
##          "BIC, nobs=100"))
## par(op)


###################################################
### chunk number 108: 
###################################################
op=par(lwd=2,bty="l",las=1,cex=1.5)
predfun = function(a.n,a.w,b.n,b.w,k.n,k.w) {
  wcode = as.numeric(X$WAVE_NON)
  a = c(a.n,a.w)[wcode]
  sl  = c(b.n,b.w)[wcode]
  k = c(k.n,k.w)[wcode]
  a*X$DBH^b
}
pred <- do.call("predfun",as.list(coef(p.abk)))
par(pty="s")
plot(pred,X$TOTCONES,xlab="Predicted cones",ylab="Actual cones",
         xlim=range(pred),
         ylim=range(pred),
         log="xy",
         pch=as.numeric(X$WAVE_NON))
par(pty="m")
abline(a=0,b=1)
legend("topleft",c("nonwave","wave"),
       pch=1:2,bty="n")
par(op)


###################################################
### chunk number 109: 
###################################################
dbh <- equal.count(X$DBH,number=4)
c1 <- coef(p.abk)
trellis.par.set(canonical.theme(color=FALSE))
print(densityplot(~TOTCONES|dbh*WAVE_NON,data=X,from=0,
                  ylim=c(0.0,0.04),
                  xlim=c(0,350),
            panel=function(x,...) {
              panel.densityplot(x,...)
              nonwave = (current.row()==1)
              elevel = current.column()
              meandbh <- mean(unlist(as.matrix(levels(dbh))[elevel,]))
              coefs <- if (nonwave) c1[c(2,4,6)] else c1[c(1,3,5)]
              m <- coefs[1]*meandbh^coefs[2]
              k <- coefs[3]
              llines(0:300,dnbinom(0:300,mu=m,size=k),col="darkgray")
            }))


###################################################
### chunk number 110:  eval=FALSE
###################################################
## ## pretty but too complex for now
## print(qqmath(~TOTCONES|dbh*WAVE_NON,data=X,
##              xlim=range(X$TOTCONES),ylim=range(X$TOTCONES),
##              panel=function(x,distribution,...) {
##                nonwave = (current.row()==1)
##                elevel = current.column()
##                meandbh <- mean(levels(dbh)[[elevel]])
##                coefs <- if (nonwave) c1[c(2,4,6)] else c1[c(1,3,5)]
##                m <- exp(coefs[1]+meandbh*coefs[2])
##                k <- coefs[3]
##                dfun <- function(q) {
##                  qnbinom(q,mu=m,size=k)
##                }
##                panel.qqmath(x,distribution=dfun,...)
##                panel.qqmathline(x,distribution=dfun,...)
##              }))
## 


