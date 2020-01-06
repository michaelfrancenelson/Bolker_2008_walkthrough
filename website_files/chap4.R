###################################################
### chunk number 1: 
###################################################
source("chapskel.R")
library(emdbook)


###################################################
### chunk number 2: 
###################################################
op <- par(cex=1.5,las=1,bty="l",lwd=2,mfrow=c(1,2))
ziunif = function(x,p,v,N) {
  ifelse(x==0,(1-v)+v/(N+1),v/(N+1))
}
zibinom = function(x,p,v,N) {
  ifelse(x==0,(1-v)+v*dbinom(0,prob=p,size=N),
         v*dbinom(x,prob=p,size=N))
}
N=5
p=0.4
v=0.7
barplot(ziunif(0:N,p,v,N),main="Zero-inflated uniform",
        xlab="Taken",ylab="Probability",names=0:5)
barplot(zibinom(0:N,p,v,N),main="Zero-inflated binomial",
        xlab="Taken",ylab="Probability",names=0:5)

par(op)


###################################################
### chunk number 3: 
###################################################
r <- 0
d <- acos(r)
scale <- c(0.5,0.3)
npoints <- 100
centre <- c(0.5,0.5)
a <- seq(0, 2 * pi, len = npoints)
m <- matrix(c(scale[1] * cos(a + d/2) + centre[1], 
              scale[2] * cos(a - d/2) + centre[2]), npoints, 2)
e <- 0.05
hyp_pts = matrix(c(0.37,1.04,
  1+e,0.8+e,
  1,-e,
  0.4,-e,
  -e,0.25),
  byrow=TRUE,ncol=2)
lab.pts = matrix(c(0.091,0.255,0.597,0.557,
  0.869,0.709,0.549,0.511,
  0.170,0.22,
  ##y
  0.865,0.613,
  0.932,0.698,0.191,0.477,
  0.087,0.277,0.077,0.31),
  ncol=2)
##hyp_pts <- hyp_pts[c(5,1:4),]
## lab.pts <- lab.pts[c(5,1:4),]
par(mar=c(0.2,0.2,0.2,0.2))
plot(1,1,type="n",xlim=c((-e),1+e),ylim=c(-e,1+e),ann=FALSE,
     xlab="",ylab="",axes=FALSE,xaxs="i",yaxs="i")
box()
polygon(m,col="lightgray",lwd=2)
polygon(c(-e,0.5,0.4,-e),c(0.25,0.5,-e,-e),density=8,angle=0,
        col="darkgray")
lines(m,lwd=2)
segments(rep(0.5,nrow(hyp_pts)),rep(0.5,nrow(hyp_pts)),
         hyp_pts[,1],hyp_pts[,2])
##text(lab.pts[,1],lab.pts[,2],1:10)
for(i in 1:5) {
  r = 2*i-1
  r2 = 2*i
  text(lab.pts[r,1],lab.pts[r,2], substitute(H[x],
                                list(x=i)),adj=0,cex=2)
  text(lab.pts[r2,1],lab.pts[r2,2], substitute(D*intersect("","","")*H[x],
                                list(x=i)),adj=0,cex=2)
}


###################################################
### chunk number 4: 
###################################################
probI=1e-6
probfp=1e-3


###################################################
### chunk number 5: 
###################################################
op <- par(cex=1.5,las=1,bty="l",lwd=2,yaxs="i")
gcols <- gray(c(0.2,0.8))
b1 <- barplot(t(matrix(c(1/3,1/3,1/3,1/4,1/4,1/2),ncol=2)),beside=TRUE,
              xlab="Predator",ylab="Probability",space=c(0.2,2),
              col=gcols,yaxs="i")
axis(side=1,at=colMeans(b1),c("raccoon","squirrel","snake"))
segments(b1[1,1],0.4,b1[2,2],0.4)
text((b1[1,1]+b1[2,2])/2,0.45,"mammalian")
par(xpd=NA)
legend(2,0.57,c("by species","by group"),
       ncol=2,fill=gcols,bty="n")
par(op)


###################################################
### chunk number 6: 
###################################################
minx <- 10
maxx <- 100
dx <- maxx-minx
dlx <- log(maxx/minx)
dlx10 <- log10(maxx/minx)
xlim <- c(0,110)
Lxlim <- c(9,110)
op <- par(cex=2,las=1,bty="l",lwd=2,mfrow=c(1,2),
          mar=c(5,4,3,0.5)+0.1,yaxs="i")
curve(ifelse(x>minx & x<maxx,1/dx,0),from=xlim[1],to=xlim[2],
      xlab="Mass",ylab="Probability density",ylim=c(0,0.04),main="linear scale",type="s",
      cex.main=1.5,xlim=c(0,100),n=400,axes=FALSE)
axis(side=1,at=c(10,100))
axis(side=2,at=c(0,0.02,0.04))
box()
curve(ifelse(x>minx & x<maxx,1/x*1/dlx,0),from=xlim[1],to=xlim[2],
      lty=2,n=400,add=TRUE)
legend("topright",c("uniform","log-uniform"),lty=1:2)
curve(ifelse(x>log(minx) & x<log(maxx),1/dlx,0),
      from=log(Lxlim[1]),to=log(Lxlim[2]),
            ylim=c(0,1.2),
      xlab="Log mass",ylab="",axes=FALSE,lty=2,main="log scale",type="s",
      cex.main=1.5,n=400)
curve(ifelse(x>log(minx) & x<log(maxx),exp(x)/dx,0),from=log(Lxlim[1]),to=log(Lxlim[2]),
      add=TRUE,n=400)
axis(side=1,
     at=log(c(10,100)),
     labels = paste("log(",c(10,100),")",sep=""))
axis(side=2,at=c(0,0.5,1))
box()
par(op)


###################################################
### chunk number 7: 
###################################################
integrate(function(p)dbinom(2,prob=p,size=3),lower=0,upper=1)


###################################################
### chunk number 8: 
###################################################
op=par(mfrow=c(2,2),mgp=c(2.5,1,0),mar=c(3.2,4.3,1,1),bty="l",las=1,
  cex=1.5,yaxs="i")
plot(0:5,dbinom(0:5,size=5,prob=0.3),
     main="",
     xlab="",
     ylab="Probability",type="h",lwd=4,
     axes=FALSE,ylim=c(0,0.4))
axis(side=1)
axis(side=2,at=c(0,0.4))

box()
plot(0:5,pbinom(0:5,size=5,prob=0.3),
     main="",type="h",lwd=4,
     xlab="",
     ylab="Cumulative\nprobability",ylim=c(0,1),axes=FALSE)
box()
axis(side=1)
axis(side=2,at=c(0,1))
curve(dexp(x,1.5),from=0,to=6,ylab="Probability density",
      main="",xaxs="i",yaxs="i",axes=FALSE,ylim=c(0,1.5))
axis(side=1,at=c(0,6))
axis(side=2,at=c(0,1.5))
box()
curve(pexp(x,1.5),from=0,to=6,ylab="Cumulative\nprobability",
      main="",xaxs="i",yaxs="i",axes=FALSE,ylim=c(0,1))
axis(side=1,at=c(0,6))
axis(side=2,at=c(0,1))
box()
par(op)


###################################################
### chunk number 9: 
###################################################
multdiscplot = function(x,parvec,parname="prob",dfun="dbinom",col,
  displ=0.2,xlab="x",ylab="Probability",lwd=7,axes=TRUE,...) {
  op = par(lend=2)
  tmpfun = function(p) {
      params = c(list(...),x=list(x),p=p)
      names(params)[length(params)]=parname
      do.call(dfun,params)
    }
  vals=sapply(parvec,tmpfun)
  n =length(parvec)
  xoff = seq(-displ*n/2,displ*n/2,length=n)
  plot(range(x)+range(xoff),range(vals),type="n",xlab=xlab,ylab=ylab,
       axes=axes)
  for (i in 1:length(parvec)) {
    points(x-xoff[i],vals[,i],col=col[i],type="h",lwd=lwd)
  }
  par(op)
}

multdiscplot2 = function(x,parmat,parnames=c("mu","size"),
  dfun="dnbinom",col,displ=0.2,xlab="x",ylab="Probability",lwd=7,...) {
  op = par(lend=2)
  n =nrow(parmat)
  vals = matrix(nrow=length(x),ncol=nrow(parmat))
  for (i in 1:n) {
    params = c(list(...),x=list(x),p1=parmat[i,1],p2=parmat[i,2])
    lp = length(params)
      names(params)[(lp-1):lp] = parnames
      vals[,i] = do.call(dfun,params)
    }
  xoff = seq(-displ*n/2,displ*n/2,length=n)
  plot(range(x)+range(xoff),range(vals),type="n",xlab=xlab,ylab=ylab)
  for (i in 1:n) {
    points(x-xoff[i],vals[,i],col=col[i],type="h",lwd=lwd)
  }
  par(op)
}



###################################################
### chunk number 10: 
###################################################
op <- par(lwd=2,cex=1.5,las=1,bty="l",mgp=c(3,1,0),yaxs="i")
## palette(gray((0:8)/8))
multdiscplot(0:10,parvec=c(0.1,0.5,0.9),size=10,col=gray((0:3)/3),lwd=5,
             xlab="# of successes",displ=0.1)
par(xpd=NA)
text(c(1,5,9),
     c(0.42,0.27,0.42),
     paste("p=",c(0.1,0.5,0.9),sep=""),
     col=gray((0:3)/3))
#legend(2,0.5,paste("N=10, p=",c(0.1,0.5,0.9),sep=""),
#       lty=1,lwd=5,col=gray((0:3)/3))
par(xpd=TRUE)
par(op)


###################################################
### chunk number 11: 
###################################################
par(lwd=2,cex=1.5,las=1,bty="l",mgp=c(2.5,1,0),yaxs="i")
pvec = c(0.8,3,12)
multdiscplot(0:20,parvec=pvec,parname="lambda",
             dfun="dpois",col=gray((0:3)/3),lwd=3,
             xlab="# of events",displ=0.1)
## construct legend by hand
## legend(10,0.4,
##        c(expression(lambda==0.8),
##          expression(lambda==3),
##          expression(lambda==12)),
##        lwd=3,col=1:3,lty=1)
text(c(1.15,2.7,10.2),
     c(0.42,0.26,0.15),
        c(expression(lambda==0.8),
          expression(lambda==3),
          expression(lambda==12)),
     col=gray((0:3)/3),adj=0)
## leg.y = seq(0.38,length=3,by=-0.02)
## segments(rep(10,3),leg.y,rep(12,3),leg.y,
##          col=gray((0:3)/3),lwd=3)
## for (i in 1:3)
##   text(rep(12.4,3),leg.y[i],
##       substitute(lambda==z,list(z=pvec[i])),adj=0)


###################################################
### chunk number 12: 
###################################################
par(lwd=2,cex=1.5,las=1,bty="l",mgp=c(2.5,1,0),yaxs="i")
xvec = 0:10
pvec = c(10,1,0.1)
multdiscplot(xvec,parvec=pvec,
             dfun="dnbinom",
             mu=2,parname=c("size"),
             col=gray((0:3)/3),lwd=4,
             displ=0.1)
##legend(4,0.6,
text(c(3.53,0.37,0.16),
     c(0.18,0.32,0.59),
       c(expression(k==10),
         expression(k==1),
         expression(k==0.1)),
     adj=0,
     col=gray((0:3)/3),lty=1,lwd=3)
##leg.y = seq(0.6,length=3,by=-0.06)
##segments(rep(4,3),leg.y,rep(6,3),leg.y,
##         col=gray((0:3)/3),lwd=3)
##for (i in 1:3)
##  text(rep(6.4,3),leg.y[i],
##       substitute(list(mu==2,k==z),list(z=pvec[i])),adj=0)


###################################################
### chunk number 13: 
###################################################
par(lwd=2,cex=1.5,las=1,bty="l",yaxs="i")
pvec=c(0.2,0.5)
multdiscplot(0:20,parvec=pvec,
             dfun="dgeom",
             parname="prob",
             col=gray((0:2)/2),lwd=3,
             displ=0.1,
             xlab="Survival time")
text(c(5.4,0.5),
     c(0.082,0.41),
     paste("p=",pvec,sep=""),
     adj=0,
     col=gray((0:2)/2))


###################################################
### chunk number 14: 
###################################################
op <- par(lwd=2,cex=1.5,las=1,bty="l",mgp=c(3,1,0),yaxs="i")
multdiscplot(0:10,parvec=c(0.5,5),parname="theta",
             dfun="dbetabinom",prob=0.5,size=10,
             col=gray((0:3)/3),lwd=5,
             xlab="# of successes",displ=0.1)
par(xpd=NA)
text(0.5,0.25,expression(theta==0.5),col=gray(0),pos=4)
text(4,0.15,  expression(theta==5),col=gray(1/3),pos=4)
#legend(2,0.5,paste("N=10, p=",c(0.1,0.5,0.9),sep=""),
#       lty=1,lwd=5,col=gray((0:3)/3))
par(xpd=TRUE)
par(op)


###################################################
### chunk number 15: 
###################################################
op <- par(lwd=2,cex=1.5,las=1,bty="l",yaxs="i")
par(mar=c(5,4,2,6))
curve(dunif(x,0,1),from=-0.3,to=2.3,
   xlab="Value",ylab="Probability density",
      axes=FALSE,type="S")
axis(side=2)
axis(side=1)
curve(dunif(x,0.5,2),add=TRUE,lty=2,type="S")
par(xpd=NA)
text(c(0.5,1.6),
     c(1.05,1/1.5+0.05),
     c("U(0,1)","U(0.5,2.5)"))
par(xpd=FALSE)
box()
par(op)


###################################################
### chunk number 16: 
###################################################
op <- par(lwd=2,cex=1.5,las=1,bty="l",yaxs="i")
par(mar=c(5,4,2,6))
curve(dnorm(x,0,1),from=-10,to=12,
   xlab="Value",ylab="Probability density",xlim=c(-14,10),
      axes=FALSE)
col2 = gray(0.7)
axis(side=2)
axis(side=1,at=seq(-10,10,by=5))
curve(dnorm(x,0,3),add=TRUE,col=1,lty=2)
curve(dnorm(x,2,1),add=TRUE,col=col2,lty=1)
curve(dnorm(x,2,3),add=TRUE,col=col2,lty=2)
par(xpd=NA)
#legend(5,0.4,
text(c(-1,-3),
     c(0.3,0.1),
     adj=1,
       c(expression(list(mu==0,sigma==1)),
         expression(list(mu==0,sigma==3))))
text(c(4,6),
     c(0.3,0.1),
     adj=0,
     c(expression(list(mu==2,sigma==1)),
       expression(list(mu==2,sigma==3))),
     col=col2)
par(op)


###################################################
### chunk number 17: 
###################################################
op <- par(lwd=2,cex=1.5,las=1,bty="l")
curve(dgamma(x,1,1),from=1e-4,to=25,
   xlab="Value",ylab="Prob. density")
cols = c("black","gray")
curve(dgamma(x,2,1),add=TRUE,col=cols[1],lty=2,from=0)
curve(dgamma(x,5,1),add=TRUE,col=cols[1],lty=3,from=0)
curve(dgamma(x,1,1/3),add=TRUE,col=cols[2],lty=1,from=0)
curve(dgamma(x,2,1/3),add=TRUE,col=cols[2],lty=2,from=0)
curve(dgamma(x,5,1/3),add=TRUE,col=cols[2],lty=3,from=0)
par(xpd=NA)
legend(5,1,
       paste("shape=",c(1,2,5,1,2,5),", scale=",
             rep(c("1","1/3"),c(3,3)),sep=""),
       col=rep(cols,each=3),lty=rep(1:3,2))
par(op)


###################################################
### chunk number 18: 
###################################################
par(lwd=2,cex=1.5,las=1,bty="l",yaxs="i",xaxs="i",mgp=c(2.5,1,0))
curve(dexp(x,1),from=1e-4,to=15,
   xlab="Value",ylab="Probability density")
curve(dexp(x,0.5),from=1e-4,add=TRUE,lty=2)
curve(dexp(x,0.2),from=1e-4,add=TRUE,lty=3)
text(c(0.73,2.17,5.7),
     c(0.67,0.2,0.11),
     c(expression(lambda==1),
       expression(lambda==1/2),
       expression(lambda==1/5)),
     col=1,adj=0)



###################################################
### chunk number 19: 
###################################################
op <- par(lwd=2,cex=1.5,las=1,bty="l",yaxs="i",mgp=c(2.5,1,0),
          mar=c(5,4,2,5)+0.1)
curve(dbeta(x,1,1),from=0,to=1,
   xlab="Value",ylab="Probability density",lty=1,
      ylim=c(0,5))
curve(dbeta(x,5,5),add=TRUE,from=0,to=1,col=1,lty=2)
curve(dbeta(x,5,1),add=TRUE,from=0,to=1,col=1,lty=3)
curve(dbeta(x,1,5),add=TRUE,from=0,to=1,col=1,lty=4)
curve(dbeta(x,0.5,0.5),add=TRUE,from=0,to=1,col=1,lty=5)
par(xpd=NA)
text(0.07,4.15,c("a=1, b=5"),adj=0)
text(0.5,2.7,c("a=5, b=5"))
text(c(1,1,0.94),
     c(4.58,2.53,0.83),
     c("a=5, b=1",
       "a=0.5, b=0.5",
       "a=1, b=1"),
     adj=0)
par(xpd=FALSE)
par(op)


###################################################
### chunk number 20: 
###################################################
par(lwd=2,cex=1.5,las=1,bty="l")
curve(dlnorm(x,0,0.5),from=0,to=12,
   xlab="Value",ylab="Probability density",n=300,lty=2,
      ylim=c(0,2))
curve(dlnorm(x,0,1),add=TRUE,col=1,lty=3,n=300,from=0)
curve(dlnorm(x,0,0.2),add=TRUE,col=1,lty=1,n=500,from=0)
curve(dlnorm(x,2,0.5),add=TRUE,col=cols[2],lty=2,from=0)
curve(dlnorm(x,2,1),add=TRUE,col=cols[2],lty=3,from=0)
curve(dlnorm(x,2,0.2),add=TRUE,col=cols[2],lty=1,from=0)
legend(3,2,
       c(expression(paste(mu[log]==0,", ",sigma[log]==0.2)),
         expression(paste(mu[log]==0,", ",sigma[log]==0.5)),
         expression(paste(mu[log]==0,", ",sigma[log]==1)),
         expression(paste(mu[log]==2,", ",sigma[log]==0.2)),
         expression(paste(mu[log]==2,", ",sigma[log]==0.5)),
         expression(paste(mu[log]==2,", ",sigma[log]==1))),
       col=rep(cols,each=3),lty=rep(1:3,2))


###################################################
### chunk number 21: 
###################################################
par(lwd=2,cex=1.5,las=1,bty="l")
curve(0.7*dnorm(x,1,2)+0.3*dnorm(x,5,1),from=-8,to=12,
      ylab="Probability density",xlab="")
mtext(side=1,"x",line=2,cex=1.5)
curve(0.7*dnorm(x,1,2),lty=2,add=TRUE)
curve(0.3*dnorm(x,5,1),lty=2,add=TRUE)


###################################################
### chunk number 22: 
###################################################
sh = 4
op=par(mfrow=c(1,2))
par(mar=c(5.1,4.5,1,0.5),cex=1.5,lwd=2,las=1,bty="l",yaxs="i")
curve(dgamma(x,shape=sh),from=0,to=20,ylab="",xlab="",lwd=2,axes=FALSE)
mtext(side=1,line=2,cex=2,"x")
axis(side=1,labels=FALSE)
axis(side=2,labels=FALSE)
box()
xvec = seq(0,5,length=100)
polygon(c(xvec,rev(xvec)),c(rep(0,length(xvec)),dgamma(rev(xvec),shape=sh)),
        col="gray",border=NA)
curve(dgamma(x,shape=sh),from=0,to=20,lwd=2,add=TRUE) ## redraw
abline(v=5,lty=3)
mtext(side=1,line=0.7,at=5,expression(x[0]),cex=1.5)
abline(h=dgamma(5,shape=sh),lty=2,lwd=2)
mtext(side=2,at=dgamma(5,shape=sh),las=1,expression(ddist(x[0])),
      line=1,cex=1.5)
text(10,0.1,expression(pdist(x[0])),cex=1)
arrows(6.5,0.1,3,0.1,angle=20,lwd=2)
mtext(side=2,at=0.07,"Probability\ndensity",cex=2,line=1,
      las=0)
set.seed(1001)
r1 <- rgamma(10,shape=sh)
points(r1,rep(0.01,10),cex=1.5)
text(11.7,0.03,"rdist(10)",adj=0,cex=1)
arrows(11,0.03,7,0.02,lwd=2,angle=20)
curve(pgamma(x,shape=sh),from=0,to=20,xlab="",ylab="",lwd=2,axes=FALSE)
axis(side=1,labels=FALSE)
axis(side=2,labels=FALSE)
mtext(side=1,line=2,cex=2,"x")
box()
abline(v=5,lty=3)
mtext(side=1,line=0.3,at=5,expression(x[0]),cex=1.5)
abline(h=pgamma(5,shape=sh),lty=2,lwd=2)
mtext(side=2,at=pgamma(5,shape=sh),las=1,expression(pdist(x[0])),
      line=1,cex=1.5)
abline(v=qgamma(0.95,shape=sh),lty=4,lwd=2)
mtext(side=2,at=0.95,las=1,0.95,
      line=par("mgp")[2],cex=1.5)
segments(par("usr")[1],0.95,qgamma(0.95,shape=sh),0.95,lty=4,lwd=2)
mtext(side=1,at=qgamma(0.95,shape=sh),text="qdist(0.95)",line=1,
      cex=1.5,adj=0.3)
mtext(side=2,at=0.4,adj=0.5,"Cumulative\ndistribution",cex=2,line=1,las=0)


###################################################
### chunk number 23: 
###################################################
set.seed(1001)


###################################################
### chunk number 24: 
###################################################
z <- rnbinom(1000,mu=10,size=0.9)


###################################################
### chunk number 25: 
###################################################
head(z)


###################################################
### chunk number 26: 
###################################################
maxz <- max(z)


###################################################
### chunk number 27: 
###################################################
f <- factor(z,levels=0:maxz)
plot(f)


###################################################
### chunk number 28: 
###################################################
obsprobs <- table(f)/1000
plot(obsprobs)


###################################################
### chunk number 29: 
###################################################
tvals <- dnbinom(0:maxz,size=0.9,mu=10)
points(0:maxz,tvals)


###################################################
### chunk number 30: 
###################################################
pnbinom(30,size=0.9,mu=10,lower.tail=FALSE)


###################################################
### chunk number 31: 
###################################################
qnbinom(0.95,size=0.9,mu=10)


###################################################
### chunk number 32: 
###################################################
qnbinom(c(0.025,0.975),size=0.9,mu=10)


###################################################
### chunk number 33: 
###################################################
mu <- 10
k <- 0.9
c(mu,mean(z))
c(mu*(1+mu/k),var(z))
c(qnbinom(0.95,size=k,mu=mu),quantile(z,0.95))


###################################################
### chunk number 34: 
###################################################
z <- rlnorm(1000,meanlog=2,sdlog=1)


###################################################
### chunk number 35: 
###################################################
hist(z,breaks=100,freq=FALSE)
lines(density(z,from=0),lwd=2)


###################################################
### chunk number 36: 
###################################################
curve(dlnorm(x,meanlog=2,sdlog=1),add=TRUE,lwd=2,from=0,
      col="darkgray")


###################################################
### chunk number 37: 
###################################################
plnorm(30,meanlog=2,sdlog=1,lower.tail=FALSE)
qlnorm(c(0.025,0.975),meanlog=2,sdlog=1)


###################################################
### chunk number 38: 
###################################################
meanlog <- 2
sdlog <- 1
c(exp(meanlog+sdlog^2/2),mean(z))
c(exp(2*meanlog+sdlog^2)*(exp(sdlog^2)-1),var(z))
c(qlnorm(0.95,meanlog=meanlog,sdlog=sdlog),quantile(z,0.95))


###################################################
### chunk number 39: 
###################################################
hist(log(z),freq=FALSE,breaks=100)
curve(dnorm(x,mean=meanlog,sd=sdlog),add=TRUE,lwd=2)


###################################################
### chunk number 40: 
###################################################
u1 <- runif(1000)
z <- ifelse(u1<0.3,rnorm(1000,mean=1,sd=2),
            rnorm(1000,mean=5,sd=1))
hist(z,breaks=100,freq=FALSE)


###################################################
### chunk number 41: 
###################################################
curve(0.3*dnorm(x,mean=1,sd=2)+0.7*dnorm(x,mean=5,sd=1),
      add=TRUE,lwd=2)


###################################################
### chunk number 42: 
###################################################
dzinbinom = function(x,mu,size,zprob) {
  ifelse(x==0,
         zprob+(1-zprob)*dnbinom(0,mu=mu,size=size),
         (1-zprob)*dnbinom(x,mu=mu,size=size))
}


###################################################
### chunk number 43: 
###################################################
rzinbinom = function(n,mu,size,zprob) {
  ifelse(runif(n)<zprob,
         0,
         rnbinom(n,mu=mu,size=size))
}


###################################################
### chunk number 44: 
###################################################
k <- 3
mu <- 10
lambda <- rgamma(1000,shape=k,scale=mu/k)
z <- rpois(1000,lambda)
P1 <- table(factor(z,levels=0:max(z)))/1000 
plot(P1)
P2 <- dnbinom(0:max(z),mu=10,size=3)
points(0:max(z),P2)


###################################################
### chunk number 45: 
###################################################
mlog <- mean(log(lambda))
sdlog <- sd(log(lambda))
lambda2 <- rlnorm(1000,meanlog=mlog,sdlog=sdlog)
z2 <- rpois(1000,lambda2)
P3 <- table(factor(z2,levels=0:max(z)))/1000
matplot(0:max(z),cbind(P1,P3),pch=1:2)
lines(0:max(z),P2)


