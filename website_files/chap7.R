###################################################
### chunk number 1: 
###################################################
library(emdbook)
library(plotrix)
library(bbmle)
library(MASS)
library(MCMCpack)
library(coda)
library(R2WinBUGS)
library(Hmisc)
source("chapskel.R")


###################################################
### chunk number 2: 
###################################################
x=c(1,1.5,2,2.5,3,4,5,6,10,15,20)
y=c(6,3.2,2,-3,3,2,0,-1,2,4,5)
y2 = spline(x,y,n=200)
plot(y2$x,y2$y,type="l",xlab="",ylab="",ylim=c(-3,15),
     axes=FALSE,xlim=c(0,20))
box()
axis(side=3,line=-11,at=seq(0,5,by=0.5),labels=rep("",11))
axis(side=3,line=-8,at=seq(0,20,by=5),labels=rep("",5))
axis(side=3,line=-5,at=seq(5,20,by=0.5),labels=rep("",31))
axis(side=3,line=-2,at=seq(0,20,by=1),labels=rep("",21))
mtext(side=3,line=-10,at=0,adj=0,"sampling grid 4")
mtext(side=3,line=-7,at=0,adj=0,"sampling grid 3")
mtext(side=3,line=-4,at=0,adj=0,"sampling grid 2")
mtext(side=3,line=-1,at=0,adj=0,"sampling grid 1")
mtext(side=3,line=-12,at=0,adj=0,expression(p[lower]))
mtext(side=3,line=-12,at=20,adj=1,expression(p[upper]))
mtext(side=3,line=-9,at=12.5,adj=0.5,expression(Delta*p))
yht <- 8.5
arrows(c(11.5,13.5),c(yht,yht),c(10,15),c(yht,yht),length=0.1,angle=20)
points(y2$x[which.min(y2$y)],min(y2$y),pch=1)
restr <- y2$x>5
xrestr <- y2$x[restr]
yrestr <- y2$y[restr]
points(xrestr[which.min(yrestr)],min(yrestr),pch=16)
#points(x,y)


###################################################
### chunk number 3: 
###################################################
data(MyxoTiter_sum)
myxdat = subset(MyxoTiter_sum,grade==1)


###################################################
### chunk number 4: 
###################################################
gm = mean(myxdat$titer)
cv = var(myxdat$titer)/mean(myxdat$titer)
shape0=gm/cv
scale0=cv


###################################################
### chunk number 5: 
###################################################
shapevec = 10:100
scalevec = seq(0.01,0.3,by=0.01)


###################################################
### chunk number 6: 
###################################################
gammaNLL1 = function(shape,scale) {
  -sum(dgamma(myxdat$titer,shape=shape,scale=scale,log=TRUE))
}


###################################################
### chunk number 7: 
###################################################
surf = matrix(nrow=length(shapevec),ncol=length(scalevec)) 
for (i in 1:length(shapevec)) {
  for (j in 1:length(scalevec)) {
    surf[i,j] = gammaNLL1(shapevec[i],scalevec[j])
  }
}


###################################################
### chunk number 8: 
###################################################
contour(shapevec,scalevec,log10(surf))


###################################################
### chunk number 9:  eval=FALSE
###################################################
## curve3d(log10(gammaNLL1(x,y)),from=c(10,0.01),to=c(100,0.3),n=c(91,30),sys3d="image")


###################################################
### chunk number 10:  eval=FALSE
###################################################
## gridsearch2d(gammaNLL1,v1min=10,v2min=0.01,v1max=100,v2max=0.3,logz=TRUE)


###################################################
### chunk number 11: 
###################################################
## run Newton's method on likelihood for
##  binomial data
data(ReedfrogPred)
x = subset(ReedfrogPred,pred=="pred" & density==10 & size=="small")
k = x$surv
N = 10
binomlik = function(p,x=k) {
  -sum(dbinom(k,prob=p,size=N,log=TRUE))
}
## D(expression(k*log(p)+(N-k)*log(1-p)),"p")
binomlik.deriv = function(p,x=k) {
  -sum(k/p-(N-k)/(1-p))
}
## D(D(expression(k*log(p)+(N-k)*log(1-p)),"p"),"p")
binomlik.deriv2 = function(p,x=k) {
  -sum(-k/p^2-(N-k)/(1-p)^2)
}
newton = function(start=0.6,graphics=TRUE,
    tol=1e-4,maxit=20,
    dfun=binomlik.deriv,
    dfun2=binomlik.deriv2,
    ...) {
  resmat = matrix(nrow=maxit,ncol=5)
  i = 1
  ## initial values
  x.guess = start
  x.dval = dfun(x.guess,...)
  x.dval2 = dfun2(x.guess,...)
  while (abs(x.dval)>tol && i<maxit) {  
    resmat[i,] = c(i,x.guess,x.dval,x.dval2,x.dval-x.dval2*x.guess)
    x.guess = x.guess - x.dval/x.dval2
    x.dval  = dfun(x.guess,...)
    x.dval2 = dfun2(x.guess,...)
    i = i+1
  }
  resmat[i,] = c(i,x.guess,x.dval,x.dval2,x.dval-x.dval2*x.guess) 
  colnames(resmat) = c("iteration","guess","dval","dval2","intercept")
  as.data.frame(resmat[1:i,])
}

plot.newton = function(n,xlab="p",xlim,ylim,npts=100,dfun=binomlik.deriv,
  dcol="darkgray",n.cex=3,pt.cex=2,
  g.lty=2,ylab="Derivative",...) {
  if (missing(xlim)) xlim=range(n$guess)
  pvec = seq(xlim[1],xlim[2],length=npts)
  dvec = sapply(pvec,dfun)
  if (missing(ylim)) ylim=range(c(n$dval,dvec))
  plot(n$guess,n$dval,type="n",
       xlab=xlab,ylab=ylab,
        xlim=xlim,ylim=ylim,...)
  lines(pvec,dvec,lwd=3,col=dcol)
  abline(h=0,lty=2)
  ng <- length(n$guess)
  apply(n[-ng,],1,function(x) {
     abline(a=x[5],b=x[4])
   })
  voff=0
  segments(n$guess,rep(0,nrow(n))-voff,n$guess,n$dval+voff,col="gray",lty=g.lty)
  points(n$guess,rep(0,nrow(n)),
         pch=21,cex=pt.cex,bg="white")
  points(n$guess,n$dval,
         cex=n.cex,pch=21,bg="white")
  points(n$guess,n$dval,pch=as.character(1:nrow(n)))
}


###################################################
### chunk number 12: 
###################################################
x0 = c(0,0.2,0.45,0.6,0.65,0.75,0.9,1)
y0 =c(-1,-0.5,0.039,0.29,0.41,0.53,0.6,0.7)
z = splinefun(x=x0,y=y0)
x1 = 0.7
y1 =z(x1)
y1d = z(x1,1)
zz = x1-y1/y1d
yint= y1-x1*y1d
dx = 0.1
z2 = y1-y1d*dx ## z(x=x1-dx)
z3 = y1+y1d*dx ## z(x=x1+dx)
##


###################################################
### chunk number 13: 
###################################################
op=par(lwd=2,bty="l",las=1,cex=1.5,mgp=c(2.5,1,0))
plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE,
     xlim=c(0.2,0.9),ylim=c(0,0.8))
abline(h=0)
curve(z(x),add=TRUE)  
abline(v=x1)
abline(a=yint,b=y1d)
points(x1,y1,pch=16)
segments(x1-0.2,y1,x1+0.1,y1,lty=3)
text(x1-0.21,y1,adj=1,expression(f(x[0])))
segments(c(x1-dx,x1+dx),c(z2,z2),c(x1+dx,x1+dx),c(z2,z3),col="gray")
text(x1+dx+0.01,z2+0.04,expression(f*minute(x[0])),adj=0)
mtext(side=1,at=x1,expression(x[0]),cex=1.5)
mtext(side=1,at=zz,expression(x[1]),cex=1.5)
points(zz,0)
par(xpd=NA)
arrows(zz,-0.02,x1,-0.02,code=3,length=0.1)
par(xpd=FALSE)
mtext(side=1,at=(zz+x1)/2,expression(f(x[0])/f*minute(x[0])),cex=1.5)
par(op)


###################################################
### chunk number 14: 
###################################################
op=par(mfrow=c(2,1),lwd=2,bty="l",las=1,cex=1.5,mgp=c(2.5,1,0))
n0 = newton(start=0.6)
n1 = n0[1:3,]
par(mar=c(0,4,2,2)+0.1)
plot.newton(n1,xlim=c(0.58,0.8),ylim=c(-28,10),
            xlab="",axes=FALSE,n.cex=2.5,pt.cex=1,
            g.lty=3,
            ylab="Derivative of -L")
axis(side=2)
axis(side=1,at=seq(0.6,0.8,by=0.05),labels=rep("",5))
text(0.6,0.5,"start",pos=3)
box()
par(mar=c(4,4,0,2)+0.1)
pvec=seq(0.6,0.8,length=100)
par(xpd=NA)
plot(pvec,sapply(pvec,binomlik),xlab="p (probability of success per trial)",
     ylab="Negative log-likelihood",type="l",lwd=3,ylim=c(7,10),
     axes=FALSE,xlim=c(0.58,0.8))
par(xpd=FALSE)
axis(side=1)
axis(side=2,at=7:10)
box()
points(n1$guess,sapply(n1$guess,binomlik),pch=21,cex=2.5,bg="white")
points(n1$guess,sapply(n1$guess,binomlik),pch=as.character(1:nrow(n1)))
par(op)


###################################################
### chunk number 15: 
###################################################
n2 <- n0[,2:4]
n2[,1] <- round(n2[,1],6)
n2[,2:3] <- round(n2[,2:3],3)
colnames(n2) <- c("Guess ($x$)","$f'(x)$","$f''(x)$")


###################################################
### chunk number 16: 
###################################################
latex(n2,file="",table.env=FALSE,title="")


###################################################
### chunk number 17: 
###################################################
data(MyxoTiter_sum)
myxdat = subset(MyxoTiter_sum,grade==1)
gammaNLL1 = function(shape,scale) {
  -sum(dgamma(myxdat$titer,shape=shape,scale=scale,log=TRUE))
}
gammaNLL2 = function(p) {
  -sum(dgamma(myxdat$titer,shape=p[1],scale=p[2],log=TRUE))
}
myxsurf <- curve3d(gammaNLL2(c(x,y)),from=c(15,0.02),to=c(60,0.3),
                   n=c(61,61),sys3d="none")
source("amoeba.R")
plotsimplex.2d <- function(x,cols=1:5,lty,...) {
  v <- matrix(x[1:9],byrow=TRUE,ncol=3)[,-3]
  code <- x[10]
  lines(rbind(v,v[1,]),col=cols[code+1],lty=ltyvec[code+1],...)
}
## 0: starting value
## 1: reflection
## 2: reflection + expansion
## 3: one-dimensional contraction
## 4: contract around lowest poin


###################################################
### chunk number 18: 
###################################################
op=par(lwd=2,bty="l",las=1,cex=1.5)
amres <- amoeba(gammaNLL2,c(20,0.05),info=TRUE)
ainfo = amres$info
##ainfo[,c(1,4,7)] = log10(ainfo[,c(1,4,7)])
contour(myxsurf$x,myxsurf$y,myxsurf$z,drawlabels=FALSE,col="gray",
        levels=seq(0,7000,by=50),xlab="",ylab="scale")
mtext(side=1,at=40,"shape",line=2,cex=1.5)
colvec <- rep(1,4) ## gray(c(0,0,0.5,0.7,0.9))
ltyvec <- c(1,1,2,3,4)
invisible(apply(ainfo,1,plotsimplex.2d,cols=colvec,lwd=2))
points(20,0.05,pch=21,bg="white")
text(20,0.05,pos=1,"start")
points(amres$estimate[1],amres$estimate[2],pch=21,bg="white")
text(amres$estimate[1],amres$estimate[2],pos=1,"end")
legend("bottomright",
       legend=c("reflection","reflect+expand","contract"),
       lty=ltyvec[2:4],col=colvec[2:4],
       bty="n",lwd=2,cex=0.7)
##       bg="white",
par(op)


###################################################
### chunk number 19:  eval=FALSE
###################################################
## contour(myxsurf$x,myxsurf$y,myxsurf$z,drawlabels=FALSE,col="gray",
##         levels=seq(0,2200,by=50),xlim=c(46,51),ylim=c(0.135,0.155))
## invisible(apply(ainfo,1,plotsimplex.2d,cols=colvec,lwd=2))


###################################################
### chunk number 20: 
###################################################
set.seed(1002)
metres  <- metropSB(gammaNLL2,c(20,0.05),retvals=TRUE,retfreq=1,
                  nmax=2500)


###################################################
### chunk number 21: 
###################################################
op=par(bty="l",las=1,cex=1.5,mgp=c(2.5,1,0))
myxsurf2 <- curve3d(gammaNLL2(c(x,y)),from=c(20,0.05),to=c(75,0.3),
                   n=c(61,61),sys3d="contour",levels=seq(0,2200,by=50),
                    drawlabels=FALSE,col="gray",
                    xlab="shape",ylab="")
mtext(side=2,at=0.17,"scale",line=3,cex=1.5,las=0)
points(metres$retvals[,"p1"],metres$retvals[,"p2"],
       col=c("gray","black")[metres$retvals[,"accept"]+1],
       cex=0.3)
tvec <- c(1,100,150,300,1000,1500,2000,2500)
n <- length(tvec)
x <- metres$retvals[tvec,"p1"]
y <- metres$retvals[tvec,"p2"]
symbols(x,y,rectangles=cbind(rep(2,n),rep(1,n)),
        inches=0.3,add=TRUE,
        bg="white")
points(metres$estimate[1],metres$estimate[2],pch=21,
       bg="white")
text(x,y,tvec,cex=0.5)


###################################################
### chunk number 22:  eval=FALSE
###################################################
## ## looking at contours of metrop output
## k1 = kde2d(log10(metres$retvals[,"p1"]),metres$retvals[,"p2"],n=50)
## contour(k1$x,k1$y,k1$z)


###################################################
### chunk number 23: 
###################################################
firstmin = which.min(metres$retvals[,"minval"])


###################################################
### chunk number 24:  eval=FALSE
###################################################
## MSBfit = metropSB(fn=gammaNLL2,start=c(20,0.05),nmax=2500)


###################################################
### chunk number 25: 
###################################################
ivec = seq(1,nrow(metres$retval),by=25)
M <- metres$retval[ivec,]
op=par(mfrow=c(2,2),cex=1.5,mgp=c(2.5,1,0),las=1)
par(mar=c(0,4.2,2.2,0))
matplot(ivec,M[,c(1,3)],type="l",col=1,lty=c(1,3),
        xlab="",ylab="",lwd=c(1,2),axes=FALSE)
box()
axis(side=2)
text(corner.loc(1,1),"shape",adj=1)
par(mar=c(0,0,2.2,4.2))
matplot(ivec,M[,c(2,4)],type="l",col=1,lty=c(1,3),ylab="",
        xlab="",lwd=c(1,2),axes=FALSE)
box()
axis(side=4)
text(corner.loc(-1,1),"scale",adj=0)
#par(mar=c(3.5,3.2,1,2)+0.1)
par(mar=c(5,4.2,0,0))
plot(ivec,M[,5]/mean(M[,5]),type="l",ylab="",
     xlab="Iterations",ylim=c(0,1.7),axes=FALSE)
axis(side=1,at=c(0,1000,2000))
axis(side=2,at=c(0,0.5,1.0))
box()
frac.accept =  filter(metres$retval[,9], rep(1, 40), sides = 1)/40
frac.accept =  frac.accept[ivec]
lines(ivec,frac.accept,lty=2)
text(corner.loc(1,1,yoff=0.12),"relative\njump size",adj=1)
text(corner.loc(1,-1,yoff=0.12),"fraction\naccepted",adj=1)
par(mar=c(5,0,0,4.2))
matplot(ivec,M[,c(7,8)],type="l",col=1,lty=c(1,3),ylab="",
        xlab="Iterations",lwd=c(1,2),axes=FALSE,ylim=c(36,50))
axis(side=1,at=c(0,1000,2000))
axis(side=4,at=seq(0,1500,by=500))
text(corner.loc(1,1,yoff=0.14),"negative\nlog likelihood",adj=1)
box()
par(op)


###################################################
### chunk number 26: 
###################################################
gammaNLL2B = function(p) {
  sum(dgamma(myxdat$titer,shape=p[1],scale=p[2],log=TRUE))
}
m3 <- MCMCmetrop1R(gammaNLL2B,theta.init=c(shape=20,scale=0.05),
                   thin=30,mcmc=30000,
                   optim.lower=rep(0.004,2),
                   optim.method="L-BFGS-B",tune=3)
colnames(m3) = c("shape","scale")


###################################################
### chunk number 27:  eval=FALSE
###################################################
## curve3d(gammaNLL2B(c(x,y)),from=c(20,0.05),to=c(70,0.3),
##         sys3d="contour",levels=seq(-200,0,by=10))
## ## since shape*scale = x*y approx constant, try with
## ##  mean = x' = x*y , var = y'=x*y^2 
## ##  so x = x'^2/y', y = y'/x'
## ## shape (20,0.05), scale (70,0.3)
## ## means mean (1,21), var (0.05,6.3)
## curve3d(sum(dgamma(myxdat$titer,shape=x^2/y,scale=y/x,log=TRUE)),
##         ##        from=c(4,0.1),to=c(10,21),
##         from=c(1,0.05),to=c(21,6.3),
##         sys3d="contour",levels=seq(-200,0,by=10))


###################################################
### chunk number 28: 
###################################################
m1 <- mcmc(metres$retvals[,1:2])
raftery.diag(m1,r=0.01)


###################################################
### chunk number 29: 
###################################################
set.seed(1002)
metres2  <- metropSB(gammaNLL2,c(70,0.2),retvals=TRUE,retfreq=1,
                  nmax=2500)
m2 <- mcmc(metres2$retvals[,1:2])


###################################################
### chunk number 30: 
###################################################
gelman.diag(mcmc.list(m1,m2))


###################################################
### chunk number 31: 
###################################################
library(R2WinBUGS)


###################################################
### chunk number 32: 
###################################################
titer <- myxdat$titer
n <- length(titer)
inits <- list(list(shape=100,rate=3),list(shape=20,rate=10))
testmyxo.bugs <- bugs(data=list("titer","n"),
                      inits,parameters.to.save=c("shape","rate"),
                      model.file="myxogamma.bug",
                      n.chains=length(inits),n.iter=5000)


###################################################
### chunk number 33: 
###################################################
testmyxo.bugs


###################################################
### chunk number 34: 
###################################################
testmyxo.coda <- as.mcmc(testmyxo.bugs)


###################################################
### chunk number 35: 
###################################################
palette(c("black","gray"))
plot(testmyxo.coda)
palette("default")


###################################################
### chunk number 36:  eval=FALSE
###################################################
## x = testmyxo.coda[[1]][,"shape"]
## plot(density(x))
## h1 = HPDinterval(testmyxo.coda)[[1]]["shape",]
## h2 = qcredint(x,verbose=TRUE,tol=0.001)
## sum(x>h1[1]&x<h1[2])/length(x)
## sum(x>h2[1]&x<h2[2])/length(x)


###################################################
### chunk number 37: 
###################################################
nr=1000
nc=2000
m <- matrix(nrow=nr,ncol=nc)
v <- numeric(nc)
t1a <- system.time({
  for (i in 1:nr) {
    for (j in 1:nc) {
      m[i,j] = rnorm(1)
    }
  }})
t1b <- system.time(m <- matrix(rnorm(nr*nc),nrow=nr,ncol=nc))
t2a <- system.time({
  for (i in 1:nr) {
    v[i] = 0
    for (j in 1:nc) {
      v[i] = v[i]+m[i,j]
    }
  }})
t2b <- system.time({
  v = numeric(nr)
  for (i in 1:nc) {
    v[i] = sum(m[,i])
  }})
t2c <- system.time(v <- apply(m,2,sum))
t2d <- system.time(v <- colSums(m))
tmpf <- function(t,d=0) { round(t["elapsed"],d)}


###################################################
### chunk number 38: 
###################################################
thrfun <- function(x,a,b,thresh) {
  ifelse(x<thresh,a,b)
}
## Simulate some ``data'' using this function.
## Establish a range of $x$ values, compute the
## deterministic $y$, and add normally distributed
## errors:
set.seed(1001)
data.x <- seq(0,5,length=20)
det.y <- thrfun(data.x,a=2,b=5,thresh=3)
data.y <- det.y+rnorm(length(data.x),0,0.5)
ssqfun <- function(pvec,x=data.x,y=data.y) {
  exp.y <- thrfun(x,pvec[1],pvec[2],pvec[3])
  sum((exp.y-data.y)^2)
}
ssqfunprof <- function(pvec,x=data.x,y=data.y,thresh) {
  exp.y <- thrfun(x,pvec[1],pvec[2],thresh)
  sum((exp.y-data.y)^2)
}
## Calculate the sum of squares as a function of 
## the threshold parameter, with the lower
## and upper values set at $p_L=2$ and $p_U=5$
## (these happen to be their true values
## and this is a \emph{slice} for the threshold
## parameter)
thrvec <- seq(0,5,length=100)
thrprof0 <- t(sapply(thrvec,
                  function(t) unlist(optim(fn=ssqfunprof,par=c(2,5),thresh=t))))
thrprof = thrprof0[,"value"]
ssqfunlogist <- function(pvec,x=data.x,y=data.y) {
  a <- pvec[1]
  b <- pvec[2]
  thresh <- pvec[3]
  slope <- pvec[4]
  exp.y <- a+(b-a)/(1+exp(-slope*(data.x-thresh)))
  sum((exp.y-data.y)^2)
}
O1 <- optim(fn=ssqfun,par=c(a=2,b=5,thresh=3))
O2 <- optim(fn=ssqfunlogist,par=c(a=2,b=5,thresh=3,slope=10))
ssqfunlogistprof <- function(pvec,x=data.x,y=data.y,thresh) {
  a <- pvec[1]
  b <- pvec[2]
  slope <- pvec[3]
  exp.y <- a+(b-a)/(1+exp(-slope*(x-thresh)))
  sum((exp.y-y)^2)
}
logistprof0 = t(sapply(thrvec,
  function(t) unlist(optim(fn=ssqfunlogistprof,par=c(a=2,b=5,slope=10),thresh=t,
                            method="BFGS",control=list(maxit=1000)))))


###################################################
### chunk number 39: 
###################################################
op <- par(mfrow=c(2,1),cex=1.5,mgp=c(2.5,1,0),las=1,lwd=2,bty="l")
par(mar=c(1,4,1,1)+0.1)
plot(data.x,data.y,ylab="y",xlab="",axes=FALSE)
axis(side=2)
box()
curve(ifelse(x<3,O1$par[1],O1$par[2]),add=TRUE,lty=1,type="s")
curve(O2$par[1]+(O2$par[2]-O2$par[1])/(1+exp(-O2$par[4]*(x-O2$par[3]))),add=TRUE,lty=2)
par(mar=c(4,4,3,1)+0.1)
plot(thrvec,thrprof,
  xlab="x",ylab="Negative log likelihood",type="s")
lines(thrvec,logistprof0[,"value"],lty=2)
text(c(0,0),c(thrprof[1],
              logistprof0[1,"value"])+5,
              c("threshold","logistic"),adj=0,xpd=NA)
par(op)


###################################################
### chunk number 40:  eval=FALSE
###################################################
## plot(data.x,data.y,xlab="x",ylab="y")
## with (as.list(logistprof0[10,]),
##       curve(par.a+(par.b-par.a)/(1+exp(-par.slope*(x-thrvec[10]))),add=TRUE))
## with (as.list(thrprof0[10,]),
##       curve(ifelse(x<thrvec[10],par1,par2),add=TRUE,col=2))


###################################################
### chunk number 41: 
###################################################
data(ReedfrogSizepred)
attach(ReedfrogSizepred,warn=FALSE)
logist2 <- function(x,sizep1,sizep2,sizep3) {
 exp(sizep1*(sizep3-x))/(1+exp(sizep2*sizep1*(sizep3-x)))
}
lik.logist2 <- function(sizep1,sizep2,sizep3) {
  p.pred <- logist2(TBL,sizep1,sizep2,sizep3)
  -sum(dbinom(Kill,size=10,prob=p.pred,log=TRUE))
}
mle.logist2 <- mle2(lik.logist2,
                  start=list(sizep1=0,sizep2=1,sizep3=12),method="Nelder-Mead")
mle.logist2B <- mle2(lik.logist2,start=list(sizep1=0,sizep2=1,sizep3=12),
                    method="BFGS",control=list(parscale=c(0.3,30,10)))
mle.logist2C <- mle2(lik.logist2,start=list(sizep1=0,sizep2=1,sizep3=12),
                    method="BFGS",control=list(maxit=1000))
mle.logist2D <- mle2(lik.logist2,start=as.list(coef(mle.logist2C)),
                    control=list(maxit=1000))
slice2 = calcslice(mle.logist2,mle.logist2C,n=1000)
detach(ReedfrogSizepred)


###################################################
### chunk number 42: 
###################################################
op <- par(cex=1.5,mgp=c(2.5,1,0),las=1,lwd=2,bty="l")
plot(slice2,type="l",ylab="Negative log-likelihood",xlab="")
## plot(slice2,type="l",ylab="Negative log-likelihood",xlab="",ylim=c(0,20))
points(0,-logLik(mle.logist2))
points(1,-logLik(mle.logist2C))
##abline(h=-logLik(mle.logist2C))
abline(h=-logLik(mle.logist2C))
abline(h=-logLik(mle.logist2C)+qchisq(0.95,1)/2,lty=2)
par(mar=c(2,2,0,0),bty="l",cex=0.7,lwd=1, mgp=c(1,0,0))
par(new=TRUE,fig=c(0.22,0.50,0.55,0.85))
sizeplot(ReedfrogSizepred$TBL,ReedfrogSizepred$Kill, xlim=c(5,30),ylim=c(-0.4,10),axes=FALSE,
         xlab="size",ylab="killed",cex.lab=1.5)
with(as.list(coef(mle.logist2)),curve(logist2(x,sizep1,sizep2,sizep3)*10,add=TRUE,lwd=2))
box()
par(new=TRUE,fig=c(0.7,0.98,0.55,0.85))
sizeplot(ReedfrogSizepred$TBL,ReedfrogSizepred$Kill, xlim=c(5,30),ylim=c(-0.4,10),axes=FALSE,
         xlab="size",ylab="killed",cex.lab=1.5)
with(as.list(coef(mle.logist2C)),curve(logist2(x,sizep1,sizep2,sizep3)*10,add=TRUE,lwd=2))
box()
par(op)


###################################################
### chunk number 43:  eval=FALSE
###################################################
## weiblikfun <- function(shape,scale) {
##   -sum(log(pweibull(dat$d2,shape,scale)-
##            pweibull(dat$d1,shape,scale)))
## }


###################################################
### chunk number 44: 
###################################################
library(emdbookx)
data(GobySurvival)
dat = subset(GobySurvival,exper==1 & density==9 & qual>median(qual))
time = (dat$d1+dat$d2)/2


###################################################
### chunk number 45: 
###################################################
weiblikfun = function(shape,scale) {
  -sum(dweibull(time,shape=shape,scale=scale,log=TRUE))
}


###################################################
### chunk number 46: 
###################################################
w1 <- mle2(weiblikfun,start=list(shape=1,scale=mean(time)))


###################################################
### chunk number 47: 
###################################################
meanfun = function(shape,scale) { scale*gamma(1+1/shape) }


###################################################
### chunk number 48: 
###################################################
weiblikfun2 <- function(shape,mu) {
  scale <- mu/gamma(1+1/shape)
  -sum(dweibull(time,shape=shape,scale=scale,log=TRUE))
}
w2 <- mle2(weiblikfun2,start=list(shape=1,mu=mean(time)))


###################################################
### chunk number 49: 
###################################################
## compute the contour surface
smat <- curve3d(weiblikfun,
                from=c(0.5,5),
                to=c(1.2,25),
                n=c(40,40),sys3d="none")
#%Now let's superimpose on this plot contours that represent
#%the mean survival time.
#%Calculate the mean values for all combinations
#%of shape and scale:
#%The \code{outer} command is a quick way to generate matrices,
#%but it only works for vectorizable functions.
mmat <- outer(smat$x,smat$y,meanfun)  ## shortcut


###################################################
### chunk number 50: 
###################################################
op=par(lwd=2,bty="l",las=1,cex=1.5,mgp=c(2.5,1,0))
alevels = c(0.8,0.9,0.95,0.99)
mlevels = seq(6,38,by=4)
contour(smat$x,smat$y,smat$z,levels=qchisq(alevels,1)/2-logLik(w1),
        labels=alevels,
        xlab="Shape",ylab="Scale",col="darkgray",labcex=1)
points(coef(w1)[1],coef(w1)[2])
contour(smat$x,smat$y,mmat,add=TRUE,col=1,levels=mlevels,labcex=1,method="edge")
p2 = profile(w2)
p2B = p2@profile$mu$par.vals
shape <- p2B[,1]
mu <- p2B[,2]
scale <- mu/gamma(1+1/shape)
lines(shape,scale,lty=3)
shvals <- approx(p2@profile$mu$z,p2B[,"shape"],xout=c(-1.96,1.96))$y
muvals <- approx(p2@profile$mu$z,p2B[,"mu"],xout=c(-1.96,1.96))$y
contour(smat$x,smat$y,mmat,add=TRUE,col=1,levels=round(muvals,1),lty=2,drawlabels=FALSE)
scvals <- muvals/gamma(1+1/shvals)
points(shvals,scvals,pch=3)


###################################################
### chunk number 51: 
###################################################
meanfun = function(shape,scale) { scale*gamma(1+1/shape) }


###################################################
### chunk number 52: 
###################################################
weiblikfun2 <- function(shape,mu) {
  scale <- mu/gamma(1+1/shape)
  -sum(dweibull(time,shape=shape,scale=scale,log=TRUE))
}


###################################################
### chunk number 53: 
###################################################
w2 <- mle2(weiblikfun2,start=list(shape=1,mu=mean(time)))
confint(w2,quietly=TRUE)


###################################################
### chunk number 54: 
###################################################
ci.prof <- confint(w2,quietly=TRUE)["mu",]


###################################################
### chunk number 55:  eval=FALSE
###################################################
## op=par(lwd=2,bty="l",las=1,cex=1.5,mgp=c(2.5,1,0))
## source("sqrprofplot.R")
##  sqrprofplot(p2,which=2,sqrVal=TRUE,conf=c(0.9,0.95),
##              col.conf="gray",plot.confstr=TRUE,
##              col.prof="black",col.minval=NA,
##              xlab="Mean",
##               ylab="Deviance")


###################################################
### chunk number 56: 
###################################################
weiblikfun3 <- function(shape,scale,target.mu,penalty=1000) {
  mu <- meanfun(shape,scale)
  NLL = -sum(dweibull(time,shape=shape,scale=scale,log=TRUE))
  pen = penalty*(mu-target.mu)^2
  NLL+pen
}
w3 <- mle2(weiblikfun3,start=list(shape=0.9,scale=13),fixed=list(target.mu=13))
##profile(w3,which="target.mu")
penlikfun3 <- function(p,meanval,penalty=1000) {
  v <-  -sum(dweibull(time,shape=p[1],scale=p[2],log=TRUE))
  pen <- penalty*(meanfun(p[1],p[2])-meanval)^2
  val <- v+pen
## cat(p[1],p[2],penalty,meanval,v,pen,val,"\n")
  return(val)
}
penlikfun2 <- function(m,penalty=1000) {
  o <- optim(fn=penlikfun3,
             par=coef(w1),
             method="BFGS",
             meanval=m,
             penalty=penalty)
  r <- c(o$value,o$par)
  names(r) <- c("NLL","shape","scale")
  r
}
critval <- c(-logLik(w1)+qchisq(0.95,1)/2)
pcritfun <- function(m,penalty=1000) {
   penlikfun2(m,penalty=penalty)[1]-critval
 }
lowx <- c(5,13)
upx <- c(14,30)
penlower <- uniroot(pcritfun,lowx)$root
penupper <- uniroot(pcritfun,upx)$root
ci.pen1 = c(penlower,penupper)
penvec <- 1:5
##v <- seq(5,13,n=101)
##plot(v,sapply(v,pcritfun,penalty=10^2),type="l")
penresults <- matrix(nrow=length(penvec),ncol=3)
for (i in 2:length(penvec)) {
  penlower <- uniroot(pcritfun,lowx,penalty=10^penvec[i])$root
  penupper <- uniroot(pcritfun,upx,penalty=10^penvec[i])$root
  penresults[i,] <- c(penvec[i],penlower,penupper)
}


###################################################
### chunk number 57:  eval=FALSE
###################################################
## shape.deriv <- -shape^2*gamma(1+1/shape)*digamma(1+1/shape)


###################################################
### chunk number 58: 
###################################################
dvar <- deltavar(fun=scale*gamma(1+1/shape),meanval=coef(w1),Sigma=vcov(w1))


###################################################
### chunk number 59: 
###################################################
sdapprox <- sqrt(dvar)
mlmean <- meanfun(coef(w1)["shape"],coef(w1)["scale"])
ci.delta <- mlmean+c(-1.96,1.96)*sdapprox


###################################################
### chunk number 60: 
###################################################
pvars <- diag(vcov(w1))
pcov <- vcov(w1)["shape","scale"]
mlshape <- coef(w1)["shape"]
mlscale <- coef(w1)["scale"]
scderiv <- gamma(1+1/mlshape)
shderiv <-mlscale*-1/mlshape^2*gamma(1+1/mlshape)*digamma(1+1/mlshape)
vapprox1 <- pvars["scale"]*scderiv^2
vapprox2 <- pvars["shape"]*shderiv^2
vapprox3 <- 2*pcov*shderiv*scderiv
vapprox <- vapprox1+vapprox2+vapprox3


###################################################
### chunk number 61: 
###################################################
vmat=mvrnorm(1000,mu=coef(w1),Sigma=vcov(w1))


###################################################
### chunk number 62: 
###################################################
dist = numeric(1000)
for (i in 1:1000) {
  dist[i] = meanfun(vmat[i,1],vmat[i,2])
}
quantile(dist,c(0.025,0.975))


###################################################
### chunk number 63: 
###################################################
ci.ppi <- quantile(dist,c(0.025,0.975))


###################################################
### chunk number 64: 
###################################################
lval <- coef(w1)["scale"]^(-coef(w1)["shape"])
n <- length(time)
inits <- list(list(shape=0.8,lambda=lval),list(shape=0.4,lambda=lval*2),
              list(shape=1.2,lambda=lval/2))
reefgoby.bugs <- bugs(data=list("time","n"),
                      inits,parameters.to.save=c("shape",
                              "scale","lambda","mean"),
                      model.file="reefgobysurv.bug",
                      n.chains=length(inits),n.iter=5000)
reefgoby.coda <- as.mcmc(reefgoby.bugs)
reefgoby.coda <- lump.mcmc.list(reefgoby.coda)
ci.bayes <- HPDinterval(reefgoby.coda)["mean",]
m <- as.matrix(reefgoby.coda)[,"mean"]
m <- reefgoby.coda[,"mean",drop=FALSE]
m <- as.matrix(reefgoby.coda)[,"mean"]
ci.bayesq <- quantile(m,c(0.025,0.975))


###################################################
### chunk number 65: 
###################################################
op=par(lwd=2,bty="l",las=1,cex=1.5,mgp=c(3,1,0))
##par(mar=c(1,4,2,2))
vals = mvrnorm(1000,mu=coef(w1),Sigma=vcov(w1))
dist = apply(vals,1,function(x)meanfun(x[1],x[2]))
textsize = 0.5
ang = 15
par(yaxs="i")
## hist(dist,breaks=30,col="gray",xlab="Mean survival time",main="",freq=FALSE,
##      axes=FALSE,xlim=c(5,35))
## axis(side=2,at=c(0,0.05,0.1),yaxs="i")
## axis(side=1)
##box()
##lines(density(dist))
conf.ppi <- quantile(dist,c(0.025,0.975))
##abline(v=conf.ppi,lty=2,lwd=2)
##par(mar=c(5,4,1,2))
par(mar=c(4,4,2,2))
d1 = density(dist)
d2 = density(lump.mcmc.list(as.mcmc(reefgoby.bugs)[,"mean"]))
plot(d1,ylim=c(0,0.12),main="",
     xlab="",axes=FALSE,xlim=c(0,30))
axis(side=2,at=c(0,0.05,0.1,0.15))
lines(d2,col="gray")
axis(side=1)
box()
mtext(side=1,line=2.5,at=15,"Mean survival time",cex=1.5)
segments(rep(conf.ppi,2),
         rep(0,4),
         rep(conf.ppi,2),
         approx(d1$x,d1$y,xout=conf.ppi)$y)
mpos = 13
labsize = 3
ht1 <- 0.006
arrows(c(mpos-labsize,mpos+labsize),
       rep(ht1,2),
       conf.ppi,
       rep(ht1,2),ang=ang)
text(mpos,ht1,"PP interval",cex=textsize)
segments(rep(ci.bayes,2),
         rep(0,4),
         rep(ci.bayes,2),
         approx(d2$x,d2$y,xout=ci.bayes)$y,
         col="gray")
ht2 <- 0.012
mpos2 <- 15
labsize2 <- 3.5
arrows(c(mpos2-labsize2,mpos2+labsize2),
       rep(ht2,2),
       ci.bayes,
       rep(ht2,2),col="gray",ang=ang)
text(mpos2,ht2,"Bayes credible",col="darkgray",cex=textsize)
## segments(rep(ci.bayesq,2),
##          rep(0,4),
##          rep(ci.bayesq,2),
##          approx(d2$x,d2$y,xout=ci.bayesq)$y,
##          col="gray",lty=3)
## ht3 <- 0.012
## mpos3 <- 16
## labsize3 <- 3.5
## arrows(c(mpos3-labsize3,mpos3+labsize3),
##        rep(ht3,2),
##        ci.bayesq,
##        rep(ht3,2),col="gray",ang=ang)
## text(mpos3,ht3,"Bayes quantile",col="darkgray",cex=textsize)
par(op)


###################################################
### chunk number 66:  eval=FALSE
###################################################
## lval <- coef(w1)["scale"]^(-coef(w1)["shape"])
## n <- length(time)
## inits <- list(list(shape=0.8,lambda=lval),list(shape=0.4,lambda=lval*2),
##               list(shape=1.2,lambda=lval/2))


###################################################
### chunk number 67:  eval=FALSE
###################################################
## reefgoby.bugs <- bugs(data=list("time","n"),
##                       inits,parameters.to.save=c("shape","scale","lambda","mean"),
##                       model.file="reefgobysurv.bug",
##                       n.chains=length(inits),n.iter=5000)


###################################################
### chunk number 68: 
###################################################
citab = rbind(ci.prof,ci.pen1,ci.delta,ci.ppi,ci.bayes)#,ci.bayesq)
cinames = c("exact profile","profile:penalty",
  "delta method","PPI",
  "Bayes credible")
#,"Bayes quantile")
dimnames(citab)=list(cinames,c("lower","upper"))


###################################################
### chunk number 69: 
###################################################
latex(round(citab,3),file="",table.env=FALSE,title="method")


###################################################
### chunk number 70:  eval=FALSE
###################################################
## tabpaste <- function(x,sep=" & ",eol="\\") {
##   x2 <- rbind(c("",colnames(x)),cbind(rownames(x),
##                                      matrix(as.character(x),nrow=nrow(x))))
##   x3 <- apply(x2,1,paste,collapse=sep)
##   x4 <- paste(x3,collapse=eol)
##   return(x4)
## }


###################################################
### chunk number 71: 
###################################################
x=c(0,0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,2,3,4,5)


###################################################
### chunk number 72: 
###################################################
fit.nb=fitdistr(x,"negative binomial")
lnb = c(-logLik(fit.nb))
fit.pois=fitdistr(x,"Poisson")
lpois = c(-logLik(fit.pois))
dev = c(-2*(lnb-lpois))
pval = pchibarsq(dev,1,lower.tail=FALSE)


###################################################
### chunk number 73: 
###################################################
devdiff = 2*(logLik(fit.nb)-logLik(fit.pois))
pchisq(devdiff,df=1,lower.tail=FALSE)


###################################################
### chunk number 74: 
###################################################
pval=pchisq(2*(logLik(fit.nb)-logLik(fit.pois)),df=1,lower.tail=FALSE)


###################################################
### chunk number 75: 
###################################################
simulated.dev = function() {
  simx = rpois(length(x),lambda=mean(x))
  simfitnb = try(fitdistr(simx,"negative binomial"))
  if (inherits(simfitnb,"try-error")) return(NA)
  simfitpois = fitdistr(simx,"Poisson")
  dev=c(2*(logLik(simfitnb)-logLik(simfitpois)))
}


###################################################
### chunk number 76: 
###################################################
set.seed(1001)
devdist = replicate(3000,simulated.dev())
devdist = na.omit(devdist)
nreps = length(devdist)


###################################################
### chunk number 77: 
###################################################
obs.dev=2*(logLik(fit.nb)-logLik(fit.pois))
sum(devdist>=obs.dev)/nreps


###################################################
### chunk number 78: 
###################################################
weiblikfun3 <- function(shape,scale,target.mu,penalty=1000) {
  mu <- meanfun(shape,scale)
  NLL = -sum(dweibull(time,shape=shape,scale=scale,log=TRUE))
  pen = penalty*(mu-target.mu)^2
  NLL+pen
}
w3 <- mle2(weiblikfun3,start=list(shape=0.9,scale=13),fixed=list(target.mu=13))


###################################################
### chunk number 79:  eval=FALSE
###################################################
##   if (shape>0) {
##     NLL = -sum(dweibull(time,shape=shape,scale=scale,log=TRUE))
##     pen = 0
##   } else {
##     NLL = -sum(dweibull(time,shape=0.0001,scale=scale,log=TRUE))
##     pen = penalty*shape^2
##   }
##   NLL+pen


###################################################
### chunk number 80: 
###################################################
critval <- -logLik(w1)+qchisq(0.95,1)/2


###################################################
### chunk number 81: 
###################################################
pcritfun <- function(target.mu,penalty=1000) {
  mfit <- mle2(weiblikfun3,
              start=list(shape=0.85,scale=12.4),
              fixed=list(target.mu=target.mu),
              data=list(penalty=penalty))
  lval <- -logLik(mfit)
  lval - critval
}


###################################################
### chunk number 82: 
###################################################
lowx <- c(5,13)
penlower <- uniroot(pcritfun,lowx)$root


###################################################
### chunk number 83: 
###################################################
upx <- c(14,30)
penupper <- uniroot(pcritfun,upx)$root


###################################################
### chunk number 84: 
###################################################
uniroot(pcritfun,lowx,penalty=1e6)$root


###################################################
### chunk number 85: 
###################################################
shape <- coef(w1)["shape"]
scale <- coef(w1)["scale"]
numericDeriv(quote(scale*gamma(1+1/shape)),c("scale","shape"))


###################################################
### chunk number 86: 
###################################################
dshape = 0.0001
x2 = scale*gamma(1+1/(shape+dshape))
x1 = scale*gamma(1+1/shape)
(x2-x1)/dshape


###################################################
### chunk number 87: 
###################################################
reefgoby.coda <- as.mcmc(reefgoby.bugs)
reefgoby.coda <- lump.mcmc.list(reefgoby.coda)
ci.bayes <- HPDinterval(reefgoby.coda)["mean",]
##ci.bayesq <- quantile(reefgoby.coda,c(0.025,0.975))


