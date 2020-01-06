###################################################
### chunk number 1: 
###################################################
## source("plotCI.R") ## ??
library(plotrix)
library(emdbook) ## for dzinbinom and rzinbinom
library(bbmle)
library(MASS)
source("chapskel.R")


###################################################
### chunk number 2: 
###################################################
op=par(mfrow=c(1,2))
par(lwd=2,bty="l",las=1,cex=1.5)
x = 1:20
a =2; b=1
y_det = a+b*x
y = y_det+rnorm(20,sd=2)
plot(x,y)
abline(lm(y~x),lty=2)
abline(a,b)
legend("topleft",c("true","best fit"),
       lty=1:2,bty="n")
##
b = 1
a = 20
k = 5
n=25
xmax=5
axmax=5
x= runif(n,min=0,max=xmax)
y_det= a*b/(b+x)
y = rnbinom(n,mu=y_det,size=k)
plot(x,y,xlim=c(0,axmax),axes=FALSE)
axis(side=1,at=0:5)
axis(side=2)
box()
curve(a*b/(b+x),add=TRUE,lty=1,to=xmax)
#mfun = function(a,b,k) {
#  -sum(dnbinom(y,mu=a/(b+x),size=k,log=TRUE))
#}
m1 = mle2(y~dnbinom(mu=a*b/(b+x),size=k),
  start=list(a=20,b=1,k=1))
curve(coef(m1)["a"]*coef(m1)["b"]/(coef(m1)["b"]+x),add=TRUE,lty=2,to=xmax)
legend("topright",
       c("true","best fit"),
       lty=1:2,
       bty="n")
par(op)


###################################################
### chunk number 3: 
###################################################
x = 1:20
a =2; b=1


###################################################
### chunk number 4: 
###################################################
y_det = a+b*x


###################################################
### chunk number 5: 
###################################################
y=rnorm(20,mean=y_det,sd=2)


###################################################
### chunk number 6: 
###################################################
a = 20
b = 1
k = 5


###################################################
### chunk number 7: 
###################################################
x= runif(50,min=0,max=5)


###################################################
### chunk number 8: 
###################################################
y_det= a*b/(b+x)
y = rnbinom(50,mu=y_det,size=k)


###################################################
### chunk number 9: 
###################################################
g = factor(rep(1:2,each=25))


###################################################
### chunk number 10: 
###################################################
a=c(20,10)
b=c(1,2)
k=5


###################################################
### chunk number 11: 
###################################################
y_det= a[g]/(b[g]+x)
y = rnbinom(50,mu=y_det,size=k)


###################################################
### chunk number 12: 
###################################################
##
## b = c(1,2)
## a = c(20,10)
## k = 5
## g = rep(1:2,each=25)
## x= runif(50,min=0,max=5)
## y_det= a[g]*b/(b+x)
## y = rnbinom(50,mu=y_det,size=k)
m2 = mle2(y~dnbinom(mu=a*b/(b+x),size=k),
  parameters=list(a~g,b~g),
  start=list(a=20,b=1,k=1))
axmax=8
xmax=5
op=par(lwd=2,las=1,mgp=c(3,1,0),bty="l",cex=1.5)
plot(x,y,pch=as.numeric(g),xlim=c(0,axmax),axes=FALSE,lwd=1)
axis(side=1,at=0:5)
axis(side=2)
box()
curve(a[1]/(b[1]+x),add=TRUE,col="gray",to=xmax)
curve(a[2]/(b[1]+x),add=TRUE,lty=2,col="gray",to=xmax)
c2 = coef(m2)
curve(c2[1]*c2[3]/(c2[3]+x),add=TRUE,lty=1,to=xmax)
curve((c2[1]+c2[2])*(c2[3]+c2[4])/(c2[3]+c2[4]+x),add=TRUE,lty=2,to=xmax)
legend("topright",
       outer(c("data","true","best fit"),
             c("(sp. 1)","(sp. 2)"),paste),
       lty=c(NA,1,1,NA,2,2),
       lwd=c(1,2,2,1,2,2),
       col=rep(c("black","gray","black"),2),
       pch=c(1,NA,NA,2,NA,NA),cex=0.8,bty="n")
par(op)


###################################################
### chunk number 13: 
###################################################
set.seed(1001)
N = 603
a = 0.696
b = 9.79
mu=25.32
zprob=0.123
k=0.932
recrprob = function(S) { a/(1+(a/b)*S) }
settlers = rzinbinom(N,mu=mu,size=k,zprob=zprob)
recr = rbinom(N,prob=recrprob(settlers),size=settlers)


###################################################
### chunk number 14: 
###################################################
op = par(mfrow=c(1,2),mar=c(5,4,2,0.2),lwd=2,las=1,cex=1.5,bty="l")
hist(settlers,breaks=40,col="gray",ylab="Frequency",xlab="Settlers",main="")
plot(settlers,recr,xlab="Settlers",ylab="Recruits")
curve(a*x/(1+(a/b)*x),add=TRUE)
par(op)


###################################################
### chunk number 15: 
###################################################
N = 603
a = 0.696
b = 9.79
mu=25.32
zprob=0.123
k=0.932


###################################################
### chunk number 16: 
###################################################
recrprob = function(S) { a/(1+(a/b)*S) }


###################################################
### chunk number 17: 
###################################################
settlers = rzinbinom(N,mu=mu,size=k,zprob=zprob)
recr = rbinom(N,prob=recrprob(settlers),size=settlers)


###################################################
### chunk number 18: 
###################################################
set.seed(1001)
L = 30
nparents = 50
offspr_per_parent = 10
noffspr = nparents*offspr_per_parent
dispdist = 2
parent_x = runif(nparents,min=0,max=L)
parent_y = runif(nparents,min=0,max=L)
angle = runif(noffspr,min=0,max=2*pi)
dist = rexp(noffspr,1/dispdist)
offspr_x = rep(parent_x,each=offspr_per_parent)+cos(angle)*dist
offspr_y = rep(parent_y,each=offspr_per_parent)+sin(angle)*dist
pos <- cbind(offspr_x,offspr_y)
ndist <- as.matrix(dist(pos,upper=TRUE,diag=TRUE))
nbrcrowd = apply(ndist<2,1,sum)-1
ci = nbrcrowd*3
M=2.3
alpha=0.49
mass_det=M/(1+ci)
mass = rgamma(length(mass_det),scale=mass_det,shape=alpha)
b = 271.6; k= 0.569
seed_det = b*mass
seed = rnbinom(length(seed_det),mu=seed_det,size=k)
f= fitdistr(nbrcrowd,"negative binomial")


###################################################
### chunk number 19: 
###################################################
op=par(mfrow=c(2,2),mar=c(5,5,0.2,0.2),lwd=1,las=1,mgp=c(3.5,1,0))
## plot 1
plot(offspr_x,offspr_y,xlab="",ylab="")
## plot 2
b1 = barplot(table(factor(nbrcrowd,levels=0:max(nbrcrowd)))/length(nbrcrowd),
xlab="Number of neighbors",ylab="Proportion")
points(dnbinom(0:max(nbrcrowd),size=f$estimate["size"],
mu=f$estimate["mu"]),pch=16,type="b")
## plot 3
plot(ci,mass,cex=0.5,log="y",xlab="Competition index",ylab="Biomass (g)")
curve(M/(1+x)*alpha,add=TRUE,from=0)
## plot 4
## x = runif(76,min=0,max=5)
seed_det = b*mass
seed = rnbinom(length(seed_det),mu=seed_det,size=k)
par(mgp=c(2.5,1,0))
plot(mass,1+seed,log="xy",xlab="Mass",ylab="1+Seed set")
curve(b*x+1,add=TRUE)
##abline(a=0,b=b,lty=1)
curve(qnbinom(0.975,mu=b*x,size=k)+1,lty=2,from=0.001,add=TRUE,type="s")
curve(qnbinom(0.025,mu=b*x,size=k)+1,lty=2,from=0.001,add=TRUE,type="s")
par(op)


###################################################
### chunk number 20: 
###################################################
set.seed(1001)
L = 30
nparents = 50
offspr_per_parent = 10
noffspr = nparents*offspr_per_parent
dispdist = 2


###################################################
### chunk number 21: 
###################################################
parent_x = runif(nparents,min=0,max=L)
parent_y = runif(nparents,min=0,max=L)


###################################################
### chunk number 22: 
###################################################
angle = runif(noffspr,min=0,max=2*pi)
dist = rexp(noffspr,1/dispdist)


###################################################
### chunk number 23: 
###################################################
offspr_x = rep(parent_x,each=offspr_per_parent)+cos(angle)*dist
offspr_y = rep(parent_y,each=offspr_per_parent)+sin(angle)*dist


###################################################
### chunk number 24: 
###################################################
pos <- cbind(offspr_x,offspr_y)
ndist <- as.matrix(dist(pos,upper=TRUE,diag=TRUE))
nbrcrowd = apply(ndist<2,1,sum)-1


###################################################
### chunk number 25: 
###################################################
ci = nbrcrowd*3
M=2.3
alpha=0.49
mass_det=M/(1+ci)
mass = rgamma(length(mass_det),scale=mass_det,shape=alpha)


###################################################
### chunk number 26: 
###################################################
b = 271.6; k= 0.569
seed_det = b*mass
seed = rnbinom(length(seed_det),mu=seed_det,size=k)


###################################################
### chunk number 27: 
###################################################
f= fitdistr(nbrcrowd,"negative binomial")


###################################################
### chunk number 28:  eval=FALSE
###################################################
## ## triangular plot
## a=2
## b=0.5
## x = runif(100,0,4)
## y = runif(100,max=a*exp(-b*x))
## plot(x,y)
## curve(a*exp(-b*x),lty=2,add=TRUE)


###################################################
### chunk number 29: 
###################################################
powplot <- function(m0=45.5,m1=54,sd=1.745,len=100,col1="gray",col2=gray(0.2)) {
  rng <- range(c(m0-4*sd,m1-4*sd,m0+4*sd,m1+4*sd))
  curve(dnorm(x,m0,sd),from=rng[1],to=rng[2],yaxs="i",
        ylim=c(0,dnorm(m0,m0,sd)*1.1),axes=FALSE,ylab="",xlab="")
  axis(side=1)
  box()
  ##  plot(pvec1,dnorm(pvec1,m0,sd),type="l",xlim=rng)
  ## lines(pvec2,dnorm(pvec2,m1,sd))
  qval <- qnorm(c(0.025,0.975),m0,sd)
  pvec <- sort(c(qval,seq(rng[1],rng[2],length=len)))
  abline(v=qval,lty=2)
  pvec.up <- pvec[pvec>=qval[2]]
  pfun <- function(x,y,...)
    polygon(c(x,rev(x)),c(rep(0,length(x)),rev(y)),...)
  pfun(pvec.up,dnorm(pvec.up,m1,sd),col=col1)
  pfun(pvec.up,dnorm(pvec.up,m0,sd),col=col2)
  pvec.down <- pvec[pvec<=qval[1]]
  pfun(pvec.down,dnorm(pvec.down,m0,sd),col=col2)
  pfun(pvec.down,dnorm(pvec.down,m1,sd),col=col1)
  ##polygon(c(pvec1b,rev(pvec1b)),
  ##          c(rep(0,length(pvec1b)),dnorm(rev(pvec1b),m1,sd)),col="gray")
  power=pnorm(qval[2],m1,sd,lower=FALSE)+pnorm(qval[1],m1,sd)
  curve(dnorm(x,m1,sd),from=rng[1],to=rng[2],add=TRUE,lwd=2)
  invisible(power)
}

powval <- function(m0=45.5,m1=54,sd=1.745,alpha=0.95) {
  1-(pnorm(qnorm(1-(1-alpha)/2,m0,sd),m1,sd)-
      pnorm(qnorm((1-alpha)/2,m0,sd),m1,sd))
}


###################################################
### chunk number 30: 
###################################################
op <- par(mfrow=c(1,2))
par(lwd=2,las=1,bty="l",cex=1.5,mgp=c(2.5,1,0))
par(mar=c(4,4,2,1))
powplot(m0=0,m1=2,sd=0.75)
mtext(side=2,las=0,"Probability",cex=1.5,line=1)
mtext(side=1,las=0,"x",cex=1.5,line=2.5)
text(c(0,2),c(0.56,0.56),
     c(expression(H[0]),expression(H[1])))
text(-2.6,0.1,expression(alpha))
arrows(c(-2.5,-2.2),
       c(0.09,0.09),
       c(-1.95,1.4),
       c(0.026,0.026),
       angle=15)
text(4,0.4,"Power")
text(4,0.37,expression(1-beta))
arrows(4,0.36,2,0.3,angle=15)
##xvec <- seq(35,55,len=100)
par(mar=c(4,4,2,2))
curve(powval(x,m0=0,sd=0.25),from=-5,to=5,ylab="Power",xlab="Effect size")
##plot(xvec,sapply(xvec,function(x)powval(m1=x)),type="l")
abline(h=0.05)
curve(powval(x,m0=0,sd=0.75),lty=2,add=TRUE)
curve(powval(x,m0=0,sd=2),lty=3,add=TRUE)
points(2,powval(2,m0=0,sd=0.75),pch=16)
par(xpd=NA,adj=0)
text(c(1.5,2.5,3.5),
     c(1.03,0.88,0.36),
     c(expression(sigma==0.25),
       expression(sigma==0.75),
       expression(sigma==2)),
     cex=1,xpd=NA)
par(op)


###################################################
### chunk number 31: 
###################################################
x = 1:20
a =2; b=1; sd=8
N = 20


###################################################
### chunk number 32: 
###################################################
y_det = a+b*x
y = rnorm(N,mean=y_det,sd=sd)
m = lm(y~x)
coef(summary(m))["x","Pr(>|t|)"]


###################################################
### chunk number 33: 
###################################################
nsim = 400


###################################################
### chunk number 34: 
###################################################
pval = numeric(nsim)


###################################################
### chunk number 35: linregpow1
###################################################
for (i in 1:nsim) {
  y_det = a+b*x
  y = rnorm(N,mean=y_det,sd=sd)
  m = lm(y~x)
  pval[i] = coef(summary(m))["x","Pr(>|t|)"]
}


###################################################
### chunk number 36: 
###################################################
sum(pval<0.05)/nsim


###################################################
### chunk number 37: linregpow2
###################################################
bvec = seq(-2,2,by=0.1)
power.b = numeric(length(bvec))
for (j in 1:length(bvec)) {
  b = bvec[j]
  for (i in 1:nsim) {
    y_det = a+b*x
    y = rnorm(N,mean=y_det,sd=sd)
    m = lm(y~x)
    pval[i] = coef(summary(m))["x","Pr(>|t|)"]
  }
  power.b[j] = sum(pval<0.05)/nsim
}


###################################################
### chunk number 38: 
###################################################
load("negbinhyp-batch2.RData")
nsim <- 200


###################################################
### chunk number 39:  eval=FALSE
###################################################
## m0 = mle2(y~dnbinom(mu=a*b*x/(b+x),size=k),start=list(a=15,b=1,k=5))
## m1 = mle2(y~dnbinom(mu=a*b*x/(b+x),size=k),
##   parameters=list(a~g,b~g),
##   start=list(a=15,b=1,k=5))
## anova(m0,m1)[2,"Pr(>Chisq)"]


###################################################
### chunk number 40: 
###################################################
op <- par(lwd=2,las=1,bty="l",cex=1.5,mgp=c(2.5,1,0))
plot(Nvec[power.N>0],power.N[power.N>0],log="x",type="b",
     xlab="Sample size",
     ylab="Power",ylim=c(0,1),axes=FALSE)
axis(side=2)
axis(side=1,cex.axis=0.8)
abline(h=1,lty=2)
box()
#grid()
par(op)


###################################################
### chunk number 41: 
###################################################
set.seed(1001)


###################################################
### chunk number 42: 
###################################################
x = rnbinom(100,mu=1,size=0.5)
f = fitdistr(x,"negative binomial")
f


###################################################
### chunk number 43: 
###################################################
Nvec = round(lseq(20,500,length=100))
estk = numeric(length(Nvec))
estksd = numeric(length(Nvec))


###################################################
### chunk number 44: negbinpow1
###################################################
set.seed(1001)
for (i in 1:length(Nvec)) {
  N = Nvec[i]
  x = rnbinom(N,mu=1,size=0.5)
  f = fitdistr(x,"negative binomial")
  estk[i] = f$estimate["size"]
  estksd[i] = f$sd["size"]
}


###################################################
### chunk number 45: 
###################################################
op <- par(mfrow=c(1,2))
par(lwd=2,las=1,bty="l",cex=1.5,mgp=c(2.5,1,0))
plot(Nvec,estk,log="x",xlab="Sample size",ylab="Estimated k",axes=FALSE)
axis(side=2)
axis(side=1,cex.axis=0.7)
box()
lines(Nvec,predict(loess(estk~Nvec)))
abline(h=0.5,lty=2)
plot(Nvec,estksd,log="xy",xlab="Sample size",ylab="Estimated sd(k)",
     mgp=c(3,1,0),axes=FALSE)
axis(side=2)
axis(side=1,cex.axis=0.7)
box()


###################################################
### chunk number 46: 
###################################################
mm = nls(recr~a*settlers/(1+(a/b)*settlers),start=c(a=0.696,b=9.79))
shep = nls(recr~a*settlers/(1+(a/b)*settlers^d),start=c(coef(mm),d=1))

recrprob = function(x,a,b,d) a/(1+(a/b)*x^d)

simdata = function(n,d=1,a=0.696,b=9.79) {
  scoefs = c(mu=25.32,k=0.932,zprob=0.123)
  settlers = rzinbinom(n,mu=scoefs["mu"],size=scoefs["k"],zprob=scoefs["zprob"])
  recr = rbinom(n,prob=recrprob(settlers,a,b,d),size=settlers)
  data.frame(settlers=settlers,recr=recr)
}

getvals = function(n=100,d=1) {
  OK = FALSE
  while (!OK) {
    z = simdata(n,d)
    mm = try(nls(recr~a*settlers/(1+(a/b)*settlers),start=c(a=0.696,b=9.79),data=z))
    shep = try(nls(recr~a*settlers/(1+(a/b)*settlers^d),start=c(coef(mm),d=1),data=z))
    OK = (class(shep)!="try-error" && class(mm)!="try-error")
    if (OK) {
      mm.ci = try(confint(mm))
      shep.ci = try(confint(shep))
      if (class(shep.ci)=="try-error") {
        s0 = summary(shep)
        ci_width = qt(0.975,s0$df[2])*s0$parameters[,"Std. Error"]
        shep.ci=c(s0$parameters[,"Estimate"]-ci_width,
          s0$parameters[,"Estimate"]+ci_width)
      }
    }
  }
  res = c(coef(mm),mm.ci,coef(shep),shep.ci)
  names(res) = c("a.mm","b.mm","a.mm.lo","b.mm.lo","a.mm.hi","b.mm.hi",
         "a.shep","b.shep","d.shep","a.shep.lo","b.shep.lo","d.shep.lo",
         "a.shep.hi","b.shep.hi","d.shep.hi")
  res
}


###################################################
### chunk number 47:  eval=FALSE
###################################################
## set.seed(1001)
## nvec = c(50,100,200,300,400,500,1000,2000)
## nsim=100
## dvec = seq(0.7,1.3,by=0.1)
## r1= getvals()
## resmat = matrix(ncol=length(r1)+3,nrow=nsim*length(nvec)*length(dvec))
## c0 = 1
## for (i in 1:length(dvec)) {
##   for (j in 1:length(nvec)) {
##      for (k in 1:nsim) {
##        cat(i,j,k,"\n")
##        resmat[c0,] = c(dvec[i],nvec[j],k,getvals(n=nvec[j],d=dvec[i]))
##      }
##    }
## }     


###################################################
### chunk number 48: 
###################################################
load("chap5-batch2.RData")
d.null = 1
true.b=9.79
## nsim=400  ## for batch2
faclist <- list(factor(resmat[,"d"]),factor(resmat[,"n"]))
a.mm.mean =   tapply(resmat[,"a.mm"],faclist,mean)
b.mm.mean =   tapply(resmat[,"b.mm"],faclist,mean)
b.mm.sd =     tapply(resmat[,"b.mm"],faclist,sd)
a.shep.mean = tapply(resmat[,"a.shep"],faclist,mean)
b.shep.mean = tapply(resmat[,"b.shep"],faclist,mean)
b.shep.sd =   tapply(resmat[,"b.shep"],faclist,sd)
d.shep.mean = tapply(resmat[,"d.shep"],faclist,mean)
d.shep.sd =   tapply(resmat[,"d.shep"],faclist,sd)
rshepOK <-    resmat[resmat[,"shep.ci.OK"]==1,]
##rshepOK <-    resmat
d.shep.OK <-  table(factor(rshepOK[,"d"]),factor(rshepOK[,"n"]))
faclist2 <- list(factor(rshepOK[,"d"]),factor(rshepOK[,"n"]))
d.shep.pow =  with(as.data.frame(rshepOK),
  tapply(d.shep.hi<d.null |  d.shep.lo>d.null,faclist2,sum))/d.shep.OK
d.shep.cov = with(as.data.frame(rshepOK),
  tapply(d.shep.hi>d & d.shep.lo<d,faclist2,sum))/d.shep.OK
d.shep.cwid = tapply(resmat[,"d.shep.hi"]-resmat[,"d.shep.lo"],
  faclist,mean)


###################################################
### chunk number 49: 
###################################################
nsim = 20
trueval=1.2
nullval=1
r1 = resmat[resmat[,"n"]==1000 & resmat[,"d"]==trueval & resmat[,"shep.ci.OK"],]
totsim=nrow(r1)
dm = mean(r1[,"d.shep"])
bias = dm-trueval
dsd = sd(r1[,"d.shep"])
dcover = sum(r1[,"d.shep.lo"]<trueval &  r1[,"d.shep.hi"]>trueval)/totsim
dpow = sum(r1[,"d.shep.lo"]>nullval |  r1[,"d.shep.hi"]<nullval)/totsim
ddens = density(r1[,"d.shep"],n=700,from=0.9,to=1.6)
ddens.l = density(r1[,"d.shep.lo"],n=700,from=0.9,to=1.6)
##ddens.u = density(r1[,"d.shep.hi"])


###################################################
### chunk number 50: 
###################################################
r=require(plotrix,quietly=TRUE,warn.conflicts=FALSE)
plotCI(1:nsim,r1[1:nsim,"d.shep"],li=r1[1:nsim,"d.shep.lo"],ui=r1[1:nsim,"d.shep.hi"],xlab="Simulation",
       ylab=expression(hat(d)),xlim=c(1,40),axes=FALSE,ylim=c(0.9,1.6))
abline(h=trueval,lwd=2)
abline(h=nullval,lwd=2,lty=2)
axis(side=2)
axis(side=1,at=c(1,10,20))
box()
axis.break(breakpos=25,axis=1,style="zigzag",brw=0.05)
n = length(ddens$x)
p1 = 35
sc=-0.8
polygon(c(p1+ddens$y*sc,p1),c(ddens$x,min(ddens$x)),col="gray")
segments(p1+ddens$y[1:(n-1)]*sc,ddens$x[1:(n-1)],
         p1+ddens$y[2:n]*sc,ddens$x[2:n])
sc2=+0.8
#polygon(c(p1+ddens.l$y*sc2,p1),c(ddens.l$x,min(ddens.l$x)),col="gray")
covvals=which(ddens.l$x<trueval)
powvals=which(ddens.l$x>nullval)
npowvals=which(ddens.l$x<nullval)
polygon(c(p1,p1+ddens.l$y[powvals]*sc2,p1),c(nullval,ddens.l$x[powvals],max(ddens.l$x)),col="lightgray")
polygon(c(p1,p1+ddens.l$y[covvals]*sc2,p1),c(min(ddens.l$x),ddens.l$x[covvals],trueval),col="darkgray")
polygon(c(p1,p1+ddens.l$y[npowvals]*sc2,p1),c(min(ddens.l$x),ddens.l$x[npowvals],nullval),col=gray(0.2))
segments(p1+ddens.l$y[1:(n-1)]*sc2,ddens.l$x[1:(n-1)],
         p1+ddens.l$y[2:n]*sc2,ddens.l$x[2:n])
rect(23,0.95,28,1.6,col="white",border=NA)
text(33,1.58,"estimates",adj=1)
text(36,1.58,"lower\nbounds",adj=0)
text(22,1.23,expression(paste("mean: ",E,"[",hat(d),"]")),adj=0)
text(22,1.18,expression(paste("true: ",d)),adj=0)
text(22,0.98,expression(paste("null: ",d[0])),adj=0)
text(rep(36,3),c(1.23,1.1,0.97),c("a","b","c"),cex=1.5,adj=0)
abline(h=dm,lty=2)
arrows(32.25,dm,32.25,dm-dsd,angle=90)
text(33,dm-dsd/2,expression(sigma[hat(d)]))


###################################################
### chunk number 51: 
###################################################
clabyoff = -0.1
op=par(mfrow=c(2,2),mar=c(4.2,4.1,2,1),mgp=c(2.5,1,0),las=1,cex=1.5,
  bty="l",lwd=2)
##mtext("a",side=2,outer=FALSE,col=2,at=0.5)
##mtext("a",side=2,col=4,las=0,cex=2,las=1,adj=0,padj=1.1)
## bias
matplot(nvec,t(d.shep.mean),type="b",col=1,pch="7890123",
       xlab="Sample size",ylab="Estimated d",axes=FALSE)
axis(side=1,cex.axis=0.8)
axis(side=2,cex.axis=0.8)
box()
par(xpd=NA)
corner.label("a",yoff=clabyoff,xoff=0,cex=2)
par(xpd=FALSE)
#plotCI(nvec[col(d.shep.sd)],d.shep.mean,d.shep.sd/5,add=TRUE)
abline(h=dvec,lty=1:length(dvec),lwd=2,col="gray")
##
##
matplot(nvec,t(d.shep.cwid),type="b",col=1,pch="7890123",
       xlab="Sample size",ylab="Confidence interval width",axes=FALSE)
axis(side=1,cex.axis=0.8)
axis(side=2,cex.axis=0.8)
box()
par(xpd=NA)
corner.label("b",yoff=clabyoff,xoff=0,cex=2)
par(xpd=FALSE)

matplot(nvec,t(d.shep.cov),type="b",pch="7890123",col=1,
        xlab="Sample size",ylab="Coverage",mgp=c(3,1,0),ylim=c(0.75,1.0),
        cex.axis=0.8)
par(xpd=NA)
corner.label("c",yoff=clabyoff,xoff=0,cex=2)
par(xpd=FALSE)
abline(h=0.95,lwd=2)
matplot(nvec,t(d.shep.pow),type="b",pch="7890123",col=1,
        ylab=expression(paste("Power or ",alpha)),
        xlab="Sample size",ylim=c(0,1),cex.axis=0.8)
par(xpd=NA)
corner.label("d",yoff=clabyoff,xoff=0,cex=2)
par(xpd=FALSE)
abline(h=0.05,lty=2)
abline(h=0.8,lty=2)
par(op)


###################################################
### chunk number 52: 
###################################################
matplot(nvec,t(b.mm.mean[3:5,]),pch="901",
        ylab="Estimate of b",xlab="Number of samples",col="darkgray",
        type="b",ylim=c(5,20))
plotCI(rep(nvec,3),t(b.mm.mean[3:5,]),t(b.mm.sd[3:5,]),add=TRUE)
matpoints(nvec,t(b.shep.mean[3:5,]),pch="901",col=1,type="b")
##
##matplot(nvec,t(b.shep.mean[3:5,]),pch="901",col=1,type="b")
tmpsd = t(b.shep.sd[3:5,])
tmpsd[tmpsd>5] = NA
##plotCI(rep(nvec,3),t(b.shep.mean[3:5,]),tmpsd,add=TRUE)
abline(h=true.b,lwd=2)
text(c(1160,1168,471,1232),
     c(6.7,9.5,12.5,15.4),
     c("B-H est.: d=1.1",
       "B-H est.: d=1.0",
       "Shepherd est.",
       "B-H est.: d=0.9"),adj=0)


###################################################
### chunk number 53: 
###################################################
## obsolete/not used
powfun = function(mu,sigma,n,alpha=0.05,two.tail=FALSE) {
  if (!two.tail) {
    critlevel = qnorm(1-alpha,mean=0,sd=sigma/sqrt(n))
    pnorm(critlevel,mean=mu,sd=sigma/sqrt(n),lower.tail=FALSE)
  } else {
    upper = qnorm(1-alpha/2,mean=0,sd=sigma/sqrt(n))
    lower = qnorm(alpha/2,mean=0,sd=sigma/sqrt(n))
    pnorm(upper,mean=mu,sd=sigma/sqrt(n),lower.tail=FALSE)+
      pnorm(lower,mean=mu,sd=sigma/sqrt(n),lower.tail=TRUE)
  }
}


###################################################
### chunk number 54:  eval=FALSE
###################################################
## xvec = seq(-5,10,by=0.1)
## m=2.75
## alpha=0.05
## upperq=1-alpha
## xvec.low = c(xvec[xvec<upperq],upperq)
## xvec.hi =  c(upperq,xvec[xvec>upperq])
## op=par(mfrow=c(1,2))
## plot(xvec,dnorm(xvec),type="l",lwd=2,xlab="x",ylab="Probability density")
## polygon(c(xvec.low,rev(xvec.low)),c(rep(0,length(xvec.low)),dnorm(rev(xvec.low),mean=m)),col="darkgray",border=NA)
## polygon(c(xvec.hi,rev(xvec.hi)),c(rep(0,length(xvec.hi)),dnorm(rev(xvec.hi),mean=m)),col="lightgray",border=NA)
## polygon(c(xvec.hi,rev(xvec.hi)),c(rep(0,length(xvec.hi)),dnorm(rev(xvec.hi),mean=0)),col=gray(0.3),border=NA)
## curve(dnorm(x),lwd=2,add=TRUE)
## upperq = qnorm(1-alpha)
## abline(v=upperq,lty=2)
## abline(v=0)
## lines(xvec,dnorm(xvec,mean=m))
## text(upperq+0.2,0.02,expression(alpha),col="white",cex=1.5,adj=0)
## text(upperq-0.2,0.02,expression(beta),cex=1.5,adj=1)
## text(m,0.2,expression(1-beta),cex=1.5,adj=0.5)
## ##
## curve(powfun(x,sigma=1,n=1,two.tail=TRUE),from=-8,to=8,lwd=2,
##       ylab=expression(paste("Power ",(1-beta))),xlab=expression(mu))
## par(op)


