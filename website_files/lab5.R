###################################################
### chunk number 1: 
###################################################
x = 1:20
a = 2; b=1


###################################################
### chunk number 2: 
###################################################
set.seed(1001)


###################################################
### chunk number 3: 
###################################################
y_det = a+b*x
y = rnorm(20,mean=y_det,sd=2)


###################################################
### chunk number 4: 
###################################################
plot(x,y)
abline(lm(y~x),lty=2)
abline(a,b)


###################################################
### chunk number 5: 
###################################################
a = 6
b = 1


###################################################
### chunk number 6: 
###################################################
x= runif(50,min=0,max=5)


###################################################
### chunk number 7: 
###################################################
y_det= a/(b+x)
y = rpois(50,y_det)


###################################################
### chunk number 8: 
###################################################
plot(x,y)
curve(a/(b+x),add=TRUE)


###################################################
### chunk number 9:  eval=FALSE
###################################################
## n=100
## x = runif(n,min=0,max=10)
## a=1
## b=0.5
## s=3
## y_det=a*x*exp(-b*x)
## y = rgamma(n,shape=s,scale=y_det/s)
## plot(x,y)
## curve(a*x*exp(-b*x),add=TRUE)


###################################################
### chunk number 10: 
###################################################
set.seed(1001)
nparents = 50
noffspr = 10
L = 30


###################################################
### chunk number 11: 
###################################################
parent_x = runif(nparents,min=0,max=L)
parent_y = runif(nparents,min=0,max=L)


###################################################
### chunk number 12: 
###################################################
angle = runif(nparents*noffspr,min=0,max=2*pi)
dist = rexp(nparents*noffspr,0.5)


###################################################
### chunk number 13: 
###################################################
offspr_x = rep(parent_x,each=noffspr)+cos(angle)*dist
offspr_y = rep(parent_y,each=noffspr)+sin(angle)*dist


###################################################
### chunk number 14: 
###################################################
dist = sqrt((outer(offspr_x,offspr_x,"-"))^2+(outer(offspr_y,offspr_y,"-"))^2)


###################################################
### chunk number 15: 
###################################################
nbrcrowd = apply(dist<2,1,sum)-1


###################################################
### chunk number 16: 
###################################################
plot(offspr_x,offspr_y,xlab="",ylab="")


###################################################
### chunk number 17: 
###################################################
b1 = barplot(table(factor(nbrcrowd,levels=0:max(nbrcrowd)))/length(nbrcrowd),
  xlab="Number of neighbors",ylab="Proportion")


###################################################
### chunk number 18: 
###################################################
ci = nbrcrowd*3


###################################################
### chunk number 19: 
###################################################
M=2.3
alpha=0.49


###################################################
### chunk number 20: 
###################################################
mass_det=M/(1+ci)


###################################################
### chunk number 21: 
###################################################
mass = rgamma(length(mass_det),scale=mass_det,shape=alpha)


###################################################
### chunk number 22: 
###################################################
plot(ci,mass,cex=0.5,xlab="Competition index",ylab="Biomass (g)")
curve(M/(1+x)*alpha,add=TRUE,from=0)


###################################################
### chunk number 23: 
###################################################
b = 271.6; k= 0.569


###################################################
### chunk number 24: 
###################################################
seed_det = b*mass
seed = rnbinom(length(seed_det),mu=seed_det,size=k)


###################################################
### chunk number 25: 
###################################################
plot(mass,1+seed,log="xy",xlab="Mass",ylab="1+Seed set")
curve(b*x+1,add=TRUE)


###################################################
### chunk number 26: 
###################################################
logxvec = seq(-7,0,length=100)
xvec = 10^logxvec
lower = qnbinom(0.025,mu=b*xvec,size=k)
upper = qnbinom(0.975,mu=b*xvec,size=k)
lines(xvec,lower+1,lty=2)
lines(xvec,upper+1,lty=2)


###################################################
### chunk number 27:  eval=FALSE
###################################################
## a = 0.696
## b = 9.79
## recrprob = function(x,a=0.696,b=9.79) a/(1+(a/b)*x)
## scoefs = c(mu=25.32,k=0.932,zprob=0.123)
## settlers = rzinbinom(603,mu=scoefs["mu"],size=scoefs["k"],zprob=scoefs["zprob"])
## rmbbinom = function(n,size,p,theta) {
##   rbinom(n,size=size,prob=rbeta(n,shape1=p*theta,shape2=(1-p)*theta))
## }
## recr = rmbbinom(603,p=recrprob(settlers),theta=10,size=settlers)
## plot(settlers,recr,xlab="Settlers",ylab="Recruits")
## curve(a*x/(1+(a/b)*x),add=TRUE)


###################################################
### chunk number 28: 
###################################################
nt=20
N0=2
dN=1
sd_process=sqrt(2)
sd_obs=sqrt(2)


###################################################
### chunk number 29: 
###################################################
Nobs = numeric(nt)
N = numeric(nt)


###################################################
### chunk number 30: 
###################################################
N[1]=N0
Nobs[1]=N[1]+rnorm(1,sd=sd_obs)


###################################################
### chunk number 31: 
###################################################
for (i in 2:nt) {
  N[i]=N[i-1]+rnorm(1,mean=dN,sd=sd_process)
  Nobs[i]=N[i]+rnorm(1,sd=sd_obs)
}


###################################################
### chunk number 32: 
###################################################
cur_N = N0
Nobs[1] = N[1]+rnorm(1,sd=sd_obs)
for (i in 2:nt) {
  cur_N = cur_N +rnorm(1,mean=dN,sd=sd_process)
  Nobs[i]=cur_N +rnorm(1,sd=sd_obs)
}


###################################################
### chunk number 33: 
###################################################
linsim = function(nt=20,N0=2,dN=1,sd_process=sqrt(2),sd_obs=sqrt(2)) {
  cur_N = N0
  Nobs[1] = N[1]+rnorm(1,sd=sd_obs)
  for (i in 2:nt) {
    cur_N = cur_N +rnorm(1,mean=dN,sd=sd_process)
    Nobs[i]=cur_N +rnorm(1,sd=sd_obs)
  }
  return(Nobs)
}


###################################################
### chunk number 34: 
###################################################
N = linsim(sd_proc=2)
tvec = 1:20
lm1 = lm(N~tvec)


###################################################
### chunk number 35: 
###################################################
plot(tvec,N,type="b")
abline(lm1)
abline(a=2,b=1,lty=2)


###################################################
### chunk number 36: 
###################################################
nsim=100
Nmat = matrix(nrow=20,ncol=100)
for (i in 1:nsim) {
  Nmat[,i] = linsim()
}


###################################################
### chunk number 37: 
###################################################
lower = apply(Nmat,1,quantile,0.025)


###################################################
### chunk number 38:  eval=FALSE
###################################################
## nsim=1000
## Nmat = matrix(nrow=20,ncol=nsim)
## for (i in 1:nsim) {
##   Nmat[,i] = linsim(sd_process=2,sd_obs=2)
## }
## matplot(1:20,Nmat,col="gray",type="l",lty=1)
## lines(1:20,rowMeans(Nmat),lwd=2)
## matlines(1:20,t(apply(Nmat,1,quantile,c(0.025,0.975))),lty=2,col=1)
## matlines(1:20,t(apply(Nmat,1,quantile,c(0.005,0.995))),lty=3,col=1)


###################################################
### chunk number 39: 
###################################################
immigsim = function(nt=20,N0=2,immig,surv) {
  N = numeric(nt)
  N[1] = N0
  for (i in 2:nt) {
    Nsurv = rbinom(1,size=N[i-1],prob=surv)
    N[i] = Nsurv+rpois(1,immig)
  }
  return(N)
}


###################################################
### chunk number 40: 
###################################################
nsim=1000
nt=30
p=0.95; N0=2; immig=10
Nmat = matrix(ncol=nsim,nrow=nt)
for (j in 1:nsim) {
 Nmat[,j] = immigsim(nt=nt,N0=N0,surv=p,immig=immig)
}
tvec=1:nt


###################################################
### chunk number 41: 
###################################################
matplot(tvec,Nmat,type="l",col="gray")
lines(tvec,rowMeans(Nmat),lwd=2)
curve(p^(x-1)*N0+(1-p^(x-1))/(1-p)*immig,add=TRUE)


###################################################
### chunk number 42: 
###################################################
library(odesolve)


###################################################
### chunk number 43: 
###################################################
derivfun = function(t,y,parms) {
 r=parms[1]
 K=parms[2]
 theta=parms[3]
 N=y[1]
 dNdt = r*N*sign(1-N/K)*abs((1-N/K))^theta
 list(dNdt,NULL)
}


###################################################
### chunk number 44: 
###################################################
tvec = seq(0,50,by=0.2)
x1 = lsoda(y=c(N=1),times=tvec,func=derivfun,
 parms=c(r=0.2,K=10,theta=1))


###################################################
### chunk number 45: 
###################################################
head(x1)


###################################################
### chunk number 46: 
###################################################
x2 = lsoda(y=c(N=1),times=tvec,func=derivfun,
 parms=c(r=0.2,K=10,theta=2))
x3 = lsoda(y=c(N=1),times=tvec,func=derivfun,
 parms=c(r=0.2,K=10,theta=0.5))


###################################################
### chunk number 47: 
###################################################
X = cbind(x1,x2[,"N"],x3[,"N"])


###################################################
### chunk number 48: 
###################################################
matplot(X[,"time"],X[,2:4],type="l",col=1,
        xlab="time",ylab="N")
r=0.2; K=10; N0=1
curve(K/((1+(K/N0-1)*exp(-r*x))),type="p",add=TRUE)
legend(30,4,c(expression(theta==1),
              expression(theta==2),
              expression(theta==0.5)),
       lty=1:3)


###################################################
### chunk number 49: 
###################################################
nt=20
sim0 = immigsim(nt=nt,N0=2,surv=0.9,immig=10)
tvec = 1:nt


###################################################
### chunk number 50: 
###################################################
lm1=lm(sim0~tvec)
slope=coef(lm1)["tvec"]
ci.slope=confint(lm1)["tvec",]


###################################################
### chunk number 51: 
###################################################
nvec=c(3,5,7,10,15,20)
nsim=500
powsimresults = matrix(nrow=length(nvec)*nsim,ncol=5)
colnames(powsimresults) = c("n","sim","slope","slope.lo","slope.hi")
ctr=1
for (i in 1:length(nvec)) {
  nt = nvec[i]
  tvec=1:nt
  cat(nt,"\n")
  for (sim in 1:nsim) {
    current.sim = immigsim(nt=nt,N0=N0,surv=p,immig=immig)
    lm1=lm(current.sim~tvec)
    slope=coef(lm1)["tvec"]
    ci.slope=confint(lm1)["tvec",]
    powsimresults[ctr,] = c(nt,sim,slope,ci.slope)
    ctr = ctr+1
  }
}


###################################################
### chunk number 52: 
###################################################
nfac = factor(powsimresults[,"n"])


###################################################
### chunk number 53: 
###################################################
slope.mean = tapply(powsimresults[,"slope"],nfac,mean)


###################################################
### chunk number 54: 
###################################################
slope.sd = tapply(powsimresults[,"slope"],nfac,sd)


###################################################
### chunk number 55: 
###################################################
ci.good = (powsimresults[,"slope.hi"]>immig) &  (powsimresults[,"slope.lo"]<immig)


###################################################
### chunk number 56: 
###################################################
nsim=500
slope.cov = tapply(ci.good,nfac,sum)/nsim


###################################################
### chunk number 57: 
###################################################
null.value=0
reject.null = (powsimresults[,"slope.hi"]<null.value) |  (powsimresults[,"slope.lo"]>null.value)


###################################################
### chunk number 58: 
###################################################
slope.pow = tapply(reject.null,nfac,sum)/nsim


###################################################
### chunk number 59: 
###################################################
lm.q=lm(sim0~tvec+I(tvec^2))


###################################################
### chunk number 60: 
###################################################
rzinbinom = function(n,mu,size,zprob) {
  ifelse(runif(n)<zprob,0,rnbinom(n,mu=mu,size=size))
}
shep = function(x,a=0.696,b=9.79,d=1) { a/(1+(a/b)*x^d) }


###################################################
### chunk number 61: 
###################################################
mu=25.32
k=0.932
zprob=0.123
a = 0.696
b = 9.79
d = 1.1
n=603


###################################################
### chunk number 62: 
###################################################
set.seed(1002)
settlers = rzinbinom(n,mu=mu,size=k,zprob=zprob)
recr = rbinom(n,prob=shep(settlers,a,b,d),size=settlers)


###################################################
### chunk number 63: 
###################################################
bh.fit = nls(recr~a*settlers/(1+(a/b)*settlers),start=c(a=0.696,b=9.79)); bh.fit


###################################################
### chunk number 64: 
###################################################
shep.fit = nls(recr~a*settlers/(1+(a/b)*settlers^d),start=c(coef(bh.fit),d=1)); shep.fit


###################################################
### chunk number 65: 
###################################################
ci = confint(shep.fit); ci


###################################################
### chunk number 66: 
###################################################
ci["d",]


###################################################
### chunk number 67: 
###################################################
getvals = function(n=100,d=1) {
  OK = FALSE
  while (!OK) {
    z = simdata(n,d)
    bh.fit = try(nls(recr~a*settlers/(1+(a/b)*settlers),start=c(a=0.696,b=9.79),data=z))
    shep.fit = try(nls(recr~a*settlers/(1+(a/b)*settlers^d),start=c(coef(bh.fit),d=1),data=z))
    OK = (class(shep.fit)!="try-error" && class(bh.fit)!="try-error")
    if (OK) {
      bh.ci = try(confint(bh.fit))
      shep.ci = try(confint(shep.ti))
      OK = (class(shep.ci)!="try-error" && class(bh.ci)!="try-error")
    }
  }
  res = c(coef(bh.fit),bh.ci,coef(shep.fit),shep.ci)
  names(res) = c("a.bh","b.bh","a.bh.lo","b.bh.lo","a.bh.hi","b.bh.hi",
         "a.shep","b.shep","d.shep","a.shep.lo","b.shep.lo","d.shep.lo",
         "a.shep.hi","b.shep.hi","d.shep.hi")
  res
}


###################################################
### chunk number 68: 
###################################################
load("chap5-batch2.RData")


###################################################
### chunk number 69: 
###################################################
faclist = list(factor(resmat[,"d"]),factor(resmat[,"n"]))


###################################################
### chunk number 70: 
###################################################
d.shep.mean = tapply(resmat[,"d.shep"],faclist,mean)


###################################################
### chunk number 71: 
###################################################
d.shep.sd = tapply(resmat[,"d.shep"],faclist,sd)


###################################################
### chunk number 72: 
###################################################
ci.good = (resmat[,"d.shep.hi"]>resmat[,"d"]) &  (resmat[,"d.shep.lo"]<resmat[,"d"])


###################################################
### chunk number 73: 
###################################################
nsim=400
d.shep.cov = tapply(ci.good,faclist,sum)/nsim


###################################################
### chunk number 74: 
###################################################
null.value=1
reject.null = (resmat[,"d.shep.hi"]<null.value) |  (resmat[,"d.shep.lo"]>null.value)


###################################################
### chunk number 75: 
###################################################
nsim=400
d.shep.pow = tapply(reject.null,faclist,sum)/nsim


###################################################
### chunk number 76: 
###################################################
nrand = 1000


###################################################
### chunk number 77: 
###################################################
pool = c(x0,x1)
n0 = length(x0)
n1 = length(x1)


###################################################
### chunk number 78: 
###################################################
rand.results = numeric(nrand)


###################################################
### chunk number 79: 
###################################################
for (i in 1:nrand) {
  x0.b = sample(pool,size=n0,replace=TRUE)
  x1.b = sample(pool,size=n1,replace=TRUE)
  rand.results[i] = t.stat(x0.b,x1.b)
}


###################################################
### chunk number 80: 
###################################################
s = sort(rand.results)
c(s[25],s[975])


###################################################
### chunk number 81: 
###################################################
clim = quantile(rand.results,c(0.025,0.975)); clim


###################################################
### chunk number 82: 
###################################################
hist(rand.results,breaks=100)
abline(v=clim,lwd=2,col="red")
abline(v=t0,lwd=2,col="blue")


