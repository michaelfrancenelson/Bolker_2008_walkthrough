###################################################
### chunk number 1: 
###################################################
set.seed(1001)


###################################################
### chunk number 2: 
###################################################
mu.true=1
k.true=0.4
x = rnbinom(50,mu=mu.true,size=k.true)


###################################################
### chunk number 3: 
###################################################
plot(table(factor(x,levels=0:max(x))),
     ylab="Frequency",xlab="x")


###################################################
### chunk number 4: 
###################################################
NLLfun1 = function(p,dat=x) {
  mu=p[1]
  k=p[2]
  -sum(dnbinom(x,mu=mu,size=k,log=TRUE))
}


###################################################
### chunk number 5: 
###################################################
nll.true=NLLfun1(c(mu=mu.true,k=k.true)); nll.true


###################################################
### chunk number 6: 
###################################################
NLLfun1(c(mu=10,k=10))


###################################################
### chunk number 7: 
###################################################
m = mean(x)
v = var(x)
mu.mom = m
k.mom = m/(v/m-1)


###################################################
### chunk number 8: 
###################################################
nll.mom=NLLfun1(c(mu=mu.mom,k=k.mom)); nll.mom


###################################################
### chunk number 9: 
###################################################
ldiff=nll.true-nll.mom;ldiff
qchisq(0.95,df=2)/2


###################################################
### chunk number 10: 
###################################################
O1 = optim(fn=NLLfun1,par=c(mu=mu.mom,k=k.mom),hessian=TRUE); O1


###################################################
### chunk number 11: 
###################################################
muvec = seq(0.4,3,by=0.05)
kvec = seq(0.01,0.7,by=0.01)


###################################################
### chunk number 12: 
###################################################
resmat = matrix(nrow=length(muvec),ncol=length(kvec))


###################################################
### chunk number 13: 
###################################################
for (i in 1:length(muvec)) {
  for (j in 1:length(kvec)) {
    resmat[i,j] = NLLfun1(c(muvec[i],kvec[j]))
  }
}


###################################################
### chunk number 14: 
###################################################
contour(muvec,kvec,resmat,xlab=expression(mu),ylab="k")
contour(muvec,kvec,resmat,levels=70:80,lty=2,add=TRUE)


###################################################
### chunk number 15: 
###################################################
alevels=c(0.5,0.9,0.95,0.99,0.999)
minval = O1$value
nll.levels = qchisq(alevels,df=2)/2+minval
contour(muvec,kvec,resmat,levels=nll.levels,labels=alevels,
        xlab=expression(mu),ylab="k")


###################################################
### chunk number 16: 
###################################################
NLLfun.mu = function(p,mu) {
  k = p[1]
  -sum(dnbinom(x,mu=mu,size=k,log=TRUE))
}


###################################################
### chunk number 17: 
###################################################
mu.profile = matrix(ncol=2,nrow=length(muvec))


###################################################
### chunk number 18: 
###################################################
NLLfun.mu2 = function(p,mu) {
  logk = p[1]
  k = exp(logk)
  -sum(dnbinom(x,mu=mu,size=k,log=TRUE))
}


###################################################
### chunk number 19: 
###################################################
for (i in 1:length(muvec)) {
  Oval = optim(fn=NLLfun.mu,par=O1$par["k"],method="L-BFGS-B",
    lower=0.002,mu=muvec[i])
  mu.profile[i,] = c(Oval$par,Oval$value)
}
colnames(mu.profile) = c("k","NLL")


###################################################
### chunk number 20: 
###################################################
NLLfun.k = function(p,k) {
  mu = p[1]
  -sum(dnbinom(x,mu=mu,size=k,log=TRUE))
}
k.profile = matrix(ncol=2,nrow=length(kvec))
for (i in 1:length(kvec)) {
  Oval = optim(fn=NLLfun.k,par=O1$par["mu"],method="L-BFGS-B",
    lower=0.002,k=kvec[i])
  k.profile[i,] = c(Oval$par,Oval$value)
}
colnames(k.profile) = c("mu","NLL")


###################################################
### chunk number 21: 
###################################################
contour(muvec,kvec,resmat,xlab=expression(mu),ylab="k")
contour(muvec,kvec,resmat,levels=70:80,lty=2,add=TRUE)
lines(muvec,mu.profile[,"k"],lwd=2)
lines(k.profile[,"mu"],kvec,lwd=2,lty=2)


###################################################
### chunk number 22: 
###################################################
plot(muvec,mu.profile[,"NLL"],type="l",
     xlab=expression(mu),ylab="Negative log-likelihood")
cutoffs = c(0,qchisq(c(0.95,0.99),1)/2)
nll.levels = O1$value+cutoffs
abline(h=nll.levels,lty=1:3)
text(rep(0.5,3),nll.levels+0.2,c("min","95%","99%"))


###################################################
### chunk number 23: 
###################################################
cutoff = O1$value+qchisq(0.95,1)/2


###################################################
### chunk number 24: 
###################################################
lowerhalf = mu.profile[muvec<1.2,"NLL"]
lowerhalf.mu = muvec[muvec<1.2]
w.lower = which.min(abs(lowerhalf-cutoff))


###################################################
### chunk number 25: 
###################################################
upperhalf = mu.profile[muvec>1.2,"NLL"]
upperhalf.mu = muvec[muvec>1.2]
w.upper = which.min(abs(upperhalf-cutoff))
ci.crude = c(lowerhalf.mu[w.lower],upperhalf.mu[w.upper])


###################################################
### chunk number 26: 
###################################################
plot(muvec,mu.profile[,"NLL"],type="l",
     xlab=expression(mu),ylab="Negative log-likelihood")
cutoffs = c(0,qchisq(c(0.95),1)/2)
nll.levels = O1$value+cutoffs
abline(h=nll.levels,lty=1:2)
abline(v=ci.crude,lty=3)


###################################################
### chunk number 27: 
###################################################
cutoff = O1$value+qchisq(c(0.95),1)/2
relheight = function(mu) {
  O2 = optim(fn=NLLfun.mu,par=O1$par["k"],method="L-BFGS-B",
        lower=0.002,mu=mu)
  O2$value-cutoff
}


###################################################
### chunk number 28: 
###################################################
relheight(mu=0.6)
relheight(mu=0.8)


###################################################
### chunk number 29: 
###################################################
lower = uniroot(relheight,interval=c(0.5,1.0))$root
upper = uniroot(relheight,interval=c(1.2,5))$root
ci.uniroot= c(lower,upper)


###################################################
### chunk number 30: 
###################################################
plot(muvec,mu.profile[,"NLL"],type="l",
     xlab=expression(mu),ylab="Negative log-likelihood")
cutoffs = c(0,qchisq(c(0.95),1)/2)
nll.levels = O1$value+cutoffs
abline(h=nll.levels,lty=1:2)
abline(v=ci.crude,lty=3)
abline(v=ci.uniroot,lty=4)


###################################################
### chunk number 31: 
###################################################
O1$hessian


###################################################
### chunk number 32: 
###################################################
s1 = solve(O1$hessian); s1


###################################################
### chunk number 33: 
###################################################
a=O1$value
b=O1$par["mu"]
c=O1$hessian["mu","mu"]/2


###################################################
### chunk number 34: 
###################################################
se.mu = sqrt(s1["mu","mu"])
ci.info = O1$par["mu"]+c(-1,1)*qnorm(0.975)*se.mu


###################################################
### chunk number 35: 
###################################################
op=par(mfrow=c(1,2))
plot(muvec,mu.profile[,"NLL"],type="l",
     xlab=expression(mu),ylab="Negative log-likelihood",ylim=c(71.8,72.2),
     xlim=c(0.7,1.7))
curve(a+c*(x-b)^2,add=TRUE,lty=2)
##
plot(muvec,mu.profile[,"NLL"],type="l",
     xlab=expression(mu),ylab="Negative log-likelihood")
cutoffs = c(0,qchisq(c(0.95),1)/2)+O1$value
curve(a+c*(x-b)^2,add=TRUE,lty=2)
abline(h=cutoffs)
abline(v=ci.info,lty=3)
abline(v=ci.uniroot,lty=1)
par(op)


###################################################
### chunk number 36: 
###################################################
library(MASS)
f=fitdistr(x,"negative binomial"); f


###################################################
### chunk number 37: 
###################################################
ci.fitdistr = f$estimate["mu"]+c(-1,1)*f$sd["mu"]*qnorm(0.975)


###################################################
### chunk number 38: 
###################################################
NLLfun2 = function(mu,k) {
  -sum(dnbinom(x,mu=mu,size=k,log=TRUE))
}


###################################################
### chunk number 39: 
###################################################
library(stats4)
m1 = mle(minuslogl=NLLfun2,start=list(mu=mu.mom,k=k.mom),
  method="L-BFGS-B",lower=0.002); m1


###################################################
### chunk number 40: 
###################################################
summary(m1)


###################################################
### chunk number 41: 
###################################################
ci.mle.all = confint(m1); ci.mle.all
ci.mle = ci.mle.all["mu",]


###################################################
### chunk number 42: 
###################################################
citab=rbind(ci.crude,ci.uniroot,ci.mle,ci.info,ci.fitdistr); citab


###################################################
### chunk number 43: 
###################################################
a = 0.696
b = 9.79
recrprob = function(x,a=0.696,b=9.79) a/(1+(a/b)*x)
scoefs = c(mu=25.32,k=0.932,zprob=0.123)
settlers = rzinbinom(603,mu=scoefs["mu"],size=scoefs["k"],zprob=scoefs["zprob"])
recr = rbinom(603,prob=recrprob(settlers),size=settlers)


###################################################
### chunk number 44: 
###################################################
NLLfun3 = function(a,b,d) {
  recrprob = a/(1+(a/b)*settlers^d)
  -sum(dbinom(recr,prob=recrprob,size=settlers,log=TRUE),na.rm=TRUE)
}


###################################################
### chunk number 45: 
###################################################
NLLfun4 = function(a,b) {
  recrprob = a/(1+(a/b)*settlers)
  -sum(dbinom(recr,prob=recrprob,size=settlers,log=TRUE),na.rm=TRUE)
}


###################################################
### chunk number 46: 
###################################################
NLLfun5 = function(a) {
  recrprob = a
  -sum(dbinom(recr,prob=recrprob,size=settlers,log=TRUE),na.rm=TRUE)
}


###################################################
### chunk number 47: 
###################################################
recr = recr[settlers>0]
settlers = settlers[settlers>0]


###################################################
### chunk number 48: 
###################################################
plot(settlers,recr)
abline(h=10)
abline(a=0,b=0.5)


###################################################
### chunk number 49: 
###################################################
m4 = mle(minuslogl=NLLfun4,start=list(a=0.5,b=10),
  method="L-BFGS-B",lower=0.003)
s3 = list(a=0.684,b=10.161,d=1)
m3 = mle(minuslogl=NLLfun3,start=s3,
  method="L-BFGS-B",lower=0.003)
m5 = mle(minuslogl=NLLfun5,start=list(a=0.5),
  method="L-BFGS-B",lower=0.003)


###################################################
### chunk number 50: 
###################################################
plot(settlers,recr)
a=coef(m5)["a"]
curve(a*x,add=TRUE,lty=3)
a=coef(m4)["a"]; b=coef(m4)["b"]
curve(a*x/(1+(a/b)*x),add=TRUE,lty=2)
a=coef(m5)["a"]; b=coef(m5)["b"]; d=coef(m5)["d"]
curve(a*x/(1+(a/b)*x^d),add=TRUE,lty=3)


###################################################
### chunk number 51: 
###################################################
nll = c(shep=-logLik(m3),BH=-logLik(m4),densind=-logLik(m5)); nll


###################################################
### chunk number 52: 
###################################################
logp=pchisq(2*nll[3]-nll[2],1,lower.tail=FALSE,log.p=TRUE)
logp


###################################################
### chunk number 53: 
###################################################
npar = c(5,4,3)
aic = nll+2*npar; aic


###################################################
### chunk number 54: 
###################################################
ndata = length(recr)
bic = nll+log(ndata)*npar;bic


###################################################
### chunk number 55: 
###################################################
confint(m3)


###################################################
### chunk number 56: 
###################################################
confint(m4)


