###################################################
### chunk number 1: 
###################################################
set.seed(1001)
mu.true=1
k.true=0.4
x = rnbinom(50,mu=mu.true,size=k.true)
NLLfun1 = function(p,dat=x) {
  mu=p[1]
  k=p[2]
  -sum(dnbinom(x,mu=mu,size=k,log=TRUE))
}


###################################################
### chunk number 2: 
###################################################
m = mean(x)
v = var(x)
mu.mom = m
k.mom = m/(v/m-1)


###################################################
### chunk number 3: 
###################################################
O1 = optim(fn=NLLfun1,par=c(mu=mu.mom,k=k.mom)); O1


###################################################
### chunk number 4: 
###################################################
muvec = seq(0.4,3,by=0.05)
kvec = seq(0.01,0.7,by=0.01)
resmat = matrix(nrow=length(muvec),ncol=length(kvec))
for (i in 1:length(muvec)) {
  for (j in 1:length(kvec)) {
    resmat[i,j] = NLLfun1(c(muvec[i],kvec[j]))
  }
}


###################################################
### chunk number 5: 
###################################################
alevels=c(0.5,0.9,0.95,0.99,0.999)
minval = O1$value
nll.levels = qchisq(alevels,df=2)/2+minval


###################################################
### chunk number 6: 
###################################################
contour(muvec,kvec,resmat,levels=nll.levels,labels=alevels)
points(O1$par["mu"],O1$par["k"],pch=16)
points(mu.true,k.true,pch=1)
points(mu.mom,k.mom,pch=2)


###################################################
### chunk number 7: 
###################################################
a = 0.696
b = 9.79
recrprob = function(x,a=0.696,b=9.79) a/(1+(a/b)*x)
scoefs = c(mu=25.32,k=0.932,zprob=0.123)
rzinbinom = function(n,mu,size,zprob) {
  ifelse(runif(n)<zprob,
         0,
         rnbinom(n,mu=mu,size=size))
}

settlers = rzinbinom(603,mu=scoefs["mu"],size=scoefs["k"],zprob=scoefs["zprob"])
recr = rbinom(603,prob=recrprob(settlers),size=settlers)


###################################################
### chunk number 8: 
###################################################
NLLfun3L = function(loga,logb,logd) {
  a = exp(loga)
  b = exp(logb)
  d = exp(logd)
  recrprob = a/(1+(a/b)*settlers^d)
  -sum(dbinom(recr,prob=recrprob,size=settlers,log=TRUE),na.rm=TRUE)
}

NLLfun4L = function(loga,logb) {
  a = exp(loga)
  b = exp(logb)
  recrprob = a/(1+(a/b)*settlers)
  -sum(dbinom(recr,prob=recrprob,size=settlers,log=TRUE),na.rm=TRUE)
}

NLLfun5L = function(loga) {
  a = exp(loga)
  recrprob = a
  r = -sum(dbinom(recr,prob=recrprob,size=settlers,log=TRUE),na.rm=TRUE)
  ##  cat(loga,a,r,"\n")
  return(r)
}


###################################################
### chunk number 9: 
###################################################
lm1 = lm(recr~settlers-1)
lm.loga = list(log(coef(lm1)))
rename=function(L,names) {
  for (i in seq(along=L)) {
    names(L[[i]]) = NULL
  }
  names(L)=names
  L
}
lm.loga=rename(lm.loga,"loga")


###################################################
### chunk number 10: 
###################################################
library(stats4)
m5L = mle(minuslogl=NLLfun5L,start=list(loga=-1),
  method="BFGS",control=list(ndeps=0.01))                              


###################################################
### chunk number 11: 
###################################################
library(stats4)
m4L = mle(minuslogl=NLLfun4L,start=list(loga=log(0.5),logb=log(10)),
  method="Nelder-Mead")


###################################################
### chunk number 12: 
###################################################
s3 = c(coef(m4L),list(logd=0))


###################################################
### chunk number 13: 
###################################################
m3L = mle(minuslogl=NLLfun3L,start=s3,
  method="Nelder-Mead")


###################################################
### chunk number 14: 
###################################################
lin = c(exp(coef(m5L)),NA,NA,-logLik(m5L))
BH =  c(exp(coef(m4L)),NA,-logLik(m4L))
shep = c(exp(coef(m3L)),-logLik(m3L))
ptab = rbind(lin,BH,shep)
colnames(ptab)=c("a","b","d","NLL")
ptab


###################################################
### chunk number 15: 
###################################################
exp(confint(m3L))
exp(confint(m4L))
exp(confint(m5L))


###################################################
### chunk number 16: 
###################################################
NLLpois = function(lambda) {
  -sum(dpois(settlers,lambda=lambda,log=TRUE))
}
NLLnb = function(mu,k) {
  -sum(dnbinom(settlers,mu=mu,size=k,log=TRUE))
}


###################################################
### chunk number 17: 
###################################################
dzinbinom = function(x,mu,size,zprob,log=FALSE) {
  v = ifelse(x==0,
         zprob+(1-zprob)*dnbinom(0,mu=mu,size=size),
         (1-zprob)*dnbinom(x,mu=mu,size=size))
  if (log) return(log(v)) else v
}


###################################################
### chunk number 18: 
###################################################
NLLzinb = function(mu,k,zprob) {
  -sum(dzinbinom(settlers,mu=mu,size=k,zprob=zprob,log=TRUE))
}


###################################################
### chunk number 19: 
###################################################
m = mean(settlers)
mpois = mle(minuslogl=NLLpois,start=list(lambda=m))
mnbinom = mle(minuslogl=NLLnb,start=list(mu=m,k=0.5))
mzinbinom = mle(minuslogl=NLLzinb,start=list(mu=m,k=0.5,zprob=0.5))


###################################################
### chunk number 20: 
###################################################
pois = c(coef(mpois),NA,NA,-logLik(mpois))
nbin =  c(coef(mnbinom),NA,-logLik(mnbinom))
zinb = c(coef(mzinbinom),-logLik(mzinbinom))
ptab = rbind(pois,nbin,zinb)
colnames(ptab)=c("lambda/mu","k","zprob","NLL")
ptab


###################################################
### chunk number 21: 
###################################################
confint(mpois)
confint(mnbinom)
confint(mzinbinom)


