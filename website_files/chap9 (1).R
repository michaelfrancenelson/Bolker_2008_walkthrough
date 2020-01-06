###################################################
### chunk number 1: 
###################################################
library(emdbook)
library(bbmle)
library(MASS)
library(nlme)
source("chapskel.R")


###################################################
### chunk number 2: 
###################################################
data(FirDBHFec)
firdata = na.omit(FirDBHFec[,c("TOTCONES","DBH","WAVE_NON")])
firdata$TOTCONES = round(firdata$TOTCONES)
cones = 1+firdata$TOTCONES
lm1 = lm(log(cones)~log(DBH),data=firdata)
lm2 = lm(cones~DBH,data=firdata)
nbNLL.0 = function(a,b,k) {
  predcones = a*firdata$DBH^b
  -sum(dnbinom(firdata$TOTCONES,mu=predcones,size=k,log=TRUE))
}


###################################################
### chunk number 3: 
###################################################
nbfit.0 = mle2(nbNLL.0,start=list(a=1,b=1,k=1))
cones2 = cones[firdata$DBH>6 & firdata$DBH<8]
ndb = fitdistr(cones2,"negative binomial")
nln = fitdistr(cones2,"lognormal")
ngm = fitdistr(cones2,"gamma")


###################################################
### chunk number 4: 
###################################################
op=par(mfrow=c(1,2),lwd=2,bty="l",las=1,cex=1.5)
plot(firdata$DBH,cones,log="xy",
     xlab="DBH",
     ylab="cones+1",col="gray",cex=0.7,axes=FALSE)
axis(side=2)
axis(side=1,at=seq(5,20,by=5))
box()
abline(v=c(6,8),col="darkgray")
curve(exp(coef(lm1)[1])*x^coef(lm1)[2],add=TRUE)
with(as.list(coef(nbfit.0)),curve(a*(x-1)^b,add=TRUE,lty=2))
curve(coef(lm2)[1]+coef(lm2)[2]*x,add=TRUE,lty=3)
par(xpd=NA)
legend(10.5,6,lty=1:3,
       c("power/LN",
         "power/NB",
         "lin/normal"),
       cex=0.9,bty="n")
par(xpd=FALSE)
plot(density(cones2,adjust=0.8,from=0),ylim=c(0,0.04),
     main="",xlab="cones+1",lty=2,xlim=c(-20,100),
     axes=FALSE)
axis(side=2)
axis(side=1,at=seq(0,100,by=25))
box()
text(40,0.04,"6 < DBH < 8")
points(cones2,runif(length(cones2),max=0.001),col="gray")
with(as.list(ndb$estimate),lines(0:100,dnbinom(0:100,size=size,mu=mu),type="s"),
     lty=1)
with(as.list(nln$estimate),curve(dlnorm(x,meanlog,sdlog),add=TRUE,lty=3))
curve(dnorm(x,mean(cones2),sd(cones2)),add=TRUE,lty=4)
abline(v=0,col="darkgray")
legend(25,0.038,
       lty=c(2,1,3,4),
       c("density","NB","LN","normal"),bty="n")

##with(as.list(ngm$estimate),curve(dgamma(x,shape=shape,rate=rate),add=TRUE,lty=4))


###################################################
### chunk number 5:  eval=FALSE
###################################################
## lm.reg = lm(y~x)


###################################################
### chunk number 6: 
###################################################
linregfun = function(a,b,sigma) {
  Y.pred = a+b*x
  -sum(dnorm(Y,mean=Y.pred,sd=sigma,log=TRUE))
}


###################################################
### chunk number 7:  eval=FALSE
###################################################
## mle2(Y~dnorm(mean=a+b*x,sd=sigma),start=...)


###################################################
### chunk number 8:  eval=FALSE
###################################################
## lm.poly = lm(y~x+I(x^2))


###################################################
### chunk number 9:  eval=FALSE
###################################################
## lm.mreg = lm(y~x1+x2+x3)


###################################################
### chunk number 10:  eval=FALSE
###################################################
## lm.1way = lm(y~f)


###################################################
### chunk number 11:  eval=FALSE
###################################################
## aov1fun = function(m1,m2,m3,sigma) {
##   Y.pred = c(m1,m2,m3)[f]
##   -sum(dnorm(Y,mean=Y.pred,sd=sigma,log=TRUE))
## }


###################################################
### chunk number 12:  eval=FALSE
###################################################
## lm.2way = lm(Y~f1*f2)


###################################################
### chunk number 13: 
###################################################
aov2fun = function(m11,m12,m21,m22,sigma) {
  intval = interaction(f1,f2)
  Y.pred = c(m11,m12,m21,m22)[intval]
  -sum(dnorm(Y,mean=Y.pred,sd=sigma,log=TRUE))
}


###################################################
### chunk number 14:  eval=FALSE
###################################################
## mle2(Y~dnorm(mean=m,sd=sigma),parameters=list(m~f1*f2))


###################################################
### chunk number 15:  eval=FALSE
###################################################
## lm(Y~f*x)


###################################################
### chunk number 16: 
###################################################
ancovafun = function(i1,i2,slope1,slope2,sigma) {
  int = c(i1,i2)[f]
  slope = c(slope1,slope2)[f]
  Y.pred =  int+ slope*x
  -sum(dnorm(Y,mean=Y.pred,sd=sigma,log=TRUE))
}


###################################################
### chunk number 17: 
###################################################
data(FirDBHFec)
firdata = na.omit(FirDBHFec[,c("TOTCONES","DBH","WAVE_NON")])
firdata$TOTCONES = round(firdata$TOTCONES)
cones = 1+firdata$TOTCONES
logDBH=log(firdata$DBH)
fir.dwi = lm(log(cones)~logDBH*WAVE_NON,data=firdata)
fir.dw = lm(log(cones)~logDBH+WAVE_NON,data=firdata)
fir.d = lm(log(cones)~logDBH,data=firdata)
fir.w = lm(log(cones)~WAVE_NON,data=firdata)
fir.0 = lm(log(cones)~1,data=firdata)


###################################################
### chunk number 18: 
###################################################
op=par(lwd=2,bty="l",las=1,cex=1.5,mar=c(5,4,2,4)+0.1,mgp=c(2.5,1,0))
plot(log(cones)~logDBH,pch=as.numeric(firdata$WAVE_NON),
     xlab="log(DBH)",ylab="log(cones+1)")
abline(a=coef(fir.dw)[1],b=coef(fir.dw)[2])
abline(a=coef(fir.dw)[1]+coef(fir.dw)[3],b=coef(fir.dw)[2],lty=2)
par(xpd=NA)
legend(2.65,1.5,lty=1:2,
       pch=1:2,c("nonwave","wave"),bty="n")
par(xpd=FALSE)
par(op)


###################################################
### chunk number 19:  eval=FALSE
###################################################
## n1 = nls(y~a*x^b,start=list(a=1,b=1))


###################################################
### chunk number 20:  eval=FALSE
###################################################
## nlsList(TOTCONES~a*DBH^b|WAVE_NON,
## data=firdata,start=list(a=0.1,b=2.7))


###################################################
### chunk number 21: 
###################################################
nlsList(TOTCONES~a*DBH^b|WAVE_NON,data=firdata,start=list(a=0.1,b=2.7))
gnls(TOTCONES~a*DBH^b,data=firdata,start=c(0.1,2.7,2.7),
     params=list(a~1,b~WAVE_NON))


###################################################
### chunk number 22:  eval=FALSE
###################################################
## gnls(TOTCONES~a*DBH^b,data=firdata,start=c(0.1,2.7,2.7),
## params=list(a~1,b~WAVE_NON))


###################################################
### chunk number 23: 
###################################################
data(ReedfrogFuncresp)
nls(Killed ~ SSmicmen(Initial, a, b), data = ReedfrogFuncresp)


###################################################
### chunk number 24: 
###################################################
op=par(lwd=2,bty="l",las=1,cex=1.5,mgp=c(2.5,1,0))
firnls = nls(TOTCONES~a*DBH^b,start=list(a=0.1,b=2.7),
              data=firdata)
plot(cones~DBH,data=firdata,ylab="Cones")
with(as.list(coef(firnls)),curve(a*x^b,add=TRUE))
curve(exp(coef(fir.d)[1])*x^(coef(fir.d)[2]),add=TRUE,lty=2)
legend("topleft",c("power/normal","power/LN"),
       lty=1:2,bty="n")
par(op)


###################################################
### chunk number 25: 
###################################################
op=par(lwd=2,bty="l",las=1,cex=1.5,mgp=c(2.5,1,0))
data(ReedfrogFuncresp)
glm1 = glm(cbind(Killed,Initial-Killed) ~ Initial, 
          family=binomial,data=ReedfrogFuncresp)
glm2 = glm(cbind(Killed,Initial-Killed) ~ Initial, 
          family=binomial(link="log"),data=ReedfrogFuncresp)
with(ReedfrogFuncresp,plot(Killed/Initial~Initial,xlab="Initial density",
                   ylab="Fraction killed"))
ivec = 1:100
lines(ivec,predict(glm1,newdata=list(Initial=ivec),type="response"))
lines(ivec,predict(glm2,newdata=list(Initial=ivec),type="response"),lty=2)
legend("topright",c("logistic regression","log-binomial model"),
       lty=1:2,bty="n")
par(op)


###################################################
### chunk number 26: 
###################################################
glm(cbind(Killed,Initial-Killed)~Initial, 
          data=ReedfrogFuncresp,family="binomial")
glm(cbind(Killed,Initial-Killed)~Initial, 
          data=ReedfrogFuncresp,family="binomial")


###################################################
### chunk number 27:  eval=FALSE
###################################################
## glm1 = glm(y~x,family="poisson")


###################################################
### chunk number 28: 
###################################################
poisregfun = function(a,b) {
  Y.pred = exp(a+b*x)
  -sum(dpois(y,lambda=Y.pred,log=TRUE))
}


###################################################
### chunk number 29:  eval=FALSE
###################################################
## glm2 = glm(cbind(y,N-y)~x,family="binomial")


###################################################
### chunk number 30: 
###################################################
logistregfun = function(a,b) {
  p.pred = exp(a+b*x)/(1+exp(a+b*x))
  -sum(dbinom(y,size=N,prob=p.pred,log=TRUE))
}


###################################################
### chunk number 31:  eval=FALSE
###################################################
## glm3 = glm(cbind(y,N-y)~x,family=binomial(link="log"))


###################################################
### chunk number 32: 
###################################################
logregfun = function(a,b) {
  p.pred = exp(a+b*x)
  -sum(dbinom(y,size=N,prob=p.pred,log=TRUE))
}


###################################################
### chunk number 33:  eval=FALSE
###################################################
## glm4 = glm(cbind(Killed,Initial-Killed) ~ Initial, 
##             family=quasi(link="inverse",variance="mu(1-mu)"),
##             data=ReedfrogFuncresp)
## x = 10:20
## p = 1/(1.2+5.5*x)
## y = rbinom(length(x),size=x,prob=p)
## testglm = glm(cbind(y,x-y)~x,
##                family=quasi(link="inverse",variance="mu(1-mu)"),
##                start=c(1.9,-0.4))
## testglm = glm(cbind(y,x-y)~x,
##                family=quasi(link="inverse",variance="constant"),
##                start=c(1.2,5.5))


###################################################
### chunk number 34: 
###################################################
a = round(coef(glm1)[1],3)
b = round(coef(glm1)[2],4)


###################################################
### chunk number 35: 
###################################################
glm(TOTCONES~log(DBH),data=firdata,family="quasipoisson")


###################################################
### chunk number 36: 
###################################################
glm.nb(TOTCONES~log(DBH),data=firdata)


###################################################
### chunk number 37: 
###################################################
logcones = log(firdata$TOTCONES+1)
lm.0 = lm(logcones~1,data=firdata)
lm.d = lm(logcones~log(DBH),data=firdata)
lm.w = lm(logcones~WAVE_NON,data=firdata)
lm.dw = lm(logcones~log(DBH)+WAVE_NON,data=firdata)
lm.dwi = lm(logcones~log(DBH)*WAVE_NON,data=firdata)


###################################################
### chunk number 38: 
###################################################
anova(lm.0,lm.d,lm.dw,lm.dwi)
AIC(lm.0,lm.d,lm.w,lm.dw,lm.dwi)


###################################################
### chunk number 39: 
###################################################
ancovafun = function(i1,i2,slope1,slope2,sigma) {
  int = c(i1,i2)[WAVE_NON]
  slope = c(slope1,slope2)[WAVE_NON]
  Y.pred =  int+ slope*log(DBH)
  -sum(dnorm(logcones,mean=Y.pred,sd=sigma,log=TRUE))
}
m1 = mle2(ancovafun,start=list(i1=-2,i2=-2,slope1=2.5,slope2=2.5,sigma=1),
  data=firdata)
AIC(m1)


###################################################
### chunk number 40: 
###################################################
coef(lm.dwi)
coef(m1)


###################################################
### chunk number 41:  eval=FALSE
###################################################
## mle2(log(TOTCONES+1)~dnorm(logDBH*WAVE_NON),data=firdata,start=...)


