###################################################
### chunk number 1: 
###################################################
sh = 4
op=par(mfrow=c(1,2),mar=c(5.1,9.1,0,0.5))
curve(dgamma(x,shape=sh),from=0,to=20,ylab="",lwd=2,axes=FALSE)
axis(side=1,labels=FALSE)
axis(side=2,labels=FALSE)
box()
xvec = seq(0,5,length=100)
polygon(c(xvec,rev(xvec)),c(rep(0,length(xvec)),dgamma(rev(xvec),shape=sh)),col="gray",border=NA)
curve(dgamma(x,shape=sh),from=0,to=20,lwd=2,add=TRUE)
abline(v=5,lty=3)
mtext(side=1,line=1.8,at=5,expression(x[0]),cex=2)
abline(h=dgamma(5,shape=sh),lty=2,lwd=2)
mtext(side=2,at=dgamma(5,shape=sh),las=1,expression(ddist(x[0])),
      line=1.8,cex=2)
text(10,0.1,expression(pdist(x[0])),cex=2,col="darkgray")
mtext(side=2,at=0.0,adj=0,"Probability density",cex=2,line=3.5)
set.seed(1001)
points(rgamma(10,shape=sh),rep(0,10),cex=1.5)
text(11.7,0.03,"rdist(10)",adj=0,cex=2)
arrows(10.8,0.023,6.6,0.008,lwd=2)
curve(pgamma(x,shape=sh),from=0,to=20,ylab="",lwd=2,axes=FALSE)
axis(side=1,labels=FALSE)
axis(side=2,labels=FALSE)
box()
abline(v=5,lty=3)
mtext(side=1,line=1.8,at=5,expression(x[0]),cex=2)
abline(h=pgamma(5,shape=sh),lty=2,lwd=2)
mtext(side=2,at=pgamma(5,shape=sh),las=1,expression(pdist(x[0])),
      line=1.8,cex=2)
abline(v=qgamma(0.95,shape=sh),lty=4,lwd=2)
mtext(side=2,at=0.95,las=1,0.95,
      line=par("mgp")[2],cex=2)
segments(par("usr")[1],0.95,qgamma(0.95,shape=sh),0.95,lty=4,lwd=2)
mtext(side=1,at=qgamma(0.95,shape=sh),text="qdist(0.95)",line=1.8,
      cex=2,adj=0.1)
mtext(side=2,at=-0.05,adj=0,"Cumulative distribution",cex=2,line=3.5)


###################################################
### chunk number 2: 
###################################################
rbinom(10,size=8,p=0.5)
rbinom(3,size=8,p=c(0.2,0.4,0.6))


###################################################
### chunk number 3:  eval=FALSE
###################################################
## plot(factor(rbinom(200,size=12,p=0.1)),xlab="# of successes",ylab="# of trials out of 200")


###################################################
### chunk number 4: 
###################################################
plot(factor(rbinom(200,size=12,p=0.1)),xlab="# of successes",ylab="# of trials out of 200")


###################################################
### chunk number 5:  eval=FALSE
###################################################
## set.seed(1001); rbinom(8,size=10,prob=0.2)


###################################################
### chunk number 6: 
###################################################
dbinom(3:5,size=10,prob=0.2)


###################################################
### chunk number 7: 
###################################################
pbinom(4,size=10,prob=0.2,lower.tail=FALSE)


###################################################
### chunk number 8: 
###################################################
set.seed(1001)
N=20; p=0.2
x = rbinom(10000,prob=p,size=N)
c(mean(x),var(x))


###################################################
### chunk number 9: 
###################################################
var_dist = replicate(1000,var(rbinom(10000,prob=p,size=N)))


###################################################
### chunk number 10: 
###################################################
summary(var_dist)
quantile(var_dist,c(0.025,0.975))


###################################################
### chunk number 11:  eval=FALSE
###################################################
## x = rbinom(10000,prob=p,size=N)


###################################################
### chunk number 12:  eval=FALSE
###################################################
## tx = table(factor(x,levels=0:12))/10000


###################################################
### chunk number 13:  eval=FALSE
###################################################
## b1 = barplot(tx,ylim=c(0,0.23),ylab="Probability")


###################################################
### chunk number 14:  eval=FALSE
###################################################
## points(b1,dbinom(0:12,prob=p,size=N),pch=16)


###################################################
### chunk number 15:  eval=FALSE
###################################################
## plot(factor(x)); points(b1,10000*dbinom(0:12,prob=p,size=N))


###################################################
### chunk number 16:  eval=FALSE
###################################################
## plot(table(x)/10000); points(0:12,dbinom(0:12,prob=p,size=N))


###################################################
### chunk number 17:  eval=FALSE
###################################################
## h= hist(x,breaks=seq(-0.5,12.5,by=1),col="gray",
##      prob=TRUE)
## points(0:12,dbinom(0:12,prob=p,size=N))


###################################################
### chunk number 18: 
###################################################
x = rbinom(10000,prob=p,size=N)
tx = table(factor(x,levels=0:12))/10000
b1 = barplot(tx,ylim=c(0,0.23),ylab="Probability")
points(b1,dbinom(0:12,prob=p,size=N),pch=16)


###################################################
### chunk number 19: 
###################################################
dat = c(5,6,5,7,5,8); dat
tabdat=table(dat); tabdat


###################################################
### chunk number 20: 
###################################################
prob=tabdat/length(dat); prob


###################################################
### chunk number 21: 
###################################################
vals = as.numeric(names(prob))
sum(prob*vals)


###################################################
### chunk number 22: 
###################################################
rep(c(1,2,3),c(2,1,5))


###################################################
### chunk number 23: 
###################################################
rep(vals,tabdat)


###################################################
### chunk number 24: 
###################################################
a = 0.696; b=9.79
dS = 0.1
S = seq(0,200,by=dS)
pS =  dexp(S,rate=1/24.5)
fS = a*S/(1+(a/b)*S)
sum(pS*fS*dS)


###################################################
### chunk number 25: 
###################################################
tmpf = function(S) { dexp(S,rate=1/24.5)*a*S/(1+(a/b)*S) }


###################################################
### chunk number 26: 
###################################################
i1 = integrate(tmpf,lower=0,upper=Inf); i1


###################################################
### chunk number 27: 
###################################################
d1 = D(expression(a*x/(1+(a/b)*x)),"x")
d2 = D(d1,"x")


###################################################
### chunk number 28: 
###################################################
Smean = 24.5
Svar = Smean^2
d2_num = eval(d2,list(a=0.696,b=9.79,x=Smean))
mval = a*Smean/(1+(a/b)*Smean)
dapprox = mval + 1/2*Svar*d2_num
merr = (mval-i1$value)/i1$value; merr
err = (dapprox-i1$value)/i1$value; err


###################################################
### chunk number 29: 
###################################################
my_dnbinom = function(x,mean,var,...) {
  mu = mean
  k = mean/(var/mean-1)
  dnbinom(x,mu=mu,size=k,...)
}

my_rnbinom = function(n,mean,var,...) {
  mu = mean
  k = mean/(var/mean-1)
  rnbinom(n,mu=mu,size=k,...)
}


###################################################
### chunk number 30: 
###################################################
x = my_rnbinom(100000,mean=1,var=4)
mean(x)
var(x)


###################################################
### chunk number 31: 
###################################################
tx = table(factor(x,levels=0:max(x)))/100000
b1 = barplot(tx,ylab="Probability")
points(b1,my_dnbinom(0:max(x),mean=1,var=4),pch=16)
abline(v=1)


###################################################
### chunk number 32: 
###################################################
dzinbinom = function(x,mu,size,zprob) {
  ifelse(x==0,
         zprob+(1-zprob)*dnbinom(0,mu=mu,size=size),
         (1-zprob)*dnbinom(x,mu=mu,size=size))
}


###################################################
### chunk number 33: 
###################################################
rzinbinom = function(n,mu,size,zprob) {
  ifelse(runif(n)<zprob,
         0,
         rnbinom(n,mu=mu,size=size))
}


###################################################
### chunk number 34: 
###################################################
rbinom(8,size=10,prob=0.8)


###################################################
### chunk number 35: 
###################################################
clutch_size = c(10,9,9,12,10,10,8,11)
rbinom(8,size=clutch_size,prob=0.8)


###################################################
### chunk number 36: 
###################################################
clutch_size = rpois(8,lambda=10)
rbinom(8,size=clutch_size,prob=0.8)


###################################################
### chunk number 37: 
###################################################
var_vals=1/rgamma(10000,shape=5,scale=1/5)


###################################################
### chunk number 38: 
###################################################
sd_vals = sqrt(var_vals)


###################################################
### chunk number 39: 
###################################################
x = rnorm(10000,mean=0,sd=sd_vals)


###################################################
### chunk number 40:  eval=FALSE
###################################################
## hist(x,prob=TRUE,breaks=100,col="gray")
## curve(dt(x,df=11),add=TRUE,lwd=2)


###################################################
### chunk number 41: 
###################################################
hist(x,prob=TRUE,breaks=100,col="gray")
curve(dt(x,df=11),add=TRUE,lwd=2)


