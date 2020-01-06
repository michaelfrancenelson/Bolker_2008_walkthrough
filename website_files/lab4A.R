###################################################
### chunk number 1: 
###################################################
set.seed(1001)
x = rbinom(n=8,size=10,prob=0.2)
sort(x)


###################################################
### chunk number 2: 
###################################################
dbinom(3:5,size=10,prob=0.2)


###################################################
### chunk number 3: 
###################################################
sum(dbinom(5:10,size=10,prob=0.2))


###################################################
### chunk number 4: 
###################################################
1-pbinom(4,size=10,prob=0.2)


###################################################
### chunk number 5: 
###################################################
pbinom(4,size=10,prob=0.2,lower.tail=FALSE)


###################################################
### chunk number 6: 
###################################################
qbinom(c(0.025,0.975),prob=0.2,size=10)


###################################################
### chunk number 7: 
###################################################
mu=2; k=0.5
x=rnbinom(10000,mu=mu,size=k)
tx = table(factor(x,levels=0:max(x)))/10000
b1 = barplot(tx,ylab="Probability")
points(b1,dnbinom(0:max(x),mu=mu,size=k),pch=1)
mean(x); var(x)
mu; mu*(1+mu/k)
p = 1/(1+mu/k)
n = k
points(b1,dnbinom(0:max(x),prob=p,size=k),pch=2)


###################################################
### chunk number 8: 
###################################################
a = 0.696; b=9.79
d1 = D(expression(a*x/(1+(a/b)*x)),"x")
d2 = D(d1,"x")
Smean = 24.5
d2_num = eval(d2,list(a=0.696,b=9.79,x=Smean))
mval = a*Smean/(1+(a/b)*Smean)


###################################################
### chunk number 9: 
###################################################
tmpf = function(S,mean=Smean,var) { dgamma(S,shape=mean^2/var,scale=var/mean)*a*S/(1+(a/b)*S) }


###################################################
### chunk number 10: 
###################################################
integrate(tmpf,lower=0,upper=Inf,var=Smean^2)


###################################################
### chunk number 11: 
###################################################
Svar_vec = c(Smean^2,100,25,1)
dapprox = mval + 1/2*Svar_vec*d2_num
exact = c(integrate(tmpf,lower=0,upper=Inf,var=Smean^2)$value,
  integrate(tmpf,lower=0,upper=Inf,var=100)$value,
  integrate(tmpf,lower=0,upper=Inf,var=25)$value,
  integrate(tmpf,lower=0,upper=Inf,var=1)$value)
merr = (mval-exact)/exact
err = (dapprox-exact)/exact
data.frame(exact=exact,mval=mval,delta=dapprox,mval.err=merr,delta.err=err)


###################################################
### chunk number 12: 
###################################################
tmpf2 = function(var) {
  integrate(tmpf,lower=0,upper=Inf,var=var)$value
}
sapply(Svar_vec,tmpf2)


###################################################
### chunk number 13: 
###################################################
my_rbeta = function(n,theta,P) {
  rbeta(n,shape1=theta*P,shape2=theta*(1-P))
}
my_dbeta = function(x,theta,P) {
  dbeta(x,shape1=theta*P,shape2=theta*(1-P))
}


###################################################
### chunk number 14: 
###################################################
x = my_rbeta(1000,theta=10,P=0.2)
hist(x,breaks=50,prob=TRUE,col="gray")
curve(my_dbeta(x,theta=10,P=0.2),add=TRUE,lwd=2)
abline(v=0.2,lwd=2,lty=3)
abline(v=mean(x),lty=2)


###################################################
### chunk number 15: 
###################################################
dzinbinom = function(x, mu, size, zprob) {
     ifelse(x == 0, zprob + (1 - zprob) * dnbinom(0, mu=mu, size=size), 
         (1 - zprob) * dnbinom(x, mu=mu, size=size))
}
rzinbinom = function(n, mu, size, zprob) {
    ifelse(runif(n) < zprob, 0, rnbinom(n, mu=mu, size=size))
}


###################################################
### chunk number 16: 
###################################################
mu=4; size=0.5; zprob=0.2
x = rzinbinom(10000,mu=mu,size=size,zprob=zprob)
tx = table(factor(x,levels=0:max(x)))/10000
b1 = barplot(tx,ylab="Probability",ylim=c(0,0.5))
points(b1,dzinbinom(0:max(x),mu=mu,size=size,zprob=zprob),pch=16)
points(b1[1],dnbinom(0,mu=mu,size=size)*(1-zprob),pch=16,col=2)


###################################################
### chunk number 17: 
###################################################
mu*(1-zprob)
mean(x)


###################################################
### chunk number 18: 
###################################################
mu=4; k=0.5
x = rpois(10000,rgamma(10000,shape=k,scale=mu/k))
plot(table(x)/10000)
points(0:max(x),dnbinom(0:max(x),mu=mu,size=k),cex=0.75)


###################################################
### chunk number 19: 
###################################################
mu=2.5; sigmasq=3
m = exp(mu+sigmasq/2)
v = exp(2*mu+sigmasq)*(exp(sigmasq)-1)
s2 = log(v/m^2+1); s2
m2 = log(m)-s2/2; m2


###################################################
### chunk number 20: 
###################################################
nsim=100000
s3 = log(32/4^2+1); s3
m3 = log(4)-s3/2; m3
lnormvals = rlnorm(nsim,meanlog=m3,sdlog=sqrt(s3))
mean(lnormvals); var(lnormvals)
poislnormvals = rpois(nsim,lnormvals)


###################################################
### chunk number 21: 
###################################################
plot(table(factor(poislnormvals,levels=0:max(poislnormvals)))/nsim,xlim=c(0,50))
points(0:50,dnbinom(0:50,mu=mu,size=k),cex=0.75,col=2)


###################################################
### chunk number 22: 
###################################################
x2 = as.numeric(table(factor(poislnormvals,levels=0:max(poislnormvals))))/nsim
plot(0:max(poislnormvals),x2,log="y")
points(0:50,dnbinom(0:50,mu=mu,size=k),cex=0.75,col=2)


###################################################
### chunk number 23: 
###################################################
var(lnormvals)
var(poislnormvals)


