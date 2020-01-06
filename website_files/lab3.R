###################################################
### chunk number 1:  eval=FALSE
###################################################
## curve(2*exp(-x/2),from=0,to=7,ylim=c(0,2),ylab="")
## curve(2*exp(-x),add=TRUE,lty=4)
## curve(x*exp(-x/2),add=TRUE,lty=2)
## curve(2*x*exp(-x/2),add=TRUE,lty=3)
## text(0.4,1.9,expression(paste("exponential: ",2*e^(-x/2))),adj=0)
## text(4,0.7,expression(paste("Ricker: ",x*e^(-x/2))))
## text(4,.25,expression(paste("Ricker: ",2*x*e^(-x/2))),adj=0)
## text(2.8,0,expression(paste("exponential: ",2*e^(-x))))


###################################################
### chunk number 2: 
###################################################
xvec = seq(0,7,length=100)
exp1_vec = 2*exp(-xvec/2)
exp2_vec = 2*exp(-xvec)
plot(xvec,exp1_vec,type="l",ylim=c(0,2),ylab="")
lines(xvec,exp2_vec,lty=4)


###################################################
### chunk number 3:  eval=FALSE
###################################################
## matplot(xvec,cbind(exp1_vec,exp2_vec),type="l",
##         ylab="")


###################################################
### chunk number 4: 
###################################################
expfun = function(x,a=1,b=1) {
   a*exp(-b*x)
 }
exp1_vec = sapply(xvec,expfun,a=2,b=1/2)
exp2_vec = sapply(xvec,expfun,a=2,b=1)


###################################################
### chunk number 5: 
###################################################
x=c(-25,-16,-9,-4,-1,0,1,4,9,16,25)
ifelse(x<0,0,sqrt(x))


###################################################
### chunk number 6: 
###################################################
op=par(mfrow=c(2,2),mgp=c(2,1,0),mar=c(4.2,3,1,1))
curve(ifelse(x<2,1,3),from=0,to=5)
curve(ifelse(x<2,2*x,4),from=0,to=5)
curve(ifelse(x<2,exp(x),exp(2)-3*(x-2)),from=0,to=5)
curve(ifelse(x<2,1,ifelse(x<4,3,5)),from=0,to=5)


###################################################
### chunk number 7: 
###################################################
d1 = D(expression(x^2),"x"); d1


###################################################
### chunk number 8: 
###################################################
eval(d1,list(x=2))


###################################################
### chunk number 9: 
###################################################
D(d1,"x")


###################################################
### chunk number 10: 
###################################################
r = 1
K = 1
n0 = 0.1  # what happens if you set n0 to 0???
curve(K/(1+(K/n0-1)*exp(-r*x)),from=0,to=10)


###################################################
### chunk number 11: 
###################################################
t_vec = seq(0,10,length=100)
logist_vec = K/(1+(K/n0-1)*exp(-r*t_vec))    
plot(t_vec,logist_vec,type="l")


###################################################
### chunk number 12: 
###################################################
logistfun = function(t,r=1,n0=0.1,K=1) {
  K/(1+(K/n0-1)*exp(-r*t))
}
logist_vec = sapply(t_vec,logistfun)


###################################################
### chunk number 13: 
###################################################
r=17
logistfun(1,r=2)
r=0
logistfun(1,r=2)


###################################################
### chunk number 14: 
###################################################
curve(logistfun(x),from=0,to=10,lwd=2)
abline(h=n0,lty=2)
abline(h=K,lty=2)
abline(h=K/2,lty=3)
abline(v=-log(n0/(K-n0))/r,lty=4)
r=1
abline(a=n0,b=r*n0*(1-n0/K),lty=5)
curve(n0*exp(r*x),from=0,lty=6,add=TRUE)


