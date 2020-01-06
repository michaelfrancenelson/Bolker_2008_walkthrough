###################################################
### chunk number 1: 
###################################################
source("chapskel.R")
library(emdbook)
library(plotrix)
library(bbmle)


###################################################
### chunk number 2: 
###################################################
data(ReedfrogSizepred)
attach(ReedfrogSizepred,warn=FALSE)
logist2 <- function(x,sizep1,sizep2,sizep3) {
 exp(sizep1*(sizep3-x))/(1+exp(sizep2*sizep1*(sizep3-x)))
}
ricker <- function(x,a,b) {
  a*x*exp(-b*x)
}
powricker <- function(x,a,b,g) {
  b*(x/a*exp(1-x/a))^g
}
tricker <- function(x,a,b,t,min=0.0001) {
  ifelse(x<t,min,b*((x-t)/a*exp(1-(x-t)/a)))
}
lik2 <- function(sizep1,sizep2,sizep3) {
  p.pred <- logist2(TBL,sizep1,sizep2,sizep3)
  -sum(dbinom(Kill,size=10,prob=p.pred,log=TRUE))
}
lik3 <- function(a,b) {
  p.pred <- ricker(TBL,a,b)
  -sum(dbinom(Kill,size=10,prob=p.pred,log=TRUE))
}
lik4 <- function(a,b,g) {
  p.pred <- powricker(TBL,a,b,g)
  -sum(dbinom(Kill,size=10,prob=p.pred,log=TRUE))
}
lik5 <- function(a,b,t) {
  p.pred <- tricker(TBL,a,b,t)
  -sum(dbinom(Kill,size=10,prob=p.pred,log=TRUE))
}
SizePredP <- mle2(lik2,start=list(sizep1=0,sizep2=1,sizep3=12),
                 method="Nelder-Mead")
rfit <- mle2(lik3,start=list(a=0.4,b=0.3))
prfit <- mle2(lik4,start=list(a=0.4,b=0.3,g=1))
trfit <- mle2(lik5,start=list(a=0.4,b=0.3,t=8))
##
op = par(lwd=2,bty="l",las=1,mgp=c(2.5,1,0),cex=1.5)
sizeplot(TBL,Kill,xlab="Tadpole Size (TBL in mm)",ylab="Number killed",
         xlim=c(0,40),scale=0.7)
with(as.list(coef(SizePredP)),curve(logist2(x,sizep1,sizep2,sizep3)*10,add=TRUE))
with(as.list(coef(rfit)),curve(ricker(x,a,b)*10,add=TRUE,lty=2))
with(as.list(coef(prfit)),curve(powricker(x,a,b,g)*10,add=TRUE,lty=3))
## with(as.list(coef(trfit)),curve(tricker(x,a,b,t)*10,add=TRUE,lty=4))
par(xpd=NA)
legend(20,5,c("Ricker","power-Ricker","modified logistic"),
       lty=c(2,3,1),lwd=2,bty="n")
par(xpd=FALSE)
par(op)
detach(ReedfrogSizepred)


###################################################
### chunk number 3: 
###################################################
op = par(lwd=2,bty="l",las=1,mgp=c(2.5,1,0),cex=1.5)
curve(2*exp(-x/2),from=0,to=7,ylim=c(0,2),ylab="",xlab="")
curve(2*exp(-x),add=TRUE,lty=4)
curve(x*exp(-x/2),add=TRUE,lty=2)
curve(2*x*exp(-x/2),add=TRUE,lty=3)
par(cex=1)
text(0.4,1.9,expression(paste(y==2*e^{-x/2})),adj=0)
text(2,1.55,expression(paste(y==x*e^{-x/2})))
text(2.4,0.84,expression(paste(y==2*x*e^{-x/2})),adj=0)
text(2.8,0,expression(paste(y==2*e^{-x})))
par(op)


###################################################
### chunk number 4: 
###################################################
xvec = seq(0,7,length=100)


###################################################
### chunk number 5: 
###################################################
ricker = function(x,a=1,b=1) {
  a*x*exp(-b*x)
}
yvals = sapply(xvec,ricker)


###################################################
### chunk number 6: 
###################################################
mortrisk = function(N,size,H=0.84) {
  a <- attackrate(size)
  a/(1+a*N*H)
}


###################################################
### chunk number 7: 
###################################################
cspp <- coef(SizePredP)
## larval density given as 125/m^2,
## no predators given as 2/tank=25/m^2 --
## says larval density was 10/tank
logist3 <- function(N,size,sizep1=cspp["sizep1"],
                    sizep2=cspp["sizep2"],
                    sizep3=cspp["sizep3"],
                    N.s=10,H=0.84,pred.s=2,days.s=3) {
  gamma.s <- logist2(size,sizep1,sizep2,sizep3)/(pred.s*days.s)
  alpha.sd <- gamma.s/(1-gamma.s*N.s*H)
  alpha.sd/(1+alpha.sd*N*H)
}                 
## for N=N.s should give size-specific pred. (gamma.s) back
## (G/(1-GNH))/(1+GNH/(1-GNH)) =
## (G/(1-GNH))/((1-GNH+GNH)/(1-GNH)) =
## (G/(1-GNH))/(1/(1-GNH)) = G.
op = par(lwd=1,bty="l",las=1,mgp=c(2.5,1,0),cex=1.5)
par(mar=c(0.5,2.5,0.5,1))
curve3d(logist3(x,y),from=c(10,0),to=c(40,30),ticktype="detailed",
        theta=50,r=1e9,xlab="Density",ylab="Size",zlab="")
par(srt=90,xpd=NA)
text(-1.8e-9,2.4e-10,"Mortality risk",srt=90)
##plot(logist3(10,1:40),type="b")
## points(logist2(1:40,sizep1=cspp["sizep1"],
##                     sizep2=cspp["sizep2"],
##                     sizep3=cspp["sizep3"])/6,type="b",col=2,lwd=2)
#curve3d(ricker.exp.2d(x,y,b2=0.5),from=c(0,0),to=c(8,8),theta=225)
par(op)


###################################################
### chunk number 8: 
###################################################
op = par(mfrow=c(1,2))
par(lwd=2,bty="l",las=1,mgp=c(2.5,1,0),cex=1.5)
curve(exp(-x),axes=FALSE,from=0,to=3,xaxs="i",
      xlab="",ylab="")
mgp = par("mgp")
jot <- 0.015
ylocs1 = c(2^(-1:-3))
ylocs2 = exp(-1:-3)
ylocs = c(ylocs1,ylocs2)
xlocs1 = -log(ylocs1)
xlocs2 = -log(ylocs2)
xlocs = -log(ylocs)
n = length(ylocs)
n1 = length(ylocs1)
n2 = length(ylocs2)
axis(side=2,at=c(1,ylocs),labels=rep("",length(ylocs)+1))
ynudge <- c(0,0,0,+0.03,0,-0.05,-0.03)
par(xpd=NA)
axis(side=2,at=c(1,ylocs)+ynudge,tck=0,
     labels=expression(a,a/2,a/4,a/8,a/e,a/e^2,a/e^3),cex.axis=0.75)
par(xpd=FALSE)
par(mgp=c(3.5,0.75,0))
xnudge = c(0,0,0.19,0,-0.06,0)
axis(side=1,at=c(0,xlocs),labels=rep("",length(xlocs)+1))
axis(side=1,at=c(0,xlocs+xnudge),tck=0,
     labels=expression(0,frac(log(2),b),2*frac(log(2),b),3*frac(log(2),b),
       frac(1,b),frac(2,b),frac(3,b)),cex.axis=0.6,padj=0.5)
segcols = c("black","black")
seglty = c(3,2)
segments(rep(0,n1),ylocs1,xlocs1,ylocs1,lty=seglty[1],col=segcols[1])
segments(xlocs1,rep(0,n1),xlocs1,ylocs1,lty=seglty[1],col=segcols[1])
segments(rep(0,n2),ylocs2,xlocs2,ylocs2,lty=seglty[2],col=segcols[2])
segments(xlocs2,rep(0,n2),xlocs2,ylocs2,lty=seglty[2],col=segcols[2])
box()
## corner.label2("a","topleft",inset=0.025)
legend("topright",
       c("half-lives",
         "e-folding\nsteps"),
       lty=seglty,col=segcols,y.intersp=0.8,bty="n")
text(2,0.45,expression(f(x)==a*e^{-b*x}),cex=1.5)
###
par(mgp=c(2.5,1,0))
curve(x/(1+x),axes=FALSE,from=0,to=8,xaxs="i",
      yaxs="i",
      xlab="",ylab="",ylim=c(0,1.02))
xlocs <- 1:5
ylocs <- xlocs/(1+xlocs)
axis(side=1,at=xlocs,
     labels=expression(b,2*b,3*b,4*b,5*b),cex.axis=0.75)
axis(side=2,at=c(ylocs,1),labels=rep("",length(ylocs)+1))
ynudge <- c(0,0,0,+0.01,+0.04)
axis(side=2,at=c(ylocs+ynudge,1),tck=0,
     labels=expression(a/2,2*a/3,3*a/4,4*a/5,5*a/6,a),cex.axis=0.75)
par(xpd=NA)
text(7.25,0.5,expression(f(x)==frac(a*x,b+x)),cex=1.5)
par(xpd=FALSE)
abline(h=1,lty=2)
n <- length(xlocs)
segments(rep(0,n),ylocs,xlocs,ylocs,lty=3,col="black")
segments(xlocs,rep(0,n),xlocs,ylocs,lty=3,col="black")
box()
## corner.label2("b","topleft",inset=0.025,bg="white")
par(op)


###################################################
### chunk number 9: 
###################################################
op = par(lwd=2,bty="l",las=1,mgp=c(2.5,1,0),cex=1.5)
tmpf <- function(x) ifelse(x<1,0,
             1+5*(x-1)/((x-1)+1))
curve(tmpf,from=1,to=8,xlim=c(0,8),
      ylim=c(0,6),ylab="",xlab="",axes=FALSE,
      xaxs="i",yaxs="i")
text(3.5,3,expression(f(x)==frac(a*(x-c),(b+(x-c)))+d),adj=0)
axis(side=2,at=c(1,3.5,6),labels=c(expression(d),
                          expression(d+frac(a,2)),
                          expression(d+a)))
axis(side=1,at=c(1,2),labels=c(expression(c),
                          expression(c+b)))
segments(2,0,2,3.5,lty=2)
segments(0,3.5,2,3.5,lty=2)
box()
abline(h=1,lty=3)
abline(v=1,lty=3)
abline(h=6,lty=3)
par(op)


###################################################
### chunk number 10:  eval=FALSE
###################################################
## plot(c(-3,2),c(-15,15),type="n")
## z <- list()
## for (i in 1:8) {
##   z[[i]] <- locator(1)
##   points(z[[i]])
## }
## vals <- matrix(unlist(z),ncol=2,byrow=TRUE)
## vals <- as.data.frame(vals)
## names(vals) <- c("x","y")
## plm <- lm(y~x+I(x^2)+I(x^3)+I(x^4),data=vals)
## xvec <- seq(-3,2,length=100)
## lines(xvec,predict(plm,newdata=data.frame(x=xvec)))


###################################################
### chunk number 11: 
###################################################
op = par(lwd=2,cex=1.5,las=1,bty="l")
## palette(gray(seq(0,.9,len=8))) ## set colors to gray scales
a = -0.5; b = -2; c = -0.75; d = 0.333; e = 0.1
a = -11; b = 0.0046; c = 4.6; d = 0.28; e = -0.36
a = 2; b = -1.9; c = -5.27; d = -0.119; e = 0.36
a = -1; b = 3.24; c = 5.23; d = 2.8; e = 0.38
curve(a+b*x+c*x^2+d*x^3+e*x^4,xlim=c(-3,2),
     ylim=c(-15,15),lwd=3,ylab="",axes=FALSE)
axis(side=1)
axis(side=2,at=seq(-15,5,by=5))
box()
## if (!badfonts) {
 ##  mtext(side=2,at=0,line=2.7,expression(a+b*x+c*x^2+d*x^3+e*x^4),las=0)
## }
val <- a
abline(h=a,lty=2)
curve(a+b*x,lty=3,add=TRUE)
curve(a+b*x+c*x^2,lty=4,add=TRUE)
curve(a+b*x+c*x^2+d*x^3,lty=5,add=TRUE)
points(0,a,pch=16,lwd=2,col="darkgray")
par(xpd=NA)
rect(par("usr")[1],8,par("usr")[2],par("usr")[4],col="white",
     border=NA)
legend("topleft",
       c("f(x)",expression(constant:~~ f(0)),
         expression(linear:~~f(0)+f*minute(0)*x),
         expression(quadratic:~~f(0)+f*minute(0)*x+(f*second(0)/2)*x^2),
         expression(cubic:~~f(0)+f*minute(0)*x+(f*second(0)/2)*x^2+(f*minute*second(0)/6)*x^3)),
       lty=1:5,lwd=c(3,rep(2,4)),cex=0.75,bty="n")
par(op)


###################################################
### chunk number 12:  eval=FALSE
###################################################
## d1 = D(expression(1/(1+exp(-b*x))),"t")
## d2 = D(D(expression(1/(1+exp(-b*x))),"t"),"t")
## t = 0
## r = 1
## eval(d1)
## eval(d2) 


###################################################
### chunk number 13:  eval=FALSE
###################################################
## op = par(lwd=2,cex=1.5,las=1,bty="l")
## curve(plogis(x),from=-4,to=4,xlab="x",ylab="logistic(x)",axes=FALSE)
## axis(side=1,at=c(-4,0,4))
## axis(side=2,at=c(0,0.5,1))
## abline(h=1/2,lty=2)
## abline(v=0,lty=3)
## abline(a=1/2,b=1/4,lty=3)
## box()
## par(op)


###################################################
### chunk number 14:  eval=FALSE
###################################################
## time = seq(1900,2000,by=10)
## pop = c(75.995,91.972,105.711,123.203,131.669,
##      150.697,179.323,203.212,226.505,249.633,281.422)
## plot(time,pop,xlim=c(1900,2020),ylim=c(0,400))
## sct = (time-1950)/50
## orders = c(1:4,8,9)
## polynoms = lapply(orders,
##    function(i) lm(pop~poly(sct,i,raw=TRUE)))
## ptime = 1900:2020
## invisible(mapply(function(m,i) lines(ptime,predict(m,newdata=data.frame(sct=(ptime-1950)/50)),
##              lty=i),polynoms,1:length(orders)))
##              
## pfun <-  function(z) {
##      predict(polynoms[[5]],
##            newdata=data.frame(sct=rep((z-1950)/50,20)))[20]
## }
## u1 = uniroot(pfun, interval=c(1990,2030))$root
## yr = round(u1)
## fracyr = u1 %% 1
## month = floor(fracyr*12)
## t2 = (fracyr-jdate/365)*365*24
## hrs = floor(t2*24)
## mins = (t2-hrs/24)*24*60


###################################################
### chunk number 15:  eval=FALSE
###################################################
## par(mfrow=c(2,2),mar=c(0,0,0,0),lwd=2,cex=1.5)
## curve(ifelse(x<3,0,3),from=0,to=10,ylim=c(0,8),axes=FALSE,ann=FALSE)
## curve(ifelse(x<5,2,8),add=TRUE,lty=2)
## corner.label("thresholds",adj=0)
## box(col="gray")
## ##
## curve(ifelse(x<3,1,1+0.75*(x-3)),from=0,to=10,ylim=c(0,8),axes=FALSE,ann=FALSE)
## curve(ifelse(x<5,x,5),add=TRUE,lty=2)
## corner.label("hockey sticks",adj=0)
## box(col="gray")
## ##
## curve(ifelse(x<4,2*x,8-3*(x-4)),from=0,to=10,ylim=c(0,10),axes=FALSE,ann=FALSE)
## curve(ifelse(x<6,0.5*x,3+1.5*(x-6)),add=TRUE,lty=2)
## corner.label("general",adj=0)
## box(col="gray")
## ##
## ## splines?
## x1 <- 1:6
## y1 <- c(0,2,4,1,2,3)
## s1 <- splinefun(x1,y1)
## curve(s1,axes=FALSE,ann=FALSE,from=0.5,to=6.5)
## points(x1,y1,pch=16)
## y2 <- c(1,1.5,2,2.1,2.2,2.3)
## points(x1,y2,pch=1)
## s2 <- splinefun(x1,y2)
## curve(s2,add=TRUE,lty=2)
## corner.label("splines",adj=0)
## box(col="gray")


###################################################
### chunk number 16: 
###################################################
par(mfrow=c(2,2),mar=c(0,0,0,0),lwd=2,cex=1.5)
eqcex=0.9
curve(ifelse(x<3,0,3),from=1,to=9,ylim=c(-1,5),
      xlim=c(0,11),axes=FALSE,ann=FALSE,type="s")
text(0,0,expression(a[1]))
text(10,3,expression(a[2]))
segments(3,-0.5,3,3.5,col="gray",lty=2)
text(3,-1,expression(s[1]))
corner.label("threshold:",adj=0)
corner.label(expression(paste(f(x)==a[1]," if ",x<s[1])),adj=0,yoff=0.13,cex=eqcex)
corner.label(expression(paste(phantom(f(x))==a[2]," if ",x>s[1])),adj=0,yoff=0.21,cex=eqcex)
box(col="gray")
##
##curve(ifelse(x<3,1,1+0.75*(x-3)),from=0,to=10,ylim=c(0,8),axes=FALSE,ann=FALSE)
curve(ifelse(x<5,x,5),from=0,to=10,xlim=c(-1,12),ylim=c(0,8),axes=FALSE,ann=FALSE)
corner.label("hockey stick:",adj=0)
corner.label(expression(paste(f(x)==a*x," if ",x<s[1])),adj=0,yoff=0.13,cex=eqcex)
corner.label(expression(paste(phantom(f(x))==a*s[1]," if ",x>s[1])),adj=0,yoff=0.21,cex=eqcex)
segments(5,0.5,5,5.5,col="gray",lty=2)
segments(c(1,2),c(1,1),c(2,2),c(1,2),col="gray")
text(2.5,1.5,"a")
text(5,0,expression(s[1]))
text(11,5,expression(a*s[1]))
box(col="gray")
##
a=2
s1=4
b=0.5
curve(ifelse(x<s1,a*x,(a*s1)-b*(x-s1)),from=0,to=20,
      ylim=c(0,12),axes=FALSE,ann=FALSE)
#curve(ifelse(x<6,0.5*x,3+1.5*(x-6)),add=TRUE,lty=2)
corner.label("general piecewise linear:",adj=0)
corner.label(expression(paste(f(x)==a*x," if ",x<s[1])),adj=0,yoff=0.13,cex=eqcex)
corner.label(expression(paste(phantom(f(x))==a*s[1]-b*(x-s[1])," if ",x>s[1])),
             adj=0,yoff=0.21,cex=eqcex)
segments(s1,0.5,s1,9,col="gray",lty=2)
segments(c(1,2),c(a,a),c(2,2),c(a,2*a),col="gray")
x1=10
x2=12
y1 = a*s1-b*(x1-s1)
y2 = a*s1-b*(x2-s1)
segments(c(x1,x1),c(y1,y2),
         c(x1,x2),c(y2,y2),col="gray")
text(x1-0.2,(y1+y2)/2,"-b",adj=1)
text(2.5,3,"a")
text(4,0,expression(s[1]))
box(col="gray")
##
## splines?
x1 <- 1:6
y1 <- c(0,2,4,1,2,3)
s1 <- splinefun(x1,y1)
curve(s1,axes=FALSE,ann=FALSE,from=0.5,to=6.5,ylim=c(0,5))
points(x1,y1,pch=16)
#y2 <- c(1,1.5,2,2.1,2.2,2.3)
#points(x1,y2,pch=1)
#s2 <- splinefun(x1,y2)
#curve(s2,add=TRUE,lty=2)
corner.label("splines:",adj=0)
corner.label("f(x) is complicated",adj=0,yoff=0.13,cex=eqcex)
box(col="gray")


###################################################
### chunk number 17: 
###################################################
par(mfrow=c(2,2),mar=c(0,0,0,0),lwd=2,cex=1.5,xaxs="i",yaxs="i")
a=2
b=1
curve(a/(b+x),from=0,to=15,ylim=c(0,3),xlim=c(-2,15),axes=FALSE,ann=FALSE)
corner.label("hyperbolic:",adj=0)
corner.label(expression(f(x)==frac(a,b+x)),adj=0,yoff=0.17,cex=eqcex)
text(-1.5,a/b,expression(frac(a,b)),adj=0)
y1=a/(2*b)
text(-1.5,y1,expression(frac(a,2*b)),adj=0)
segments(0.2,y1,3,y1,lty=2,col="gray")
x1=b
segments(b,0.25,b,y1+0.25,lty=2,col="gray")
text(x1,0.1,"b")
box(col="gray")
##
a=1
b=1
curve(a*x/(b+x),from=0,to=10,ylim=c(0,b*1.4),
      xlim=c(0,11),axes=FALSE,ann=FALSE)
corner.label("Michaelis-Menten:",adj=0)
corner.label(expression(f(x)==frac(a*x,b+x)),adj=0,yoff=0.17,cex=eqcex)
text(10.9,a,"a",adj=1)
segments(0,a,10,a,lty=2,col="gray")
x1=0.02
x2=0.4
y1=a*x1/(b+x1)
y2=a*x2/(b+x2)
segments(c(x1,x2),c(y1,y1),c(x2,x2),c(y1,y2),col="gray")
text(x2+0.04,(y1+y2)/2,expression(frac(a,b)),adj=0,cex=0.7)
x1=b
y1=a/2
segments(0,y1,4,y1,lty=2,col="gray")
text(4.1,y1,adj=0,expression(frac(a,2)))
segments(x1,0.15,x1,y1,lty=2,col="gray")
text(x1+0.1,0.1,"b",adj=0)
box(col="gray")
##
a=1
b=1
curve(a*x^2/(b^2+x^2),from=0,to=5.5,ylim=c(0,b*1.4),
      xlim=c(0,6),axes=FALSE,ann=FALSE)
corner.label("Holling type III:",adj=0)
corner.label(expression(f(x)==frac(a*x^2,b^2+x^2)),adj=0,yoff=0.18,cex=eqcex)
segments(0,a,5.5,a,lty=2,col="gray")
text(5.9,a,"a",adj=1)
text(b,0.05,"b")
segments(b,0.1,b,a/2,lty=2,col="gray")
segments(0,a/2,b+0.5,a/2,lty=2,col="gray")
text(b+0.55,a/2,expression(frac(a,2)),adj=0)
box(col="gray")
##
a=1
b=1
c=-1.5
curve(a*x^2/(b+c*x+x^2),from=0,to=5.5,ylim=c(0,3.5),
      xlim=c(0,6),axes=FALSE,ann=FALSE)
corner.label("Holling type IV (c<0):",adj=0)
corner.label(expression(f(x)==frac(a*x^2,b+c*x+x^2)),adj=0,yoff=0.18,cex=eqcex)
segments(0,a,5.5,a,lty=2,col="gray")
text(5.9,a,"a",adj=1)
segments(-2*b/c,0.9,-2*b/c,2.8,col="gray",lty=2)
text(-2*b/c,0.35,expression(frac(-2*b,c)))
## text(b,0.05,"b")
## segments(b,0.1,b,a/2,lty=2,col="gray")
## segments(0,a/2,b+0.5,a/2,lty=2,col="gray")
## text(b+0.55,a/2,expression(frac(a,2)),adj=0)
box(col="gray")


###################################################
### chunk number 18: 
###################################################
op <- par(mfrow=c(2,2),mar=c(0,0,0,0),lwd=2,cex=1.5,xaxs="i",yaxs="i")
a=1
b=1
curve(a*exp(-b*x),from=0,to=3,ylim=c(0,1.3),axes=FALSE,ann=FALSE,
      xlim=c(-0.5,3))
corner.label("negative exponential:",adj=0)
corner.label(expression(paste(f(x)==a*e^{-b*x})),adj=0,yoff=0.13,cex=eqcex)
text(-0.3,1,adj=1,"a")
segments(-0.25,a*exp(-1),1.4,a*exp(-1),col="gray",lty=2)
text(-0.3,a*exp(-1),adj=1,expression(frac(a,e)))
segments(1/b,0.26,1/b,a*exp(-1),col="gray",lty=2)
text(1/b,0.15,expression(frac(1,b)))
box(col="gray")
##
a=1
b=1
xmax=7
curve(a*(1-exp(-b*x)),from=0,to=xmax,ylim=c(0,1.5),axes=FALSE,ann=FALSE)
curve(a*x/(b+x),from=0,to=xmax,add=TRUE,col="gray")
text(xmax*0.8,0.75,"M-M",col="gray")
corner.label("monomolecular:",adj=0)
corner.label(expression(f(x)==a*(1-e^{-b*x})),
             adj=0,yoff=0.13,cex=eqcex)
segments(0.3,a,10,a,col="gray",lty=2)
text(0,a,adj=0,"a")
x1=0.05
x2=0.3
y1=a*(1-exp(-b*x1))
y2=a*(1-exp(-b*x2))
segments(c(x1,x2),c(y1,y1),
         c(x2,x2),c(y1,y2),col="gray")
text(x2+0.1,(y1+y2)/2,adj=0,expression(a*b))
box(col="gray")
##
a=1
b=1
curve(a*x*exp(-b*x),from=0,to=6,ylim=c(0,0.5),axes=FALSE,ann=FALSE)
corner.label("Ricker:",adj=0)
corner.label(expression(f(x)==a*x*e^{-b*x}),
             adj=0,yoff=0.13,cex=eqcex)
x1=0.01
x2=0.2
y1=a*x1*exp(-b*x1)
y2=a*x2*exp(-b*x2)
segments(c(x1,x2),c(y1,y1),
         c(x2,x2),c(y1,y2),col="gray")
text(x2+0.1,(y1+y2)/2,adj=0,expression(a))
y1 = a/b*exp(-1)
segments(1/b,0.1,1/b,y1,col="gray",lty=2)
text(1/b,0.05,expression(frac(1,b)))
segments(0,y1,1/b+1,y1,col="gray",lty=2)
text(1/b+1+0.1,y1,adj=0,expression(frac(a,b)*e^{-1}))
box(col="gray")
##
a=-5
b=1
curve(plogis(a+b*x),from=0,to=10,ylim=c(0,1.2),axes=FALSE,ann=FALSE)
corner.label("logistic:",adj=0)
corner.label(expression(f(x)==frac(e^{a+b*x},1+e^{a+b*x})),
             adj=0,yoff=0.17,cex=eqcex)
x1=-a-0.5
x2=-a+0.5
y1=plogis(a+b*x1)
y2=plogis(a+b*x2)
segments(c(x1,x2),c(y1,y1),
         c(x2,x2),c(y1,y2),col="gray")
text(x2+0.1,(y1+y2)/2,adj=0,expression(frac(b,4)))
segments(-a/b,0.25,-a/b,0.5,col="gray",lty=2)
text(-a/b-0.2,0.11,expression(-frac(a,b)),adj=0.5)
segments(-a/b-1,0.5,-a/b,0.5,col="gray",lty=2)
text(-a/b-1.2,0.5,adj=1,expression(frac(1,2)))
segments(10,1,7,1,col="gray",lty=2)
text(6.9,1,adj=1,1)
box(col="gray")


###################################################
### chunk number 19: 
###################################################
op <- par(mfrow=c(2,2),mar=c(0,0,0,0),lwd=2,cex=1.5,xaxs="i",yaxs="i")
a=1
b=2.5
curve(a*x^b,from=0,to=3,ylim=c(0,20),axes=FALSE,ann=FALSE)
a2=10
b2=0.5
curve(a2*x^b2,from=0,to=3,add=TRUE,lty=2)
a3=2
b3=-1
curve(a3*x^b3,from=0,to=3,add=TRUE,lty=3)
corner.label("power laws:",adj=0)
corner.label(expression(f(x)==a*x^b),
             adj=0,yoff=0.13,cex=eqcex)
text(1.95,15,expression(paste(0<b,""<1)),adj=1)
text(2.3,7.1,expression(b>1),adj=0)
text(2,2.4, expression(b<0))
box(col="gray")
##
a <- 1
k <- 1
d <- 2/3
xmax=15
curve(a*(1 - exp(-k*(a-d)*x))^(1/(1-d)),from=0,to=xmax,ylim=c(0,1.3),
      ann=FALSE,axes=FALSE)
corner.label("von Bertalanffy:",adj=0,yoff=0.06)
corner.label(expression(f(x)==a*(1-e^{-k*(a-d)*x})^(1/(1-d))),adj=0,yoff=0.145,cex=eqcex)
segments(xmax-5,a,xmax,a,lty=2,col="gray")
text(xmax-5.5,a,adj=1,"a")
box(col="gray")
## Shepherd, Hassell
##
curve(x/(1+x)^1.8,from=0,to=3,ylim=c(0,0.45),axes=FALSE,ann=FALSE,lty=2)
curve(0.6*x/(1+x^1.8),add=TRUE,lty=1)
curve(0.8*x*exp(-x),add=TRUE,col="gray")
text(2,0.13,"Ricker",col="gray",adj=0)
corner.label("Shepherd, Hassell:",adj=0)
corner.label(expression(list(f(x)==frac(a*x,b+x^c),f(x)==frac(a*x,(b+x)^c))),adj=0,yoff=0.17,cex=eqcex)
text(2.8,0.3,"H")
text(2.8,0.2,"S")
box(col="gray")
##
xi <- 0.9
alpha <- 1
pmax <- 20
curve(1/(2*xi)*(alpha*x+pmax-sqrt((alpha*x+pmax)^2-4*xi*alpha*x*pmax)),from=0,to=200,ylim=c(0,25),
      ann=FALSE,axes=FALSE)
curve(20*x/(x+10),add=TRUE,col="gray")
corner.label("non-rectangular\nhyperbola:",adj=0,yoff=0.09)
text(154,17.7,"M-M",col="gray")
box(col="gray")


###################################################
### chunk number 20:  eval=FALSE
###################################################
## library(odesolve)
## tlfun = function(parms) { lsoda(y=0.1,times=seq(0,10,by=0.1),parms=parms,
##        func=function(t,y,parms) { r=parms[1]; K=parms[2]; theta=parms[3]; 
##                          list(r*y*(1-(y/K)^theta),NULL) }) }
## x1 = tlfun(c(r=1,K=1,theta=1))
## x2 = tlfun(c(r=1,K=1,theta=2))
## x3 = tlfun(c(r=1,K=1,theta=0.5))
## matplot(x1[,1],cbind(x1[,2],x2[,2],x3[,2]),lty=1:3,col=1,type="l",
##     xlab="Time",ylab="Pop. density")
## text(c(2.25,3.38,5.5),c(0.8,0.84,0.75),
##      c(expression(theta==2),expression(theta==1),expression(theta==0.5)))


###################################################
### chunk number 21: 
###################################################
## Richards curve:
richards <- function(x,a,k,delta,gamma) {
  if (delta==1) {
    k*exp(-exp(-k*(x-gamma)))
  } else {
    a*(1+(delta-1)*exp(-k*(x-gamma)))^(1/(1-delta))
  }
}


###################################################
### chunk number 22: 
###################################################
curve(2*x/(1+x))


###################################################
### chunk number 23: 
###################################################
micmen <- function(x,a=2,b=1) {
  a*x/(b+x)
}


###################################################
### chunk number 24: 
###################################################
curve(micmen(x),from=0,to=8,ylim=c(0,10))
curve(micmen(x,b=3),add=TRUE,col=2)
curve(micmen(x,a=8),add=TRUE,col=3)
abline(h=8)


###################################################
### chunk number 25: 
###################################################
xvec <- seq(0,10,by=0.1)


###################################################
### chunk number 26: 
###################################################
curve(ifelse(x<5,1,2),from=0,to=10)


###################################################
### chunk number 27: 
###################################################
curve(ifelse(x<5,1+x,6-3*(x-5)),from=0,to=10)


###################################################
### chunk number 28: 
###################################################
curve(ifelse(x<5,1+x,
             ifelse(x<8,
                    6-3*(x-5),
                    -3+2*(x-8))),from=0,to=10)


###################################################
### chunk number 29: 
###################################################
D(expression(log(x)),"x")
D(expression(x^2),"x")
logist <- expression(exp(x)/(1+exp(x)))
dfun <- deriv(logist,"x",function.arg=TRUE)
xvec <- seq(-4,4,length=40)
y <- dfun(xvec)
plot(xvec,y)
lines(xvec,attr(y,"grad"))


###################################################
### chunk number 30: 
###################################################
d1 <- D(expression(a*x/(b+x)),"x")
d1
eval(d1,list(a=2,b=1,x=3))


