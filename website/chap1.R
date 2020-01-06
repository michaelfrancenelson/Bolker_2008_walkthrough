###################################################
### chunk number 1: 
###################################################
## badfonts <- .Platform$OS.type=="unix"
badfonts <- FALSE
library(Hmisc,warn.conflicts=FALSE)
source("chapskel.R")


###################################################
### chunk number 2: 
###################################################
## read in and attach seed predation data
library(emdbook)
data(SeedPred)
seedsize <- c(0.11,0.03,0.21,0.16,0.35,0.01,0.52,0.27)
ss.levels <- levels(SeedPred$species)[order(seedsize)]
SeedPred$species <- factor(SeedPred$species,levels=ss.levels)
seedsub = subset(SeedPred,seeds>0 & !is.na(taken))
attach(seedsub,warn.conflicts=FALSE)
## grab useful? libraries
invisible(require(lattice,quietly=TRUE))
invisible(require(plotrix,quietly=TRUE,warn.conflicts=FALSE))
## compute none/any values for number of seeds taken by species
t1 <- table(taken>0,species)
## total number of obs per species
tot <- table(species)
## number of "successes" and proportion for all species
succ <- t1[2,]
props <- succ/tot
## restrict to just the first two species
spp <- c("pol","psd")
succ <- succ[spp]
tot <- tot[spp]
vals <- t1[2:1,spp] ## reverse order
## MLE estimates of binomial p
p.mle <- succ/tot
ptot.mle <- sum(succ)/sum(tot)
obs.pratio <- p.mle["pol"]/p.mle["psd"]
## let's be frequentist: Fisher exact test
F1 <- fisher.test(vals)
F1.1 <- fisher.test(vals,alternative="greater")
## likelihood: restricted and full likelihood
pr <- sum(succ)/sum(tot)
L0 <- dbinom(sum(succ),size=sum(tot),prob=pr,log=TRUE)
L0 <- sum(dbinom(succ,size=tot,prob=pr,log=TRUE))  ## ??? sum(log-likelihoods) neq log(pooled data)???
L0x <- dbinom(succ,size=tot,prob=pr,log=TRUE)
L0bx <- lchoose(tot,succ)+succ*log(pr)+log(1-pr)*(tot-succ)
L0b <- lchoose(sum(tot),sum(succ))+sum(succ)*log(pr)+log(1-pr)*(sum(tot)-sum(succ))
## correct!
L1 <- sum(dbinom(succ,size=tot,prob=succ/tot,log=TRUE))
g1 <- glm(t(vals)~factor(1:2),family="binomial")
g2 <- glm(t(vals)~1,family="binomial")
## anova(g1,g2,test="Chisq")
## ??? logLik(g1) matches L1
## logLik(g2) = -22.23, L0 = -2.84
## plogis(coef(g2))
## plogis(coef(g1)[1])
## plogis(sum(coef(g1)))
## likelihood ratio test:
LRT.p <- pchisq(2*(L1-L0),df=1,lower.tail=FALSE)
OR.to.PR <- function(OR,p1) {
  O1 <- p1/(1-p1)
  O2 <- O1*OR
  p2 <- O2/(1+O2)
  p2/p1
}
pr.conf <- OR.to.PR(F1$conf,p.mle[2])

## calculate num taken from odds ratio,
## given margins v1,v2,v3
OR.to.nums <- function(OR,vals) {
  v1 <- colSums(vals)[1] ## total pol
  v2 <- colSums(vals)[2] ## total psd
  v3 <- rowSums(vals)[1] ## total num taken
}

## calculate odds ratio given number in [1,1]
x.to.OR <- function(x,vals) {
  v1 <- colSums(vals)[1] ## total pol
  v2 <- colSums(vals)[2] ## total psd
  v3 <- rowSums(vals)[1] ## total num taken
  m <- matrix(nrow=2,ncol=2)
  m[1,1] <- x
  m[2,1] <- v1-x
  m[1,2] <- v3-x
  m[2,2] <- v2-m[1,2]
  (m[1,1]/m[2,1])/(m[1,2]/m[2,2])
}

## PR = P1/P2 = (m[1,1]/v1)/(m[1,2]/v2)
## m[1,2] = v3-m[1,1]
## PR = (m[1,1]/v1)/((v3-m[1,1])/v2)
## ((v3-m[1,1])/v2)*PR = m[1,1]/v1
## v3/v2*PR = m[1,1]*(1/v1+PR/v2)
## m[1,1] = v3/v2*PR/(1/v1+PR/v2)
## calculate number in [1,1] given prob. ratio,
## marginals
PR.to.x <- function(PR,vals) {
  v1 <- colSums(vals)[1] ## total pol
  v2 <- colSums(vals)[2] ## total psd
  v3 <- rowSums(vals)[1] ## total num taken
  v3/v2*PR/(1/v1+PR/v2)
}

ors <- sapply(0:40,x.to.OR,vals=vals)
prs <- sapply(ors,OR.to.PR,p1=vals[1,1]/colSums(vals)[1])
or2 <- pmax(ors,1/ors)
smvals <- (0:40)[ors<1 & or2>F1$estimate]
lgvals <- vals[1,1]:40


###################################################
### chunk number 3: 
###################################################
vals2 = rbind(vals,colSums(vals))
dimnames(vals2)[[1]] <- c("any taken ($t$)","none taken","total ($N$)")
latex(vals2,file="",title="",table.env=FALSE)


###################################################
### chunk number 4: 
###################################################
op <- par(cex=1.5,lwd=2,bty="l",yaxs="i")
par(mar=c(5,4,4,2)+0.1)
## stuff for Fisher's test: one-sided test
## vals
maxv <- 40
xvec <- 0:maxv
v1 <- colSums(vals)[1] ## total pol
v2 <- colSums(vals)[2] ## total psd
v3 <- rowSums(vals)[1] ## total num taken
hvec <- dhyper(xvec,v1,v2,v3,log=TRUE)/log(10)
## convert to log-10 units
hvec2 <- dhyper(xvec,v1,v2,v3)
##phyper(29,125,439,52,lower.tail=FALSE)
pval <- fisher.test(vals,alternative="greater")$p.val
pval2 <- fisher.test(vals,alternative="two.sided")$p.val
pval3 <- pval2-pval
op <- par(cex=1.5,lwd=2,bty="l",mgp=c(2.5,1,0),yaxs="i")
## barplot(hvec)
plot(xvec,hvec,type="n",xlab="Number of pol stations with any seeds taken",
     ylab="",ylim=c(min(hvec),min(hvec)+diff(range(hvec))*1.1),
      axes=FALSE,xlim=c(-5,40))
rect(xvec-0.5,rep(min(hvec),length(xvec)),xvec+0.5,hvec,
     border="black",lwd=1,
     col=rep(c("darkgray","white","darkgray"),
       c(length(smvals),
         length(xvec)-length(smvals)-length(lgvals), 
         length(lgvals))))
## plot(xvec,hvec2,type="l",log="y")
if (!badfonts) mtext("Probability",side=2,at=-10,line=3,las=0,cex=1.5)
axis(side=1,at=seq(0,40,by=10))
axis(side=2,at=c(-5,-10,-15),labels=NA)
prvals <- seq(0,10,by=2)
axis(side=3,at=PR.to.x(prvals,vals=vals),labels=prvals,
     line=0.5,cex.axis=0.8)
v0 <- vals[1,1]
obs=v0
obsPR=3.62
par(xpd=NA)
arrows(26,4.34,26,1.8,length=0.15,angle=15) ## 1.585
text(26,5.0, ## 4.82,
     paste("obs=",obsPR,sep=""),cex=0.5)
par(xpd=TRUE)
## kluge
mtext(side=3,at=PR.to.x(10,vals=vals),text="10",
      line=1.5,cex=0.8*1.5)
mtext(side=3,at=20,line=3,"Probability ratio")
par(las=1)
invisible(lapply(seq(-5,-15,by=-5),
                 function(x) { mtext(text=substitute(10^x,list(x=x)),
                                     side=2,at=x,line=par("mgp")[2],cex=1.5) }))
v0 <- vals[1,1]
x1 <- lgvals
##polygon(c(x1,rev(x1)),c(hvec[xvec>=v0],rep(par("usr")[3],length(x1))),
##        col="gray",border=NA)
x2 <- smvals
v5 <- max(smvals)
## polygon(c(x2,rev(x2)),c(hvec[xvec<=v5],rep(par("usr")[3],length(x2))),
##        col="gray",border=NA)
## lines(xvec,hvec)
box(lwd=2)
u1=par("usr")[3]
u2=par("usr")[4]
segments(v0,u1,v0,-2.7,lty=2)
segments(v5,u1,v5,u2,lty=2)
par(xpd=NA)
text(v0,1,paste("observed\n=",v0,sep=""),pos=1)
do.call("text",list(40,-8,scinot(pval,"expression",digits=2,pref="p="),
                    cex=0.5))
do.call("text",list(-2,-2,scinot(pval3,"expression",digits=2,pref="p="),
                    cex=0.5))
arrows(40,-9,36,-15,angle=15)
arrows(-2,-3,1,-6,angle=15)
par(op)


###################################################
### chunk number 5: 
###################################################
or = vals[1,1]*vals[2,2]/(vals[1,2]*vals[2,1])


###################################################
### chunk number 6: 
###################################################
tprobs <- c(0.05,0.04)
bprobs <- dbinom(sum(succ),prob=tprobs,size=sum(tot))
ml <- dbinom(sum(succ),prob=sum(succ)/sum(tot),size=sum(tot))


###################################################
### chunk number 7: 
###################################################
op <- par(mfrow=c(2,1),las=1,lwd=2,cex=1.5,bty="l")
par(mar=c(0,4,2,2)+0.2)
curve(dbinom(sum(succ),size=sum(tot),prob=x),from=0.01,to=0.2,ylim=c(0,0.06),
      axes=FALSE,
      ylab="Likelihood",xlab="")
##      xlab="Prob. any seeds\ntaken in one trial",
axis(side=2,at=seq(0,0.06,by=0.03))
box()
abline(v=sum(succ)/sum(tot),lty=2)
par(mar=c(4,4,1,2)+0.2)
corner.label2("a","topleft",inset=0.025)
curve(dbinom(sum(succ),size=sum(tot),prob=x,log=TRUE),from=0.01,to=0.2,
      xlab="",
      axes=FALSE,
      ylab="Log-likelihood",ylim=c(-66,0))
mtext(side=1,expression(paste("P(seeds taken), ",p[s])),
    at=0.1,line=2.5,cex=1.5)
axis(side=2,at=seq(0,-60,by=-20))
box()
par(xpd=NA)
corner.label2("b","topleft",inset=0.025)
par(xpd=FALSE)
axis(side=1)
par(las=1)
## invisible(lapply(seq(-5,-25,by=-10),
##        function(x) { mtext(text=substitute(10^x,list(x=x)),
##                            side=2,at=10^(x),line=par("mgp")[2],cex=1.5) }))
#mtext(side=2,at=1e-15,"Likelihood",line=3.2,las=0,cex=1.5)
#box()
abline(v=sum(succ)/sum(tot),lty=2)
par(op)


###################################################
### chunk number 8: 
###################################################
clip <- function(x) { min(1,max(0,x)) }
plik1 <- function(r) {
  optimize(f=function(x) {-sum(dbinom(x=succ,
             size=tot,prob=c(x,
                        clip(x*r)),log=TRUE))},
           interval=c(0.01,0.4))$objective
}

## use odds ratio instead of probability ratio 
plik1.OR <- function(r) {
  optimize(f=function(x) {
    or2 <- (x/(1-x))*r; p2 <- or2/(1+or2);
    -sum(dbinom(x=succ,size=tot,prob=c(x,p2),log=TRUE))},
           interval=c(0.05,0.4))$objective
}
## range of probability ratios to try
## changed range -- species 1 now has higher number taken
rvec <- seq(0.01,1,length=300)
rvec2 <- seq(1.5,7,length=100)
## negative log-likelihood
lprof <- sapply(rvec,plik1)
lprof2 <- sapply(1/rvec2,plik1)
## likelihood profile on odds ratio
lprof.OR <- sapply(rvec,plik1.OR)

fit.mle <- optim(fn=function(p) { -sum(dbinom(x=succ,size=tot,prob=c(p[1],clip(p[1]*p[2])),
                   log=TRUE)) },par=c(0.2,2))

##library(stats4)
##invisible(require(nlme,quietly=TRUE))
##library(methods)
invisible(require(bbmle,quietly=TRUE))
fit2.mle <- mle2(minuslogl=function(p,r) { -sum(dbinom(x=succ,size=tot,prob=c(p,clip(p*r)),
                   log=TRUE)) },start=list(p=0.2,r=2))
ci.mle <- confint(fit2.mle,quietly=TRUE)


###################################################
### chunk number 9: 
###################################################
lprior <- function(prob,ratio) { 1 } ## improper prior in both directions
plik2 <- function(prob,ratio) {
  exp(sum(dbinom(x=succ,size=tot,prob=c(prob,clip(prob*ratio)),log=TRUE)))
}
## odds ratio:
## prob1 + odds ratio -> new odds ratio = prob1/(1-prob1)*oddsratio
##  new probability = OR/(1+OR)
plik2.OR <- function(prob,ratio) {
  or2 <- (prob/(1-prob))*ratio  ## changed from *ratio to /ratio
  p2 <- or2/(1+or2)
  exp(sum(dbinom(x=succ,size=tot,prob=c(prob,p2),log=TRUE)))
}
pnum <- function(P) lprior(P[1],P[2])*plik2(P[1],P[2])
library(adapt)
denom <- adapt(ndim=2,lower=c(0,0),upper=c(1,10),functn=pnum,eps=0.005)$value
plik2B <- function(prob,ratio) {
  sapply(prob,plik2,ratio=ratio)
}
plik2B.OR <- function(prob,ratio) {
  sapply(prob,plik2.OR,ratio=ratio)
}
plik3 <- function(r) {
  integrate(plik2B,lower=0,upper=1,ratio=r)$value
}
plik3B <- function(p) {
  sapply(p,plik3)
}

plik3.OR <- function(r) {
  integrate(plik2B.OR,lower=0,upper=1,ratio=r)$value
}
## what probability ratio corresponds to a given odds ratio???
ratiopost <- sapply(rvec,plik3)
ratiopost.OR <- sapply(rvec,plik3.OR)


###################################################
### chunk number 10:  eval=FALSE
###################################################
## plot(rvec,ratiopost)
## par(new=TRUE)
## plot(rvec,exp(-lprof),col=2,type="l")
## abline(v=p.mle[2]/p.mle[1],lty=2)
## 
## plot(rvec,ratiopost.OR)
## par(new=TRUE)
## plot(rvec,exp(-lprofC),col=2,type="l")
## abline(v=1.9535)


###################################################
### chunk number 11: 
###################################################
densfac <- 1/sum(ratiopost)/diff(rvec)[1]
## find credible interval ...
## find lower and upper values for which prob dens = target value
lims <- function(pdens) {
  lower <- uniroot(function(p) {plik3(p)*densfac-pdens},
                   interval=c(0.0001,rvec[which.max(ratiopost)]))$root
  upper <- uniroot(function(p) {plik3(p)*densfac-pdens},
                   interval=c(rvec[which.max(ratiopost)],3))$root
  c(lower,upper)
}

## find area between target values
limarea <- function(pdens) {
  intlim <- lims(pdens)
  (integrate(plik3B,intlim[1],intlim[2])$value)*densfac
}

## find credible interval
credint <- function(level) {
  u <- uniroot(function(x) {limarea(x)-level},
               interval=c(1,1e-5))
  intlim <- lims(u$root)
  c(intlim,plik3(intlim[1])*densfac,limarea(u$root))
}

cred1 <- credint(0.95)

bayesmean <- integrate(function(x)plik3B(x)*x,0.001,0.6)$value*densfac  ## range covers plausible range
bayesmode <- optimize(plik3B,interval=c(0.001,6),maximum=TRUE)$maximum  ##
bayeshyp <- integrate(plik3B,1,Inf)$value*densfac


###################################################
### chunk number 12: 
###################################################
op <- par(las=1,lwd=2,cex=1.5,bty="l",mar=c(5,4,4,4)+0.1)
plot(rvec2,-lprof2,type="l",xlab="",ylab="Log-likelihood")
mtext(side=1,"Ratio of (pol prob)/(psd prob)",line=2,cex=1.5)
abline(h=-fit.mle$value,lty=2)
## plot(rvec2,exp(-lprof),type="l",col=2)
abline(v=p.mle["pol"]/p.mle["psd"],lty=2)
abline(v=1/ci.mle[2,],lty=3)
abline(h=logLik(fit2.mle)+c(0,-qchisq(0.95,1)/2),lty=3)
##abline(v=1,lty=4)
par(xpd=NA)
text(1/coef(fit2.mle)[2],-4.5,"MLE")
text(12,logLik(fit2.mle),"maximum\nlikelihood",adj=1)
text(12,logLik(fit2.mle)-qchisq(0.95,1)/2-0.2,"95% cutoff",adj=1)
text(1/ci.mle[2,],-4.5,c("upper","lower"))
##text(1,-17,"null",cex=1.5)
##arrows(c(1.3,2.0),-12,ci.mle[2,],-12,lwd=2)
par(xpd=FALSE)


###################################################
### chunk number 13: 
###################################################
## shade interior of credible interval?
load("ch1bayes.RData")
op <- par(cex=1.5,lwd=2,bty="l",mgp=c(2.5,1,0),las=1,yaxs="i")
plot(rvec,ratiopost.norm,type="l",
     xlab="Ratio of (pol prob)/(psd prob)",
     ylab="Probability density",xlim=c(0.9,10),ylim=c(0,0.45))
abline(h=cred1["p"])
##abline(v=cred1[1:2],lty=2)
##abline(v=1,lty=3)
abline(v=bayesmode,lty=1)
abline(v=bayesmean,lty=3)
labht <- 0.45
## text(bayesmean,1.0,"mean",cex=1.5)
## text(bayesmode,0.9,"mode",cex=1.5)
rcred <- cred1["p"]
x1 <- c(rvec[rvec<cred1["lwr"]],cred1["lwr"])
polygon(c(x1,rev(x1)),c(ratiopost.norm[rvec<cred1[1]],cred1["p"],rep(0,length(x1))),
        col="gray",border=NA)
x2 <- c(cred1["upr"],rvec[rvec>cred1["upr"]])
polygon(c(x2,rev(x2)),c(cred1["p"],ratiopost.norm[rvec>cred1["upr"]],rep(0,length(x2))),
        col="gray",border=NA)
lines(rvec,ratiopost.norm)  ## redraw
##text(1,1,"null",cex=1.5)
##text(1.7,0.1,"credible\ninterval",cex=1.5)
##arrows(c(1.4,2.0),0.08,cred1[1:2],0.08,lwd=2)
arrht <- 0.03
arrows(cred1[1],arrht,cred1[2],arrht,angle=15,length=0.18,code=3)
text(mean(cred1[1:2]),0.015,"credible interval",cex=0.65)
par(xpd=NA)
labgap <- 0.5
text(bayesmode-labgap,labht,"mode",adj=1)
arrows(bayesmode-labgap*0.8,labht,bayesmode,labht,angle=15,length=0.15)
text(bayesmean+labgap,labht,"mean",adj=0)
arrows(bayesmean+labgap*0.8,labht,bayesmean,labht,angle=15,length=0.15)
box()
par(op)


###################################################
### chunk number 14: 
###################################################
detach(seedsub)


###################################################
### chunk number 15: 
###################################################
2*8
sqrt(25)


###################################################
### chunk number 16: 
###################################################
x = sqrt(36)
x


###################################################
### chunk number 17:  eval=FALSE
###################################################
## install.packages("plotrix")


###################################################
### chunk number 18:  eval=FALSE
###################################################
## install.packages(c("ellipse","plotrix"))


###################################################
### chunk number 19:  eval=FALSE
###################################################
## install.packages("plotrix",repos=NULL)


###################################################
### chunk number 20:  eval=FALSE
###################################################
## mypkgdir="c:/Documents and Settings/Bolker/Desktop/Rpkgs"
## # options(repos="http://cran.us.r-project.org")
## install.packages("plotrix",destdir=mypkgdir,lib=mypkgdir)


###################################################
### chunk number 21:  eval=FALSE
###################################################
## library(plotrix,lib=mypkgdir)


###################################################
### chunk number 22: 
###################################################
## TO DO: ignore packages? (weaver)
r1 <- system("grep \"require(\" chap[0-9].Rnw chap[0-9][0-9].Rnw | sed -e 's/^.*require(\\([^),]*\\).*/\\1/'",intern=TRUE)
r2 <- system("grep \"library(\" chap[0-9].Rnw chap[0-9][0-9].Rnw | sed -e 's/^.*library(\\([^),]*\\).*/\\1/'",intern=TRUE)
r3 <- c(r1,r2) 
r3 <- unique(sort(r3[-grep("\\\\",r3)]))
i1 <- installed.packages()
omitpkgs <- rownames(i1)[!is.na(i1[,"Priority"])]
omitpkgs <- c(omitpkgs,"Hmisc","weaver","pkg")
r4 <- r3[!(r3 %in% omitpkgs)]
## OVERRIDE: add chron and gtools/gdata
r4 <- c("adapt","bbmle","chron","coda","ellipse","emdbook","gplots","gtools","gdata",
     "MCMCpack","odesolve","plotrix","R2WinBUGS","reshape","rgl","scatterplot3d")
cat(r4,"\n",fill=50)


###################################################
### chunk number 23:  eval=FALSE
###################################################
## help.start()


###################################################
### chunk number 24: 
###################################################
set.seed(101)


###################################################
### chunk number 25: 
###################################################
frogs = c(1.1,1.3,1.7,1.8,1.9,2.1,2.3,2.4,2.5,2.8,
      3.1,3.3,3.6,3.7,3.9,4.1,4.5,4.8,5.1,5.3)
tadpoles = rnorm(n=20,mean=2*frogs,sd=0.5)


###################################################
### chunk number 26: 
###################################################
tadpoles


###################################################
### chunk number 27: 
###################################################
plot(frogs,tadpoles)
abline(a=0,b=2)


###################################################
### chunk number 28: 
###################################################
plot(frogs,tadpoles)
abline(a=0,b=2)


###################################################
### chunk number 29: 
###################################################
log_tadpoles = log(tadpoles)
plot(frogs,log_tadpoles)


###################################################
### chunk number 30: 
###################################################
n=1:length(frogs)
plot(n,frogs)


###################################################
### chunk number 31: 
###################################################
h = rev((1:6)/6)
xr = seq(0,0.9,length=20)
xpos = 1.2
op = par(xpd=TRUE,mar=c(0,0,0,5))
plot(0:1,range(h),type="n",axes=FALSE,ann=FALSE)
points(xr,rep(h[1],20),pch=1:20)
text(xpos,h[1],"pch: point type",adj=1)
palette(gray((0:19)/19))
points(xr,rep(h[2],20),col=1:20,pch=16)
palette("default")
text(xpos,h[2],"col: point color",adj=1)
points(xr,rep(h[3],20),cex=seq(0.01,2,length=20))
text(xpos,h[3],"cex: point size",adj=1)
text(xr,rep(h[4],20),c(letters[1:5],LETTERS[1:5],0:9))
text(xpos,h[4],"text",adj=1)
##text(seq(0.1,0.9,length=20),rep(0.6,20),letters[1:20],cex=seq(0.01,2,length=20))
nl = 6
nlsp = 0.1
lend = 0.95
##segments(seq(0,0.9,length=nl),rep(h[5],2*nl),
##         seq(1/nl,1,length=nl),rep(h[5],2*nl),
## want to go from 0 to lend, in nl segments,
## using a fraction fr of the space for interspaces
## each segment + space takes up lend/nl
## 
fr = 0.7
seg1 = seq(0,by=lend/nl,length=nl)
seg2 = seq(lend/nl,to=lend,by=lend/nl)-(lend*(1-fr)/nl)
segments(seg1,rep(h[5],nl),
         seg2,
         rep(h[5],nl),
         lty=1:nl,lwd=3)
text(xpos,h[5],"lty: line type",adj=1)
nv = 10
segments(seq(0,1/(1+nlsp),length=nv),rep(h[6],nv),
         seq(0.05,0.95,length=nv),rep(h[6],nv),
         lwd=1:10)
text(xpos,h[6],"lwd: line width",adj=1)
par(op)


###################################################
### chunk number 32: 
###################################################
mean(tadpoles)
sd(tadpoles)
summary(tadpoles)


###################################################
### chunk number 33: 
###################################################
cor(frogs,tadpoles)


###################################################
### chunk number 34: 
###################################################
cor.test(frogs,tadpoles)


###################################################
### chunk number 35:  eval=FALSE
###################################################
## help.search("correlation")


