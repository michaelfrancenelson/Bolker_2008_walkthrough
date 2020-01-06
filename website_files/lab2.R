###################################################
### chunk number 1: 
###################################################
data = read.table("seedpred.dat",header=TRUE)


###################################################
### chunk number 2: 
###################################################
data$available=data$taken+data$remaining


###################################################
### chunk number 3:  eval=FALSE
###################################################
## count.fields("myfile.dat",sep=",")


###################################################
### chunk number 4:  eval=FALSE
###################################################
## cf = count.fields("myfile.dat",sep=",")
## which(cf!=cf[1])


###################################################
### chunk number 5:  eval=FALSE
###################################################
## mydata <- read.csv("myfile.dat",fill=FALSE)


###################################################
### chunk number 6: 
###################################################
sapply(data,class)


###################################################
### chunk number 7:  eval=FALSE
###################################################
## mydata <- read.table("mydata.dat",na.strings="*")


###################################################
### chunk number 8: 
###################################################
loc = factor(rep(LETTERS[1:3],2))
day = factor(rep(1:2,each=3))
val = round(runif(6),3)
d = data.frame(loc,day,val)


###################################################
### chunk number 9: 
###################################################
d2 = reshape(d,direction="wide",idvar="loc",timevar="day"); d2


###################################################
### chunk number 10: 
###################################################
reshape(d2,direction="long",varying=c("val.1","val.2"),
        timevar="day",idvar="loc")


###################################################
### chunk number 11: 
###################################################
data2 = read.table("seedpred.dat",header=TRUE,as.is="Species")
data2 = read.table("seedpred.dat",header=TRUE,as.is=1)
sapply(data2,class)


###################################################
### chunk number 12: 
###################################################
data2 = read.table("seedpred.dat",header=TRUE,colClasses=c("character",rep("numeric",4)))


###################################################
### chunk number 13: 
###################################################
data2 = read.table("seedpred.dat",header=TRUE)
sapply(data2,class)
data2$Species = as.character(data2$Species)
sapply(data2,class)


###################################################
### chunk number 14: 
###################################################
data2 = read.table("seedpred.dat",header=TRUE,colClasses=c(rep("factor",2),rep("numeric",3)))
sapply(data2,class)


###################################################
### chunk number 15: 
###################################################
f = factor(1:10); levels(f)


###################################################
### chunk number 16: 
###################################################
f = factor(as.character(1:10)); levels(f)


###################################################
### chunk number 17: 
###################################################
f = factor(as.character(1:10),levels=1:10)
x = c("north","middle","south")
f = factor(x,levels=c("far_north","north","middle","south"))


###################################################
### chunk number 18: 
###################################################
f = factor(c(3,3,5,6,7,8,10),levels=3:10)


###################################################
### chunk number 19: 
###################################################
f = factor(c("a","b","c","d"))
f2 = f[1:2]
levels(f2)
f2 = factor(as.character(f2))
levels(f2)


###################################################
### chunk number 20: 
###################################################
as.Date(c("1jan1960", "2jan1960", "31mar1960", "30jul1960"),
        format="%d%b%Y")
as.Date(c("02/27/92", "02/27/92", "01/14/92", "02/28/92", "02/01/92"),
        format="%m/%d/%y")


###################################################
### chunk number 21: 
###################################################
year = c(2004,2004,2004,2005)
month = c(10,11,12,1)
day = c(20,18,28,17)
datestr = paste(year,month,day,sep="/")
date = as.Date(datestr)
date


###################################################
### chunk number 22:  eval=FALSE
###################################################
## data = read.table("datafile",sep="|",quote="") 


###################################################
### chunk number 23: 
###################################################
attach(data)


###################################################
### chunk number 24:  eval=FALSE
###################################################
## data(dataset)


###################################################
### chunk number 25:  eval=FALSE
###################################################
## install.packages("plotrix")


###################################################
### chunk number 26: 
###################################################
library(plotrix)


###################################################
### chunk number 27:  eval=FALSE
###################################################
## sizeplot(available,taken,xlab="Available",ylab="Taken")


###################################################
### chunk number 28: 
###################################################
t1 = table(available,taken)


###################################################
### chunk number 29:  eval=FALSE
###################################################
## r = row(t1)
## c = col(t1)-1
## text(r[t1>0],c[t1>0],t1[t1>0])


###################################################
### chunk number 30:  eval=FALSE
###################################################
## barplot(t(log10(t1+1)),beside=TRUE,legend=TRUE,xlab="Available",
##         ylab="log10(1+# observations)")
## op = par(xpd=TRUE)
## text(34.5,3.05,"Number taken")
## par(op)


###################################################
### chunk number 31: 
###################################################
x = 1:10
col_vec = rep(1:2,length=10)
pch_vec = rep(1:2,each=5)
plot(x,col=col_vec,pch=pch_vec)


###################################################
### chunk number 32:  eval=FALSE
###################################################
## v = as.numeric(log10(1+t1))
## plot(sort(v))
## r = row(t1)
## c = col(t1)
## plot(v,col=r,pch=c)
## plot(v[order(v)],col=r[order(v)],pch=c[order(v)])


###################################################
### chunk number 33:  eval=FALSE
###################################################
## library(lattice)


###################################################
### chunk number 34:  eval=FALSE
###################################################
## library(lattice)
## barchart(log10(1+table(available,taken)),stack=FALSE,
##          auto.key=TRUE)


###################################################
### chunk number 35:  eval=FALSE
###################################################
## barchart(log10(1+table(available,Species,taken)),stack=FALSE,
##          auto.key=TRUE)


###################################################
### chunk number 36: 
###################################################
frac_taken = taken/available


###################################################
### chunk number 37: 
###################################################
mean_frac_by_avail = tapply(frac_taken,available,mean)


###################################################
### chunk number 38: 
###################################################
n_by_avail = table(available)
se_by_avail = tapply(frac_taken,available,sd)/sqrt(n_by_avail)


###################################################
### chunk number 39:  eval=FALSE
###################################################
## library(gplots)
## lower_lim = mean_frac_by_avail-se_by_avail
## upper_lim = mean_frac_by_avail+se_by_avail
## b = barplot2(mean_frac_by_avail,plot.ci=TRUE,
##   ci.l=lower_lim,ci.u=upper_lim,
##   xlab="Number available",
##   ylab="Mean number taken")


###################################################
### chunk number 40:  eval=FALSE
###################################################
## histogram(~frac_taken|Species,xlab="Fraction taken")


###################################################
### chunk number 41: 
###################################################
splitdat = split(frac_taken,Species)


###################################################
### chunk number 42:  eval=FALSE
###################################################
## op=par(mfrow=c(3,3),mar=c(2,2,1,1))


###################################################
### chunk number 43:  eval=FALSE
###################################################
## h=lapply(splitdat,hist,xlab="",ylab="",main="",col="gray")


###################################################
### chunk number 44: 
###################################################
detach(data)
rm(list=ls())
data = read.table("ewcitmeas.dat",header=TRUE,na.strings="*")
attach(data)


###################################################
### chunk number 45: 
###################################################
date = as.Date(paste(year+1900,mon,day,sep="/"))
city_names = colnames(data)[4:10]


###################################################
### chunk number 46: 
###################################################
data = cbind(data,date)
data_long=reshape(data,direction="long",
  varying=list(city_names),v.name="incidence",
  drop=c("day","mon","year"),times=factor(city_names),
  timevar="city")


###################################################
### chunk number 47:  eval=FALSE
###################################################
## matplot(date,data[,4:10],type="l",col=1:7,lty=1:7,axes=FALSE,
##         ylab="Weekly incidence",xlab="Date")
## axis(side=2)
## axis.Date(side=1,x=date)
## vacc.date = as.Date("1968/1/1")
## biennial =  seq.Date(as.Date("1948/9/1"),as.Date("1986/9/1"),by="2 years")
## abline(v=biennial,col="gray",lty=2)
## abline(v=vacc.date,lty=2,lwd=2)
## legend(x=1970,y=5000,city_names,col=1:7,lty=1:7,lwd=2,
##        bg="white")
## box()


###################################################
### chunk number 48:  eval=FALSE
###################################################
## xyplot(incidence~date,groups=city,data=data_long,type="l",auto.key=TRUE)


###################################################
### chunk number 49: 
###################################################
allvals = na.omit(c(as.matrix(data[,4:10])))
logvals = log10(1+allvals)


###################################################
### chunk number 50:  eval=FALSE
###################################################
## hist(logvals,col="gray",
##      main="",xlab="Log weekly incidence",ylab="Density",freq=FALSE,
##      ylim=c(0,0.6))


###################################################
### chunk number 51:  eval=FALSE
###################################################
## lines(density(logvals),lwd=2)
## lines(density(logvals,adjust=0.5),lwd=2,lty=2)


###################################################
### chunk number 52:  eval=FALSE
###################################################
## curve(dnorm(x,mean=mean(logvals),sd=sd(logvals)),lty=3,lwd=2,add=TRUE)


###################################################
### chunk number 53:  eval=FALSE
###################################################
## legend(x=2.1,y=0.62,
##        legend=c("density, default",
##          "density, adjust=0.5","normal"),
##        lwd=2,lty=c(1,2,3))


###################################################
### chunk number 54:  eval=FALSE
###################################################
## logscaledat = as.data.frame(log10(scale(1+data[,4:10],
##                            center=FALSE,
##                            scale=colMeans(1+data[,4:10],na.rm=TRUE))))


###################################################
### chunk number 55: 
###################################################
city_means <- tapply(1+data_long$incidence,data_long$city,mean,na.rm=TRUE)


###################################################
### chunk number 56: 
###################################################
scdat <- (1+data_long$incidence)/city_means[data_long$city]


###################################################
### chunk number 57:  eval=FALSE
###################################################
## plot(density(na.omit(logscaledat[,1])),
##      type="n",main="",xlab="Log scaled incidence")


###################################################
### chunk number 58: 
###################################################
tmpfun = function(x,i) {
  lines(density(na.omit(x)),lwd=2,col=i,lty=i)
}


###################################################
### chunk number 59:  eval=FALSE
###################################################
## m = mapply(tmpfun,logscaledat,1:7)


###################################################
### chunk number 60: 
###################################################
legend(-2.6,0.65,city_names,lwd=2,col=1:7,lty=1:7)


###################################################
### chunk number 61:  eval=FALSE
###################################################
## densityplot(~log10(scdat),groups=data_long$city,plot.points=FALSE,auto.key=TRUE,
##             lty=1:7)


###################################################
### chunk number 62: 
###################################################
city_abbr = substr(city_names,1,3)


###################################################
### chunk number 63:  eval=FALSE
###################################################
## boxplot(log10(1+incidence)~city,data=data_long,ylab="Log(incidence+1)",
##         names=city_abbr)


###################################################
### chunk number 64:  eval=FALSE
###################################################
## bwplot(log10(1+incidence)~city,data=data_long,
##              panel=panel.violin,
##              horizontal=FALSE,
##              scales=list(abbreviate=TRUE))


###################################################
### chunk number 65: 
###################################################
data(quakes)


###################################################
### chunk number 66:  eval=FALSE
###################################################
## pairs(quakes,pch=".")
## splom(quakes,pch=".")


###################################################
### chunk number 67:  eval=FALSE
###################################################
## coplot(lat ~ long | depth, data = quakes)


###################################################
### chunk number 68: 
###################################################
tmpdat = quakes[quakes$long>175,]


###################################################
### chunk number 69:  eval=FALSE
###################################################
## plot(tmpdat$long,tmpdat$depth,xlab="Longitude",ylab="Depth",
##      col="darkgray",pch=".")


###################################################
### chunk number 70:  eval=FALSE
###################################################
## lines(lowess(tmpdat$long,tmpdat$depth),lwd=2)


###################################################
### chunk number 71: 
###################################################
lines(smooth.spline(tmpdat$long,tmpdat$depth),lwd=2,lty=2)
lines(smooth.spline(tmpdat$long,tmpdat$depth,df=4),lwd=2,lty=3)


###################################################
### chunk number 72: 
###################################################
abline(lm(depth ~long,data=tmpdat),lwd=2,col="gray")


###################################################
### chunk number 73: 
###################################################
quad.lm = lm(depth ~long+I(long^2),data=tmpdat)


###################################################
### chunk number 74: 
###################################################
lvec = seq(176,188,length=100)
quadvals = predict(quad.lm,newdata=data.frame(long=lvec))


###################################################
### chunk number 75: 
###################################################
lines(lvec,quadvals,lwd=2,lty=2,col="gray")
legend(183.2,690,c("lowess","spline (default)","spline (df=4)","regression",
                   "quad. regression"),
       lwd=2,lty=c(1,2,3,1,2),col=c(rep("black",3),rep("gray",2)))


