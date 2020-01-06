###################################################
### chunk number 1: 
###################################################
loc = factor(rep(LETTERS[1:3],2))
day = factor(rep(1:2,each=3))
set.seed(1001)
val = round(runif(6),3)
d = data.frame(loc,day,val); d


###################################################
### chunk number 2: 
###################################################
unstack(d,val~day)


###################################################
### chunk number 3: 
###################################################
unstack(d,val~loc)


###################################################
### chunk number 4: 
###################################################
f=factor(c(3,3,5,6,7,8,10))
op=par(mfrow=c(1,2))
plot(f)
f=factor(c(3,3,5,6,7,8,10),levels=3:10)
plot(f)
par(op)


###################################################
### chunk number 5: 
###################################################
data = read.table("seedpred.dat",header=TRUE)
data$available=data$remaining+data$taken
t1 = table(data$available,data$taken)
v = as.numeric(log10(1+t1))
r = row(t1)
c = col(t1)


###################################################
### chunk number 6: 
###################################################
v_sorted = v[order(v)]
r_sorted = r[order(v)]
c_sorted = c[order(v)]


###################################################
### chunk number 7: 
###################################################
op=par(mfrow=c(2,2),mgp=c(2,1,0),mar=c(4.2,3,1,1))
plot(sort(v))
plot(v,col=r,pch=c)
plot(v_sorted,col=r_sorted,pch=c_sorted)
legend(0,2.8,pch=1,col=1:5,legend=1:5)
legend(6,2.8,pch=1:6,col=1,legend=0:5)
text(0,3,"available",adj=0)
text(8,3,"taken",adj=0)
par(op)


###################################################
### chunk number 8: 
###################################################
data = read.table("seedpred.dat",header=TRUE)
data2 = data
data2$available=data2$remaining+data2$taken
data2 = data2[data2$available==5,]
t1 = table(data2$taken,data2$Species)


###################################################
### chunk number 9: 
###################################################
op=par(mfrow=c(2,1),mgp=c(2.5,1,0),mar=c(4.1,3.5,1.1,1.1))
logt1=log10(1+t1)
barplot(logt1,beside=TRUE,ylab="log10(1+taken)")
library(gplots)
barplot2(t1+1,beside=TRUE,log="y",ylab="taken+1")
par(op)


###################################################
### chunk number 10: 
###################################################
data = read.table("ewcitmeas.dat",header=TRUE,na.strings="*")


###################################################
### chunk number 11: 
###################################################
incidence = data[,4:10]
imin = apply(incidence,2,min,na.rm=TRUE)
imax = apply(incidence,2,max,na.rm=TRUE)
irange = imax-imin


###################################################
### chunk number 12: 
###################################################
iranges = apply(incidence,2,range,na.rm=TRUE); iranges
irange = iranges[2,]-iranges[1,]


###################################################
### chunk number 13: 
###################################################
rangediff = function(x) {
  diff(range(x,na.rm=TRUE))
}
irange = apply(incidence,2,rangediff)


###################################################
### chunk number 14: 
###################################################
scaled_incidence = scale(incidence,center=imin,scale=irange)


###################################################
### chunk number 15: 
###################################################
summary(scaled_incidence)
apply(scaled_incidence,2,range,na.rm=TRUE)


###################################################
### chunk number 16: 
###################################################
imean = colMeans(incidence,na.rm=TRUE)
scaled_incidence = sweep(incidence,2,imean,"-")


###################################################
### chunk number 17: 
###################################################
c1 = colMeans(scaled_incidence,na.rm=TRUE); c1


###################################################
### chunk number 18: 
###################################################
all(abs(c1)<1e-11)


###################################################
### chunk number 19: 
###################################################
date = as.Date(paste(data$year+1900,data$mon,data$day,sep="/"))
city_names = colnames(data)[4:10]
data = cbind(data,date)
data_long=reshape(data,direction="long",
  varying=list(city_names),v.name="incidence",
  drop=c("day","mon","year"),times=factor(city_names),
  timevar="city")


###################################################
### chunk number 20: 
###################################################
city_max = tapply(data_long$incidence,data_long$city,max,na.rm=TRUE)
city_min = tapply(data_long$incidence,data_long$city,min,na.rm=TRUE)
range1 = city_max-city_min


###################################################
### chunk number 21: 
###################################################
scdat1 = data_long$incidence-city_min[data_long$city]
scdat = scdat1/range1[data_long$city]


###################################################
### chunk number 22: 
###################################################
tapply(scdat,data_long$city,range,na.rm=TRUE)


