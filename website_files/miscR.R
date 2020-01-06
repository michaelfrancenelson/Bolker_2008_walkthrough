###################################################
### chunk number 1:  eval=FALSE
###################################################
## desktop=getwd()
## options(repos="http://cran.us.r-project.org")
## install.packages("plotrix",destdir=desktop,lib=desktop)
## library(plotrix,lib=desktop)
## install.packages("gplots",destdir=desktop,lib=desktop)
## install.packages("gtools",destdir=desktop,lib=desktop)
## install.packages("gdata",destdir=desktop,lib=desktop)
## library(gtools,lib=desktop)
## library(gdata,lib=desktop)
## library(gplots,lib=desktop)


###################################################
### chunk number 2: 
###################################################
timevec1 = c("11:00:00","11:25:30","15:30:20")
times1 = times(timevec1)


###################################################
### chunk number 3: 
###################################################
set.seed(1001)
mydata = data.frame(indiv=rep(1:3,c(3,4,5)),
  sex=factor(c(rep("F",7),rep("M",5))),
  day=c(1:3,1:4,1:5),dist=runif(12))


###################################################
### chunk number 4: 
###################################################
r1 = reshape(mydata,direction="wide",idvar="indiv",timevar="day",
        v.names="dist"); r1


###################################################
### chunk number 5: 
###################################################
table(r1$sex)


###################################################
### chunk number 6: 
###################################################
splitdata = split.data.frame(mydata,mydata$indiv)
firstlines = lapply(splitdata,function(x)x[1,])
recombined = do.call("rbind",firstlines)


