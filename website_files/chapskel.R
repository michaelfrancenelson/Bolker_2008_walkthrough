rm(list=ls())
badfonts = FALSE
.detach.stuff = function() {
  s1 = grep("package|Autoloads",search())
  nattach = length(search())
  xattach = search()[c(-1,-s1)]
  for (i in xattach) 
    eval(substitute(detach(i),list(i=i)))
}
.detach.stuff()
## ps.options(family="NimbusSan")
ps.options(colormodel="cmyk")
options(width=65)
## palette(gray((0:8)/8))
library(lattice)
trellis.par.set(canonical.theme(color=FALSE))

## corner label
corner.label2 <- function(label,x,y,inset=0,pos,cex=par("cex"),
                          bg=NULL,border=NA,...) {
  w = strwidth(label,cex=cex)
  h = strheight(label,cex=cex)
  auto <- if (is.character(x)) 
    match.arg(x, c("bottomright", "bottom", "bottomleft", 
                   "left", "topleft", "top", "topright", "right", "center"))
  else NA
  if (is.na(auto)) {
    xy <- xy.coords(x, y)
    x <- xy$x
    y <- xy$y
    nx <- length(x)
    if (nx < 1 || nx > 2) 
      stop("invalid coordinate lengths")
  }
  else nx <- 0
  if (is.na(auto)) {
    left <- x - xjust * w
    top <- y + (1 - yjust) * h
  }
  else {
    usr <- par("usr")
    inset <- rep(inset, length.out = 2)
    insetx <- inset[1] * (usr[2] - usr[1])
    left <- switch(auto, bottomright = , topright = , 
                   right = usr[2] - w - insetx, bottomleft = , left = , 
                   topleft = usr[1] + insetx, bottom = , top = , 
                   center = (usr[1] + usr[2] - w)/2)
    insety <- inset[2] * (usr[4] - usr[3])
    top <- switch(auto, bottomright = , bottom = , bottomleft = usr[3] + 
                  h + insety, topleft = , top = , topright = usr[4] - 
                  insety, left = , right = , center = (usr[3] + 
                                                       usr[4] + h)/2)
  }
  if (par("xlog")) { left <- 10^left }
  if (par("ylog")) { top <- 10^top }
  if (!is.null(bg)) {
    if (bg==TRUE) bg <- "white"
    ## doesn't really work for log scales yet
    if (!is.na(border) && border==TRUE) border <- par("fg")
    rect(left,top-h/2,left+w,top+h/2,col=bg,border=border)
  }                    
  text(left, top, label, adj=0, ...)
}

## rectangle behind text
textrec = function(x,y,labels,cex,border=NA,bg=NULL,expand=1.1) {
  if (!is.null(bg)) {
    if (bg==TRUE) bg <- "white"
  }
  if (!is.na(border) && border==TRUE) border <- par("fg")
  w = strwidth(labels,cex=cex)
  h = strheight(labels,cex=cex)
  rect(x-w/2*expand,y-h/2*expand,x+w/2*expand,y+h/2*expand,
       col=bg,border=border)
}

