##
## Pawel Gajer miscellaneous routines
##

library(Hmisc) # for ldesc()
library(dirichletprocess)

intersect.list <- function(y.l)
{
    cn <- names(y.l[[1]])
    for ( i in 2:length(y.l) )
    {
        cn <- intersect(cn, names(y.l[[i]]))
    }
    cn
}

plot.dirichletprocess.1d.posterior <- function(obj,
                                              with.hist=TRUE,
                                              br=100,
                                              hist.col="red",
                                              ci.size = 0.05,
                                              n.xgrid.pts=100,
                                              n.quant.pts=100,
                                              xlab="", ylab="Density",
                                              pt.pch=20, pt.cex=1, with.rug=T,
                                              med.col="blue", bci.col="gray90", ...)
{
    x.grid <- pretty(obj$data, n = n.xgrid.pts)
    posteriorFit <- replicate(n.quant.pts, PosteriorFunction(obj)(x.grid))
    posteriorCI <- apply(posteriorFit, 1, quantile, probs = c(ci.size/2, 0.5, 1 - ci.size/2), na.rm = TRUE)

    y.mean <- posteriorCI[2, ] # mean
    y.l <- posteriorCI[1, ] # lower boundary of 95% CI
    y.u <- posteriorCI[3, ] # upper boundary of 95% CI

    if ( with.hist ){
        hist(obj$data, br=br, col=2, las=1, ylim=c(0, max(y.u)), freq=FALSE,
             xlab=xlab, ylab=ylab, ...)
        polygon(c(x.grid, rev(x.grid)), c(y.l, rev(y.u)),  border = NA, col=bci.col)
        lines(x.grid, y.mean, col=med.col, lwd=1)
    } else {
        plot(x.grid, y.mean, type = "n", ylim=range(c(y.l,y.u)), las=1,  xlab=xlab, ylab=ylab, ...)
        polygon(c(x.grid, rev(x.grid)), c(y.l, rev(y.u)),  border = NA, col=bci.col)
        lines(x.grid, y.mean, col=med.col, lwd=1)
    }
}


plot.logit <- function(x, y,
                      y.med, y.l, y.u,
                      ylim=c(0,1),
                      title="",
                      alpha=0.05,
                      xlab="x", ylab="pr( y=1 )",
                      pt.pch=20, pt.cex=1, with.rug=T,
                      med.col="blue", bci.col="gray90")
{
    idx <- !is.na(y.med) & !is.na(y.l) & !is.na(y.u)
    x <- x[idx]
    y <- y[idx]
    y.med <- y.med[idx]
    y.l <- y.l[idx]
    y.u <- y.u[idx]

    plot(x, y.med, type = "n", ylim=ylim, las=1,  xlab=xlab, ylab=ylab, main=title) # ylim=range(c(theta.l, theta.u))
    polygon(c(x, rev(x)), c(y.l, rev(y.u)),  border = NA, col=bci.col)
    ##points(x, y, pch=pt.pch, cex=pt.cex)
    lines(x, y.med, col=med.col, lwd=1)
    if ( with.rug )
    {
        rug(x[y==0])
        rug(x[y==1], col='red', lwd=1)
    }
}

fit.and.plot.logit <- function(x, y,
                              title="",
                              xlab="x", ylab="pr( y=1 )",
                              pt.pch=20, pt.cex=1,
                              med.col="blue", bci.col="gray90", ...)
{
    o <- order(x)
    x <- x[o]
    y <- y[o]

    m <- glm(y~x, family=binomial)

    n <- length(x)
    t.val <- qt(0.975, n - 2) # Calculate critical t-value

    pred <- predict(m, type = "response", se.fit=TRUE)
    y.lower <- pred$fit - t.val*pred$se.fit
    y.hat <- pred$fit
    y.upper <- pred$fit + t.val*pred$se.fit

    plot.logit(x, y, y.hat, y.lower, y.upper,
               title=title, ylab=ylab, xlab=xlab, bci.col=bci.col, med.col=med.col, ...)

    invisible(m)
}

## helper function for bp.plot()
plot.vertical.ticks <- function(x, y, dx=0.1, col="gray")
{
    for ( i in seq(y) )
    {
        segments(x-dx, y[i], x+dx, y[i], col=col)
    }
}

## A modification of Hmisc library bbplot() routine
## with number of samples on the 3rd axis and ticks for sample values
bp.plot <- function (..., name = TRUE, main = "Box-Percentile Plot", xlab = "",
    ylab = "", srtx = 0, plotopts = NULL, dx=0.05, v.tick.col="gray")
{
    all.x <- list(...)
    nam <- character(0)
    if (is.list(all.x[[1]])) {
        all.x <- all.x[[1]]
        if (is.logical(name) && name)
            name <- names(...)
    }
    n <- length(all.x)
    n.samples <- sapply(all.x, length)
    centers <- seq(from = 0, by = 1.2, length = n)
    ymax <- max(sapply(all.x, max, na.rm = TRUE))
    ymin <- min(sapply(all.x, min, na.rm = TRUE))
    xmax <- max(centers) + 0.5
    xmin <- -0.5
    pargs <- c(list(c(xmin, xmax), c(ymin, ymax), type = "n", las=1,
        main = main, xlab = "", ylab = ylab, xaxt = "n"), plotopts)
    do.call("plot", pargs)
    for (i in 1:n) {
        plot.values <- bpx(all.x[[i]], centers[i])
        lines(plot.values$x1, plot.values$y1)
        lines(plot.values$x2, plot.values$y2)
        lines(plot.values$q1.x, plot.values$q1.y)
        lines(plot.values$q3.x, plot.values$q3.y)
        lines(plot.values$med.x, plot.values$med.y)
        plot.vertical.ticks(centers[i], all.x[[i]], dx=dx, col=v.tick.col)
    }
    mgp.axis(3, centers, n.samples)
    abline(v=centers, col="gray90")
    if (is.logical(name)) {
        if (name)
            mgp.axis(1, centers, sapply(substitute(list(...)),
                deparse)[2:(n + 1)], srt = srtx, adj = if (srtx ==
                0)
                0.5
            else 1, axistitle = xlab)
    }
    else mgp.axis(1, centers, name, srt = srtx,
                  adj = if (srtx==0) 0.5 else 1,
                  axistitle = xlab)
    invisible(centers)
}

##  a stick plot
stick.plot <- function(x, y, init=TRUE, ...)
{
    if ( init )
    {
        plot(x, y, type='n', xlab="", ylab="", las=1, ...)
    }
    for ( i in seq(x) )
    {
        segments(x[i], 0, x[i], y[i], ...)
    }
}

plot.loess <- function( x, y,
                     span=0.75,
                     no.ci=NULL,
                     title="",
                     xlab="log10( Relative Abundance )",
                     ylab="",
                     pt.col="gray20",
                     med.col='blue',
                     pch=20,
                     ylim=NULL,
                     pt.cex=0.5,
                     show.h.line=TRUE,
                     h.col="gray50",
                     h.lty=1, ...)
{
    o <- order(x)
    x <- x[o]
    y <- y[o]

    if ( is.null(ylim) )
    {
        ylim <- range(y)
    }

    plot(x, y, ylim=ylim, xlab=xlab, ylab=ylab, main=title, las=1, type='n', ...)
    ##plot(x, y)
    plx <- predict(loess(y ~ x, span=span), se=T)

    if ( is.null(no.ci) )
    {
        z <- qt(0.975,plx$df)
        yu <- plx$fit + z*plx$se
        yl <- plx$fit - z*plx$se
        polygon( c(x, rev(x)), c(yl, rev(yu)), col='gray90', border=NA )
        if ( show.h.line )
        {
            abline(h=plx$fit[1], col=h.col, lty=h.lty)
        }
    }

    lines( x, plx$fit, col=med.col )
    points( x, y, col=pt.col, pch=pch, cex=pt.cex )
    ## if ( is.null(no.ci) )
    ## {
    ##     ##abline(h=mean(y), col=h.col, lty=h.lty)

    ## }
    ## nStr <- paste("(n=",length(y),")", sep="")
    ## legend("bottomright", legend=c(nStr), inset=0.05)
    ## lines(x, plx$fit - qt(0.975,plx$df)*plx$se, lty=2)
    ## lines(x, plx$fit + qt(0.975,plx$df)*plx$se, lty=2)
}

# compute the mode of a density function 'd'
# or any object with x and y components
denMode <- function(d) d$x[which.max(d$y)]

## mode of a vector of values
mode <- function(x)
{
    d <- density(x)
    denMode(d)
}


loess.gEff <- function(ph, ref.ph, idx1, span=0.75)
{
    ph.nReads <- ct2[,ph]
    ref.ph.nReads <- ct2[,ref.ph]
    y <- mt$delStatus.n
    idx <- idx1 & ph.nReads > 0 & ref.ph.nReads > 0
    ph.nReads <- log10(ph.nReads[idx])
    ref.ph.nReads <- log10(ref.ph.nReads[idx])
    y <- y[idx]
    x <- ph.nReads - ref.ph.nReads

    d <- density(x)
    mode.log.rat <- denMode(d)

    plx <- predict(loess(y ~ x, span=span), se=T)
    z <- qt(0.975,plx$df)
    yu <- plx$fit + z*plx$se
    yl <- plx$fit - z*plx$se

    y <- plx$fit

    i.min <- which.min(y)
    x.min <- x[i.min]
    y.min <- y[i.min]

    i.max <- which.max(y)
    x.max <- x[i.max]
    y.max <- y[i.max]

    pval.max <- pnorm(y[i.min], y[i.max], plx$se[i.max])[[1]]
    pval.min <- 1 - pnorm(y[i.max], y[i.min], plx$se[i.min])[[1]]

    if ( pval.max > pval.min )
    {
        pval <- pval.max
    } else {
        pval <- pval.min
    }

    n <- length(y)
    eff <- NA

    ## percentage of points over which first derivative of y.fit is positive
    ## (y.fit is increasing)
    dy.fit <- diff(y)
    perc.incr <- 100*sum(dy.fit>0)/n
    incr.thld <- 60

    if ( (perc.incr > incr.thld) || (i.min!=1 && i.min!=n) ) # risk is increasing or U-shaped
    {
        eff <- y.max - y.min
    } else {
        eff <- -(y.max - y.min)
    }

    list(pval=pval, pval.max=pval.max, pval.min=pval.min,
         eff=eff, x.min=x.min, y.min=y.min, x.max=x.max, y.max=y.max,
         mode.log.rat=mode.log.rat)
}


n.stats2 <- function(ph, ref.ph, idx1)
{
    x1 <- ct2[,ph]
    x2 <- ct2[,ref.ph]
    ##n <- rowSums(ct2)
    y <- mt$delStatus.n
    s <- mt$subjID
    idx <- idx1 & x1 > 0 & x2 > 0
    x1 <- x1[idx]
    x2 <- x2[idx]
    y <- y[idx]
    s <- s[idx]
    ## x <- log10(x1) - log10(x2)
    ## q <- quantile(x, probs=c(0.005, 0.995)) # removing potential outliers
    ## idx <- x > q[1] & x/n < q[2]
    ## x <- x[idx]
    ## y <- y[idx]
    ## s <- s[idx]
    nCaseSubjs <- length(unique(s[y==1]))
    nCtrkSubjs <- length(unique(s[y==0]))

    list(nCaseSubjs=nCaseSubjs,
         nCtrkSubjs=nCtrkSubjs,
         n=length(y), n.sPTB=sum(y),
         n.TERM=length(y)-sum(y))
}


##
## R interface to varInfo() C routine in
## /Users/pgajer/.Rlocal/lib/Pmisc.c
##

varInfo <- function(cl1, cl2){

  if ( length(cl1) != length(cl2) )
    stop("the arguments of varInfo must be of the same length")

  .C("varInfo",
     as.integer(cl1),
     as.integer(cl2),
     as.integer(length(cl1)),
     as.double(0.0))[[4]]
}

normVarInfo <- function(cl1, cl2){

  if ( length(cl1) != length(cl2) )
    stop("the arguments of varInfo must be of the same length")

  .C("normVarInfo",
     as.integer(cl1),
     as.integer(cl2),
     as.integer(length(cl1)),
     as.double(0.0))[[4]]
}

# cross validation of clustering parameter (number of clusters)
# x        - data table,
# nClrs    - number of clusters
# method   - hclust method
# f        - fraction of input data (rows)
# nSamples - number of random samples
#
crossValidate <- function(x, nClrs, f=0.9, nSamples=100, method="average"){

  n <- nrow(x)
  sampleSize = as.integer(f * n)
  v <- numeric(nSamples)
  timeStr <- format(Sys.time(), "%X")
  timeStr <- gsub(":",".",timeStr)

  for ( j in 1:nSamples ){

    s1 <- sample(1:n,sampleSize)
    s2 <- sample(1:n,sampleSize)

    x1 <- x[s1,]
    x2 <- x[s2,]

    hc1 <- hclust(dist(x1),method=method)
    memb1 <- cutree(hc1,k=nClrs)
    m1 <- numeric(n)
    for ( i in 1:sampleSize )
      m1[s1[i]] <- memb1[i]

    hc2 <- hclust(dist(x2),method=method)
    memb2 <- cutree(hc2,k=nClrs)
    m2 <- numeric(n)
    for ( i in 1:sampleSize )
      m2[s2[i]] <- memb2[i]

    s <- intersect(s1,s2)

    cl1 <- m1[s]
    cl2 <- m2[s]

    v[j] <- varInfo(cl1,cl2)
  }

  v
}


prob <- function(cl)
{
  nc <- 0
  for ( i in 1:length(cl) )
    if ( cl[i] > nc )
      nc <- cl[i]

  nc <- nc+1
  pp <- rep(0,nc)
  h <- 0

  .C("prob2",
     as.integer(cl),
     as.integer(length(cl)),
     as.double(h),
     p=as.double(pp))$p
}


entropy <- function(p)
{
  h <- 0

  .C("entropy2",
     as.double(p),
     as.integer(length(p)),
     out=as.double(h))$out
}

relEntropy <- function(p,q){

  if ( length(p) != length(q) )
    stop("Error in relEntropy() lengths of the arguments have to be the same.")

  h <- 0

  .C("relEntropy2",
     as.double(p),
     as.double(q),
     as.integer(length(p)),
     out=as.double(h))$out
}


JSdist <- function(p,q)
{
    p <- p/sum(p)
    q <- q/sum(q)
    a <- (p+q)/2
    (relEntropy(p,a) + relEntropy(q,a))/2
}

## compute Jensen-Shannon distance between rows of matrix x assuming that each
## row is a vector of counts for given two rows, rarefy the one with larger sum
## the sum of the other vector and then compute Jensen-Shannon distance on the
## proportion vectors of these count vectors
ssJSdist <- function(x)
{
    x <- as.matrix(x)
    nr <- nrow(x)
    nr1 <- nr-1
    nc <- ncol(x)

    d <- matrix(0, nrow=nr, ncol=nr)
    for ( i in 1:nr1 )
    {
        for ( j in (i+1):nr )
        {
            if ( (n1 <- sum(x[i,])) > (n2 <- sum(x[j,])) )
            {
                r1 <- rrarefy(x[i,], n2)
                d[i,j] <- JSdist(r1, x[j,])
                d[j,i] <- d[i,j]
            } else {
                r1 <- rrarefy(x[j,], n1)
                d[i,j] <- JSdist(r1, x[i,])
                d[j,i] <- d[i,j]
            }
        }
    }
    d
}


symmRelEntropy <- function(p,q){

  (relEntropy(p,q) + relEntropy(q,p))/2
}


matr.shrink <- function(x)
{
    nr <- nrow(x)
    nc <- ncol(x)
    ##xout <- .C("shrink_mat", as.double(x), nr, nc, sx = double(nr*nc))$sx
    xout <- .C("matr_shrink", x=as.double(x), nr, nc)$x
    matrix(xout,nrow=nr,ncol=nc,byrow=F)
}

vect.shrink <- function(x)
{
    n <- length(x)
    .C("vect_shrink", as.double(x), n, sx = double(n))$sx
}

freqs.shrink2 <- function(x)
{
    n <- length(x)
    .C("freqs_shrink", xout=double(x), n)$xout
}


ldesc <- function(data, nomargins=T, width=8.5, height=11) {
    options(xdvicmd='xdvi')
    d <- describe(data, desc=deparse(substitute(data)))
    dvi(latex(d, file='/tmp/z.tex'),
        nomargins=nomargins, width=width, height=height)
}


cexFn <- function(x){ -0.75*(x-2) + 2}

volcano.plot <- function(signal.mat, q.thld=0.05, e.thld=0.1,
                        legend.title="Median Relative Abundance",
                        xlab=expression(paste("max ", Delta, " pr( sPTB )")), ylab="-log10( q-values )",
                        xlim=c(-1.1,1.1), ylim=c(0,16), lfactor=1, rfactor=1.4,
                        yjust=-0.15, xjust=-0.15, adj=c(-0.07, -0.6), text.cex=0.8)
{
    op <- par(mar=c(4, 4, 4, 0.5), mgp=c(2.25,0.6,0),tcl = -0.3) #cex.axis=0.7)
    x <- signal.mat[,1]
    idx <- signal.mat[,3] < 1e-14
    if ( sum(idx) )
    {
        signal.mat[idx,3] <- runif(sum(idx), 1e-16, 1e-14)
    }
    y <- -log10(signal.mat[,3])
    if ( any(is.na(xlim)) )
    {
        r <- range(x)
        xlim <- c(lfactor*r[1], rfactor*r[2])
    }
    if ( any(is.na(ylim)) )
    {
        ylim <- range(y[is.finite(y)])
    }
    plot(x,y, las=1, xlab=xlab, ylab=ylab, pch=19, xlim=xlim, ylim=ylim, type='n')
    abline(h=-log10(q.thld), col='gray80')
    abline(v=0, col='gray80')
    for ( i in seq(nrow(signal.mat)) )
    {
        if ( signal.mat[i,3] < q.thld & abs(signal.mat[i,1]) > e.thld )
        {
            if ( x[i] > 0 ){
                col <- "red"
            } else {
                col <- "green"
            }
            text(x[i], y[i], labels=gsub("_"," ",rownames(signal.mat)[i]), adj=adj, cex=text.cex, font=3)
            points(x[i], y[i], pch=19, col=col, cex=cexFn(signal.mat[i,ncol(signal.mat)]))
        } else {
            points(x[i], y[i], pch=20)
        }
    }
    legend(par('usr')[1], par('usr')[4], yjust=yjust, xjust=xjust, xpd=NA,
           legend=c("0.001","0.01","0.1"), pch=19, pt.cex=cexFn(c(3,2,1)),ncol=3,
           title=legend.title)
    par(op)
}


poly2.volcano.plot <- function(signal.mat, q.thld=0.05,
                              xlab="log OR of pr( sPTB )", ylab="-log10( q-values )",
                              xlim1=NA, xlim2=NA, ylim=c(0,16), lfactor=1, rfactor=2.4,
                              yjust=-0.15, xjust=-0.15, adj=c(-0.07, -0.6), text.cex=0.8)
{
    op <- par(mfrow=c(1,2), mar=c(4, 4, 4, 0.5), mgp=c(2.25,0.6,0),tcl = -0.3) #cex.axis=0.7)
    ## linear term
    x <- signal.mat[,1]
    idx <- signal.mat[,3] < 1e-14
    if ( sum(idx) )
    {
        signal.mat[idx,3] <- runif(sum(idx), 1e-16, 1e-14)
    }
    y <- -log10(signal.mat[,3])
    if ( any(is.na(xlim1)) )
    {
        r <- range(x)
        xlim1 <- c(lfactor*r[1], rfactor*r[2])
    }
    if ( any(is.na(ylim)) )
    {
        ylim <- range(y[is.finite(y)])
    }
    plot(x,y, las=1, xlab=xlab, ylab=ylab, pch=19, xlim=xlim1, ylim=ylim, type='n')
    abline(h=-log10(q.thld), col='gray80')
    abline(v=0, col='gray80')
    for ( i in seq(nrow(signal.mat)) )
    {
        if ( signal.mat[i,3] < q.thld )
        {
            if ( x[i] > 0 ){
                col <- "red"
            } else {
                col <- "green"
            }
            text(x[i], y[i], labels=gsub("_"," ",rownames(signal.mat)[i]), adj=adj, cex=text.cex, font=3)
            points(x[i], y[i], pch=19, col=col, cex=cexFn(signal.mat[i,ncol(signal.mat)]))
        } else {
            points(x[i], y[i], pch=20)
        }
    }
    legend(par('usr')[1], par('usr')[4], yjust=yjust, xjust=xjust, xpd=NA,
           legend=c("0.001","0.01","0.1"), pch=19, pt.cex=cexFn(c(3,2,1)),ncol=3,
           title="Median Relative Abundance")
    ## quadratic term
    x <- signal.mat[,1+3]
    idx <- signal.mat[,3+3] < 1e-14
    if ( sum(idx) )
    {
        signal.mat[idx,3+3] <- runif(sum(idx), 1e-16, 1e-14)
    }
    y <- -log10(signal.mat[,3+3])
    if ( any(is.na(xlim2)) )
    {
        r <- range(x)
        xlim2 <- c(lfactor*r[1], rfactor*r[2])
    }
    if ( any(is.na(ylim)) )
    {
        ylim <- range(y[is.finite(y)])
    }
    plot(x,y, las=1, xlab=xlab, ylab=ylab, pch=19, xlim=xlim2, ylim=ylim, type='n')
    abline(h=-log10(q.thld), col='gray80')
    abline(v=0, col='gray80')
    for ( i in seq(nrow(signal.mat)) )
    {
        if ( signal.mat[i,3+3] < q.thld )
        {
            if ( x[i] > 0 ){
                col <- "red"
            } else {
                col <- "green"
            }
            text(x[i], y[i], labels=gsub("_"," ",rownames(signal.mat)[i]), adj=adj, cex=text.cex, font=3)
            points(x[i], y[i], pch=19, col=col, cex=cexFn(signal.mat[i,ncol(signal.mat)]))
        } else {
            points(x[i], y[i], pch=20)
        }
    }
    legend(par('usr')[1], par('usr')[4], yjust=yjust, xjust=xjust, xpd=NA,
           legend=c("0.001","0.01","0.1"), pch=19, pt.cex=cexFn(c(3,2,1)),ncol=3,
           title="Median Relative Abundance")
    par(op)
}
