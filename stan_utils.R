##
## Stan utilities
##

## in /Users/pgajer/projects/Statistics/MCMC_Models
## cp spmrf_binomial_horseshoe_order2.stan spmrf_bernoulli_horseshoe_order2.stan
## manually edited spmrf_bernoulli_horseshoe_order2.stan to be bernoulli model
## spmrf.bernoulli.o2.model <- stan_model(file="/Users/pgajer/projects/Statistics/MCMC_Models/spmrf_bernoulli_horseshoe_order2.stan")
## save(spmrf.bernoulli.o2.model, file="/Users/pgajer/projects/Statistics/MCMC_Models/spmrf_bernoulli_horseshoe_order2.rda")
load("/Users/pgajer/projects/Statistics/MCMC_Models/spmrf_bernoulli_horseshoe_order2.rda")

## spmrf.bernoulli.o3.model <- stan_model(file="/Users/pgajer/projects/Statistics/MCMC_Models/spmrf_bernoulli_horseshoe_order3.stan")
## save(spmrf.bernoulli.o3.model, file="/Users/pgajer/projects/Statistics/MCMC_Models/spmrf_bernoulli_horseshoe_order3.rda")
load("/Users/pgajer/projects/Statistics/MCMC_Models/spmrf_bernoulli_horseshoe_order3.rda")

load("/Users/pgajer/projects/M_and_M/data/mad_log10RA_vs_medianlog10nReads_plot_exp_model.rda") # m.dat, m.fit, nexp.m, nexpFn

## $ scp /Users/pgajer/organizer/programming/R/libs/stan_utils.R pawel@OPEE-iMac.local:organizer/programming/R/libs
## $ scp /Users/pgajer/organizer/programming/R/libs/stan_utils.R cadbane.igs.umaryland.edu:/home/pgajer/projects/Stat_Utilities

## b2p <- stan_model(file="/Users/pgajer/projects/Statistics/MCMC_Models/bin_2prop.stan")
## save(b2p, file="/Users/pgajer/projects/Statistics/MCMC_Models/bin_2prop.rda")
load("/Users/pgajer/projects/Statistics/MCMC_Models/bin_2prop.rda") # b2p


## extracting quantiles
get.qq <- function(y)
{
    q <- quantile(y, probs = c(0.25,0.5, 0.75))
    q1 <- ifelse(y<=q[1], 1, 0)
    q2 <- ifelse(y<=q[2] & y>q[1], 1, 0)
    q3 <- ifelse(y<=q[3] & y>q[2], 1, 0)
    q4 <- ifelse(y>q[3], 1, 0)
    qq <- y
    qq[q1==1] <- "q1"
    qq[q2==1] <- "q2"
    qq[q3==1] <- "q3"
    qq[q4==1] <- "q4"

    qq
}

denMode <- function(d) d$x[which.max(d$y)]

expit <- function (x) 1/(1 + exp(-x))

## computing log ratio, 95% CI and the corresponding p-value of two proportions
## from m.dat
m.2prop.fn <- function(m.dat, nItr=50000, thin=5, nChains=4, nCores=4)
{
    b2p.fit <- sampling(b2p, data=m.dat, iter=nItr, thin=thin, chains=nChains, cores=nCores)
    m <- summary(b2p.fit, pars=c("logRat"), probs=c(0.025,0.975))$summary
    logRat.ci <- m[,c(1,4,5)]

    m <- summary(b2p.fit, pars=c("p"), probs=c(0.025,0.975))$summary
    p.ci <- m[,c(1,4,5)]
    rownames(p.ci) <- c("sPTB","TERM")

    logRat <- as.vector(extract(b2p.fit, pars=c("logRat"), permuted = TRUE)$logRat)
    ## str(logRat)
    ## qqnorm(logRat)
    ## qqline(logRat, col=2)

    logRat.bar <- median(logRat)

    pFn <- function(p)
    {
        ret <- 0
        ##p2 <- p/2
        p2 <- p
        mu.q <- quantile(logRat, probs=c(p2,1-p2))
        if ( logRat.bar > 0 )
        {
            ret <- mu.q[1]
        }
        else
        {
            ret <- mu.q[2]
        }
        ret
    }

    fn.pval <- NA
    if ( sign(pFn(0)*pFn(1))==-1 )
    {
        fn.pval <- uniroot(pFn, c(0, 1))$root
    }


    norm.pval <- NA
    if (  logRat.bar < 0 )
    {
        norm.pval <- 1 - pnorm(0, logRat.bar, mad(logRat))
    } else {
        norm.pval <- pnorm(0, logRat.bar, mad(logRat))
    }
    ##norm.pval

    if(is.na(fn.pval))
    {
        fn.pval <- norm.pval
    }

    logRat.ci <- c(logRat.ci, max(c(norm.pval,fn.pval)), norm.pval, fn.pval)
    names(logRat.ci)[1] <- "logRat"
    names(logRat.ci)[4] <- "p-value"
    names(logRat.ci)[5] <- "norm p-value"
    names(logRat.ci)[6] <- "fn p-value"

    list(p.ci=p.ci, logRat=logRat, logRat.ci=logRat.ci, m.dat=m.dat) # , m.fit=m.fit
}

sig.thld.fn <- function(r)
{
    thld <- NA
    if (  r$gEff[1] > 0 )
    {
        y.min <- r$r2$y.min
        ci.lower <- r$y.l
        yy <- approxfun(r$x, ci.lower-y.min)
        rg <- range(r$x[is.finite(r$y.l)])
    } else {
        y.max <- r$r2$y.max
        ci.upper <- r$y.u
        yy <- approxfun(r$x, ci.upper-y.max)
        rg <- range(r$x[is.finite(r$y.u)])
    }

    if ( sign(yy(rg[1])*yy(rg[2]))==-1 )
    {
        thld <- uniroot(yy, rg)$root
    }

    thld
}


## here the significance thld is defined as the intersection of y.max with y.u
## when y.med is increasing and y.min with y.l otherwise
sig.thld.fn2 <- function(r)
{
    thld <- NA
    if (  r$gEff[1] > 0 )
    {
        y.max <- r$r2$y.max
        x.min <- r$r2$x.min
        ci.upper <- r$y.u
        yy <- approxfun(r$x, ci.upper-y.max)
        idx <- r$x > x.min
    } else {
        y.min <- r$r2$y.min
        x.max <- r$r2$x.max
        ci.lower <- r$y.u
        yy <- approxfun(r$x, ci.lower-y.min)
        idx <- r$x < x.min
    }

    rg <- range(r$x[idx])
    if ( sign(yy(rg[1])*yy(rg[2]))==-1 )
    {
        thld <- uniroot(yy, rg)$root
    }

    thld
}


## here the significance thld is defined as the intersection of y.med with y1 = y.min + p*gEff
## where p is set to 0.1
sig.thld.fn3 <- function(r, p=0.1)
{
    thld <- NA
    if (  r$gEff[1] > 0 )
    {
        y1 <- expit(r$y.med) - (expit(r$r2$y.min) + p*r$gEff[1])
        yy <- approxfun(r$x, y1)
        idx <- r$x > r$r2$x.min
    } else {
        y.min <- r$r2$y.min
        x.max <- r$r2$x.max
        ci.lower <- r$y.u
        yy <- approxfun(r$x, ci.lower-y.min)
        idx <- r$x < x.min
    }

    rg <- range(r$x[idx])
    ##rg[2] <- rg[2]-0.1
    if ( sign(yy(rg[1])*yy(rg[2]))==-1 )
    {
        thld <- uniroot(yy, rg)$root
    }

    thld
}



old.sig.thld.fn <- function(r)
{
    if (  r$gEff[1] > 0 )
    {
        y.min <- r$r2$y.min
        x.min <- r$r2$x.min
        ci.lower <- r$y.l
        yy <- abs(ci.lower-y.min)
        idx <- r$x > x.min
    } else {
        y.max <- r$r2$y.max
        x.max <- r$r2$x.max
        ci.upper <- r$y.u
        yy <- abs(ci.upper-y.max)
        idx <- r$x > x.max
    }

    i.h <- which.min(yy[idx])
    xx <- r$x[idx]
    ##abline(v=xx[i.h], col="gray80")
    xx[i.h][[1]]
}


aa.sig.thld.fn <- function(r)
{
    x <- r$x + r$bLoad
    y <- r$y.med
    y.l <- r$y.l
    y.u <- r$y.u
    o <- order(x)
    x <- x[o]
    y <- y[o]
    y.l <- y.l[o]
    y.u <- y.u[o]
    plot(x,expit(y), type='l')

    if (  r$gEff[1] > 0 )
    {
        i.min <- which.min(y)
        y.min <- y[i.min]
        x.min <- x[i.min]
        ci.lower <- y.l
        yy <- abs(ci.lower-y.min)
        idx <- x > x.min
    } else {
        i.max <- which.max(y)
        y.max <- y[i.max]
        x.max <- x[i.max]
        ci.upper <- y.u
        yy <- abs(ci.upper-y.max)
        idx <- x > x.max
    }

    i.h <- which.min(yy[idx])
    xx <- x[idx]
    ##abline(v=xx[i.h], col="gray80")
    xx[i.h][[1]]
}

plot.logit.incr <- function(ph, r, sig.thld)
{
    op <- par(mar=c(3.5, 3.5, 2.5, 0.5), mgp=c(2.25,0.6,0),tcl = -0.3)
    title <- paste0(ph," (n=",length(r$y),")")
    plot.logit(r$x, r$y, expit(r$y.med), expit(r$y.l), expit(r$y.u),
               title=title, ylab="pr( sPTB )", xlab="log10( Relative Abundance )")
    abline(v=sig.thld, col="gray", lty=2)
    abline(h=expit(r$r2$y.min), col="gray", lty=2)
    par(op)
}

plot.logit.incr2 <- function(ph, r, sig.thld)
{
    op <- par(mar=c(3.5, 3.5, 2.5, 0.5), mgp=c(2.25,0.6,0),tcl = -0.3)
    title <- paste0(ph," (n=",length(r$y),")")
    plot.logit(r$x, r$y, expit(r$y.med), expit(r$y.l), expit(r$y.u),
               title=title, ylab="pr( sPTB )", xlab="log10( Relative Abundance )")
    abline(v=sig.thld, col="gray", lty=2)
    abline(h=expit(r$r2$y.max), col="gray", lty=2)
    par(op)
}


n.stats <- function(ph, idx1)
{
    x <- pt2[,ph]
    y <- mt$delStatus.n
    s <- mt$subjID
    idx <- idx1 & x > 0
    x <- -log10(x[idx])
    y <- y[idx]
    s <- s[idx]

    d <- density(x)
    mode.log10.ra <- denMode(d)

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
         n.TERM=length(y)-sum(y),
         mode.log10.ra=mode.log10.ra)
}

n.stats.comb <- function(phs, idx1)
{
    x <- pt2[,phs[1]]
    for ( i in 2:length(phs) )
    {
        x <- x + pt2[,phs[i]]
    }
    y <- mt$delStatus.n
    s <- mt$subjID
    idx <- idx1 & x > 0
    x <- -log10(x[idx])
    y <- y[idx]
    s <- s[idx]

    d <- density(x)
    mode.log10.ra <- denMode(d)

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
         n.TERM=length(y)-sum(y),
         mode.log10.ra=mode.log10.ra)
}

n.stats.silva <- function(ph, idx1)
{
    x <- pt.silva[,ph]
    y <- mt$delStatus.n
    s <- mt$subjID
    idx <- idx1 & x > 0
    x <- -log10(x[idx])
    y <- y[idx]
    s <- s[idx]

    d <- density(x)
    mode.log10.ra <- denMode(d)

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
         n.TERM=length(y)-sum(y),
         mode.log10.ra=mode.log10.ra)
}

n.stats.relman <- function(ph)
{
    x <- pt.relman[,ph]
    y <- mt.relman$preterm.n
    s <- mt.relman$subjID
    idx <- x > 0
    x <- -log10(x[idx])
    y <- y[idx]
    s <- s[idx]

    d <- density(x)
    mode.log10.ra <- denMode(d)

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
         n.TERM=length(y)-sum(y),
         mode.log10.ra=mode.log10.ra)
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
    x <- log10(x1) - log10(x2)
    d <- density(x)
    mode.log.rat <- denMode(d)

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
         n.TERM=length(y)-sum(y),
         mode.log.rat=mode.log.rat)
}


n.stats2.no <- function(ph, ref.ph, idx1)
{
    x1 <- ct2[,ph]
    x2 <- ct2[,ref.ph]
    y <- mt$delStatus.n
    s <- mt$subjID
    idx <- idx1 & x1 > 0 & x2 == 0
    x1 <- x1[idx]
    y <- y[idx]
    s <- s[idx]
    x <- log10(x1)
    d <- density(x)
    mode.log.rat <- denMode(d)

    nCaseSubjs <- length(unique(s[y==1]))
    nCtrkSubjs <- length(unique(s[y==0]))

    list(nCaseSubjs=nCaseSubjs,
         nCtrkSubjs=nCtrkSubjs,
         n=length(y), n.sPTB=sum(y),
         n.TERM=length(y)-sum(y),
         mode.log.rat=mode.log.rat)
}


spmrf.get.data <- function(data)
{
    ## Check if y is missing
    if (is.null(data$y)) stop("Data list must contain observation variable named y")

    ## Check if N is missing
    if (is.null(data$N)){
        if (!is.null(data$J)) {
            data$N <- length(data$y)
            data <- data[names(data)!="J"]  #get rid of J input
        } else if (is.null(data$J)) {
            data$N <- length(data$y)
        }
    }
    if (length(data$y)!=data$N) stop("number of observations N must equal length of y")
    ## Check if a xvar1 is specified
    ## if missing, assume equal spaced grid of length N
    if (is.null(data$xvar1)) {
        data$xvar1 <- 1:data$N
    }
    ## Find unique values of xvar1
    uxv1 <- unique(data$xvar1)
    ## Get ranks of locations
    ruxv1 <- rank(uxv1)
    ## order by ranks
    m.xv <- cbind(1:length(uxv1), uxv1, ruxv1)
    m.xv <- m.xv[order(m.xv[,2]),]
    ## get grid cell widths for ordered cells
    duxv1 <- diff(m.xv[,2])
    suxv1 <- sort(uxv1)
    ## create mapping of obs to xvar ranks
    rnk.xv <- integer(data$N)
    for (ii in 1:data$N){
        rnk.xv[ii] <- ruxv1[which(uxv1==data$xvar1[ii])]
    }
    ## add ranks and cell widths to data
    data$J <- length(uxv1)  #this is number of grid cells
    data$duxvar1 <- duxv1  #length
    data$xrank1 <- rnk.xv

    data
}

## nItr=1000
## nModels=3
## thin=1
## nChains=1
## nCores=1
## alpha=0.05
## picsDir="../pics/coeffects_dir/"

## logistic regression spline model for risk of sPTB as a function of log10
## absolute abundances of a given phylotype
aa.spmrf.bernoulli.o2.fn <- function(ph, nModels=5, nItr=2000, thin=1, nChains=1, nCores=1,
                                    alpha=0.05,
                                    aa.min.lim=NA,
                                    aa.max.lim=NA,
                                    picsDir="../pics/aa_spmrf_o2_gEff_dir/target_dir/")
{
    plow <- alpha/2
    phigh <- 1 - alpha/2

    ##j0 <- grep(ph, selPhs)
    x <- pt2[,ph]
    bLoad <- mt$bLoad
    y <- mt$delStatus.n
    ## idx1 <- mt$delStatus.n != 2 #& mt$nullip==1 & mt$visit=="V3" & !is.na(mt$race) & mt$raceChar=="black"
    x <- x[idx1]
    y <- y[idx1]
    bLoad <- bLoad[idx1]
    idx <- x > 0
    x <- x[idx]
    y <- y[idx]
    bLoad <- bLoad[idx]
    length(x)
    q <- quantile(x, probs=c(0.005, 0.995)) # removing potential outliers
    idx <- x > q[1] & x < q[2]
    x <- x[idx]
    y <- y[idx]
    bLoad <- bLoad[idx]
    nx <- length(x)
    sigma <- nexpFn(log10(x))
    x <- log10(x)
    bLoad <- log10(bLoad)
    aa <- x + bLoad
    if ( is.numeric(aa.min.lim) ){
        idx <- aa > aa.min.lim
        x <- x[idx]
        nx <- length(x)
        y <- y[idx]
        bLoad <- bLoad[idx]
        sigma <- sigma[idx]
        aa <- aa[idx]
    }
    if ( is.numeric(aa.max.lim) ){
        idx <- aa < aa.max.lim
        x <- x[idx]
        nx <- length(x)
        y <- y[idx]
        bLoad <- bLoad[idx]
        sigma <- sigma[idx]
        aa <- aa[idx]
    }

    ## o <- order(x)
    ## y <- y[o]
    ## x <- x[o]
    ## bLoad <- bLoad[o]
    ## sigma <- sigma[o]

    m.o2 <- list()
    ##m.o2.dat <- list()
    y.med <- list()
    y.mad <- list()
    y.l <- list()
    y.u <- list()
    for ( i in seq(nModels) )
    {
        print(paste(ph,i))
        z <- numeric(length(x))
        for ( j in seq(nx) )
        {
            z[j] <- rnorm(1, mean=x[j], sd=sigma[j])
        }
        z <- z + bLoad
        o <- order(z)
        ##m.o2.dat[[i]] <- list(o=o, z=z)
        sy <- y[o]
        z <- z[o]
        ##bLoad <- bLoad[o]
        m.dat <- list(y=sy, xvar1=z, J=length(y))
        m.dat2 <- spmrf.get.data(m.dat)
        m.o2[[i]] <- sampling(spmrf.bernoulli.o2.model, data=m.dat2, control=list(adapt_delta=0.96, max_treedepth=12), iter=nItr, thin=thin, chains=nChains, cores=nCores)
        theta <- rstan::extract(m.o2[[i]], "theta")[[1]]
        y.med[[i]] <- approx(z, apply(theta, 2, median), xout=aa)
        y.mad[[i]] <- approx(z, apply(theta, 2, mad), xout=aa)
        y.l[[i]] <- approx(z, apply(theta, 2, quantile, probs = plow), xout=aa)
        y.u[[i]] <- approx(z, apply(theta, 2, quantile, probs = phigh), xout=aa)
    }

    y.med.a <- c()
    y.mad.a <- c()
    y.l.a <- c()
    y.u.a <- c()
    for ( i in seq(nModels) )
    {
        y.med.a <- rbind(y.med.a, y.med[[i]]$y)
        y.mad.a <- rbind(y.mad.a, y.mad[[i]]$y)
        y.l.a <- rbind(y.l.a, y.l[[i]]$y)
        y.u.a <- rbind(y.u.a, y.u[[i]]$y)
    }

    non.na.columns <- apply(y.med.a, 2, function(x) !any(is.na(x)))
    y.med.a <- y.med.a[, non.na.columns]
    y.mad.a <- y.mad.a[, non.na.columns]
    y.l.a <- y.l.a[, non.na.columns]
    y.u.a <- y.u.a[, non.na.columns]
    x <- x[non.na.columns]
    y <- y[non.na.columns]
    bLoad <- bLoad[non.na.columns]
    aa <- aa[non.na.columns]

    y.med <- apply(y.med.a, 2, function(x) mean(x, na.rm=T))
    y.mad <- apply(y.mad.a, 2, function(x) mean(x, na.rm=T))
    y.l <- apply(y.l.a, 2, function(x) mean(x, na.rm=T))
    y.u <- apply(y.u.a, 2, function(x) mean(x, na.rm=T))

    r2 <- gEff.pval(x, y.med, y.mad)

    gEff <- numeric(8)
    names(gEff) <- c("gEff", "p-val","q-val","n(TERM)","n(sPTB)", "n(subj TERM)","n(subj sPTB)","mode($\\log_{10}$(rat))")
    r3 <- n.stats(ph, idx1)
    gEff[c(1,2,4:8)] <- c(r2$eff, r2$pval, r3$n.TERM, r3$n.sPTB,
                         r3$nCtrkSubjs, r3$nCaseSubjs, r3$mode.log10.ra)

    o <- order(aa)

    file <- paste0(picsDir, ph, ".pdf")
    pdf(file, width=6, height=6)
    op <- par(mar=c(3.5, 3.5, 2.5, 0.5), mgp=c(2.25,0.6,0),tcl = -0.3)
    title <- paste0(ph," (n=",length(y),")")
    plot.logit(aa[o], y[o], expit(y.med[o]), expit(y.l[o]), expit(y.u[o]), title=title, ylab="pr( sPTB )", xlab="log10( Absolute Abundance )")
    par(op)
    dev.off()

    list(ph=ph,
         x=x, y=y,
         bLoad=bLoad,
         gEff=gEff, r2=r2,
         y.med=y.med,
         y.mad=y.mad,
         y.l=y.l,
         y.u=y.u,
         file=file)
}


## logistic regression spline model for risk of sPTB as a function of log10
## absolute abundances of a group of phylotypes
aa.spmrf.bernoulli.o2.comb.fn <- function(phs, phs.name="", nModels=5, nItr=2000, thin=1, nChains=1, nCores=1,
                                    alpha=0.05,
                                    aa.min.lim=NA,
                                    aa.max.lim=NA,
                                    picsDir="../pics/aa_spmrf_o2_gEff_dir/target_dir/")
{
    plow <- alpha/2
    phigh <- 1 - alpha/2

    x <- pt2[,phs[1]]
    for ( i in 2:length(phs) )
    {
        x <- x + pt2[,phs[i]]
    }
    bLoad <- mt$bLoad
    y <- mt$delStatus.n
    ## idx1 <- mt$delStatus.n != 2 #& mt$nullip==1 & mt$visit=="V3" & !is.na(mt$race) & mt$raceChar=="black"
    x <- x[idx1]
    y <- y[idx1]
    bLoad <- bLoad[idx1]
    idx <- !is.na(bLoad) & x > 0
    x <- x[idx]
    y <- y[idx]
    bLoad <- bLoad[idx]
    length(x)
    q <- quantile(x, probs=c(0.005, 0.995)) # removing potential outliers
    idx <- x > q[1] & x < q[2]
    x <- x[idx]
    y <- y[idx]
    bLoad <- bLoad[idx]
    nx <- length(x)
    sigma <- nexpFn(log10(x))
    x <- log10(x)
    bLoad <- log10(bLoad)
    aa <- x + bLoad
    if ( is.numeric(aa.min.lim) ){
        idx <- aa > aa.min.lim
        x <- x[idx]
        nx <- length(x)
        y <- y[idx]
        bLoad <- bLoad[idx]
        sigma <- sigma[idx]
        aa <- aa[idx]
    }
    if ( is.numeric(aa.max.lim) ){
        idx <- aa < aa.max.lim
        x <- x[idx]
        nx <- length(x)
        y <- y[idx]
        bLoad <- bLoad[idx]
        sigma <- sigma[idx]
        aa <- aa[idx]
    }

    m.o2 <- list()
    y.med <- list()
    y.mad <- list()
    y.l <- list()
    y.u <- list()
    for ( i in seq(nModels) )
    {
        print(paste(phs.name,i))
        z <- numeric(length(x))
        for ( j in seq(nx) )
        {
            z[j] <- rnorm(1, mean=x[j], sd=sigma[j])
        }
        z <- z + bLoad
        o <- order(z)
        ##m.o2.dat[[i]] <- list(o=o, z=z)
        sy <- y[o]
        z <- z[o]
        ##bLoad <- bLoad[o]
        m.dat <- list(y=sy, xvar1=z, J=length(y))
        m.dat2 <- spmrf.get.data(m.dat)
        m.o2[[i]] <- sampling(spmrf.bernoulli.o2.model, data=m.dat2, control=list(adapt_delta=0.96, max_treedepth=12), iter=nItr, thin=thin, chains=nChains, cores=nCores)
        theta <- rstan::extract(m.o2[[i]], "theta")[[1]]
        y.med[[i]] <- approx(z, apply(theta, 2, median), xout=aa)
        y.mad[[i]] <- approx(z, apply(theta, 2, mad), xout=aa)
        y.l[[i]] <- approx(z, apply(theta, 2, quantile, probs = plow), xout=aa)
        y.u[[i]] <- approx(z, apply(theta, 2, quantile, probs = phigh), xout=aa)
    }

    y.med.a <- c()
    y.mad.a <- c()
    y.l.a <- c()
    y.u.a <- c()
    for ( i in seq(nModels) )
    {
        y.med.a <- rbind(y.med.a, y.med[[i]]$y)
        y.mad.a <- rbind(y.mad.a, y.mad[[i]]$y)
        y.l.a <- rbind(y.l.a, y.l[[i]]$y)
        y.u.a <- rbind(y.u.a, y.u[[i]]$y)
    }

    non.na.columns <- apply(y.med.a, 2, function(x) !any(is.na(x)))
    y.med.a <- y.med.a[, non.na.columns]
    y.mad.a <- y.mad.a[, non.na.columns]
    y.l.a <- y.l.a[, non.na.columns]
    y.u.a <- y.u.a[, non.na.columns]
    x <- x[non.na.columns]
    y <- y[non.na.columns]
    bLoad <- bLoad[non.na.columns]
    aa <- aa[non.na.columns]

    y.med <- apply(y.med.a, 2, function(x) mean(x, na.rm=T))
    y.mad <- apply(y.mad.a, 2, function(x) mean(x, na.rm=T))
    y.l <- apply(y.l.a, 2, function(x) mean(x, na.rm=T))
    y.u <- apply(y.u.a, 2, function(x) mean(x, na.rm=T))

    r2 <- gEff.pval(x, y.med, y.mad)

    gEff <- numeric(8)
    names(gEff) <- c("gEff", "p-val","q-val","n(TERM)","n(sPTB)", "n(subj TERM)","n(subj sPTB)","mode($\\log_{10}$(rat))")
    r3 <- n.stats.comb(phs, idx1)
    gEff[c(1,2,4:8)] <- c(r2$eff, r2$pval, r3$n.TERM, r3$n.sPTB,
                         r3$nCtrkSubjs, r3$nCaseSubjs, r3$mode.log10.ra)

    o <- order(aa)

    file <- paste0(picsDir, gsub(" ","_",phs.name), ".pdf")
    pdf(file, width=6, height=6)
    op <- par(mar=c(3.5, 3.5, 2.5, 0.5), mgp=c(2.25,0.6,0),tcl = -0.3)
    title <- paste0(phs.name," (n=",length(y),")")
    plot.logit(aa[o], y[o], expit(y.med[o]), expit(y.l[o]), expit(y.u[o]), title=title, ylab="pr( sPTB )", xlab="log10( Absolute Abundance )")
    par(op)
    dev.off()

    list(phs=phs,
         phs.name=phs.name,
         x=x, y=y,
         bLoad=bLoad,
         gEff=gEff, r2=r2,
         y.med=y.med,
         y.mad=y.mad,
         y.l=y.l,
         y.u=y.u,
         file=file)
}



spmrf.bernoulli.o2.fn.xy <- function(x, y, nModels=5, nItr=2000, thin=1, nChains=1, nCores=1, alpha=0.05, showMessages=showMessages) # picsDir="../pics/spmrf_o2_RA_gEff_dir/target_dir/"
{
    plow <- alpha/2
    phigh <- 1 - alpha/2

    sigma <- nexpFn(x) # x is assumed to be log10 relAb
    o <- order(x)
    y <- y[o]
    x <- x[o]
    sigma <- sigma[o]

    if ( nModels > 1 )
    {
        m.o2 <- list()
        m.o2.dat <- list()
        y.med <- list()
        y.mad <- list()
        y.l <- list()
        y.u <- list()
        for ( i in seq(nModels) )
        {
            print(paste(i))
            z <- numeric(length(x))
            for ( j in seq(x) )
            {
                z[j] <- rnorm(1, mean=x[j], sd=sigma[j])
            }
            o <- order(z)
            m.o2.dat[[i]] <- list(o=o, z=z)
            sy <- y[o]
            z <- z[o]
            m.dat <- list(y=sy, xvar1=z, J=length(y))
            m.dat2 <- spmrf.get.data(m.dat)
            m.o2[[i]] <- sampling(spmrf.bernoulli.o2.model, data=m.dat2, iter=nItr, thin=thin, chains=nChains, cores=nCores, show_messages=showMessages) # control=list(adapt_delta=0.96, max_treedepth=12),
            theta <- rstan::extract(m.o2[[i]], "theta")[[1]]
            y.med[[i]] <- approx(z, apply(theta, 2, median), xout=x)
            y.mad[[i]] <- approx(z, apply(theta, 2, mad), xout=x)
            y.l[[i]] <- approx(z, apply(theta, 2, quantile, probs = plow), xout=x)
            y.u[[i]] <- approx(z, apply(theta, 2, quantile, probs = phigh), xout=x)
        }

        y.med.a <- c()
        y.mad.a <- c()
        y.l.a <- c()
        y.u.a <- c()
        for ( i in seq(nModels) )
        {
            y.med.a <- rbind(y.med.a, y.med[[i]]$y)
            y.mad.a <- rbind(y.mad.a, y.mad[[i]]$y)
            y.l.a <- rbind(y.l.a, y.l[[i]]$y)
            y.u.a <- rbind(y.u.a, y.u[[i]]$y)
        }

        non.na.columns <- apply(y.med.a, 2, function(x) !any(is.na(x)))
        y.med.a <- y.med.a[, non.na.columns]
        y.mad.a <- y.mad.a[, non.na.columns]
        y.l.a <- y.l.a[, non.na.columns]
        y.u.a <- y.u.a[, non.na.columns]
        x <- x[non.na.columns]
        y <- y[non.na.columns]

        y.med <- apply(y.med.a, 2, function(x) mean(x, na.rm=T))
        y.mad <- apply(y.mad.a, 2, function(x) mean(x, na.rm=T))
        y.l <- apply(y.l.a, 2, function(x) mean(x, na.rm=T))
        y.u <- apply(y.u.a, 2, function(x) mean(x, na.rm=T))
    } else {
        m.dat <- list(y=y, xvar1=x, J=length(y))
        m.dat2 <- spmrf.get.data(m.dat)
        m <- sampling(spmrf.bernoulli.o2.model, data=m.dat2, iter=nItr, thin=thin, chains=nChains, cores=nCores, show_messages=showMessages)  # control=list(adapt_delta=0.96, max_treedepth=12)
        theta <- rstan::extract(m, "theta")[[1]]
        y.med <- apply(theta, 2, median)
        y.mad <- apply(theta, 2, mad)
        y.l <- apply(theta, 2, quantile, probs = plow)
        y.u <- apply(theta, 2, quantile, probs = phigh)

    }

    idx <- seq(y.med)
    r2 <- gEff.pval(x[idx], y.med, y.mad)

    ## gEff <- numeric(8)
    ## names(gEff) <- c("gEff", "p-val","q-val","n(TERM)","n(sPTB)", "n(subj TERM)","n(subj sPTB)","mode($\\log_{10}$(rat))")
    ## r3 <- n.stats(ph, idx1)
    ## gEff[c(1,2,4:8)] <- c(r2$eff, r2$pval.min, r3$n.TERM, r3$n.sPTB,
    ##                      r3$nCtrkSubjs, r3$nCaseSubjs, r3$mode.log10.ra)

    list(x=x,
         y=y,
         ##gEff=gEff,
         r2=r2,
         y.med=y.med,
         y.mad=y.mad,
         y.l=y.l,
         y.u=y.u)
}

spmrf.bernoulli.o3.fn.xy <- function(x, y, nModels=5, nItr=2000, thin=1, nChains=1, nCores=1, alpha=0.05, showMessages=showMessages) # picsDir="../pics/spmrf_o2_RA_gEff_dir/target_dir/"
{
    plow <- alpha/2
    phigh <- 1 - alpha/2

    sigma <- nexpFn(x) # x is assumed to be log10 relAb
    o <- order(x)
    y <- y[o]
    x <- x[o]
    sigma <- sigma[o]

    if ( nModels > 1 )
    {
        m.o3 <- list()
        m.o3.dat <- list()
        y.med <- list()
        y.mad <- list()
        y.l <- list()
        y.u <- list()
        for ( i in seq(nModels) )
        {
            print(paste(i))
            z <- numeric(length(x))
            for ( j in seq(x) )
            {
                z[j] <- rnorm(1, mean=x[j], sd=sigma[j])
            }
            o <- order(z)
            m.o3.dat[[i]] <- list(o=o, z=z)
            sy <- y[o]
            z <- z[o]
            m.dat <- list(y=sy, xvar1=z, J=length(y))
            m.dat2 <- spmrf.get.data(m.dat)
            m.o3[[i]] <- sampling(spmrf.bernoulli.o3.model, data=m.dat2, iter=nItr, thin=thin, chains=nChains, cores=nCores, show_messages=showMessages) # control=list(adapt_delta=0.96, max_treedepth=12),
            theta <- rstan::extract(m.o3[[i]], "theta")[[1]]
            y.med[[i]] <- approx(z, apply(theta, 2, median), xout=x)
            y.mad[[i]] <- approx(z, apply(theta, 2, mad), xout=x)
            y.l[[i]] <- approx(z, apply(theta, 2, quantile, probs = plow), xout=x)
            y.u[[i]] <- approx(z, apply(theta, 2, quantile, probs = phigh), xout=x)
        }

        y.med.a <- c()
        y.mad.a <- c()
        y.l.a <- c()
        y.u.a <- c()
        for ( i in seq(nModels) )
        {
            y.med.a <- rbind(y.med.a, y.med[[i]]$y)
            y.mad.a <- rbind(y.mad.a, y.mad[[i]]$y)
            y.l.a <- rbind(y.l.a, y.l[[i]]$y)
            y.u.a <- rbind(y.u.a, y.u[[i]]$y)
        }

        non.na.columns <- apply(y.med.a, 2, function(x) !any(is.na(x)))
        y.med.a <- y.med.a[, non.na.columns]
        y.mad.a <- y.mad.a[, non.na.columns]
        y.l.a <- y.l.a[, non.na.columns]
        y.u.a <- y.u.a[, non.na.columns]
        x <- x[non.na.columns]
        y <- y[non.na.columns]

        y.med <- apply(y.med.a, 2, function(x) mean(x, na.rm=T))
        y.mad <- apply(y.mad.a, 2, function(x) mean(x, na.rm=T))
        y.l <- apply(y.l.a, 2, function(x) mean(x, na.rm=T))
        y.u <- apply(y.u.a, 2, function(x) mean(x, na.rm=T))
    } else {
        m.dat <- list(y=y, xvar1=x, J=length(y))
        m.dat2 <- spmrf.get.data(m.dat)
        m <- sampling(spmrf.bernoulli.o3.model, data=m.dat2, iter=nItr, thin=thin, chains=nChains, cores=nCores, show_messages=showMessages)  # control=list(adapt_delta=0.96, max_treedepth=12)
        theta <- rstan::extract(m, "theta")[[1]]
        y.med <- apply(theta, 2, median)
        y.mad <- apply(theta, 2, mad)
        y.l <- apply(theta, 2, quantile, probs = plow)
        y.u <- apply(theta, 2, quantile, probs = phigh)

    }

    idx <- seq(y.med)
    r2 <- gEff.pval(x[idx], y.med, y.mad)

    ## gEff <- numeric(8)
    ## names(gEff) <- c("gEff", "p-val","q-val","n(TERM)","n(sPTB)", "n(subj TERM)","n(subj sPTB)","mode($\\log_{10}$(rat))")
    ## r3 <- n.stats(ph, idx1)
    ## gEff[c(1,2,4:8)] <- c(r2$eff, r2$pval.min, r3$n.TERM, r3$n.sPTB,
    ##                      r3$nCtrkSubjs, r3$nCaseSubjs, r3$mode.log10.ra)

    list(x=x,
         y=y,
         ##gEff=gEff,
         r2=r2,
         y.med=y.med,
         y.mad=y.mad,
         y.l=y.l,
         y.u=y.u)
}



spmrf.bernoulli.o2.fn <- function(ph, nModels=5, nItr=2000, thin=1, nChains=1, nCores=1,
                                    alpha=0.05, picsDir="../pics/spmrf_o2_RA_gEff_dir/target_dir/")
{
    plow <- alpha/2
    phigh <- 1 - alpha/2

    ##j0 <- grep(ph, selPhs)
    x <- ct2[,ph]
    n <- rowSums(ct2)
    y <- mt$delStatus.n
    ## idx1 <- mt$delStatus.n != 2 #& mt$nullip==1 & mt$visit=="V3" & !is.na(mt$race) & mt$raceChar=="black"
    x <- x[idx1]
    y <- y[idx1]
    n <- n[idx1]
    idx <- x > 0
    x <- x[idx]
    y <- y[idx]
    n <- n[idx]
    length(x)
    q <- quantile(x/n, probs=c(0.005, 0.995)) # removing potential outliers
    idx <- x/n > q[1] & x/n < q[2]
    x <- x[idx]
    y <- y[idx]
    n <- n[idx]
    nx <- length(x)
    sigma <- nexpFn(log10(x/n))
    x <- log10(x/n)
    o <- order(x)
    y <- y[o]
    x <- x[o]
    sigma <- sigma[o]

    m.o2 <- list()
    m.o2.dat <- list()
    y.med <- list()
    y.mad <- list()
    y.l <- list()
    y.u <- list()
    for ( i in seq(nModels) )
    {
        print(paste(ph,i))
        z <- numeric(length(x))
        for ( j in seq(nx) )
        {
            z[j] <- rnorm(1, mean=x[j], sd=sigma[j])
        }
        o <- order(z)
        m.o2.dat[[i]] <- list(o=o, z=z)
        sy <- y[o]
        z <- z[o]
        m.dat <- list(y=sy, xvar1=z, J=length(y))
        m.dat2 <- spmrf.get.data(m.dat)
        m.o2[[i]] <- sampling(spmrf.bernoulli.o2.model, data=m.dat2, control=list(adapt_delta=0.96, max_treedepth=12), iter=nItr, thin=thin, chains=nChains, cores=nCores)
        theta <- rstan::extract(m.o2[[i]], "theta")[[1]]
        y.med[[i]] <- approx(z, apply(theta, 2, median), xout=x)
        y.mad[[i]] <- approx(z, apply(theta, 2, mad), xout=x)
        y.l[[i]] <- approx(z, apply(theta, 2, quantile, probs = plow), xout=x)
        y.u[[i]] <- approx(z, apply(theta, 2, quantile, probs = phigh), xout=x)
    }

    y.med.a <- c()
    y.mad.a <- c()
    y.l.a <- c()
    y.u.a <- c()
    for ( i in seq(nModels) )
    {
        y.med.a <- rbind(y.med.a, y.med[[i]]$y)
        y.mad.a <- rbind(y.mad.a, y.mad[[i]]$y)
        y.l.a <- rbind(y.l.a, y.l[[i]]$y)
        y.u.a <- rbind(y.u.a, y.u[[i]]$y)
    }

    non.na.columns <- apply(y.med.a, 2, function(x) !any(is.na(x)))
    y.med.a <- y.med.a[, non.na.columns]
    y.mad.a <- y.mad.a[, non.na.columns]
    y.l.a <- y.l.a[, non.na.columns]
    y.u.a <- y.u.a[, non.na.columns]
    x <- x[non.na.columns]
    y <- y[non.na.columns]

    y.med <- apply(y.med.a, 2, function(x) mean(x, na.rm=T))
    y.mad <- apply(y.mad.a, 2, function(x) mean(x, na.rm=T))
    y.l <- apply(y.l.a, 2, function(x) mean(x, na.rm=T))
    y.u <- apply(y.u.a, 2, function(x) mean(x, na.rm=T))

    r2 <- gEff.pval(x, y.med, y.mad)

    gEff <- numeric(8)
    names(gEff) <- c("gEff", "p-val","q-val","n(TERM)","n(sPTB)", "n(subj TERM)","n(subj sPTB)","mode($\\log_{10}$(rat))")
    r3 <- n.stats(ph, idx1)
    gEff[c(1,2,4:8)] <- c(r2$eff, r2$pval, r3$n.TERM, r3$n.sPTB,
                         r3$nCtrkSubjs, r3$nCaseSubjs, r3$mode.log10.ra)

    file <- paste0(picsDir, ph, ".pdf")
    pdf(file, width=6, height=6)
    op <- par(mar=c(3.5, 3.5, 2.5, 0.5), mgp=c(2.25,0.6,0),tcl = -0.3)
    title <- paste0(ph," (n=",length(y),")")
    plot.logit(x, y, expit(y.med), expit(y.l), expit(y.u), title=title, ylab="pr( sPTB )", xlab="log10( Relative Abundance )")
    par(op)
    dev.off()

    list(ph=ph,
         x=x, y=y,
         gEff=gEff, r2=r2,
         y.med=y.med,
         y.mad=y.mad,
         y.l=y.l,
         y.u=y.u,
         file=file)
}


## version of spmrf.bernoulli.o2.fn()
## accepting a vector of phylotypes
spmrf.bernoulli.o2.comb.fn <- function(phs, phs.name="", nModels=5, nItr=2000, thin=1, nChains=1, nCores=1,
                                      alpha=0.05, picsDir="../pics/spmrf_o2_RA_gEff_dir/target_dir/")
{
    plow <- alpha/2
    phigh <- 1 - alpha/2

    ##j0 <- grep(ph, selPhs)
    x <- ct2[,phs[1]]
    for ( i in 2:length(phs) )
    {
        x <- x + ct2[,phs[i]]
    }
    n <- rowSums(ct2)
    y <- mt$delStatus.n
    ## idx1 <- mt$delStatus.n != 2 #& mt$nullip==1 & mt$visit=="V3" & !is.na(mt$race) & mt$raceChar=="black"
    x <- x[idx1]
    y <- y[idx1]
    n <- n[idx1]
    idx <- x > 0
    x <- x[idx]
    y <- y[idx]
    n <- n[idx]
    length(x)
    q <- quantile(x/n, probs=c(0.005, 0.995)) # removing potential outliers
    idx <- x/n > q[1] & x/n < q[2]
    x <- x[idx]
    y <- y[idx]
    n <- n[idx]
    nx <- length(x)
    sigma <- nexpFn(log10(x/n))
    x <- log10(x/n)
    o <- order(x)
    y <- y[o]
    x <- x[o]
    sigma <- sigma[o]

    m.o2 <- list()
    m.o2.dat <- list()
    y.med <- list()
    y.mad <- list()
    y.l <- list()
    y.u <- list()
    for ( i in seq(nModels) )
    {
        print(paste(phs.name,i))
        z <- numeric(length(x))
        for ( j in seq(nx) )
        {
            z[j] <- rnorm(1, mean=x[j], sd=sigma[j])
        }
        o <- order(z)
        m.o2.dat[[i]] <- list(o=o, z=z)
        sy <- y[o]
        z <- z[o]
        m.dat <- list(y=sy, xvar1=z, J=length(y))
        m.dat2 <- spmrf.get.data(m.dat)
        m.o2[[i]] <- sampling(spmrf.bernoulli.o2.model, data=m.dat2, control=list(adapt_delta=0.96, max_treedepth=12), iter=nItr, thin=thin, chains=nChains, cores=nCores)
        theta <- rstan::extract(m.o2[[i]], "theta")[[1]]
        y.med[[i]] <- approx(z, apply(theta, 2, median), xout=x)
        y.mad[[i]] <- approx(z, apply(theta, 2, mad), xout=x)
        y.l[[i]] <- approx(z, apply(theta, 2, quantile, probs = plow), xout=x)
        y.u[[i]] <- approx(z, apply(theta, 2, quantile, probs = phigh), xout=x)
    }

    y.med.a <- c()
    y.mad.a <- c()
    y.l.a <- c()
    y.u.a <- c()
    for ( i in seq(nModels) )
    {
        y.med.a <- rbind(y.med.a, y.med[[i]]$y)
        y.mad.a <- rbind(y.mad.a, y.mad[[i]]$y)
        y.l.a <- rbind(y.l.a, y.l[[i]]$y)
        y.u.a <- rbind(y.u.a, y.u[[i]]$y)
    }

    non.na.columns <- apply(y.med.a, 2, function(x) !any(is.na(x)))
    y.med.a <- y.med.a[, non.na.columns]
    y.mad.a <- y.mad.a[, non.na.columns]
    y.l.a <- y.l.a[, non.na.columns]
    y.u.a <- y.u.a[, non.na.columns]
    x <- x[non.na.columns]
    y <- y[non.na.columns]

    y.med <- apply(y.med.a, 2, function(x) mean(x, na.rm=T))
    y.mad <- apply(y.mad.a, 2, function(x) mean(x, na.rm=T))
    y.l <- apply(y.l.a, 2, function(x) mean(x, na.rm=T))
    y.u <- apply(y.u.a, 2, function(x) mean(x, na.rm=T))

    r2 <- gEff.pval(x, y.med, y.mad)

    gEff <- numeric(8)
    names(gEff) <- c("gEff", "p-val","q-val","n(TERM)","n(sPTB)", "n(subj TERM)","n(subj sPTB)","mode($\\log_{10}$(rat))")
    r3 <- n.stats.comb(phs, idx1)
    gEff[c(1,2,4:8)] <- c(r2$eff, r2$pval, r3$n.TERM, r3$n.sPTB,
                         r3$nCtrkSubjs, r3$nCaseSubjs, r3$mode.log10.ra)

    file <- paste0(picsDir, gsub(" ","_",phs.name), ".pdf")
    pdf(file, width=6, height=6)
    op <- par(mar=c(3.5, 3.5, 2.5, 0.5), mgp=c(2.25,0.6,0),tcl = -0.3)
    title <- paste0(phs.name," (n=",length(y),")")
    plot.logit(x, y, expit(y.med), expit(y.l), expit(y.u), title=title, ylab="pr( sPTB )", xlab="log10( Relative Abundance )")
    par(op)
    dev.off()

    list(phs=phs,
         x=x, y=y,
         gEff=gEff, r2=r2,
         y.med=y.med,
         y.mad=y.mad,
         y.l=y.l,
         y.u=y.u,
         file=file)
}


spmrf.bernoulli.o2.silva.fn <- function(ph, nModels=5, nItr=2000, thin=1, nChains=1, nCores=1,
                                       alpha=0.05, picsDir="../pics/spmrf_o2_RA_gEff_dir/target_dir/")
{
    plow <- alpha/2
    phigh <- 1 - alpha/2

    ##j0 <- grep(ph, selPhs)
    x <- ct.silva[,ph]
    n <- rowSums(ct.silva)
    y <- mt$delStatus.n
    ## idx1 <- mt$delStatus.n != 2 #& mt$nullip==1 & mt$visit=="V3" & !is.na(mt$race) & mt$raceChar=="black"
    x <- x[idx1]
    y <- y[idx1]
    n <- n[idx1]
    idx <- x > 0
    x <- x[idx]
    y <- y[idx]
    n <- n[idx]
    length(x)
    q <- quantile(x/n, probs=c(0.005, 0.995)) # removing potential outliers
    idx <- x/n > q[1] & x/n < q[2]
    x <- x[idx]
    y <- y[idx]
    n <- n[idx]
    nx <- length(x)
    sigma <- nexpFn(log10(x/n))
    x <- log10(x/n)
    o <- order(x)
    y <- y[o]
    x <- x[o]
    sigma <- sigma[o]

    m.o2 <- list()
    m.o2.dat <- list()
    y.med <- list()
    y.mad <- list()
    y.l <- list()
    y.u <- list()
    for ( i in seq(nModels) )
    {
        print(paste(ph,i))
        z <- numeric(length(x))
        for ( j in seq(nx) )
        {
            z[j] <- rnorm(1, mean=x[j], sd=sigma[j])
        }
        o <- order(z)
        m.o2.dat[[i]] <- list(o=o, z=z)
        sy <- y[o]
        z <- z[o]
        m.dat <- list(y=sy, xvar1=z, J=length(y))
        m.dat2 <- spmrf.get.data(m.dat)
        m.o2[[i]] <- sampling(spmrf.bernoulli.o2.model, data=m.dat2, control=list(adapt_delta=0.96, max_treedepth=12), iter=nItr, thin=thin, chains=nChains, cores=nCores)
        theta <- rstan::extract(m.o2[[i]], "theta")[[1]]
        y.med[[i]] <- approx(z, apply(theta, 2, median), xout=x)
        y.mad[[i]] <- approx(z, apply(theta, 2, mad), xout=x)
        y.l[[i]] <- approx(z, apply(theta, 2, quantile, probs = plow), xout=x)
        y.u[[i]] <- approx(z, apply(theta, 2, quantile, probs = phigh), xout=x)
    }

    y.med.a <- c()
    y.mad.a <- c()
    y.l.a <- c()
    y.u.a <- c()
    for ( i in seq(nModels) )
    {
        y.med.a <- rbind(y.med.a, y.med[[i]]$y)
        y.mad.a <- rbind(y.mad.a, y.mad[[i]]$y)
        y.l.a <- rbind(y.l.a, y.l[[i]]$y)
        y.u.a <- rbind(y.u.a, y.u[[i]]$y)
    }

    non.na.columns <- apply(y.med.a, 2, function(x) !any(is.na(x)))
    y.med.a <- y.med.a[, non.na.columns]
    y.mad.a <- y.mad.a[, non.na.columns]
    y.l.a <- y.l.a[, non.na.columns]
    y.u.a <- y.u.a[, non.na.columns]
    x <- x[non.na.columns]
    y <- y[non.na.columns]

    y.med <- apply(y.med.a, 2, function(x) mean(x, na.rm=T))
    y.mad <- apply(y.mad.a, 2, function(x) mean(x, na.rm=T))
    y.l <- apply(y.l.a, 2, function(x) mean(x, na.rm=T))
    y.u <- apply(y.u.a, 2, function(x) mean(x, na.rm=T))

    r2 <- gEff.pval(x, y.med, y.mad)

    gEff <- numeric(8)
    names(gEff) <- c("gEff", "p-val","q-val","n(TERM)","n(sPTB)", "n(subj TERM)","n(subj sPTB)","mode($\\log_{10}$(rat))")
    r3 <- n.stats.silva(ph, idx1)
    gEff[c(1,2,4:8)] <- c(r2$eff, r2$pval, r3$n.TERM, r3$n.sPTB,
                         r3$nCtrkSubjs, r3$nCaseSubjs, r3$mode.log10.ra)

    file <- paste0(picsDir, ph, ".pdf")
    pdf(file, width=6, height=6)
    op <- par(mar=c(3.5, 3.5, 2.5, 0.5), mgp=c(2.25,0.6,0),tcl = -0.3)
    title <- paste0(ph," (n=",length(y),")")
    plot.logit(x, y, expit(y.med), expit(y.l), expit(y.u), title=title, ylab="pr( sPTB )", xlab="log10( Relative Abundance )")
    par(op)
    dev.off()

    list(ph=ph,
         x=x, y=y,
         gEff=gEff, r2=r2,
         y.med=y.med,
         y.mad=y.mad,
         y.l=y.l,
         y.u=y.u,
         file=file)
}


## nModels=3; nItr=500; thin=1; nChains=1; nCores=1; alpha=0.05; picsDir=all.picsDir

spmrf.bernoulli.o2.relman.fn <- function(ph, nModels=5, nItr=2000, thin=1, nChains=1, nCores=1,
                                       alpha=0.05, picsDir="../pics/spmrf_o2_RA_gEff_dir/target_dir/")
{
    plow <- alpha/2
    phigh <- 1 - alpha/2

    ##j0 <- grep(ph, selPhs)
    x <- ct.relman[,ph]
    n <- rowSums(ct.relman)
    y <- mt.relman$preterm.n
    idx <- x > 0
    x <- x[idx]
    y <- y[idx]
    n <- n[idx]
    length(x)
    q <- quantile(x/n, probs=c(0.005, 0.995)) # removing potential outliers
    idx <- x/n > q[1] & x/n < q[2]
    x <- x[idx]
    y <- y[idx]
    n <- n[idx]
    nx <- length(x)
    sigma <- nexpFn(log10(x/n))
    x <- log10(x/n)
    o <- order(x)
    y <- y[o]
    x <- x[o]
    sigma <- sigma[o]

    m.o2 <- list()
    m.o2.dat <- list()
    y.med <- list()
    y.mad <- list()
    y.l <- list()
    y.u <- list()
    for ( i in seq(nModels) )
    {
        print(paste(ph,i))
        z <- numeric(length(x))
        for ( j in seq(nx) )
        {
            z[j] <- rnorm(1, mean=x[j], sd=sigma[j])
        }
        o <- order(z)
        m.o2.dat[[i]] <- list(o=o, z=z)
        sy <- y[o]
        z <- z[o]
        m.dat <- list(y=sy, xvar1=z, J=length(y))
        m.dat2 <- spmrf.get.data(m.dat)
        m.o2[[i]] <- sampling(spmrf.bernoulli.o2.model, data=m.dat2, control=list(adapt_delta=0.96, max_treedepth=12), iter=nItr, thin=thin, chains=nChains, cores=nCores)
        theta <- rstan::extract(m.o2[[i]], "theta")[[1]]
        y.med[[i]] <- approx(z, apply(theta, 2, median), xout=x)
        y.mad[[i]] <- approx(z, apply(theta, 2, mad), xout=x)
        y.l[[i]] <- approx(z, apply(theta, 2, quantile, probs = plow), xout=x)
        y.u[[i]] <- approx(z, apply(theta, 2, quantile, probs = phigh), xout=x)
    }

    y.med.a <- c()
    y.mad.a <- c()
    y.l.a <- c()
    y.u.a <- c()
    for ( i in seq(nModels) )
    {
        y.med.a <- rbind(y.med.a, y.med[[i]]$y)
        y.mad.a <- rbind(y.mad.a, y.mad[[i]]$y)
        y.l.a <- rbind(y.l.a, y.l[[i]]$y)
        y.u.a <- rbind(y.u.a, y.u[[i]]$y)
    }

    non.na.columns <- apply(y.med.a, 2, function(x) !any(is.na(x)))
    y.med.a <- y.med.a[, non.na.columns]
    y.mad.a <- y.mad.a[, non.na.columns]
    y.l.a <- y.l.a[, non.na.columns]
    y.u.a <- y.u.a[, non.na.columns]
    x <- x[non.na.columns]
    y <- y[non.na.columns]

    y.med <- apply(y.med.a, 2, function(x) mean(x, na.rm=T))
    y.mad <- apply(y.mad.a, 2, function(x) mean(x, na.rm=T))
    y.l <- apply(y.l.a, 2, function(x) mean(x, na.rm=T))
    y.u <- apply(y.u.a, 2, function(x) mean(x, na.rm=T))

    r2 <- gEff.pval(x, y.med, y.mad)

    gEff <- numeric(8)
    names(gEff) <- c("gEff", "p-val","q-val","n(TERM)","n(sPTB)", "n(subj TERM)","n(subj sPTB)","mode($\\log_{10}$(rat))")
    r3 <- n.stats.relman(ph)
    gEff[c(1,2,4:8)] <- c(r2$eff, r2$pval, r3$n.TERM, r3$n.sPTB,
                         r3$nCtrkSubjs, r3$nCaseSubjs, r3$mode.log10.ra)

    file <- paste0(picsDir, gsub("/","",ph), ".pdf")
    pdf(file, width=6, height=6)
    op <- par(mar=c(3.5, 3.5, 2.5, 0.5), mgp=c(2.25,0.6,0),tcl = -0.3)
    title <- paste0(ph," (n=",length(y),")")
    plot.logit(x, y, expit(y.med), expit(y.l), expit(y.u), title=title, ylab="pr( sPTB )", xlab="log10( Relative Abundance )")
    par(op)
    dev.off()

    list(ph=ph,
         x=x, y=y,
         gEff=gEff, r2=r2,
         y.med=y.med,
         y.mad=y.mad,
         y.l=y.l,
         y.u=y.u,
         file=file)
}




aa.spmrf.bernoulli.o2.coeffects.fn <- function(ph, nModels=5, nItr=2000, thin=1, nChains=1, nCores=1,
                                              alpha=0.05,
                                              picsDir="../pics/coeffects_dir/target_dir/")
{
    plow <- alpha/2
    phigh <- 1 - alpha/2

    j0 <- grep(ph, other.phs)
    ph.x <- pt2[,ph]
    ref.ph.x <- pt2[,ref.ph]
    bLoad <- mt$bLoad
    y <- mt$delStatus.n
    idx <- idx1 & ph.x > 0 & ref.ph.x > 0
    y <- y[idx]
    bLoad <- log10(bLoad[idx])
    ph.x <- log10(ph.x[idx])
    ref.ph.x <- log10(ref.ph.x[idx])
    ph.sigma <- nexpFn(ph.x)
    ref.ph.sigma <- nexpFn(ref.ph.x)
    x.log.rat <- ph.x - ref.ph.x
    aa.log.rat <- x.log.rat + bLoad
    ## o <- order(x.log.rat)
    ## x <- x.log.rat[o]
    ## y <- y[o]
    nx <- length(x.log.rat)

    y.med <- list()
    y.mad <- list()
    y.l <- list()
    y.u <- list()
    m.o2 <- list()
    for ( i in seq(nModels) )
    {
        print(paste(ph, "vs", ref.ph, j0, i))
        ph.z <- numeric(nx)
        ref.ph.z <- numeric(nx)
        for ( j in seq(nx) )
        {
            ph.z[j] <- rnorm(1, mean=ph.x[j], sd=ph.sigma[j])
            ref.ph.z[j] <- rnorm(1, mean=ref.ph.x[j], sd=ref.ph.sigma[j])
        }
        z.log.rat <- ph.z - ref.ph.z + bLoad
        o <- order(z.log.rat)
        z.log.rat <- z.log.rat[o]
        z <- z.log.rat
        sy <- y[o]
        m.dat <- list(y=sy, xvar1=z.log.rat, J=length(y))
        m.dat2 <- spmrf.get.data(m.dat)
        m.o2[[i]] <- sampling(spmrf.bernoulli.o2.model, data=m.dat2, control=list(adapt_delta=0.96, max_treedepth=12), iter=nItr, thin=thin, chains=nChains, cores=nCores)
        theta <- rstan::extract(m.o2[[i]], "theta")[[1]]
        y.med[[i]] <- approx(z, apply(theta, 2, median), xout=aa.log.rat)
        y.mad[[i]] <- approx(z, apply(theta, 2, mad), xout=aa.log.rat)
        y.l[[i]] <- approx(z, apply(theta, 2, quantile, probs = plow), xout=aa.log.rat)
        y.u[[i]] <- approx(z, apply(theta, 2, quantile, probs = phigh), xout=aa.log.rat)
    }

    y.med.a <- c()
    y.mad.a <- c()
    y.l.a <- c()
    y.u.a <- c()
    for ( i in seq(nModels) )
    {
        y.med.a <- rbind(y.med.a, y.med[[i]]$y)
        y.mad.a <- rbind(y.mad.a, y.mad[[i]]$y)
        y.l.a <- rbind(y.l.a, y.l[[i]]$y)
        y.u.a <- rbind(y.u.a, y.u[[i]]$y)
    }

    non.na.columns <- apply(y.med.a, 2, function(x) !any(is.na(x)))
    y.med.a <- y.med.a[, non.na.columns]
    y.mad.a <- y.mad.a[, non.na.columns]
    y.l.a <- y.l.a[, non.na.columns]
    y.u.a <- y.u.a[, non.na.columns]
    x <- x[non.na.columns]
    y <- y[non.na.columns]
    bLoad <- bLoad[non.na.columns]
    aa.log.rat <- aa.log.rat[non.na.columns]

    y.med <- apply(y.med.a, 2, function(x) mean(x, na.rm=T))
    y.mad <- apply(y.mad.a, 2, function(x) mean(x, na.rm=T))
    y.l <- apply(y.l.a, 2, function(x) mean(x, na.rm=T))
    y.u <- apply(y.u.a, 2, function(x) mean(x, na.rm=T))

    r2 <- gEff.pval(x, y.med, y.mad)
    r3 <- n.stats2(ph, ref.ph, idx1)

    gEff <- numeric(8)
    names(gEff) <- c("gEff", "p-val","q-val","n(TERM)","n(sPTB)", "n(subj TERM)","n(subj sPTB)","mode($\\log_{10}$(rat))")
    gEff[c(1,2,4:8)] <- c(r2$eff, r2$pval, r3$n.TERM, r3$n.sPTB,
                             r3$nCtrkSubjs, r3$nCaseSubjs, r3$mode.log.rat)

    gEff.min <- numeric(8)
    names(gEff.min) <- c("gEff", "p-val","q-val","n(TERM)","n(sPTB)", "n(subj TERM)","n(subj sPTB)","mode($\\log_{10}$(rat))")
    gEff.min[c(1,2,4:8)] <- c(r2$eff, r2$pval.min, r3$n.TERM, r3$n.sPTB,
                             r3$nCtrkSubjs, r3$nCaseSubjs, r3$mode.log.rat)

    o <- order(aa.log.rat)

    file <- paste0(picsDir,ph,"_",ref.ph,".pdf")
    pdf(file, width=6, height=6)
    op <- par(mar=c(3.5, 3.5, 2.5, 0.5), mgp=c(2.25,0.6,0),tcl = -0.3)
    title <- paste0(ph," vs ",ref.ph)
    plot.logit(aa.log.rat[o], y[o], expit(y.med[o]), expit(y.l[o]), expit(y.u[o]), title=title, ylab="pr( sPTB )", xlab="log10( absAb ratio )")
    par(op)
    dev.off()

    list(ph=ph,
         ref.ph=ref.ph,
         x=x, y=y,
         r2=r2,
         gEff=gEff,
         gEff.min=gEff.min,
         y.med=y.med,
         y.mad=y.mad,
         y.l=y.l,
         y.u=y.u,
         y.med.a=y.med.a,
         y.mad.a=y.mad.a,
         y.l.a=y.l.a,
         y.u.a=y.u.a,
         file=file)
}

spmrf.bernoulli.o2.coeffects.fn <- function(ph, nModels=5, nItr=5000, thin=1, nChains=1, nCores=1,
                                           alpha=0.05,
                                           thld.min=NA,
                                           thld.max=NA,
                                           picsDir="../pics/coeffects_dir/target_dir/")
{
    plow <- alpha/2
    phigh <- 1 - alpha/2

    j0 <- grep(ph, other.phs)
    ph.x <- pt2[,ph]
    ref.ph.x <- pt2[,ref.ph]
    y <- mt$delStatus.n
    idx <- idx1 & ph.x > 0 & ref.ph.x > 0
    y <- y[idx]
    ph.x <- log10(ph.x[idx])
    ref.ph.x <- log10(ref.ph.x[idx])
    if ( !is.na(thld.min) )
    {
        idx <- ref.ph.x > thld.min
        y <- y[idx]
        ph.x <- ph.x[idx]
        ref.ph.x <- ref.ph.x[idx]
    }
    if ( !is.na(thld.max) )
    {
        idx <- ref.ph.x < thld.max
        y <- y[idx]
        ph.x <- ph.x[idx]
        ref.ph.x <- ref.ph.x[idx]
    }
    ph.sigma <- nexpFn(ph.x)
    ref.ph.sigma <- nexpFn(ref.ph.x)
    x.log.rat <- ph.x - ref.ph.x
    o <- order(x.log.rat)
    x <- x.log.rat[o]
    ## y <- y[o]
    nx <- length(x.log.rat)

    y.med <- list()
    y.mad <- list()
    y.l <- list()
    y.u <- list()
    m.o2 <- list()
    for ( i in seq(nModels) )
    {
        print(paste(ph, "vs", ref.ph, j0, i))
        ph.z <- numeric(nx)
        ref.ph.z <- numeric(nx)
        for ( j in seq(nx) )
        {
            ph.z[j] <- rnorm(1, mean=ph.x[j], sd=ph.sigma[j])
            ref.ph.z[j] <- rnorm(1, mean=ref.ph.x[j], sd=ref.ph.sigma[j])
        }
        z.log.rat <- ph.z - ref.ph.z
        o <- order(z.log.rat)
        z.log.rat <- z.log.rat[o]
        z <- z.log.rat
        sy <- y[o]
        m.dat <- list(y=sy, xvar1=z.log.rat, J=length(y))
        m.dat2 <- spmrf.get.data(m.dat)
        m.o2[[i]] <- sampling(spmrf.bernoulli.o2.model, data=m.dat2, control=list(adapt_delta=0.96, max_treedepth=12), iter=nItr, thin=thin, chains=nChains, cores=nCores)
        theta <- rstan::extract(m.o2[[i]], "theta")[[1]]
        y.med[[i]] <- approx(z, apply(theta, 2, median), xout=x)
        y.mad[[i]] <- approx(z, apply(theta, 2, mad), xout=x)
        y.l[[i]] <- approx(z, apply(theta, 2, quantile, probs = plow), xout=x)
        y.u[[i]] <- approx(z, apply(theta, 2, quantile, probs = phigh), xout=x)
    }

    y.med.a <- c()
    y.mad.a <- c()
    y.l.a <- c()
    y.u.a <- c()
    for ( i in seq(nModels) )
    {
        y.med.a <- rbind(y.med.a, y.med[[i]]$y)
        y.mad.a <- rbind(y.mad.a, y.mad[[i]]$y)
        y.l.a <- rbind(y.l.a, y.l[[i]]$y)
        y.u.a <- rbind(y.u.a, y.u[[i]]$y)
    }

    non.na.columns <- apply(y.med.a, 2, function(x) !any(is.na(x)))
    y.med.a <- y.med.a[, non.na.columns]
    y.mad.a <- y.mad.a[, non.na.columns]
    y.l.a <- y.l.a[, non.na.columns]
    y.u.a <- y.u.a[, non.na.columns]
    x <- x[non.na.columns]
    y <- y[non.na.columns]

    y.med <- apply(y.med.a, 2, function(x) mean(x, na.rm=T))
    y.mad <- apply(y.mad.a, 2, function(x) mean(x, na.rm=T))
    y.l <- apply(y.l.a, 2, function(x) mean(x, na.rm=T))
    y.u <- apply(y.u.a, 2, function(x) mean(x, na.rm=T))

    r2 <- gEff.pval(x, y.med, y.mad)
    r3 <- n.stats2(ph, ref.ph, idx1)

    gEff <- numeric(8)
    names(gEff) <- c("gEff", "p-val","q-val","n(TERM)","n(sPTB)", "n(subj TERM)","n(subj sPTB)","mode($\\log_{10}$(rat))")
    gEff[c(1,2,4:8)] <- c(r2$eff, r2$pval, r3$n.TERM, r3$n.sPTB,
                             r3$nCtrkSubjs, r3$nCaseSubjs, r3$mode.log.rat)

    gEff.min <- numeric(8)
    names(gEff.min) <- c("gEff", "p-val","q-val","n(TERM)","n(sPTB)", "n(subj TERM)","n(subj sPTB)","mode($\\log_{10}$(rat))")
    gEff.min[c(1,2,4:8)] <- c(r2$eff, r2$pval.min, r3$n.TERM, r3$n.sPTB,
                             r3$nCtrkSubjs, r3$nCaseSubjs, r3$mode.log.rat)

    file <- paste0(picsDir,ph,"_",ref.ph,".pdf")
    pdf(file, width=6, height=6)
    op <- par(mar=c(3.5, 3.5, 2.5, 0.5), mgp=c(2.25,0.6,0),tcl = -0.3)
    title <- paste0(ph," vs ",ref.ph)
    plot.logit(x, y, expit(y.med), expit(y.l), expit(y.u), title=title, ylab="pr( sPTB )", xlab="log10( absAb ratio )")
    par(op)
    dev.off()

    list(ph=ph,
         ref.ph=ref.ph,
         x=x, y=y,
         r2=r2,
         gEff=gEff,
         gEff.min=gEff.min,
         y.med=y.med,
         y.mad=y.mad,
         y.l=y.l,
         y.u=y.u,
         y.med.a=y.med.a,
         y.mad.a=y.mad.a,
         y.l.a=y.l.a,
         y.u.a=y.u.a,
         file=file)
}


spmrf.bernoulli.o2.no.coeffects.fn <- function(ph, nModels=5, nItr=5000, thin=1, nChains=1, nCores=1,
                                              alpha=0.05,
                                              picsDir="../pics/coeffects_dir/target_dir/")
{
    plow <- alpha/2
    phigh <- 1 - alpha/2

    j0 <- grep(ph, other.phs)
    ph.x <- pt2[,ph]
    ref.ph.x <- pt2[,ref.ph]
    y <- mt$delStatus.n
    idx <- idx1 & ph.x > 0 & ref.ph.x == 0
    y <- y[idx]
    ph.x <- log10(ph.x[idx])
    ph.sigma <- nexpFn(ph.x)
    x <- ph.x
    o <- order(x)
    x <- x[o]
    nx <- length(x)

    y.med <- list()
    y.mad <- list()
    y.l <- list()
    y.u <- list()
    m.o2 <- list()
    for ( i in seq(nModels) )
    {
        print(paste(ph, "vs", ref.ph, j0, i))
        ph.z <- numeric(nx)
        for ( j in seq(nx) )
        {
            ph.z[j] <- rnorm(1, mean=ph.x[j], sd=ph.sigma[j])
        }
        z <- ph.z
        o <- order(z)
        z <- z[o]
        sy <- y[o]
        m.dat <- list(y=sy, xvar1=z, J=length(y))
        m.dat2 <- spmrf.get.data(m.dat)
        m.o2[[i]] <- sampling(spmrf.bernoulli.o2.model, data=m.dat2, control=list(adapt_delta=0.96, max_treedepth=12), iter=nItr, thin=thin, chains=nChains, cores=nCores)
        theta <- rstan::extract(m.o2[[i]], "theta")[[1]]
        y.med[[i]] <- approx(z, apply(theta, 2, median), xout=x)
        y.mad[[i]] <- approx(z, apply(theta, 2, mad), xout=x)
        y.l[[i]] <- approx(z, apply(theta, 2, quantile, probs = plow), xout=x)
        y.u[[i]] <- approx(z, apply(theta, 2, quantile, probs = phigh), xout=x)
    }

    y.med.a <- c()
    y.mad.a <- c()
    y.l.a <- c()
    y.u.a <- c()
    for ( i in seq(nModels) )
    {
        y.med.a <- rbind(y.med.a, y.med[[i]]$y)
        y.mad.a <- rbind(y.mad.a, y.mad[[i]]$y)
        y.l.a <- rbind(y.l.a, y.l[[i]]$y)
        y.u.a <- rbind(y.u.a, y.u[[i]]$y)
    }

    non.na.columns <- apply(y.med.a, 2, function(x) !any(is.na(x)))
    y.med.a <- y.med.a[, non.na.columns]
    y.mad.a <- y.mad.a[, non.na.columns]
    y.l.a <- y.l.a[, non.na.columns]
    y.u.a <- y.u.a[, non.na.columns]
    x <- x[non.na.columns]
    y <- y[non.na.columns]

    y.med <- apply(y.med.a, 2, function(x) mean(x, na.rm=T))
    y.mad <- apply(y.mad.a, 2, function(x) mean(x, na.rm=T))
    y.l <- apply(y.l.a, 2, function(x) mean(x, na.rm=T))
    y.u <- apply(y.u.a, 2, function(x) mean(x, na.rm=T))

    r2 <- gEff.pval(x, y.med, y.mad)
    r3 <- n.stats2.no(ph, ref.ph, idx1)

    gEff <- numeric(8)
    names(gEff) <- c("gEff", "p-val","q-val","n(TERM)","n(sPTB)", "n(subj TERM)","n(subj sPTB)","mode($\\log_{10}$(rat))")
    gEff[c(1,2,4:8)] <- c(r2$eff, r2$pval, r3$n.TERM, r3$n.sPTB,
                             r3$nCtrkSubjs, r3$nCaseSubjs, r3$mode.log.rat)

    gEff.min <- numeric(8)
    names(gEff.min) <- c("gEff", "p-val","q-val","n(TERM)","n(sPTB)", "n(subj TERM)","n(subj sPTB)","mode($\\log_{10}$(rat))")
    gEff.min[c(1,2,4:8)] <- c(r2$eff, r2$pval.min, r3$n.TERM, r3$n.sPTB,
                             r3$nCtrkSubjs, r3$nCaseSubjs, r3$mode.log.rat)

    file <- paste0(picsDir,ph,"__no_",ref.ph,".pdf")
    pdf(file, width=6, height=6)
    op <- par(mar=c(3.5, 3.5, 2.5, 0.5), mgp=c(2.25,0.6,0),tcl = -0.3)
    title <- paste0(ph," vs ",ref.ph)
    plot.logit(x, y, expit(y.med), expit(y.l), expit(y.u), title=title, ylab="pr( sPTB )", xlab="log10( absAb )")
    par(op)
    dev.off()

    list(ph=ph,
         ref.ph=ref.ph,
         x=x, y=y,
         r2=r2,
         gEff=gEff,
         gEff.min=gEff.min,
         y.med=y.med,
         y.mad=y.mad,
         y.l=y.l,
         y.u=y.u,
         y.med.a=y.med.a,
         y.mad.a=y.mad.a,
         y.l.a=y.l.a,
         y.u.a=y.u.a,
         file=file)
}


## norm.2gr.diff.model <- stan_model(file="/Users/pgajer/projects/Statistics/MCMC_Models/norm_2gr_diff.stan")
## save(norm.2gr.diff.model, file="/Users/pgajer/projects/Statistics/MCMC_Models/norm_2gr_diff.rda")
load("/Users/pgajer/projects/Statistics/MCMC_Models/norm_2gr_diff.rda")

## nItr=5000; thin=5; nChains=3; nCores=3

dEff.pval <- function(y.min, se.min, y.max, se.max, eff, n=100,
                     nItr=5000, thin=5, nChains=3, nCores=3)
{
    y1 <- rnorm(n, mean=y.min, sd=se.min)
    y2 <- rnorm(n, mean=y.max, sd=se.max)

    m.dat <- list(y=c(y1,y2),
                 N=2*n,
                 gr=c(rep(1,n),rep(2,n)),
                 diff=diff)

    fit <- sampling(norm.2gr.diff.model, data=m.dat, iter=nItr, thin=thin, chains=nChains, cores=nCores)
    delta <- rstan::extract(fit, "delta")[[1]]

    delta.mean <- mean(delta)
    delta.l <- quantile(delta, probs = 0.025)
    ##delta.u <- apply(delta, 2, quantile, probs = phigh)
    sig <- delta.l > 0

    list(sig=sig, delta.l=delta.l)
}

gEff.pval <- function(x, y.med, y.mad)
{
    idx <- !is.na(y.med) & !is.na(y.mad)
    y.med <- y.med[idx]
    y.mad <- y.mad[idx]
    x <- x[idx]

    y <- y.med
    se <- y.mad

    i.min <- which.min(y)
    x.min <- x[i.min]
    y.min <- y[i.min]

    i.max <- which.max(y)
    x.max <- x[i.max]
    y.max <- y[i.max]

    pval.max <- pnorm(y[i.min], y[i.max], se[i.max])[[1]]
    pval.min <- 1 - pnorm(y[i.max], y[i.min], se[i.min])[[1]]

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
        eff <- expit(y.max) - expit(y.min)
    } else {
        eff <- -(expit(y.max) - expit(y.min))
    }

    list(pval=pval, pval.max=pval.max, pval.min=pval.min,
         eff=eff, x.min=x.min, y.min=y.min, x.max=x.max, y.max=y.max)
}


spmrf.logit.gEff <- function(fit.list, x, alpha=0.05)
{
    theta <- c()
    for ( i in seq(fit.list) )
    {
        new.theta <- rstan::extract(fit.list[[i]], "theta")[[1]]
        ##new.theta <- new.theta[,o.inv.list[[i]]]
        theta <- rbind(theta, new.theta)
    }

    plow <- alpha/2
    phigh <- 1 - alpha/2
    y <- apply(theta, 2, median)
    se <- apply(theta, 2, mad)

    i.min <- which.min(y)
    x.min <- x[i.min]
    y.min <- y[i.min]

    i.max <- which.max(y)
    x.max <- x[i.max]
    y.max <- y[i.max]

    pval.max <- pnorm(y[i.min], y[i.max], se[i.max])[[1]]
    pval.min <- 1 - pnorm(y[i.max], y[i.min], se[i.min])[[1]]

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
        eff <- expit(y.max) - expit(y.min)
    } else {
        eff <- -(expit(y.max) - expit(y.min))
    }

    list(pval=pval, pval.max=pval.max, pval.min=pval.min, theta=theta,
         theta.md=y, theta.mad=se,
         eff=eff, x.min=x.min, y.min=y.min, x.max=x.max, y.max=y.max)
}

spmrf.logit.spline <- function(fit)
{
    theta <- rstan::extract(fit, "theta")[[1]]
    theta.md <- expit(apply(theta, 2, median))
    theta.md
}

plot.spmrf.logit.list <- function(fit.list,
                                 x,
                                 ##o.inv.list,
                                 title="",
                                 alpha=0.05,
                                 xlab="x", ylab="pr( y=1 )",
                                 pt.pch=20, pt.cex=1,
                                 med.col="blue", bci.col="gray90")
{
    theta <- c()
    for ( i in seq(fit.list) )
    {
        new.theta <- rstan::extract(fit.list[[i]], "theta")[[1]]
        ##new.theta <- new.theta[,o.inv.list[[i]]]
        theta <- rbind(theta, new.theta)
    }

    plow <- alpha/2
    phigh <- 1 - alpha/2
    theta.md <- expit(apply(theta, 2, median))
    theta.l <- expit(apply(theta, 2, quantile, probs = plow))
    theta.u <- expit(apply(theta, 2, quantile, probs = phigh))

    plot(x, theta.md, type = "n", ylim=c(0,1), las=1,  xlab=xlab, ylab=ylab, main=title) # ylim=range(c(theta.l, theta.u))
    polygon(c(x, rev(x)), c(theta.l, rev(theta.u)),  border = NA, col=bci.col)
    ##points(x, y, pch=pt.pch, cex=pt.cex)
    lines(x, theta.md, col=med.col, lwd=1)
    rug(x[y==0])
    rug(x[y==1], col='red', lwd=2)
}

plot.logit <- function(x, y,
                      y.med, y.l, y.u,
                      ylim=c(0,1),
                      title="",
                      alpha=0.05,
                      xlab="x", ylab="pr( y=1 )",
                      pt.pch=20, pt.cex=1,
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
    rug(x[y==0])
    rug(x[y==1], col='red', lwd=2)
}

ylim=c(0,1)
title=""
alpha=0.05
xlab="x"
ylab="pr( y=1 )"
pt.pch=20
pt.cex=1
med.col="blue"
bci.col="gray90"

fit.and.plot.logit <- function(x, y,
                              ylim=c(0,1),
                              title="",
                              alpha=0.05,
                              xlab="x", ylab="pr( y=1 )",
                              pt.pch=20, pt.cex=1,
                              med.col="blue", bci.col="gray90")
{
    o <- order(x, decreasing = T)
    x <- x[o]
    y <- y[o]

    m <- glm(y~x, family=binomial)

    n <- length(x)
    t.val <- qt(0.975, n - 2) # Calculate critical t-value

    pred <- predict(m, type = "response", se.fit=TRUE)
    y.l <- pred$fit - t.val*pred$se.fit
    y.med <- pred$fit
    y.u <- pred$fit + t.val*pred$se.fit

    plot(x, y.med, type = "n", ylim=ylim, las=1,  xlab=xlab, ylab=ylab, main=title) # ylim=range(c(theta.l, theta.u))
    polygon(c(x, rev(x)), c(y.l, rev(y.u)),  border = NA, col=bci.col)
    lines(x, y.med, col=med.col, lwd=1)

    invisible(m)
}


plot.spmrf.logit <- function(fit, x, y,
                            title="",
                            alpha=0.05,
                            xlab="x", ylab="pr( y=1 )",
                            pt.pch=20, pt.cex=1,
                            med.col="blue", bci.col="gray90")
{
    theta <- rstan::extract(fit, "theta")[[1]]
    theta.md <- expit(apply(theta, 2, median))

    plow <- alpha/2
    phigh <- 1 - alpha/2
    theta.l <- expit(apply(theta, 2, quantile, probs = plow))
    theta.u <- expit(apply(theta, 2, quantile, probs = phigh))

    ## thp <- list(postmed = theta.md, bci.lower = theta.l, bci.upper = theta.u)

    x <- x[seq(length(theta.md))]

    plot(x, theta.md, type = "n", ylim=c(0,1), las=1,  xlab=xlab, ylab=ylab, main=title) # ylim=range(c(theta.l, theta.u))
    polygon(c(x, rev(x)), c(theta.l, rev(theta.u)),  border = NA, col=bci.col)
    ##points(x, y, pch=pt.pch, cex=pt.cex)
    lines(x, theta.md, col=med.col, lwd=1)
    rug(x[y==0])
    rug(x[y==1], col='red', lwd=2)
}

plot.spmrf.logit.vs.baseline <- function(fit, x, y,
                                        title="",
                                        alpha=0.05,
                                        xlab="x", ylab="pr( y=1 )",
                                        pt.pch=20, pt.cex=1,
                                        med.col="blue", bci.col="gray90")
{
    delta <- rstan::extract(fit, "delta")[[1]]
    delta.md <- apply(delta, 2, median)

    plow <- alpha/2
    phigh <- 1 - alpha/2
    delta.l <- apply(delta, 2, quantile, probs = plow)
    delta.u <- apply(delta, 2, quantile, probs = phigh)

    x <- x[seq(length(delta.md))]

    plot(x, delta.md, ylim=range(c(delta.l, delta.u)), type = "n", las=1,  xlab=xlab, ylab=ylab, main=title) # ylim=c(0,1) ylim=range(c(delta.l, delta.u))
    polygon(c(x, rev(x)), c(delta.l, rev(delta.u)),  border = NA, col=bci.col)
    ##points(x, y, pch=pt.pch, cex=pt.cex)
    lines(x, delta.md, col=med.col, lwd=1)
    abline(h=0, col='gray70', lwd=0.1)
    rug(x[y==0])
    rug(x[y==1], col='red', lwd=2)
}


Lb.cstColTbl <- c("darkorange", "royalblue")
names(Lb.cstColTbl) <- c("Lb", "IV")

plot.spmrf.logit.Lb.cst <- function(fit, x, y, cst,
                                title="",
                                alpha=0.05,
                                xlab="x", ylab="pr( y=1 )",
                                pt.pch=20, pt.cex=1,
                                med.col="blue", bci.col="gray90")
{
    theta <- rstan::extract(fit, "theta")[[1]]

    plow <- alpha/2
    phigh <- 1 - alpha/2
    theta.md <- expit(apply(theta, 2, median))
    theta.l <- expit(apply(theta, 2, quantile, probs = plow))
    theta.u <- expit(apply(theta, 2, quantile, probs = phigh))

    ## thp <- list(postmed = theta.md, bci.lower = theta.l, bci.upper = theta.u)

    x <- x[seq(length(theta.md))]
    cst <- cst[seq(length(theta.md))]
    cst <- ifelse(cst=="IV-A" | cst=="IV-B", "IV", "Lb")

    op <- par(mar=c(3.75, 3.75, 2.5, 3.5), mgp=c(2.5,0.5,0),tcl = -0.3)
    plot(x, theta.md, type = "n", ylim=c(0,1), las=1,  xlab=xlab, ylab=ylab, main=title) # ylim=range(c(theta.l, theta.u))
    polygon(c(x, rev(x)), c(theta.l, rev(theta.u)),  border = NA, col=bci.col)
    ##points(x, y, pch=pt.pch, cex=pt.cex)
    lines(x, theta.md, col=med.col, lwd=1)
    rug(x[y==0])
    rug(x[y==1], col='red', lwd=2)
    csts <- unique(cst)
    for ( one.cst in csts )
    {
        idx <- cst==one.cst
        rug(x[idx], side=3, col=Lb.cstColTbl[one.cst])
    }
    legend(par('usr')[1], par('usr')[4], xjust=-0.1, yjust=2.25 ,
           xpd=NA, cex=0.8, legend=names(Lb.cstColTbl), ncol=length(Lb.cstColTbl), fill=Lb.cstColTbl)
    par(op)
}


cstColTbl <- c("red1", "chartreuse", "darkorange", "aquamarine4", "royalblue", "yellow")
names(cstColTbl) <- c("I", "II", "III", "IV-A", "IV-B", "V")

plot.spmrf.logit.cst <- function(fit, x, y, cst,
                                title="",
                                alpha=0.05,
                                xlab="x", ylab="pr( y=1 )",
                                pt.pch=20, pt.cex=1,
                                med.col="blue", bci.col="gray90")
{
    theta <- rstan::extract(fit, "theta")[[1]]

    plow <- alpha/2
    phigh <- 1 - alpha/2
    theta.md <- expit(apply(theta, 2, median))
    theta.l <- expit(apply(theta, 2, quantile, probs = plow))
    theta.u <- expit(apply(theta, 2, quantile, probs = phigh))

    ## thp <- list(postmed = theta.md, bci.lower = theta.l, bci.upper = theta.u)

    x <- x[seq(length(theta.md))]
    cst <- cst[seq(length(theta.md))]

    op <- par(mar=c(3.75, 3.75, 2.5, 3.5), mgp=c(2.5,0.5,0),tcl = -0.3)
    plot(x, theta.md, type = "n", ylim=c(0,1), las=1,  xlab=xlab, ylab=ylab, main=title) # ylim=range(c(theta.l, theta.u))
    polygon(c(x, rev(x)), c(theta.l, rev(theta.u)),  border = NA, col=bci.col)
    ##points(x, y, pch=pt.pch, cex=pt.cex)
    lines(x, theta.md, col=med.col, lwd=1)
    rug(x[y==0])
    rug(x[y==1], col='red', lwd=2)
    csts <- unique(cst)
    for ( one.cst in csts )
    {
        idx <- cst==one.cst
        rug(x[idx], side=3, col=cstColTbl[one.cst])
    }
    legend(par('usr')[1], par('usr')[4], xjust=-0.1, yjust=2.25 , xpd=NA, cex=0.8, legend=names(cstColTbl), ncol=length(cstColTbl), fill=cstColTbl)  #title="Community State Type", cex=1
    par(op)
}


b.spline <- function(fit)
{
    ff <- extract(fit)
    mu.med <- array(NA, length(y))
    for (i in 1:length(y))
    {
        mu.med[i] <- expit(median(ff$mu[,i]))
    }
    mu.med
}

plot.b.spline <- function(fit, x, y, med.col="blue", bci.col="gray90", title="", ... )
{
    ff <- extract(fit)
    mu.med <- array(NA, length(y))
    mu.ub <- array(NA, length(y))
    mu.lb <- array(NA, length(y))
    for (i in 1:length(y)) {
        mu.med[i] <- expit(median(ff$mu[,i]))
        mu.lb[i] <- expit(quantile(ff$mu[,i],probs = 0.025))
        mu.ub[i] <- expit(quantile(ff$mu[,i],probs = 0.975))
    }

    plot(x, y, las=1, main=title, ...) # col="azure4"
    polygon(c(rev(x), x), c(rev(mu.lb), mu.ub), col=bci.col, border = NA)
    lines(x, mu.med, col=med.col, lty=1)
}

b.spline.gEff <- function(fit, x, y)
{
    ff <- extract(fit)
    mu.med <- array(NA, length(y))
    mu.ub <- array(NA, length(y))
    mu.lb <- array(NA, length(y))
    se <- array(NA, length(y))
    for (i in 1:length(y)) {
        mu.med[i] <- median(ff$mu[,i])
        mu.lb[i] <- quantile(ff$mu[,i],probs = 0.025)
        mu.ub[i] <- quantile(ff$mu[,i],probs = 0.975)
        se[i] <- mad(ff$mu[,i])
    }

    y <- mu.med
    i.min <- which.min(y)
    x.min <- x[i.min]
    y.min <- y[i.min]

    i.max <- which.max(y)
    x.max <- x[i.max]
    y.max <- y[i.max]

    pval.max <- pnorm(y[i.min], y[i.max], se[i.max])[[1]]
    pval.min <- 1 - pnorm(y[i.max], y[i.min], se[i.min])[[1]]

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
        eff <- expit(y.max) - expit(y.min)
    } else {
        eff <- -(expit(y.max) - expit(y.min))
    }

    list(pval=pval, pval.max=pval.max, pval.min=pval.min,
         eff=eff, x.min=x.min, y.min=y.min, x.max=x.max, y.max=y.max)
}

## compute p-value of a coefficient of m.fit
## see ~/projects/louisville_cdc/R/models.R
## for examples
get.stan.pval <- function(fit, coefName)
{
    mu <- as.vector(extract(fit, pars=coefName, permuted = TRUE)[[1]])

    ## qqnorm(mu)
    ## qqline(mu, col=2)

    mu.bar <- median(mu)

    mu.pFn <- function(p)
    {
        ret <- 0
        p2 <- p/2
        mu.q <- quantile(mu, probs=c(p2,1-p2))
        if ( mu.bar > 0 )
        {
            ret <- mu.q[1]
        }
        else
        {
            ret <- mu.q[2]
        }
        ret
    }

    mu.pval.Fn <- NA
    if ( sign(mu.pFn(0)*mu.pFn(1))==-1 )
    {
        mu.pval.Fn <- uniroot(mu.pFn, c(0, 1))$root
    }
    ##mu.pval.Fn

    mu.pval <- NA
    if (  mu.bar < 0 )
    {
        mu.pval <- 1 - pnorm(0, mu.bar, mad(mu))
    } else {
        mu.pval <- pnorm(0, mu.bar, mad(mu))
    }
    ##mu.pval

    mu.pval2 <- NA
    if (  mu.bar < 0 )
    {
        mu.pval2 <- 1 - pnorm(0, mu.bar, sd(mu))
    } else {
        mu.pval2 <- pnorm(0, mu.bar, sd(mu))
    }

    if ( !is.na(mu.pval.Fn) ) #& mu.pval.Fn > mu.pval )
    {
        pval <- mu.pval.Fn
    }
    else
    {
        pval <- mu.pval
    }

    list(pval=pval, mu.pval=mu.pval, mu.pval2=mu.pval2, mu.pval.Fn=mu.pval.Fn, mu=mu, mu.bar=mu.bar)
}



bbp2me.gEff.tbl <- function(res, selPhs, q.thld=0.05, e.thld=0.1)
{
    res.phs <- c()
    for ( i in seq(res) )
    {
        r <- res[[i]]
        res.phs <- c(res.phs, r$ph)
    }

    gEff <- matrix(nrow=length(selPhs), ncol=8)
    colnames(gEff) <- c("gEff", "p-val","q-val","n(TERM)","n(sPTB)", "n(subj TERM)","n(subj sPTB)","-med($\\log_{10}$(RA))")
    rownames(gEff) <- selPhs
    for ( ph in selPhs )
    {
        i <- grep(ph, res.phs)
        r <- res[[i]]
        ## ph <- r$ph
        f <- stan.bbp2me.eff.pval(r$fit)
        gEff[ph, c(1,2,4:8)] <- c(f$eff, f$pval, (r$n - r$n.sPTB), r$n.sPTB,
                                         r$nCtrkSubjs, r$nCaseSubjs, r$med.x)
    }

    gEff[,3] <- p.adjust(gEff[,2], method="fdr")

    idx <- gEff[,3] < q.thld & abs(gEff[,1]) > e.thld
    gEff.sig <- as.data.frame(gEff)[idx,]

    o <- order(gEff[,2])
    gEff <- gEff[o,]

    list(gEff=gEff, gEff.sig=gEff.sig)
}

bbp2me.gEff.tbl2 <- function(res, selPhs, idx2, q.thld=0.05, e.thld=0.1)
{
    res.phs <- c()
    for ( i in seq(res) )
    {
        r <- res[[i]]
        res.phs <- c(res.phs, r$ph)
    }

    selPhs.res.phs.diff <- setdiff(selPhs, res.phs)
    res.phs.selPhs.diff <- setdiff(res.phs, selPhs)
    selPhs <- intersect(selPhs, res.phs)

    gEff <- matrix(nrow=length(selPhs), ncol=8)
    colnames(gEff) <- c("gEff", "p-val","q-val","n(TERM)","n(sPTB)", "n(subj TERM)","n(subj sPTB)","-med($\\log_{10}$(RA))")
    rownames(gEff) <- selPhs
    for ( ph in selPhs )
    {
        i <- grep(ph, res.phs)
        if ( length(i)>0 )
        {
            r <- res[[i]]
            ph <- r$ph
            rr <- n.stats(ph, idx1)
            f <- stan.bbp2me.eff.pval(r$fit)
            gEff[ph, c(1,2,4:8)] <- c(f$eff, f$pval,
                                     rr$n.TERM, rr$n.sPTB,
                                     rr$nCtrkSubjs, rr$nCaseSubjs, rr$med.x)
        }
    }
    gEff[,3] <- p.adjust(gEff[,2], method="fdr")

    idx <- gEff[,3] < q.thld & abs(gEff[,1]) > e.thld
    gEff.sig <- as.data.frame(gEff)[idx,]

    o <- order(gEff[,2])
    gEff <- gEff[o,]

    list(gEff=gEff, gEff.sig=gEff.sig,
         selPhs.res.phs.diff=selPhs.res.phs.diff,
         res.phs.selPhs.diff=res.phs.selPhs.diff)
}


bbp2me.eff.tbl <- function(res, selPh, q.thld=0.1)
{
    res.phs <- c()
    for ( i in seq(res) )
    {
        r <- res[[i]]
        res.phs <- c(res.phs, r$ph)
    }

    mat <- matrix(nrow=length(selPh), ncol=11)
    colnames(mat) <- c("lin", "lin p-val","lin q-val","quad", "quad p-val","quad q-val","n(TERM)","n(sPTB)","n(subj TERM)", "n(subj sPTB)","-med($\\log_{10}$(RA))")
    rownames(mat) <- selPh
    for ( ph in selPh )
    {
        i <- grep(ph, res.phs)
        r <- res[[i]]
        fit <- r$fit
        m <- r$m
        mat[ph, c(1,2,4,5,7:11)] <- c(m["beta[1]",1],
                                     get.stan.pval(fit, "beta[1]")$mu.pval,
                                     m["beta[2]",1],
                                     get.stan.pval(fit, "beta[2]")$mu.pval,
                                     (r$n - r$n.sPTB), r$n.sPTB,
                                     r$nCtrkSubjs, r$nCaseSubjs, r$med.x)
    }

    mat[,3] <- p.adjust(mat[,2], method="fdr")
    mat[,6] <- p.adjust(mat[,5], method="fdr")

    o <- order(mat[,2])
    mat <- mat[o,]

    idx <- mat[,3] < q.thld | mat[,6] < q.thld
    (mat.sig <- as.data.frame(mat)[idx,])

    list(mat=mat, mat.sig=mat.sig)
}


bbp2me.eff.tbl2 <- function(res, selPh, idx1, q.thld=0.1)
{
    res.phs <- c()
    for ( i in seq(res) )
    {
        r <- res[[i]]
        res.phs <- c(res.phs, r$ph)
    }

    selPhs.res.phs.diff <- setdiff(selPh, res.phs)
    res.phs.selPhs.diff <- setdiff(res.phs, selPh)
    selPh <- intersect(selPh, res.phs)

    mat <- matrix(nrow=length(selPh), ncol=11)
    colnames(mat) <- c("lin", "lin p-val","lin q-val","quad", "quad p-val","quad q-val","n(TERM)","n(sPTB)","n(subj TERM)", "n(subj sPTB)","-med($\\log_{10}$(RA))")
    rownames(mat) <- selPh
    for ( ph in selPh )
    {
        i <- grep(ph, res.phs)
        if ( length(i)>0 )
        {
            r <- res[[i]]
            ph <- r$ph
            rr <- n.stats(ph, idx1)
            fit <- r$fit
            m <- r$m
            mat[ph, c(1,2,4,5,7:11)] <- c(m["beta[1]",1],
                                         get.stan.pval(fit, "beta[1]")$mu.pval,
                                         m["beta[2]",1],
                                         get.stan.pval(fit, "beta[2]")$mu.pval,
                                         rr$n.TERM, rr$n.sPTB,
                                         rr$nCtrkSubjs, rr$nCaseSubjs, rr$med.x)
        }
    }
    mat[,3] <- p.adjust(mat[,2], method="fdr")
    mat[,6] <- p.adjust(mat[,5], method="fdr")

    o <- order(mat[,2])
    mat <- mat[o,]

    idx <- mat[,3] < q.thld | mat[,6] < q.thld
    (mat.sig <- as.data.frame(mat)[idx,])

    list(mat=mat, mat.sig=mat.sig,
         selPhs.res.phs.diff=selPhs.res.phs.diff,
         res.phs.selPhs.diff=res.phs.selPhs.diff)
}


## computing maximum change of sPTB over the whole range of log10 Relative
## Abundances of a given phylotype in the bbp2me model
stan.bbp2me.eff.pval <- function(fit)
{
    m.theta <- summary(fit, pars=c("theta"), probs=c(0.025,0.975))$summary
    m.tlp <- summary(fit, pars=c("tlp"), probs=c(0.025,0.975))$summary

    xx <- m.tlp[,1]
    yy <- m.theta[,1]
    o <- order(xx)
    xx <- xx[o]
    yy <- yy[o]
    lm <- loess(yy ~ xx)
    ## plot(xx,yy)
    ## lines(xx,predict(lm), col=2, lwd=2)

    x.fit <- xx
    y.fit <- expit(predict(lm))

    x <- x.fit
    y <- y.fit

    i.min <- which.min(y)
    x.min <- x[i.min]
    y.min <- y[i.min]

    i.max <- which.max(y)
    x.max <- x[i.max]
    y.max <- y[i.max]

    theta <- extract(fit, pars=c("theta"), permuted = TRUE)$theta
    theta <- theta[,o]

    theta.max <- theta[,i.max]
    ## myHist(theta.max)
    ## library(car) # for qqPlot()
    ## qqPlot(theta.max)
    pval.max <- pnorm(yy[i.min], yy[i.max], mad(theta.max))

    theta.min <- theta[,i.min]
    ## myHist(theta.min)
    ## library(car) # for qqPlot()
    ## qqPlot(theta.min)
    pval.min <- 1 - pnorm(yy[i.max], yy[i.min], mad(theta.min))

    if ( pval.max > pval.min )
    {
        pval <- pval.max
    }
    else
    {
        pval <- pval.min
    }

    n <- length(y)
    eff <- NA

    ## percentage of points over which first derivative of y.fit is positive
    ## (y.fit is increasing)
    dy.fit <- diff(y.fit)
    perc.incr <- 100*sum(dy.fit>0)/n
    incr.thld <- 60

    ##if ( (i.min==1 && i.max==n) || (i.min!=1 && i.min!=n) ) # risk is increasing or U-shaped
    if ( (perc.incr > incr.thld) || (i.min!=1 && i.min!=n) ) # risk is increasing or U-shaped
    {
        eff <- y.max - y.min
    } else {
        eff <- -(y.max - y.min)
    }

    list(pval=pval, eff=eff, x.min=x.min, y.min=y.min, x.max=x.max, y.max=y.max )
}

## generate plot of sPTB risk
plot.bbp2me <- function(fit, ph)
{
    m.theta <- summary(fit, pars=c("theta"), probs=c(0.025,0.975))$summary
    m.tlp <- summary(fit, pars=c("tlp"), probs=c(0.025,0.975))$summary

    xx <- m.tlp[,1]
    yy <- m.theta[,4]
    o <- order(xx)
    xx <- xx[o]
    yy <- yy[o]
    ll <- loess(yy ~ xx)
    ## plot(xx,yy)
    ## lines(xx,predict(ll), col=2, lwd=2)

    xx <- m.tlp[,1]
    yy <- m.theta[,5]
    o <- order(xx)
    xx <- xx[o]
    yy <- yy[o]
    lu <- loess(yy ~ xx)
    ## plot(xx,yy)
    ## lines(xx,predict(lu), col=2, lwd=2)

    xx <- m.tlp[,1]
    yy <- m.theta[,1]
    o <- order(xx)
    xx <- xx[o]
    yy <- yy[o]
    lm <- loess(yy ~ xx)
    ## plot(xx,yy)
    ## lines(xx,predict(lm), col=2, lwd=2)

    x.fit <- xx
    y.fit <- expit(predict(lm))
    ci.lower <- expit(predict(ll))
    ci.upper <- expit(predict(lu))

    x <- x.fit
    y <- y.fit

    i.min <- which.min(y)
    x.min <- x[i.min]
    y.min <- y[i.min]

    i.max <- which.max(y)
    x.max <- x[i.max]
    y.max <- y[i.max]

    h.y <- NA
    n <- length(y)
    ## thld <- 10
    ## if ( (i.min < thld && i.max > (n-thld)) || (i.min > thld && i.min < (n-thld)) ) # risk is increasing or U-shaped
    ## {
    ##     h.y <- y.min
    ## }
    ## else
    ## {
    ##     h.y <- y.max
    ## }

    ## h.y <- y.min

    ## percentage of points over which first derivative of y.fit is positive
    ## (y.fit is increasing)
    dy.fit <- diff(y.fit)
    perc.incr <- 100*sum(dy.fit>0)/n
    incr.thld <- 60

    incr.or.U.shape <- (perc.incr > incr.thld) || (i.min!=1 && i.min!=n)
    if ( incr.or.U.shape ) # risk is increasing or U-shaped
    {
        h.y <- y.min
    } else {
        h.y <- y.max
    }

    plot(x.fit, y.fit, las=1, ylim=c(0,1), type='n', xlab="log10( Relative Abundance )", ylab="pr( sPTB )", main=ph)
    polygon( c(x.fit, rev(x.fit)), c(ci.lower, rev(ci.upper)), col='gray90', border=NA )
    lines(x.fit, y.fit, col=2, lwd=1)
    abline(h=h.y, col="gray80")

    invisible(list(h.y=h.y, ci.lower=ci.lower, ci.upper=ci.upper,
                   x.min=x.min, x.max=x.max,
                   incr.or.U.shape=incr.or.U.shape,
                   x.fit=x.fit, y.fit=y.fit))
}



## generate plot of sPTB risk vs log10 absAb's
plot.aa.bbp2me <- function(fit, lbl, ph)
{
    m.theta <- summary(fit, pars=c("theta"), probs=c(0.025,0.975))$summary
    m.tlp <- summary(fit, pars=c("tlp"), probs=c(0.025,0.975))$summary

    xx <- m.tlp[,1] + lbl
    yy <- m.theta[,4]
    o <- order(xx)
    xx <- xx[o]
    yy <- yy[o]
    ll <- loess(yy ~ xx)
    ## plot(xx,yy)
    ## lines(xx,predict(ll), col=2, lwd=2)

    xx <- m.tlp[,1] + lbl
    yy <- m.theta[,5]
    o <- order(xx)
    xx <- xx[o]
    yy <- yy[o]
    lu <- loess(yy ~ xx)
    ## plot(xx,yy)
    ## lines(xx,predict(lu), col=2, lwd=2)

    xx <- m.tlp[,1] + lbl
    yy <- m.theta[,1]
    o <- order(xx)
    xx <- xx[o]
    yy <- yy[o]
    lm <- loess(yy ~ xx)
    ## plot(xx,yy)
    ## lines(xx,predict(lm), col=2, lwd=2)

    x.fit <- xx
    y.fit <- expit(predict(lm))
    ci.lower <- expit(predict(ll))
    ci.upper <- expit(predict(lu))

    x <- x.fit
    y <- y.fit

    i.min <- which.min(y)
    x.min <- x[i.min]
    y.min <- y[i.min]

    i.max <- which.max(y)
    x.max <- x[i.max]
    y.max <- y[i.max]

    h.y <- NA
    n <- length(y)
    ## thld <- 10
    ## if ( (i.min < thld && i.max > (n-thld)) || (i.min > thld && i.min < (n-thld)) ) # risk is increasing or U-shaped
    ## {
    ##     h.y <- y.min
    ## }
    ## else
    ## {
    ##     h.y <- y.max
    ## }

    ## h.y <- y.min

    ## percentage of points over which first derivative of y.fit is positive
    ## (y.fit is increasing)
    dy.fit <- diff(y.fit)
    perc.incr <- 100*sum(dy.fit>0)/n
    incr.thld <- 60

    incr.or.U.shape <- (perc.incr > incr.thld) || (i.min!=1 && i.min!=n)
    if ( incr.or.U.shape ) # risk is increasing or U-shaped
    {
        h.y <- y.min
    } else {
        h.y <- y.max
    }

    plot(x.fit, y.fit, las=1, ylim=c(0,1), type='n', xlab="log10( Absolute Abundance )", ylab="pr( sPTB )", main=ph)
    polygon( c(x.fit, rev(x.fit)), c(ci.lower, rev(ci.upper)), col='gray90', border=NA )
    lines(x.fit, y.fit, col=2, lwd=1)
    abline(h=h.y, col="gray80")

    invisible(list(h.y=h.y, ci.lower=ci.lower, ci.upper=ci.upper,
                   x.min=x.min, x.max=x.max,
                   incr.or.U.shape=incr.or.U.shape,
                   x.fit=x.fit, y.fit=y.fit))
}


## check_div
## check_treedepth
## check_energy
## partition_div

## are from https://github.com/betanalpha/knitr_case_studies/tree/master/rstan_workflow

check_div <- function(fit)
{
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
  divergent <- do.call(rbind, sampler_params)[,'divergent__']
  n = sum(divergent)
  ##N = length(divergent)
  ##c(n, N, 100 * n / N)
  n
}

#' Check transitions that ended with a divergence
orig_check_div <- function(fit) {
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
  divergent <- do.call(rbind, sampler_params)[,'divergent__']
  n = sum(divergent)
  N = length(divergent)

  print(sprintf('%s of %s iterations ended with a divergence (%s%%)',
                n, N, 100 * n / N))
  if (n > 0)
    print('Try running with larger adapt_delta to remove the divergences')
}


# Check transitions that ended prematurely due to maximum tree depth limit
orig_check_treedepth <- function(fit, max_depth = 10) {
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
  treedepths <- do.call(rbind, sampler_params)[,'treedepth__']
  n = length(treedepths[sapply(treedepths, function(x) x == max_depth)])
  N = length(treedepths)

  print(sprintf('%s of %s iterations saturated the maximum tree depth of %s (%s%%)',
                n, N, max_depth, 100 * n / N))
  if (n > 0)
    print('Run again with max_depth set to a larger value to avoid saturation')
}

check_treedepth <- function(fit, max_depth = 10)
{
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
  treedepths <- do.call(rbind, sampler_params)[,'treedepth__']
  n = length(treedepths[sapply(treedepths, function(x) x == max_depth)])
  n
}

# Checks the energy Bayesian fraction of missing information (E-BFMI)
check_energy <- function(fit) {
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
  for (n in 1:length(sampler_params)) {
    energies = sampler_params[n][[1]][,'energy__']
    numer = sum(diff(energies)**2) / length(energies)
    denom = var(energies)
    if (numer / denom < 0.2) {
        print(sprintf('Chain %s: E-BFMI = %s', n, numer / denom))
        print('E-BFMI below 0.2 indicates you may need to reparameterize your model')
    }
  }
}

# Returns parameter arrays separated into divergent and non-divergent transitions
partition_div <- function(fit) {
  nom_params <- extract(fit, permuted=FALSE)
  n_chains <- dim(nom_params)[2]
  params <- as.data.frame(do.call(rbind, lapply(1:n_chains, function(n) nom_params[,n,])))

  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
  divergent <- do.call(rbind, sampler_params)[,'divergent__']
  params$divergent <- divergent

  div_params <- params[params$divergent == 1,]
  nondiv_params <- params[params$divergent == 0,]

  return(list(div_params, nondiv_params))
}


##
##
##
spmrf.xy.fn <- function(x, y, nItr=2000, thin=1, nChains=3, nCores=3, alpha=0.05)
{
    m.dat <- list(y=y, xvar1=x, J=length(y))
    m.dat2 <- spmrf.get.data(m.dat)

    r <- sampling(spmrf.bernoulli.o2.model, data=m.dat2, iter=nItr, thin=thin, chains=nChains, cores=nCores)

    theta <- rstan::extract(r, "theta")[[1]]
    plow <- alpha/2
    phigh <- 1 - alpha/2
    y.med <- apply(theta, 2, median)
    y.mad <- apply(theta, 2, mad)
    y.l <- apply(theta, 2, quantile, probs = plow)
    y.u <- apply(theta, 2, quantile, probs = phigh)

    idx <- seq(y.med)
    r2 <- gEff.pval(x[idx], y.med, y.mad)

    list(x=x,
         y=y,
         r2=r2,
         y.med=y.med,
         y.mad=y.mad,
         y.l=y.l,
         y.u=y.u)
}

spmrf.ph2.T.csts.AA.fn <- function(T, phs, ph="Mobiluncus_curtisii_mulieris", nItr=2000, thin=1, nChains=3, nCores=3)
{
    x <- log10(pt.AA[, ph])

    if ( length(phs)==1 ){
        xT <- log10(pt.AA[, phs])
    } else {
        xT <- rowSums(ct.AA[, phs])
        n <- rowSums(ct.AA)
        xT <- log10(xT/n)
    }

    idx <- is.finite(x) & is.finite(xT)
    xT <- xT[idx]
    x <- x[idx]
    mt.AA.R <- mt.AA[idx,]

    qT <- quantile(xT, probs = c(0.3333,0.6666))

    tT <- xT
    tT[xT<=qT[1]] <- "T1"
    tT[xT<=qT[2] & xT>qT[1]] <- "T2"
    tT[xT>qT[2]] <- "T3"

    idx <- tT==T
    mt.AA.R <- mt.AA.R[idx,]
    x <- x[idx]
    o <- order(x)
    mt.AA.R <- mt.AA.R[o,]
    x <- x[o]

    cst.med.mat <- matrix(nrow=sum(idx), ncol=length(CSTs))
    colnames(cst.med.mat) <- CSTs

    for ( cst in CSTs )
    {
        print(cst)
        y <- ifelse(mt.AA.R$cst==cst, 1, 0)
        r <- spmrf.xy.fn(x, y, nItr=nItr, thin=thin, nChains=nChains, nCores=nCores)
        lo <- loess(r$y.med~r$x)
        p <- predict(lo)
        cst.med.mat[seq(p),cst] <- p
    }

    cst.med.mat <- cst.med.mat[!is.na(cst.med.mat[,1]),]

    rs <- rowSums(cst.med.mat)
    ##plot(rs)
    cst.med.mat <- cst.med.mat/rs
    ## rs <- rowSums(cst.med.mat)
    ## plot(rs)

    list(cst.med.mat=cst.med.mat, x=r$x)
}

## T="T1"; phs=c("Lactobacillus_iners"); ph="Mobiluncus_curtisii_mulieris"

spmrf.ph.T.AA.fn <- function(T, phs, ph="Mobiluncus_curtisii_mulieris", nItr=2000, thin=1, nChains=3, nCores=3)
{
    x <- log10(pt.AA[, ph])
    y <- mt.AA$sPTB
    s <- mt.AA$subjID

    if ( length(phs)==1 ){
        xT <- log10(pt.AA[, phs])
    } else {
        xT <- rowSums(ct.AA[, phs])
        n <- rowSums(ct.AA)
        xT <- log10(xT/n)
    }

    idx <- is.finite(x) & is.finite(xT)
    y <- y[idx]
    x <- x[idx]
    s <- s[idx]
    xT <- xT[idx]

    qT <- quantile(xT, probs = c(0.3333,0.6666))

    tT <- xT
    tT[xT<=qT[1]] <- "T1"
    tT[xT<=qT[2] & xT>qT[1]] <- "T2"
    tT[xT>qT[2]] <- "T3"

    idx <- tT==T
    y <- y[idx]
    x <- x[idx]
    s <- s[idx]
    o <- order(x)
    y <- y[o]
    x <- x[o]
    s <- s[o]

    m.dat <- list(y=y, xvar1=x, J=length(y))
    m.dat2 <- spmrf.get.data(m.dat)
    m <- sampling(spmrf.bernoulli.o2.model, data=m.dat2, control=list(adapt_delta=0.96, max_treedepth=12), iter=nItr, thin=thin, chains=nChains, cores=nCores)

    list(m=m, x=x, y=y, s=s, ph=ph, phs=phs, T=T)
}




spmrf.ph2.T.csts.fn <- function(T, phs, ph="Mobiluncus_curtisii_mulieris", nItr=2000, thin=1, nChains=3, nCores=3)
{
    x <- log10(pt.no.mPTB[, ph])

    if ( length(phs)==1 ){
        xT <- log10(pt.no.mPTB[, phs])
    } else {
        xT <- rowSums(ct.no.mPTB[, phs])
        n <- rowSums(ct.no.mPTB)
        xT <- log10(xT/n)
    }

    idx <- is.finite(x) & is.finite(xT)
    xT <- xT[idx]
    x <- x[idx]
    mt.no.mPTB.R <- mt.no.mPTB[idx,]

    qT <- quantile(xT, probs = c(0.3333,0.6666))

    tT <- xT
    tT[xT<=qT[1]] <- "T1"
    tT[xT<=qT[2] & xT>qT[1]] <- "T2"
    tT[xT>qT[2]] <- "T3"

    idx <- tT==T
    mt.no.mPTB.R <- mt.no.mPTB.R[idx,]
    x <- x[idx]
    o <- order(x)
    mt.no.mPTB.R <- mt.no.mPTB.R[o,]
    x <- x[o]

    cst.med.mat <- matrix(nrow=sum(idx), ncol=length(CSTs))
    colnames(cst.med.mat) <- CSTs

    for ( cst in CSTs )
    {
        print(cst)
        y <- ifelse(mt.no.mPTB.R$cst==cst, 1, 0)
        r <- spmrf.xy.fn(x, y, nItr=nItr, thin=thin, nChains=nChains, nCores=nCores)
        lo <- loess(r$y.med~r$x)
        p <- predict(lo)
        cst.med.mat[seq(p),cst] <- p
    }

    cst.med.mat <- cst.med.mat[!is.na(cst.med.mat[,1]),]

    rs <- rowSums(cst.med.mat)
    ##plot(rs)
    cst.med.mat <- cst.med.mat/rs
    ## rs <- rowSums(cst.med.mat)
    ## plot(rs)

    list(cst.med.mat=cst.med.mat, x=r$x)
}

## T="T1"; phs=c("Lactobacillus_iners"); ph="Mobiluncus_curtisii_mulieris"

spmrf.ph.T.fn <- function(T, phs, ph="Mobiluncus_curtisii_mulieris", nItr=2000, thin=1, nChains=3, nCores=3)
{
    x <- log10(pt.no.mPTB[, ph])
    y <- mt.no.mPTB$sPTB

    if ( length(phs)==1 ){
        xT <- log10(pt.no.mPTB[, phs])
    } else {
        xT <- rowSums(ct.no.mPTB[, phs])
        n <- rowSums(ct.no.mPTB)
        xT <- log10(xT/n)
    }

    idx <- is.finite(x) & is.finite(xT)
    y <- y[idx]
    x <- x[idx]
    xT <- xT[idx]

    qT <- quantile(xT, probs = c(0.3333,0.6666))

    tT <- xT
    tT[xT<=qT[1]] <- "T1"
    tT[xT<=qT[2] & xT>qT[1]] <- "T2"
    tT[xT>qT[2]] <- "T3"

    idx <- tT==T
    y <- y[idx]
    x <- x[idx]
    o <- order(x)
    y <- y[o]
    x <- x[o]

    m.dat <- list(y=y, xvar1=x, J=length(y))
    m.dat2 <- spmrf.get.data(m.dat)
    m <- sampling(spmrf.bernoulli.o2.model, data=m.dat2, control=list(adapt_delta=0.96, max_treedepth=12), iter=nItr, thin=thin, chains=nChains, cores=nCores)

    list(m=m, x=x, y=y, ph=ph, phs=phs, T=T)
}



spmrf.ph.T.AA.fn <- function(T, phs, ph="Mobiluncus_curtisii_mulieris", nItr=2000, thin=1, nChains=3, nCores=3)
{
    x <- log10(pt.AA[, ph])
    y <- mt.AA$sPTB

    if ( length(phs)==1 ){
        xT <- log10(pt.AA[, phs])
    } else {
        xT <- rowSums(ct.AA[, phs])
        n <- rowSums(ct.AA)
        xT <- log10(xT/n)
    }

    idx <- is.finite(x) & is.finite(xT)
    y <- y[idx]
    x <- x[idx]
    xT <- xT[idx]

    qT <- quantile(xT, probs = c(0.3333,0.6666))

    tT <- xT
    tT[xT<=qT[1]] <- "T1"
    tT[xT<=qT[2] & xT>qT[1]] <- "T2"
    tT[xT>qT[2]] <- "T3"

    idx <- tT==T
    y <- y[idx]
    x <- x[idx]
    o <- order(x)
    y <- y[o]
    x <- x[o]

    m.dat <- list(y=y, xvar1=x, J=length(y))
    m.dat2 <- spmrf.get.data(m.dat)
    m <- sampling(spmrf.bernoulli.o2.model, data=m.dat2, control=list(adapt_delta=0.96, max_treedepth=12), iter=nItr, thin=thin, chains=nChains, cores=nCores)

    list(m=m, x=x, y=y, ph=ph, phs=phs, T=T)
}




spmrf.ph2.T.csts.fn <- function(T, phs, ph="Mobiluncus_curtisii_mulieris", nItr=2000, thin=1, nChains=3, nCores=3)
{
    x <- log10(pt.no.mPTB[, ph])

    if ( length(phs)==1 ){
        xT <- log10(pt.no.mPTB[, phs])
    } else {
        xT <- rowSums(ct.no.mPTB[, phs])
        n <- rowSums(ct.no.mPTB)
        xT <- log10(xT/n)
    }

    idx <- is.finite(x) & is.finite(xT)
    xT <- xT[idx]
    x <- x[idx]
    mt.no.mPTB.R <- mt.no.mPTB[idx,]

    qT <- quantile(xT, probs = c(0.3333,0.6666))

    tT <- xT
    tT[xT<=qT[1]] <- "T1"
    tT[xT<=qT[2] & xT>qT[1]] <- "T2"
    tT[xT>qT[2]] <- "T3"

    idx <- tT==T
    mt.no.mPTB.R <- mt.no.mPTB.R[idx,]
    x <- x[idx]
    o <- order(x)
    mt.no.mPTB.R <- mt.no.mPTB.R[o,]
    x <- x[o]

    cst.med.mat <- matrix(nrow=sum(idx), ncol=length(CSTs))
    colnames(cst.med.mat) <- CSTs

    for ( cst in CSTs )
    {
        print(cst)
        y <- ifelse(mt.no.mPTB.R$cst==cst, 1, 0)
        r <- spmrf.xy.fn(x, y, nItr=nItr, thin=thin, nChains=nChains, nCores=nCores)
        lo <- loess(r$y.med~r$x)
        p <- predict(lo)
        cst.med.mat[seq(p),cst] <- p
    }

    cst.med.mat <- cst.med.mat[!is.na(cst.med.mat[,1]),]

    rs <- rowSums(cst.med.mat)
    ##plot(rs)
    cst.med.mat <- cst.med.mat/rs
    ## rs <- rowSums(cst.med.mat)
    ## plot(rs)

    list(cst.med.mat=cst.med.mat, x=r$x)
}

## T="T1"; phs=c("Lactobacillus_iners"); ph="Mobiluncus_curtisii_mulieris"

spmrf.ph.fn <- function(ph="Mobiluncus_curtisii_mulieris", nItr=2000, thin=1, nChains=3, nCores=3)
{
    x <- log10(pt.no.mPTB[, ph])
    y <- mt.no.mPTB$sPTB

    idx <- is.finite(x)
    y <- y[idx]
    x <- x[idx]
    o <- order(x)
    y <- y[o]
    x <- x[o]

    m.dat <- list(y=y, xvar1=x, J=length(y))
    m.dat2 <- spmrf.get.data(m.dat)
    m <- sampling(spmrf.bernoulli.o2.model, data=m.dat2, control=list(adapt_delta=0.96, max_treedepth=12), iter=nItr, thin=thin, chains=nChains, cores=nCores)

    list(m=m, x=x, y=y, ph=ph, phs=phs)
}


spmrf.ph.AA.fn <- function(ph="Mobiluncus_curtisii_mulieris", nItr=2000, thin=1, nChains=3, nCores=3)
{
    x <- log10(pt.AA[, ph])
    y <- mt.AA$sPTB
    s <- mt.AA$subjID

    idx <- is.finite(x)
    y <- y[idx]
    x <- x[idx]
    s <- s[idx]
    o <- order(x)
    y <- y[o]
    x <- x[o]
    s <- s[o]

    m.dat <- list(y=y, xvar1=x, J=length(y))
    m.dat2 <- spmrf.get.data(m.dat)
    m <- sampling(spmrf.bernoulli.o2.model, data=m.dat2, control=list(adapt_delta=0.96, max_treedepth=12), iter=nItr, thin=thin, chains=nChains, cores=nCores)

    list(m=m, x=x, y=y, s=s, ph=ph, phs=phs)
}
