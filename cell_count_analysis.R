##
## cell count analysis
##

## Are cell counts significantly different between aBV, ND and sBV samples after
## subsampling to 30 women per group?

## 1. stratify sBV, aBV and noBV by CSTs/Lb-CSTs. Are these groups different from the point of view of CSTs?

## 2. compare log( cell counts ) b/w sBV, aBV and noBV groups stratified by CSTs/Lb-CSTs.

## 3. sci tabl: perform (superficial cells)/(total cells) between sBV, aBV and noBV groups stratified by CSTs/Lb-CSTs.

## In all the above analysis 3 comparisons: sBV-aBV, sBV-noBV and aBV-noBV, are
## to be done for CST (with CSTs I, II, III, IIIA, IV, IVB, V, VI) and for
## Lb-CST: (Lb-CSTs, CST-IV).

## Thus, we are looking at 3x3x2 = 18 analyses.


## save(gt, cst, cat, cc, file="hmp_shedding.rda")
load("hmp_shedding.rda")


table(cat)
##  aBV noBV  sBV
##   71  104   17

##
## 1. Stratify sBV, aBV and noBV by CSTs/Lb-CSTs. Are these groups different
##    from the point of view of CSTs?
##




## aBV-noBV
pval.aBV.noBV <- c()
for ( sel.cst in CSTs )
{
    print(sel.cst)
    y <- ifelse(cst==sel.cst, 1, 0)
    x <- mt$HC.type
    s <- mt$subjID
    m.dat <- list(y=y, type=as.integer(factor(x, levels=c("noHC","HC"))), subjID=s, N=length(y), nTypes=2, nSubj=length(unique(s)))
    ## str(m.dat)
    fit <- sampling(btr.m, data=m.dat, iter=nItr, thin=thin, chains=nChains, cores=nChains, control = list(adapt_delta = 0.99, max_treedepth = 15))# init=ini)
    ## m <- summary(fit, pars=c("logRat"), probs=c(0.025,0.975))$summary
    r <- get.stan.pval(fit,"delta")
    if ( !is.na(r$mu.pval.Fn) )
    {
        pvals[sel.cst] <- r$mu.pval.Fn
    } else {
        pvals[sel.cst] <- r$mu.pval
    }
}
pvals


## aBV-sBV


## sBV-noBV
