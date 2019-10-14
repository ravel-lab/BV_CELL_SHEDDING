##
## Compare log( cell counts ) b/w sBV, aBV and noBV groups stratified by CSTs/Lb-CSTs
##

library(car)
library(chron)
library(rstan)
rstan_options(auto_write = TRUE)

source("stan_utils.R")
source("Pmisc.R") # for bp.plot()

## dexp.3gr.ri.m <- stan_model(file="dexp_3groups_dexp_ri.stan")
load("dexp_3groups_dexp_ri.rds")

## dexp.gr.dexp.ri.m <- stan_model(file="dexp_gr_dexp_ri.stan")
load("dexp_gr_dexp_ri.rds")

## save(gt, cst, cat, cc, CSTs, file="hmp_shedding.rda")
load("hmp_shedding.rda")

## save(gt2, cst2, Lb.cst2, cat2, gt.noBV, gt.sBV, gt.aBV, file="hmp_shedding_v2.rda")
load("hmp_shedding_v2.rda")

cc <- gt$count
cc2 <- gt2$count

## checking out distribution structure of cc2

myHist(cc)
abline(v=50)

myHist(log10(cc))
abline(v=log10(50))

q <- qqPlot(log10(cc), distribution="norm")

qqPlot(log10(cc[cc<=50]), distribution="norm")

qqPlot(log10(cc[cc>50]), distribution="norm")


##
## plotting cell counts within BV groups stratified by Lb-CSTs
##

file <- "../pics/cell_count_by_gr_LbCST_hists.pdf"
pdf(file, width=9, height=6)
op <- par(mfrow=c(2,3), mar=c(3.5, 3.4, 1.5, 1.5), mgp=c(2.25,0.6,0),tcl = -0.3)
for ( cst.type in c("IV","Lb") )
{
    for ( BV.type in c("noBV", "aBV", "sBV") )
    {
        idx <- Lb.cst2==cst.type & cat2==BV.type
        myHist(cc2[idx], main=paste(cst.type,BV.type), xlab="")
    }
}
par(op)
dev.off()
file

cc2.list <- list()
for ( cst.type in c("IV","Lb") )
{
    for ( BV.type in c("noBV", "aBV", "sBV") )
    {
        idx <- Lb.cst2==cst.type & cat2==BV.type
        id <- paste0(cst.type,"::",BV.type)
        cc2.list[[id]] <- cc2[idx]
    }
}



file <- "../pics/cell_count_by_gr_LbCST_bpplots.pdf"
pdf(file, width=9, height=6)
op <- par(mar=c(3.5, 3.4, 1.5, 1.5), mgp=c(2.25,0.6,0),tcl = -0.3)
bp.plot(cc2.list, main="")
par(op)
dev.off()
file


##
## no stratification by CST
##

cc.list <- list()
cc.subj.list <- list()
log10.cc.list <- list()
log10.cc.subj.list <- list()
for ( BV.type in c("noBV", "aBV", "sBV") )
{
    idx <- cat==BV.type
    cc.list[[BV.type]] <- cc[idx]
    log10.cc.list[[BV.type]] <- log10(cc[idx])
    y <- aggregate(x = cc[idx], by = list(gt$subjID[idx]), FUN = "median")
    cc.subj.list[[BV.type]] <- y$x
    log10.cc.subj.list[[BV.type]] <- log10(y$x)
}


file <- "../pics/cell_count_by_gr_bpplots.pdf"
pdf(file, width=9, height=6)
op <- par(mar=c(3.5, 3.4, 1.5, 1.5), mgp=c(2.25,0.6,0),tcl = -0.3)
bp.plot(cc.list, main="")
par(op)
dev.off()
file

file <- "../pics/subj_cell_count_by_gr_bpplots.pdf"
pdf(file, width=9, height=6)
op <- par(mar=c(3.5, 3.4, 1.5, 1.5), mgp=c(2.25,0.6,0),tcl = -0.3)
bp.plot(cc.subj.list, main="")
par(op)
dev.off()
file

## after log transform

file <- "../pics/log10_cell_count_by_gr_bpplots.pdf"
pdf(file, width=9, height=6)
op <- par(mar=c(3.5, 3.4, 1.5, 1.5), mgp=c(2.25,0.6,0),tcl = -0.3)
bp.plot(log10.cc.list, main="", dx=0.05)
par(op)
dev.off()
file

file <- "../pics/log10_subj_cell_count_by_gr_bpplots.pdf"
pdf(file, width=9, height=6)
op <- par(mar=c(3.5, 3.4, 1.5, 1.5), mgp=c(2.25,0.6,0),tcl = -0.3)
bp.plot(log10.cc.subj.list, main="")
par(op)
dev.off()
file




##
## comparison of median cell counts between BV groups - 3 way comparison
##

m.dat <- list(y=log10(cc),
             N=length(cc),
             gr=cat.n,
             nGr=3,
             subjID=as.integer(factor(gt$subjID)),
             nSubj=length(unique(gt$subjID)))
str(m.dat)

nItr <- 5000
thin <- 5
nChains <- 4

## fit.raw <- fit

fit <- sampling(dexp.3gr.ri.m, data=m.dat, iter=nItr, thin=thin, chains=nChains, cores=nChains)#, control = list(adapt_delta = 0.99, max_treedepth = 15))# init=ini)

d1.ci <- summary(fit, pars=c("delta1"), probs=c(0.025,0.975))$summary
d2.ci <- summary(fit, pars=c("delta2"), probs=c(0.025,0.975))$summary
d3.ci <- summary(fit, pars=c("delta3"), probs=c(0.025,0.975))$summary

r1.1 <- get.stan.pval(fit,"delta1[1]")
r1.2 <- get.stan.pval(fit,"delta1[2]")
pvals1 <- c(r1.1$pval, r1.2$pval)

r2.1 <- get.stan.pval(fit,"delta2[1]")
r2.2 <- get.stan.pval(fit,"delta2[2]")
pvals2 <- c(r2.1$pval, r2.2$pval)

r3.1 <- get.stan.pval(fit,"delta3[1]")
r3.2 <- get.stan.pval(fit,"delta3[2]")
pvals3 <- c(r3.1$pval, r3.2$pval)


(t1 <- cbind(d1.ci[,c(1,4,5)], pvals1))
##                 mean       2.5%      97.5%       pvals1
## delta1[1] -0.4046618 -0.5005609 -0.3100666 2.220446e-16
## delta1[2]  0.2647050  0.1756154  0.3579456 7.205353e-09

tblTexFile <- "../docs/Tables/log10_cell_count_noBVref.tex"
tblLabel <- "cc.noBV.ref:tbl"
lt <- latex(formatC(t1, digits=2), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)



(t2 <- cbind(d2.ci[,c(1,4,5)], pvals2))
##                mean      2.5%     97.5%       pvals2
## delta2[1] 0.4046618 0.3100666 0.5005609 2.385188e-16
## delta2[2] 0.6693668 0.5515170 0.7828199 4.505485e-33

tblTexFile <- "../docs/Tables/log10_cell_count_aBVref.tex"
tblLabel <- "cc.aBV.ref:tbl"
lt <- latex(formatC(t2, digits=2), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)



(t3 <- cbind(d3.ci[,c(1,4,5)], pvals3))
##                 mean       2.5%      97.5%       pvals3
## delta3[1] -0.2647050 -0.3579456 -0.1756154 7.205353e-09
## delta3[2] -0.6693668 -0.7828199 -0.5515170 0.000000e+00

tblTexFile <- "../docs/Tables/log10_cell_count_sBVref.tex"
tblLabel <- "cc.sBV.ref:tbl"
lt <- latex(formatC(t3, digits=2), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)



summary(fit, pars=c("mu"), probs=c(0.025,0.975))$summary
##           mean   se_mean       sd      2.5%    97.5%     n_eff      Rhat
## mu[1] 29.73164 0.1190146 4.614073 20.637573 38.41938 1503.0331 1.0012572
## mu[2] 15.55229 0.1516282 4.055113  7.279989 23.04669  715.2301 1.0017392
## mu[3]  3.37011 0.1119280 5.022299 -6.900129 13.43530 2013.3885 0.9996933

##           mean      se_mean         sd     2.5%    97.5%    n_eff      Rhat
## mu[1] 1.936378 0.0006457579 0.02709602 1.882536 1.987478 1760.644 1.0010769
## mu[2] 1.531716 0.0009624521 0.04217905 1.450585 1.614276 1920.592 0.9989891
## mu[3] 2.201083 0.0009714801 0.04161036 2.116588 2.280788 1834.573 1.0011382

describe(log10(cc.list[["noBV"]]))
describe(log10(cc.list[["aBV"]]))
describe(log10(cc.list[["sBV"]]))
##        n  missing distinct     Info     Mean      Gmd      .05      .10      .25      .50
##      104        0       78        1     90.2    52.34    19.30    27.30    57.75    90.50
##      .75      .90      .95
##   126.00   146.00   152.70

## lowest :  11  12  13  19  21, highest: 153 158 186 205 232

describe(cc.subj.list[["noBV"]])
##        n  missing distinct     Info     Mean      Gmd      .05      .10      .25      .50
##       97        0       76        1    90.84    53.66     19.0     26.6     57.0     91.0
##      .75      .90      .95
##    127.0    147.2    153.0

## lowest :  11  12  13  19  21, highest: 153 158 186 205 232



##
## The same but in CST IV samples
##

idx <- Lb.cst2=="IV"

log10.cc.IV.list <- list()
log10.cc.IV.subj.list <- list()
for ( BV.type in c("noBV", "aBV", "sBV") )
{
    idx <- Lb.cst2=="IV" & cat2==BV.type
    log10.cc.IV.list[[BV.type]] <- log10(cc2[idx])
    y <- aggregate(x = cc2[idx], by = list(gt2$subjID[idx]), FUN = "median")
    log10.cc.IV.subj.list[[BV.type]] <- log10(y$x)
}

file <- "../pics/log10_cell_count_in_CST_IV_by_gr_bpplots.pdf"
pdf(file, width=9, height=6)
op <- par(mar=c(3.5, 3.4, 1.5, 1.5), mgp=c(2.25,0.6,0),tcl = -0.3)
bp.plot(log10.cc.IV.list, main="")
par(op)
dev.off()
file

file <- "../pics/log10_subj_cell_count_in_CST_IV_by_gr_bpplots.pdf"
pdf(file, width=9, height=6)
op <- par(mar=c(3.5, 3.4, 1.5, 1.5), mgp=c(2.25,0.6,0),tcl = -0.3)
bp.plot(log10.cc.IV.subj.list, main="")
par(op)
dev.off()
file


##
## CST IV
##

gt2.IV <- gt2[Lb.cst2=="IV",]

m.dat <- list(y=log10(gt2.IV$count),
             N=length(gt2.IV$count),
             gr=cat2.n[Lb.cst2=="IV"],
             nGr=3,
             subjID=as.integer(factor(gt2.IV$subjID)),
             nSubj=length(unique(gt2.IV$subjID)))
str(m.dat)

nItr <- 5000
thin <- 5
nChains <- 4

## fit.raw <- fit

fit.IV <- sampling(dexp.3gr.ri.m, data=m.dat, iter=nItr, thin=thin, chains=nChains, cores=nChains)#, control = list(adapt_delta = 0.99, max_treedepth = 15))# init=ini)

d1.IV.ci <- summary(fit.IV, pars=c("delta1"), probs=c(0.025,0.975))$summary
d2.IV.ci <- summary(fit.IV, pars=c("delta2"), probs=c(0.025,0.975))$summary
d3.IV.ci <- summary(fit.IV, pars=c("delta3"), probs=c(0.025,0.975))$summary

r1.1 <- get.stan.pval(fit.IV,"delta1[1]")
r1.2 <- get.stan.pval(fit.IV,"delta1[2]")
pvals1 <- c(r1.1$pval, r1.2$pval)

r2.1 <- get.stan.pval(fit.IV,"delta2[1]")
r2.2 <- get.stan.pval(fit.IV,"delta2[2]")
pvals2 <- c(r2.1$pval, r2.2$pval)

r3.1 <- get.stan.pval(fit.IV,"delta3[1]")
r3.2 <- get.stan.pval(fit.IV,"delta3[2]")
pvals3 <- c(r3.1$pval, r3.2$pval)


(t1 <- cbind(d1.IV.ci[,c(1,4,5)], pvals1))

tblTexFile <- "../docs/Tables/log10_cell_count_CST_IV_noBVref.tex"
tblLabel <- "cc.noBV.ref:tbl"
lt <- latex(formatC(t1, digits=2), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)



(t2 <- cbind(d2.IV.ci[,c(1,4,5)], pvals2))

tblTexFile <- "../docs/Tables/log10_cell_count_CST_IV_aBVref.tex"
tblLabel <- "cc.aBV.ref:tbl"
lt <- latex(formatC(t2, digits=2), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)



(t3 <- cbind(d3.IV.ci[,c(1,4,5)], pvals3))

tblTexFile <- "../docs/Tables/log10_cell_count_CST_IV_sBVref.tex"
tblLabel <- "cc.sBV.ref:tbl"
lt <- latex(formatC(t3, digits=2), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)


##
## Lb CSTs
##

log10.cc.Lb.list <- list()
log10.cc.Lb.subj.list <- list()
for ( BV.type in c("noBV", "aBV", "sBV") )
{
    idx <- Lb.cst2=="Lb" & cat2==BV.type
    log10.cc.Lb.list[[BV.type]] <- log10(cc2[idx])
    y <- aggregate(x = cc2[idx], by = list(gt2$subjID[idx]), FUN = "median")
    log10.cc.Lb.subj.list[[BV.type]] <- log10(y$x)
}

file <- "../pics/log10_cell_count_in_CST_Lb_by_gr_bpplots.pdf"
pdf(file, width=9, height=6)
op <- par(mar=c(3.5, 3.4, 1.5, 1.5), mgp=c(2.25,0.6,0),tcl = -0.3)
bp.plot(log10.cc.Lb.list, main="")
par(op)
dev.off()
file

file <- "../pics/log10_subj_cell_count_in_CST_Lb_by_gr_bpplots.pdf"
pdf(file, width=9, height=6)
op <- par(mar=c(3.5, 3.4, 1.5, 1.5), mgp=c(2.25,0.6,0),tcl = -0.3)
bp.plot(log10.cc.Lb.subj.list, main="")
par(op)
dev.off()
file


gt2.Lb <- gt2[Lb.cst2=="Lb",]

m.dat <- list(y=log10(gt2.Lb$count),
             N=length(gt2.Lb$count),
             gr=cat2.n[Lb.cst2=="Lb"],
             nGr=3,
             subjID=as.integer(factor(gt2.Lb$subjID)),
             nSubj=length(unique(gt2.Lb$subjID)))
str(m.dat)

nItr <- 5000
thin <- 5
nChains <- 4

fit.Lb <- sampling(dexp.3gr.ri.m, data=m.dat, iter=nItr, thin=thin, chains=nChains, cores=nChains)#, control = list(adapt_delta = 0.99, max_treedepth = 15))# init=ini)

d1.Lb.ci <- summary(fit.Lb, pars=c("delta1"), probs=c(0.025,0.975))$summary
d2.Lb.ci <- summary(fit.Lb, pars=c("delta2"), probs=c(0.025,0.975))$summary
d3.Lb.ci <- summary(fit.Lb, pars=c("delta3"), probs=c(0.025,0.975))$summary

r1.1 <- get.stan.pval(fit.Lb,"delta1[1]")
r1.2 <- get.stan.pval(fit.Lb,"delta1[2]")
pvals1 <- c(r1.1$pval, r1.2$pval)

r2.1 <- get.stan.pval(fit.Lb,"delta2[1]")
r2.2 <- get.stan.pval(fit.Lb,"delta2[2]")
pvals2 <- c(r2.1$pval, r2.2$pval)

r3.1 <- get.stan.pval(fit.Lb,"delta3[1]")
r3.2 <- get.stan.pval(fit.Lb,"delta3[2]")
pvals3 <- c(r3.1$pval, r3.2$pval)


(t1 <- cbind(d1.Lb.ci[,c(1,4,5)], pvals1))

tblTexFile <- "../docs/Tables/log10_cell_count_CST_Lb_noBVref.tex"
tblLabel <- "cc.noBV.ref:tbl"
lt <- latex(formatC(t1, digits=2), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)



(t2 <- cbind(d2.Lb.ci[,c(1,4,5)], pvals2))

tblTexFile <- "../docs/Tables/log10_cell_count_CST_Lb_aBVref.tex"
tblLabel <- "cc.aBV.ref:tbl"
lt <- latex(formatC(t2, digits=2), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)



(t3 <- cbind(d3.Lb.ci[,c(1,4,5)], pvals3))

tblTexFile <- "../docs/Tables/log10_cell_count_CST_Lb_sBVref.tex"
tblLabel <- "cc.sBV.ref:tbl"
lt <- latex(formatC(t3, digits=2), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)


##
## within group CST IV vs Lb comparisons of log cell counts
##


##
## noBV
##

table(gt.noBV$Lb.cst)
## IV Lb
## 31 70

log10.cc.noBV.list <- list()
for ( cst.type in c("IV", "Lb") )
{
    idx <- gt.noBV$Lb.cst==cst.type
    log10.cc.noBV.list[[cst.type]] <- log10(gt.noBV$count[idx])
}

file <- "../pics/log10_cell_count_in_noBV_bpplots.pdf"
pdf(file, width=9, height=6)
op <- par(mar=c(3.5, 3.4, 1.5, 1.5), mgp=c(2.25,0.6,0),tcl = -0.3)
bp.plot(log10.cc.noBV.list, main="")
par(op)
dev.off()
file


m.dat <- list(y=log10(gt.noBV$count),
             N=length(gt.noBV$count),
             gr=ifelse(gt.noBV$Lb.cst=="Lb",2,1),# Lb is the reference class
             nGr=2,
             subjID=as.integer(factor(gt.noBV$subjID)),
             nSubj=length(unique(gt.noBV$subjID)))
str(m.dat)

nItr <- 5000
thin <- 5
nChains <- 4

fit <- sampling(dexp.gr.dexp.ri.m, data=m.dat, iter=nItr, thin=thin, chains=nChains, cores=nChains)

(d.ci <- summary(fit, pars=c("delta"), probs=c(0.025,0.975))$summary)

r <- get.stan.pval(fit,"delta[1]")
(pval <- r$pval)

(et <- as.table(c(d.ci[1,c(1,4,5)], pval)))
names(et)[4] <- "p-val"
names(et)[1] <- "$\\Delta$"

et <- rbind(et)

tblTexFile <- "../docs/Tables/cmp_median_log10_cell_counts_in_noBV_bw_Lb_dom_CST_type.tex"
tblLabel <- "cc.noBV:tbl"
lt <- latex(formatC(et, digits=2), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowname=NULL, caption="",
            where='!htp',vbar=FALSE, latexVerbatim=T)


##
## aBV
##

table(gt.aBV$Lb.cst)
## IV Lb
## 62  4

log10.cc.aBV.list <- list()
for ( cst.type in c("IV", "Lb") )
{
    idx <- gt.aBV$Lb.cst==cst.type
    log10.cc.aBV.list[[cst.type]] <- log10(gt.aBV$count[idx])
}

file <- "../pics/log10_cell_count_in_aBV_bpplots.pdf"
pdf(file, width=9, height=6)
op <- par(mar=c(3.5, 3.4, 1.5, 1.5), mgp=c(2.25,0.6,0),tcl = -0.3)
bp.plot(log10.cc.aBV.list, main="")
par(op)
dev.off()
file


m.dat <- list(y=log10(gt.aBV$count),
             N=length(gt.aBV$count),
             gr=ifelse(gt.aBV$Lb.cst=="Lb",2,1),# Lb is the reference class
             nGr=2,
             subjID=as.integer(factor(gt.aBV$subjID)),
             nSubj=length(unique(gt.aBV$subjID)))
str(m.dat)

nItr <- 5000
thin <- 5
nChains <- 4

fit <- sampling(dexp.gr.dexp.ri.m, data=m.dat, iter=nItr, thin=thin, chains=nChains, cores=nChains)

(d.ci <- summary(fit, pars=c("delta"), probs=c(0.025,0.975))$summary)

r <- get.stan.pval(fit,"delta[1]")
(pval <- r$pval)

(et <- as.table(c(d.ci[1,c(1,4,5)], pval)))
names(et)[4] <- "p-val"
names(et)[1] <- "$\\Delta$"

et <- rbind(et)

tblTexFile <- "../docs/Tables/cmp_median_log10_cell_counts_in_aBV_bw_Lb_dom_CST_type.tex"
tblLabel <- "cc.aBV:tbl"
lt <- latex(formatC(et, digits=2), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowname=NULL, caption="",
            where='!htp',vbar=FALSE, latexVerbatim=T)


##
## sBV
##

table(gt.sBV$Lb.cst)
## IV Lb
## 13  4

log10.cc.sBV.list <- list()
for ( cst.type in c("IV", "Lb") )
{
    idx <- gt.sBV$Lb.cst==cst.type
    log10.cc.sBV.list[[cst.type]] <- log10(gt.sBV$count[idx])
}

file <- "../pics/log10_cell_count_in_sBV_bpplots.pdf"
pdf(file, width=9, height=6)
op <- par(mar=c(3.5, 3.4, 1.5, 1.5), mgp=c(2.25,0.6,0),tcl = -0.3)
bp.plot(log10.cc.sBV.list, main="")
par(op)
dev.off()
file


m.dat <- list(y=log10(gt.sBV$count),
             N=length(gt.sBV$count),
             gr=ifelse(gt.sBV$Lb.cst=="Lb",2,1),# Lb is the reference class
             nGr=2,
             subjID=as.integer(factor(gt.sBV$subjID)),
             nSubj=length(unique(gt.sBV$subjID)))
str(m.dat)

nItr <- 5000
thin <- 5
nChains <- 4

fit <- sampling(dexp.gr.dexp.ri.m, data=m.dat, iter=nItr, thin=thin, chains=nChains, cores=nChains)

(d.ci <- summary(fit, pars=c("delta"), probs=c(0.025,0.975))$summary)

r <- get.stan.pval(fit,"delta[1]")
(pval <- r$pval)

(et <- as.table(c(d.ci[1,c(1,4,5)], pval)))
names(et)[4] <- "p-val"
names(et)[1] <- "$\\Delta$"

et <- rbind(et)

tblTexFile <- "../docs/Tables/cmp_median_log10_cell_counts_in_sBV_bw_Lb_dom_CST_type.tex"
tblLabel <- "cc.sBV:tbl"
lt <- latex(formatC(et, digits=2), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowname=NULL, caption="",
            where='!htp',vbar=FALSE, latexVerbatim=T)
