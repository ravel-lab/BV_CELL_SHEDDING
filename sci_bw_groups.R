##
## perform (superficial cells)/(total cells) = sc2tc.rat, between sBV, aBV and noBV groups stratified by CSTs/Lb-CSTs.
##

library(car)
library(chron)

library(rstan)
rstan_options(auto_write = TRUE)

source("stan_utils.R")
source("Pmisc.R") # for bp.plot()

## bin.tr.ri.m <- stan_model(file="bin_type_ri.stan")
load("bin_type_ri.rds")
## beta.tr.ri.m <- stan_model(file="beta_tr_ri.stan")
load("beta_tr_ri.rds")

## save(gt, cst, cat, cc, CSTs, file="hmp_shedding.rda")
load("hmp_shedding.rda")

## save(gt2, cst2, Lb.cst2, cat2, gt.noBV, gt.sBV, gt.aBV, file="hmp_shedding_v2.rda")
load("hmp_shedding_v2.rda")

## save(st, file="scit.rda")
load("scit.rda")


## checking out distribution structure of sci

myHist(st$sci)

myHist(log10(st$sci))

qqPlot(log10(st$sci[st$sci>0]), distribution="norm")

sum(st$sci==0) #  6
nrow(st)       # 76

table(st$cat)
## aBV noBV
##  53   23

sci.list <- list()
log10.sci.list <- list()
for ( BV.type in c("noBV", "aBV") )
{
    idx <- st$cat==BV.type
    sci.list[[BV.type]] <- st$sci[idx]
    log10.sci.list[[BV.type]] <- log10(st$sci[idx])
}


file <- "../pics/sci_by_gr_bpplots.pdf"
pdf(file, width=9, height=6)
op <- par(mar=c(3.5, 3.4, 1.5, 1.5), mgp=c(2.25,0.6,0),tcl = -0.3)
bp.plot(sci.list, main="")
par(op)
dev.off()
file


file <- "../pics/log10_sci_by_gr_bpplots.pdf"
pdf(file, width=9, height=6)
op <- par(mar=c(3.5, 3.4, 1.5, 1.5), mgp=c(2.25,0.6,0),tcl = -0.3)
bp.plot(log10.sci.list, main="", dx=0.05)
par(op)
dev.off()
file


##
## comparison of median sc2tc.rat(s) between aBV and noBV groups using binomial model
##

## creating 'number of superficial cells' and 'number of total cells' variables
## together with subject ID and category

st$subjID.n <- as.integer(factor(st$subjID))

sc <- c()
tc <- c()
sID <- c()
bvType <- c()
for ( i in seq(nrow(st)) )
{
    sc <- c(sc, st[i,"sp.1"], st[i,"sp.2"], st[i,"sp.3"])
    tc <- c(tc, st[i,"ct.1"], st[i,"ct.2"], st[i,"ct.3"])
    sID <- c(sID, rep(st[i,"subjID.n"],3))
    bvType <- c(bvType, rep(st[i,"cat.n"],3))
}

idx <- tc<sc
sum(idx)

tc[idx] <- 2*sc[idx]

m.dat <- list(y=sc,
             total16S=tc,
             N=length(sc),
             type=bvType,
             nTypes=2,
             subjID=sID,
             nSubj=length(unique(sID)))
str(m.dat)

nItr <- 5000
thin <- 5
nChains <- 4

fit <- sampling(bin.tr.ri.m, data=m.dat, iter=nItr, thin=thin, chains=nChains, cores=nChains)

d.ci <- summary(fit, pars=c("delta"), probs=c(0.025,0.975))$summary

r <- get.stan.pval(fit,"delta[1]")
(pval <- r$pval)

qqPlot(r$mu)

(et <- c(d.ci[,c(1,4,5)], pval))

names(et)[4] <- "p-val"
names(et)[1] <- "$\\Delta$"


tblTexFile <- "../docs/Tables/sci_delta_bin_type_ri_model.tex"
tblLabel <- "sci.delta:tbl"
lt <- latex(formatC(rbind(et), digits=2), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowname=NULL, rowlabel="", caption="",
            where='!htp',vbar=FALSE)


##
## comparison of median sc2tc.rat(s) between aBV and noBV groups using beta model
##

y <- st$sci
ymin <- min(y[y>0])
y[y==0] <- runif(sum(y==0), min=0, max=ymin)

m.dat <- list(y=y,
             N=length(st$sci),
             tr=st$cat.n,
             nTr=2,
             subjID=as.integer(factor(st$subjID)),
             nSubj=length(unique(st$subjID)))
str(m.dat)

nItr <- 5000
thin <- 5
nChains <- 4

fit <- sampling(beta.tr.ri.m, data=m.dat, iter=nItr, thin=thin, chains=nChains, cores=nChains)

d.ci <- summary(fit, pars=c("delta"), probs=c(0.025,0.975))$summary

r <- get.stan.pval(fit,"delta[1]")
(pval <- r$pval)

qqPlot(r$mu)

(et <- c(d.ci[,c(1,4,5)], pval))

names(et)[4] <- "p-val"
names(et)[1] <- "$\\Delta$"


tblTexFile <- "../docs/Tables/sci_delta.tex"
tblLabel <- "sci.delta:tbl"
lt <- latex(formatC(rbind(et), digits=2), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowname=NULL, rowlabel="", caption="",
            where='!htp',vbar=FALSE)



(t2 <- cbind(d2.ci[,c(1,4,5)], pvals2))
##                mean      2.5%     97.5%       pvals2
## delta2[1] 0.4046618 0.3100666 0.5005609 2.385188e-16
## delta2[2] 0.6693668 0.5515170 0.7828199 4.505485e-33

tblTexFile <- "../docs/Tables/log10_sci_aBVref.tex"
tblLabel <- "sci.aBV.ref:tbl"
lt <- latex(formatC(t2, digits=2), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)



(t3 <- cbind(d3.ci[,c(1,4,5)], pvals3))
##                 mean       2.5%      97.5%       pvals3
## delta3[1] -0.2647050 -0.3579456 -0.1756154 7.205353e-09
## delta3[2] -0.6693668 -0.7828199 -0.5515170 0.000000e+00

tblTexFile <- "../docs/Tables/log10_sci_sBVref.tex"
tblLabel <- "sci.sBV.ref:tbl"
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

describe(log10(sci.list[["noBV"]]))
describe(log10(sci.list[["aBV"]]))
describe(log10(sci.list[["sBV"]]))
##        n  missing distinct     Info     Mean      Gmd      .05      .10      .25      .50
##      104        0       78        1     90.2    52.34    19.30    27.30    57.75    90.50
##      .75      .90      .95
##   126.00   146.00   152.70

## lowest :  11  12  13  19  21, highest: 153 158 186 205 232

describe(sci.subj.list[["noBV"]])
##        n  missing distinct     Info     Mean      Gmd      .05      .10      .25      .50
##       97        0       76        1    90.84    53.66     19.0     26.6     57.0     91.0
##      .75      .90      .95
##    127.0    147.2    153.0

## lowest :  11  12  13  19  21, highest: 153 158 186 205 232



##
## The same but in CST IV samples
##

idx <- Lb.cst2=="IV"

log10.sci.IV.list <- list()
log10.sci.IV.subj.list <- list()
for ( BV.type in c("noBV", "aBV", "sBV") )
{
    idx <- Lb.cst2=="IV" & cat2==BV.type
    log10.sci.IV.list[[BV.type]] <- log10(sci2[idx])
    y <- aggregate(x = sci2[idx], by = list(gt2$subjID[idx]), FUN = "median")
    log10.sci.IV.subj.list[[BV.type]] <- log10(y$x)
}

file <- "../pics/log10_sci_in_CST_IV_by_gr_bpplots.pdf"
pdf(file, width=9, height=6)
op <- par(mar=c(3.5, 3.4, 1.5, 1.5), mgp=c(2.25,0.6,0),tcl = -0.3)
bp.plot(log10.sci.IV.list, main="")
par(op)
dev.off()
file

file <- "../pics/log10_subj_sci_in_CST_IV_by_gr_bpplots.pdf"
pdf(file, width=9, height=6)
op <- par(mar=c(3.5, 3.4, 1.5, 1.5), mgp=c(2.25,0.6,0),tcl = -0.3)
bp.plot(log10.sci.IV.subj.list, main="")
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

tblTexFile <- "../docs/Tables/log10_sci_CST_IV_noBVref.tex"
tblLabel <- "sci.noBV.ref:tbl"
lt <- latex(formatC(t1, digits=2), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)



(t2 <- cbind(d2.IV.ci[,c(1,4,5)], pvals2))

tblTexFile <- "../docs/Tables/log10_sci_CST_IV_aBVref.tex"
tblLabel <- "sci.aBV.ref:tbl"
lt <- latex(formatC(t2, digits=2), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)



(t3 <- cbind(d3.IV.ci[,c(1,4,5)], pvals3))

tblTexFile <- "../docs/Tables/log10_sci_CST_IV_sBVref.tex"
tblLabel <- "sci.sBV.ref:tbl"
lt <- latex(formatC(t3, digits=2), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)


##
## Lb CSTs
##

log10.sci.Lb.list <- list()
log10.sci.Lb.subj.list <- list()
for ( BV.type in c("noBV", "aBV", "sBV") )
{
    idx <- Lb.cst2=="Lb" & cat2==BV.type
    log10.sci.Lb.list[[BV.type]] <- log10(sci2[idx])
    y <- aggregate(x = sci2[idx], by = list(gt2$subjID[idx]), FUN = "median")
    log10.sci.Lb.subj.list[[BV.type]] <- log10(y$x)
}

file <- "../pics/log10_sci_in_CST_Lb_by_gr_bpplots.pdf"
pdf(file, width=9, height=6)
op <- par(mar=c(3.5, 3.4, 1.5, 1.5), mgp=c(2.25,0.6,0),tcl = -0.3)
bp.plot(log10.sci.Lb.list, main="")
par(op)
dev.off()
file

file <- "../pics/log10_subj_sci_in_CST_Lb_by_gr_bpplots.pdf"
pdf(file, width=9, height=6)
op <- par(mar=c(3.5, 3.4, 1.5, 1.5), mgp=c(2.25,0.6,0),tcl = -0.3)
bp.plot(log10.sci.Lb.subj.list, main="")
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

tblTexFile <- "../docs/Tables/log10_sci_CST_Lb_noBVref.tex"
tblLabel <- "sci.noBV.ref:tbl"
lt <- latex(formatC(t1, digits=2), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)



(t2 <- cbind(d2.Lb.ci[,c(1,4,5)], pvals2))

tblTexFile <- "../docs/Tables/log10_sci_CST_Lb_aBVref.tex"
tblLabel <- "sci.aBV.ref:tbl"
lt <- latex(formatC(t2, digits=2), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowlabel="", caption="",
            where='!htp',vbar=FALSE)



(t3 <- cbind(d3.Lb.ci[,c(1,4,5)], pvals3))

tblTexFile <- "../docs/Tables/log10_sci_CST_Lb_sBVref.tex"
tblLabel <- "sci.sBV.ref:tbl"
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

log10.sci.noBV.list <- list()
for ( cst.type in c("IV", "Lb") )
{
    idx <- gt.noBV$Lb.cst==cst.type
    log10.sci.noBV.list[[cst.type]] <- log10(gt.noBV$count[idx])
}

file <- "../pics/log10_sci_in_noBV_bpplots.pdf"
pdf(file, width=9, height=6)
op <- par(mar=c(3.5, 3.4, 1.5, 1.5), mgp=c(2.25,0.6,0),tcl = -0.3)
bp.plot(log10.sci.noBV.list, main="")
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

tblTexFile <- "../docs/Tables/cmp_median_log10_scis_in_noBV_bw_Lb_dom_CST_type.tex"
tblLabel <- "sci.noBV:tbl"
lt <- latex(formatC(et, digits=2), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowname=NULL, caption="",
            where='!htp',vbar=FALSE, latexVerbatim=T)


##
## aBV
##

table(gt.aBV$Lb.cst)
## IV Lb
## 62  4

log10.sci.aBV.list <- list()
for ( cst.type in c("IV", "Lb") )
{
    idx <- gt.aBV$Lb.cst==cst.type
    log10.sci.aBV.list[[cst.type]] <- log10(gt.aBV$count[idx])
}

file <- "../pics/log10_sci_in_aBV_bpplots.pdf"
pdf(file, width=9, height=6)
op <- par(mar=c(3.5, 3.4, 1.5, 1.5), mgp=c(2.25,0.6,0),tcl = -0.3)
bp.plot(log10.sci.aBV.list, main="")
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

tblTexFile <- "../docs/Tables/cmp_median_log10_scis_in_aBV_bw_Lb_dom_CST_type.tex"
tblLabel <- "sci.aBV:tbl"
lt <- latex(formatC(et, digits=2), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowname=NULL, caption="",
            where='!htp',vbar=FALSE, latexVerbatim=T)


##
## sBV
##

table(gt.sBV$Lb.cst)
## IV Lb
## 13  4

log10.sci.sBV.list <- list()
for ( cst.type in c("IV", "Lb") )
{
    idx <- gt.sBV$Lb.cst==cst.type
    log10.sci.sBV.list[[cst.type]] <- log10(gt.sBV$count[idx])
}

file <- "../pics/log10_sci_in_sBV_bpplots.pdf"
pdf(file, width=9, height=6)
op <- par(mar=c(3.5, 3.4, 1.5, 1.5), mgp=c(2.25,0.6,0),tcl = -0.3)
bp.plot(log10.sci.sBV.list, main="")
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

tblTexFile <- "../docs/Tables/cmp_median_log10_scis_in_sBV_bw_Lb_dom_CST_type.tex"
tblLabel <- "sci.sBV:tbl"
lt <- latex(formatC(et, digits=2), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowname=NULL, caption="",
            where='!htp',vbar=FALSE, latexVerbatim=T)



##
## restricting to CST IV or Lb dom + comparisons b/w them
##

idx <- st$cst!=""
sum(idx)
st2 <- st[idx,]

st2$Lb.cst <- ifelse(st2$cst=="IV","IV","Lb")

(f <- table(st2$cat, st2$Lb.cst))
##      IV Lb
## aBV  48  2
## noBV 11 12

tblTexFile <- "../docs/Tables/sci_cat_vs_LbCSTs_freq.tex"
tblLabel <- "sci.cst.CSTs:tbl"
lt <- latex(f, file=tblTexFile, label=tblLabel,
           caption.loc='bottom',longtable=FALSE, caption="",
           where='!htp',vbar=FALSE, latexVerbatim=T)


##
## Superficial cells proportions within CST IV
##

sc <- c()
tc <- c()
sID <- c()
Lb.cst <- c()
bvType <- c()
for ( i in seq(nrow(st2)) )
{
    sc <- c(sc, st2[i,"sp.1"], st2[i,"sp.2"], st2[i,"sp.3"])
    tc <- c(tc, st2[i,"ct.1"], st2[i,"ct.2"], st2[i,"ct.3"])
    sID <- c(sID, rep(st2[i,"subjID"],3))
    bvType <- c(bvType, rep(st2[i,"cat.n"],3))
    Lb.cst <- c(Lb.cst, rep(st2[i,"Lb.cst"],3))
}

idx <- tc<sc
sum(idx)
tc[idx] <- 2*sc[idx]

idx <- Lb.cst=="IV"
sum(idx) # 177


m.dat <- list(y=sc[idx],
             total16S=tc[idx],
             N=length(sc[idx]),
             type=bvType[idx],
             nTypes=2,
             subjID=as.integer(factor(sID[idx])),
             nSubj=length(unique(sID[idx])))
str(m.dat)

nItr <- 10000
thin <- 10
nChains <- 4

fit <- sampling(bin.tr.ri.m, data=m.dat, iter=nItr, thin=thin, chains=nChains, cores=nChains)

d.ci <- summary(fit, pars=c("delta"), probs=c(0.025,0.975))$summary

r <- get.stan.pval(fit,"delta[1]")
(pval <- r$pval)

qqPlot(r$mu)

(et <- c(d.ci[,c(1,4,5)], pval))

names(et)[4] <- "p-val"
names(et)[1] <- "$\\Delta$"


tblTexFile <- "../docs/Tables/sci_CSTIV_delta_bin_type_ri_model.tex"
tblLabel <- "sci.delta.cstIV:tbl"
lt <- latex(formatC(rbind(et), digits=2), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowname=NULL, rowlabel="", caption="",
            where='!htp',vbar=FALSE)

sci.cstIV.list <- list()
for ( BV.type in c("noBV", "aBV") )
{
    idx <- st2$cat==BV.type & st2$Lb.cst=="IV"
    sci.cstIV.list[[BV.type]] <- st2$sci[idx]
}


file <- "../pics/sci_CSTIV_by_gr_bpplots.pdf"
pdf(file, width=9, height=6)
op <- par(mar=c(3.5, 3.4, 1.5, 1.5), mgp=c(2.25,0.6,0),tcl = -0.3)
bp.plot(sci.cstIV.list, main="")
par(op)
dev.off()
file


##
## Within noBV group Lactobacillus dominated vs CTS IV superficial
## cell proportions comparison
##

sc <- c()
tc <- c()
sID <- c()
Lb.cst <- c()
bvType <- c()
for ( i in seq(nrow(st2)) )
{
    sc <- c(sc, st2[i,"sp.1"], st2[i,"sp.2"], st2[i,"sp.3"])
    tc <- c(tc, st2[i,"ct.1"], st2[i,"ct.2"], st2[i,"ct.3"])
    sID <- c(sID, rep(st2[i,"subjID"],3))
    bvType <- c(bvType, rep(st2[i,"cat.n"],3))
    Lb.cst <- c(Lb.cst, rep(st2[i,"Lb.cst"],3))
}

idx <- tc<sc
sum(idx)
tc[idx] <- 2*sc[idx]

##
## within noBV analysis
##
idx <- bvType==2
sum(idx) # 69

m.dat <- list(y=sc[idx],
             total16S=tc[idx],
             N=length(sc[idx]),
             type=ifelse(Lb.cst[idx]=="Lb",1,2),
             nTypes=2,
             subjID=as.integer(factor(sID[idx])),
             nSubj=length(unique(sID[idx])))
str(m.dat)

nItr <- 10000
thin <- 10
nChains <- 4

fit <- sampling(bin.tr.ri.m, data=m.dat, iter=nItr, thin=thin, chains=nChains, cores=nChains)

d.ci <- summary(fit, pars=c("delta"), probs=c(0.025,0.975))$summary

r <- get.stan.pval(fit,"delta[1]")
(pval <- r$pval)

qqPlot(r$mu)

(et <- c(d.ci[,c(1,4,5)], pval))

names(et)[4] <- "p-val"
names(et)[1] <- "$\\Delta$"


tblTexFile <- "../docs/Tables/within_noBV_between_CSTs_sci_comp.tex"
tblLabel <- "within.noBV.bw.CSTs.sci.comp:tbl"
lt <- latex(formatC(rbind(et), digits=2), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowname=NULL, rowlabel="", caption="",
            where='!htp',vbar=FALSE)

sci.cstIV.list <- list()
for ( BV.type in c("noBV", "aBV") )
{
    idx <- st2$cat==BV.type & st2$Lb.cst=="IV"
    sci.cstIV.list[[BV.type]] <- st2$sci[idx]
}


within.noBV.sci.list <- list()
for ( cst.type in c("IV", "Lb") )
{
    idx <- st2$cat=="noBV" & st2$Lb.cst==cst.type
    within.noBV.sci.list[[cst.type]] <- st2$sci[idx]
}


file <- "../pics/within_noBV_sci_cst_cmp_bpplots.pdf"
pdf(file, width=9, height=6)
op <- par(mar=c(3.5, 3.4, 1.5, 1.5), mgp=c(2.25,0.6,0),tcl = -0.3)
bp.plot(within.noBV.sci.list, main="")
par(op)
dev.off()
file

##
## within aBV analysis
##
idx <- bvType==1
sum(idx) # 150

m.dat <- list(y=sc[idx],
             total16S=tc[idx],
             N=length(sc[idx]),
             type=ifelse(Lb.cst[idx]=="Lb",1,2),
             nTypes=2,
             subjID=as.integer(factor(sID[idx])),
             nSubj=length(unique(sID[idx])))
str(m.dat)

nItr <- 10000
thin <- 10
nChains <- 4

fit <- sampling(bin.tr.ri.m, data=m.dat, iter=nItr, thin=thin, chains=nChains, cores=nChains)

d.ci <- summary(fit, pars=c("delta"), probs=c(0.025,0.975))$summary

r <- get.stan.pval(fit,"delta[1]")
(pval <- r$pval)

qqPlot(r$mu)

(et <- c(d.ci[,c(1,4,5)], pval))

names(et)[4] <- "p-val"
names(et)[1] <- "$\\Delta$"


tblTexFile <- "../docs/Tables/within_aBV_between_CSTs_sci_comp.tex"
tblLabel <- "within.aBV.bw.CSTs.sci.comp:tbl"
lt <- latex(formatC(rbind(et), digits=2), file=tblTexFile, label=tblLabel,
            caption.loc='bottom',longtable=FALSE, rowname=NULL, rowlabel="", caption="",
            where='!htp',vbar=FALSE)

within.aBV.sci.list <- list()
for ( cst.type in c("IV", "Lb") )
{
    idx <- st2$cat=="aBV" & st2$Lb.cst==cst.type
    within.aBV.sci.list[[cst.type]] <- st2$sci[idx]
}


file <- "../pics/within_aBV_sci_cst_cmp_bpplots.pdf"
pdf(file, width=9, height=6)
op <- par(mar=c(3.5, 3.4, 1.5, 1.5), mgp=c(2.25,0.6,0),tcl = -0.3)
bp.plot(within.aBV.sci.list, main="")
par(op)
dev.off()
file
