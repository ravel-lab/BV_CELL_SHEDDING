##
## cell_count_noBV_bw_CSTIV_vs_LbCSTs.R
##

library(car)
library(chron)
library(rstan)
source("stan_utils.R")

source("Pmisc.R") # for bp.plot()

rstan_options(auto_write = TRUE)

dexp.gr.ri.m <- stan_model(file="/Users/pgajer/projects/Statistics/MCMC_Models/dexp_gr_ri.stan")

## gt2$Lb.cst <- ifelse(gt2$cst=="IV","IV","Lb")
## save(gt, gt2, file="hmp_shedding_final.rda")
load("hmp_shedding_final.rda")


idx <- gt2$cat=="noBV" ## gt2$cst=="III" &
sum(idx)
y <- log10(gt2$Mean.cell.count[idx])
gr <- as.integer(factor(gt2$Lb.cst[idx], levels=c("Lb","IV")))

m.dat <- list(y=y,
             N=length(y),
             gr=gr,
             nGr=2,
             subjID=as.integer(factor(gt2$subjID[idx])),
             nSubj=length(unique(gt2$subjID[idx])))
str(m.dat)

nItr <- 5000
thin <- 5
nChains <- 4

fit <- sampling(dexp.gr.ri.m, data=m.dat, iter=nItr, thin=thin, chains=nChains, cores=nChains)#, control = list(adapt_delta = 0.99, max_treedepth = 15))# init=ini)

(d.ci <- summary(fit, pars=c("delta"), probs=c(0.025,0.975))$summary)
##               mean     se_mean         sd        2.5%     97.5%    n_eff     Rhat
## delta[1] 0.1204438 0.003995643 0.06301363 0.006822789 0.2438405 248.7113 1.026427

r <- get.stan.pval(fit,"delta[1]")
r$pval

(t1 <- c(d.ci[,c(1,4,5)], r$pval))
##        mean        2.5%       97.5%
## 0.120443831 0.006822789 0.243840550 0.038374433

(mu.ci <- summary(fit, pars=c("mu"), probs=c(0.025,0.975))$summary)
##           mean     se_mean         sd     2.5%    97.5%    n_eff     Rhat
## mu[1] 1.968987 0.001868389 0.03250048 1.905262 2.029656 302.5833 1.011961
## mu[2] 1.848544 0.002427117 0.05673719 1.734883 1.952309 546.4549 1.015492


file <- "../pics/log10_cell_count_in_noBV_bs_CSTIV_and_LbCST_bpplots.pdf"
pdf(file, width=6, height=6)
op <- par(mar=c(3, 4.4, 2.5, 1.5), mgp=c(2.75,0.6,0),tcl = -0.)
bp.plot.f(factor(gt2$Lb.cst[idx]), y, name = TRUE, main="", ylab="log10( cell count )")
par(op)
dev.off()
file
