## This file is where I do the analysis and explore the MCMC fits from the
## already run models
source("startup.R")

### Make plots for each model to examine what's going on.
hake <- readRDS('results/pilot_rwm_hake.RDS')
plot.slow(hake)
## Look at which parameter MLE vs posterior variances are different
plot.sds(hake)
## Compare estimates of DQs
xlims <- list(c(0, 7e6), c(0, 1.5), c(0, 3e6))
ylims <- list(c(0, 6e-7), c(0, 3),c(0, 2e-6))
plot.uncertainties(hake, xlims=xlims, ylims=ylims)
plot.marginal(hake, save=TRUE)

## For models which need regularization (all but hake), only make the SD
## and management plots for the regularized version
halibut <- readRDS('results/pilot_rwm_halibut.RDS')
plot.slow(halibut)
plot.marginal(halibut, save=TRUE)
halibut2 <- readRDS('results/pilot_rwm_halibut2.RDS')
plot.slow(halibut2)
plot.sds(halibut2)
xlims <- list(c(350000, 650000), c(100000, 300000), c(100000, 300000))
ylims <- list(c(0, 2.5e-5), c(0, 5e-5),c(0, 4e-5))
plot.uncertainties(regularized=halibut2, original=halibut, xlims=xlims, ylims=ylims)
plot.improvement(halibut, halibut2)
plot.marginal(halibut2, save=TRUE)

canary <- readRDS('results/pilot_rwm_canary.RDS')
plot.slow(canary)
plot.marginal(canary, save=TRUE)
canary2 <- readRDS('results/pilot_rwm_canary2.RDS')
plot.slow(canary2)
plot.sds(canary2)
plot.marginal(canary2, save=TRUE)
xlims <- list(c(0, 1.1), c(2000, 5000), c(0, 3500))
ylims <- list(c(0, 7), c(0, 0.002),c(0, 2e-03))
plot.uncertainties(original=canary, regularized=canary2,
                   xlims=xlims, ylims=ylims)
plot.improvement(canary, canary2)

snowcrab <- readRDS('results/pilot_rwm_snowcrab.RDS')
plot.slow(snowcrab)
plot.marginal(snowcrab, save=TRUE)
snowcrab2 <- readRDS('results/pilot_rwm_snowcrab2.RDS')
plot.slow(snowcrab2, n.slow=12)
plot.sds(snowcrab2)
plot.marginal(snowcrab2, save=TRUE)
xlims <- list(c(250, 350), c(.8, 2), c(15, 40))
ylims <- list(c(0, .06), c(0, 4),c(0, .2))
plot.uncertainties(snowcrab, snowcrab2, xlims=xlims, ylims=ylims)
plot.improvement(snowcrab, snowcrab2)
### End of looking at the pilot chains


### Examine other things ( in development)

## halibut2 NUTs chains having issues with some pars
fit <- readRDS(file='results/halibut2_fits.RDS')[[1]]

pars <- c("SR_parm[1]", "SR_parm[5]", paste0("recdev2[",22:26, "]"), "selparm[53]")
png('plots/halibut2_divergences.png', width=5, height=5, units='in', res=500)
pairs_admb(fit, pars=pars)
dev.off()

## look at what divergent trajectories SSB look like
post <- extract_samples(fit)
divs <- extract_sampler_params(fit)$divergent__
## psv and .sso files probably overwritten so ready them
adnuts:::.write_psv('halibut2', post, 'halibut2')
## setwd('halibut2');system('halibut2 -mceval');setwd('..')
library(r4ss)
xx <- SSgetMCMC(dir='halibut2')[[1]]
ssb <- xx[,grep("SPB_", x=names(xx))[-(1:2)]]
years <- 1996:2019
png("plots/halibut2_divergences_ssb.png", width=7, height=5, units='in', res=500)
plot(0,0, xlim=range(years), ylim=c(0, 1.1*max(ssb)), type='n',
     xlab='year', ylab='SSB')
ssb2 <- ssb[order(divs),]
divs2 <- sort(divs)
for(i in 1:nrow(ssb2)){
  lines(years, ssb2[i,], col=ifelse(divs2[i]==0, rgb(0,0,0,.1),2),
        lwd=ifelse(divs2[i]==0, 1, 3))
}
dev.off()

## get all the single variables, I manually found the bad ones
## ind <- grep('\\[', xx, invert=TRUE)
## plot.slow(snowcrab, pars=pars, save=FALSE)
## pars <- xx[ind[c(1,2,3,4, 12,13, 29,34,37)]]
xx <- dimnames(snowcrab2$samples)[3][[1]]
post <- extract_samples(snowcrab2)
plot(post$Fem_F_prop_constant)

plot.slow(snowcrab2, fast=FALSE, save=FALSE)
pars <- c("af", "am", "bf", "bm", "log_avg_fmortdf", "log_avg_fmortf", "Fem_F_prop_constant",
"srv3_q_f", "srvind_sel95_f", "srv10ind_q_f")
plot.slow(snowcrab2, pars=pars, save=FALSE)

pars <- grep('matestfe', x=xx)
plot.slow(snowcrab2, pars=pars, save=FALSE)
pars <- grep('mateste', x=xx)
plot.slow(snowcrab2, pars=pars, save=FALSE)

## pars <- grep('recdev', x=xx)
## plot.slow(snowcrab2, pars=pars, save=FALSE)
pars <- grep('rec_devf', x=xx)
plot.slow(snowcrab2, pars=pars[1:18], save=FALSE)
plot.slow(snowcrab2, pars=pars[18:36], save=FALSE)
pars <- grep('fmort_dev', x=xx)
plot.slow(snowcrab2, pars=pars[1:18], save=FALSE)
plot.slow(snowcrab2, pars=pars[18:35], save=FALSE)
pars <- grep('fmortdf_dev', x=xx)
plot.slow(snowcrab2, pars=pars[1:18], save=FALSE)
plot.slow(snowcrab2, pars=pars[18:35], save=FALSE)

pars <- grep('fmortt_dev_era1', x=xx)
plot.slow(snowcrab2, pars=pars, save=FALSE)
pars <- grep('fmortt_dev_era2', x=xx)
plot.slow(snowcrab2, pars=pars[1:13], save=FALSE)
plot.slow(snowcrab2, pars=pars[13:25], save=FALSE)


##pars <- grep('log_sel50_dev_mo', x=xx)
## plot.slow(snowcrab2, pars=pars, save=FALSE)
pars <- grep('selsmo10ind', x=xx)
plot.slow(snowcrab2, pars=pars, save=FALSE)
pars <- grep('selsmo09ind', x=xx)
plot.slow(snowcrab2, pars=pars, save=FALSE)
pars <- xx[1:10]
plot.slow(snowcrab2, pars=pars, save=FALSE)



## maybe look at NUTS vs RWM from the comparison chains?
n.slow <- 10
## Look at NUTS fit
fit.nuts <- readRDS('results/canary2_fits.RDS')[[1]]
fit.nuts$ess <- monitor(fit.nuts$samples, warmup=fit.nuts$warmup, print=FALSE)[,'n_eff']
slow <- names(sort(fit.nuts$ess))[1:n.slow]
pairs_admb(fit.nuts, diag='acf', pars=slow)
fit.rwm <- readRDS('results/canary2_fits.RDS')[[2]]
fit.rwm$ess <- monitor(fit.rwm$samples, warmup=fit.rwm$warmup, print=FALSE)[,'n_eff']
#slow <- names(sort(fit.rwm$ess))[1:n.slow]
pairs_admb(fit.rwm, diag='acf', pars=slow)

n.slow <- 6

halibut <- readRDS('results/pilot_rwm_halibut.RDS')
plot.slow(halibut)

## This is old  code to look at specfic halibut parameters, which now can
## be down with plot.slow.
halibut.post <- extract_samples(halibut.rwm, inc_lp=TRUE)
slow <- names(sort(halibut.rwm$ess))[1:n.slow]
png('plots/pairs.halibut.rwm.png', width=7, height=5, units='in', res=500)
pairs_admb(halibut.post, mle=halibut.rwm$mle, pars=slow);dev.off()
halibut2.rwm <- readRDS('results/pilot_rwm_halibut2.RDS')
chain <- rep(1:dim(halibut.rwm$samples)[2], each=dim(halibut.rwm$samples)[1]-halibut.rwm$warmup)
halibut2.post <- extract_samples(halibut2.rwm, inc_lp=TRUE)
slow <- names(sort(halibut2.rwm$ess))[1:n.slow]
png('plots/pairs.halibut2.rwm.png', width=7, height=5, units='in', res=500)
pairs_admb(halibut2.post, mle=halibut2.rwm$mle, diag='trace', pars=slow);dev.off()
recdev2 <- names(halibut2.post)[4:34][1:7]
png('plots/pairs.halibut2.rwm.recdev2.png', width=7, height=5, units='in', res=500)
pairs_admb(halibut2.post, mle=halibut2.rwm$mle, diag='trace',chain=chain, pars=recdev2);dev.off()
var.post <- apply(extract_samples(halibut2.rwm),2, sd)
var.mle <- halibut2.rwm$mle$se[1:halibut2.rwm$mle$nopar]
vars <- data.frame(post=var.post, mle=var.mle)
g <- ggplot(vars, aes(x=log10(mle), log10(post))) + geom_point(alpha=.7) +
  geom_abline(slope=1) + xlab("MLE Variance") + ylab("Posterior Variance")
ggsave(paste0('plots/vars.halibut2.png'), g, width=7, height=5)
halibut3.rwm <- readRDS('results/pilot_rwm_halibut3.RDS')
chain <- rep(1:dim(halibut.rwm$samples)[2], each=dim(halibut.rwm$samples)[1]-halibut.rwm$warmup)
halibut3.post <- extract_samples(halibut3.rwm, inc_lp=TRUE)
slow <- names(sort(halibut3.rwm$ess))[1:n.slow]
png('plots/pairs.halibut3.rwm.png', width=7, height=5, units='in', res=500)
pairs_admb(halibut3.post, mle=halibut3.rwm$mle, diag='trace', pars=slow);dev.off()
recdev2 <- names(halibut3.post)[4:34][1:7]
png('plots/pairs.halibut3.rwm.recdev2.png', width=7, height=5, units='in', res=500)
pairs_admb(halibut3.post, mle=halibut3.rwm$mle, diag='trace',chain=chain, pars=recdev2);dev.off()
var.post <- apply(extract_samples(halibut3.rwm),2, sd)
var.mle <- halibut3.rwm$mle$se[1:halibut3.rwm$mle$nopar]
vars <- data.frame(post=var.post, mle=var.mle)
g <- ggplot(vars, aes(x=log10(mle), log10(post))) + geom_point(alpha=.7) +
  geom_abline(slope=1) + xlab("MLE Variance") + ylab("Posterior Variance")
ggsave(paste0('plots/vars.halibut3.png'), g, width=7, height=5)
library(vioplot)
png(paste0('plots/ess_improvement_halibut.png'), width=5, height=5,
    units='in', res=500)
vioplot(log10(halibut.rwm$ess), log10(halibut2.rwm$ess),
        log10(halibut3.rwm$ess), names=c("Original", "Fixed", "SigmaR Small"))
mtext('halibut', line=1, cex=1.5)
mtext("log10(ESS)", side=2, line=2.5, cex=1.25)
dev.off()


n.slow <- 16
fit <- readRDS('results/pilot_rwm_snowcrab.RDS')
fit$ess <- monitor(fit$samples, warmup=fit$warmup, print=FALSE)[,'n_eff']
slow <- names(sort(fit$ess, FALSE))[1:n.slow]
png('plots/pairs.snowcrab.slow.png', width=7, height=5, units='in', res=500)
pairs_admb(fit,  diag='trace', pars=slow)
dev.off()
## I found this manually by looking at par file
hitbounds <- sort(c(293,309, 310:312, 323, 324,331,6,268, 269, 284, 286, 287))
png('plots/pairs.snowcrab.hitbounds.png', width=7, height=5, units='in', res=500)
pairs_admb(snowcrab.post, mle=fit$mle, chain=chain, diag='trace', pars=hitbounds);dev.off()
snowcrab2.rwm <- readRDS('results/pilot_rwm_snowcrab2.RDS')
snowcrab2.post <- extract_samples(snowcrab2.rwm, TRUE, inc_lp=TRUE)
chain <- rep(1:dim(snowcrab2.rwm$samples)[2], each=dim(snowcrab2.rwm$samples)[1]-snowcrab2.rwm$warmup)
slow <- names(sort(snowcrab2.rwm$ess, FALSE))[1:n.slow]
png('plots/pairs.snowcrab2.slow.png', width=7, height=5, units='in', res=500)
pairs_admb(snowcrab2.post, mle=snowcrab2.rwm$mle, chain=chain, diag='trace', pars=slow);dev.off()
plot.improvement(fit, snowcrab2.rwm)
## launch_shinyadmb(fit)
## launch_shinyadmb(snowcrab.nuts)




all.fits <- list(cod.rwm, cod.nuts, halibut.rwm, halibut.nuts, hake.rwm,
                 hake.nuts)
perf.wide <- ldply(all.fits, function(x){
  data.frame(model=x$model, alg=x$algorithm, runtime=sum(x$time.total),
        minESS=min(x$ess), maxRhat=max(x$Rhat))})
g <- ggplot(perf.wide, aes(model, maxRhat, color=alg)) + geom_point()
ggsave('plots/maxRhat_comparison.png', g, width=7, height=5)
perf.wide$perf <- with(perf.wide, minESS/runtime)
g <- ggplot(perf.wide, aes(model, y=log(perf), color=alg)) + geom_point()
ggsave('plots/efficienty_comparison.png', g, width=7, height=5)
temp <- reshape2::melt(subset(perf.wide, select=-c(maxRhat)), c('model', 'alg'))
perf.long <- reshape2::dcast(temp, model+variable~alg)
g <- ggplot(perf.long, aes(x=log10(RWM), y=log10(NUTS), color=model)) + geom_point() +
  geom_abline(slope=1)+ facet_wrap('variable', scales='free') + coord_equal()
ggsave('plots/perf_comparison.png', g, width=7, height=5)

