## This file is where I do the analysis and explore the MCMC fits from the
## already run models
source("startup.R")
library(adnuts)

hake.rwm <- readRDS('results/pilot_rwm_hake.RDS')
hake.post <- extract_samples(hake.rwm, inc_lp=TRUE)
hake.post <- cbind(hake.post, hake.rwm$dq.post)
mle <- hake.rwm$mle
## Run mcsave and get generated quantities
chain <- rep(1:dim(hake.rwm$samples)[2], each=dim(hake.rwm$samples)[1]-hake.rwm$warmup)
slow <- c(names(sort(hake.rwm$ess))[1:n.slow], names(hake.rwm$dq.post), 'lp__')
png('plots/pairs.hake.slow.png', width=7, height=5, units='in', res=500)
pairs_admb(hake.post, mle=mle, chains=chain, pars=slow, diag='acf')
dev.off()
## Look at which parameter MLE vs posterior variances are different
var.post <- apply(extract_samples(hake.rwm),2, var)
var.mle <- diag(hake.rwm$mle$cov)[1:hake.rwm$mle$nopar]
vars <- data.frame(post=var.post, mle=var.mle)
g <- ggplot(vars, aes(x=log10(mle), log10(post))) + geom_point(alpha=.7) +
  geom_abline(slope=1) + xlab("MLE Variance") + ylab("Posterior Variance")
ggsave(paste0('plots/vars.', m, '.png'), g, width=7, height=5)
## Compare estimates of DQs
xlims <- list(c(0, 7e6), c(0, 1.5), c(0, 3e6))
ylims <- list(c(0, 6e-7), c(0, 3),c(0, 2e-6))
plot.uncertainties(hake.rwm, xlims=xlims, ylims=ylims)
hake.nuts <- readRDS('results/pilot_nuts_hake.RDS')
divs <- extract_sampler_params(hake.nuts)$divergent__
hake.post <- extract_samples(hake.nuts, inc_lp=TRUE)
pairs_admb(hake.post, mle=mle, chains=chain, pars=slow, diag='acf', divergences=divs)
plot.ess(rwm=hake.rwm, nuts=hake.nuts)

n.slow <- 10
canary.rwm <- readRDS('results/pilot_rwm_canary.RDS')
canary.post <- extract_samples(canary.rwm, inc_lp=TRUE)
#canary.post <- cbind(canary.post[-1,], canary.rwm$dq.post)
mle <- canary.rwm$mle
## Run mcsave and get generated quantities
chain <- rep(1:dim(canary.rwm$samples)[2], each=dim(canary.rwm$samples)[1]-canary.rwm$warmup)
slow <- c(names(sort(canary.rwm$ess))[1:n.slow])
png('plots/pairs.canary.slow.png', width=7, height=5, units='in', res=500)
pairs_admb(canary.post, mle=mle, chains=chain, pars=slow, diag='acf')
dev.off()
## Look at which parameter MLE vs posterior variances are different
var.post <- apply(extract_samples(canary.rwm),2, var)
var.mle <- canary.rwm$mle$se[1:canary.rwm$mle$nopar]
vars <- data.frame(post=var.post, mle=var.mle)
g <- ggplot(vars, aes(x=log10(mle), log10(post))) + geom_point(alpha=.7) +
  geom_abline(slope=1) + xlab("MLE Variance") + ylab("Posterior Variance")
ggsave(paste0('plots/vars.canary.png'), g, width=7, height=5)
## ## Compare estimates of DQs
## xlims <- list(c(0, 7e6), c(0, 1.5), c(0, 3e6))
## ylims <- list(c(0, 6e-7), c(0, 3),c(0, 2e-6))
## plot.uncertainties(canary.rwm, xlims=xlims, ylims=ylims)
canary2.rwm <- readRDS('results/pilot_rwm_canary2.RDS')
canary2.post <- extract_samples(canary2.rwm, inc_lp=F)
#canary2.post <- cbind(canary2.post[-1,], canary2.rwm$dq.post)
mle <- canary2.rwm$mle
## Run mcsave and get generated quantities
chain <- rep(1:dim(canary2.rwm$samples)[2], each=dim(canary2.rwm$samples)[1]-canary2.rwm$warmup)
slow <- c(names(sort(canary2.rwm$ess))[1:n.slow])
png('plots/pairs.canary2.slow.png', width=7, height=5, units='in', res=500)
pairs_admb(canary2.post, mle=mle, chains=chain, pars=slow, diag='trace')
dev.off()
## Look at which parameter MLE vs posterior variances are different
var.post <- apply(extract_samples(canary2.rwm),2, sd)
var.mle <- canary2.rwm$mle$se[1:canary2.rwm$mle$nopar]
vars <- data.frame(post=var.post, mle=var.mle)
g <- ggplot(vars, aes(x=log10(mle), log10(post))) + geom_point(alpha=.7) +
  geom_abline(slope=1) + xlab("MLE Variance") + ylab("Posterior Variance")
ggsave(paste0('plots/vars.canary2.png'), g, width=7, height=5)
plot.improvement(canary.rwm, canary2.rwm)
## ## Compare estimates of DQs
## canary.rwm$dq
## canary2.rwm$dq

## Look at NUTS fit
n.slow <- 10
fit.nuts <- readRDS('results/canary2_fits.RDS')[[1]]
fit.nuts$ess <- monitor(fit.nuts$samples, warmup=fit.nuts$warmup, print=FALSE)[,'n_eff']
slow <- names(sort(fit.nuts$ess))[1:n.slow]
pairs_admb(fit.nuts, diag='acf', pars=slow)
fit.rwm <- readRDS('results/canary2_fits.RDS')[[2]]
fit.rwm$ess <- monitor(fit.rwm$samples, warmup=fit.rwm$warmup, print=FALSE)[,'n_eff']
#slow <- names(sort(fit.rwm$ess))[1:n.slow]
pairs_admb(fit.rwm, diag='acf', pars=slow)



n.slow <- 6
halibut.rwm <- readRDS('results/pilot_rwm_halibut.RDS')
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

n.slow <- 15
tanner.rwm <- readRDS('results/pilot_rwm_tanner.RDS')
tanner.post <- extract_samples(tanner.rwm, inc_lp=TRUE)
chain <- rep(1:dim(tanner.rwm$samples)[2], each=dim(tanner.rwm$samples)[1]-tanner.rwm$warmup)
slow <- names(sort(tanner.rwm$ess, FALSE))[1:n.slow]
png('plots/pairs.tanner.slow.png', width=7, height=5, units='in', res=500)
pairs_admb(tanner.post, mle=tanner.rwm$mle, chain=chain, diag='trace', pars=slow);dev.off()
## I found this manually by looking at par file
hitbounds <- c(1,7,33, 67,70)
png('plots/pairs.tanner.rwm.hitbounds.png', width=7, height=5, units='in', res=500)
pairs_admb(tanner.post, mle=tanner.rwm$mle, chain=chain, diag='trace',
           pars=hitbounds)
tanner2.rwm <- readRDS('results/pilot_rwm_tanner2.RDS')
tanner2.post <- extract_samples(tanner2.rwm, inc_lp=TRUE)
chain <- rep(1:dim(tanner2.rwm$samples)[2], each=dim(tanner2.rwm$samples)[1]-tanner2.rwm$warmup)
slow <- names(sort(tanner2.rwm$ess, FALSE))[1:n.slow]
png('plots/pairs.tanner2.slow.png', width=7, height=5, units='in', res=500)
pairs_admb(tanner2.post, mle=tanner2.rwm$mle, chain=chain, diag='trace', pars=slow);dev.off()
plot.improvement(tanner.rwm, tanner2.rwm)


## Look at which parameter MLE vs posterior variances are different
var.poster <- apply(extract_samples(snowcrab.rwm),2, var)
var.mle <- diag(snowcrab.rwm$mle$cov)[1:snowcrab.rwm$mle$nopar]
plot(log10(var.mle), log10(var.poster)); abline(0,1)

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

