## This file is where I do the analysis and explore the MCMC fits from the
## already run models

message("Loading result fits into workspace..")
hake <- readRDS('results/pilot_hake.RDS')
hake2 <- hake##readRDS('results/pilot_hake2.RDS')
halibut <- readRDS('results/pilot_halibut.RDS')
halibut2 <- readRDS('results/pilot_halibut2.RDS')
canary <- readRDS('results/pilot_canary.RDS')
canary2 <- readRDS('results/pilot_canary2.RDS')
snowcrab <- readRDS('results/pilot_snowcrab.RDS')
snowcrab2 <- readRDS('results/pilot_snowcrab2.RDS')
## Now the comparison fits
hakefits <- readRDS('results/hake_fits.RDS')
halibutfits <- readRDS('results/halibut2_fits.RDS')
canaryfits <- readRDS('results/canary2_fits.RDS')
snowcrabfits <- readRDS('results/snowcrab2_fits.RDS')

## Look at convergence and ESS for pilot chains
conv <- ldply(list(hake, hake2, halibut, halibut2, snowcrab, snowcrab2,
                                    canary, canary2), function(x) x$monitor)
conv <- melt(conv, measure.vars=c('Rhat', 'n_eff'))
write.csv(conv, file='results/table_convergence_pilot.csv')
g <- ggplot(conv, aes(version, y=value)) + geom_jitter() +
  facet_grid(variable~model, scales='free') + geom_hline(yintercept=1.1) +
  scale_y_log10()
ggsave('plots/pilot_convergence.png', g, width=7, height=5)

## Look at convergence and ESS for pilot chains
conv <- ldply(c(hakefits, halibutfits, snowcrabfits, canaryfits),
              function(x) x$monitor)
conv <- melt(conv, measure.vars=c('Rhat', 'n_eff'))
conv$version2 <- paste0(conv$alg, "_", conv$metric)
write.csv(conv, file='results/table_convergence_updated.csv')
g <- ggplot(conv, aes(version2, y=value)) + geom_jitter() +
  facet_grid(variable~model, scales='free') + geom_hline(yintercept=1.1) +
  scale_y_log10()
ggsave('plots/comparisons_convergence.png', g, width=7, height=5)

## Quick check that divergences aren't a problem. < 1% seems very
## reasonable
ff <- function(fit){
  x <- extract_sampler_params(fit)
   data.frame(model=fit$model, pct.divs=100*sum(x$divergent__)/nrow(x))
}
rbind(ff(hakefits[[2]]), ff(halibutfits[[2]]), ff(snowcrabfits[[2]]), ff(canaryfits[[2]]))

message("Make model specific plots and diagnostics...")
### Make plots for each model to examine what's going on.
plot.slow(hake)
plot.marginal(hake, save=TRUE)
plot.slow(hake2)
plot.sds(hake2)
xlims <- list(c(0, 9e6), c(0, 2.4), c(0, 2e6))
ylims <- list(c(0, 6e-7), c(0, 3),c(0, 2e-6))
inputs <- list(model='hake', posterior=hakefits[[2]]$dq.post, mle=hake2$dq,
               mle0=hake$dq)
plot.uncertainties(inputs, xlims=xlims, ylims=ylims)
plot.improvement(hake, hake2)
plot.marginal(hake2, save=TRUE)

## For models which need regularization (all but hake), only make the SD
## and management plots for the regularized version
plot.slow(halibut)
plot.marginal(halibut, save=TRUE)
plot.slow(halibut2)
plot.sds(halibut2)
xlims <- list(c(350000, 650000), c(100000, 300000), c(100000, 300000))
ylims <- list(c(0, 2.5e-5), c(0, 5e-5),c(0, 4e-5))
inputs <- list(model='halibut', posterior=halibutfits[[2]]$dq.post, mle=halibut2$dq,
               mle0=halibut$dq)
plot.uncertainties(inputs, xlims=xlims, ylims=ylims)
plot.improvement(halibut, halibut2)
plot.marginal(halibut2, save=TRUE)

plot.slow(canary)
plot.marginal(canary, save=TRUE)
plot.slow(canary2)
plot.sds(canary2)
plot.marginal(canary2, save=TRUE)
xlims <- list(c(0, 1.1), c(2000, 5000), c(0, 3500))
ylims <- list(c(0, 7), c(0, 0.002),c(0, 2e-03))
inputs <- list(model='canary', posterior=canaryfits[[2]]$dq.post, mle=canary2$dq,
               mle0=canary$dq)
plot.uncertainties(inputs, xlims=xlims, ylims=ylims)
plot.improvement(canary, canary2)

plot.slow(snowcrab)
plot.marginal(snowcrab, save=TRUE)
plot.slow(snowcrab2)
plot.sds(snowcrab2)
plot.marginal(snowcrab2, save=TRUE)
xlims <- list(c(250, 350), c(.8, 2), c(15, 40))
ylims <- list(c(0, .06), c(0, 4),c(0, .2))
inputs <- list(model='snowcrab', posterior=snowcrabfits[[2]]$dq.post, mle=snowcrab2$dq,
               mle0=snowcrab$dq)
plot.uncertainties(inputs, xlims=xlims, ylims=ylims)
plot.improvement(snowcrab, snowcrab2)
### End of looking at the pilot chains


### End of analysis
