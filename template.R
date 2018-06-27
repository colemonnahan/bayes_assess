set.seed(seed)
seeds <- sample(1:1e4, size=reps)
d <- m

## Run model with hbf=1 to get right covariance matrix and MLEs for NUTS
setwd(d)
if(recompile)
  system(paste('admb', m, '-f')) # -f does optimized mode
## Add the -mcmc so SS turns off bias adjustment and we get a covariance
## that matches this to use for sampling. Otherwise there is a mismatch.
if(m!='snowcrab2'){
  system(paste(m, '-hbf 1 -mcmc 15 '))
} else {
  system(paste(m, '-hbf 1 -ainp snowcrab2.par -phase 50'))
}
setwd('..')

### Run NUTS for different mass matrices
## Mass matrix is MLE covariance
fit.nuts.mle <-
  sample_admb(m, iter=iter, init=inits, algorithm='NUTS',  seeds=seeds,
               parallel=TRUE, chains=reps, warmup=warmup, path=d, cores=reps,
              control=list(max_treedepth=td, metric="mle", adapt_delta=ad))
## Mass matrix is dense one estimated from previous run
fit.nuts.dense <-
  sample_admb(m, iter=iter, init=inits, algorithm='NUTS', seeds=seeds,
               parallel=TRUE, chains=reps, warmup=warmup, path=d, cores=reps,
              control=list(max_treedepth=td, metric=fit.nuts.mle$covar.est, adapt_delta=ad))

## Now run RWM but using a thinning rate similar to NUTS so the time is
## roughly equivalent.
## Rerun model with hbf=0
setwd(m)
if(m!='snowcrab2'){
  system(paste(m, '-hbf 0 -mcmc 15 -nox'))
} else {
  system(paste0(m, ' -nox -phase 50 -ainp ', m,'.par'))
#  system(paste(m, '-hbf 0 -ainp snowcrab2.par -phase 50 -nox'))
}
setwd('..')

tt <- floor(4*mean(extract_sampler_params( fit.nuts.mle)$n_leapfrog__))
fit.rwm.mle <-
  sample_admb(m, iter=tt*iter, init=inits, thin=tt, seeds=seeds,
              parallel=TRUE, chains=reps, warmup=tt*warmup,
              path=d, cores=reps, control=list(metric=NULL),
              algorithm='RWM')
## tt <- floor(4*mean(extract_sampler_params( fit.nuts.dense)$n_leapfrog__))
## fit.rwm.dense <-
##   sample_admb(m, iter=tt*iter, init=inits, thin=tt,  seeds=seeds,
##               parallel=TRUE, chains=reps, warmup=tt*warmup,
##               path=d, cores=reps, control=list(metric=fit.rwm.mle$covar.est),
##               algorithm='RWM')
fit.rwm.dense <- fit.rwm.mle
## Save fits
saveRDS(list(fit.nuts.mle=fit.nuts.mle, fit.nuts.dense=fit.nuts.dense,
             fit.rwm.mle=fit.rwm.mle, fit.rwm.dense=fit.rwm.dense),
        file=paste0('results/', d,'_fits.RDS'))


## Gather adaptation and performance metrics
ff <- function(labels, ...){
  fits <- list(...)
  if(length(fits)!=length(labels)) stop("aa")
  do.call(rbind,  lapply(1:length(fits), function(i){
    x <- fits[[i]]$sampler_params
    data.frame(m=labels[i], chain=1:length(x),
               eps=as.numeric(do.call(rbind, lapply(x, function(l) tail(l[,2],1)))),
               divergences=as.numeric(do.call(rbind, lapply(x, function(l) sum(l[-(1:fits[[i]]$warmup),5])))),
               accept_prob=as.numeric(do.call(rbind, lapply(x,
                          function(l) mean(l[-(1:fits[[i]]$warmup),1])))),
               nsteps=as.numeric(do.call(rbind, lapply(x,
                          function(l) mean(l[-(1:fits[[i]]$warmup),4])))))
  }))}
adaptation <- ff(c("mle", "dense"), fit.nuts.mle, fit.nuts.dense)
stats.nuts.mle <- with(fit.nuts.mle, data.frame(alg='nuts', m='mle', time.total=sum(time.total), rstan::monitor(samples, warmup=warmup, probs=.5, print=FALSE)))
perf.nuts.mle <- data.frame(alg='nuts', m='mle', efficiency=min(stats.nuts.mle$n_eff)/sum(fit.nuts.mle$time.total))
stats.nuts.dense <- with(fit.nuts.dense, data.frame(alg='nuts', m='dense', time.total=sum(time.total), rstan::monitor(samples, warmup=warmup, probs=.5, print=FALSE)))
perf.nuts.dense <- data.frame(alg='nuts', m='dense', efficiency=min(stats.nuts.dense$n_eff)/sum(fit.nuts.dense$time.total))
stats.rwm.mle <- with(fit.rwm.mle, data.frame(alg='rwm', m='mle', time.total=sum(time.total), rstan::monitor(samples, warmup=warmup, probs=.5, print=FALSE)))
perf.rwm.mle <- data.frame(alg='rwm', m='mle', efficiency=min(stats.rwm.mle$n_eff)/sum(fit.rwm.mle$time.total))
stats.rwm.dense <- with(fit.rwm.dense, data.frame(alg='rwm', m='dense', time.total=sum(time.total), rstan::monitor(samples, warmup=warmup, probs=.5, print=FALSE)))
perf.rwm.dense <- data.frame(alg='rwm', m='dense', efficiency=min(stats.rwm.dense$n_eff)/sum(fit.rwm.dense$time.total))
stats.all <- rbind(stats.nuts.mle, stats.nuts.dense, stats.rwm.mle, stats.rwm.dense)
stats.all[,c('mean', 'se_mean', 'sd', 'X50.')] <- NULL
stats.all <- ddply(stats.all, .(alg, m), mutate, perf=(n_eff)/time.total)
stats.long <- reshape2::melt(stats.all, c('alg', 'm'))
perf.all <- rbind(perf.nuts.mle, perf.nuts.dense, perf.rwm.mle,
                  perf.rwm.dense)
adaptation.long <- reshape2::melt(adaptation, c('m', 'chain'))

## Quick plots
ggwidth <- 7
ggheight <- 5
pdf(paste0('plots/', d, '_comparison.pdf'), width=7, height=5)
g <- ggplot(adaptation.long, aes(y=value, x=m)) + geom_point(alpha=.5) +
  facet_wrap('variable', scales='free')
print(g)
g <- ggplot(stats.long, aes(y=value, x=m, color=alg)) + geom_jitter(alpha=.5) +
  facet_wrap('variable', scales='free')
print(g)
g <- ggplot(perf.all, aes(m, efficiency, color=alg)) + geom_point()
print(g)
plot.ess(rwm=fit.rwm.mle, nuts=fit.nuts.mle)
dev.off()

