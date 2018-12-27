## This script is used to test recdev behavior for the halibut model
## between RWM and different NUTS settings.

## Run model with hbf=1 to get right covariance matrix and MLEs
d <- m <- 'halibut3'
setwd(d)
system(paste(m, '-hbf 1 -nox -mcmc 10'))
setwd('..')

## Fit with two levels of adapt_delta.
mat <- 'mle' # use mle matrix
fit.nuts1 <-
  sample_admb(m, iter=iter, init=inits, parallel=TRUE, chains=reps,
              warmup=warmup, path=d, cores=reps, seeds=1:reps,
              control=list(max_treedepth=td, metric=mat,
                           adapt_delta=.8))
fit.nuts2 <-
  sample_admb(m, iter=iter, init=inits, parallel=TRUE, chains=reps,
              warmup=warmup, path=d, cores=reps, seeds=seeds,
              control=list(max_treedepth=td, metric=mat,
                           adapt_delta=.98))

## Now run RWM but using a thinning rate similar to NUTS so the time is
## roughly equivalent.
tt <- 4*floor(mean(extract_sampler_params( fit.nuts2)$n_leapfrog__))
## Rerun model with hbf=0
setwd(d)
system(paste(m, '-nox -mcmc 10 -phase 10'))
setwd('..')
fit.rwm <-
  sample_admb(m, iter=tt*iter, init=inits, thin=tt, seeds=1:reps,
              parallel=TRUE, chains=reps, warmup=tt*warmup,
              path=d, cores=reps, control=list(metric=NULL),
              algorithm='RWM')

post.nuts1 <-
  cbind(alg='NUTS, 0.8',
        div=extract_sampler_params(fit.nuts1)$divergent__,
        extract_samples(fit.nuts1))
post.nuts2 <-
  cbind(alg='NUTS, 0.99',
        div=extract_sampler_params(fit.nuts2)$divergent__,
        extract_samples(fit.nuts2))
post.rwm <- cbind(alg='RWM', div=0, extract_samples(fit.rwm))
df <- rbind(post.nuts1, post.nuts2, post.rwm)
## Trick to get ggplot to read par names with brackets
p1 <- "recdev2[7]"
p2 <- "recdev2[8]"
df$r1 <- df[,p1]
df$r2 <- df[,p2]
g <- ggplot(df, aes(x=r1, y=r2, color=factor(div))) + geom_point(alpha=.5, size=.1) +
  facet_wrap('alg', nrow=3)
g
## ggsave("plots/hake_tests_recdevs.png", g, width=7, height=5)

## Save fits
saveRDS(df, file=paste0('results/halibut_bias_fits.RDS'))

## pairs_admb(fit=fit.nuts2, pars=grep('recdev2', names(post.nuts1))[4:8])
## pairs_admb(fit=fit.nuts2, pars=grep('recdev2', names(post.nuts1))[11:20])
## pairs_admb(fit=fit.nuts2, pars=grep('recdev2', names(post.nuts1))[21:30])
