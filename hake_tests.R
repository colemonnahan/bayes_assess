## This script is used to test recdev behavior for the hake model between
## RWM and different NUTS settings
source("startup.R")
iter <- 1000
warmup <- iter/4
chains <- 5
td <- 10
inits <- get.inits('hake', chains, 234)

## Run model with hbf=1 to get right covariance matrix and MLEs
d <- m <- 'hake'
setwd(d)
system(paste(m, '-hbf 1 -nox -mcmc 10'))
setwd('..')

## Fit with a low acceptance rate
mat <- 'mle' # for now using mle matrix
fit.nuts1 <-
  sample_admb(m, iter=iter, init=inits, parallel=TRUE, chains=chains,
              warmup=warmup, path=d, cores=chains,
              control=list(max_treedepth=td, metric=mat,
                           adapt_delta=.8)
fit.nuts2 <-
  sample_admb(m, iter=iter, init=inits, parallel=TRUE, chains=chains,
              warmup=warmup, path=d, cores=chains,
              control=list(max_treedepth=td, metric=mat,
                           adapt_delta=.9)
fit.nuts3 <-
  sample_admb(m, iter=iter, init=inits, parallel=TRUE, chains=chains,
              warmup=warmup, path=d, cores=chains,
              control=list(max_treedepth=td, metric=mat,
                           adapt_delta=.98)

## Now run RWM but using a thinning rate similar to NUTS so the time is
## roughly equivalent.
tt <- 4*floor(mean(extract_sampler_params( fit.nuts3)$n_leapfrog__))
## Rerun model with hbf=0
setwd(d)
system(paste(m, '-nox -mcmc 10 -phase 10'))
setwd('..')
fit.rwm <-
  sample_admb(m, iter=tt*iter, init=inits, thin=tt,
              parallel=TRUE, chains=chains, warmup=tt*warmup,
              path=d, cores=chains, control=list(metric=NULL),
              algorithm='RWM')

post.nuts1 <- cbind(alg='NUTS, 0.8', extract_samples(fit.nuts1))
post.nuts2 <- cbind(alg='NUTS, 0.9', extract_samples(fit.nuts2))
post.nuts3 <- cbind(alg='NUTS, 0.98', extract_samples(fit.nuts3))
post.rwm <- cbind(alg='RWM', extract_samples(fit.rwm))
df <- rbind(post.nuts1, post.nuts2, post.nuts3, post.rwm)
## Trick to get ggplot to read par names with brackets
p1="recdev2[22]"
p2="recdev2[23]"
df$r1 <- df[,p1]
df$r2 <- df[,p2]
g <- ggplot(df, aes(x=r1, y=r2)) + geom_point(alpha=.5, size=.1) +
  facet_wrap('alg', nrow=2)
ggsave("plots/hake_tests_recdevs.png", g, width=7, height=5)

## Save fits
saveRDS(list(fit.nuts1=fit.nuts1, fit.nuts2=fit.nuts2, fit.nuts3=fit.nuts3,
             fit.rwm=fit.rwm),
             file=paste0('results/hake_test_fits.RDS'))

## pairs_admb(fit=fit.rwm, pars=grep('recdev_early', names(post.nuts1))[1:10])
