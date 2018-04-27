

source("startup.R")
iter <- 1000
warmup <- iter/4
td <- 12
ad <- .8

## Run model with hbf=1 to get right covariance matrix and MLEs
d <- m <- 'hake'
setwd(d)
system(paste('admb', m))
system(paste(m, '-hbf 1 -nox'))
setwd('..')

sfStop()
## Draw inits from MVT using MLE and covar
inits <- NULL
eps <- NULL
sfInit(parallel=TRUE, cpus=reps)
sfExportAll()

fit.temp <- readRDS('results/hake_fits.RDS')[[2]]
covar.dense <- fit.temp$covar.est
fit.nuts.dense <-
  sample_admb(m, iter=iter, init=inits, algorithm='NUTS',
               parallel=TRUE, chains=reps, warmup=warmup, dir=d, cores=reps,
              control=list(max_treedepth=td, stepsize=eps, extra.args='-noest',
                           metric=covar.dense, adapt_delta=ad))

## Now run RWM but using a thinning rate similar to NUTS so the time is
## roughly equivalent.
tt <- 4*floor(mean(extract_sampler_params( fit.nuts.dense)$n_leapfrog__))
## Rerun model with hbf=0
setwd(d)
system(paste(m, '-nox'))
setwd('..')
fit.rwm.mle <-
  sample_admb(m, iter=tt*iter, init=inits, thin=tt,
              parallel=TRUE, chains=reps, warmup=tt*warmup,
              dir=d, cores=reps, control=list(metric=NULL),
              algorithm='RWM')

nuts <- data.frame(alg='NUTS', par=names(fit.nuts.dense$Rhat),
                        Rhat=as.numeric(fit.nuts.dense$Rhat),
                        ess=as.numeric(fit.nuts.dense$ess),
                        efficiency=as.numeric(fit.nuts.dense$ess)/sum(fit.nuts.dense$time.total))
rwm <- data.frame(alg='RWM', par=names(fit.rwm.mle$Rhat),
                        Rhat=as.numeric(fit.rwm.mle$Rhat),
                        ess=as.numeric(fit.rwm.mle$ess),
                        efficiency=as.numeric(fit.rwm.mle$ess)/sum(fit.rwm.mle$time.total))
df <- rbind(nuts,rwm)
df.long <- reshape2::melt(df, c("alg", "par"))
g <- ggplot(df.long, aes(alg, value)) + geom_violin() + facet_wrap('variable', scales='free_y')
ggsave("plots/hake_tests_stats.png", g, width=7, height=5)

post.nuts <- cbind(alg='NUTS', extract_samples(fit.nuts.dense))
post.rwm <- cbind(alg='RWM', extract_samples(fit.rwm.mle))
df <- rbind(post.nuts, post.rwm)
p1="recdev1[22]"
p2="recdev1[23]"
df$r1 <- df[,p1]
df$r2 <- df[,p2]
g <- ggplot(df, aes(x=r1, y=r2)) + geom_point(alpha=.5, size=.1) +
  facet_wrap('alg', nrow=2)
ggsave("plots/hake_tests_recdevs.png", g, width=7, height=5)

## Save fits
saveRDS(list(fit.nuts.dense, fit.rwm.mle),
        file=paste0('results/hake_test_fits.RDS'))
