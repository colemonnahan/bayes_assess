## devtools::install('C:/Users/Cole/adnuts')
library(shinystan)
library(rstan)
library(adnuts)
library(snowfall)
recompile <- TRUE # recomile .tpl files?

reps <- 5                        # chains/reps to run
iter <- 2000; warmup <- (iter/10)
m <- d <- 'halibut3'
ad <- .8                                # adapt_delta
td <- 12
inits <- NULL
source('template.R')

reps <- 5                        # chains/reps to run
iter <- 2000; warmup <- (iter/10)
m <- d <- 'hake'
ad <- .8                                # adapt_delta
td <- 12
inits <- NULL
source('template.R')

reps <- 5                        # chains/reps to run
iter <- 500; warmup <- (iter/10)
m <- d <- 'canary2'
ad <- .8                                # adapt_delta
td <- 10
source('template.R')

reps <- 3                        # chains/reps to run
iter <- 1000; warmup <- (iter/10)
m <- d <- 'tanner2'
ad <- .8                                # adapt_delta
td <- 10
source('template.R')

## Old code to look at unbounded space
rotated <- read.csv(file.path(d,'rotated.csv'), head=FALSE)
unbounded <- read.csv(file.path(d,'unbounded.csv'), head=FALSE)
## Examine pairs in bounded space.
posterior <- extract_samples(fit.rwm.mle, inc_lp=TRUE)
ess <- fit.rwm.mle$ess
pars <- names(sort(ess))[1:8]
pairs_admb(posterior=posterior, mle=fit.rwm.mle$mle, pars=pars)
launch_shinyadmb(fit.rwm.mle)
## Now convert to unbounded and check
mle <- fit.rwm.mle$mle
unbounded$lp__ <- posterior$lp__
names(unbounded) <- names(posterior)
mle$cov <- fit.rwm.mle$covar.est
mle$se <- sqrt(diag(mle$cov))
mle$cor <- cov2cor(mle$cov)
mle$est <- apply(unbounded, 2, mean)
pairs_admb(posterior=unbounded, mle=mle, pars=pars, diag='trace')


