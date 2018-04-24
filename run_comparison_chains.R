## devtools::install('C:/Users/Cole/adnuts')
library(shinystan)
library(rstan)
library(adnuts)
library(snowfall)
recompile <- TRUE # recomile .tpl files?


## Investigate performance differences between algorithms and settings
reps <- 3                        # chains/reps to run
iter <- 2000; warmup <- (iter/2)
m <- d <- 'cod'
ad <- .8                                # adapt_delt
td <- 12
set.seed(252352)
#inits <- function() 128.966, 0.212482, 0.10938, 0.0980505, 19, -0.485899, -0.467345, 0.204968, 0.177974, -0.457728, -0.019547, -0.646809, -0.0316159, 0.247758, 0.307354, 0.433508, 0.312063, -0.5211, 0.139359, -0.276049, -0.0195735, 0.231638, -0.0404556, -0.35596, 0.462223, 0.135116, -0.408404, 0.344326, 0.313476, 0.0455946, 0.0060742, 0.00217732, -0.319812, 0.10312, -0.1992, -0.0980967, 0.146143, -0.0411011, -0.209175, 0.350999, -0.660133, -0.0507754, 0.184323, -0.0961166, 0.191711, -0.715587, 0.245363, -0.0715448, 0.741792, 0.294105, 0.0474952, 0.172519, 0.328332, 0.00412565, 0.0183907, -0.0253599, 52.6217, 5.26589, 37.3857, 4.4318)
#inits <- NULL
source('template.R')

reps <- 3                        # chains/reps to run
iter <- 1000; warmup <- (iter/10)
m <- d <- 'halibut3'
ad <- .8                                # adapt_delta
td <- 12
set.seed(23525)
source('template.R')

reps <- 4                        # chains/reps to run
iter <- 1000; warmup <- (iter/10)
m <- d <- 'hake'
ad <- .9                                # adapt_delta
td <- 10
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


