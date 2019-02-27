### This file demonstrates how to run Bayesian inference on ADMB stock
### assessments using the adnuts R package. We demonstrate a Stock
### Synthesis (SS) model called 'hake'. It is the supplemental material for
### the paper:

### Monnahan, C.C., T.A. Branch, J.T. Thorson, I.J Stewart, C.S. Szuwalski
### (2019). Overcoming long Bayesian run times in integrated fisheries stock
### assessments. In review at ICES Journal of Marine Science.

### The use of SS necessitates slightly different workflow for technical
### reasons. First, when optimizing before MCMC initiate -mcmc 50 to tell
### SS to turn off bias adjustment for recdevs. Otherwise the estimated
### mass matrix will be mismatched when executing the real MCMC chains. Be
### careful not to use MLE estimates from these runs for inference. To save
### time we recommend setting SS to read from the .par file to speed up the
### optimizations below.

### 2/2019 Cole Monnahan | monnahc@uw.edu

library(adnuts)
library(snowfall)
library(rstan)
library(shinystan)
reps <- parallel::detectCores()-1 # chains to run in parallel
## Reproducible seeds are passed to ADMB
set.seed(352)
seeds <- sample(1:1e4, size=reps)

## Here we assume the hake.exe model is in a folder called 'hake'
## as well. This folder gets copied during parallel runs.
m <- 'hake'
## First optimize the model to make sure the Hessian is good.
setwd(m); system('hake -nox -iprint 200 -mcmc 15'); setwd('..')

## Then run parallel RWM chains as a first test
thin <- 100
iter <- 1000*thin; warmup <- iter/4
inits <- NULL ## start chains from MLE
pilot <- sample_admb(m, iter=iter, thin=thin, seeds=seeds, init=inits,
                     parallel=TRUE, chains=reps, warmup=warmup,
                     path=m, cores=reps, algorithm='RWM')

## Check convergence
mon <- monitor(pilot$samples, warmup=pilot$warmup, print=FALSE)
max(mon[,'Rhat'])
min(mon[,'n_eff'])
## Examine the slowest mixing parameters
slow <- names(sort(mon[,'n_eff']))[1:8]
pairs_admb(fit=pilot, pars=slow)
## Or can specify them by name
pairs_admb(fit=pilot, pars=c('MGparm[1]', 'SR_parm[1]', 'SR_parm[2]'))

## After regularizing we can run NUTS chains. First reoptimize to get the
## correct mass matrix for NUTS. Note the -hbf 1 argument. This is a
## technical requirement b/c NUTS uses a different set of bounding
## functions and thus the mass matrix will be different.
setwd(m); system(paste(m, '-hbf 1 -nox -iprint 200 -mcmc 15')); setwd('..')
## Use default MLE covariance (mass matrix) and short parallel NUTS chains
## started from the MLE.
nuts.mle <-
  sample_admb(model=m, iter=500, init=NULL, algorithm='NUTS',  seeds=seeds,
               parallel=TRUE, chains=reps, warmup=100, path=m, cores=reps,
              control=list(metric="mle", adapt_delta=0.8))
## Check for issues like slow mixing, divergences, max treedepths with
## ShinyStan and pairs_admb as above. Fix and rerun this part as needed.
launch_shinyadmb(nuts.mle)

## If good, run again for inference using updated mass matrix. Increase
## adapt_delta toward 1 if you have divergences (runs will take longer).
mass <- nuts.mle$covar.est # note this is in unbounded parameter space
inits <- sample_inits(nuts.mle, reps) ## use inits from pilot run
nuts.updated <-
  sample_admb(model=m, iter=1000, init=inits, algorithm='NUTS',  seeds=seeds,
               parallel=TRUE, chains=reps, warmup=100, path=m, cores=reps,
              mceval=TRUE, control=list(metric=mass, adapt_delta=0.9))
## Again check for issues of nonconvergence and other standard checks. Then
## use for inference.
mon <- monitor(nuts.updated$samples, warmup=nuts.updated$warmup, print=FALSE)
max(mon[,'Rhat'])
min(mon[,'n_eff'])

## We can calculate efficiecy as ess/time. Since there's multiple chains
## add the time together because the ESS is summed across chains too.
(eff <- min(mon[,'n_eff'])/sum(nuts.updated$time.total))
## Or how long to get 1000 effective samples
1000/eff                                # in seconds
1000/eff/60                             # in minutes

## NOTE: the mceval=TRUE argument tells ADMB to run -mceval on ALL chains
## combined AFTER discarding warmup period and thinning. Thus whatever your
## model outputs during mceval is ready for use in
## management. Alternatively you can run -mceval from the command
## line. sample_admb will merge samples into the .psv file in the main
## folder so either way works.
