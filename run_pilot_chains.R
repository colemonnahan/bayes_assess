## Run pilot chains for individual models using RWM, and then the "fixed"
## version of the model. Each model saves a few key management parameters
## in MLE and MCMC estimates, for both versions of models if they exist.

## This file assumes the models have already been run and produced .hes
## files (which contain the MLEs from which the chains start). Note that SS
## turns off the bias adjustment for recdevs during MCMC which can cause
## big problems if starting from the MLE. So make sure to optimize the
## model with the -mcmc flag so that SS turns it off for both.

source('startup.R')
reps <- 5 # chains to run in parallel
set.seed(352)
seeds <- sample(1:1e4, size=reps)

m <- 'hake'
setwd(m); system(paste(m,"-mcmc 100 -nox")); setwd('..')
thin <- 100
iter <- 2000
warmup <- iter/4
inits <- get.inits(m, reps, seed=12)
fit.rwm <-
  sample_admb(model=m, iter=iter*thin, thin=thin, seeds=seeds, init=inits,
              parallel=TRUE, chains=reps, warmup=warmup*thin,
              path=m, cores=reps, algorithm='RWM')
## Get posterior draws of dqs to cbind onto parameter draws later. Need to
## rerun model so r4ss works right
setwd(m);system(m); system(paste(m, '-mceval')); setwd('..')
dq.names <- c("SSB_MSY", "SPB_2013", "Bratio_2013")
fit.rwm$dq.post <- r4ss::SSgetMCMC(dir=m)[[1]][,dq.names]
xx <- SS_output(m, model=m, verbose=F, covar=TRUE)
## Get estimates for derived quantitiesd
dq <- subset(xx$derived_quants, LABEL %in% dq.names)[,1:3]
names(dq) <- c('dq','mle', 'se'); rownames(dq) <- NULL
fit.rwm$dq <- dq
saveRDS(fit.rwm, file=paste0("results/pilot_rwm_", m, ".RDS"))

m <- 'halibut'
m <- 'halibut2'
setwd(m); system(paste(m,"-mcmc 100 -nox")); setwd('..')
thin <- 100
iter <- 2000
warmup <- iter/4
inits <- get.inits(m, reps, seed=121)
fit.rwm <-
  sample_admb(m, iter=iter*thin, thin=thin, seeds=seeds, init=inits,
              parallel=TRUE, chains=reps, warmup=warmup*thin,
              path=m, cores=reps, algorithm='RWM')
## Get posterior draws of dqs to cbind onto parameter draws later. Need to
## rerun model so r4ss works right
setwd(m);system(m); system(paste(m, '-mceval')); setwd('..')
dq.names <- c("SPB_2000", "SPB_2010", "SPB_2015")
fit.rwm$dq.post <- r4ss::SSgetMCMC(dir=m)[[1]][,dq.names]
xx <- SS_output(m, model=m, verbose=FALSE, covar=TRUE)
## Get estimates for derived quantitiesd
dq <- subset(xx$derived_quants, LABEL %in% dq.names)[,1:3]
names(dq) <- c('dq','mle', 'se'); rownames(dq) <- NULL
fit.rwm$dq <- dq
saveRDS(fit.rwm, file=paste0("results/pilot_rwm_", m, ".RDS"))

m <- 'canary'
m <- 'canary2'
setwd(m); system(paste(m,"-mcmc 100 -nox")); setwd('..')
thin <- 100
iter <- 2000
warmup <- iter/4
inits <- get.inits(m, reps, seed=12)
fit.rwm <-
  sample_admb(m, iter=iter*thin, thin=thin, seeds=seeds, init=inits,
              parallel=TRUE, chains=reps, warmup=warmup*thin,
              path=m, cores=reps, algorithm='RWM')
## Get posterior draws of dqs to cbind onto parameter draws later
setwd(m);system(m);
##system(paste(m, '-mceval'));
setwd('..')
dq.names <- c("SSB_MSY", "OFLCatch_2015", "Bratio_2015")
fit.rwm$dq.post <- r4ss::SSgetMCMC(dir=m)[[1]][,dq.names]
xx <- SS_output(m, model=m, verbose=FALSE, covar=T, ncols=500)
## Get estimates for derived quantitiesd
dq <- subset(xx$derived_quants, LABEL %in% dq.names)[,1:3]
names(dq) <- c('dq','mle', 'se'); rownames(dq) <- NULL
fit.rwm$dq <- dq
saveRDS(fit.rwm, file=paste0("results/pilot_rwm_", m, ".RDS"))


m <- 'snowcrab';
m <- 'snowcrab2';
setwd(m); system(m); setwd('..')
thin <- 100
iter <- 2000
warmup <- iter/4
inits <- get.inits(m, reps, seed=678)
fit.rwm <-
  sample_admb(m, iter=iter*thin, thin=thin, seeds=seeds, init=inits,
              parallel=TRUE, chains=reps, warmup=warmup*thin,
              path=m, cores=reps, algorithm='RWM')
## Get posterior draws of dqs to cbind onto parameter draws later
dq.names <- c("SSB_2015", "F35sd", "OFL_main")
fit.rwm$dq.post <- read.csv(file.path(m, "posterior.csv"))
## Get estimates for derived quantities
xx <- R2admb::read_admb(file.path(m,m))
fit.rwm$dq <- data.frame(dq=dq.names, mle=xx$coefficients[dq.names], se=xx$se[dq.names])
saveRDS(fit.rwm, file=paste0("results/pilot_rwm_", m, ".RDS"))

