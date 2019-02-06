
## Run pilot chains for individual models using RWM, and then the
## "regularized" version of the model which has "2" appended. Each model
## saves a few key management parameters in MLE and MCMC estimates are
## added to the saved file later, after running NUTS

## This file assumes the models have already been run and produced .hes
## files (which contain the MLEs from which the chains start). Note that SS
## turns off the bias adjustment for recdevs during MCMC which can cause
## big problems if starting from the MLE. So make sure to optimize the
## model with the -mcmc flag so that SS turns it off for both. Also make
## sure to reoptimize without the flag for grabbing the MLE derived
## quantities and parameters.

## The inits are taken from a previous run and saved to file so that each
## model starts from a reproducible initial value that is not the MLE. See
## startup.R for how this was done.

for(m in c('hake', 'hake2')[1]){
setwd(m); system(paste(m,"-mcmc 100 -nox")); setwd('..')
inits <- pilot.inits[[m]][1:reps]
fit.rwm <-
  sample_admb(model=m, iter=iter, thin=thin, seeds=seeds, init=inits,
              parallel=TRUE, chains=reps, warmup=warmup,
              path=m, cores=reps, algorithm='RWM')
## Calculate and save monitor information from rstan
fit.rwm <- add.monitor(fit.rwm)
## Need to rerun model so r4ss works right
setwd(m);system(m);  setwd('..')
## Get estimates for derived quantitiesd
fit.rwm$dq <- get.dq(m)
saveRDS(fit.rwm, file=paste0("results/pilot_", m, ".RDS"))
}

for(m in c('halibut', 'halibut2')){
setwd(m); system(paste(m,"-mcmc 100 -nox")); setwd('..')
inits <- pilot.inits[[m]][1:reps]
fit.rwm <-
  sample_admb(m, iter=iter, thin=thin, seeds=seeds, init=inits,
              parallel=TRUE, chains=reps, warmup=warmup,
              path=m, cores=reps, algorithm='RWM')
fit.rwm <- add.monitor(fit.rwm)
## Need to rerun model so r4ss works right
setwd(m);system(m); setwd('..')
fit.rwm$dq <- get.dq(m)
saveRDS(fit.rwm, file=paste0("results/pilot_", m, ".RDS"))
}

## this is the only non-SS model so ignore the advice about the -mcmc flag above
for(m in c('snowcrab', 'snowcrab2')){
setwd(m); system(paste0(m, ' -nox -phase 50 -ainp ', m,'.par')); setwd('..')
inits <- pilot.inits[[m]][1:reps]
fit.rwm <-
  sample_admb(m, iter=iter, thin=thin, seeds=seeds, init=inits,
              parallel=TRUE, chains=reps, warmup=warmup,
              mceval=FALSE,
              path=m, cores=reps, algorithm='RWM')
fit.rwm <- add.monitor(fit.rwm)
fit.rwm$dq <- get.dq(m)
saveRDS(fit.rwm, file=paste0("results/pilot_", m, ".RDS"))
}

## this is by far the slowest model
for(m in c('canary', 'canary2')){
setwd(m); system(paste(m,"-mcmc 100 -nox")); setwd('..')
inits <- pilot.inits[[m]][1:reps]
fit.rwm <-
  sample_admb(m, iter=iter, thin=thin, seeds=seeds, init=inits,
              parallel=TRUE, chains=reps, warmup=warmup,
              path=m, cores=reps, algorithm='RWM')
fit.rwm <- add.monitor(fit.rwm)
## Need to rerun model so r4ss works right
setwd(m);system(m); setwd('..')
fit.rwm$dq <- get.dq(m)
saveRDS(fit.rwm, file=paste0("results/pilot_", m, ".RDS"))
}
