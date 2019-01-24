
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
## Need to rerun model so r4ss works right
setwd(m);system(m);  setwd('..')
dq.names <- c("SSB_MSY", "SPB_2013", "Bratio_2013")
fit.rwm$dq.post <- r4ss::SSgetMCMC(dir=m)[[1]][,dq.names]
xx <- SS_output(m, model=m, verbose=F, covar=TRUE)
## Get estimates for derived quantitiesd
dq <- subset(xx$derived_quants, LABEL %in% dq.names)[,1:3]
names(dq) <- c('dq','mle', 'se'); rownames(dq) <- NULL
fit.rwm$dq <- dq
saveRDS(fit.rwm, file=paste0("results/pilot_", m, ".RDS"))
}

for(m in c('halibut', 'halibut2')){
setwd(m); system(paste(m,"-mcmc 100 -nox")); setwd('..')
inits <- pilot.inits[[m]][1:reps]
fit.rwm <-
  sample_admb(m, iter=iter, thin=thin, seeds=seeds, init=inits,
              parallel=TRUE, chains=reps, warmup=warmup,
              path=m, cores=reps, algorithm='RWM')
## Need to rerun model so r4ss works right
setwd(m);system(m); setwd('..')
dq.names <- c("SPB_2000", "SPB_2010", "SPB_2015")
fit.rwm$dq.post <- r4ss::SSgetMCMC(dir=m)[[1]][,dq.names]
xx <- SS_output(m, model=m, verbose=FALSE, covar=TRUE)
dq <- subset(xx$derived_quants, LABEL %in% dq.names)[,1:3]
names(dq) <- c('dq','mle', 'se'); rownames(dq) <- NULL
fit.rwm$dq <- dq
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
## No need to rerun it, the DQs are ready
dq.names <- c("SSB_2015", "F35sd", "OFL_main")
fit.rwm$dq.post <- read.csv(file.path(m, "posterior.csv"))
xx <- R2admb::read_admb(file.path(m,m))
fit.rwm$dq <- data.frame(dq=dq.names, mle=xx$coefficients[dq.names], se=xx$se[dq.names])
saveRDS(fit.rwm, file=paste0("results/pilot_", m, ".RDS"))
}

## this is the only non-SS model so ignore the advice about the -mcmc flag above
for(m in c('pollock', 'pollock2')[1]){
setwd(m); system(paste0(m, ' -nox -phase 50 -ainp ', m,'.par')); setwd('..')
inits <- pilot.inits[[m]][1:reps]
fit.rwm <-
  sample_admb(m, iter=iter, thin=thin, seeds=seeds, init=inits,
              parallel=TRUE, chains=reps, warmup=warmup,
              path=m, cores=reps, algorithm='RWM')
## No need to rerun it, the DQs are ready
## dq.names <- c("SSB_2015", "F35sd", "OFL_main")
## fit.rwm$dq.post <- read.csv(file.path(m, "posterior.csv"))
## xx <- R2admb::read_admb(file.path(m,m))
## fit.rwm$dq <- data.frame(dq=dq.names, mle=xx$coefficients[dq.names], se=xx$se[dq.names])
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
## Need to rerun model so r4ss works right
setwd(m);system(m); setwd('..')
dq.names <- c("SSB_MSY", "OFLCatch_2015", "Bratio_2015")
fit.rwm$dq.post <- r4ss::SSgetMCMC(dir=m)[[1]][,dq.names]
xx <- SS_output(m, model=m, verbose=FALSE, covar=T, ncols=500)
dq <- subset(xx$derived_quants, LABEL %in% dq.names)[,1:3]
names(dq) <- c('dq','mle', 'se'); rownames(dq) <- NULL
fit.rwm$dq <- dq
saveRDS(fit.rwm, file=paste0("results/pilot_", m, ".RDS"))
}
