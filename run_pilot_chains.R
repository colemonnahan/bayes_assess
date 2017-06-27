## Run pilot chains for individual models using RWM, and then the "fixed"
## version of the model
## devtools::install('C:/Users/Cole/adnuts')
library(shinystan)
library(adnuts)
library(snowfall)
library(r4ss)
reps <- 5 # chains to run in parallel

sfStop()
d <- m <- 'cod'
d <- m <- 'cod2'
thin <- 100
iter <- 1000
warmup <- iter/4
inits <- NULL
sfInit(parallel=TRUE, cpus=reps)
sfExportAll()
fit.rwm <-
  sample_admb(m, iter=iter*thin, init=inits, thin=thin, mceval=TRUE,
              parallel=TRUE, chains=reps, warmup=warmup*thin,
              dir=d, cores=reps, algorithm='RWM')
## Get posterior draws of dqs to cbind onto parameter draws later
dq.names <- c("SSB_MSY", "SPB_50", "Bratio_50")
fit.rwm$dq.post <- r4ss::SSgetMCMC(dir=m)[[1]] [,dq.names]
xx <- SS_output(m, model=m, verbose=TRUE, covar=T)
## Get estimates for derived quantitiesd
dq <- subset(xx$derived_quants, LABEL %in% dq.names)[,1:3]
names(dq) <- c('dq','mle', 'se'); rownames(dq) <- NULL
fit.rwm$dq <- dq
saveRDS(fit.rwm, file=paste0("results/pilot_rwm_", m, ".RDS"))


d <- m <- 'hake'
thin <- 1000
iter <- 1000
warmup <- iter/4
inits <- NULL
sfStop()
sfInit(parallel=TRUE, cpus=reps)
sfExportAll()
fit.rwm <- sample_admb(m, iter=iter*thin, init=inits,  thin=thin,
              parallel=TRUE, chains=reps, warmup=warmup*thin, mceval=TRUE,
              dir=d, cores=reps, algorithm='RWM')
## Get posterior draws of dqs to cbind onto parameter draws later
dq.names <- c("SSB_MSY", "SPB_2013", "Bratio_2013")
fit.rwm$dq.post <- r4ss::SSgetMCMC(dir=m)[[1]][,dq.names]
xx <- SS_output(m, model=m, verbose=F, covar=T)
## Get estimates for derived quantitiesd
dq <- subset(xx$derived_quants, LABEL %in% dq.names)[,1:3]
names(dq) <- c('dq','mle', 'se'); rownames(dq) <- NULL
fit.rwm$dq <- dq
saveRDS(fit.rwm, file=paste0("results/pilot_rwm_", m, ".RDS"))


sfStop()
d <- m <- 'tanner2'
d <- m <- 'tanner'
thin <- 100
iter <- 1000
warmup <- iter/4
inits <- NULL
sfInit(parallel=TRUE, cpus=reps)
sfExportAll()
fit.rwm <- sample_admb(m, iter=iter*thin, init=inits, thin=thin,
             parallel=TRUE, chains=reps, warmup=warmup*thin, mceval=TRUE,
              dir=d, cores=reps, algorithm='RWM')
## Get posterior draws of dqs to cbind onto parameter draws later
dq.names <- c("SSB_2015", "depletion_2015", "OFL_main")
post <- read.csv(file.path(m, 'posterior.csv'), header=FALSE)
names(post) <- dq.names
fit.rwm$dq.post <- post
## Get estimates for derived quantitiesd
ind <- which(fit.rwm$mle$names.all %in% dq.names)
fit.rwm$dq <- data.frame(dq=fit.rwm$mle$names.all[ind], mle=fit.rwm$mle$est[ind], se=fit.rwm$mle$se[ind])
saveRDS(fit.rwm, file=paste0("results/pilot_rwm_", m, ".RDS"))


sfStop()
d <- m <- 'halibut'
d <- m <- 'halibut3'
d <- m <- 'halibut2'
thin <- 100
iter <- 1000
warmup <- iter/4
mle <- r4ss::read.admbFit(paste0(d,'/',m))
N <- mle$nopar
inits <- lapply(1:reps, function(i) mle$est[1:N])
sfInit(parallel=TRUE, cpus=reps)
sfExportAll()
fit.rwm <- sample_admb(m, iter=iter*thin, init=inits, thin=thin,
              parallel=TRUE, chains=reps, warmup=warmup*thin, mceval=TRUE,
              dir=d, cores=reps, algorithm='RWM')
## Get posterior draws of dqs to cbind onto parameter draws later
dq.names <- c("SSB_MSY", "SPB_2015", "Bratio_2015")
fit.rwm$dq.post <- r4ss::SSgetMCMC(dir=m)[[1]][,dq.names]
xx <- SS_output(m, model=m, verbose=TRUE, covar=T)
## Get estimates for derived quantitiesd
dq <- subset(xx$derived_quants, LABEL %in% dq.names)[,1:3]
names(dq) <- c('dq','mle', 'se'); rownames(dq) <- NULL
fit.rwm$dq <- dq
saveRDS(fit.rwm, file=paste0("results/pilot_rwm_", m, ".RDS"))


sfStop()
d <- 'snowcrab'; m <- 'snowcrab'
d <- 'snowcrab2'; m <- 'snowcrab2'
thin <- 100
iter <- 1000
warmup <- iter/4
inits <- NULL
sfInit(parallel=TRUE, cpus=reps)
sfExportAll()
fit.rwm <- sample_admb(m, iter=iter*thin, init=inits,  thin=thin,
              parallel=TRUE, chains=reps, warmup=warmup*thin, mceval=TRUE,
              dir=d, cores=reps, algorithm='RWM')
## Get posterior draws of dqs to cbind onto parameter draws later
dq.names <- c("SSB_2015", "depletion_2015", "OFL_main")
fit.rwm$dq.post <- read.csv(file.path(m, "posterior.csv"))
## Get estimates for derived quantities
xx <- read_mle_fit(m,d); ind <- which(xx$names.all %in% dq.names)
fit.rwm$dq <- data.frame(dq=dq.names, mle=xx$est[ind], se=xx$se[ind])
saveRDS(fit.rwm, file=paste0("results/pilot_rwm_", m, ".RDS"))


sfStop()
d <- m <- 'canary'
d <- m <- 'canary2'
thin <- 100
iter <- 1000
warmup <- iter/4
inits <- NULL
sfInit(parallel=TRUE, cpus=reps)
sfExportAll()
fit.rwm <- sample_admb(m, iter=iter*thin, init=inits, thin=thin,
              parallel=TRUE, chains=reps, warmup=warmup*thin, mceval=TRUE,
              dir=d, cores=reps, algorithm='RWM')
## Get posterior draws of dqs to cbind onto parameter draws later
dq.names <- c("SSB_MSY", "OFLCatch_2015", "Bratio_2015")
fit.rwm$dq.post <- r4ss::SSgetMCMC(dir=m)[[1]][,dq.names]
xx <- SS_output(m, model=m, verbose=TRUE, covar=T, ncols=500)
## Get estimates for derived quantitiesd
dq <- subset(xx$derived_quants, LABEL %in% dq.names)[,1:3]
names(dq) <- c('dq','mle', 'se'); rownames(dq) <- NULL
fit.rwm$dq <- dq
saveRDS(fit.rwm, file=paste0("results/pilot_rwm_", m, ".RDS"))

launch_shinyadmb(fit.rwm)

