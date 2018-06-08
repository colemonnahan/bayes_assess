## For each model look at how long it takes to calculate derivatives vs
## objective function.

## I think it's best to do this with mcdiag so it's not doing any of those
## calculations
source("startup.R")

compare.speed <- function(m, iter, tt=4){
  setwd(m)
  on.exit(setwd('..'))
  system(paste(m, '-hbf 0 -nox'))
  ts <- Sys.time()
  system(paste(m, '-mcmc', iter*tt,' -rwm -mcdiag -noest'))
  t1 <- Sys.time()-ts
  ts <- Sys.time()
  system(paste(m, '-mcmc', iter*tt,' -rwm -mcdiag -noest -mcdiag'))
  t4 <- Sys.time()-ts
  system(paste(m, '-hbf 1 -nox'))
  ts <- Sys.time()
  system(paste(m, '-mcmc', iter,' -hmc -hynstep 1 -hyeps .0001 -noest'))
  t2 <- Sys.time()-ts
  ts <- Sys.time()
  system(paste(m, '-mcmc', iter,' -nuts -max_treedepth 0 -hyeps .0001 -noest'))
  t3 <- Sys.time()-ts
  ts <- Sys.time()
  system(paste(m, '-mcmc', iter,' -hmc -hynstep 1 -hyeps .0001 -noest -mcdiag'))
  t5 <- Sys.time()-ts
  ts <- Sys.time()
  system(paste(m, '-mcmc', iter,' -nuts -max_treedepth 0 -hyeps .0001 -noest -mcdiag'))
  t6 <- Sys.time()-ts
  t1 <- 1000*as.double(t1, units='secs')
  t2 <- 1000*as.double(t2, units='secs')
  t3 <- 1000*as.double(t3, units='secs')
  t4 <- 1000*as.double(t4, units='secs')
  t5 <- 1000*as.double(t5, units='secs')
  t6 <- 1000*as.double(t6, units='secs')
  ## With dense matrix
  x1 <- round(c(rwm=t1/(iter*tt),hmc=t2/iter,nuts=t3/iter), 4)
  ## With mcdiag
  x2 <- round(c(rwm=t4/(iter*tt),hmc=t5/iter,nuts=t6/iter), 4)
  return(rbind(dense=x1,diag=x2))
}

hake <- compare.speed('hake', 2000)
halibut <- compare.speed('halibut', 2000)
canary <- compare.speed('canary', 2000)
