## This file runs longer chains for comparing efficiency (speed) between
## RWM and NUTS.
source("startup.R")

recompile <- FALSE # recomile .tpl files?

## I'm using a smaller warmup since we're not using mass matrix adaptation,
## and starting from the mode so the chains need less warmup time.

## adapt_delta (ad) and max_treedepth (td) are set at defaults unless
## needed to be tweaked. inits NULL means to use the mode.

reps <- 5                        # chains/reps to run
iter <- 1000; warmup <- (iter/5)
m <- 'halibut2'
ad <- .9; td <- 10
inits <- get.inits(m, reps, 12)
source('template.R')

reps <- 5                        # chains/reps to run
iter <- 1000; warmup <- (iter/5)
m <- 'hake2'
ad <- .9; td <- 12
inits <- get.inits(m, reps, 12)
source('template.R')

reps <- 10                        # chains/reps to run
iter <- 500; warmup <- (iter/5)
m <- 'canary2'
ad <- .8; td <- 10
inits <- get.inits(m, reps, 12)
source('template.R')

reps <- 5                        # chains/reps to run
iter <- 1000; warmup <- (iter/4)
m <- 'snowcrab2'
ad <- .8; td <- 10
inits <- get.inits(m, reps, 12)
source('template.R')
