## This file runs longer chains for comparing efficiency (speed) between
## RWM and NUTS.
source("startup.R")

recompile <- TRUE # recomile .tpl files?

## I'm using a smaller warmup since we're not using mass matrix adaptation,
## and starting from the mode so the chains need less warmup time.

## adapt_delta (ad) and max_treedepth (td) are set at defaults unless
## needed to be tweaked. inits NULL means to use the mode.

reps <- 5                        # chains/reps to run
iter <- 2000; warmup <- (iter/10)
m <- 'halibut3'
ad <- .8; td <- 12
inits <- NULL
source('template.R')

reps <- 5                        # chains/reps to run
iter <- 2000; warmup <- (iter/10)
m <- 'hake'
ad <- .8; td <- 12
inits <- NULL
source('template.R')

reps <- 5                        # chains/reps to run
iter <- 2000; warmup <- (iter/10)
m <- 'canary2'
ad <- .8; td <- 12
source('template.R')

reps <- 5                        # chains/reps to run
iter <- 2000; warmup <- (iter/10)
m <- 'tanner2'
ad <- .8; td <- 10
source('template.R')
