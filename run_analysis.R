#### This file runs the analysis for the paper:

## Monnahan et al. (in review). Overcoming long run times in Bayesian stock
## assessments. ICES JMS.

## You can recreate the analysis by installing necessary packages and
## running the code below. You need to obtain .exe files directly from the
## author (monnahc@uw.edu) as the .tpl files are proprietary.

source('startup.R')
packageVersion('adnuts')
packageVersion('r4ss')
reps <- 5 # chains to run in parallel
set.seed(352)
seeds <- sample(1:1e4, size=reps)

## The pilot runs
thin <- 200
iter <- 2000
warmup <- iter/4
source('run_pilot_chains.R')

## Now do the comparisons of RWM vs NUTS and the different mass matrix
## versions on the regularized models. I'm using a smaller warmup since
## we're not using mass matrix adaptation, and starting from the mode so
## the chains need less warmup time. adapt_delta (ad) and max_treedepth
## (td) are set at defaults unless needed to be tweaked. inits is from the
## previous run.
seed <- 2352 # used to draw seeds in the template file
iter <- 2000; warmup <- (iter/5)
m <- 'halibut2'
ad <- .9; td <- 12
source('template.R')
m <- 'hake'
ad <- .9; td <- 12
source('template.R')
m <- 'snowcrab2'
ad <- .8; td <- 12
source('template.R')
m <- 'canary2' # this one takes way longer
ad <- .92; td <- 12
source('template.R')


## If desired, run the bias checks which looks at the recdevs for the
## halibut model. Compares RWM vs NUTS and explores adapt_delta.
source("startup.R")
iter <- 4000
warmup <- iter/2
chains <- 5
td <- 15
inits <- NULL # start from MLE
seeds <- 1:chains

## Now run the analysis script
source("analysis.R")

