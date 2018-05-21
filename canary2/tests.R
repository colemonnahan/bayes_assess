## Some quick code to identify why the regularization changes biomass so
## much

library(r4ss)

dq.names <- c("SSB_MSY", "OFLCatch_2015", "Bratio_2015")
xx <- SS_output('../canary', model='canary', verbose=FALSE, covar=TRUE,
                ncols=300)
## Get estimates for derived quantitiesd
dq <- subset(xx$derived_quants, LABEL %in% dq.names)[,1:3]
names(dq) <- c('dq','mle', 'se'); rownames(dq) <- NULL

dq.original <- dq

getdq <- function(){
  m <- 'canary2'
  xx <- SS_output(getwd(), model=m, verbose=FALSE, covar=TRUE, ncols=300)
  ## Get estimates for derived quantitiesd
  dq <- subset(xx$derived_quants, LABEL %in% dq.names)[,1:3]
  names(dq) <- c('dq','mle', 'se'); rownames(dq) <- NULL
  return(dq)
}

## This is the base canary2 I left off with.. seems too big of changes
dq0 <- getdq()

## turn back up the sd devs from .25 back to .5
system('canary2')
dq1 <- getdq()
## turn back on a bunch of stuff and turn off priors
system('canary2')
dq2 <- getdq() ## MSY has no SE??
## undo selparm[10] and selparm[12] changes
system('canary2')
dq3 <- getdq()
## no ch ange so putting back on, also udno changes to selparm[30]
system('canary2')
dq3 <- getdq()

