## Some quick code to identify why the regularization changes biomass so
## much

library(r4ss)

setwd('../canary')
system("canary -sdonly")
setwd('../canary2')

dq.names <- c("SSB_MSY", "OFLCatch_2015", "Bratio_2015")
xx <- SS_output('../canary', model='canary', verbose=FALSE, covar=TRUE,
                ncols=300)
## Get estimates for derived quantitiesd
dq <- subset(xx$derived_quants, LABEL %in% dq.names)[,1:3]
names(dq) <- c('dq','mle', 'se'); rownames(dq) <- NULL
dq.original <- dq
dq

getdq <- function(covar=TRUE){
  m <- 'canary2'
  xx <- SS_output(getwd(), model=m, verbose=FALSE, covar=TRUE, ncols=300)
  ## Get estimates for derived quantitiesd
  dq <- subset(xx$derived_quants, LABEL %in% dq.names)[,1:3]
  names(dq) <- c('dq','mle', 'se'); rownames(dq) <- NULL
  return(dq)
}

## canary2 with no control file changes; ensure matches canary
system('canary2 -nox')
dq0 <- getdq()

## regularize selparms 4,6, 10, 12
system('canary2')
dq1 <- getdq()


