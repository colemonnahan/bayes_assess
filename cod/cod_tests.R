library(r4ss)

## Verify MLE is good
d <- getwd()
system('cod')
out <- SS_output(dir=d, model='cod')
SS_plots(out, html=TRUE, uncertainty=TRUE)

## Run RWM algorithm. Note that burnin and thinning are done in the
## command line, first 500 are burnin second 500 are samples
system("cod -nox -mcmc 10000 -rwm -mcscale 5000 -chain 1 -mcseed 3229621 -mcsave 10")
system("cod -mceval")

mcmc <- SSgetMCMC(dir=d, writecsv=TRUE)
spb <- mcmc$model1[,grep('SPB_', x= names(mcmc$model1))]
## Plot the end year biomass vs negative log likelihood
plot(spb[,'SPB_50'], mcmc$model1$Objective_function)
## It is crashing!
plot(0,0, xlim=c(1,52), ylim=c(0, max(spb)*1.05), type='n')
trash <- apply(spb, 1, function(i)
  lines(1:52, y=i))

## The RWM spits out the current NLL during the -mcmc phase. So we can
## compare this to the mceval phase calculated by SS
lp <- read.table('rwm_lp.txt', header=TRUE)[,1]
plot(-lp, mcmc$model1$Objective_function, xlab='NLL during MCMC',
     ylab='NLL during mceval')
abline(0,1)

## These are very different. SS is calculating a different NLL during the
## MCMC phase than the mceval phase. This happens with all three
## algorithms, so I think it must be on the ss.tpl side of things.
