library(r4ss)
library(R2admb)
library(adnuts)
library(shinystan)
library(rstan)
library(plyr)
library(snowfall)
library(vioplot)
library(reshape2)


## Get initial values for the pilot chains instead of starting from the
## MLE, using a previously run model and then saving these to file so they
## are reproducible.
## pilot.inits <- list(halibut=get.inits('halibut', 50, seed=1),
##               hake=get.inits('hake', 50, seed=1),
##               snowcrab=get.inits('snowcrab', 50, seed=1),
##               canary=get.inits('canary', 50, seed=1),
##               halibut2=get.inits('halibut2', 50, seed=1),
##               ## hake2=get.inits('hake2', 50, seed=1),
##               snowcrab2=get.inits('snowcrab2', 50, seed=1),
##               canary2=get.inits('canary2', 50, seed=1))
## saveRDS(pilot.inits, 'results/pilot.inits.RDS')
pilot.inits <- readRDS('results/pilot.inits.RDS')

plot.marginal <- function(fit, nrow=5, ncol=5, save=FALSE){
  post <- extract_samples(fit=fit)
  stds <- fit$mle$se[1:ncol(post)]
  mles <- fit$mle$est[1:ncol(post)]
  par.names <- names(post)
  names(mles) <- names(stds) <- par.names
  if(save)
    ## png(paste0('plots/marginal_fits_', fit$model,'_%02d.png'), units='in', res=500,
    ##     width=7,height=5)
    pdf(paste0('plots/marginal_fits_', fit$model,'.pdf'),
        width=7,height=5)
  par(mfrow=c(nrow,ncol), mar=c(1.5,0,0,0), mgp=c(1,.01, 0), tck=-.01)
  for(pp in par.names){
    temp <- hist(post[,pp], plot=FALSE)
    std <- stds[pp]
    est <- mles[pp]
    x <- seq(-3*std+est, to=3*std+est, len=1000)
    y <- dnorm(x, mean=est, sd=std)
    xlim <- range(c(x, temp$mids))
    ylim <- c(0, 1.3*max(c(y, temp$density)))
    plot(x,y, xlim=xlim, ylim=ylim, type='n', ann=FALSE, axes=FALSE)
    hist(post[,pp], col=gray(.8), plot=TRUE, add=TRUE, freq=FALSE)
    lines(x,y, lwd=2, col='red')
    axis(1, col=gray(.5))
    box(col=gray(.5))
    mtext(pp, line=-1.5)
  }
  if(save) dev.off()
}

plot.slow <- function(fit, n.slow=10, pars=NULL, fast=FALSE, save=TRUE){
  if(is.null(pars)){
    ## if pars nott specified use slowest mixing parameters
    ess <- monitor(fit$samples, warmup=fit$warmup, print=FALSE)[,'n_eff']
    pars <- names(sort(ess, decreasing=fast))[1:n.slow]
  }
  if(save)
    png(paste0('plots/pairs.', fit$model,'.fit.slow.png'), width=7, height=5, units='in', res=500)
  pairs_admb(fit, pars=pars, diag='trace')
  if(save)
    dev.off()
}

plot.sds <- function(fit){
  m <- fit$model
  sd.post <- apply(extract_samples(fit),2, sd)
  sd.mle <- fit$mle$se[1:fit$mle$nopar]
  sds <- data.frame(post=sd.post, mle=sd.mle)
  g <- ggplot(sds, aes(x=log10(mle), log10(post))) + geom_point(alpha=.7) +
    geom_abline(slope=1) + xlab("log10 MLE SE") + ylab("log10 Posterior SD")
  ggsave(paste0('plots/sds.', m, '.png'), g, width=7, height=5)
}

plot.uncertainties <- function(regularized, original=NULL, xlims, ylims){
  fit2 <- original; fit1 <- regularized
  n <- NROW(fit1$dq)
  png(paste0('plots/uncertainties_', fit1$model, '.png'), units='in', width=7,
             height=3, res=300)
  par(mfrow=c(1,n), mar=c(5,2,1,1 ), oma=c(0,0, 0, 0))
  for(i in 1:n){
    ii <- as.character(fit1$dq$dq[i])
    xx <- fit1$dq.post[,ii]
    hist(xx, freq=FALSE, xlim=xlims[[i]], ylim=ylims[[i]], col=gray(.8),
         border=gray(.8), breaks=50, xlab=ii, main=NA, ylab=NA)
    abline(v=mean(xx), col=2)
    lines(x <- seq(min(xlims[[i]]), max(xlims[[i]]), len=1000),
          y=dnorm(x, fit1$dq[i, 'mle'], fit1$dq[i,'se']))
    abline(v=fit1$dq[i, 'mle'], col=1)
    if(!is.null(fit2)){
      lines(x <- seq(min(xlims[[i]]), max(xlims[[i]]), len=1000),
            y=dnorm(x, fit2$dq[i, 'mle'], fit2$dq[i,'se']), col='blue')
      abline(v=fit2$dq[i, 'mle'], col='blue')
    }
  }
  dev.off()
}
plot.ess <- function(rwm, nuts){
  model <- rwm$model
  x <- monitor(rwm$samples, print=FALSE)[,'n_eff']/sum(rwm$time.total);
  y <-  monitor(nuts$samples, print=FALSE)[,'n_eff']/sum(nuts$time.total)
  temp <- range(c(x, y,0))
  ## png(paste0('plots/ess_comparison_',model, '.png'), width=7, height=3,
  ##     units='in', res=500)
  par(mfrow=c(1,3))
  plot(x=x, y=y, xlim=temp, ylim=temp, xlab='RWM', ylab='NUTS')
  abline(0,1)
  col1 <- gray(.7)
  barplot(sort(x), ylim=temp, main='RWM', col=col1, border=col1)
  barplot(sort(y), ylim=temp, main='NUTS', col=col1, border=col1)
#  dev.off()
}
plot.improvement <- function(fit1, fit2){
  ## xx <- rbind(data.frame(model=fit1$model, ess=fit1$ess,
  ##                  perf=fit1$ess/sum(fit1$time.total)),
  ##       data.frame(model=fit2$model, ess=fit2$ess,
  ##                  perf=fit2$ess/sum(fit2$time.total)) )
  ## xx$par <- row.names(xx)
  ## levels(xx$model) <- c("Original", "Fixed")
  ## ggplot(xx, aes(x=model, y=ess)) +geom_violin() + scale_y_log10()
  x <- monitor(fit1$samples, print=FALSE)[,'n_eff'];
  y <- monitor(fit2$samples, print=FALSE)[,'n_eff'];
  png(paste0('plots/ess_improvement_',fit1$model, '.png'), width=3, height=5,
      units='in', res=500)
  vioplot(log10(x), log10(y), names=c("Original", "Fixed"))
  mtext(fit1$model, line=1, cex=1.5)
  mtext("log10(ESS)", side=2, line=2.5, cex=1.25)
  dev.off()
}

get.inits <- function(model, chains, seed){
  set.seed(seed)
  file <- file.path(getwd(),'results', paste0('pilot_',model, '.RDS'))
  if(!file.exists(file)) {
    warning("pilot file not found for inits so using NULL")
    print(file)
    return(NULL)
  }
  fit <- readRDS(file)
  post <- extract_samples(fit)
  ## select random rows
  ind <- sample(1:nrow(post), size=chains)
  inits <- lapply(ind, function(i) as.numeric(post[i,]))
  return(inits)
}


check.identifiable <- function(path, model=path){
  ## Check eigendecomposition
  fit <- R2admb::read_pars(file.path(path, model))
  hes <- r4ss::getADMBHessian(file.path(path,'admodel.hes'),NULL)$hes
  ev  <-  eigen(hes)
  WhichBad <-  which( ev$values < sqrt(.Machine$double.eps) )
  if(length(WhichBad)==0){
    message( "All parameters are identifiable" )
    return(NULL)
  }
  ## Check for parameters
  RowMax  <-  apply(ev$vectors[, WhichBad], MARGIN=1, FUN=function(vec){max(abs(vec))} )
  bad <- data.frame(ParNum=1:nrow(hes), Param=names(fit$coefficients)[1:nrow(hes)], MLE=fit$coefficients[1:nrow(hes)], Param_check=ifelse(RowMax>0.1, "Bad","OK"))
  row.names(bad) <- NULL
  bad <- subset(bad, Param_check=='Bad')
  print(bad)
  return(invisible(bad))
}
