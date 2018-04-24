
# cod
inits <- lapply(1:reps, function(i)
  c(runif(1, 18, 20),
    runif(1, 100,130),
    runif(1, .2,.3),
    runif(1, .01, .2),
    runif(1, .01, .3),
    runif(1, 20.5, 21.5),
    runif(50, 0, 2),
    runif(1, -1, 0),
    runif(1, 40, 60),
    runif(1, 3, 6),
    runif(1, 20, 50),
    runif(1, 2, 6)))

## halibut
setwd('halibut3')
file.copy('halibut3.par', 'halibut3.parold', overwrite=TRUE)

ff <- function(n)
  boxplot(post[,n], sapply(1:1000, function(i) inits()[n]))
inits.fn <- function()
  c(runif(1, .1, .2),
    runif(1, 10, 11),
    runif(1, -.5,.5),
    runif(21, -3,3),
    runif(1, 12,15),
    runif(11,-3,3),
    runif(1, -7, -4), #37
    runif(18, -.1,.1), #55
    runif(1, 12,18),
    runif(1, 0,5),
    runif(1, -6,0), #58
    runif(1, -4,5),
    runif(1, 7, 14),
    runif(1, -4,5),
    runif(1, -4,5), #62
    runif(1, -4,5),
    runif(1, -4,5), #64
    runif(1, -4, 0),
    runif(1, .1,1),
    runif(1, 0,25), #67
    runif(1, 0,4),
    runif(1, 0, 3),
    runif(1, 0,1),
    runif(105, -1,1))
set.seed(2352)

for(i in 1:3){
  ## overwrite par file
  i0 <- c(0,inits.fn())
  write.table(x=i0, file='init.pin', row.names=FALSE, col.names=FALSE)
  ## run model
  system('halibut3 -nohess -ainp init.pin')
  out <- SSget
}


set.seed(325)
inits <- lapply(1:5, function(i) inits.fn())
post <- sample_admb('halibut3', iter=2000, algorithm='RWM', init=inits, chains=5)


## hake
inits <- NULL
inits <- lapply(1:reps, function(i)
  c(runif(1, .1, .4),
    runif(1, 15, 17),
    runif(1, .1, 1),
    runif(24, -4,4),
    runif(45, -4,4),
    runif(5, -6,6),
    runif(1, .03, 2),
    runif(9, -2, 6),
    runif(130, -1,1)))
