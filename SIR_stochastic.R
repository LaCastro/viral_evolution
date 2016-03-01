## Specify model parameters use within() to make assignments *inside* an 
## empty (or existing) list. Yhis is a handy R trick that allows you to 
## refer to existing list elements on right hand side (RHS)
##
## Note the braces, <-, and and no commas here:  everything in braces is a
## regular code block, except that assignments happen *inside* the list
params <- list()
params <- within(params, {
  
  ## set rng state
  seed <- 0
  tau <- 0.001 # in years
  nyears <- 10
  
  ## total number of steps
  nsteps <- nyears/tau
  
  mu <- 1/70 #death rate
  gamma <- 365/10 #recovery rate
  R0 <- 10
  ## refers to R0 above
  beta <- R0*(gamma+mu) #transmission rate
  nu <- mu #birth rate
  
  ## initial conditions, list within list
  ## use within() to modify empty list, as above
  init <- within(list(), {
    pop <- 1e6
    S <- round(pop/R0)
    I <- round(pop*mu*(1-1/R0)/(gamma+mu))
    ## refers to S,I above
    R <- pop-S-I
  })
})

set.seed(params$seed)

## run the model once
result.df <- tauleapCpp(params)

library(plyr)
nsim <- 12

## run many sims, combine all results into one data.frame
## plyr will combine results for us
result.rep <- ldply(1:nsim, function(.nn) {
  set.seed(.nn)
  ## run the model
  result <- tauleapCpp(params)
  ## this wastes space, but is very simple and aids plotting
  result$nsim <- .nn
  return(result)
})

library(lattice)

## lattice plot of results
plot(
  xyplot(I ~ time | sprintf("Simulation %02d",nsim), 
         data=result.rep, type=c('l','g'), as.table=T,
         ylab='Infected', xlab='Year',
         scales=list(y=list(alternating=F))
  )
)


