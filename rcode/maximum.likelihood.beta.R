
new.model.time <- time_records_all(newmodel)
new.model.time <- combine_time_records(new.model.time)


ggplot(new.model.time, aes(x = vtime, y = vI, group = iter, color = iter)) + geom_line()

require(deSolve)

sir.model.closed <- function (t, x, params) { #here we begin a function with three arguments
  S <- x[1] #create local variable S, the first element of x
  I <- x[2] #create local variable I
  R <- x[3] #create local variable R
  with( #we can simplify code using "with"
    as.list(params), #this argument to "with" lets us use the variable names
    { #the system of rate equations
      dS <- -beta*S*I
      dI <- beta*S*I-gamma*I
      dR <- gamma*I
      dx <- c(dS,dI,dR) #combine results into a single vector dx
      list(dx) #return result as a list
    }
  )
}

# Return a vector of the model predictions
vt = seq(tbeg,tend,delta_t)
vparameters = c(gamma=gamma,beta=beta)
inits = c(S=S_0,I=I_0,R=R_0)

sirmodel = as.data.frame(lsoda(inits, vt, SIRfunc, vparameters))

prediction <- function (inits, times, vparameters) {
  sirmodel = as.data.frame(lsoda(inits, vt, SIRfunc, vparameters))
  return(sirmodel[-1, 3])
  
  #xstart <- params[c("S.0","I.0", "R.0")]
  #browser()
  #out <- ode(
    
   # y=xstart,
  #  times=c(0,times),
  #  parms=params
  #)

  #out[-1,3] # return the I variable only
}

model.predict = prediction(inits = params[c("S_0", "I_0", "R_0")], 
                           times = times,
                           vparameters = params[c("gamma", "beta")])


loglik <- function (params, data) {
  #times <- data$biweek/26
  times <- data$vtime
  pred <- prediction(inits = params[c("S_0", "I_0", "R_0")], 
                     times = times,
                     vparameters = params[c("gamma", "beta")])
 # browser()
  sum(dnorm(x=data$vI,mean=pred,sd=params["sigma"],log=TRUE))
}

#dat <- data.frame(biweek=seq(1:dim(niamey)[1]), measles=niamey[,1])
new.model.time %>% filter(iter == 1) %>% select(vtime, vI) -> trial.1

#params <- c(S.0=1000,I.0=10,gamma=365/13 ,beta=NA,sigma=1)
params <- c(S_0 = 1000, I_0 = 10, R_0 = 0, gamma = 1/3, beta = NA, sigma = 1)

f <- function (beta, data) {
  par <- params
  #browser()
  par["beta"] <- beta
  loglik(par,data)
}

#beta <- seq(from=0,to=0.02,by=0.0001)
beta <- seq(from=.3, to=.8, by=0.01)


ll.beta <- daply(new.model.time, .variables = "iter", function(x) {  
  ll <- sapply(beta,f, x)
  beta.hat <- beta[which.max(ll)]
  return(beta.hat)
})

summary(ll.beta)
hist(ll.beta)
summary(ll.beta)/gamma

#plot(beta,-ll,type='l',ylab="-log(L)")
beta.hat <- beta[which.max(ll)]
abline(v=beta.hat,lty=2)

load('data.RData')
