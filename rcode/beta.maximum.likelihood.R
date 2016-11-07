### R code from vignette source 'likelihood.rnw'


###################################################
### code chunk number 9: closed-sir-model-defn
###################################################
require(deSolve)

SIRfunc=function(t, x, vparameters){
  S = x[1]  
  I = x[2]  
  R = x[3]  
  if (I<0) I=0 
  
  with(as.list(vparameters),{
    npop = S+I+R   
    dS = -beta*S*I/npop            
    dI = +beta*S*I/npop - gamma*I  
    dR = +gamma*I                  
    out = c(dS,dI,dR)
    list(out)
  })
}

###################################################
### code chunk number 10: closed-sir-predictions
###################################################

prediction <- function (inits, vparameters, times) {
  sirmodel = as.data.frame(lsoda(inits, times, SIRfunc, vparameters))
  sirmodel[-1,3] # return the I variable only
}


###################################################
### code chunk number 11: closed-sir-negloglik
###################################################
loglik <- function (params, data) {
  # browser()
  times <- data$vtime
  vparameters <- params[c("gamma", "beta")]
  inits <- unname(params[c("S_0", "I_0", "R_0")])
  pred <- prediction(inits, vparameters, times)
  sum(dnorm(x=data$vI,mean=pred,sd=params["sigma"],log=TRUE))
}

f <- function (beta, data, params) {
  #browser()
  params["beta"] <- beta
  loglik(params,data)
}


beta.hat.estimate <- daply(.data = newmodel, .variables = "iter", function(x) {
  ll <- sapply(beta,f, x)
  beta.hat <- beta[which.max(ll)]
  return(beta.hat)
})


trajectories <- adply(.data = data.files.trial, .margins = 1,.id=NULL, .expand = F, function(x) {
  load(as.character(x$file.list))
  trial.params <- get_params(x)
  #epidemic.index <- x$epidemic.trials[[1]]$epidemic.trials
  #if (epidemic.index[1] > 0) {
  #time.records.a <- align_time_series_all(time.records = time.records[epidemic.index])
  vI.trajectory <- combine_time_records(time.records)
  vI.trajectory <- cbind(trial.params, vI.trajectory)
  #vI.trajectory$iter = epidemic.index[vI.trajectory$iter]
  return(vI.trajectory)
  #}
})


trajectories %>% group_by(rnott, pop.size)  %>% 
  filter(rnott == 1 & pop.size == 1000) -> trial.data

sample.iter = sample(1:100, 20)

trajectories %>% filter(pop.size == 10000 & iter %in% sample.iter) -> trial.10000


beta.hat.estimate.10000 <- daply(.data = trial.10000,
                                .variables = c("rnott", "pop.size", "iter"), function(x) {
  N = unique(x$pop.size) 
  I_0 = .01 * N
  r0.value = unique(x$rnott)
  gamma = 1/3
  beta.expected = r0.value * gamma
  beta <- seq(from= (beta.expected - .3),to=(beta.expected + .3), by=0.001)
  params <- c(S_0=N-I_0, I_0=I_0, R_0 = 0, gamma = 1/3, beta=NA, sigma=1)
  ll <- sapply(beta,f,x, params)
  beta.hat <- beta[which.max(ll)]
  print(paste0(r0.value, "_iter_", unique(x$iter)))
  cbind(r0.value, beta.hat)
})


save(beta.hat.estimate.10000, file = "beta.hat.esimate.10000.RData")

beta.hat.estimate.10000 = data.frame(beta.hat.estimate.10000)
beta.hat.estimate.10000 = beta.hat.estimate.10000[,20:40]
colnames(beta.hat.estimate.10000)[1] = "rnott"
colnames(beta.hat.estimate.10000)[2:21] = seq(1:20)
beta.hat.estimate.10000 = gather(beta.hat.estimate.10000, iter, beta.hat, -1)

beta.hat.master = join(x = beta.hat.estimate.100, y = beta.hat.estimate.1000)

beta.hat.master = join(x = beta.hat.master, y = beta.hat.estimate.10000, by = "rnott")

beta.hat.master = beta.hat.master[, -2]


## Graphs
 
# Graphs based on final proportion

ggplot(final.size.test, aes(x = as.factor(rnott), y = effective.r0, 
                            color = as.factor(pop.size))) +
  geom_boxplot() + 
  scale_y_continuous(limits = c(0.9, 4), breaks = seq(0.9, 4, .2)) +
  
  ggplot(final.size.test, aes(prop.infected, fill = as.factor(pop.size)))+
  geom_density(alpha = .5)+facet_wrap(~rnott) + 
  geom_vline(aes(xintercept = final.size))


ggplot(beta.hat.master, aes(x = as.factor(rnott), y = rnott.hat, 
                            color = as.factor(beta.hat))) +
  geom_boxplot() + 
  scale_y_continuous(limits = c(0.9, 4), breaks = seq(0.9, 4, .2)) 
  