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


trajectories <- adply(.data = data.files, .margins = 1,.id=NULL, .expand = F, function(x) {
  load(as.character(x$file.trial))
  trial.params <- get_params(x)
  #epidemic.index <- x$epidemic.trials[[1]]$epidemic.trials
  #if (epidemic.index[1] > 0) {
  #time.records.a <- align_time_series_all(time.records = time.records[epidemic.index])
  time.records.a <- align_time_series_all(time.records = time.records)
  vI.trajectory <- combine_time_records(time.records.a)
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

final.size.test = join(final.infected.calculatedr0, deterministic.final.size, by = "rnott")
colnames(deterministic.final.size)[1] = "rnott"

final.prop <- ggplot(final.size.test, aes(x = as.factor(rnott), y = effective.r0, 
                            fill = as.factor(pop.size))) +
  geom_boxplot() + 
  scale_y_continuous(limits = c(0.9, 4), breaks = seq(0.9, 4, .2)) +
  scale_fill_discrete(name="Population Size",
                      breaks=c("beta.hat.100", "beta.hat.1000", "beta.hat.10000"),
                      labels=c("100", "1000", "10000")) +
  labs(x = "Simulation Beta/Gamma", y = "effective.rnott", title = "Calculated From Final Proportion")
 
save_plot(filename = paste0(fig_path, "effective.rnott.pdf"), final.prop, base_aspect_ratio = 1.5)


final.size.test$final.size[which(final.size.test$rnott == 1.7)] = rnott1.7
compare.prop <- ggplot(final.size.test, aes(prop.infected, fill = as.factor(pop.size)))+
  geom_density(alpha = .5)+facet_wrap(~rnott, scales = "free") + 
  geom_vline(aes(xintercept = final.size), size = 1.5, color = "purple")+
  scale_fill_discrete(name="Population\n Size") +#,
                      #breaks=c("beta.hat.100", "beta.hat.1000", "beta.hat.10000"),
                      #labels=c("100", "1000", "10000")) +
  labs(x = "Final Prop. Infected", y = "Density")
save_plot(paste0(fig_path, "compare.final.prop.pdf"), compare.prop, base_height = 8, base_aspect_ratio = 1.5)

head(final.size.test)
mle.beta <- ggplot(beta.hat.master, aes(x = as.factor(rnott), y = rnott.hat, 
                            fill = as.factor(beta.hat))) +
  geom_boxplot() + 
  scale_y_continuous(limits = c(0.9, 3.1), breaks = seq(0.9, 4, .2)) +
  scale_fill_discrete(name="Population Size",
                      breaks=c("beta.hat.100", "beta.hat.1000", "beta.hat.10000"),
                      labels=c("100", "1000", "10000")) +
  labs(x = "Simulation Beta/Gamma", y = "rnott.hat", title = "MLE of Beta")

save_plot(filename = paste0(fig_path, "mle.beta.pdf"), mle.beta, base_aspect_ratio = 1.5)
  