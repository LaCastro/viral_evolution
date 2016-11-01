##################################################################################
##################################################################################
# An R script to perform a stochastic epidemic simulation using an
# Agent Based Model and homogeneous mixing in the population.
#
# Author:  Sherry Towers
#          smtowers@asu.edu
# Created: Feb 12, 2016
#
# Copyright Sherry Towers, 2016
#
# This script is not guaranteed to be free of bugs and/or errors
#
# This script can be freely used and shared as long as the author and 
# copyright information in this header remain intact.
##################################################################################
set.seed(578194)

require("deSolve")
source("sir_agent_func.R")
nrealisations = 10

##################################################################################
##################################################################################
# this is a function which, given a value of S,I and R at time t
# calculates the time derivatives of S I and R
# vparameters contains the parameters of the model, like the
# recovery period, gamma, and the transmission rate, beta
# this function gets passed to the deSolve package
##################################################################################
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



##################################################################################
##################################################################################
# Set up initial conditions
##################################################################################
N = 10000    # population size
I_0 = 10     # number intially infected people in the population
S_0 = N-I_0
R_0 = 0      # assume no one has recovered at first

delta_t = 0.1       # nominal time step
tbeg  = 0           # begin day
tend  = 120         # end day
gamma = 1/3         # recovery period of influenza in days^{-1}
R0    = 2     # R0 of a hypothetical strain of pandemic influenza
beta = R0*gamma     # "reverse engineer" beta from R0 and gamma

##################################################################################
# first simulate the model with deterministic ODE's, so that we have something
# to compare our stochastic simulation to.
##################################################################################
vt = seq(tbeg,tend,delta_t)
vparameters = c(gamma=gamma,beta=beta)
inits = c(S=S_0,I=I_0,R=R_0)

sirmodel = as.data.frame(lsoda(inits, vt, SIRfunc, vparameters))

##################################################################################
# now plot the results of the deterministic model
##################################################################################
par(mfrow=c(2,1))  # divides the page into two plotting areas 
plot(sirmodel$time,sirmodel$I/N,ylim=c(0,1.5*max(sirmodel$I/N)),type="l",col=1,lwd=5,xlab="Time, in days",ylab="Fraction infected (prevalence)",main=paste("Influenza pandemic in population of ",N,sep=""))
cat("The final size of epidemic from the deterministic model is ",max(sirmodel$R/N),"\n")
deterministic = max(sirmodel$R/N)

##################################################################################
# now do several simulations using the agent based model, and overlay the
# results on those from the deterministic model
##################################################################################
time_seq <- seq(0.1, 1, .1)
deterministic = max(sirmodel$R/N)
nrealisations = 100
R0 = 1.2

mean.infect.moved <- adply(.data = time_seq, .margins = 1, .expand = FALSE, function(x) {
  delta_t = x
  vfinal = numeric(0)
  for (iter in 1:nrealisations) {
    myagent = SIR_agent(N,I_0,S_0,gamma,R0,tbeg,tend,delta_t)
    vfinal = append(vfinal,myagent$final_size) 
  }
  final.size = cbind(delta_t, as.data.frame(vfinal))
  return(final.size)
})

exponential.plot.moved.lowr0 <- ggplot(mean.infect.moved, aes(x = as.factor(delta_t), vfinal))+
  geom_jitter(alpha = .5) + 
  geom_violin(alpha = .5 )+
  geom_hline(yintercept = deterministic, color = "purple") +
  theme(strip.background = element_blank()) +
  labs(x = "Delta_t (Days)", y = "Final Epidemic Size", title = "Exponential r0 = .95 v2")

mean.infect %>% group_by(delta_t) %>%
  dplyr::summarize(average = mean(vfinal)) %>% 
  ggplot2::ggplot(aes(x = delta_t, y = average)) + geom_point() + 
      geom_hline(yintercept = deterministic, color = "purple") +
      labs(x = "delta t", y = "Average Final Epidemic Size") -> average.infect.plot



lines(sirmodel$time,sirmodel$I/N,lwd=5)
cat("The final size of epidemic from the deterministic model is ",max(sirmodel$R/N),"\n")
legend("topright",legend=c("Deterministic","Agent Based simulation"),col=c(1,2),lwd=3,lty=c(1,3),bty="n")

hist(vfinal,xlab="Distribution of epidemic final size",main="") # histogram the final size
lines(c(max(sirmodel$R/N),max(sirmodel$R/N)),c(-1000,1000),col=4,lwd=2,lty=2)
legend("topleft",legend=c("Agent Based simulation final size","Deterministic model final size"),col=c(1,4),lwd=2,lty=c(1,2),bty="n")



########### Exponentially distributed 
set.seed(578194)
SIRfunc_erlang=function(t, x, vparameters){
  k = length(x)-2
  
  S = x[1]  
  I = x[2:(length(x)-1)]  
  R = x[length(x)]  
  
  #cat(length(S),length(I),length(R),"\n")
  with(as.list(vparameters),{
    npop = S+sum(I)+R   
    dS = -beta*S*sum(I)/npop            
    dI = rep(0,length(I))
    dI[1] = +beta*S*sum(I)/npop - k*gamma*I[1]  
    if (k>1){
      dI[2:k] = +k*gamma*I[1:(k-1)] - k*gamma*I[2:k]
    }
    dR = +k*gamma*I[k]                 
    out = c(dS,dI,dR)
    list(out)
  })
}



##################################################################################
# Set up initial conditions
##################################################################################
set.seed(578194)
N = 10000    # population size
delta_t = 0.1       # nominal time step
tbeg  = 0           # begin day
tend  = 120         # end day
gamma = 1/3         # recovery period of influenza in days^{-1}
R0    = 1.2       # R0 of a hypothetical strain of pandemic influenza
beta = R0*gamma     # "reverse engineer" beta from R0 and gamma
k = 5

I_0 = rep(0,k)
I_0[1] = 10
S_0 = N-sum(I_0)
R_0 = 0      # assume no one has recovered at first

##################################################################################
# first simulate the model with deterministic ODE's, so that we have something
# to compare our stochastic simulation to.
##################################################################################
vt = seq(tbeg,tend,delta_t)
vparameters = c(gamma=gamma,beta=beta)
inits = c(S=S_0,I=I_0,R=R_0)

sirmodel = as.data.frame(lsoda(inits, vt, SIRfunc_erlang, vparameters))
sirmodel$I = rowSums(sirmodel[3:(3+k-1)])
deterministic =  max(sirmodel$R/N)
##################################################################################
# now plot the results of the deterministic model
##################################################################################
par(mfrow=c(2,1))  # divides the page into two plotting areas 
plot(sirmodel$time,sirmodel$I/N,ylim=c(0,1.5*max(sirmodel$I/N)),type="l",col=1,lwd=5,xlab="Time, in days",ylab="Fraction infected (prevalence)",main=paste("Influenza pandemic in population of ",N,sep=""))
cat("The final size of epidemic from the deterministic model is ", max(sirmodel$R/N),"\n")

##################################################################################
# now do several simulations using the agent based model, and overlay the
# results on those from the deterministic model
##################################################################################
nrealisations = 100

vfinal = numeric(0)
for (iter in 1:nrealisations) {
  myagent = SIR_agent_erlang(N,sum(I_0),S_0,gamma,k,R0,tbeg,tend,delta_t)
  vfinal = append(vfinal,myagent$final_size) 
}
vfinal.moved = vfinal

vfinal.original = numeric(0)
for (iter in 1:nrealisations) {
  myagent = SIR_agent_erlang(N,sum(I_0),S_0,gamma,k,R0,tbeg,tend,delta_t)
  vfinal.original= append(vfinal.original,myagent$final_size) 
}
 
mean.infect.erlang <- adply(.data = time_seq, .margins = 1, .expand = FALSE, function(x) {
  delta_t = x
  vfinal = numeric(0)
  for (iter in 1:nrealisations) {
      myagent = SIR_agent_erlang(N,sum(I_0),S_0,gamma,k,R0,tbeg,tend,delta_t)
      vfinal = append(vfinal,myagent$final_size) 
  }
  final.size = cbind(delta_t, as.data.frame(vfinal))
  return(final.size)
})

erlang.plot.r1.2 <- ggplot(mean.infect.erlang, aes(x = as.factor(delta_t), vfinal))+
  geom_jitter(alpha = .5) + 
  geom_violin(alpha = .5) +
  geom_hline(yintercept = deterministic, color = "purple") +
  theme(strip.background = element_blank()) +
  labs(x = "Delta_t (Days)", y = "Final Epidemic Size", title = "ErLang r0 = 1.2")


mean.infect.erlang %>% group_by(delta_t) %>%
  dplyr::summarize(average = mean(vfinal)) %>% 
  ggplot2::ggplot(aes(x = delta_t, y = average)) + geom_point() + 
  geom_hline(yintercept = deterministic, color = "purple") +
  labs(x = "delta t", y = "Average Final Epidemic Size") -> average.infect.plot

mean.infect.erlang %>% filter(delta_t == 0.1) %>%
  ggplot(aes(vfinal)) + geom_histogram(bins = 5)+theme(strip.background = element_blank())


meta.data %>% gather(method, final.size, -delta_t) %>%
  group_by(delta_t, method) %>%
  ggplot(aes(x = as.factor(delta_t), y = final.size, color = method))+
  geom_jitter(alpha = .3) 




lines(sirmodel$time,sirmodel$I/N,lwd=5)
cat("The final size of epidemic from the deterministic model is ",max(sirmodel$R/N),"\n")
legend("topright",legend=c("Deterministic","Agent Based simulation"),col=c(1,2),lwd=3,lty=c(1,3),bty="n")

hist(vfinal,xlab="Distribution of epidemic final size",main="") # histogram the final size
lines(c(max(sirmodel$R/N),max(sirmodel$R/N)),c(-1000,1000),col=4,lwd=2,lty=2)
legend("topleft",legend=c("Agent Based simulation final size","Deterministic model final size"),col=c(1,4),lwd=2,lty=c(1,2),bty="n")
