##################################################################################
##################################################################################
# An R script to perform a stochastic epidemic simulation using an
# Agent Based Model and homogeneous mixing in the population.
##################################################################################

rm(list=ls())
if(grepl('meyerslab', Sys.info()['login'])) setwd('~/Documents/projects/viral_evolution/viral_evolution_repo/rcode/')
if(grepl('laurencastro', Sys.info()['login'])) setwd('~/Documents/projects/viral_evolution/rcode/')

require("deSolve")
source("sir_agent_func.R")
source("samping_agent_func.R")
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(reshape2)

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

delta_t = .1      # nominal time step
tbeg  = 0           # begin day
tend  = 120         # end day
gamma = 1/5         # recovery period of influenza in days^{-1}
R0    = 1.5         # R0 of a hypothetical strain of pandemic influenza
beta = R0*gamma     # "reverse engineer" beta from R0 and gamma

##################################################################################
# first simulate the model with deterministic ODE's, so that we have something
# to compare our stochastic simulation to.
##################################################################################
vt = seq(tbeg,tend,delta_t)
vparameters = c(gamma=gamma,beta=beta)
inits = c(S=S_0,I=I_0,R=R_0)

sirmodel = as.data.frame(lsoda(y = inits, times = vt, func = SIRfunc, parms = vparameters))
sirmodel$I_N = sirmodel$I/N #Getting Percentage Infected (Prevalence)

##################################################################################
# now plot the results of the deterministic model
##################################################################################
#par(mfrow=c(2,1))  # divides the page into two plotting areas 

plot_deterministic <- ggplot(data = sirmodel, aes(x = time, y = I_N)) + 
  geom_line(size = 3) + ylim(0, 1.5*max(sirmodel$I_N)) +  
  labs(title = paste("Influenza pandemic in population of", N, sep = " "),
       x = "Time, in days", y = "Fraction infected (prevalence)") + 
  theme_cowplot() %+replace% theme(strip.background=element_blank(), strip.text.x = element_blank(),
 legend.title.align = 0.5,legend.position = c(0.2, 0.17)) 

plot_deterministic
#cat("The final size of epidemic from the deterministic model is ",max(sirmodel$R/N),"\n")

##################################################################################
# now do several simulations using the agent based model, and overlay the
# results on those from the deterministic model
# myagent returns list of time, I, S, and final size 
##################################################################################
set.seed(578194)
nrealisations = 25


vfinal = numeric(0) # we will fill this vector with the epidemic final size estimates from the simulations
runs <- list()



for (iter in 1:nrealisations) {
  #myagent = SIR_agent(N = N,I_0 = I_0,S_0 = S_0,gamma = gamma,R0,tbeg,tend,delta_t)
  myagent = SIR_agent_detection(N = N, I_0=I_0, S_0=S_0, gamma = gamma, R0, tbeg, tend, delta_t) 
  I_sim = myagent$I/N
  runs[[iter]] <- data.frame(trial = rep(iter, length(I_sim)), time = myagent$time, I = I_sim, D = myagent$D, C = myagent$C)
  #lines(myagent$time,myagent$I/N,lwd=2,col=(iter+1),lty=3)

  #cat(iter,nrealisations,"The final size of the epidemic from the agent based stochastic model is ",myagent$final_size,"\n")
  
  vfinal = append(vfinal,myagent$final_size) 
}

vfinal <- data.frame(vfinal)
runs.master.df <- do.call("rbind", runs) 

stochastic <- ggplot(runs.master.df, aes(x = time, y = I, color = trial,  group = trial)) +   geom_line() + guides(color = FALSE) 
stochastic


combined <- ggplot() + 
  geom_line(data = runs.master.df, aes(x = time, y = I, color = trial, group = trial)) + 
  guides(color = FALSE) + 
  geom_line(data = sirmodel, aes(x = time, y = I_N), linetype = 5, color = "red", size = 1.5)  + 
  labs(title = paste("Influenza pandemic in population of", N, sep = " "),
       x = "Time, in days", y = "Fraction infected (prevalence)") 

combined
## Need to figure out how to put legend right
#legend("topright",legend=c("Deterministic","Agent Based simulation"),col=c(1,2),lwd=3,lty=c(1,3),bty="n")

final.size <- ggplot(data = vfinal, aes(vfinal)) + geom_histogram(binwidth = .025) + 
  geom_vline(xintercept = max(sirmodel$R/N), size = 1.5, colour="red", linetype = "longdash") +
  labs(x = "Distribution of epidemic final size", main = "", y = "Count")

final.size
plot_grid(combined, final.size, ncol = 1)
vfinal


# Detected vs Cumulative

plot_prevalences <- function(df){
  
  df$disc_prob <- paste0(calculate.discover(df$disc_prob), "%")
  
  ggplot(df, aes(detected, median, color=as.factor(r_not), linetype=as.factor(disc_prob), fill=as.factor(r_not))) + 
    geom_line(size=1)+
    geom_ribbon(aes(ymax=max, ymin=min), alpha=0.1, color=NA)+
    scale_y_log10(expand=c(0.01,0.01), limits=c(1,100), breaks = c(5,10,25,50,100))+
    scale_x_continuous(expand=c(0.01,0.01), limits=c(0,30))+
    scale_color_brewer(palette="Set1", direction = -1)+
    scale_fill_brewer(palette="Set1", direction=-1)+
    guides(linetype=FALSE)+
    labs(x = "Cumulative Reported Cases", 
         y = "Prevalence (log scale)", 
         color = expression("R"[0]), 
         fill = expression("R"[0]))
}

runs.master.m  <- melt(data = runs.master.df, id.vars = c( "trial","time"), measure.vars = c("D", "C"))


sampling <- ggplot(runs.master.m, aes(x = time, y = value, color = as.factor(variable), linetype = as.factor(variable), fill = as.factor(variable))) +
  geom_line(size = 1, alpha = 0.2) +
  stat_summary(fun.y = "mean", color = "black", size = 1, geom = "line") + 
  guides(linetype = FALSE) +
  labs(x = "Time", y = "Number of Cases", color = "Type of Case")
sampling


############## Erlang Distribution
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








N = 10000    # population size

delta_t = 0.1       # nominal time step
tbeg  = 0           # begin day
tend  = 120         # end day
gamma = 1/3         # recovery period of influenza in days^{-1}
R0    = 1.5         # R0 of a hypothetical strain of pandemic influenza
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

##################################################################################
# now plot the results of the deterministic model
##################################################################################
par(mfrow=c(2,1))  # divides the page into two plotting areas 
plot(sirmodel$time,sirmodel$I/N,ylim=c(0,1.5*max(sirmodel$I/N)),type="l",col=1,lwd=5,xlab="Time, in days",ylab="Fraction infected (prevalence)",main=paste("Influenza pandemic in population of ",N,sep=""))
cat("The final size of epidemic from the deterministic model is ",max(sirmodel$R/N),"\n")

##################################################################################
# now do several simulations using the agent based model, and overlay the
# results on those from the deterministic model
##################################################################################
vfinal = numeric(0) # we will fill this vector with the epidemic final size estimates from the simulations
for (iter in 1:nrealisations){
  myagent = SIR_agent_erlang(N,sum(I_0),S_0,gamma,k,R0,tbeg,tend,delta_t)
  lines(myagent$time,myagent$I/N,lwd=2,col=(iter+1),lty=3)
  cat(iter,nrealisations,"The final size of the epidemic from the agent based stochastic model is ",myagent$final_size,"\n")
  vfinal = append(vfinal,myagent$final_size) 
}
lines(sirmodel$time,sirmodel$I/N,lwd=5)
cat("The final size of epidemic from the deterministic model is ",max(sirmodel$R/N),"\n")
legend("topright",legend=c("Deterministic","Agent Based simulation"),col=c(1,2),lwd=3,lty=c(1,3),bty="n")

hist(vfinal,xlab="Distribution of epidemic final size",main="") # histogram the final size
lines(c(max(sirmodel$R/N),max(sirmodel$R/N)),c(-1000,1000),col=4,lwd=2,lty=2)
legend("topleft",legend=c("Agent Based simulation final size","Deterministic model final size"),col=c(1,4),lwd=2,lty=c(1,2),bty="n")

