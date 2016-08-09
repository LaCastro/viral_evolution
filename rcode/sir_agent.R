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
rm(list=ls())
if(grepl('meyerslab', Sys.info()['login'])) setwd('~/Documents/projects/viral_evolution/rcode/')
if(grepl('laurencastro', Sys.info()['login'])) setwd('~/Documents/projects/viral_evolution/rcode/')

require("deSolve")
source("sir_agent_func.R")
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)

set.seed(578194)
nrealisations = 50

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
gamma = 1/3         # recovery period of influenza in days^{-1}
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

#cat("The final size of epidemic from the deterministic model is ",max(sirmodel$R/N),"\n")

##################################################################################
# now do several simulations using the agent based model, and overlay the
# results on those from the deterministic model
# myagent returns list of time, I, S, and final size 
##################################################################################
vfinal = numeric(0) # we will fill this vector with the epidemic final size estimates from the simulations
runs <- list()

for (iter in 1:nrealisations){
  myagent = SIR_agent(N = N,I_0 = I_0,S_0 = S_0,gamma = gamma,R0,tbeg,tend,delta_t)
  
  I_sim = myagent$I/N
  runs[[iter]] <- data.frame(trial = rep(iter, length(I_sim)), time = myagent$time, I = I_sim)
  #lines(myagent$time,myagent$I/N,lwd=2,col=(iter+1),lty=3)

  #cat(iter,nrealisations,"The final size of the epidemic from the agent based stochastic model is ",myagent$final_size,"\n")
  
  vfinal = append(vfinal,myagent$final_size) 
}
vfinal <- data.frame(vfinal)
runs.master.df <- do.call("rbind", runs) 

stochastic <- ggplot(runs.master.df, aes(x = time, y = I, color = trial,  group = trial)) +   geom_line() + guides(color = FALSE) 

combined <- ggplot() + 
  geom_line(data = runs.master.df, aes(x = time, y = I, color = trial, group = trial)) + 
  guides(color = FALSE) + 
  geom_line(data = sirmodel, aes(x = time, y = I_N), linetype = 5, color = "red", size = 1.5)  + 
  labs(title = paste("Influenza pandemic in population of", N, sep = " "),
       x = "Time, in days", y = "Fraction infected (prevalence)") 


## Need to figure out how to put legend right
#legend("topright",legend=c("Deterministic","Agent Based simulation"),col=c(1,2),lwd=3,lty=c(1,3),bty="n")

final.size <- ggplot(data = vfinal, aes(vfinal)) + geom_histogram(binwidth = .025) + 
  geom_vline(xintercept = max(sirmodel$R/N), size = 1.5, colour="red", linetype = "longdash") +
  labs(x = "Distribution of epidemic final size", main = "", y = "Count")

plot_grid(combined, final.size, ncol = 1)
