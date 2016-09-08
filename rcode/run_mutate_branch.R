##################################################################################
##################################################################################
# An R script to perform a stochastic epidemic simulation using an
# Agent Based Model and homogeneous mixing in the population, while mutation is occuring 
###################################################################################
rm(list=ls())
if(grepl('meyerslab', Sys.info()['login'])) setwd('~/Documents/projects/viral_evolution/viral_evolution_repo/rcode/')
if(grepl('laurencastro', Sys.info()['login'])) setwd('~/Documents/projects/viral_evolution/rcode/')

sapply(c('sir_agent_func.R','evo_functions.R', 'sir_mutation_func.R'), source)

library(deSolve)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(reshape2)
library(plyr)

##################################################################################
# Set up initial conditions and parameters 
##################################################################################
set.seed(578194)


epi_mut_params <- function(N = 1000,
                           I_0 = 10,
                           S_0 = N-I_0,
                           delta_t = 1, 
                           tbeg = 1, 
                           tend = 121,
                           gamma = 1/3,
                           R0 = 1.5, 
                           beta = R0*gamma, 
                           seq_len = 100,
                           alphabet = c(1, 2, 3, 4),
                           year_mut_rate = 1.8*10^-2, #1.8*10^-3
                           mut_rate = year_mut_rate / 365 / delta_t) #per site per day per delta t)
return(as.list(environment()))



##################################################################################
# Set up meta-analysis frames to keep track off 
##################################################################################
params <- epi_mut_params(N = 1000, I_0 = 10, delta_t = 1, R0 = 2)

epi_runs <- list()
epi_size_final = numeric(0) # vector with the epidemic final size estimates from the simulations
nrealisations = 100



for (iter in 1:nrealisations) {
  myagent = sir_mutation_agent(params)
  epi_runs[[iter]] <- data.frame(trial = rep(iter, nrow(myagent$time_record)), myagent$time_record) 
  epi_size_final = append(epi_size_final, myagent$final_size) 
}

epi_size_final <- data.frame(epi_size_final)
runs.master.df <- do.call("rbind", epi_runs) 

stochastic.vI <- ggplot(runs.master.df, aes(x = vtime, y = vI, color = trial,  group = trial)) +  
  geom_line() + guides(color = FALSE) 
stochastic.vI

stochastic.diverge <- ggplot(runs.master.df, aes(x = vtime, y = diverge, color = trial,  group = trial)) +  
  geom_line() + guides(color = FALSE) 
stochastic.diverge

number.circulating.strains <- ggplot(runs.master.df, aes(x = vtime, y = cir.strains, color = trial, group = trial)) + 
  geom_line() + guides(color = FALSE)
number.circulating.strains


stochastic.diversity <- ggplot(runs.master.df, aes(x = vtime, y = diversity, color = trial, group = trial)) + 
  geom_line() + guides(color = FALSE)
stochastic.diversity


final.size <- ggplot(data = epi_size_final, aes(epi_size_final)) + geom_histogram(binwidth = .025) + 
  #geom_vline(xintercept = max(sirmodel$R/N), size = 1.5, colour="red", linetype = "longdash") +
  labs(x = "Distribution of epidemic final size", main = "", y = "Count")
final.size

# Detected vs Cumulative
runs.master.m  <- melt(data = runs.master.df, id.vars = c( "trial","time"), measure.vars = c("D", "C"))


sampling <- ggplot(runs.master.m, aes(x = time, y = value, color = as.factor(variable), linetype = as.factor(variable), fill = as.factor(variable))) +
  geom_line(size = 1, alpha = 0.2) +
  stat_summary(fun.y = "mean", color = "black", size = 1, geom = "line") + 
  guides(linetype = FALSE) +
  labs(x = "Time", y = "Number of Cases", color = "Type of Case")
sampling
