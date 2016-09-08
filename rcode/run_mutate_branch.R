##################################################################################
##################################################################################
# An R script to perform a stochastic epidemic simulation using an
# Agent Based Model and homogeneous mixing in the population, while mutation is occuring 
###################################################################################
rm(list=ls())
if(grepl('meyerslab', Sys.info()['login'])) setwd('~/Documents/projects/viral_evolution/viral_evolution_repo/rcode/')
if(grepl('laurencastro', Sys.info()['login'])) setwd('~/Documents/projects/viral_evolution/rcode/')

sapply(c('sir_agent_func.R','evo_functions.R', 'sir_mutation_func.R', 'analyze_saved_sims.R', 'plotting_functions.R'), source)

library(deSolve)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(reshape2)
library(plyr)

##################################################################################
# Set up initial functions to designate parameters and multiple runs 
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

run_mutate_branches_inc <- function(num_reps, ...) {
  rlply(.n = num_reps, .expr = sir_mutation_agent(...) ) 
}

##################################################################################
# Set up meta-analysis frames to keep track off 
##################################################################################
params <- epi_mut_params(N = 1000, I_0 = 10, delta_t = 1)
nrealisations = 100

trial <- run_mutate_branches_inc(num_reps = nrealisations, params)


#Getting Final Sizes (Proportions of each outbreak)
epi.size.all <- epi_size_all(trial = trial) 

## Getting Combined Time Records
time.records.all <- time_records_all(trial)


## Combining 
strain.records.all <- strain_freq_all(trial)


### Plotting 

time.max.diversity <- all_time_max_diversity(time.records.all)
time.max.infected <- all_time_max_infected(time.records.all)


debug(plot_max_times)

plot_max_times(time.records.all, type = "divergence")



## Histograms 
plot_final_sizes(epi.size.all)
plot_max_divergence(time.records.all)
plot_max_diversity(time.records.all)


### Line Plots 
combined.time.records <- combine_time_records(time.records.all)
stochastic.vI <- ggplot(combined.time.records, aes(x = vtime, y = vI, color = iter,  group = iter)) +  
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




### Left over from sampling 
sampling <- ggplot(runs.master.m, aes(x = time, y = value, color = as.factor(variable), linetype = as.factor(variable), fill = as.factor(variable))) +
  geom_line(size = 1, alpha = 0.2) +
  stat_summary(fun.y = "mean", color = "black", size = 1, geom = "line") + 
  guides(linetype = FALSE) +
  labs(x = "Time", y = "Number of Cases", color = "Type of Case")
sampling
