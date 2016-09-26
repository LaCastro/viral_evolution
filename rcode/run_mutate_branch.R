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
library(dplyr)

fig_path = "~/Documents/projects/viral_evolution/viral_evolution_repo/figs/"

##################################################################################
# Set up initial functions to designate parameters and multiple runs 
##################################################################################
set.seed(578194)


epi_mut_params <- function(N = 1000,
                           I_0 = 10,
                           S_0 = N-I_0,
                           delta_t = 1, 
                           tbeg = 1, 
                           tend = 500,
                           gamma = 1/3,
                           R0 = 1.5, 
                           beta = R0*gamma, 
                           seq_len = 100,
                           alphabet = c(1, 2, 3, 4),
                           year_mut_rate = 1.8*10^-2, #1.8*10^-3
                           mut_rate = year_mut_rate / 365 / delta_t) #per site per day per delta t)
return(as.list(environment()))


##################################################################################
## Analysis 
##################################################################################
params <- epi_mut_params(N = 1000, I_0 = 1, delta_t = 1, R0 = 1.5)
nrealisations = 100

N = c(100, 1000, 10000)
r0_seq = seq(.95,1.2, 0.05)

if(grepl('meyerslab', Sys.info()['login'])) data_path <- "~/Documents/projects/viral_evolution/viral_evolution_repo/data/"
if(grepl('laurencastro', Sys.info()['login'])) data_path <- "~/Documents/projects/viral_evolution/data/"

## multiple parameter combinations, could feed this a_ply 


### need to do this but for the strain frequencies 
 
for (size in 1:length(N)) {
  for (r0 in 1:length(r0_seq)) {
    trial <- run_mutate_branches_inc(num_reps = nrealisations,
                                     params = epi_mut_params(N = N[size], R0 = r0_seq[r0], delta_t = 1))
    filename.time <- paste0(data_path, "/trial/trial_rnott", r0_seq[r0], "_N", N[size])
    time.records <- time_records_all(trial)
    save(time.records, file = paste0(filename.time, ".RData"))
    
    filename.strain <- paste0(data_path,"/strain/strain_rnott", r0_seq[r0], "_N", N[size])
    strain.records <- strain_freq_all(trial)
    save(strain.records, file = paste0(filename.strain, ".RData"))
  }
}

trial.1000 <- run_mutate_branches_inc(num_reps = nrealisations,
                                 params = epi_mut_params(N = 10000))
