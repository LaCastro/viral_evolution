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
params <- epi_mut_params(N = 1000, I_0 = 10, delta_t = 1, R0 = 4)
nrealisations = 100

r0_seq <- seq(1, 4, 0.5)

## Figure out how to write this so you don't have to set it up 

trial_1 <- run_mutate_branches_inc(num_reps = nrealisations, params)
trail_1.5 <- run_mutate_branches_inc(num_reps = nrealisations, params)
trial_2 <- run_mutate_branches_inc(num_reps = nrealisations, params)
trial_2.5 <- run_mutate_branches_inc(num_reps = nrealisations, params)
trial_3 <- run_mutate_branches_inc(num_reps = nrealisations, params)
trial_3.5 <- run_mutate_branches_inc(num_reps = nrealisations, params)
trial_4 <- run_mutate_branches_inc(num_reps = nrealisations, params)


#Getting Final Sizes (Proportions of each outbreak)
epi.size.all <- epi_size_all(trial = trial_1) 


## Creating Max Metric 
master.data <- data.frame(matrix(ncol = 3))
colnames(master.data) <- c("rnott", "divergence", "diversity")
master.data <- master.data[-1,]

master.time.data <- data.frame(matrix(ncol = 3))
colnames(master.time.data) <- c("rnott", "divergence", "diversity")



time.records <- time_records_all(trial_4)
divergence.trial <- all_time_max_diverge(time.records)
diversity.trial <- all_time_max_diversity(time.records)
data <- cbind(divergence.trial, diversity.trial)
data <- data.frame(cbind(rep(4, length(diversity.trial)), data))
colnames(data) <- c("rnott", "divergence", "diversity")
master.time.data <- rbind(master.time.data, data)
master.time.data <- master.time.data[-1,]


master.data.threshold <- data.frame(matrix(ncol = 4))
colnames(master.data.threshold) <- c("rnott", "divergence", "diversity", "time")


time.records <- time_records_all(trial_4)
data <- threshold_metrics(threshold = .05, N, records.list = time.records)
data <- data.frame(cbind(rep(4, nrow(data)), data))
colnames(data) <- c("rnott", "divergence", "diversity", "time")
master.data.threshold <-  rbind(master.data.threshold , data)
master.data.threshold <- master.data.threshold[-1,]



### boxplots of Max diversity and divergence with R0s
master.data.threshold.m <- melt(data = master.data.threshold, id.vars = c("rnott"), measure.vars = c("divergence", "diversity", "time"))
colnames(master.data.threshold.m) <- c("rnott", "type", "value") 

tail(master.data.threshold.m)

threshold.plot <- ggplot(master.data.threshold.m, aes(factor(rnott), value)) + facet_grid(type~., scales = "free") +
  geom_boxplot()+
  labs(x = expression("R"[0]))

threshold.plot

save_plot(paste0(fig_path, "threshold.plot.pdf"), threshold.plot, base_height = 8, base_aspect_ratio = 1.2)



### Get time of max infected, will compare against time of max diversity 

master.max.time <- data.frame(matrix(ncol = 2))
colnames(master.max.time ) <- c("rnott", "time")


time.records <- time_records_all(trial_4)
data <- all_time_max_infected(x = time.records)
data <- data.frame(cbind(rep(4,length(data)), data))
colnames(data) <- c("rnott", "time")
master.max.time <-  rbind(master.max.time , data)

master.max.time <- master.max.time[-1,]
# need to get times of max diversity/divergence 
colnames(master.max.time) <- c("rnott", "infected_time")

head(master.max.time)
head(master.time.data)
merged.max.data <- merge(x = master.max.time, y = master.time.data, by = "rnott")
head(merged.max.data)


merged.data.m <- melt(data = merged.max.data, id.vars = c("rnott", "infected_time"), measure.vars = c("divergence", "diversity"))

head(merged.data.m)

max.data.plot <- ggplot(merged.data.m, aes(x = infected_time, y = value)) + geom_point() + facet_grid(variable ~ rnott)
                
max.data.plot


                            aes(color = factor(rnott))) + facet_wrap(rnott ~ rnott., nrow = 1)

which(is.na(merged.data.m))

factor(merged.data.m$rnott)


### boxplots of Max diversity and divergence with R0s
master.data.m <- melt(data = master.data, id.vars = c("rnott"), measure.vars = c("divergence", "diversity"))
colnames(master.data.m) <- c("rnott", "type", "value") 
  

max.metric.plot <- ggplot(master.data.m, aes(factor(rnott), value)) + facet_wrap(~type, nrow = 1) +
  geom_boxplot()+
  labs(x = expression("R"[0]), y = "Distance")

save_plot(paste0(fig_path, "max.metric.plot.pdf"), max.metric.plot, base_height = 4, base_aspect_ratio = 2.2)




## Combining 

strain.records.all <- strain_freq_all(trial)


### Plotting 

time.max.diversity <- all_time_max_diversity(time.records.all)
time.max.infected <- all_time_max_infected(time.records.all)


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
