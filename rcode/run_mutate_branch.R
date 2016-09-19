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
## Analysis 
##################################################################################
params <- epi_mut_params(N = 100, I_0 = 1, delta_t = 1, R0 = 1.5)
nrealisations = 100

N = c(100, 1000, 10000)
r0_seq = seq(1,4, 0.5)

data_path <- "~/Documents/projects/viral_evolution/viral_evolution_repo/data/"

## multiple parameter combinations, could feed this a_ply 


### need to do this but for the strain frequencies 
 
for (size in 1:length(N)) {
  for (r0 in 1:length(r0_seq)) {
    trial <- run_mutate_branches_inc(num_reps = nrealisations,
                                     params = epi_mut_params(N = N[size], R0 = r0_seq[r0], delta_t = 1))
     filename <- paste0(data_path,"trial_rnott", r0_seq[r0], "_N", N[size])
     
     time.records <- time_records_all(trial)
     
     save(time.records, file = paste0(filename, ".RData"))
  }
}


#Getting Final Sizes (Proportions of each outbreak)
epi.size.all <- epi_size_all(trial = trail_1.5_10e4) 


### Reading all desired data and combining for analysis 
fig_path <- "~/Documents/projects/viral_evolution/figs/"
data_path <- "~/Documents/projects/viral_evolution/data/"

file.list <- get_vec_of_files(dir_path = data_path, r0_seq, N)
combos <- sort(apply(expand.grid(r0_seq, N), 1, paste, collapse = "_", sep = "")) 
data.files <- data.frame(cbind(combos, file.list))

###### Code for producing boxplots of max genetic variation metrics 
## Want to load the data and calculate the max metrics for all 
max.metrics.master <-  adply(.data = data.files, .margins = 1, function(x) {
  load(as.character(x$file.list))
  rnott = as.numeric(sub("_.*", "", x$combos))
  pop.size = as.numeric(sub(".*_", "", x$combos))
  genetic.metrics <- get_max_genetic_metrics(time.records = time.records, rnott = rnott)
  genetic.metrics <-  cbind(rep(pop.size, nrow(genetic.metrics)), genetic.metrics)
  colnames(genetic.metrics)[1] <- "pop.size"
  return(genetic.metrics)
})

# plotting code 
max.metrics.m <- melt(data = max.metrics.master, id.vars = c("rnott", "pop.size"),
                      measure.vars = c("max.divergence", "max.diversity"))
colnames(max.metrics.m) <- c("rnott", "pop.size", "type", "value") 
max.metric.plot <- ggplot(max.metrics.m, aes(factor(rnott), value)) + 
  facet_grid(pop.size~type) +
  geom_boxplot()+
  labs(x = expression("R"[0]), y = "Distance")
save_plot(paste0(fig_path, "max.metric.plot.pdf"), max.metric.plot, base_height = 8, base_aspect_ratio = 1.2)




######## Code for Combing threshold metrics and plotting 
threshold.metrics.master <- adply(.data = data.files, .margins = 1,.id=NULL, .expand = F, function(x) {
  load(as.character(x$file.list))
  rnott = as.numeric(sub("_.*", "", x$combos))
  pop.size = as.numeric(sub(".*_", "", x$combos))
  threshold.metrics = threshold_metrics(threshold = .1, pop.size, records.list = time.records, rnott = rnott)
  threshold.metrics = cbind(rep(pop.size, nrow(threshold.metrics)), threshold.metrics)
  colnames(threshold.metrics)[1] <- "pop.size"
  return(threshold.metrics)
})

threshold.metrics.m <- melt(data = threshold.metrics.master, id.vars = c("rnott", "pop.size"), 
                            measure.vars = c("diverge", "diversity"))
colnames(threshold.metrics.m) <- c("rnott", "pop.size", "type", "value") 
threshold.plot <-ggplot(threshold.metrics.m, aes(factor(rnott), value)) + 
  facet_grid(pop.size~type) +
  geom_boxplot()+
  labs(x = expression("R"[0]), y = "Distance")
save_plot(paste0(fig_path, "threshold.metric.plot.1v2.pdf"), threshold.plot, base_height = 8, base_aspect_ratio = 1.2)



######## Code for combining beginning threshold metrics and plotting 
beg.threshold.metrics.master <- adply(.data = data.files, .margins = 1, .id=NULL, .expand = F, function(x) {
  load(as.character(x$file.list))
  rnott = as.numeric(sub("_.*", "", x$combos))
  pop.size = as.numeric(sub(".*_", "", x$combos))
  threshold.metrics = beginning_threshold_metrics(threshold = .1, pop.size, records.list = time.records, rnott = rnott)
  threshold.metrics = cbind(rep(pop.size, nrow(threshold.metrics)), threshold.metrics)
  colnames(threshold.metrics)[1] <- "pop.size"
  return(threshold.metrics)
})


beg.threshold.metrics.m <- melt(data = beg.threshold.metrics.master, id.vars = c("rnott", "pop.size"), 
                            measure.vars = c("diverge", "diversity"))
colnames(beg.threshold.metrics.m) <- c("rnott", "pop.size", "type", "value") 
threshold.plot <-ggplot(beg.threshold.metrics.m, aes(factor(rnott), value)) + 
  facet_grid(pop.size~type) +
  geom_boxplot()+
  labs(x = expression("R"[0]), y = "Distance")
save_plot(paste0(fig_path, "beg.threshold.metric.plot.pdf"), threshold.plot, base_height = 8, base_aspect_ratio = 1.2)



###### Code to compare the average rate at which mutations accumulate depending on R0 

cum.strains.master <- adply(.data = data.files, .margins = 1, .id=NULL, .expand = F, function(x) {
  load(as.character(x$file.list))
  rnott = as.numeric(sub("_.*", "", x$combos))
  pop.size = as.numeric(sub(".*_", "", x$combos))
  meta.data = data.frame(cbind(rnott, pop.size))
  cum.strains.trial = combine_mutations(time.records = time.records)
  cum.strains.trial = cbind(meta.data, cum.strains.trial)
  return(cum.strains.trial)
})

## This works! 
cum.strains.groups <- group_by(cum.strains.master, rnott, pop.size, vtime)
cum.strains.avg <- summarise(cum.strains.groups, avg.strains = mean(cum.strains), sd.strains = sd(cum.strains))


# plot the average increase for each rnott, pop.size 
increase.mutations <- ggplot(cum.strains.avg, aes(x = vtime, y = avg.strains)) +
  geom_ribbon(aes(ymin = (avg.strains -sd.strains), ymax = (avg.strains + sd.strains)),
                                  alpha = .5, fill="steelblue2", color="steelblue2") +
  geom_line(color = "black", lwd = 1) +
  facet_grid(pop.size~rnott) +
  theme(panel.margin = unit(1.5, "lines")) + 
  labs(x = "Time", y = "Number of Cumulative Mutations")
  
save_plot(paste0(fig_path, "increase.mutations.plot.pdf"), increase.mutations, base_height = 8, base_aspect_ratio = 1.75)







######### Code to  Get time of max infected and compare against time  max diversity 
### Example will be a low R0  (5) and a high R0 (20) for N = 1000
## Want a two row figure for max diversity and max divergence 

low.r0.path <- file.list[5]
high.r0.path <- file.list[20]

load(low.r0.path)
times.low.r0 <- all_time_max_infected(time.records)
diversity.low.r0 <- all_time_max_diversity(time.records)
diverge.low.r0 <- all_time_max_diverge(time.records)

low.r0.data <- data.frame(cbind(times.low.r0, diversity.low.r0, diverge.low.r0))
low.r0.data <- cbind(rep("r0 = 1.5", nrow(low.r0.data)), low.r0.data)
colnames(low.r0.data) <- c("type", "infected", "diversity", "divergence")

load(high.r0.path)
times.high.r0 <- all_time_max_infected(time.records)
diversity.high.r0 <- all_time_max_diversity(time.records)
diverge.high.r0 <- all_time_max_diverge(time.records)
high.r0.data <- data.frame(cbind(times.high.r0, diversity.high.r0, diverge.high.r0))
high.r0.data <- cbind(rep("r0 = 4", nrow(high.r0.data)), high.r0.data)
colnames(high.r0.data) <- c("type", "infected", "diversity", "divergence")

time.data.master <- rbind(low.r0.data, high.r0.data)
head(time.data.master)
time.data.m <- melt(data = time.data.master, id.vars = c("type", "infected"),
                    measure.vars = c("diversity", "divergence"))


combine.low <- combine_time_records(time.records)
head(combine.low)

stochastic <- ggplot(combined.time.records.m, aes(x = vtime, y = value, color = iter,  group = iter)) +  
  geom_line() + guides(color = FALSE) +
  facet_grid(type~variable) + 
  theme_cowplot() %+replace% theme(strip.background=element_rect(fill=NULL,color="black", size=1, linetype=1)) +
  labs(x = "Time", y = "Distance")

stochastic <- ggplot(combined.time.records.m, aes(x = vtime, y = value, color = iter,  group = iter)) +  
  geom_line() + guides(color = FALSE) +
  facet_grid(type~variable) + 
  theme_cowplot() %+replace% theme(strip.background=element_rect(fill=NULL,color="black", size=1, linetype=1)) +
  labs(x = "Time", y = "Distance")

head(combine.low)

cir.strains.plot <- ggplot(combine.low, aes(x = vtime, y = cir.strains, color = iter,  group = iter)) +  
  geom_point() + guides(color = FALSE) +
  #facet_grid(type~variable) + 
  theme_cowplot() %+replace% theme(strip.background=element_rect(fill=NULL,color="black", size=1, linetype=1)) +
  labs(x = "Time", y = "Number of STrains")


time.max.low <- all_time_max_infected(time.records)
time.max.strains.low <- all_time_max_numstrains(time.records)
head(strains.time)

strains.time <- data.frame(cbind(time.max.low, time.max.strains.low))
strains.time  <- cbind(rep("r0 = 1.5", nrow(strains.time)), strains.time)
colnames(strains.time) <- c("type", "infected", "num.strains")



load(high.r0.path)
times.high.r0 <- all_time_max_infected(time.records)
times.high.strains <- all_time_max_numstrains(time.records)

strains.high <- data.frame(cbind(times.high.r0, times.high.strains))
strains.high <- cbind(rep("r0 = 4", nrow(strains.high)), strains.high)
colnames(strains.high) <- c("type", "infected", "num.strains")

strains.data.master <- rbind(strains.time, strains.high)
head(strains.data.master)
#time.data.m <- melt(data = time.data.master, id.vars = c("type", "infected"),
#                    measure.vars = c("diversity", "divergence"))



time.scatter.strains <- ggplot(data = strains.data.master, aes(x = infected, y = num.strains)) + 
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(facets = "type", nrow = 1) +
  theme_cowplot() %+replace% theme(strip.background=element_rect(fill=NULL,color="black", size=1, linetype=1)) +
  labs(x = "Time of Max Infected", y = "Time of Max Circulating Strains")


save_plot(paste0(fig_path, "strains.time.scatter.pdf"), time.scatter.strains, base_height = 4, base_aspect_ratio = 1.8)



### Scatterplot 
time.scatter <- ggplot(data = time.data.m, aes(x = infected, y = value)) + 
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  facet_grid(type~variable) +
  theme_cowplot() %+replace% theme(strip.background=element_rect(fill=NULL,color="black", size=1, linetype=1)) +
  labs(x = "Time of Max Infected", y = "Time of Max Distance")

save_plot(paste0(fig_path, "time.scatter.pdf"), time.scatter, base_height = 8, base_aspect_ratio = 1.2)


## Combining 
strain.records.all <- strain_freq_all(trial)


### Line Plots 
combined.high <- combine_time_records(time.records) 
combined.high <- cbind(rep("r0 = 4", nrow(combined.high)), combined.high)
colnames(combined.high)[1] <- "type"

combined.low <- combine_time_records(time.records)
combined.low <- cbind(rep("r0 = 1.5", nrow(combined.low)), combined.low)
colnames(combined.low)[1] <- "type"

combined.time.records <- rbind(combined.low, combined.high)
combined.time.records.m <- melt(data = combined.time.records, id.vars = c("type", "iter", "vtime"),
                                measure.vars = c("diverge", "diversity"))


stochastic <- ggplot(combined.time.records.m, aes(x = vtime, y = value, color = iter,  group = iter)) +  
  geom_line() + guides(color = FALSE) +
  facet_grid(type~variable) + 
  theme_cowplot() %+replace% theme(strip.background=element_rect(fill=NULL,color="black", size=1, linetype=1)) +
  labs(x = "Time", y = "Distance")

save_plot(filename = paste0(fig_path, "stochastic.line.pdf"), plot = stochastic, base_height = 4, base_width = 7)


stochastic.diverge <- ggplot(runs.master.df, aes(x = vtime, y = diverge, color = trial,  group = trial)) +  
  geom_line() + guides(color = FALSE) 
stochastic.diverge

number.circulating.strains <- ggplot(runs.master.df, aes(x = vtime, y = cir.strains, color = trial, group = trial)) + 
  geom_line() + guides(color = FALSE)
number.circulating.strains

stochastic.diversity <- ggplot(runs.master.df, aes(x = vtime, y = diversity, color = trial, group = trial)) + 
  geom_line() + guides(color = FALSE)
stochastic.diversity




## Histograms 
plot_final_sizes(epi.size.all)
plot_max_divergence(time.records)
plot_max_diversity(time.records)



### Left over from sampling 
sampling <- ggplot(runs.master.m, aes(x = time, y = value, color = as.factor(variable), linetype = as.factor(variable), fill = as.factor(variable))) +
  geom_line(size = 1, alpha = 0.2) +
  stat_summary(fun.y = "mean", color = "black", size = 1, geom = "line") + 
  guides(linetype = FALSE) +
  labs(x = "Time", y = "Number of Cases", color = "Type of Case")
sampling
