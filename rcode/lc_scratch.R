## Scratch Script for  working out functions
## Exploratory Analysis Figures 

rm(list=ls())
if(grepl('meyerslab', Sys.info()['login'])) setwd('~/Documents/projects/viral_evolution/viral_evolution_repo/rcode/')
if(grepl('laurencastro', Sys.info()['login'])) setwd('~/Documents/projects/viral_evolution/rcode/')

sapply(c('analyze_saved_sims.R', 'run_mutate_branch.R'), source)
library(deSolve)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(data.table)

## Setting Up Figure and Path directions
# Figure Path
if(grepl('meyerslab', Sys.info()['login'])) fig_path = "~/Documents/projects/viral_evolution/viral_evolution_repo/figs/"
if(grepl('laurencastro', Sys.info()['login'])) fig_path <- "~/Documents/projects/viral_evolution/figs/"

# Data Path 
if(grepl('meyerslab', Sys.info()['login'])) data_path <- "~/Documents/projects/viral_evolution/viral_evolution_repo/data/"
if(grepl('laurencastro', Sys.info()['login'])) data_path <- "~/Documents/projects/viral_evolution/data/"


# Reading in files from saved data runs for chosen R0s and popsizes
N = c(100, 1000, 10000)
r0_seq = c(seq(0.9, 2, .2), seq(2,5, 0.5))

file.list <- get_vec_of_files(dir_path = data_path, type = "trial", r0_seq, N)
combos <- sort(apply(expand.grid(r0_seq, N), 1, paste, collapse = "_", sep = "")) 
data.files <- data.frame(cbind(combos, file.list))

# Deciding which trials are epidemics based on a threshold 
# for prevalance and total proportion infected
epidemic.trials <- set_epidemic_criteria(time.records, threshold.prev = .025, threshold.prop = .25)
data.files$epidemic.trials <- epidemic.trials


#########################################################
# TIME SERIES ANALYSIS
#########################################################

# 1. Plotting infected time-series/genetic metrics 
trajectories <- adply(.data = data.files, .margins = 1,.id=NULL, .expand = F, function(x) {
  load(as.character(x$file.list))
  trial.params <- get_params(x)
  epidemic.index <- x$epidemic.trials[[1]]$epidemic.trials
  time.records.a <- align_time_series_all(time.records = time.records[epidemic.index])
  vI.trajectory <- combine_time_records(time.records.a)
  vI.trajectory <- cbind(trial.params, vI.trajectory)
  return(vI.trajectory)
})
colnames(trajectories) <- c("rnott", "pop.size", "iter", "vtime", "vI")

# plot a sample 
#random.iter <- sample(x = seq(1:1000), size = 50)
#trajectories %>% filter(iter %in% random.iter) -> trajectories.sample
#ggplot(data = trajectories.sample, aes(x = vtime, y = vI, group = iter, color = iter))+
#  geom_smooth()+guides(color = FALSE) +
#  facet_grid(rnott~pop.size, scales = "free")

average.vi <- group_by(trajectory.master, rnott, pop.size, shifted.time) %>%
  filter(rnott %in% desired.rnotts) %>%
  filter(iter %in% random.iter)  %>%
  summarise(infect = mean(vI),
            infect.sd = sd(vI)) 

stochastic.entropy <- ggplot(trajectory.long, aes(x=vtime, y = entropy)) +
  geom_line(size = 1.25) + 
  geom_ribbon(aes(ymin = (entropy -2*entropy.sd), ymax = (entropy + 2*entropy.sd)),
              alpha = .5, fill="indianred1", color="indianred1") +
  facet_grid(pop.size~rnott) + 
  labs(x = "Time", y = "Entropy") +
  theme(axis.text.x  = element_text(angle=45, vjust=0.5)) +
  theme(strip.background = element_blank()) + 
  theme(panel.margin = unit(1.5, "lines")) 


average.vi.plot <- ggplot(average.vi, aes(x = vtime, y = infect)) + geom_line(color = "red", linetype  = 1, size = 1.5)+
  facet_grid(pop.size~rnott, scales = "free_y") + 
  geom_ribbon(aes(ymin = (infect -2*infect.sd), ymax = (infect + 2*infect.sd)),
              alpha = .5, fill="indianred1", color="indianred1") +
  theme(strip.background = element_blank()) +
  theme(axis.text.x  = element_text(angle=45, vjust=0.5)) +
  labs(x = "Time", y = "Infected Individuals")  + 
  theme(panel.margin = unit(1.5, "lines")) 


# 2. Align Time Series and Summarise 
time.records.a <- align_time_series_all(time.records)
aligned.time.series <- combine_time_records(time.records.a)
aligned.time.series %>% group_by(shifted.time) %>%
  summarise(mean.diversity = mean(diversity),
            sd.diversity = sd(diversity),
            mean.entropy = mean(entropy),
            sd.entropy = sd(entropy),
            mean.cir.strains = mean(cir.strains),
            sd.cir.strains = sd(cir.strains),
            mean.cum.strains = mean(cum.strains),
            sd.cum.strains = sd(cum.strains)) -> aligned.time.series.average

# Plot by time 
align.ts.long <- gather(aligned.time.series, metric, value, -iter, -vtime, -shifted.time) %>%
  group_by(vtime) 

ggplot(align.ts.long, aes(x = shifted.time, y = value, group = iter, color = iter))+guides(color = FALSE)+
  facet_wrap(~metric, scales = "free") + geom_line()

#3. Compare the average rate at which mutations accumulate depending on R0 
cum.strains.master <- adply(.data = data.files, .margins = 1, .id=NULL, .expand = F, function(x) {
  load(as.character(x$file.list))
  trial.params <- get_params(x)
  epidemic.index <- x$epidemic.trials[[1]]$epidemic.trials
  if (epidemic.index[1] > 0) {
    time.records.a <- align_time_series_all(time.records[epidemic.index])
    cum.strains.trial = combine_mutations(time.records = time.records.a)
    cum.strains.trial = cbind(trial.params, cum.strains.trial)
    return(cum.strains.trial)
  }
})

cum.strains.master %>% group_by(rnott, pop.size, shifted.time) %>%
  summarise(avg.strains = mean(cum.strains), 
            sd.strains = sd(cum.strains))  %>%
  filter(rnott %in% desired.rnotts) -> cum.strains.avg

# 4. Trajectories Plots of Entropy - should I align these? 








##########################################
# Metrics Analysis 
##########################################

# 1. Get the length of time each simulation takes
final.time.master <-  adply(.data = data.files, .margins = 1,.id=NULL, .expand = F, function(x) {
  load(as.character(x$file.list))
  trial.params <- get_params(x)
  epidemic.index <- x$epidemic.trials[[1]]$epidemic.trials
  final.time.trial <- data.frame(all_epi_time(time.records[epidemic.index])); colnames(final.time.trial) = "days"
  final.time.trial <-  cbind(trial.params, final.time.trial)
  return(final.time.trial)
})

# Scatter Plot 
epi.time.plot <- ggplot(final.time.master, aes(factor(rnott), days)) + 
  facet_wrap(facets = ~pop.size, nrow = 1) +
  geom_jitter(alpha = .5) +
  labs(x = expression("R"[0]), y = "Time (days)")

# Histogram 
epi.time.plot.hist <- ggplot(final.time.master, aes(days)) + 
  facet_grid(rnott~pop.size) +
  geom_histogram() +
  labs(x = "Time (days)", y = "Count")
save_plot(paste0(fig_path, "epi.time.r0around1.pdf"), epi.time.plot, base_height = 8, base_aspect_ratio = 1.7)


# 2. Look at the proportion of the population that was infected 
final.infected.master <-  adply(.data = data.files, .margins = 1,.id=NULL, .expand = F, function(x) {
  load(as.character(x$file.list))
  trial.params <- get_params(x)
  epidemic.index <- x$epidemic.trials[[1]]$epidemic.trials
  final.infected.trial <- data.frame(all_total_infected(time.records = time.records[epidemic.trials], trial.N = trial.params$pop.size))
  final.infected.trial <- cbind(trial.params, final.infected.trial)
  return(final.infected.trial)
})
colnames(final.infected.master)[3] <- "prop.infected"

ggplot(final.infected.master, aes(factor(rnott), prop.infected)) + geom_jitter(alpha = .5) + 
  geom_violin(alpha = .5) + facet_wrap(~pop.size, nrow =1 ) + geom_hline(yintercept = .13, color = "purple") +
  geom_hline(yintercept = .25, color = "orange")

# 3. Maximum number of people infected at a single time
maximum.prevalance <- adply(.data = data.files, .margins = 1, .id= NULL, .expand = F, function(x) {
  load(as.character(x$file.list))
  trial.params <- get_params(x)
  max.prev <- data.frame(all_max_infected(x = time.records))
  histogram.data <- cbind(trial.params,max.prev)
  return(histogram.data)
})
colnames(maximum.prevalance)[3] <- "max.prev"

ggplot(data = maximum.prevalance, aes(max.prev))+geom_histogram(bins = 20) + 
  facet_grid(rnott~pop.size, scales = "free_x") 

# 3. Comparing the relationship between proportion and maximum prevalance-used to help determine threshold 
meta.data.master <- cbind(final.infected.master, final.time.master[,3]); colnames(meta.data.master)[4] <- "days"
colnames(meta.data.master) <- c("rnott", "pop.size", "prop.infected", "days", "max.prev", "ratio")

# Trying Different Relationships 
ggplot(data = meta.data.master, aes(x = ratio, y = prop.infected)) + geom_jitter(alpha = .5)+facet_grid(rnott~pop.size, scales = "free_x")
ggplot(data = meta.data.prop25, aes(x = days, y = prop.infected)) + geom_jitter(alpha = .5)+facet_grid(rnott~pop.size, scales = "free_x")
ggplot(meta.data.two.critera, aes(x = max.prev, y = prop.infected)) + geom_jitter(alpha = .5) + facet_grid(rnott~pop.size, scales = "free_x")

# 4. Calculate max metrics for all-at moximum of the metric itself 
max.metrics.master <-  adply(.data = data.files, .margins = 1, .id=NULL, .expand = F,function(x) {
  load(as.character(x$file.list))
  trial.params <- get_params(x)
  epidemic.index <- x$epidemic.trials[[1]]$epidemic.trials
  genetic.metrics <- get_max_genetic_metrics(time.records = time.records[epidemic.index], rnott = trial.params$rnott)
  genetic.metrics <-  cbind(rep(trial.params$pop.size, nrow(genetic.metrics)), genetic.metrics)
  colnames(genetic.metrics)[1] <- "pop.size"
  return(genetic.metrics)
})
max.metrics.m <- melt(data = max.metrics.master, id.vars = c("rnott", "pop.size"),
                      measure.vars = c("max.entropy", "max.diversity"))
colnames(max.metrics.m) <- c("rnott", "pop.size", "type", "value") 

max.metric.plot <- ggplot(max.metrics.m, aes(factor(rnott), value)) + 
  facet_grid(type~pop.size, scales = "free_y") +
  geom_boxplot()+
  labs(x = expression("R"[0]), y = "Distance")
#save_plot(paste0(fig_path, "max.metric.entropy.diversity.epidemics.pdf"),max.metric.plot, base_height = 8, base_aspect_ratio = 1.2)


# 5. Code for getting genetic metrics at time of max infection 
max.infect.metrics.master <-  adply(.data = data.files, .margins = 1, .id=NULL, .expand = F,function(x) {
  load(as.character(x$file.list))
  trial.params <- get_params(x)
  epidemic.index <- x$epidemic.trials[[1]]$epidemic.trials
  time.records.a <- align_time_series_all(time.records[epidemic.index])
  metrics <- metrics_at_max_infect(time.records = time.records.a)
  metrics <-  cbind(rep(trial.params, nrow(metrics)), metrics)
  return(metrics)
})

plot.metric.at.maxinfect(max.infect.metrics.master, desired.rnotts = r0_seq, desired.metric = "diversity")


# 6. Code for combining beginning threshold metrics and plotting 
beg.threshold.metrics.master <- adply(.data = data.files, .margins = 1, .id=NULL, .expand = F, function(x) {
  load(as.character(x$file.list))
  trial.params <- get_params(x)
  epidemic.index <- x$epidemic.trials[[1]]$epidemic.trials
  threshold.metrics = beginning_threshold_metrics(threshold = .1,trial.params$pop.size, 
                                                  records.list = time.records[epidemic.index], rnott = trial.params$rnott)
  threshold.metrics = cbind(rep(trial.params$pop.size, nrow(threshold.metrics)), threshold.metrics)
  colnames(threshold.metrics)[1] <- "pop.size"
  return(threshold.metrics)
})
beg.threshold.metrics.m <- melt(data = beg.threshold.metrics.master, id.vars = c("rnott", "pop.size"), 
                                measure.vars = c("entropy", "diversity"))
colnames(beg.threshold.metrics.m) <- c("rnott", "pop.size", "type", "value") 

threshold.plot <-ggplot(beg.threshold.metrics.m, aes(factor(rnott), value)) + 
  facet_grid(type~pop.size, scales = "free_y") +
  geom_boxplot()+
  labs(x = expression("R"[0]), y = "Distance")
#save_plot(paste0(fig_path, "beg.threshold.entropy.short.pdf"), threshold.plot, base_height = 8, base_aspect_ratio = 1.2)

# Combine the time points of max genetic metrics and infections 
time.max.master <- adply(.data = data.files, .margins = 1, .id=NULL, .expand = F, function(x) {
  load(as.character(x$file.list))
  trial.params <- get_params(x)
  epidemic.index <- x$epidemic.trials[[1]]$epidemic.trials
  max.infect <- all_time_max_infected(time.records[epidemic.index])
  max.diversity <- all_time_max_diversity(time.records[epidemic.index])
  max.entropy <- all_time_max_entropy(time.records[epidemic.index])
  max.times.trial <- cbind(trial.params,max.infect, max.diversity, max.entropy)
  return(max.times.trial)
})

# filter based on desired rnott 
time.max.groups <- group_by(time.max.master, rnott, pop.size) %>%
  summarise(infect = mean(max.infect),
            diversity = mean(max.diversity),
            entropy = mean(max.entropy)) %>%
  #gather(type, day, -rnott, -pop.size) %>%
  filter(rnott <= 1 | rnott >= 2)


# plot the average increase for each rnott, pop.size 
increase.mutations.plot <- ggplot(cum.strains.avg, aes(x = shifted.time, y = avg.strains)) +
  geom_ribbon(aes(ymin = (avg.strains -2*sd.strains), ymax = (avg.strains + 2*sd.strains)),
              alpha = .5, fill="steelblue2", color="steelblue2") +
  geom_line(color = "black", lwd = 1) +
  facet_grid(pop.size~rnott, scales = "free_y") +
  #theme(panel.margin = unit(1.5, "lines")) + 
  labs(x = "Time", y = "Number of Cumulative Mutations") +
  theme(axis.text.x  = element_text(angle=45, vjust=0.5, size=12)) +
  geom_vline(aes(xintercept = 0), linetype = "dotdash") +
  scale_x_continuous(breaks=seq(0, 150, 50)) +
  
  geom_vline(data = time.max.groups, aes(xintercept = diversity), linetype = 2, color = "red") +
  
  geom_vline(data = time.max.groups, aes(xintercept = entropy), linetype = "twodash", color = "purple")

increase.mutations.plot
save_plot(paste0(fig_path, "increase.mutations.long.range.pdf"), increase.mutations.plot, base_height = 8, base_aspect_ratio = 2)



# Looking at the derivative of increase mutations
cum.strains.master <- adply(.data = data.files, .margins = 1, .id=NULL, .expand = F, function(x) {
  load(as.character(x$file.list))
  trial.params <- get_params(x)
  cum.strains.trial = combine_mutations(time.records = time.records)
  cum.strains.trial = cbind(trial.params, cum.strains.trial)
  return(cum.strains.trial)
})


desired.rnotts <- seq(0.9, 1.9, .2)
desired.rnotts <- c(.9, seq(1,4,1))


time.max.groups <- group_by(time.max.master, rnott, pop.size) %>%
  summarise(infect = mean(max.infect),
            diversity = mean(max.diversity),
            entropy = mean(max.entropy)) %>%
  filter(rnott %in% desired.rnotts)

### Need to take the derivative of all first, and then average 
cum.strains.avg <- group_by(cum.strains.master, rnott, pop.size, vtime) %>%
  filter(rnott %in% desired.rnotts) %>%
  summarise(avg.strains = mean(cum.strains),
            sd.strains = sd(cum.strains)) # standard deviation, not of the mean 

cum.strains.avg %>%
  group_by(rnott, pop.size) %>%
  mutate(derivative = avg.strains - lag(avg.strains, default = avg.strains[1])) -> cum.strains.derivative

mutations.derivative <- ggplot(cum.strains.derivative, aes(x = vtime, y = derivative)) +
  geom_line(color = "black", lwd = 1, size = 1.5) +
  facet_grid(pop.size~rnott, scales = "free_y") +
  labs(x = "Time", y = "Number of Mutations, Prime") +
  theme(axis.text.x  = element_text(angle=45, vjust=0.5, size=12)) +
  #scale_x_continuous(breaks=seq(0, 150, 50)) +
  theme(strip.background = element_blank()) + 
  geom_vline(data = time.max.groups, aes(xintercept = diversity), linetype = 2, color = "red", size = 1) +
  geom_vline(data = time.max.groups, aes(xintercept = infect), linetype = "dotdash", size = 1) +
  geom_vline(data = time.max.groups, aes(xintercept = entropy), linetype = "twodash", color = "purple", size = 1)

save_plot(filename = paste0(fig_path, "mutations.derivative.long.pdf"), mutations.derivative, base_height = 7, base_aspect_ratio = 1.5)

### Look at Derivative of Diversity and of Entropy with confidence intervals 
desired.rnotts <- seq(0.9, 1.9, .2)
desired.rnotts <- c(.9, seq(1,4,1))

time.max.groups <- group_by(time.max.master, rnott, pop.size) %>%
  summarise(infect = mean(max.infect),
            diversity = mean(max.diversity),
            entropy = mean(max.entropy)) %>%
  filter(rnott %in% desired.rnotts)


trajectory.master %>% 
  group_by(rnott, pop.size, vtime) %>%
  filter(rnott %in% desired.rnotts) %>%
  summarise(mean.entropy = mean(entropy),
            sd.entropy = sd(entropy),
            mean.diversity = mean(diversity),
            sd.diversity = sd(diversity)) -> trajectory.metrics


trajectory.metrics %>%
  group_by(rnott, pop.size) %>%
  mutate(derivative.entropy = mean.entropy - lag(mean.entropy, default = mean.entropy[1]),
         derivative.diversity = mean.diversity - lag(mean.diversity, default = mean.diversity[1])) %>%
  gather(Type, Value, derivative.entropy, derivative.diversity)  -> trajectory.derivative


trajectory.metrics <- gather(trajectory.metrics, Type, Value, 4:7) # %>%
trajectory.metrics  %>%  filter(Type == "mean.diversity" | Type == "sd.diversity") %>%
  spread(Type, Value) -> trajectory.metrics


timeseries.diveristy <- ggplot(data = trajectory.metrics, aes(x = vtime, y = mean.diversity)) +
  geom_line() +
  geom_ribbon(aes(ymin = (mean.diversity -2*sd.diversity), ymax = (mean.diversity + 2*sd.diversity)),
              alpha = .5, fill="indianred1", color="indianred1") +
  facet_grid(pop.size~rnott) +
  labs(x = "Time", y = "Entropy") +
  theme(axis.text.x  = element_text(angle=45, vjust=0.5, size=12)) +
  #scale_x_continuous(breaks=seq(0, 150, 50)) +
  theme(strip.background = element_blank()) +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(data = time.max.groups, aes(xintercept = diversity), linetype = 2, color = "red", size = 1) +
  geom_vline(data = time.max.groups, aes(xintercept = infect), linetype = "dotdash", size = 1) +
  geom_vline(data = time.max.groups, aes(xintercept = entropy), linetype = "twodash", color = "purple", size = 1)


timeseries.diveristy

save_plot(paste0(fig_path, "timeseries.diversity.long.pdf"), timeseries.diveristy, base_height = 8, base_aspect_ratio = 1.5)





trajectory.entropy <- filter(.data = trajectory.derivative, Type == "derivative.entropy")
trajectory.diversity <- filter(.data = trajectory.derivative, Type == "derivative.diversity")


#save(trajectory.metrics, file = paste0(data_path, "trajectory.metrics.RData"))
#save(trajectory.derivative, file = paste0(data_path, "trajectory.derivative.RData"))

metrics.entropy <- ggplot(data = trajectory.entropy, aes(x = vtime, y = Value))  + 
  geom_line() +  #, color = Type, group = Type)) +
  facet_grid(pop.size~rnott) +
  labs(x = "Time", y = "Metrics") +
  theme(axis.text.x  = element_text(angle=45, vjust=0.5, size=12)) +
  #scale_x_continuous(breaks=seq(0, 150, 50)) +
  theme(strip.background = element_blank()) +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(data = time.max.groups, aes(xintercept = diversity), linetype = 2, color = "red", size = 1) +
  geom_vline(data = time.max.groups, aes(xintercept = infect), linetype = "dotdash", size = 1) +
  geom_vline(data = time.max.groups, aes(xintercept = entropy), linetype = "twodash", color = "purple", size = 1)

metrics.entropy

save_plot(filename = paste0(fig_path,"trajectory.entropy.derivative.long.pdf"), metrics.entropy, base_height = 8, base_aspect_ratio = 1.5)

save_plot(filename = paste0(fig_path,"trajectory.diversity.derivative.long.pdf"), metrics.diversity, base_height = 8, base_aspect_ratio = 1.5)

### Scatterplot of Time of Max Events (Diversity vs Entropy, Entropy vs Infect, etc )
time.max.master <- adply(.data = data.files, .margins = 1, .id=NULL, .expand = F, function(x) {
  load(as.character(x$file.list))
  trial.params <- get_params(x)
  max.infect <- all_time_max_infected(time.records)
  max.diversity <- all_time_max_diversity(time.records)
  max.entropy <- all_time_max_entropy(time.records)
  max.times.trial <- cbind(trial.params,max.infect, max.diversity, max.entropy)
  return(max.times.trial)
})

max.times.trial.m <- melt(data = time.max.master, id.vars = c("rnott", "pop.size", "max.infect"), 
                          measure.vars = c("max.diversity", "max.entropy"))

max.times.trial.long <- filter(max.times.trial.m, rnott <= 1 | rnott >= 2)
# Having Some Trouble with it getting masked 
save(max.times.trial.m, file = paste0(data_path, "max.times.trial.long.RData"))

mini.times = sample_n(max.times.trial.long, 5000)
time.scatter <- ggplot(data = max.times.trial.long) + 
  geom_jitter(aes(x = max.infect, y = value, color = variable, shape = variable), alpha = .75 )+
  geom_abline(intercept = 0, slope = 1) +
  facet_grid(rnott~pop.size, scales = "free_x") +
  #scale_x_continuous(breaks=seq(0, 150, 50)) +
  scale_y_continuous(breaks=seq(0, 150, 50)) +
  scale_colour_manual(values = c("purple","orange")) +
  labs(x = "Time of Max Infected", y = "Time of Max Metric")


save_plot(paste0(fig_path, "time.scatter.pdf"), time.scatter, base_height = 8, base_aspect_ratio = 1.2)






