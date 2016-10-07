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


## Reading in files from saved data runs
# Whatever R0s and popsizes we are interested in 
N = c(100, 1000, 10000)
r0_seq = c(seq(0.9, 2, .1), seq(2.5,5, 0.5))



file.list <- get_vec_of_files(dir_path = data_path, type = "trial", r0_seq, N)
combos <- sort(apply(expand.grid(r0_seq, N), 1, paste, collapse = "_", sep = "")) 
data.files <- data.frame(cbind(combos, file.list))

######## Getting Final Times across R0s and popsizes 
final.time.master <-  adply(.data = data.files, .margins = 1,.id=NULL, .expand = F, function(x) {
  load(as.character(x$file.list))
  trial.params <- get_params(x)
  final.time.trial <- data.frame(all_epi_time(time.records)); colnames(final.time.trial) = "days"
  final.time.trial <-  cbind(trial.params, final.time.trial)
  return(final.time.trial)
})

epi.time.plot <- ggplot(final.time.master, aes(factor(rnott), days)) + 
  facet_wrap(facets = ~pop.size, nrow = 1) +
  geom_boxplot()+
  labs(x = expression("R"[0]), y = "Time (days)")
save_plot(paste0(fig_path, "epi.time.r0around1.pdf"), epi.time.plot, base_height = 8, base_aspect_ratio = 1.7)



####### Boxplots of lifespan of wildtype strains across R0s and popsizes 
wildtype.life.master <- adply(.data = data.files, .margins = 1, .id = NULL, .expand = F, function(x) {
  load(as.character(x$file.list))
  trial.params <- get_params(x)
  trial.lifespan <- data.frame(all_time_life_wildtype(strain.records)); colnames(trial.lifespan) = "days"
  trial.lifespan <- cbind(trial.params, trial.lifespan)
  return(trial.lifespan)
})

time.wildtype <- ggplot(wildtype.life.master, aes(factor(rnott), days)) + 
  facet_wrap(~pop.size) +
  geom_boxplot()+
  labs(x = expression("R"[0]), y = "Life Span of Wildtype")
save_plot(paste0(fig_path, "time.wildtype.aroundR0.pdf"), time.wildtype, base_height = 8, base_aspect_ratio = 1.2)



##### Histograms of Lifespans of Mutants
## This takes awhile 
mutant.life.master <- adply(.data = data.files, .margins = 1, .id = NULL, .expand = F, function(x) {
  load(as.character(x$file.list))
  trial.params <- get_params(x)
  mutant.lifespan <- data.frame(all_mutant_lifespan(strain.records))
  colnames(mutant.lifespan) = "lifespan"
  mutant.lifespan <- cbind(trial.params, mutant.lifespan)
  return(mutant.lifespan)
})


## Plotting Code
mutant.lifespan.plot <- ggplot(mutant.life.master, aes(lifespan)) + 
  facet_grid(pop.size~rnott, scales = "free_y") +
  geom_histogram() +
  theme(panel.margin = unit(1.5, "lines")) +
  labs(x = "Lifespan (Days)", y = "Count")
save_plot(paste0(fig_path, "mutant.life.span.aroundr01.pdf"), mutant.lifespan.plot, base_height = 8, base_aspect_ratio = 1.2)

save(mutant.life.master, file = "mutant.life.big.master.RData")
mini.mutant.life = sample_n(mutant.life.master, 10000)


mutant.lifespan.plot <- ggplot(mini.mutant.life, aes(x = rnott, y  = lifespan)) +
  facet_wrap(~pop.size, nrow = 3, scales = "free_y") +
  geom_point() +
  geom_violin(aes(group = rnott)) +
  geom_jitter(alpha = 0.3) +
  labs(x = labs(x =expression("R"[0]), y = "Lifespan(Days)"))+
  theme(strip.background = element_blank())



mutant.lifespan.plot <- ggplot(mini.mutant.life, aes(factor(rnott), lifespan)) +
  facet_wrap(~ pop.size, nrow = 3,scales = "free_y") +
  geom_boxplot() +
 # geom_violin(aes(group = rnott)) +
  #geom_jitter(alpha = 0.3) +
 # labs(x = labs(x =expression("R"[0]), y = "Lifespan(Days)"))+
  theme(strip.background = element_blank())


### Mutant 
mutant.lifespan.plot.density <- ggplot(aes(lifespan, color = factor(rnott), fill = factor(rnott)), data = mini.mutant.life) +  
  geom_density(adjust = 5, alpha = 0.1, position = "stack") +
  facet_wrap(~pop.size, nrow = 1) +
  xlim(0, 50) + 
  theme(strip.background = element_blank())
mutant.lifespan.plot.density
save_plot(paste0(fig_path, "mutant.life.big.span.pdf"), mutant.lifespan.plot.density, base_height = 8, base_aspect_ratio = 1.2)

## Means and Proportions of Mutant Life Spans
mutant.life.master <- group_by(mutant.life.master, rnott, pop.size)
mean.mutant.life <- summarise(mutant.life.master, mean.life = mean(lifespan))

total.entires <- mutant.life.master %>%
                    group_by(rnott, pop.size) %>%
                    summarise(total = n())

dead.end.entries <- mutant.life.master %>%
                      group_by(rnott, pop.size) %>%
                      filter(lifespan == 0) %>%
                      count(rnott, pop.size)
mutant.life.proportion <- inner_join(total.entires, dead.end.entries) %>%
                             mutate(freq.dead = n/total)
                            #inner_join(mean.mutant.life)  
mutant.life.proportion <- inner_join(mutant.life.proportion, mean.mutant.life)

plot.mean.life <- ggplot(mutant.life.proportion, aes(y = mean.life, x = factor(rnott), width = .75)) +
  facet_wrap(~pop.size, nrow = 1) + 
  geom_bar(stat="identity") + 
  labs(x =expression("R"[0]), y = "Days")

plot.dead.end.freq <- ggplot(mutant.life.proportion, aes(y = freq.dead, x = factor(rnott), width = .75)) +
  facet_wrap(~pop.size, nrow = 1) + 
  geom_bar(stat="identity") + 
  labs(x =expression("R"[0]), y = "Dead-End Proportion")





###### Code for producing boxplots of max genetic variation metrics 
## Want to load the data and calculate the max metrics for all 
max.metrics.master <-  adply(.data = data.files, .margins = 1, .id=NULL, .expand = F,function(x) {
  load(as.character(x$file.list))
  trial.params <- get_params(x)
  genetic.metrics <- get_max_genetic_metrics(time.records = time.records, rnott = trial.params$rnott)
  genetic.metrics <-  cbind(rep(trial.params$pop.size, nrow(genetic.metrics)), genetic.metrics)
  colnames(genetic.metrics)[1] <- "pop.size"
  return(genetic.metrics)
})

# plotting code 
max.metrics.m <- melt(data = max.metrics.master, id.vars = c("rnott", "pop.size"),
                      measure.vars = c("max.entropy", "max.diversity"))
colnames(max.metrics.m) <- c("rnott", "pop.size", "type", "value") 

max.metrics.long.range <- filter(max.metrics.m, rnott == c(.9, seq(1,5,.5)))
max.metrics.short.range <- filter(max.metrics.m, rnott == c(.9, seq(1,2,.1)))

max.metric.plot <- ggplot(max.metrics.short.range, aes(factor(rnott), value)) + 
  facet_grid(type~pop.size, scales = "free_y") +
  geom_boxplot()+
  labs(x = expression("R"[0]), y = "Distance")
max.metric.plot
save_plot(paste0(fig_path, "max.metric.entropy.diversity.short.pdf"), max.metric.plot, base_height = 8, base_aspect_ratio = 1.2)


############## Code for getting metrics at infection max 
max.infect.metrics.master <-  adply(.data = data.files, .margins = 1, .id=NULL, .expand = F,function(x) {
  load(as.character(x$file.list))
  trial.params <- get_params(x)
  time.records.a <- align_time_series_all(time.records)
  metrics <- metrics_at_max_infect(time.records = time.records.a)
  metrics <-  cbind(rep(trial.params, nrow(metrics)), metrics)
  return(metrics)
})


desired.rnotts.long <- c(.9, seq(1,4,1))
desired.rnotts.short <- seq(0.9, 1.9, .2)

plot.metric.at.maxinfect(max.infect.metrics.master, desired.rnotts = desired.rnotts.long, desired.metric = "entropy")


######## Code for Combing threshold metrics and plotting 
# This isn't working 
threshold.metrics.master <- adply(.data = data.files, .margins = 1,.id=NULL, .expand = F, function(x) {
  load(as.character(x$file.list))
  trial.params <- get_params(x)
  threshold.metrics = end_threshold_metrics(threshold = .1, trial.params$pop.size, records.list = time.records, rnott = trial.params$rnott)
  threshold.metrics = cbind(rep(trial.params$pop.size, nrow(threshold.metrics)), threshold.metrics)
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
  trial.params <- get_params(x)
  threshold.metrics = beginning_threshold_metrics(threshold = .1,trial.params$pop.size, records.list = time.records, rnott = trial.params$rnott)
  threshold.metrics = cbind(rep(trial.params$pop.size, nrow(threshold.metrics)), threshold.metrics)
  colnames(threshold.metrics)[1] <- "pop.size"
  return(threshold.metrics)
})

beg.threshold.metrics.m <- melt(data = beg.threshold.metrics.master, id.vars = c("rnott", "pop.size"), 
                                measure.vars = c("entropy", "diversity"))
colnames(beg.threshold.metrics.m) <- c("rnott", "pop.size", "type", "value") 

beg.threshold.long.range <- filter(beg.threshold.metrics.m, rnott == c(.9, seq(1,5,.5)))
beg.threshold.short.range <- filter(beg.threshold.metrics.m, rnott == c(.9, seq(1,2,.1)))

threshold.plot <-ggplot(beg.threshold.short.range, aes(factor(rnott), value)) + 
  facet_grid(type~pop.size, scales = "free_y") +
  geom_boxplot()+
  labs(x = expression("R"[0]), y = "Distance")
save_plot(paste0(fig_path, "beg.threshold.entropy.short.pdf"), threshold.plot, base_height = 8, base_aspect_ratio = 1.2)


### beginning threshold metrics for different thresholds
### 
thresholds = seq(0.01, .1, .01)


multiple.beg.thresholds <- adply(.data = thresholds, .margins = 1, .id = NULL, .expand = F, function(x) {
  threshold = x
  beg.metrics <- beg_threshold_metrics(data.files,threshold)
  beg.metrics <- cbind(threshold, beg.metrics)
  return(beg.metrics)
})

edited.m.beg.thres <- filter(multiple.beg.thresholds, rnott == 0.9 | rnott == 1.1 | rnott == 1.2, pop.size == 10000) %>%
                                  group_by(threshold, pop.size, rnott) %>%
                                  summarise(mean.diversity = mean(diversity), 
                                            mean.entropy = mean(entropy))
                                            #mean.time = mean(time))  

m.thresholds.long <- gather(edited.m.beg.thres, "type", "value", 4:5)

beg.range.thresholds.plot <- ggplot(m.thresholds.long, aes(x = threshold, y = value)) + 
  facet_grid(type~rnott, scale = "free_y") +
  geom_line(size = 2) 
save_plot(paste0(fig_path, "beg.thresholds.around1.pdf"), beg.range.thresholds.plot, base_height = 8, base_aspect_ratio = 1.75)


###### Code to compare the average rate at which mutations accumulate depending on R0 
cum.strains.master <- adply(.data = data.files, .margins = 1, .id=NULL, .expand = F, function(x) {
  load(as.character(x$file.list))
  trial.params <- get_params(x)
  cum.strains.trial = combine_mutations(time.records = time.records)
  cum.strains.trial = cbind(trial.params, cum.strains.trial)
  return(cum.strains.trial)
})

## This works! 
cum.strains.groups <- group_by(cum.strains.master, rnott, pop.size, vtime)
  summarise(avg.strains = mean(cum.strains), 
            sd.strains = sd(cum.strains)) %>%
  filter(rnott %in% desired.rnotts) -> cum.strains.avg
  


time.max.master <- adply(.data = data.files, .margins = 1, .id=NULL, .expand = F, function(x) {
  load(as.character(x$file.list))
  trial.params <- get_params(x)
  max.infect <- all_time_max_infected(time.records)
  max.diversity <- all_time_max_diversity(time.records)
  max.entropy <- all_time_max_entropy(time.records)
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
increase.mutations.plot <- ggplot(cum.strains.avg, aes(x = vtime, y = avg.strains)) +
  geom_ribbon(aes(ymin = (avg.strains -2*sd.strains), ymax = (avg.strains + 2*sd.strains)),
              alpha = .5, fill="steelblue2", color="steelblue2") +
  geom_line(color = "black", lwd = 1) +
  facet_grid(pop.size~rnott, scales = "free_y") +
  #theme(panel.margin = unit(1.5, "lines")) + 
  labs(x = "Time", y = "Number of Cumulative Mutations") +
  theme(axis.text.x  = element_text(angle=45, vjust=0.5, size=12)) +
  scale_x_continuous(breaks=seq(0, 150, 50)) +
  geom_vline(data = time.max.groups, aes(xintercept = diversity), linetype = 2, color = "red") +
  geom_vline(data = time.max.groups, aes(xintercept = infect), linetype = "dotdash") +
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

head(time.max.master)

time.max.groups <- group_by(time.max.master, rnott, pop.size) %>%
  summarise(infect = mean(max.infect),
            diversity = mean(max.diversity),
            entropy = mean(max.entropy)) %>%
  filter(rnott %in% desired.rnotts)


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
  scale_x_continuous(breaks=seq(0, 150, 50)) +
  theme(strip.background = element_blank()) + 
  geom_vline(data = time.max.groups, aes(xintercept = diversity), linetype = 2, color = "red", size = 1) +
  geom_vline(data = time.max.groups, aes(xintercept = infect), linetype = "dotdash", size = 1) +
  geom_vline(data = time.max.groups, aes(xintercept = entropy), linetype = "twodash", color = "purple", size = 1)

mutations.derivative
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

head(trajectory.metrics)

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


### Time Plots of Entropy

trajectory.master <- adply(.data = data.files, .margins = 1, .id=NULL, .expand = F, function(x) {
  load(as.character(x$file.list))
  trial.params <- get_params(x)
  trajectories <- combine_time_records(time.records)
  trajectories <- cbind(trial.params, trajectories)
  return(trajectories)
})
save(trajectory.master, file = paste0(data_path,"trajectory.master.RData"))

#Sampling 20 random iterations from 100 to plot 
random.iter <- sample(x = seq(0,100), size = 20)
desired.rnotts <- c(0.9, seq(1, 5, 1))
head(trajectory.master)

trajectory.long <- group_by(trajectory.master, rnott, pop.size, vtime) %>%
  filter(rnott %in% desired.rnotts) %>%
  filter(iter %in% random.iter)  %>%
  summarise(entropy = mean(entropy),
            entropy.sd = sd(entropy)) 

average.vi <- group_by(trajectory.master, rnott, pop.size, vtime) %>%
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


combined.short.trajectories <- plot_grid(stochastic.entropy, average.vi.plot, labels = c("A", "B"))
save_plot(filename = paste0(fig_path, "combined.short.traj.pdf"), plot = combined.short.trajectories, base_height = 8, base_aspect_ratio = 2)



###### 


## Sampling 20/100 strain. records. 
sample <- sample(x = seq(1:100), size = 20)

# function to go into strain.records and put frequency into long form
long.strains <- function(index, strain.records) {
  trial.strain.record <- strain.records[[index]]
  trial.strain.record$time <- seq(1:nrow(trial.strain.record))
  strain.record.long <- gather(trial.strain.record, strain.name, frequency, -iter, -time)
  strain.record.long[is.na(strain.record.long)] <- 0
  return(strain.record.long)
}


long.strains.master <- adply(.data = sample, .margins = 1, function(x) {
  long.strains.m <- long.strains(x, strain.records)
  return(long.strains.m)
})


#### Getting Sample Time for Samples

time.max.samples <- adply(.data = sample, .margins = 1, function(x) {
  time.record <- time.records[[x]]
  max.infect <- time_max_infected(time.record)
  max.diversity < time_max_diversity(time.record)
  max.entropy <- time_max_entropy(time.record)
  max.strains <- time_max_strains(time.record)
  max.times.trial <- cbind(max.infect, max.diversity, max.entropy, max.strains)
  return(max.times.trial)
})

time.max.samples.long <- gather(time.max.samples, type, value, -X1)


rnott.1.5.pop.1000 <- ggplot(data = long.strains.master, aes(x = time, y = frequency, fill = strain.name)) +
  facet_wrap(~ X1, scales = "free") +
  geom_area(color = 'black', size = .2, alpha = .4) + guides(fill = FALSE) +
  geom_vline(data = time.max.samples.long, aes(xintercept = jitter(value, .15), color = type), size = 1.25, alpha = .85) +
  scale_color_manual(values = c("darkorange1", "mediumblue", "purple", "black")) +
  theme(strip.background = element_blank(),  strip.text.x = element_blank()) +
  labs(x = "Time", y = "Frequency")
rnott.1.5.pop.1000

  

time.metrics <- data.frame(rbind(max.entropy, max.diversity, max.infect)); colnames(time.metrics)[1] <- "value"
time.metrics$type <- NA
time.metrics$type <- c("max.entropy", "max.diversity", "max.infect")

max.entropy <- time_max_entropy(x = time.record.1)
max.diversity <- time_max_diversity(x=time.record.1)
max.infect <- time_max_infected(x = time.record.1)


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


head(aligned.time.series.average)
align.ts.long <- gather(aligned.time.series, metric, value, -iter, -vtime, -shifted.time) %>%
  group_by(vtime) %>%
  



ggplot(align.ts.long, aes(x = shifted.time, y = value, group = iter, color = iter))+guides(color = FALSE)+
  facet_wrap(~metric, scales = "free") + geom_line()
