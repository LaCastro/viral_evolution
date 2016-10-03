## Scratch Script for  working out functions
## Exploratory Analysis Figures 

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
library(tidyr)


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

r0_seq = c(.9, seq(1,4,.5))

r0_seq = seq(.9,1.2, 0.05)


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

threshold.plot


### beginning threshold metrics for different thresholds
### 
thresholds = seq(0.01, .05, .01)


multiple.beg.thresholds <- adply(.data = thresholds, .margins = 1, .id = NULL, .expand = F, function(x) {
  threshold = x
  beg.metrics <- beg_threshold_metrics(data.files,threshold)
  beg.metrics <- cbind(threshold, beg.metrics)
  return(beg.metrics)
})

head(multiple.beg.thresholds)
edited.m.beg.thres <- filter(multiple.beg.thresholds, rnott == 0.9 | rnott == 1.1 | rnott == 1.2, pop.size == 10000) %>%
                                  group_by(threshold, pop.size, rnott) %>%
                                  summarise(mean.diversity = mean(diversity), 
                                            mean.entropy = mean(entropy))
                                            #mean.time = mean(time))  

#### 
m.thresholds.long <- gather(edited.m.beg.thres, "type", "value", 4:5)

beg.range.thresholds.plot <- ggplot(m.thresholds.long, aes(x = threshold, y = value)) + 
  facet_grid(type~rnott, scale = "free_y") +
  geom_line(size = 2) 
save_plot(paste0(fig_path, "increase.mutations.aroundr01.pdf"), increase.mutations, base_height = 8, base_aspect_ratio = 1.75)



########## A function but can also stand on its own 
beg_threshold_metrics <- function(data.files, threshold) { 
    beg.threshold.master <-  adply(.data = data.files, .margins = 1, .id=NULL, .expand = F, function(x) {
      load(as.character(x$file.list))
      trial.params <- get_params(x)
      threshold.metrics = beginning_threshold_metrics(threshold = threshold,trial.params$pop.size, records.list = time.records, rnott = trial.params$rnott)
      threshold.metrics = cbind(rep(trial.params$pop.size, nrow(threshold.metrics)), threshold.metrics)
      colnames(threshold.metrics)[1] <- "pop.size"
      return(threshold.metrics)
    })
    return(beg.threshold.master)
}
  



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
#cum.strains.groups <- group_by(cum.strains.master, vtime)
cum.strains.avg <- summarise(cum.strains.groups, avg.strains = mean(cum.strains), sd.strains = sd(cum.strains)) # standard deviation, not of the mean 


time.max.master <- adply(.data = data.files, .margins = 1, .id=NULL, .expand = F, function(x) {
  load(as.character(x$file.list))
  trial.params <- get_params(x)
  max.infect <- all_time_max_infected(time.records)
  max.diversity <- all_time_max_diversity(time.records)
  max.entropy <- all_time_max_entropy(time.records)
  max.times.trial <- cbind(trial.params,max.infect, max.diversity, max.entropy)
  return(max.times.trial)
})


time.max.groups <- group_by(time.max.master, rnott, pop.size) %>%
                    summarise(infect = mean(max.infect),
                              diversity = mean(max.diversity),
                              entropy = mean(max.entropy)) %>%
                    #gather(type, day, -rnott, -pop.size) %>%
                    filter(rnott <= 2)


#cum.strains.avg$rnott <- factor(cum.strains.avg$rnott)
cum.strains.avg$pop.size <- factor(cum.strains.avg$pop.size)

#time.max.avg$rnott <- factor(time.max.avg$rnott)
#time.max.avg$pop.size <- factor(time.max.avg$pop.size)

# plot the average increase for each rnott, pop.size 
increase.mutations.long.range <- filter(cum.strains.avg, rnott == c(.9, seq(1,5,.5)))
increase.mutations.short.range <- filter(cum.strains.avg, rnott <= 2)

#time.max.avg$rnott = as.numeric(time.max.avg$rnott)
#time.max.long.range <- filter(time.max.long, rnott == c(.9, seq(1,5,.5)))

str(time.max.avg)

increase.mutations.plot <- ggplot(increase.mutations.short.range, aes(x = vtime, y = avg.strains)) +
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
save_plot(paste0(fig_path, "increase.mutations.short.range.pdf"), increase.mutations, base_height = 8, base_aspect_ratio = 1.75)


######### Code to  Get time of max infected and compare against time  max diversity 
### Example will be a low R0  (5) and a high R0 (20) for N = 1000
## Want a two row figure for max diversity and max divergence 

#low.r0.path <- file.list[5]
#high.r0.path <- file.list[20]

#load(low.r0.path)
#times.low.r0 <- all_time_max_infected(time.records)
#diversity.low.r0 <- all_time_max_diversity(time.records)
#diverge.low.r0 <- all_time_max_diverge(time.records)

#low.r0.data <- data.frame(cbind(times.low.r0, diversity.low.r0, diverge.low.r0))
#low.r0.data <- cbind(rep("r0 = 1.5", nrow(low.r0.data)), low.r0.data)
#colnames(low.r0.data) <- c("type", "infected", "diversity", "divergence")

#load(high.r0.path)
#times.high.r0 <- all_time_max_infected(time.records)
#diversity.high.r0 <- all_time_max_diversity(time.records)
#diverge.high.r0 <- all_time_max_diverge(time.records)
#high.r0.data <- data.frame(cbind(times.high.r0, diversity.high.r0, diverge.high.r0))
#high.r0.data <- cbind(rep("r0 = 4", nrow(high.r0.data)), high.r0.data)
#colnames(high.r0.data) <- c("type", "infected", "diversity", "divergence")

#time.data.master <- rbind(low.r0.data, high.r0.data)
#time.data.m <- melt(data = time.data.master, id.vars = c("type", "infected"),
#                    measure.vars = c("diversity", "divergence"))


#combine.low <- combine_time_records(time.records)

#stochastic <- ggplot(combined.time.records.m, aes(x = vtime, y = value, color = iter,  group = iter)) +  
#  geom_line() + guides(color = FALSE) +
#  facet_grid(type~variable) + 
#  theme_cowplot() %+replace% theme(strip.background=element_rect(fill=NULL,color="black", size=1, linetype=1)) +
#  labs(x = "Time", y = "Distance")

#cir.strains.plot <- ggplot(combine.low, aes(x = vtime, y = cir.strains, color = iter,  group = iter)) +  
#  geom_point() + guides(color = FALSE) +
#  #facet_grid(type~variable) + 
#  theme_cowplot() %+replace% theme(strip.background=element_rect(fill=NULL,color="black", size=1, linetype=1)) +
#  labs(x = "Time", y = "Number of STrains")


#time.max.low <- all_time_max_infected(time.records)
#time.max.strains.low <- all_time_max_numstrains(time.records)


#strains.time <- data.frame(cbind(time.max.low, time.max.strains.low))
#strains.time  <- cbind(rep("r0 = 1.5", nrow(strains.time)), strains.time)
#colnames(strains.time) <- c("type", "infected", "num.strains")



#load(high.r0.path)
#times.high.r0 <- all_time_max_infected(time.records)
#times.high.strains <- all_time_max_numstrains(time.records)

#strains.high <- data.frame(cbind(times.high.r0, times.high.strains))
#strains.high <- cbind(rep("r0 = 4", nrow(strains.high)), strains.high)
#colnames(strains.high) <- c("type", "infected", "num.strains")

#strains.data.master <- rbind(strains.time, strains.high)

#time.scatter.strains <- ggplot(data = strains.data.master, aes(x = infected, y = num.strains)) + 
#  geom_point() +
#  geom_abline(intercept = 0, slope = 1) +
#  facet_wrap(facets = "type", nrow = 1) +
#  theme_cowplot() %+replace% theme(strip.background=element_rect(fill=NULL,color="black", size=1, linetype=1)) +
#  labs(x = "Time of Max Infected", y = "Time of Max Circulating Strains")
#save_plot(paste0(fig_path, "strains.time.scatter.pdf"), time.scatter.strains, base_height = 4, base_aspect_ratio = 1.8)



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

save(max.times.trial.m, file = paste0(data_path, "max.times.trial.long.RData"))

mini.times = sample_n(max.times.trial.long, 5000)

time.scatter <- ggplot(data = max.times.trial.long) + 
  geom_jitter(aes(x = max.infect, y = value, color = variable, shape = variable), alpha = .75 )+
 #geom_jitter() + 
  #geom_point(aes(x = max.infect, y = max.entropy), color = "orange") + 
  geom_abline(intercept = 0, slope = 1) +
  facet_grid(rnott~pop.size, scales = "free_x") +
  #scale_x_continuous(breaks=seq(0, 150, 50)) +
  scale_y_continuous(breaks=seq(0, 150, 50)) +
  scale_colour_manual(values = c("purple","orange")) +
#  theme_cowplot() %+replace% theme(strip.background=element_rect(fill=NULL,color="black", size=1, linetype=1)) +
 labs(x = "Time of Max Infected", y = "Time of Max Metric")

time.scatter
#save_plot(paste0(fig_path, "time.scatter.pdf"), time.scatter, base_height = 8, base_aspect_ratio = 1.2)


## Combining 
#strain.records.all <- strain_freq_all(trial)

### Line Plots 
#combined.high <- combine_time_records(time.records) 
#combined.high <- cbind(rep("r0 = 4", nrow(combined.high)), combined.high)
#colnames(combined.high)[1] <- "type"

#combined.low <- combine_time_records(time.records)
#combined.low <- cbind(rep("r0 = 1.5", nrow(combined.low)), combined.low)
#colnames(combined.low)[1] <- "type"

#combined.time.records <- rbind(combined.low, combined.high)
#combined.time.records.m <- melt(data = combined.time.records, id.vars = c("type", "iter", "vtime"),
#                                measure.vars = c("diverge", "diversity"))


#stochastic <- ggplot(combined.time.records.m, aes(x = vtime, y = value, color = iter,  group = iter)) +  
#  geom_line() + guides(color = FALSE) +
#  facet_grid(type~variable) + 
#  theme_cowplot() %+replace% theme(strip.background=element_rect(fill=NULL,color="black", size=1, linetype=1)) +
#  labs(x = "Time", y = "Distance")

#save_plot(filename = paste0(fig_path, "stochastic.line.pdf"), plot = stochastic, base_height = 4, base_width = 7)


#stochastic.diverge <- ggplot(runs.master.df, aes(x = vtime, y = diverge, color = trial,  group = trial)) +  
#  geom_line() + guides(color = FALSE) 
#stochastic.diverge

#number.circulating.strains <- ggplot(runs.master.df, aes(x = vtime, y = cir.strains, color = trial, group = trial)) + 
#  geom_line() + guides(color = FALSE)
#number.circulating.strains

#stochastic.diversity <- ggplot(runs.master.df, aes(x = vtime, y = diversity, color = trial, group = trial)) + 
#  geom_line() + guides(color = FALSE)
#stochastic.diversity

stochastic.infectives <- ggplot(time.records.combined, aes(x=vtime, y = vI, color = iter, group = iter)) +
  geom_line() + guides(color = FALSE)
stochastic.infectives


stochastic.diverge <- ggplot(time.records.combined, aes(x=vtime, y = diverge, color = iter, group = iter)) +
  geom_line() + guides(color = FALSE)


stochastic.entropy <- ggplot(time.records.combined, aes(x=vtime, y = entropy, color = iter, group = iter)) +
  geom_line() + guides(color = FALSE)


time.records <- time_records_all(trial.1000)
time.records.combined <- combine_time_records(time.records)

time.records.combined %>%
  filter(iter == 3) %>%
  ggplot(aes(x=vtime, y = vI)) + geom_line() 

time.records.combined %>%
  filter(iter == 3) %>%
  ggplot(aes(x=vtime, y = entropy)) + geom_line() 

