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
r0_seq = c(seq(0.9, 2, .2), seq(2,5, 0.5))

#r0_seq = c(seq(0.9, 2, .1))
#r0_seq = c(.9, 1, 1.1)

file.list <- get_vec_of_files(dir_path = data_path, type = "trial", r0_seq, N)
combos <- sort(apply(expand.grid(r0_seq, N), 1, paste, collapse = "_", sep = "")) 
data.files <- data.frame(cbind(combos, file.list))

######## Getting Final Times across R0s and popsizes 
# Get trajectories of infected to see what is occuring when R0 ~ 1
epidemic.trials <- alply(.data = data.files, .margins = 1, function(x) {
  load(as.character(x$file.list))
  trial.params <- get_params(x)
  epidemic.trials <- get_epidemic_index(time.records, trial.N = trial.params$pop.size, threshold.prev = .025,
                                        threshold.prop = .25)
  if (length(epidemic.trials) == 0) epidemic.trials <- 0
  epidemic.trials <- cbind(trial.params, epidemic.trials)
  return(epidemic.trials)
})

data.files$epidemic.trials <- epidemic.trials


trajectories <- adply(.data = data.files, .margins = 1,.id=NULL, .expand = F, function(x) {
  load(as.character(x$file.list))
  trial.params <- get_params(x)
  epidemic.index <- x$epidemic.trials[[1]]$epidemic.trials
  vI.trajectory <- combine_time_records(time.records[epidemic.index])
  vI.trajectory <- cbind(trial.params, vI.trajectory$iter, vI.trajectory$vtime, vI.trajectory$vI)
  return(vI.trajectory)
})
colnames(trajectories) <- c("rnott", "pop.size", "iter", "vtime", "vI")

random.iter <- sample(x = seq(1:1000), size = 50)
trajectories %>% filter(iter %in% random.iter) -> trajectories.sample

ggplot(data = trajectories.sample, aes(x = vtime, y = vI, group = iter, color = iter))+geom_smooth()+guides(color = FALSE) +
  facet_grid(rnott~pop.size, scales = "free")

###### Get trajectories of outbreaks that will be considered epidemcis 
epidemic.trajectories <- adply(.data = data.files, .margins = 1, .id = NULL, .expand = F, function(x) {
  load(as.character(x$file.list))
  trial.params <- get_params(x)
  epidemic.index <- x$epidemic.trials[[1]]$epidemic.trials
  vI.trajectory <- combine_time_records(time.records[epidemic.index])
  vI.trajectory <- cbind(trial.params, vI.trajectory$iter, vI.trajectory$vtime, vI.trajectory$vI)
  return(vI.trajectory)
})




colnames(epidemic.trajectories) <- c("rnott", "pop.size", "iter", "vtime", "vI")
random.iter <- sample(x = seq(1:800), size = 200)
epidemic.trajectories %>% filter(iter %in% random.iter) -> epidemic.trajectories.sample

ggplot(data = epidemic.trajectories.sample, aes(x = vtime, y = vI, group = iter, color = iter))+geom_smooth()+guides(color = FALSE) +
  facet_grid(pop.size~rnott, scales = "free")


##############################
# Get the length of time each simulation takes
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

#Histogram 
epi.time.plot.hist <- ggplot(final.time.master, aes(days)) + 
  facet_grid(rnott~pop.size) +
  geom_histogram() +
  labs(x = "Time (days)", y = "Count")
save_plot(paste0(fig_path, "epi.time.r0around1.pdf"), epi.time.plot, base_height = 8, base_aspect_ratio = 1.7)


# Look at the proportion of the population that was infected 
final.infected.master <-  adply(.data = data.files, .margins = 1,.id=NULL, .expand = F, function(x) {
  load(as.character(x$file.list))
  trial.params <- get_params(x)
 # browser()
  final.infected.trial <- data.frame(all_total_infected(time.records = time.records, trial.N = trial.params$pop.size))
  final.infected.trial <- cbind(trial.params, final.infected.trial)
  return(final.infected.trial)
})
colnames(final.infected.master)[3] <- "prop.infected"

final.infected.master %>% filter(pop.size == 1000) -> final.infected.1000
final.infected.1000 %>% filter(rnott == 0.9) %>%
  summarise(max.infect = max(prop.infected))

ggplot(final.infected.master, aes(factor(rnott), prop.infected)) + geom_jitter(alpha = .5) + 
  geom_violin(alpha = .5) + facet_wrap(~pop.size, nrow =1 ) + geom_hline(yintercept = .13, color = "purple") +
  geom_hline(yintercept = .25, color = "orange")

final.infected.master %>% filter(prop.infected > 0.25) -> epidemics

final.infected.master %>% group_by(rnott, pop.size) %>%
                                 filter(prop.infected < 0.35) %>%
                                 tally() %>% 
                                 mutate(prop.excluded = n/9000) -> prop.excluded

final.infected.master %>% filter(prop.infected >= 0.35) -> final.infected.35

ggplot(final.infected.35, aes(factor(rnott), prop.infected)) + geom_jitter(alpha = .5) + 
  geom_violin(alpha = .5) + facet_wrap(~pop.size, nrow =1 )


# Look at the distribution of maxinum number of people infected in a time
histogram.prevalance <- adply(.data = data.files, .margins = 1, .id= NULL, .expand = F, function(x) {
  load(as.character(x$file.list))
  trial.params <- get_params(x)
  max.prev <- data.frame(all_max_infected(x = time.records))
  histogram.data <- cbind(trial.params,max.prev)
  return(histogram.data)
})
colnames(histogram.prevalance)[3] <- "max.prev"
head(histogram.prevalance)
histogram.prevalance <- mutate(histogram.prevalance, ratio = max.prev/(.01*pop.size))

# Filter to look at area right around 1 
# Testing Thresholds 
deterministic.prop %>% filter(pop.size == 10000) -> deterministic.prop.10000

time.thresholds <- data.frame(small = .025*100, medium = .03*1000, large = .03*10000)

histogram.prevalance %>% filter(rnott < 1.3) -> max.prev.around.1

ggplot(data = histogram.prevalance, aes(max.prev))+geom_histogram(bins = 20) +facet_grid(rnott~pop.size, scales = "free_x") #+
  geom_vline(xintercept = 2.5)+geom_vline(xintercept = 25, color = "blue")+geom_vline(xintercept = 250, color = "green")

total.infected.10000 <- ggplot(final.infected.10000, aes(x = factor(rnott), prop.infected)) +
  geom_violin()+
  geom_jitter(alpha = .5) +
  geom_hline(yintercept = .2,  color = "purple", size = 2)+
  geom_hline(yintercept = .36, color = "blue", size = 2) + 
  labs(x =  expression("R"[0]), y = "Proportion of Population Infected-10000")

total.infected.10000.density <- ggplot(final.infected.10000, aes(prop.infected)) +
  geom_density(adjust = 3) +
  facet_wrap(~rnott) +
  geom_vline(xintercept = .2, color = "purple", size = 1.5) +
  geom_vline(xintercept = .36, color = "blue", size = 1.5) + 
  geom_vline(data = deterministic.prop.10000, aes(xintercept = prop), color = "orange", size = 1.5, linetype = "longdash") + 
  labs(y = "Density", x = "Proportion of Population Infected-10000")

total.infected.10000.density  

total.infected.plot <- ggplot(final.infected.master, aes(prop.infected)) + 
  facet_grid(rnott~pop.size) +
  geom_histogram() +
  labs(x = expression("R"[0]), y = "Proportion of Population Infected")



meta.data.master <- cbind(final.infected.master, final.time.master[,3]); colnames(meta.data.master)[4] <- "days"
meta.data.master <- cbind(meta.data.master, histogram.prevalance$max.prev, histogram.prevalance$ratio)
colnames(meta.data.master) <- c("rnott", "pop.size", "prop.infected", "days", "max.prev", "ratio")
head(meta.data.master)
meta.data.master %>% filter(prop.infected >=  .25) -> meta.data.prop25
meta.data.prop35 %>% filter(max.prev >= .02*pop.size) -> meta.data.two.critera 



# Trying Different Relationships 
ggplot(data = meta.data.master, aes(x = ratio, y = prop.infected)) + geom_jitter(alpha = .5)+facet_grid(rnott~pop.size, scales = "free_x")

ggplot(data = meta.data.prop25, aes(x = days, y = prop.infected)) + geom_jitter(alpha = .5)+facet_grid(rnott~pop.size, scales = "free_x")

ggplot(meta.data.two.critera, aes(x = max.prev, y = prop.infected)) + geom_jitter(alpha = .5) + facet_grid(rnott~pop.size, scales = "free_x")






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

data.files$epidemic.indicies <- epidemic.indices

max.metrics.master <-  adply(.data = data.files, .margins = 1, .id=NULL, .expand = F,function(x) {
  load(as.character(x$file.list))
  trial.params <- get_params(x)
  genetic.metrics <- get_max_genetic_metrics(time.records = time.records, rnott = trial.params$rnott)
  genetic.metrics <-  cbind(rep(trial.params$pop.size, nrow(genetic.metrics)), genetic.metrics)
  colnames(genetic.metrics)[1] <- "pop.size"
  return(genetic.metrics)
})

data.files

max.metrics.master <-  adply(.data = data.files, .margins = 1, .id=NULL, .expand = F,function(x) {
  load(as.character(x$file.list))
  trial.params <- get_params(x)
  epidemic.index <- x$epidemic.trials[[1]]$epidemic.trials
  genetic.metrics <- get_max_genetic_metrics(time.records = time.records[epidemic.index], rnott = trial.params$rnott)
  genetic.metrics <-  cbind(rep(trial.params$pop.size, nrow(genetic.metrics)), genetic.metrics)
  colnames(genetic.metrics)[1] <- "pop.size"
  return(genetic.metrics)
})

head(max.metrics.master)
# plotting code 
max.metrics.m <- melt(data = max.metrics.master, id.vars = c("rnott", "pop.size"),
                      measure.vars = c("max.entropy", "max.diversity"))
colnames(max.metrics.m) <- c("rnott", "pop.size", "type", "value") 
#head()
max.metrics.long.range <- filter(max.metrics.m, rnott == c(.9, seq(1,5,.5)))
max.metrics.short.range <- filter(max.metrics.m, rnott == c(.9, seq(1,2,.1)))

max.metric.plot <- ggplot(max.metrics.m, aes(factor(rnott), value)) + 
  facet_grid(type~pop.size, scales = "free_y") +
  geom_boxplot()+
  labs(x = expression("R"[0]), y = "Distance")
max.metric.plot
save_plot(paste0(fig_path, "max.metric.entropy.diversity.epidemics.pdf"), max.metric.plot, base_height = 8, base_aspect_ratio = 1.2)


############## Code for getting metrics at infection max 


max.infect.metrics.master <-  adply(.data = data.files, .margins = 1, .id=NULL, .expand = F,function(x) {
  load(as.character(x$file.list))
  trial.params <- get_params(x)
  epidemic.index <- x$epidemic.trials[[1]]$epidemic.trials
  time.records.a <- align_time_series_all(time.records[epidemic.index])
  metrics <- metrics_at_max_infect(time.records = time.records.a)
  metrics <-  cbind(rep(trial.params, nrow(metrics)), metrics)
  return(metrics)
})

desired.rnotts <- r0_seq
desired.rnotts.long <- c(.9, seq(1,4,1))
desired.rnotts.short <- seq(0.9, 1.9, .2)

plot.metric.at.maxinfect(max.infect.metrics.master, desired.rnotts = desired.rnotts, desired.metric = "diversity")


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
  epidemic.index <- x$epidemic.trials[[1]]$epidemic.trials
  if (epidemic.index[1] > 0) {
    time.records.a <- align_time_series_all(time.records[epidemic.index])
    cum.strains.trial = combine_mutations(time.records = time.records.a)
    cum.strains.trial = cbind(trial.params, cum.strains.trial)
    return(cum.strains.trial)
  }
})

## This works! 
head(cum.strains.master)

cum.strains.master %>% group_by(rnott, pop.size, shifted.time) %>%
  summarise(avg.strains = mean(cum.strains), 
            sd.strains = sd(cum.strains))  %>%
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
  #scale_x_continuous(breaks=seq(0, 150, 50)) +
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
  epidemic.index <- x$epidemic.trials[[1]]$epidemic.trials
  time.records.a <- align_time_series_all(time.records[epidemic.index])
  trajectories <- combine_time_records(time.records.a)
  trajectories <- cbind(trial.params, trajectories)
  browser()
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
  group_by(vtime) 
  


ggplot(align.ts.long, aes(x = shifted.time, y = value, group = iter, color = iter))+guides(color = FALSE)+
  facet_wrap(~metric, scales = "free") + geom_line()
