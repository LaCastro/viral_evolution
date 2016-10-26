### Script for Analyzing Based on strains records
if(grepl('meyerslab', Sys.info()['login'])) setwd('~/Documents/projects/viral_evolution/viral_evolution_repo/rcode/')
if(grepl('laurencastro', Sys.info()['login'])) setwd('~/Documents/projects/viral_evolution/rcode/')

sapply(c('analyze_saved_sims.R', 'run_mutate_branch.R'), source)

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


file.list <- get_vec_of_files(dir_path = data_path, type = "trial", r0_seq, N)
combos <- sort(apply(expand.grid(r0_seq, N), 1, paste, collapse = "_", sep = "")) 
data.files.trial <- data.frame(cbind(combos, file.list))

epidemic.trials <- set_epidemic_criteria(time.records, threshold.prev = .025, threshold.prop = .25)
data.files.trial$epidemic.trials <- epidemic.trials


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



# function to go into strain.records and put frequency into long form
desired.rnotts <- c(0.9, 1.5, 2.5)
sample <- sample(1:100, 5)


epidemic.trials <- set_epidemic_criteria(data.files.trial, threshold.prev = .025, threshold.prop = .25)
data.files.trial$epidemic.trials <- epidemic.trials

## First get time at max for each trial
# Build in max number of circulating strains to see what the clonal interference chart looks lik 

trajectories <- adply(.data = data.files.trial, .margins = 1,.id=NULL, .expand = F, function(x) {
  load(as.character(x$file.list))
  trial.params <- get_params(x)
  epidemic.index <- x$epidemic.trials[[1]]$epidemic.trials
  if (epidemic.index[1] > 0) {
  time.records.a <- align_time_series_all(time.records = time.records[epidemic.index])
  vI.trajectory <- combine_time_records(time.records.a)
  vI.trajectory <- cbind(trial.params, vI.trajectory)
  vI.trajectory$iter = epidemic.index[vI.trajectory$iter]
  return(vI.trajectory)
  }
})

trajectories %>% filter(iter == 31 & pop.size == 1000) %>%
  group_by(rnott, pop.size) %>%
  select(vtime, shifted.time, diversity, entropy) -> trajectories.sample



time.max.master <- adply(.data = data.files.trial, .margins = 1, .id=NULL, .expand = F, function(x) {
  load(as.character(x$file.list))
  trial.params <- get_params(x)
  #epidemic.index <- x$epidemic.trials[[1]]$epidemic.trials
  max.infect <- all_time_max_infected(time.records)
  max.diversity <- all_time_max_diversity(time.records)
  max.entropy <- all_time_max_entropy(time.records)
  iter <- seq(1:100)
  max.times.trial <- cbind(trial.params,iter, max.infect, max.diversity, max.entropy)
  return(max.times.trial)
})

time.max.master %>% filter(iter == 31 & pop.size == 1000) %>%
  gather(metric, value, 4:6) -> time.max.sample

file.list <- get_vec_of_files(dir_path = data_path, type = "strain", r0_seq, N)
combos <- sort(apply(expand.grid(r0_seq, N), 1, paste, collapse = "_", sep = "")) 
data.files.strain <- data.frame(cbind(combos, file.list))

## Sampling 20/100 strain records
sample = 31

long.strains.master <- adply(.data = data.files.strain, .margins = 1,.id = NULL, .expand = FALSE, function(x) {
  load(as.character(x$file.list))
  trial.params <- get_params(x)
  long.strain.df <- adply(.data = sample, .margins = 1, function(x) {
    long.strains.m <- long.strains(x, strain.records)
    return(long.strains.m)
  })
  long.strain.df <- cbind(trial.params, long.strain.df)
  return(long.strain.df)
})

long.strains.master %>% filter(pop.size == 1000) -> long.strains.sample

# Clonal Interference Graph 
#time.max.sample %>% gather(metric, value, 4:6) -> time.max.sample.long


ci.pop.1000 <- ggplot2::ggplot(data = long.strains.sample, aes(x = time, y = frequency, fill = strain.name)) +
  facet_wrap(~rnott, nrow = 3) +
  geom_area(color = "black", size = .2, alpha = .3) + guides(fill = FALSE) +
  geom_vline(data = time.max.sample, aes(xintercept = jitter(value, .25), color = metric), size = 1.25, alpha = .85) +
  scale_color_manual(values = c("black", "mediumblue", "purple"), guide = FALSE) + 
  #theme(strip.background = element_blank(),  strip.text.x = element_blank()) +
  labs(x = "Time", y = "Frequency")

diversity.plot <- ggplot(trajectories.sample, aes(x = vtime, y = diversity)) + facet_wrap(~rnott, nrow = 3) +
  geom_line(color = "black", size = 2)

entropy.plot <- ggplot(trajectories.sample, aes(x = vtime, y = entropy)) + facet_wrap(~rnott, nrow = 3) +
  geom_line(color = "blue4", size = 2)

combined <- plot_grid(ci.pop.1000, diversity.plot, entropy.plot, nrow = 1)
save_plot(filename = paste0(fig_path, "combined.CI2.pdf"), combined, base_height = 8, base_aspect_ratio = 1.9)
