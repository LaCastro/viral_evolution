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


file.list <- get_vec_of_files(dir_path = data_path, type = "strain", r0_seq, N)
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
long.strains <- function(index, strain.records) {
  trial.strain.record <- strain.records[[index]]
  trial.strain.record$time <- seq(1:nrow(trial.strain.record))
  strain.record.long <- gather(trial.strain.record, strain.name, frequency, -iter, -time)
  strain.record.long[is.na(strain.record.long)] <- 0
  return(strain.record.long)
}

## Sampling 20/100 strain records
sample <- sample(x = seq(1:100), size = 20)

long.strains.master <- adply(.data = sample, .margins = 1, function(x) {
  long.strains.m <- long.strains(x, strain.records)
  return(long.strains.m)
})


