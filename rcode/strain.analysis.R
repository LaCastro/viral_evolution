## Scratch Strain Sampling Analysis 
library(grid)
library(gridExtra)
library(cowplot)
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)

## Code Path
if(grepl('meyerslab', Sys.info()['login'])) setwd('~/Documents/projects/viral_evolution/viral_evolution_repo/rcode/')

# Figure Path
if(grepl('meyerslab', Sys.info()['login'])) fig_path = "~/Documents/projects/viral_evolution/viral_evolution_repo/figs/"
if(grepl('laurencastro', Sys.info()['login'])) fig_path <- "~/Documents/projects/viral_evolution/figs/"

# Data Path 
if(grepl('meyerslab', Sys.info()['login'])) data_path <- "~/Documents/projects/viral_evolution/viral_evolution_repo/data/"
if(grepl('laurencastro', Sys.info()['login'])) data_path <- "~/Documents/projects/viral_evolution/data/"

sapply(c('analyze_saved_sims.R', 'samplers.R'), source)


N = c(100, 1000, 10000)
r0_seq = c(seq(0.9, 2, .1), seq(2.5,5, 0.5))

file.list <- get_vec_of_files(dir_path = data_path, type = "strain", r0_seq, N)
combos <- sort(apply(expand.grid(r0_seq, N), 1, paste, collapse = "_", sep = "")) 
data.files <- data.frame(cbind(combos, file.list))

load(as.character(x.trail$file.list))
load(as.character(x.strain$file.list))

## Plot to compare the number of strains sampled under each method 
## Versus actual number of strains circulating 

cases.detected.all <- sample_uniform_all(time.records = time.records, n.samples = 9, low.detection.threshold = .1)
strains.sampled.all <- sampler_strains_all(strain.records = strain.records, cases.sample.all = cases.detected.all)
unique.strains <- get_all_unique_samples(samples.strains.all = strains.sampled.all)


### Turn into function to figure outh how many should be sampled per time step 
## to keep the proportional even with the uniform
cases.sampled.mean <- group_by(cases.sampled.list, iter) %>%
  summarise(length.time = length(detected.time),
            total.detected = sum(detected.cases)) %>%
  summarise(mean.time = mean(length.time),
            mean.detected = mean(total.detected)) %>%
  mutate(detected.per.time = mean.detected/mean.time)



sample.strain.plot <- ggplot2::ggplot(data = unique.strains, aes(x = time.steps, y = unique.strains, group = iter, color = iter)) +
  geom_smooth() +
  guides(color = FALSE) +
  labs(x = "Time", y = "Unique Strains Sampled")
sample.strain.plot

circulating.strains.trial = combine_circulating(time.records = time.records)
time.max.infect = data.frame(all_time_max_infected(time.records)); colnames(time.max.infect) <- "infect.max"
time.max.infect <- summarise(time.max.infect,
                             mean.time = mean(infect.max),
                             sd.time = sd(infect.max))

time.max.diversity <- data.frame(all_time_max_diversity(time.records)); colnames(time.max.diversity) <- "diversity.max"
time.max.diversity <- summarise(time.max.diversity,
                             mean.time = mean(diversity.max),
                             sd.time = sd(diversity.max))

time.max.entropy <- data.frame(all_time_max_entropy(time.records)); colnames(time.max.entropy) <- "entropy.max"
time.max.entropy <- summarise(time.max.entropy,
                                mean.time = mean(entropy.max),
                                sd.time = sd(entropy.max))


time.metrics <- rbind(time.max.diversity, time.max.entropy, time.max.infect)
time.metrics <- cbind((seq(1:3)), time.metrics); colnames(time.metrics)[1] <- "id"

sample.vs.real <- cbind(sample.vs.real, unique.strains)
sample.vs.real<- mutate(sample.vs.real, ratio.uniform = unique.strains.uniform/cir.strains)


strain.capture.plot.uniform <- ggplot(data = sample.vs.real, aes(x = unique.strains.uniform, y = cir.strains)) +
  geom_jitter() + 
  geom_abline(intercept = 0, slope = 1)+
  labs(x = "Number of Sampled Strains-Uniform", y = "Number of Circulating Strains") +
  scale_x_continuous(breaks = 1:7) +
  scale_y_continuous(breaks = 1:15)
strain.capture.plot.uniform

strain.capture.scatter.comparison <- ggplot(data = sample.vs.real.long, aes(x =  cir.strains, y = num.strains))+
  geom_jitter()+
  geom_abline(intercept = 0, slope = 1) + 
  facet_wrap(~ratio.method) +
  labs(x = "Number of Sampled Strains", y = "Number of Circulating Strains")
strain.capture.scatter.comparison
save_plot(paste0(fig_path, "capture.method.jitter.pdf"), strain.capture.scatter.comparison, base_height = 7, base_aspect_ratio = 1.5)


strain.capture.hist.comparison <- ggplot(data = sample.vs.real.long, aes(value)) +
  geom_histogram(bins = 20) +
  facet_wrap(~ratio.method)
  strain.capture.hist.uniform
  strain.capture.hist.comparison



ratio.time.uniform<- ggplot(data = sample.vs.real, aes(x = time.steps, y = ratio.uniform, group = iter, color = iter))+
  geom_smooth() + guides(color = FALSE) + 
  ylim(0, 1.5) +
  geom_vline(data = time.metrics, aes(xintercept = mean.time, color = id))



sample.vs.real %>% group_by(iter) %>%
  gather(sample.method, num.strains, unique.strains.prop, unique.strains.uniform) %>%
  gather(ratio.method, value, ratio, ratio.uniform) -> sample.vs.real.long
View(sample.vs.real.long)

levels(sample.vs.real.long$ratio.method) 

ratio.comparison <- ggplot(sample.vs.real.long, aes(x=time.steps, y = value, group = iter, color = iter)) +
  facet_wrap(~ratio.method, nrow = 2) + geom_smooth()+
  guides(color = FALSE)+
  ylim(0, 1.5) +
  labs(x = "Time", y = "Ratio of Sampled to Circulating Strains") +
  theme(strip.background = element_blank()) 


save_plot(paste0(fig_path, "ratio.comparison.pdf"), ratio.comparison, base_height = 8, base_aspect_ratio = 1.1)
