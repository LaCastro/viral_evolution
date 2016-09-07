##################################################################################
##################################################################################
# An R script to perform a stochastic epidemic simulation using an
# Agent Based Model and homogeneous mixing in the population, while mutation is occuring 
###################################################################################
rm(list=ls())
if(grepl('meyerslab', Sys.info()['login'])) setwd('~/Documents/projects/viral_evolution/viral_evolution_repo/rcode/')
if(grepl('laurencastro', Sys.info()['login'])) setwd('~/Documents/projects/viral_evolution/rcode/')

sapply(c('sir_agent_func.R','evo_functions.R', 'sir_mutation_func.R'), source)

library(deSolve)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(reshape2)
library(plyr)

##################################################################################
# Set up initial conditions and parameters 
##################################################################################
set.seed(578194)
nrealisations = 25

epi_mut_params <- function(N = 1000,
                           I_0 = 10,
                           S_0 = N-I_0,
                           R_0 = 0,
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



##################################################################################
# Set up meta-analysis frames to keep track off 
##################################################################################
epi_runs <- list()
epi_size_final = numeric(0) # vector with the epidemic final size estimates from the simulations

# Creating dataframes to hold mutation data
genetic.diversity.master <- data.frame(matrix("numeric", nrow = length(times), ncol = iterations))
genetic.divergence.master <- data.frame(matrix("numeric", nrow = length(times), ncol = iterations))
number.strains.master <- data.frame(matrix("numeric", nrow =length(times), ncol = iterations))
total.cumulative.strains.master <- data.frame(matrix("numeric", nrow = length(times), ncol = iterations))
number.mutations.master <- data.frame(matrix("numeric", nrow = length(times), ncol = iterations))

params <- epi_mut_params(N = 10000, I_0 = 10, delta_t = 1)
trial1 <- sir_mutation_agent(params)
trial1$time_record

strain.freq <- trial1$strain.freq


for (iter in 1:nrealisations) {
  myagent = SIR_agent(N = N,I_0 = I_0,S_0 = S_0,gamma = gamma,R0,tbeg,tend,delta_t)
  #myagent = SIR_agent_detection(N = N, I_0=I_0, S_0=S_0, gamma = gamma, R0, tbeg, tend, delta_t) 
  #myagent = sir_mutation_agent(N = N, I_0=I_0, S_0=S_0, gamma = gamma, R0, tbeg, tend, delta_t) 
 # myagent <- sir_mutation_agent(N=N, I_0, S_0, gamma, R0 = R0 ,tbeg =tbeg , tend = tend, delta_t = delta_t , seq_len = seq_len, alphabet = alphabet, mut_rate = mut_rate)
  I_sim = myagent$I/N
  epi_runs[[iter]] <- data.frame(trial = rep(iter, length(I_sim)), time = myagent$time, I = I_sim) #, D = myagent$D, C = myagent$C)
  #lines(myagent$time,myagent$I/N,lwd=2,col=(iter+1),lty=3)
  
  #cat(iter,nrealisations,"The final size of the epidemic from the agent based stochastic model is ",myagent$final_size,"\n")
  
  epi_size_final = append(epi_size_final, myagent$final_size) 
}

epi_size_final <- data.frame(epi_size_final)
runs.master.df <- do.call("rbind", epi_runs) 

stochastic <- ggplot(runs.master.df, aes(x = time, y = I, color = trial,  group = trial)) +   geom_line() + guides(color = FALSE) 
stochastic

combined <- ggplot() + 
  geom_line(data = runs.master.df, aes(x = time, y = I, color = trial, group = trial)) + 
  guides(color = FALSE) + 
  geom_line(data = sirmodel, aes(x = time, y = I_N), linetype = 5, color = "red", size = 1.5)  + 
  labs(title = paste("Influenza pandemic in population of", N, sep = " "),
       x = "Time, in days", y = "Fraction infected (prevalence)") 


## Need to figure out how to put legend right
#legend("topright",legend=c("Deterministic","Agent Based simulation"),col=c(1,2),lwd=3,lty=c(1,3),bty="n")

final.size <- ggplot(data = epi_size_final, aes(epi_size_final)) + geom_histogram(binwidth = .025) + 
  #geom_vline(xintercept = max(sirmodel$R/N), size = 1.5, colour="red", linetype = "longdash") +
  labs(x = "Distribution of epidemic final size", main = "", y = "Count")
final.size
plot_grid(combined, final.size, ncol = 1)
vfinal


# Detected vs Cumulative

plot_prevalences <- function(df){
  
  df$disc_prob <- paste0(calculate.discover(df$disc_prob), "%")
  
  ggplot(df, aes(detected, median, color=as.factor(r_not), linetype=as.factor(disc_prob), fill=as.factor(r_not))) + 
    geom_line(size=1)+
    geom_ribbon(aes(ymax=max, ymin=min), alpha=0.1, color=NA)+
    scale_y_log10(expand=c(0.01,0.01), limits=c(1,100), breaks = c(5,10,25,50,100))+
    scale_x_continuous(expand=c(0.01,0.01), limits=c(0,30))+
    scale_color_brewer(palette="Set1", direction = -1)+
    scale_fill_brewer(palette="Set1", direction=-1)+
    guides(linetype=FALSE)+
    labs(x = "Cumulative Reported Cases", 
         y = "Prevalence (log scale)", 
         color = expression("R"[0]), 
         fill = expression("R"[0]))
}

runs.master.m  <- melt(data = runs.master.df, id.vars = c( "trial","time"), measure.vars = c("D", "C"))


sampling <- ggplot(runs.master.m, aes(x = time, y = value, color = as.factor(variable), linetype = as.factor(variable), fill = as.factor(variable))) +
  geom_line(size = 1, alpha = 0.2) +
  stat_summary(fun.y = "mean", color = "black", size = 1, geom = "line") + 
  guides(linetype = FALSE) +
  labs(x = "Time", y = "Number of Cases", color = "Type of Case")
sampling
