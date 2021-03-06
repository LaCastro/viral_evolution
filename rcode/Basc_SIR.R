if(grepl('meyerslab', Sys.info()['login'])) setwd('~/Documents/projects/viral_evolution/viral_evolution_repo/rcode/')
if(grepl('laurencastro', Sys.info()['login'])) setwd('~/Documents/zika_alarm/zika_code/')


#Simulation for SIR 
library(deSolve)
source("disease_models.R")
source("evo_functions.R")


init <- c(S = 1-1e-5, I = 1e-5, 0.0) #1e-5
beta = 0.999
gamma  = 0.14286
times <- seq(1, 80, by = 1)

parameters <- c(beta, gamma)
out.regular <- as.data.frame(ode(y = init, times = times, func = sir, parms = parameters))
out.regular$time <- NULL

matplot(times, out.regular, type = "l", xlab = "Time", ylab = "Susceptibles and Recovereds", main = "SIR Model", lwd = 1, lty = 1, bty = "l", col = 2:4)
legend("topright", c("Susceptibles", "Infecteds", "Recovereds"), pch = 1, col = 2:4)

population.out <- 10^5*out.regular #population x proportion to get numbers 

#For use in the mutation simulation 
infected <- population.out$I


#Parameters for mutation model
seq_len <- 100
alphabet = c(1, 2, 3, 4)
year_mut_rate <- 1.8*10^-3
day_mut_rate <- year_mut_rate / 365 #per site per day

# Number of times to run model 
iterations <- 20


## Creating dataframes to hold analysis 
genetic.diversity.master <- data.frame(matrix("numeric", nrow = length(times), ncol = iterations))
genetic.divergence.master <- data.frame(matrix("numeric", nrow = length(times), ncol = iterations))
number.strains.master <- data.frame(matrix("numeric", nrow =length(times), ncol = iterations))
total.cumulative.strains.master <- data.frame(matrix("numeric", nrow = length(times), ncol = iterations))
number.mutations.master <- data.frame(matrix("numeric", nrow = length(times), ncol = iterations))


#########################
#Initialization of Population 
#initalize diveristy and divergence vectors
tic()

for(d in 1:iterations) {
  print(paste("Beginning Iteration", d, sep = "-"))
  diversity <- vector("numeric", length(times))
  divergence <- vector("numeric", length(times))
  total.strains <- vector("numeric", length(times)) # For each time step number of new strains
  circulating.strains <- vector("numeric", length(times)) 
  number.mutations <- vector("numeric", length(times))
  
  
  #initialize sequence 
  
  base_haplotype = rep(1, seq_len)
  
  #Population is defined for each time step by: 
  #h index (keeping track of the number of total strains)
  # frequency of each strain 
  # strain sequences
  
  population <- init_population(times = times, base_haplotype = base_haplotype)
  


  #Beginning Mutation and Drift 
  for (i in 1:(max(times)-1)) {
    
    #marking the time step 
    timestep <- times[i]
    #print(paste("timestep", timestep, sep = "-"))
    
    #calculate diversity and divergence for the timestep before the mutation and replication event
    total.strains[i] <- length(population$hindex)
    
    
    current.haplotypes <- get_current(population, timestep)  
    #num.current.strains <- length(current.haplotypes$hindex)
    circulating.strains[i] <- length(current.haplotypes$hindex)
    diversity[timestep] <- get_diversity(current.haplotypes)
    divergence[timestep] <- get_divergence(current.haplotypes, base_haplotype) 
    
    
    
    #get expected number of mutations assuming one sequence per individual
    mean.mutation <- day_mut_rate*seq_len*infected[timestep]
    num.mutations <- rpois(1, mean.mutation)
    number.mutations[i] <- num.mutations
    
    
    #If there are mutations....   
    if (num.mutations > 0) {    
      
      #choose equal number of haplotypes as number of mutations from infected by their frequency in the population
      chosen.strains <- as.vector(rmultinom(n = 1, size = num.mutations, prob = current.haplotypes$frequencies))
      original.current.strains <- length(chosen.strains) #how many strains are there in the population to mutate
      
      for (j in 1:original.current.strains) { #for each of the original strains in the population
        num.current.strains <- length(current.haplotypes$hindex) # how many are there-this will be added to 
        num.times.specific.strain <- chosen.strains[j] # how many times did we sample this particular strains
        
        #if we sampled this particular strain at least once 
        if (num.times.specific.strain > 0) {
          for (m in 1:num.times.specific.strain) { #for however many times we are going to mutate strain j
            
            # Mutate the strain 
            new.strain <- mutate.strain(strain = current.haplotypes$strain[[j]]) 
            
            #check if strain is the same as any current strain 
            is.it.current <- check.with.current(new.strain, current.haplotypes) 
            
            
            if(is.it.current$distance > 0) {
              index.same.strain <- is.it.current$index
              current.haplotypes$frequencies[index.same.strain] <- ((current.haplotypes$frequencies[index.same.strain]*infected[timestep])+1)/infected[timestep]
            } else {
              
              #Index the next 
              num.current.strains <- num.current.strains+1
              current.haplotypes$strain[[num.current.strains]] <- new.strain
              current.haplotypes$hindex[num.current.strains] <- "temporary"
              
              # Calculate and change the frequencies
              #calulate.mutated.frequencies(current.haplotypes, infected, og.strain.index, strain.index, timestep)
              new.strain.frequency <- 1/infected[timestep]
              old.strain.frequency <- (current.haplotypes$frequencies[j]*infected[timestep]-1)/infected[timestep] # new frequency of old strain
              if(old.strain.frequency < 0) {
                old.strain.frequency <- 0
                # print(paste("old strain less than 0"))
              }
              current.haplotypes$frequencies[j] <- old.strain.frequency
              
              current.haplotypes$frequencies[num.current.strains] <- new.strain.frequency
              #print(paste("frequency.after.mutation",sum(current.haplotypes$frequencies), sep = "-"))
            }
          }
        }
      } 
    } 
    
    
    
    
    #Drift - Replication Part 
    num.new.infecteds <- ceiling(infected[timestep+1])
    #print(paste("num.infections.in next", num.new.infecteds, sep = "-"))
    new.strains <- rmultinom(n = 1, size = num.new.infecteds, prob = current.haplotypes$frequencies)
    
    #list of new haplotypes
    nex.gen.hap <- list(hindex = current.haplotypes$hindex, frequencies = as.vector(new.strains)/num.new.infecteds, strain = current.haplotypes$strain)
    #print(paste("sum of freq after rep",sum(nex.gen.hap$frequencies), sep = "-"))
    #check if any have strains that were mutated have 0 frequency - i.e. mutated and didn't replicate 
    
    if (num.mutations > 0) { 
      temporary.strains <- which(nex.gen.hap$hindex == "temporary") #get what strains are temporary 
      if(length(temporary.strains) > 0) {
        erase.strains <- c()
        
        for(j in 1:length(temporary.strains)) { #for how many temporary strains there are 
          if (nex.gen.hap$frequencies[temporary.strains[j]] == 0) { # if the don't have a frequency, store them to delete later
            erase.strains <- c(erase.strains, j)
          } 
        }
        
        #deleteing 0 frequencies 
        if(length(erase.strains) > 0) {
          nex.gen.hap$hindex <- nex.gen.hap$hindex[-temporary.strains[erase.strains]]
          nex.gen.hap$frequencies <-nex.gen.hap$frequencies[-temporary.strains[erase.strains]]
          nex.gen.hap$strain <- nex.gen.hap$strain[-temporary.strains[erase.strains]]
        }
      }
    }
    
    #(paste("sum of freq in next gen",sum(nex.gen.hap$frequencies), sep = "-"))
    #update historical record 
    #print(paste("updatinghistorical.record", timestep, sep = "."))
    
    new.hindex <- as.character(nex.gen.hap$hindex)
    pop.hindex <- as.character(population$hindex)
    
    same.hindex <- which(new.hindex%in%pop.hindex) # These indicies of nex.gen.hap that have the same as population and can add them 
    same.indices <- as.numeric(new.hindex[same.hindex]) # The hindex in the population 
    
    #For those that are the same, set their frequencies in next time step
    for(l in 1:length(same.hindex)) {
      population$frequency[timestep+1,same.indices] <- nex.gen.hap$frequencies[same.hindex]
    }
    
    different.hindex <- which(!new.hindex%in%pop.hindex) 
    #print("split differint index")
    
    different.strains <- split.list(nex.gen.hap, different.hindex) #The different nex.gen.hap strain indcices 
    
    #older.hindex <- which(!pop.hindex%in%new.hindex) #which hindex were once present but not in current
    #older.strains <- split.list(population, older.hindex)
    #print("split older index")
    
    
    #If temporary, check if strain matches previous that have gone extinct
    #Loops through all those strains that weren't in the previous generation
    
    #rewriting this - 
    
    if (length(different.hindex) > 0 ) {
      population <- update.ancestor.pop(population, different.strains)
      #print("entered this loop")
      #print(paste("sum of frequencies in next.gen", sum(population$frequency[timestep+1,]), sep  = "-"))
    }
    #print(paste("sum of frequencies in next.gen", sum(population$frequency[timestep+1,]), sep  = "-"))
  }
  
  #Last Time Step 
  genetic.metrics <- calculate_last.time.step(times = times, population,
                                              diversity, divergence, total.strains = total.strains, 
                                              circulating.strains = circulating.strains)
  
  par(mfrow = c(1,1))
  plot_results(out, genetic.metrics) 
  
  
  
  
  number.mutations.master[,d] <- number.mutations
  genetic.divergence.master[,d] <- genetic.metrics$divergence
  genetic.diversity.master[,d] <- genetic.metrics$diversity
  total.cumulative.strains.master[,d] <- genetic.metrics$total.strains
  number.strains.master[,d] <- genetic.metrics$circulating.strains
  par(mfrow = c(1,1))
  matplot(times, population$frequency, type = "l", xlab = "Time", ylab = "Frequency", main = "One Simulation R0 = 7", lwd = 1, lty = 1, bty = "l", col = 2:5)
}


simulation.metrics.master <- list()
simulation.metrics.master[[1]] <- genetic.divergence.master
simulation.metrics.master[[2]] <- genetic.diversity.master
simulation.metrics.master[[3]] <- total.cumulative.strains.master
simulation.metrics.master[[4]] <- number.strains.master
simulation.metrics.master[[5]] <- number.mutations.master
simulation.metrics.master[[6]] <- infected


plot.final(simulation.metrics.master, times, beta = beta, gamma = gamma)

par(mfrow = c(2,2))
divergence.range <- boxplot(genetic.divergence.master)
peak.diversity <- boxplot(colMax(genetic.diversity.master), main = "Peak Diversity")
peak.circulating.strains <- boxplot(colMax(number.strains.master), main = "Number of Circulating Strains at Peak")
total.strains <- boxplot(colMax(total.cumulative.strains.master), main = "Cumulative Novel Strains")

toc()


