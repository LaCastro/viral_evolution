get_current <- function(population,timestep) {
  num.strains.pop <- length(population$hindex)
  
  #get the population of present strains and their frequencies 
  present.sequences <- list()
  present.frequences <- c()
  present.hindex <- c()
  
  total.present <- 0
  for (m in 1:num.strains.pop) {
    strain <- population$strain[[m]]
    frequency <- population$frequency[timestep,m] 
    hindex <- population$hindex[m]
    
    #Is this sequence present in the present timestep? If it is, record it 
    if (frequency > 0) {
      total.present <- total.present+ 1
      present.sequences[[total.present]] <- strain
      present.frequences[total.present] <- frequency
      present.hindex[total.present] <- hindex
    } else {   
    }
  } 
  current.haplotypes <- list(hindex = present.hindex,frequencies =  present.frequences, strain = present.sequences)
  return(current.haplotypes)
}



#Mutate Strain and Recalculate frequencies 
mutate.strain <- function(strain)  {
  strain <- current.haplotypes$strain[[j]] #Get selected strain from list of current haplotypes
  baseindex <- round(runif(1, min = 1, max = seq_len)) # select an index to mutate
  nucleotide <- strain[baseindex]
  possible.mutations <- alphabet[which(alphabet != nucleotide)]
  strain[baseindex] <- sample(possible.mutations, 1) #Mutate the strain to any of the other options
  return(strain)
}
  

#calculate frequencies -can't get this to work yet 
#calulate.mutated.frequencies <- function(current.haplotypes, og.strain.index, infected, strain.index, timestep) {
  
#  new.strain.frequency <- 1/infected[timestep]
#  old.strain.frequency <- (current.haplotypes$frequencies[og.strain.index]*infected[timestep]-1)/infected[timestep] # new frequency of old strain
#  current.haplotypes$frequencies[og.strain.index] <- old.strain.frequency
#  current.haplotypes$frequencies[strain.index] <- new.strain.frequency




#calculate the diversity between all strains - #working 
get_diversity <- function(current.haplotypes) {
  num.of.strains <- length(current.haplotypes$frequencies)
  diversity = 0
  for (i in 1:num.of.strains) {
    for (j in 1:num.of.strains) {
      strain_a <- current.haplotypes$strain[[i]]
      strain_b <- current.haplotypes$strain[[j]]
      frequency_a <- current.haplotypes$frequencies[i]
      frequency_b <- current.haplotypes$frequencies[j]
      frequency_ab <- frequency_a* frequency_b
      diversity = diversity + frequency_ab * get_distance(strain_a, strain_b)
    }
  }
  return(diversity)
}



#calculate the distance between two particular strains 
#This is working 
get_distance <- function(strain_a, strain_b) {
  differences = 0
    for(i in 1:length(strain_a)) {
      if((strain_a[i]-strain_b[i]) != 0) {
        differences = differences + 1
      } else {
        }
    }
  return(differences/length(strain_a))
}


#calculate the divergence - working! 
get_divergence <- function(current.haplotypes, base_haplotype) {
  divergence = 0
  for (i in 1:length(current.haplotypes$frequencies)) {
    strain_evolved <- current.haplotypes$strain[[i]]
    frequency = current.haplotypes$frequencies[i]
    divergence = divergence + frequency*get_distance(base_haplotype, strain_evolved)
  }
  return(divergence)
}



#initialize population
init_population <- function(times, base_haplotype) {
  #Setting up the population list to store haplotypes through time 
  frequency.df <- data.frame(rep(0, length(times))) #storing the frequencies of each strain
  colnames(frequency.df) <- paste("strain",1, sep = ".")
  
  #storing number of unique strains in the population 
  hindex <- c()
  hindex[1] <- 1 #initializing the first 
  
  #List of strains 
  strain <- list()
  strain[[1]] <- base_haplotype
  
  population <- list(hindex = hindex, frequency = frequency.df, strain = strain )
  population$frequency$strain.1[1] <- 1 #When there is only one strain introduced   
  return(population)
}

calculate_last.time.step <- function(times, population, diversity, divergence, total.strains, circulating.strains) {
  last.time.step <- length(times)
  total.strains[last.time.step] <- length(population$hindex)
  current.haplotypes <- get_current(population, last.time.step)  
  circulating.strains[last.time.step] <- length(current.haplotypes$hindex)
  diversity[last.time.step] <- get_diversity(current.haplotypes)
  divergence[last.time.step] <- get_divergence(current.haplotypes, base_haplotype)
  return(list(diversity = diversity, divergence = divergence, total.strains = total.strains, circulating.strains = circulating.strains))
}

plot_results <- function(out, genetic.metrics) {
  par(mfrow = c(2,2))
  matplot(times, out, type = "l", xlab = "Time", ylab = "Susceptibles and Recovereds", main = "SIR Model", lwd = 1, lty = 1, bty = "l", col = 2:4)
  legend(40, 0.7, c("Susceptibles", "Infecteds", "Recovereds"), pch = 1, col = 2:4)
  
  
  plot(times, genetic.metrics$diversity, type = "l", xlab = "Time", ylab = "Diversity", main = 'SIR Model', col = "red")
  plot(times, genetic.metrics$divergence, type = "l", xlab = "Time", ylab = "Divergence", main = "SIR Model", col = "blue")
  plot(times, genetic.metrics$total.strains, type = 'l', xlab = "Time", ylab = "Strain Numbers", main = "Number of Circulating Strains", col = "2")
  legend("topleft", c("Total", "Current"), pch = 1, col = 2:3)
  lines(genetic.metrics$circulating.strains, col = 3)
  
}