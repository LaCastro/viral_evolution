
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
  

#calculate the diversity between all strains 
get_diversity <- function(current.haplotypes) {
  num.of.strains <- length(current.haplotypes$frequencies)
  diversity = 0
  if (num.of.strains > 1) {
    for (i in 1:(num.of.strains-1)) {
      for (j in (i+1):num.of.strains) {
        strain_a <- current.haplotypes$strain[[i]]
        strain_b <- current.haplotypes$strain[[j]]
        frequency_a <- current.haplotypes$frequencies[i]
        frequency_b <- current.haplotypes$frequencies[j]
        frequency_ab <- frequency_a* frequency_b
        diversity = diversity + frequency_ab * get_distance(strain_a, strain_b)
      }
    }
  } else {
    diversity = 0
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
    } 
  }
  return(differences/length(strain_a))
}



#calculate the divergence
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
  matplot(times, population.out, type = "l", xlab = "Time", ylab = "Susceptibles and Recovereds", main = "SIR Model", lwd = 1, lty = 1, bty = "l", col = 2:4)
  legend("topright", c("Susceptibles", "Infecteds", "Recovereds"), pch = 1, col = 2:4)
  
  
  plot(times, genetic.metrics$diversity, type = "l", xlab = "Time", ylab = "Diversity", main = 'SIR Model', col = "red")
  plot(times, genetic.metrics$divergence, type = "l", xlab = "Time", ylab = "Divergence", main = "SIR Model", col = "blue")
  plot(times, genetic.metrics$total.strains, type = 'l', xlab = "Time", ylab = "Strain Numbers", main = "Number of Circulating Strains", col = "2")
  legend("topleft", c("Total", "Current"), pch = 1, col = 2:3)
  lines(genetic.metrics$circulating.strains, col = 3)
  
}

check.with.current <- function(new.strain, current.haplotypes) {
  num.strains <- length(current.haplotypes$hindex)
  distance <- 0
  for (i in 1:num.strains) {
    if (distance > 0) { 
      return
      } else { 
      calculated.distance <- get_distance(new.strain, current.haplotypes$strain[[i]])
      if(calculated.distance == 0) {
      distance = distance + 1
      }
    }
  return(list(distance = distance, index = i))
  }
}

add.newstrain.population <- function(population, new.strain, nex.gen.hap,index) {
  new.strain.hindex <- length(population$hindex)+1
  population$hindex[new.strain.hindex] <- new.strain.hindex
  population$strain[[new.strain.hindex]] <- new.strain
  new.strain.name <- paste("strain", new.strain.hindex, sep = ".")
  population$frequency$new.strain.hindex <- rep(0, length(times))
  names(population$frequency)[names(population$frequency) == 'new.strain.hindex'] <- new.strain.name
  population$frequency[timestep+1,new.strain.hindex] <- nex.gen.hap$frequencies[index]    
} 

#function for collecting an index of a list
split.list <- function(population, index) {
  desired.hindex <- population$hindex[index]
  desired.freq <- population$frequencies[index]
  desired.strains <- population$strain[index]
  return(list(hindex = desired.hindex, frequencies = desired.freq, strains = desired.strains))
}


#returns the indices of two sublists for strains that match in each-the indices that are returned are the strain indices in older population
compare.two.list <- function(list1_new, list2_old) {
  strains_1 <- list1_new$strains
  strains_2 <- list2_old$strains
  
  num.unplaced.strains <- rep(0,length(list1_new$hindex))
  
  if (length(strains_2$hindex)==0) {
    return
  } else {
    
    for (m in 1:length(list1_new$hindex)) { #for how many different strains there are 
      temp_strain <- strains_1[[m]]
      for (k in 1:length(list2_old$hindex)) {
        ancestor_strain <- strains_2[[k]]
        distance <- get_distance(temp_strain, ancestor_strain)
        if (distance == 0) num.unplaced.strains[m] <- k
      }
    }
  }
  return(num.unplaced.strains)
}

#function to add frequencies to ancestror strains that have not been in the population but have returned-why am I doing this, wouldn't it give an underestimation 
#Just change it so that if it''s different you add a new strain
# will have a new function at the end to compare which strains were the same, and compare 
#update.ancestor.pop <- function(population, comparing.new.old, timestep, older.hindex, different.hindex, different.strains) {
 # for (i in 1:length(comparing.new.old)) {
  #  if (comparing.new.old[i] == 0) {
      
   #   new.strain.hindex <- length(population$hindex)+1
    #  population$hindex[new.strain.hindex] <- new.strain.hindex
     # population$strain[new.strain.hindex] <- different.strains$strains[i]
      #new.strain.name <- paste("strain", new.strain.hindex, sep = ".")
      #population$frequency$new.strain.hindex <- rep(0, length(times))
      #names(population$frequency)[names(population$frequency) == 'new.strain.hindex'] <- new.strain.name
     # population$frequency[timestep+1,new.strain.hindex] <- different.strains$frequencies[i]    
      
    #} else {
    #  older.hindex <- comparing.new.old[i]
   #   population$frequency[timestep+1, older.hindex] <- nex.gen.hap$frequencies[different.hindex[i]]
  #  }
 # }
#  return(population)
#}

update.ancestor.pop <- function(population, different.strains) {
  for (i in 1:length(different.strains$hindex)) {
      new.strain.hindex <- length(population$hindex)+1
      population$hindex[new.strain.hindex] <- new.strain.hindex
      population$strain[new.strain.hindex] <- different.strains$strains[i]
      new.strain.name <- paste("strain", new.strain.hindex, sep = ".")
      population$frequency$new.strain.hindex <- rep(0, length(times))
      names(population$frequency)[names(population$frequency) == 'new.strain.hindex'] <- new.strain.name
      population$frequency[timestep+1,new.strain.hindex] <- different.strains$frequencies[i]  
  }
  return(population)
}

plot.final <- function(simulation.metrics.master, times, beta, gamma) {
  oma=c(0,0,3,0)
  simulation.divergence <- simulation.metrics.master[[1]]
  simulation.diversity <- simulation.metrics.master[[2]]
  simulation.total.strains <- simulation.metrics.master[[3]]
  simulation.num.strains <- simulation.metrics.master[[4]]
  
  par(mfrow = c(2,2))

  matplot(times, simulation.divergence, type = "l", xlab = "Time", ylab = "Divergence", lwd = 1, lty = 1, bty = "l", 
          col = 1:10, main =  paste("SIR Model: R0", round(beta/gamma, digits = 1), sep = " "))
  matplot(times, simulation.diversity, type = "l", xlab = "Time", ylab = "Diversity", lwd = 1, lty = 1, bty = "l", col = 1:10)
  matplot(times, simulation.total.strains, type = "l", xlab = "Time", ylab = "Total Number of Strains", lwd = 1, lty = 1, bty = "l", col = 1:10)
  matplot(times, simulation.num.strains, type = "l", xlab = "Time", ylab = "Number of Circulating Strains", lwd = 1, lty = 1, bty = "l", col = 1:10)
}

colMax <- function(data) sapply(data, max, na.rm = TRUE)

