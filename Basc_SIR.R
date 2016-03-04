#Simulation for SIR 
init <- c(S = 1-1e-5, I = 1e-5, 0.0) #1e-4
parameters <- c(beta = 0.599, gamma = 0.14286)
times <- seq(1, 70, by = 1)
out <- as.data.frame(ode(y = init, times = times, func = sir, parms = parameters))
out$time <- NULL
population.out <- 10^5*out #population x proportion to get numbers 


#For use in the mutation simulation 
infected <- population.out$I[1:70]
times <- c(1:70)


#########################
#Initialization of Population 
#initalize diveristy and divergence vectors
diversity <- rep(0, length(times))
divergence <- rep(0, length(times))
total.strains <- rep(0, length(times))
circulating.strains <- rep(0, length(times))



#Parameters for mutation model
seq_len <- 100
alphabet = c(1, 2, 3, 4)
year_mut_rate <- 1.8*10^-3
day_mut_rate <- year_mut_rate / 365 #per site per day

#initialize sequence 
base_haplotype = rep(1, seq_len)

population <- init_population(times = times, base_haplotype = base_haplotype)


#Beginning Mutation and Drift 
for (i in 1:(max(times)-1)) {
  
  #marking the time step 
  timestep <- times[i]
  print(paste("timestep", timestep, sep = "-"))
  
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
  print(paste("number of mutations", num.mutations, sep = "-"))
  
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
        for (m in 1:num.times.specific.strain) { #for however many times we did 
          
          # Mutate the strain 
          new.strain <- mutate.strain(strain = current.haplotypes$strain[[j]]) 
          
          #Index the next 
          num.current.strains <- num.current.strains+1
          
          current.haplotypes$strain[[num.current.strains]] <- new.strain
          current.haplotypes$hindex[num.current.strains] <- "temporary"
          
          # Calculate and change the frequencies
          #calulate.mutated.frequencies(current.haplotypes, infected, og.strain.index, strain.index, timestep)
          new.strain.frequency <- 1/infected[timestep]
          old.strain.frequency <- (current.haplotypes$frequencies[j]*infected[timestep]-1)/infected[timestep] # new frequency of old strain
          current.haplotypes$frequencies[j] <- old.strain.frequency
          current.haplotypes$frequencies[num.current.strains] <- new.strain.frequency
          print(paste("frequency.after.mutation",sum(current.haplotypes$frequencies), sep = "-"))
        }
      }
    } 
  } 
  
  
  
  
  #Drift - Replication Part 
  num.new.infecteds <- ceiling(infected[timestep+1])
  print(paste("num.infections.in next", num.new.infecteds, sep = "-"))
  new.strains <- rmultinom(n = 1, size = num.new.infecteds, prob = current.haplotypes$frequencies)
  
  #list of new haplotypes
  nex.gen.hap <- list(hindex = current.haplotypes$hindex, frequencies = as.vector(new.strains)/num.new.infecteds, strain = current.haplotypes$strain)
  print(paste("sum of freq in next gen",sum(nex.gen.hap$frequencies), sep = "-"))
  
  #check if any have strains that were mutated have 0 frequency - i.e. mutated and didn't replicate 
  
  if (num.mutations > 0) { 
    temporary.strains <- which(nex.gen.hap$hindex == "temporary") #get what strains are temporary 
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
    
    
    #update historical record 
    #update_population <- function(nex.gen.hap, population, timestep) {
    
    print(paste("updatinghistorical.record", timestep, sep = "."))
    
    new.hindex <- as.character(nex.gen.hap$hindex)
    pop.hindex <- as.character(population$hindex)
    same.hindex <- which(new.hindex%in%pop.hindex)
    different.hindex <- which(!new.hindex%in%pop.hindex) #which new strain indexes aren't in the population
    older.strains <- which(!pop.hindex%in%new.hindex) #which hindex were once present but not in current
    num.old.strains <- length(older.strains)
    
    #For those that are the same, update their frequencies in next time step
    for(l in 1:length(same.hindex)) {
      population$frequency[timestep+1,same.hindex[l]] <- nex.gen.hap$frequencies[same.hindex[l]]
    }
    
    #If temporary, check if strain matches previous that have gone extinct
    #Loops through all those strains that weren't in the previous generation
    
    #rewriting this - 
    
    if (length(different.hindex) > 0 ) {
      print("there are temporary strains")
      for (m in 1:length(different.hindex)) { #for how many different strains there are 
        temp_strain <- nex.gen.hap$strain[[different.hindex[m]]]
        
        #Loops through all the strains that have been in the population but weren't in the previous
        #If not any, go straight to add new 
        if (num.old.strains > 0 ) {
          matched = 0
          for(k in 1:num.old.strains) {
            ancestor_strain <- population$strain[[older.strains[k]]]      
            #checks the distance between the two, if its 0 (i.e., same frequency) puts that frequency in the timestep 
            if(get_distance(temp_strain, ancestor_strain) == 0) {
              print(paste("temporary matched previous", k, sep = "-"))
              population$frequency[timestep+1,older.strains[k]] <- nex.gen.hap$frequencies[different.hindex[m]]
              matched = 1
              print("match")
            }
          }
          if(matched == 1) return #stop looking through older strains
        } 
        
        
        #Adding to population
        new.strain.hindex <- length(population$hindex)+1
        population$hindex[new.strain.hindex] <- new.strain.hindex
        population$strain[[new.strain.hindex]] <- temp_strain
        new.strain.name <- paste("strain", new.strain.hindex, sep = ".")
        population$frequency$new.strain.hindex <- rep(0, length(times))
        names(population$frequency)[names(population$frequency) == 'new.strain.hindex'] <- new.strain.name
        population$frequency[timestep+1,new.strain.hindex] <- nex.gen.hap$frequencies[different.hindex[m]]    
      } #Loop if there aren't any strains that were once present, but not in the most recent generation
    }
    print(paste("finished.updating.records", timestep, sep  = "-"))
  }

  #Last Time Step 
genetic.metrics <- calculate_last.time.step(times, population,
                                              diversity, divergence, total.strains = total.strains, 
                                              circulating.strains = circulating.strains)

plot_results(out, genetic.metrics) 


